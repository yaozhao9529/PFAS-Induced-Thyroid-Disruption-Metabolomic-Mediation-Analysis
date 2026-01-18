## 0) Packages ----
pkgs <- c("readr","dplyr","tidyr","stringr","lme4","broom.mixed","truncnorm","purrr")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install) > 0) install.packages(to_install)
lapply(pkgs, library, character.only = TRUE)

library(dplyr)
library(readr)
library(lme4)
library(broom.mixed)
library(truncnorm)


## 1) Read & merge ----
pfas <- read_csv("thyroid-PFAS-imputed.csv")
covs <- read_csv("sample_info.csv")

df_all <- pfas %>%
  left_join(covs, by = "ID") %>%
  mutate(
    ID          = as.factor(ID),
    period      = as.factor(period),
    education   = factor(education),
    employment  = factor(employment),
    gravidity   = factor(gravidity),
    parity      = factor(parity),
    income      = factor(income),
    smoking     = factor(smoking),
    alcohol     = factor(alcohol),
    coffee      = factor(coffee),
    fetal_sex   = factor(fetal_sex)
  )

## 2) PFAS list + LOQ ----
pollutant_vars <- c("C9_HFPO","PFMOAA","PFO5DoDA","PFOA","PFNA","PFDA",
                    "PFUnDA","PFDoDA","PFHxS","PFHpS","PFOS","PFESA62")

loq <- c(C9_HFPO = 0.05,
         PFMOAA = 0.01,
         PFO5DoDA = 0.02,
         PFOA = 0.02,
         PFNA = 0.02,
         PFDA = 0.02,
         PFUnDA = 0.02,
         PFDoDA = 0.02,
         PFHxS = 0.05,
         PFHpS = 0.02,
         PFOS = 0.02,
         PFESA62 = 0.05)

## 3) Specify outcomes + covariates
thyroid_vars <- c("TSH","FT3","FT4","T3","T4","T4_T3")

covariates <- c("age", "BMI", "UIC", 
                "education", "employment", "gravidity", "parity", 
                "income", "smoking", "alcohol", "coffee", "fetal_sex")
covars_imp <- c("period", covariates)

## 4) Clean PFAS: "< LOQ" -> NA
clean_pfas_numeric <- function(df, pollutant_vars, loq) {
  df2 <- df
  for (p in pollutant_vars) {
    x <- df2[[p]]
    x <- ifelse(as.character(x) == "< LOQ", NA, x)
    x <- suppressWarnings(as.numeric(x))
    df2[[p]] <- x
  }
  df2
}
df_all_num <- clean_pfas_numeric(df_all, pollutant_vars, loq)

## 5) MI function: lognormal + truncated normal draw
mi_leftcensor_lognormal <- function(df, y, loq_value, covars, m = 20, seed = 2026) {
  set.seed(seed)
  
  needed_cols <- unique(c(y, covars))
  d <- df[, needed_cols, drop = FALSE]
  
  y_vec <- d[[y]]
  cens_idx <- which(is.na(y_vec))
  det_idx  <- which(!is.na(y_vec) & y_vec >= loq_value)
  
  det_cc <- det_idx[complete.cases(d[det_idx, , drop = FALSE])]
  
  if (length(det_cc) < 30) {
    stop(paste0("[", y, "] detected complete cases too few (", length(det_cc),
                "). Consider reducing covariates in covars_imp."))
  }
  
  fml <- as.formula(paste0("log(", y, ") ~ ", paste(covars, collapse = " + ")))
  fit <- lm(fml, data = d[det_cc, , drop = FALSE])
  sigma <- summary(fit)$sigma
  
  mu_hat <- predict(fit, newdata = d)
  
  out_list <- vector("list", m)
  for (k in 1:m) {
    dfk <- df
    if (length(cens_idx) > 0) {
      usable_cens <- cens_idx[complete.cases(d[cens_idx, covars, drop = FALSE])]
      if (length(usable_cens) > 0) {
        log_draw <- truncnorm::rtruncnorm(
          n = length(usable_cens),
          a = -Inf,
          b = log(loq_value),
          mean = mu_hat[usable_cens],
          sd = sigma
        )
        dfk[[y]][usable_cens] <- exp(log_draw)
      }
    }
    out_list[[k]] <- dfk
  }
  out_list
}

## 6) Build m imputed datasets with all PFAS imputed ----
make_imputed_datasets <- function(df, pollutant_vars, loq, covars_imp, m = 20, seed = 2026) {
  imp_list <- mi_leftcensor_lognormal(df, pollutant_vars[1], loq[[pollutant_vars[1]]],
                                      covars = covars_imp, m = m, seed = seed)
  if (length(pollutant_vars) > 1) {
    for (p in pollutant_vars[-1]) {
      tmp_list <- mi_leftcensor_lognormal(df, p, loq[[p]],
                                          covars = covars_imp, m = m, seed = seed + 13)
      for (k in 1:m) {
        imp_list[[k]][[p]] <- tmp_list[[k]][[p]]
      }
    }
  }
  imp_list
}

## 7) Add log columns after imputation ----
add_log_cols <- function(df, thyroid_vars, pollutant_vars) {
  df2 <- df
  for (t in thyroid_vars) {
    df2[[paste0(t, "_log")]] <- log(df2[[t]])
  }
  for (p in pollutant_vars) {
    df2[[paste0(p, "_log")]] <- log(df2[[p]])
  }
  df2
}

## 8) Fit lmer for one dataset and extract exposure term ----
run_models_one_dataset <- function(df, thyroid_vars, pollutant_vars, covariates) {
  res_all <- list()
  
  for (thyroid in thyroid_vars) {
    outcome_col <- paste0(thyroid, "_log")
    
    for (pollutant in pollutant_vars) {
      exposure_col <- paste0(pollutant, "_log")
      
      fml <- as.formula(
        paste0(outcome_col, " ~ ", exposure_col,
               " + period + ",
               paste(covariates, collapse = " + "),
               " + (1|ID)")
      )
      
      fit <- lme4::lmer(fml, data = df, REML = FALSE)
      
      tab <- broom.mixed::tidy(fit, conf.int = FALSE) %>%
        dplyr::filter(effect == "fixed", term == exposure_col) %>%
        dplyr::mutate(
          thyroid = thyroid,
          pollutant = pollutant
        )
      
      model_df <- df %>%
        dplyr::select(ID, period, all_of(c(outcome_col, exposure_col, covariates))) %>%
        dplyr::filter(complete.cases(.))
      n_obs <- nrow(model_df)
      n_id  <- dplyr::n_distinct(model_df$ID)
      
      tab <- tab %>%
        dplyr::mutate(n_obs = n_obs, n_id = n_id)
      
      res_all[[length(res_all) + 1]] <- tab
    }
  }
  
  dplyr::bind_rows(res_all)
}

## 9) Rubin pooling with df, 95% CI, p-value ----
pool_rubin_full <- function(est, se) {
  m <- length(est)
  qbar <- mean(est)
  W <- mean(se^2)
  B <- stats::var(est)
  
  if (is.na(B) || B == 0) B <- 1e-12
  
  Tvar <- W + (1 + 1/m) * B
  se_pool <- sqrt(Tvar)
  
  # Barnard-Rubin type df approximation
  r <- ((1 + 1/m) * B) / W
  df <- (m - 1) * (1 + 1/r)^2
  if (is.na(df) || df < 1) df <- m - 1
  
  tval <- qbar / se_pool
  pval <- 2 * stats::pt(abs(tval), df = df, lower.tail = FALSE)
  
  crit <- stats::qt(0.975, df = df)
  ci_low  <- qbar - crit * se_pool
  ci_high <- qbar + crit * se_pool
  
  data.frame(
    estimate = qbar,
    std.error = se_pool,
    statistic = tval,
    df = df,
    p.value = pval,
    conf.low = ci_low,
    conf.high = ci_high
  )
}

## 10) Run MI + models + pool ----
m <- 20
imp_list <- make_imputed_datasets(df_all_num, pollutant_vars, loq, covars_imp, m = m, seed = 2026)

res_list <- vector("list", m)
for (k in 1:m) {
  dfk <- add_log_cols(imp_list[[k]], thyroid_vars, pollutant_vars)
  rk <- run_models_one_dataset(dfk, thyroid_vars, pollutant_vars, covariates) %>%
    mutate(imputation = k)
  res_list[[k]] <- rk
}
res_m <- bind_rows(res_list)

# pool
pooled <- res_m %>%
  group_by(thyroid, pollutant) %>%
  summarise(
    pooled = list(pool_rubin_full(estimate, std.error)),
    n_obs = first(n_obs),
    n_id  = first(n_id),
    .groups = "drop"
  ) %>%
  tidyr::unnest(pooled) %>%

  mutate(FDR_BH = p.adjust(p.value, method = "BH")) %>%
  group_by(thyroid) %>%
  mutate(FDR_within_thyroid = p.adjust(p.value, method = "BH")) %>%
  ungroup()

## 11) Format as Table Sx ----
Table_Sx <- pooled %>%
  transmute(
    Thyroid = thyroid,
    PFAS = pollutant,
    Beta = estimate,
    SE = std.error,
    CI95_L = conf.low,
    CI95_U = conf.high,
    t = statistic,
    df = df,
    P = p.value,
    FDR = FDR_BH,
    FDR_within = FDR_within_thyroid,
    N_obs = n_obs,
    N_ID = n_id
  ) %>%
  arrange(Thyroid, FDR)

Table_Sx <- pooled %>%
  mutate(
    PFAS_order = factor(pollutant, levels = pollutant_vars),
    Thyroid_order = factor(thyroid, levels = thyroid_vars)
  ) %>%
  transmute(
    Thyroid = thyroid,
    PFAS = pollutant,
    Beta = estimate,
    SE = std.error,
    CI95_L = conf.low,
    CI95_U = conf.high,
    t = statistic,
    df = df,
    P = p.value,
    FDR = FDR_BH,
    FDR_within = FDR_within_thyroid,
    N_obs = n_obs,
    N_ID = n_id,
    PFAS_order = PFAS_order,
    Thyroid_order = Thyroid_order
  ) %>%
  arrange(Thyroid_order, PFAS_order) %>%
  select(-PFAS_order, -Thyroid_order)


## 12) Write out ----
out_path <- "model-results-imputed.csv"
readr::write_csv(Table_Sx, out_path)
message("Saved: ", out_path)



