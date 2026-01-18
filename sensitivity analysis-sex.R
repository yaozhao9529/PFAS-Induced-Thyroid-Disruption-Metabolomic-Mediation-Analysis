library(tidyverse)
library(lmerTest)
library(broom.mixed)
library(openxlsx)


df_main <- read_csv("thyroid-PFAS.csv")
df_covariate <- read_csv("sample-info.csv")

df_all <- df_main %>%
  left_join(df_covariate, by = "ID")

pollutant_vars <- c("C9_HFPO", "PFMOAA", "PFO5DoDA", "PFOA",
                    "PFNA", "PFDA", "PFUnDA", "PFDoDA",
                    "PFHxS", "PFHpS", "PFOS", "PFESA62")

thyroid_vars   <- c("TSH", "FT3", "FT4", "T3", "T4", "T4_T3")

covariates <- c("age", "BMI", "UIC",
                "education", "employment", "gravidity", "parity",
                "income", "smoking", "alcohol", "coffee", "fetal_sex")

ln_vars <- c(thyroid_vars, pollutant_vars)

df_all <- df_all %>%
  mutate(across(all_of(ln_vars), ~ log(.x), .names = "{.col}_log")) %>%
  mutate(
    period      = factor(period, levels = c(1, 2), labels = c("T1", "T3")),
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

run_lme_by_sex_adjusted <- function(dat, sex_value,
                                    pollutant_vars, thyroid_vars, covariates,
                                    conf_method = "Wald") {
  
  dat_sex <- dat %>% filter(fetal_sex == sex_value)
  
  covars_adj <- setdiff(covariates, "fetal_sex")
  
  fit_one <- function(outcome, pollutant) {
    outcome_col  <- paste0(outcome, "_log")
    exposure_col <- paste0(pollutant, "_log")
    
    formula_adjusted <- as.formula(
      paste0(outcome_col, " ~ ", exposure_col, " + period + ",
             paste(covars_adj, collapse = " + "),
             " + (1|ID)")
    )
    
    out <- tryCatch({
      m2 <- lmer(formula_adjusted, data = dat_sex, REML = TRUE)
      
      r2 <- broom.mixed::tidy(m2, conf.int = TRUE, conf.method = conf_method) %>%
        filter(effect == "fixed") %>%
        select(term, estimate, std.error, statistic, p.value, conf.low, conf.high) %>%
        mutate(model = "Adjusted")
      
      r2 %>%
        mutate(
          outcome  = outcome,
          exposure = pollutant,
          sex      = as.character(sex_value),
          exposure_term = exposure_col
        )
    }, error = function(e) {
      tibble(
        term = NA_character_, estimate = NA_real_, std.error = NA_real_,
        statistic = NA_real_, p.value = NA_real_, conf.low = NA_real_, conf.high = NA_real_,
        model = "Adjusted",
        outcome = outcome, exposure = pollutant, sex = as.character(sex_value),
        exposure_term = exposure_col,
        error_msg = conditionMessage(e)
      )
    })
    
    out
  }
  
  res <- tidyr::expand_grid(outcome = thyroid_vars, pollutant = pollutant_vars) %>%
    mutate(result = purrr::map2(outcome, pollutant, fit_one)) %>%
    select(result) %>%
    tidyr::unnest(result)
  
  res_exposure <- res %>%
    group_by(sex) %>% 
    mutate(p.value.FDR = p.adjust(p.value, method = "BH")) %>%
    ungroup()
  
  list(
    exposure_only = res_exposure,
    all_fixed = res
  )
}


res_sex1 <- run_lme_by_sex_adjusted(df_all, sex_value = 1,
                                    pollutant_vars = pollutant_vars,
                                    thyroid_vars = thyroid_vars,
                                    covariates = covariates)

res_sex2 <- run_lme_by_sex_adjusted(df_all, sex_value = 2,
                                    pollutant_vars = pollutant_vars,
                                    thyroid_vars = thyroid_vars,
                                    covariates = covariates)


final_sex_strat <- bind_rows(res_sex1$exposure_only, res_sex2$exposure_only) %>%
  mutate(
    exposure = factor(exposure, levels = pollutant_vars),
    outcome  = factor(outcome, levels = thyroid_vars),
    sex      = factor(sex, levels = c("1", "2"))
  ) %>%
  arrange(outcome, exposure, sex) %>%
  mutate(
    pct_change = (exp(estimate) - 1) * 100,
    pct_low    = (exp(conf.low) - 1) * 100,
    pct_high   = (exp(conf.high) - 1) * 100
  )


iqr_by_sex_ln <- df_all %>%
  select(fetal_sex, all_of(pollutant_vars)) %>%
  pivot_longer(
    cols = all_of(pollutant_vars),
    names_to = "exposure",
    values_to = "x_raw"
  ) %>%
  mutate(
    x_ln = log(x_raw),
    x_ln = ifelse(is.finite(x_ln), x_ln, NA_real_)
  ) %>%
  group_by(fetal_sex, exposure) %>%
  summarise(
    n = sum(!is.na(x_ln)),
    q1_ln = as.numeric(quantile(x_ln, 0.25, na.rm = TRUE, type = 7)),
    q3_ln = as.numeric(quantile(x_ln, 0.75, na.rm = TRUE, type = 7)),
    iqr_ln = q3_ln - q1_ln,
    .groups = "drop"
  ) %>%
  mutate(
    fetal_sex = factor(fetal_sex, levels = c(1, 2)),
    exposure  = factor(exposure, levels = pollutant_vars)
  ) %>%
  arrange(exposure, fetal_sex)


final_sex_strat_IQRln <- final_sex_strat %>%
  left_join(
    iqr_by_sex_ln %>% transmute(sex = as.character(fetal_sex), exposure, iqr_ln, q1_ln, q3_ln, n),
    by = c("sex", "exposure")
  ) %>%
  mutate(
    estimate_IQR    = estimate * iqr_ln,
    conf.low_IQR    = conf.low * iqr_ln,
    conf.high_IQR   = conf.high * iqr_ln,
    pct_change_IQR  = (exp(estimate_IQR) - 1) * 100,
    pct_low_IQR     = (exp(conf.low_IQR) - 1) * 100,
    pct_high_IQR    = (exp(conf.high_IQR) - 1) * 100
  )


write.xlsx(final_sex_strat_IQRln, "model_results_sex_stratified.xlsx")


