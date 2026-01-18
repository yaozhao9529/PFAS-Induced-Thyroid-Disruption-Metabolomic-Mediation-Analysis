library(dplyr)
library(tidyr)
library(stats)

thyroid_PFAS <- read.csv("thyroid-PFAS.csv")

df_covariate <- read.csv("sample-info.csv")

df_all <- thyroid_PFAS %>%
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
  mutate(across(all_of(ln_vars), ~ log(.x), .names = "{.col}_log"))

df_all <- df_all %>%
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

df_all <- cbind(1:nrow(df_all), df_all)
colnames(df_all)[1] <- c("sum-ID")

df_met <- read.csv("Significant_Metabolites.csv")
colnames(df_met)[1:4] <- c("sum-ID", "period", "group", "Sample-ID")
metabolites_ID <- colnames(df_met)[5:2226]

data_metabolites <- cbind(df_met[1],df_met[metabolites_ID])

df_all <- df_all %>%
  inner_join(data_metabolites, by = "sum-ID")

data_T1 <- df_all %>% filter(period == "T1")
data_T3 <- df_all %>% filter(period == "T3")

perform_metabolite_outcome_glm <- function(data, period_label, outcomes, metabolites, covariates) {
  result_list <- list()
  
  for (outcome in outcomes) {
    for (metabolite in metabolites) {

      formula_str <- paste(paste0(outcome, "_log"), "~", metabolite, "+", paste(covariates, collapse = " + "))
      model <- try(glm(as.formula(formula_str), data = data, family = gaussian()), silent = TRUE)
      
      if (!inherits(model, "try-error")) {
        coef_summary <- summary(model)$coefficients
        
        if (metabolite %in% rownames(coef_summary)) {
          estimate <- coef_summary[metabolite, "Estimate"]
          p_val <- coef_summary[metabolite, "Pr(>|t|)"]
          
          conf <- try(confint(model, parm = metabolite, level = 0.95), silent = TRUE)
          if (!inherits(conf, "try-error")) {
            ci_lower <- conf[1]
            ci_upper <- conf[2]
          } else {
            ci_lower <- NA
            ci_upper <- NA
          }
          
          result_list[[paste(outcome, metabolite, sep = "_")]] <- data.frame(
            outcome = outcome,
            metabolite = metabolite,
            estimate = estimate,
            CI_lower = ci_lower,
            CI_upper = ci_upper,
            p = p_val,
            period = period_label
          )
        }
      }
    }
  }
  
  result_df <- do.call(rbind, result_list)
  result_df$FDR <- p.adjust(result_df$p, method = "BH")
  return(result_df)
}

res_T1 <- perform_metabolite_outcome_glm(data_T1, "T1", thyroid_vars, metabolites_ID, covariates)
res_T3 <- perform_metabolite_outcome_glm(data_T3, "T3", thyroid_vars, metabolites_ID, covariates)

all_results_df <- rbind(res_T1, res_T3)

sig_T1 <- res_T1 %>% filter(FDR < 0.2)
sig_T3 <- res_T3 %>% filter(FDR < 0.2)

sig_metabolites_both <- inner_join(sig_T1, sig_T3, by = c("outcome", "metabolite"))

write.csv(sig_metabolites_both, "Key_metabolites_with_outcomes.csv", row.names = FALSE)


all_results_df$sig_flag <- ifelse(all_results_df$FDR > 0.2, 0, 1)


all_results_df$key <- paste(all_results_df$outcome, all_results_df$metabolite, sep = "_")
sig_metabolites_both$key <- paste(sig_metabolites_both$outcome, sig_metabolites_both$metabolite, sep = "_")

all_results_df$shared_flag <- ifelse(
  all_results_df$key %in% sig_metabolites_both$key,
  2,
  all_results_df$sig_flag
)

all_results_df$key <- NULL
sig_metabolites_both$key <- NULL

write.csv(all_results_df, "MWAS_results_with_flag_outcomes.csv", row.names = FALSE)



library(tidyverse)  
library(glmnet)
library(dplyr)
library(ggplot2)
library(caret)
library(glmmLasso)

thyroid_outcome <- thyroid_PFAS[,3:8]
thyroid_outcome <- cbind(1:nrow(thyroid_outcome), thyroid_outcome) 
colnames(thyroid_outcome)[1] <- c("sum-ID")

thyroid_outcome$"lnTSH" <- log(thyroid_outcome$"TSH")
thyroid_outcome$"lnFT3" <- log(thyroid_outcome$"FT3")
thyroid_outcome$"lnFT4" <- log(thyroid_outcome$"FT4")
thyroid_outcome$"lnT3" <- log(thyroid_outcome$"T3")
thyroid_outcome$"lnT4" <- log(thyroid_outcome$"T4")
thyroid_outcome$"lnT4_T3" <- log(thyroid_outcome$"T4_T3")

df_lasso <- df_met %>%
  inner_join(thyroid_outcome, by = "sum-ID")

colnames(df_lasso)[1:4] <- c("sum-ID", "period", "group", "Sample-ID")
df_lasso$group <- as.factor(df_lasso$group)


find_best_lambda <- function(outcome_var, df_lasso, metabolites_ID, batch_size = 500, K = 10) {
  
  y <- df_lasso[[outcome_var]]
  
  lambda_values <- exp(seq(log(0.01), log(100), length.out = 20))
  
  folds <- createFolds(y, k = K, list = TRUE)
  
  VIP_features <- c()
  
  for (start in seq(1, length(metabolites_ID), by = batch_size)) {
    batch_vars <- metabolites_ID[start:min(start + batch_size - 1, length(metabolites_ID))]
    formula <- reformulate(batch_vars, response = outcome_var)
    
    lasso_model <- glmmLasso(
      fix = formula,
      rnd = list(group = ~1),
      lambda = 10,
      data = df_lasso,
      family = gaussian()
    )
    
    coef_df <- as.data.frame(lasso_model$coefficients[-1], stringsAsFactors = FALSE)
    colnames(coef_df) <- c("value")
    coef_df$feature <- names(lasso_model$coefficients[-1])

    selected_batch_features <- coef_df %>%
      filter(value != 0) %>%
      arrange(desc(abs(value))) %>%
      pull(feature)
    
    VIP_features <- unique(c(VIP_features, selected_batch_features))
  }
  
  VIP_features <- head(VIP_features, 500)
  
  cv_results <- sapply(lambda_values, function(lambda) {
    mse_list <- numeric(K)
    
    for (i in seq_len(K)) {
      test_idx <- folds[[i]]
      train_data <- df_lasso[-test_idx, ]
      test_data <- df_lasso[test_idx, ]
      
      final_formula <- reformulate(VIP_features, response = outcome_var)
      
      final_model <- glmmLasso(
        fix = final_formula,
        rnd = list(group = ~1),
        lambda = lambda,
        data = train_data,
        family = gaussian()
      )

      y_pred <- predict(final_model, newdata = test_data)
      mse_list[i] <- mean((test_data[[outcome_var]] - y_pred)^2)
    }
    
    return(mean(mse_list))
  })
  
  best_lambda <- lambda_values[which.min(cv_results)]
  
  print(paste(outcome_var, "best lambda:", best_lambda))
  return(best_lambda)
}

metabolites_ID <- colnames(df_lasso)[5:2226]

# TSH
best_lambda_TSH <- find_best_lambda("lnTSH", df_lasso, metabolites_ID)

# FT3
best_lambda_FT3 <- find_best_lambda("lnFT3", df_lasso, metabolites_ID)

# FT4
best_lambda_FT4 <- find_best_lambda("lnFT4", df_lasso, metabolites_ID)

# T3
best_lambda_T3 <- find_best_lambda("lnT3", df_lasso, metabolites_ID)

# T4
best_lambda_T4 <- find_best_lambda("lnT4", df_lasso, metabolites_ID)

# T4/T3
best_lambda_T4_T3 <- find_best_lambda("lnT4_T3", df_lasso, metabolites_ID)


perform_glmmLasso <- function(outcome_var, best_lambda, df_lasso, metabolites_ID, batch_size = 500) {
  library(glmmLasso)
  library(dplyr)
  
  coef_df <- data.frame()
  
  for (i in seq(1, length(metabolites_ID), by = batch_size)) {
    predictors <- paste(metabolites_ID[i:min(i+batch_size-1, length(metabolites_ID))], collapse = " + ")
    formula <- as.formula(paste(outcome_var, "~", predictors))
    
    lasso_model <- glmmLasso(
      fix = formula,
      rnd = list("group" = ~1),
      lambda = best_lambda,
      data = df_lasso,
      family = gaussian()
    )
    
    coef_values <- lasso_model$coefficients[-1]
    coef_df <- rbind(coef_df, data.frame(feature = names(coef_values), value = coef_values))
  }
  
  coef_final <- coef_df %>%
    group_by(feature) %>%
    slice_max(order_by = abs(value), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    filter(value != 0)
  
  return(coef_final)
}

important_features_TSH <- perform_glmmLasso("lnTSH", 8.85866790410083, df_lasso, metabolites_ID)
important_features_FT3 <- perform_glmmLasso("lnFT3", 0.01, df_lasso, metabolites_ID)
important_features_FT4 <- perform_glmmLasso("lnFT4", 23.3572146909013, df_lasso, metabolites_ID)
important_features_T3 <- perform_glmmLasso("lnT3", 0.0263665089873036, df_lasso, metabolites_ID)
important_features_T4 <- perform_glmmLasso("lnT4", 10, df_lasso, metabolites_ID)
important_features_T4_T3 <- perform_glmmLasso("lnT4_T3", 37.9269019073225, df_lasso, metabolites_ID)


feature_list <- list(
  important_features_TSH = important_features_TSH,
  important_features_FT3 = important_features_FT3,
  important_features_FT4 = important_features_FT4,
  important_features_T3 = important_features_T3,
  important_features_T4 = important_features_T4,
  important_features_T4_T3 = important_features_T4_T3
)

important_features_df <- bind_rows(feature_list, .id = "related_outcome")
important_features_df$related_outcome <- gsub("important_features_", "", important_features_df$related_outcome)

important_features_summary <- important_features_df %>%
  group_by(feature) %>%                     
  summarise(
    value = value[which.max(abs(value))],
    counts = n()  
  ) %>%
  ungroup()

write.csv(important_features_df, "important_features_6 thyroid outcomes.csv")
write.csv(important_features_summary, "important_features_summary.csv")



