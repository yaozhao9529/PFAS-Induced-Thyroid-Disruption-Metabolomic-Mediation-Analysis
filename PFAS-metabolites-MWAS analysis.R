library(dplyr)
library(tidyr)
library(stats)


thyroid_PFAS <- read.csv("thyroid-PFAS.csv")
df_covariate <- read.csv("sample-info.csv")

df_all <- thyroid_PFAS %>%
  left_join(df_covariate, by = "ID")

colSums(is.na(df_all))


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


perform_mwas <- function(data, period_label, pollutants, medium_met_id, covariates) {
  result_list <- list()
  
  for (pollutant in pollutants) {
    
    pollutant <- paste0(pollutant, "_log")
    
    for (metabolite in medium_met_id) {
      
      formula_str <- paste(metabolite, "~", pollutant, "+", paste(covariates, collapse = " + "))
      model <- try(glm(as.formula(formula_str), data = data), silent = TRUE)
      

      if (!inherits(model, "try-error")) {
        
        coef_summary <- summary(model)$coefficients
        
        if (pollutant %in% rownames(coef_summary)) {
          estimate <- coef_summary[pollutant, "Estimate"]
          p_val <- coef_summary[pollutant, "Pr(>|t|)"]
          
          conf <- try(confint(model, parm = pollutant, level = 0.95), silent = TRUE)
          if (!inherits(conf, "try-error")) {
            ci_lower <- conf[1]
            ci_upper <- conf[2]
          } else {
            ci_lower <- NA
            ci_upper <- NA
          }
        
          result_list[[paste(pollutant, metabolite, sep = "_")]] <- data.frame(
            pollutant = pollutant,
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

res_T1 <- perform_mwas(data_T1, "T1", pollutant_vars, metabolites_ID, covariates)
res_T3 <- perform_mwas(data_T3, "T3", pollutant_vars, metabolites_ID, covariates)

all_results_df <- rbind(res_T1, res_T3)

sig_T1 <- res_T1 %>% filter(FDR < 0.05)
sig_T3 <- res_T3 %>% filter(FDR < 0.05)

sig_metabolites_both <- inner_join(sig_T1, sig_T3, by = c("pollutant", "metabolite"))

write.csv(sig_metabolites_both, "Key_metabolites_with_PFAS.csv", row.names = FALSE)

all_results_df$sig_flag <- ifelse(all_results_df$FDR > 0.05, 0, 1)

all_results_df$key <- paste(all_results_df$pollutant, all_results_df$metabolite, sep = "_")
sig_metabolites_both$key <- paste(sig_metabolites_both$pollutant, sig_metabolites_both$metabolite, sep = "_")

all_results_df$shared_flag <- ifelse(
  all_results_df$key %in% sig_metabolites_both$key,
  2,
  all_results_df$sig_flag 
)

all_results_df$key <- NULL
sig_metabolites_both$key <- NULL


write.csv(all_results_df, "MWAS_results_with_flag_PFAS.csv", row.names = FALSE)





