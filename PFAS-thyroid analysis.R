library(tidyverse)     
library(lmerTest)       
library(broom.mixed)
library(dplyr)
library(openxlsx)          
library(psych)


df_main <- read_csv("thyroid-PFAS.csv")
df_covariate <- read_csv("sample_info.csv")

df_all <- df_main %>%
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


all_results <- list()

counter <- 1

for (thyroid in thyroid_vars) {
  
  outcome_col <- paste0(thyroid, "_log")
  
  for (pollutant in pollutant_vars) {
    
    exposure_col <- paste0(pollutant, "_log")
    
    # 1) crude
    formula_crude <- as.formula(
      paste0(outcome_col, " ~ ", exposure_col, " + period + (1|ID)")
    )
    
    model_crude <- lmer(formula_crude, data = df_all)
    
    # 2) adjusted
    formula_adjusted <- as.formula(
      paste0(outcome_col, " ~ ", exposure_col, " + period + ",
             paste(covariates, collapse = " + "), 
             " + (1|ID)")
    )
    
    model_adjusted <- lmer(formula_adjusted, data = df_all)
    
    # 3) results
    res_crude    <- broom.mixed::tidy(model_crude, conf.int = TRUE, conf.method = "Wald")
    res_adjusted <- broom.mixed::tidy(model_adjusted, conf.int = TRUE, conf.method = "Wald")
    
    res_crude_fixed <- res_crude %>%
      filter(effect == "fixed") %>%
      select(term, estimate, std.error, statistic, p.value, conf.low, conf.high) %>%
      mutate(model = "Crude")
    
    res_adjusted_fixed <- res_adjusted %>%
      filter(effect == "fixed") %>%
      select(term, estimate, std.error, statistic, p.value, conf.low, conf.high) %>%
      mutate(model = "Adjusted")
    
    res_current <- bind_rows(res_crude_fixed, res_adjusted_fixed) %>%
      mutate(
        outcome   = thyroid,       
        exposure  = pollutant   
      )
    
    all_results[[counter]] <- res_current
    counter <- counter + 1
  }
}

final_results <- bind_rows(all_results)
final_results <- final_results %>%
  mutate(p.value.FDR = p.adjust(p.value, method = "BH"))

write.xlsx(final_results, "model_results.xlsx")

