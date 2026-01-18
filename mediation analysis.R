library(dplyr)
library(readxl)
library(mediation)
library(lme4) 
library(stats)
library(openxlsx)  

set.seed(123)

thyroid_PFAS <- read.csv("thyroid-PFAS.csv")
thyroid_PFAS <- data.frame(thyroid_PFAS[,c(1,2,24)], log(thyroid_PFAS[,3:20]))
df_covariate <- read.csv("sample-info.csv")
thyroid_PFAS <- thyroid_PFAS %>%
  left_join(df_covariate, by = "ID")

met_data <- read.csv("metabolites.csv")
mid_columns <- grep("^MID", colnames(met_data), value = TRUE)
selected_data <- met_data[ , metabolite_list, drop = FALSE]
selected_data <- log(selected_data)

colnames(selected_data) <- paste0("M", seq_len(ncol(selected_data)))
df_met <- cbind(met_data[, c(1:1)], met_data[, c(4:5)], met_data[, c(7:7)], selected_data)
colnames(df_met)[1:4] <- c("sum-ID", "period", "group", "Sample-ID")

outcomes <- c("TSH", "FT3", "FT4", "T3", "T4", "T4_T3")  

pollutants <- c("C9_HFPO", "PFMOAA", "PFO5DoDA", "PFOA",
                "PFNA", "PFDA", "PFUnDA", "PFDoDA", 
                "PFHxS", "PFHpS", "PFOS", "PFESA62")

covariates <- c("age", "BMI", "UIC", 
                "education", "employment", "gravidity", "parity", 
                "income", "smoking", "alcohol", "coffee", "fetal_sex")


mediation_analysis <- function(outcome = outcome, pollutants = pollutants, covariates = covariates,
                               medium_met_id = medium_met_id, 
                               df_met = df_met, df_outcomes_PFAS = thyroid_PFAS) {
  
  extracted_met <- df_met[, medium_met_id, drop = FALSE]
  
  mediation_data <- data.frame(df_outcomes_PFAS, extracted_met)

  results <- list()
  results_df <- data.frame()
  
  for (pollutant in pollutants) {
    for (metabolite in medium_met_id) {
      try({
        covariate_formula <- paste(covariates, collapse = " + ")
        
        mediator_model <- lmer(as.formula(paste(metabolite, "~", pollutant, "+", covariate_formula, "+ period + (1|ID)")),
                               data = mediation_data, REML = FALSE)
        
        outcome_model <- lmer(as.formula(paste(outcome, "~", pollutant, "+", metabolite, "+", covariate_formula, "+ period + (1|ID)")),
                              data = mediation_data, REML = FALSE)
        
        med_out <- mediate(mediator_model, outcome_model, treat = pollutant, mediator = metabolite, sims = 1000)
        
        summary_med <- summary(med_out)
        results[[paste(pollutant, metabolite, outcome, sep = "_")]] <- summary_med
        
        temp_df <- data.frame(
          Pollutant = pollutant,
          Metabolite = metabolite,
          Outcome = outcome,
          ACME = summary_med$d0,               
          ACME_Lower = summary_med$d0.ci[1],  
          ACME_Upper = summary_med$d0.ci[2],
          ACME_p = summary_med$d0.p,          
          ADE = summary_med$z0,             
          ADE_Lower = summary_med$z0.ci[1],   
          ADE_Upper = summary_med$z0.ci[2], 
          ADE_p = summary_med$z0.p,           
          Total_Effect = summary_med$tau.coef, 
          Total_Lower = summary_med$tau.ci[1], 
          Total_Upper = summary_med$tau.ci[2],
          Total_p = summary_med$tau.p,       
          Prop_Mediated = summary_med$n0,    
          Prop_Lower = summary_med$n0.ci[1],  
          Prop_Upper = summary_med$n0.ci[2],  
          Prop_p = summary_med$n0.p     
        )
        
        results_df <- rbind(results_df, temp_df)
        
        print(paste("mediation analysis:", pollutant, "->", metabolite, "->", outcome))
      }, silent = TRUE)
    }
  }
  
  results_df$ACME_FDR <- p.adjust(results_df$ACME_p, method = "BH")
  results_df$ADE_FDR <- p.adjust(results_df$ADE_p, method = "BH")
  results_df$Total_FDR <- p.adjust(results_df$Total_p, method = "BH")
  results_df$Prop_FDR <- p.adjust(results_df$Prop_p, method = "BH")
  
  return(results_df)
}


wb <- createWorkbook()

medium_met_id <- colnames(selected_data)

all_results <- data.frame()

for (outcome in outcomes) {
  print(paste("analysis:", outcome))
  result <- mediation_analysis(outcome, pollutants, covariates, medium_met_id, df_met, thyroid_PFAS)
  
  write.csv(result, paste0(outcome, "_Mediation_Analysis_Results.csv"), row.names = FALSE)
  
  sheet_name <- substr(paste0(outcome, "_Mediation"), 1, 31)
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet = sheet_name, x = result)  
  
  all_results <- rbind(all_results, result)
}

addWorksheet(wb, "All_Mediation")  
writeData(wb, "All_Mediation", all_results)

saveWorkbook(wb,  
             file      = "Mediation_Analysis_Results.xlsx",
             overwrite = TRUE)  


