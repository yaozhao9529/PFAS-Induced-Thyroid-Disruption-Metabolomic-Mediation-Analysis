library("hdmed")
library(readxl)
library(dplyr)

set.seed(123)

thyroid_PFAS <- read.csv("thyroid-PFAS.csv")
cols_keep_id <- c("ID", "period", "UIC")
thyroid_PFAS <- data.frame(thyroid_PFAS[,cols_keep_id], log(thyroid_PFAS[,3:20]))

df_met <- read.csv("Significant_Metabolites.csv")
colnames(df_met)[1:4] <- c("sum_ID", "period", "ID", "Sample_ID")

df_pathway <- read.csv("metabolites for pathway.csv")
pathway_names <- colnames(df_pathway)[12:21]

pathway_list <- setNames(
  lapply(pathway_names, function(p) {
    unique(na.omit(df_pathway$Metabolites_ID[!is.na(df_pathway[[p]])]))
  }),
  gsub("\\.", "_", pathway_names)
)

df_cov <- read.csv("sample-info.csv")
colnames(df_cov)[1] <- "ID"


df_all <- thyroid_PFAS %>%
  left_join(df_cov, by = "ID")

df_all <- df_all %>%
  left_join(df_met, by = c("ID","period"))

any(is.na(df_all))

covariates <- c("age", "BMI", "UIC", "period",
                "education", "employment", "gravidity", "parity",
                "income", "smoking", "alcohol", "coffee", "fetal_sex")

outcomes <- c("TSH", "FT3", "FT4", "T3", "T4", "T4_T3")

pollutants <- c("C9_HFPO", "PFMOAA", "PFO5DoDA", "PFOA",
                "PFNA", "PFDA", "PFUnDA", "PFDoDA", 
                "PFHxS", "PFHpS", "PFOS", "PFESA62")

hdma_results_list  <- list()
i <- 1

for (outcome in outcomes) {
  
  Y_all <- df_all[[outcome]]
  
  for (pollutant in pollutants) {
    
    A_all <- df_all[[pollutant]]
    
    for (p in seq_along(pathway_list)) {
      
      cat("\n analyzing：", pollutant, "→", names(pathway_list)[p], "→", outcome, "\n")
      
      M_all <- df_all[, pathway_list[[p]], drop = FALSE]
      
      if (ncol(M_all) < 2) {
        message("metabolites< 2：", names(pathway_list)[p])
        next
      }
      
      cov_df <- df_all[, covariates, drop = FALSE]
      
      cov_df$education  <- factor(cov_df$education)
      cov_df$employment <- factor(cov_df$employment)
      cov_df$gravidity  <- factor(cov_df$gravidity)
      cov_df$parity     <- factor(cov_df$parity)
      cov_df$income     <- factor(cov_df$income)
      cov_df$period     <- factor(cov_df$period)
      cov_df$fetal_sex  <- factor(cov_df$fetal_sex)
      cov_df$smoking    <- factor(cov_df$smoking)
      cov_df$alcohol    <- factor(cov_df$alcohol)
      cov_df$coffee     <- factor(cov_df$coffee)
      
      C_all <- model.matrix(~ . - 1, data = cov_df)
      C_all <- as.matrix(C_all)
      
      cc <- complete.cases(A_all, Y_all, M_all, C_all)
      A <- A_all[cc]
      Y <- Y_all[cc]
      M <- M_all[cc, , drop = FALSE]
      C <- C_all[cc, , drop = FALSE]

      hdma_out <- tryCatch({
        mediate_hdma(A, M, Y, C1 = C, C2 = C)
      }, error = function(e) {
        message("HDMA error：", conditionMessage(e))
        return(NULL)
      })
      
      if (!is.null(hdma_out) &&
          !is.null(hdma_out$contributions) &&
          nrow(hdma_out$contributions) > 0) {
        
        hdma_results_list[[i]] <- list(
          outcome   = outcome,
          pollutant = pollutant,
          pathway_id = p,
          pathway_name = names(pathway_list)[p],
          contributions = hdma_out$contributions,
          effects       = hdma_out$effects
        )
      } else {
        message("HDMA：no mediator")
      }
      
      i <- i + 1
    }
  }
}

contrib_df_hdma <- bind_rows(
  lapply(hdma_results_list, function(res) {
    cbind(outcome   = res$outcome,
          pollutant = res$pollutant,
          pathway   = res$pathway_name,
          res$contributions)
  })
)

effects_df_hdma <- bind_rows(
  lapply(hdma_results_list, function(res){
    ef <- res$effects
    data.frame(
      outcome   = res$outcome,
      pollutant = res$pollutant,
      pathway   = res$pathway_name,
      indirect_estimate = ef$estimate[ef$effect == "indirect"],
      direct_estimate   = ef$estimate[ef$effect == "direct"],
      total_estimate    = ef$estimate[ef$effect == "total"]
    )
  })
)

pathway_pval_df <- contrib_df_hdma %>%
  group_by(outcome, pollutant, pathway) %>%
  summarise(
    n_mediators = n(),
    fisher_X2 = -2 * sum(log(pmax(ab_pv, .Machine$double.xmin)), na.rm = TRUE),
    df = 2 * n(),
    p_value = pchisq(fisher_X2, df, lower.tail = FALSE),
    .groups = "drop"
  ) %>%
  mutate(FDR = p.adjust(p_value, method = "BH"))

effects_df_hdma <- left_join(effects_df_hdma, pathway_pval_df,
                             by = c("outcome", "pollutant", "pathway"))

hdma_contrib_ci <- contrib_df_hdma %>%
  mutate(
    z_alpha = qnorm(1 - alpha_pv / 2),
    z_beta  = qnorm(1 - beta_pv  / 2),
    se_alpha = abs(alpha) / z_alpha,
    se_beta  = abs(beta)  / z_beta,
    
    se_ab = sqrt((beta^2) * (se_alpha^2) + (alpha^2) * (se_beta^2)),
    
    ab_ci_lower = (alpha * beta) - 1.96 * se_ab,
    ab_ci_upper = (alpha * beta) + 1.96 * se_ab
  )

merged_data <- contrib_df_hdma %>%
  left_join(effects_df_hdma, 
            by = c("outcome", "pollutant", "pathway")) %>%
  select(outcome, pollutant, pathway, mediator, 
         alpha, alpha_pv, beta, beta_pv, alpha_beta, ab_pv,
         indirect_estimate, direct_estimate, total_estimate, 
         n_mediators, fisher_X2, df, p_value, FDR)

write.csv(contrib_df_hdma, "hdma_contributions_all.csv", row.names = FALSE)  
write.csv(effects_df_hdma, "hdma_effects_all.csv", row.names = FALSE)
write.csv(merged_data, "hdma_merged_data.csv", row.names = FALSE)  


