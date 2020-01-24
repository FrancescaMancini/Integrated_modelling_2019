#########################################################
## Script for integrated modelling of hoverfly data
## Author: Francesca
## Date created: 2019/02/14
## Date modified: 2020-01-24
########################################################

library(R2jags)
library(dplyr)
library(sparta)
library(ggplot2)

## Assess convergence across all chains ----

source("./Functions/combine_chains_modified.R")
source("./Functions/extract_metrics.R")

dir.create("./Outputs/HRS/csv_results")

dir.create("./Outputs/Integrated/csv_results")

combine_chains(data_dir = "Outputs/HRS",
               output_dir = "Outputs/HRS/csv_results",
               target = 16666)

combine_chains(data_dir = "Outputs/Integrated",
               output_dir = "Outputs/Integrated/csv_results", model_type = "integrated",
               target = 16666)


## Extract model results and plot ----

load("./Data/HRS_formatted.rdata")
load("./Data/POMS_formatted.rdata")

hrs_2017_sp_names <- read.csv("./Data/HRS_2017_spNames.csv", 
                              header = TRUE, stringsAsFactors = FALSE)

hov_POMS_2017 <- read.csv("./Data/POMS_2017_5Km.csv", 
                          header = TRUE, stringsAsFactors = FALSE)


species_list <- read.csv("./Data/common_spNames.csv", 
                         header = TRUE, stringsAsFactors = FALSE)


model_outputs <- lapply(species_list$Species[-1], extract_metrics, 
                        path = "./Outputs/")

model_outputs_df <- do.call(rbind, model_outputs)

str(model_outputs_df)

# write.csv(model_outputs_df, "./Outputs/hoverflies_model_metrics.csv", row.names = F)

model_outputs_df <- read.csv("./Output/hoverflies_model_metrics_best_species.csv",
                             header = T)

str(model_outputs_df)

# subset the results to only keep those that converged
# and plot some main stats


Prec_a <- ggplot(model_outputs_df%>%
                   filter_at(vars(contains("Rhat_a")), all_vars(.<=1.15))) +
  geom_point(aes(x = Prec_a_HRS, y = Prec_a_integrated)) +
  geom_smooth(aes(x = Prec_a_HRS, y = Prec_a_integrated),method = "lm", colour = "black") +
  geom_abline(intercept = 0, slope = 1,linetype = "dashed") +
  coord_equal()

Prec_a_log <- ggplot(model_outputs_df%>%
                       filter_at(vars(contains("Rhat_a")), all_vars(.<=1.15))) +
  geom_point(aes(x = log(Prec_a_HRS), y = log(Prec_a_integrated))) +
  stat_smooth(aes(x = log(Prec_a_HRS), y = log(Prec_a_integrated)),method = "lm", colour = "black") +
  geom_abline(intercept = 0, slope = 1,linetype = "dashed") 


Prec_psi <- ggplot(model_outputs_df%>%
                     filter_at(vars(contains("Rhat_psi")), all_vars(.<=1.15))) +
  geom_point(aes(x = Prec_psi_HRS, y = Prec_psi_integrated)) +
  stat_smooth(aes(x = Prec_psi_HRS, y = Prec_psi_integrated),method = "lm") +
  geom_abline(intercept = 0, slope = 1,linetype = "dashed") +
  coord_equal()

Prec_psi_log <- ggplot(model_outputs_df%>%
                         filter_at(vars(contains("Rhat_psi")), all_vars(.<=1.15))) +
  geom_point(aes(x = log(Prec_psi_HRS), y = log(Prec_psi_integrated))) +
  stat_smooth(aes(x = log(Prec_psi_HRS), y = log(Prec_psi_integrated)),method = "lm") +
  geom_abline(intercept = 0, slope = 1,linetype = "dashed") +
  coord_equal()


Prec_beta1 <- ggplot(model_outputs_df%>%
                       filter_at(vars(contains("Rhat_beta1")), all_vars(.<=1.15))) +
  geom_point(aes(x = Prec_beta1_HRS, y = Prec_beta1_integrated)) +
  stat_smooth(aes(x = Prec_beta1_HRS, y = Prec_beta1_integrated),method = "lm") +
  geom_abline(intercept = 0, slope = 1,linetype = "dashed")


Prec_beta2 <- ggplot(model_outputs_df%>%
                       filter_at(vars(contains("Rhat_beta2")), all_vars(.<=1.15))) +
  geom_point(aes(x = Prec_beta2_HRS, y = Prec_beta2_integrated)) +
  stat_smooth(aes(x = Prec_beta2_HRS, y = Prec_beta2_integrated),method = "lm") +
  geom_abline(intercept = 0, slope = 1,linetype = "dashed")


Prec_beta3 <- ggplot(model_outputs_df%>%
                       filter_at(vars(contains("Rhat_beta3")), all_vars(.<=1.15))) +
  geom_point(aes(x = Prec_beta3_HRS, y = Prec_beta3_integrated)) +
  stat_smooth(aes(x = Prec_beta3_HRS, y = Prec_beta3_integrated),method = "lm") +
  geom_abline(intercept = 0, slope = 1,linetype = "dashed")


source("../Descriptive Stats/multiplot.R")

convergence_psi <- ggplot(model_outputs_df) +
  geom_point(aes(x = Rhat_psi_HRS, y = Rhat_psi_integrated)) +
  stat_smooth(aes(x = Rhat_psi_HRS, y = Rhat_psi_integrated),method = "lm") +
  geom_abline(intercept = 0, slope = 1,linetype = "dashed")

convergence_a <- ggplot(model_outputs_df) +
  geom_point(aes(x = Rhat_a_HRS, y = Rhat_a_integrated)) +
  stat_smooth(aes(x = Rhat_a_HRS, y = Rhat_a_integrated),method = "lm") +
  geom_abline(intercept = 0, slope = 1,linetype = "dashed") 



conv_a_records_integrated <- ggplot(model_outputs_df) +
  geom_point(aes(x = n_records_PoMS, y = Rhat_a_integrated)) +
  stat_smooth(aes(x = n_records_PoMS, y = Rhat_a_integrated),method = "lm") 



conv_a_sites_integrated <- ggplot(model_outputs_df) +
  geom_point(aes(x = n_sites_PoMS, y = Rhat_a_integrated)) +
  stat_smooth(aes(x = n_sites_PoMS, y = Rhat_a_integrated),method = "lm") 

conv_a_records_HRS <- ggplot(model_outputs_df) +
  geom_point(aes(x = n_records_HRS, y = Rhat_a_integrated)) +
  stat_smooth(aes(x = n_records_HRS, y = Rhat_a_integrated),method = "lm") 

conv_a_sites_HRS <- ggplot(model_outputs_df) +
  geom_point(aes(x = n_sites_HRS, y = Rhat_a_integrated)) +
  stat_smooth(aes(x = n_sites_HRS, y = Rhat_a_integrated),method = "lm") 


conv_a_records_total <- ggplot(model_outputs_df) +
  geom_point(aes(x = (n_records_PoMS + n_records_HRS), y = Rhat_a_integrated)) +
  stat_smooth(aes(x = (n_records_PoMS + n_records_HRS), y = Rhat_a_integrated),method = "lm")

conv_a_sites_total <- ggplot(model_outputs_df) +
  geom_point(aes(x = (n_sites_PoMS + n_sites_HRS), y = Rhat_a_integrated)) +
  stat_smooth(aes(x = (n_sites_PoMS + n_sites_HRS), y = Rhat_a_integrated),method = "lm")


Prec_a_records_PoMS <- ggplot(model_outputs_df) +
  geom_point(aes(x = n_records_PoMS, y = Prec_a_integrated)) +
  stat_smooth(aes(x = n_records_PoMS, y = Prec_a_integrated),method = "lm") 

Prec_a_sites_PoMS <- ggplot(model_outputs_df) +
  geom_point(aes(x = n_sites_PoMS, y = Prec_a_integrated)) +
  stat_smooth(aes(x = n_sites_PoMS, y = Prec_a_integrated),method = "lm") 


Prec_a_records_HRS <- ggplot(model_outputs_df) +
  geom_point(aes(x = n_records_HRS, y = Prec_a_integrated)) +
  stat_smooth(aes(x = n_records_HRS, y = Prec_a_integrated),method = "lm") 

Prec_a_sites_HRS <- ggplot(model_outputs_df) +
  geom_point(aes(x = n_sites_HRS, y = Prec_a_integrated)) +
  stat_smooth(aes(x = n_sites_HRS, y = Prec_a_integrated),method = "lm") 

Prec_a_records_total <- ggplot(model_outputs_df) +
  geom_point(aes(x = (n_records_PoMS + n_records_HRS), y = Prec_a_integrated)) +
  stat_smooth(aes(x = (n_records_PoMS + n_records_HRS), y = Prec_a_integrated),method = "lm")

Prec_a_sites_total <- ggplot(model_outputs_df) +
  geom_point(aes(x = (n_sites_PoMS + n_sites_HRS), y = Prec_a_integrated)) +
  stat_smooth(aes(x = (n_sites_PoMS + n_sites_HRS), y = Prec_a_integrated),method = "lm")

