## Function to extract metrics from data and model csv summaries
## species - a vector containing all the species names
## path - the path to the directory containing the csv files
## the function is hard -coded 
## the directory needs to have two folders (Integrated and HRS)
## containing another folder (csv_results)
## lines 11 to 21 can be modified to suit a different folder structure
## also data needs to be already loaded and named 
## PoMS_formatted and HRS_formatted
## hov_POMS_2017 and hrs_2017_sp_names
## returns a dataframe with a set of parameters etracted from the 
## model summaries and datasets

extract_metrics <- function(species, path = getwd()){
  
  HRS_model_file <- list.files(paste0(path, "HRS/csv_results/"), 
                               pattern = paste0(species, "_it16666_ep2e+05.csv"))
  
  integrated_model_summary <- read.csv(paste0(path, "Integrated/csv_results/", 
                                              species, "_it16666_ep2e+05.csv"),
                                       header = T, row.names = 1)
  
  HRS_model_summary <- read.csv(paste0(path, "HRS/csv_results/", 
                                       species, "_it16666_ep2e+05.csv"),
                                header = T, row.names = 1)
  
  n_records_PoMS <- sum(as.numeric(PoMS_formatted$spp_vis[,species]))
  
  n_records_HRS <- sum(as.numeric(HRS_formatted$spp_vis[,species]))
  
  n_sites_PoMS <- n_distinct(hov_POMS_2017[which(hov_POMS_2017$taxon_standardised == species), "GRIDREF_5KM_PREC"])
  
  n_sites_HRS <- n_distinct(hrs_2017_sp_names[which(hrs_2017_sp_names$NAME == species), "GRIDREF_5KM_PREC"])
  
  psi_integrated <- integrated_model_summary["psi.fs", "Mean"]
  psi_HRS <- HRS_model_summary["psi.fs", "Mean"]
  
  a_integrated <- integrated_model_summary["a", "Mean"]
  a_HRS <- HRS_model_summary["a", "Mean"]
  
  dtype1_integrated <- integrated_model_summary["dtype1.p", "Mean"]
  dtype1_HRS <- HRS_model_summary["dtype1.p", "Mean"]
  
  dtype2_integrated <- integrated_model_summary["dtype2.p", "Mean"]
  dtype2_HRS <- HRS_model_summary["dtype2.p", "Mean"]
  
  dtype3_integrated <- integrated_model_summary["dtype3.p", "Mean"]
  dtype3_HRS <- HRS_model_summary["dtype3.p", "Mean"]
  
  beta1_integrated <- integrated_model_summary["beta1", "Mean"]
  beta1_HRS <- HRS_model_summary["beta1", "Mean"]
  
  beta2_integrated <- integrated_model_summary["beta2", "Mean"]
  beta2_HRS <- HRS_model_summary["beta2", "Mean"]
  
  beta3_integrated <- integrated_model_summary["beta3", "Mean"]
  beta3_HRS <- HRS_model_summary["beta3", "Mean"]
  
  Prec_psi_integrated <- 1/integrated_model_summary["psi.fs", "SD"]^2
  Prec_psi_HRS <- 1/HRS_model_summary["psi.fs", "SD"]^2
  
  Prec_a_integrated <- 1/integrated_model_summary["a", "SD"]^2
  Prec_a_HRS <- 1/HRS_model_summary["a", "SD"]^2
  
  Prec_dtype1_integrated <- 1/integrated_model_summary["dtype1.p", "SD"]^2
  Prec_dtype1_HRS <- 1/HRS_model_summary["dtype1.p", "SD"]^2
  
  Prec_dtype2_integrated <- 1/integrated_model_summary["dtype2.p", "SD"]^2
  Prec_dtype2_HRS <- 1/HRS_model_summary["dtype2.p", "SD"]^2
  
  Prec_dtype3_integrated <- 1/integrated_model_summary["dtype3.p", "SD"]^2
  Prec_dtype3_HRS <- 1/HRS_model_summary["dtype3.p", "SD"]^2
  
  Prec_beta1_integrated <- 1/integrated_model_summary["beta1", "SD"]^2
  Prec_beta1_HRS <- 1/HRS_model_summary["beta1", "SD"]^2
  
  Prec_beta2_integrated <- 1/integrated_model_summary["beta2", "SD"]^2
  Prec_beta2_HRS <- 1/HRS_model_summary["beta2", "SD"]^2
  
  Prec_beta3_integrated <- 1/integrated_model_summary["beta3", "SD"]^2
  Prec_beta3_HRS <- 1/HRS_model_summary["beta3", "SD"]^2
  
  
  Rhat_psi_integrated <- integrated_model_summary["psi.fs", "PSRF"]
  Rhat_psi_HRS <- HRS_model_summary["psi.fs", "PSRF"]
  
  Rhat_a_integrated <- integrated_model_summary["a", "PSRF"]
  Rhat_a_HRS <- HRS_model_summary["a", "PSRF"]
  
  Rhat_dtype1_integrated <- integrated_model_summary["dtype1.p", "PSRF"]
  Rhat_dtype1_HRS <- HRS_model_summary["dtype1.p", "PSRF"]
  
  Rhat_dtype2_integrated <- integrated_model_summary["dtype2.p", "PSRF"]
  Rhat_dtype2_HRS <- HRS_model_summary["dtype2.p", "PSRF"]
  
  Rhat_dtype3_integrated <- integrated_model_summary["dtype3.p", "PSRF"]
  Rhat_dtype3_HRS <- HRS_model_summary["dtype3.p", "PSRF"]
  
  Rhat_beta1_integrated <- integrated_model_summary["beta1", "PSRF"]
  Rhat_beta1_HRS <- HRS_model_summary["beta1", "PSRF"]
  
  Rhat_beta2_integrated <- integrated_model_summary["beta2", "PSRF"]
  Rhat_beta2_HRS <- HRS_model_summary["beta2", "PSRF"]
  
  Rhat_beta3_integrated <- integrated_model_summary["beta3", "PSRF"]
  Rhat_beta3_HRS <- HRS_model_summary["beta3", "PSRF"]
  
  rm(list = c("integrated_model_summary", "HRS_model_summary"))
  
  model_outputs <- data.frame(species = species, n_records_PoMS = n_records_PoMS, 
                              n_records_HRS = n_records_HRS, n_sites_PoMS = n_sites_PoMS, 
                              n_sites_HRS = n_sites_HRS, psi_integrated = psi_integrated,
                              psi_HRS = psi_HRS, a_integrated = a_integrated,
                              a_HRS = a_HRS, dtype1_integrated = dtype1_integrated,
                              dtype1_HRS = dtype1_HRS, dtype2_integrated = dtype2_integrated,
                              dtype2_HRS = dtype2_HRS, dtype3_integrated = dtype3_integrated,
                              dtype3_HRS = dtype3_HRS, beta1_integrated = beta1_integrated,
                              beta1_HRS = beta1_HRS, beta2_integrated = beta2_integrated,
                              beta2_HRS = beta2_HRS, beta3_integrated = beta3_integrated,
                              beta3_HRS = beta3_HRS, Prec_psi_integrated = Prec_psi_integrated, 
                              Prec_psi_HRS = Prec_psi_HRS, Prec_a_integrated = Prec_a_integrated, 
                              Prec_a_HRS = Prec_a_HRS, Prec_dtype1_integrated = Prec_dtype1_integrated, 
                              Prec_dtype1_HRS = Prec_dtype1_HRS, Prec_dtype2_integrated = Prec_dtype2_integrated, 
                              Prec_dtype2_HRS = Prec_dtype2_HRS, Prec_dtype3_integrated = Prec_dtype3_integrated, 
                              Prec_dtype3_HRS = Prec_dtype3_HRS, Prec_beta1_integrated = Prec_beta1_integrated,
                              Prec_beta1_HRS = Prec_beta1_HRS, Prec_beta2_integrated = Prec_beta2_integrated,
                              Prec_beta2_HRS = Prec_beta2_HRS, Prec_beta3_integrated = Prec_beta3_integrated,
                              Prec_beta3_HRS = Prec_beta3_HRS, Rhat_psi_integrated = Rhat_psi_integrated, 
                              Rhat_psi_HRS = Rhat_psi_HRS, Rhat_a_integrated = Rhat_a_integrated,
                              Rhat_a_HRS = Rhat_a_HRS, Rhat_dtype1_integrated = Rhat_dtype1_integrated,
                              Rhat_dtype1_HRS = Rhat_dtype1_HRS, Rhat_dtype2_integrated = Rhat_dtype2_integrated,
                              Rhat_dtype2_HRS = Rhat_dtype2_HRS, Rhat_dtype3_integrated = Rhat_dtype3_integrated,
                              Rhat_dtype3_HRS = Rhat_dtype3_HRS, Rhat_beta1_integrated = Rhat_beta1_integrated,
                              Rhat_beta1_HRS = Rhat_beta1_HRS, Rhat_beta2_integrated = Rhat_beta2_integrated,
                              Rhat_beta2_HRS = Rhat_beta2_HRS, Rhat_beta3_integrated = Rhat_beta3_integrated,
                              Rhat_beta3_HRS = Rhat_beta3_HRS)
  
  return(model_outputs)
  
}
