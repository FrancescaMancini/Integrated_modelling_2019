#########################################################
## Script for integrated modelling of hoverfly data
## Author: Francesca
## Date created: 2019/02/14
## Date modified: 2019-12-05
########################################################

library(BRCmap)
library(RODBC)
library(R2jags)
library(dplyr)
library(sparta)
library(ggplot2)
#library(snowfall)


# load the hoverfly recording scheme data
hrs <- read.csv("..\\Data\\HRS\\RD_HRS.csv", header = T, stringsAsFactors = F)
str(hrs)

summary(hrs$DT_ID)
# I should delete everything that has DT_ID = 8, no date provided

hrs <- subset(hrs, DT_ID != 8)

# convert to date format
hrs$startdate <- as.Date(hrs$TO_STARTDATE, format = "%d/%m/%Y")
hrs$enddate <- as.Date(hrs$TO_ENDDATE, format = "%d/%m/%Y")

# retain only 2017 records
hrs_2017 <- subset(hrs, format(hrs$enddate,"%Y")==2017)
# retain only records with precision 1Km or less
hrs_2017 <- subset(hrs_2017, TO_PRECISION <= 1000)

POMS <- read.csv("P:\\NEC06214_UK Pollinator Monitoring and Research Partnership\\Data and analysis\\data outputs current versions\\tblEXPORT_1kmPanTrap_insects.csv", 
                 header = T, stringsAsFactors = F)
str(POMS)

# select only flies
hov_POMS <- subset(POMS, taxon_group == "insect - true fly (Diptera)")

# select only 2017
hov_POMS$dates <- as.Date(hov_POMS$date, format = "%d/%m/%Y %H:%M:%S")
hov_POMS_2017 <- subset(hov_POMS, format(hov_POMS$dates,"%Y")==2017)

# at 1Km resolution, there is no overlap between HRS and POMS sites
# change precision of grid references to be 5Km
hrs_2017$GRIDREF_5KM_PREC <- reformat_gr(hrs_2017$TO_GRIDREF, prec_out = 5000, precision = hrs_2017$TO_PRECISION)
hov_POMS_2017$GRIDREF_5KM_PREC <- reformat_gr(hov_POMS_2017$sample_gridref, prec_out = 5000)

# write.csv(hrs_2017, "HRS_2017_5Km.csv", row.names = FALSE)
# write.csv(hov_POMS_2017, "POMS_2017_5Km.csv", row.names = FALSE)


length(which(unique(hov_POMS_2017$GRIDREF_5KM_PREC) %in% unique(hrs_2017$GRIDREF_5KM_PREC)))
# at 5Km resolution 22 sites overlap

# extract the overlapping sites
overlap <- unique(hov_POMS_2017[which(hov_POMS_2017$GRIDREF_5KM_PREC %in% hrs_2017$GRIDREF_5KM_PREC), "GRIDREF_5KM_PREC"])


# visualise spatial overlap between the sites from the two data sources
#plot first all the HRS sites and then highlight those overlapping with PoMS in red
par(mar = c(0.1,0.1,1,0.1))
plot_GIS(UK$britain, new.window = FALSE, show.axis = FALSE, 
         show.grid = FALSE, xlab = "", ylab = "", main = "HRS sites overlapping with PoMS")

plotUK_gr(hrs_2017$GRIDREF_5KM_PREC, gr_prec = 5000, border = "grey")
plotUK_gr(overlap, gr_prec = 5000, border = "red")

# and vice versa
par(mar = c(0.1,0.1,1,0.1))
plot_GIS(UK$britain, new.window = FALSE, show.axis = FALSE, 
         show.grid = FALSE, xlab = "", ylab = "", main = "PoMS sites overlapping with HRS")

plotUK_gr(hov_POMS_2017$GRIDREF_5KM_PREC, gr_prec = 5000, border = "grey")
plotUK_gr(overlap, gr_prec = 5000, border = "red")


# obtain species names from concept codes -----

# if(.Platform$GUI =="RStudio"){
#   brc_pass = winDialogString("Enter database password","")
#   channel = odbcConnect("FRAMAN",pwd=brc_pass, believeNRows = FALSE)
#   rm(brc_pass)
# }else{
#   channel = odbcConnect("FRAMAN", believeNRows = FALSE)
# }
# 
# sp_names <- sqlQuery(channel, paste("select name, concept from brc.taxa_taxon_register where concept in (", 
#                                    paste(shQuote(unique(hrs_2017[, "CONCEPT"]), type="sh"), collapse = ","), ") and valid in ('V', 'X')"))
#
# write.csv(sp_names, "../Data/HRS/Sp_Names.csv", row.names = FALSE)

sp_names <- read.csv("../Data/HRS/Sp_Names.csv")

str(sp_names)

# returns 240 species names, when we only have 238 species codes
which(duplicated(sp_names$CONCEPT) == TRUE)

# species Dip_3430 and Dip_3434 have two associated names: Eristalis abusivus/abusiva and Eristalis intricarius/intricaria
# in PoMS these species are referred as abusivus and intricarius, therefore we are using same names here

sp_names <- sp_names[-c(94, 98), ]

hrs_2017_sp_names <- left_join(hrs_2017, sp_names, by = "CONCEPT")


## Run the model for all species with records in both datasets -----

# How many species have records in both datasets?

common_species <- intersect(hov_POMS_2017$taxon_standardised, hrs_2017_sp_names$NAME)
length(common_species)

# there are 63 species that have records in both datasets

# Data reformatting ----
# y1 is the detection history for one species from the HRS dataset
# y2 is the detection history fort he same species from the PoMS dataset

# Site1 is a list of sites visited in HRS dataset (no unique)
# Site2 is a list of sites visited in PoMS (no unique)
# nsite is the total number of sites visited in both datasets (unique)

# nvisit1 is the numbe of visits in the HRS dataset
# nvisit2 is the number of visits in the PoMS dataset

# DATATYPE2 is a column that indicates if that visit produced a medium list length (1) or not (0) (for HRS dataset only)
# DATATYPE3 same as above for long list length

# can probably use formatOccData from sparta

# create year to supply as value for closure_period to formatOccData
# this is because we cannot use date as to survey = 
# because having decreased te precision of the gridref
# we might have multiple visits in the same 5 Km site on the same day
# and we want to keep them as replicates
hrs_2017_sp_names$YEAR <- format(hrs_2017_sp_names$enddate,"%Y")

HRS_formatted <- formatOccData(taxa = hrs_2017_sp_names$NAME, 
                               survey = paste(hrs_2017_sp_names$enddate,hrs_2017_sp_names$TO_GRIDREF, sep = "-"), 
                               # create survey value from date and 1Km gridref to keep all visits as replicates
                               site = hrs_2017_sp_names$GRIDREF_5KM_PREC,
                               closure_period = hrs_2017_sp_names$YEAR)

# create Julian day variable from the date in visit
# The variable visit at the moment is a combination of the 5Km grid ref
# the date and the 1Km grid ref
# the code extracts the date from these character strings (characters number 7 to 16)
# formats them as dates and converts them to numbers
JulDate1 <- as.numeric(format(as.POSIXlt(substr(HRS_formatted$spp_vis$visit, 7, 16),
                                         format = "%Y-%m-%d"), "%j"))

Site1 <- HRS_formatted$occDetdata$site
summary(HRS_formatted$occDetdata$L)
DATATYPE2 <- HRS_formatted$occDetdata %>%
  mutate(L = case_when(L < 2 ~ 0,
                       L < 4 ~ 1,
                       L >= 4 ~ 0)) %>%
  select(L)


DATATYPE3 <- HRS_formatted$occDetdata %>%
  mutate(L = case_when(L < 4 ~ 0,
                       L >= 4 ~ 1)) %>%
  select(L)

# for now just change the dates by a number of day that is equal to the pantrap station number (1 to 5)
# eventually we want to modify formatOccData to allow for another identifier for survey
hov_POMS_2017$dates_new <- hov_POMS_2017$dates + hov_POMS_2017$pan_trap_station

# create year to supply as value for closure_period to formatOccData
# this is because we cannot use date as to survey = 
# because having decreased te precision of the gridref
# we might have multiple visits in the same 5 Km site on the same day
# and we want to keep them as replicates

hov_POMS_2017$YEAR <- format(hov_POMS_2017$dates_new,"%Y")

PoMS_formatted <- formatOccData(taxa = hov_POMS_2017$taxon_standardised,
                                survey = paste(hov_POMS_2017$dates_new, hov_POMS_2017$sample_gridref, sep = "-"),
                                # create survey value from date and 1Km gridref to keep all visits as replicates
                                site = hov_POMS_2017$GRIDREF_5KM_PREC,
                                closure_period = hov_POMS_2017$YEAR)

# create Julian day variable from the date in visit
# The variable visit at the moment is a combination of the 5Km grid ref
# the date and the 1Km grid ref
# the code extracts the date from these character strings (characters number 7 to 16)
# formats them as dates and converts them to numbers
JulDate2 <- as.numeric(format(as.POSIXlt(substr(PoMS_formatted$spp_vis$visit, 7, 16),
                                         format = "%Y-%m-%d"), "%j")) 


Site2 <- PoMS_formatted$occDetdata$site


nsite <- length(unique(c(as.character(Site1), as.character(Site2))))

# MCMC settings
ni <- 25000
nt <- 3
nb <- 10000
nc <- 3

integrated_model <- function(species) { 
  
  write.table(paste("Model for species",species,"has begun at ",date(),sep=" "),
              paste("./Output/progress.txt",sep=""),
              append=TRUE,col.names=FALSE,row.names=FALSE,quote=FALSE)
  
  # format the data - HRS
  
  y1 <- as.integer(HRS_formatted$spp_vis[, species]) 
  
  nvisit1 <- length(y1)
  
  # format the data - POMS
  y2 <- as.integer(PoMS_formatted$spp_vis[, species])
  
  nvisit2 <- length(y2)
  
  
  ## Integrated model - Julian data as gaussian density function ----
  
  jags_integrated_data <- list(y1 = y1, y2 = y2, Site1 = Site1, Site2 = Site2, nvisit1 = nvisit1,
                               nvisit2 = nvisit2, DATATYPE2 = DATATYPE2, DATATYPE3 = DATATYPE3, 
                               nsite = nsite, JulDate1 = JulDate1, JulDate2 = JulDate2)
  
  inits_integrated <- function() list(a = 2, dtype1.p = -0.5, dtype2.p = 1.1, 
                                      dtype3.p = 2, pan.p = 0.2, beta1 = 180, 
                                      beta2 = 16, beta3 = 0.1)
  
  params_integrated <- c("a", "dtype1.p", "dtype2.p", "dtype3.p", "pan.p", "f_x",
                         "psi.fs", "p1", "p2", "beta1", "beta2", "beta3")
  
  
  out_integrated_model <- jags(data = jags_integrated_data, parameters.to.save = params_integrated,
                               model.file = "integrated_model_3.txt", inits = inits_integrated,
                               n.chains=nc, n.iter=ni, n.thin=nt, n.burnin=nb)
  
  
  save(out_integrated_model, file = paste0("./Output/", species, "_integrated_model.rdata"))
  
  recompile(out_integrated_model)
  out_integrated_model_updated <- update(out_integrated_model, n.iter = 1000, n.thin = 3)
  
  save(out_integrated_model_updated, file = paste0("./Output/", species, "_integrated_model_1000.rdata"))
}

# ## Run on the cluster ##
# hosts<-as.character(read.table(Sys.getenv('PBS_NODEFILE'),header=FALSE)[,1]) # read the nodes to use
# 
# sfSetMaxCPUs(length(hosts)) # ensure that snowfall can cope with this many hosts
# 
# sfInit(parallel=TRUE,type="MPI",cpus=length(hosts),useRscript=TRUE) # initialise the connection
# 
# sfExportAll()
# sfLibrary(R2jags)
# 
# sfClusterApplyLB(common_species, integrated_model)
# 
# # Stop snowfall connection
# sfStop()


HRS_model <- function(species) { 
  
  write.table(paste("Model for species",species,"has begun at ",date(),sep=" "),
              paste("./Output/progress.txt",sep=""),
              append=TRUE,col.names=FALSE,row.names=FALSE,quote=FALSE)
  
  # format the data - HRS
  
  y1 <- as.integer(HRS_formatted$spp_vis[, species]) 
  
  nvisit1 <- length(y1)
  
  ## HRS-only model - Julian data as gaussian density function ----
  jags_HRS_data <- list(y1 = y1, Site1 = Site1, nvisit1 = nvisit1,
                        DATATYPE2 = DATATYPE2, DATATYPE3 = DATATYPE3,
                        nsite = length(unique(Site1)), JulDate1 = JulDate1)
  
  inits_HRS <- function() list(a = 2, dtype1.p = -0.5, dtype2.p = 1.1,
                               dtype3.p = 2, beta1 = 180,
                               beta2 = 16, beta3 = 0.1)
  
  params_HRS <- c("a", "dtype1.p", "dtype2.p", "dtype3.p", "psi.fs",
                  "p1", "beta1", "beta2", "beta3", "f_x")
  
  
  # Call JAGS from R
  
  
  out_HRS_model <- jags(data = jags_HRS_data, parameters.to.save = params_HRS,
                        model.file = "HRS_model_3.txt", inits = inits_HRS,
                        n.chains=nc, n.iter=ni, n.thin=nt, n.burnin=nb)
  
  
  save(out_HRS_model, file = paste0("./Output/", species, "_HRS_model.rdata"))
  
  recompile(out_HRS_model)
  out_HRS_model_updated <- update(out_HRS_model, n.iter = 1000, n.thin = 3)
  
  save(out_HRS_model_updated, file = paste0("./Output/", species, "_HRS_model_1000.rdata"))
}


# ## Run on the cluster ##
# hosts<-as.character(read.table(Sys.getenv('PBS_NODEFILE'),header=FALSE)[,1]) # read the nodes to use
# 
# sfSetMaxCPUs(length(hosts)) # ensure that snowfall can cope with this many hosts
# 
# sfInit(parallel=TRUE,type="MPI",cpus=length(hosts),useRscript=TRUE) # initialise the connection
# 
# sfExportAll()
# sfLibrary(R2jags)
# 
# sfClusterApplyLB(common_species, HRS_model)
# 
# # Stop snowfall connection
# sfStop()

## Assess convergence across all chains ----

source("combine_chains_HRS.R")

combine_chains(data_dir = "../Integrated_modelling/Output/Rerun/HRS",
               target = 2000)

source("combine_chains_integrated.R")


combine_chains(data_dir = "../Integrated_modelling/Output/Rerun/Integrated",
               target = 2000)


## Extract model results and plot ----

best_species <- read.table("../Integrated_modelling/20_species_most_data.txt", stringsAsFactors = F)

extract_metrics <- function(species, path = getwd()){
  
  integrated_model_file <- list.files(paste0(path, "Integrated/"), 
                                      pattern = paste0(species, "_it2000_ep201000.csv"))
  

  HRS_model_file <- list.files(paste0(path, "HRS/"), 
                               pattern = paste0(species, "_it2000_ep201000.csv"))
  
  integrated_model_summary <- read.csv(paste0(path, "Integrated/", integrated_model_file),
                                       header = T, row.names = 1)
  HRS_model_summary <- read.csv(paste0(path, "HRS/", HRS_model_file),
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

model_outputs <- lapply(best_species$x, extract_metrics, 
                        path = "../Integrated_modelling/Output/Rerun/")

model_outputs_df <- do.call(rbind, model_outputs)

str(model_outputs_df)

# write.csv(model_outputs_df, "./Output/hoverflies_model_metrics_best_species.csv", row.names = F)

model_outputs_df <- read.csv("./Output/hoverflies_model_metrics_best_species.csv",
                             header = T)

str(model_outputs_df)

# model_outputs_df_sub <- subset(model_outputs_df, species != "Volucella inanis")

# Prec_a_int_converged <- ggplot(model_outputs_df) +
#   geom_point(aes(x = Prec_a_HRS, y = Prec_a_integrated, 
#                  color = cut(Rhat_a_integrated, c(0, 1.1,100)))) +
#   stat_smooth(aes(x = Prec_a_HRS, y = Prec_a_integrated),method = "lm", colour = "black") +
#   geom_abline(intercept = 0, slope = 1,linetype = "dashed") +
#   scale_color_manual(name = "Convergence", 
#                      values = c("cadetblue", "darkred"),
#                      labels = c("Converged (Rhat < 1.1)", "Not converged (Rhat > 1.1")) 

Prec_a <- ggplot(model_outputs_df) +
  geom_point(aes(x = Prec_a_HRS, y = Prec_a_integrated)) +
  stat_smooth(aes(x = Prec_a_HRS, y = Prec_a_integrated),method = "lm", colour = "black") +
  geom_abline(intercept = 0, slope = 1,linetype = "dashed") +
  coord_equal()

Prec_a_log <- ggplot(model_outputs_df) +
  geom_point(aes(x = log(Prec_a_HRS), y = log(Prec_a_integrated))) +
  stat_smooth(aes(x = log(Prec_a_HRS), y = log(Prec_a_integrated)),method = "lm", colour = "black") +
  geom_abline(intercept = 0, slope = 1,linetype = "dashed") 


Prec_psi <- ggplot(model_outputs_df) +
  geom_point(aes(x = Prec_psi_HRS, y = Prec_psi_integrated)) +
  stat_smooth(aes(x = Prec_psi_HRS, y = Prec_psi_integrated),method = "lm") +
  geom_abline(intercept = 0, slope = 1,linetype = "dashed") +
  coord_equal()

Prec_psi_log <- ggplot(model_outputs_df) +
  geom_point(aes(x = log(Prec_psi_HRS), y = log(Prec_psi_integrated))) +
  stat_smooth(aes(x = log(Prec_psi_HRS), y = log(Prec_psi_integrated)),method = "lm") +
  geom_abline(intercept = 0, slope = 1,linetype = "dashed") +
  coord_equal()

Prec_dtype1 <- ggplot(model_outputs_df) +
  geom_point(aes(x = Prec_dtype1_HRS, y = Prec_dtype1_integrated)) +
  stat_smooth(aes(x = Prec_dtype1_HRS, y = Prec_dtype1_integrated),method = "lm") +
  geom_abline(intercept = 0, slope = 1,linetype = "dashed")

Prec_dtype2 <- ggplot(model_outputs_df) +
  geom_point(aes(x = Prec_dtype2_HRS, y = Prec_dtype2_integrated)) +
  stat_smooth(aes(x = Prec_dtype2_HRS, y = Prec_dtype2_integrated),method = "lm") +
  geom_abline(intercept = 0, slope = 1,linetype = "dashed")

Prec_dtype3 <- ggplot(model_outputs_df) +
  geom_point(aes(x = Prec_dtype3_HRS, y = Prec_dtype3_integrated)) +
  stat_smooth(aes(x = Prec_dtype3_HRS, y = Prec_dtype3_integrated),method = "lm") +
  geom_abline(intercept = 0, slope = 1,linetype = "dashed")

Prec_beta1 <- ggplot(model_outputs_df) +
  geom_point(aes(x = Prec_beta1_HRS, y = Prec_beta1_integrated)) +
  stat_smooth(aes(x = Prec_beta1_HRS, y = Prec_beta1_integrated),method = "lm") +
  geom_abline(intercept = 0, slope = 1,linetype = "dashed")

Prec_beta1_log <- ggplot(model_outputs_df) +
  geom_point(aes(x = log(Prec_beta1_HRS), y = log(Prec_beta1_integrated))) +
  stat_smooth(aes(x = log(Prec_beta1_HRS), y = log(Prec_beta1_integrated)),method = "lm") +
  geom_abline(intercept = 0, slope = 1,linetype = "dashed")

Prec_beta2 <- ggplot(model_outputs_df) +
  geom_point(aes(x = Prec_beta2_HRS, y = Prec_beta2_integrated)) +
  stat_smooth(aes(x = Prec_beta2_HRS, y = Prec_beta2_integrated),method = "lm") +
  geom_abline(intercept = 0, slope = 1,linetype = "dashed")

Prec_beta2_log <- ggplot(model_outputs_df) +
  geom_point(aes(x = log(Prec_beta2_HRS), y = log(Prec_beta2_integrated))) +
  stat_smooth(aes(x = log(Prec_beta2_HRS), y = log(Prec_beta2_integrated)),method = "lm") +
  geom_abline(intercept = 0, slope = 1,linetype = "dashed")

Prec_beta3 <- ggplot(model_outputs_df) +
  geom_point(aes(x = Prec_beta3_HRS, y = Prec_beta3_integrated)) +
  stat_smooth(aes(x = Prec_beta3_HRS, y = Prec_beta3_integrated),method = "lm") +
  geom_abline(intercept = 0, slope = 1,linetype = "dashed")

Prec_beta3_log <- ggplot(model_outputs_df) +
  geom_point(aes(x = log(Prec_beta3_HRS), y = log(Prec_beta3_integrated))) +
  stat_smooth(aes(x = log(Prec_beta3_HRS), y = log(Prec_beta3_integrated)),method = "lm") +
  geom_abline(intercept = 0, slope = 1,linetype = "dashed")

source("../Descriptive Stats/multiplot.R")

conv_psi_HRS <- ggplot(model_outputs_df) +
  geom_histogram(aes(Rhat_psi_HRS, fill = Rhat_psi_HRS > 1.1), 
                 binwidth = 0.05) +
  xlim(0,5)

conv_psi_integrated <- ggplot(model_outputs_df) +
  geom_histogram(aes(Rhat_psi_integrated, fill = Rhat_psi_integrated > 1.1), 
                 binwidth = 0.05) +
  xlim(0,5)

multiplot(conv_psi_integrated, conv_psi_HRS)

conv_a_HRS <- ggplot(model_outputs_df) +
  geom_histogram(aes(Rhat_a_HRS, fill = Rhat_a_HRS > 1.1), 
                 binwidth = 0.05) 

conv_a_integrated <- ggplot(model_outputs_df) +
  geom_histogram(aes(Rhat_a_integrated, fill = Rhat_a_integrated > 1.1), 
                 binwidth = 0.05) 

multiplot(conv_a_HRS, conv_a_integrated)

convergence_psi <- ggplot(model_outputs_df) +
  geom_point(aes(x = Rhat_psi_HRS, y = Rhat_psi_integrated)) +
  stat_smooth(aes(x = Rhat_psi_HRS, y = Rhat_psi_integrated),method = "lm") +
  geom_abline(intercept = 0, slope = 1,linetype = "dashed")

convergence_a <- ggplot(model_outputs_df) +
  geom_point(aes(x = Rhat_a_HRS, y = Rhat_a_integrated)) +
  stat_smooth(aes(x = Rhat_a_HRS, y = Rhat_a_integrated),method = "lm") +
  geom_abline(intercept = 0, slope = 1,linetype = "dashed") 


Prec_a_converged <- ggplot(converged_models) +
  geom_point(aes(x = Prec_a_HRS, y = Prec_a_integrated)) +
  stat_smooth(aes(x = Prec_a_HRS, y = Prec_a_integrated),method = "lm", colour = "black") +
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

# Look at the species for which there is the biggest advantage in precision for a

model_outputs_df[which(model_outputs_df$Prec_a_integrated>40), "species"]
# [1] Eupeodes luniger

load("./Output/Rerun/Eupeodes luniger_HRS_model_it161000.rdata")

load("./Output/Rerun/Eupeodes luniger_integrated_model_it161000.rdata")


out_integrated_model$BUGSoutput$summary

# # check convergence visually

traceplot(out_integrated_model, varname = c("a", "pan.p", "dtype1.p", "dtype2.p", "dtype3.p",
                                            "psi.fs", "beta1", "beta2", "beta3"),
          mfrow = c(3,3), ask = FALSE)


# # plot posterior distribution of occupancy and compare to
# # prop of sites in which species was detected in the 2 datasets

par(mfrow = c(1,1))
hist(out_integrated_model$BUGSoutput$sims.list$psi.fs, xlim = c(0,1),
     col = "gray", main = "", xlab = "Occupancy")

abline(v = sum(aggregate(as.integer(HRS_formatted$spp_vis$'Eupeodes luniger'),
                         by = list(HRS_formatted$occDetdata$site), max)[,2])/length(unique(Site1)),
       lwd = 2, col = "red")
abline(v = sum(aggregate(as.integer(PoMS_formatted$spp_vis$'Eupeodes luniger'),
                         by = list(PoMS_formatted$occDetdata$site), max)[,2])/length(unique(Site2)),
       lwd = 2, col = "blue")
legend("topright", legend = c("Reporting rate HRS", "Reporting rate PoMS"), col = c("red", "blue"), lty = 1, lwd = 2)



out_HRS_model$BUGSoutput$summary

# # check convergence visually

traceplot(out_HRS_model, varname = c("a", "dtype1.p", "dtype2.p", "dtype3.p", "psi.fs",
                                     "beta1", "beta2", "beta3"),
          mfrow = c(3,3), ask = FALSE)


# # plot posterior distribution of occupancy and compare to
# # prop of sites in which species was detected in the 2 datasets

hist(out_HRS_model$BUGSoutput$sims.list$psi.fs, xlim = c(0,1),
     col = "gray", main = "", xlab = "Occupancy")
abline(v = sum(aggregate(as.integer(HRS_formatted$spp_vis$'Eupeodes luniger'),
                         by = list(HRS_formatted$occDetdata$site), max)[,2])/length(unique(Site1)),
       lwd = 2, col = "red")
legend("topleft", legend = c("Reporting rate HRS"), col = c("red"), lty = 1, lwd = 2)






# # Compare the precison of the occupancy

var.psi.hrs <- var(out_HRS_model$BUGSoutput$sims.list$psi.fs)
var.psi.int <- var(out_integrated_model$BUGSoutput$sims.list$psi.fs)

prec.psi.hrs <- 1/var.psi.hrs
prec.psi.int <- 1/var.psi.int

var.a.hrs <- var(out_HRS_model$BUGSoutput$sims.list$a)
var.a.int <- var(out_integrated_model$BUGSoutput$sims.list$a)

prec.a.hrs <- 1/var.a.hrs
prec.a.int <- 1/var.a.int

precisions <- data.frame(HRS_only = c(prec.psi.hrs,prec.a.hrs),
                         integrated = c(prec.psi.int, prec.a.int),
                         row.names = c("precision of psi", "precision of a"))
#
# write.csv(precisions, "./Output/Precisions.csv", row.names = T)
#
width.CI.psi.hrs <- out_HRS_model$BUGSoutput$summary["psi.fs", "97.5%"] -
  out_HRS_model$BUGSoutput$summary["psi.fs", "2.5%"]
width.CI.psi.int <- out_integrated_model$BUGSoutput$summary["psi.fs", "97.5%"] -
  out_integrated_model$BUGSoutput$summary["psi.fs", "2.5%"]
width.CI.a.hrs <- out_HRS_model$BUGSoutput$summary["a", "97.5%"] -
  out_HRS_model$BUGSoutput$summary["a", "2.5%"]
width.CI.a.int <- out_integrated_model$BUGSoutput$summary["a", "97.5%"] -
  out_integrated_model$BUGSoutput$summary["a", "2.5%"]

width.CI <- data.frame(HRS_only = c(width.CI.psi.hrs, width.CI.a.hrs),
                       integrated = c(width.CI.psi.int, width.CI.a.int),
                       row.names = c("width of 95CI for psi", "width of 95CI for a"))

# write.csv(width.CI, "./Output/Width.CI.csv", row.names = T)
#


# Look at the species for which there is the biggest disadvantage 

model_outputs_df[which(model_outputs_df$Prec_a_HRS>10 & model_outputs_df$Prec_a_HRS<20), "species"]
# [1] Eupeodes corollae

load("./Output/Rerun/Eupeodes corollae_HRS_model_it161000.rdata")

load("./Output/Rerun/Eupeodes corollae_integrated_model_it161000.rdata")


out_integrated_model$BUGSoutput$summary

# # check convergence visually

traceplot(out_integrated_model, varname = c("a", "pan.p", "dtype1.p", "dtype2.p", "dtype3.p",
                                            "psi.fs", "beta1", "beta2", "beta3"),
          mfrow = c(3,3), ask = FALSE)


# # plot posterior distribution of occupancy and compare to
# # prop of sites in which species was detected in the 2 datasets


hist(out_integrated_model$BUGSoutput$sims.list$psi.fs, xlim = c(0,1),
     col = "gray", main = "", xlab = "Occupancy")

abline(v = sum(aggregate(as.integer(HRS_formatted$spp_vis$'Eupeodes corollae'),
                         by = list(HRS_formatted$occDetdata$site), max)[,2])/length(unique(Site1)),
       lwd = 2, col = "red")
abline(v = sum(aggregate(as.integer(PoMS_formatted$spp_vis$'Eupeodes corollae'),
                         by = list(PoMS_formatted$occDetdata$site), max)[,2])/length(unique(Site2)),
       lwd = 2, col = "blue")
legend("topright", legend = c("Reporting rate HRS", "Reporting rate PoMS"), col = c("red", "blue"), lty = 1, lwd = 2)



out_HRS_model$BUGSoutput$summary



# Look at the species for which there is no difference in precision for a 

model_outputs_df[which(model_outputs_df$Prec_a_HRS>10 & model_outputs_df$Prec_a_HRS<20), "species"]
# [1] Eupeodes corollae

load("./Output/Rerun/Eupeodes corollae_HRS_model_it161000.rdata")

load("./Output/Rerun/Eupeodes corollae_integrated_model_it161000.rdata")


out_integrated_model$BUGSoutput$summary

# # check convergence visually

traceplot(out_integrated_model, varname = c("a", "pan.p", "dtype1.p", "dtype2.p", "dtype3.p",
                                            "psi.fs", "beta1", "beta2", "beta3"),
          mfrow = c(3,3), ask = FALSE)


# # plot posterior distribution of occupancy and compare to
# # prop of sites in which species was detected in the 2 datasets

par(mfrow = c(1,1))
hist(out_integrated_model$BUGSoutput$sims.list$psi.fs, xlim = c(0,1),
     col = "gray", main = "", xlab = "Occupancy")

abline(v = sum(aggregate(as.integer(HRS_formatted$spp_vis$'Eupeodes corollae'),
                         by = list(HRS_formatted$occDetdata$site), max)[,2])/length(unique(Site1)),
       lwd = 2, col = "red")
abline(v = sum(aggregate(as.integer(PoMS_formatted$spp_vis$'Eupeodes corollae'),
                         by = list(PoMS_formatted$occDetdata$site), max)[,2])/length(unique(Site2)),
       lwd = 2, col = "blue")
legend("topright", legend = c("Reporting rate HRS", "Reporting rate PoMS"), col = c("red", "blue"), lty = 1, lwd = 2)



out_HRS_model$BUGSoutput$summary

traceplot(out_HRS_model, varname = c("a", "dtype1.p", "dtype2.p", "dtype3.p",
                                     "psi.fs", "beta1", "beta2", "beta3"),
          mfrow = c(3,3), ask = FALSE)

par(mfrow = c(1,1))
hist(out_HRS_model$BUGSoutput$sims.list$psi.fs, xlim = c(0,1),
     col = "gray", main = "", xlab = "Occupancy")

abline(v = sum(aggregate(as.integer(HRS_formatted$spp_vis$'Eupeodes corollae'),
                         by = list(HRS_formatted$occDetdata$site), max)[,2])/length(unique(Site1)),
       lwd = 2, col = "red")




# need to find one that converged

model_outputs_df[which(model_outputs_df$Rhat_a_integrated <= 1.1 & 
                         abs(model_outputs_df$Prec_a_HRS - model_outputs_df$Prec_a_integrated) < 2), "species"]
# [1] Eristalis pertinax    Helophilus pendulus   Syritta pipiens       Meliscaeva auricollis Sphaerophoria scripta

load("./Output/Rerun/Helophilus pendulus_HRS_model_it161000.rdata")

load("./Output/Rerun/Helophilus pendulus_integrated_model_it161000.rdata")


out_integrated_model$BUGSoutput$summary

# # check convergence visually

traceplot(out_integrated_model, varname = c("a", "pan.p", "dtype1.p", "dtype2.p", "dtype3.p",
                                            "psi.fs", "beta1", "beta2", "beta3"),
          mfrow = c(3,3), ask = FALSE)


# # plot posterior distribution of occupancy and compare to
# # prop of sites in which species was detected in the 2 datasets

par(mfrow = c(1,1))
hist(out_integrated_model$BUGSoutput$sims.list$psi.fs, xlim = c(0,1),
     col = "gray", main = "", xlab = "Occupancy")

abline(v = sum(aggregate(as.integer(HRS_formatted$spp_vis$'Meliscaeva auricollis'),
                         by = list(HRS_formatted$occDetdata$site), max)[,2])/length(unique(Site1)),
       lwd = 2, col = "red")
abline(v = sum(aggregate(as.integer(PoMS_formatted$spp_vis$'Meliscaeva auricollis'),
                         by = list(PoMS_formatted$occDetdata$site), max)[,2])/length(unique(Site2)),
       lwd = 2, col = "blue")
legend("topright", legend = c("Reporting rate HRS", "Reporting rate PoMS"), col = c("red", "blue"), lty = 1, lwd = 2)


out_HRS_model$BUGSoutput$summary

# # check convergence visually

traceplot(out_HRS_model, varname = c("a", "dtype1.p", "dtype2.p", "dtype3.p", "psi.fs",
                                     "beta1", "beta2", "beta3"),
          mfrow = c(3,3), ask = FALSE)


# # plot posterior distribution of occupancy and compare to
# # prop of sites in which species was detected in the 2 datasets
par(mfrow = c(1,1))
hist(out_HRS_model$BUGSoutput$sims.list$psi.fs, xlim = c(0,1),
     col = "gray", main = "", xlab = "Occupancy")
abline(v = sum(aggregate(as.integer(HRS_formatted$spp_vis$'Sphaerophoria scripta'),
                         by = list(HRS_formatted$occDetdata$site), max)[,2])/length(unique(Site1)),
       lwd = 2, col = "red")
legend("topleft", legend = c("Reporting rate HRS"), col = c("red"), lty = 1, lwd = 2)


##  Episyrphus balteatus ----
load(paste0("./Output/Rerun/", best_species$x[1], "_HRS_model_it161000.rdata"))

load(paste0("./Output/Rerun/", best_species$x[1], "_integrated_model_it161000.rdata"))


out_integrated_model$BUGSoutput$summary

# # check convergence visually

traceplot(out_integrated_model, varname = c("a", "pan.p", "dtype1.p", "dtype2.p", "dtype3.p",
                                            "psi.fs", "beta1", "beta2", "beta3"),
          mfrow = c(3,3), ask = FALSE)


# # plot posterior distribution of occupancy and compare to
# # prop of sites in which species was detected in the 2 datasets

par(mfrow = c(1,1))
hist(out_integrated_model$BUGSoutput$sims.list$psi.fs, xlim = c(0,1),
     col = "gray", main = "", xlab = "Occupancy")

abline(v = sum(aggregate(as.integer(HRS_formatted$spp_vis$'Eupeodes luniger'),
                         by = list(HRS_formatted$occDetdata$site), max)[,2])/length(unique(Site1)),
       lwd = 2, col = "red")
abline(v = sum(aggregate(as.integer(PoMS_formatted$spp_vis$'Eupeodes luniger'),
                         by = list(PoMS_formatted$occDetdata$site), max)[,2])/length(unique(Site2)),
       lwd = 2, col = "blue")
legend("topright", legend = c("Reporting rate HRS", "Reporting rate PoMS"), col = c("red", "blue"), lty = 1, lwd = 2)



out_HRS_model$BUGSoutput$summary

# # check convergence visually

traceplot(out_HRS_model, varname = c("a", "dtype1.p", "dtype2.p", "dtype3.p", "psi.fs",
                                     "beta1", "beta2", "beta3"),
          mfrow = c(3,3), ask = FALSE)


# # plot posterior distribution of occupancy and compare to
# # prop of sites in which species was detected in the 2 datasets

hist(out_HRS_model$BUGSoutput$sims.list$psi.fs, xlim = c(0,1),
     col = "gray", main = "", xlab = "Occupancy")
abline(v = sum(aggregate(as.integer(HRS_formatted$spp_vis$'Eupeodes luniger'),
                         by = list(HRS_formatted$occDetdata$site), max)[,2])/length(unique(Site1)),
       lwd = 2, col = "red")
legend("topleft", legend = c("Reporting rate HRS"), col = c("red"), lty = 1, lwd = 2)






# # Compare the precison of the occupancy

var.psi.hrs <- var(out_HRS_model$BUGSoutput$sims.list$psi.fs)
var.psi.int <- var(out_integrated_model$BUGSoutput$sims.list$psi.fs)

prec.psi.hrs <- 1/var.psi.hrs
prec.psi.int <- 1/var.psi.int

var.a.hrs <- var(out_HRS_model$BUGSoutput$sims.list$a)
var.a.int <- var(out_integrated_model$BUGSoutput$sims.list$a)

prec.a.hrs <- 1/var.a.hrs
prec.a.int <- 1/var.a.int

precisions <- data.frame(HRS_only = c(prec.psi.hrs,prec.a.hrs),
                         integrated = c(prec.psi.int, prec.a.int),
                         row.names = c("precision of psi", "precision of a"))
#
# write.csv(precisions, "./Output/Precisions.csv", row.names = T)
#
width.CI.psi.hrs <- out_HRS_model$BUGSoutput$summary["psi.fs", "97.5%"] -
  out_HRS_model$BUGSoutput$summary["psi.fs", "2.5%"]
width.CI.psi.int <- out_integrated_model$BUGSoutput$summary["psi.fs", "97.5%"] -
  out_integrated_model$BUGSoutput$summary["psi.fs", "2.5%"]
width.CI.a.hrs <- out_HRS_model$BUGSoutput$summary["a", "97.5%"] -
  out_HRS_model$BUGSoutput$summary["a", "2.5%"]
width.CI.a.int <- out_integrated_model$BUGSoutput$summary["a", "97.5%"] -
  out_integrated_model$BUGSoutput$summary["a", "2.5%"]

width.CI <- data.frame(HRS_only = c(width.CI.psi.hrs, width.CI.a.hrs),
                       integrated = c(width.CI.psi.int, width.CI.a.int),
                       row.names = c("width of 95CI for psi", "width of 95CI for a"))

# write.csv(width.CI, "./Output/Width.CI.csv", row.names = T)
#

# Eristalis tenax ----
load(paste0("./Output/Rerun/", best_species$x[2], "_HRS_model_it161000.rdata"))

load(paste0("./Output/Rerun/", best_species$x[2], "_integrated_model_it161000.rdata"))


out_integrated_model$BUGSoutput$summary

# # check convergence visually

traceplot(out_integrated_model, varname = c("a", "pan.p", "dtype1.p", "dtype2.p", "dtype3.p",
                                            "psi.fs", "beta1", "beta2", "beta3"),
          mfrow = c(3,3), ask = FALSE)


# # plot posterior distribution of occupancy and compare to
# # prop of sites in which species was detected in the 2 datasets

par(mfrow = c(1,1))
hist(out_integrated_model$BUGSoutput$sims.list$psi.fs, xlim = c(0,1),
     col = "gray", main = "", xlab = "Occupancy")

abline(v = sum(aggregate(as.integer(HRS_formatted$spp_vis$'Eristalis tenax'),
                         by = list(HRS_formatted$occDetdata$site), max)[,2])/length(unique(Site1)),
       lwd = 2, col = "red")
abline(v = sum(aggregate(as.integer(PoMS_formatted$spp_vis$'Eristalis tenax'),
                         by = list(PoMS_formatted$occDetdata$site), max)[,2])/length(unique(Site2)),
       lwd = 2, col = "blue")
legend("topright", legend = c("Reporting rate HRS", "Reporting rate PoMS"), col = c("red", "blue"), lty = 1, lwd = 2)



out_HRS_model$BUGSoutput$summary

# # check convergence visually

traceplot(out_HRS_model, varname = c("a", "dtype1.p", "dtype2.p", "dtype3.p", "psi.fs",
                                     "beta1", "beta2", "beta3"),
          mfrow = c(3,3), ask = FALSE)


# # plot posterior distribution of occupancy and compare to
# # prop of sites in which species was detected in the 2 datasets

hist(out_HRS_model$BUGSoutput$sims.list$psi.fs, xlim = c(0,1),
     col = "gray", main = "", xlab = "Occupancy")
abline(v = sum(aggregate(as.integer(HRS_formatted$spp_vis$'Eupeodes luniger'),
                         by = list(HRS_formatted$occDetdata$site), max)[,2])/length(unique(Site1)),
       lwd = 2, col = "red")
legend("topleft", legend = c("Reporting rate HRS"), col = c("red"), lty = 1, lwd = 2)






# # Compare the precison of the occupancy

var.psi.hrs <- var(out_HRS_model$BUGSoutput$sims.list$psi.fs)
var.psi.int <- var(out_integrated_model$BUGSoutput$sims.list$psi.fs)

prec.psi.hrs <- 1/var.psi.hrs
prec.psi.int <- 1/var.psi.int

var.a.hrs <- var(out_HRS_model$BUGSoutput$sims.list$a)
var.a.int <- var(out_integrated_model$BUGSoutput$sims.list$a)

prec.a.hrs <- 1/var.a.hrs
prec.a.int <- 1/var.a.int

precisions <- data.frame(HRS_only = c(prec.psi.hrs,prec.a.hrs),
                         integrated = c(prec.psi.int, prec.a.int),
                         row.names = c("precision of psi", "precision of a"))


# Eristalis pertinax ----
load(paste0("./Output/Rerun/", best_species$x[3], "_HRS_model_it161000.rdata"))

load(paste0("./Output/Rerun/", best_species$x[3], "_integrated_model_it161000.rdata"))


out_integrated_model$BUGSoutput$summary

# # check convergence visually

traceplot(out_integrated_model, varname = c("a", "pan.p", "dtype1.p", "dtype2.p", "dtype3.p",
                                            "psi.fs", "beta1", "beta2", "beta3"),
          mfrow = c(3,3), ask = FALSE)


# # plot posterior distribution of occupancy and compare to
# # prop of sites in which species was detected in the 2 datasets

par(mfrow = c(1,1))
hist(out_integrated_model$BUGSoutput$sims.list$psi.fs, xlim = c(0,1),
     col = "gray", main = "", xlab = "Occupancy")

abline(v = sum(aggregate(as.integer(HRS_formatted$spp_vis$'Eristalis pertinax'),
                         by = list(HRS_formatted$occDetdata$site), max)[,2])/length(unique(Site1)),
       lwd = 2, col = "red")
abline(v = sum(aggregate(as.integer(PoMS_formatted$spp_vis$'Eristalis pertinax'),
                         by = list(PoMS_formatted$occDetdata$site), max)[,2])/length(unique(Site2)),
       lwd = 2, col = "blue")
legend("topright", legend = c("Reporting rate HRS", "Reporting rate PoMS"), col = c("red", "blue"), lty = 1, lwd = 2)



out_HRS_model$BUGSoutput$summary

# # check convergence visually

traceplot(out_HRS_model, varname = c("a", "dtype1.p", "dtype2.p", "dtype3.p", "psi.fs",
                                     "beta1", "beta2", "beta3"),
          mfrow = c(3,3), ask = FALSE)


# # plot posterior distribution of occupancy and compare to
# # prop of sites in which species was detected in the 2 datasets

hist(out_HRS_model$BUGSoutput$sims.list$psi.fs, xlim = c(0,1),
     col = "gray", main = "", xlab = "Occupancy")
abline(v = sum(aggregate(as.integer(HRS_formatted$spp_vis$'Eristalis pertinax'),
                         by = list(HRS_formatted$occDetdata$site), max)[,2])/length(unique(Site1)),
       lwd = 2, col = "red")
legend("topleft", legend = c("Reporting rate HRS"), col = c("red"), lty = 1, lwd = 2)






# # Compare the precison of the occupancy

var.psi.hrs <- var(out_HRS_model$BUGSoutput$sims.list$psi.fs)
var.psi.int <- var(out_integrated_model$BUGSoutput$sims.list$psi.fs)

prec.psi.hrs <- 1/var.psi.hrs
prec.psi.int <- 1/var.psi.int

var.a.hrs <- var(out_HRS_model$BUGSoutput$sims.list$a)
var.a.int <- var(out_integrated_model$BUGSoutput$sims.list$a)

prec.a.hrs <- 1/var.a.hrs
prec.a.int <- 1/var.a.int

precisions <- data.frame(HRS_only = c(prec.psi.hrs,prec.a.hrs),
                         integrated = c(prec.psi.int, prec.a.int),
                         row.names = c("precision of psi", "precision of a"))


# Melanostoma scalare ----
load(paste0("./Output/Rerun/", best_species$x[4], "_HRS_model_it161000.rdata"))

load(paste0("./Output/Rerun/", best_species$x[4], "_integrated_model_it161000.rdata"))


out_integrated_model$BUGSoutput$summary

# # check convergence visually

traceplot(out_integrated_model, varname = c("a", "pan.p", "dtype1.p", "dtype2.p", "dtype3.p",
                                            "psi.fs", "beta1", "beta2", "beta3"),
          mfrow = c(3,3), ask = FALSE)


# # plot posterior distribution of occupancy and compare to
# # prop of sites in which species was detected in the 2 datasets

par(mfrow = c(1,1))
hist(out_integrated_model$BUGSoutput$sims.list$psi.fs, xlim = c(0,1),
     col = "gray", main = "", xlab = "Occupancy")

abline(v = sum(aggregate(as.integer(HRS_formatted$spp_vis$'Melanostoma scalare'),
                         by = list(HRS_formatted$occDetdata$site), max)[,2])/length(unique(Site1)),
       lwd = 2, col = "red")
abline(v = sum(aggregate(as.integer(PoMS_formatted$spp_vis$'Melanostoma scalare'),
                         by = list(PoMS_formatted$occDetdata$site), max)[,2])/length(unique(Site2)),
       lwd = 2, col = "blue")
legend("topright", legend = c("Reporting rate HRS", "Reporting rate PoMS"), col = c("red", "blue"), lty = 1, lwd = 2)



out_HRS_model$BUGSoutput$summary

# # check convergence visually

traceplot(out_HRS_model, varname = c("a", "dtype1.p", "dtype2.p", "dtype3.p", "psi.fs",
                                     "beta1", "beta2", "beta3"),
          mfrow = c(3,3), ask = FALSE)


# # plot posterior distribution of occupancy and compare to
# # prop of sites in which species was detected in the 2 datasets

hist(out_HRS_model$BUGSoutput$sims.list$psi.fs, xlim = c(0,1),
     col = "gray", main = "", xlab = "Occupancy")
abline(v = sum(aggregate(as.integer(HRS_formatted$spp_vis$'Melanostoma scalare'),
                         by = list(HRS_formatted$occDetdata$site), max)[,2])/length(unique(Site1)),
       lwd = 2, col = "red")
legend("topleft", legend = c("Reporting rate HRS"), col = c("red"), lty = 1, lwd = 2)






# # Compare the precison of the occupancy

var.psi.hrs <- var(out_HRS_model$BUGSoutput$sims.list$psi.fs)
var.psi.int <- var(out_integrated_model$BUGSoutput$sims.list$psi.fs)

prec.psi.hrs <- 1/var.psi.hrs
prec.psi.int <- 1/var.psi.int

var.a.hrs <- var(out_HRS_model$BUGSoutput$sims.list$a)
var.a.int <- var(out_integrated_model$BUGSoutput$sims.list$a)

prec.a.hrs <- 1/var.a.hrs
prec.a.int <- 1/var.a.int

precisions <- data.frame(HRS_only = c(prec.psi.hrs,prec.a.hrs),
                         integrated = c(prec.psi.int, prec.a.int),
                         row.names = c("precision of psi", "precision of a"))


# Platycheirus albimanus ----
load(paste0("./Output/Rerun/", best_species$x[5], "_HRS_model_it161000.rdata"))

load(paste0("./Output/Rerun/", best_species$x[5], "_integrated_model_it161000.rdata"))


out_integrated_model$BUGSoutput$summary

# # check convergence visually

traceplot(out_integrated_model, varname = c("a", "pan.p", "dtype1.p", "dtype2.p", "dtype3.p",
                                            "psi.fs", "beta1", "beta2", "beta3"),
          mfrow = c(3,3), ask = FALSE)


# # plot posterior distribution of occupancy and compare to
# # prop of sites in which species was detected in the 2 datasets

par(mfrow = c(1,1))
hist(out_integrated_model$BUGSoutput$sims.list$psi.fs, xlim = c(0,1),
     col = "gray", main = "", xlab = "Occupancy")

abline(v = sum(aggregate(as.integer(HRS_formatted$spp_vis$'Platycheirus albimanus'),
                         by = list(HRS_formatted$occDetdata$site), max)[,2])/length(unique(Site1)),
       lwd = 2, col = "red")
abline(v = sum(aggregate(as.integer(PoMS_formatted$spp_vis$'Platycheirus albimanus'),
                         by = list(PoMS_formatted$occDetdata$site), max)[,2])/length(unique(Site2)),
       lwd = 2, col = "blue")
legend("topright", legend = c("Reporting rate HRS", "Reporting rate PoMS"), col = c("red", "blue"), lty = 1, lwd = 2)



out_HRS_model$BUGSoutput$summary

# # check convergence visually

traceplot(out_HRS_model, varname = c("a", "dtype1.p", "dtype2.p", "dtype3.p", "psi.fs",
                                     "beta1", "beta2", "beta3"),
          mfrow = c(3,3), ask = FALSE)


# # plot posterior distribution of occupancy and compare to
# # prop of sites in which species was detected in the 2 datasets

hist(out_HRS_model$BUGSoutput$sims.list$psi.fs, xlim = c(0,1),
     col = "gray", main = "", xlab = "Occupancy")
abline(v = sum(aggregate(as.integer(HRS_formatted$spp_vis$'Platycheirus albimanus'),
                         by = list(HRS_formatted$occDetdata$site), max)[,2])/length(unique(Site1)),
       lwd = 2, col = "red")
legend("topleft", legend = c("Reporting rate HRS"), col = c("red"), lty = 1, lwd = 2)






# # Compare the precison of the occupancy

var.psi.hrs <- var(out_HRS_model$BUGSoutput$sims.list$psi.fs)
var.psi.int <- var(out_integrated_model$BUGSoutput$sims.list$psi.fs)

prec.psi.hrs <- 1/var.psi.hrs
prec.psi.int <- 1/var.psi.int

var.a.hrs <- var(out_HRS_model$BUGSoutput$sims.list$a)
var.a.int <- var(out_integrated_model$BUGSoutput$sims.list$a)

prec.a.hrs <- 1/var.a.hrs
prec.a.int <- 1/var.a.int

precisions <- data.frame(HRS_only = c(prec.psi.hrs,prec.a.hrs),
                         integrated = c(prec.psi.int, prec.a.int),
                         row.names = c("precision of psi", "precision of a"))


# Helophilus pendulus ----
load(paste0("./Output/Rerun/", best_species$x[6], "_HRS_model_it161000.rdata"))

load(paste0("./Output/Rerun/", best_species$x[6], "_integrated_model_it161000.rdata"))


out_integrated_model$BUGSoutput$summary

# # check convergence visually

traceplot(out_integrated_model, varname = c("a", "pan.p", "dtype1.p", "dtype2.p", "dtype3.p",
                                            "psi.fs", "beta1", "beta2", "beta3"),
          mfrow = c(3,3), ask = FALSE)


# # plot posterior distribution of occupancy and compare to
# # prop of sites in which species was detected in the 2 datasets

par(mfrow = c(1,1))
hist(out_integrated_model$BUGSoutput$sims.list$psi.fs, xlim = c(0,1),
     col = "gray", main = "", xlab = "Occupancy")

abline(v = sum(aggregate(as.integer(HRS_formatted$spp_vis$'Helophilus pendulus'),
                         by = list(HRS_formatted$occDetdata$site), max)[,2])/length(unique(Site1)),
       lwd = 2, col = "red")
abline(v = sum(aggregate(as.integer(PoMS_formatted$spp_vis$'Helophilus pendulus'),
                         by = list(PoMS_formatted$occDetdata$site), max)[,2])/length(unique(Site2)),
       lwd = 2, col = "blue")
legend("topright", legend = c("Reporting rate HRS", "Reporting rate PoMS"), col = c("red", "blue"), lty = 1, lwd = 2)



out_HRS_model$BUGSoutput$summary

# # check convergence visually

traceplot(out_HRS_model, varname = c("a", "dtype1.p", "dtype2.p", "dtype3.p", "psi.fs",
                                     "beta1", "beta2", "beta3"),
          mfrow = c(3,3), ask = FALSE)


# # plot posterior distribution of occupancy and compare to
# # prop of sites in which species was detected in the 2 datasets
par(mfrow = c(1,1))
hist(out_HRS_model$BUGSoutput$sims.list$psi.fs, xlim = c(0,1),
     col = "gray", main = "", xlab = "Occupancy")
abline(v = sum(aggregate(as.integer(HRS_formatted$spp_vis$'Helophilus pendulus'),
                         by = list(HRS_formatted$occDetdata$site), max)[,2])/length(unique(Site1)),
       lwd = 2, col = "red")
legend("topleft", legend = c("Reporting rate HRS"), col = c("red"), lty = 1, lwd = 2)






# # Compare the precison of the occupancy

var.psi.hrs <- var(out_HRS_model$BUGSoutput$sims.list$psi.fs)
var.psi.int <- var(out_integrated_model$BUGSoutput$sims.list$psi.fs)

prec.psi.hrs <- 1/var.psi.hrs
prec.psi.int <- 1/var.psi.int

var.a.hrs <- var(out_HRS_model$BUGSoutput$sims.list$a)
var.a.int <- var(out_integrated_model$BUGSoutput$sims.list$a)

prec.a.hrs <- 1/var.a.hrs
prec.a.int <- 1/var.a.int

precisions <- data.frame(HRS_only = c(prec.psi.hrs,prec.a.hrs),
                         integrated = c(prec.psi.int, prec.a.int),
                         row.names = c("precision of psi", "precision of a"))


# Myathropa florea ----
load(paste0("./Output/Rerun/", best_species$x[7], "_HRS_model_it161000.rdata"))

load(paste0("./Output/Rerun/", best_species$x[7], "_integrated_model_it161000.rdata"))


out_integrated_model$BUGSoutput$summary

# # check convergence visually

traceplot(out_integrated_model, varname = c("a", "pan.p", "dtype1.p", "dtype2.p", "dtype3.p",
                                            "psi.fs", "beta1", "beta2", "beta3"),
          mfrow = c(3,3), ask = FALSE)


# # plot posterior distribution of occupancy and compare to
# # prop of sites in which species was detected in the 2 datasets

par(mfrow = c(1,1))
hist(out_integrated_model$BUGSoutput$sims.list$psi.fs, xlim = c(0,1),
     col = "gray", main = "", xlab = "Occupancy")

abline(v = sum(aggregate(as.integer(HRS_formatted$spp_vis$'Myathropa florea'),
                         by = list(HRS_formatted$occDetdata$site), max)[,2])/length(unique(Site1)),
       lwd = 2, col = "red")
abline(v = sum(aggregate(as.integer(PoMS_formatted$spp_vis$'Myathropa florea'),
                         by = list(PoMS_formatted$occDetdata$site), max)[,2])/length(unique(Site2)),
       lwd = 2, col = "blue")
legend("topright", legend = c("Reporting rate HRS", "Reporting rate PoMS"), col = c("red", "blue"), lty = 1, lwd = 2)



out_HRS_model$BUGSoutput$summary

# # check convergence visually

traceplot(out_HRS_model, varname = c("a", "dtype1.p", "dtype2.p", "dtype3.p", "psi.fs",
                                     "beta1", "beta2", "beta3"),
          mfrow = c(3,3), ask = FALSE)


# # plot posterior distribution of occupancy and compare to
# # prop of sites in which species was detected in the 2 datasets
par(mfrow = c(1,1))
hist(out_HRS_model$BUGSoutput$sims.list$psi.fs, xlim = c(0,1),
     col = "gray", main = "", xlab = "Occupancy")
abline(v = sum(aggregate(as.integer(HRS_formatted$spp_vis$'Myathropa florea'),
                         by = list(HRS_formatted$occDetdata$site), max)[,2])/length(unique(Site1)),
       lwd = 2, col = "red")
legend("topleft", legend = c("Reporting rate HRS"), col = c("red"), lty = 1, lwd = 2)






# # Compare the precison of the occupancy

var.psi.hrs <- var(out_HRS_model$BUGSoutput$sims.list$psi.fs)
var.psi.int <- var(out_integrated_model$BUGSoutput$sims.list$psi.fs)

prec.psi.hrs <- 1/var.psi.hrs
prec.psi.int <- 1/var.psi.int

var.a.hrs <- var(out_HRS_model$BUGSoutput$sims.list$a)
var.a.int <- var(out_integrated_model$BUGSoutput$sims.list$a)

prec.a.hrs <- 1/var.a.hrs
prec.a.int <- 1/var.a.int

precisions <- data.frame(HRS_only = c(prec.psi.hrs,prec.a.hrs),
                         integrated = c(prec.psi.int, prec.a.int),
                         row.names = c("precision of psi", "precision of a"))



# Syritta pipiens ----
load(paste0("./Output/Rerun/", best_species$x[8], "_HRS_model_it161000.rdata"))

load(paste0("./Output/Rerun/", best_species$x[8], "_integrated_model_it161000.rdata"))


out_integrated_model$BUGSoutput$summary

# # check convergence visually

traceplot(out_integrated_model, varname = c("a", "pan.p", "dtype1.p", "dtype2.p", "dtype3.p",
                                            "psi.fs", "beta1", "beta2", "beta3"),
          mfrow = c(3,3), ask = FALSE)


# # plot posterior distribution of occupancy and compare to
# # prop of sites in which species was detected in the 2 datasets

par(mfrow = c(1,1))
hist(out_integrated_model$BUGSoutput$sims.list$psi.fs, xlim = c(0,1),
     col = "gray", main = "", xlab = "Occupancy")

abline(v = sum(aggregate(as.integer(HRS_formatted$spp_vis$'Syritta pipiens'),
                         by = list(HRS_formatted$occDetdata$site), max)[,2])/length(unique(Site1)),
       lwd = 2, col = "red")
abline(v = sum(aggregate(as.integer(PoMS_formatted$spp_vis$'Syritta pipiens'),
                         by = list(PoMS_formatted$occDetdata$site), max)[,2])/length(unique(Site2)),
       lwd = 2, col = "blue")
legend("topright", legend = c("Reporting rate HRS", "Reporting rate PoMS"), col = c("red", "blue"), lty = 1, lwd = 2)



out_HRS_model$BUGSoutput$summary

# # check convergence visually

traceplot(out_HRS_model, varname = c("a", "dtype1.p", "dtype2.p", "dtype3.p", "psi.fs",
                                     "beta1", "beta2", "beta3"),
          mfrow = c(3,3), ask = FALSE)


# # plot posterior distribution of occupancy and compare to
# # prop of sites in which species was detected in the 2 datasets
par(mfrow = c(1,1))
hist(out_HRS_model$BUGSoutput$sims.list$psi.fs, xlim = c(0,1),
     col = "gray", main = "", xlab = "Occupancy")
abline(v = sum(aggregate(as.integer(HRS_formatted$spp_vis$'Syritta pipiens'),
                         by = list(HRS_formatted$occDetdata$site), max)[,2])/length(unique(Site1)),
       lwd = 2, col = "red")
legend("topleft", legend = c("Reporting rate HRS"), col = c("red"), lty = 1, lwd = 2)






# # Compare the precison of the occupancy

var.psi.hrs <- var(out_HRS_model$BUGSoutput$sims.list$psi.fs)
var.psi.int <- var(out_integrated_model$BUGSoutput$sims.list$psi.fs)

prec.psi.hrs <- 1/var.psi.hrs
prec.psi.int <- 1/var.psi.int

var.a.hrs <- var(out_HRS_model$BUGSoutput$sims.list$a)
var.a.int <- var(out_integrated_model$BUGSoutput$sims.list$a)

prec.a.hrs <- 1/var.a.hrs
prec.a.int <- 1/var.a.int

precisions <- data.frame(HRS_only = c(prec.psi.hrs,prec.a.hrs),
                         integrated = c(prec.psi.int, prec.a.int),
                         row.names = c("precision of psi", "precision of a"))


# Syrphus ribesii ----
load(paste0("./Output/Rerun/", best_species$x[9], "_HRS_model_it161000.rdata"))

load(paste0("./Output/Rerun/", best_species$x[9], "_integrated_model_it161000.rdata"))


out_integrated_model$BUGSoutput$summary

# # check convergence visually

traceplot(out_integrated_model, varname = c("a", "pan.p", "dtype1.p", "dtype2.p", "dtype3.p",
                                            "psi.fs", "beta1", "beta2", "beta3"),
          mfrow = c(3,3), ask = FALSE)


# # plot posterior distribution of occupancy and compare to
# # prop of sites in which species was detected in the 2 datasets

par(mfrow = c(1,1))
hist(out_integrated_model$BUGSoutput$sims.list$psi.fs, xlim = c(0,1),
     col = "gray", main = "", xlab = "Occupancy")

abline(v = sum(aggregate(as.integer(HRS_formatted$spp_vis$'Syrphus ribesii'),
                         by = list(HRS_formatted$occDetdata$site), max)[,2])/length(unique(Site1)),
       lwd = 2, col = "red")
abline(v = sum(aggregate(as.integer(PoMS_formatted$spp_vis$'Syrphus ribesii'),
                         by = list(PoMS_formatted$occDetdata$site), max)[,2])/length(unique(Site2)),
       lwd = 2, col = "blue")
legend("topright", legend = c("Reporting rate HRS", "Reporting rate PoMS"), col = c("red", "blue"), lty = 1, lwd = 2)



out_HRS_model$BUGSoutput$summary

# # check convergence visually

traceplot(out_HRS_model, varname = c("a", "dtype1.p", "dtype2.p", "dtype3.p", "psi.fs",
                                     "beta1", "beta2", "beta3"),
          mfrow = c(3,3), ask = FALSE)


# # plot posterior distribution of occupancy and compare to
# # prop of sites in which species was detected in the 2 datasets
par(mfrow = c(1,1))
hist(out_HRS_model$BUGSoutput$sims.list$psi.fs, xlim = c(0,1),
     col = "gray", main = "", xlab = "Occupancy")
abline(v = sum(aggregate(as.integer(HRS_formatted$spp_vis$'Syrphus ribesii'),
                         by = list(HRS_formatted$occDetdata$site), max)[,2])/length(unique(Site1)),
       lwd = 2, col = "red")
legend("topleft", legend = c("Reporting rate HRS"), col = c("red"), lty = 1, lwd = 2)






# # Compare the precison of the occupancy

var.psi.hrs <- var(out_HRS_model$BUGSoutput$sims.list$psi.fs)
var.psi.int <- var(out_integrated_model$BUGSoutput$sims.list$psi.fs)

prec.psi.hrs <- 1/var.psi.hrs
prec.psi.int <- 1/var.psi.int

var.a.hrs <- var(out_HRS_model$BUGSoutput$sims.list$a)
var.a.int <- var(out_integrated_model$BUGSoutput$sims.list$a)

prec.a.hrs <- 1/var.a.hrs
prec.a.int <- 1/var.a.int

precisions <- data.frame(HRS_only = c(prec.psi.hrs,prec.a.hrs),
                         integrated = c(prec.psi.int, prec.a.int),
                         row.names = c("precision of psi", "precision of a"))


# Xylota segnis ----
load(paste0("./Output/Rerun/", best_species$x[11], "_HRS_model_it161000.rdata"))

load(paste0("./Output/Rerun/", best_species$x[11], "_integrated_model_it161000.rdata"))


out_integrated_model$BUGSoutput$summary

# # check convergence visually

traceplot(out_integrated_model, varname = c("a", "pan.p", "dtype1.p", "dtype2.p", "dtype3.p",
                                            "psi.fs", "beta1", "beta2", "beta3"),
          mfrow = c(3,3), ask = FALSE)


# # plot posterior distribution of occupancy and compare to
# # prop of sites in which species was detected in the 2 datasets

par(mfrow = c(1,1))
hist(out_integrated_model$BUGSoutput$sims.list$psi.fs, xlim = c(0,1),
     col = "gray", main = "", xlab = "Occupancy")

abline(v = sum(aggregate(as.integer(HRS_formatted$spp_vis$'Xylota segnis'),
                         by = list(HRS_formatted$occDetdata$site), max)[,2])/length(unique(Site1)),
       lwd = 2, col = "red")
abline(v = sum(aggregate(as.integer(PoMS_formatted$spp_vis$'Xylota segnis'),
                         by = list(PoMS_formatted$occDetdata$site), max)[,2])/length(unique(Site2)),
       lwd = 2, col = "blue")
legend("topright", legend = c("Reporting rate HRS", "Reporting rate PoMS"), col = c("red", "blue"), lty = 1, lwd = 2)



out_HRS_model$BUGSoutput$summary

# # check convergence visually

traceplot(out_HRS_model, varname = c("a", "dtype1.p", "dtype2.p", "dtype3.p", "psi.fs",
                                     "beta1", "beta2", "beta3"),
          mfrow = c(3,3), ask = FALSE)


# # plot posterior distribution of occupancy and compare to
# # prop of sites in which species was detected in the 2 datasets
par(mfrow = c(1,1))
hist(out_HRS_model$BUGSoutput$sims.list$psi.fs, xlim = c(0,1),
     col = "gray", main = "", xlab = "Occupancy")
abline(v = sum(aggregate(as.integer(HRS_formatted$spp_vis$'Xylota segnis'),
                         by = list(HRS_formatted$occDetdata$site), max)[,2])/length(unique(Site1)),
       lwd = 2, col = "red")
legend("topleft", legend = c("Reporting rate HRS"), col = c("red"), lty = 1, lwd = 2)






# # Compare the precison of the occupancy

var.psi.hrs <- var(out_HRS_model$BUGSoutput$sims.list$psi.fs)
var.psi.int <- var(out_integrated_model$BUGSoutput$sims.list$psi.fs)

prec.psi.hrs <- 1/var.psi.hrs
prec.psi.int <- 1/var.psi.int

var.a.hrs <- var(out_HRS_model$BUGSoutput$sims.list$a)
var.a.int <- var(out_integrated_model$BUGSoutput$sims.list$a)

prec.a.hrs <- 1/var.a.hrs
prec.a.int <- 1/var.a.int

precisions <- data.frame(HRS_only = c(prec.psi.hrs,prec.a.hrs),
                         integrated = c(prec.psi.int, prec.a.int),
                         row.names = c("precision of psi", "precision of a"))



# Rhingia campestris ----
load(paste0("./Output/Rerun/", best_species$x[12], "_HRS_model_it161000.rdata"))

load(paste0("./Output/Rerun/", best_species$x[12], "_integrated_model_it161000.rdata"))


out_integrated_model$BUGSoutput$summary

# # check convergence visually

traceplot(out_integrated_model, varname = c("a", "pan.p", "dtype1.p", "dtype2.p", "dtype3.p",
                                            "psi.fs", "beta1", "beta2", "beta3"),
          mfrow = c(3,3), ask = FALSE)


# # plot posterior distribution of occupancy and compare to
# # prop of sites in which species was detected in the 2 datasets

par(mfrow = c(1,1))
hist(out_integrated_model$BUGSoutput$sims.list$psi.fs, xlim = c(0,1),
     col = "gray", main = "", xlab = "Occupancy")

abline(v = sum(aggregate(as.integer(HRS_formatted$spp_vis$'Rhingia campestris'),
                         by = list(HRS_formatted$occDetdata$site), max)[,2])/length(unique(Site1)),
       lwd = 2, col = "red")
abline(v = sum(aggregate(as.integer(PoMS_formatted$spp_vis$'Rhingia campestris'),
                         by = list(PoMS_formatted$occDetdata$site), max)[,2])/length(unique(Site2)),
       lwd = 2, col = "blue")
legend("topright", legend = c("Reporting rate HRS", "Reporting rate PoMS"), col = c("red", "blue"), lty = 1, lwd = 2)



out_HRS_model$BUGSoutput$summary

# # check convergence visually

traceplot(out_HRS_model, varname = c("a", "dtype1.p", "dtype2.p", "dtype3.p", "psi.fs",
                                     "beta1", "beta2", "beta3"),
          mfrow = c(3,3), ask = FALSE)


# # plot posterior distribution of occupancy and compare to
# # prop of sites in which species was detected in the 2 datasets
par(mfrow = c(1,1))
hist(out_HRS_model$BUGSoutput$sims.list$psi.fs, xlim = c(0,1),
     col = "gray", main = "", xlab = "Occupancy")
abline(v = sum(aggregate(as.integer(HRS_formatted$spp_vis$'Rhingia campestris'),
                         by = list(HRS_formatted$occDetdata$site), max)[,2])/length(unique(Site1)),
       lwd = 2, col = "red")
legend("topleft", legend = c("Reporting rate HRS"), col = c("red"), lty = 1, lwd = 2)






# # Compare the precison of the occupancy

var.psi.hrs <- var(out_HRS_model$BUGSoutput$sims.list$psi.fs)
var.psi.int <- var(out_integrated_model$BUGSoutput$sims.list$psi.fs)

prec.psi.hrs <- 1/var.psi.hrs
prec.psi.int <- 1/var.psi.int

var.a.hrs <- var(out_HRS_model$BUGSoutput$sims.list$a)
var.a.int <- var(out_integrated_model$BUGSoutput$sims.list$a)

prec.a.hrs <- 1/var.a.hrs
prec.a.int <- 1/var.a.int

precisions <- data.frame(HRS_only = c(prec.psi.hrs,prec.a.hrs),
                         integrated = c(prec.psi.int, prec.a.int),
                         row.names = c("precision of psi", "precision of a"))


# Eupeodes corollae ----
load(paste0("./Output/Rerun/", best_species$x[13], "_HRS_model_it161000.rdata"))

load(paste0("./Output/Rerun/", best_species$x[13], "_integrated_model_it161000.rdata"))


out_integrated_model$BUGSoutput$summary

# # check convergence visually

traceplot(out_integrated_model, varname = c("a", "pan.p", "dtype1.p", "dtype2.p", "dtype3.p",
                                            "psi.fs", "beta1", "beta2", "beta3"),
          mfrow = c(3,3), ask = FALSE)


# # plot posterior distribution of occupancy and compare to
# # prop of sites in which species was detected in the 2 datasets

par(mfrow = c(1,1))
hist(out_integrated_model$BUGSoutput$sims.list$psi.fs, xlim = c(0,1),
     col = "gray", main = "", xlab = "Occupancy")

abline(v = sum(aggregate(as.integer(HRS_formatted$spp_vis$'Eupeodes corollae'),
                         by = list(HRS_formatted$occDetdata$site), max)[,2])/length(unique(Site1)),
       lwd = 2, col = "red")
abline(v = sum(aggregate(as.integer(PoMS_formatted$spp_vis$'Eupeodes corollae'),
                         by = list(PoMS_formatted$occDetdata$site), max)[,2])/length(unique(Site2)),
       lwd = 2, col = "blue")
legend("topright", legend = c("Reporting rate HRS", "Reporting rate PoMS"), col = c("red", "blue"), lty = 1, lwd = 2)



out_HRS_model$BUGSoutput$summary

# # check convergence visually

traceplot(out_HRS_model, varname = c("a", "dtype1.p", "dtype2.p", "dtype3.p", "psi.fs",
                                     "beta1", "beta2", "beta3"),
          mfrow = c(3,3), ask = FALSE)


# # plot posterior distribution of occupancy and compare to
# # prop of sites in which species was detected in the 2 datasets
par(mfrow = c(1,1))
hist(out_HRS_model$BUGSoutput$sims.list$psi.fs, xlim = c(0,1),
     col = "gray", main = "", xlab = "Occupancy")
abline(v = sum(aggregate(as.integer(HRS_formatted$spp_vis$'Eupeodes corollae'),
                         by = list(HRS_formatted$occDetdata$site), max)[,2])/length(unique(Site1)),
       lwd = 2, col = "red")
legend("topleft", legend = c("Reporting rate HRS"), col = c("red"), lty = 1, lwd = 2)






# # Compare the precison of the occupancy

var.psi.hrs <- var(out_HRS_model$BUGSoutput$sims.list$psi.fs)
var.psi.int <- var(out_integrated_model$BUGSoutput$sims.list$psi.fs)

prec.psi.hrs <- 1/var.psi.hrs
prec.psi.int <- 1/var.psi.int

var.a.hrs <- var(out_HRS_model$BUGSoutput$sims.list$a)
var.a.int <- var(out_integrated_model$BUGSoutput$sims.list$a)

prec.a.hrs <- 1/var.a.hrs
prec.a.int <- 1/var.a.int

precisions <- data.frame(HRS_only = c(prec.psi.hrs,prec.a.hrs),
                         integrated = c(prec.psi.int, prec.a.int),
                         row.names = c("precision of psi", "precision of a"))




# Eupeodes luniger ----
load(paste0("./Output/Rerun/", best_species$x[14], "_HRS_model_it161000.rdata"))

load(paste0("./Output/Rerun/", best_species$x[14], "_integrated_model_it161000.rdata"))


out_integrated_model$BUGSoutput$summary

# # check convergence visually

traceplot(out_integrated_model, varname = c("a", "pan.p", "dtype1.p", "dtype2.p", "dtype3.p",
                                            "psi.fs", "beta1", "beta2", "beta3"),
          mfrow = c(3,3), ask = FALSE)


# # plot posterior distribution of occupancy and compare to
# # prop of sites in which species was detected in the 2 datasets

par(mfrow = c(1,1))
hist(out_integrated_model$BUGSoutput$sims.list$psi.fs, xlim = c(0,1),
     col = "gray", main = "", xlab = "Occupancy")

abline(v = sum(aggregate(as.integer(HRS_formatted$spp_vis$'Eupeodes luniger'),
                         by = list(HRS_formatted$occDetdata$site), max)[,2])/length(unique(Site1)),
       lwd = 2, col = "red")
abline(v = sum(aggregate(as.integer(PoMS_formatted$spp_vis$'Eupeodes luniger'),
                         by = list(PoMS_formatted$occDetdata$site), max)[,2])/length(unique(Site2)),
       lwd = 2, col = "blue")
legend("topright", legend = c("Reporting rate HRS", "Reporting rate PoMS"), col = c("red", "blue"), lty = 1, lwd = 2)



out_HRS_model$BUGSoutput$summary

# # check convergence visually

traceplot(out_HRS_model, varname = c("a", "dtype1.p", "dtype2.p", "dtype3.p", "psi.fs",
                                     "beta1", "beta2", "beta3"),
          mfrow = c(3,3), ask = FALSE)


# # plot posterior distribution of occupancy and compare to
# # prop of sites in which species was detected in the 2 datasets
par(mfrow = c(1,1))
hist(out_HRS_model$BUGSoutput$sims.list$psi.fs, xlim = c(0,1),
     col = "gray", main = "", xlab = "Occupancy")
abline(v = sum(aggregate(as.integer(HRS_formatted$spp_vis$'Eupeodes luniger'),
                         by = list(HRS_formatted$occDetdata$site), max)[,2])/length(unique(Site1)),
       lwd = 2, col = "red")
legend("topleft", legend = c("Reporting rate HRS"), col = c("red"), lty = 1, lwd = 2)






# # Compare the precison of the occupancy

var.psi.hrs <- var(out_HRS_model$BUGSoutput$sims.list$psi.fs)
var.psi.int <- var(out_integrated_model$BUGSoutput$sims.list$psi.fs)

prec.psi.hrs <- 1/var.psi.hrs
prec.psi.int <- 1/var.psi.int

var.a.hrs <- var(out_HRS_model$BUGSoutput$sims.list$a)
var.a.int <- var(out_integrated_model$BUGSoutput$sims.list$a)

prec.a.hrs <- 1/var.a.hrs
prec.a.int <- 1/var.a.int

precisions <- data.frame(HRS_only = c(prec.psi.hrs,prec.a.hrs),
                         integrated = c(prec.psi.int, prec.a.int),
                         row.names = c("precision of psi", "precision of a"))



# Eristalis arbustorum ----
load(paste0("./Output/Rerun/", best_species$x[15], "_HRS_model_it161000.rdata"))

load(paste0("./Output/Rerun/", best_species$x[15], "_integrated_model_it161000.rdata"))


out_integrated_model$BUGSoutput$summary

# # check convergence visually

traceplot(out_integrated_model, varname = c("a", "pan.p", "dtype1.p", "dtype2.p", "dtype3.p",
                                            "psi.fs", "beta1", "beta2", "beta3"),
          mfrow = c(3,3), ask = FALSE)


# # plot posterior distribution of occupancy and compare to
# # prop of sites in which species was detected in the 2 datasets

par(mfrow = c(1,1))
hist(out_integrated_model$BUGSoutput$sims.list$psi.fs, xlim = c(0,1),
     col = "gray", main = "", xlab = "Occupancy")

abline(v = sum(aggregate(as.integer(HRS_formatted$spp_vis$'Eristalis arbustorum'),
                         by = list(HRS_formatted$occDetdata$site), max)[,2])/length(unique(Site1)),
       lwd = 2, col = "red")
abline(v = sum(aggregate(as.integer(PoMS_formatted$spp_vis$'Eristalis arbustorum'),
                         by = list(PoMS_formatted$occDetdata$site), max)[,2])/length(unique(Site2)),
       lwd = 2, col = "blue")
legend("topright", legend = c("Reporting rate HRS", "Reporting rate PoMS"), col = c("red", "blue"), lty = 1, lwd = 2)



out_HRS_model$BUGSoutput$summary

# # check convergence visually

traceplot(out_HRS_model, varname = c("a", "dtype1.p", "dtype2.p", "dtype3.p", "psi.fs",
                                     "beta1", "beta2", "beta3"),
          mfrow = c(3,3), ask = FALSE)


# # plot posterior distribution of occupancy and compare to
# # prop of sites in which species was detected in the 2 datasets
par(mfrow = c(1,1))
hist(out_HRS_model$BUGSoutput$sims.list$psi.fs, xlim = c(0,1),
     col = "gray", main = "", xlab = "Occupancy")
abline(v = sum(aggregate(as.integer(HRS_formatted$spp_vis$'Eristalis arbustorum'),
                         by = list(HRS_formatted$occDetdata$site), max)[,2])/length(unique(Site1)),
       lwd = 2, col = "red")
legend("topleft", legend = c("Reporting rate HRS"), col = c("red"), lty = 1, lwd = 2)






# # Compare the precison of the occupancy

var.psi.hrs <- var(out_HRS_model$BUGSoutput$sims.list$psi.fs)
var.psi.int <- var(out_integrated_model$BUGSoutput$sims.list$psi.fs)

prec.psi.hrs <- 1/var.psi.hrs
prec.psi.int <- 1/var.psi.int

var.a.hrs <- var(out_HRS_model$BUGSoutput$sims.list$a)
var.a.int <- var(out_integrated_model$BUGSoutput$sims.list$a)

prec.a.hrs <- 1/var.a.hrs
prec.a.int <- 1/var.a.int

precisions <- data.frame(HRS_only = c(prec.psi.hrs,prec.a.hrs),
                         integrated = c(prec.psi.int, prec.a.int),
                         row.names = c("precision of psi", "precision of a"))



# Eristalis nemorum ----
load(paste0("./Output/Rerun/", best_species$x[16], "_HRS_model_it161000.rdata"))

load(paste0("./Output/Rerun/", best_species$x[16], "_integrated_model_it161000.rdata"))


out_integrated_model$BUGSoutput$summary

# # check convergence visually

traceplot(out_integrated_model, varname = c("a", "pan.p", "dtype1.p", "dtype2.p", "dtype3.p",
                                            "psi.fs", "beta1", "beta2", "beta3"),
          mfrow = c(3,3), ask = FALSE)


# # plot posterior distribution of occupancy and compare to
# # prop of sites in which species was detected in the 2 datasets

par(mfrow = c(1,1))
hist(out_integrated_model$BUGSoutput$sims.list$psi.fs, xlim = c(0,1),
     col = "gray", main = "", xlab = "Occupancy")

abline(v = sum(aggregate(as.integer(HRS_formatted$spp_vis$'Eristalis nemorum'),
                         by = list(HRS_formatted$occDetdata$site), max)[,2])/length(unique(Site1)),
       lwd = 2, col = "red")
abline(v = sum(aggregate(as.integer(PoMS_formatted$spp_vis$'Eristalis nemorum'),
                         by = list(PoMS_formatted$occDetdata$site), max)[,2])/length(unique(Site2)),
       lwd = 2, col = "blue")
legend("topright", legend = c("Reporting rate HRS", "Reporting rate PoMS"), col = c("red", "blue"), lty = 1, lwd = 2)



out_HRS_model$BUGSoutput$summary

# # check convergence visually

traceplot(out_HRS_model, varname = c("a", "dtype1.p", "dtype2.p", "dtype3.p", "psi.fs",
                                     "beta1", "beta2", "beta3"),
          mfrow = c(3,3), ask = FALSE)


# # plot posterior distribution of occupancy and compare to
# # prop of sites in which species was detected in the 2 datasets
par(mfrow = c(1,1))
hist(out_HRS_model$BUGSoutput$sims.list$psi.fs, xlim = c(0,1),
     col = "gray", main = "", xlab = "Occupancy")
abline(v = sum(aggregate(as.integer(HRS_formatted$spp_vis$'Eristalis nemorum'),
                         by = list(HRS_formatted$occDetdata$site), max)[,2])/length(unique(Site1)),
       lwd = 2, col = "red")
legend("topleft", legend = c("Reporting rate HRS"), col = c("red"), lty = 1, lwd = 2)






# # Compare the precison of the occupancy

var.psi.hrs <- var(out_HRS_model$BUGSoutput$sims.list$psi.fs)
var.psi.int <- var(out_integrated_model$BUGSoutput$sims.list$psi.fs)

prec.psi.hrs <- 1/var.psi.hrs
prec.psi.int <- 1/var.psi.int

var.a.hrs <- var(out_HRS_model$BUGSoutput$sims.list$a)
var.a.int <- var(out_integrated_model$BUGSoutput$sims.list$a)

prec.a.hrs <- 1/var.a.hrs
prec.a.int <- 1/var.a.int

precisions <- data.frame(HRS_only = c(prec.psi.hrs,prec.a.hrs),
                         integrated = c(prec.psi.int, prec.a.int),
                         row.names = c("precision of psi", "precision of a"))


# Cheilosia illustrata ----
load(paste0("./Output/Rerun/", best_species$x[17], "_HRS_model_it161000.rdata"))

load(paste0("./Output/Rerun/", best_species$x[17], "_integrated_model_it161000.rdata"))


out_integrated_model$BUGSoutput$summary

# # check convergence visually

traceplot(out_integrated_model, varname = c("a", "pan.p", "dtype1.p", "dtype2.p", "dtype3.p",
                                            "psi.fs", "beta1", "beta2", "beta3"),
          mfrow = c(3,3), ask = FALSE)


# # plot posterior distribution of occupancy and compare to
# # prop of sites in which species was detected in the 2 datasets

par(mfrow = c(1,1))
hist(out_integrated_model$BUGSoutput$sims.list$psi.fs, xlim = c(0,1),
     col = "gray", main = "", xlab = "Occupancy")

abline(v = sum(aggregate(as.integer(HRS_formatted$spp_vis$'Cheilosia illustrata'),
                         by = list(HRS_formatted$occDetdata$site), max)[,2])/length(unique(Site1)),
       lwd = 2, col = "red")
abline(v = sum(aggregate(as.integer(PoMS_formatted$spp_vis$'Cheilosia illustrata'),
                         by = list(PoMS_formatted$occDetdata$site), max)[,2])/length(unique(Site2)),
       lwd = 2, col = "blue")
legend("topright", legend = c("Reporting rate HRS", "Reporting rate PoMS"), col = c("red", "blue"), lty = 1, lwd = 2)



out_HRS_model$BUGSoutput$summary

# # check convergence visually

traceplot(out_HRS_model, varname = c("a", "dtype1.p", "dtype2.p", "dtype3.p", "psi.fs",
                                     "beta1", "beta2", "beta3"),
          mfrow = c(3,3), ask = FALSE)


# # plot posterior distribution of occupancy and compare to
# # prop of sites in which species was detected in the 2 datasets
par(mfrow = c(1,1))
hist(out_HRS_model$BUGSoutput$sims.list$psi.fs, xlim = c(0,1),
     col = "gray", main = "", xlab = "Occupancy")
abline(v = sum(aggregate(as.integer(HRS_formatted$spp_vis$'Cheilosia illustrata'),
                         by = list(HRS_formatted$occDetdata$site), max)[,2])/length(unique(Site1)),
       lwd = 2, col = "red")
legend("topleft", legend = c("Reporting rate HRS"), col = c("red"), lty = 1, lwd = 2)






# # Compare the precison of the occupancy

var.psi.hrs <- var(out_HRS_model$BUGSoutput$sims.list$psi.fs)
var.psi.int <- var(out_integrated_model$BUGSoutput$sims.list$psi.fs)

prec.psi.hrs <- 1/var.psi.hrs
prec.psi.int <- 1/var.psi.int

var.a.hrs <- var(out_HRS_model$BUGSoutput$sims.list$a)
var.a.int <- var(out_integrated_model$BUGSoutput$sims.list$a)

prec.a.hrs <- 1/var.a.hrs
prec.a.int <- 1/var.a.int

precisions <- data.frame(HRS_only = c(prec.psi.hrs,prec.a.hrs),
                         integrated = c(prec.psi.int, prec.a.int),
                         row.names = c("precision of psi", "precision of a"))



# Sphaerophoria scripta ----
load(paste0("./Output/Rerun/", best_species$x[18], "_HRS_model_it161000.rdata"))

load(paste0("./Output/Rerun/", best_species$x[18], "_integrated_model_it161000.rdata"))


out_integrated_model$BUGSoutput$summary

# # check convergence visually

traceplot(out_integrated_model, varname = c("a", "pan.p", "dtype1.p", "dtype2.p", "dtype3.p",
                                            "psi.fs", "beta1", "beta2", "beta3"),
          mfrow = c(3,3), ask = FALSE)


# # plot posterior distribution of occupancy and compare to
# # prop of sites in which species was detected in the 2 datasets

par(mfrow = c(1,1))
hist(out_integrated_model$BUGSoutput$sims.list$psi.fs, xlim = c(0,1),
     col = "gray", main = "", xlab = "Occupancy")

abline(v = sum(aggregate(as.integer(HRS_formatted$spp_vis$'Sphaerophoria scripta'),
                         by = list(HRS_formatted$occDetdata$site), max)[,2])/length(unique(Site1)),
       lwd = 2, col = "red")
abline(v = sum(aggregate(as.integer(PoMS_formatted$spp_vis$'Sphaerophoria scripta'),
                         by = list(PoMS_formatted$occDetdata$site), max)[,2])/length(unique(Site2)),
       lwd = 2, col = "blue")
legend("topright", legend = c("Reporting rate HRS", "Reporting rate PoMS"), col = c("red", "blue"), lty = 1, lwd = 2)



out_HRS_model$BUGSoutput$summary

# # check convergence visually

traceplot(out_HRS_model, varname = c("a", "dtype1.p", "dtype2.p", "dtype3.p", "psi.fs",
                                     "beta1", "beta2", "beta3"),
          mfrow = c(3,3), ask = FALSE)


# # plot posterior distribution of occupancy and compare to
# # prop of sites in which species was detected in the 2 datasets
par(mfrow = c(1,1))
hist(out_HRS_model$BUGSoutput$sims.list$psi.fs, xlim = c(0,1),
     col = "gray", main = "", xlab = "Occupancy")
abline(v = sum(aggregate(as.integer(HRS_formatted$spp_vis$'Sphaerophoria scripta'),
                         by = list(HRS_formatted$occDetdata$site), max)[,2])/length(unique(Site1)),
       lwd = 2, col = "red")
legend("topleft", legend = c("Reporting rate HRS"), col = c("red"), lty = 1, lwd = 2)






# # Compare the precison of the occupancy

var.psi.hrs <- var(out_HRS_model$BUGSoutput$sims.list$psi.fs)
var.psi.int <- var(out_integrated_model$BUGSoutput$sims.list$psi.fs)

prec.psi.hrs <- 1/var.psi.hrs
prec.psi.int <- 1/var.psi.int

var.a.hrs <- var(out_HRS_model$BUGSoutput$sims.list$a)
var.a.int <- var(out_integrated_model$BUGSoutput$sims.list$a)

prec.a.hrs <- 1/var.a.hrs
prec.a.int <- 1/var.a.int

precisions <- data.frame(HRS_only = c(prec.psi.hrs,prec.a.hrs),
                         integrated = c(prec.psi.int, prec.a.int),
                         row.names = c("precision of psi", "precision of a"))



# Scaeva pyrastri ----
load(paste0("./Output/Rerun/", best_species$x[19], "_HRS_model_it161000.rdata"))

load(paste0("./Output/Rerun/", best_species$x[19], "_integrated_model_it161000.rdata"))


out_integrated_model$BUGSoutput$summary

# # check convergence visually

traceplot(out_integrated_model, varname = c("a", "pan.p", "dtype1.p", "dtype2.p", "dtype3.p",
                                            "psi.fs", "beta1", "beta2", "beta3"),
          mfrow = c(3,3), ask = FALSE)


# # plot posterior distribution of occupancy and compare to
# # prop of sites in which species was detected in the 2 datasets

par(mfrow = c(1,1))
hist(out_integrated_model$BUGSoutput$sims.list$psi.fs, xlim = c(0,1),
     col = "gray", main = "", xlab = "Occupancy")

abline(v = sum(aggregate(as.integer(HRS_formatted$spp_vis$'Scaeva pyrastri'),
                         by = list(HRS_formatted$occDetdata$site), max)[,2])/length(unique(Site1)),
       lwd = 2, col = "red")
abline(v = sum(aggregate(as.integer(PoMS_formatted$spp_vis$'Scaeva pyrastri'),
                         by = list(PoMS_formatted$occDetdata$site), max)[,2])/length(unique(Site2)),
       lwd = 2, col = "blue")
legend("topright", legend = c("Reporting rate HRS", "Reporting rate PoMS"), col = c("red", "blue"), lty = 1, lwd = 2)



out_HRS_model$BUGSoutput$summary

# # check convergence visually

traceplot(out_HRS_model, varname = c("a", "dtype1.p", "dtype2.p", "dtype3.p", "psi.fs",
                                     "beta1", "beta2", "beta3"),
          mfrow = c(3,3), ask = FALSE)


# # plot posterior distribution of occupancy and compare to
# # prop of sites in which species was detected in the 2 datasets
par(mfrow = c(1,1))
hist(out_HRS_model$BUGSoutput$sims.list$psi.fs, xlim = c(0,1),
     col = "gray", main = "", xlab = "Occupancy")
abline(v = sum(aggregate(as.integer(HRS_formatted$spp_vis$'Scaeva pyrastri'),
                         by = list(HRS_formatted$occDetdata$site), max)[,2])/length(unique(Site1)),
       lwd = 2, col = "red")
legend("topleft", legend = c("Reporting rate HRS"), col = c("red"), lty = 1, lwd = 2)






# # Compare the precison of the occupancy

var.psi.hrs <- var(out_HRS_model$BUGSoutput$sims.list$psi.fs)
var.psi.int <- var(out_integrated_model$BUGSoutput$sims.list$psi.fs)

prec.psi.hrs <- 1/var.psi.hrs
prec.psi.int <- 1/var.psi.int

var.a.hrs <- var(out_HRS_model$BUGSoutput$sims.list$a)
var.a.int <- var(out_integrated_model$BUGSoutput$sims.list$a)

prec.a.hrs <- 1/var.a.hrs
prec.a.int <- 1/var.a.int

precisions <- data.frame(HRS_only = c(prec.psi.hrs,prec.a.hrs),
                         integrated = c(prec.psi.int, prec.a.int),
                         row.names = c("precision of psi", "precision of a"))



# Sericomyia silentis ----
load(paste0("./Output/Rerun/", best_species$x[20], "_HRS_model_it161000.rdata"))

load(paste0("./Output/Rerun/", best_species$x[20], "_integrated_model_it161000.rdata"))


out_integrated_model$BUGSoutput$summary

# # check convergence visually

traceplot(out_integrated_model, varname = c("a", "pan.p", "dtype1.p", "dtype2.p", "dtype3.p",
                                            "psi.fs", "beta1", "beta2", "beta3"),
          mfrow = c(3,3), ask = FALSE)


# # plot posterior distribution of occupancy and compare to
# # prop of sites in which species was detected in the 2 datasets

par(mfrow = c(1,1))
hist(out_integrated_model$BUGSoutput$sims.list$psi.fs, xlim = c(0,1),
     col = "gray", main = "", xlab = "Occupancy")

abline(v = sum(aggregate(as.integer(HRS_formatted$spp_vis$'Sericomyia silentis'),
                         by = list(HRS_formatted$occDetdata$site), max)[,2])/length(unique(Site1)),
       lwd = 2, col = "red")
abline(v = sum(aggregate(as.integer(PoMS_formatted$spp_vis$'Sericomyia silentis'),
                         by = list(PoMS_formatted$occDetdata$site), max)[,2])/length(unique(Site2)),
       lwd = 2, col = "blue")
legend("topright", legend = c("Reporting rate HRS", "Reporting rate PoMS"), col = c("red", "blue"), lty = 1, lwd = 2)



out_HRS_model$BUGSoutput$summary

# # check convergence visually

traceplot(out_HRS_model, varname = c("a", "dtype1.p", "dtype2.p", "dtype3.p", "psi.fs",
                                     "beta1", "beta2", "beta3"),
          mfrow = c(3,3), ask = FALSE)


# # plot posterior distribution of occupancy and compare to
# # prop of sites in which species was detected in the 2 datasets
par(mfrow = c(1,1))
hist(out_HRS_model$BUGSoutput$sims.list$psi.fs, xlim = c(0,1),
     col = "gray", main = "", xlab = "Occupancy")
abline(v = sum(aggregate(as.integer(HRS_formatted$spp_vis$'Sericomyia silentis'),
                         by = list(HRS_formatted$occDetdata$site), max)[,2])/length(unique(Site1)),
       lwd = 2, col = "red")
legend("topleft", legend = c("Reporting rate HRS"), col = c("red"), lty = 1, lwd = 2)






# # Compare the precison of the occupancy

var.psi.hrs <- var(out_HRS_model$BUGSoutput$sims.list$psi.fs)
var.psi.int <- var(out_integrated_model$BUGSoutput$sims.list$psi.fs)

prec.psi.hrs <- 1/var.psi.hrs
prec.psi.int <- 1/var.psi.int

var.a.hrs <- var(out_HRS_model$BUGSoutput$sims.list$a)
var.a.int <- var(out_integrated_model$BUGSoutput$sims.list$a)

prec.a.hrs <- 1/var.a.hrs
prec.a.int <- 1/var.a.int

precisions <- data.frame(HRS_only = c(prec.psi.hrs,prec.a.hrs),
                         integrated = c(prec.psi.int, prec.a.int),
                         row.names = c("precision of psi", "precision of a"))



