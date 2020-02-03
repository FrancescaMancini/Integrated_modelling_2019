#########################################
## Script to create list of parameters 
## and job script for JASMIN
## Adapted from: Mark Logie
## Date created: 2019-12-17
## Date modified:
#########################################

library(dplyr)
library(sparta)
rm(list = ls())

# Load the data -----

hrs_2017 <- read.csv("../Data/HRS/HRS_2017_spNames.csv", header = T)

POMS_2017 <- read.csv("../Data/PoMS/POMS_2017_5Km.csv", header = T)

# format the data for jags
# HRS
# create year to supply as value for closure_period to formatOccData
# this is because we cannot use date as to survey = 
# because having decreased te precision of the gridref
# we might have multiple visits in the same 5 Km site on the same day
# and we want to keep them as replicates

hrs_2017$YEAR <- format(as.Date(hrs_2017$enddate, format = "%Y-%m-%d"),"%Y")

HRS_formatted <- formatOccData(taxa = hrs_2017$NAME, 
                               survey = paste(hrs_2017$enddate,hrs_2017$TO_GRIDREF, sep = "-"), 
                               # create survey value from date and 1Km gridref to keep all visits as replicates
                               site = hrs_2017$GRIDREF_5KM_PREC,
                               closure_period = hrs_2017$YEAR)


# save(HRS_formatted, file = "../Data/HRS/HRS_formatted.rdata")

# POMS
# for now just change the dates by a number of day that is equal to the pantrap station number (1 to 5)
# eventually we want to modify formatOccData to allow for another identifier for survey
POMS_2017$dates_new <- as.Date(POMS_2017$dates, format = "%Y-%m-%d") + POMS_2017$pan_trap_station

# create year to supply as value for closure_period to formatOccData
# this is because we cannot use date as to survey = 
# because having decreased te precision of the gridref
# we might have multiple visits in the same 5 Km site on the same day
# and we want to keep them as replicates

POMS_2017$YEAR <- format(POMS_2017$dates_new,"%Y")

PoMS_formatted <- formatOccData(taxa = POMS_2017$taxon_standardised,
                                survey = paste(POMS_2017$dates_new, POMS_2017$sample_gridref, sep = "-"),
                                # create survey value from date and 1Km gridref to keep all visits as replicates
                                site = POMS_2017$GRIDREF_5KM_PREC,
                                closure_period = POMS_2017$YEAR)

# save(PoMS_formatted, file = "../Data/PoMS/POMS_formatted.rdata")


# find the commen species in the two datasets
taxa <- intersect(POMS_2017$taxon_standardised, hrs_2017$NAME)
length(taxa)


chain <- 200
start_it <- 200000
max_it <- 400000

write.csv(data.frame(species = lapply(taxa,
                                      FUN = function(x){rep(x,chain)}) %>% unlist(),
                     start = rep(seq(start_it,((max_it/chain)*(chain-1)+1),(max_it/chain)),length(taxa)),
                     end = rep(seq((max_it/chain),max_it,(max_it/chain)),length(taxa)),
                     chain = rep(200:400,length(taxa)),
                     seed = round(runif(n = chain*length(taxa), min = 0, max = 2147483647),0)),
          'parameters.csv',row.names = FALSE)

job_num <- paste0(seq(1,(chain*length(taxa)-chain+1),chain),'-',
                  seq(chain,(chain*length(taxa)),chain))
bsub_lines <- paste0('bsub -q short-serial -n 1 -R "rusage[mem=20000]" -M 20000',
                     ' -W 23:59 -J HRS_modelling_Job[', job_num,
                     ']%1 -oo Outputs/Log_files/R-%J-%I.o -eo Outputs/Log_files/R-%J-%I.e < multi_array.job')
write.table(bsub_lines, file = 'bsub_lines.sh', sep = '',
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# Run create_visit_data.R
# chmod +x bsub_lines.sh (may need to open in WinSCP, then save, to remove windows linebreaks)
# ./bsub_lines.sh

# Restart failures
# From ./functions
# source('restart_failures.R')
# restart_failures(dir = '../../dragonflies',originalJob = '../../dragonflies/bsub_lines.sh')

# Combine chains
# From ./functions
# source('combine_chains.R')
# combine_chains(data_dir = '../../dragonflies', target = 3340)

# Get posteriors
# From ./functions
# source('extract_posteriors.R')
# extract_posteriors(data_dir = '../../dragonflies/output',output_dir = '../../dragonflies/output/posteriors, target = 3340)

# Plot occupancy
# source('plot_occupancy.R')
# summary_occ <- get_summary_occ(output_dir = '../../dragonflies', source_data = '../../dragonflies/data/dragonflies_outputs.rdata')
# NOT RUN create_summary_plots(summary_occ, plots_dir = '../../dragonflies/plots')

# Calculate and plot trend estimates
# Run full calculate_trend.R script

