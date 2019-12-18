rm(list = ls())

# MCMC settings
ni <- 1000
nt <- 3
nb <- 500
nc <- 3



library(R2jags)
library(dplyr)
library(sparta)

#location of data input and output

data_dir <- "Data/"

output_folder <- "./Outputs/Integrated"


# Load the data -----

load(paste0(data_dir,"HRS_formatted.rdata"))

load(paste0(data_dir,"POMS_formatted.rdata"))


# Get job index
index <- as.numeric(commandArgs(trailingOnly = TRUE))[1]
str(index)


parameters <- read.csv(paste0(data_dir,"parameters.csv"))[index,]
print(parameters)

set.seed(parameters$seed)


if(parameters$start == 1){
  
  # Data reformatting ----
  # y1 is the detection history for one species from the HRS dataset
  
  # Site1 is a list of sites visited in HRS dataset (no unique)
  # nsite is the total number of sites visited (unique)
  
  # nvisit1 is the numbe of visits in the HRS dataset
  
  # DATATYPE2 is a column that indicates if that visit produced a medium list length (1) or not (0) (for HRS dataset only)
  # DATATYPE3 same as above for long list length
  
  # create Julian day variable from the date in visit (HRS)
  # The variable visit at the moment is a combination of the 5Km grid ref
  # the date and the 1Km grid ref
  # the code extracts the date from these character strings (characters number 7 to 16)
  # formats them as dates and converts them to numbers
  JulDate1 <- as.numeric(format(as.POSIXlt(substr(HRS_formatted$spp_vis$visit, 7, 16),
                                           format = "%Y-%m-%d"), "%j"))
  
  Site1 <- HRS_formatted$occDetdata$site
  DATATYPE2 <- HRS_formatted$occDetdata %>%
    mutate(L = case_when(L < 2 ~ 0,
                         L < 4 ~ 1,
                         L >= 4 ~ 0)) %>%
    select(L)
  
  
  DATATYPE3 <- HRS_formatted$occDetdata %>%
    mutate(L = case_when(L < 4 ~ 0,
                         L >= 4 ~ 1)) %>%
    select(L)
  
  # create Julian day variable from the date in visit (POMS)
  # The variable visit at the moment is a combination of the 5Km grid ref
  # the date and the 1Km grid ref
  # the code extracts the date from these character strings (characters number 7 to 16)
  # formats them as dates and converts them to numbers
  JulDate2 <- as.numeric(format(as.POSIXlt(substr(PoMS_formatted$spp_vis$visit, 7, 16),
                                           format = "%Y-%m-%d"), "%j")) 
  
  
  Site2 <- PoMS_formatted$occDetdata$site
  
  
  nsite <- length(unique(c(as.character(Site1), as.character(Site2))))
  
  
  
  ## HRS-only model - Julian data as gaussian density function ----
  
  
  write.table(paste("Integrated model for species",parameters$species,"has begun at ",date(),sep=" "),
              paste("./Outputs/progress_int.txt",sep=""),
              append=TRUE,col.names=FALSE,row.names=FALSE,quote=FALSE)
  
  # format the data - HRS
  
  y1 <- as.integer(HRS_formatted$spp_vis[, parameters$species]) 
  
  nvisit1 <- length(y1)
  
  # format the data - POMS
  y2 <- as.integer(PoMS_formatted$spp_vis[, parameters$species])
  
  nvisit2 <- length(y2)
  
  
  
  ## Integrated model - Julian data as gaussian density function ----
  
  jags_integrated_data <- list(y1 = y1, y2 = y2, Site1 = Site1, Site2 = Site2, nvisit1 = nvisit1,
                               nvisit2 = nvisit2, DATATYPE2 = DATATYPE2, DATATYPE3 = DATATYPE3, 
                               nsite = nsite, JulDate1 = JulDate1, JulDate2 = JulDate2)
  
  inits_integrated <- function() list(a = 2, dtype1.p = -0.5, dtype2.p = 1.1, 
                                      dtype3.p = 2, pan.p = 0.2, beta1 = 180, 
                                      beta2 = 16, beta3 = 0.1)
  
  params_integrated <- c("a", "dtype1.p", "dtype2.p", "dtype3.p", "pan.p", 
                         "psi.fs", "beta1", "beta2", "beta3")
  
  
  out_integrated_model <- jags(data = jags_integrated_data, parameters.to.save = params_integrated,
                               model.file = "integrated_model_3.txt", inits = inits_integrated,
                               n.chains=nc, n.iter=ni, n.thin=nt, n.burnin=nb)
  


  save(out_integrated_model, file = file.path(output_folder,
                                       paste0(paste(parameters$species,
                                                    paste0("it",
                                                           parameters$end),
                                                    sep = '_'),
                                              '.rdata')))
  
} else {
  
  source('daisy_chain_functions.R')
  
  # get the path to the previous output
  filename <- file.path(output_folder,
                        paste0(paste(parameters$species,
                                     paste0("it", 
                                            parameters$start - 1),
                                     sep = '_'),
                               '.rdata'))
  # run the daisy
  out_integrated_model <- daisyChain(rdataFile = filename,
                              outDir = output_folder,
                              nDaisy = 1,
                              by.it = length(parameters$start:parameters$end),
                              n.thin = 3,
                              quiet = TRUE, 
                              model = "integrated", 
                              update = TRUE)
  
  save(out_integrated_model, file = file.path(output_folder,
                                       paste0(paste(parameters$species,
                                                    paste0("it",
                                                           parameters$end),
                                                    sep = '_'),
                                              '.rdata')))
  
}


quit(save = 'no', runLast = FALSE)
