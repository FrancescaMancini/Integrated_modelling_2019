#######################################
## Data cleaning
## Author: Francesca Mancini
## Date created: 2019-08-08
## Date modified: 2019-08-19
#######################################

library(BRCmap)
library(RODBC)
library(dplyr)
library(ggplot2)
library(sparta)


# load anc clean the hoverfly recording scheme data ----
hrs <- read.csv("..\\Data\\HRS\\RD_HRS.csv", header = T, stringsAsFactors = F)
str(hrs)

summary(hrs$DT_ID)
# I should delete everything that has DT_ID = 8, no date provided

# hrs <- subset(hrs, DT_ID != 8)
hrs <- subset(hrs, DT_ID == 1)


# convert to date format
hrs$startdate <- as.Date(hrs$TO_STARTDATE, format = "%d/%m/%Y")
hrs$enddate <- as.Date(hrs$TO_ENDDATE, format = "%d/%m/%Y")

# retain only 2017 records
hrs <- subset(hrs, format(hrs$enddate,"%Y")==2017)
# retain only records with precision 1Km or less
hrs <- subset(hrs, TO_PRECISION <= 1000)

# change precision of gridref to 5 Km to increase overlap with PoMS data
hrs$GRIDREF_5KM_PREC <- reformat_gr(hrs$TO_GRIDREF, prec_out = 5000, 
                                         precision = hrs$TO_PRECISION)


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
#                                    paste(shQuote(unique(hrs[, "CONCEPT"]), type="sh"), collapse = ","), ") and valid in ('V', 'X')"))
# 
# write.csv(sp_names, "../Data/HRS/Sp_Names_all.csv", row.names = FALSE)

sp_names <- read.csv("../Data/HRS/Sp_Names_all.csv")

str(sp_names)

# returns 240 species names, when we only have 238 species codes
which(duplicated(sp_names$CONCEPT) == TRUE)

# species Dip_3430 and Dip_3434 have two associated names: Eristalis abusivus/abusiva and Eristalis intricarius/intricaria
# in PoMS these species are referred as abusivus and intricarius, therefore we are using same names here

sp_names <- sp_names[-c(117, 122), ]

hrs_sp_names <- left_join(hrs, sp_names, by = "CONCEPT")



## find the mistakes ----

## the only possible erroneous combinations are 
## those where the day is reported as 1 to 12
## Also, the errors that are likely to cause problems are 
## those in the winter months, because that will 
## make the Julian date estimation difficult to converge
## If the 6th of July has been reported as the 7Th of June
## it won't make much of a difference because in both months
## we expect similar numbers of hoverflies recorded


# extract the observation with those characteristics

poss_errors <- subset(hrs_sp_names,
                        format(hrs_sp_names$startdate,"%d") < 13)

summary(poss_errors)

summary(hrs_sp_names)

# group the observations by site and date
# to calculate a list length for every 
# possibly wrong visit and plot the LL by data
poss_errors_grouped <- poss_errors %>%
  group_by(TO_GRIDREF, startdate) %>%
  summarise(list_length = n_distinct(NAME)) %>%
  ungroup()

plot(poss_errors_grouped$list_length ~ poss_errors_grouped$startdate)

# from the plot it is clear that there are some potentially wrong observations
# in January, February and December, we can extract them from the plot

poss_outliers <- identify(poss_errors_grouped$startdate, poss_errors_grouped$list_length)

poss_outliers <- poss_errors_grouped[outliers,]

# Following Kath's observation that some of these records
# may be duplicates of correct records, where the date has been transposed
# we can try to extract from the original dataset all the records
# from the sites in outliers and the same and trasnposed date

poss_duplicates <- hrs_sp_names %>%
  filter(TO_GRIDREF %in% outliers$TO_GRIDREF) %>%
  filter(startdate %in% outliers$startdate | 
         startdate %in% as.Date(outliers$startdate, format = "%Y-d%-m%")) #%>%
  # group_by(TO_GRIDREF, CONCEPT, ABUNDANCE_COMMENT) %>%
  # filter(n() > 1) #%>%
#  do(filter(as.Date(TO_STARTDATE, format = "%d/%m/%Y")))



# write.csv(poss_duplicates, "./Data/HRS_2017_possible_duplicates.csv",
#           row.names = FALSE)

# after manuallly looking at the csv file
# it appears that the records are not duplicated, 
# therefore the best thing to do at the moment
# is to exclude all those records that are very unlikely
# from the dataset
# these dates will have to be revisited and 
# the accuracy f the dates verified from the original data files

plot(poss_errors_grouped$list_length ~ poss_errors_grouped$startdate)

outliers <- identify(poss_errors_grouped$startdate, poss_errors_grouped$list_length)

outliers <- poss_errors_grouped[outliers,]

# filter the original data using site and date from outliers

hrs_sp_names_filtered <- hrs_sp_names %>%
  anti_join(outliers, by = c("TO_GRIDREF", "startdate"))

# make sure all outliers have been removed

hrs_sp_names_filtered_grouped <- hrs_sp_names_filtered %>%
  group_by(TO_GRIDREF, startdate) %>%
  summarise(list_length = n_distinct(NAME)) %>%
  ungroup()

plot(hrs_sp_names_filtered_grouped$list_length ~ hrs_sp_names_filtered_grouped$startdate)


write.csv(hrs_sp_names_filtered, "./Data/HRS_2017_filtered.csv",
          row.names = FALSE)
