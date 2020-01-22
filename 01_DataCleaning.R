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

# poss_errors_grouped <- poss_errors %>%
#   group_by(TO_GRIDREF, startdate) %>%
#   summarise(list_length = n_distinct(NAME)) %>%
#   ungroup()
# 
# plot(poss_errors_grouped$list_length ~ poss_errors_grouped$startdate)


poss_duplicates <- hrs_sp_names %>%
  group_by(TO_GRIDREF, CONCEPT, ABUNDANCE_COMMENT) %>%
  filter(startdate %in% poss_errors$startdate | 
         startdate %in% as.Date(poss_errors$startdate, format = "%Y-d%-m%")) %>%
  filter(n() > 1) #%>%
#  do(filter(as.Date(TO_STARTDATE, format = "%d/%m/%Y")))






hrs_sp_names$JulDate <- as.numeric(format(hrs_sp_names$startdate, "%j"))
hist(hrs_sp_names$JulDate)



hrs_sp_names$YEAR <- format(hrs_sp_names$startdate,"%Y")


HRS_formatted <- formatOccData(taxa = hrs_sp_names$NAME, 
                                      survey = paste(hrs_sp_names$startdate,
                                                     hrs_sp_names$TO_GRIDREF, sep = "-"), 
                                      # create survey value from date and 1Km gridref to keep all visits as replicates
                                      site = hrs_sp_names$GRIDREF_5KM_PREC,
                                      closure_period = hrs_sp_names$YEAR)


HRS_formatted$occDetdata$date <- as.numeric(format(as.POSIXlt(substr(HRS_formatted$occDetdata$visit, 7, 16),
                                                                     format = "%Y-%m-%d"), "%j"))


plot(HRS_formatted$occDetdata$L ~ HRS_formatted$occDetdata$date)

# find the outliers
identify(HRS_formatted$occDetdata$L ~ HRS_formatted$occDetdata$date)
outliers <- HRS_formatted$occDetdata[c(2434, 2435, 2437, 2754, 3196, 3358, 11741),]
# extract them from the dataset
outliers_df <- filter(hrs_sp_names, GRIDREF_5KM_PREC %in% outliers$site & JulDate %in% outliers$date)



HRS_formatted$occDetdata$DAY <- as.numeric(format(as.POSIXlt(substr(HRS_formatted$occDetdata$visit, 7, 16),
                                                                    format = "%Y-%m-%d"), "%d"))

HRS_formatted$occDetdata$MONTH <- as.numeric(format(as.POSIXlt(substr(HRS_formatted$occDetdata$visit, 7, 16),
                                                                      format = "%Y-%m-%d"), "%m"))

HRS_formatted$occDetdata$YEAR <- as.numeric(format(as.POSIXlt(substr(HRS_formatted$occDetdata$visit, 7, 16),
                                                               format = "%Y-%m-%d"), "%Y"))

ggplot(HRS_formatted$occDetdata, aes(x = DAY, y = L)) +
  geom_point() +
  facet_wrap(~ MONTH)



hrs_winter <- subset(hrs_sp_names, JulDate < 80 | JulDate > 300)
hist(hrs_winter$JulDate)


hrs_winter$YEAR <- format(hrs_winter$startdate,"%Y")
hrs_winter$MONTH <- format(hrs_winter$startdate,"%m")
hrs_winter$DAY <- format(hrs_winter$startdate,"%d")



HRS_winter_formatted <- formatOccData(taxa = hrs_winter$NAME,
                                      survey = paste(hrs_winter$startdate,hrs_winter$TO_GRIDREF, sep = "-"),
                                      # create survey value from date and 1Km gridref to keep all visits as replicates
                                      site = hrs_winter$GRIDREF_5KM_PREC,
                                      closure_period = hrs_winter$YEAR)

HRS_winter_formatted$occDetdata$date <- as.numeric(format(as.POSIXlt(substr(HRS_winter_formatted$occDetdata$visit, 7, 16),
                                                                     format = "%Y-%m-%d"), "%j"))


plot(HRS_winter_formatted$occDetdata$L ~ HRS_winter_formatted$occDetdata$date)


HRS_winter_formatted$occDetdata$DAY <- as.numeric(format(as.POSIXlt(substr(HRS_winter_formatted$occDetdata$visit, 7, 16),
                                                                     format = "%Y-%m-%d"), "%d"))

HRS_winter_formatted$occDetdata$MONTH <- as.numeric(format(as.POSIXlt(substr(HRS_winter_formatted$occDetdata$visit, 7, 16),
                                                                    format = "%Y-%m-%d"), "%m"))


# HRS_winter_coordinates <- gr2gps_latlon(HRS_winter_formatted$occDetdata$site, precision = 5000)
# 
# HRS_winter_formatted$occDetdata$LATITUDE <- HRS_winter_coordinates$LATITUDE
# HRS_winter_formatted$occDetdata$LONGITUDE <- HRS_winter_coordinates$LONGITUDE
# 
# UK_df <- fortify(UK$britain)
# hrs_new_coordinates <- LatLong_Cartesian(HRS_winter_formatted$occDetdata$LONGITUDE,
#                                          HRS_winter_coordinates$LATITUDE)
# HRS_winter_formatted$occDetdata$LATITUDE <- hrs_new_coordinates$y
# HRS_winter_formatted$occDetdata$LONGITUDE <- hrs_new_coordinates$x
# 
# ggplot() +
#   geom_polygon(data = UK_df, aes(x = long, y = lat, group = id)) +
#   geom_point(data = HRS_winter_formatted$occDetdata, 
#              aes(x = LONGITUDE, y = LATITUDE, size = L), alpha = 0.5)


color <- rgb(0,0,0, alpha = 0.3)

par(mar = c(0.1,0.1,1,0.1))

plot_GIS(UK$britain, new.window = FALSE, show.axis = FALSE, 
         show.grid = FALSE, xlab = "", ylab = "")

plotUK_gr(HRS_winter_formatted$occDetdata$site, gr_prec = 5000, border = "red")

plotUK_gr_points(HRS_winter_formatted$occDetdata$site, 
                 cex = HRS_winter_formatted$occDetdata$L/5,
                 col = color, pch = 19)




# transform grid references to 10Km resolution
# this is just to be able to calculate statistics for the next step

HRS_formatted$occDetdata$region <- reformat_gr(HRS_formatted$occDetdata$site, prec_out = 10000, 
                                    precision = 5000)

HRS_L_summary <- HRS_formatted$occDetdata %>%
  group_by(region, MONTH) %>%
  summarise(median_L = median(L), 
            min_L = min(L),
            max_L = max(L),
            quantile_95 = quantile(L, probs = 0.95))
