
library(dplyr)   
library(tidyr)
library(sp) 
library(ncdf4) 
library(terra)
library(raster) 
library(sf)
library(dplyr)   
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(tmap)
library(RColorBrewer)
library(scales)
library("gridExtra") 
library(sp) 
library(ncdf4) 
library(terra)
library(raster) 
library(sf)
library(dplyr)   
library(stars)
library(viridis)
library(ggplot2)
library(scales) 
library(forcats)
library(stringr)
library(pracma)
library(ggrepel)
library(extrafont)
library(spdep)
library(showtext)
library(sysfonts)
library(grid)
library(trend)
library(ggtext)
library(tibble)
library(ggh4x)
library(patchwork)
library(maps)
options(warn = -1)
options(scipen = 999)
font_add_google("Lato", "lato")
showtext_auto()



#################################################################################
#################################################################################
#  Data preparation  ############################################################
#################################################################################
#################################################################################

#Surface NO2
no2_tract <- st_sf(state_id = integer(), state_name = character(), county_id = integer(), tract_id = integer(), GEOID = numeric(), year = numeric(), surface_no2 = numeric(), geometry = st_sfc(), crs = 4269)
for (i in c(2019:2024)) {
  no2_raster <- terra::rast(paste0("/Users/k_sy_n_imac/Library/CloudStorage/Dropbox/GWU/Lab/FINESST/Research/Data/NO2/annual_mean_tropomi_lur_conus_surface_no2_", i, ".v1.02.nc"))
  
  no2_tract_ <- terra::extract(no2_raster, ustract, fun=mean, exact=TRUE, weights=TRUE, na.rm=TRUE)
  no2_tract_$year <- i
  no2_tract_ <- cbind(ustract_df, no2_tract_)
  no2_tract_ <- no2_tract_[, c("state_id", "state_name","county_id", "tract_id", "GEOID", "year", "surface_no2", "geometry")]
  no2_tract_ <- st_as_sf(no2_tract_)
  
  no2_tract <- rbind(no2_tract, no2_tract_)
}
no2_tract_df <- st_drop_geometry(no2_tract)
no2_tract_df <- no2_tract_df[, c("GEOID", "state_name", "year", "surface_no2")]
no2_tract_df$GEOID <- as.numeric(no2_tract_df$GEOID)

xx <- subset(no2_tract_df, is.na(surface_no2))
unique(xx$GEOID) # 44 tracts have NA values

write.csv(no2_tract_df, "/Users/k_sy_n_imac/Library/CloudStorage/Dropbox/GWU/Lab/FINESST/Research/Results/no2_tract_df.csv") #save


#Mortality burdens attributable to NO2 (reading csv, Calculation: DB_FINESST_tract_level_mor.R)
final3 <- read.csv("/Users/k_sy_n_imac/Library/CloudStorage/Dropbox/GWU/Lab/FINESST/Research/Data/final3.csv")
final3$rate_pe <- ifelse(final3$pop_over20==0, 0, final3$rate_pe)
final3$rate_low <- ifelse(final3$pop_over20==0, 0, final3$rate_low)
final3$rate_up <- ifelse(final3$pop_over20==0, 0, final3$rate_up)


#State & census tract shapefile
usstate <- st_read("/Users/k_sy_n_imac/Library/CloudStorage/Dropbox/GWU/Lab/FINESST/Research/Data/Boundaries/cb_2019_us_state_500k/cb_2019_us_state_500k.shp")
usstate <- st_crop(usstate, c(xmin = -125.0017, ymin = 20, xmax = -66.9346, ymax = 49.3844)) #CONUS
usstate <- vect(usstate)

ustract <- st_read("/Users/k_sy_n_imac/Library/CloudStorage/Dropbox/GWU/Lab/FINESST/Research/Data/Boundaries/nhgis0001_shape/nhgis0001_shapefile_tl2020_us_tract_2020/US_tract_2020.shp")
ustract <- st_transform(ustract, crs = 4269)
ustract <- st_make_valid(ustract)
ustract <- st_crop(ustract, c(xmin = -125.0017, ymin = 20, xmax = -66.9346, ymax = 49.3844)) #CONUS
ustract <- vect(ustract)

ustract_df <- st_as_sf(ustract)
ustract_df <- as.data.frame(ustract_df)
usstate_df <- as.data.frame(usstate)
usstate_df$state_name <- usstate_df$NAME
usstate_df <- usstate_df[, c("state_name", "STATEFP")]
ustract_df <- merge(ustract_df, usstate_df, by='STATEFP')
ustract_df$state_id <- ustract_df$STATEFP
ustract_df$county_id <- ustract_df$COUNTYFP
ustract_df$tract_id <- ustract_df$TRACTCE
ustract_df <- ustract_df[, c("state_id", "state_name","county_id", "tract_id", "GEOID", "geometry")]
ustract_df$GEOID <- as.numeric(ustract_df$GEOID)
ustract_df<- ustract_df %>% mutate(county_id_2020 = paste0(state_id, county_id))
ustract_df$county_id_2020 <- as.numeric(ustract_df$county_id_2020)

saveRDS(ustract_df, "/Users/k_sy_n_imac/Library/CloudStorage/Dropbox/GWU/Lab/FINESST/Research/Data/ustract_df.rds")

#SES factors
#** Correction for Connecticut GEOID  (see the section below for more info)
# race & tract file (TIGER): 2020 values
# the other SES files: 2022 values
# Using "con" below, I correct Connecticut GEOIDs in the SES data other than race/ethnicity
con <- read.csv("/Users/k_sy_n_imac/Library/CloudStorage/Dropbox/GWU/Lab/FINESST/Research/Data/2022tractcrosswalk.csv")
con <- con[,c(1,2)]
con <- rename(con, GEOID = tract_fips_2020, GEOID_tmp = Tract_fips_2022)

#Race/ethnicity
race <- read.csv("/Users/k_sy_n_imac/Library/CloudStorage/Dropbox/GWU/Lab/FINESST/Research/Data/ACS_2023/race_ethnicity.csv") #ACS2023
race <- race %>%
  mutate(GEOID = str_extract(GEOID, "(?<=US)\\d+"))
race$GEOID_tmp <- as.numeric(race$GEOID)
race$pct_white <- as.numeric(race$percent_non_hispanic_white)
race <- race %>%
  mutate(white_decile = ntile(pct_white, 10))
race <- race[, c("GEOID_tmp", "pct_white", "white_decile")]
race <- race %>%
  dplyr::left_join(con, by = "GEOID_tmp") %>%
  dplyr::mutate(GEOID_tmp = coalesce(GEOID, GEOID_tmp)) %>%  # replace if match
  dplyr::select(-GEOID)  # drop extra column if not needed
ses <- rename(race, GEOID=GEOID_tmp)

#Educational level
edu <- read.csv("/Users/k_sy_n_imac/Library/CloudStorage/Dropbox/GWU/Lab/FINESST/Research/Data/ACS_2023/edu.csv") #ACS2023
edu <- edu %>%
  mutate(GEOID = str_extract(GEOID, "(?<=US)\\d+"))
edu$GEOID_tmp <- as.numeric(edu$GEOID)
edu$percent_over25_less_than_high_school <- as.numeric(edu$percent_over25_less_than_high_school)
edu$pct_high_over <- 100-edu$percent_over25_less_than_high_school
edu <- edu %>%
  mutate(edu_decile = ntile(pct_high_over, 10))
xxxx <- edu[, c("GEOID_tmp","pct_high_over", "edu_decile")]
xxxx <- xxxx %>%
  dplyr::left_join(con, by = "GEOID_tmp") %>%
  dplyr::mutate(GEOID_tmp = coalesce(GEOID, GEOID_tmp)) %>%  # replace if match
  dplyr::select(-GEOID)  # drop extra column if not needed
xxxx <- rename(xxxx, GEOID=GEOID_tmp)
ses <- merge(ses, xxxx, by="GEOID", all=TRUE)

#Employment rate
employ <- read.csv("/Users/k_sy_n_imac/Library/CloudStorage/Dropbox/GWU/Lab/FINESST/Research/Data/ACS_2023/employment.csv") #ACS2023
employ <- employ %>%
  mutate(GEOID = str_extract(GEOID, "(?<=US)\\d+"))
employ$GEOID_tmp <- as.numeric(employ$GEOID)
employ$employ_rate <- employ$Pop_employed/employ$Pop_over16
employ <- employ %>%
  mutate(employ_decile = ntile(employ_rate, 10)) # % of population 16 years and over who worked in the last 5 years
xxxx <- employ[, c("GEOID_tmp","employ_rate", "employ_decile")]
xxxx <- xxxx %>%
  dplyr::left_join(con, by = "GEOID_tmp") %>%
  dplyr::mutate(GEOID_tmp = coalesce(GEOID, GEOID_tmp)) %>%  # replace if match
  dplyr::select(-GEOID)  # drop extra column if not needed
xxxx <- rename(xxxx, GEOID=GEOID_tmp)
ses <- merge(ses, xxxx, by="GEOID", all=TRUE)

#Poverty
pov <- read.csv("/Users/k_sy_n_imac/Library/CloudStorage/Dropbox/GWU/Lab/FINESST/Research/Data/ACS_2023/poverty.csv") #ACS2023
pov$pop_below_poverty_line <- as.numeric(pov$pop_below_poverty_line)
pov$total_pop <- as.numeric(pov$total_pop)
pov$pov_rate <- pov$pop_below_poverty_line/pov$total_pop
pov <- pov %>%
  mutate(GEOID = str_extract(GEOID, "(?<=US)\\d+"))
pov$GEOID_tmp <- as.numeric(pov$GEOID)
pov <- pov %>% mutate(pov_decile = ntile(pov_rate, 10))
xxxx <- pov[, c("GEOID_tmp","pov_rate", "pov_decile")]
xxxx <- xxxx %>%
  dplyr::left_join(con, by = "GEOID_tmp") %>%
  dplyr::mutate(GEOID_tmp = coalesce(GEOID, GEOID_tmp)) %>%  # replace if match
  dplyr::select(-GEOID)  # drop extra column if not needed
xxxx <- rename(xxxx, GEOID=GEOID_tmp)
ses <- merge(ses, xxxx, by="GEOID", all=TRUE)

#Mean household income
income <- read.csv("/Users/k_sy_n_imac/Library/CloudStorage/Dropbox/GWU/Lab/FINESST/Research/Data/ACS_2023/income.csv") #ACS2023
income <- income %>%
  mutate(GEOID = str_extract(GEOID, "(?<=US)\\d+"))
income$GEOID_tmp <- as.numeric(income$GEOID)
income$mean_household_income <- as.numeric(income$mean_household_income)
income <- income %>%
  mutate(house_income_decile = ntile(mean_household_income, 10))
xxxx <- income[, c("GEOID_tmp","mean_household_income", "house_income_decile")]
xxxx <- xxxx %>%
  dplyr::left_join(con, by = "GEOID_tmp") %>%
  dplyr::mutate(GEOID_tmp = coalesce(GEOID, GEOID_tmp)) %>%  # replace if match
  dplyr::select(-GEOID)  # drop extra column if not needed
xxxx <- rename(xxxx, GEOID=GEOID_tmp)
ses <- merge(ses, xxxx, by="GEOID", all=TRUE)

write.csv(ses, "/Users/k_sy_n_imac/Library/CloudStorage/Dropbox/GWU/Lab/FINESST/Research/Results/ses.csv") 


#Urbanicity (MSA/non-MSA)
#MSA categorization by country
urb <- read.csv("/Users/k_sy_n_imac/Library/CloudStorage/Dropbox/GWU/Lab/FINESST/Research/Data/MSA.csv")
urb$county_id_2022 <- urb$County.Code
urb$metro <- "MSA"
urb <- urb[,c("county_id_2022", "metro", "MSA_title")]

#Correction for Connecticut
#** GEOID of tracts in Connecticut changes recently.
#** Changes should be addressed using the csv file for matching previous/current GEOID from 
#** https://www.ctdata.org/blog/geographic-resources-for-connecticuts-new-county-equivalent-geography
#** 879/883 have changed recently
ct <- read.csv("/Users/k_sy_n_imac/Library/CloudStorage/Dropbox/GWU/Lab/FINESST/Research/Data/2022tractcrosswalk.csv")
ct$county_id_2020 <- ct$county_fips_2020
ct$county_id_2022 <- ct$ce_fips_2022
ct <- ct[,c("county_id_2020", "county_id_2022")]
urb <- merge(urb, ct, by="county_id_2022", all.x=TRUE)
urb <- urb %>% mutate(county_id_2020 = coalesce(county_id_2020, county_id_2022))
urb <- urb[, c("county_id_2020", "metro", "MSA_title")]

#MSA & MicroSA 
urb$metro <- ifelse(grepl("MSA$", urb$MSA_title), "MSA",
                        ifelse(grepl("oSA$", urb$MSA_title), "MicroSA", NA))
xx <- ustract_df
xx <- merge(xx, urb, by="county_id_2020", all.x=TRUE)
xx <- xx %>% replace_na(list(metro = "Non-urban"))
urb <- xx[,c("GEOID", "geometry", "metro", "MSA_title")]
urb <- urb[!duplicated(urb$GEOID), ]

#EPA_region
no2_ind <- subset(no2_tract_df, year==2019)
urb1 <- merge(urb, no2_ind, by="GEOID", all.x=TRUE)
urb1 <- subset(urb1, select = -c(X, surface_no2, year))
urb1 <- urb1 %>%
  mutate(
    EPA_region = case_when(
      state_name %in% c("Connecticut", "Maine", "Massachusetts", "New Hampshire", "Rhode Island", "Vermont") ~ "Region 1",
      state_name %in% c("New Jersey", "New York", "Puerto Rico", "U.S. Virgin Islands") ~ "Region 2",
      state_name %in% c("Delaware", "District of Columbia", "Maryland", "Pennsylvania", "Virginia", "West Virginia") ~ "Region 3",
      state_name %in% c("Alabama", "Florida", "Georgia", "Kentucky", "Mississippi", "North Carolina", "South Carolina", "Tennessee") ~ "Region 4",
      state_name %in% c("Illinois", "Indiana", "Michigan", "Minnesota", "Ohio", "Wisconsin") ~ "Region 5",
      state_name %in% c("Arkansas", "Louisiana", "New Mexico", "Oklahoma", "Texas") ~ "Region 6",
      state_name %in% c("Iowa", "Kansas", "Missouri", "Nebraska") ~ "Region 7",
      state_name %in% c("Colorado", "Montana", "North Dakota", "South Dakota", "Utah", "Wyoming") ~ "Region 8",
      state_name %in% c("Arizona", "California", "Hawaii", "Nevada", "American Samoa", 
                        "Commonwealth of the Northern Mariana Islands", "Federated States of Micronesia",
                        "Guam", "Marshall Islands", "Republic of Palau") ~ "Region 9",
      state_name %in% c("Alaska", "Idaho", "Oregon", "Washington") ~ "Region 10",
      TRUE ~ NA_character_
    )
  )

#Name for non-urban areas
urb1 <- urb1 %>%
  mutate(
    MSA_title = case_when(
      metro == "Non-urban" & state_name %in% c("Connecticut", "Maine", "Massachusetts", "New Hampshire", "Rhode Island", "Vermont") ~ "Non-urban (Region 1)",
      metro == "Non-urban" & state_name %in% c("New Jersey", "New York", "Puerto Rico", "U.S. Virgin Islands") ~ "Non-urban (Region 2)",
      metro == "Non-urban" & state_name %in% c("Delaware", "District of Columbia", "Maryland", "Pennsylvania", "Virginia", "West Virginia") ~ "Non-urban (Region 3)",
      metro == "Non-urban" & state_name %in% c("Alabama", "Florida", "Georgia", "Kentucky", "Mississippi", "North Carolina", "South Carolina", "Tennessee") ~ "Non-urban (Region 4)",
      metro == "Non-urban" & state_name %in% c("Illinois", "Indiana", "Michigan", "Minnesota", "Ohio", "Wisconsin") ~ "Non-urban (Region 5)",
      metro == "Non-urban" & state_name %in% c("Arkansas", "Louisiana", "New Mexico", "Oklahoma", "Texas") ~ "Non-urban (Region 6)",
      metro == "Non-urban" & state_name %in% c("Iowa", "Kansas", "Missouri", "Nebraska") ~ "Non-urban (Region 7)",
      metro == "Non-urban" & state_name %in% c("Colorado", "Montana", "North Dakota", "South Dakota", "Utah", "Wyoming") ~ "Non-urban (Region 8)",
      metro == "Non-urban" & state_name %in% c("Arizona", "California", "Hawaii", "Nevada", "American Samoa", "Commonwealth of the Northern Mariana Islands",
                                               "Federated States of Micronesia", "Guam", "Marshall Islands", "Republic of Palau") ~ "Non-urban (Region 9)",
      metro == "Non-urban" & state_name %in% c("Alaska", "Idaho", "Oregon", "Washington") ~ "Non-urban (Region 10)",
      TRUE ~ MSA_title  # leave other rows unchanged
    )
  )

#total pop (tract level)
pop <- read.csv("/Users/k_sy_n_imac/Library/CloudStorage/Dropbox/GWU/Lab/FINESST/Research/Data/DECENNIALDP2020/pop.csv")
pop$total_dc2020 <- pop$total
pop <- pop %>%
  mutate(GEOID = str_extract(GEOID, "(?<=US)\\d+"))
pop$GEOID <- as.numeric(pop$GEOID)
tmp <- pop[,c("GEOID", "total_dc2020")]
urb1 <- merge(urb1, tmp, by="GEOID", all.x=TRUE)

#total pop (MSA level)
pop1 <- urb1 %>% group_by(MSA_title) %>% summarize(total_pop_MSA=sum(total_dc2020))

#Adding another metro cateogory with 1M criteria
urb5 <- merge(urb1, pop1, by="MSA_title")
urb5$metro4 <- urb5$metro
urb5$metro4[urb5$total_pop_MSA >= 1000000 & urb5$metro == "MSA"] <- "Large MSA"
urb5$metro4 <- ifelse(urb5$metro4=="MSA", "Small MSA", urb5$metro4)

#Merge SES data with urbancity data
xxxx <- merge(ses, urb5, by="GEOID", all=TRUE)
saveRDS(xxxx, "/Users/k_sy_n_imac/Library/CloudStorage/Dropbox/GWU/Lab/FINESST/Research/Data/xxxx.rds")











#################################################################################
#################################################################################
#  Quick data read  #############################################################
#################################################################################
#################################################################################
# NO2 data
no2_tract_df <- read.csv("/Users/k_sy_n_imac/Library/CloudStorage/Dropbox/GWU/Lab/FINESST/Research/Results/no2_tract_df.csv") 
# Mortality burdens 
final3 <- read.csv("/Users/k_sy_n_imac/Library/CloudStorage/Dropbox/GWU/Lab/FINESST/Research/Data/final3.csv")
final3$rate_pe <- ifelse(final3$pop_over20==0, 0, final3$rate_pe)
final3$rate_low <- ifelse(final3$pop_over20==0, 0, final3$rate_low)
final3$rate_up <- ifelse(final3$pop_over20==0, 0, final3$rate_up)
# tracts
ustract_df <- readRDS("/Users/k_sy_n_imac/Library/CloudStorage/Dropbox/GWU/Lab/FINESST/Research/Data/ustract_df.rds")
# GEOID, SES, & geometry of census tracts
xxxx <- readRDS("/Users/k_sy_n_imac/Library/CloudStorage/Dropbox/GWU/Lab/FINESST/Research/Data/xxxx.rds")
# ** Will merge these in each block for every use ** ##










#################################################################################
#################################################################################
#  Data checking  ###############################################################
#################################################################################
#################################################################################

#Checking NAs
df_with_na <- xxxx %>% 
  filter(if_any(everything(), is.na))
summary(df_with_na)
df_with_na1 <- merge(df_with_na, ustract_df, by="GEOID")
df_with_na2 <- merge(df_with_na, ustract_df, by="GEOID", all.x=TRUE)

table(df_with_na1$state_name) 
table(df_with_na2$state_name)
sum(is.na(df_with_na1$state_name))
sum(is.na(df_with_na2$state_name))

test <- urb %>% 
  filter(if_any(everything(), is.na))

#Comparison between estimates from state-level and tract-level mor
no2_tract_df <- read.csv("/Users/k_sy_n_imac/Library/CloudStorage/Dropbox/GWU/Lab/FINESST/Research/Results/no2_tract_df.csv")

final3 <- read.csv("/Users/k_sy_n_imac/Library/CloudStorage/Dropbox/GWU/Lab/FINESST/Research/Data/final3.csv")
final3$rate_pe <- ifelse(final3$pop_over20==0, 0, final3$rate_pe)
final3$rate_low <- ifelse(final3$pop_over20==0, 0, final3$rate_low)
final3$rate_up <- ifelse(final3$pop_over20==0, 0, final3$rate_up)
list_year <- c(2019, 2020, 2021, 2022, 2023, 2024) #data is too large to handle all the years and causes at the same time. Select each to proceed
list_cause <- c("All causes", "Respiratory diseases") 
i <- 6
j <- 1
no2_ind <- subset(no2_tract_df, year==list_year[i])
db_ind <- subset(final3, year==list_year[i] & cause_name==list_cause[j])
final4 <- merge(no2_ind, db_ind, by="GEOID")
final4$rate_state_mor <- final4$rate_pe 
final4$case_state_mor <- final4$case_pe 


final3 <- read.csv("/Users/k_sy_n_imac/Library/CloudStorage/Dropbox/GWU/Lab/FINESST/Research/Data/final3_new.csv")
final3$rate_pe <- ifelse(final3$pop_over20==0, 0, final3$rate_pe)
final3$rate_low <- ifelse(final3$pop_over20==0, 0, final3$rate_low)
final3$rate_up <- ifelse(final3$pop_over20==0, 0, final3$rate_up)
list_year <- c(2019, 2020, 2021, 2022, 2023, 2024) #data is too large to handle all the years and causes at the same time. Select each to proceed
list_cause <- c("all", "res") 
i <- 6
j <- 2
no2_ind <- subset(no2_tract_df, year==list_year[i])
db_ind <- subset(final3, year==list_year[i] & cause_name==list_cause[j])
final4_new <- merge(no2_ind, db_ind, by="GEOID")
final4_new$rate_tract_mor <- final4_new$rate_pe 
final4_new$case_tract_mor <- final4_new$case_pe 

merge <- merge(final4, final4_new, by="GEOID")
ggplot(merge, aes(x = rate_state_mor, y = rate_tract_mor)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") 

ggplot(merge, aes(x = case_state_mor, y = case_tract_mor)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") 












#################################################################################
#################################################################################
#  Descriptive analysis  ########################################################
#################################################################################
#################################################################################
# Data
list_year <- c(2019, 2020, 2021, 2022, 2023, 2024) #data is too large to handle all the years and causes at the same time. Select each to proceed
list_cause <- c("All causes", "Cardiovascular diseases", "Respiratory diseases", "Ischemic heart disease", "Tracheal, bronchus, and lung cancer") 
i <- 6
j <- 1
no2_ind <- subset(no2_tract_df, year==list_year[i])
db_ind <- subset(final3, year==list_year[i] & cause_name==list_cause[j])
final4 <- merge(no2_ind, db_ind, by="GEOID")
final4 <- merge(final4, xxxx, by="GEOID")
final4$year <- final4$year.x
# National summary
# Total N
sum(!is.na(final4$surface_no2))
sum(!is.na(final4$case_pe))
# Summary statistics in 2024
summary(final4$surface_no2)
# Summary statistics (full period)
xx <- final4 %>% group_by(year.x) %>% 
  summarize(mean=mean(surface_no2, na.rm=TRUE), 
            pw_mean=sum(surface_no2 * total_dc2020, na.rm = TRUE) / sum(total_dc2020, na.rm = TRUE),
            median=median(surface_no2, na.rm=TRUE), min=min(surface_no2, na.rm=TRUE), max=max(surface_no2, na.rm=TRUE), n=n())
xx
# n/% of tracts above WHO AQG in 2024
sum(final4$surface_no2 > 5.31 & final4$year == 2024, na.rm=TRUE)
sum(final4$surface_no2 > 5.31 & final4$year == 2024, na.rm=TRUE)/sum(!is.na(no2_tract_df$surface_no2) & no2_tract_df$year == 2024)*100
# Urban group of tracts above AQG
x <- subset(final4, final4$surface_no2 > 5.31 & final4$year == 2024, na.rm=TRUE)
table(x$metro4)
# Summary of premature deaths (case)
xx <- final4 %>% group_by(year) %>% summarize(case=sum(case_pe,na.rm=TRUE), case_low=sum(case_low,na.rm=TRUE), case_up=sum(case_up,na.rm=TRUE), pop=sum(pop_over20,na.rm=TRUE))
xx
# Summary of premature deaths (rate: case per 100,000)
xx$case/xx$pop*100000
xx$case_low/xx$pop*100000
xx$case_up/xx$pop*100000

# Sub-group summary
# Number of tracts in each sub-group
table_counts <- final4 %>%
  count(EPA_region, metro4) %>%              # count rows by both factors
  pivot_wider(names_from = metro4,           # make metro2 values become columns
              values_from = n, 
              values_fill = 0)               # fill missing combos with 0
table_counts
table_counts <- final4 %>%
  count(metro4) %>%              # count rows by both factors
  pivot_wider(names_from = metro4,           # make metro2 values become columns
              values_from = n, 
              values_fill = 0)               # fill missing combos with 0
table_counts
# Number of statistical areas in each sub-group
final4 <- final4 %>%
  mutate(
    MSA_title = if_else(
      metro4 == "Non-urban",
      paste0("Non_urban_", state_name.x),
      MSA_title
    )
  )

group_df <- final4 %>% group_by(MSA_title, metro4, EPA_region) %>%
  summarize(n_2024 = n(), .groups = "drop")
msa_regions <- group_df %>% count(MSA_title) %>% filter(n > 1)
majority_regions <- final4 %>%
  group_by(MSA_title, EPA_region) %>%
  summarise(total_pop = sum(pop_over20, na.rm = TRUE), .groups = "drop") %>%
  group_by(MSA_title) %>%
  slice_max(total_pop, n = 1) %>%  # ← Gets region with HIGHEST population
  ungroup() %>%
  dplyr::select(MSA_title, EPA_region_new = EPA_region)
final4 <- final4 %>%
  left_join(majority_regions, by = "MSA_title") %>%
  mutate(EPA_region = if_else(!is.na(EPA_region_new), 
                              EPA_region_new, 
                              EPA_region))

msa_counts_wide <- final4 %>%
  group_by(EPA_region, metro4) %>%
  summarise(n_unique_MSA = n_distinct(MSA_title), .groups = "drop") %>%
  pivot_wider(names_from = metro4,
              values_from = n_unique_MSA,
              values_fill = 0)
msa_counts_wide
# By urbanicity 
xx <- final4 %>% group_by(metro4) %>% summarize(w_mean = weighted.mean(surface_no2, total_dc2020, na.rm = TRUE),case=sum(case_pe,na.rm=TRUE), case_low=sum(case_low,na.rm=TRUE), case_up=sum(case_up,na.rm=TRUE), pop=sum(pop_over20,na.rm=TRUE))
xx
x <- subset(final4, final4$case_pe == 0 & final4$year == 2024, na.rm=TRUE)
table(x$metro4)
# By EPA region
xx <- final4 %>% group_by(EPA_region) %>% summarize(w_mean = weighted.mean(surface_no2, total_dc2020, na.rm = TRUE), case=sum(case_pe,na.rm=TRUE), case_low=sum(case_low,na.rm=TRUE), case_up=sum(case_up,na.rm=TRUE), pop=sum(pop_over20,na.rm=TRUE))
xx


# Figure S1: EPA REGIONS 
# Merge & Choose specific years/cause of death
list_year <- c(2019, 2020, 2021, 2022, 2023, 2024) #data is too large to handle all the years and causes at the same time. Select each to proceed
list_cause <- c("All causes", "Cardiovascular diseases", "Respiratory diseases", "Ischemic heart disease", "Tracheal, bronchus, and lung cancer") 
i <- 6
j <- 1
no2_ind <- subset(no2_tract_df, year==list_year[i])
db_ind <- subset(final3, year==list_year[i] & cause_name==list_cause[j])
final4 <- merge(no2_ind, db_ind, by="GEOID")
final4 <- merge(final4, xxxx, by="GEOID", all.x=TRUE)

region_labels <- c(
  "Region 1 (New England)",
  "Region 2 (Northeast/Caribbean)",
  "Region 3 (Mid-Atlantic)",
  "Region 4 (Southeast)",
  "Region 5 (Great Lakes)",
  "Region 6 (South Central)",
  "Region 7 (Heartland)",
  "Region 8 (Mountains/Plains)",
  "Region 9 (Pacific Southwest)",
  "Region 10 (Northwest)"
)

region_colors <- c(
  "Region 1 (New England)"        = "#1B2F70",
  "Region 2 (Northeast/Caribbean)"= "#2C4FD7",
  "Region 3 (Mid-Atlantic)"       = "#5A8AC6",
  "Region 4 (Southeast)"          = "#8BC8C8",
  "Region 5 (Great Lakes)"        = "#A1D99B",
  "Region 6 (South Central)"      = "#D9C89B",
  "Region 7 (Heartland)"          = "#DDAA66",
  "Region 8 (Mountains/Plains)"   = "#F2A6A0",
  "Region 9 (Pacific Southwest)"  = "#7A5AC9",
  "Region 10 (Northwest)"         = "#4A0066"
)

# ── Recode EPA_region to named labels before plotting ────────────────────────
final4 <- final4 %>%
  mutate(
    EPA_region = recode(EPA_region,
                        "Region 1"  = "Region 1 (New England)",
                        "Region 2"  = "Region 2 (Northeast/Caribbean)",
                        "Region 3"  = "Region 3 (Mid-Atlantic)",
                        "Region 4"  = "Region 4 (Southeast)",
                        "Region 5"  = "Region 5 (Great Lakes)",
                        "Region 6"  = "Region 6 (South Central)",
                        "Region 7"  = "Region 7 (Heartland)",
                        "Region 8"  = "Region 8 (Mountains/Plains)",
                        "Region 9"  = "Region 9 (Pacific Southwest)",
                        "Region 10" = "Region 10 (Northwest)"
    ),
    EPA_region = factor(EPA_region, levels = region_labels)  # enforce order
  )

final4 <- st_as_sf(final4)
st_crs(final4) <- 4326
states_sf <- st_as_sf(maps::map("state", plot = FALSE, fill = TRUE))
states_sf <- st_transform(states_sf, st_crs(final4))

pdf(file = paste0("/Users/k_sy_n_imac/Library/CloudStorage/Dropbox/GWU/Lab/FINESST/Research/Results/EPA_region_map_final.pdf"), 
    width = 9, 
    height = 4)
par(family = "Lato")
ggplot(final4) +
  geom_sf(aes(geometry = geometry, fill = EPA_region), color = NA) +
  geom_sf(data = states_sf, fill = NA, color = "#FAF9F6", linewidth = 0.15) +
  scale_fill_manual(values = region_colors) +
  labs(
    fill = "EPA Region"
    #title = "Census Tracts by EPA Region",
   # subtitle = "Custom Colorblind-Safe Region Palette"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    plot.title = element_text(face = "bold"),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank()
  )
dev.off()







########################################################################################
########################################################################################
# Spatial & Temporal trends  ###########################################################
########################################################################################
########################################################################################
# Figure 1 ###############################################################################
# Spatial: NO2
# Merge & Choose specific years/cause of death
list_year <- c(2019, 2020, 2021, 2022, 2023, 2024) #data is too large to handle all the years and causes at the same time. Select each to proceed
list_cause <- c("All causes", "Cardiovascular diseases", "Respiratory diseases", "Ischemic heart disease", "Tracheal, bronchus, and lung cancer") 
i <- 6
j <- 1
no2_ind <- subset(no2_tract_df, year==list_year[i])
db_ind <- subset(final3, year==list_year[i] & cause_name==list_cause[j])
final4 <- merge(no2_ind, db_ind, by="GEOID")
final4 <- merge(final4, xxxx, by="GEOID")

# Plotting
xxx <- merge(final4, ustract_df[,c("GEOID", "geometry")], by="GEOID")
q <- quantile(xxx$surface_no2, 0.95, na.rm=TRUE)
min_val <- min(xxx$surface_no2, na.rm = TRUE)
max_val <- max(xxx$surface_no2, na.rm = TRUE)

brks <- c(seq(0, q, length.out = 10), max_val)
bin_labels <- scales::number_format(accuracy = 0.1)(brks)
bin_labels[length(bin_labels) - 1] <- paste0(bin_labels[length(bin_labels) - 1])
bin_labels[length(bin_labels)] <- ""

xxx <- st_as_sf(xxx)
st_crs(xxx) <- 4326
states_sf <- st_as_sf(map("state", plot = FALSE, fill = TRUE))
states_sf <- st_transform(states_sf, st_crs(xxx))


pdf(file = paste0("/Users/k_sy_n_imac/Library/CloudStorage/Dropbox/GWU/Lab/FINESST/Research/Results/NO2_", list_year[i],"_final.pdf"), 
    width = 11, 
    height = 5)
par(family = "Lato")
ggplot() +  
  geom_sf(data = xxx, aes(geometry = geometry.x,fill = surface_no2), color = NA) +  
  geom_sf(data = states_sf, fill = NA, color = "#FAF9F6", linewidth = 0.15) +
  scale_fill_stepsn(
    colors = viridis::viridis(10, option = "plasma", begin = 0, end = 1),
    values = scales::rescale(brks),
    breaks = brks,                              # ✅ 11 breakpoint ticks
    labels = bin_labels,                        # ✅ labels at breakpoints
    limits = c(0, max_val),
    oob = scales::squish,
    na.value = "grey80",
    name = expression("NO"[2] ~ "(ppb)"),
    guide = guide_colorsteps(
      show.limits = TRUE,
      ticks = TRUE,
      even.steps = TRUE,                       # last bin can be wider
      label.position = "right",
      barwidth = unit(0.8, "cm"),
      barheight = unit(10, "cm")
    )
  ) +
  theme_minimal(base_size = 13) +
  #labs(title = bquote("Surface NO"[2] ~ " concentrations (CONUS census tracts, " * .(list_year[i]) * ")"),
  #     fill = expression("NO"[2] ~ "(ppb)")) +  
  theme( 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),  
    plot.title = element_text(size = 20, face = "bold"),  
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),  
    legend.text = element_text(size = 15),  
    legend.title = element_text(size = 20)  
  )  
dev.off()

#Descriptive
yy <- subset(xxx, xxx$surface_no2>=q)
yy <- subset(xxx, xxx$surface_no2>=5.31)
table(yy$metro3)
sum(!is.na(xxx$surface_no2))



# Spatial: Mortality burdens 
#Merge & Choose specific years/cause of death
list_year <- c(2019, 2020, 2021, 2022, 2023, 2024) #data is too large to handle all the years and causes at the same time. Select each to proceed
list_cause <- c("All causes", "Cardiovascular diseases", "Respiratory diseases", "Ischemic heart disease", "Tracheal, bronchus, and lung cancer") 
i <- 6
j <- 1

# Plotting
for (j in c(1:5)) {
  
  no2_ind <- subset(no2_tract_df, year==list_year[i])
  db_ind <- subset(final3, year==list_year[i] & cause_name==list_cause[j])
  final4 <- merge(no2_ind, db_ind, by="GEOID")
  final4 <- merge(final4, xxxx, by="GEOID")
  
  xxx <- merge(final4, ustract_df[,c("GEOID", "geometry")], by="GEOID")
  q <- quantile(xxx$rate_pe, 0.95, na.rm=TRUE)
  minv <- min(xxx$rate_pe, na.rm = TRUE)
  maxv <- max(xxx$rate_pe, na.rm = TRUE)
  brks <- c(seq(minv, q, length.out = 10), maxv)
  
  bin_labels <- scales::number_format(accuracy = 0.1)(brks)
  bin_labels[length(bin_labels) - 1] <- paste0(bin_labels[length(bin_labels) - 1])
  bin_labels[length(bin_labels)] <- ""
  
  xxx <- st_as_sf(xxx)
  st_crs(xxx) <- 4326
  states_sf <- st_as_sf(map("state", plot = FALSE, fill = TRUE))
  states_sf <- st_transform(states_sf, st_crs(xxx))
  
  pdf(file = paste0("/Users/k_sy_n_imac/Library/CloudStorage/Dropbox/GWU/Lab/FINESST/Research/Results/",
                    list_cause[j], "_", list_year[i], "_final.pdf"),
      width = 11, height = 5)
  par(family = "Lato")
  print(
    ggplot(xxx) +
      geom_sf(data = xxx, aes(geometry = geometry.x, fill = rate_pe), color = NA) +
      geom_sf(data = states_sf, fill = NA, color = "#FAF9F6", linewidth = 0.15) +
      scale_fill_stepsn(
        colors = viridis::viridis(10, option = "plasma", begin = 0, end = 1),
        values = scales::rescale(brks),
        breaks = brks,                              # ✅ 11 breakpoint ticks
        labels = bin_labels,                        # ✅ labels at breakpoints
        limits = c(minv, maxv),
        oob = scales::squish,
        na.value = "grey80",
        name = "Deaths per\n100,000",
        guide = guide_colorsteps(
          show.limits = TRUE,
          ticks = TRUE,
          even.steps = TRUE,                       # last bin can be wider
          label.position = "right",
          barwidth = unit(0.8, "cm"),
          barheight = unit(10, "cm")
        )
      ) +
      theme_minimal() +
      #labs(
     #   title = bquote(atop("NO"[2] * "-attributable mortality rates ",
      #                      "(" * .(list_cause[j]) * ", CONUS census tracts, " * .(list_year[i]) * ")"))
     # ) +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title   = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 15),  
        legend.title = element_text(size = 17),
        legend.spacing.y = unit(0, "cm")
      )
  )
  dev.off()
}

#Descriptive
a <- xxx %>% group_by(EPA_region) %>% summarize(x = sum(case_pe))
a <- xxx  %>% summarize(x = sum(case_pe), pop=sum(pop_over20))
a$x/a$pop * 100000

yy <- subset(xxx, xxx$metro3 %in% c("MSA_high", "MSA_low"))
sum(yy$case_pe)

sum(!is.na(xxx$rate_pe))


# Temporal: NO2
#Merge & Choose specific years/cause of death   **FOR TEMPORAL CHANGE PLOTS**
list_year <- c(2019, 2020, 2021, 2022, 2023, 2024) #data is too large to handle all the years and causes at the same time. Select each to proceed
list_cause <- c("All causes", "Cardiovascular diseases", "Respiratory diseases", "Ischemic heart disease", "Tracheal, bronchus, and lung cancer") 
i1 <- 1
i2 <- 6
j <- 1
no2_ind <- subset(no2_tract_df, year %in% c(list_year[i1], list_year[i2]))
db_ind <- subset(final3, year %in% c(list_year[i1], list_year[i2]) & cause_name==list_cause[j])
final4 <- merge(no2_ind, db_ind, by=c("GEOID", "year"))
final4 <- merge(final4, xxxx, by="GEOID", all.x=TRUE)

#Differences calculation
xxx <- merge(final4, ustract_df[,c("GEOID", "geometry")], by="GEOID")
x <- subset(final4, final4$year==list_year[i1])
x <- rename(x, no2_pre=surface_no2)
xx <- subset(final4, final4$year==list_year[i2])
xx <- rename(xx, no2_post=surface_no2)
xx <- xx[,c("GEOID", "no2_post")]

xxx <- merge(x, xx, by="GEOID", all=TRUE)
xxx$abs_dif <- - (xxx$no2_pre - xxx$no2_post)
xxx$rel_dif <- - (xxx$no2_pre - xxx$no2_post) / xxx$no2_pre * 100

# Descriptive
summary(xxx$abs_dif)
summary(xxx$rel_dif)

#Absolute dif
xxx <- merge(xxx, ustract_df[,c("GEOID", "geometry")], by="GEOID")
q99 <- quantile(xxx$abs_dif, 0.99, na.rm=TRUE)
q1 <- quantile(xxx$abs_dif, 0.01, na.rm=TRUE)
xxx$abs_dif <- ifelse(xxx$abs_dif>=q99, q99, xxx$abs_dif)
xxx$abs_dif <- ifelse(xxx$abs_dif<=q1, q1, xxx$abs_dif)
min_val <- min(xxx$abs_dif, na.rm = TRUE)
max_val <- max(xxx$abs_dif, na.rm = TRUE)
mid_val <- (min_val + max_val) / 2

xxx <- st_as_sf(xxx)
st_crs(xxx) <- 4326
states_sf <- st_as_sf(map("state", plot = FALSE, fill = TRUE))
states_sf <- st_transform(states_sf, st_crs(xxx))

pdf(file = paste0("/Users/k_sy_n_imac/Library/CloudStorage/Dropbox/GWU/Lab/FINESST/Research/Results/NO2_abs_dif_", list_year[i1], "_", list_year[i2],"_final.pdf"), 
    width = 11, 
    height = 5)
par(family = "Lato")
ggplot() + 
  geom_sf(data = xxx, aes(geometry = geometry.x, fill = abs_dif), color = NA) + 
  geom_sf(data = states_sf, fill = NA, color = "black", linewidth = 0.1) +
  scale_fill_gradientn(
    colors = c("#1034A6", "grey90", "#C21807"),  # Blue → Grey → Red
    values = scales::rescale(c(min_val, 0, max_val)),  # Ensure grey at 0
    name = expression("Δ" ~ NO[2] ~ "(ppb)"),
    na.value = "grey10",
    limits = c(min_val, max_val),
    breaks = c(min_val, 0, max_val),
    labels = function(x) {
      ifelse(x == max_val, paste0(">= ", number_format(accuracy = 0.1)(x)), 
             ifelse(x == min_val, paste0("<= ", number_format(accuracy = 0.1)(x)), 
                    number_format(accuracy = 0.1)(x))) }) +  
  theme_minimal() +  
  labs(
    #title = bquote("Changes in surface NO"[2] * " concentrations from " * .(list_year[i1]) * "-" * .(list_year[i2]) * " (CONUS census tracts)"),
    fill = expression("NO"[2] ~ "(ppb)")
  ) + 
  guides(fill = guide_colorbar(
    barwidth = unit(0.8, "cm"), barheight = unit(10, "cm"),
    label.position = "right",
    title.theme = element_text(size = 20, face = "bold", margin = margin(b = 20)),
    label.theme = element_text(size = 16)
  )) +   
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    legend.text = element_text(size = 17),
    legend.title = element_text(size = 20)
  )
dev.off()

#Descriptive
sum(!is.na(xxx$abs_dif))
pct_positive_by_region <- aggregate(
  abs_dif ~ EPA_region,
  data = xxx,
  FUN = function(x) mean(x > 0, na.rm = TRUE) * 100
)
pct_positive_by_region
mean_abs_dif_by_region <- aggregate(
  abs_dif ~ EPA_region,
  data = xxx,
  FUN = function(x) mean(x, na.rm = TRUE)
)
mean_abs_dif_by_region

# Temporal: mortality burdens
#Merge & Choose specific years/cause of death   
list_year <- c(2019, 2020, 2021, 2022, 2023, 2024) #data is too large to handle all the years and causes at the same time. Select each to proceed
list_cause <- c("All causes", "Cardiovascular diseases", "Respiratory diseases", "Ischemic heart disease", "Tracheal, bronchus, and lung cancer") 
i1 <- 1
i2 <- 6
j <- 1
no2_ind <- subset(no2_tract_df, year %in% c(list_year[i1], list_year[i2]))
db_ind <- subset(final3, year %in% c(list_year[i1], list_year[i2]) & cause_name==list_cause[j])
final4 <- merge(no2_ind, db_ind, by=c("GEOID", "year"))
final4 <- merge(final4, xxxx, by="GEOID", all.x=TRUE)

#Differences calculation
xxx <- merge(final4, ustract_df[,c("GEOID", "geometry")], by="GEOID")
x <- subset(final4, final4$year==list_year[i1])
x <- rename(x, no2_pre=rate_pe)
xx <- subset(final4, final4$year==list_year[i2])
xx <- rename(xx, no2_post=rate_pe)
xx <- xx[,c("GEOID", "no2_post")]

xxx <- merge(x, xx, by="GEOID", all=TRUE)
xxx$abs_dif <- - (xxx$no2_pre - xxx$no2_post)
xxx$rel_dif <- - (xxx$no2_pre - xxx$no2_post) / xxx$no2_pre * 100

# Descriptive
summary(xxx$abs_dif)
summary(xxx$rel_dif)

#Absolute dif
xxx <- merge(xxx, ustract_df[,c("GEOID", "geometry")], by="GEOID")
q99 <- quantile(xxx$abs_dif, 0.99, na.rm=TRUE)
q1 <- quantile(xxx$abs_dif, 0.01, na.rm=TRUE)
xxx$abs_dif <- ifelse(xxx$abs_dif>=q99, q99, xxx$abs_dif)
xxx$abs_dif <- ifelse(xxx$abs_dif<=q1, q1, xxx$abs_dif)
min_val <- min(xxx$abs_dif, na.rm = TRUE)
max_val <- max(xxx$abs_dif, na.rm = TRUE)
mid_val <- (min_val + max_val) / 2

xxx <- st_as_sf(xxx)
st_crs(xxx) <- 4326
states_sf <- st_as_sf(map("state", plot = FALSE, fill = TRUE))
states_sf <- st_transform(states_sf, st_crs(xxx))

pdf(file = paste0("/Users/k_sy_n_imac/Library/CloudStorage/Dropbox/GWU/Lab/FINESST/Research/Results/NO2burden_abs_dif_", list_year[i1], "_", list_year[i2],"_final.pdf"), 
    width = 11, 
    height = 5)
par(family = "Lato")
ggplot() + 
  geom_sf(data = xxx, aes(geometry = geometry.x, fill = abs_dif), color = NA) + 
  geom_sf(data = states_sf, fill = NA, color = "black", linewidth = 0.1) +
  scale_fill_gradientn(
    colors = c("#1034A6", "grey90", "#C21807"),  # Blue → Grey → Red
    values = scales::rescale(c(min_val, 0, max_val)),  # Ensure grey at 0
    name = "Δ Deaths\nper 100,000",
    na.value = "grey10",
    limits = c(min_val, max_val),
    breaks = c(min_val, 0, max_val),
    labels = function(x) {
      ifelse(x == max_val, paste0(">= ", number_format(accuracy = 0.1)(x)), 
             ifelse(x == min_val, paste0("<= ", number_format(accuracy = 0.1)(x)), 
                    number_format(accuracy = 0.1)(x)))
    }
  ) +  
  theme_minimal() +  
#  labs(
  #  title = bquote("Changes in NO"[2] * "-attributable premature deaths from " * .(list_year[i1]) * "-" * .(list_year[i2]) * " (CONUS census tracts)"),
  #  fill = expression("NO"[2] ~ "(ppb)")
#  ) + 
  guides(fill = guide_colorbar(
    barwidth = unit(0.8, "cm"), barheight = unit(9.5, "cm"),
    label.position = "right", title.hjust = 0,
    title.theme = element_text(size = 20, margin = margin(b = 20)),
    label.theme = element_text(size = 16)
  )) +   
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(size = 20, face = "bold"),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 17)
  )
dev.off()

#Descriptive
sum(!is.na(xxx$abs_dif))
summary(xxx$abs_dif)








# Sensitivity analysis ######################################################

# Spatial: Mortality burdens 
#Merge & Choose specific years/cause of death
final3_new <- read.csv("/Users/k_sy_n_imac/Library/CloudStorage/Dropbox/GWU/Lab/FINESST/Research/Data/final3_new.csv")
final3_new$rate_pe <- ifelse(is.na(final3_new$rate_pe), 0, final3_new$rate_pe)
list_year <- c(2019, 2020, 2021, 2022, 2023, 2024) #data is too large to handle all the years and causes at the same time. Select each to proceed
list_cause <- c("All causes", "Cardiovascular diseases", "Respiratory diseases", "Ischemic heart disease", "Tracheal, bronchus, and lung cancer") 
i <- 6
j <- 1

# Plotting
for (j in c(1:5)) {
  
  no2_ind <- subset(no2_tract_df, year==list_year[i])
  db_ind <- subset(final3_new, year==list_year[i] & cause_name=="all")
  final4 <- merge(no2_ind, db_ind, by="GEOID")
  final4 <- merge(final4, xxxx, by="GEOID")
  
  xxx <- merge(final4, ustract_df[,c("GEOID", "geometry")], by="GEOID")
  q <- quantile(xxx$rate_pe, 0.95, na.rm=TRUE)
  minv <- min(xxx$rate_pe, na.rm = TRUE)
  maxv <- max(xxx$rate_pe, na.rm = TRUE)
  brks <- c(seq(minv, q, length.out = 10), maxv)
  
  bin_labels <- scales::number_format(accuracy = 0.1)(brks)
  bin_labels[length(bin_labels) - 1] <- paste0(bin_labels[length(bin_labels) - 1])
  bin_labels[length(bin_labels)] <- ""
  
  pdf(file = paste0("/Users/k_sy_n_imac/Library/CloudStorage/Dropbox/GWU/Lab/FINESST/Research/Results/",
                    list_cause[j], "_", list_year[i], "_final.pdf"),
      width = 11, height = 5)
  par(family = "Lato")
  print(
    ggplot(xxx) +
      geom_sf(data = xxx, aes(geometry = geometry.x, fill = rate_pe), color = NA) +
      scale_fill_stepsn(
        colors = viridis::viridis(10, option = "plasma", begin = 0, end = 1),
        values = scales::rescale(brks),
        breaks = brks,                              # ✅ 11 breakpoint ticks
        labels = bin_labels,                        # ✅ labels at breakpoints
        limits = c(minv, maxv),
        oob = scales::squish,
        na.value = "grey80",
        name = "Deaths per\n100,000",
        guide = guide_colorsteps(
          show.limits = TRUE,
          ticks = TRUE,
          even.steps = TRUE,                       # last bin can be wider
          label.position = "right",
          barwidth = unit(0.8, "cm"),
          barheight = unit(10, "cm")
        )
      ) +
      theme_minimal() +
      #labs(
      #   title = bquote(atop("NO"[2] * "-attributable mortality rates ",
      #                      "(" * .(list_cause[j]) * ", CONUS census tracts, " * .(list_year[i]) * ")"))
      # ) +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title   = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 15),  
        legend.title = element_text(size = 17),
        legend.spacing.y = unit(0, "cm")
      )
  )
  dev.off()
}

#Descriptive
a <- xxx %>% group_by(EPA_region) %>% summarize(x = sum(case_pe))
a <- xxx  %>% summarize(x = sum(case_pe), pop=sum(pop_over20))
a$x/a$pop * 100000

yy <- subset(xxx, xxx$metro3 %in% c("MSA_high", "MSA_low"))
sum(yy$case_pe)

sum(!is.na(xxx$rate_pe))


# Temporal: mortality burdens
#Merge & Choose specific years/cause of death   
list_year <- c(2019, 2020, 2021, 2022, 2023, 2024) #data is too large to handle all the years and causes at the same time. Select each to proceed
list_cause <- c("All causes", "Cardiovascular diseases", "Respiratory diseases", "Ischemic heart disease", "Tracheal, bronchus, and lung cancer") 
i1 <- 1
i2 <- 6
j <- 1
no2_ind <- subset(no2_tract_df, year %in% c(list_year[i1], list_year[i2]))
db_ind <- subset(final3, year %in% c(list_year[i1], list_year[i2]) & cause_name==list_cause[j])
final4 <- merge(no2_ind, db_ind, by=c("GEOID", "year"))
final4 <- merge(final4, xxxx, by="GEOID", all.x=TRUE)

#Differences calculation
xxx <- merge(final4, ustract_df[,c("GEOID", "geometry")], by="GEOID")
x <- subset(final4, final4$year==list_year[i1])
x <- rename(x, no2_pre=rate_pe)
xx <- subset(final4, final4$year==list_year[i2])
xx <- rename(xx, no2_post=rate_pe)
xx <- xx[,c("GEOID", "no2_post")]

xxx <- merge(x, xx, by="GEOID", all=TRUE)
xxx$abs_dif <- - (xxx$no2_pre - xxx$no2_post)
xxx$rel_dif <- - (xxx$no2_pre - xxx$no2_post) / xxx$no2_pre * 100

# Descriptive
summary(xxx$abs_dif)
summary(xxx$rel_dif)

#Absolute dif
xxx <- merge(xxx, ustract_df[,c("GEOID", "geometry")], by="GEOID")
q99 <- quantile(xxx$abs_dif, 0.99, na.rm=TRUE)
q1 <- quantile(xxx$abs_dif, 0.01, na.rm=TRUE)
xxx$abs_dif <- ifelse(xxx$abs_dif>=q99, q99, xxx$abs_dif)
xxx$abs_dif <- ifelse(xxx$abs_dif<=q1, q1, xxx$abs_dif)
min_val <- min(xxx$abs_dif, na.rm = TRUE)
max_val <- max(xxx$abs_dif, na.rm = TRUE)
mid_val <- (min_val + max_val) / 2

pdf(file = paste0("/Users/k_sy_n_imac/Library/CloudStorage/Dropbox/GWU/Lab/FINESST/Research/Results/NO2burden_abs_dif_", list_year[i1], "_", list_year[i2],"_final.pdf"), 
    width = 11, 
    height = 5)
par(family = "Lato")
ggplot() + 
  geom_sf(data = xxx, aes(geometry = geometry.x, fill = abs_dif), color = NA) + 
  scale_fill_gradientn(
    colors = c("#1034A6", "grey90", "#C21807"),  # Blue → Grey → Red
    values = scales::rescale(c(min_val, 0, max_val)),  # Ensure grey at 0
    name = "Δ Deaths\nper 100,000",
    na.value = "grey10",
    limits = c(min_val, max_val),
    breaks = c(min_val, 0, max_val),
    labels = function(x) {
      ifelse(x == max_val, paste0(">= ", number_format(accuracy = 0.1)(x)), 
             ifelse(x == min_val, paste0("<= ", number_format(accuracy = 0.1)(x)), 
                    number_format(accuracy = 0.1)(x)))
    }
  ) +  
  theme_minimal() +  
  #  labs(
  #  title = bquote("Changes in NO"[2] * "-attributable premature deaths from " * .(list_year[i1]) * "-" * .(list_year[i2]) * " (CONUS census tracts)"),
  #  fill = expression("NO"[2] ~ "(ppb)")
  #  ) + 
  guides(fill = guide_colorbar(
    barwidth = unit(0.8, "cm"), barheight = unit(9.5, "cm"),
    label.position = "right", title.hjust = 0,
    title.theme = element_text(size = 20, margin = margin(b = 20)),
    label.theme = element_text(size = 16)
  )) +   
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(size = 20, face = "bold"),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 17)
  )
dev.off()

#Descriptive
sum(!is.na(xxx$abs_dif))
summary(xxx$abs_dif)












# Figure S2: time-series plots 
# Merge & Choose specific years/cause of death
list_year <- c(2019, 2020, 2021, 2022, 2023, 2024) #data is too large to handle all the years and causes at the same time. Select each to proceed
list_cause <- c("All causes", "Cardiovascular diseases", "Respiratory diseases", "Ischemic heart disease", "Tracheal, bronchus, and lung cancer") 
i <- 6
j <- 1
db_ind <- subset(final3, cause_name==list_cause[j])
final4 <- merge(no2_tract_df, db_ind, by=c("year","GEOID"))
final4 <- merge(final4, xxxx, by="GEOID", all.x=TRUE)

region_labels <- c(
  "Region 1 (New England)",
  "Region 2 (Northeast/Caribbean)",
  "Region 3 (Mid-Atlantic)",
  "Region 4 (Southeast)",
  "Region 5 (Great Lakes)",
  "Region 6 (South Central)",
  "Region 7 (Heartland)",
  "Region 8 (Mountains/Plains)",
  "Region 9 (Pacific Southwest)",
  "Region 10 (Northwest)"
)

recode_region <- function(x) {
  recode(x,
         "Region 1"  = "Region 1 (New England)",
         "Region 2"  = "Region 2 (Northeast/Caribbean)",
         "Region 3"  = "Region 3 (Mid-Atlantic)",
         "Region 4"  = "Region 4 (Southeast)",
         "Region 5"  = "Region 5 (Great Lakes)",
         "Region 6"  = "Region 6 (South Central)",
         "Region 7"  = "Region 7 (Heartland)",
         "Region 8"  = "Region 8 (Mountains/Plains)",
         "Region 9"  = "Region 9 (Pacific Southwest)",
         "Region 10" = "Region 10 (Northwest)"
  )
}
final4 <- final4 %>%
  mutate(EPA_region = factor(recode_region(EPA_region), levels = region_labels))

# NO2
nation_data <- final4 %>%
  group_by(year) %>%
  summarise(
    pop_weighted_no2 = sum(surface_no2 * total_dc2020, na.rm = TRUE) /
      sum(total_dc2020, na.rm = TRUE),
    .groups = "drop"
  )

mean_data <- final4 %>%
  group_by(EPA_region, metro4, year) %>%
  summarise(
    pop_weighted_no2 = sum(surface_no2 * total_dc2020, na.rm = TRUE) /
      sum(total_dc2020, na.rm = TRUE),
    case_pe = sum(case_pe, na.rm = TRUE),
    pop     = sum(total_dc2020, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(metro4 = factor(metro4, levels = c("Large MSA", "Small MSA",
                                            "MicroSA", "Non-urban")))
# EPA_region already a named factor from final4 — no need to re-mutate

metro_colors <- c(
  "Large MSA" = "#2E004F",
  "Small MSA" = "#54278F",
  "MicroSA"   = "#756BB1",
  "Non-urban" = "#CBC9E2"
)


pdf(file = paste0("/Users/k_sy_n_imac/Library/CloudStorage/Dropbox/GWU/Lab/FINESST/Research/Results/time_series_NO2_final.pdf"), 
    width = 20, 
    height = 7)
par(family = "Lato")
ggplot(mean_data, aes(x = year, y = pop_weighted_no2,
                      group = metro4, color = metro4)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 1) +
  geom_line(
    data = nation_data,
    aes(x = year, y = pop_weighted_no2, color = "National"),
    inherit.aes = FALSE,
    linewidth = 1.1
  ) +
  facet_wrap(~ EPA_region, ncol = 5) +
  scale_color_manual(
    values = c(
      metro_colors,
      "National" = "#B11226"
    ),
    name = "Urbanicity"
  ) +
  labs(
    x = "Year",
    y = bquote("Population-weighted average NO"[2] ~ "(ppb)")
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 16, margin = margin(t = 15)),
    axis.title.y = element_text(size = 16,margin = margin(r = 15)),    
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    legend.position = "bottom"
  )
dev.off()

# Urbanicity only
urbanicity_data <- final4 %>%
  group_by(metro4, year) %>%
  summarise(
    pop_weighted_no2 = sum(surface_no2 * total_dc2020, na.rm = TRUE) /
      sum(total_dc2020, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    metro4 = factor(metro4, levels = c("Large MSA", "Small MSA",
                                       "MicroSA", "Non-urban"))
  )

pdf(
  file = "/Users/k_sy_n_imac/Library/CloudStorage/Dropbox/GWU/Lab/FINESST/Research/Results/time_series_NO2_by_urbanicity.pdf",
  width = 10,
  height = 7
)

par(family = "Lato")

ggplot(urbanicity_data, aes(x = year, y = pop_weighted_no2,
                            group = metro4, color = metro4)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  geom_line(
    data = nation_data,
    aes(x = year, y = pop_weighted_no2, color = "National"),
    inherit.aes = FALSE,
    linewidth = 1.1
  ) +
  scale_color_manual(
    values = c(
      metro_colors,
      "National" = "#B11226"
    ),
    name = "Urbanicity"
  ) +
  labs(
    x = "Year",
    y = bquote("Population-weighted average NO"[2] ~ "(ppb)")
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 18, margin = margin(t = 15)),
    axis.title.y = element_text(size = 18, margin = margin(r = 15)),
    axis.text = element_text(size = 17),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 18),
    legend.position = "bottom"
  )
dev.off()

# Region only
region_colors <- c(
  "Region 1 (New England)"         = "#1B2F70",
  "Region 2 (Northeast/Caribbean)" = "#2C4FD7",
  "Region 3 (Mid-Atlantic)"        = "#5A8AC6",
  "Region 4 (Southeast)"           = "#8BC8C8",
  "Region 5 (Great Lakes)"         = "#A1D99B",
  "Region 6 (South Central)"       = "#D9C89B",
  "Region 7 (Heartland)"           = "#DDAA66",
  "Region 8 (Mountains/Plains)"    = "#F2A6A0",
  "Region 9 (Pacific Southwest)"   = "#7A5AC9",
  "Region 10 (Northwest)"          = "#4A0066"
)

region_data <- final4 %>%
  group_by(EPA_region, year) %>%
  summarise(
    pop_weighted_no2 = sum(surface_no2 * total_dc2020, na.rm = TRUE) /
      sum(total_dc2020, na.rm = TRUE),
    .groups = "drop"
  )

pdf(
  file = "/Users/k_sy_n_imac/Library/CloudStorage/Dropbox/GWU/Lab/FINESST/Research/Results/time_series_NO2_by_region.pdf",
  width = 15,
  height = 8)
par(family = "Lato")

ggplot(
  region_data,
  aes(x = year,
      y = pop_weighted_no2,
      group = EPA_region,
      color = EPA_region)
) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  
  # National summary
  geom_line(
    data = nation_data,
    aes(x = year,
        y = pop_weighted_no2,
        color = "National"),
    inherit.aes = FALSE,
    linewidth = 1.3
  ) +
  
  scale_color_manual(
    values = c(
      region_colors,
      "National" = "#B11226"
    ),
    name = "EPA region"
  ) +
  
  labs(
    x = "Year",
    y = bquote("Population-weighted average NO"[2] ~ "(ppb)")
  ) +
  
  theme_minimal() +
  
  theme(
    axis.title.x = element_text(size = 18, margin = margin(t = 15)),
    axis.title.y = element_text(size = 18, margin = margin(r = 15)),
    axis.text = element_text(size = 17),
    legend.title = element_text(size = 17),
    legend.text = element_text(size = 17),
    legend.position = "bottom"
  )
dev.off()







##########################################################################################
# Figure S3 ##############################################################################
##########################################################################################
# Mortality burden
mean_data <- final4 %>%
  group_by(EPA_region, metro4, year) %>%
  summarise(
    case_pe = sum(case_pe, na.rm = TRUE),
    pop     = sum(pop_over20, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    metro4 = factor(
      metro4,
      levels = c("Large MSA", "Small MSA",
                 "MicroSA", "Non-urban")
    )
  )

regional_pop <- final4 %>%
  group_by(EPA_region, year) %>%
  summarise(
    regional_pop = sum(pop_over20, na.rm = TRUE),
    .groups = "drop"
  )

mean_data <- merge(
  mean_data,
  regional_pop,
  by = c("EPA_region", "year")
) %>%
  mutate(rate = case_pe / regional_pop * 100000)


# National mean
national_data <- final4 %>%
  group_by(year) %>%
  summarise(
    case_pe = sum(case_pe, na.rm = TRUE),
    pop     = sum(pop_over20, na.rm = TRUE),
    .groups = "drop") %>%
  mutate(rate = case_pe / pop * 100000)


# Regional baseline mortality rate
mor <- read.csv(
  "/Users/k_sy_n_imac/Library/CloudStorage/Dropbox/GWU/Lab/FINESST/Research/Data/IHME-GBD_2023/IHME-GBD_2023_DATA-58df5c3a-1/IHME-GBD_2023_DATA-58df5c3a-1.csv",
  header = TRUE)
mor <- rename(mor,
              state_name = location_name,
              mor_rate = val)
mor <- subset(mor, mor$metric_name == "Rate")
mor <- mor[c(
  "state_name",
  "age_name",
  "cause_name",
  "year",
  "mor_rate"
)]
xx <- subset(mor, mor$year == 2023)
xx$year <- 2024
mor <- rbind(mor, xx)

# Population
pop <- final4 %>%
  filter(year == 2019) %>%
  distinct(state_name.x, .keep_all = TRUE)

pop <- pop[, c(
  "state_name.x",
  "pop_over20",
  "EPA_region"
)]

pop <- rename(pop,
              state_name = state_name.x)

merge <- merge(mor, pop, by = "state_name")

mor_region <- merge %>%
  group_by(EPA_region, cause_name, year) %>%
  summarise(
    mor_rate = sum(mor_rate * pop_over20, na.rm = TRUE) /
      sum(pop_over20, na.rm = TRUE),
    .groups = "drop"
  )


pdf(file = "/Users/k_sy_n_imac/Library/CloudStorage/Dropbox/GWU/Lab/FINESST/Research/Results/time_series_NO2burden_final.pdf",
  width = 20, height = 7)

ggplot(
  mean_data,
  aes(
    x = year,
    y = rate,
    group = metro4,
    color = metro4)) +
  
  # Transparent blue bars
  geom_col(
    data = mor_region,
    aes(
      x = year,
      y = mor_rate,
      fill = "Regional baseline mortality rate (over 20)"
    ),
    inherit.aes = FALSE,
    alpha = 0.22,
    width = 0.5
  ) +
  
  geom_line(linewidth = 1.2) +
  geom_point(size = 1.5) +
  facet_wrap(~ EPA_region, ncol = 5) +
  
  geom_line(
    data = national_data,
    aes(
      x = year,
      y = rate,
      group = 1,
      color = "National"
    ),
    linewidth = 1.1,
    inherit.aes = FALSE
  ) +
  
  scale_color_manual(
    values = c(
      metro_colors,
      "National" = "#B11226"
    ),
    breaks = c(
      names(metro_colors),
      "National"
    ),
    name = NULL
  ) +
  
  scale_fill_manual(
    values = c(
      "Regional baseline mortality rate (over 20)" = "grey"
    ),
    breaks = "Regional baseline mortality rate (over 20)",
    name = NULL
  ) +
  scale_x_continuous(
    breaks = 2019:2024
  ) +
  labs(
    x = "Year",
    y = bquote(
      "Premature cases attributable to NO"[2] ~
        "(per 100,000)"
    )
  ) +
  theme_minimal(base_family = "Lato") +
  theme(
    strip.text = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 16, margin = margin(t = 14)),
    axis.title.y = element_text(size = 16, margin = margin(r = 14)),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    legend.position = "bottom"
  )
dev.off()



# Urbanicity only
urbanicity_burden <- final4 %>%
  group_by(metro4, year) %>%
  summarise(
    case_pe = sum(case_pe, na.rm = TRUE),
    pop     = sum(pop_over20, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    rate = case_pe / pop * 100000,
    metro4 = factor(
      metro4,
      levels = c("Large MSA", "Small MSA",
                 "MicroSA", "Non-urban")
    )
  )
# National mean
national_data <- final4 %>%
  group_by(year) %>%
  summarise(
    case_pe = sum(case_pe, na.rm = TRUE),
    pop     = sum(pop_over20, na.rm = TRUE),
    .groups = "drop") %>%
  mutate(rate = case_pe / pop * 100000)

# Plot
pdf(file = "/Users/k_sy_n_imac/Library/CloudStorage/Dropbox/GWU/Lab/FINESST/Research/Results/time_series_NO2burden_by_urbanicity.pdf",
  width = 10, height = 7)

ggplot(
  urbanicity_burden,
  aes(
    x = year,
    y = rate,
    group = metro4,
    color = metro4
  )
) +
  
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  geom_line(
    data = national_data,
    aes(
      x = year,
      y = rate,
      group = 1,
      color = "National"
    ),
    linewidth = 1.1,
    inherit.aes = FALSE
  ) +
  
  scale_color_manual(
    values = c(
      metro_colors,
      "National" = "#B11226"
    ),
    name = "Urbanicity"
  ) +
  
  labs(
    x = "Year",
    y = bquote(
      "Premature cases attributable to NO"[2] ~
        "(per 100,000)"
    )
  ) +
  
  theme_minimal(base_family = "Lato") +
  
  theme(
    axis.title.x = element_text(size = 18, margin = margin(t = 15)),
    axis.title.y = element_text(size = 18, margin = margin(r = 15)),
    axis.text = element_text(size = 17),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 18),
    legend.position = "bottom"
  )
dev.off()


# Region only
region_burden <- final4 %>%
  group_by(EPA_region, year) %>%
  summarise(
    case_pe = sum(case_pe, na.rm = TRUE),
    pop     = sum(total_dc2020, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(rate = case_pe / pop * 100000)

pdf(
  file = "/Users/k_sy_n_imac/Library/CloudStorage/Dropbox/GWU/Lab/FINESST/Research/Results/time_series_NO2burden_by_region.pdf",
  width = 15,
  height = 7
)

ggplot(
  region_burden,
  aes(
    x = year,
    y = rate,
    group = EPA_region,
    color = EPA_region
  )
) +
  
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  
  geom_line(
    data = national_data,
    aes(
      x = year,
      y = rate,
      group = 1,
      color = "National"
    ),
    linewidth = 1.1,
    inherit.aes = FALSE
  ) +
  
  scale_color_manual(
    values = c(
      region_colors,
      "National" = "#B11226"
    ),
    name = "EPA region"
  ) +
  
  labs(
    x = "Year",
    y = bquote(
      "Premature cases attributable to NO"[2] ~
        "(per 100,000)"
    )
  ) +
  
  theme_minimal(base_family = "Lato") +
  
  theme(
    axis.title.x = element_text(size = 17, margin = margin(t = 15)),
    axis.title.y = element_text(size = 15, margin = margin(r = 15)),
    axis.text = element_text(size = 17),
    legend.title = element_text(size = 17),
    legend.text = element_text(size = 17),
    legend.position = "bottom"
  )
dev.off()



# Descriptive
zero_both <- final4 %>%
  sf::st_drop_geometry() %>%
  filter(year %in% c(2019, 2024)) %>%
  group_by(GEOID, metro4) %>%
  summarise(
    zero_2019 = any(year == 2019 & case_pe == 0, na.rm = TRUE),
    zero_2024 = any(year == 2024 & case_pe == 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(zero_both = zero_2019 & zero_2024)
zero_both_by_metro <- zero_both %>%
  group_by(metro4) %>%
  summarise(
    n_total = n(),
    n_zero_both = sum(zero_both),
    pct_zero_both = 100 * n_zero_both / n_total,
    .groups = "drop"
  )

















########################################################################################
########################################################################################
# Sub-group analysis for trends  #######################################################
########################################################################################
########################################################################################
# Figure 2 ###############################################################################
list_year <- c(2019, 2020, 2021, 2022, 2023, 2024) #data is too large to handle all the years and causes at the same time. Select each to proceed
list_cause <- c("All causes", "Cardiovascular diseases", "Respiratory diseases", "Ischemic heart disease", "Tracheal, bronchus, and lung cancer") 
i1 <- 1
i2 <- 6
j <- 1
no2_ind <- subset(no2_tract_df, year %in% c(list_year[i1], list_year[i2]))
db_ind <- subset(final3, year %in% c(list_year[i1], list_year[i2]) & cause_name==list_cause[j])
final4 <- merge(no2_ind, db_ind, by=c("GEOID", "year"))
final4 <- merge(final4, xxxx, by="GEOID", all.x=TRUE)
final4 <- final4 %>%
  mutate(
    MSA_title = if_else(
      metro4 == "Non-urban",
      paste0("Non_urban_", state_name.x),
      MSA_title
    )
  )
group_df <- final4 %>% group_by(MSA_title, metro4, EPA_region) %>%
  summarize(n_2024 = n(), .groups = "drop")
msa_regions <- group_df %>% count(MSA_title) %>% filter(n > 1)
majority_regions <- final4 %>%
  group_by(MSA_title, EPA_region) %>%
  summarise(total_pop = sum(pop_over20, na.rm = TRUE), .groups = "drop") %>%
  group_by(MSA_title) %>%
  slice_max(total_pop, n = 1) %>%  # ← Gets region with HIGHEST population
  ungroup() %>%
  dplyr::select(MSA_title, EPA_region_new = EPA_region)
final4 <- final4 %>%
  left_join(majority_regions, by = "MSA_title") %>%
  mutate(EPA_region = if_else(!is.na(EPA_region_new), 
                              EPA_region_new, 
                              EPA_region))

# Descriptive
x <- subset(final4, final4$year==list_year[i1])
test <- x %>%
  group_by(EPA_region, metro4) %>%
  summarise(n_MSA = n_distinct(MSA_title), .groups = "drop") %>%
  pivot_wider(names_from = metro4, values_from = n_MSA, values_fill = 0)

# Calculation
# Sub-group
# 2019
x <- subset(final4, final4$year==list_year[i1])
group_df_2019 <- x %>%
  group_by(MSA_title, EPA_region, metro4) %>%
  summarize(
    pw_no2_2019 = weighted.mean(surface_no2, total_dc2020, na.rm = TRUE),
    .groups = "drop"
  )

# 2024
xx <- subset(final4, final4$year==list_year[i2])
group_df_2024 <- xx %>%
  group_by(MSA_title, EPA_region, metro4) %>%
  summarize(
    pw_no2_2024 = weighted.mean(surface_no2, total_dc2020, na.rm = TRUE),
    .groups = "drop"
  )

group_df_merge <- merge(
  group_df_2019,
  group_df_2024,
  by = c("MSA_title", "EPA_region", "metro4")
)

group_df_merge$pct_change <-
  (group_df_merge$pw_no2_2024 - group_df_merge$pw_no2_2019) /
  group_df_merge$pw_no2_2024 * 100


################################################################################
# Prepare raw (MSA-level) data
################################################################################
region_labels <- c(
  " Region 1\n (New England)",
  " Region 2\n (Northeast/Caribbean)",
  " Region 3\n (Mid-Atlantic)",
  " Region 4\n (Southeast)",
  " Region 5\n (Great Lakes)",
  " Region 6\n (South Central)",
  " Region 7\n (Heartland)",
  " Region 8\n (Mountains/Plains)",
  " Region 9\n (Pacific Southwest)",
  " Region 10\n (Northwest)")

recode_region <- function(x) {
  recode(x,
         "Region 1"  = " Region 1\n (New England)",
         "Region 2"  = " Region 2\n (Northeast/Caribbean)",
         "Region 3"  = " Region 3\n (Mid-Atlantic)",
         "Region 4"  = " Region 4\n (Southeast)",
         "Region 5"  = " Region 5\n (Great Lakes)",
         "Region 6"  = " Region 6\n (South Central)",
         "Region 7"  = " Region 7\n (Heartland)",
         "Region 8"  = " Region 8\n (Mountains/Plains)",
         "Region 9"  = " Region 9\n (Pacific Southwest)",
         "Region 10" = " Region 10\n (Northwest)")
}

region_colors <- c(
  " Region 1\n (New England)"         = "#1B2F70",
  " Region 2\n (Northeast/Caribbean)" = "#2C4FD7",
  " Region 3\n (Mid-Atlantic)"        = "#5A8AC6",
  " Region 4\n (Southeast)"           = "#8BC8C8",
  " Region 5\n (Great Lakes)"         = "#A1D99B",
  " Region 6\n (South Central)"       = "#D9C89B",
  " Region 7\n (Heartland)"           = "#DDAA66",
  " Region 8\n (Mountains/Plains)"    = "#F2A6A0",
  " Region 9\n (Pacific Southwest)"   = "#7A5AC9",
  " Region 10\n (Northwest)"          = "#4A0066")

plot_df <- group_df_merge %>%
  mutate(
    metro4     = factor(metro4, levels = c("Large MSA", "Small MSA",
                                           "MicroSA", "Non-urban")),
    EPA_region = factor(recode_region(EPA_region), levels = region_labels)
  )

################################################################################
# Plot: raw scatter
################################################################################

pdf(
  file = "/Users/k_sy_n_imac/Library/CloudStorage/Dropbox/GWU/Lab/FINESST/Research/Results/Quadrant_NO2_raw.pdf",
  width = 13,
  height = 10
)
par(family = "Lato")

ggplot(plot_df, aes(x = pct_change, y = pw_no2_2024,
                    color = EPA_region, shape = metro4)) +
  
  # Reference lines
  geom_hline(yintercept = 5.31, linetype = "solid", color = "black", linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.8) +
  
  geom_vline(xintercept = -5, linetype = "dashed", color = "gray40", linewidth = 0.5) +
  geom_vline(xintercept = 5, linetype = "dashed", color = "gray40", linewidth = 0.5) +
  
  # Label
  annotate("text",
           x = 5, y = 5.1,
           label = "5",
           size = 5,
           #fontface = "bold",
           color = "black") +
  annotate("text",
           x = -5, y = 5.1,
           label = "-5",
           size = 5,
           #fontface = "bold",
           color = "black") +
  annotate("text",
           x = -31, y = 5.6,
           label = "WHO AQG = 5.3 ppb",
           size = 5,
           #fontface = "bold",
           color = "black") +
  
  # Raw points (one per MSA)
  geom_point(size = 2.5, alpha = 0.8) +
  
  scale_color_manual(values = region_colors) +
  scale_shape_manual(values = c(
    "Large MSA" = 16,
    "Small MSA" = 17,
    "MicroSA" = 15,
    "Non-urban" = 8
  )) +
  labs(
    x = bquote(bold("% Change in population-weighted NO"[bolditalic(2)] ~ "between 2019 and 2024")),
    y = bquote(bold("Population-weighted NO"[bolditalic(2)] ~ "in 2024 (ppb)")),
    shape = "Urbanicity",
    color = "EPA Region"
  ) + 
  guides(
    shape = guide_legend(
      order = 1,
      nrow = 4,
      override.aes = list(size = 3)
    ),
    color = guide_legend(
      order = 2,
      nrow = 10,
      theme = theme(
        legend.text = element_text(
          size = 13,
          margin = margin(b = 8)
        )
      )
    )
  ) +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 13),
    axis.text = element_text(size = 13),
    axis.title.x = element_text(size = 15, margin = margin(t = 15, b = 15)),
    axis.title.y = element_text(size = 15, margin = margin(r = 15, l = 15)),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 13),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray95"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )
dev.off()


# Descriptive
sum(plot_df$pw_no2_2024 > 5.3, na.rm = TRUE)
mean(plot_df$pw_no2_2024 < 5.3, na.rm = TRUE) * 100

test <- subset(plot_df, plot_df$pw_no2_2024 >= 5.3)
table(test$metro4)

test <- subset(plot_df, plot_df$metro4=="Non-urban")

plot_df %>% group_by(metro4) %>% summarize(mean=mean(pw_no2_2024))
plot_df %>% group_by(metro4) %>% summarize(mean=mean(pct_change))
plot_df %>% group_by(EPA_region) %>% summarize(mean=mean(pw_no2_2024))
plot_df %>% group_by(EPA_region) %>% summarize(mean=mean(pct_change))

sum(abs(plot_df$pct_change) > 5, na.rm = TRUE)
mean(abs(plot_df$pct_change) > 5, na.rm = TRUE) * 100
sum(plot_df$pct_change < -5, na.rm = TRUE)
sum(plot_df$pct_change > 5, na.rm = TRUE)

test <- subset(plot_df, plot_df$pct_change < -5)
table(test$metro4)
test <- subset(plot_df, plot_df$pct_change > 5)
table(test$metro4)
table(test$EPA_region)
table(test$EPA_region, test$metro4)



################################################################################# 
# Summary by EPA REGION ONLY (across all urbanicity types)
#############################################################################
group_summary_by_region <- group_df_merge %>%
  group_by(EPA_region) %>%
  summarize(
    n_msa = n(),
    pw_no2 = mean(pw_no2_2024, na.rm = TRUE),
    var_pw_no2 = var(pw_no2_2024, na.rm = TRUE),
    pct_chg = mean(pct_change, na.rm = TRUE),
    var_pct_change = var(pct_change, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    se_pw_no2 = sqrt(var_pw_no2 / n_msa),
    se_pct_chg = sqrt(var_pct_change / n_msa),
    x0 = pct_chg,
    y0 = pw_no2,
    rx = 1.96 * se_pct_chg,
    ry = 1.96 * se_pw_no2
  ) %>%
  mutate(
    EPA_region = factor(
      EPA_region,
      levels = c("Region 1", "Region 2", "Region 3", "Region 4", "Region 5",
                 "Region 6", "Region 7", "Region 8", "Region 9", "Region 10")
    )
  )

nation_df <- group_df_merge %>%
  summarize(
    n_msa = n(),
    pw_no2 = mean(pw_no2_2024, na.rm = TRUE),
    var_pw_no2 = var(pw_no2_2024, na.rm = TRUE),
    pct_chg = mean(pct_change, na.rm = TRUE),
    var_pct_change = var(pct_change, na.rm = TRUE),
    .groups = "drop"
  )  %>%
  mutate(
    se_pw_no2 = sqrt(var_pw_no2 / n_msa),
    se_pct_chg = sqrt(var_pct_change / n_msa),
    x0 = pct_chg,
    y0 = pw_no2,
    rx = 1.96 * se_pct_chg,
    ry = 1.96 * se_pw_no2
  )
  
  
# Horizontal thin diamonds (x-axis 95% CI - replacing horizontal error bars)
diamond_region_x <- bind_rows(lapply(seq_len(nrow(group_summary_by_region)), function(i) {
  tibble(
    EPA_region = group_summary_by_region$EPA_region[i],
    diamond_id = paste0(group_summary_by_region$EPA_region[i], "_x"),
    x = c(group_summary_by_region$x0[i] - group_summary_by_region$rx[i],
          group_summary_by_region$x0[i],
          group_summary_by_region$x0[i] + group_summary_by_region$rx[i],
          group_summary_by_region$x0[i]),
    y = c(group_summary_by_region$y0[i],
          group_summary_by_region$y0[i] + 0.025,
          group_summary_by_region$y0[i],
          group_summary_by_region$y0[i] - 0.025)
  )
}))

# Vertical thin diamonds (y-axis 95% CI - replacing vertical error bars)
diamond_region_y <- bind_rows(lapply(seq_len(nrow(group_summary_by_region)), function(i) {
  tibble(
    EPA_region = group_summary_by_region$EPA_region[i],
    diamond_id = paste0(group_summary_by_region$EPA_region[i], "_y"),
    x = c(group_summary_by_region$x0[i],
          group_summary_by_region$x0[i] + 0.09,
          group_summary_by_region$x0[i],
          group_summary_by_region$x0[i] - 0.09),
    y = c(group_summary_by_region$y0[i] - group_summary_by_region$ry[i],
          group_summary_by_region$y0[i],
          group_summary_by_region$y0[i] + group_summary_by_region$ry[i],
          group_summary_by_region$y0[i])
  )
}))

diamond_region <- bind_rows(diamond_region_x, diamond_region_y)

region_colors <- c(
  "Region 1"  = "#1B2F70",
  "Region 2"  = "#2C4FD7",
  "Region 3"  = "#5A8AC6",
  "Region 4"  = "#8BC8C8",
  "Region 5"  = "#A1D99B",
  "Region 6"  = "#D9C89B",
  "Region 7"  = "#DDAA66",
  "Region 8"  = "#F2A6A0",
  "Region 9"  = "#7A5AC9",
  "Region 10" = "#4A0066"
)

pdf(file = "/Users/k_sy_n_imac/Library/CloudStorage/Dropbox/GWU/Lab/FINESST/Research/Results/Quadrant_NO2_EPA_Region.pdf",
    width = 9, height = 9)

ggplot() +
  # Reference lines
  geom_hline(yintercept = 5.31, linetype = "solid", color = "black", linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.5) +
  
  geom_vline(xintercept = -5, linetype = "dashed", color = "gray40", linewidth = 0.5) +
  geom_vline(xintercept = 5, linetype = "dashed", color = "gray40", linewidth = 0.5) +
  annotate("text",
           x = -9.8, y = 5.1,
           label = "WHO AQG = 5.3 ppb",
           size = 6.5,
           #fontface = "bold",
           color = "black") +
  
  # Thin diamond uncertainty regions (replacing error bars) - x and y axis separately
  geom_polygon(data = diamond_region, 
               aes(x = x, y = y, fill = EPA_region, group = diamond_id),
               alpha = 1, linewidth = 1, color = NA) +
  
  # Center points
  geom_point(data = group_summary_by_region,
             aes(x = x0, y = y0, color = EPA_region),
             size = 5, stroke = 1.2) +
  
  # National summary diamonds (uncertainty)
  geom_polygon(data = {
    h_diamond <- tibble(
      x = c(nation_df$pct_chg - nation_df$rx,
            nation_df$pct_chg,
            nation_df$pct_chg + nation_df$rx,
            nation_df$pct_chg),
      y = c(nation_df$pw_no2,
            nation_df$pw_no2 + 0.03,
            nation_df$pw_no2,
            nation_df$pw_no2 - 0.03),
      diamond_id = "national_x"
    )
    v_diamond <- tibble(
      x = c(nation_df$pct_chg,
            nation_df$pct_chg + 0.08,
            nation_df$pct_chg,
            nation_df$pct_chg - 0.08),
      y = c(nation_df$pw_no2 - nation_df$ry,
            nation_df$pw_no2,
            nation_df$pw_no2 + nation_df$ry,
            nation_df$pw_no2),
      diamond_id = "national_y"
    )
    bind_rows(h_diamond, v_diamond)
  },
  aes(x = x, y = y, group = diamond_id),
  alpha = 0.7, linewidth = 1, color = NA, fill = "#B11226") +
  
  # National summary point
  geom_point(data = nation_df,
             aes(x = pct_chg, y = pw_no2),
             color = "#B11226", size = 7, shape = 18) +
  geom_text(data = nation_df,
            aes(x = pct_chg, y = pw_no2),
            label = "National Summary",
            color = "#B11226", size = 7.5, fontface = "bold",
            vjust = -1.6, hjust = 0.05) +
  
  scale_color_manual(values = region_colors) +
  scale_fill_manual(values = region_colors) +
  
  labs(
    x = bquote(bold("% Change from 2019")),
    y = bquote(bold("NO"[bolditalic(2)] ~ "in 2024 (ppb)")),
    color = "EPA Region",
    fill = "EPA Region"
  ) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 22),
    axis.title.x = element_text(size = 25, margin = margin(t = 15, b = 15)),
    axis.title.y = element_text(size = 25, margin = margin(r = 15, l = 15)),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray95"),
    legend.position = "none"
  ) + 
  xlim(-12,5) + 
  ylim(1.5, 7.1)

dev.off()



################################################################################# 
# Summary by URBANICITY ONLY (across all EPA regions)
#############################################################################

group_summary_by_metro <- group_df_merge %>%
  group_by(metro4) %>%
  summarize(
    n_msa = n(),
    pw_no2 = mean(pw_no2_2024, na.rm = TRUE),
    var_pw_no2 = var(pw_no2_2024, na.rm = TRUE),
    pct_chg = mean(pct_change, na.rm = TRUE),
    var_pct_change = var(pct_change, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    se_pw_no2 = sqrt(var_pw_no2 / n_msa),
    se_pct_chg = sqrt(var_pct_change / n_msa),
    x0 = pct_chg,
    y0 = pw_no2,
    rx = 1.96 * se_pct_chg,
    ry = 1.96 * se_pw_no2
  ) %>%
  mutate(
    metro4 = factor(metro4, levels = c("Large MSA", "Small MSA", "MicroSA", "Non-urban"))
  )

nation_df <- group_df_merge %>%
  summarize(
    n_msa = n(),
    pw_no2 = mean(pw_no2_2024, na.rm = TRUE),
    var_pw_no2 = var(pw_no2_2024, na.rm = TRUE),
    pct_chg = mean(pct_change, na.rm = TRUE),
    var_pct_change = var(pct_change, na.rm = TRUE),
    .groups = "drop"
  )  %>%
  mutate(
    se_pw_no2 = sqrt(var_pw_no2 / n_msa),
    se_pct_chg = sqrt(var_pct_change / n_msa),
    x0 = pct_chg,
    y0 = pw_no2,
    rx = 1.96 * se_pct_chg,
    ry = 1.96 * se_pw_no2
  )

# Horizontal thin diamonds (x-axis 95% CI - replacing horizontal error bars)
diamond_metro_x <- bind_rows(lapply(seq_len(nrow(group_summary_by_metro)), function(i) {
  tibble(
    metro4 = group_summary_by_metro$metro4[i],
    diamond_id = paste0(group_summary_by_metro$metro4[i], "_x"),
    x = c(group_summary_by_metro$x0[i] - group_summary_by_metro$rx[i],
          group_summary_by_metro$x0[i],
          group_summary_by_metro$x0[i] + group_summary_by_metro$rx[i],
          group_summary_by_metro$x0[i]),
    y = c(group_summary_by_metro$y0[i],
          group_summary_by_metro$y0[i] + 0.025,
          group_summary_by_metro$y0[i],
          group_summary_by_metro$y0[i] - 0.025)
  )
}))

# Vertical thin diamonds (y-axis 95% CI - replacing vertical error bars)
diamond_metro_y <- bind_rows(lapply(seq_len(nrow(group_summary_by_metro)), function(i) {
  tibble(
    metro4 = group_summary_by_metro$metro4[i],
    diamond_id = paste0(group_summary_by_metro$metro4[i], "_y"),
    x = c(group_summary_by_metro$x0[i],
          group_summary_by_metro$x0[i] + 0.09,
          group_summary_by_metro$x0[i],
          group_summary_by_metro$x0[i] - 0.09),
    y = c(group_summary_by_metro$y0[i] - group_summary_by_metro$ry[i],
          group_summary_by_metro$y0[i],
          group_summary_by_metro$y0[i] + group_summary_by_metro$ry[i],
          group_summary_by_metro$y0[i])
  )
}))

diamond_metro <- bind_rows(diamond_metro_x, diamond_metro_y)

################################################################################# 
# Plot 2: URBANICITY STRATIFICATION
#############################################################################

pdf(file = "/Users/k_sy_n_imac/Library/CloudStorage/Dropbox/GWU/Lab/FINESST/Research/Results/Quadrant_NO2_Urbanicity.pdf",
    width = 9, height = 9)

ggplot() +
  # Reference lines
  geom_hline(yintercept = 5.31, linetype = "solid", color = "black", linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.5) +
  
  geom_vline(xintercept = -5, linetype = "dashed", color = "gray40", linewidth = 0.5) +
  geom_vline(xintercept = 5, linetype = "dashed", color = "gray40", linewidth = 0.5) +
  annotate("text",
           x = -9.8, y = 5.1,
           label = "WHO AQG = 5.3 ppb",
           size = 6.5,
           color = "black") +
  
  # Thin diamond uncertainty regions (all black)
  geom_polygon(data = diamond_metro, 
               aes(x = x, y = y, group = diamond_id),
               fill = "black", alpha = 0.7, linewidth = 1, color = NA) +
  
  # Center points (black, distinguished by shape only)
  geom_point(data = group_summary_by_metro,
             aes(x = x0, y = y0, shape = metro4),
             color = "black", size = 5, stroke = 1.2) +
  
  # National summary diamonds (uncertainty)
  geom_polygon(data = {
    h_diamond <- tibble(
      x = c(nation_df$pct_chg - nation_df$rx,
            nation_df$pct_chg,
            nation_df$pct_chg + nation_df$rx,
            nation_df$pct_chg),
      y = c(nation_df$pw_no2,
            nation_df$pw_no2 + 0.02,
            nation_df$pw_no2,
            nation_df$pw_no2 - 0.02),
      diamond_id = "national_x"
    )
    v_diamond <- tibble(
      x = c(nation_df$pct_chg,
            nation_df$pct_chg + 0.05,
            nation_df$pct_chg,
            nation_df$pct_chg - 0.05),
      y = c(nation_df$pw_no2 - nation_df$ry,
            nation_df$pw_no2,
            nation_df$pw_no2 + nation_df$ry,
            nation_df$pw_no2),
      diamond_id = "national_y"
    )
    bind_rows(h_diamond, v_diamond)
  },
  aes(x = x, y = y, group = diamond_id),
  alpha =  0.7, linewidth = 1, color = NA, fill = "#B11226") +
  
  # National summary point
  geom_point(data = nation_df,
             aes(x = pct_chg, y = pw_no2),
             color = "#B11226", size = 7, shape = 18) +
  geom_text(data = nation_df,
            aes(x = pct_chg, y = pw_no2),
            label = "National Summary",
            color = "#B11226", size = 7.5, fontface = "bold",
            vjust = -1.5, hjust = 0.3) +
  
  scale_shape_manual(values = c(
    "Large MSA" = 16,
    "Small MSA" = 17,
    "MicroSA" = 15,
    "Non-urban" = 8
  )) +
  
  labs(
    x = bquote(bold("% Change from 2019")),
    y = bquote(bold("NO"[bolditalic(2)] ~ "in 2024 (ppb)")),
    shape = "Urbanicity"
  ) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 22),
    axis.title.x = element_text(size = 25, margin = margin(t = 15, b = 15)),
    axis.title.y = element_text(size = 25, margin = margin(r = 15, l = 15)),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray95"),
    legend.position = "none"
  ) + 
  xlim(-12,5) + 
  ylim(1.5, 7.1)

dev.off()





###############################################################################################
# Two-way ANOVA ###############################################################################
###############################################################################################
# NO2 level in 2024
model_no2 <- lm(pw_no2_2024 ~ EPA_region + metro4, data = group_df_merge)
anova(model_no2)

# NO2 percent change
model_change <- lm(pct_change ~ EPA_region + metro4, data = group_df_merge)
anova(model_change)








###################################################################################################
# Figure 3: TREEMAP ###############################################################################
###################################################################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

# ── Named region labels ───────────────────────────────────────────────────────
region_labels <- c(
  "Region 1\n(New England)",
  "Region 2\n(Northeast/Caribbean)",
  "Region 3\n(Mid-Atlantic)",
  "Region 4\n(Southeast)",
  "Region 5\n(Great Lakes)",
  "Region 6\n(South Central)",
  "Region 7\n(Heartland)",
  "Region 8\n(Mountains/Plains)",
  "Region 9\n(Pacific Southwest)",
  "Region 10\n(Northwest)"
)

recode_region <- function(x) {
  recode(
    x,
    "Region 1"  = "Region 1\n(New England)",
    "Region 2"  = "Region 2\n(Northeast/Caribbean)",
    "Region 3"  = "Region 3\n(Mid-Atlantic)",
    "Region 4"  = "Region 4\n(Southeast)",
    "Region 5"  = "Region 5\n(Great Lakes)",
    "Region 6"  = "Region 6\n(South Central)",
    "Region 7"  = "Region 7\n(Heartland)",
    "Region 8"  = "Region 8\n(Mountains/Plains)",
    "Region 9"  = "Region 9\n(Pacific Southwest)",
    "Region 10" = "Region 10\n(Northwest)"
  )
}

# ─────────────────────────────────────────────────────────────────────────────
# Data prep
# ─────────────────────────────────────────────────────────────────────────────
list_year  <- c(2019, 2020, 2021, 2022, 2023, 2024)
list_cause <- c(
  "All causes",
  "Cardiovascular diseases",
  "Respiratory diseases",
  "Ischemic heart disease",
  "Tracheal, bronchus, and lung cancer"
)

i1 <- 1
i2 <- 6
j  <- 1

no2_ind <- subset(
  no2_tract_df,
  year %in% c(list_year[i1], list_year[i2]))
db_ind <- subset(
  final3,
  year %in% c(list_year[i1], list_year[i2]) &
    cause_name == list_cause[j])

final4 <- merge(no2_ind, db_ind, by = c("GEOID", "year"))
final4 <- merge(final4, xxxx, by = "GEOID", all.x = TRUE)
final4 <- final4 %>%
  mutate(
    MSA_title = if_else(
      metro4 == "Non-urban",
      paste0("Non_urban_", state_name.x),
      MSA_title),
    EPA_region = recode_region(EPA_region))

group_df <- final4 %>% group_by(MSA_title, metro4, EPA_region) %>%
  summarize(n_2024 = n(), .groups = "drop")
msa_regions <- group_df %>% count(MSA_title) %>% filter(n > 1)
majority_regions <- final4 %>%
  group_by(MSA_title, EPA_region) %>%
  summarise(total_pop = sum(pop_over20, na.rm = TRUE), .groups = "drop") %>%
  group_by(MSA_title) %>%
  slice_max(total_pop, n = 1) %>%  # ← Gets region with HIGHEST population
  ungroup() %>%
  dplyr::select(MSA_title, EPA_region_new = EPA_region)
final4 <- final4 %>%
  left_join(majority_regions, by = "MSA_title") %>%
  mutate(EPA_region = if_else(!is.na(EPA_region_new), 
                              EPA_region_new, 
                              EPA_region))

# ─────────────────────────────────────────────────────────────────────────────
# Calculations
# ─────────────────────────────────────────────────────────────────────────────
x <- subset(final4, year == list_year[i1])

group_df_2019 <- x %>%
  group_by(MSA_title, metro4, EPA_region) %>%
  summarize(
    pw_no2_2019 = weighted.mean(surface_no2, total_dc2020, na.rm = TRUE),
    n_2019 = n(),
    .groups = "drop"
  )

xx <- subset(final4, year == list_year[i2])

group_df_2024 <- xx %>%
  group_by(MSA_title, metro4, EPA_region) %>%
  summarize(
    pw_no2_2024 = weighted.mean(surface_no2, total_dc2020, na.rm = TRUE),
    n_2024 = n(),
    .groups = "drop"
  )

group_df <- merge(group_df_2019, group_df_2024, by = c("MSA_title", "EPA_region", "metro4")) %>%
  mutate(
    pct_change = (pw_no2_2024 - pw_no2_2019) / pw_no2_2019 * 100,
    quadrant = case_when(
      pct_change >  5 & pw_no2_2024 >= 5.31 ~ "High & Increasing",
      pct_change < -5 & pw_no2_2024 >= 5.31 ~ "High & Decreasing",
      pct_change >  5 & pw_no2_2024 <  5.31 ~ "Low & Increasing",
      pct_change < -5 & pw_no2_2024 <  5.31 ~ "Low & Decreasing",
      pw_no2_2024 >= 5.31                    ~ "High & Stable",
      pw_no2_2024 <  5.31                    ~ "Low & Stable"
    )
  )

# ─────────────────────────────────────────────────────────────────────────────
# Levels
# ─────────────────────────────────────────────────────────────────────────────
metro_levels   <- c("Large MSA", "Small MSA", "MicroSA", "Non-urban")
region_levels  <- region_labels
metro_levels2  <- c(metro_levels, "Total")
region_levels2 <- c(region_levels, "Total")

# ─────────────────────────────────────────────────────────────────────────────
# Aggregate cells
# ─────────────────────────────────────────────────────────────────────────────
obs_df_main <- group_df %>%
  group_by(metro4, EPA_region, quadrant) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(metro4, EPA_region) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

obs_df_row_total <- group_df %>%
  group_by(metro4, quadrant) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(metro4) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  mutate(EPA_region = "Total")

obs_df_col_total <- group_df %>%
  group_by(EPA_region, quadrant) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(EPA_region) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  mutate(metro4 = "Total")

obs_df_grand_total <- group_df %>%
  group_by(quadrant) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(
    prop = n / sum(n),
    metro4 = "Total",
    EPA_region = "Total"
  )

obs_df <- bind_rows(
  obs_df_main,
  obs_df_row_total,
  obs_df_col_total,
  obs_df_grand_total
)

# ─────────────────────────────────────────────────────────────────────────────
# Wide format
# ─────────────────────────────────────────────────────────────────────────────
plot_base <- obs_df %>%
  dplyr::select(metro4, EPA_region, quadrant, prop) %>%
  pivot_wider(
    names_from = quadrant,
    values_from = prop,
    values_fill = 0
  ) %>%
  transmute(
    metro4,
    EPA_region,
    Q1 = `High & Increasing`,
    Q2 = `High & Decreasing`,
    Q3 = `Low & Decreasing`,
    Q4 = `Low & Increasing`,
    Q5 = `High & Stable`,
    Q6 = `Low & Stable`
  ) %>%
  mutate(
    metro4     = factor(metro4, levels = metro_levels2),
    EPA_region = factor(EPA_region, levels = region_levels2),
    x = as.numeric(EPA_region),
    y = rev(seq_along(metro_levels2))[match(metro4, metro_levels2)]
  )

# ─────────────────────────────────────────────────────────────────────────────
# Nudge zeros
# ─────────────────────────────────────────────────────────────────────────────
epsilon <- 0.02

plot_base_nudged <- plot_base %>%
  mutate(
    Q1 = if_else(Q1 == 0, epsilon, Q1),
    Q2 = if_else(Q2 == 0, epsilon, Q2),
    Q3 = if_else(Q3 == 0, epsilon, Q3),
    Q4 = if_else(Q4 == 0, epsilon, Q4),
    Q5 = if_else(Q5 == 0, epsilon, Q5),
    Q6 = if_else(Q6 == 0, epsilon, Q6)
  ) %>%
  mutate(
    total = Q1 + Q2 + Q3 + Q4 + Q5 + Q6,
    Q1 = Q1 / total,
    Q2 = Q2 / total,
    Q3 = Q3 / total,
    Q4 = Q4 / total,
    Q5 = Q5 / total,
    Q6 = Q6 / total
  ) %>%
  dplyr::select(-total)

# ─────────────────────────────────────────────────────────────────────────────
# Build rect_df
# ─────────────────────────────────────────────────────────────────────────────
pad   <- 0.09
min_w <- 0.03
min_h <- 0.03

rect_df <- plot_base %>%
  rowwise() %>%
  do({
    row <- .
    
    x_left  <- row$x - 0.5 + pad
    x_right <- row$x + 0.5 - pad
    y_bot   <- row$y - 0.5 + pad
    y_top   <- row$y + 0.5 - pad
    
    total_w <- x_right - x_left
    total_h <- y_top - y_bot
    
    alloc <- function(vals, total_px, min_px, tiny_prop = 0.1) {
      
      vals <- replace_na(vals, 0)
      
      # Cells that are zero OR very close to zero
      # will receive exactly the same size as NA/zero cells.
      tiny_cells <- vals <= tiny_prop
      
      n_tiny <- sum(tiny_cells)
      reserved <- n_tiny * min_px
      remaining <- max(total_px - reserved, 0)
      
      non_tiny_sum <- sum(vals[!tiny_cells])
      
      out <- numeric(length(vals))
      
      # Tiny or zero cells get the same fixed size
      out[tiny_cells] <- min_px
      
      # Other cells divide the remaining space proportionally
      if (non_tiny_sum > 0) {
        out[!tiny_cells] <- remaining * vals[!tiny_cells] / non_tiny_sum
      } else {
        out[] <- total_px / length(vals)
      }
      
      out
    }
    
    left_raw  <- row$Q2 + row$Q3
    mid_raw   <- row$Q5 + row$Q6
    right_raw <- row$Q1 + row$Q4
    
    x_widths <- alloc(
      c(left_raw, mid_raw, right_raw),
      total_w,
      min_w
    )
    
    x_mid1 <- x_left + x_widths[1]
    x_mid2 <- x_mid1 + x_widths[2]
    
    q2 <- row$Q2
    q3 <- row$Q3
    
    y_heights_l <- alloc(
      c(q3, q2),
      total_h,
      min_h
    )
    
    y_left_split <- y_bot + y_heights_l[1]
    
    q5 <- row$Q5
    q6 <- row$Q6
    
    y_heights_m <- alloc(
      c(q6, q5),
      total_h,
      min_h
    )
    
    y_mid_split <- y_bot + y_heights_m[1]
    
    q1 <- row$Q1
    q4 <- row$Q4
    
    y_heights_r <- alloc(
      c(q4, q1),
      total_h,
      min_h
    )
    
    y_right_split <- y_bot + y_heights_r[1]
    
    tibble(
      metro4     = row$metro4,
      EPA_region = row$EPA_region,
      quadrant   = c("Q1", "Q2", "Q3", "Q4", "Q5", "Q6"),
      prob       = c(row$Q1, row$Q2, row$Q3, row$Q4, row$Q5, row$Q6),
      xmin       = c(x_mid2, x_left, x_left, x_mid2, x_mid1, x_mid1),
      xmax       = c(x_right, x_mid1, x_mid1, x_right, x_mid2, x_mid2),
      ymin       = c(y_right_split, y_left_split, y_bot, y_bot, y_mid_split, y_bot),
      ymax       = c(y_top, y_top, y_left_split, y_right_split, y_top, y_mid_split)
    )
  }) %>%
  ungroup() %>%
  mutate(
    is_stable = quadrant %in% c("Q5", "Q6")
  )

rect_df <- rect_df %>%
  left_join(
    obs_df %>%
      mutate(
        quadrant = recode(
          quadrant,
          "High & Increasing" = "Q1",
          "High & Decreasing" = "Q2",
          "Low & Decreasing"  = "Q3",
          "Low & Increasing"  = "Q4",
          "High & Stable"     = "Q5",
          "Low & Stable"      = "Q6"
        ),
        metro4     = factor(metro4, levels = metro_levels2),
        EPA_region = factor(EPA_region, levels = region_levels2)
      ) %>%
      dplyr::select(metro4, EPA_region, quadrant, n),
    by = c("metro4", "EPA_region", "quadrant")
  ) %>%
  mutate(
    n = replace_na(n, 0),
    fill_col = if_else(n == 0, "black", quadrant)
  )

# ─────────────────────────────────────────────────────────────────────────────
# label_df
# ─────────────────────────────────────────────────────────────────────────────

narrow_width_threshold <- 1 / 6
bottom_text_pad <- 0.06

label_df <- obs_df %>%
  mutate(
    metro4     = factor(metro4, levels = metro_levels2),
    EPA_region = factor(EPA_region, levels = region_levels2),
    quadrant = recode(
      quadrant,
      "High & Increasing" = "Q1",
      "High & Decreasing" = "Q2",
      "Low & Decreasing"  = "Q3",
      "Low & Increasing"  = "Q4",
      "High & Stable"     = "Q5",
      "Low & Stable"      = "Q6"
    )
  ) %>%
  left_join(
    rect_df %>%
      mutate(
        x = (xmin + xmax) / 2,
        y = (ymin + ymax) / 2,
        area = (xmax - xmin) * (ymax - ymin),
        cell_width = xmax - xmin,
        cell_height = ymax - ymin
      ) %>%
      dplyr::select(
        metro4,
        EPA_region,
        quadrant,
        x,
        y,
        area,
        xmin,
        xmax,
        ymin,
        ymax,
        cell_width,
        cell_height
      ),
    by = c("metro4", "EPA_region", "quadrant")
  ) %>%
  mutate(
    label = case_when(
      n == 0        ~ NA_character_,
      area >= 0.002 ~ paste0(n, "\n(", round(prop * 100), "%)"),
      TRUE          ~ NA_character_
    ),
    
    txt_size = case_when(
      area >= 0.008 ~ 5.2,
      area >= 0.002 ~ 5.2,
      TRUE          ~ 5.2
    ),
    
    narrow_cell = cell_width < narrow_width_threshold,
    short_cell = cell_height < (min_h + 0.1),
    # Only align at bottom if narrow AND not short
    push_to_bottom = narrow_cell & !short_cell,
    y = if_else(
      push_to_bottom,
      ymin + bottom_text_pad,
      y
    ),
    text_vjust = if_else(
      push_to_bottom,
      0,
      0.5
    )
  ) %>%
  filter(!is.na(label))
# ─────────────────────────────────────────────────────────────────────────────
# Manual label adjustment function
# ─────────────────────────────────────────────────────────────────────────────
label_df <- label_df %>%
  mutate(
    x_original = x,
    y_original = y,
    label_length = nchar(gsub("\n", "", label))
  )

two_line_shift <- 0.16

move_shorter_label <- function(
    df,
    cell1,
    cell2,
    overlap_type = c("horizontal", "vertical"),
    shift = two_line_shift
) {
  
  overlap_type <- match.arg(overlap_type)
  
  idx1 <- which(
    as.character(df$metro4) == cell1$metro4 &
      as.character(df$EPA_region) == cell1$EPA_region &
      as.character(df$quadrant) == cell1$quadrant
  )
  
  idx2 <- which(
    as.character(df$metro4) == cell2$metro4 &
      as.character(df$EPA_region) == cell2$EPA_region &
      as.character(df$quadrant) == cell2$quadrant
  )
  
  if (length(idx1) == 0 | length(idx2) == 0) {
    return(df)
  }
  
  shorter_idx <- ifelse(
    df$label_length[idx1] <= df$label_length[idx2],
    idx1,
    idx2
  )
  
  if (overlap_type == "horizontal") {
    
    # If two labels overlap horizontally,
    # move the shorter one downward.
    df$y[shorter_idx] <- df$y[shorter_idx] - shift
    
    # Keep the label inside its rectangle.
    df$y[shorter_idx] <- pmax(
      df$ymin[shorter_idx] + 0.06,
      df$y[shorter_idx]
    )
  }
  
  if (overlap_type == "vertical") {
    
    # If two labels overlap vertically,
    # move the shorter one toward the center of its own rectangle.
    cell_center_x <- (df$xmin[shorter_idx] + df$xmax[shorter_idx]) / 2
    cell_center_y <- (df$ymin[shorter_idx] + df$ymax[shorter_idx]) / 2
    
    df$x[shorter_idx] <- cell_center_x
    df$y[shorter_idx] <- cell_center_y
    
    # Optional slight inward/downward adjustment.
    df$y[shorter_idx] <- df$y[shorter_idx] - shift
    
    # Keep label inside its rectangle.
    df$x[shorter_idx] <- pmin(
      df$xmax[shorter_idx] - 0.04,
      pmax(df$xmin[shorter_idx] + 0.04, df$x[shorter_idx])
    )
    
    df$y[shorter_idx] <- pmin(
      df$ymax[shorter_idx] - 0.06,
      pmax(df$ymin[shorter_idx] + 0.06, df$y[shorter_idx])
    )
  }
  
  df
}

# ─────────────────────────────────────────────────────────────────────────────
# Manual label adjustments
# Edit this section based on the overlapping labels you see.
# ─────────────────────────────────────────────────────────────────────────────

# Example 1:
# If Q1 and Q5 overlap horizontally within the same cell,
# the shorter label will be moved downward.
label_df <- move_shorter_label(
  label_df,
  cell1 = list(
    metro4 = "Large MSA",
    EPA_region = "Region 6\n(South Central)",
    quadrant = "Q1"
  ),
  cell2 = list(
    metro4 = "Large MSA",
    EPA_region = "Region 6\n(South Central)",
    quadrant = "Q5"
  ),
  overlap_type = "horizontal"
)

# Example 2:
# If two labels overlap vertically, the shorter label will be moved
# toward the center of its own rectangle.
label_df <- move_shorter_label(
  label_df,
  cell1 = list(
    metro4 = "Small MSA",
    EPA_region = "Region 7\n(Heartland)",
    quadrant = "Q4"
  ),
  cell2 = list(
    metro4 = "MicroSA",
    EPA_region = "Region 7\n(Heartland)",
    quadrant = "Q4"
  ),
  overlap_type = "vertical"
)

# Add more manual adjustment blocks as needed:
#
# label_df <- move_shorter_label(
#   label_df,
#   cell1 = list(
#     metro4 = "YOUR_METRO4",
#     EPA_region = "YOUR_REGION",
#     quadrant = "Q?"
#   ),
#   cell2 = list(
#     metro4 = "YOUR_METRO4",
#     EPA_region = "YOUR_REGION",
#     quadrant = "Q?"
#   ),
#   overlap_type = "horizontal"
# )

# ─────────────────────────────────────────────────────────────────────────────
# Colors
# ─────────────────────────────────────────────────────────────────────────────
#quad_colors <- c(
#  "Q1" = "#d73027",
#  "Q2" = "#fee08b",
#  "Q3" = "#6baed6",
#  "Q4" = "#9E8FCC",
#  "Q5" = "#fdb863",
#  "Q6" = "#2196a8"
#)

quad_colors <- c(
  "Q1" = "#c0152a",  # dark red   (top-right)
  "Q2" = "#fcc5c0",  # light red  (top-left)
  "Q3" = "#c6dbef",  # light blue (bottom-left)
  "Q4" = "#1a6bb5",  # dark blue  (bottom-right)
  "Q5" = "#f1717e",  # mid red    (top-middle)
  "Q6" = "#6baed6"   # mid blue   (bottom-middle)
)


# ─────────────────────────────────────────────────────────────────────────────
# Legend
# ─────────────────────────────────────────────────────────────────────────────
legend_plot <- ggplot() +
  annotate(
    "rect",
    xmin = 0.007, xmax = 1,
    ymin = 0.007, ymax = 1,
    fill = quad_colors["Q1"]
  ) +
  annotate(
    "rect",
    xmin = -1, xmax = -0.007,
    ymin = 0.007, ymax = 1,
    fill = quad_colors["Q2"]
  ) +
  annotate(
    "rect",
    xmin = -1, xmax = -0.007,
    ymin = -1, ymax = -0.007,
    fill = quad_colors["Q3"]
  ) +
  annotate(
    "rect",
    xmin = 0.007, xmax = 1,
    ymin = -1, ymax = -0.007,
    fill = quad_colors["Q4"]
  ) +
  annotate(
    "rect",
    xmin = -0.2, xmax = 0.2,
    ymin = 0, ymax = 1,
    fill = quad_colors["Q5"]
  ) +
  annotate(
    "rect",
    xmin = -0.2, xmax = 0.2,
    ymin = -1, ymax = 0,
    fill = quad_colors["Q6"]
  ) +
  annotate(
    "rect",
    xmin = -0.009, xmax = 0.009,
    ymin = -1, ymax = 1,
    fill = "grey70"
  ) +
  annotate(
    "rect",
    xmin = -1, xmax = 1,
    ymin = -0.009, ymax = 0.009,
    fill = "grey70"
  ) +
  annotate(
    "text",
    x = 0.6, y = 0.3,
    label = "S1\nHigh &\nIncreasing",
    size = 4,
    fontface = "bold",
    color = "white",
    lineheight = 0.95
  ) +
  annotate(
    "text",
    x = -0.6, y = 0.3,
    label = "S2\nHigh &\nDecreasing",
    size = 4,
    fontface = "bold",
    color = "grey20",
    lineheight = 0.95
  ) +
  annotate(
    "text",
    x = -0.6, y = -0.3,
    label = "S3\nLow &\nDecreasing",
    size = 4,
    fontface = "bold",
    color = "grey20",
    lineheight = 0.95
  ) +
  annotate(
    "text",
    x = 0.6, y = -0.3,
    label = "S4\nLow &\nIncreasing",
    size = 4,
    fontface = "bold",
    color = "white",
    lineheight = 0.95
  ) +
  annotate(
    "text",
    x = 0.00, y = 0.7,
    label = "S5\nHigh &\nStable",
    size = 3.8,
    fontface = "bold",
    color = "grey20",
    lineheight = 0.95
  ) +
  annotate(
    "text",
    x = 0.00, y = -0.7,
    label = "S6\nLow &\nStable",
    size = 3.8,
    fontface = "bold",
    color = "grey20",
    lineheight = 0.95
  ) +
  annotate(
    "text",
    x = 0, y = -1.12,
    label = "% Change from 2019",
    size = 4
  ) +
  annotate(
    "text",
    x = -1.2, y = -0.05,
    label = "NO\u2082 in 2024",
    size = 4,
    angle = 90
  ) +
  theme_void()

# ─────────────────────────────────────────────────────────────────────────────
# Treemap
# ─────────────────────────────────────────────────────────────────────────────
x_faces <- ifelse(region_levels2 == "Total", "bold", "plain")
y_faces <- ifelse(metro_levels2  == "Total", "bold", "plain")


treemap_plot <- ggplot(rect_df) +
  geom_rect(
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax,
      fill = fill_col
    ),
    color = "white",
    linewidth = 0.5
  ) +
  scale_fill_manual(
    values = c(quad_colors, "black" = "grey40")
  ) +
  guides(fill = "none") +
  geom_text(
    data = label_df,
    aes(
      x = x,
      y = y,
      label = label,
      size = txt_size,
      vjust = text_vjust
    ),
    lineheight = 0.9
  ) +
  scale_size_identity() +
  scale_x_continuous(
    breaks = seq_along(region_levels2),
    labels = gsub(" \\(", "\n(", region_levels2),
    expand = c(0.02, 0.02),
    position = "top"
  ) +
  scale_y_continuous(
    breaks = rev(seq_along(metro_levels2)),
    labels = metro_levels2,
    expand = c(0, 0)
  ) +
  coord_fixed(clip = "off") +
  labs(x = "", y = "") +
  theme_minimal(base_family = "Lato", base_size = 13) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid       = element_blank(),
    panel.border     = element_rect(color = "white", fill = NA, linewidth = 0.8),
    axis.text.x.top      = element_text(
      angle = 30, color="black",
      hjust = 0.5,
      vjust = 0.5,
      face = x_faces,
      size = 20, margin = margin(b = 5)
    ),
    axis.text.y      = element_text(
      size = 20, 
      color="black", 
      face = y_faces,
      angle = 30,        # ← Tilt y-axis text
      hjust = 0.5,
      vjust = 1,
      margin = margin(r = 0)
    ),
    axis.ticks       = element_blank(),
    plot.margin      = margin(t = 15, r = 15, b = 15, l = 5) 
  )

# ─────────────────────────────────────────────────────────────────────────────
# Save
# ─────────────────────────────────────────────────────────────────────────────
pdf(file = "/Users/k_sy_n_imac/Library/CloudStorage/Dropbox/GWU/Lab/FINESST/Research/Results/Quadrant_NO2_treemap_with_totals.pdf",
    width = 26, height = 12)
print(treemap_plot)
dev.off()

pdf(file = "/Users/k_sy_n_imac/Library/CloudStorage/Dropbox/GWU/Lab/FINESST/Research/Results/Quadrant_NO2_treemap_with_totals(legend).pdf",
    width = 3.5, height = 3.5)
print(legend_plot)
dev.off()










# Supplementary texts & : sensitivity analysis

###################################################################################################
# Figure ST1: TREEMAP ###############################################################################
###################################################################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

# ── Named region labels ───────────────────────────────────────────────────────
region_labels <- c(
  "Region 1\n(New England)",
  "Region 2\n(Northeast/Caribbean)",
  "Region 3\n(Mid-Atlantic)",
  "Region 4\n(Southeast)",
  "Region 5\n(Great Lakes)",
  "Region 6\n(South Central)",
  "Region 7\n(Heartland)",
  "Region 8\n(Mountains/Plains)",
  "Region 9\n(Pacific Southwest)",
  "Region 10\n(Northwest)"
)

recode_region <- function(x) {
  recode(
    x,
    "Region 1"  = "Region 1\n(New England)",
    "Region 2"  = "Region 2\n(Northeast/Caribbean)",
    "Region 3"  = "Region 3\n(Mid-Atlantic)",
    "Region 4"  = "Region 4\n(Southeast)",
    "Region 5"  = "Region 5\n(Great Lakes)",
    "Region 6"  = "Region 6\n(South Central)",
    "Region 7"  = "Region 7\n(Heartland)",
    "Region 8"  = "Region 8\n(Mountains/Plains)",
    "Region 9"  = "Region 9\n(Pacific Southwest)",
    "Region 10" = "Region 10\n(Northwest)"
  )
}

# ─────────────────────────────────────────────────────────────────────────────
# Data prep
# ─────────────────────────────────────────────────────────────────────────────
list_year  <- c(2019, 2020, 2021, 2022, 2023, 2024)
list_cause <- c(
  "All causes",
  "Cardiovascular diseases",
  "Respiratory diseases",
  "Ischemic heart disease",
  "Tracheal, bronchus, and lung cancer"
)

i1 <- 1
i2 <- 6
j  <- 1

no2_ind <- subset(
  no2_tract_df,
  year %in% c(list_year[i1], list_year[i2]))
db_ind <- subset(
  final3,
  year %in% c(list_year[i1], list_year[i2]) &
    cause_name == list_cause[j])

final4 <- merge(no2_ind, db_ind, by = c("GEOID", "year"))
final4 <- merge(final4, xxxx, by = "GEOID", all.x = TRUE)
final4 <- final4 %>%
  mutate(
    MSA_title = if_else(
      metro4 == "Non-urban",
      paste0("Non_urban_", state_name.x),
      MSA_title),
    EPA_region = recode_region(EPA_region))

group_df <- final4 %>% group_by(MSA_title, metro4, EPA_region) %>%
  summarize(n_2024 = n(), .groups = "drop")
msa_regions <- group_df %>% count(MSA_title) %>% filter(n > 1)
majority_regions <- final4 %>%
  group_by(MSA_title, EPA_region) %>%
  summarise(total_pop = sum(pop_over20, na.rm = TRUE), .groups = "drop") %>%
  group_by(MSA_title) %>%
  slice_max(total_pop, n = 1) %>%  # ← Gets region with HIGHEST population
  ungroup() %>%
  dplyr::select(MSA_title, EPA_region_new = EPA_region)
final4 <- final4 %>%
  left_join(majority_regions, by = "MSA_title") %>%
  mutate(EPA_region = if_else(!is.na(EPA_region_new), 
                              EPA_region_new, 
                              EPA_region))

# ─────────────────────────────────────────────────────────────────────────────
# Calculations
# ─────────────────────────────────────────────────────────────────────────────
x <- subset(final4, year == list_year[i1])

group_df_2019 <- x %>%
  group_by(MSA_title, metro4, EPA_region) %>%
  summarize(
    pw_no2_2019 = weighted.mean(surface_no2, total_dc2020, na.rm = TRUE),
    n_2019 = n(),
    .groups = "drop"
  )

xx <- subset(final4, year == list_year[i2])

group_df_2024 <- xx %>%
  group_by(MSA_title, metro4, EPA_region) %>%
  summarize(
    pw_no2_2024 = weighted.mean(surface_no2, total_dc2020, na.rm = TRUE),
    n_2024 = n(),
    .groups = "drop"
  )

group_df <- merge(group_df_2019, group_df_2024, by = c("MSA_title", "EPA_region", "metro4")) %>%
  mutate(
    pct_change = (pw_no2_2024 - pw_no2_2019) / pw_no2_2019 * 100,
    quadrant = case_when(
      pct_change >  10 & pw_no2_2024 >= 5.31 ~ "High & Increasing",
      pct_change < -10 & pw_no2_2024 >= 5.31 ~ "High & Decreasing",
      pct_change >  10 & pw_no2_2024 <  5.31 ~ "Low & Increasing",
      pct_change < -10 & pw_no2_2024 <  5.31 ~ "Low & Decreasing",
      pw_no2_2024 >= 5.31                    ~ "High & Stable",
      pw_no2_2024 <  5.31                    ~ "Low & Stable"
    )
  )

# ─────────────────────────────────────────────────────────────────────────────
# Levels
# ─────────────────────────────────────────────────────────────────────────────
metro_levels   <- c("Large MSA", "Small MSA", "MicroSA", "Non-urban")
region_levels  <- region_labels
metro_levels2  <- c(metro_levels, "Total")
region_levels2 <- c(region_levels, "Total")

# ─────────────────────────────────────────────────────────────────────────────
# Aggregate cells
# ─────────────────────────────────────────────────────────────────────────────
obs_df_main <- group_df %>%
  group_by(metro4, EPA_region, quadrant) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(metro4, EPA_region) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

obs_df_row_total <- group_df %>%
  group_by(metro4, quadrant) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(metro4) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  mutate(EPA_region = "Total")

obs_df_col_total <- group_df %>%
  group_by(EPA_region, quadrant) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(EPA_region) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  mutate(metro4 = "Total")

obs_df_grand_total <- group_df %>%
  group_by(quadrant) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(
    prop = n / sum(n),
    metro4 = "Total",
    EPA_region = "Total"
  )

obs_df <- bind_rows(
  obs_df_main,
  obs_df_row_total,
  obs_df_col_total,
  obs_df_grand_total
)

# ─────────────────────────────────────────────────────────────────────────────
# Wide format
# ─────────────────────────────────────────────────────────────────────────────
plot_base <- obs_df %>%
  dplyr::select(metro4, EPA_region, quadrant, prop) %>%
  pivot_wider(
    names_from  = quadrant,
    values_from = prop,
    values_fill = 0
  ) %>%
  # Add any missing quadrant columns as 0 (happens when no MSA falls in that category)
  mutate(
    `High & Increasing` = if ("High & Increasing" %in% names(.)) `High & Increasing` else 0,
    `High & Decreasing` = if ("High & Decreasing" %in% names(.)) `High & Decreasing` else 0,
    `Low & Decreasing`  = if ("Low & Decreasing"  %in% names(.)) `Low & Decreasing`  else 0,
    `Low & Increasing`  = if ("Low & Increasing"  %in% names(.)) `Low & Increasing`  else 0,
    `High & Stable`     = if ("High & Stable"     %in% names(.)) `High & Stable`     else 0,
    `Low & Stable`      = if ("Low & Stable"      %in% names(.)) `Low & Stable`      else 0
  ) %>%
  transmute(
    metro4,
    EPA_region,
    Q1 = `High & Increasing`,
    Q2 = `High & Decreasing`,
    Q3 = `Low & Decreasing`,
    Q4 = `Low & Increasing`,
    Q5 = `High & Stable`,
    Q6 = `Low & Stable`
  ) %>%
  mutate(
    metro4     = factor(metro4, levels = metro_levels2),
    EPA_region = factor(EPA_region, levels = region_levels2),
    x = as.numeric(EPA_region),
    y = rev(seq_along(metro_levels2))[match(metro4, metro_levels2)]
  )

# ─────────────────────────────────────────────────────────────────────────────
# Nudge zeros
# ─────────────────────────────────────────────────────────────────────────────
epsilon <- 0.02

plot_base_nudged <- plot_base %>%
  mutate(
    Q1 = if_else(Q1 == 0, epsilon, Q1),
    Q2 = if_else(Q2 == 0, epsilon, Q2),
    Q3 = if_else(Q3 == 0, epsilon, Q3),
    Q4 = if_else(Q4 == 0, epsilon, Q4),
    Q5 = if_else(Q5 == 0, epsilon, Q5),
    Q6 = if_else(Q6 == 0, epsilon, Q6)
  ) %>%
  mutate(
    total = Q1 + Q2 + Q3 + Q4 + Q5 + Q6,
    Q1 = Q1 / total,
    Q2 = Q2 / total,
    Q3 = Q3 / total,
    Q4 = Q4 / total,
    Q5 = Q5 / total,
    Q6 = Q6 / total
  ) %>%
  dplyr::select(-total)

# ─────────────────────────────────────────────────────────────────────────────
# Build rect_df
# ─────────────────────────────────────────────────────────────────────────────
pad   <- 0.09
min_w <- 0.03
min_h <- 0.03

rect_df <- plot_base %>%
  rowwise() %>%
  do({
    row <- .
    
    x_left  <- row$x - 0.5 + pad
    x_right <- row$x + 0.5 - pad
    y_bot   <- row$y - 0.5 + pad
    y_top   <- row$y + 0.5 - pad
    
    total_w <- x_right - x_left
    total_h <- y_top - y_bot
    
    alloc <- function(vals, total_px, min_px, tiny_prop = 0.1) {
      
      vals <- replace_na(vals, 0)
      
      # Cells that are zero OR very close to zero
      # will receive exactly the same size as NA/zero cells.
      tiny_cells <- vals <= tiny_prop
      
      n_tiny <- sum(tiny_cells)
      reserved <- n_tiny * min_px
      remaining <- max(total_px - reserved, 0)
      
      non_tiny_sum <- sum(vals[!tiny_cells])
      
      out <- numeric(length(vals))
      
      # Tiny or zero cells get the same fixed size
      out[tiny_cells] <- min_px
      
      # Other cells divide the remaining space proportionally
      if (non_tiny_sum > 0) {
        out[!tiny_cells] <- remaining * vals[!tiny_cells] / non_tiny_sum
      } else {
        out[] <- total_px / length(vals)
      }
      
      out
    }
    
    left_raw  <- row$Q2 + row$Q3
    mid_raw   <- row$Q5 + row$Q6
    right_raw <- row$Q1 + row$Q4
    
    x_widths <- alloc(
      c(left_raw, mid_raw, right_raw),
      total_w,
      min_w
    )
    
    x_mid1 <- x_left + x_widths[1]
    x_mid2 <- x_mid1 + x_widths[2]
    
    q2 <- row$Q2
    q3 <- row$Q3
    
    y_heights_l <- alloc(
      c(q3, q2),
      total_h,
      min_h
    )
    
    y_left_split <- y_bot + y_heights_l[1]
    
    q5 <- row$Q5
    q6 <- row$Q6
    
    y_heights_m <- alloc(
      c(q6, q5),
      total_h,
      min_h
    )
    
    y_mid_split <- y_bot + y_heights_m[1]
    
    q1 <- row$Q1
    q4 <- row$Q4
    
    y_heights_r <- alloc(
      c(q4, q1),
      total_h,
      min_h
    )
    
    y_right_split <- y_bot + y_heights_r[1]
    
    tibble(
      metro4     = row$metro4,
      EPA_region = row$EPA_region,
      quadrant   = c("Q1", "Q2", "Q3", "Q4", "Q5", "Q6"),
      prob       = c(row$Q1, row$Q2, row$Q3, row$Q4, row$Q5, row$Q6),
      xmin       = c(x_mid2, x_left, x_left, x_mid2, x_mid1, x_mid1),
      xmax       = c(x_right, x_mid1, x_mid1, x_right, x_mid2, x_mid2),
      ymin       = c(y_right_split, y_left_split, y_bot, y_bot, y_mid_split, y_bot),
      ymax       = c(y_top, y_top, y_left_split, y_right_split, y_top, y_mid_split)
    )
  }) %>%
  ungroup() %>%
  mutate(
    is_stable = quadrant %in% c("Q5", "Q6")
  )

rect_df <- rect_df %>%
  left_join(
    obs_df %>%
      mutate(
        quadrant = recode(
          quadrant,
          "High & Increasing" = "Q1",
          "High & Decreasing" = "Q2",
          "Low & Decreasing"  = "Q3",
          "Low & Increasing"  = "Q4",
          "High & Stable"     = "Q5",
          "Low & Stable"      = "Q6"
        ),
        metro4     = factor(metro4, levels = metro_levels2),
        EPA_region = factor(EPA_region, levels = region_levels2)
      ) %>%
      dplyr::select(metro4, EPA_region, quadrant, n),
    by = c("metro4", "EPA_region", "quadrant")
  ) %>%
  mutate(
    n = replace_na(n, 0),
    fill_col = if_else(n == 0, "black", quadrant)
  )

# ─────────────────────────────────────────────────────────────────────────────
# label_df
# ─────────────────────────────────────────────────────────────────────────────

narrow_width_threshold <- 1 / 6
bottom_text_pad <- 0.06

label_df <- obs_df %>%
  mutate(
    metro4     = factor(metro4, levels = metro_levels2),
    EPA_region = factor(EPA_region, levels = region_levels2),
    quadrant = recode(
      quadrant,
      "High & Increasing" = "Q1",
      "High & Decreasing" = "Q2",
      "Low & Decreasing"  = "Q3",
      "Low & Increasing"  = "Q4",
      "High & Stable"     = "Q5",
      "Low & Stable"      = "Q6"
    )
  ) %>%
  left_join(
    rect_df %>%
      mutate(
        x = (xmin + xmax) / 2,
        y = (ymin + ymax) / 2,
        area = (xmax - xmin) * (ymax - ymin),
        cell_width = xmax - xmin,
        cell_height = ymax - ymin
      ) %>%
      dplyr::select(
        metro4,
        EPA_region,
        quadrant,
        x,
        y,
        area,
        xmin,
        xmax,
        ymin,
        ymax,
        cell_width,
        cell_height
      ),
    by = c("metro4", "EPA_region", "quadrant")
  ) %>%
  mutate(
    label = case_when(
      n == 0        ~ NA_character_,
      area >= 0.002 ~ paste0(n, "\n(", round(prop * 100), "%)"),
      TRUE          ~ NA_character_
    ),
    
    txt_size = case_when(
      area >= 0.008 ~ 5.2,
      area >= 0.002 ~ 5.2,
      TRUE          ~ 5.2
    ),
    
    narrow_cell = cell_width < narrow_width_threshold,
    short_cell = cell_height < (min_h + 0.1),
    # Only align at bottom if narrow AND not short
    push_to_bottom = narrow_cell & !short_cell,
    y = if_else(
      push_to_bottom,
      ymin + bottom_text_pad,
      y
    ),
    text_vjust = if_else(
      push_to_bottom,
      0,
      0.5
    )
  ) %>%
  filter(!is.na(label))
# ─────────────────────────────────────────────────────────────────────────────
# Manual label adjustment function
# ─────────────────────────────────────────────────────────────────────────────
label_df <- label_df %>%
  mutate(
    x_original = x,
    y_original = y,
    label_length = nchar(gsub("\n", "", label))
  )

two_line_shift <- 0.16

move_shorter_label <- function(
    df,
    cell1,
    cell2,
    overlap_type = c("horizontal", "vertical"),
    shift = two_line_shift
) {
  
  overlap_type <- match.arg(overlap_type)
  
  idx1 <- which(
    as.character(df$metro4) == cell1$metro4 &
      as.character(df$EPA_region) == cell1$EPA_region &
      as.character(df$quadrant) == cell1$quadrant
  )
  
  idx2 <- which(
    as.character(df$metro4) == cell2$metro4 &
      as.character(df$EPA_region) == cell2$EPA_region &
      as.character(df$quadrant) == cell2$quadrant
  )
  
  if (length(idx1) == 0 | length(idx2) == 0) {
    return(df)
  }
  
  shorter_idx <- ifelse(
    df$label_length[idx1] <= df$label_length[idx2],
    idx1,
    idx2
  )
  
  if (overlap_type == "horizontal") {
    
    # If two labels overlap horizontally,
    # move the shorter one downward.
    df$y[shorter_idx] <- df$y[shorter_idx] - shift
    
    # Keep the label inside its rectangle.
    df$y[shorter_idx] <- pmax(
      df$ymin[shorter_idx] + 0.06,
      df$y[shorter_idx]
    )
  }
  
  if (overlap_type == "vertical") {
    
    # If two labels overlap vertically,
    # move the shorter one toward the center of its own rectangle.
    cell_center_x <- (df$xmin[shorter_idx] + df$xmax[shorter_idx]) / 2
    cell_center_y <- (df$ymin[shorter_idx] + df$ymax[shorter_idx]) / 2
    
    df$x[shorter_idx] <- cell_center_x
    df$y[shorter_idx] <- cell_center_y
    
    # Optional slight inward/downward adjustment.
    df$y[shorter_idx] <- df$y[shorter_idx] - shift
    
    # Keep label inside its rectangle.
    df$x[shorter_idx] <- pmin(
      df$xmax[shorter_idx] - 0.04,
      pmax(df$xmin[shorter_idx] + 0.04, df$x[shorter_idx])
    )
    
    df$y[shorter_idx] <- pmin(
      df$ymax[shorter_idx] - 0.06,
      pmax(df$ymin[shorter_idx] + 0.06, df$y[shorter_idx])
    )
  }
  
  df
}

# ─────────────────────────────────────────────────────────────────────────────
# Manual label adjustments
# Edit this section based on the overlapping labels you see.
# ─────────────────────────────────────────────────────────────────────────────

# Example 1:
# If Q1 and Q5 overlap horizontally within the same cell,
# the shorter label will be moved downward.
label_df <- move_shorter_label(
  label_df,
  cell1 = list(
    metro4 = "Large MSA",
    EPA_region = "Region 6\n(South Central)",
    quadrant = "Q1"
  ),
  cell2 = list(
    metro4 = "Large MSA",
    EPA_region = "Region 6\n(South Central)",
    quadrant = "Q5"
  ),
  overlap_type = "horizontal"
)

# Example 2:
# If two labels overlap vertically, the shorter label will be moved
# toward the center of its own rectangle.
label_df <- move_shorter_label(
  label_df,
  cell1 = list(
    metro4 = "Small MSA",
    EPA_region = "Region 7\n(Heartland)",
    quadrant = "Q4"
  ),
  cell2 = list(
    metro4 = "MicroSA",
    EPA_region = "Region 7\n(Heartland)",
    quadrant = "Q4"
  ),
  overlap_type = "vertical"
)

# Add more manual adjustment blocks as needed:
#
# label_df <- move_shorter_label(
#   label_df,
#   cell1 = list(
#     metro4 = "YOUR_METRO4",
#     EPA_region = "YOUR_REGION",
#     quadrant = "Q?"
#   ),
#   cell2 = list(
#     metro4 = "YOUR_METRO4",
#     EPA_region = "YOUR_REGION",
#     quadrant = "Q?"
#   ),
#   overlap_type = "horizontal"
# )

# ─────────────────────────────────────────────────────────────────────────────
# Colors
# ─────────────────────────────────────────────────────────────────────────────
quad_colors <- c(
  "Q1" = "#c0152a",  # dark red   (top-right)
  "Q2" = "#fcc5c0",  # light red  (top-left)
  "Q3" = "#c6dbef",  # light blue (bottom-left)
  "Q4" = "#1a6bb5",  # dark blue  (bottom-right)
  "Q5" = "#f1717e",  # mid red    (top-middle)
  "Q6" = "#6baed6"   # mid blue   (bottom-middle)
)

# ─────────────────────────────────────────────────────────────────────────────
# Treemap
# ─────────────────────────────────────────────────────────────────────────────
x_faces <- ifelse(region_levels2 == "Total", "bold", "plain")
y_faces <- ifelse(metro_levels2  == "Total", "bold", "plain")


treemap_plot <- ggplot(rect_df) +
  geom_rect(
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax,
      fill = fill_col
    ),
    color = "white",
    linewidth = 0.5
  ) +
  scale_fill_manual(
    values = c(quad_colors, "black" = "grey40")
  ) +
  guides(fill = "none") +
  geom_text(
    data = label_df,
    aes(
      x = x,
      y = y,
      label = label,
      size = txt_size,
      vjust = text_vjust
    ),
    lineheight = 0.9
  ) +
  scale_size_identity() +
  scale_x_continuous(
    breaks = seq_along(region_levels2),
    labels = gsub(" \\(", "\n(", region_levels2),
    expand = c(0.02, 0.02),
    position = "top"
  ) +
  scale_y_continuous(
    breaks = rev(seq_along(metro_levels2)),
    labels = metro_levels2,
    expand = c(0, 0)
  ) +
  coord_fixed(clip = "off") +
  labs(x = "", y = "") +
  theme_minimal(base_family = "Lato", base_size = 13) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid       = element_blank(),
    panel.border     = element_rect(color = "white", fill = NA, linewidth = 0.8),
    axis.text.x.top      = element_text(
      angle = 30, color="black",
      hjust = 0.5,
      vjust = 0.5,
      face = x_faces,
      size = 20, margin = margin(b = 5)
    ),
    axis.text.y      = element_text(
      size = 20, 
      color="black", 
      face = y_faces,
      angle = 30,        # ← Tilt y-axis text
      hjust = 0.5,
      vjust = 1,
      margin = margin(r = 0)
    ),
    axis.ticks       = element_blank(),
    plot.margin      = margin(t = 15, r = 15, b = 15, l = 5) 
  )

# ─────────────────────────────────────────────────────────────────────────────
# Save
# ─────────────────────────────────────────────────────────────────────────────
pdf(file = "/Users/k_sy_n_imac/Library/CloudStorage/Dropbox/GWU/Lab/FINESST/Research/Results/Quadrant_NO2_treemap_with_totals_sens.pdf",
    width = 26, height = 12)
print(treemap_plot)
dev.off()
















###############################################################
# Figure S5 ###############################################################
###############################################################
library(dplyr)
library(tidyr)
library(stringr)

base_dir <- "/Users/k_sy_n_imac/Library/CloudStorage/Dropbox/GWU/Lab/FINESST/Research/Data/ACS_5Y_Pop"
years <- 2019:2024

age_pop_all <- data.frame()

for (yr in years) {
  
  fp <- file.path(base_dir,
                  paste0("ACSDP5Y", yr),
                  paste0("ACSDP5Y", yr, ".DP05-Data.csv"))
  
  a <- read.csv(fp)
  a <- a[, 1:37]
  
  age_pop <- a %>%
    slice(-1) %>%
    transmute(
      GEOID = str_remove(GEO_ID, "^1400000US"),
      NAME,
      DP05_0009E, DP05_0010E, DP05_0011E, DP05_0012E,
      DP05_0013E, DP05_0014E, DP05_0015E,
      DP05_0016E, DP05_0017E
    ) %>%
    pivot_longer(
      cols = starts_with("DP05_"),
      names_to = "var",
      values_to = "pop"
    ) %>%
    mutate(
      pop = as.integer(pop),
      new_age_id_pop = recode(var,
                              "DP05_0009E" = 20,
                              "DP05_0010E" = 25,
                              "DP05_0011E" = 35,
                              "DP05_0012E" = 45,
                              "DP05_0013E" = 55,
                              "DP05_0014E" = 60,
                              "DP05_0015E" = 65,
                              "DP05_0016E" = 75,
                              "DP05_0017E" = 85
      ),
      year = yr
    ) %>%
    dplyr::select(GEOID, NAME, pop, new_age_id_pop, year)
  
  age_pop_all <- bind_rows(age_pop_all, age_pop)
}

age_pop_all <- age_pop_all %>%
  mutate(
    state = str_trim(str_extract(NAME, "[^,;]+$"))
  )

test <- age_pop_all %>%
  group_by(new_age_id_pop, year, state) %>%
  summarise(pop = sum(pop), .groups = "drop")

test_pct <- test %>%
  group_by(state, year) %>%
  mutate(pop_pct = pop / sum(pop) * 100) %>%
  ungroup()

test_pct <- subset(test_pct, !state %in% c( "Alaska","Hawaii", "Puerto Rico") )
test_pct <- test_pct %>%
  mutate(age_group = recode(new_age_id_pop,
                            `20` = "20–24",
                            `25` = "25–34",
                            `35` = "35–44",
                            `45` = "45–54",
                            `55` = "55–59",
                            `60` = "60–64",
                            `65` = "65–74",
                            `75` = "75–84",
                            `85` = "85+"
  ))

pdf(
  file = "/Users/k_sy_n_imac/Library/CloudStorage/Dropbox/GWU/Lab/FINESST/Research/Results/age_structure.pdf",
  width = 14,
  height = 8
)
par(family = "Lato")
ggplot(test_pct,
       aes(x = factor(year),
           y = pop_pct,
           fill = factor(age_group))) +
  geom_col(width = 0.8) +
  facet_wrap(~ state, nrow = 5, ncol = 10) +
  labs(
    x = "Year",
    y = "Population (%)",
    fill = "Age group"
  ) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold", size = 8),
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6),
    axis.title = element_text(size = 8),
    legend.position = "bottom"
  )
dev.off()













# Figure S4 ###############################################################
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggthemes)
library(scales)

# ── 1. Load & filter NOx ─────────────────────────────────────────────────────
df <- read_excel(
  path  = "/Users/k_sy_n_imac/Library/CloudStorage/Dropbox/GWU/Lab/FINESST/Research/Data/national_state_sector_2002_2024_caps_21feb2025_tons.xlsx",
  sheet = "State"
)
region_labels <- c(
  "Region 1\n(New England)",
  "Region 2\n(Northeast/Caribbean)",
  "Region 3\n(Mid-Atlantic)",
  "Region 4\n(Southeast)",
  "Region 5\n(Great Lakes)",
  "Region 6\n(South Central)",
  "Region 7\n(Heartland)",
  "Region 8\n(Mountains/Plains)",
  "Region 9\n(Pacific Southwest)",
  "Region 10\n(Northwest)"
)

recode_region <- function(x) {
  recode(x,
         "Region 1"  = "Region 1\n(New England)",
         "Region 2"  = "Region 2\n(Northeast/Caribbean)",
         "Region 3"  = "Region 3\n(Mid-Atlantic)",
         "Region 4"  = "Region 4\n(Southeast)",
         "Region 5"  = "Region 5\n(Great Lakes)",
         "Region 6"  = "Region 6\n(South Central)",
         "Region 7"  = "Region 7\n(Heartland)",
         "Region 8"  = "Region 8\n(Mountains/Plains)",
         "Region 9"  = "Region 9\n(Pacific Southwest)",
         "Region 10" = "Region 10\n(Northwest)"
  )
}

df_nox <- df %>%
  filter(Pollutant == "NOX") %>%
  dplyr::select(-c(5:21)) %>%
  mutate(across(starts_with("emissions"), as.numeric))

# ── 2. State → EPA region ─────────────────────────────────────────────────────
state_to_epa <- function(st) {
  case_when(
    st %in% c("CT","ME","MA","NH","RI","VT")           ~ 1L,
    st %in% c("NJ","NY","PR","VI")                     ~ 2L,
    st %in% c("DE","DC","MD","PA","VA","WV")           ~ 3L,
    st %in% c("AL","FL","GA","KY","MS","NC","SC","TN") ~ 4L,
    st %in% c("IL","IN","MI","MN","OH","WI")           ~ 5L,
    st %in% c("AR","LA","NM","OK","TX")                ~ 6L,
    st %in% c("IA","KS","MO","NE")                     ~ 7L,
    st %in% c("CO","MT","ND","SD","UT","WY")           ~ 8L,
    st %in% c("AZ","CA","HI","NV","AS","GU","MP")      ~ 9L,
    st %in% c("AK","ID","OR","WA")                     ~ 10L,
    TRUE ~ NA_integer_
  )
}


# ── 3. Sector mapping ─────────────────────────────────────────────────────────
df_mapped <- df_nox %>%
  mutate(
    sector = case_when(
      Sector %in% c("Agriculture - Livestock Waste",
                    "Agriculture - Crops & Livestock Dust",
                    "Agriculture - Fertilizer Application")
      ~ "Agriculture",
      Sector %in% c("Fuel Comb - Electric Generation - Biomass",
                    "Fuel Comb - Electric Generation - Coal",
                    "Fuel Comb - Electric Generation - Natural Gas",
                    "Fuel Comb - Electric Generation - Oil",
                    "Fuel Comb - Electric Generation - Other")
      ~ "Electric Generation",
      Sector %in% c("Fires - Prescribed Fires",
                    "Fires - Agricultural Field Burning")
      ~ "Prescribed Fires",
      Sector == "Fires - Wildfires"
      ~ "Wildfires",
      Sector %in% c("Mobile - Non-Road Equipment - Diesel",
                    "Mobile - Non-Road Equipment - Gasoline",
                    "Mobile - Non-Road Equipment - Other")
      ~ "Non-road Vehicles",
      Sector %in% c("Mobile - On-Road Diesel Heavy Duty Vehicles",
                    "Mobile - On-Road Diesel Light Duty Vehicles",
                    "Mobile - On-Road non-Diesel Heavy Duty Vehicles",
                    "Mobile - On-Road non-Diesel Light Duty Vehicles")
      ~ "On-road Vehicles",
      Sector %in% c("Mobile - Aircraft",
                    "Mobile - Commercial Marine Vessels",
                    "Mobile - Locomotives")
      ~ "Other Transport (Aircraft, Vessel, Locomotive)",
      Sector %in% c("Fuel Comb - Residential - Natural Gas",
                    "Fuel Comb - Residential - Oil",
                    "Fuel Comb - Residential - Other",
                    "Fuel Comb - Residential - Wood")
      ~ "Residential Combustion",
      Sector == "Industrial Processes - Oil & Gas Production"
      ~ "Oil & Gas Production",
      Sector %in% c("Fuel Comb - Comm/Institutional - Biomass",
                    "Fuel Comb - Comm/Institutional - Coal",
                    "Fuel Comb - Comm/Institutional - Natural Gas",
                    "Fuel Comb - Comm/Institutional - Oil",
                    "Fuel Comb - Comm/Institutional - Other",
                    "Fuel Comb - Industrial Boilers, ICEs - Biomass",
                    "Fuel Comb - Industrial Boilers, ICEs - Coal",
                    "Fuel Comb - Industrial Boilers, ICEs - Natural Gas",
                    "Fuel Comb - Industrial Boilers, ICEs - Oil",
                    "Fuel Comb - Industrial Boilers, ICEs - Other",
                    "Industrial Processes - Cement Manuf",
                    "Industrial Processes - Chemical Manuf",
                    "Industrial Processes - Ferrous Metals",
                    "Industrial Processes - Mining",
                    "Industrial Processes - NEC",
                    "Industrial Processes - Non-ferrous Metals",
                    "Industrial Processes - Petroleum Refineries",
                    "Industrial Processes - Pulp & Paper",
                    "Industrial Processes - Storage and Transfer")
      ~ "Industrial & Commercial Combustion",
      Sector %in% c("Miscellaneous Non-Industrial NEC",
                    "Bulk Gasoline Terminals", "Commercial Cooking",
                    "Gas Stations", "Dust - Construction Dust",
                    "Dust - Paved Road Dust", "Dust - Unpaved Road Dust")
      ~ "Other Sources (Miscellaneous)",
      Sector == "Waste Disposal" ~ "Waste",
      Sector %in% c("Solvent - Consumer & Commercial Solvent Use",
                    "Solvent - Degreasing", "Solvent - Dry Cleaning",
                    "Solvent - Graphic Arts",
                    "Solvent - Industrial Surface Coating & Solvent Use",
                    "Solvent - Non-Industrial Surface Coating")
      ~ "Solvent Use",
      TRUE ~ NA_character_
    ),
    EPA_region       = state_to_epa(State),
    EPA_region_label = factor(
      recode_region(paste0("Region ", EPA_region)),  # ← named labels
      levels = region_labels                          # ← enforces correct order
    )
  ) %>%
  filter(!is.na(EPA_region), !is.na(sector))


# ── 4. Pivot long + tons → ktons ─────────────────────────────────────────────
df_long <- df_mapped %>%
  pivot_longer(
    cols      = starts_with("emissions"),
    names_to  = "year",
    values_to = "emissions_kton"
  ) %>%
  mutate(
    year           = as.integer(gsub("emissions", "", year)),
    emissions_kton = emissions_kton / 1000
  )

# ── 5. Aggregate ──────────────────────────────────────────────────────────────
plot_df <- df_long %>%
  group_by(EPA_region, EPA_region_label, year, sector) %>%
  summarise(emissions = sum(emissions_kton, na.rm = TRUE), .groups = "drop")

total_df <- plot_df %>%
  group_by(EPA_region, EPA_region_label, year) %>%
  summarise(total = sum(emissions, na.rm = TRUE), .groups = "drop")

# ── 6. ONE global scale factor ────────────────────────────────────────────────
# Both axes share the same unit (ktons); we just compress total bars
# to fit visually within the sector line range.
global_max_sector <- max(plot_df$emissions, na.rm = TRUE)
global_max_total  <- max(total_df$total,    na.rm = TRUE)

# total_scaled = total / scale_factor  → lives on same axis as sector lines
# sec_axis transform = ~ . * scale_factor  → converts back to real ktons
scale_factor <- global_max_total / global_max_sector

total_df <- total_df %>%
  mutate(total_scaled = total / scale_factor)

# plot setup
all_years <- sort(unique(plot_df$year)) 

sector_levels <- c(
  # Mobile sources
  "On-road Vehicles",
  "Non-road Vehicles",
  "Other Transport (Aircraft, Vessel, Locomotive)",
  
  # Stationary combustion 
  "Electric Generation",
  "Industrial & Commercial Combustion",
  "Oil & Gas Production",
  "Residential Combustion",
  
  # Fire-related
  "Wildfires",
  "Prescribed Fires",
  
  # Agriculture
  "Agriculture",
  
  # Minor / miscellaneous
  "Waste",
  "Solvent Use",
  "Other Sources (Miscellaneous)"
)

# ── Colors matched to group (similar hues within group) ───────────────────────
sector_colors <- c(
  # Mobile → blues
  "On-road Vehicles"                               = "#08519c",  # dark blue
  "Non-road Vehicles"                              = "#2171b5",  # mid blue
  "Other Transport (Aircraft, Vessel, Locomotive)" = "#6baed6",  # light blue
  
  # Stationary combustion → reds/oranges
  "Electric Generation"                            = "#e41a1c",  # dark red
  "Industrial & Commercial Combustion"             = "#a65628",  # red
  "Oil & Gas Production"                           = "#fb6a4a",  # orange-red
  "Residential Combustion"                         = "#f781bf",  # light salmon
  
  # Fires → browns/amber
  "Wildfires"                                      = "#e6ab02",  # dark brown
  "Prescribed Fires"                               = "#878500",  # amber-brown
  
  # Agriculture → green
  "Agriculture"                                    = "#4daf4a",  # green
  
  # Miscellaneous → greys/neutral
  "Waste"                                          = "#737373",  # dark grey
  "Solvent Use"                                    = "#984ea3",   # light grey
  "Other Sources (Miscellaneous)"                  = "#252525"  # mid grey
)


plot_df <- plot_df %>%
  mutate(sector = factor(sector, levels = sector_levels))


# ── 7. Plot ───────────────────────────────────────────────────────────────────
p <- ggplot() +
  geom_col(
    data  = total_df,
    aes(x = year, y = total_scaled),
    fill  = "#e31a1c",
    alpha = 0.1,
    width = 0.5
  ) +
  geom_line(
    data      = plot_df,
    aes(x = year, y = emissions,
        color = sector,
        group = sector),
    linewidth = 0.8
  ) +
  facet_wrap(~ EPA_region_label, ncol = 5) +   # no free_y → shared scale
  scale_y_continuous(
    name   = "Sector NOx emissions (ktons)",
    labels = label_comma(),
    sec.axis = sec_axis(
      transform = ~ . * scale_factor,           # one global transform
      name      = "Total NOx emissions (ktons)",
      labels    = label_comma(accuracy = 1)
    )
  ) +
  scale_x_continuous(
    breaks = all_years,                              # ← every year
    labels = as.character(all_years)                 # ← plain numbers, no abbrev
  ) +
  scale_color_manual(
    values = sector_colors,
    breaks = sector_levels,
    name   = "Sector"
  ) +  labs(
    #title    = "NOx Emissions by Sector (lines) + Total (bars) across EPA Regions",
    #subtitle = "Left axis = sector emissions  |  Right axis = total emissions  |  Both in ktons, shared scale across all regions",
    x        = "Year",
    color    = "Sector"
  ) +
  theme_minimal(base_family = "Lato") +
  theme(
    legend.position    = "bottom",
    legend.title       = element_text(face = "bold"),
    strip.text         = element_text(face = "bold", size = 9),
    axis.title.y.left  = element_text(size = 9, margin = margin(r = 6)),
    axis.title.y.right = element_text(color = "#e31a1c", size = 9,
                                      margin = margin(l = 6)),
    axis.text.y.right  = element_text(color = "#e31a1c"),
    plot.title         = element_text(face = "bold", size = 13),
    plot.subtitle      = element_text(size = 8, color = "grey45")
  )

pdf(
  file   = "/Users/k_sy_n_imac/Library/CloudStorage/Dropbox/GWU/Lab/FINESST/Research/Results/nox_emission.pdf",
  width  = 16,
  height = 9
)
print(p)
dev.off()












#########################################################################################################################################
# Figure S5  #######################################################################################################################################
######################################################################################################################################################
list_year <- c(2019, 2020, 2021, 2022, 2023, 2024) #data is too large to handle all the years and causes at the same time. Select each to proceed
list_cause <- c("All causes", "Cardiovascular diseases", "Respiratory diseases", "Ischemic heart disease", "Tracheal, bronchus, and lung cancer") 
i1 <- 1
i2 <- 6
j <- 1
no2_ind <- subset(no2_tract_df, year %in% c(list_year[i1], list_year[i2]))
db_ind <- subset(final3, year %in% c(list_year[i1], list_year[i2]) & cause_name==list_cause[j])
final4 <- merge(no2_ind, db_ind, by=c("GEOID", "year"))
final4 <- merge(final4, xxxx, by="GEOID", all.x=TRUE)
final4 <- final4 %>%
  mutate(
    MSA_title = if_else(
      metro4 == "Non-urban",
      paste0("Non_urban_", state_name.x),
      MSA_title
    )
  )

group_df <- final4 %>% group_by(MSA_title, metro4, EPA_region) %>%
  summarize(n_2024 = n(), .groups = "drop")
msa_regions <- group_df %>% count(MSA_title) %>% filter(n > 1)
majority_regions <- final4 %>%
  group_by(MSA_title, EPA_region) %>%
  summarise(total_pop = sum(pop_over20, na.rm = TRUE), .groups = "drop") %>%
  group_by(MSA_title) %>%
  slice_max(total_pop, n = 1) %>%  # ← Gets region with HIGHEST population
  ungroup() %>%
  dplyr::select(MSA_title, EPA_region_new = EPA_region)
final4 <- final4 %>%
  left_join(majority_regions, by = "MSA_title") %>%
  mutate(EPA_region = if_else(!is.na(EPA_region_new), 
                              EPA_region_new, 
                              EPA_region))

# Calculation
# Sub-group
# 2019
x <- subset(final4, final4$year==list_year[i1])
group_df_2019 <- x %>%
  group_by(MSA_title, EPA_region, metro4) %>%
  summarize(case_pe_2019 = sum(case_pe, na.rm = TRUE), pop_2019 = sum(total_dc2020, na.rm = TRUE), 
            rate_pe_2019 = case_pe_2019/pop_2019 * 100000, n_2019 = n())
xx <- subset(final4, final4$year==list_year[i2])
group_df_2024 <- xx %>%
  group_by(MSA_title, EPA_region, metro4) %>%
  summarize(case_pe_2024 = sum(case_pe, na.rm = TRUE), pop_2024 = sum(total_dc2020, na.rm = TRUE), 
            rate_pe_2024 = case_pe_2024/pop_2024 * 100000, n_2024 = n())
group_df <- merge(
  group_df_2019,
  group_df_2024,
  by = c("MSA_title", "EPA_region", "metro4")
)

group_df$pct_change <- (group_df$rate_pe_2024-group_df$rate_pe_2019)/group_df$rate_pe_2019*100
group_df$pct_change <- ifelse(is.na(group_df$pct_change), 0, group_df$pct_change) # Non-urban cases from zero to zero 
group_df$pct_change <- ifelse(group_df$pct_change>=100, 100, group_df$pct_change) # Unify values over 100 as 100



################################################################################
# Prepare raw (MSA-level) data
################################################################################
region_labels <- c(
  " Region 1\n (New England)",
  " Region 2\n (Northeast/Caribbean)",
  " Region 3\n (Mid-Atlantic)",
  " Region 4\n (Southeast)",
  " Region 5\n (Great Lakes)",
  " Region 6\n (South Central)",
  " Region 7\n (Heartland)",
  " Region 8\n (Mountains/Plains)",
  " Region 9\n (Pacific Southwest)",
  " Region 10\n (Northwest)")

recode_region <- function(x) {
  recode(x,
         "Region 1"  = " Region 1\n (New England)",
         "Region 2"  = " Region 2\n (Northeast/Caribbean)",
         "Region 3"  = " Region 3\n (Mid-Atlantic)",
         "Region 4"  = " Region 4\n (Southeast)",
         "Region 5"  = " Region 5\n (Great Lakes)",
         "Region 6"  = " Region 6\n (South Central)",
         "Region 7"  = " Region 7\n (Heartland)",
         "Region 8"  = " Region 8\n (Mountains/Plains)",
         "Region 9"  = " Region 9\n (Pacific Southwest)",
         "Region 10" = " Region 10\n (Northwest)")
}

region_colors <- c(
  " Region 1\n (New England)"         = "#1B2F70",
  " Region 2\n (Northeast/Caribbean)" = "#2C4FD7",
  " Region 3\n (Mid-Atlantic)"        = "#5A8AC6",
  " Region 4\n (Southeast)"           = "#8BC8C8",
  " Region 5\n (Great Lakes)"         = "#A1D99B",
  " Region 6\n (South Central)"       = "#D9C89B",
  " Region 7\n (Heartland)"           = "#DDAA66",
  " Region 8\n (Mountains/Plains)"    = "#F2A6A0",
  " Region 9\n (Pacific Southwest)"   = "#7A5AC9",
  " Region 10\n (Northwest)"          = "#4A0066")

plot_df <- group_df %>%
  mutate(
    metro4     = factor(metro4, levels = c("Large MSA", "Small MSA",
                                           "MicroSA", "Non-urban")),
    EPA_region = factor(recode_region(EPA_region), levels = region_labels)
  )

################################################################################
# Plot: raw scatter
################################################################################

pdf(
  file = "/Users/k_sy_n_imac/Library/CloudStorage/Dropbox/GWU/Lab/FINESST/Research/Results/Quadrant_NO2burden_raw.pdf",
  width = 13,
  height = 10
)
par(family = "Lato")

ggplot(plot_df, aes(x = pct_change, y = rate_pe_2024,
                    color = EPA_region, shape = metro4)) +
  
  # Reference lines
  geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.8) +   # threshold
  geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.8) +
  
  geom_vline(xintercept = -5, linetype = "dashed", color = "gray40", linewidth = 0.5) +
  geom_vline(xintercept = 5, linetype = "dashed", color = "gray40", linewidth = 0.5) +
  
  # Label
  annotate("text",
           x = 8.5, y = -0.8,
           label = "5",
           size = 5,
           #fontface = "bold",
           color = "black") +
  annotate("text",
           x = -8.5, y = -0.8,
           label = "-5",
           size = 5,
           #fontface = "bold",
           color = "black") +
  
  # Raw points (one per MSA)
  geom_point(size = 2.5, alpha = 0.8) +
  
  scale_color_manual(values = region_colors) +
  scale_shape_manual(values = c(
    "Large MSA" = 16,
    "Small MSA" = 17,
    "MicroSA" = 15,
    "Non-urban" = 8
  )) +
  labs(
    x = bquote(bold("% Change NO"[bolditalic(2)] ~ "-attributable deaths per 100,000 between 2019 and 2024")),
    y = bquote(bold("NO"[bolditalic(2)] ~ "-attributable deaths in 2024 (per 100,000)")),
    shape = "Urbanicity",
    color = "EPA Region"
  ) + 
  guides(
    shape = guide_legend(
      order = 1,
      nrow = 4,
      override.aes = list(size = 3)
    ),
    color = guide_legend(
      order = 2,
      nrow = 10,
      theme = theme(
        legend.text = element_text(
          size = 13,
          margin = margin(b = 8)
        )
      )
    )
  ) +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 13),
    axis.text = element_text(size = 13),
    axis.title.x = element_text(size = 15, margin = margin(t = 15, b = 15)),
    axis.title.y = element_text(size = 15, margin = margin(r = 15, l = 15)),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 13),
    panel.grid.minor = element_blank()  ,
    panel.grid.major = element_line(color = "gray95"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  ) 
dev.off()



# ============================================================================
# Section 1: EPA REGION STRATIFICATION
# ============================================================================

# Calculate summary statistics by EPA region
group_summary_by_region_burden <- group_df %>%
  group_by(EPA_region) %>%
  summarize(
    n_msa = n(),
    rate_pe = mean(rate_pe_2024, na.rm = TRUE),
    var_rate = var(rate_pe_2024, na.rm = TRUE),
    pct_chg = mean(pct_change, na.rm = TRUE),
    var_pct_change = var(pct_change, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    se_rate = sqrt(var_rate / n_msa),
    se_pct_chg = sqrt(var_pct_change / n_msa),
    x0 = pct_chg,
    y0 = rate_pe,
    rx = 1.96 * se_pct_chg,
    ry = 1.96 * se_rate
  ) %>%
  mutate(
    EPA_region = factor(
      EPA_region,
      levels = c("Region 1", "Region 2", "Region 3", "Region 4", "Region 5",
                 "Region 6", "Region 7", "Region 8", "Region 9", "Region 10")
    )
  )

# Generate horizontal thin diamonds for x-axis uncertainty
diamond_region_burden_x <- bind_rows(lapply(seq_len(nrow(group_summary_by_region_burden)), function(i) {
  tibble(
    EPA_region = group_summary_by_region_burden$EPA_region[i],
    diamond_id = paste0(group_summary_by_region_burden$EPA_region[i], "_x"),
    x = c(group_summary_by_region_burden$x0[i] - group_summary_by_region_burden$rx[i],
          group_summary_by_region_burden$x0[i],
          group_summary_by_region_burden$x0[i] + group_summary_by_region_burden$rx[i],
          group_summary_by_region_burden$x0[i]),
    y = c(group_summary_by_region_burden$y0[i],
          group_summary_by_region_burden$y0[i] + 0.02,
          group_summary_by_region_burden$y0[i],
          group_summary_by_region_burden$y0[i] - 0.02)
  )
}))

# Generate vertical thin diamonds for y-axis uncertainty
diamond_region_burden_y <- bind_rows(lapply(seq_len(nrow(group_summary_by_region_burden)), function(i) {
  tibble(
    EPA_region = group_summary_by_region_burden$EPA_region[i],
    diamond_id = paste0(group_summary_by_region_burden$EPA_region[i], "_y"),
    x = c(group_summary_by_region_burden$x0[i],
          group_summary_by_region_burden$x0[i] + 0.2,
          group_summary_by_region_burden$x0[i],
          group_summary_by_region_burden$x0[i] - 0.2),
    y = c(group_summary_by_region_burden$y0[i] - group_summary_by_region_burden$ry[i],
          group_summary_by_region_burden$y0[i],
          group_summary_by_region_burden$y0[i] + group_summary_by_region_burden$ry[i],
          group_summary_by_region_burden$y0[i])
  )
}))

diamond_region_burden <- bind_rows(diamond_region_burden_x, diamond_region_burden_y)

# Calculate national summary
nation_burden_df <- group_df %>%
  summarize(
    n_msa = n(),
    rate_pe = mean(rate_pe_2024, na.rm = TRUE),
    var_rate = var(rate_pe_2024, na.rm = TRUE),
    pct_chg = mean(pct_change, na.rm = TRUE),
    var_pct_change = var(pct_change, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    se_rate = sqrt(var_rate / n_msa),
    se_pct_chg = sqrt(var_pct_change / n_msa),
    x0 = pct_chg,
    y0 = rate_pe,
    rx = 1.96 * se_pct_chg,
    ry = 1.96 * se_rate
  )

# Define region colors
region_colors <- c(
  "Region 1"  = "#1B2F70",
  "Region 2"  = "#2C4FD7",
  "Region 3"  = "#5A8AC6",
  "Region 4"  = "#8BC8C8",
  "Region 5"  = "#A1D99B",
  "Region 6"  = "#D9C89B",
  "Region 7"  = "#DDAA66",
  "Region 8"  = "#F2A6A0",
  "Region 9"  = "#7A5AC9",
  "Region 10" = "#4A0066"
)

# Create EPA Region plot
pdf(file = "/Users/k_sy_n_imac/Library/CloudStorage/Dropbox/GWU/Lab/FINESST/Research/Results/Quadrant_NO2burden_EPA_Region.pdf", 
    width = 9, height = 9)

ggplot() +
  # Reference lines
  geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.5) +
  geom_vline(xintercept = -5, linetype = "dashed", color = "gray40", linewidth = 0.5) +
  geom_vline(xintercept = 5, linetype = "dashed", color = "gray40", linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.5) +
  
  # Thin diamond uncertainty regions
  geom_polygon(data = diamond_region_burden,
               aes(x = x, y = y, fill = EPA_region, group = diamond_id),
               alpha = 1, linewidth = 1, color = NA) +
  
  # Center points
  geom_point(data = group_summary_by_region_burden,
             aes(x = x0, y = y0, color = EPA_region),
             size = 5, stroke = 1.2) +
  
  # National summary diamonds
  geom_polygon(data = {
    h_diamond <- tibble(
      x = c(nation_burden_df$x0 - nation_burden_df$rx,
            nation_burden_df$x0,
            nation_burden_df$x0 + nation_burden_df$rx,
            nation_burden_df$x0),
      y = c(nation_burden_df$y0,
            nation_burden_df$y0 + 0.03,
            nation_burden_df$y0,
            nation_burden_df$y0 - 0.03),
      diamond_id = "national_x"
    )
    v_diamond <- tibble(
      x = c(nation_burden_df$x0,
            nation_burden_df$x0 + 0.3,
            nation_burden_df$x0,
            nation_burden_df$x0 - 0.3),
      y = c(nation_burden_df$y0 - nation_burden_df$ry,
            nation_burden_df$y0,
            nation_burden_df$y0 + nation_burden_df$ry,
            nation_burden_df$y0),
      diamond_id = "national_y"
    )
    bind_rows(h_diamond, v_diamond)
  },
  aes(x = x, y = y, group = diamond_id),
  alpha = 0.7, linewidth = 1, color = NA, fill = "#B11226") +
  
  # National summary point
  geom_point(data = nation_burden_df,
             aes(x = x0, y = y0),
             color = "#B11226", size = 7, shape = 18) +
  geom_text(data = nation_burden_df,
            aes(x = x0, y = y0),
            label = "National Summary",
            color = "#B11226", size = 7.5, fontface = "bold",
            vjust = 1.8, hjust = 0.1) +
  
  scale_color_manual(values = region_colors) +
  scale_fill_manual(values = region_colors) +
  
  labs(
    x = bquote(bold("% Change from 2019")),
    y = bquote(bold("NO"[bolditalic(2)] ~ "-attributable deaths in 2024")),
    color = "EPA Region",
    fill = "EPA Region"
  ) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 22),
    axis.title.x = element_text(size = 20, margin = margin(t = 15, b = 15)),
    axis.title.y = element_text(size = 25, margin = margin(r = 15, l = 15)),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray95"),
    legend.position = "none"
  ) +
  xlim(-45, 25) +
  ylim(0, 14)
dev.off()





# ============================================================================
# Section 2: URBANICITY STRATIFICATION
# ============================================================================

# Calculate summary statistics by urbanicity
group_summary_by_metro_burden <- group_df %>%
  group_by(metro4) %>%
  summarize(
    n_msa = n(),
    rate_pe = mean(rate_pe_2024, na.rm = TRUE),
    var_rate = var(rate_pe_2024, na.rm = TRUE),
    pct_chg = mean(pct_change, na.rm = TRUE),
    var_pct_change = var(pct_change, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    se_rate = sqrt(var_rate / n_msa),
    se_pct_chg = sqrt(var_pct_change / n_msa),
    x0 = pct_chg,
    y0 = rate_pe,
    rx = 1.96 * se_pct_chg,
    ry = 1.96 * se_rate
  ) %>%
  mutate(
    metro4 = factor(metro4, levels = c("Large MSA", "Small MSA", "MicroSA", "Non-urban"))
  )

# Generate horizontal thin diamonds for x-axis uncertainty
diamond_metro_burden_x <- bind_rows(lapply(seq_len(nrow(group_summary_by_metro_burden)), function(i) {
  tibble(
    metro4 = group_summary_by_metro_burden$metro4[i],
    diamond_id = paste0(group_summary_by_metro_burden$metro4[i], "_x"),
    x = c(group_summary_by_metro_burden$x0[i] - group_summary_by_metro_burden$rx[i],
          group_summary_by_metro_burden$x0[i],
          group_summary_by_metro_burden$x0[i] + group_summary_by_metro_burden$rx[i],
          group_summary_by_metro_burden$x0[i]),
    y = c(group_summary_by_metro_burden$y0[i],
          group_summary_by_metro_burden$y0[i] + 0.07,
          group_summary_by_metro_burden$y0[i],
          group_summary_by_metro_burden$y0[i] - 0.07)
  )
}))

# Generate vertical thin diamonds for y-axis uncertainty
diamond_metro_burden_y <- bind_rows(lapply(seq_len(nrow(group_summary_by_metro_burden)), function(i) {
  tibble(
    metro4 = group_summary_by_metro_burden$metro4[i],
    diamond_id = paste0(group_summary_by_metro_burden$metro4[i], "_y"),
    x = c(group_summary_by_metro_burden$x0[i],
          group_summary_by_metro_burden$x0[i] + 0.16,
          group_summary_by_metro_burden$x0[i],
          group_summary_by_metro_burden$x0[i] - 0.16),
    y = c(group_summary_by_metro_burden$y0[i] - group_summary_by_metro_burden$ry[i],
          group_summary_by_metro_burden$y0[i],
          group_summary_by_metro_burden$y0[i] + group_summary_by_metro_burden$ry[i],
          group_summary_by_metro_burden$y0[i])
  )
}))

diamond_metro_burden <- bind_rows(diamond_metro_burden_x, diamond_metro_burden_y)


# Create Urbanicity plot
pdf(file = "/Users/k_sy_n_imac/Library/CloudStorage/Dropbox/GWU/Lab/FINESST/Research/Results/Quadrant_NO2burden_urbanicity.pdf", 
    width = 9, height = 9)

ggplot() +
  # Reference lines
  geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.5) +
  geom_vline(xintercept = -5, linetype = "dashed", color = "gray40", linewidth = 0.5) +
  geom_vline(xintercept = 5, linetype = "dashed", color = "gray40", linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.5) +
  
  # Thin diamond uncertainty regions (all black)
  geom_polygon(data = diamond_metro_burden,
               aes(x = x, y = y, group = diamond_id),
               fill = "black", alpha = 0.7, linewidth = 1, color = NA) +
  
  # Center points
  geom_point(data = group_summary_by_metro_burden,
             aes(x = x0, y = y0, shape = metro4),
             color = "black", size = 5, stroke = 1.2) +
  
  # National summary diamonds
  geom_polygon(data = {
    h_diamond <- tibble(
      x = c(nation_burden_df$x0 - nation_burden_df$rx,
            nation_burden_df$x0,
            nation_burden_df$x0 + nation_burden_df$rx,
            nation_burden_df$x0),
      y = c(nation_burden_df$y0,
            nation_burden_df$y0 + 0.08,
            nation_burden_df$y0,
            nation_burden_df$y0 - 0.08),
      diamond_id = "national_x"
    )
    v_diamond <- tibble(
      x = c(nation_burden_df$x0,
            nation_burden_df$x0 + 0.18,
            nation_burden_df$x0,
            nation_burden_df$x0 - 0.18),
      y = c(nation_burden_df$y0 - nation_burden_df$ry,
            nation_burden_df$y0,
            nation_burden_df$y0 + nation_burden_df$ry,
            nation_burden_df$y0),
      diamond_id = "national_y"
    )
    bind_rows(h_diamond, v_diamond)
  },
  aes(x = x, y = y, group = diamond_id),
  alpha = 0.7, linewidth = 1, color = NA, fill = "#B11226") +
  
  # National summary point
  geom_point(data = nation_burden_df,
             aes(x = x0, y = y0),
             color = "#B11226", size = 7, shape = 18) +
  geom_text(data = nation_burden_df,
            aes(x = x0, y = y0),
            label = "National Summary",
            color = "#B11226", size = 7.5, fontface = "bold",
            vjust = -1.5, hjust = 0.1) +
  
  scale_shape_manual(values = c(
    "Large MSA" = 16,
    "Small MSA" = 17,
    "MicroSA" = 15,
    "Non-urban" = 8
  )) +
  
  labs(
    x = bquote(bold("% Change from 2019")),
    y = bquote(bold("NO"[bolditalic(2)] ~ "-attributable deaths in 2024")),
    shape = "Urbanicity"
  ) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 22),
    axis.title.x = element_text(size = 20, margin = margin(t = 15, b = 15)),
    axis.title.y = element_text(size = 25, margin = margin(r = 15, l = 15)),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray95"),
    legend.position = "none"
  ) +
  xlim(-45, 25) +
  ylim(0, 14)
dev.off()




# Descriptive
plot_df %>%
  group_by(metro4) %>%
  summarise(n_zero = sum(rate_pe_2024 < 1, na.rm = TRUE), n = n(), pct = n_zero/n*100)
plot_df %>%
  group_by(metro4) %>%
  summarise(n_zero = sum(pct_change >= 100, na.rm = TRUE), n = n(), pct = n_zero/n*100)



















########################################################################################
########################################################################################
# Disparity analysis  ##################################################################
########################################################################################
########################################################################################

# ============================================================
# Figure 4 & Figure S6: Multi-panel PNG:
#   Top row:    LA, Chicago, Atlanta, DMV
#   Middle row: full map + inset bivar legend
#   Bottom row: Dallas, Houston, NYC, Miami
# ============================================================
# -----------------------------
# 0) Your data prep (same logic)
# -----------------------------
list_year  <- c(2019, 2020, 2021, 2022, 2023, 2024)
list_cause <- c("All causes", "Cardiovascular diseases", "Respiratory diseases",
                "Ischemic heart disease", "Tracheal, bronchus, and lung cancer")
i <- 6
j <- 1
no2_ind <- subset(no2_tract_df, year == list_year[i])
db_ind  <- subset(final3,       year == list_year[i] & cause_name == list_cause[j])
final4 <- merge(no2_ind, db_ind, by = "GEOID")
final4 <- merge(final4, xxxx,   by = "GEOID", all.x = TRUE)
# If merge() stripped sf class, put it back (assuming geometry column exists)
if (!inherits(final4, "sf")) {
  final4 <- sf::st_as_sf(final4)}  # works if final4 already has an sf geometry column

# 5 x 5 bivariate variables
final4 <- final4 %>%
  mutate(
    white_quintile = dplyr::if_else(!is.na(pct_white) & !is.na(surface_no2),
                                    dplyr::ntile(pct_white, 5), NA_integer_),
    no2_quintile_race = dplyr::if_else(!is.na(pct_white) & !is.na(surface_no2),
                                       dplyr::ntile(surface_no2, 5), NA_integer_),
    white_no2_group = dplyr::if_else(!is.na(white_quintile) & !is.na(no2_quintile_race),
                                     paste0("W", white_quintile, "_N", no2_quintile_race),NA_character_),
    
    edu_quintile = dplyr::if_else(!is.na(pct_high_over) & !is.na(surface_no2),
                                  dplyr::ntile(pct_high_over, 5), NA_integer_),
    no2_quintile_edu = dplyr::if_else(!is.na(pct_high_over) & !is.na(surface_no2),
                                      dplyr::ntile(surface_no2, 5), NA_integer_),
    edu_no2_group = dplyr::if_else(!is.na(edu_quintile) & !is.na(no2_quintile_edu),
                                   paste0("W", edu_quintile, "_N", no2_quintile_edu),NA_character_),
    
    inc_quintile = dplyr::if_else(!is.na(mean_household_income) & !is.na(surface_no2),
                                  dplyr::ntile(mean_household_income, 5), NA_integer_),
    no2_quintile_inc = dplyr::if_else(!is.na(mean_household_income) & !is.na(surface_no2),
                                      dplyr::ntile(surface_no2, 5), NA_integer_),
    inc_no2_group = dplyr::if_else(!is.na(inc_quintile) & !is.na(no2_quintile_inc),
                                   paste0("W", inc_quintile, "_N", no2_quintile_inc),NA_character_),
    
    emp_quintile = dplyr::if_else(!is.na(employ_rate) & !is.na(surface_no2),
                                  dplyr::ntile(employ_rate, 5), NA_integer_),
    no2_quintile_emp = dplyr::if_else(!is.na(employ_rate) & !is.na(surface_no2),
                                      dplyr::ntile(surface_no2, 5), NA_integer_),
    emp_no2_group = dplyr::if_else(!is.na(emp_quintile) & !is.na(no2_quintile_emp),
                                   paste0("W", emp_quintile, "_N", no2_quintile_emp),NA_character_),
    
    pov_rate_reverse = 100-pov_rate,
    pov_quintile = dplyr::if_else(!is.na(pov_rate_reverse) & !is.na(surface_no2),
                                  dplyr::ntile(pov_rate_reverse, 5), NA_integer_),
    no2_quintile_pov = dplyr::if_else(!is.na(pov_rate_reverse) & !is.na(surface_no2),
                                      dplyr::ntile(surface_no2, 5), NA_integer_),
    pov_no2_group = dplyr::if_else(!is.na(pov_quintile) & !is.na(no2_quintile_pov),
                                   paste0("W", pov_quintile, "_N", no2_quintile_pov),NA_character_),
  )

# -----------------------------
# 1) Functions
# -----------------------------
# color palette
make_pal_WN <- function(dim = 5,
                        W_lo = "#F28E2B",
                        W_hi = "#008080",
                        N_lo = "#FFFFFF",
                        N_hi = "#000000",
                        w_W  = 0.45) {
  
  ramp_W <- grDevices::colorRampPalette(c(W_lo, W_hi))(dim)
  ramp_N <- grDevices::colorRampPalette(c(N_lo, N_hi))(dim)
  
  blend <- function(c1, c2, a = 0.5) {
    c1 <- grDevices::col2rgb(c1); c2 <- grDevices::col2rgb(c2)
    grDevices::rgb((1-a)*c1[1,] + a*c2[1,],
                   (1-a)*c1[2,] + a*c2[2,],
                   (1-a)*c1[3,] + a*c2[3,],
                   maxColorValue = 255)}
  
  cols <- character(dim * dim)
  nms  <- character(dim * dim)
  k <- 1
  for (w in 1:dim) {
    for (n in 1:dim) {
      cols[k] <- blend(ramp_N[n], ramp_W[w], a = w_W)
      nms[k]  <- sprintf("W%d_N%d", w, n)
      k <- k + 1}}
  stats::setNames(cols, nms)
}
pal_vec     <- make_pal_WN(dim = 5)
# bivariate legend
bivar_legend <- function(pal_vec, dim = 5, ylab,
                         reverse_y_for_minority = TRUE) {
  
  df <- expand.grid(N = 1:dim, W = 1:dim)
  df$key <- sprintf("W%d_N%d", df$W, df$N)
  df$Wy  <- if (reverse_y_for_minority) dim - df$W + 1 else df$W
  
  ggplot(df, aes(x = N, y = Wy, fill = key)) +
    geom_tile() +
    scale_fill_manual(values = pal_vec, guide = "none") +
    scale_x_continuous(expand = c(0,0), breaks = NULL) +
    scale_y_continuous(expand = c(0,0), breaks = NULL) +
    coord_fixed(clip = "off") +
    theme_void() +
    annotate("segment",
             x = 0.1, xend = dim + 0.7,
             y = 0.1, yend = 0.1,
             arrow = arrow(type = "closed", length = unit(5, "pt")),
             size = 2) +
    annotate("segment",
             x = 0.1, xend = 0.1,
             y = 0.1, yend = dim + 0.7,
             arrow = arrow(type = "closed", length = unit(5, "pt")),
             size = 2) +
    annotate("text",
             x = (dim + 1)/2, y = -0.7,
             label = expression(bold(NO[2]~concentration)), vjust = 0.5, size = 14, fontface="bold") +
    annotate("text",
             x = -0.7, y = (dim + 1)/2,
             label = ylab, angle = 90, hjust = 0.5, size = 14, fontface="bold")
}
# zoom in panels
crop_to_bbox <- function(sf_df, xlim, ylim) {
  bb_mat <- matrix(
    c(xlim[1], ylim[1],
      xlim[2], ylim[1],
      xlim[2], ylim[2],
      xlim[1], ylim[2],
      xlim[1], ylim[1]),
    ncol = 2, byrow = TRUE
  )
  bb_sfc <- st_sfc(st_polygon(list(bb_mat)), crs = st_crs(sf_df))
  st_intersection(sf_df, bb_sfc)
}

make_zoom_panel <- function(sf_df, pal_vec, title, xlim, ylim) {
  sf_crop <- crop_to_bbox(sf_df, xlim, ylim)
  
  ggplot() +
    geom_sf(data = sf_crop, aes(fill = white_no2_group), color = NA) +
    scale_fill_manual(values = pal_vec, guide = "none") +
    theme_void() +
    ggtitle(paste0("  ", title)) +     theme(
      plot.margin = margin(2, 2, 2, 2),
      text       = element_text(size = 34, face="bold"),
      legend.text = element_text(size = 34, face="bold"),
      legend.title = element_text(size = 34, face="bold")
    )
}

metro_bbox <- list(
  LA           = list(xlim = c(-118.95, -117.55), ylim = c(33.55, 34.45)),
  Chicago      = list(xlim = c( -88.35,  -87.25), ylim = c(41.45, 42.25)),
  NYC          = list(xlim = c( -74.35,  -73.55), ylim = c(40.45, 41.05)),
  Philadelphia = list(xlim = c( -75.40,  -74.90), ylim = c(39.80, 40.20)),
  DMV          = list(xlim = c( -77.60,  -76.70), ylim = c(38.55, 39.20)),
  Dallas       = list(xlim = c( -97.15,  -96.35), ylim = c(32.55, 33.05)),
  Houston      = list(xlim = c( -95.95,  -94.95), ylim = c(29.45, 30.15)),
  Miami        = list(xlim = c( -80.55,  -80.00), ylim = c(25.55, 26.15))
)
metro_names <- c("LA", "Chicago", "NYC", "Philadelphia", "DMV", "Dallas", "Houston", "Miami")
metro_full <- c(
  LA           = "Los Angeles",
  Chicago      = "Chicago",
  NYC          = "New York City",
  Philadelphia = "Philadelphia",
  DMV          = "DC-Maryland-Virginia",
  Dallas       = "Dallas",
  Houston      = "Houston",
  Miami        = "Miami"
)
panel_letters <- LETTERS[1:8]
names(panel_letters) <- metro_names

rect_df <- do.call(rbind, lapply(metro_names, function(nm) {
  bb <- metro_bbox[[nm]]
  data.frame(
    metro = nm,
    xmin  = bb$xlim[1],
    xmax  = bb$xlim[2],
    ymin  = bb$ylim[1],
    ymax  = bb$ylim[2],
    label = panel_letters[nm]
  )
}))

rect_sf <- st_as_sf(rect_df, coords = c("xmin", "ymin"), crs = st_crs(final4), remove = FALSE)
rect_sf$geometry <- st_sfc(lapply(seq_len(nrow(rect_df)), function(k) {
  with(rect_df[k, ], st_polygon(list(matrix(
    c(xmin, ymin,
      xmax, ymin,
      xmax, ymax,
      xmin, ymax,
      xmin, ymin),
    ncol = 2, byrow = TRUE
  ))))
}), crs = st_crs(final4))
rect_sf <- st_as_sf(rect_sf)

label_out_df <- rect_df %>%
  mutate(
    label = panel_letters[metro],
    x_lab = xmin + 0.2,
    y_lab = ymax + 0.4)
label_out_pts <- st_as_sf(label_out_df, coords = c("x_lab", "y_lab"),
                          crs = st_crs(final4))

label_df <- do.call(rbind, lapply(metro_names, function(nm) {
  bb <- metro_bbox[[nm]]
  data.frame(
    metro = nm,
    label = panel_letters[nm],
    x = mean(bb$xlim),
    y = mean(bb$ylim)
  )
}))
label_pts <- st_as_sf(label_df, coords = c("x", "y"), crs = st_crs(final4))

# Creating actual plots
ses_list <- list(
  race = list(group = "white_no2_group", ylab = "% Non-white", suffix = "race"),
  edu  = list(group = "edu_no2_group",   ylab = "% Low education", suffix = "edu"),
  inc  = list(group = "inc_no2_group",   ylab = "Low income", suffix = "income"),
  emp  = list(group = "emp_no2_group",   ylab = "% Unemployed", suffix = "employment"),
  pov  = list(group = "pov_no2_group",   ylab = "Poverty rate", suffix = "poverty")
)

for (ses_name in names(ses_list)) {
  
  cat("Making map for:", ses_name, "\n")
  
  # Set mapping variables
  group_col <- ses_list[[ses_name]]$group
  ylab_text <- ses_list[[ses_name]]$ylab
  file_tag  <- ses_list[[ses_name]]$suffix
  
  # Legend
  legend_plot <- bivar_legend(pal_vec, dim = 5, ylab = ylab_text)
  
  # Zoom panels
  zoom_panels <- lapply(metro_names, function(nm) {
    bb <- metro_bbox[[nm]]
    make_zoom_panel(
      final4,
      pal_vec = pal_vec,
      title = paste0(panel_letters[nm], ". ", metro_full[nm]),
      xlim = bb$xlim,
      ylim = bb$ylim
    )
  })
  
  top_row    <- wrap_plots(zoom_panels[1:4], nrow = 1)
  bottom_row <- wrap_plots(zoom_panels[5:8], nrow = 1)
  
  # State boundaries
  final4 <- st_as_sf(final4)
  st_crs(final4) <- 4326
  states_sf <- st_as_sf(map("state", plot = FALSE, fill = TRUE))
  states_sf <- st_transform(states_sf, st_crs(final4))
  
  # Main map
  p_map <- ggplot() +
    geom_sf(data = final4, aes(fill = .data[[group_col]]), color = NA) +
    geom_sf(data = states_sf, fill = NA, color = "#FAF9F6", linewidth = 0.15) +
    scale_fill_manual(values = pal_vec, guide = "none", na.value = "grey70") +
    geom_sf(data = rect_sf,
            inherit.aes = FALSE,
            fill = NA,
            color = "black",
            linewidth = 0.6) +
    geom_sf_text(data = label_out_pts,
                 inherit.aes = FALSE,
                 aes(label = label),
                 size = 12, fontface = "bold", color = "black") +
    theme_void() +
    theme(
      plot.margin = margin(2,2,2,2),
      text = element_text(size = 34)
    )
  
  p_main_with_legend <- p_map +
    inset_element(legend_plot,
                  left = 0.75, bottom = 0.10,
                  right = 1.00, top = 0.35,
                  align_to = "panel")
  
  final_panel <- top_row / p_main_with_legend / bottom_row +
    plot_layout(heights = c(1.2, 3.0, 1.2))
  
  # Output file name
  out_png <- paste0(
    "/Users/k_sy_n_imac/Library/CloudStorage/Dropbox/GWU/Lab/FINESST/Research/Results/",
    "disparity_map_multipanel_", file_tag, ".png"
  )
  
  # Save
  ggsave(
    filename = out_png,
    plot     = final_panel,
    width    = 14,
    height   = 12,
    units    = "in",
    dpi      = 300,
    bg       = "white"
  )
}

# Descriptive 
summary_table <- final4 %>%
  sf::st_drop_geometry() %>%
  filter(white_no2_group %in% c("W1_N5", "W5_N1")) %>%
  count(white_no2_group) %>%
  mutate(percent = 100 * n / sum(n))
summary_table <- final4 %>%
  sf::st_drop_geometry() %>%
  filter(white_no2_group %in% c("W1_N5", "W5_N1")) %>%
  summarise(n = n(), percent = 100 * n() / nrow(sf::st_drop_geometry(final4)))
summary_table

summary_table <- final4 %>%
  sf::st_drop_geometry() %>%
  count(edu_no2_group) %>%
  mutate(percent = 100 * n / sum(n))
summary_table <- final4 %>%
  sf::st_drop_geometry() %>%
  filter(edu_no2_group %in% c("W1_N5", "W5_N1")) %>%
  summarise(n = n(), percent = 100 * n() / nrow(sf::st_drop_geometry(final4)))
summary_table

summary_table <- final4 %>%
  sf::st_drop_geometry() %>%
  count(emp_no2_group) %>%
  mutate(percent = 100 * n / sum(n))
summary_table <- final4 %>%
  sf::st_drop_geometry() %>%
  filter(emp_no2_group %in% c("W1_N5", "W5_N1")) %>%
  summarise(n = n(), percent = 100 * n() / nrow(sf::st_drop_geometry(final4)))
summary_table

summary_table <- final4 %>%
  sf::st_drop_geometry() %>%
  count(inc_no2_group) %>%
  mutate(percent = 100 * n / sum(n))
summary_table <- final4 %>%
  sf::st_drop_geometry() %>%
  filter(inc_no2_group %in% c("W1_N5", "W5_N1")) %>%
  summarise(n = n(), percent = 100 * n() / nrow(sf::st_drop_geometry(final4)))
summary_table

summary_table <- final4 %>%
  sf::st_drop_geometry() %>%
  count(pov_no2_group) %>%
  mutate(percent = 100 * n / sum(n))
summary_table <- final4 %>%
  sf::st_drop_geometry() %>%
  filter(pov_no2_group %in% c("W1_N5", "W5_N1")) %>%
  summarise(n = n(), percent = 100 * n() / nrow(sf::st_drop_geometry(final4)))
summary_table





# Figure 5: concentration curve ####################################################################
# Functions for plotting
make_ccurve_plot_bivar_stacked <- function(df,
                                           varname,
                                           group_col,      # e.g., "white_no2_group"
                                           pal_vec,
                                           panel_title,
                                           xlabel,
                                           xbin_width = 1,     # percentile-bin width; try 0.5, 1, 2
                                           strip_height = 12,  # y-units reserved below 0
                                           strip_scale = c("global", "panel")) {
  
  strip_scale <- match.arg(strip_scale)
  if (inherits(df, "sf")) df <- sf::st_drop_geometry(df)
  
  d <- df %>%
    dplyr::select(surface_no2, dplyr::all_of(varname), dplyr::all_of(group_col)) %>%
    dplyr::filter(!is.na(surface_no2), !is.na(.data[[varname]]), !is.na(.data[[group_col]])) 
  ecdf_function <- stats::ecdf(d[[varname]])
  d <- d[order(d[[varname]]), ]
  d$x_pct <- ecdf_function(d[[varname]]) * 100
  d <- d[order(d$x_pct), ]
  d$cumsum_y <- cumsum(d$surface_no2) / sum(d$surface_no2) * 100
  
  auc <- pracma::trapz(d$x_pct, d$cumsum_y)
  CI  <- (5000 - auc) / 5000
  
  # --- Build stacked strip data ---
  breaks <- seq(0, 100, by = xbin_width)
  
  strip_long <- d %>%
    dplyr::mutate(
      x_bin = cut(x_pct, breaks = breaks, include.lowest = TRUE, right = FALSE)
    ) %>%
    dplyr::group_by(x_bin, key = .data[[group_col]]) %>%
    dplyr::summarise(
      no2_sum = sum(surface_no2, na.rm = TRUE),
      x0 = min(x_pct),
      x1 = max(x_pct),
      .groups = "drop"
    ) %>%
    dplyr::filter(!is.na(key))
  
  # total NO2 per bin (for scaling)
  bin_tot <- strip_long %>%
    dplyr::group_by(x_bin) %>%
    dplyr::summarise(bin_no2 = sum(no2_sum), x0 = min(x0), x1 = max(x1), .groups = "drop")
  
  # choose scaling: same scale within this panel ("panel") or (later) across all panels ("global")
  # Here "panel" means max bin total within this df panel determines scaling.
  scale_denom <- if (strip_scale == "panel") max(bin_tot$bin_no2, na.rm = TRUE) else 1
  
  # join totals and compute stacked ymin/ymax in scaled units
  strip_rect <- strip_long %>%
    dplyr::left_join(bin_tot %>% dplyr::select(x_bin, bin_no2, x0_bin = x0, x1_bin = x1),
                     by = "x_bin") %>%
    dplyr::arrange(x_bin, key) %>%
    dplyr::group_by(x_bin) %>%
    dplyr::mutate(
      # stack in original ppb first
      y0_ppb = cumsum(dplyr::lag(no2_sum, default = 0)),
      y1_ppb = y0_ppb + no2_sum
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      # scale ppb stack into [0, strip_height]
      # if global scaling, set scale_denom yourself outside this function (see note below)
      y0_sc = (y0_ppb / scale_denom) * strip_height,
      y1_sc = (y1_ppb / scale_denom) * strip_height,
      # place strip below 0
      ymin = -strip_height + y0_sc,
      ymax = -strip_height + y1_sc,
      xmin = x0_bin,
      xmax = x1_bin
    )
  
  ggplot() +
    # --- stacked strip under the curve ---
    geom_rect(
      data = strip_rect,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = key),
      linewidth = 0
    ) +
    scale_fill_manual(values = pal_vec, guide = "none", na.value = "grey70") +
    
    # --- concentration curve ---
    geom_ribbon(
      data = d,
      aes(x = x_pct, ymin = 0, ymax = cumsum_y),
      fill = "royalblue3", alpha = 0.25
    ) +
    geom_line(
      data = d,
      aes(x = x_pct, y = cumsum_y),
      color = "royalblue3", linewidth = 1.4
    ) +
    geom_abline(intercept = 0, slope = 1, color = "grey40", linewidth = 1.2) +
    geom_ribbon(
      data = d,
      aes(x = x_pct, ymin = 0, ymax = x_pct),
      fill = "grey50", alpha = 0.18
    ) +
    
    scale_x_continuous(limits = c(0, 100), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-strip_height, 100), expand = c(0, 0)) +
    
    labs(
      title = panel_title,
      x = xlabel,
      y = bquote("Cumulative " * NO[2] * " concentrations (%)")
    ) +
    annotate("text", x = 17, y = 90,
             label = paste0("CIX = ", round(CI, 3)),
             size = 7, color = "black") +
    coord_cartesian(clip = "off") +
    theme_minimal() +
    theme(
      plot.title   = element_text(size = 22, face = "bold", hjust = 0.5),
      axis.title.x = element_text(size = 18, margin = margin(t = 10)),
      axis.title.y = element_text(size = 18, margin = margin(r = 10)),
      axis.text    = element_text(size = 18),
      plot.margin  = margin(t = 10, r = 10, b = 30, l = 10)
    )
}

# Actual plotting for five different factors
vars2 <- tibble::tribble(
  ~key,                    ~var,                    ~xlabel,                                                     ~group_col,
  "Race/ethnicity",        "pct_white",             "Census tract percentile\n(0 = least white to 100 = most white)", "white_no2_group",
  "Educational attainment","pct_high_over",         "Census tract percentile\n(0 = lowest educational attainment to 100 = highest)", "edu_no2_group",
  "Employment rate",       "employ_rate",           "Census tract percentile\n(0 = lowest employment rate to 100 = highest)", "emp_no2_group",
  "Mean household income", "mean_household_income", "Census tract percentile\n(0 = lowest mean household income to 100 = highest)", "inc_no2_group",
  "Poverty rate",          "pov_rate_reverse",      "Census tract percentile\n(0 = highest poverty to 100 = lowest)", "pov_no2_group"
)

final4_df <- if (inherits(final4, "sf")) sf::st_drop_geometry(final4) else final4

plots <- lapply(seq_len(nrow(vars2)), function(k) {
  make_ccurve_plot_bivar_stacked(
    df = final4_df,
    varname = vars2$var[k],
    group_col = vars2$group_col[k],
    pal_vec = pal_vec,
    panel_title = vars2$key[k],
    xlabel = vars2$xlabel[k],
    xbin_width = 1,       # try 0.5 if you want more granular x slices
    strip_height = 12,
    strip_scale = "panel" # each panel rescales its strip to fit nicely
  )
})

blank <- ggplot() + theme_void()

pdf("/Users/k_sy_n_imac/Library/CloudStorage/Dropbox/GWU/Lab/FINESST/Research/Results/concentration_curve_bivar_NO2_2024_final.pdf", width = 20, height = 14)
patchwork::wrap_plots(c(plots, list(blank)), ncol = 3)
dev.off()




# Table 2: CI in 2019 and 2024 ####################################################################
list_year <- c(2019, 2020, 2021, 2022, 2023, 2024)
list_cause <- c("All causes", "Cardiovascular diseases", "Respiratory diseases", "Ischemic heart disease", "Tracheal, bronchus, and lung cancer") 

for (j in 1) {
  tmp_ci <- data.frame(year=numeric(), var=character(), ci=numeric())
  
  for (i in c(1:6)) {
    
    no2_ind <- subset(no2_tract_df, year==list_year[i])
    db_ind <- subset(final3, year==list_year[i] & cause_name==list_cause[j])
    final4 <- merge(no2_ind, db_ind, by="GEOID")
    final4 <- merge(final4, xxxx, by="GEOID", all.x=TRUE)
    
    #race/ethnicity
    final_clean_filtered <- final4 %>% dplyr::filter(!is.na(pct_white), !is.na(surface_no2)) 
    ecdf_function <- stats::ecdf(final_clean_filtered[["pct_white"]])
    final_clean_filtered <- final_clean_filtered[order(final_clean_filtered[["pct_white"]]), ]
    final_clean_filtered$x_pct <- ecdf_function(final_clean_filtered[["pct_white"]]) * 100
    final_clean_filtered <- final_clean_filtered[order(final_clean_filtered$x_pct), ]
    final_clean_filtered$cumsum_y <- cumsum(final_clean_filtered$surface_no2) / sum(final_clean_filtered$surface_no2) * 100
    
    auc <- pracma::trapz(final_clean_filtered$x_pct, final_clean_filtered$cumsum_y)
    CI  <- (5000 - auc) / 5000
    
    bb <- data.frame(year = list_year[i], var  = "% white", ci = CI, stringsAsFactors = FALSE)
    tmp_ci <- rbind(tmp_ci, bb)
    
    #education
    final_clean_filtered <- final4 %>% dplyr::filter(!is.na(pct_high_over), !is.na(surface_no2)) 
    ecdf_function <- stats::ecdf(final_clean_filtered[["pct_high_over"]])
    final_clean_filtered <- final_clean_filtered[order(final_clean_filtered[["pct_high_over"]]), ]
    final_clean_filtered$x_pct <- ecdf_function(final_clean_filtered[["pct_high_over"]]) * 100
    final_clean_filtered <- final_clean_filtered[order(final_clean_filtered$x_pct), ]
    final_clean_filtered$cumsum_y <- cumsum(final_clean_filtered$surface_no2) / sum(final_clean_filtered$surface_no2) * 100
    
    auc <- pracma::trapz(final_clean_filtered$x_pct, final_clean_filtered$cumsum_y)
    CI  <- (5000 - auc) / 5000
    
    bb <- data.frame(year = list_year[i], var  = "Educational attainment", ci = CI, stringsAsFactors = FALSE)
    tmp_ci <- rbind(tmp_ci, bb)
    
    #employment
    final_clean_filtered <- final4 %>% dplyr::filter(!is.na(employ_rate), !is.na(surface_no2)) 
    ecdf_function <- stats::ecdf(final_clean_filtered[["employ_rate"]])
    final_clean_filtered <- final_clean_filtered[order(final_clean_filtered[["employ_rate"]]), ]
    final_clean_filtered$x_pct <- ecdf_function(final_clean_filtered[["employ_rate"]]) * 100
    final_clean_filtered <- final_clean_filtered[order(final_clean_filtered$x_pct), ]
    final_clean_filtered$cumsum_y <- cumsum(final_clean_filtered$surface_no2) / sum(final_clean_filtered$surface_no2) * 100
    
    auc <- pracma::trapz(final_clean_filtered$x_pct, final_clean_filtered$cumsum_y)
    CI  <- (5000 - auc) / 5000
    
    bb <- data.frame(year = list_year[i], var  = "Employment rate", ci = CI, stringsAsFactors = FALSE)
    tmp_ci <- rbind(tmp_ci, bb)
    
    #household income
    final_clean_filtered <- final4 %>% dplyr::filter(!is.na(mean_household_income), !is.na(surface_no2)) 
    ecdf_function <- stats::ecdf(final_clean_filtered[['mean_household_income']])
    final_clean_filtered <- final_clean_filtered[order(final_clean_filtered[['mean_household_income']]), ]
    final_clean_filtered$x_pct <- ecdf_function(final_clean_filtered[['mean_household_income']]) * 100
    final_clean_filtered <- final_clean_filtered[order(final_clean_filtered$x_pct), ]
    final_clean_filtered$cumsum_y <- cumsum(final_clean_filtered$surface_no2) / sum(final_clean_filtered$surface_no2) * 100
    
    auc <- pracma::trapz(final_clean_filtered$x_pct, final_clean_filtered$cumsum_y)
    CI  <- (5000 - auc) / 5000
    
    bb <- data.frame(year = list_year[i], var  = "Mean household income", ci = CI, stringsAsFactors = FALSE)
    tmp_ci <- rbind(tmp_ci, bb)
    
    #poverty rate
    final_clean_filtered <- final4 %>% dplyr::filter(!is.na(pov_rate), !is.na(surface_no2)) 
    final_clean_filtered$pov_rate_reverse <- 100-final_clean_filtered$pov_rate
    ecdf_function <- stats::ecdf(final_clean_filtered[['pov_rate_reverse']])
    final_clean_filtered <- final_clean_filtered[order(final_clean_filtered[['pov_rate_reverse']]), ]
    final_clean_filtered$x_pct <- ecdf_function(final_clean_filtered[['pov_rate_reverse']]) * 100
    final_clean_filtered <- final_clean_filtered[order(final_clean_filtered$x_pct), ]
    final_clean_filtered$cumsum_y <- cumsum(final_clean_filtered$surface_no2) / sum(final_clean_filtered$surface_no2) * 100
    
    auc <- pracma::trapz(final_clean_filtered$x_pct, final_clean_filtered$cumsum_y)
    CI  <- (5000 - auc) / 5000
    
    bb <- data.frame(year = list_year[i], var  = "Poverty rate", ci = CI, stringsAsFactors = FALSE)
    tmp_ci <- rbind(tmp_ci, bb)
  }
}
tmp_ci














##############################################################################################
# Figure S7, S8, S9, S10, S11: Time series of CI ###############################################################
##############################################################################################




library(pracma)
library(dplyr)
library(ggplot2)

# Settings
list_year  <- c(2019, 2020, 2021, 2022, 2023, 2024)
list_cause <- c(
  "All causes",
  "Cardiovascular diseases",
  "Chronic respiratory diseases",
  "Ischemic heart disease",
  "Tracheal, bronchus, and lung cancer"
)
j <- 1

region_labels <- c(
  "Region 1\n(New England)",
  "Region 2\n(Northeast/Caribbean)",
  "Region 3\n(Mid-Atlantic)",
  "Region 4\n(Southeast)",
  "Region 5\n(Great Lakes)",
  "Region 6\n(South Central)",
  "Region 7\n(Heartland)",
  "Region 8\n(Mountains/Plains)",
  "Region 9\n(Pacific Southwest)",
  "Region 10\n(Northwest)"
)

recode_region <- function(x) {
  recode(x,
         "Region 1"  = "Region 1\n(New England)",
         "Region 2"  = "Region 2\n(Northeast/Caribbean)",
         "Region 3"  = "Region 3\n(Mid-Atlantic)",
         "Region 4"  = "Region 4\n(Southeast)",
         "Region 5"  = "Region 5\n(Great Lakes)",
         "Region 6"  = "Region 6\n(South Central)",
         "Region 7"  = "Region 7\n(Heartland)",
         "Region 8"  = "Region 8\n(Mountains/Plains)",
         "Region 9"  = "Region 9\n(Pacific Southwest)",
         "Region 10" = "Region 10\n(Northwest)"
  )
}

# -------------------------------------------------------------------
# Factor definitions
# Each entry: list(var = column name, label = display label, reverse = T/F)
# reverse = TRUE for poverty rate (higher poverty = worse, so flip)
# -------------------------------------------------------------------
factor_list <- list(
  list(var = "pct_white",            label = "% white",                reverse = FALSE),
  list(var = "pct_high_over",        label = "Educational attainment",  reverse = FALSE),
  list(var = "employ_rate",          label = "Employment rate",         reverse = FALSE),
  list(var = "mean_household_income",label = "Mean household income",   reverse = FALSE),
  list(var = "pov_rate",             label = "Poverty rate",            reverse = TRUE)
)

# -------------------------------------------------------------------
# Generalized CI function
# -------------------------------------------------------------------
calc_ci <- function(df, var, reverse = FALSE) {
  
  df <- df %>% filter(!is.na(.data[[var]]), !is.na(surface_no2))
  if (nrow(df) < 2) return(NA)
  
  # Reverse the variable if needed (e.g., poverty rate → higher = better)
  if (reverse) {
    df[[var]] <- 100 - df[[var]]
  }
  
  ecdf_fn   <- ecdf(df[[var]])
  df        <- df %>% arrange(.data[[var]])
  df$x_pct  <- ecdf_fn(df[[var]]) * 100
  df        <- df %>% arrange(x_pct)
  df$cumsum <- cumsum(df$surface_no2) / sum(df$surface_no2) * 100
  
  auc <- trapz(df$x_pct, df$cumsum)
  CI  <- (5000 - auc) / 5000
  return(CI)
}

# -------------------------------------------------------------------
# Output containers
# -------------------------------------------------------------------
ci_group_timeseries  <- data.frame()
ci_region_timeseries <- data.frame()
ci_metro_timeseries  <- data.frame()

# -------------------------------------------------------------------
# Main loop: year × factor
# -------------------------------------------------------------------
for (yr in list_year) {
  cat("Processing year:", yr, "\n")
  
  no2_ind  <- subset(no2_tract_df, year == yr)
  db_ind   <- subset(final3, year == yr & cause_name == list_cause[j])
  
  final4_yr <- merge(no2_ind, db_ind, by = "GEOID")
  final4_yr <- merge(final4_yr, xxxx, by = "GEOID", all.x = TRUE)
  
  if ("year.x" %in% names(final4_yr)) {
    final4_yr$year <- final4_yr$year.x
  }
  
  final4_yr <- final4_yr %>%
    mutate(EPA_region = recode_region(EPA_region))
  
  for (fac in factor_list) {
    
    var_name  <- fac$var
    var_label <- fac$label
    rev_flag  <- fac$reverse
    
    cat("  Factor:", var_label, "\n")
    
    # ================================================================
    # 1. Urbanicity × EPA region
    # ================================================================
    combos <- final4_yr %>% distinct(EPA_region, metro4)
    
    for (ii in 1:nrow(combos)) {
      reg_i   <- combos$EPA_region[ii]
      metro_i <- combos$metro4[ii]
      
      sub_df <- final4_yr %>%
        filter(EPA_region == reg_i, metro4 == metro_i)
      
      CI_val <- calc_ci(sub_df, var_name, rev_flag)
      
      ci_group_timeseries <- rbind(
        ci_group_timeseries,
        data.frame(
          year       = yr,
          factor     = var_label,
          EPA_region = reg_i,
          metro4     = metro_i,
          CI_no2     = CI_val,
          n_tract    = nrow(sub_df)
        )
      )
    }
    
    # ================================================================
    # 2. EPA region only
    # ================================================================
    for (reg_i in unique(final4_yr$EPA_region)) {
      
      sub_df <- final4_yr %>% filter(EPA_region == reg_i)
      CI_val <- calc_ci(sub_df, var_name, rev_flag)
      
      ci_region_timeseries <- rbind(
        ci_region_timeseries,
        data.frame(
          year       = yr,
          factor     = var_label,
          EPA_region = reg_i,
          CI_no2     = CI_val,
          n_tract    = nrow(sub_df)
        )
      )
    }
    
    # ================================================================
    # 3. Urbanicity only
    # ================================================================
    for (metro_i in unique(final4_yr$metro4)) {
      
      sub_df <- final4_yr %>% filter(metro4 == metro_i)
      CI_val <- calc_ci(sub_df, var_name, rev_flag)
      
      ci_metro_timeseries <- rbind(
        ci_metro_timeseries,
        data.frame(
          year    = yr,
          factor  = var_label,
          metro4  = metro_i,
          CI_no2  = CI_val,
          n_tract = nrow(sub_df)
        )
      )
    }
  }
  
  gc()
}

# -------------------------------------------------------------------
# Factor ordering (consistent across all three data frames)
# -------------------------------------------------------------------
factor_levels <- c(
  "% white",
  "Educational attainment",
  "Employment rate",
  "Mean household income",
  "Poverty rate"
)

metro_levels <- c("Large MSA", "Small MSA", "MicroSA", "Non-urban")

ci_group_timeseries <- ci_group_timeseries %>%
  mutate(
    factor     = factor(factor, levels = factor_levels),
    metro4     = factor(metro4, levels = metro_levels),
    EPA_region = factor(EPA_region, levels = region_labels)
  )

ci_region_timeseries <- ci_region_timeseries %>%
  mutate(
    factor     = factor(factor, levels = factor_levels),
    EPA_region = factor(EPA_region, levels = region_labels)
  )

ci_metro_timeseries <- ci_metro_timeseries %>%
  mutate(
    factor = factor(factor, levels = factor_levels),
    metro4 = factor(metro4, levels = metro_levels)
  )

# -------------------------------------------------------------------
# Color palettes
# -------------------------------------------------------------------
metro_colors <- c(
  "Large MSA" = "#2E004F",
  "Small MSA" = "#54278F",
  "MicroSA"   = "#756BB1",
  "Non-urban" = "#CBC9E2"
)

region_colors <- c(
  "Region 1\n(New England)"         = "#1B2F70",
  "Region 2\n(Northeast/Caribbean)" = "#2C4FD7",
  "Region 3\n(Mid-Atlantic)"        = "#5A8AC6",
  "Region 4\n(Southeast)"           = "#8BC8C8",
  "Region 5\n(Great Lakes)"         = "#A1D99B",
  "Region 6\n(South Central)"       = "#D9C89B",
  "Region 7\n(Heartland)"           = "#DDAA66",
  "Region 8\n(Mountains/Plains)"    = "#F2A6A0",
  "Region 9\n(Pacific Southwest)"   = "#7A5AC9",
  "Region 10\n(Northwest)"          = "#4A0066"
)

# -------------------------------------------------------------------
# National reference lines — one per factor (from tmp_ci)
# -------------------------------------------------------------------
nation_data_list <- list(
  "% white"                = tmp_ci %>% filter(var == "% white")                %>% dplyr::select(year, CI_mean = ci),
  "Educational attainment" = tmp_ci %>% filter(var == "Educational attainment") %>% dplyr::select(year, CI_mean = ci),
  "Employment rate"        = tmp_ci %>% filter(var == "Employment rate")        %>% dplyr::select(year, CI_mean = ci),
  "Mean household income"  = tmp_ci %>% filter(var == "Mean household income")  %>% dplyr::select(year, CI_mean = ci),
  "Poverty rate"           = tmp_ci %>% filter(var == "Poverty rate")           %>% dplyr::select(year, CI_mean = ci)
)

# -------------------------------------------------------------------
# Plot helper: shared theme
# -------------------------------------------------------------------
base_theme <- function(base_size = 16) {
  theme_minimal() +
    theme(
      strip.text    = element_text(size = base_size, face = "bold"),
      axis.title.x  = element_text(size = base_size, margin = margin(t = 15)),
      axis.title.y  = element_text(size = base_size, margin = margin(r = 15)),
      axis.text     = element_text(size = base_size - 4),
      legend.title  = element_text(size = base_size),
      legend.text   = element_text(size = base_size),
      legend.position = "bottom"
    )
}

# -------------------------------------------------------------------
# 1. Urbanicity × EPA region — one PDF per factor (faceted by region)
# -------------------------------------------------------------------
for (fac_label in factor_levels) {
  
  safe_name <- gsub("[^A-Za-z0-9]", "_", fac_label)
  
  pdf(
    file  = paste0(
      "/Users/k_sy_n_imac/Library/CloudStorage/Dropbox/GWU/Lab/FINESST/Research/Results/test/",
      "time_series_CI_final_", safe_name, ".pdf"
    ),
    width = 20, height = 7
  )
  par(family = "Lato")
  
  plot_data <- ci_group_timeseries %>% filter(factor == fac_label)
  
  p <- ggplot(
    plot_data,
    aes(x = year, y = CI_no2, color = metro4, group = metro4)
  ) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 1) +
    geom_line(
      data = nation_data_list[[fac_label]], 
      aes(x = year, y = CI_mean, color = "National"),
      inherit.aes = FALSE, linewidth = 1.1
    ) +
    facet_wrap(~ EPA_region, ncol = 5) +
    scale_color_manual(
      values = c(metro_colors, "National" = "#B11226"),
      name   = "Urbanicity"
    ) +
    labs(
      x     = "Year",
      y     = "Concentration index"
    ) +
    theme_minimal() +
    theme(
      strip.text = element_text(size = 16, face = "bold"),
      axis.title.x = element_text(size = 16, margin = margin(t = 15)),
      axis.title.y = element_text(size = 16,margin = margin(r = 15)),    axis.title = element_text(size = 16),
      axis.text = element_text(size = 12),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 16),
      legend.position = "bottom"
    )
  print(p)
  dev.off()
}

# -------------------------------------------------------------------
# 2. EPA region only — one PDF per factor
# -------------------------------------------------------------------
for (fac_label in factor_levels) {
  
  safe_name <- gsub("[^A-Za-z0-9]", "_", fac_label)
  
  pdf(
    file  = paste0(
      "/Users/k_sy_n_imac/Library/CloudStorage/Dropbox/GWU/Lab/FINESST/Research/Results/test/",
      "time_series_CI_by_region_", safe_name, ".pdf"
    ),
    width = 15, height = 7
  )
  par(family = "Lato")
  
  plot_data <- ci_region_timeseries %>% filter(factor == fac_label)
  
  p <- ggplot(
    plot_data,
    aes(x = year, y = CI_no2, color = EPA_region, group = EPA_region)
  ) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2) +
    geom_line(
      data = nation_data_list[[fac_label]], 
      aes(x = year, y = CI_mean, color = "National"),
      inherit.aes = FALSE, linewidth = 1.1
    ) +
    scale_color_manual(
      values = c(region_colors, "National" = "#B11226"),
      breaks = c(names(region_colors), "National"),
      labels = c(
        "Region 1 (New England)", "Region 2 (Northeast/Caribbean)",
        "Region 3 (Mid-Atlantic)", "Region 4 (Southeast)",
        "Region 5 (Great Lakes)", "Region 6 (South Central)",
        "Region 7 (Heartland)", "Region 8 (Mountains/Plains)",
        "Region 9 (Pacific Southwest)", "Region 10 (Northwest)",
        "National"
      ),
      name = "EPA region"
    ) +
    labs(
      x     = "Year",
      y     = "Concentration index"
    ) +
    theme_minimal() + 
    theme(
      strip.text = element_text(size = 18, face = "bold"),
      axis.title.x = element_text(size = 18, margin = margin(t = 15)),
      axis.title.y = element_text(size = 18,margin = margin(r = 15)),  
      axis.text = element_text(size = 17),
      legend.title = element_text(size = 15),
      legend.text = element_text(size = 15),
      legend.position = "bottom")
  
  print(p)
  dev.off()
}

# -------------------------------------------------------------------
# 3. Urbanicity only — one PDF per factor
# -------------------------------------------------------------------
for (fac_label in factor_levels) {
  
  safe_name <- gsub("[^A-Za-z0-9]", "_", fac_label)
  
  pdf(
    file  = paste0(
      "/Users/k_sy_n_imac/Library/CloudStorage/Dropbox/GWU/Lab/FINESST/Research/Results/test/",
      "time_series_CI_by_urbanicity_", safe_name, ".pdf"
    ),
    width = 10, height = 7
  )
  par(family = "Lato")
  
  plot_data <- ci_metro_timeseries %>% filter(factor == fac_label)
  
  p <- ggplot(
    plot_data,
    aes(x = year, y = CI_no2, color = metro4, group = metro4)
  ) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 1.2) +
    geom_line(
      data = nation_data_list[[fac_label]], 
      aes(x = year, y = CI_mean, color = "National"),
      inherit.aes = FALSE, linewidth = 1.1
    ) +
    scale_color_manual(
      values = c(metro_colors, "National" = "#B11226"),
      name   = "Urbanicity"
    ) +
    labs(
      x     = "Year",
      y     = "Concentration index"
    ) +
    theme_minimal() + 
    theme(
      strip.text = element_text(size = 18, face = "bold"),
      axis.title.x = element_text(size = 18, margin = margin(t = 15)),
      axis.title.y = element_text(size = 18,margin = margin(r = 15)),  
      axis.text = element_text(size = 17),
      legend.title = element_text(size = 15),
      legend.text = element_text(size = 15),
      legend.position = "bottom")
  
  print(p)
  dev.off()
}





# -------------------------------------------------------------------
# Maximum temporal swing (max - min CI across years) within each group
# -------------------------------------------------------------------

# 1. Urbanicity × EPA region
swing_group <- ci_group_timeseries %>%
  group_by(factor, EPA_region, metro4) %>%
  summarise(
    CI_min    = min(CI_no2,  na.rm = TRUE),
    CI_max    = max(CI_no2,  na.rm = TRUE),
    year_min  = year[which.min(CI_no2)],
    year_max  = year[which.max(CI_no2)],
    swing     = CI_max - CI_min,
    .groups   = "drop"
  ) %>%
  arrange(desc(swing))

# 2. EPA region only
swing_region <- ci_region_timeseries %>%
  group_by(factor, EPA_region) %>%
  summarise(
    CI_min   = min(CI_no2,  na.rm = TRUE),
    CI_max   = max(CI_no2,  na.rm = TRUE),
    year_min = year[which.min(CI_no2)],
    year_max = year[which.max(CI_no2)],
    swing    = CI_max - CI_min,
    .groups  = "drop"
  ) %>%
  arrange(desc(swing))

# 3. Urbanicity only
swing_metro <- ci_metro_timeseries %>%
  group_by(factor, metro4) %>%
  summarise(
    CI_min   = min(CI_no2,  na.rm = TRUE),
    CI_max   = max(CI_no2,  na.rm = TRUE),
    year_min = year[which.min(CI_no2)],
    year_max = year[which.max(CI_no2)],
    swing    = CI_max - CI_min,
    .groups  = "drop"
  ) %>%
  arrange(desc(swing))

# -------------------------------------------------------------------
# Overall maximum across all three stratifications
# -------------------------------------------------------------------
all_swings <- bind_rows(
  swing_group  %>% mutate(stratification = "Urbanicity x EPA region",
                          group = paste(factor, EPA_region, metro4, sep = " | ")),
  swing_region %>% mutate(stratification = "EPA region only",
                          group = paste(factor, EPA_region, sep = " | ")),
  swing_metro  %>% mutate(stratification = "Urbanicity only",
                          group = paste(factor, metro4, sep = " | "))
) %>%
  dplyr::select(stratification, group, factor, CI_min, year_min,
                CI_max, year_max, swing) %>%
  arrange(desc(swing))

# -------------------------------------------------------------------
# Print top results
# -------------------------------------------------------------------
cat("=== Top 10 largest temporal swings (across all stratifications) ===\n")
print(head(all_swings, 10))

cat("\n=== OVERALL MAXIMUM ===\n")
print(all_swings[1, ])

cat("\n=== Top 5: Urbanicity x EPA region ===\n")
print(head(swing_group, 5))

cat("\n=== Top 5: EPA region only ===\n")
print(head(swing_region, 5))

cat("\n=== Top 5: Urbanicity only ===\n")
print(head(swing_metro, 5))












################################################################################
# Figure 6 ################################################################################
################################################################################
var_names_2024 <- c("MSA_title", "metro4", "EPA_region",
                    "CI_no2_2024", "abs_no2_2024", "ratio_no2_2024",
                    "n_tract_total_2024", "n_tract_valid_2024")
est <- data.frame(matrix(ncol = length(var_names_2024), nrow = 0))
colnames(est) <- var_names_2024

var_names_2019 <- c("MSA_title", "metro4", "EPA_region",
                    "CI_no2_2019", "abs_no2_2019", "ratio_no2_2019",
                    "n_tract_total_2019", "n_tract_valid_2019")
est1 <- data.frame(matrix(ncol = length(var_names_2019), nrow = 0))
colnames(est1) <- var_names_2019

# Calculating CI in 2024
no2_2024 <- subset(no2_tract_df, year == 2024)
db_2024  <- subset(final3, year == 2024 & cause_name == list_cause[j])

final4_2024 <- merge(no2_2024, db_2024, by = "GEOID")
final4_2024 <- merge(final4_2024, xxxx, by = "GEOID", all.x = TRUE)

final4_2024 <- final4_2024 %>%
  mutate(
    MSA_title = if_else(
      metro4 == "Non-urban",
      paste0("Non_urban_", state_name.x),
      MSA_title
    )
  )
group_df <- final4_2024 %>% group_by(MSA_title, metro4, EPA_region) %>%
  summarize(n_2024 = n(), .groups = "drop")
msa_regions <- group_df %>% count(MSA_title) %>% filter(n > 1)
majority_regions <- final4_2024 %>%
  group_by(MSA_title, EPA_region) %>%
  summarise(total_pop = sum(pop_over20, na.rm = TRUE), .groups = "drop") %>%
  group_by(MSA_title) %>%
  slice_max(total_pop, n = 1) %>%  # ← Gets region with HIGHEST population
  ungroup() %>%
  dplyr::select(MSA_title, EPA_region_new = EPA_region)
final4_2024 <- final4_2024 %>%
  left_join(majority_regions, by = "MSA_title") %>%
  mutate(EPA_region = if_else(!is.na(EPA_region_new), 
                              EPA_region_new, 
                              EPA_region))

# Use all three variables as unique MSA key
MSA_list_2024 <- final4_2024 %>%
  dplyr::select(MSA_title, EPA_region, metro4) %>%
  distinct() %>%
  filter(!is.na(MSA_title), !is.na(EPA_region), !is.na(metro4))

for (ii in seq_len(nrow(MSA_list_2024))) {
  
  msa_i   <- MSA_list_2024$MSA_title[ii]
  epa_i   <- MSA_list_2024$EPA_region[ii]
  metro_i <- MSA_list_2024$metro4[ii]
  
  final4_metro <- final4_2024 %>%
    filter(MSA_title == msa_i, EPA_region == epa_i, metro4 == metro_i)
  
  n_total <- nrow(final4_metro)
  
  final_clean_filtered <- final4_metro %>%
    filter(!is.na(pct_white), !is.na(surface_no2))
  
  n_valid <- nrow(final_clean_filtered)
  
  if (n_valid < 5) next
  
  ecdf_function <- stats::ecdf(final_clean_filtered[["pct_white"]])
  final_clean_filtered <- final_clean_filtered[order(final_clean_filtered[["pct_white"]]), ]
  final_clean_filtered$x_pct <- ecdf_function(final_clean_filtered[["pct_white"]]) * 100
  final_clean_filtered <- final_clean_filtered[order(final_clean_filtered$x_pct), ]
  final_clean_filtered$cumsum_y <- cumsum(final_clean_filtered$surface_no2) /
    sum(final_clean_filtered$surface_no2) * 100
  
  auc <- pracma::trapz(final_clean_filtered$x_pct, final_clean_filtered$cumsum_y)
  CI_no2 <- (5000 - auc) / 5000
  
  final_dec <- final_clean_filtered %>%
    mutate(white_decile = ntile(pct_white, 10)) %>%
    group_by(white_decile) %>%
    summarise(mean_no2 = mean(surface_no2, na.rm = TRUE), .groups = "drop") %>%
    right_join(data.frame(white_decile = 1:10), by = "white_decile") %>%
    arrange(white_decile)
  
  abs_no2   <- unname(-final_dec$mean_no2[1] + final_dec$mean_no2[10])
  ratio_no2 <- unname(final_dec$mean_no2[1] / final_dec$mean_no2[10])
  
  est <- rbind(est, data.frame(
    MSA_title          = msa_i,
    metro4             = metro_i,
    EPA_region         = epa_i,
    CI_no2_2024        = CI_no2,
    abs_no2_2024       = abs_no2,
    ratio_no2_2024     = ratio_no2,
    n_tract_total_2024 = n_total,
    n_tract_valid_2024 = n_valid,
    stringsAsFactors   = FALSE
  ))
}

rm(no2_2024, db_2024, final4_2024, final4_metro, final_clean_filtered, ecdf_function, final_dec)
gc()

# Calculating CI in 2019
no2_2019 <- subset(no2_tract_df, year == 2019)
db_2019  <- subset(final3, year == 2019 & cause_name == list_cause[j])

final4_2019 <- merge(no2_2019, db_2019, by = "GEOID")
final4_2019 <- merge(final4_2019, xxxx, by = "GEOID", all.x = TRUE)

final4_2019 <- final4_2019 %>%
  mutate(
    MSA_title = if_else(
      metro4 == "Non-urban",
      paste0("Non_urban_", state_name.x),
      MSA_title
    )
  )
group_df <- final4_2019 %>% group_by(MSA_title, metro4, EPA_region) %>%
  summarize(n_2024 = n(), .groups = "drop")
msa_regions <- group_df %>% count(MSA_title) %>% filter(n > 1)
majority_regions <- final4_2019 %>%
  group_by(MSA_title, EPA_region) %>%
  summarise(total_pop = sum(pop_over20, na.rm = TRUE), .groups = "drop") %>%
  group_by(MSA_title) %>%
  slice_max(total_pop, n = 1) %>%  # ← Gets region with HIGHEST population
  ungroup() %>%
  dplyr::select(MSA_title, EPA_region_new = EPA_region)
final4_2019 <- final4_2019 %>%
  left_join(majority_regions, by = "MSA_title") %>%
  mutate(EPA_region = if_else(!is.na(EPA_region_new), 
                              EPA_region_new, 
                              EPA_region))

# Use all three variables as unique MSA key
MSA_list_2019 <- final4_2019 %>%
  dplyr::select(MSA_title, EPA_region, metro4) %>%
  distinct() %>%
  filter(!is.na(MSA_title), !is.na(EPA_region), !is.na(metro4))

for (ii in seq_len(nrow(MSA_list_2019))) {
  
  msa_i   <- MSA_list_2019$MSA_title[ii]
  epa_i   <- MSA_list_2019$EPA_region[ii]
  metro_i <- MSA_list_2019$metro4[ii]
  
  final4_metro <- final4_2019 %>%
    filter(MSA_title == msa_i, EPA_region == epa_i, metro4 == metro_i)
  
  n_total <- nrow(final4_metro)
  
  final_clean_filtered <- final4_metro %>%
    filter(!is.na(pct_white), !is.na(surface_no2))
  
  n_valid <- nrow(final_clean_filtered)
  
  if (n_valid < 5) next
  
  ecdf_function <- stats::ecdf(final_clean_filtered[["pct_white"]])
  final_clean_filtered <- final_clean_filtered[order(final_clean_filtered[["pct_white"]]), ]
  final_clean_filtered$x_pct <- ecdf_function(final_clean_filtered[["pct_white"]]) * 100
  final_clean_filtered <- final_clean_filtered[order(final_clean_filtered$x_pct), ]
  final_clean_filtered$cumsum_y <- cumsum(final_clean_filtered$surface_no2) /
    sum(final_clean_filtered$surface_no2) * 100
  
  auc <- pracma::trapz(final_clean_filtered$x_pct, final_clean_filtered$cumsum_y)
  CI_no2 <- (5000 - auc) / 5000
  
  final_dec <- final_clean_filtered %>%
    mutate(white_decile = ntile(pct_white, 10)) %>%
    group_by(white_decile) %>%
    summarise(mean_no2 = mean(surface_no2, na.rm = TRUE), .groups = "drop") %>%
    right_join(data.frame(white_decile = 1:10), by = "white_decile") %>%
    arrange(white_decile)
  
  abs_no2   <- unname(-final_dec$mean_no2[1] + final_dec$mean_no2[10])
  ratio_no2 <- unname(final_dec$mean_no2[1] / final_dec$mean_no2[10])
  
  est1 <- rbind(est1, data.frame(
    MSA_title          = msa_i,
    metro4             = metro_i,
    EPA_region         = epa_i,
    CI_no2_2019        = CI_no2,
    abs_no2_2019       = abs_no2,
    ratio_no2_2019     = ratio_no2,
    n_tract_total_2019 = n_total,
    n_tract_valid_2019 = n_valid,
    stringsAsFactors   = FALSE
  ))
}

rm(no2_2019, db_2019, final4_2019, final4_metro, final_clean_filtered, ecdf_function, final_dec)
gc()

# Merge on all three keys
final_est <- full_join(est, est1, by = c("MSA_title", "EPA_region", "metro4"))
final_est$CI_dif <- final_est$CI_no2_2024 - final_est$CI_no2_2019

n_msa_total_2024 <- nrow(MSA_list_2024)
n_msa_valid_2024 <- nrow(est)
n_msa_total_2019 <- nrow(MSA_list_2019)
n_msa_valid_2019 <- nrow(est1)


################################################################################
# Prepare raw (MSA-level) data
################################################################################

region_labels <- c(
  " Region 1\n (New England)",
  " Region 2\n (Northeast/Caribbean)",
  " Region 3\n (Mid-Atlantic)",
  " Region 4\n (Southeast)",
  " Region 5\n (Great Lakes)",
  " Region 6\n (South Central)",
  " Region 7\n (Heartland)",
  " Region 8\n (Mountains/Plains)",
  " Region 9\n (Pacific Southwest)",
  " Region 10\n (Northwest)")

recode_region <- function(x) {
  recode(x,
         "Region 1"  = " Region 1\n (New England)",
         "Region 2"  = " Region 2\n (Northeast/Caribbean)",
         "Region 3"  = " Region 3\n (Mid-Atlantic)",
         "Region 4"  = " Region 4\n (Southeast)",
         "Region 5"  = " Region 5\n (Great Lakes)",
         "Region 6"  = " Region 6\n (South Central)",
         "Region 7"  = " Region 7\n (Heartland)",
         "Region 8"  = " Region 8\n (Mountains/Plains)",
         "Region 9"  = " Region 9\n (Pacific Southwest)",
         "Region 10" = " Region 10\n (Northwest)")
}

region_colors <- c(
  " Region 1\n (New England)"         = "#1B2F70",
  " Region 2\n (Northeast/Caribbean)" = "#2C4FD7",
  " Region 3\n (Mid-Atlantic)"        = "#5A8AC6",
  " Region 4\n (Southeast)"           = "#8BC8C8",
  " Region 5\n (Great Lakes)"         = "#A1D99B",
  " Region 6\n (South Central)"       = "#D9C89B",
  " Region 7\n (Heartland)"           = "#DDAA66",
  " Region 8\n (Mountains/Plains)"    = "#F2A6A0",
  " Region 9\n (Pacific Southwest)"   = "#7A5AC9",
  " Region 10\n (Northwest)"          = "#4A0066")

plot_df <- final_est %>%
  mutate(
    metro4     = factor(metro4, levels = c("Large MSA", "Small MSA",
                                           "MicroSA", "Non-urban")),
    EPA_region = factor(recode_region(EPA_region), levels = region_labels)
  )


################################################################################
# Plot: raw scatter
################################################################################

pdf(
  file = "/Users/k_sy_n_imac/Library/CloudStorage/Dropbox/GWU/Lab/FINESST/Research/Results/Quadrant_NO2disparities_raw.pdf",
  width = 13,
  height = 10
)
par(family = "Lato")

ggplot(plot_df, aes(x = CI_dif, y = CI_no2_2024,
                    color = EPA_region, shape = metro4)) +
  
  # Reference lines
  geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.8) +
  
  geom_vline(xintercept =  0.03, linetype = "dashed", color = "gray40", linewidth = 0.5) +
  geom_vline(xintercept = -0.03, linetype = "dashed", color = "gray40", linewidth = 0.5) +
  geom_hline(yintercept =  0.13, linetype = "dashed", color = "gray40", linewidth = 0.5) +
  geom_hline(yintercept = -0.13, linetype = "dashed", color = "gray40", linewidth = 0.5) +
  
  annotate("text", x =  0.03,  y = -0.02, label =  "0.03", size = 5, color = "black") +
  annotate("text", x = -0.03,  y = -0.02, label = "-0.03", size = 5, color = "black") +
  annotate("text", x =  0.005, y =  0.13, label =  "0.13", size = 5, color = "black") +
  annotate("text", x =  0.005, y = -0.13, label = "-0.13", size = 5, color = "black") +
  
  # Raw points (one per unique MSA x EPA region x urbanicity)
  geom_point(size = 2.5, alpha = 0.8) +
  
  scale_color_manual(values = region_colors) +
  scale_shape_manual(values = c(
    "Large MSA" = 16,
    "Small MSA" = 17,
    "MicroSA"   = 15,
    "Non-urban" = 8
  )) +
  
  labs(
    x      = bquote(bold("Change in concentration index (race/ethnicity) between 2019 and 2024")),
    y      = bquote(bold("Concentration index (race/ethnicity) in 2024")),
    shape = "Urbanicity",
    color = "EPA Region"
  ) +
  guides(
    shape = guide_legend(
      order = 1,
      nrow = 4,
      override.aes = list(size = 3)
    ),
    color = guide_legend(
      order = 2,
      nrow = 10,
      theme = theme(
        legend.text = element_text(
          size = 13,
          margin = margin(b = 8)
        )
      )
    )
  ) +
  theme(
    plot.title    = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 13),
    axis.text     = element_text(size = 13),
    axis.title.x  = element_text(size = 15, face = "bold", margin = margin(t = 15, b = 15)),
    axis.title.y  = element_text(size = 15, face = "bold", margin = margin(r = 15, l = 15)),
    legend.title  = element_text(size = 13, face = "bold"),
    legend.text   = element_text(size = 13),
    panel.grid.minor  = element_blank(),
    panel.grid.major  = element_line(color = "gray95"),
    panel.background  = element_rect(fill = "white", color = NA),
    plot.background   = element_rect(fill = "white", color = NA)
  ) +
  
  coord_cartesian(
    xlim = c(-0.05, 0.05),
    ylim = c(-0.22,  0.22)
  )
dev.off()


################################################################################
# Summary by EPA REGION ONLY (CI)
################################################################################

group_summary_ci_region <- final_est %>%
  filter(!is.na(CI_no2_2024), !is.na(CI_dif)) %>%
  mutate(EPA_region = factor(EPA_region,
                             levels = c("Region 1", "Region 2", "Region 3", "Region 4", "Region 5",
                                        "Region 6", "Region 7", "Region 8", "Region 9", "Region 10"))) %>%
  group_by(EPA_region) %>%
  summarize(
    n_msa        = n(),
    ci_2024      = mean(CI_no2_2024, na.rm = TRUE),
    var_ci_2024  = var(CI_no2_2024, na.rm = TRUE),
    ci_dif       = mean(CI_dif, na.rm = TRUE),
    var_ci_dif   = var(CI_dif, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    se_ci_2024 = sqrt(var_ci_2024 / n_msa),
    se_ci_dif  = sqrt(var_ci_dif  / n_msa),
    x0 = ci_dif,
    y0 = ci_2024,
    rx = 1.96 * se_ci_dif,
    ry = 1.96 * se_ci_2024
  )

################################################################################
# National summary (CI)
################################################################################

nation_ci <- final_est %>%
  filter(!is.na(CI_no2_2024), !is.na(CI_dif)) %>%
  summarize(
    n_msa       = n(),
    ci_2024     = mean(CI_no2_2024, na.rm = TRUE),
    var_ci_2024 = var(CI_no2_2024, na.rm = TRUE),
    ci_dif      = mean(CI_dif, na.rm = TRUE),
    var_ci_dif  = var(CI_dif, na.rm = TRUE)
  ) %>%
  mutate(
    se_ci_2024 = sqrt(var_ci_2024 / n_msa),
    se_ci_dif  = sqrt(var_ci_dif  / n_msa),
    x0 = ci_dif,
    y0 = ci_2024,
    rx = 1.96 * se_ci_dif,
    ry = 1.96 * se_ci_2024
  )

################################################################################
# Diamond shapes for EPA region (CI)
################################################################################

diamond_ci_region_x <- bind_rows(lapply(seq_len(nrow(group_summary_ci_region)), function(i) {
  tibble(
    EPA_region = group_summary_ci_region$EPA_region[i],
    diamond_id = paste0(group_summary_ci_region$EPA_region[i], "_x"),
    x = c(group_summary_ci_region$x0[i] - group_summary_ci_region$rx[i],
          group_summary_ci_region$x0[i],
          group_summary_ci_region$x0[i] + group_summary_ci_region$rx[i],
          group_summary_ci_region$x0[i]),
    y = c(group_summary_ci_region$y0[i],
          group_summary_ci_region$y0[i] + 0.0005,
          group_summary_ci_region$y0[i],
          group_summary_ci_region$y0[i] - 0.0005)
  )
}))

diamond_ci_region_y <- bind_rows(lapply(seq_len(nrow(group_summary_ci_region)), function(i) {
  tibble(
    EPA_region = group_summary_ci_region$EPA_region[i],
    diamond_id = paste0(group_summary_ci_region$EPA_region[i], "_y"),
    x = c(group_summary_ci_region$x0[i],
          group_summary_ci_region$x0[i] + 0.0002,
          group_summary_ci_region$x0[i],
          group_summary_ci_region$x0[i] - 0.0002),
    y = c(group_summary_ci_region$y0[i] - group_summary_ci_region$ry[i],
          group_summary_ci_region$y0[i],
          group_summary_ci_region$y0[i] + group_summary_ci_region$ry[i],
          group_summary_ci_region$y0[i])
  )
}))

diamond_ci_region <- bind_rows(diamond_ci_region_x, diamond_ci_region_y)

################################################################################
# Plot 1: EPA REGION STRATIFICATION (CI)
################################################################################

region_colors <- c(
  "Region 1"  = "#1B2F70",
  "Region 2"  = "#2C4FD7",
  "Region 3"  = "#5A8AC6",
  "Region 4"  = "#8BC8C8",
  "Region 5"  = "#A1D99B",
  "Region 6"  = "#D9C89B",
  "Region 7"  = "#DDAA66",
  "Region 8"  = "#F2A6A0",
  "Region 9"  = "#7A5AC9",
  "Region 10" = "#4A0066"
)

pdf(
  file = "/Users/k_sy_n_imac/Library/CloudStorage/Dropbox/GWU/Lab/FINESST/Research/Results/Quadrant_CI_EPA_Region.pdf",
  width = 9, height = 9)

ggplot() +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.5) +
  geom_hline(yintercept = 0.13, linetype = "dashed", color = "gray40", linewidth = 0.5) +
  geom_hline(yintercept = -0.13, linetype = "dashed", color = "gray40", linewidth = 0.5) +
  
  geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.5) +
  geom_vline(xintercept = 0.03, linetype = "dashed", color = "gray40", linewidth = 0.5) +
  geom_vline(xintercept = -0.03, linetype = "dashed", color = "gray40", linewidth = 0.5) +
  
  geom_polygon(data = diamond_ci_region,
               aes(x = x, y = y, fill = EPA_region, group = diamond_id),
               alpha = 1, linewidth = 1, color = NA) +
  
  geom_point(data = group_summary_ci_region,
             aes(x = x0, y = y0, color = EPA_region),
             size = 5, stroke = 1.2) +
  
  # National summary diamonds
  geom_polygon(data = {
    h_diamond <- tibble(
      x = c(nation_ci$x0 - nation_ci$rx,
            nation_ci$x0,
            nation_ci$x0 + nation_ci$rx,
            nation_ci$x0),
      y = c(nation_ci$y0,
            nation_ci$y0 + 0.0002,
            nation_ci$y0,
            nation_ci$y0 - 0.0002),
      diamond_id = "national_x"
    )
    v_diamond <- tibble(
      x = c(nation_ci$x0,
            nation_ci$x0 + 0.0005,
            nation_ci$x0,
            nation_ci$x0 - 0.0005),
      y = c(nation_ci$y0 - nation_ci$ry,
            nation_ci$y0,
            nation_ci$y0 + nation_ci$ry,
            nation_ci$y0),
      diamond_id = "national_y"
    )
    bind_rows(h_diamond, v_diamond)
  },
  aes(x = x, y = y, group = diamond_id),
  alpha = 0.7, linewidth = 1, color = NA, fill = "#B11226") +
  
  geom_point(data = nation_ci,
             aes(x = x0, y = y0),
             color = "#B11226", size = 7, shape = 18) +
  geom_text(data = nation_ci,
            aes(x = x0, y = y0),
            label = "National Summary",
            color = "#B11226", size = 7.5, fontface = "bold",
            vjust = 0.4, hjust = -0.1) +
  
  scale_color_manual(values = region_colors) +
  scale_fill_manual(values = region_colors) +
  
  labs(
    x     = bquote(bold("Change from 2019")),
    y     = bquote(bold("Concentration index in 2024")),
    color = "EPA Region",
    fill  = "EPA Region"
  ) +
  theme_minimal() +
  theme(
    axis.text     = element_text(size = 22),
    axis.title.x  = element_text(size = 25, margin = margin(t = 15, b = 15)),
    axis.title.y  = element_text(size = 25, margin = margin(r = 15, l = 15)),
    legend.title  = element_text(size = 13, face = "bold"),
    legend.text   = element_text(size = 11),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray95"),
    legend.position  = "none"
  ) +
  xlim(-0.03, 0.03) +
  ylim(-0.13, 0.01)
dev.off()

################################################################################
# Summary by URBANICITY ONLY (CI)
################################################################################

group_summary_ci_metro <- final_est %>%
  filter(!is.na(CI_no2_2024), !is.na(CI_dif)) %>%
  mutate(metro4 = factor(metro4, levels = c("Large MSA", "Small MSA",
                                            "MicroSA", "Non-urban"))) %>%
  group_by(metro4) %>%
  summarize(
    n_msa        = n(),
    ci_2024      = mean(CI_no2_2024, na.rm = TRUE),
    var_ci_2024  = var(CI_no2_2024, na.rm = TRUE),
    ci_dif       = mean(CI_dif, na.rm = TRUE),
    var_ci_dif   = var(CI_dif, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    se_ci_2024 = sqrt(var_ci_2024 / n_msa),
    se_ci_dif  = sqrt(var_ci_dif  / n_msa),
    x0 = ci_dif,
    y0 = ci_2024,
    rx = 1.96 * se_ci_dif,
    ry = 1.96 * se_ci_2024
  )

################################################################################
# Diamond shapes for urbanicity (CI)
################################################################################

diamond_ci_metro_x <- bind_rows(lapply(seq_len(nrow(group_summary_ci_metro)), function(i) {
  tibble(
    metro4     = group_summary_ci_metro$metro4[i],
    diamond_id = paste0(group_summary_ci_metro$metro4[i], "_x"),
    x = c(group_summary_ci_metro$x0[i] - group_summary_ci_metro$rx[i],
          group_summary_ci_metro$x0[i],
          group_summary_ci_metro$x0[i] + group_summary_ci_metro$rx[i],
          group_summary_ci_metro$x0[i]),
    y = c(group_summary_ci_metro$y0[i],
          group_summary_ci_metro$y0[i] + 0.0005,
          group_summary_ci_metro$y0[i],
          group_summary_ci_metro$y0[i] - 0.0005)
  )
}))

diamond_ci_metro_y <- bind_rows(lapply(seq_len(nrow(group_summary_ci_metro)), function(i) {
  tibble(
    metro4     = group_summary_ci_metro$metro4[i],
    diamond_id = paste0(group_summary_ci_metro$metro4[i], "_y"),
    x = c(group_summary_ci_metro$x0[i],
          group_summary_ci_metro$x0[i] + 0.0002,
          group_summary_ci_metro$x0[i],
          group_summary_ci_metro$x0[i] - 0.0002),
    y = c(group_summary_ci_metro$y0[i] - group_summary_ci_metro$ry[i],
          group_summary_ci_metro$y0[i],
          group_summary_ci_metro$y0[i] + group_summary_ci_metro$ry[i],
          group_summary_ci_metro$y0[i])
  )
}))

diamond_ci_metro <- bind_rows(diamond_ci_metro_x, diamond_ci_metro_y)

################################################################################
# Plot 2: URBANICITY STRATIFICATION (CI)
################################################################################

pdf(file = "/Users/k_sy_n_imac/Library/CloudStorage/Dropbox/GWU/Lab/FINESST/Research/Results/Quadrant_CI_Urbanicity.pdf",
  width = 9, height = 9)

ggplot() +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.5) +
  geom_hline(yintercept = 0.13, linetype = "dashed", color = "gray40", linewidth = 0.5) +
  geom_hline(yintercept = -0.13, linetype = "dashed", color = "gray40", linewidth = 0.5) +
  
  geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.5) +
  geom_vline(xintercept = 0.03, linetype = "dashed", color = "gray40", linewidth = 0.5) +
  geom_vline(xintercept = -0.03, linetype = "dashed", color = "gray40", linewidth = 0.5) +
  
  geom_polygon(data = diamond_ci_metro,
               aes(x = x, y = y, group = diamond_id),
               fill = "black", alpha = 0.7, linewidth = 1, color = NA) +
  
  geom_point(data = group_summary_ci_metro,
             aes(x = x0, y = y0, shape = metro4),
             color = "black", size = 4, stroke = 1.2) +
  
  # National summary diamonds
  geom_polygon(data = {
    h_diamond <- tibble(
      x = c(nation_ci$x0 - nation_ci$rx,
            nation_ci$x0,
            nation_ci$x0 + nation_ci$rx,
            nation_ci$x0),
      y = c(nation_ci$y0,
            nation_ci$y0 + 0.0003,
            nation_ci$y0,
            nation_ci$y0 - 0.0003),
      diamond_id = "national_x"
    )
    v_diamond <- tibble(
      x = c(nation_ci$x0,
            nation_ci$x0 + 0.0004,
            nation_ci$x0,
            nation_ci$x0 - 0.0004),
      y = c(nation_ci$y0 - nation_ci$ry,
            nation_ci$y0,
            nation_ci$y0 + nation_ci$ry,
            nation_ci$y0),
      diamond_id = "national_y"
    )
    bind_rows(h_diamond, v_diamond)
  },
  aes(x = x, y = y, group = diamond_id),
  alpha = 0.7, linewidth = 1, color = NA, fill = "#B11226") +
  
  geom_point(data = nation_ci,
             aes(x = x0, y = y0),
             color = "#B11226", size = 7, shape = 18) +
  geom_text(data = nation_ci,
            aes(x = x0, y = y0),
            label = "National Summary",
            color = "#B11226", size = 7.5, fontface = "bold",
            vjust = 0.5, hjust = -0.1) +
  
  scale_shape_manual(values = c(
    "Large MSA" = 16,
    "Small MSA" = 17,
    "MicroSA"   = 15,
    "Non-urban" = 8
  )) +
  
  labs(
    x     = bquote(bold("Change from 2019")),
    y     = bquote(bold("Concentration index in 2024")),
    shape = "Urbanicity"
  ) +
  theme_minimal() +
  theme(
    axis.text     = element_text(size = 22),
    axis.title.x  = element_text(size = 25, margin = margin(t = 15, b = 15)),
    axis.title.y  = element_text(size = 25, margin = margin(r = 15, l = 15)),
    legend.title  = element_text(size = 13, face = "bold"),
    legend.text   = element_text(size = 11),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray95"),
    legend.position  = "none"
  ) +
  xlim(-0.03,0.03) +
  ylim(-0.13, 0.01)

dev.off()












###############################################################################
# Two-way ANOVA ###############################################################################
###############################################################################

# NO2 level in 2024
model_no2 <- lm(CI_no2_2024 ~ EPA_region + metro4, data = plot_df)
anova(model_no2)

# NO2 percent change
model_change <- lm(CI_dif ~ EPA_region + metro4, data = plot_df)
anova(model_change)









#######################################################################################################################################
# Figure 7: TREEMAPS ##################################################################################################################
#######################################################################################################################################
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

# ── Named region labels ───────────────────────────────────────────────────────
region_labels <- c(
  "Region 1\n(New England)",
  "Region 2\n(Northeast/Caribbean)",
  "Region 3\n(Mid-Atlantic)",
  "Region 4\n(Southeast)",
  "Region 5\n(Great Lakes)",
  "Region 6\n(South Central)",
  "Region 7\n(Heartland)",
  "Region 8\n(Mountains/Plains)",
  "Region 9\n(Pacific Southwest)",
  "Region 10\n(Northwest)"
)

recode_region <- function(x) {
  recode(x,
         "Region 1"  = "Region 1\n(New England)",
         "Region 2"  = "Region 2\n(Northeast/Caribbean)",
         "Region 3"  = "Region 3\n(Mid-Atlantic)",
         "Region 4"  = "Region 4\n(Southeast)",
         "Region 5"  = "Region 5\n(Great Lakes)",
         "Region 6"  = "Region 6\n(South Central)",
         "Region 7"  = "Region 7\n(Heartland)",
         "Region 8"  = "Region 8\n(Mountains/Plains)",
         "Region 9"  = "Region 9\n(Pacific Southwest)",
         "Region 10" = "Region 10\n(Northwest)"
  )
}

# ── Quadrant colors ───────────────────────────────────────────────────────────
#quad_colors <- c(
#  "Q1" = "#1a9896",
#  "Q2" = "#8b78c0",
#  "Q3" = "#e63232",
#  "Q4" = "#ffd166",
#  "Q5" = "#90c4dd",
#  "Q6" = "#8b5e72",
#  "Q7" = "#f4a261",
#  "Q8" = "#5a9e8a",
#  "Q9" = "#cccccc"
#)

quad_colors <- c(
  "Q1" = "#fcc5c0",  # light red    (top-right)
  "Q2" = "#c0152a",  # dark red     (top-left)
  "Q3" = "#1a6bb5",  # dark blue    (bottom-left)
  "Q4" = "#c6dbef",  # light blue   (bottom-right)
  "Q5" = "#f1717e",  # mid red      (top-middle)
  "Q6" = "#8e3aab",  # dark purple  (mid-left)
  "Q7" = "#6baed6",  # mid blue     (bottom-middle)
  "Q8" = "#e8cfe8",  # light purple (mid-right)
  "Q9" = "#b07cc6"   # mid purple   (mid-middle)
)

# ─────────────────────────────────────────────────────────────────────────────
# Data prep
# ─────────────────────────────────────────────────────────────────────────────
list_year  <- c(2019, 2020, 2021, 2022, 2023, 2024)
list_cause <- c("All causes", "Cardiovascular diseases",
                "Chronic respiratory diseases", "Ischemic heart disease",
                "Tracheal, bronchus, and lung cancer")
j <- 1

var_names_2024 <- c("MSA_title", "metro4", "EPA_region",
                    "CI_no2_2024", "abs_no2_2024", "ratio_no2_2024",
                    "n_tract_total_2024", "n_tract_valid_2024")
est  <- data.frame(matrix(ncol = length(var_names_2024), nrow = 0))
colnames(est) <- var_names_2024

var_names_2019 <- c("MSA_title",
                    "CI_no2_2019", "abs_no2_2019", "ratio_no2_2019",
                    "n_tract_total_2019", "n_tract_valid_2019")
est1 <- data.frame(matrix(ncol = length(var_names_2019), nrow = 0))
colnames(est1) <- var_names_2019

# ── CI 2024 ───────────────────────────────────────────────────────────────────
no2_2024    <- subset(no2_tract_df, year == 2024)
db_2024     <- subset(final3, year == 2024 & cause_name == list_cause[j])
final4_2024 <- merge(no2_2024, db_2024, by = "GEOID")
final4_2024 <- merge(final4_2024, xxxx, by = "GEOID", all.x = TRUE)
final4_2024 <- final4_2024 %>%
  mutate(
    MSA_title  = if_else(metro4 == "Non-urban",
                         paste0("Non_urban_", state_name.x), MSA_title),
    EPA_region = recode_region(EPA_region)
  )

group_df <- final4_2024 %>% group_by(MSA_title, metro4, EPA_region) %>%
  summarize(n_2024 = n(), .groups = "drop")
msa_regions <- group_df %>% count(MSA_title) %>% filter(n > 1)
majority_regions <- final4_2024 %>%
  group_by(MSA_title, EPA_region) %>%
  summarise(total_pop = sum(pop_over20, na.rm = TRUE), .groups = "drop") %>%
  group_by(MSA_title) %>%
  slice_max(total_pop, n = 1) %>%  # ← Gets region with HIGHEST population
  ungroup() %>%
  dplyr::select(MSA_title, EPA_region_new = EPA_region)
final4_2024 <- final4_2024 %>%
  left_join(majority_regions, by = "MSA_title") %>%
  mutate(EPA_region = if_else(!is.na(EPA_region_new), 
                              EPA_region_new, 
                              EPA_region))

MSA_list_2024 <- final4_2024 %>% distinct(MSA_title, EPA_region, metro4)
for (ii in 1:nrow(MSA_list_2024)) {
  msa_i   <- MSA_list_2024$MSA_title[ii]
  epa_i   <- MSA_list_2024$EPA_region[ii]
  metro_i <- MSA_list_2024$metro4[ii]
  
  final4_metro <- final4_2024 %>%
    filter(MSA_title == msa_i, 
           EPA_region == epa_i, 
           metro4 == metro_i)
  
  n_total      <- nrow(final4_metro)
  final_clean  <- final4_metro %>% filter(!is.na(pct_white), !is.na(surface_no2))
  n_valid      <- nrow(final_clean)
  if (n_valid < 5) next
  ecdf_fn     <- stats::ecdf(final_clean[["pct_white"]])
  final_clean <- final_clean[order(final_clean[["pct_white"]]), ]
  final_clean$x_pct    <- ecdf_fn(final_clean[["pct_white"]]) * 100
  final_clean <- final_clean[order(final_clean$x_pct), ]
  final_clean$cumsum_y <- cumsum(final_clean$surface_no2) /
    sum(final_clean$surface_no2) * 100
  auc    <- pracma::trapz(final_clean$x_pct, final_clean$cumsum_y)
  CI_no2 <- (5000 - auc) / 5000
  final_dec <- final_clean %>%
    mutate(white_decile = ntile(pct_white, 10)) %>%
    group_by(white_decile) %>%
    summarise(mean_no2 = mean(surface_no2, na.rm = TRUE), .groups = "drop") %>%
    right_join(data.frame(white_decile = 1:10), by = "white_decile") %>%
    arrange(white_decile)
  abs_no2   <- unname(-final_dec$mean_no2[1] + final_dec$mean_no2[10])
  ratio_no2 <- unname(final_dec$mean_no2[1]  / final_dec$mean_no2[10])
  EPA_i     <- unique(final4_metro$EPA_region); EPA_i   <- EPA_i[!is.na(EPA_i)][1]
  metro_i   <- unique(final4_metro$metro4);     metro_i <- metro_i[!is.na(metro_i)][1]
  est <- rbind(est, data.frame(
    MSA_title = msa_i, metro4 = metro_i, EPA_region = EPA_i,
    CI_no2_2024 = CI_no2, abs_no2_2024 = abs_no2, ratio_no2_2024 = ratio_no2,
    n_tract_total_2024 = n_total, n_tract_valid_2024 = n_valid,
    stringsAsFactors = FALSE
  ))
}
rm(no2_2024, db_2024, final4_2024, final4_metro, final_clean, ecdf_fn, final_dec)
gc()

# ── CI 2019 ───────────────────────────────────────────────────────────────────
no2_2019    <- subset(no2_tract_df, year == 2019)
db_2019     <- subset(final3, year == 2019 & cause_name == list_cause[j])
final4_2019 <- merge(no2_2019, db_2019, by = "GEOID")
final4_2019 <- merge(final4_2019, xxxx, by = "GEOID", all.x = TRUE)
final4_2019 <- final4_2019 %>%
  mutate(
    MSA_title  = if_else(metro4 == "Non-urban",
                         paste0("Non_urban_", state_name.x), MSA_title),
    EPA_region = recode_region(EPA_region)
  )

group_df <- final4_2019 %>% group_by(MSA_title, metro4, EPA_region) %>%
  summarize(n_2024 = n(), .groups = "drop")
msa_regions <- group_df %>% count(MSA_title) %>% filter(n > 1)
majority_regions <- final4_2019 %>%
  group_by(MSA_title, EPA_region) %>%
  summarise(total_pop = sum(pop_over20, na.rm = TRUE), .groups = "drop") %>%
  group_by(MSA_title) %>%
  slice_max(total_pop, n = 1) %>%  # ← Gets region with HIGHEST population
  ungroup() %>%
  dplyr::select(MSA_title, EPA_region_new = EPA_region)
final4_2019 <- final4_2019 %>%
  left_join(majority_regions, by = "MSA_title") %>%
  mutate(EPA_region = if_else(!is.na(EPA_region_new), 
                              EPA_region_new, 
                              EPA_region))

MSA_list_2019 <- final4_2019 %>% distinct(MSA_title, EPA_region, metro4)
for (ii in 1:nrow(MSA_list_2024)) {
  msa_i   <- MSA_list_2019$MSA_title[ii]
  epa_i   <- MSA_list_2019$EPA_region[ii]
  metro_i <- MSA_list_2019$metro4[ii]
  
  final4_metro <- final4_2019 %>%
    filter(MSA_title == msa_i, 
           EPA_region == epa_i, 
           metro4 == metro_i) 
  n_total      <- nrow(final4_metro)
  final_clean  <- final4_metro %>% filter(!is.na(pct_white), !is.na(surface_no2))
  n_valid      <- nrow(final_clean)
  if (n_valid < 5) next
  ecdf_fn     <- stats::ecdf(final_clean[["pct_white"]])
  final_clean <- final_clean[order(final_clean[["pct_white"]]), ]
  final_clean$x_pct    <- ecdf_fn(final_clean[["pct_white"]]) * 100
  final_clean <- final_clean[order(final_clean$x_pct), ]
  final_clean$cumsum_y <- cumsum(final_clean$surface_no2) /
    sum(final_clean$surface_no2) * 100
  auc    <- pracma::trapz(final_clean$x_pct, final_clean$cumsum_y)
  CI_no2 <- (5000 - auc) / 5000
  final_dec <- final_clean %>%
    mutate(white_decile = ntile(pct_white, 10)) %>%
    group_by(white_decile) %>%
    summarise(mean_no2 = mean(surface_no2, na.rm = TRUE), .groups = "drop") %>%
    right_join(data.frame(white_decile = 1:10), by = "white_decile") %>%
    arrange(white_decile)
  abs_no2   <- unname(-final_dec$mean_no2[1] + final_dec$mean_no2[10])
  ratio_no2 <- unname(final_dec$mean_no2[1]  / final_dec$mean_no2[10])
  est1 <- rbind(est1, data.frame(
    MSA_title = msa_i, CI_no2_2019 = CI_no2, abs_no2_2019 = abs_no2,
    ratio_no2_2019 = ratio_no2, n_tract_total_2019 = n_total,
    n_tract_valid_2019 = n_valid, stringsAsFactors = FALSE
  ))
}
rm(no2_2019, db_2019, final4_2019, final4_metro, final_clean, ecdf_fn, final_dec)
gc()

# ── Merge & quadrant assignment ───────────────────────────────────────────────
final_est        <- full_join(est, est1, by = "MSA_title")
final_est$CI_dif <- final_est$CI_no2_2024 - final_est$CI_no2_2019

group_df <- final_est %>%
  mutate(
    quadrant = case_when(
      CI_no2_2024 <  -0.13 & CI_dif <  -0.03 ~ "Marginalized overburdened & Burden to marginalized",
      CI_no2_2024 <  -0.13 & CI_dif >   0.03 ~ "Marginalized overburdened & Burden to non-marginalized",
      CI_no2_2024 >   0.13 & CI_dif >   0.03 ~ "Non-marginalized overburdened & Burden to non-marginalized",
      CI_no2_2024 >   0.13 & CI_dif <  -0.03 ~ "Non-marginalized overburdened & Burden to marginalized",
      CI_no2_2024 <  -0.13 & between(CI_dif, -0.03, 0.03) ~ "Marginalized overburdened & stable",
      CI_no2_2024 >   0.13 & between(CI_dif, -0.03, 0.03) ~ "Non-marginalized overburdened & stable",
      between(CI_no2_2024, -0.13, 0.13) & CI_dif <  -0.03 ~ "Neutral & Burden to marginalized",
      between(CI_no2_2024, -0.13, 0.13) & CI_dif >   0.03 ~ "Neutral & Burden to non-marginalized",
      between(CI_no2_2024, -0.13, 0.13) & between(CI_dif, -0.03, 0.03) ~ "Neutral & stable"
    )
  )

# ─────────────────────────────────────────────────────────────────────────────
# Levels
# ─────────────────────────────────────────────────────────────────────────────
metro_levels   <- c("Large MSA", "Small MSA", "MicroSA", "Non-urban")
region_levels  <- region_labels
metro_levels2  <- c(metro_levels,  "Total")
region_levels2 <- c(region_levels, "Total")

# ─────────────────────────────────────────────────────────────────────────────
# Aggregate cells
# ─────────────────────────────────────────────────────────────────────────────
obs_df_main <- group_df %>%
  group_by(metro4, EPA_region, quadrant) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(metro4, EPA_region) %>%
  mutate(prop = n / sum(n)) %>% ungroup()

obs_df_row_total <- group_df %>%
  group_by(metro4, quadrant) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(metro4) %>%
  mutate(prop = n / sum(n)) %>% ungroup() %>%
  mutate(EPA_region = "Total")

obs_df_col_total <- group_df %>%
  group_by(EPA_region, quadrant) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(EPA_region) %>%
  mutate(prop = n / sum(n)) %>% ungroup() %>%
  mutate(metro4 = "Total")

obs_df_grand_total <- group_df %>%
  group_by(quadrant) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(prop = n / sum(n), metro4 = "Total", EPA_region = "Total")

obs_df <- bind_rows(obs_df_main, obs_df_row_total,
                    obs_df_col_total, obs_df_grand_total)

# ─────────────────────────────────────────────────────────────────────────────
# Wide format
# ─────────────────────────────────────────────────────────────────────────────
plot_base <- obs_df %>%
  dplyr::select(metro4, EPA_region, quadrant, prop) %>%
  pivot_wider(names_from = quadrant, values_from = prop, values_fill = 0) %>%
  {
    missing_cols <- setdiff(c(
      "Non-marginalized overburdened & Burden to non-marginalized",
      "Non-marginalized overburdened & Burden to marginalized",
      "Non-marginalized overburdened & stable",
      "Marginalized overburdened & Burden to marginalized",
      "Marginalized overburdened & Burden to non-marginalized",
      "Marginalized overburdened & stable",
      "Neutral & Burden to marginalized",
      "Neutral & Burden to non-marginalized",
      "Neutral & stable"
    ), names(.))
    for (col in missing_cols) { . <- mutate(., !!col := 0) }
    .
  } %>%
  transmute(
    metro4, EPA_region,
    Q1 = `Non-marginalized overburdened & Burden to non-marginalized`,
    Q2 = `Non-marginalized overburdened & Burden to marginalized`,
    Q3 = `Marginalized overburdened & Burden to marginalized`,
    Q4 = `Marginalized overburdened & Burden to non-marginalized`,
    Q5 = `Non-marginalized overburdened & stable`,
    Q6 = `Neutral & Burden to marginalized`,
    Q7 = `Marginalized overburdened & stable`,
    Q8 = `Neutral & Burden to non-marginalized`,
    Q9 = `Neutral & stable`
  ) %>%
  mutate(
    metro4     = factor(metro4,     levels = metro_levels2),
    EPA_region = factor(EPA_region, levels = region_levels2),
    x = as.numeric(EPA_region),
    y = rev(seq_along(metro_levels2))[match(metro4, metro_levels2)]
  )

# ─────────────────────────────────────────────────────────────────────────────
# Nudge zeros
# ─────────────────────────────────────────────────────────────────────────────
epsilon <- 0.02

plot_base_nudged <- plot_base %>%
  mutate(
    Q1 = if_else(Q1 == 0, epsilon, Q1),
    Q2 = if_else(Q2 == 0, epsilon, Q2),
    Q3 = if_else(Q3 == 0, epsilon, Q3),
    Q4 = if_else(Q4 == 0, epsilon, Q4),
    Q5 = if_else(Q5 == 0, epsilon, Q5),
    Q6 = if_else(Q6 == 0, epsilon, Q6)
  ) %>%
  mutate(
    total = Q1 + Q2 + Q3 + Q4 + Q5 + Q6,
    Q1 = Q1 / total,
    Q2 = Q2 / total,
    Q3 = Q3 / total,
    Q4 = Q4 / total,
    Q5 = Q5 / total,
    Q6 = Q6 / total
  ) %>%
  dplyr::select(-total)

# ─────────────────────────────────────────────────────────────────────────────
# Build rect_df
# ─────────────────────────────────────────────────────────────────────────────
pad   <- 0.09
min_w <- 0.03
min_h <- 0.03

alloc <- function(vals, total_px, min_px, tiny_prop = 0.1) {
  
  vals <- replace_na(vals, 0)
  
  # Cells that are zero OR very close to zero
  # will receive exactly the same size as NA/zero cells.
  tiny_cells <- vals <= tiny_prop
  
  n_tiny <- sum(tiny_cells)
  reserved <- n_tiny * min_px
  remaining <- max(total_px - reserved, 0)
  
  non_tiny_sum <- sum(vals[!tiny_cells])
  
  out <- numeric(length(vals))
  
  # Tiny or zero cells get the same fixed size
  out[tiny_cells] <- min_px
  
  # Other cells divide the remaining space proportionally
  if (non_tiny_sum > 0) {
    out[!tiny_cells] <- remaining * vals[!tiny_cells] / non_tiny_sum
  } else {
    out[] <- total_px / length(vals)
  }
  
  out
}

# ─────────────────────────────────────────────────────────────────────────────
# Build rect_df
# ─────────────────────────────────────────────────────────────────────────────
rect_df <- plot_base %>%
  rowwise() %>%
  do({
    row <- .
    x_left  <- row$x - 0.5 + pad
    x_right <- row$x + 0.5 - pad
    y_bot   <- row$y - 0.5 + pad
    y_top   <- row$y + 0.5 - pad
    total_w <- x_right - x_left
    total_h <- y_top   - y_bot
    
    # ── X: 3 columns ─────────────────────────────────────────────────────────
    left_raw  <- row$Q2 + row$Q6 + row$Q3
    mid_raw   <- row$Q5 + row$Q9 + row$Q7
    right_raw <- row$Q1 + row$Q8 + row$Q4
    x_widths  <- alloc(c(left_raw, mid_raw, right_raw), total_w, min_w)
    x_mid1    <- x_left + x_widths[1]
    x_mid2    <- x_mid1 + x_widths[2]
    
    # ── Y within left col: Q3 (bot), Q6 (mid), Q2 (top) ─────────────────────
    yl         <- alloc(c(row$Q3, row$Q6, row$Q2), total_h, min_h)
    y_left_bot <- y_bot      + yl[1]
    y_left_top <- y_left_bot + yl[2]
    
    # ── Y within mid col: Q7 (bot), Q9 (mid), Q5 (top) ──────────────────────
    ym        <- alloc(c(row$Q7, row$Q9, row$Q5), total_h, min_h)
    y_mid_bot <- y_bot     + ym[1]
    y_mid_top <- y_mid_bot + ym[2]
    
    # ── Y within right col: Q4 (bot), Q8 (mid), Q1 (top) ────────────────────
    yr          <- alloc(c(row$Q4, row$Q8, row$Q1), total_h, min_h)
    y_right_bot <- y_bot       + yr[1]
    y_right_top <- y_right_bot + yr[2]
    
    tibble(
      metro4     = row$metro4,
      EPA_region = row$EPA_region,
      quadrant   = c("Q1",        "Q2",       "Q3",       "Q4",         "Q5",      "Q6",       "Q7",     "Q8",         "Q9"),
      prob       = c(row$Q1,      row$Q2,     row$Q3,     row$Q4,       row$Q5,    row$Q6,     row$Q7,   row$Q8,       row$Q9),
      xmin       = c(x_mid2,      x_left,     x_left,     x_mid2,       x_mid1,    x_left,     x_mid1,   x_mid2,       x_mid1),
      xmax       = c(x_right,     x_mid1,     x_mid1,     x_right,      x_mid2,    x_mid1,     x_mid2,   x_right,      x_mid2),
      ymin       = c(y_right_top, y_left_top, y_bot,      y_bot,        y_mid_top, y_left_bot, y_bot,    y_right_bot,  y_mid_bot),
      ymax       = c(y_top,       y_top,      y_left_bot, y_right_bot,  y_top,     y_left_top, y_mid_bot, y_right_top, y_mid_top)
    )
  }) %>%
  ungroup()

# ── Join counts ───────────────────────────────────────────────────────────────
rect_df <- rect_df %>%
  left_join(
    obs_df %>%
      mutate(
        quadrant = recode(quadrant,
                          "Non-marginalized overburdened & Burden to non-marginalized" = "Q1",
                          "Non-marginalized overburdened & Burden to marginalized"     = "Q2",
                          "Marginalized overburdened & Burden to marginalized"         = "Q3",
                          "Marginalized overburdened & Burden to non-marginalized"     = "Q4",
                          "Non-marginalized overburdened & stable"                     = "Q5",
                          "Neutral & Burden to marginalized"                           = "Q6",
                          "Marginalized overburdened & stable"                         = "Q7",
                          "Neutral & Burden to non-marginalized"                       = "Q8",
                          "Neutral & stable"                                           = "Q9"
        ),
        metro4     = factor(metro4,     levels = metro_levels2),
        EPA_region = factor(EPA_region, levels = region_levels2)
      ) %>%
      dplyr::select(metro4, EPA_region, quadrant, n),
    by = c("metro4", "EPA_region", "quadrant")
  ) %>%
  mutate(
    n        = replace_na(n, 0),
    fill_col = if_else(n == 0, "empty", quadrant)
  )

# ─────────────────────────────────────────────────────────────────────────────
# label_df
# ─────────────────────────────────────────────────────────────────────────────
narrow_width_threshold <- 1 / 6
bottom_text_pad <- 0.06

label_df <- obs_df %>%
  mutate(
    metro4     = factor(metro4,     levels = metro_levels2),
    EPA_region = factor(EPA_region, levels = region_levels2),
    quadrant   = recode(quadrant,
                        "Non-marginalized overburdened & Burden to non-marginalized" = "Q1",
                        "Non-marginalized overburdened & Burden to marginalized"     = "Q2",
                        "Marginalized overburdened & Burden to marginalized"         = "Q3",
                        "Marginalized overburdened & Burden to non-marginalized"     = "Q4",
                        "Non-marginalized overburdened & stable"                     = "Q5",
                        "Neutral & Burden to marginalized"                           = "Q6",
                        "Marginalized overburdened & stable"                         = "Q7",
                        "Neutral & Burden to non-marginalized"                       = "Q8",
                        "Neutral & stable"                                           = "Q9"
    )
  ) %>%
  left_join(
    rect_df %>%
      mutate(
        x = (xmin + xmax) / 2,
        y = (ymin + ymax) / 2,
        area = (xmax - xmin) * (ymax - ymin),
        cell_width = xmax - xmin,
        cell_height = ymax - ymin
      ) %>%
      dplyr::select(
        metro4,
        EPA_region,
        quadrant,
        x,
        y,
        area,
        xmin,
        xmax,
        ymin,
        ymax,
        cell_width,
        cell_height
      ),
    by = c("metro4", "EPA_region", "quadrant")
  ) %>%
  mutate(
    label = case_when(
      n == 0        ~ NA_character_,
      area >= 0.002 ~ paste0(n, "\n(", round(prop * 100), "%)"),
      TRUE          ~ NA_character_
    ),
    
    txt_size = case_when(
      area >= 0.008 ~ 5.2,
      area >= 0.002 ~ 5.2,
      TRUE          ~ 5.2
    ),
    
    narrow_cell = cell_width < narrow_width_threshold,
    short_cell = cell_height < (min_h + 0.1),
    # Only align at bottom if narrow AND not short
    push_to_bottom = narrow_cell & !short_cell,
    y = if_else(
      push_to_bottom,
      ymin + bottom_text_pad,
      y
    ),
    text_vjust = if_else(
      push_to_bottom,
      0,
      0.5
    )
  ) %>%
  filter(!is.na(label))

label_df <- label_df %>%
  mutate(
    x_original = x,
    y_original = y,
    label_length = nchar(gsub("\n", "", label))
  )

two_line_shift <- 0.16

move_shorter_label <- function(
    df,
    cell1,
    cell2,
    overlap_type = c("horizontal", "vertical"),
    shift = two_line_shift
) {
  
  overlap_type <- match.arg(overlap_type)
  
  idx1 <- which(
    as.character(df$metro4) == cell1$metro4 &
      as.character(df$EPA_region) == cell1$EPA_region &
      as.character(df$quadrant) == cell1$quadrant
  )
  
  idx2 <- which(
    as.character(df$metro4) == cell2$metro4 &
      as.character(df$EPA_region) == cell2$EPA_region &
      as.character(df$quadrant) == cell2$quadrant
  )
  
  if (length(idx1) == 0 | length(idx2) == 0) {
    return(df)
  }
  
  shorter_idx <- ifelse(
    df$label_length[idx1] <= df$label_length[idx2],
    idx1,
    idx2
  )
  
  if (overlap_type == "horizontal") {
    
    # If two labels overlap horizontally,
    # move the shorter one downward.
    df$y[shorter_idx] <- df$y[shorter_idx] - shift
    
    # Keep the label inside its rectangle.
    df$y[shorter_idx] <- pmax(
      df$ymin[shorter_idx] + 0.06,
      df$y[shorter_idx]
    )
  }
  
  if (overlap_type == "vertical") {
    
    # If two labels overlap vertically,
    # move the shorter one toward the center of its own rectangle.
    cell_center_x <- (df$xmin[shorter_idx] + df$xmax[shorter_idx]) / 2
    cell_center_y <- (df$ymin[shorter_idx] + df$ymax[shorter_idx]) / 2
    
    df$x[shorter_idx] <- cell_center_x
    df$y[shorter_idx] <- cell_center_y
    
    # Optional slight inward/downward adjustment.
    df$y[shorter_idx] <- df$y[shorter_idx] - shift
    
    # Keep label inside its rectangle.
    df$x[shorter_idx] <- pmin(
      df$xmax[shorter_idx] - 0.04,
      pmax(df$xmin[shorter_idx] + 0.04, df$x[shorter_idx])
    )
    
    df$y[shorter_idx] <- pmin(
      df$ymax[shorter_idx] - 0.06,
      pmax(df$ymin[shorter_idx] + 0.06, df$y[shorter_idx])
    )
  }
  
  df
}

# ─────────────────────────────────────────────────────────────────────────────
# Legend
# ─────────────────────────────────────────────────────────────────────────────
library(ggtext)

x1 <- -1; x2 <- -0.28; x3 <- 0.28; x4 <- 1
y1 <- -1; y2 <- -0.28; y3 <- 0.28; y4 <- 1   
gap <- 0.01

legend_plot <- ggplot() +
  annotate("rect", xmin = x1,     xmax = x2-gap, ymin = y3+gap, ymax = y4,     fill = quad_colors["Q2"]) +
  annotate("rect", xmin = x2+gap, xmax = x3-gap, ymin = y3+gap, ymax = y4,     fill = quad_colors["Q5"]) +
  annotate("rect", xmin = x3+gap, xmax = x4,     ymin = y3+gap, ymax = y4,     fill = quad_colors["Q1"]) +
  annotate("rect", xmin = x1,     xmax = x2-gap, ymin = y2+gap, ymax = y3-gap, fill = quad_colors["Q6"]) +
  annotate("rect", xmin = x2+gap, xmax = x3-gap, ymin = y2+gap, ymax = y3-gap, fill = quad_colors["Q9"]) +
  annotate("rect", xmin = x3+gap, xmax = x4,     ymin = y2+gap, ymax = y3-gap, fill = quad_colors["Q8"]) +
  annotate("rect", xmin = x1,     xmax = x2-gap, ymin = y1,     ymax = y2-gap, fill = quad_colors["Q3"]) +
  annotate("rect", xmin = x2+gap, xmax = x3-gap, ymin = y1,     ymax = y2-gap, fill = quad_colors["Q7"]) +
  annotate("rect", xmin = x3+gap, xmax = x4,     ymin = y1,     ymax = y2-gap, fill = quad_colors["Q4"]) +
  ggtext::geom_richtext(
    aes(x = (x1+x2)/2, y = (y3+y4)/2),
    label = "<b>S2<br>Non-marg NO<sub>2</sub>↑<br> &amp; Exposure to<br>marg</b>",
    size = 3.5,
    color = "white",
    fill = NA,
    label.color = NA,
    lineheight = 0.85) +
  ggtext::geom_richtext(
    aes(x = 0, y = (y3+y4)/2),
    label = "<b>S5<br>Non-marg<br>NO<sub>2</sub>↑ &amp; Stable</b>",
    size = 3.5,
    color = "grey20",
    fill = NA,
    label.color = NA,
    lineheight = 0.85) +
  ggtext::geom_richtext(
    aes(x = (x3+x4)/2, y = (y3+y4)/2),
    label = "<b>S1<br>Non-marg NO<sub>2</sub>↑<br>&amp; Exposure to<br>non-marg</b>",
    size = 3.5,
    color = "grey20",
    fill = NA,
    label.color = NA,
    lineheight = 0.85) +
  ggtext::geom_richtext(
    aes(x = (x1+x2)/2, y = (y1+y2)/2),
    label = "<b>S3<br>Marg NO<sub>2</sub>↑<br>&amp; Exposure to<br>marg</b>",
    size = 3.5,
    color = "white",
    fill = NA,
    label.color = NA,
    lineheight = 0.85) +
  ggtext::geom_richtext(
    aes(x = 0, y = (y1+y2)/2),
    label = "<b>S7<br>Marg NO<sub>2</sub>↑<br>&amp; Stable</b>",
    size = 3.5,
    color = "grey20",
    fill = NA,
    label.color = NA,
    lineheight = 0.85) +
  ggtext::geom_richtext(
    aes(x = (x3+x4)/2, y = (y1+y2)/2),
    label = "<b>S4<br>Marg NO<sub>2</sub>↑<br>&amp; Exposure to<br>non-marg</b>",
    size = 3.5,
    color = "grey20",
    fill = NA,
    label.color = NA,
    lineheight = 0.85) + 
  annotate("text", x = (x1+x2)/2, y = 0,
           label = "S6\nNeutral &\nExposure to\nmarg",
           size = 3.5, fontface = "bold", color = "white", lineheight = 0.85) +
  annotate("text", x = 0,         y = 0,
           label = "S9\nNeutral\n& Stable",
           size = 3.5, fontface = "bold", color = "grey20", lineheight = 0.85) +
  annotate("text", x = (x3+x4)/2, y = 0,
           label = "S8\nNeutral &\nExposure to\nnon-marg",
           size = 3.5, fontface = "bold", color = "grey20", lineheight = 0.85) +
  
  annotate("text", x = 0, y = -1.1, label = "Change from 2019", size = 3.5) +
  annotate("text", x = -1.1, y = 0, label = "Concentration index in 2024", size = 3.5, angle = 90, lineheight = 0.85) +
  coord_fixed(xlim = c(x1-0.55, x4+0.2), ylim = c(y1-0.3, y4+0.2), clip = "off") +
  theme_void() +
  theme(plot.margin = margin(20, 10, 20, 0))


# ─────────────────────────────────────────────────────────────────────────────
# Treemap
# ─────────────────────────────────────────────────────────────────────────────
x_faces <- ifelse(region_levels2 == "Total", "bold", "plain")
y_faces <- ifelse(metro_levels2  == "Total", "bold", "plain")

treemap_plot <- ggplot(rect_df) +
  geom_rect(
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill_col),
    color = "white", linewidth = 0.5
  ) +
  scale_fill_manual(values = c(quad_colors, "empty" = "grey40")) +
  guides(fill = "none") +
  geom_text(
    data      = label_df,
    aes(x = x, y = y, label = label, size = txt_size),
    lineheight = 0.9
  ) +
  scale_size_identity() +
  scale_x_continuous(
    breaks   = seq_along(region_levels2),
    labels   = region_levels2,
    expand   = c(0.02, 0.02),
    position = "top"
  ) +
  scale_y_continuous(
    breaks = rev(seq_along(metro_levels2)),
    labels = metro_levels2,
    expand = c(0.0, 0.00)
  ) +
  coord_fixed(clip = "off") +
  labs(x = "", y = "") +
  theme_minimal(base_family = "Lato", base_size = 13) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid       = element_blank(),
    panel.border     = element_rect(color = NA, fill = NA, linewidth = 0.8),
    axis.text.x.top      = element_text(
      angle = 30, color="black",
      hjust = 0.5,
      vjust = 0.5,
      face = x_faces,
      size = 20, margin = margin(b = 5)
    ),
    axis.text.y      = element_text(
      size = 20, 
      color="black", 
      face = y_faces,
      angle = 30,        # ← Tilt y-axis text
      hjust = 0.5,
      vjust = 1,
      margin = margin(r = 0)
    ),
    plot.margin      = margin(t = 15, r = 15, b = 30, l = 5) ,
    axis.ticks       = element_blank(),
  )

# ─────────────────────────────────────────────────────────────────────────────
# Save
# ─────────────────────────────────────────────────────────────────────────────
pdf(file = "/Users/k_sy_n_imac/Library/CloudStorage/Dropbox/GWU/Lab/FINESST/Research/Results/Quadrant_CI_treemap_with_totals.pdf",
    width = 26, height = 12)
print(treemap_plot)
dev.off()

pdf(file = "/Users/k_sy_n_imac/Library/CloudStorage/Dropbox/GWU/Lab/FINESST/Research/Results/Quadrant_CI_treemap_with_totals(legend).pdf",
    width = 5.5, height = 5.5)
print(legend_plot)
dev.off()


































#####################################################################################################
# Additional test/summary for discussion ############################################################
#####################################################################################################
# Checking Population age-structure
pop <- read.csv("/Users/k_sy_n_imac/Library/CloudStorage/Dropbox/GWU/Lab/FINESST/Research/Data/DECENNIALDP2020/pop.csv")
pop_long <- pop %>%
  pivot_longer(cols = 3:21, names_to = "age_id", values_to = "pop") %>%
  dplyr::select(GEOID, NAME, age_id, pop) 
pop_long$new_age_id <- ifelse(pop_long$age_id == "age_20", 20,
                              ifelse(pop_long$age_id == "age_25", 25,
                                     ifelse(pop_long$age_id == "age_30", 30,
                                            ifelse(pop_long$age_id == "age_35", 35,
                                                   ifelse(pop_long$age_id == "age_40", 40,        
                                                          ifelse(pop_long$age_id == "age_45", 45,
                                                                 ifelse(pop_long$age_id == "age_50", 50,
                                                                        ifelse(pop_long$age_id == "age_55", 55,
                                                                               ifelse(pop_long$age_id == "age_60", 60,
                                                                                      ifelse(pop_long$age_id == "age_65", 65,
                                                                                             ifelse(pop_long$age_id == "age_70", 70,
                                                                                                    ifelse(pop_long$age_id == "age_75", 75,
                                                                                                           ifelse(pop_long$age_id == "age_80", 80, 
                                                                                                                  ifelse(pop_long$age_id == "over_85", 85,  NA))))))))))))))
pop_long <- subset(pop_long, pop_long$new_age_id != "NA")
pop_long <- pop_long %>% dplyr::select(-age_id)
pop_long <- pop_long %>%
  mutate(GEOID = str_extract(GEOID, "(?<=US)\\d+"))
pop_long$GEOID <- as.numeric(pop_long$GEOID)
pop_long <- pop_long %>% mutate(state = trimws(sub(".*[;,]\\s*", "", NAME)))

summary <- pop_long %>% group_by(state, new_age_id) %>% summarize(pop=sum(pop))


# Checking Baseline mortality rate
mor <- read.csv("/Users/k_sy_n_imac/Library/CloudStorage/Dropbox/GWU/Lab/FINESST/Research/Data/IHME-GBD_2023/IHME-GBD_2023_DATA-e9f6b249-1.csv", header=TRUE)
mor$new_age_id <- ifelse(mor$age_id == 9, 20,
                         ifelse(mor$age_id == 10, 25,
                                ifelse(mor$age_id == 11, 30,
                                       ifelse(mor$age_id == 12, 35,
                                              ifelse(mor$age_id == 13, 40,        
                                                     ifelse(mor$age_id == 14, 45,
                                                            ifelse(mor$age_id == 15, 50,
                                                                   ifelse(mor$age_id == 16, 55,
                                                                          ifelse(mor$age_id == 17, 60,
                                                                                 ifelse(mor$age_id == 18, 65,
                                                                                        ifelse(mor$age_id == 19, 70,
                                                                                               ifelse(mor$age_id == 20, 75,
                                                                                                      ifelse(mor$age_id == 30, 80, 
                                                                                                             ifelse(mor$age_id == 160, 85, NA))))))))))))))
mor <- rename(mor, state_name=location_name, mor_rate=val)
mor <- subset(mor, mor$metric_name=="Rate")
mor <- mor[c("state_name", "new_age_id", "age_name", "cause_name", "year", "mor_rate")]
mor <- mor %>% filter(!is.na(new_age_id))

mor$cause_name1 <- ifelse(mor$cause_name %in% c("Lower respiratory infections", "Upper respiratory infections", "Chronic respiratory diseases"), "Respiratory diseases", mor$cause_name)
mor <- mor %>%
  group_by(state_name, new_age_id, age_name, cause_name1, year) %>%
  summarize(mor_rate = sum(mor_rate), .groups = "drop")

xx <- subset(mor, mor$year==2023) #no data for 2022-. Use the same mor rate in 2021
xx$year <- 2024
mor <- rbind(mor, xx) 
mor <- subset(mor, cause_name1 == "All causes")
summary <- mor %>% group_by(state_name, year) %>% summarize(mor=sum(mor_rate))
summary <- mor %>% group_by(age_name, year) %>% summarize(mor=sum(mor_rate))
summary



# Temporal: NO2
#Merge & Choose specific years/cause of death   **FOR TEMPORAL CHANGE PLOTS**
list_year <- c(2019, 2020, 2021, 2022, 2023, 2024) #data is too large to handle all the years and causes at the same time. Select each to proceed
list_cause <- c("All causes", "Cardiovascular diseases", "Respiratory diseases", "Ischemic heart disease", "Tracheal, bronchus, and lung cancer") 
i1 <- 1
i2 <- 6
j <- 1
no2_ind <- subset(no2_tract_df, year %in% c(list_year[i1], list_year[i2]))
db_ind <- subset(final3, year %in% c(list_year[i1], list_year[i2]) & cause_name==list_cause[j])
final4 <- merge(no2_ind, db_ind, by=c("GEOID", "year"))
final4 <- merge(final4, xxxx, by="GEOID", all.x=TRUE)

#Differences calculation
xxx <- merge(final4, ustract_df[,c("GEOID", "geometry")], by="GEOID")
x <- subset(final4, final4$year==list_year[i1])
x <- rename(x, no2_pre=surface_no2, db_pre=rate_pe)
xx <- subset(final4, final4$year==list_year[i2])
xx <- rename(xx, no2_post=surface_no2, db_post=rate_pe)
xx <- xx[,c("GEOID", "no2_post", "db_post")]

xxx <- merge(x, xx, by="GEOID", all=TRUE)
xxx$abs_dif_no2 <- - (xxx$no2_pre - xxx$no2_post)
xxx$abs_dif_db <- - (xxx$db_pre - xxx$db_post)

keep <- !is.na(xxx$abs_dif_no2) & !is.na(xxx$abs_dif_db) &
  xxx$abs_dif_no2 != 0 & xxx$abs_dif_db != 0 &
  sign(xxx$abs_dif_no2) != sign(xxx$abs_dif_db)
xxx <- xxx[keep, ]

table(xxx$state_name.x)
table(xxx$metro4)










# Merge & Choose specific years/cause of death
list_year <- c(2019, 2020, 2021, 2022, 2023, 2024) #data is too large to handle all the years and causes at the same time. Select each to proceed
list_cause <- c("All causes", "Cardiovascular diseases", "Respiratory diseases", "Ischemic heart disease", "Tracheal, bronchus, and lung cancer") 
i1 <- 1
i2 <- 6
j <- 1
no2_ind <- subset(no2_tract_df, year %in% c(list_year[i1], list_year[i2]))
db_ind <- subset(final3, year %in% c(list_year[i1], list_year[i2]) & cause_name==list_cause[j])
final4 <- merge(no2_ind, db_ind, by=c("GEOID", "year"))
final4 <- merge(final4, xxxx, by="GEOID", all.x=TRUE)

# Calculation
# Nationwide
nation_df <- final4 %>%
  group_by(year) %>%
  summarize(pw_no2 = weighted.mean(surface_no2, total_dc2020, na.rm = TRUE),
            n = n())
pw_no2_2024_nation <- as.numeric(nation_df[2,2])
pct_change_nation  <- as.numeric((nation_df[2,2]-nation_df[1,2])/nation_df[1,2]*100)
# Sub-group
x <- subset(final4, final4$year==list_year[i1])
group_df_2019 <- x %>%
  group_by(EPA_region, metro4) %>%
  summarize(pw_no2_2019 = weighted.mean(surface_no2, total_dc2020, na.rm = TRUE),
            n_2019 = n())
xx <- subset(final4, final4$year==list_year[i2])
group_df_2024 <- xx %>%
  group_by(EPA_region, metro4) %>%
  summarize(pw_no2_2024 = weighted.mean(surface_no2, total_dc2020, na.rm = TRUE),
            n_2024 = n())
group_df <- merge(group_df_2019, group_df_2024, by=c("EPA_region", "metro4"))
group_df$pct_change <- (group_df$pw_no2_2024-group_df$pw_no2_2019)/group_df$pw_no2_2019*100

# Plotting

region_labels <- c(
  "Region 1 (New England)",
  "Region 2 (Northeast/Caribbean)",
  "Region 3 (Mid-Atlantic)",
  "Region 4 (Southeast)",
  "Region 5 (Great Lakes)",
  "Region 6 (South Central)",
  "Region 7 (Heartland)",
  "Region 8 (Mountains/Plains)",
  "Region 9 (Pacific Southwest)",
  "Region 10 (Northwest)"
)

recode_region <- function(x) {
  recode(x,
         "Region 1"  = "Region 1 (New England)",
         "Region 2"  = "Region 2 (Northeast/Caribbean)",
         "Region 3"  = "Region 3 (Mid-Atlantic)",
         "Region 4"  = "Region 4 (Southeast)",
         "Region 5"  = "Region 5 (Great Lakes)",
         "Region 6"  = "Region 6 (South Central)",
         "Region 7"  = "Region 7 (Heartland)",
         "Region 8"  = "Region 8 (Mountains/Plains)",
         "Region 9"  = "Region 9 (Pacific Southwest)",
         "Region 10" = "Region 10 (Northwest)"
  )
}

region_colors <- c(
  "Region 1 (New England)"         = "#1B2F70",
  "Region 2 (Northeast/Caribbean)" = "#2C4FD7",
  "Region 3 (Mid-Atlantic)"        = "#5A8AC6",
  "Region 4 (Southeast)"           = "#8BC8C8",
  "Region 5 (Great Lakes)"         = "#A1D99B",
  "Region 6 (South Central)"       = "#D9C89B",
  "Region 7 (Heartland)"           = "#DDAA66",
  "Region 8 (Mountains/Plains)"    = "#F2A6A0",
  "Region 9 (Pacific Southwest)"   = "#7A5AC9",
  "Region 10 (Northwest)"          = "#4A0066"
)

group_df <- group_df %>%
  mutate(
    metro4     = factor(metro4, levels = c("Large MSA", "Small MSA",
                                           "MicroSA", "Non-urban")),
    EPA_region = factor(recode_region(EPA_region), levels = region_labels)
  )


pdf(file = paste0("/Users/k_sy_n_imac/Library/CloudStorage/Dropbox/GWU/Lab/FINESST/Research/Results/Quadrant_NO2_final.pdf"), 
    width = 10, 
    height = 10)
par(family = "Lato")
ggplot(group_df, aes(x = pct_change, y = pw_no2_2024,
                     color = EPA_region, shape = metro4)) +
  geom_hline(yintercept = 5.31, linetype = "dashed", color = "gray40", linewidth = 0.7) +   # threshold
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.7) +         # no-change line
  geom_point(size = 5, stroke=1.5, alpha = 0.8) +
  scale_color_manual(values = region_colors) +
  scale_shape_manual(
    values = c(
      "Large MSA" = 16,           # circle
      "Small MSA" = 17,       # triangle
      "MicroSA" = 15,       # square
      "Non-urban" = 3     # reverse triangle
    )) +
  geom_point(x = pct_change_nation, y = pw_no2_2024_nation, color = "#B11226", size = 6, shape = 18, inherit.aes = FALSE) +
  geom_text(
    x = pct_change_nation,
    y = pw_no2_2024_nation,
    label = "National Summary",
    color = "#B11226",
    size = 5,
    vjust = -1.5,   # move text slightly above the point
    hjust = 0.5,
    fontface = "bold",
    inherit.aes = FALSE) + 
  theme_minimal() +
  labs(
    x = bquote(bold("% Change in population-weighted NO"[bolditalic(2)] ~ "between 2019 and 2024")),
    y = bquote(bold("Population-weighted NO"[bolditalic(2)] ~ "in 2024 (ppb)")),
    color = "EPA Region",
    shape = "Urbanicity"
  ) +
  annotate("text",
           x = max(group_df$pct_change, na.rm = TRUE),
           y = 5.31,
           label = "WHO AQG (5.3 ppb)",
           hjust = 0.8,          # left aligned
           vjust = -0.8,       # slightly above the line
           size = 4.5,
           color = "gray40",
           fontface = "bold"
  ) +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 13),
    axis.text = element_text(size = 13),
    axis.title.x = element_text(size = 15, margin = margin(t = 15, b = 15)),
    axis.title.y = element_text(size = 15, margin = margin(r = 15, l = 15)),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 13)
  )
dev.off()

# Descriptive
sum(!is.na(xxx$rel_dif))
sum(!is.na(xxx$metro4))









# Figure 6: Quadrant plots of CI ##############################################################################
# Settings
list_year <- c(2019, 2020, 2021, 2022, 2023, 2024)
list_cause <- c("All causes", "Cardiovascular diseases", "Chronic respiratory diseases", "Ischemic heart disease", "Tracheal, bronchus, and lung cancer") 
j <- 1

var_names_2024 <- c("MSA_EPA_combo", "metro4", "EPA_region",
                    "CI_no2_2024", "abs_no2_2024", "ratio_no2_2024",
                    "n_tract_total_2024", "n_tract_valid_2024")
est <- data.frame(matrix(ncol = length(var_names_2024), nrow = 0))
colnames(est) <- var_names_2024

var_names_2019 <- c("MSA_EPA_combo",
                    "CI_no2_2019", "abs_no2_2019", "ratio_no2_2019",
                    "n_tract_total_2019", "n_tract_valid_2019")
est1 <- data.frame(matrix(ncol = length(var_names_2019), nrow = 0))
colnames(est1) <- var_names_2019

# Calculating CI in 2024
no2_2024 <- subset(no2_tract_df, year == 2024)
db_2024  <- subset(final3,       year == 2024 & cause_name == list_cause[j])

final4_2024 <- merge(no2_2024, db_2024, by = "GEOID")
final4_2024 <- merge(final4_2024, xxxx, by = "GEOID", all.x = TRUE)

final4_2024 <- final4_2024 %>%
  mutate(MSA_EPA_combo = paste(MSA_title, EPA_region, sep = " | "))

MSA_list_2024 <- unique(final4_2024$MSA_EPA_combo)

for (ii in seq_along(MSA_list_2024)) {
  
  msa_i <- MSA_list_2024[ii]
  final4_metro <- subset(final4_2024, MSA_EPA_combo == msa_i)
  
  n_total <- nrow(final4_metro)
  
  final_clean_filtered <- final4_metro %>%
    filter(!is.na(pct_white), !is.na(surface_no2))
  
  n_valid <- nrow(final_clean_filtered)
  
  # 5-tract minimum
  if (n_valid < 5) next
  
  # ---- CI-NO2 ----
  ecdf_function <- stats::ecdf(final_clean_filtered[["pct_white"]])
  final_clean_filtered <- final_clean_filtered[order(final_clean_filtered[["pct_white"]]), ]
  final_clean_filtered$x_pct <- ecdf_function(final_clean_filtered[["pct_white"]]) * 100
  final_clean_filtered <- final_clean_filtered[order(final_clean_filtered$x_pct), ]
  final_clean_filtered$cumsum_y <- cumsum(final_clean_filtered$surface_no2) / sum(final_clean_filtered$surface_no2) * 100
  
  auc <- pracma::trapz(final_clean_filtered$x_pct, final_clean_filtered$cumsum_y)
  CI_no2  <- (5000 - auc) / 5000
  
  # ---- Absolute & Ratio NO2 (deciles) ----
  final_dec <- final_clean_filtered %>%
    mutate(white_decile = ntile(pct_white, 10)) %>%
    group_by(white_decile) %>%
    summarise(mean_no2 = mean(surface_no2, na.rm = TRUE), .groups = "drop") %>%
    # force all deciles so we can safely pick 1 and 10
    right_join(data.frame(white_decile = 1:10), by = "white_decile") %>%
    arrange(white_decile)
  
  abs_no2   <- unname(-final_dec$mean_no2[1] + final_dec$mean_no2[10])
  ratio_no2 <- unname(final_dec$mean_no2[1] / final_dec$mean_no2[10])
  
  EPA_i   <- unique(final4_metro$EPA_region); EPA_i <- EPA_i[!is.na(EPA_i)][1]
  metro_i <- unique(final4_metro$metro4);     metro_i <- metro_i[!is.na(metro_i)][1]
  
  est_tmp <- data.frame(
    MSA_EPA_combo      = msa_i,
    metro4             = metro_i,
    EPA_region         = EPA_i,
    CI_no2_2024        = CI_no2,
    abs_no2_2024       = abs_no2,
    ratio_no2_2024     = ratio_no2,
    n_tract_total_2024 = n_total,
    n_tract_valid_2024 = n_valid,
    stringsAsFactors   = FALSE
  )
  
  est <- rbind(est, est_tmp)
}
rm(no2_2024, db_2024, final4_2024, final4_metro, final_clean_filtered, ecdf_function, final_dec)
gc()

# Calculating CI in 2019
no2_2019 <- subset(no2_tract_df, year == 2019)
db_2019  <- subset(final3,       year == 2019 & cause_name == list_cause[j])

final4_2019 <- merge(no2_2019, db_2019, by = "GEOID")
final4_2019 <- merge(final4_2019, xxxx, by = "GEOID", all.x = TRUE)

final4_2019 <- final4_2019 %>%
  mutate(MSA_EPA_combo = paste(MSA_title, EPA_region, sep = " | "))

MSA_list_2019 <- unique(final4_2019$MSA_EPA_combo)

for (ii in seq_along(MSA_list_2019)) {
  
  msa_i <- MSA_list_2019[ii]
  final4_metro <- subset(final4_2019, MSA_EPA_combo == msa_i)
  
  n_total <- nrow(final4_metro)
  
  final_clean_filtered <- final4_metro %>%
    filter(!is.na(pct_white), !is.na(surface_no2))
  
  n_valid <- nrow(final_clean_filtered)
  
  # Apply 5-tract minimum
  if (n_valid < 5) next
  
  # ---- CI-NO2  ----
  ecdf_function <- stats::ecdf(final_clean_filtered[["pct_white"]])
  final_clean_filtered <- final_clean_filtered[order(final_clean_filtered[["pct_white"]]), ]
  final_clean_filtered$x_pct <- ecdf_function(final_clean_filtered[["pct_white"]]) * 100
  final_clean_filtered <- final_clean_filtered[order(final_clean_filtered$x_pct), ]
  final_clean_filtered$cumsum_y <- cumsum(final_clean_filtered$surface_no2) / sum(final_clean_filtered$surface_no2) * 100
  
  auc <- pracma::trapz(final_clean_filtered$x_pct, final_clean_filtered$cumsum_y)
  CI_no2  <- (5000 - auc) / 5000
  
  # ---- Absolute & Ratio NO2 (deciles) ----
  final_dec <- final_clean_filtered %>%
    mutate(white_decile = ntile(pct_white, 10)) %>%
    group_by(white_decile) %>%
    summarise(mean_no2 = mean(surface_no2, na.rm = TRUE), .groups = "drop") %>%
    right_join(data.frame(white_decile = 1:10), by = "white_decile") %>%
    arrange(white_decile)
  
  abs_no2   <- unname(-final_dec$mean_no2[1] + final_dec$mean_no2[10])
  ratio_no2 <- unname(final_dec$mean_no2[1] / final_dec$mean_no2[10])
  
  est_tmp <- data.frame(
    MSA_EPA_combo      = msa_i,
    CI_no2_2019        = CI_no2,
    abs_no2_2019       = abs_no2,
    ratio_no2_2019     = ratio_no2,
    n_tract_total_2019 = n_total,
    n_tract_valid_2019 = n_valid,
    stringsAsFactors   = FALSE
  )
  
  est1 <- rbind(est1, est_tmp)
}
rm(no2_2019, db_2019, final4_2019, final4_metro, final_clean_filtered, ecdf_function, final_dec)
gc()

# Merge
final_est <- full_join(est, est1, by = "MSA_EPA_combo")
final_est$CI_dif <- final_est$CI_no2_2024 - final_est$CI_no2_2019

# Check: overall totals (how many MSAs made it through the min_n filter)
n_msa_total_2024 <- length(unique(MSA_list_2024))
n_msa_valid_2024 <- nrow(est)
n_msa_total_2019 <- length(unique(MSA_list_2019))
n_msa_valid_2019 <- nrow(est1)

# Group
agg_df <- final_est %>% group_by(metro4, EPA_region) %>% summarize(CI_2024 = mean(CI_no2_2024), CI_dif = mean(CI_dif))
agg_df <- agg_df %>%
  mutate(
    metro4 = factor(metro4,
                    levels = c("Large MSA", "Small MSA",
                               "MicroSA",
                               "Non-urban")),
    EPA_region = factor(EPA_region,
                        levels = c("Region 1",
                                   "Region 2",
                                   "Region 3",
                                   "Region 4",
                                   "Region 5",
                                   "Region 6",
                                   "Region 7",
                                   "Region 8",
                                   "Region 9",
                                   "Region 10")))

region_colors <- c(
  "#1B2F70",  # Region 1 - darker blue
  "#2C4FD7",  # Region 2 - deep blue
  "#5A8AC6",  # Region 3 - slightly blue
  "#8BC8C8",  # Region 4 - teal-blue
  "#A1D99B",  # Region 5 - greenish
  "#D9C89B",  # Region 6 - yellow-green
  "#DDAA66",  # Region 7 - warm tan
  "#F2A6A0",  # Region 8 - rose
  "#7A5AC9",  # Region 9 - purple-rose
  "#4A0066"   # Region 10 - deep purple
)

# National summary from the previous block of codes
national_ci_2024 <- as.numeric(tmp_ci[26, 3])
national_ci_change <- as.numeric(tmp_ci[26, 3]-tmp_ci[1, 3])

pdf(file = paste0("/Users/k_sy_n_imac/Library/CloudStorage/Dropbox/GWU/Lab/FINESST/Research/Results/Quadrant_CI_final.pdf"), 
    width = 10, 
    height = 10)
par(family = "Lato")
ggplot(agg_df, aes(x = CI_dif, y = CI_2024,
                   color = EPA_region, shape = metro4)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.7) +   # threshold
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.7) +         # no-change line
  geom_point(size = 5, stroke=1.5, alpha = 0.8) +
  scale_color_manual(values = region_colors) +
  scale_shape_manual(
    values = c(
      "Large MSA" = 16,           # circle
      "Small MSA" = 17,       # triangle
      "MicroSA" = 15,       # square
      "Non-urban" = 3     # reverse triangle
    )) +
  geom_point(x = national_ci_change, y = national_ci_2024, color = "#B11226", size = 6, shape = 18, inherit.aes = FALSE) +
  geom_text(
    x = national_ci_change,
    y = national_ci_2024,
    label = "National Summary",
    color = "#B11226",
    size = 5,
    vjust = -1.5,   # move text slightly above the point
    hjust = 0.5,
    fontface = "bold",
    inherit.aes = FALSE) + 
  theme_minimal() +
  coord_cartesian(
    xlim = range(c(agg_df$CI_dif,  national_ci_change), na.rm = TRUE),
    ylim = range(c(agg_df$CI_2024, national_ci_2024),   na.rm = TRUE)
  ) +
  labs(
    #  title = bquote("Mean % change vs. 2024 CI by EPA Region and Urbanity"),
    #  subtitle = "Negative values indicate higher NO2 in marginalized tracts",
    x = "Changes in concentration index (race/ethnicity) between 2019 and 2024",
    y = bquote("Concentration index in 2024"),
    color = "EPA Region",
    shape = "Urbanicity"
  ) +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 13),
    axis.text = element_text(size = 13),
    axis.title.x = element_text(size = 15, face = "bold", margin = margin(t = 15, b = 15)),
    axis.title.y = element_text(size = 15, face = "bold", margin = margin(r = 15, l = 15)),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 13)
  )
dev.off()




