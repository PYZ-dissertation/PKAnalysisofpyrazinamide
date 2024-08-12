###############
#An R code file to clean the patient concentration data, add in relevant dosing events and add in demographic,
#physical and biochemical data
#Adds:
# ~ BLQ (0 if DV >0, 1 if <LRL)
# ~ EVID (0 if non-dosing event, 1 otherwise)
# ~ HIV (0 if HIV - negative, 1 if HIV- positive)
# ~ MDV - missing concentrations (0 if non-missing DV, 1 if DV <LRL or missing)
# ~ CMT - compartment (1 if dosing event, 2 otherwise)
# ~ RATE (-2 if dosing event, 0 otherwise)
# ~ Sex (0 if male, 1 if female
# ~ Age (years)
# ~ Ethnicity
# ~ Village
# ~ District
# ~ WT - weight (kg)
# ~ CREATRES - serum creatinine levels (mg/dL)
# ~ TOTBRES - total bilirubin level (mg/dL)
# ~ ALTRES - alanine aminotransferase level (units/L)
# ~ AMT - doing event column (0 if no dose given, otherwise amount of pyrazinamide recieved (ng/))
###############

#Load packages
library(readxl)
library(lubridate)
library(dplyr)
library(openxlsx)
library(tidyr)
library(zoo)
library(stringr)

#Loading required data
tbdata <- as.data.frame(read_excel()) #Change here for own concentration data file 
dem_data <- as.data.frame(read_excel()) #Change here for own demographic data file 
cov_data <- as.data.frame(read_excel()) #Change here for own physical data file
bio_data <- as.data.frame(read_excel()) #Change here for own biochemical data file
daily_data <- as.data.frame(read_excel())#Change here for own dosing file

#Adding BLQ column
tbdata <- mutate(tbdata, BLQ = ifelse(DV == "<LRL", 1, 0))
# Convert TIMEC to POSIXct format
tbdata$TIMEC <- as.POSIXct(tbdata$TIMEC, format = "%H:%M:%S")
time_only <- format(tbdata$TIMEC, "%H:%M:%S")
# Combine DATE with extracted time
new_datetime <- paste(tbdata$DATE, time_only)
# Convert the new combined datetime to POSIXct format
new_datetime <- as.POSIXct(new_datetime)
# Update TIMEC with the new datetime
tbdata$TIMEC <- new_datetime
tbdata$DV <- gsub("[^0-9]", "", tbdata$DV)

########COMBING DAILY DOSE WITH TBDATA##########################################
# Extract just the ID from daily data
daily_data$ID <- as.numeric(str_extract(daily_data$Code, "\\d+$"))

#Make the TIME column the same format as tbdata
combined_datetime <- paste(daily_data$DATE, daily_data$TIMEC, sep=" ")
combined_datetime <- as.POSIXct(combined_datetime, format="%Y-%m-%d %H:%M")
daily_data$TIMEC <- combined_datetime

#Take just the unique rows
daily_data <- unique(daily_data)

#Adding in columns to both dataframes so they both have the same types
tbdata$dtbdrug <-NA
daily_data$Code<-NULL
daily_data$'Laboratory ID'<-NA
daily_data$'Study arm'<-NA
daily_data$'Time.point'<-NA
daily_data$'DV'<-NA
daily_data$'Hemolysis'<-NA
daily_data$'BLQ'<-NA

# Convert 'DATE' columns to Date class
tbdata$DATE <- as.Date(tbdata$DATE)
daily_data$DATE <- as.Date(daily_data$DATE)

#Check which rowxs are in Daily data and not tb data
missing_rows <- daily_data %>%
  anti_join(tbdata, by = c("ID", "DATE", "TIMEC")) %>%
  arrange(ID, DATE)

#Combine the missing rows with tbdata
tbdata <- rbind(tbdata, missing_rows)
#order them by ID, date and then time
tbdata <- tbdata[order(tbdata$ID, tbdata$DATE, tbdata$TIMEC), ]
#Remove rows where pyrazinamide wasn't given as a drug
tbdata <- tbdata[!(tbdata$dtbdrug %in% c("FDC2", "FDC3")), ]
#Get rid of the type pf drug given column
tbdata$dtbdrug  <-NULL

# Function to calculate the difference in days between two dates to fill in all dosing events
date_diff <- function(date1, date2) {
  as.numeric(difftime(date1, date2, units = "days"))
}

for (i in 1:nrow(tbdata)) {
  # Check if Time.point = NA
  if (is.na(tbdata[i, "Time.point"])) {
    # Find corresponding D0 H0 date to each ID
    d0_h0_date <- subset(tbdata, ID == tbdata$ID[i] & Time.point == "D0 H0", select = DATE)
    # Check if d0_h0_date is not empty
    if (!is.null(d0_h0_date) && nrow(d0_h0_date) > 0) {
      # Extract the first row of d0_h0_date
      d0_h0_date <- d0_h0_date[1, , drop = FALSE]
      # Calculate the difference in days between the current date and D0 H0 date
      days_diff <- date_diff(tbdata$DATE[i], d0_h0_date$DATE)
      # Assign appropriate Time.point based on the conditions
      if (days_diff < 14) {
        tbdata$Time.point[i] <- "P1 H0"#these get w1
      } else if (days_diff >= 14 & days_diff < 28) {
        tbdata$Time.point[i] <- "P2 H0"#these get w2
      } else if (days_diff >= 28 & days_diff < 60) {
        tbdata$Time.point[i] <- "P3 H0"#these get w4
      } else if (days_diff >= 60 & days_diff < 90) {
        tbdata$Time.point[i] <- "P4 H0"#these get m2
      } else if (days_diff >= 90 & days_diff < 120) {
        tbdata$Time.point[i] <- "P5 H0"#these get m3
      } else if (days_diff >= 120 & days_diff < 150) {
        tbdata$Time.point[i] <- "P6 H0"#these get m4
      } else if (days_diff >= 150 & days_diff < 180) {
        tbdata$Time.point[i] <- "P7 H0"#these get m5
      } else if (days_diff >= 180) {
        tbdata$Time.point[i] <- "P8 H0"#these get m6
      }
    }
  }
}
#Set EVID =1 when there is a H0 (i.e. a dose is given)
tbdata$EVID <- ifelse(grepl("H0($|\\s)", tbdata$'Time.point'), 1, 0)#]

#Set HIV status to all rows
tbdata <- tbdata %>%
  group_by(ID) %>%
  fill(`Study arm`, .direction = "downup") %>%
  ungroup()
tbdata<-as.data.frame(tbdata)

# Convert 'DV' column to numeric
tbdata$DV <- as.numeric(tbdata$DV)

# Change the format of DATE to "YYYY-MM-DD"
tbdata$DATE <- format(tbdata$DATE, "%Y-%m-%d")

#Transform HIV covariate to binary
tbdata$HIV <- ifelse(tbdata$`Study arm` == "HIV pos", 1, 0)#1 IF POSITIVE, 0 IF NEGATIVE

#Create missing concentration column
tbdata$MDV <- ifelse(is.na(tbdata$DV), 1, 0)

#Create CMT column
tbdata$CMT <- ifelse(grepl("H0($|\\s)", tbdata$'Time.point'), 1, 2)

#Change any NAs in BLQ column to 0 
tbdata$BLQ[is.na(tbdata$BLQ)] <- 0

#Creatre a rate column
tbdata$RATE <- ifelse(grepl("H0($|\\s)", tbdata$'Time.point'), -2, 0)

# Extract subject number from Subject Label
dem_data$Subject_Number <- as.numeric(str_extract(dem_data$'Subject Label', "\\d+$"))
#Creating a merged data frame by ID and Subject Number
merged_data <- merge(tbdata, dem_data, by.x = "ID", by.y = "Subject_Number")

#################################DEMOGRAPHICS###################################
# Adding demographic data
#Sex
tbdata$SEX <- merged_data$Sex
#Age
tbdata$AGE <- merged_data$'Age (year)'
#Ethnicity
tbdata$Ethnicity <- merged_data$Ethnicity
#Village
tbdata$Village <- merged_data$Village
#District
tbdata$District <- merged_data$District

#Converting Sex to binayr
tbdata$SEX <- ifelse(tbdata$SEX == "Male", 0, 1)#1 if FEMALE, 0 IF MALE


#################################COVARIATES#####################################
#Adding covariate data, extract subject number from label
cov_data$Subject_Number <- as.numeric(str_extract(cov_data$'label', "\\d+$"))
#Creating a merged data frame by ID and Subject Number for covariate data
merged_data <- merge(tbdata, cov_data, by.x = "ID", by.y = "Subject_Number")

#WEIGHT
merged_data$WT <- NA
for (i in 1:nrow(merged_data)) {
  time_point <- merged_data[i, 'Time.point']
  id <- merged_data[i, "ID"]
  
  # Assign weight based on time point
  if (grepl("D0", time_point)) {
    d0_subset <- merged_data[merged_data$ID == id & grepl("D0", merged_data$'Time.point'), "weight"]
    merged_data[i, "WT"] <- d0_subset[1]
  } else if (grepl("W6", time_point)) {
    w6_subset <- merged_data[merged_data$ID == id & grepl("W6", merged_data$'Time.point'), "weight_w4"]
    merged_data[i, "WT"] <- w6_subset[1]
  } else if (grepl("M2", time_point)) {
    M2_subset <- merged_data[merged_data$ID == id & grepl("M2", merged_data$'Time.point'), "weight_m2"]
    merged_data[i, "WT"] <- M2_subset[1]
  } else if (grepl("M3", time_point)) {  
    M3_subset <- merged_data[merged_data$ID == id & grepl("M3", merged_data$'Time.point'), "weight_m3"]
    merged_data[i, "WT"] <- M3_subset[1]
  } else if (grepl("M4", time_point)) { 
    M4_subset <- merged_data[merged_data$ID == id & grepl("M4", merged_data$'Time.point'), "weight_m4"]
    merged_data[i, "WT"] <- M4_subset[1]
  } else if (grepl("M5", time_point)) {   
    M5_subset <- merged_data[merged_data$ID == id & grepl("M5", merged_data$'Time.point'), "weight_m5"]
    merged_data[i, "WT"] <- M5_subset[1]
  } else if (grepl("M6", time_point)) {     
    M6_subset <- merged_data[merged_data$ID == id & grepl("M6", merged_data$'Time.point'), "weight_m6"]
    merged_data[i, "WT"] <- M6_subset[1]
  } else if (grepl("P1", time_point)) {     
    M6_subset <- merged_data[merged_data$ID == id & grepl("P1", merged_data$'Time.point'), "weight"]
    merged_data[i, "WT"] <- M6_subset[1]
  } else if (grepl("P2", time_point)) {     
    M6_subset <- merged_data[merged_data$ID == id & grepl("P2", merged_data$'Time.point'), "weight_w2"]
    merged_data[i, "WT"] <- M6_subset[1]
  } else if (grepl("P3", time_point)) {     
    M6_subset <- merged_data[merged_data$ID == id & grepl("P3", merged_data$'Time.point'), "weight_w4"]
    merged_data[i, "WT"] <- M6_subset[1]
  } else if (grepl("P4", time_point)) {     
    M6_subset <- merged_data[merged_data$ID == id & grepl("P4", merged_data$'Time.point'), "weight_m2"]
    merged_data[i, "WT"] <- M6_subset[1]
  } else if (grepl("P5", time_point)) {     
    M6_subset <- merged_data[merged_data$ID == id & grepl("P5", merged_data$'Time.point'), "weight_m3"]
    merged_data[i, "WT"] <- M6_subset[1]
  } else if (grepl("P6", time_point)) {     
    M6_subset <- merged_data[merged_data$ID == id & grepl("P6", merged_data$'Time.point'), "weight_m4"]
    merged_data[i, "WT"] <- M6_subset[1]
  } else if (grepl("P7", time_point)) {     
    M6_subset <- merged_data[merged_data$ID == id & grepl("P7", merged_data$'Time.point'), "weight_m5"]
    merged_data[i, "WT"] <- M6_subset[1]
  } else if (grepl("P8", time_point)) {     
    M6_subset <- merged_data[merged_data$ID == id & grepl("P8", merged_data$'Time.point'), "weight_m6"]
    merged_data[i, "WT"] <- M6_subset[1]
  }
  
  }
# Assign WT back to tbdata
tbdata$WT <- merged_data$WT

#BIOCHEMISTRY####################################################################
bio_data$Subject_Number <- as.numeric(str_extract(bio_data$'label', "\\d+$"))
#Creating a merged data frame by ID and Subject Number for covariate data
merged_data <- merge(tbdata, bio_data, by.x = "ID", by.y = "Subject_Number")

#CREATRES
merged_data$CREATRES<- NA
for (i in 1:nrow(merged_data)) {
  time_point <- merged_data[i, 'Time.point']
  id <- merged_data[i, "ID"]
  
  # Assign creatres based on time point
  if (grepl("D0", time_point)) {
    d0_subset <- merged_data[merged_data$ID == id & grepl("D0", merged_data$'Time.point'), "creatres"]
    merged_data[i, "CREATRES"] <- d0_subset[1]
  } else if (grepl("W6", time_point)) {
    w6_subset <- merged_data[merged_data$ID == id & grepl("W6", merged_data$'Time.point'), "creatres_w4"]
    merged_data[i, "CREATRES"] <- w6_subset[1]
  } else if (grepl("M2", time_point)) {
    M2_subset <- merged_data[merged_data$ID == id & grepl("M2", merged_data$'Time.point'), "creatres_m2"]
    merged_data[i, "CREATRES"] <- M2_subset[1]
  } else if (grepl("M3", time_point)) {  
    M3_subset <- merged_data[merged_data$ID == id & grepl("M3", merged_data$'Time.point'), "creatres_m2"]
    merged_data[i, "CREATRES"] <- M3_subset[1]
  } else if (grepl("M4", time_point)) { 
    M4_subset <- merged_data[merged_data$ID == id & grepl("M4", merged_data$'Time.point'), "creatres_m2"]
    merged_data[i, "CREATRES"] <- M4_subset[1]
  } else if (grepl("M5", time_point)) {   
    M5_subset <- merged_data[merged_data$ID == id & grepl("M5", merged_data$'Time.point'), "creatres_m2"]
    merged_data[i, "CREATRES"] <- M5_subset[1]
  } else if (grepl("M6", time_point)) {     
    M6_subset <- merged_data[merged_data$ID == id & grepl("M6", merged_data$'Time.point'), "creatres_m2"]
    merged_data[i, "CREATRES"] <- M6_subset[1]
  } else if (grepl("P1", time_point)) {     
    M6_subset <- merged_data[merged_data$ID == id & grepl("P1", merged_data$'Time.point'), "creatres"]
    merged_data[i, "CREATRES"] <- M6_subset[1]
  } else if (grepl("P2", time_point)) {     
    M6_subset <- merged_data[merged_data$ID == id & grepl("P2", merged_data$'Time.point'), "creatres_w2"]
    merged_data[i, "CREATRES"] <- M6_subset[1]
  } else if (grepl("P3", time_point)) {     
    M6_subset <- merged_data[merged_data$ID == id & grepl("P3", merged_data$'Time.point'), "creatres_w4"]
    merged_data[i, "CREATRES"] <- M6_subset[1]
  } else if (grepl("P4", time_point)) {     
    M6_subset <- merged_data[merged_data$ID == id & grepl("P4", merged_data$'Time.point'), "creatres_m2"]
    merged_data[i, "CREATRES"] <- M6_subset[1]
  } else if (grepl("P5", time_point)) {     
    M6_subset <- merged_data[merged_data$ID == id & grepl("P5", merged_data$'Time.point'), "creatres_m2"]
    merged_data[i, "CREATRES"] <- M6_subset[1]
  } else if (grepl("P6", time_point)) {     
    M6_subset <- merged_data[merged_data$ID == id & grepl("P6", merged_data$'Time.point'), "creatres_m2"]
    merged_data[i, "CREATRES"] <- M6_subset[1]
  } else if (grepl("P7", time_point)) {     
    M6_subset <- merged_data[merged_data$ID == id & grepl("P7", merged_data$'Time.point'), "creatres_m2"]
    merged_data[i, "CREATRES"] <- M6_subset[1]
  } else if (grepl("P8", time_point)) {     
    M6_subset <- merged_data[merged_data$ID == id & grepl("P8", merged_data$'Time.point'), "creatres_m2"]
    merged_data[i, "CREATRES"] <- M6_subset[1]
  }
  
}
# Assign CREATRES back to tbdata
tbdata$CREATRES <- merged_data$CREATRES

#ALTRES
merged_data$ALTRES<- NA
for (i in 1:nrow(merged_data)) {
  time_point <- merged_data[i, 'Time.point']
  id <- merged_data[i, "ID"]
  
  # Assign altres based on time point
  if (grepl("D0", time_point)) {
    d0_subset <- merged_data[merged_data$ID == id & grepl("D0", merged_data$'Time.point'), "altres"]
    merged_data[i, "ALTRES"] <- d0_subset[1]
  } else if (grepl("W6", time_point)) {
    w6_subset <- merged_data[merged_data$ID == id & grepl("W6", merged_data$'Time.point'), "altres_w4"]
    merged_data[i, "ALTRES"] <- w6_subset[1]
  } else if (grepl("M2", time_point)) {
    M2_subset <- merged_data[merged_data$ID == id & grepl("M2", merged_data$'Time.point'), "altres_m2"]
    merged_data[i, "ALTRES"] <- M2_subset[1]
  } else if (grepl("M3", time_point)) {  
    M3_subset <- merged_data[merged_data$ID == id & grepl("M3", merged_data$'Time.point'), "altres_m2"]
    merged_data[i, "ALTRES"] <- M3_subset[1]
  } else if (grepl("M4", time_point)) { 
    M4_subset <- merged_data[merged_data$ID == id & grepl("M4", merged_data$'Time.point'), "altres_m2"]
    merged_data[i, "ALTRES"] <- M4_subset[1]
  } else if (grepl("M5", time_point)) {   
    M5_subset <- merged_data[merged_data$ID == id & grepl("M5", merged_data$'Time.point'), "altres_m2"]
    merged_data[i, "ALTRES"] <- M5_subset[1]
  } else if (grepl("M6", time_point)) {     
    M6_subset <- merged_data[merged_data$ID == id & grepl("M6", merged_data$'Time.point'), "altres_m2"]
    merged_data[i, "ALTRES"] <- M6_subset[1]
  }else if (grepl("P1", time_point)) {     
    M6_subset <- merged_data[merged_data$ID == id & grepl("P1", merged_data$'Time.point'), "altres"]
    merged_data[i, "ALTRES"] <- M6_subset[1]
  } else if (grepl("P2", time_point)) {     
    M6_subset <- merged_data[merged_data$ID == id & grepl("P2", merged_data$'Time.point'), "altres_w2"]
    merged_data[i, "ALTRES"] <- M6_subset[1]
  } else if (grepl("P3", time_point)) {     
    M6_subset <- merged_data[merged_data$ID == id & grepl("P3", merged_data$'Time.point'), "altres_w4"]
    merged_data[i, "ALTRES"] <- M6_subset[1]
  } else if (grepl("P4", time_point)) {     
    M6_subset <- merged_data[merged_data$ID == id & grepl("P4", merged_data$'Time.point'), "altres_m2"]
    merged_data[i, "ALTRES"] <- M6_subset[1]
  } else if (grepl("P5", time_point)) {     
    M6_subset <- merged_data[merged_data$ID == id & grepl("P5", merged_data$'Time.point'), "altres_m2"]
    merged_data[i, "ALTRES"] <- M6_subset[1]
  } else if (grepl("P6", time_point)) {     
    M6_subset <- merged_data[merged_data$ID == id & grepl("P6", merged_data$'Time.point'), "altres_m2"]
    merged_data[i, "ALTRES"] <- M6_subset[1]
  } else if (grepl("P7", time_point)) {     
    M6_subset <- merged_data[merged_data$ID == id & grepl("P7", merged_data$'Time.point'), "altres_m2"]
    merged_data[i, "ALTRES"] <- M6_subset[1]
  } else if (grepl("P8", time_point)) {     
    M6_subset <- merged_data[merged_data$ID == id & grepl("P8", merged_data$'Time.point'), "altres_m2"]
    merged_data[i, "ALTRES"] <- M6_subset[1]
  }
  
}
# Assign ALTRES back to tbdata
tbdata$ALTRES <- merged_data$ALTRES

#TOTBRES
merged_data$TOTBRES<- NA
for (i in 1:nrow(merged_data)) {
  time_point <- merged_data[i, 'Time.point']
  id <- merged_data[i, "ID"]
  
  # Assign totbres based on time point
  if (grepl("D0", time_point)) {
    d0_subset <- merged_data[merged_data$ID == id & grepl("D0", merged_data$'Time.point'), "totbres"]
    merged_data[i, "TOTBRES"] <- d0_subset[1]
  } else if (grepl("W6", time_point)) {
    w6_subset <- merged_data[merged_data$ID == id & grepl("W6", merged_data$'Time.point'), "totbres_w4"]
    merged_data[i, "TOTBRES"] <- w6_subset[1]
  } else if (grepl("M2", time_point)) {
    M2_subset <- merged_data[merged_data$ID == id & grepl("M2", merged_data$'Time.point'), "totbres_m2"]
    merged_data[i, "TOTBRES"] <- M2_subset[1]
  } else if (grepl("M3", time_point)) {  
    M3_subset <- merged_data[merged_data$ID == id & grepl("M3", merged_data$'Time.point'), "totbres_m2"]
    merged_data[i, "TOTBRES"] <- M3_subset[1]
  } else if (grepl("M4", time_point)) { 
    M4_subset <- merged_data[merged_data$ID == id & grepl("M4", merged_data$'Time.point'), "totbres_m2"]
    merged_data[i, "TOTBRES"] <- M4_subset[1]
  } else if (grepl("M5", time_point)) {   
    M5_subset <- merged_data[merged_data$ID == id & grepl("M5", merged_data$'Time.point'), "totbres_m2"]
    merged_data[i, "TOTBRES"] <- M5_subset[1]
  } else if (grepl("M6", time_point)) {     
    M6_subset <- merged_data[merged_data$ID == id & grepl("M6", merged_data$'Time.point'), "totbres_m2"]
    merged_data[i, "TOTBRES"] <- M6_subset[1]
  }else if (grepl("P1", time_point)) {     
    M6_subset <- merged_data[merged_data$ID == id & grepl("P1", merged_data$'Time.point'), "totbres"]
    merged_data[i, "TOTBRES"] <- M6_subset[1]
  } else if (grepl("P2", time_point)) {     
    M6_subset <- merged_data[merged_data$ID == id & grepl("P2", merged_data$'Time.point'), "totbres_w2"]
    merged_data[i, "TOTBRES"] <- M6_subset[1]
  } else if (grepl("P3", time_point)) {     
    M6_subset <- merged_data[merged_data$ID == id & grepl("P3", merged_data$'Time.point'), "totbres_w4"]
    merged_data[i, "TOTBRES"] <- M6_subset[1]
  } else if (grepl("P4", time_point)) {     
    M6_subset <- merged_data[merged_data$ID == id & grepl("P4", merged_data$'Time.point'), "totbres_m2"]
    merged_data[i, "TOTBRES"] <- M6_subset[1]
  } else if (grepl("P5", time_point)) {     
    M6_subset <- merged_data[merged_data$ID == id & grepl("P5", merged_data$'Time.point'), "totbres_m2"]
    merged_data[i, "TOTBRES"] <- M6_subset[1]
  } else if (grepl("P6", time_point)) {     
    M6_subset <- merged_data[merged_data$ID == id & grepl("P6", merged_data$'Time.point'), "totbres_m2"]
    merged_data[i, "TOTBRES"] <- M6_subset[1]
  } else if (grepl("P7", time_point)) {     
    M6_subset <- merged_data[merged_data$ID == id & grepl("P7", merged_data$'Time.point'), "totbres_m2"]
    merged_data[i, "TOTBRES"] <- M6_subset[1]
  } else if (grepl("P8", time_point)) {     
    M6_subset <- merged_data[merged_data$ID == id & grepl("P8", merged_data$'Time.point'), "totbres_m2"]
    merged_data[i, "TOTBRES"] <- M6_subset[1]
  }
  
}
# Assign TOTBRES back to tbdata
tbdata$TOTBRES <- merged_data$TOTBRES

#Ensuring nothing is NA and any missing covariates are filled in from last known point 
# Specify the columns you want to fill missing values for
columns_to_fill <- c( "SEX", "AGE", "Ethnicity", "Village", "District", "WT", "CREATRES", "ALTRES", "TOTBRES")

# Group by ID and fill missing values in specified columns
tbdata[columns_to_fill] <- ave(tbdata[columns_to_fill], tbdata$ID, FUN = function(x) na.locf(x, na.rm = FALSE))

#Filling in missing TOTBRES data with NS
tbdata <- tbdata %>%
  group_by(ID) %>%
  mutate(TOTBRES = ifelse(is.na(TOTBRES), zoo::na.locf(TOTBRES, fromLast = TRUE), TOTBRES))


#ADDING IN DOSING ADJUSTED FOR WEIGHT###########################################
#400mg per tablet = 4E+08 ng
tbdata <- tbdata %>%
  mutate(AMT = case_when(
    grepl("H0($|\\s)", Time.point) & WT >= 21 & WT < 30 ~ 2*4E+08,
    grepl("H0($|\\s)", Time.point) & WT >= 30 & WT < 35 ~ 2*4E+08,
    grepl("H0($|\\s)", Time.point) & WT >= 35 & WT < 40 ~ 2.5*4E+08,
    grepl("H0($|\\s)", Time.point) & WT >= 40 & WT < 50 ~ 3*4E+08,
    grepl("H0($|\\s)", Time.point) & WT >= 50 & WT < 55 ~ 3*4E+08,
    grepl("H0($|\\s)", Time.point) & WT >= 55 & WT < 65 ~ 4*4E+08,
    grepl("H0($|\\s)", Time.point) & WT >= 65 & WT <= 70 ~ 4*4E+08,
    grepl("H0($|\\s)", Time.point) & WT > 70 ~ 5*4E+08,
    TRUE ~ 0
  ))

# Convert TIMEC to POSIXct format
tbdata$TIMEC <- as.POSIXct(tbdata$TIMEC, format = "%Y-%m-%d %H:%M:%S")

# Extract only the time part from the POSIXct object
tbdata$TIMEC <- format(tbdata$TIMEC, format = "%H:%M")
View(tbdata)
#Writing an excel file
write.xlsx(tbdata, "PYZfulldata.xlsx")
