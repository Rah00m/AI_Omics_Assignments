# -------------------------------------------------------------------
# Author: Rahma Fathy
# Purpose: Clean and prepare patient_info.csv dataset for analysis
# -------------------------------------------------------------------

# Clear workspace
rm(list = ls())

# Load the data
patient_data <- read.csv("C:/Users/LENOVO/Documents/AI_Omics_Internship_2025/data/patient_info.csv", 
                         stringsAsFactors = FALSE)

# Convert appropriate variables to correct data types
patient_data$age <- as.numeric(patient_data$age)
patient_data$bmi <- as.numeric(patient_data$bmi)

patient_data$gender <- as.factor(patient_data$gender)
patient_data$smoker <- as.factor(patient_data$smoker)
patient_data$diagnosis <- as.factor(patient_data$diagnosis)

# Create a new binary variable 
patient_data$smoker_binary <- ifelse(patient_data$smoker == "Yes", 1, 0)

# Reorder columns 
patient_data <- patient_data[, c("patient_id", "age", "gender", "bmi", 
                                 "smoker", "smoker_binary", "diagnosis")]

#  Check the structure and summary to confirm changes
str(patient_data)
summary(patient_data)

#  Save cleaned data 
dir.create("C:/Users/LENOVO/Documents/AI_Omics_Internship_2025/clean_data", showWarnings = FALSE)
write.csv(patient_data, 
          "C:/Users/LENOVO/Documents/AI_Omics_Internship_2025/clean_data/patient_info_clean.csv", 
          row.names = FALSE)

save.image(file = "Rahma_Class_Ib_Assignment.RData")
