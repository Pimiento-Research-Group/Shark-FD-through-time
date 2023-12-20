###############################################################################################################
# 01. Data treatment
## This R code prepares the shark trait dataset for the analyses performed in all other codes
## it produces one Rdata file to be loaded into all subsequent codes
###############################################################################################################

## Import packages
library(readxl)
library(tidyverse)

# Import data
data <- read_xlsx("~/Data S1.xlsx")
attach(data)

# Prepare epoch non-matches and filter
data <- data %>% 
  mutate(Epoch = Epoch_earliest) %>%
  rowwise() %>%
  mutate(Epoch = replace(Epoch, !Epoch_earliest %in% Epoch_latest, "nonmatch"))

# Filter epochs for seven consistent time bins (only 2 teeth dated to "Holocene" and both taxa included in rest of data)
data <- data %>% 
  filter(Epoch!="Holocene")

# Factorise dental characters
data$Epoch <- ordered(data$Epoch, levels=c("Palaeocene","Eocene","Oligocene","Miocene",
                                           "Pliocene","Pleistocene","Recent","nonmatch"))
data$CH <- ordered(data$CH, levels=c("<5mm","5-20mm","20-50mm",">50mm"))
data$CW <- ordered(data$CW, levels=c("<10mm","10-35mm",">35mm"))
data$CE <- as.factor(data$CE)
data$LC <- as.factor(data$LC)
data$XO <- as.factor(data$XO)
data$LO <- as.factor(data$LO)

# Save data to load into other codes
save(data,file="~/Cleaned data.RData")
