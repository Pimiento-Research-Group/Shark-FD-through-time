###############################################################################################################
# 02. Data exploration
## This R code provides exploration of the shark trait data, including taxonomic comparisons to FINS
## it produces Figure S1 and Table S3
###############################################################################################################

## Import packages
library(readxl)
library(tidyverse)
library(doBy)
library(reshape2)
library(DescTools)
library(corrr)
library(ggsci)

# Load data
load(file="~/Cleaned data.RData")

# Number of genera
Gen_data <- data %>% 
  count(Genus)      # 163 genera
# Number of taxa (total)
Taxa_data <- data %>% 
  count(Taxon_corrected)    # 590 total taxa (species & genus-level)

Sp_data <- data %>% 
  filter(Species!="sp.") %>% 
  count(Taxon_corrected)    # 507 species-level taxa

# Number of species (fossil)
Fossil_data <- data %>% 
  filter(Epoch!="Recent") %>% 
  filter(Species!="sp.") %>% 
  count(Taxon_corrected)

# Geographic occurrences of shark teeth - produces Figure S1
world <- map_data("world")
gmap <- ggplot() +
  geom_map(data = world, map = world,
           aes(long,lat,map_id=region),
           color="black",fill="lightgray",size=0.5) +
  xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_minimal() +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent"))

data_filter <- data %>% 
  filter(Epoch!="nonmatch")

dfmap <- data_filter %>% 
  filter(Latitude!="NA") %>% 
  mutate(Latitude = as.numeric(Latitude)) %>%
  mutate(Longitude = as.numeric(Longitude))

tdata<-summaryBy(dfmap~Formation+Epoch+Longitude+Latitude, data=dfmap, FUN=length)
tdata <- tdata %>% 
  rename(Teeth = Formation.length)

gmap + 
  geom_point(
    data = tdata,
    aes(Longitude,Latitude, color = Epoch, size = Teeth),
    alpha = 0.8,
    inherit.aes = FALSE) +
  scale_color_manual(values=c("#FED439FF","#709AE1FF","#8A9197FF","#D2AF81FF","#FD7446FF","#D5E4A2FF","orchid4")) +
  guides(fill = guide_legend(title = "Status")) +
  guides(fill = guide_legend(title = "Teeth")) +
  theme(plot.title = element_text(hjust = 0.5))

# Dental character correlations
## Check correlations between dental characters sharing biological information
### Numerical tooth sizes
data %>% 
  select(CH_mm,CW_mm) %>% 
  correlate()

# Dental character polychoric correlations - ordinal tooth size and nominal tooth shape
CorPolychor(data$CH,data$CW)
CorPolychor(data$XO,data$LO)

# Tooth-position-dental character correlations - produces Table S3
CorPolychor(data$Tooth_position,data$CH)
CorPolychor(data$Tooth_position,data$CW)
CorPolychor(data$Tooth_position,data$CE)
CorPolychor(data$Tooth_position,data$LC)
CorPolychor(data$Tooth_position,data$XO)
CorPolychor(data$Tooth_position,data$LO)
