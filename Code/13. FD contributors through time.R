###################################################################################################################
# 13. FD contributors through time 
## This R code provides analyses on species contribution to functional diversity through time
# it produces the highest FOri & FSpe taxa per epoch; included in Figures 1 and S4
###################################################################################################################

## Import packages
library(tidyverse)
library(tibble)
library(scales)
library(cowplot)
library(doBy)
library(reshape2)
library(RColorBrewer)
library(mFD)
library(ggsci)
library(data.table)

# Load data
load(file="C:/Users/2022207/Dropbox/Jack's PhD/Chapter 2. FD changes over time/Analyses/Current Analyses/R codes/Cleaned data.RData")

# Construct per-epoch functional space and highlight highest FOri & FSpe contributors
# Occurrence matrix
Occ_data <- data %>%
  group_by(Taxon_corrected) %>% 
  distinct(Epoch,Epoch_earliest,Epoch_latest,Current_status)
Occ_data <- Occ_data[order(Occ_data$Taxon_corrected),]
Occ_data$Taxa_dup<-paste(Occ_data$Taxon_corrected, Occ_data$Epoch, sep="+")

# Cast into wide-format data
wp.abun <- Occ_data %>% 
  pivot_wider(names_from = Epoch, values_from = Taxa_dup)

wp.abun <- wp.abun[, c("Taxon_corrected","Epoch_earliest","Epoch_latest","Current_status",
                       "Palaeocene","Eocene","Oligocene","Miocene","Pliocene","Pleistocene","Recent","nonmatch")]


# Mark 0 and 1 for absence and presence data in all epochs
wp.abun$Palaeocene[is.na(wp.abun$Palaeocene)] <- "0" 
wp.abun$Eocene[is.na(wp.abun$Eocene)] <- "0"
wp.abun$Oligocene[is.na(wp.abun$Oligocene)] <- "0"
wp.abun$Miocene[is.na(wp.abun$Miocene)] <- "0"
wp.abun$Pliocene[is.na(wp.abun$Pliocene)] <- "0"
wp.abun$Pleistocene[is.na(wp.abun$Pleistocene)] <- "0"
wp.abun$Recent[is.na(wp.abun$Recent)] <- "0"
wp.abun$nonmatch[is.na(wp.abun$nonmatch)] <- "0"

wp.abun <- wp.abun %>% 
  mutate(Palaeocene = replace(Palaeocene,Palaeocene!="0", "1")) %>%
  mutate(Eocene = replace(Eocene,Eocene!="0", "1")) %>% 
  mutate(Oligocene = replace(Oligocene,Oligocene!="0", "1")) %>% 
  mutate(Miocene = replace(Miocene,Miocene!="0", "1")) %>%  
  mutate(Pliocene = replace(Pliocene,Pliocene!="0", "1")) %>%
  mutate(Pleistocene = replace(Pleistocene,Pleistocene!="0", "1")) %>% 
  mutate(Recent = replace(Recent,Recent!="0", "1")) %>%
  mutate(nonmatch = replace(nonmatch,nonmatch!="0", "1"))

## Fill in non-matches
wp.abun <- wp.abun %>% 
  rowwise() %>% 
  mutate(Palaeocene = case_when(Epoch_earliest=="Palaeocene" && nonmatch=="1"~"1",TRUE~Palaeocene)) %>% 
  mutate(Eocene = case_when(Epoch_latest=="Eocene" && nonmatch=="1"~"1",TRUE~Eocene)) %>% 
  mutate(Eocene = case_when(Epoch_earliest=="Eocene" && nonmatch=="1"~"1",TRUE~Eocene)) %>%
  mutate(Oligocene = case_when(Epoch_latest=="Oligocene" && nonmatch=="1"~"1",TRUE~Oligocene)) %>% 
  mutate(Oligocene = case_when(Epoch_earliest=="Oligocene" && nonmatch=="1"~"1",TRUE~Oligocene)) %>% 
  mutate(Miocene = case_when(Epoch_latest=="Miocene" && nonmatch=="1"~"1",TRUE~Miocene)) %>% 
  mutate(Miocene = case_when(Epoch_earliest=="Miocene" && nonmatch=="1"~"1",TRUE~Miocene)) %>% 
  mutate(Pliocene = case_when(Epoch_latest=="Pliocene" && nonmatch=="1"~"1",TRUE~Pliocene))

# Range through via "extant" status
# Automatically range-through extant taxa
wp.abun <- wp.abun %>% 
  rowwise() %>% 
  mutate(Eocene = case_when(Palaeocene == "1" && Current_status == "Extant"~"1", TRUE ~ Eocene)) %>% 
  mutate(Oligocene = case_when(Eocene == "1" && Current_status == "Extant"~"1", TRUE ~ Oligocene)) %>% 
  mutate(Miocene = case_when(Oligocene == "1" && Current_status == "Extant"~"1", TRUE ~ Miocene)) %>% 
  mutate(Pliocene = case_when(Miocene == "1" && Current_status == "Extant"~"1", TRUE ~ Pliocene)) %>% 
  mutate(Pleistocene = case_when(Pliocene == "1" && Current_status == "Extant"~"1", TRUE ~ Pleistocene)) %>% 
  mutate(Recent = case_when(Pleistocene == "1" && Current_status == "Extant"~"1", TRUE ~ Recent))

baskets.range <- wp.abun %>% 
  select(-c(Epoch_earliest,Epoch_latest,nonmatch))

# Remove all genus-level taxa from Recent where all species are accounted for (i.e., make 0)
baskets.range <- baskets.range %>% 
  rowwise() %>% 
  mutate(Recent = case_when(Recent == "1" && Taxon_corrected == "Alopias sp."~"0", TRUE ~ Recent)) %>%
  mutate(Recent = case_when(Recent == "1" && Taxon_corrected == "Carcharias sp."~"0", TRUE ~ Recent)) %>% 
  mutate(Recent = case_when(Recent == "1" && Taxon_corrected == "Carcharodon sp."~"0", TRUE ~ Recent)) %>% 
  mutate(Recent = case_when(Recent == "1" && Taxon_corrected == "Dalatias sp."~"0", TRUE ~ Recent)) %>% 
  mutate(Recent = case_when(Recent == "1" && Taxon_corrected == "Galeocerdo sp."~"0", TRUE ~ Recent)) %>% 
  mutate(Recent = case_when(Recent == "1" && Taxon_corrected == "Galeorhinus sp."~"0", TRUE ~ Recent)) %>% 
  mutate(Recent = case_when(Recent == "1" && Taxon_corrected == "Galeus sp."~"0", TRUE ~ Recent)) %>% 
  mutate(Recent = case_when(Recent == "1" && Taxon_corrected == "Hemipristis sp."~"0", TRUE ~ Recent)) %>% 
  mutate(Recent = case_when(Recent == "1" && Taxon_corrected == "Heptranchias sp."~"0", TRUE ~ Recent)) %>% 
  mutate(Recent = case_when(Recent == "1" && Taxon_corrected == "Isurus sp."~"0", TRUE ~ Recent)) %>%
  mutate(Recent = case_when(Recent == "1" && Taxon_corrected == "Megachasma sp."~"0", TRUE ~ Recent)) %>% 
  mutate(Recent = case_when(Recent == "1" && Taxon_corrected == "Nebrius sp."~"0", TRUE ~ Recent)) %>% 
  mutate(Recent = case_when(Recent == "1" && Taxon_corrected == "Notorynchus sp."~"0", TRUE ~ Recent))

# Make 0 NAs, merge duplicates, reform NAs as 0
baskets.range[baskets.range == 0] <- NA
coalesce_all_columns <- function(df) {
  return(coalesce(!!! as.list(df)))
}

baskets.range <- baskets.range %>%
  group_by(Taxon_corrected) %>%
  summarise_all(coalesce_all_columns)

baskets.range$Palaeocene[is.na(baskets.range$Palaeocene)] <- "0" 
baskets.range$Eocene[is.na(baskets.range$Eocene)] <- "0"
baskets.range$Oligocene[is.na(baskets.range$Oligocene)] <- "0"
baskets.range$Miocene[is.na(baskets.range$Miocene)] <- "0"
baskets.range$Pliocene[is.na(baskets.range$Pliocene)] <- "0"
baskets.range$Pleistocene[is.na(baskets.range$Pleistocene)] <- "0"
baskets.range$Recent[is.na(baskets.range$Recent)] <- "0"

# Fill in any gaps in range automatically
baskets.range.taxon <- baskets.range %>% 
  rowwise() %>% 
  mutate(Eocene = case_when(Palaeocene == "1" && Miocene == "1"~ "1", TRUE ~ Eocene)) %>%
  mutate(Oligocene = case_when(Palaeocene == "1" && Miocene == "1"~ "1", TRUE ~ Oligocene)) %>%
  mutate(Oligocene = case_when(Eocene == "1" && Miocene == "1"~ "1", TRUE ~ Oligocene)) %>%
  mutate(Miocene = case_when(Oligocene == "1" && Pliocene == "1"~ "1", TRUE ~ Miocene)) %>%
  mutate(Pliocene = case_when(Miocene == "1" && Pleistocene == "1"~ "1", TRUE ~ Pliocene)) %>%
  mutate(Oligocene = case_when(Eocene == "1" && Pliocene == "1"~ "1", TRUE ~ Oligocene)) %>%
  mutate(Miocene = case_when(Eocene == "1" && Pliocene == "1"~ "1", TRUE ~ Miocene))

# Mark and remove species that only occur in Recent
selected_taxa <- baskets.range %>%
  filter(
    Recent == "1" &
      Palaeocene == "0" &
      Eocene == "0" &
      Oligocene == "0" &
      Miocene == "0" &
      Pliocene == "0" &
      Pleistocene == "0"
  ) %>%
  select(Taxon_corrected)
# Filter taxa that are NOT present in "selected_taxa"
filtered_baskets.range <- baskets.range %>%
  anti_join(selected_taxa, by = "Taxon_corrected")

# Make taxa row names & tidy dataframe to form final matrix
baskets <- filtered_baskets.range %>% 
  remove_rownames %>% 
  column_to_rownames(var="Taxon_corrected") %>% 
  select(-Current_status)

# Transpose to form final occurrence matrix
baskets_sharks_weights <- t(baskets)

# Function to calculate mode
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

# Calculate average CH & mode CE, LC, XO, LO per species & genus
Av_data <- data %>%
  group_by(Taxon_corrected) %>% 
  reframe(CH = mean(CH_mm),
          CW = mean(CW_mm),
          CE = getmode(CE),
          LC = getmode(LC),
          XO = getmode(XO),
          LO = getmode(LO))

Av_sharks <- Av_data %>% 
  rowwise() %>% 
  mutate(CH = cut(CH,
                  breaks = c(-Inf,5,20,50,Inf),
                  labels = c("<5mm","5-20mm","20-50mm",">50mm"),
                  right = FALSE)) %>% 
  mutate(CW = cut(CW,
                  breaks = c(-Inf,10,35,Inf),
                  labels = c("<10mm","10-35mm",">35mm"),
                  right = FALSE))
Av_sharks$CH <- ordered(Av_sharks$CH, levels=c("<5mm","5-20mm","20-50mm",">50mm"))
Av_sharks$CW <- ordered(Av_sharks$CW, levels=c("<10mm","10-35mm",">35mm"))

# Filter taxa that are NOT present in "selected_taxa"
Av_sharks_filtered <- Av_sharks %>%
  anti_join(selected_taxa, by = "Taxon_corrected")

# Make species-trait matrix
sharks_traits <- Av_sharks_filtered %>%
  remove_rownames %>% 
  column_to_rownames(var="Taxon_corrected")

# Form trait category matrix
sharks_traits_cat <- tibble(trait_name = c("CH", "CW", "CE", "LC", "XO", "LO"),
                            trait_type = c("O", "O", "N", "N", "N", "N"),
                            trait_weight = c(0.5,0.5,1,1,0.33,0.67))

# Replicate FD analyses
# make assemblage data into a matrix & ensure all matrix values are numeric
baskets_sharks_weights <- data.matrix(baskets_sharks_weights, rownames.force = NA)
class(baskets_sharks_weights) <- "numeric"

# summarise dataset
asb_sp_sharks_summ <- mFD::asb.sp.summary(asb_sp_w = baskets_sharks_weights)

# retrieve species occurrences for all assemblages (epochs):
asb_sp_sharks_occ <- asb_sp_sharks_summ$asb_sp_occ

# Form the functional entities
sp_to_fe_sharks <- mFD::sp.to.fe(
  sp_tr       = sharks_traits, 
  tr_cat      = sharks_traits_cat, 
  fe_nm_type  = "fe_rank", 
  check_input = TRUE)

# Number of functional entities - 61 possible functional entities (2.82%)
sp_to_fe_sharks$"fe_nm"

# Calculate and display metrics
alpha_fd_fe_sharks <- alpha.fd.fe(asb_sp_sharks_occ,
                                  sp_to_fe_sharks,
                                  ind_nm = c("fred", "fored", "fvuln"),
                                  check_input = TRUE,
                                  details_returned = TRUE)

# Form dataframe of results for plotting
FEmetrics_taxon <- as.data.frame(alpha_fd_fe_sharks$asb_fdfe) %>% 
  tibble::rownames_to_column("Epoch")
FEmetrics_taxon$Epoch <- ordered(FEmetrics_taxon$Epoch, levels=c("Palaeocene","Eocene","Oligocene","Miocene",
                                                                 "Pliocene","Pleistocene"))

# Construct trait distance matrix using species (using FEs not needed since all species unique here)
sp_dist_sharks <- mFD::funct.dist(
  sp_tr         = sharks_traits,
  tr_cat        = sharks_traits_cat,
  metric        = "gower",
  scale_euclid  = "scale_center",
  ordinal_var   = "classic",
  weight_type   = "equal",
  stop_if_NA    = TRUE)

# Assess quality
fspaces_quality_sharks <- mFD::quality.fspaces(
  sp_dist             = sp_dist_sharks,
  maxdim_pcoa         = 10,
  deviation_weighting = c("absolute","squared"),
  fdist_scaling       = FALSE,
  fdendro             = "average")

# Return coordinates of each axis
sp_faxes_coord_sharks <- fspaces_quality_sharks$"details_fspaces"$"sp_pc_coord"

# Compute the range of functional axes:
range_sp_coord  <- range(sp_faxes_coord_sharks)

range_faxes <- range_sp_coord +
  c(-1, 1) * (range_sp_coord[2] - range_sp_coord[1]) * 0.05
range_faxes

# get species coordinates along the axes:
sp_faxes_coord_xy <- sp_faxes_coord_sharks[, c("PC1", "PC2")]
sp_faxes_coord_xz <- sp_faxes_coord_sharks[, c("PC1", "PC3")]
sp_faxes_coord_yz <- sp_faxes_coord_sharks[, c("PC2", "PC3")]

# Calculate FD metrics - needed for functional space plotting
alpha_fd_indices_sharks <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = sp_faxes_coord_sharks[ , c("PC1", "PC2", "PC3")],
  asb_sp_w         = baskets_sharks_weights,
  ind_vect         = c("fric", "fori", "fspe"),
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = TRUE)

## Filter by epochs (baskets)
## Palaeocene:
sp_filter_pal_12 <- mFD::sp.filter(asb_nm = c("Palaeocene"),
                                sp_faxes_coord = sp_faxes_coord_xy,
                                asb_sp_w = baskets_sharks_weights)
sp_filter_pal_13 <- mFD::sp.filter(asb_nm = c("Palaeocene"),
                                   sp_faxes_coord = sp_faxes_coord_xz,
                                   asb_sp_w = baskets_sharks_weights)
sp_filter_pal_23 <- mFD::sp.filter(asb_nm = c("Palaeocene"),
                                   sp_faxes_coord = sp_faxes_coord_yz,
                                   asb_sp_w = baskets_sharks_weights)
## get species coordinates (Palaeocene):
sp_faxes_coord_pal_12 <- sp_filter_pal_12$`species coordinates`
sp_faxes_coord_pal_13 <- sp_filter_pal_13$`species coordinates`
sp_faxes_coord_pal_23 <- sp_filter_pal_23$`species coordinates`

# Eocene:
sp_filter_eo_12 <- mFD::sp.filter(asb_nm = c("Eocene"),
                               sp_faxes_coord = sp_faxes_coord_xy,
                               asb_sp_w = baskets_sharks_weights)
sp_filter_eo_13 <- mFD::sp.filter(asb_nm = c("Eocene"),
                               sp_faxes_coord = sp_faxes_coord_xz,
                               asb_sp_w = baskets_sharks_weights)
sp_filter_eo_23 <- mFD::sp.filter(asb_nm = c("Eocene"),
                               sp_faxes_coord = sp_faxes_coord_yz,
                               asb_sp_w = baskets_sharks_weights)
## get species coordinates (Eocene):
sp_faxes_coord_eo_12 <- sp_filter_eo_12$`species coordinates`
sp_faxes_coord_eo_13 <- sp_filter_eo_13$`species coordinates`
sp_faxes_coord_eo_23 <- sp_filter_eo_23$`species coordinates`

# Oligocene:
sp_filter_oli_12 <- mFD::sp.filter(asb_nm = c("Oligocene"),
                                sp_faxes_coord = sp_faxes_coord_xy,
                                asb_sp_w = baskets_sharks_weights)
sp_filter_oli_13 <- mFD::sp.filter(asb_nm = c("Oligocene"),
                                sp_faxes_coord = sp_faxes_coord_xz,
                                asb_sp_w = baskets_sharks_weights)
sp_filter_oli_23 <- mFD::sp.filter(asb_nm = c("Oligocene"),
                                sp_faxes_coord = sp_faxes_coord_yz,
                                asb_sp_w = baskets_sharks_weights)
## get species coordinates (Oligocene):
sp_faxes_coord_oli_12 <- sp_filter_oli_12$`species coordinates`
sp_faxes_coord_oli_13 <- sp_filter_oli_13$`species coordinates`
sp_faxes_coord_oli_23 <- sp_filter_oli_23$`species coordinates`

# Miocene:
sp_filter_mio_12 <- mFD::sp.filter(asb_nm = c("Miocene"),
                                sp_faxes_coord = sp_faxes_coord_xy,
                                asb_sp_w = baskets_sharks_weights)
sp_filter_mio_13 <- mFD::sp.filter(asb_nm = c("Miocene"),
                                sp_faxes_coord = sp_faxes_coord_xz,
                                asb_sp_w = baskets_sharks_weights)
sp_filter_mio_23 <- mFD::sp.filter(asb_nm = c("Miocene"),
                                sp_faxes_coord = sp_faxes_coord_yz,
                                asb_sp_w = baskets_sharks_weights)
## get species coordinates (Miocene):
sp_faxes_coord_mio_12 <- sp_filter_mio_12$`species coordinates`
sp_faxes_coord_mio_13 <- sp_filter_mio_13$`species coordinates`
sp_faxes_coord_mio_23 <- sp_filter_mio_23$`species coordinates`

# Pliocene:
sp_filter_plio_12 <- mFD::sp.filter(asb_nm = c("Pliocene"),
                                 sp_faxes_coord = sp_faxes_coord_xy,
                                 asb_sp_w = baskets_sharks_weights)
sp_filter_plio_13 <- mFD::sp.filter(asb_nm = c("Pliocene"),
                                 sp_faxes_coord = sp_faxes_coord_xz,
                                 asb_sp_w = baskets_sharks_weights)
sp_filter_plio_23 <- mFD::sp.filter(asb_nm = c("Pliocene"),
                                 sp_faxes_coord = sp_faxes_coord_yz,
                                 asb_sp_w = baskets_sharks_weights)
## get species coordinates (Pliocene):
sp_faxes_coord_plio_12 <- sp_filter_plio_12$`species coordinates`
sp_faxes_coord_plio_13 <- sp_filter_plio_13$`species coordinates`
sp_faxes_coord_plio_23 <- sp_filter_plio_23$`species coordinates`

# Pleistocene:
sp_filter_ple_12 <- mFD::sp.filter(asb_nm = c("Pleistocene"),
                                sp_faxes_coord = sp_faxes_coord_xy,
                                asb_sp_w = baskets_sharks_weights)
sp_filter_ple_13 <- mFD::sp.filter(asb_nm = c("Pleistocene"),
                                sp_faxes_coord = sp_faxes_coord_xz,
                                asb_sp_w = baskets_sharks_weights)
sp_filter_ple_23 <- mFD::sp.filter(asb_nm = c("Pleistocene"),
                                sp_faxes_coord = sp_faxes_coord_yz,
                                asb_sp_w = baskets_sharks_weights)
## get species coordinates (Pleistocene):
sp_faxes_coord_ple_12 <- sp_filter_ple_12$`species coordinates`
sp_faxes_coord_ple_13 <- sp_filter_ple_13$`species coordinates`
sp_faxes_coord_ple_23 <- sp_filter_ple_23$`species coordinates`

# Recent:
sp_filter_rec_12 <- mFD::sp.filter(asb_nm = c("Recent"),
                                sp_faxes_coord = sp_faxes_coord_xy,
                                asb_sp_w = baskets_sharks_weights)
sp_filter_rec_13 <- mFD::sp.filter(asb_nm = c("Recent"),
                                   sp_faxes_coord = sp_faxes_coord_xz,
                                   asb_sp_w = baskets_sharks_weights)
sp_filter_rec_23 <- mFD::sp.filter(asb_nm = c("Recent"),
                                   sp_faxes_coord = sp_faxes_coord_yz,
                                   asb_sp_w = baskets_sharks_weights)
## get species coordinates (Recent):
sp_faxes_coord_rec_12 <- sp_filter_rec_12$`species coordinates`
sp_faxes_coord_rec_13 <- sp_filter_rec_13$`species coordinates`
sp_faxes_coord_rec_23 <- sp_filter_rec_23$`species coordinates`

### Add epochs to fuse-formed dataframe and filter to identify highest contributors of FOri & FSpe
GE <- filtered_baskets.range %>% 
  select(Taxon_corrected, Current_status) %>% 
  distinct(Taxon_corrected,Current_status) %>%  
  mutate(Current_status = case_when(Current_status == "Extant" ~ "1", TRUE ~ Current_status)) %>% 
  mutate(Current_status = case_when(Current_status == "Extinct" ~ "0", TRUE ~ Current_status))

GE <- GE[order(GE$Taxon_corrected),]
GE$Current_status <- as.numeric(GE$Current_status)

load("C:/Users/2022207/Dropbox/Jack's PhD/Chapter 2. FD changes over time/Analyses/Current Analyses/R codes/Taxon code/Data/FUn_Fsp metrics.RData")

# Form dataset with mean metrics, filtered taxa, status, taxonomic & epoch information
FUSE_sum <- FUSE_spp %>% 
  select(Taxon,FUSE,FUGE,FSGE,FUn_std,FSp_std) %>% 
  group_by(Taxon) %>% 
  summarise(FUSE = mean(FUSE),
            FUGE = mean(FUGE),
            FSGE = mean(FSGE),
            Fori = mean(FUn_std),
            Fspe = mean(FSp_std))
Taxa_info <- data %>% 
  group_by(Taxon_corrected) %>% 
  distinct(Taxon_corrected,Order,Family,Genus,Species) %>% 
  as.data.frame()
Taxa_filtered <- Taxa_info %>%
  anti_join(selected_taxa, by = "Taxon_corrected")
Taxa_filtered <- Taxa_filtered[order(Taxa_filtered$Taxon_corrected),]

FUSE_sum$Status <- GE$Current_status
FUSE_sum$Status <- as.factor(FUSE_sum$Status)
FUSE_sum <- FUSE_sum %>% 
  mutate(Status = case_when(Status == "1" ~ "Extant", TRUE ~ Status)) %>% 
  mutate(Status = case_when(Status == "0" ~ "Extinct", TRUE ~ Status))
FUSE_sum$Order <- Taxa_filtered$Order
FUSE_sum$Palaeocene <- filtered_baskets.range$Palaeocene
FUSE_sum$Eocene <- filtered_baskets.range$Eocene
FUSE_sum$Oligocene <- filtered_baskets.range$Oligocene
FUSE_sum$Miocene <- filtered_baskets.range$Miocene
FUSE_sum$Pliocene <- filtered_baskets.range$Pliocene
FUSE_sum$Pleistocene <- filtered_baskets.range$Pleistocene
FUSE_sum$Recent <- filtered_baskets.range$Recent

# Identify highest contributors through time by filtering by epoch
# Palaeocene
FUSE_pal <- FUSE_sum %>%
  select(Taxon,Fori,Fspe,Status,Order,Palaeocene) %>% 
  filter(Palaeocene == "1")

pal_fori <- FUSE_pal %>%
  arrange(desc(Fori)) 
head(pal_fori, 1)

pal_fspe <- FUSE_pal %>%
  arrange(desc(Fspe)) 
head(pal_fspe, 1)
## Highest contributors: Squaliodalatias sp. (FOri) & Palaeocarcharodon orientalis (FSpe)

# Eocene
FUSE_eo <- FUSE_sum %>%
  select(Taxon,Fori,Fspe,Status,Order,Eocene) %>% 
  filter(Eocene == "1")

eo_fori <- FUSE_eo %>%
  arrange(desc(Fori)) 
head(eo_fori, 1)

eo_fspe <- FUSE_eo %>%
  arrange(desc(Fspe)) 
head(eo_fspe, 1)
## Highest contributors: Heterodontus woodwardi (FOri) & Otodus nodai (FSpe)

# Oligocene
FUSE_oli <- FUSE_sum %>%
  select(Taxon,Fori,Fspe,Status,Order,Oligocene) %>% 
  filter(Oligocene == "1")

oli_fori <- FUSE_oli %>%
  arrange(desc(Fori)) 
head(oli_fori, 1)

oli_fspe <- FUSE_oli %>%
  arrange(desc(Fspe)) 
head(oli_fspe, 1)
## Highest contributors: Dalatias sp. (FOri) & Otodus angustidens (FSpe)

# Miocene
FUSE_mio <- FUSE_sum %>%
  select(Taxon,Fori,Fspe,Status,Order,Miocene) %>% 
  filter(Miocene == "1")

mio_fori <- FUSE_mio %>%
  arrange(desc(Fori)) 
head(mio_fori, 1)

mio_fspe <- FUSE_mio %>%
  arrange(desc(Fspe)) 
head(mio_fspe, 1)
## Highest contributors: Echinorhinus blakei (FOri) & Otodus angustidens (FSpe)

# Pliocene
FUSE_plio <- FUSE_sum %>%
  select(Taxon,Fori,Fspe,Status,Order,Pliocene) %>% 
  filter(Pliocene == "1")

plio_fori <- FUSE_plio %>%
  arrange(desc(Fori)) 
head(plio_fori, 1)

plio_fspe <- FUSE_plio %>%
  arrange(desc(Fspe)) 
head(plio_fspe, 1)
## Highest contributors: Echinorhinus blakei (FOri) & Otodus megalodon (FSpe)

# Pleistocene
FUSE_ple <- FUSE_sum %>%
  select(Taxon,Fori,Fspe,Status,Order,Pleistocene) %>% 
  filter(Pleistocene == "1")

ple_fori <- FUSE_ple %>%
  arrange(desc(Fori)) 
head(ple_fori, 1)

ple_fspe <- FUSE_ple %>%
  arrange(desc(Fspe)) 
head(ple_fspe, 1)
## Highest contributors: Megachasma sp. (FOri) & Dalatias sp. (FSpe)

# Recent
FUSE_rec <- FUSE_sum %>%
  select(Taxon,Fori,Fspe,Status,Order,Recent) %>% 
  filter(Recent == "1")

rec_fori <- FUSE_rec %>%
  arrange(desc(Fori)) 
head(rec_fori, 1)

rec_fspe <- FUSE_rec %>%
  arrange(desc(Fspe)) 
head(rec_fspe, 1)
## Highest contributors: Sphyrna tiburo (FOri) & Galeocerdo cuvier (FSpe)
