#########################################################################################################################
# 09. Recent resampling analyses
## This R code provides how we account for Recent sampling
## it produces 4 Rdata files used to plot the final results
#########################################################################################################################

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
library(foreach)
library(parallel)
library(doParallel)
library(deeptime)
library(grid)
library(gridExtra)

# Load data
load(file="C:/Users/2022207/Dropbox/Jack's PhD/Chapter 2. FD changes over time/Analyses/Current Analyses/R codes/Cleaned data.RData")

# Functional spaces per epoch - do not use "selected taxa" so as to maintain *all* Recent taxa
## Occurrence matrix - need to group by Taxa by epoch so we can mark their occurrences through time 
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
  mutate(Recent = case_when(Recent == "1" && Taxon_corrected == "Echinorhinus sp."~"0", TRUE ~ Recent)) %>%
  mutate(Recent = case_when(Recent == "1" && Taxon_corrected == "Galeocerdo sp."~"0", TRUE ~ Recent)) %>% 
  mutate(Recent = case_when(Recent == "1" && Taxon_corrected == "Galeorhinus sp."~"0", TRUE ~ Recent)) %>% 
  mutate(Recent = case_when(Recent == "1" && Taxon_corrected == "Galeus sp."~"0", TRUE ~ Recent)) %>% 
  mutate(Recent = case_when(Recent == "1" && Taxon_corrected == "Hemipristis sp."~"0", TRUE ~ Recent)) %>% 
  mutate(Recent = case_when(Recent == "1" && Taxon_corrected == "Heptranchias sp."~"0", TRUE ~ Recent)) %>% 
  mutate(Recent = case_when(Recent == "1" && Taxon_corrected == "Isurus sp."~"0", TRUE ~ Recent)) %>%
  mutate(Recent = case_when(Recent == "1" && Taxon_corrected == "Lamna sp."~"0", TRUE ~ Recent)) %>%
  mutate(Recent = case_when(Recent == "1" && Taxon_corrected == "Megachasma sp."~"0", TRUE ~ Recent)) %>%
  mutate(Recent = case_when(Recent == "1" && Taxon_corrected == "Mitsukurina sp."~"0", TRUE ~ Recent)) %>%
  mutate(Recent = case_when(Recent == "1" && Taxon_corrected == "Nebrius sp."~"0", TRUE ~ Recent)) %>% 
  mutate(Recent = case_when(Recent == "1" && Taxon_corrected == "Negaprion sp."~"0", TRUE ~ Recent)) %>%
  mutate(Recent = case_when(Recent == "1" && Taxon_corrected == "Notorynchus sp."~"0", TRUE ~ Recent)) %>% 
  mutate(Recent = case_when(Recent == "1" && Taxon_corrected == "Odontaspis sp."~"0", TRUE ~ Recent))

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

# Make taxa row names & tidy dataframe to form final matrix
baskets <- baskets.range.taxon %>% 
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

# Calculate average CH, CW & mode CE, LC, XO, LO per species & genus
## Used to form functional space due to being most common tooth morphology per species; see Code S4 for iteration-based analyses
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

# Make species-trait matrix
sharks_traits <- Av_sharks %>%
  remove_rownames %>% 
  column_to_rownames(var="Taxon_corrected")

# Form trait category matrix
sharks_traits_cat <- tibble(trait_name = c("CH", "CW", "CE", "LC", "XO", "LO"),
                            trait_type = c("O", "O", "N", "N", "N", "N"),
                            trait_weight = c(0.5,0.5,1,1,0.33,0.67))

# Summarise dataset
sharks_traits_summ <- mFD::sp.tr.summary(
  tr_cat     = sharks_traits_cat,   
  sp_tr      = sharks_traits, 
  stop_if_NA = TRUE)

# Reveal number of states for each dental character
sharks_traits_summ$tr_summary_list

# Use states to check number of possible combinations - 2160 possible combinations
4*3*3*2*6*5

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

# Number of functional entities - 65 possible functional entities per iteration (3%); 1 extra entity added
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
                                                                 "Pliocene","Pleistocene","Recent"))

# Construct trait distance matrix using species
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

# Check pocket - data best represented by 4 dimensions/axes; but difference is negilgible from 3 axes
round(fspaces_quality_sharks$"quality_fspaces", 5)

# Plot quality - produces Figure S2
mFD::quality.fspaces.plot(
  fspaces_quality            = fspaces_quality_sharks,
  quality_metric             = "mad",
  fspaces_plot               = c("tree_average", "pcoa_2d", "pcoa_3d", 
                                 "pcoa_4d", "pcoa_5d"),
  name_file                  = NULL,
  range_dist                 = NULL,
  range_dev                  = NULL,
  range_qdev                 = NULL,
  gradient_deviation         = c(neg = "darkblue", nul = "grey80", pos = "darkred"),
  gradient_deviation_quality = c(low = "yellow", high = "red"),
  x_lab                      = "Trait-based distance")

# Check % of inertia accounted for by axes
eig <- fspaces_quality_sharks$details_fspaces$pc_eigenvalues %>% 
  as.data.frame()

eig$variance <- (eig$Eigenvalues/68.09946)*100  
# 3 axes accounts for 83.5% of total inertia; even more than when removing non-fossil spp.

# Return coordinates of each axis
sp_faxes_coord_sharks <- fspaces_quality_sharks$"details_fspaces"$"sp_pc_coord"

# Correlate each axis against traits
sharks_tr_faxes <- mFD::traits.faxes.cor(
  sp_tr          = sharks_traits, 
  sp_faxes_coord = sp_faxes_coord_sharks[ , c("PC1", "PC2", "PC3")], 
  plot           = TRUE)

# Print traits with significant effect - produces Table S3
sharks_tr_faxes$"tr_faxes_stat"[which(sharks_tr_faxes$"tr_faxes_stat"$"p.value" < 0.05), ]

# Return plots - plots in blue indicate which traits drive each axis - produces Figure S3
sharks_tr_faxes$"tr_faxes_plot"
# Most correlated DCs per axis the same as before

## Form functional space of first 3 axes
big_plot <- mFD::funct.space.plot(
  sp_faxes_coord  = sp_faxes_coord_sharks[ , c("PC1", "PC2", "PC3")],
  faxes           = c("PC1", "PC2", "PC3"),
  name_file       = NULL,
  faxes_nm        = NULL,
  range_faxes     = c(NA, NA),
  color_bg        = NA,
  color_pool      = "darkgreen",
  fill_pool       = "white",
  shape_pool      = 21,
  size_pool       = 1,
  plot_ch         = TRUE,
  color_ch        = "#000000",
  fill_ch         = "white",
  alpha_ch        = 1,
  plot_vertices   = TRUE,
  color_vert      = "blueviolet",
  fill_vert       = "blueviolet",
  shape_vert      = 23,
  size_vert       = 1,
  plot_sp_nm      = NULL,
  nm_size         = 3,
  nm_color        = "black",
  nm_fontface     = "plain",
  check_input     = TRUE)
big_plot$patchwork

# Calculate FD metrics
alpha_fd_indices_sharks <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = sp_faxes_coord_sharks[ , c("PC1", "PC2", "PC3")],
  asb_sp_w         = baskets_sharks_weights,
  ind_vect         = c("fric", "fori", "fspe"),
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = TRUE)

fd_ind_values_sharks <- alpha_fd_indices_sharks$"functional_diversity_indices"

## Form dataframe - combine with FE metrics
FDindices_taxon <- as.data.frame(fd_ind_values_sharks) %>% 
  tibble::rownames_to_column("Epoch")

# Form spaces for each epoch, with highest contributors highlighted - Figure 1b-h
## Set axis limits
x_limits_PC1 <- c(-0.4, 0.7)
y_limits_PC12 <- c(-0.3, 1) 
y_limits_PC13 <- c(-0.4, 0.3)

## Palaeocene
plots_Pal <- mFD::alpha.multidim.plot(
  output_alpha_fd_multidim = alpha_fd_indices_sharks,
  plot_asb_nm              = "Palaeocene",
  ind_nm                   = c("fric", "fori", "fspe"),
  faxes                    = NULL,
  faxes_nm                 = NULL,
  range_faxes              = c(NA, NA),
  color_bg                 = NA,
  shape_sp                 = c(pool = 21, asb1 = 21),
  size_sp                  = c(pool = 1, asb1 = 1),
  color_sp                 = c(pool = "darkgrey", asb1 = "#FBA75F"),
  color_vert               = c(pool = "darkgrey", asb1 = "#FBA75F"),
  fill_sp                  = c(pool = "darkgrey", asb1 = "#FBA75F"),
  fill_vert                = c(pool = "darkgrey", asb1 = "#FBA75F"),
  color_ch                 = c(pool = 'black', asb1 = "#FBA75F"),
  fill_ch                  = c(pool = "white", asb1 = "#FBA75F"),
  alpha_ch                 = c(pool = 1, asb1 = 0.3),
  shape_centroid_fdis      = c(asb1 = 22),
  shape_centroid_fdiv      = c(asb1 = 22),
  shape_centroid_fspe      = 23,
  color_centroid_fspe      = "black",
  size_sp_nm               = 3, 
  color_sp_nm              = "black",
  plot_sp_nm               = NULL,
  fontface_sp_nm           = "plain",
  save_file                = FALSE,
  check_input              = TRUE)

Pal <- plots_Pal$fric$PC1_PC2 +
  labs(x = "PCoA1", y = "PCoA2")+
  coord_fixed(ratio = 1, xlim = x_limits_PC1, ylim = y_limits_PC12)

## Eocene
plots_Eo <- mFD::alpha.multidim.plot(
  output_alpha_fd_multidim = alpha_fd_indices_sharks,
  plot_asb_nm              = "Eocene",
  ind_nm                   = c("fric", "fori", "fspe"),
  faxes                    = NULL,
  faxes_nm                 = NULL,
  range_faxes              = c(NA, NA),
  color_bg                 = NA,
  shape_sp                 = c(pool = 21, asb1 = 21),
  size_sp                  = c(pool = 1, asb1 = 1),
  color_sp                 = c(pool = "darkgrey", asb1 = "#FDB46C"),
  color_vert               = c(pool = "darkgrey", asb1 = "#FDB46C"),
  fill_sp                  = c(pool = "darkgrey", asb1 = "#FDB46C"),
  fill_vert                = c(pool = "darkgrey", asb1 = "#FDB46C"),
  color_ch                 = c(pool = 'black', asb1 = "#FDB46C"),
  fill_ch                  = c(pool = "white", asb1 = "#FDB46C"),
  alpha_ch                 = c(pool = 1, asb1 = 0.3),
  shape_centroid_fdis      = c(asb1 = 22),
  shape_centroid_fdiv      = c(asb1 = 22),
  shape_centroid_fspe      = 23,
  color_centroid_fspe      = "black",
  size_sp_nm               = 3, 
  color_sp_nm              = "black",
  plot_sp_nm               = NULL,
  fontface_sp_nm           = "plain",
  save_file                = FALSE,
  check_input              = TRUE)

Eo <- plots_Eo$fric$PC1_PC2 +
  labs(x = "PCoA1", y = "PCoA2")+
  coord_fixed(ratio = 1, xlim = x_limits_PC1, ylim = y_limits_PC12)

## Oligocene
plots_Oli <- mFD::alpha.multidim.plot(
  output_alpha_fd_multidim = alpha_fd_indices_sharks,
  plot_asb_nm              = "Oligocene",
  ind_nm                   = c("fric", "fori", "fspe"),
  faxes                    = NULL,
  faxes_nm                 = NULL,
  range_faxes              = c(NA, NA),
  color_bg                 = NA,
  shape_sp                 = c(pool = 21, asb1 = 21),
  size_sp                  = c(pool = 1, asb1 = 1),
  color_sp                 = c(pool = "darkgrey", asb1 = "#FDC07A"),
  color_vert               = c(pool = "darkgrey", asb1 = "#FDC07A"),
  fill_sp                  = c(pool = "darkgrey", asb1 = "#FDC07A"),
  fill_vert                = c(pool = "darkgrey", asb1 = "#FDC07A"),
  color_ch                 = c(pool = 'black', asb1 = "#FDC07A"),
  fill_ch                  = c(pool = "white", asb1 = "#FDC07A"),
  alpha_ch                 = c(pool = 1, asb1 = 0.3),
  shape_centroid_fdis      = c(asb1 = 22),
  shape_centroid_fdiv      = c(asb1 = 22),
  shape_centroid_fspe      = 23,
  color_centroid_fspe      = "black",
  size_sp_nm               = 3, 
  color_sp_nm              = "black",
  plot_sp_nm               = NULL,
  fontface_sp_nm           = "plain",
  save_file                = FALSE,
  check_input              = TRUE)

Oli <- plots_Oli$fric$PC1_PC2 +
  labs(x = "PCoA1", y = "PCoA2")+
  coord_fixed(ratio = 1, xlim = x_limits_PC1, ylim = y_limits_PC12)

## Miocene
plots_Mio <- mFD::alpha.multidim.plot(
  output_alpha_fd_multidim = alpha_fd_indices_sharks,
  plot_asb_nm              = "Miocene",
  ind_nm                   = c("fric", "fori", "fspe"),
  faxes                    = NULL,
  faxes_nm                 = NULL,
  range_faxes              = c(NA, NA),
  color_bg                 = NA,
  shape_sp                 = c(pool = 21, asb1 = 21),
  size_sp                  = c(pool = 1, asb1 = 1),
  color_sp                 = c(pool = "darkgrey", asb1 = "#FFFF90"),
  color_vert               = c(pool = "darkgrey", asb1 = "#FFFF90"),
  fill_sp                  = c(pool = "darkgrey", asb1 = "#FFFF90"),
  fill_vert                = c(pool = "darkgrey", asb1 = "#FFFF90"),
  color_ch                 = c(pool = 'black', asb1 = "#FFFF90"),
  fill_ch                  = c(pool = "white", asb1 = "#FFFF90"),
  alpha_ch                 = c(pool = 1, asb1 = 0.3),
  shape_centroid_fdis      = c(asb1 = 22),
  shape_centroid_fdiv      = c(asb1 = 22),
  shape_centroid_fspe      = 23,
  color_centroid_fspe      = "black",
  size_sp_nm               = 3, 
  color_sp_nm              = "black",
  plot_sp_nm               = NULL,
  fontface_sp_nm           = "plain",
  save_file                = FALSE,
  check_input              = TRUE)

Mio <- plots_Mio$fric$PC1_PC2 +
  labs(x = "PCoA1", y = "PCoA2")+
  coord_fixed(ratio = 1, xlim = x_limits_PC1, ylim = y_limits_PC12)

## Pliocene
plots_Plio <- mFD::alpha.multidim.plot(
  output_alpha_fd_multidim = alpha_fd_indices_sharks,
  plot_asb_nm              = "Pliocene",
  ind_nm                   = c("fric", "fori", "fspe"),
  faxes                    = NULL,
  faxes_nm                 = NULL,
  range_faxes              = c(NA, NA),
  color_bg                 = NA,
  shape_sp                 = c(pool = 21, asb1 = 21),
  size_sp                  = c(pool = 1, asb1 = 1),
  color_sp                 = c(pool = "darkgrey", asb1 = "#FFFF99"),
  color_vert               = c(pool = "darkgrey", asb1 = "#FFFF99"),
  fill_sp                  = c(pool = "darkgrey", asb1 = "#FFFF99"),
  fill_vert                = c(pool = "darkgrey", asb1 = "#FFFF99"),
  color_ch                 = c(pool = 'black', asb1 = "#FFFF99"),
  fill_ch                  = c(pool = "white", asb1 = "#FFFF99"),
  alpha_ch                 = c(pool = 1, asb1 = 0.3),
  shape_centroid_fdis      = c(asb1 = 22),
  shape_centroid_fdiv      = c(asb1 = 22),
  shape_centroid_fspe      = 23,
  color_centroid_fspe      = "black",
  size_sp_nm               = 3, 
  color_sp_nm              = "black",
  plot_sp_nm               = NULL,
  fontface_sp_nm           = "plain",
  save_file                = FALSE,
  check_input              = TRUE)

Plio <- plots_Plio$fric$PC1_PC2 +
  labs(x = "PCoA1", y = "PCoA2")+
  coord_fixed(ratio = 1, xlim = x_limits_PC1, ylim = y_limits_PC12)

## Pleistocene
plots_Ple <- mFD::alpha.multidim.plot(
  output_alpha_fd_multidim = alpha_fd_indices_sharks,
  plot_asb_nm              = "Pleistocene",
  ind_nm                   = c("fric", "fori", "fspe"),
  faxes                    = NULL,
  faxes_nm                 = NULL,
  range_faxes              = c(NA, NA),
  color_bg                 = NA,
  shape_sp                 = c(pool = 21, asb1 = 21),
  size_sp                  = c(pool = 1, asb1 = 1),
  color_sp                 = c(pool = "darkgrey", asb1 = "#FFF2AE"),
  color_vert               = c(pool = "darkgrey", asb1 = "#FFF2AE"),
  fill_sp                  = c(pool = "darkgrey", asb1 = "#FFF2AE"),
  fill_vert                = c(pool = "darkgrey", asb1 = "#FFF2AE"),
  color_ch                 = c(pool = 'black', asb1 = "#FFF2AE"),
  fill_ch                  = c(pool = "white", asb1 = "#FFF2AE"),
  alpha_ch                 = c(pool = 1, asb1 = 0.3),
  shape_centroid_fdis      = c(asb1 = 22),
  shape_centroid_fdiv      = c(asb1 = 22),
  shape_centroid_fspe      = 23,
  color_centroid_fspe      = "black",
  size_sp_nm               = 3, 
  color_sp_nm              = "black",
  plot_sp_nm               = NULL,
  fontface_sp_nm           = "plain",
  save_file                = FALSE,
  check_input              = TRUE)

Ple <- plots_Ple$fric$PC1_PC2 +
  labs(x = "PCoA1", y = "PCoA2")+
  coord_fixed(ratio = 1, xlim = x_limits_PC1, ylim = y_limits_PC12)

## Recent
plots_Rec <- mFD::alpha.multidim.plot(
  output_alpha_fd_multidim = alpha_fd_indices_sharks,
  plot_asb_nm              = "Recent",
  ind_nm                   = c("fric", "fori", "fspe"),
  faxes                    = NULL,
  faxes_nm                 = NULL,
  range_faxes              = c(NA, NA),
  color_bg                 = NA,
  shape_sp                 = c(pool = 21, asb1 = 21),
  size_sp                  = c(pool = 1, asb1 = 1),
  color_sp                 = c(pool = "darkgrey", asb1 = "#FEF2E0"),
  color_vert               = c(pool = "darkgrey", asb1 = "#FEF2E0"),
  fill_sp                  = c(pool = "darkgrey", asb1 = "#FEF2E0"),
  fill_vert                = c(pool = "darkgrey", asb1 = "#FEF2E0"),
  color_ch                 = c(pool = 'black', asb1 = "#FEF2E0"),
  fill_ch                  = c(pool = "white", asb1 = "#FEF2E0"),
  alpha_ch                 = c(pool = 1, asb1 = 0.3),
  shape_centroid_fdis      = c(asb1 = 22),
  shape_centroid_fdiv      = c(asb1 = 22),
  shape_centroid_fspe      = 23,
  color_centroid_fspe      = "black",
  size_sp_nm               = 3, 
  color_sp_nm              = "black",
  plot_sp_nm               = NULL,
  fontface_sp_nm           = "plain",
  save_file                = FALSE,
  check_input              = TRUE)

Rec <- plots_Rec$fric$PC1_PC2 +
  labs(x = "PCoA1", y = "PCoA2")+
  coord_fixed(ratio = 1, xlim = x_limits_PC1, ylim = y_limits_PC12)

Spaces <- plot_grid(Pal,Eo,Oli,Mio,Plio,Ple,Rec,
                        labels = c("Palaeocene","Eocene","Oligocene",
                                   "Miocene","Pliocene","Pleistocene","Recent"), 
                        label_size = 12,align = "hv", label_fontface = "bold", hjust = -0.15,  nrow=2)

## Now loop based on FTUs to see if results are similar (prior to resampling) - on 1000 iterations for this test
## Form dataframe of traits and occurrences; to be split into two input matrices in the loop
# Occurrence matrix
# Form FTU unique combinations, with epoch included
abun <- data %>% 
  group_by(Taxon_corrected) %>% 
  distinct(CH,CW,CE,LC,XO,LO,Epoch,Epoch_earliest,Epoch_latest,Current_status)

# Order taxa alphabetically and order Epochs
abun <- abun[order(abun$Taxon_corrected),]
abun$Epoch <- ordered(abun$Epoch, levels=c("Palaeocene","Eocene","Oligocene","Miocene",
                                           "Pliocene","Pleistocene", "Recent","nonmatch"))

# Form FTUs column
abun$FTU<-paste(abun$Taxon_corrected, abun$CH, abun$CW, abun$CE, abun$LC, abun$XO, abun$LO, sep="+")

# Cast into wide-format data
wp.abun <- abun %>% 
  pivot_wider(names_from = Epoch, values_from = FTU)

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
  mutate(Recent = case_when(Recent == "1" && Taxon_corrected == "Echinorhinus sp."~"0", TRUE ~ Recent)) %>%
  mutate(Recent = case_when(Recent == "1" && Taxon_corrected == "Galeocerdo sp."~"0", TRUE ~ Recent)) %>% 
  mutate(Recent = case_when(Recent == "1" && Taxon_corrected == "Galeorhinus sp."~"0", TRUE ~ Recent)) %>% 
  mutate(Recent = case_when(Recent == "1" && Taxon_corrected == "Galeus sp."~"0", TRUE ~ Recent)) %>% 
  mutate(Recent = case_when(Recent == "1" && Taxon_corrected == "Hemipristis sp."~"0", TRUE ~ Recent)) %>% 
  mutate(Recent = case_when(Recent == "1" && Taxon_corrected == "Heptranchias sp."~"0", TRUE ~ Recent)) %>% 
  mutate(Recent = case_when(Recent == "1" && Taxon_corrected == "Isurus sp."~"0", TRUE ~ Recent)) %>%
  mutate(Recent = case_when(Recent == "1" && Taxon_corrected == "Lamna sp."~"0", TRUE ~ Recent)) %>%
  mutate(Recent = case_when(Recent == "1" && Taxon_corrected == "Megachasma sp."~"0", TRUE ~ Recent)) %>%
  mutate(Recent = case_when(Recent == "1" && Taxon_corrected == "Mitsukurina sp."~"0", TRUE ~ Recent)) %>%
  mutate(Recent = case_when(Recent == "1" && Taxon_corrected == "Nebrius sp."~"0", TRUE ~ Recent)) %>% 
  mutate(Recent = case_when(Recent == "1" && Taxon_corrected == "Negaprion sp."~"0", TRUE ~ Recent)) %>%
  mutate(Recent = case_when(Recent == "1" && Taxon_corrected == "Notorynchus sp."~"0", TRUE ~ Recent)) %>% 
  mutate(Recent = case_when(Recent == "1" && Taxon_corrected == "Odontaspis sp."~"0", TRUE ~ Recent))

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


# Make taxa row names & tidy dataframe to form final matrix
baskets <- baskets.range.taxon %>% 
  remove_rownames %>% 
  column_to_rownames(var="Taxon_corrected") %>% 
  select(-Current_status)

# Make species richness object (as this will be constant regardless of FTU variation) to add later
baskets.epoch <- baskets.range.taxon %>% 
  select(c(Taxon_corrected,Palaeocene,Eocene,Oligocene,Miocene,Pliocene,Pleistocene,Recent)) %>% 
  rowwise() %>% 
  mutate(Palaeocene = replace(Palaeocene,Palaeocene=="1",Taxon_corrected)) %>% 
  mutate(Eocene = replace(Eocene,Eocene=="1",Taxon_corrected)) %>% 
  mutate(Oligocene = replace(Oligocene,Oligocene=="1",Taxon_corrected)) %>% 
  mutate(Miocene = replace(Miocene,Miocene=="1",Taxon_corrected)) %>% 
  mutate(Pliocene = replace(Pliocene,Pliocene=="1",Taxon_corrected)) %>% 
  mutate(Pleistocene = replace(Pleistocene,Pleistocene=="1",Taxon_corrected)) %>% 
  mutate(Recent = replace(Recent,Recent=="1",Taxon_corrected))

baskets.epoch[baskets.epoch == 0]<-NA

Spp_richness <- c(length(unique(baskets.epoch$Palaeocene)),
                  length(unique(baskets.epoch$Eocene)),
                  length(unique(baskets.epoch$Oligocene)),
                  length(unique(baskets.epoch$Miocene)),
                  length(unique(baskets.epoch$Pliocene)),
                  length(unique(baskets.epoch$Pleistocene)),
                  length(unique(baskets.epoch$Recent)))

# Form functional taxonomic units (FTUs), grouping by species to account for intraspecific variability
FTU <- data %>% 
  group_by(Taxon_corrected) %>% 
  distinct(CH,CW,CE,LC,XO,LO)

# Order FTUs alphabetically to help ensure matches between species-trait matrix and occurrence matrix later
FTU <- FTU[order(FTU$Taxon_corrected),]

## Run iterations & calculate mean by taxon
# Set up model to run on n iterations
registerDoParallel(cores = 10)
getDoParWorkers()
res.TaxonVar<-NULL
res.TaxonVar<- lapply(1:1000, function(x){
  # Form sliced species-trait matrix
  sharks_traits.var <- FTU %>%
    group_by(Taxon_corrected) %>%
    slice_sample(n=1, replace = FALSE) %>% 
    remove_rownames %>% 
    column_to_rownames(var="Taxon_corrected")
  
  # Transpose to form final occurrence matrix
  baskets_sharks_weights.var <- t(baskets)
  
  # Make matrix and numeric
  baskets_sharks_weights.var <- data.matrix(baskets_sharks_weights.var, rownames.force = NA)
  class(baskets_sharks_weights.var) <- "numeric"
  
  # Summarise & extract occurrences
  asb_sp_sharks_summ.var <- mFD::asb.sp.summary(asb_sp_w = baskets_sharks_weights.var)
  
  asb_sp_sharks_occ.var <- asb_sp_sharks_summ.var$asb_sp_occ
  
  ############################################################################################
  # Make trait_cat matrix
  sharks_traits_cat.var <- tibble(trait_name = c("CH", "CW", "CE", "LC", "XO", "LO"),
                                  trait_type = c("O", "O", "N", "N", "N", "N"),
                                  trait_weight = c(0.5,0.5,1,1,0.33,0.67))
  
  ############################################################################################
  # Form the functional entities
  sp_to_fe_sharks.var <- mFD::sp.to.fe(
    sp_tr       = sharks_traits.var, 
    tr_cat      = sharks_traits_cat.var, 
    fe_nm_type  = "fe_rank", 
    check_input = TRUE)
  
  # Calculate and display metrics
  alpha_fd_fe_sharks.var <- alpha.fd.fe(asb_sp_sharks_occ.var,
                                        sp_to_fe_sharks.var,
                                        ind_nm = c("fred", "fored", "fvuln"),
                                        check_input = TRUE,
                                        details_returned = TRUE)
  
  # Form dataframe
  FEmetrics.var <- as.data.frame(alpha_fd_fe_sharks.var$asb_fdfe) %>% 
    tibble::rownames_to_column("Epoch")
  FEmetrics.var$Epoch <- ordered(FEmetrics.var$Epoch, levels=c("Palaeocene","Eocene","Oligocene","Miocene",
                                                               "Pliocene","Pleistocene","Recent"))
  
  # Construct trait distance matrix
  sp_dist_sharks.var <- mFD::funct.dist(
    sp_tr         = sharks_traits.var,
    tr_cat        = sharks_traits_cat.var,
    metric        = "gower",
    scale_euclid  = "scale_center",
    ordinal_var   = "classic",
    weight_type   = "equal",
    stop_if_NA    = TRUE)
  
  # Assess quality
  fspaces_quality_sharks.var <- mFD::quality.fspaces(
    sp_dist             = sp_dist_sharks.var,
    maxdim_pcoa         = 10,
    deviation_weighting = c("absolute","squared"),
    fdist_scaling       = FALSE,
    fdendro             = "average")
  
  # Return coordinates of each axis
  sp_faxes_coord_sharks.var <- fspaces_quality_sharks.var$"details_fspaces"$"sp_pc_coord"
  
  # Form alpha diversity indices
  alpha_fd_indices_sharks.var <- mFD::alpha.fd.multidim(
    sp_faxes_coord   = sp_faxes_coord_sharks.var[ , c("PC1", "PC2", "PC3")],
    asb_sp_w         = baskets_sharks_weights.var,
    ind_vect         = c("fric", "fori", "fspe"),
    scaling          = TRUE,
    check_input      = TRUE,
    details_returned = TRUE)
  
  fd_ind_values_sharks.var <- alpha_fd_indices_sharks.var$"functional_diversity_indices"
  
  ## Form dataframe
  FD.var <- as.data.frame(fd_ind_values_sharks.var) %>% 
    tibble::rownames_to_column("Epoch")
  FD.var$Epoch <- ordered(FD.var$Epoch, levels=c("Palaeocene","Eocene","Oligocene","Miocene",
                                                 "Pliocene","Pleistocene","Recent"))
  
  # Form list to merge datasets
  FDind.var = list(FEmetrics.var,FD.var)
  
  # Output
  FDind.var %>% 
    reduce(inner_join, by = "Epoch")
  
  #close loop
})

stopImplicitCluster()

# Merge lists into 1 dataframe

res_df_TaxonVar <- res.TaxonVar %>% 
  bind_rows()

# Remove sp_richn, which is the same as nb_sp - redundant
res_df_TaxonVar$sp_richn<-NULL

# Format dataframe to be loaded for comparison plots
FullRecent_TaxonVar<- res_df_TaxonVar %>% 
  select(Epoch:fspe) # remember to save as something
save(FullRecent_TaxonVar, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 2. FD changes over time/Analyses/Current Analyses/R codes/Taxon code/Data/Full_Recent.RData")

# Melt data - remember to save as something before clearing console
FullRecent_long<- melt(FullRecent_TaxonVar, id.vars= "Epoch")
save(FullRecent_long, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 2. FD changes over time/Analyses/Current Analyses/R codes/Taxon code/Data/Full_Recent_long.RData")


# Form dataframe of mean, median and standard deviation of all FD metrics
## Spp_richness + Median and SD values produce Table 1
Recent_var <- FullRecent_TaxonVar %>% 
  group_by(Epoch) %>%
  summarise(Sp_mean = mean(nb_sp),
            Sp_med = median(nb_sp),
            Sp_sd = sd(nb_sp),
            FE_mean = mean(nb_fe),
            FE_med = median(nb_fe),
            FE_sd = sd(nb_fe),
            Red_mean = mean(fred),
            Red_med = median(fred),
            Red_sd = sd(fred),
            Ored_mean = mean(fored),
            Ored_med = median(fored),
            Ored_sd = sd(fored),
            Vul_mean = mean(fvuln),
            Vul_med = median(fvuln),
            Vul_sd = sd(fvuln),
            FRic_mean = mean(fric),
            FRic_med = median(fric),
            FRic_sd = sd(fric),
            Fori_mean = mean(fori),
            Fori_med = median(fori),
            Fori_sd = sd(fori),
            Fspe_mean = mean(fspe),
            Fspe_med = median(fspe),
            Fspe_sd = sd(fspe))

Recent_var$Epoch <- ordered(Recent_var$Epoch, levels=c("Palaeocene","Eocene","Oligocene","Miocene",
                                                     "Pliocene","Pleistocene","Recent"))

# Plot
epochs_custom <- data.frame(
  name = c("Palaeocene","Eocene","Oligocene","Miocene","Pliocene","Pleistocene","Recent"), 
  max_age = c(1.5,2.5,3.5,4.5,5.5,6.5,7.6),
  min_age = c(0.4,1.5,2.5,3.5,4.5,5.5,6.5), 
  color = c("#FBA75F","#FDB46C","#FDC07A","#FFFF90","#FFFF99","#FFF2AE","#FEF2E0")
)
## Functional entities
FullRecent_FE <- FullRecent_long %>% 
  filter(variable == "nb_fe")
FE_full_Recent <- ggplot()+
  scale_x_discrete(limits=c("Palaeocene","Eocene","Oligocene","Miocene", "Pliocene", "Pleistocene", "Recent"))+
  geom_boxplot(data = FullRecent_FE, aes(x = Epoch, y = value),
               fill=c("#FBA75F","#FDB46C","#FDC07A","#FFFF90","#FFFF99","#FFF2AE","#FEF2E0"))+
  labs(x = "", y = "# FEs")+
  theme_minimal() +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size= 6.5, color= "black"),
        axis.title= element_text(size= 8), panel.background= element_rect(fill= "white"))+
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent")) +
  theme(axis.text.x = element_blank())

## Functional redundancy
FullRecent_Fred <- FullRecent_long %>% 
  filter(variable == "fred")
FRed_full_Recent <- ggplot()+
  scale_x_discrete(limits=c("Palaeocene","Eocene","Oligocene","Miocene", "Pliocene","Pleistocene","Recent"))+
  geom_boxplot(data = FullRecent_Fred, aes(x = Epoch, y = value), 
               fill=c("#FBA75F","#FDB46C","#FDC07A","#FFFF90","#FFFF99","#FFF2AE","#FEF2E0"))+
  labs(x = "", y = "FRed")+
  #ylim(0,1) +
  #scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  theme_minimal() +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size= 6.5, color= "black"),
        axis.title= element_text(size= 8), panel.background= element_rect(fill= "white"))+
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent")) +
  theme(axis.text.x = element_blank())

## Functional over-redundancy
FullRecent_Fored <- FullRecent_long %>% 
  filter(variable == "fored")
FOred_full_Recent <- ggplot()+
  scale_x_discrete(limits=c("Palaeocene","Eocene","Oligocene","Miocene", "Pliocene", "Pleistocene","Recent"))+
  geom_boxplot(data = FullRecent_Fored, aes(x = Epoch, y = value), 
               fill=c("#FBA75F","#FDB46C","#FDC07A","#FFFF90","#FFFF99","#FFF2AE","#FEF2E0"))+
  labs(x = "", y = "FOred")+
  #ylim(0.3,0.6) +
  #scale_y_continuous(limits = c(0.3, 0.6), breaks = seq(0.3, 0.6, by = 0.1)) +
  theme_minimal() +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size= 6.5, color= "black"),
        axis.title= element_text(size= 8), panel.background= element_rect(fill= "white"))+
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent"))+
  theme(axis.text.x = element_blank())

## Functional vulnerability
FullRecent_Fvuln <- FullRecent_long %>% 
  filter(variable == "fvuln")
FVuln_full_Recent <- ggplot()+
  scale_x_discrete(limits=c("Palaeocene","Eocene","Oligocene","Miocene", "Pliocene", "Pleistocene","Recent"))+
  geom_boxplot(data = FullRecent_Fvuln, aes(x = Epoch, y = value), 
               fill=c("#FBA75F","#FDB46C","#FDC07A","#FFFF90","#FFFF99","#FFF2AE","#FEF2E0"))+
  labs(x = "", y = "FVuln")+
  theme_minimal() +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size= 6.5, color= "black"),
        axis.title= element_text(size= 8), panel.background= element_rect(fill= "white"))+
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent"))+
  theme(axis.text.x = element_blank())

## Functional richness
FullRecent_FRic <- FullRecent_long %>% 
  filter(variable == "fric")
FRic_full_Recent <- ggplot()+
  scale_x_discrete(limits=c("Palaeocene","Eocene","Oligocene","Miocene", "Pliocene", "Pleistocene","Recent"))+
  geom_boxplot(data = FullRecent_FRic, aes(x = Epoch, y = value), 
               fill=c("#FBA75F","#FDB46C","#FDC07A","#FFFF90","#FFFF99","#FFF2AE","#FEF2E0"))+
  labs(x = "", y = "FRic")+
  theme_minimal() +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size= 6.5, color= "black"),
        axis.title= element_text(size= 8), panel.background= element_rect(fill= "white"))+
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent")) +
  theme(axis.text.x = element_blank())

## Functional originality
FullRecent_FOri <- FullRecent_long %>% 
  filter(variable == "fori")
FOri_full_Recent <- ggplot()+
  scale_x_discrete(limits=c("Palaeocene","Eocene","Oligocene","Miocene", "Pliocene", "Pleistocene","Recent"))+
  geom_boxplot(data = FullRecent_FOri, aes(x = Epoch, y = value), 
               fill=c("#FBA75F","#FDB46C","#FDC07A","#FFFF90","#FFFF99","#FFF2AE","#FEF2E0"))+
  labs(x = "", y = "FOri")+
  theme_minimal() +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size= 6.5, color= "black"),
        axis.title= element_text(size= 8), panel.background= element_rect(fill= "white"))+
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent")) +
  theme(axis.text.x = element_blank())


## Functional specialisation
FullRecent_FSpe <- FullRecent_long %>% 
  filter(variable == "fspe")
FSpe_full_Recent <- ggplot()+
  scale_x_discrete(limits=c("Palaeocene","Eocene","Oligocene","Miocene", "Pliocene", "Pleistocene","Recent"))+
  geom_boxplot(data = FullRecent_FSpe, aes(x = Epoch, y = value), 
               fill=c("#FBA75F","#FDB46C","#FDC07A","#FFFF90","#FFFF99","#FFF2AE","#FEF2E0"))+
  labs(x = "Epoch", y = "FSpe")+
  theme_minimal() +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size= 6.5, color= "black"),
        axis.title= element_text(size= 8), panel.background= element_rect(fill= "white"))+
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent")) +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.ticks.x = element_blank()) +
  theme(axis.text.x = element_blank())+
  coord_geo(
    dat = epochs_custom, pos = "bottom", expand = TRUE,
    skip = NULL, abbrv = TRUE, dat_is_discrete = TRUE,
    size = 3
  )


# Plot everything together
Fig_X_full <- plot_grid(FE_full_Recent,
                   FRed_full_Recent,FOred_full_Recent,
                   FRic_full_Recent,
                   FOri_full_Recent,FSpe_full_Recent,
                   labels= c("(a)","(b)","(c)","(d)","(e)","(f)"), 
                   label_size = 10,align = "hv", label_fontface = "bold",  nrow=6)

## Now do loop analyses based on resampling Recent sample
# Set up model to run on n iterations
registerDoParallel(cores = 10)
getDoParWorkers()
Taxon.res<-NULL
Taxon.res<- lapply(1:1000, function(x){
  # set sample of Recent
  baskets.varRes <- baskets.range.taxon %>%
    mutate(Recent = ifelse(Recent == "1", sample(c("0", "1"), size = 114, replace = TRUE), Recent))
  
  baskets.varRes <- baskets.varRes %>%
    remove_rownames %>% 
    column_to_rownames(var="Taxon_corrected") %>% 
    select(-Current_status) %>% 
    filter(rowSums(. != "0") > 0)
  
  # Transpose to form final occurrence matrix
  baskets_sharks_weights.varRes <- t(baskets.varRes)
  
  # Filter rows in FTU where Taxon_corrected is present in baskets
  valid_taxa <- unique(rownames(baskets.varRes))
  FTU.varRes <- FTU %>%
    filter(Taxon_corrected %in% valid_taxa)
  
  # Order FTUs alphabetically
  FTU.varRes <- FTU.varRes[order(FTU.varRes$Taxon_corrected),]
  
  # Form sliced species-trait matrix
  sharks_traits.varRes <- FTU.varRes %>%
    group_by(Taxon_corrected) %>%
    slice_sample(n=1, replace = FALSE) %>% 
    remove_rownames %>% 
    column_to_rownames(var="Taxon_corrected")
  
  # Transpose to form final occurrence matrix
  baskets_sharks_weights.varRes <- t(baskets.varRes)
  
  # Make matrix and numeric
  baskets_sharks_weights.varRes <- data.matrix(baskets_sharks_weights.varRes, rownames.force = NA)
  class(baskets_sharks_weights.varRes) <- "numeric"
  
  # Summarise & extract occurrences
  asb_sp_sharks_summ.varRes <- mFD::asb.sp.summary(asb_sp_w = baskets_sharks_weights.varRes)
  
  asb_sp_sharks_occ.varRes <- asb_sp_sharks_summ.varRes$asb_sp_occ
  
  ############################################################################################
  # Make trait_cat matrix
  sharks_traits_cat.varRes <- tibble(trait_name = c("CH", "CW", "CE", "LC", "XO", "LO"),
                                     trait_type = c("O", "O", "N", "N", "N", "N"),
                                     trait_weight = c(0.5,0.5,1,1,0.33,0.67))
  
  ############################################################################################
  # Form the functional entities
  sp_to_fe_sharks.varRes <- mFD::sp.to.fe(
    sp_tr       = sharks_traits.varRes, 
    tr_cat      = sharks_traits_cat.varRes, 
    fe_nm_type  = "fe_rank", 
    check_input = TRUE)
  
  # Calculate and display metrics
  alpha_fd_fe_sharks.varRes <- alpha.fd.fe(asb_sp_sharks_occ.varRes,
                                           sp_to_fe_sharks.varRes,
                                           ind_nm = c("fred", "fored", "fvuln"),
                                           check_input = TRUE,
                                           details_returned = TRUE)
  
  # Form dataframe
  FEmetrics.varRes <- as.data.frame(alpha_fd_fe_sharks.varRes$asb_fdfe) %>% 
    tibble::rownames_to_column("Epoch")
  FEmetrics.varRes$Epoch <- ordered(FEmetrics.varRes$Epoch, levels=c("Palaeocene","Eocene","Oligocene","Miocene",
                                                                     "Pliocene","Pleistocene","Recent"))
  
  # Construct trait distance matrix
  sp_dist_sharks.varRes <- mFD::funct.dist(
    sp_tr         = sharks_traits.varRes,
    tr_cat        = sharks_traits_cat.varRes,
    metric        = "gower",
    scale_euclid  = "scale_center",
    ordinal_var   = "classic",
    weight_type   = "equal",
    stop_if_NA    = TRUE)
  
  # Assess quality
  fspaces_quality_sharks.varRes <- mFD::quality.fspaces(
    sp_dist             = sp_dist_sharks.varRes,
    maxdim_pcoa         = 10,
    deviation_weighting = c("absolute","squared"),
    fdist_scaling       = FALSE,
    fdendro             = "average")
  
  # Return coordinates of each axis
  sp_faxes_coord_sharks.varRes <- fspaces_quality_sharks.varRes$"details_fspaces"$"sp_pc_coord"
  
  # Form alpha diversity indices
  alpha_fd_indices_sharks.varRes <- mFD::alpha.fd.multidim(
    sp_faxes_coord   = sp_faxes_coord_sharks.varRes[ , c("PC1", "PC2", "PC3")],
    asb_sp_w         = baskets_sharks_weights.varRes,
    ind_vect         = c("fric", "fori", "fspe"),
    scaling          = TRUE,
    check_input      = TRUE,
    details_returned = TRUE)
  
  fd_ind_values_sharks.varRes <- alpha_fd_indices_sharks.varRes$"functional_diversity_indices"
  
  ## Form dataframe
  FD.varRes <- as.data.frame(fd_ind_values_sharks.varRes) %>% 
    tibble::rownames_to_column("Epoch")
  FD.varRes$Epoch <- ordered(FD.varRes$Epoch, levels=c("Palaeocene","Eocene","Oligocene","Miocene",
                                                       "Pliocene","Pleistocene","Recent"))
  
  # Form list to merge datasets
  FDind.varRes = list(FEmetrics.varRes,FD.varRes)
  
  # Output
  FDind.varRes %>% 
    reduce(inner_join, by = "Epoch")
  
  #close loop
})

# Merge lists into 1 dataframe
res_Taxonvar_resamp_list <- list(Taxon.res)

res_Taxonvar_resamp_df <- res_Taxonvar_resamp_list %>% 
  bind_rows()

# Remove sp_richn, which is the same as nb_sp - redundant
res_Taxonvar_resamp_df$sp_richn<-NULL

# Format res for FD metrics
Recent_resampled<- res_Taxonvar_resamp_df %>% 
  select(Epoch:fspe)
save(Recent_resampled, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 2. FD changes over time/Analyses/Current Analyses/R codes/Taxon code/Data/Recent_resampling.RData")

# Melt data & save
Recent_resampled_long<- melt(Recent_resampled, id.vars= "Epoch")
save(Recent_resampled_long, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 2. FD changes over time/Analyses/Current Analyses/R codes/Taxon code/Data/Recent_resampling_long.RData")

Recent_resamp_var <-Recent_resampled %>% 
  group_by(Epoch) %>%
  summarise(FE_med = median(nb_fe),
            FE_sd = sd(nb_fe),
            Red_med = median(fred),
            Red_sd = sd(fred),
            Ored_med = median(fored),
            Ored_sd = sd(fored),
            FRic_med = median(fric),
            FRic_sd = sd(fric),
            Fori_med = median(fori),
            Fori_sd = sd(fori),
            Fspe_med = median(fspe),
            Fspe_sd = sd(fspe))

Recent_resamp_var$Epoch <- ordered(Recent_resamp_var$Epoch, levels=c("Palaeocene","Eocene","Oligocene","Miocene",
                                                     "Pliocene","Pleistocene","Recent"))

# Plot resampling in violin plots
TR_Recent_Taxonvar <- Recent_resampled_long %>% 
  filter(variable == "nb_sp")
FEmetrics_Recent_Taxonvar <- Recent_resampled_long %>% 
  filter(variable == "nb_fe")
FRedmetrics_Recent_Taxonvar <- Recent_resampled_long %>% 
  filter(variable == "fred")
FOredmetrics_Recent_Taxonvar <- Recent_resampled_long %>% 
  filter(variable == "fored")
FVulnmetrics_Recent_Taxonvar <- Recent_resampled_long %>% 
  filter(variable == "fvuln")
FRicmetrics_Recent_Taxonvar <- Recent_resampled_long %>% 
  filter(variable == "fric")
FOrimetrics_Recent_Taxonvar <- Recent_resampled_long %>% 
  filter(variable == "fori")
FSpemetrics_Recent_Taxonvar <- Recent_resampled_long %>% 
  filter(variable == "fspe")

## Produce geological time scale (for FSpe plots)
epochs_custom <- data.frame(
  name = c("Palaeocene","Eocene","Oligocene","Miocene","Pliocene","Pleistocene","Recent"), 
  max_age = c(1.5,2.5,3.5,4.5,5.5,6.5,7.6),
  min_age = c(0.4,1.5,2.5,3.5,4.5,5.5,6.5), 
  color = c("#FBA75F","#FDB46C","#FDC07A","#FFFF90","#FFFF99","#FFF2AE","#FEF2E0")
)

## Functional entities
FE_Recent_resamp <- ggplot()+
  scale_x_discrete(limits=c("Palaeocene","Eocene","Oligocene","Miocene", "Pliocene", "Pleistocene", "Recent"))+
  geom_boxplot(data = FEmetrics_Recent_Taxonvar, aes(x = Epoch, y = value),
               fill=c("#FBA75F","#FDB46C","#FDC07A","#FFFF90","#FFFF99","#FFF2AE","#FEF2E0"))+
  labs(x = "", y = "# FEs")+
  theme_minimal() +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size= 6.5, color= "black"),
        axis.title= element_text(size= 8), panel.background= element_rect(fill= "white"))+
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent")) +
  theme(axis.text.x = element_blank())

## Functional redundancy
FRed_Recent_resamp <- ggplot()+
  scale_x_discrete(limits=c("Palaeocene","Eocene","Oligocene","Miocene", "Pliocene", "Pleistocene", "Recent"))+
  geom_boxplot(data = FRedmetrics_Recent_Taxonvar, aes(x = Epoch, y = value),
               fill=c("#FBA75F","#FDB46C","#FDC07A","#FFFF90","#FFFF99","#FFF2AE","#FEF2E0"))+
  labs(x = "", y = "FRed")+
  theme_minimal() +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size= 6.5, color= "black"),
        axis.title= element_text(size= 8), panel.background= element_rect(fill= "white"))+
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent")) +
  theme(axis.text.x = element_blank())

## Functional over-redundancy
FOred_Recent_resamp <- ggplot()+
  scale_x_discrete(limits=c("Palaeocene","Eocene","Oligocene","Miocene", "Pliocene", "Pleistocene", "Recent"))+
  geom_boxplot(data = FOredmetrics_Recent_Taxonvar, aes(x = Epoch, y = value),
               fill=c("#FBA75F","#FDB46C","#FDC07A","#FFFF90","#FFFF99","#FFF2AE","#FEF2E0"))+
  labs(x = "", y = "FOred")+
  theme_minimal() +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size= 6.5, color= "black"),
        axis.title= element_text(size= 8), panel.background= element_rect(fill= "white"))+
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent")) +
  theme(axis.text.x = element_blank())

## Functional vulnerability
FVuln_Recent_resamp <- ggplot()+
  scale_x_discrete(limits=c("Palaeocene","Eocene","Oligocene","Miocene", "Pliocene", "Pleistocene", "Recent"))+
  geom_boxplot(data = FVulnmetrics_Recent_Taxonvar, aes(x = Epoch, y = value),
               fill=c("#FBA75F","#FDB46C","#FDC07A","#FFFF90","#FFFF99","#FFF2AE","#FEF2E0"))+
  labs(x = "", y = "FVuln")+
  theme_minimal() +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size= 6.5, color= "black"),
        axis.title= element_text(size= 8), panel.background= element_rect(fill= "white"))+
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent")) +
  theme(axis.text.x = element_blank())

## Functional richness
FRic_Recent_resamp <- ggplot()+
  scale_x_discrete(limits=c("Palaeocene","Eocene","Oligocene","Miocene", "Pliocene", "Pleistocene", "Recent"))+
  geom_boxplot(data = FRicmetrics_Recent_Taxonvar, aes(x = Epoch, y = value),
               fill=c("#FBA75F","#FDB46C","#FDC07A","#FFFF90","#FFFF99","#FFF2AE","#FEF2E0"))+
  labs(x = "", y = "FRic")+
  theme_minimal() +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size= 6.5, color= "black"),
        axis.title= element_text(size= 8), panel.background= element_rect(fill= "white"))+
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent")) +
  theme(axis.text.x = element_blank())

## Functional originality
FOri_Recent_resamp <- ggplot()+
  scale_x_discrete(limits=c("Palaeocene","Eocene","Oligocene","Miocene", "Pliocene", "Pleistocene", "Recent"))+
  geom_boxplot(data = FOrimetrics_Recent_Taxonvar, aes(x = Epoch, y = value),
               fill=c("#FBA75F","#FDB46C","#FDC07A","#FFFF90","#FFFF99","#FFF2AE","#FEF2E0"))+
  labs(x = "", y = "FOri")+
  theme_minimal() +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size= 6.5, color= "black"),
        axis.title= element_text(size= 8), panel.background= element_rect(fill= "white"))+
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent")) +
  theme(axis.text.x = element_blank())

## Functional specialisation
FSpe_Recent_resamp <- ggplot()+
  scale_x_discrete(limits=c("Palaeocene","Eocene","Oligocene","Miocene", "Pliocene", "Pleistocene", "Recent"))+
  geom_boxplot(data = FSpemetrics_Recent_Taxonvar, aes(x = Epoch, y = value),
               fill=c("#FBA75F","#FDB46C","#FDC07A","#FFFF90","#FFFF99","#FFF2AE","#FEF2E0"))+
  labs(x = "", y = "FSpe")+
  theme_minimal() +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size= 6.5, color= "black"),
        axis.title= element_text(size= 8), panel.background= element_rect(fill= "white"))+
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent")) +
  theme(axis.ticks.x = element_blank())+
  theme(axis.text.x = element_blank())+
  coord_geo(
    dat = epochs_custom, pos = "bottom", expand = TRUE,
    skip = NULL, abbrv = TRUE, dat_is_discrete = TRUE,
    size = 3
  )


# Plot everything together  - produces Figure 2
Fig_X <- plot_grid(FE_Recent_resamp,
                   FRed_Recent_resamp,FOred_Recent_resamp,
                   FRic_Recent_resamp,
                   FOri_Recent_resamp,FSpe_Recent_resamp,
                   labels= c("(a)","(b)","(c)","(d)","(e)","(f)"), 
                   label_size = 10,align = "hv", label_fontface = "bold",  nrow=6)
