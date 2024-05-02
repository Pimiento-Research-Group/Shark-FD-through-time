#########################################################################################################################
# 03. Functional space
## This R code provides how functional space is formed
## it produces Figure 1, S2-S4; and Table S3
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
library(deeptime)
library(grid)
library(gridExtra)

# Load data
load(file="~/Cleaned data.RData")

# Functional spaces per epoch
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
  mutate(Recent = case_when(Recent == "1" && Taxon_corrected == "Galeocerdo sp."~"0", TRUE ~ Recent)) %>% 
  mutate(Recent = case_when(Recent == "1" && Taxon_corrected == "Galeorhinus sp."~"0", TRUE ~ Recent)) %>% 
  mutate(Recent = case_when(Recent == "1" && Taxon_corrected == "Galeus sp."~"0", TRUE ~ Recent)) %>% 
  mutate(Recent = case_when(Recent == "1" && Taxon_corrected == "Hemipristis sp."~"0", TRUE ~ Recent)) %>% 
  mutate(Recent = case_when(Recent == "1" && Taxon_corrected == "Heptranchias sp."~"0", TRUE ~ Recent)) %>% 
  mutate(Recent = case_when(Recent == "1" && Taxon_corrected == "Isurus sp."~"0", TRUE ~ Recent)) %>%
  mutate(Recent = case_when(Recent == "1" && Taxon_corrected == "Megachasma sp."~"0", TRUE ~ Recent)) %>% 
  mutate(Recent = case_when(Recent == "1" && Taxon_corrected == "Nebrius sp."~"0", TRUE ~ Recent)) %>% 
  mutate(Recent = case_when(Recent == "1" && Taxon_corrected == "Notorynchus sp."~"0", TRUE ~ Recent))

baskets.range <- baskets.range %>% 
  rowwise() %>% 
  mutate(Current_status = case_when(Taxon_corrected == "Alopias sp."~"Extinct", TRUE ~ Current_status)) %>%
  mutate(Current_status = case_when(Taxon_corrected == "Carcharias sp."~"Extinct", TRUE ~ Current_status)) %>% 
  mutate(Current_status = case_when(Taxon_corrected == "Carcharodon sp."~"Extinct", TRUE ~ Current_status)) %>% 
  mutate(Current_status = case_when(Taxon_corrected == "Dalatias sp."~"Extinct", TRUE ~ Current_status)) %>% 
  mutate(Current_status = case_when(Taxon_corrected == "Galeocerdo sp."~"Extinct", TRUE ~ Current_status)) %>% 
  mutate(Current_status = case_when(Taxon_corrected == "Galeorhinus sp."~"Extinct", TRUE ~ Current_status)) %>% 
  mutate(Current_status = case_when(Taxon_corrected == "Galeus sp."~"Extinct", TRUE ~ Current_status)) %>% 
  mutate(Current_status = case_when(Taxon_corrected == "Hemipristis sp."~"Extinct", TRUE ~ Current_status)) %>% 
  mutate(Current_status = case_when(Taxon_corrected == "Heptranchias sp."~"Extinct", TRUE ~ Current_status)) %>% 
  mutate(Current_status = case_when(Taxon_corrected == "Isurus sp."~"Extinct", TRUE ~ Current_status)) %>%
  mutate(Current_status = case_when(Taxon_corrected == "Megachasma sp."~"Extinct", TRUE ~ Current_status)) %>% 
  mutate(Current_status = case_when(Taxon_corrected == "Nebrius sp."~"Extinct", TRUE ~ Current_status)) %>% 
  mutate(Current_status = case_when(Taxon_corrected == "Notorynchus sp."~"Extinct", TRUE ~ Current_status))

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

# Save excluded taxa
save(selected_taxa,file = "~/Filtered taxa.RData")

# Filter taxa that are NOT present in "selected_taxa" from sample
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

# Filter taxa that are NOT present in "selected_taxa" from sample
Av_sharks_filtered <- Av_sharks %>%
  anti_join(selected_taxa, by = "Taxon_corrected")

# Save species-trait matrix for use in randomisation simulations (Code S6)
save(Av_sharks_filtered,file = "~/Species-trait matrix.RData")


# Make species-trait matrix
sharks_traits <- Av_sharks_filtered %>%
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

# Number of functional entities - 64 possible functional entities per iteration (2.96%)
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

mad_df <- round(fspaces_quality_sharks$"quality_fspaces", 5) %>% 
  as.data.frame()
plot(mad_df$mad, xlab = "Axes", ylab = "mad")

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

# Check % of inertia accounted for by axes - reveals that 3 axes accounts for 78.4% inertia
eig <- fspaces_quality_sharks$details_fspaces$pc_eigenvalues %>% 
  as.data.frame()

eig$variance <- (eig$Eigenvalues/68.09946)*100

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

FD_df <- left_join(FEmetrics_taxon, FDindices_taxon, by = "Epoch")


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

## Plot space with highest FOri & FSpe (see Code S13 for FOri/FSpe calculations)
Pal_contributors12 <- plots_Pal$fric$PC1_PC2 +
  labs(x = "PCoA1", y = "PCoA2")+
  annotate("point", x = -0.236591721, y = 0.04938241, color = "#00AFBB", size = 3)+
  annotate("point", x = 0.071143081, y = 0.75269101, color = "#FC4E07", size = 3)+
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

## Plot space with highest FOri & FSpe (see Code S13 for FOri/FSpe calculations)
Eo_contributors12 <- plots_Eo$fric$PC1_PC2 +
  labs(x = "PCoA1", y = "PCoA2")+
  annotate("point", x = 0.330376859, y = -0.24171603, color = "#00AFBB", size = 3)+
  annotate("point", x = 0.075912198, y = 0.871759751, color = "#FC4E07", size = 3)+
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

## Plot space with highest FOri & FSpe (see Code S13 for FOri/FSpe calculations)
Oli_contributors12 <- plots_Oli$fric$PC1_PC2 +
  labs(x = "PCoA1", y = "PCoA2")+
  annotate("point", x = 0.038887820, y = 0.528810970, color = "#00AFBB", size = 3)+
  annotate("point", x = 0.074004551, y = 0.824132255, color = "#FC4E07", size = 3)+
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

## Plot space with highest FOri & FSpe (see Code S13 for FOri/FSpe calculations)
Mio_contributors12 <- plots_Mio$fric$PC1_PC2 +
  labs(x = "PCoA1", y = "PCoA2")+
  annotate("point", x = -0.235718809, y = 0.22022928, color = "#00AFBB", size = 3)+
  annotate("point", x = 0.074004551, y = 0.824132255, color = "#FC4E07", size = 3)+
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

## Plot space with highest FOri & FSpe (see Code S13 for FOri/FSpe calculations)
Plio_contributors12 <- plots_Plio$fric$PC1_PC2 +
  labs(x = "PCoA1", y = "PCoA2")+
  annotate("point", x = -0.235718809, y = 0.22022928, color = "#00AFBB", size = 3)+
  annotate("point", x = 0.573863865, y = 0.300282886, color = "#FC4E07", size = 3)+
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

## Plot space with highest FOri & FSpe (see Code S13 for FOri/FSpe calculations)
Ple_contributors12 <- plots_Ple$fric$PC1_PC2 +
  labs(x = "PCoA1", y = "PCoA2")+
  annotate("point", x = 0.039169586, y = -0.118377674, color = "#00AFBB", size = 3)+
  annotate("point", x = 0.038887820, y = 0.52881097, color = "#FC4E07", size = 3)+
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

## Plot space with highest FOri & FSpe (see Code S13 for FOri/FSpe calculations)
Rec_contributors12 <- plots_Rec$fric$PC1_PC2 +
  labs(x = "PCoA1", y = "PCoA2")+
  annotate("point", x = 0.163877685, y = -0.237613103, color = "#00AFBB", size = 3)+
  annotate("point", x = 0.456052215, y = 0.161540500, color = "#FC4E07", size = 3)+
  coord_fixed(ratio = 1, xlim = x_limits_PC1, ylim = y_limits_PC12)

# Functional space extinct vs extant sharks - Figure 1i
Av_data <- data %>%
  group_by(Taxon_corrected,Current_status) %>% 
  reframe(CH = mean(CH_mm),
          CW = mean(CW_mm),
          CE = getmode(CE),
          LC = getmode(LC),
          XO = getmode(XO),
          LO = getmode(LO))

Av_sharks <- Av_data %>% 
  rowwise %>% 
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
Av_filtered <- Av_sharks %>%
  anti_join(selected_taxa, by = "Taxon_corrected")

sharks_traits <- Av_filtered %>%
  remove_rownames %>% 
  column_to_rownames(var="Taxon_corrected") %>% 
  select(-Current_status)

# Occurrence matrix
baskets.all.taxa <- Av_filtered %>% 
  select(Taxon_corrected,Current_status)

filtered_baskets.taxa <- baskets.all.taxa %>%
  anti_join(selected_taxa, by = "Taxon_corrected")

filtered_baskets.taxa$Extinct <- 0
filtered_baskets.taxa$Extant <- 0

filtered_baskets.taxa <- filtered_baskets.taxa %>% 
  rowwise() %>% 
  mutate(Extinct = case_when(Current_status == "Extinct"~1,TRUE~Extinct),
         Extant = case_when(Current_status == "Extant"~1,TRUE~Extant))

baskets.space <- filtered_baskets.taxa %>% 
  column_to_rownames(var="Taxon_corrected") %>% 
  select(-Current_status)
baskets.taxa <- t(baskets.space)

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
baskets_weights <- data.matrix(baskets.taxa, rownames.force = NA)
class(baskets_weights) <- "numeric"

# summarise dataset
asb_sp_sharks_summ <- mFD::asb.sp.summary(asb_sp_w = baskets_weights)

# retrieve species occurrences for all assemblages (epochs):
asb_sp_sharks_occ <- asb_sp_sharks_summ$asb_sp_occ

# Form the functional entities
sp_to_fe_sharks <- mFD::sp.to.fe(
  sp_tr       = sharks_traits, 
  tr_cat      = sharks_traits_cat, 
  fe_nm_type  = "fe_rank", 
  check_input = TRUE)

# Calculate and display metrics
alpha_fd_fe_sharks <- alpha.fd.fe(asb_sp_sharks_occ,
                                  sp_to_fe_sharks,
                                  ind_nm = c("fred", "fored", "fvuln"),
                                  check_input = TRUE,
                                  details_returned = TRUE)

# Form dataframe of results for plotting
FEmetrics_taxon <- as.data.frame(alpha_fd_fe_sharks$asb_fdfe) %>% 
  tibble::rownames_to_column("Status")

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

# Return coordinates of each axis
sp_faxes_coord_sharks <- fspaces_quality_sharks$"details_fspaces"$"sp_pc_coord"

# Compute the range of functional axes:
range_sp_coord  <- range(sp_faxes_coord_sharks)

# Based on the range of species coordinates values, compute a nice range ...
# ... for functional axes:
range_faxes <- range_sp_coord +
  c(-1, 1) * (range_sp_coord[2] - range_sp_coord[1]) * 0.05
range_faxes

# get species coordinates:
sp_faxes_coord_xy <- sp_faxes_coord_sharks[, c("PC1", "PC2")]
sp_faxes_coord_xz <- sp_faxes_coord_sharks[, c("PC1", "PC3")]

# Form indices
alpha_fd_indices_sharks_space <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = sp_faxes_coord_sharks[ , c("PC1", "PC2", "PC3")],
  asb_sp_w         = baskets_weights,
  ind_vect         = c("fric", "fori", "fspe"),
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = TRUE)

# Plot space
plots_space <- mFD::alpha.multidim.plot(
  output_alpha_fd_multidim = alpha_fd_indices_sharks_space,
  plot_asb_nm              = c("Extinct", "Extant"),
  ind_nm                   = c("fric","fori", "fspe"),
  faxes                    = NULL,
  faxes_nm                 = NULL,
  range_faxes              = c(NA, NA),
  color_bg                 = "white",
  shape_sp                 = c(pool = 21, asb1 = 21, asb2 = 21),
  size_sp                  = c(pool = 0.7, asb1 = 1, asb2 = 1),
  color_sp                 = c(pool = "black", asb1 = "#0073C2FF", asb2 = "#868686FF"),
  color_vert               = c(pool = "black", asb1 = "#0073C2FF", asb2 = "#868686FF"),
  fill_sp                  = c(pool = NA, asb1 = "#0073C2FF", asb2 = "#868686FF"),
  fill_vert                = c(pool = NA, asb1 = "#0073C2FF", asb2 = "#868686FF"),
  color_ch                 = c(pool = NA, asb1 = "#0073C2FF", asb2 = "#868686FF"),
  fill_ch                  = c(pool = "white", asb1 = "#0073C2FF", asb2 = "#868686FF"),
  alpha_ch                 = c(pool = 1, asb1 = 0.3, asb2 = 0.3),
  shape_centroid_fdis      = c(asb1 = 22,  asb2 = 24),
  shape_centroid_fdiv      = c(asb1 = 22,  asb2 = 24),
  shape_centroid_fspe      = 23,
  color_centroid_fspe      = "black",
  size_sp_nm               = 3, 
  color_sp_nm              = "black",
  plot_sp_nm               = NULL,
  fontface_sp_nm           = "plain",
  save_file                = FALSE,
  check_input              = TRUE) 

# Whole spaces
whole_space_a <- plots_space$"fric"$PC1_PC2+
  labs(x = "PCoA1", y = "PCoA2")+
  coord_fixed(ratio = 1, xlim = x_limits_PC1, ylim = y_limits_PC12)

whole_space_b <- plots_space$"fric"$PC1_PC3+
  labs(x = "PCoA1", y = "PCoA3")+
  coord_fixed(ratio = 1, xlim = x_limits_PC1, ylim = y_limits_PC13)

# Plot all spaces on first 2 axes - produces Figure 1b-i
Spaces <- plot_grid(Pal_contributors12,Eo_contributors12,Oli_contributors12,Mio_contributors12,
                    Plio_contributors12,Ple_contributors12,Rec_contributors12,whole_space_a,
                    labels = c("(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)"), 
                    label_size = 12,align = "hv", label_fontface = "bold", hjust = -0.15,  nrow=3)

# Supplement PC1-PC3 spaces (See Code S13 for FOri/FSpe values & coordinates) - produces Figure S4
## PCoA1_PCoA3
# Palaeocene
Pal_contributors_13 <- plots_Pal$fric$PC1_PC3 +
  labs(x = "PCoA1", y = "PCoA3")+
  annotate("point", x = -0.236591721, y = -0.15123018, color = "#00AFBB", size = 4)+
  annotate("point", x = 0.071143081, y = -0.127133020, color = "#FC4E07", size = 4)+
  coord_fixed(ratio = 1, xlim = x_limits_PC1, ylim = y_limits_PC13)

# Eocene
Eo_contributors_13 <- plots_Eo$fric$PC1_PC3 +
  labs(x = "PCoA1", y = "PCoA3")+
  annotate("point", x = 0.330376859, y = -0.149695127, color = "#00AFBB", size = 4)+
  annotate("point", x = 0.075912198, y = -0.136009956, color = "#FC4E07", size = 4)+
  coord_fixed(ratio = 1, xlim = x_limits_PC1, ylim = y_limits_PC13)

# Oligocene
Oli_contributors_13 <- plots_Oli$fric$PC1_PC3 +
  labs(x = "PCoA1", y = "PCoA3")+
  annotate("point", x = 0.038887820, y = -0.369841873, color = "#00AFBB", size = 4)+
  annotate("point", x = 0.074004551, y = -0.132459181, color = "#FC4E07", size = 4)+
  coord_fixed(ratio = 1, xlim = x_limits_PC1, ylim = y_limits_PC13)

# Miocene
Mio_contributors_13 <- plots_Mio$fric$PC1_PC3 +
  labs(x = "PCoA1", y = "PCoA3")+
  annotate("point", x = -0.235718809, y = -0.007746932, color = "#00AFBB", size = 4)+
  annotate("point", x = 0.074004551, y = -0.132459181, color = "#FC4E07", size = 4)+
  coord_fixed(ratio = 1, xlim = x_limits_PC1, ylim = y_limits_PC13)

# Pliocene
Plio_contributors_13 <- plots_Plio$fric$PC1_PC3 +
  labs(x = "PCoA1", y = "PCoA3")+
  annotate("point", x = -0.235718809, y = -0.007746932, color = "#00AFBB", size = 4)+
  annotate("point", x = 0.573863865, y = 0.132176195, color = "#FC4E07", size = 4)+
  coord_fixed(ratio = 1, xlim = x_limits_PC1, ylim = y_limits_PC13)

# Pleistocene
Ple_contributors_13 <- plots_Ple$fric$PC1_PC3 +
  labs(x = "PCoA1", y = "PCoA3")+
  annotate("point", x = 0.039169586, y = 0.204724630, color = "#00AFBB", size = 4)+
  annotate("point", x = 0.038887820, y = -0.369841873, color = "#FC4E07", size = 4)+
  coord_fixed(ratio = 1, xlim = x_limits_PC1, ylim = y_limits_PC13)

# Recent
Rec_contributors_13 <- plots_Rec$fric$PC1_PC3 +
  labs(x = "PCoA1", y = "PCoA3")+
  annotate("point", x = 0.163877685, y = 0.04022443, color = "#00AFBB", size = 4)+
  annotate("point", x = 0.456052215, y = 0.235554632, color = "#FC4E07", size = 4)+
  coord_fixed(ratio = 1, xlim = x_limits_PC1, ylim = y_limits_PC13)

Spaces_sup <- plot_grid(whole_space_b,Pal_contributors_13,Eo_contributors_13,Oli_contributors_13,
                    Mio_contributors_13,Plio_contributors_13,Ple_contributors_13,Rec_contributors_13,
                    labels = c("(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)"), 
                    label_size = 12,align = "hv", label_fontface = "bold", hjust = -0.15,  nrow=2)


## Make functional space with vertices species - Figure 1a
# get species coordinates along axes:
sp_faxes_coord_xy <- sp_faxes_coord_sharks[, c("PC1", "PC2")]
sp_faxes_coord_xz <- sp_faxes_coord_sharks[, c("PC1", "PC3")]

# Plot background:
plot_12 <- mFD::background.plot(range_faxes = range_faxes, 
                               faxes_nm = c("PC1", "PC2"),
                               color_bg = "white")
plot_13 <- mFD::background.plot(range_faxes = range_faxes, 
                                faxes_nm = c("PC1", "PC3"),
                                color_bg = "white")

# Retrieve vertices coordinates along the functional axes:
vert12 <- mFD::vertices(sp_faxes_coord = sp_faxes_coord_xy,  
                      order_2D = FALSE, 
                      check_input = TRUE)
vert13 <- mFD::vertices(sp_faxes_coord = sp_faxes_coord_xz,  
                        order_2D = FALSE, 
                        check_input = TRUE)

# Plot vertices species
plot_sp_vert12 <- mFD::pool.plot(ggplot_bg = plot_12,
                               sp_coord2D = sp_faxes_coord_xy,
                               vertices_nD = vert12,
                               plot_pool = TRUE,
                               color_pool = "darkgrey",
                               fill_pool = "darkgrey",
                               alpha_ch =  0.8,
                               color_ch = "black",
                               fill_ch = "white",
                               shape_pool = 21,
                               size_pool = 1.5,
                               shape_vert = 16,
                               size_vert = 1,
                               color_vert = "darkgrey",
                               fill_vert = "darkgrey")+
  coord_fixed(ratio = 1, xlim = x_limits_PC1, ylim = y_limits_PC12)+
annotate("point", x = 0.118653287, y = -0.26891626, color = "black", size = 3) +
annotate("point", x = -0.341272802, y = -0.04882540, color = "black", size = 3) +
annotate("point", x = 0.570100534, y = 0.299540649, color = "black", size = 3) +
annotate("point", x = 0.078050620, y = 0.866069711, color = "black", size = 3)+
  labs(x = "PCoA1", y = "PCoA2")

plot_sp_vert13 <- mFD::pool.plot(ggplot_bg = plot_13,
                                 sp_coord2D = sp_faxes_coord_xz,
                                 vertices_nD = vert13,
                                 plot_pool = TRUE,
                                 color_pool = "darkgrey",
                                 fill_pool = "darkgrey",
                                 alpha_ch =  1.5,
                                 color_ch = "black",
                                 fill_ch = "white",
                                 shape_pool = 21,
                                 size_pool = 1.5,
                                 shape_vert = 16,
                                 size_vert = 1,
                                 color_vert = "darkgrey",
                                 fill_vert = "darkgrey")+
  labs(x = "PCoA1", y = "PCoA3")+
  coord_fixed(ratio = 1, xlim = x_limits_PC1, ylim = y_limits_PC13)+
annotate("point", x = 0.027526015, y = 0.249759754, color = "black", size = 3) +
annotate("point", x = -0.338562290, y = -0.118183001, color = "black", size = 3) +
annotate("point", x = 0.573863865, y = 0.132176195, color = "black", size = 3) +
annotate("point", x = -0.177841550, y = -0.378911567, color = "black", size = 3)
  

# Produce Figure 1a
Fig_1a <- plot_grid(plot_sp_vert12,plot_sp_vert13,
          labels = NA,
          label_size = 12,align = "hv", label_fontface = "bold", hjust = -0.15,  nrow=2)

# Produce near-complete Figure 1
Fig_1a_j <- plot_grid(Fig_1a,Spaces,
                    labels = c("(a)",""),
                    label_size = 12,align = "hv", label_fontface = "bold", hjust = -0.15,  nrow=1)

# Produce complete Figure 1
Fig_1 <- grid.arrange(
  Fig_1a_j,
  Spp_occ,
  ncol = 1,
  heights = c(2, 1)
)
