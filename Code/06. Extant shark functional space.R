#########################################################################################################################
# 06. Extant shark functional space
## This R code provides how functional space is formed for all living sharks to test range of extant shark data
## it produces Figure S9
#########################################################################################################################

## Import packages
library(readxl)
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

# Load tooth data & isolate extant shark species
load(file="~/Cleaned data.RData")

extant_sharks <- data %>%
  select(Taxon_corrected,Current_status,Order,Family,Genus,Species) %>% 
  group_by(Taxon_corrected) %>% 
  filter(Current_status == "Extant")

extant_sharks <- extant_sharks %>% 
  distinct(Taxon_corrected,Current_status,Order,Family,Genus,Species)
extant_sharks <- extant_sharks[order(extant_sharks$Taxon_corrected),]

# Count Orders, families, genera 
n_distinct(extant_sharks$Order)
## 9 - 100% of 9 (Ebert et al. 2021)
n_distinct(extant_sharks$Family)
## 32 - 84.2% of 38 (Ebert et al. 2021)
n_distinct(extant_sharks$Genus)
## 68 - 64.2% of 106 (Ebert et al. 2021)

extant_species <- extant_sharks %>%
  ungroup() %>% 
  filter(Species!="sp.")

# Load all living species & traits
living_sharks <- read_xlsx("~/Data S2.xlsx")

living_sharks$Body_size <- ordered(living_sharks$Body_size, 
                                   levels=c("Small","Medium","Large","Giant"))

living_sharks$Prey_preference <- as.factor(living_sharks$Prey_preference)
living_sharks$Feeding_mechanism <- as.factor(living_sharks$Feeding_mechanism)
living_sharks$size <- as.numeric(living_sharks$size)
# Form species-trait matrix from living_sharks
sharks_traits <- living_sharks %>% 
  select(Species,Body_size,Prey_preference,Feeding_mechanism) %>% 
  remove_rownames %>% 
  column_to_rownames(var="Species")

# Form occurrence matrix
Occ_data <- living_sharks %>% 
  select(Species)

# Mark presence or absence in tooth data
Occ_data$Data <- as.numeric(Occ_data$Species %in% extant_species$Taxon_corrected)
Occ_data$FullFRic <- 1 #Done to ensure every species is present in an assemblage for coding; & to mark all as extant for fuse function

baskets <- Occ_data %>% 
  remove_rownames %>% 
  column_to_rownames(var="Species")

baskets_sharks_weights <- t(baskets)

# Form trait category matrix
sharks_traits_cat <- tibble(trait_name = c("Body_size","Prey_preference","Feeding_mechanism"),
                            trait_type = c("O", "N", "N"))
# Summarise dataset - mark stop_if_NA as false due to 13% of species missing feeding mechanism data
sharks_traits_summ <- mFD::sp.tr.summary(
  tr_cat     = sharks_traits_cat,   
  sp_tr      = sharks_traits, 
  stop_if_NA = FALSE)

# Reveal number of states for each dental character
sharks_traits_summ$tr_summary_list

# make assemblage data into a matrix & ensure all matrix values are numeric
baskets_sharks_weights <- data.matrix(baskets_sharks_weights, rownames.force = NA)
class(baskets_sharks_weights) <- "numeric"

# summarise dataset
asb_sp_sharks_summ <- mFD::asb.sp.summary(asb_sp_w = baskets_sharks_weights)

# retrieve species occurrences for all assemblages - FE calculations not possible with NAs
asb_sp_sharks_occ <- asb_sp_sharks_summ$asb_sp_occ

# Construct trait distance matrix using species
sp_dist_sharks <- mFD::funct.dist(
  sp_tr         = sharks_traits,
  tr_cat        = sharks_traits_cat,
  metric        = "gower",
  scale_euclid  = "scale_center",
  ordinal_var   = "classic",
  weight_type   = "equal",
  stop_if_NA    = FALSE)

# Assess quality
fspaces_quality_sharks <- mFD::quality.fspaces(
  sp_dist             = sp_dist_sharks,
  maxdim_pcoa         = 10,
  deviation_weighting = c("absolute","squared"),
  fdist_scaling       = FALSE,
  fdendro             = "average")

# Check pocket - data best represented by 5 dimensions/axes; but what about reduced to 4 & 3 dimensions?
round(fspaces_quality_sharks$"quality_fspaces", 5)

# Plot quality
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

# Check % of inertia accounted for by axes - reveals that 4 axes accounts for 73.25% inertia
eig <- fspaces_quality_sharks$details_fspaces$pc_eigenvalues %>% 
  as.data.frame()

eig$variance <- (eig$Eigenvalues/sum(eig$Eigenvalues))*100

# Return coordinates of each axis - correlation checks not possible with NAs
sp_faxes_coord_sharks <- fspaces_quality_sharks$"details_fspaces"$"sp_pc_coord"

## Form functional space of first 4 axes
big_plot <- mFD::funct.space.plot(
  sp_faxes_coord  = sp_faxes_coord_sharks[ , c("PC1", "PC2", "PC3", "PC4")],
  faxes           = c("PC1", "PC2", "PC3", "PC4"),
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
  sp_faxes_coord   = sp_faxes_coord_sharks[ , c("PC1", "PC2", "PC3", "PC4")],
  asb_sp_w         = baskets_sharks_weights,
  ind_vect         = c("fric", "fori", "fspe"),
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = TRUE)

fd_ind_values_sharks <- alpha_fd_indices_sharks$"functional_diversity_indices"

FDmetrics <- fd_ind_values_sharks %>% 
  as.data.frame()

x_limits_PC1 <- c(-0.5, 0.5)
y_limits_PC12 <- c(-0.8, 0.25) 

## Species in tooth data
plots_data <- mFD::alpha.multidim.plot(
  output_alpha_fd_multidim = alpha_fd_indices_sharks,
  plot_asb_nm              = "Data",
  ind_nm                   = c("fric","fori", "fspe"),
  faxes                    = NULL,
  faxes_nm                 = NULL,
  range_faxes              = c(NA, NA),
  color_bg                 = "white",
  shape_sp                 = c(pool = 21, asb1 = 21),
  size_sp                  = c(pool = 1.75, asb1 = 1),
  color_sp                 = c(pool = "black", asb1 = "#00AFBB"),
  color_vert               = c(pool = "black", asb1 = "#00AFBB"),
  fill_sp                  = c(pool = "black", asb1 = "#00AFBB"),
  fill_vert                = c(pool = "black", asb1 = "#00AFBB"),
  color_ch                 = c(pool = "black", asb1 = "#00AFBB"),
  fill_ch                  = c(pool = "white", asb1 = "#00AFBB"),
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

data_space <- plots_data$"fric"$PC1_PC2+
  labs(x = "PCoA1", y = "PCoA2")+
  coord_fixed(ratio = 1, xlim = x_limits_PC1, ylim = y_limits_PC12)

# Calculate highest ranked FOri & FSpe species to see if in the data
# Calculate FUSE for each species - set neighbours to 1 to calculate originality instead of uniqueness
FUSE.extant <- fuse(sp_dist        = sp_dist_sharks, 
                    sp_faxes_coord = as.matrix(sp_faxes_coord_sharks), 
                    nb_NN          = 1,  
                    GE             = Occ_data$FullFRic,
                    standGE        = TRUE)

# Add presence/absence in tooth data
FUSE.extant$Data <- Occ_data$Data
FUSE.extant$Data <- as.factor(FUSE.extant$Data)

status_colors <- c("0" = "black", "1" = "#00AFBB")

# Top 5% contributors - FSpe; species that drive FRic
## Rank contributors (top 5% = 27 spp.)
FUSE.extant <- FUSE.extant %>% 
  rownames_to_column(var = "Species")

FUSE_fspe_rank <- FUSE.extant %>%
  select(Species, FSp_std, Data) %>%
  arrange(desc(FSp_std)) %>%
  head(50)
FUSE_fspe_rank$Species <- factor(FUSE_fspe_rank$Species, levels = FUSE_fspe_rank$Species[order(FUSE_fspe_rank$FSp_std, decreasing = TRUE)])

# Create lollipop plot for FSpe
plot_FSpe_rank_lollipop <- ggplot(FUSE_fspe_rank, aes(x = Species, y = FSp_std, fill = Data, color = Data)) +
  geom_segment(aes(xend = Species, yend = 0), size = 1.5) +
  geom_point(size = 3) +
  coord_flip() +
  labs(x = "", y = "FSpe") +
  scale_fill_manual(values = status_colors, guide = guide_legend(title = "Status")) +
  scale_color_manual(values = status_colors, guide = guide_legend(title = "Status")) +
  theme_minimal() +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text.y = element_text(size = 10, color = "black", face = "italic"),  # Set y-axis labels to italic
        axis.title = element_text(size = 12), panel.background = element_rect(fill = "white")) +
  theme(panel.grid = element_blank())

# Plot everything - Produces Figure S9
plot_grid(
  data_space,
  plot_FSpe_rank_lollipop,
  labels = c("(a)", "(b)"), 
  label_size = 12, align = "hv", label_fontface = "bold", nrow = 1)
