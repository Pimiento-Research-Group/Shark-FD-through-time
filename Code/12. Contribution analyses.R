###################################################################################################################
# 12. Contribution analyses 
## This R code provides analyses on species contribution to functional diversity
## it produces Figures 3 and S12
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
library(foreach)
library(parallel)
library(doParallel)
library(grid)

# Load data
load(file="~/Cleaned data.RData")

# Run functional diversity inputs 
# Occurrence matrix
# Form FTU unique combinations, with epoch included
abun <- data %>% 
  group_by(Taxon_corrected) %>% 
  distinct(CH,CW,CE,LC,XO,LO,Epoch,Epoch_earliest,Epoch_latest,Current_status)

# Order taxa alphabetically like before and order Epochs
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
selected_taxa <- baskets.range.taxon %>%
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
filtered_baskets.range <- baskets.range.taxon %>%
  anti_join(selected_taxa, by = "Taxon_corrected")

# Make taxa row names & tidy dataframe to form final matrix
baskets <- filtered_baskets.range %>% 
  remove_rownames %>% 
  column_to_rownames(var="Taxon_corrected") %>% 
  select(-Current_status)

# Make species richness object (as this will be constant regardless of FTU variation) to add later
baskets.epoch <- filtered_baskets.range %>% 
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

# Form functional taxonomic units (FTUs), grouping by species to account for intraspecific variability
FTU <- data %>% 
  group_by(Taxon_corrected) %>% 
  distinct(CH,CW,CE,LC,XO,LO)

# Order FTUs alphabetically to help ensure matches between species-trait matrix and occurrence matrix later
FTU <- FTU[order(FTU$Taxon_corrected),]

# Filter taxa that are NOT present in "selected_taxa"
FTU_filtered <- FTU %>%
  anti_join(selected_taxa, by = "Taxon_corrected")

## Calculate FOri and FSpe contributions based on 1000 iterations of species
# Set up GE
GE <- filtered_baskets.range %>% 
  select(Taxon_corrected, Current_status) %>% 
  distinct(Taxon_corrected,Current_status) %>%  
  mutate(Current_status = case_when(Current_status == "Extant" ~ "1", TRUE ~ Current_status)) %>% 
  mutate(Current_status = case_when(Current_status == "Extinct" ~ "0", TRUE ~ Current_status))

GE <- GE[order(GE$Taxon_corrected),]
GE$Current_status <- as.numeric(GE$Current_status)

# Set up model to run on n iterations
registerDoParallel(cores = 10)
getDoParWorkers()
FUSE_res<-NULL
FUSE_res<- lapply(1:1000, function(x){
  # Slice by Taxon in the larger dataframe
  # Form sliced species-trait matrix
  sharks_traits.var <- FTU_filtered %>%
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
  
  fspaces_quality_sharks.var <- mFD::quality.fspaces(
    sp_dist             = sp_dist_sharks.var,
    maxdim_pcoa         = 10,
    deviation_weighting = c("absolute","squared"),
    fdist_scaling       = FALSE,
    fdendro             = "average")
  
  # Return coordinates of each axis
  sp_faxes_coord_sharks.var <- fspaces_quality_sharks.var$"details_fspaces"$"sp_pc_coord"
  
  # Extract taxa to new df to add to FUSE matrix
  sp_faxes_coord_sharks.taxa <- sp_faxes_coord_sharks.var %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "Taxon")
  
  # Calculate FUSE for each species - set neighbours to 1 to calculate originality instead of uniqueness
  FUSE_res <- fuse(sp_dist        = sp_dist_sharks.var, 
                   sp_faxes_coord = as.matrix(sp_faxes_coord_sharks.var), 
                   nb_NN          = 1,  
                   GE             = GE$Current_status,
                   standGE        = TRUE)
  # Add Taxon column for mean values after loop
  FUSE_res$Taxon <- sp_faxes_coord_sharks.taxa$Taxon
  
  # Output - dataframe of FUSE values for each taxon per iteration
  FUSE_res
  
  #close loop
})

# Merge lists into 1 dataframe
res_FUSE_list <- list(FUSE_res)

res_FUSE_df <- res_FUSE_list %>% 
  bind_rows()

# House-keeping
FUSE_spp <- res_FUSE_df %>% 
  rownames_to_column(var = "Iteration")
save(FUSE_spp, file = "~/FUn_Fsp metrics.RData")

# Set up dataset of mean metric values, filtered taxa, status, & taxonomic information 
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
FUSE_sum$Family <- Taxa_filtered$Family
FUSE_sum$Genus <- Taxa_filtered$Genus
FUSE_sum$Fori_log <- log(FUSE_sum$Fori)
FUSE_sum$Fspe_log <- log(FUSE_sum$Fspe)

status_colors <- c("Extant" = "#868686FF", "Extinct" = "#0073C2FF")

# Top contributors by order and status
## Mean FOri and FSpe across orders
FUSE_sum %>% 
  ungroup() %>% 
  select(Order,Fori,Fspe) %>% 
  group_by(Order) %>% 
  summarise(Ori = mean(Fori),
            Spe = mean(Fspe))

## Mean FOri and FSpe across extinct and extant taxa
FUSE_sum %>% 
  ungroup() %>% 
  select(Status,Fori,Fspe) %>% 
  group_by(Status) %>% 
  summarise(Ori = mean(Fori),
            Spe = mean(Fspe))

## Plots - by Order
FUSE_fori <- FUSE_sum %>%
  select(Taxon, Order, Fori, Status)
FUSE_fori$Taxon <- factor(FUSE_fori$Taxon, levels = FUSE_fori$Taxon[order(FUSE_fori$Fori, decreasing = TRUE)])
FUSE_fori <- FUSE_fori[FUSE_fori$Order != "incert.fam", ]
FUSE_fori$Order <- factor(FUSE_fori$Order, levels = unique(FUSE_fori$Order)[order(unique(FUSE_fori$Fori), decreasing = TRUE)])

plot_FOri_order <- ggplot(FUSE_fori, aes(x = Order, y = log(Fori))) +
  geom_boxplot(data = FUSE_fori, aes(x = reorder(Order,Fori), y = log(Fori)), 
               fill=c("#FAFD7C","#5A9599FF","#FF6F00FF","#84D7E1FF","#008EA0FF",
                               "#8A4198FF","#C71000FF","#ADE2D0FF","#FF95A8FF")) +
                                 coord_flip()+
  labs(x = "", y = "log(FOri)") +
  theme_minimal() +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size= 10, color= "black"),
        axis.title= element_text(size= 12), panel.background= element_rect(fill= "white"))+
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent"))

FUSE_fspe <- FUSE_sum %>%
  select(Taxon, Order, Fspe, Status)
FUSE_fspe$Taxon <- factor(FUSE_fspe$Taxon, levels = FUSE_fspe$Taxon[order(FUSE_fspe$Fspe, decreasing = TRUE)])
FUSE_fspe <- FUSE_fspe[FUSE_fspe$Order != "incert.fam", ]
FUSE_fspe$Order <- factor(FUSE_fspe$Order, levels = unique(FUSE_fspe$Order)[order(unique(FUSE_fspe$Fspe), decreasing = TRUE)])
plot_FSpe_order <- ggplot(FUSE_fspe, aes(x = Order, y = Fspe)) +
  geom_boxplot(data = FUSE_fspe, aes(x = reorder(Order,Fspe), y = Fspe), 
               fill=c("#FF6348FF","#008EA0FF","#8A4198FF","#84D7E1FF","#5A9599FF",
                                 "#FF95A8FF","#FF6F00FF","#FAFD7C","#C71000FF","#ADE2D0FF")) +
                                   coord_flip()+
  labs(x = "", y = "FSpe") +
  theme_minimal() +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size= 10, color= "black"),
        axis.title= element_text(size= 12), panel.background= element_rect(fill= "white"))+
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent"))

# Plot both panels - produces Figure S12
library(patchwork)
combined_plot_order <- plot_FOri_order + plot_FSpe_order +
  plot_layout(widths = c(1, 1))

## Across extinct vs extant taxa
t.test(Fori ~ Status, data = FUSE_fori)   # P = 0.37
t.test(Fspe ~ Status, data = FUSE_fspe)   # P = 0.83

# Produce Figure 3a
plot_FOri_status <- ggplot(FUSE_fori, aes(x = Status, y = log(Fori))) +
  geom_boxplot(data = FUSE_fori, aes(x = Status, y = log(Fori)), 
               fill=c(status_colors)) +
  coord_flip()+
  labs(x = "", y = "log(FOri)") +
  theme_minimal() +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size= 10, color= "black"),
        axis.title= element_text(size= 12), panel.background= element_rect(fill= "white"))+
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent"))

plot_grid(
  plot_FOri_status,
  labels = c("(a)"), 
  label_size = 12, align = "hv", label_fontface = "bold", nrow = 1)

# Produce Figure 3c
plot_FSpe_status <- ggplot(FUSE_fspe, aes(x = Status, y = Fspe)) +
  geom_boxplot(data = FUSE_fspe, aes(x = Status, y = Fspe), 
               fill=c(status_colors)) +
  coord_flip()+
  labs(x = "", y = "FSpe") +
  theme_minimal() +
  theme_bw() +
  theme(legend.position = "right") +
  theme(axis.text = element_text(size= 10, color= "black"),
        axis.title= element_text(size= 12), panel.background= element_rect(fill= "white"))+
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent"))

plot_grid(
  plot_FSpe_status,
  labels = c("(c)"), 
  label_size = 12, align = "hv", label_fontface = "bold", nrow = 1)


# Top 5% contributors
## Rank contributors (top 5% = 27 spp.)
FUSE_fori_rank <- FUSE_sum %>%
  select(Taxon, Fori, Status) %>%
  arrange(desc(Fori)) %>%
  head(27)
FUSE_fori_rank$Taxon <- factor(FUSE_fori_rank$Taxon, levels = FUSE_fori_rank$Taxon[order(FUSE_fori_rank$Fori, decreasing = TRUE)])

FUSE_fspe_rank <- FUSE_sum %>%
  select(Taxon, Fspe, Status) %>%
  arrange(desc(Fspe)) %>%
  head(27)
FUSE_fspe_rank$Taxon <- factor(FUSE_fspe_rank$Taxon, levels = FUSE_fspe_rank$Taxon[order(FUSE_fspe_rank$Fspe, decreasing = TRUE)])

# Create lollipop plot for FOri
plot_FOri_rank_lollipop <- ggplot(FUSE_fori_rank, aes(x = Taxon, y = Fori, color = Status, fill = Status)) +
  geom_segment(aes(xend = Taxon, yend = 0), size = 1.5) +
  geom_point(size = 3) +
  coord_flip() +
  labs(x = "", y = "FOri") +
  scale_fill_manual(values = status_colors, guide = guide_legend(title = "Status")) +
  scale_color_manual(values = status_colors, guide = guide_legend(title = "Status")) +
  theme_minimal() +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text.y = element_text(size = 10, color = "black", face = "italic"),  # Set y-axis labels to italic
        axis.title = element_text(size = 12), panel.background = element_rect(fill = "white")) +
  theme(panel.grid = element_blank(), 
        panel.border = element_blank()  # Remove the panel border
  )

# Rotate - produces Figure 3b
plot_FOri_rank_lollipop <- print(plot_FOri_rank_lollipop, vp=viewport(width = 0.9,
                                                                      height = 0.7,
                                                                      angle=90))

# Create the plot with the reversed x-axis
plot_FSpe_rank_lollipop <- ggplot(FUSE_fspe_rank, aes(x = Taxon, y = Fspe, fill = Status, color = Status)) +
  geom_segment(aes(xend = Taxon, yend = 0), size = 1.5) +
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
  theme(panel.grid = element_blank(), 
        panel.border = element_blank()  # Remove the panel border
  )

# Rotate - produces Figure 3d
plot_FSpe_rank_lollipop <- print(plot_FSpe_rank_lollipop, vp=viewport(width = 0.9,
                                                                      height = 0.7,
                                                                      angle=90))
# Produce Figure 3
plot_grid(
  plot_FOri_status, plot_FOri_rank_lollipop,
  plot_FSpe_status, plot_FSpe_rank_lollipop,
  labels = c("(a)","(b)","(c)","(d)"), 
  label_size = 12, align = "hv", label_fontface = "bold", nrow = 2
)

## Test that extinct FSpe outliers significantly differ from same distribution of extant sharks
# Extract FSpe values for extinct and extant species
fspe_extinct <- FUSE_fspe_rank$Fspe[FUSE_fspe_rank$Status == "Extinct"]
fspe_extant <- FUSE_fspe_rank$Fspe[FUSE_fspe_rank$Status == "Extant"]

# Compute 90th quantiles for each group
quantile_extinct <- quantile(fspe_extinct, probs = 0.9)
quantile_extant <- quantile(fspe_extant, probs = 0.9)

# Define a function to compute the difference in 90th quantiles
compute_diff <- function(x, y) {
  quantile(x, probs = 0.9) - quantile(y, probs = 0.9)
}

# Compute observed difference
observed_diff <- compute_diff(fspe_extinct, fspe_extant)

# Generate permutations of Status variable
n_permutations <- 5000
permuted_diffs <- replicate(n_permutations, {
  status_permuted <- sample(FUSE_fspe_rank$Status)
  fspe_extinct_permuted <- FUSE_fspe_rank$Fspe[status_permuted == "Extinct"]
  fspe_extant_permuted <- FUSE_fspe_rank$Fspe[status_permuted == "Extant"]
  compute_diff(fspe_extinct_permuted, fspe_extant_permuted)
})

# Compute p-value
p_value <- mean(abs(permuted_diffs) >= abs(observed_diff))

# Output results
cat("Observed Difference:", observed_diff, "\n")
cat("P-value:", p_value, "\n")
