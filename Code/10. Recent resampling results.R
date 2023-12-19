#########################################################################################################################
# 10. Recent resampling results
## This R code plots the results of our Recent resampling analyses
## it produces Figure S7
#########################################################################################################################

# Import packages
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
library(xfun)
library(deeptime)

# Load data
load(file="C:/Users/2022207/Dropbox/Jack's PhD/Chapter 2. FD changes over time/Analyses/Current Analyses/R codes/Cleaned data.RData")

# Load original results
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 2. FD changes over time/Analyses/Current Analyses/R codes/Taxon code/Data/Taxon_variation_metrics.RData")
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 2. FD changes over time/Analyses/Current Analyses/R codes/Taxon code/Data/Taxon_variation_long_metrics.RData")
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 2. FD changes over time/Analyses/Current Analyses/R codes/Taxon code/Data/Mean_Taxon_metrics.RData")

# Load 'full Recent' data - called "FullRecent"
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 2. FD changes over time/Analyses/Current Analyses/R codes/Taxon code/Data/Full_Recent.RData")
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 2. FD changes over time/Analyses/Current Analyses/R codes/Taxon code/Data/Full_Recent_long.RData")

FullRecent_var <- FullRecent_TaxonVar %>% 
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

# Load 'Recent resampled' data - called "Recent_resampled"
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 2. FD changes over time/Analyses/Current Analyses/R codes/Taxon code/Data/Recent_resampling.RData")
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 2. FD changes over time/Analyses/Current Analyses/R codes/Taxon code/Data/Recent_resampling_long.RData")

Resampled_var <- Recent_resampled %>% 
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

# Load null model - for original & resampled data - called "Null_FDmetrics_taxonvar"
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 2. FD changes over time/Analyses/Current Analyses/R codes/Taxon code/Data/Taxon_variation_null.RData")
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 2. FD changes over time/Analyses/Current Analyses/R codes/Taxon code/Data/Taxon_variation_null_long.RData")


# Null model for 'full Recent'
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
  sharks_traits.nullvar <- FTU %>%
    group_by(Taxon_corrected) %>%
    slice_sample(n=1, replace = FALSE)
  
  # Form species-trait matrix
  sharks_traits.nullvar <- sharks_traits.nullvar %>% 
    ungroup() %>%
    mutate(ID = c(1:590))
  sharks_traits.nullvar$ID <- numbers_to_words(sharks_traits.nullvar$ID)
  
  sharks_traits_null <- sharks_traits.nullvar %>% 
    column_to_rownames(var = "ID") %>% 
    select(-c(Taxon_corrected))
  
  # Form randomised number column in baskets
  baskets.null <- baskets.range.taxon %>%
    ungroup() %>% 
    mutate(vector = sample(1:590, 590, replace= FALSE))
  
  # Order 'vector' and set to column names
  baskets.null <- baskets.null[order(baskets.null$vector),]
  baskets.null$vector <- numbers_to_words(baskets.null$vector)
  
  baskets.null <- baskets.null %>%
    column_to_rownames(var = "vector") %>% 
    select(-c(Taxon_corrected,Current_status))
  
  # Transpose to form final occurrence matrix
  baskets_sharks_null <- t(baskets.null)
  
  
  # Make matrix and numeric
  baskets_sharks_weights.null <- data.matrix(baskets_sharks_null, rownames.force = NA)
  class(baskets_sharks_weights.null) <- "numeric"
  
  # Summarise & extract occurrences
  asb_sp_sharks_summ_null <- mFD::asb.sp.summary(asb_sp_w = baskets_sharks_weights.null)
  
  asb_sp_sharks_occ_null <- asb_sp_sharks_summ_null$asb_sp_occ
  
  ############################################################################################
  # Make trait_cat matrix
  sharks_traits_cat_null <- tibble(trait_name = c("CH", "CW", "CE", "LC", "XO", "LO"),
                                  trait_type = c("O", "O", "N", "N", "N", "N"),
                                  trait_weight = c(0.5,0.5,1,1,0.33,0.67))
  
  ############################################################################################
  # Form the functional entities
  sp_to_fe_sharks_null <- mFD::sp.to.fe(
    sp_tr       = sharks_traits_null, 
    tr_cat      = sharks_traits_cat_null, 
    fe_nm_type  = "fe_rank", 
    check_input = TRUE)
  
  # Calculate and display metrics
  alpha_fd_fe_sharks_null <- alpha.fd.fe(asb_sp_sharks_occ_null,
                                        sp_to_fe_sharks_null,
                                        ind_nm = c("fred", "fored", "fvuln"),
                                        check_input = TRUE,
                                        details_returned = TRUE)
  
  # Form dataframe
  FEmetrics_null <- as.data.frame(alpha_fd_fe_sharks_null$asb_fdfe) %>% 
    tibble::rownames_to_column("Epoch")
  FEmetrics_null$Epoch <- ordered(FEmetrics_null$Epoch, levels=c("Palaeocene","Eocene","Oligocene","Miocene",
                                                               "Pliocene","Pleistocene","Recent"))
  
  # Construct trait distance matrix
  sp_dist_sharks_null <- mFD::funct.dist(
    sp_tr         = sharks_traits_null,
    tr_cat        = sharks_traits_cat_null,
    metric        = "gower",
    scale_euclid  = "scale_center",
    ordinal_var   = "classic",
    weight_type   = "equal",
    stop_if_NA    = TRUE)
  
  # Assess quality
  fspaces_quality_sharks_null <- mFD::quality.fspaces(
    sp_dist             = sp_dist_sharks_null,
    maxdim_pcoa         = 10,
    deviation_weighting = c("absolute","squared"),
    fdist_scaling       = FALSE,
    fdendro             = "average")
  
  # Return coordinates of each axis
  sp_faxes_coord_sharks_null <- fspaces_quality_sharks_null$"details_fspaces"$"sp_pc_coord"
  
  # Form alpha diversity indices
  alpha_fd_indices_sharks_null <- mFD::alpha.fd.multidim(
    sp_faxes_coord   = sp_faxes_coord_sharks_null[ , c("PC1", "PC2", "PC3")],
    asb_sp_w         = baskets_sharks_weights.null,
    ind_vect         = c("fric", "fori", "fspe"),
    scaling          = TRUE,
    check_input      = TRUE,
    details_returned = TRUE)
  
  fd_ind_values_sharks_null <- alpha_fd_indices_sharks_null$"functional_diversity_indices"
  
  ## Form dataframe
  FD_null <- as.data.frame(fd_ind_values_sharks_null) %>% 
    tibble::rownames_to_column("Epoch")
  FD_null$Epoch <- ordered(FD_null$Epoch, levels=c("Palaeocene","Eocene","Oligocene","Miocene",
                                                 "Pliocene","Pleistocene","Recent"))
  
  # Form list to merge datasets
  FDind_null = list(FEmetrics_null,FD_null)
  
  # Output
  FDind_null %>% 
    reduce(inner_join, by = "Epoch")
  
  #close loop
})

stopImplicitCluster()

# Merge lists into 1 dataframe

Recent_null_df <- res.TaxonVar %>% 
  bind_rows()

# Remove sp_richn, which is the same as nb_sp - redundant
Recent_null_df$sp_richn<-NULL

# Format dataframe to be loaded for comparison plots
FullRecent_null<- Recent_null_df %>% 
  select(Epoch:fspe) # remember to save as something
save(FullRecent_null, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 2. FD changes over time/Analyses/Current Analyses/R codes/Taxon code/Data/Full_Recent_null.RData")

# Melt data - remember to save as something before clearing console
FullRecentNull_long<- melt(FullRecent_null, id.vars= "Epoch")
save(FullRecentNull_long, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 2. FD changes over time/Analyses/Current Analyses/R codes/Taxon code/Data/Full_Recent_null_long.RData")


# Plot - FRic x 3
Recent_colour <- data.frame(
  name = "Recent", 
  max_age = 1.6,
  min_age = 0.4, 
  color = "#FEF2E0"
)

# Original plot
FDmetrics_FRic <- FDmetrics_long_TaxonVar %>%
  filter(variable == "fric" & Epoch == "Recent")
FRicmetrics_null_Taxonvar <- FDmetrics_null_long_taxonvar %>% 
  filter(variable == "fric" & Epoch == "Recent")
FRic_null_variation <- ggplot(data=FRicmetrics_null_Taxonvar, aes(x= Epoch, y= value))+
  geom_violin(trim=TRUE, fill='grey', color="grey", alpha = 0.25)+
  scale_x_discrete(limits="Recent")+
  geom_boxplot(data = FDmetrics_FRic, aes(x = Epoch, y = value), 
               fill="#FEF2E0")+
  labs(x = "", y = "FRic")+
  ylim(0,1) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  theme_minimal() +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size= 6.5, color= "black"),
        axis.title= element_text(size= 8), panel.background= element_rect(fill= "white"))+
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent")) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.text.x = element_blank())+
  coord_geo(
    dat = Recent_colour, pos = "bottom", expand = TRUE,
    skip = NULL, abbrv = TRUE, dat_is_discrete = TRUE,
    size = 3
  )


# Full Recent plot
FullRecent_FRic <- FullRecent_long %>% 
  filter(variable == "fric" & Epoch == "Recent")
FRic_null_FullRecent <- FullRecentNull_long  %>% 
  filter(variable == "fric" & Epoch == "Recent")

FRic_Recent_full <- ggplot(data=FRic_null_FullRecent, aes(x= Epoch, y= value))+
  geom_violin(trim=TRUE, fill='grey', color="grey", alpha = 0.25)+
  scale_x_discrete(limits="Recent")+
  geom_boxplot(data = FullRecent_FRic, aes(x = Epoch, y = value), 
               fill="#FEF2E0")+
  labs(x = "", y = "FRic")+
  ylim(0,1) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  theme_minimal() +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size= 6.5, color= "black"),
        axis.title= element_text(size= 8), panel.background= element_rect(fill= "white"))+
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent")) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.text.x = element_blank())+
  coord_geo(
    dat = Recent_colour, pos = "bottom", expand = TRUE,
    skip = NULL, abbrv = TRUE, dat_is_discrete = TRUE,
    size = 3
  )

# Full Recent resampled plot
RecentResamp_FRic <- Recent_resampled_long %>% 
  filter(variable == "fric" & Epoch == "Recent")

FRic_Recent_resampled <- ggplot(data=FRicmetrics_null_Taxonvar, aes(x= Epoch, y= value))+
  geom_violin(trim=TRUE, fill='grey', color="grey", alpha = 0.25)+
  scale_x_discrete(limits="Recent")+
  geom_boxplot(data = RecentResamp_FRic, aes(x = Epoch, y = value), 
               fill="#FEF2E0")+
  labs(x = "", y = "FRic")+
  ylim(0,1) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  theme_minimal() +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size= 6.5, color= "black"),
        axis.title= element_text(size= 8), panel.background= element_rect(fill= "white"))+
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent")) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.text.x = element_blank())+
  coord_geo(
    dat = Recent_colour, pos = "bottom", expand = TRUE,
    skip = NULL, abbrv = TRUE, dat_is_discrete = TRUE,
    size = 3
  )

# Make functional spaces
## Original
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

# Return coordinates of each axis
sp_faxes_coord_sharks <- fspaces_quality_sharks$"details_fspaces"$"sp_pc_coord"

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

## Plot
epochs_custom <- data.frame(
  name = c("Palaeocene","Eocene","Oligocene","Miocene","Pliocene","Pleistocene","Recent"), 
  max_age = c(1.5,2.5,3.5,4.5,5.5,6.5,7.6),
  min_age = c(0.4,1.5,2.5,3.5,4.5,5.5,6.5), 
  color = c("#FBA75F","#FDB46C","#FDC07A","#FFFF90","#FFFF99","#FFF2AE","#FEF2E0")
)


# Form spaces for each epoch, with highest contributors highlighted - Figure 1b-h
## Set axis limits
x_limits_PC1 <- c(-0.4, 0.7)
y_limits_PC12 <- c(-0.3, 1) 
y_limits_PC13 <- c(-0.4, 0.3)

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

## Plot space with highest FOri & FSpe (see Code S10 for FOri/FSpe calculations)
Rec_contributors12 <- plots_Rec$fric$PC1_PC2 +
  labs(x = "PCoA1", y = "PCoA2")+
  coord_fixed(ratio = 1, xlim = x_limits_PC1, ylim = y_limits_PC12)

## Full Recent
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

# Return coordinates of each axis
sp_faxes_coord_sharks <- fspaces_quality_sharks$"details_fspaces"$"sp_pc_coord"

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

## Resampled Recent
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

plots_Rec <- mFD::alpha.multidim.plot(
  output_alpha_fd_multidim = alpha_fd_indices_sharks.var,
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

Rec_resamp <- plots_Rec$fric$PC1_PC2 +
  labs(x = "PCoA1", y = "PCoA2")+
  coord_fixed(ratio = 1, xlim = x_limits_PC1, ylim = y_limits_PC12)
# Z scores for Full Recent & Recent_resampled
# Form empty dataframe
fullRecent.null<- as.data.frame(matrix(data= NA,nrow= 7, ncol= 16, dimnames= list(c("Palaeocene", "Eocene", "Oligocene", "Miocene", "Pliocene", "Pleistocene","Recent"),
                                                                           c("nb_sp","nb_fe", "fred", "fored", "fvuln", "fric", "fori","fspe",
                                                                             "nb_sp_Z","nb_fe_Z", "fred_Z", "fored_Z", "fvuln_Z", "fric_Z", "fori_Z","fspe_Z"))))
# Separate by epoch for each analysis
## Empirical 
FullRecent_TaxonVar_Pal<- FullRecent_TaxonVar[grep("Palaeocene", FullRecent_TaxonVar$Epoch), (2:9)]
FullRecent_TaxonVar_Eo<- FullRecent_TaxonVar[grep("Eocene", FullRecent_TaxonVar$Epoch), (2:9)]
FullRecent_TaxonVar_Oli<- FullRecent_TaxonVar[grep("Oligocene", FullRecent_TaxonVar$Epoch), (2:9)]
FullRecent_TaxonVar_Mio<- FullRecent_TaxonVar[grep("Miocene", FullRecent_TaxonVar$Epoch), (2:9)]
FullRecent_TaxonVar_Plio<- FullRecent_TaxonVar[grep("Pliocene", FullRecent_TaxonVar$Epoch), (2:9)]
FullRecent_TaxonVar_Ple<- FullRecent_TaxonVar[grep("Pleistocene", FullRecent_TaxonVar$Epoch), (2:9)]
FullRecent_TaxonVar_Rec<- FullRecent_TaxonVar[grep("Recent", FullRecent_TaxonVar$Epoch), (2:9)]
## Null model
FullRecent_null_Pal<- FullRecent_null[grep("Palaeocene", FullRecent_null$Epoch), (2:9)]
FullRecent_null_Eo<- FullRecent_null[grep("Eocene", FullRecent_null$Epoch), (2:9)]
FullRecent_null_Oli<- FullRecent_null[grep("Oligocene", FullRecent_null$Epoch), (2:9)]
FullRecent_null_Mio<- FullRecent_null[grep("Miocene", FullRecent_null$Epoch), (2:9)]
FullRecent_null_Plio<- FullRecent_null[grep("Pliocene", FullRecent_null$Epoch), (2:9)]
FullRecent_null_Ple<- FullRecent_null[grep("Pleistocene", FullRecent_null$Epoch), (2:9)]
FullRecent_null_Rec<- FullRecent_null[grep("Recent", FullRecent_null$Epoch), (2:9)]

# Make empirical dataframes - means or medians of empirical & null models
FDempFullRecent <- FullRecent_TaxonVar %>% 
  group_by(Epoch) %>%
  summarise(Sp_med = median(nb_sp),
            FE_med = median(nb_fe),
            Red_med = median(fred),
            Ored_med = median(fored),
            Vul_med = median(fvuln),
            FRic_med = median(fric),
            Fori_med = median(fori),
            Fspe_med = median(fspe))
FDempFullRecent <- FDempFullRecent %>% 
  column_to_rownames(var = "Epoch")

FDnull_FullRecent <- FullRecent_null %>% 
  group_by(Epoch) %>%
  summarise(Sp_med = median(nb_sp),
            FE_med = median(nb_fe),
            Red_med = median(fred),
            Ored_med = median(fored),
            Vul_med = median(fvuln),
            FRic_med = median(fric),
            Fori_med = median(fori),
            Fspe_med = median(fspe))
FDnull_FullRecent <- FDnull_FullRecent %>% 
  column_to_rownames(var = "Epoch")

#Calculate empirical differences; use FDemp as epochs are rownames rather than a column
for(e in 1:8) { 
  for (g in 1:7){
    fullRecent.null[g,e]<- as.numeric(FDempFullRecent[g,e])-as.numeric(FDnull_FullRecent[g,e])
    
  }
}


#write function to calculate slopes for null model values
slopes_fun<- function(data1, data2){
  output<- matrix(data= NA, nrow= dim(data1)[1], ncol= dim(data1)[2])
  for (k in 1:8) {
    output[,k]<- data2[,k]- data1[,k]
  }
  return(output)
}

N_Pal_resamp_Slopes<- slopes_fun(FullRecent_TaxonVar_Pal, FullRecent_null_Pal)
N_Eo_resamp_Slopes<- slopes_fun(FullRecent_TaxonVar_Eo, FullRecent_null_Eo)
N_Oli_resamp_Slopes<- slopes_fun(FullRecent_TaxonVar_Oli, FullRecent_null_Oli)
N_Mio_resamp_Slopes<- slopes_fun(FullRecent_TaxonVar_Mio, FullRecent_null_Mio)
N_Plio_resamp_Slopes<- slopes_fun(FullRecent_TaxonVar_Plio, FullRecent_null_Plio)
N_Ple_resamp_Slopes<- slopes_fun(FullRecent_TaxonVar_Ple, FullRecent_null_Ple)
N_Rec_resamp_Slopes<- slopes_fun(FullRecent_TaxonVar_Rec, FullRecent_null_Rec)

## Calculate Z scores
for (i in 1:8){
  fullRecent.null[1,i+8]<- (fullRecent.null[1,i]- median(N_Pal_resamp_Slopes[,i], na.rm = TRUE))/sd(N_Pal_resamp_Slopes[,i], na.rm= TRUE)}

for (i in 1:8){
  fullRecent.null[2,i+8]<- (fullRecent.null[2,i]- median(N_Eo_resamp_Slopes[,i],  na.rm = TRUE))/sd(N_Eo_resamp_Slopes[,i],  na.rm = TRUE)}

for (i in 1:8){
  fullRecent.null[3,i+8]<- (fullRecent.null[3,i]- median(N_Oli_resamp_Slopes[,i], na.rm = TRUE))/sd(N_Oli_resamp_Slopes[,i],  na.rm = TRUE)}

for (i in 1:8){
  fullRecent.null[4,i+8]<- (fullRecent.null[4,i]- median(N_Mio_resamp_Slopes[,i],  na.rm = TRUE))/sd(N_Mio_resamp_Slopes[,i],  na.rm = TRUE)}

for (i in 1:8){
  fullRecent.null[5,i+8]<- (fullRecent.null[5,i]- median(N_Plio_resamp_Slopes[,i],  na.rm = TRUE))/sd(N_Plio_resamp_Slopes[,i],  na.rm = TRUE)}

for (i in 1:8){
  fullRecent.null[6,i+8]<- (fullRecent.null[6,i]- median(N_Ple_resamp_Slopes[,i],  na.rm = TRUE))/sd(N_Ple_resamp_Slopes[,i],  na.rm = TRUE)}

for (i in 1:8){
  fullRecent.null[7,i+8]<- (fullRecent.null[7,i]- median(N_Rec_resamp_Slopes[,i],  na.rm = TRUE))/sd(N_Rec_resamp_Slopes[,i],  na.rm = TRUE)}


fullRecent.null
## FRic is significantly lower than expected (Z = -4.56)

# Recent_resampled
emp.resamp<- as.data.frame(matrix(data= NA,nrow= 7, ncol= 16, dimnames= list(c("Palaeocene", "Eocene", "Oligocene", "Miocene", "Pliocene", "Pleistocene","Recent"),
                                                                           c("nb_sp","nb_fe", "fred", "fored", "fvuln", "fric", "fori","fspe",
                                                                             "nb_sp_Z","nb_fe_Z", "fred_Z", "fored_Z", "fvuln_Z", "fric_Z", "fori_Z","fspe_Z"))))
# Separate by epoch for each analysis
## Empirical_Recent resampled
FDmetrics_taxonvar_Pal<- Recent_resampled[grep("Palaeocene", Recent_resampled$Epoch), (2:9)]
FDmetrics_taxonvar_Eo<- Recent_resampled[grep("Eocene", Recent_resampled$Epoch), (2:9)]
FDmetrics_taxonvar_Oli<- Recent_resampled[grep("Oligocene", Recent_resampled$Epoch), (2:9)]
FDmetrics_taxonvar_Mio<- Recent_resampled[grep("Miocene", Recent_resampled$Epoch), (2:9)]
FDmetrics_taxonvar_Plio<- Recent_resampled[grep("Pliocene", Recent_resampled$Epoch), (2:9)]
FDmetrics_taxonvar_Ple<- Recent_resampled[grep("Pleistocene", Recent_resampled$Epoch), (2:9)]
FDmetrics_taxonvar_Rec<- Recent_resampled[grep("Recent", Recent_resampled$Epoch), (2:9)]
## Null model
Null_FDmetrics_taxonvar_Pal<- Null_FDmetrics_taxonvar[grep("Palaeocene", Null_FDmetrics_taxonvar$Epoch), (2:9)]
Null_FDmetrics_taxonvar_Eo<- Null_FDmetrics_taxonvar[grep("Eocene", Null_FDmetrics_taxonvar$Epoch), (2:9)]
Null_FDmetrics_taxonvar_Oli<- Null_FDmetrics_taxonvar[grep("Oligocene", Null_FDmetrics_taxonvar$Epoch), (2:9)]
Null_FDmetrics_taxonvar_Mio<- Null_FDmetrics_taxonvar[grep("Miocene", Null_FDmetrics_taxonvar$Epoch), (2:9)]
Null_FDmetrics_taxonvar_Plio<- Null_FDmetrics_taxonvar[grep("Pliocene", Null_FDmetrics_taxonvar$Epoch), (2:9)]
Null_FDmetrics_taxonvar_Ple<- Null_FDmetrics_taxonvar[grep("Pleistocene", Null_FDmetrics_taxonvar$Epoch), (2:9)]
Null_FDmetrics_taxonvar_Rec<- Null_FDmetrics_taxonvar[grep("Recent", Null_FDmetrics_taxonvar$Epoch), (2:9)]


# Make empirical dataframes - means or medians of empirical & null models
FDresamp <- Recent_resampled %>% 
  group_by(Epoch) %>%
  summarise(Sp_med = median(nb_sp),
            FE_med = median(nb_fe),
            Red_med = median(fred),
            Ored_med = median(fored),
            Vul_med = median(fvuln),
            FRic_med = median(fric),
            Fori_med = median(fori),
            Fspe_med = median(fspe))
FDresamp <- FDresamp %>% 
  column_to_rownames(var = "Epoch")

FDnull <- Null_FDmetrics_taxonvar %>% 
  group_by(Epoch) %>%
  summarise(Sp_med = median(nb_sp),
            FE_med = median(nb_fe),
            Red_med = median(fred),
            Ored_med = median(fored),
            Vul_med = median(fvuln),
            FRic_med = median(fric),
            Fori_med = median(fori),
            Fspe_med = median(fspe))
FDnull <- FDnull %>% 
  column_to_rownames(var = "Epoch")

#Calculate empirical differences; use FDemp as epochs are rownames rather than a column
for(e in 1:8) { 
  for (g in 1:7){
    emp.resamp[g,e]<- as.numeric(FDresamp[g,e])-as.numeric(FDnull[g,e])
    
  }
}


#write function to calculate slopes for null model values
slopes_fun<- function(data1, data2){
  output<- matrix(data= NA, nrow= dim(data1)[1], ncol= dim(data1)[2])
  for (k in 1:8) {
    output[,k]<- data2[,k]- data1[,k]
  }
  return(output)
}

N_Pal_null_Slopes<- slopes_fun(FDmetrics_taxonvar_Pal, Null_FDmetrics_taxonvar_Pal)
N_Eo_null_Slopes<- slopes_fun(FDmetrics_taxonvar_Eo, Null_FDmetrics_taxonvar_Eo)
N_Oli_null_Slopes<- slopes_fun(FDmetrics_taxonvar_Oli, Null_FDmetrics_taxonvar_Oli)
N_Mio_null_Slopes<- slopes_fun(FDmetrics_taxonvar_Mio, Null_FDmetrics_taxonvar_Mio)
N_Plio_null_Slopes<- slopes_fun(FDmetrics_taxonvar_Plio, Null_FDmetrics_taxonvar_Plio)
N_Ple_null_Slopes<- slopes_fun(FDmetrics_taxonvar_Ple, Null_FDmetrics_taxonvar_Ple)
N_Rec_null_Slopes<- slopes_fun(FDmetrics_taxonvar_Rec, Null_FDmetrics_taxonvar_Rec)

## Calculate Z scores
for (i in 1:8){
  emp.resamp[1,i+8]<- (emp.resamp[1,i]- median(N_Pal_null_Slopes[,i], na.rm = TRUE))/sd(N_Pal_null_Slopes[,i], na.rm= TRUE)}

for (i in 1:8){
  emp.resamp[2,i+8]<- (emp.resamp[2,i]- median(N_Eo_null_Slopes[,i],  na.rm = TRUE))/sd(N_Eo_null_Slopes[,i],  na.rm = TRUE)}

for (i in 1:8){
  emp.resamp[3,i+8]<- (emp.resamp[3,i]- median(N_Oli_null_Slopes[,i], na.rm = TRUE))/sd(N_Oli_null_Slopes[,i],  na.rm = TRUE)}

for (i in 1:8){
  emp.resamp[4,i+8]<- (emp.resamp[4,i]- median(N_Mio_null_Slopes[,i],  na.rm = TRUE))/sd(N_Mio_null_Slopes[,i],  na.rm = TRUE)}

for (i in 1:8){
  emp.resamp[5,i+8]<- (emp.resamp[5,i]- median(N_Plio_null_Slopes[,i],  na.rm = TRUE))/sd(N_Plio_null_Slopes[,i],  na.rm = TRUE)}

for (i in 1:8){
  emp.resamp[6,i+8]<- (emp.resamp[6,i]- median(N_Ple_null_Slopes[,i],  na.rm = TRUE))/sd(N_Ple_null_Slopes[,i],  na.rm = TRUE)}

for (i in 1:8){
  emp.resamp[7,i+8]<- (emp.resamp[7,i]- median(N_Rec_null_Slopes[,i],  na.rm = TRUE))/sd(N_Rec_null_Slopes[,i],  na.rm = TRUE)}


emp.resamp 
## FRic does not significantly deviate from expectations in Recent (Z = -1.69)

# Plot spaces & box/violin plots together
plot_grid(Rec_contributors12, Rec, Rec_resamp,
          FRic_null_variation,
          FRic_Recent_full, FRic_Recent_resampled,
                   labels= c("(a)","(c)","(e)","(b)","(d)","(f)"), 
                   label_size = 10,align = "hv", label_fontface = "bold",  nrow=2)
