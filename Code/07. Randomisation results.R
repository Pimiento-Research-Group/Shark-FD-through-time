###################################################################################################################
# 07. Randomisation results 
## This R code provides the random simulations on shark function diversity through time
## randomisation tests are on resampling and the relationship between taxon richness and FRic
## it produces Figure S5 and S6
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
library(deeptime)

# Load data
load(file="C:/Users/2022207/Dropbox/Jack's PhD/Chapter 2. FD changes over time/Analyses/Current Analyses/R codes/Cleaned data.RData")

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

# Filter taxa that are NOT present in "selected_taxa"
FTU_filtered <- FTU %>%
  anti_join(selected_taxa, by = "Taxon_corrected")

### Resampling model - based on variation
# Separate by epoch
Pal <- data %>% 
  filter(Epoch=="Palaeocene")
Eo <- data %>% 
  filter(Epoch=="Eocene")
Oli <- data %>% 
  filter(Epoch=="Oligocene")
Mio <- data %>% 
  filter(Epoch=="Miocene")
Plio <- data %>% 
  filter(Epoch=="Pliocene")
Ple <- data %>% 
  filter(Epoch=="Pleistocene")
Rec <- data %>% 
  filter(Epoch=="Recent")
NM <- data %>% 
  filter(Epoch=="nonmatch")

# Resample - make big dataframe & then split like before using sampled data within the loop
# Set up model to run on n iterations
registerDoParallel(cores = 10)
getDoParWorkers()
Taxon.res<-NULL
Taxon.res<- lapply(1:1000, function(x){
  # set sample to lowest & randomise (lowest in Pleistocene)
  small.pal<-Pal[sample(1:nrow(Pal), size=309, replace=F),]
  small.eo<-Eo[sample(1:nrow(Eo), size=309, replace=F),]
  small.oli<-Oli[sample(1:nrow(Oli), size=309, replace=F),]
  small.mio<-Mio[sample(1:nrow(Mio), size=309, replace=F),]
  small.plio<-Plio[sample(1:nrow(Plio), size=309, replace=F),]
  small.rec<-Rec[sample(1:nrow(Rec), size=309, replace=F),]
  #rbind to form new dataset
  data.new<-rbind(small.pal,small.eo,small.oli,small.mio,small.plio,Ple,small.rec,NM)
  data.new$Epoch <- ordered(data.new$Epoch, levels=c("Palaeocene","Eocene","Oligocene","Miocene",
                                                     "Pliocene","Pleistocene","Recent","nonmatch"))
  
  abun.res <- data.new %>% 
    group_by(Taxon_corrected) %>% 
    distinct(CH,CW,CE,LC,XO,LO,Epoch,Epoch_earliest,Epoch_latest,Current_status)
  
  # Order taxa alphabetically and order Epochs
  abun.res <- abun.res[order(abun.res$Taxon_corrected),]
  abun.res$Epoch <- ordered(abun.res$Epoch, levels=c("Palaeocene","Eocene","Oligocene","Miocene",
                                                     "Pliocene","Pleistocene", "Recent","nonmatch"))
  
  # Form FTUs column
  abun.res$FTU<-paste(abun.res$Taxon_corrected, abun.res$CH, abun.res$CW, abun.res$CE, abun.res$LC, abun.res$XO, abun.res$LO, sep="+")
  
  # Cast into wide-format data
  wp.abun.res <- abun.res %>% 
    pivot_wider(names_from = Epoch, values_from = FTU)
  
  wp.abun.res <- wp.abun.res[, c("Taxon_corrected","Epoch_earliest","Epoch_latest","Current_status",
                                 "Palaeocene","Eocene","Oligocene","Miocene","Pliocene","Pleistocene","Recent","nonmatch")]
  
  
  # Mark 0 and 1 for absence and presence data in all epochs
  wp.abun.res$Palaeocene[is.na(wp.abun.res$Palaeocene)] <- "0" 
  wp.abun.res$Eocene[is.na(wp.abun.res$Eocene)] <- "0"
  wp.abun.res$Oligocene[is.na(wp.abun.res$Oligocene)] <- "0"
  wp.abun.res$Miocene[is.na(wp.abun.res$Miocene)] <- "0"
  wp.abun.res$Pliocene[is.na(wp.abun.res$Pliocene)] <- "0"
  wp.abun.res$Pleistocene[is.na(wp.abun.res$Pleistocene)] <- "0"
  wp.abun.res$Recent[is.na(wp.abun.res$Recent)] <- "0"
  wp.abun.res$nonmatch[is.na(wp.abun.res$nonmatch)] <- "0"
  
  wp.abun.res <- wp.abun.res %>% 
    mutate(Palaeocene = replace(Palaeocene,Palaeocene!="0", "1")) %>%
    mutate(Eocene = replace(Eocene,Eocene!="0", "1")) %>% 
    mutate(Oligocene = replace(Oligocene,Oligocene!="0", "1")) %>% 
    mutate(Miocene = replace(Miocene,Miocene!="0", "1")) %>%  
    mutate(Pliocene = replace(Pliocene,Pliocene!="0", "1")) %>%
    mutate(Pleistocene = replace(Pleistocene,Pleistocene!="0", "1")) %>% 
    mutate(Recent = replace(Recent,Recent!="0", "1")) %>%
    mutate(nonmatch = replace(nonmatch,nonmatch!="0", "1"))
  
  ## Fill in non-matches
  wp.abun.res <- wp.abun.res %>% 
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
  wp.abun.res <- wp.abun.res %>% 
    rowwise() %>% 
    mutate(Eocene = case_when(Palaeocene == "1" && Current_status == "Extant"~"1", TRUE ~ Eocene)) %>% 
    mutate(Oligocene = case_when(Eocene == "1" && Current_status == "Extant"~"1", TRUE ~ Oligocene)) %>% 
    mutate(Miocene = case_when(Oligocene == "1" && Current_status == "Extant"~"1", TRUE ~ Miocene)) %>% 
    mutate(Pliocene = case_when(Miocene == "1" && Current_status == "Extant"~"1", TRUE ~ Pliocene)) %>% 
    mutate(Pleistocene = case_when(Pliocene == "1" && Current_status == "Extant"~"1", TRUE ~ Pleistocene)) %>% 
    mutate(Recent = case_when(Pleistocene == "1" && Current_status == "Extant"~"1", TRUE ~ Recent))
  
  baskets.range.res <- wp.abun.res %>% 
    select(-c(Epoch_earliest,Epoch_latest,nonmatch))
  
  # Remove all genus-level taxa from Recent where all species are accounted for (i.e., make 0)
  baskets.range.res <- baskets.range.res %>% 
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
  baskets.range.res[baskets.range.res == 0] <- NA
  coalesce_all_columns <- function(df) {
    return(coalesce(!!! as.list(df)))
  }
  
  baskets.range.res <- baskets.range.res %>%
    group_by(Taxon_corrected) %>%
    summarise_all(coalesce_all_columns)
  
  baskets.range.res$Palaeocene[is.na(baskets.range.res$Palaeocene)] <- "0" 
  baskets.range.res$Eocene[is.na(baskets.range.res$Eocene)] <- "0"
  baskets.range.res$Oligocene[is.na(baskets.range.res$Oligocene)] <- "0"
  baskets.range.res$Miocene[is.na(baskets.range.res$Miocene)] <- "0"
  baskets.range.res$Pliocene[is.na(baskets.range.res$Pliocene)] <- "0"
  baskets.range.res$Pleistocene[is.na(baskets.range.res$Pleistocene)] <- "0"
  baskets.range.res$Recent[is.na(baskets.range.res$Recent)] <- "0"
  
  # Fill in any gaps in range automatically
  baskets.range.taxon.res <- baskets.range.res %>% 
    rowwise() %>% 
    mutate(Eocene = case_when(Palaeocene == "1" && Miocene == "1"~ "1", TRUE ~ Eocene)) %>%
    mutate(Oligocene = case_when(Palaeocene == "1" && Miocene == "1"~ "1", TRUE ~ Oligocene)) %>%
    mutate(Oligocene = case_when(Eocene == "1" && Miocene == "1"~ "1", TRUE ~ Oligocene)) %>%
    mutate(Miocene = case_when(Oligocene == "1" && Pliocene == "1"~ "1", TRUE ~ Miocene)) %>%
    mutate(Pliocene = case_when(Miocene == "1" && Pleistocene == "1"~ "1", TRUE ~ Pliocene)) %>%
    mutate(Oligocene = case_when(Eocene == "1" && Pliocene == "1"~ "1", TRUE ~ Oligocene)) %>%
    mutate(Miocene = case_when(Eocene == "1" && Pliocene == "1"~ "1", TRUE ~ Miocene))
  
  # Mark and remove species that only occur in Recent
  selected_taxa.res <- baskets.range.taxon.res %>%
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
  filtered_baskets.range.res <- baskets.range.taxon.res %>%
    anti_join(selected_taxa.res, by = "Taxon_corrected")
  
  # Make taxa row names & tidy dataframe to form final matrix
  baskets.varRes <- filtered_baskets.range.res %>% 
    remove_rownames %>% 
    column_to_rownames(var="Taxon_corrected") %>% 
    select(-Current_status)
  
  # Form functional taxonomic units (FTUs), grouping by species to account for intraspecific variability
  FTU.res <- data.new %>% 
    group_by(Taxon_corrected) %>% 
    distinct(CH,CW,CE,LC,XO,LO)
  
  # Order FTUs alphabetically to help ensure matches between species-trait matrix and occurrence matrix later
  FTU.res <- FTU.res[order(FTU.res$Taxon_corrected),]
  
  # Filter taxa that are NOT present in "selected_taxa"
  FTU_filtered.res <- FTU.res %>%
    anti_join(selected_taxa.res, by = "Taxon_corrected")
  
  # Form sliced species-trait matrix
  sharks_traits.varRes <- FTU_filtered.res %>%
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
resamp_FDmetrics_taxonvar<- res_Taxonvar_resamp_df %>% 
  select(Epoch:fspe)
save(resamp_FDmetrics_taxonvar, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 2. FD changes over time/Analyses/Current Analyses/R codes/Taxon code/Data/Taxon_variation_resampling_metrics.RData")

# Melt data & save
FDmetrics_resamp_long_taxonvar<- melt(resamp_FDmetrics_taxonvar, id.vars= "Epoch")
save(FDmetrics_resamp_long_taxonvar, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 2. FD changes over time/Analyses/Current Analyses/R codes/Taxon code/Data/Taxon_variation_resampling_long.RData")

# Plot resampling in violin plots
TR_resampling_Taxonvar <- FDmetrics_resamp_long_taxonvar %>% 
  filter(variable == "nb_sp")
FEmetrics_resampling_Taxonvar <- FDmetrics_resamp_long_taxonvar %>% 
  filter(variable == "nb_fe")
FRedmetrics_resampling_Taxonvar <- FDmetrics_resamp_long_taxonvar %>% 
  filter(variable == "fred")
FOredmetrics_resampling_Taxonvar <- FDmetrics_resamp_long_taxonvar %>% 
  filter(variable == "fored")
FVulnmetrics_resampling_Taxonvar <- FDmetrics_resamp_long_taxonvar %>% 
  filter(variable == "fvuln")
FRicmetrics_resampling_Taxonvar <- FDmetrics_resamp_long_taxonvar %>% 
  filter(variable == "fric")
FOrimetrics_resampling_Taxonvar <- FDmetrics_resamp_long_taxonvar %>% 
  filter(variable == "fori")
FSpemetrics_resampling_Taxonvar <- FDmetrics_resamp_long_taxonvar %>% 
  filter(variable == "fspe")

# Load empirical results saved from code S4
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 2. FD changes over time/Analyses/Current Analyses/R codes/Taxon code/Data/Taxon_variation_metrics.RData")
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 2. FD changes over time/Analyses/Current Analyses/R codes/Taxon code/Data/Taxon_variation_long_metrics.RData")
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 2. FD changes over time/Analyses/Current Analyses/R codes/Taxon code/Data/Mean_Taxon_metrics.RData")

# Plot empirical and resampling results
epochs_custom <- data.frame(
  name = c("Palaeocene","Eocene","Oligocene","Miocene","Pliocene","Pleistocene","Recent"), 
  max_age = c(1.5,2.5,3.5,4.5,5.5,6.5,7.6),
  min_age = c(0.4,1.5,2.5,3.5,4.5,5.5,6.5), 
  color = c("#FBA75F","#FDB46C","#FDC07A","#FFFF90","#FFFF99","#FFF2AE","#FEF2E0")
)

## Taxon richness - median extracted into Table 1
FDmetrics_TR <- FDmetrics_long_TaxonVar %>% 
  filter(variable == "nb_sp")
epoch_colors <- c("#FBA75F", "#FDB46C", "#FDC07A", "#FFFF90", "#FFFF99", "#FFF2AE", "#FEF2E0")


TR_resampling_variation <- ggplot(data = TR_resampling_Taxonvar, aes(x = Epoch, y = value)) +
  geom_violin(trim=TRUE, fill='#868686', color="#868686", alpha = 0.5) +
  scale_x_discrete(limits = c("Palaeocene", "Eocene", "Oligocene", "Miocene", "Pliocene", "Pleistocene", "Recent")) +
  labs(x = "", y = "TR") +
  scale_fill_manual(values = c("#FBA75F", "#FDB46C", "#FDC07A", "#FFFF90", "#FFFF99", "#FFF2AE", "#FEF2E0")) +
  theme_minimal() +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size = 6.5, color = "black"),
        axis.title = element_text(size = 8),
        panel.background = element_rect(fill = "white")) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(fill = "transparent")) +
  theme(axis.text.x = element_blank()) +
  ylim(70, 280)+
  geom_point(x = 1,y = 105, pch = 21, size = 5, colour = "black", fill = "#FBA75F", stroke = 1)+
  geom_point(x = 2,y = 268, pch = 21, size = 5, colour = "black", fill = "#FDB46C", stroke = 1)+
  geom_point(x = 3,y = 122, pch = 21, size = 5, colour = "black", fill = "#FDC07A", stroke = 1)+
  geom_point(x = 4,y = 235, pch = 21, size = 5, colour = "black", fill = "#FFFF90", stroke = 1)+
  geom_point(x = 5,y = 147, pch = 21, size = 5, colour = "black", fill = "#FFFF99", stroke = 1)+
  geom_point(x = 6,y = 134, pch = 21, size = 5, colour = "black", fill = "#FFF2AE", stroke = 1)+
  geom_point(x = 7,y = 115, pch = 21, size = 5, colour = "black", fill = "#FEF2E0", stroke = 1)
save(TR_resampling_variation, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 2. FD changes over time/Analyses/Current Analyses/R codes/Taxon code/Data/Taxon richness_resampling.RData")

# Plot functional diversity metrics - empirical and resampling; for Figure S7
## Functional entities
FDmetrics_FE <- FDmetrics_long_TaxonVar %>% 
  filter(variable == "nb_fe")
FE_resampling_variation <- ggplot(data=FEmetrics_resampling_Taxonvar, aes(x= Epoch, y= value))+
  geom_violin(trim=TRUE, fill='grey', color="grey", alpha = 0.25)+
  scale_x_discrete(limits=c("Palaeocene","Eocene","Oligocene","Miocene", "Pliocene", "Pleistocene","Recent"))+
  geom_boxplot(data = FDmetrics_FE, aes(x = Epoch, y = value), 
               fill=c("#FBA75F","#FDB46C","#FDC07A","#FFFF90","#FFFF99","#FFF2AE","#FEF2E0"))+
  labs(x = "", y = "# FEs")+
  #ylim(0,55) +
  #scale_y_continuous(limits = c(0, 55), breaks = seq(0, 55, by = 10)) +
  #ggtitle("Functional entities - Resampled taxa with variation")+
  theme_minimal() +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size= 10, color= "black"),
        axis.title= element_text(size= 12), panel.background= element_rect(fill= "white"))+
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent")) +
  theme(axis.text.x = element_blank())

## Functional redundancy
FDmetrics_Fred <- FDmetrics_long_TaxonVar %>% 
  filter(variable == "fred")
FRed_resampling_variation <- ggplot(data=FRedmetrics_resampling_Taxonvar, aes(x= Epoch, y= value))+
  geom_violin(trim=TRUE, fill='grey', color="grey", alpha = 0.25)+
  scale_x_discrete(limits=c("Palaeocene","Eocene","Oligocene","Miocene", "Pliocene", "Pleistocene","Recent"))+
  geom_boxplot(data = FDmetrics_Fred, aes(x = Epoch, y = value), 
               fill=c("#FBA75F","#FDB46C","#FDC07A","#FFFF90","#FFFF99","#FFF2AE","#FEF2E0"))+
  labs(x = "", y = "FRed")+
  #ylim(0,1) +
  #scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  #ggtitle("Functional redundancy - Resampled taxa with variation")+
  theme_minimal() +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size= 10, color= "black"),
        axis.title= element_text(size= 12), panel.background= element_rect(fill= "white"))+
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent")) +
  theme(axis.text.x = element_blank())

## Functional over-redundancy
FDmetrics_Fored <- FDmetrics_long_TaxonVar %>% 
  filter(variable == "fored")
FOred_resampling_variation <- ggplot(data=FOredmetrics_resampling_Taxonvar, aes(x= Epoch, y= value))+
  geom_violin(trim=TRUE, fill='grey', color="grey", alpha = 0.25)+
  scale_x_discrete(limits=c("Palaeocene","Eocene","Oligocene","Miocene", "Pliocene", "Pleistocene","Recent"))+
  geom_boxplot(data = FDmetrics_Fored, aes(x = Epoch, y = value), 
               fill=c("#FBA75F","#FDB46C","#FDC07A","#FFFF90","#FFFF99","#FFF2AE","#FEF2E0"))+
  labs(x = "", y = "FOred")+
  #ylim(0.3,0.6) +
  #scale_y_continuous(limits = c(0.3, 0.6), breaks = seq(0.3, 0.6, by = 0.1)) +
  theme_minimal() +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size= 10, color= "black"),
        axis.title= element_text(size= 12), panel.background= element_rect(fill= "white"))+
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent")) +
  theme(axis.text.x = element_blank())

## Functional vulnerability
FDmetrics_Fvuln <- FDmetrics_long_TaxonVar %>% 
  filter(variable == "fvuln")
FVuln_resampling_variation <- ggplot(data=FVulnmetrics_resampling_Taxonvar, aes(x= Epoch, y= value))+
  geom_violin(trim=TRUE, fill='grey', color="grey", alpha = 0.25)+
  scale_x_discrete(limits=c("Palaeocene","Eocene","Oligocene","Miocene", "Pliocene", "Pleistocene","Recent"))+
  geom_boxplot(data = FDmetrics_Fvuln, aes(x = Epoch, y = value), 
               fill=c("#FBA75F","#FDB46C","#FDC07A","#FFFF90","#FFFF99","#FFF2AE","#FEF2E0"))+
  labs(x = "", y = "FVuln")+
  #ylim(0,1) +
  #scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  theme_minimal() +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size= 10, color= "black"),
        axis.title= element_text(size= 12), panel.background= element_rect(fill= "white"))+
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent")) +
  theme(axis.text.x = element_blank())

## Functional richness
FDmetrics_FRic <- FDmetrics_long_TaxonVar %>% 
  filter(variable == "fric")
FRic_resampling_variation <- ggplot(data=FRicmetrics_resampling_Taxonvar, aes(x= Epoch, y= value))+
  geom_violin(trim=TRUE, fill='grey', color="grey", alpha = 0.25)+
  scale_x_discrete(limits=c("Palaeocene","Eocene","Oligocene","Miocene", "Pliocene", "Pleistocene","Recent"))+
  geom_boxplot(data = FDmetrics_FRic, aes(x = Epoch, y = value), 
               fill=c("#FBA75F","#FDB46C","#FDC07A","#FFFF90","#FFFF99","#FFF2AE","#FEF2E0"))+
  labs(x = "", y = "FRic")+
  ylim(0,1) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  #ggtitle("Functional richness - Resampled taxa with variation")+
  theme_minimal() +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size= 10, color= "black"),
        axis.title= element_text(size= 12), panel.background= element_rect(fill= "white"))+
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent")) +
  theme(axis.text.x = element_blank())

## Functional originality
FDmetrics_Fori <- FDmetrics_long_TaxonVar %>% 
  filter(variable == "fori")
FOri_resampling_variation <- ggplot(data=FOrimetrics_resampling_Taxonvar, aes(x= Epoch, y= value))+
  geom_violin(trim=TRUE, fill='grey', color="grey", alpha = 0.25)+
  scale_x_discrete(limits=c("Palaeocene","Eocene","Oligocene","Miocene", "Pliocene", "Pleistocene","Recent"))+
  geom_boxplot(data = FDmetrics_Fori, aes(x = Epoch, y = value), 
               fill=c("#FBA75F","#FDB46C","#FDC07A","#FFFF90","#FFFF99","#FFF2AE","#FEF2E0"))+
  labs(x = "", y = "FOri")+
  #ylim(0,1) +
  #scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  theme_minimal() +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size= 10, color= "black"),
        axis.title= element_text(size= 12), panel.background= element_rect(fill= "white"))+
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent")) +
  theme(axis.text.x = element_blank())

## Functional specialisation
FDmetrics_Fspe <- FDmetrics_long_TaxonVar %>% 
  filter(variable == "fspe")
FSpe_resampling_variation <- ggplot(data=FSpemetrics_resampling_Taxonvar, aes(x= Epoch, y= value))+
  geom_violin(trim=TRUE, fill='grey', color="grey", alpha = 0.25)+
  scale_x_discrete(limits=c("Palaeocene","Eocene","Oligocene","Miocene", "Pliocene", "Pleistocene","Recent"))+
  geom_boxplot(data = FDmetrics_Fspe, aes(x = Epoch, y = value), 
               fill=c("#FBA75F","#FDB46C","#FDC07A","#FFFF90","#FFFF99","#FFF2AE","#FEF2E0"))+
  labs(x = "Epoch", y = "FSpe")+
  #ylim(0,1) +
  #scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  theme_minimal() +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size= 10, color= "black"),
        axis.title= element_text(size= 12), panel.background= element_rect(fill= "white"))+
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

# Plot all results in a single grid - produces Figure S7
plot_grid(FE_resampling_variation,
          FRed_resampling_variation,
          FOred_resampling_variation,
          FRic_resampling_variation,
          FOri_resampling_variation,
          FSpe_resampling_variation,
          labels= c("(a)","(b)","(c)","(d)","(e)","(f)"), 
          label_size = 12,align = "hv", label_fontface = "bold",  nrow=6)

# Calculate slopes in resampling data and save for statistical analyses
## Separate by epoch
Resampled_Pal<- resamp_FDmetrics_taxonvar[grep("Palaeocene", resamp_FDmetrics_taxonvar$Epoch), (2:9)]
Resampled_Eo<- resamp_FDmetrics_taxonvar[grep("Eocene", resamp_FDmetrics_taxonvar$Epoch), (2:9)]
Resampled_Oli<- resamp_FDmetrics_taxonvar[grep("Oligocene", resamp_FDmetrics_taxonvar$Epoch), (2:9)]
Resampled_Mio<- resamp_FDmetrics_taxonvar[grep("Miocene", resamp_FDmetrics_taxonvar$Epoch), (2:9)]
Resampled_Plio<- resamp_FDmetrics_taxonvar[grep("Pliocene", resamp_FDmetrics_taxonvar$Epoch), (2:9)]
Resampled_Ple<- resamp_FDmetrics_taxonvar[grep("Pleistocene", resamp_FDmetrics_taxonvar$Epoch), (2:9)]
Resampled_Rec<- resamp_FDmetrics_taxonvar[grep("Recent", resamp_FDmetrics_taxonvar$Epoch), (2:9)]

## Calculate & save resampling slopes
slopes_fun<- function(data1, data2){
  output<- matrix(data= NA, nrow= dim(data1)[1], ncol= dim(data1)[2])
  for (k in 1:8) {
    output[,k]<- data2[,k]- data1[,k]
  }
  return(output)
}

Resamp_Pal_Eo_Slopes<- slopes_fun(Resampled_Pal, Resampled_Eo)
save(Resamp_Pal_Eo_Slopes, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 2. FD changes over time/Analyses/Current Analyses/R codes/Taxon code/Data/Pal-Eo resampled slope.RData")
Resamp_Eo_Oli_Slopes<- slopes_fun(Resampled_Eo, Resampled_Oli)
save(Resamp_Eo_Oli_Slopes, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 2. FD changes over time/Analyses/Current Analyses/R codes/Taxon code/Data/Eo-Oli resampled slope.RData")
Resamp_Oli_Mio_Slopes<- slopes_fun(Resampled_Oli, Resampled_Mio)
save(Resamp_Oli_Mio_Slopes, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 2. FD changes over time/Analyses/Current Analyses/R codes/Taxon code/Data/Oli-Mio resampled slope.RData")
Resamp_Mio_Plio_Slopes<- slopes_fun(Resampled_Mio, Resampled_Plio)
save(Resamp_Mio_Plio_Slopes, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 2. FD changes over time/Analyses/Current Analyses/R codes/Taxon code/Data/Mio-Plio resampled slope.RData")
Resamp_Plio_Ple_Slopes<- slopes_fun(Resampled_Plio, Resampled_Ple)
save(Resamp_Plio_Ple_Slopes, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 2. FD changes over time/Analyses/Current Analyses/R codes/Taxon code/Data/Plio-Ple resampled slope.RData")
Resamp_Ple_Rec_Slopes<- slopes_fun(Resampled_Ple, Resampled_Rec)
save(Resamp_Ple_Rec_Slopes, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 2. FD changes over time/Analyses/Current Analyses/R codes/Taxon code/Data/Ple-Rec resampled slope.RData")

# FRic-species richness - produced based on functional space produced in Code S3
N=round(seq(10, length(unique(data$Taxon_corrected)), length.out=50))
reps<-1:100
sims.dat<-expand.grid(N, reps);names(sims.dat)<-c("N", "Rep")
sims.dat.list <- split(sims.dat, seq(nrow(sims.dat)))

species.names <- filtered_baskets.range %>% distinct(Taxon_corrected) %>% as.data.frame()

# Load species-trait matrix from Code S3
load(file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 2. FD changes over time/Analyses/Current Analyses/R codes/Taxon code//Data/Average species-trait matrix.RData")
# Make species-trait matrix
sharks_traits <- Av_sharks_filtered %>%
  remove_rownames %>% 
  column_to_rownames(var="Taxon_corrected")

# Function to calculate mode
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

# Randomisation of species richness-FRic relationship
registerDoParallel(cores = 10)
getDoParWorkers()
res<-NULL
res<-lapply(sims.dat.list, function(x){
  sample.size<-x$N   
  small.unique.taxa <- species.names %>% 
    slice_sample(n=sample.size)
  
  small.taxa <- data %>% 
    filter(Taxon_corrected %in% small.unique.taxa$Taxon_corrected)
  
  Av_data.res <- small.taxa %>%
    group_by(Taxon_corrected) %>% 
    reframe(CH = mean(CH_mm),
            CW = mean(CW_mm),
            CE = getmode(CE),
            LC = getmode(LC),
            XO = getmode(XO),
            LO = getmode(LO))
  
  Av_sharks.res <- Av_data.res %>% 
    rowwise %>% 
    mutate(CH = cut(CH,
                    breaks = c(-Inf,5,20,50,Inf),
                    labels = c("<5mm","5-20mm","20-50mm",">50mm"),
                    right = FALSE)) %>% 
    mutate(CW = cut(CW,
                    breaks = c(-Inf,10,35,Inf),
                    labels = c("<10mm","10-35mm",">35mm"),
                    right = FALSE))
  Av_sharks.res$CH <- ordered(Av_sharks.res$CH, levels=c("<5mm","5-20mm","20-50mm",">50mm"))
  Av_sharks.res$CW <- ordered(Av_sharks.res$CW, levels=c("<10mm","10-35mm",">35mm"))
  
  # Make species-trait matrix
  sharks_traits.res <- Av_sharks.res %>%
    remove_rownames %>% 
    column_to_rownames(var="Taxon_corrected")
  
  ####################################################################
  # Occurrence matrix
  Av_sharks_filtered$original<-1
  Av_sharks_filtered$resampled<-1
  Av_sharks_filtered$resampled[!Av_sharks_filtered$Taxon_corrected %in% Av_sharks.res$Taxon_corrected]<-0  # Ensures species & all associated FTUs as 0
  
  abun<-unique(Av_sharks_filtered[,c("Taxon_corrected","original","resampled")])
  abun.res <- abun %>%
    remove_rownames %>% 
    column_to_rownames(var="Taxon_corrected")
  baskets_sharks_weights.res<-t(abun.res)
  
  # Make matrix and numeric; set NAs to 0
  baskets_sharks_weights.res <- data.matrix(baskets_sharks_weights.res, rownames.force = NA)
  baskets_sharks_weights.res[is.na(baskets_sharks_weights.res)]<-0
  class(baskets_sharks_weights.res) <- "numeric"
  
  # Summarise & extract occurrences
  asb_sp_sharks_summ.res <- mFD::asb.sp.summary(asb_sp_w = baskets_sharks_weights.res)
  
  asb_sp_sharks_occ.res <- asb_sp_sharks_summ.res$asb_sp_occ
  
  ############################################################################################
  # Make trait_cat matrix
  sharks_traits_cat.res <- tibble(trait_name = c("CH", "CW", "CE", "LC", "XO", "LO"),
                                  trait_type = c("O", "O", "N", "N", "N", "N"),
                                  trait_weight = c(0.5,0.5,1,1,0.33,0.67))
  
  ############################################################################################
  # Form the functional entities
  sp_to_fe_sharks.res <- mFD::sp.to.fe(
    sp_tr       = sharks_traits, 
    tr_cat      = sharks_traits_cat.res, 
    fe_nm_type  = "fe_rank", 
    check_input = TRUE)
  
  # Calculate and display metrics
  alpha_fd_fe_sharks.res <- alpha.fd.fe(asb_sp_sharks_occ.res,
                                        sp_to_fe_sharks.res,
                                        ind_nm = c("fred", "fored", "fvuln"),
                                        check_input = TRUE,
                                        details_returned = TRUE)
  
  # Form dataframe
  FEmetrics.res <- as.data.frame(alpha_fd_fe_sharks.res$asb_fdfe)
  
  # Construct trait distance matrix
  sp_dist_sharks.res <- mFD::funct.dist(
    sp_tr         = sharks_traits,
    tr_cat        = sharks_traits_cat.res,
    metric        = "gower",
    scale_euclid  = "scale_center",
    ordinal_var   = "classic",
    weight_type   = "equal",
    stop_if_NA    = TRUE)
  
  # Assess quality
  fspaces_quality_sharks.res <- mFD::quality.fspaces(
    sp_dist             = sp_dist_sharks.res,
    maxdim_pcoa         = 10,
    deviation_weighting = c("absolute","squared"),
    fdist_scaling       = FALSE,
    fdendro             = "average")
  
  # Return coordinates of each axis
  sp_faxes_coord_sharks.res <- fspaces_quality_sharks.res$"details_fspaces"$"sp_pc_coord"
  
  # Form alpha diversity indices
  alpha_fd_indices_sharks.res <- mFD::alpha.fd.multidim(
    sp_faxes_coord   = sp_faxes_coord_sharks.res[ , c("PC1", "PC2", "PC3")],
    asb_sp_w         = baskets_sharks_weights.res,
    ind_vect         = c("fric", "fori", "fspe"),
    scaling          = TRUE,
    check_input      = TRUE,
    details_returned = TRUE)
  
  fd_ind_values_sharks.res <- alpha_fd_indices_sharks.res$"functional_diversity_indices"
  
  ## Form dataframe
  FD.res <- as.data.frame(fd_ind_values_sharks.res)
  
  # Form list to merge datasets
  FDind.res = list(FEmetrics.res,FD.res)
  
  # Output
  FDind.res %>% 
    bind_cols() %>% 
    as.data.frame()
  
  #close loop
})

# Merge res into 1 dataframe
classes.in.list=lapply(res,class)          # find classes in your list
idx=which(classes.in.list=="data.frame") # find indices of data.frames
mydfs=res[idx]                             # subset your list
resamp_df <- data.table::rbindlist(mydfs,idcol='id')     # fast rbind with ID column (install data.table library)

# Focus on even numbered rows, which are all resamples
row_odd <- seq_len(nrow(resamp_df)) %% 2              # Create row indicator
row_odd

resampled_cor <- resamp_df[row_odd == 0, ]            # Subset even rows

save(resampled_cor,file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 2. FD changes over time/Analyses/Current Analyses/R codes/Taxon code//Data/Species richness sequence.RData")

# Correlation test to FRic
cor.test(resampled_cor$nb_sp,resampled_cor$fric)      # Strong positive association (rho = 0.83; P < 2.2e-16)

# Get mean values to plot lines - FRic as key focus with CIs (can we get CI this way?)
resamp_means <- resampled_cor %>%
  group_by(nb_sp) %>%
  summarise(
    mean_fric = mean(fric),
    ci_fric_lower = ifelse(sd(fric) == 0, mean(fric), t.test(fric)$conf.int[1]),
    ci_fric_upper = ifelse(sd(fric) == 0, mean(fric), t.test(fric)$conf.int[2]),
    ci_fric_50_lower = quantile(fric, 0.25),
    ci_fric_50_upper = quantile(fric, 0.75),
    pi_fric_lower = quantile(fric, 0.025),
    pi_fric_upper = quantile(fric, 0.975),
    mean_fe = mean(nb_fe),
    mean_fred = mean(fred),
    mean_fored = mean(fored),
    mean_fvuln = mean(fvuln),
    mean_fori = mean(fori),
    mean_fspe = mean(fspe)
  )

cor.test(resamp_means$nb_sp,resamp_means$mean_fric) # Strong positive association (rho = 0.88; P < 0.0001)

# Plot results - produces Figure S6
annotations <- data.frame(
  nb_sp = c(104, 265, 115, 234, 143, 133, 114),
  fric = c(0.66, 0.86, 0.78, 0.87, 0.69, 0.60, 0.43),
  label = c("Pal", "Eo", "Oli", "Mio", "Plio", "Ple", "Rec")
)

Fig_S6 <- ggplot(data = resampled_cor, aes(x = nb_sp, y = fric)) +
  geom_point(color = "grey80", pch = 1, size = 3.5) +
  geom_ribbon(data = resamp_means,
              aes(x = nb_sp, ymin = pi_fric_lower, ymax = pi_fric_upper),
              fill = "darkorchid3", alpha = .40,
              inherit.aes = FALSE) +
  geom_ribbon(data = resamp_means,
              aes(x = nb_sp, ymin = ci_fric_50_lower, ymax = ci_fric_50_upper),
              fill = "darkorchid3", alpha = .70,
              inherit.aes = FALSE) +
  geom_line(data = resamp_means, aes(x = nb_sp, y = mean_fric), color = "black", linewidth = 1) +
  geom_point(data = annotations,
             aes(x = nb_sp, y = fric),
             pch = 21, size = 5, color = "black", fill = c("#FBA75F", "#FDB46C", "#FDC07A", "#FFFF90", "#FFFF99", "#FFF2AE", "#FEF2E0"),
             stroke = 1) +
  geom_text(data = annotations,
            aes(x = nb_sp, y = fric, label = label),
            vjust = -0.5, hjust = 0.5, size = 3.5, color = "black", position = position_nudge(y = 0.01)) +
  labs(x = "Taxonomic richness", y = "FRic") +
  xlim(0, 537) +
  scale_x_continuous(limits = c(0, 537), breaks = seq(0, 537, by = 20)) +
  theme_minimal() +
  theme_bw() +
  theme(legend.position = "none") +
  theme(
    axis.text = element_text(size = 6.5, color = "black"),
    axis.title = element_text(size = 8),
    panel.background = element_rect(fill = "white")
  ) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(fill = "transparent")
  )
