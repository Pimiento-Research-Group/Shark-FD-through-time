########################################################################################################################################
# 14. Stage-based FD null model
## This R code provides the null mode for functional diversity analyses on stage-level analyses
## it produces data necessary to produce supplementary figure S6
#######################################################################################################################################

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
library(xfun)

# Load data
load(file="~/Cleaned data.RData")

# Form stage column with non-matches
data <- data %>% 
  mutate(Stage = Stage_earliest) %>%
  rowwise() %>%
  mutate(Stage = replace(Stage, !Stage_earliest %in% Stage_latest, "nonmatch"))

## Form dataframe of traits and occurrences; to be split into two input matrices in the loop
# Occurrence matrix
# Form FTU unique combinations, with stage included
abun <- data %>% 
  group_by(Taxon_corrected) %>% 
  distinct(CH,CW,CE,LC,XO,LO,Stage,Stage_earliest,Stage_latest,Current_status)

# Order taxa alphabetically and order Epochs
abun <- abun[order(abun$Taxon_corrected),]
abun$Stage <- ordered(abun$Stage, 
                      levels=c("Danian","Selandian","Thanetian",
                               "Ypresian","Lutetian","Bartonian","Priabonian",
                               "Rupelian","Chattian","Aquitanian","Burdigalian",
                               "Langhian","Serravallian","Tortonian","Messinian",
                               "Zanclean","Piacenzian","Gelasian","Calabrian",
                               "Chibanian","Late/Upper","Recent","nonmatch"))

# Form FTUs column
abun$FTU<-paste(abun$Taxon_corrected, abun$CH, abun$CW, abun$CE, abun$LC, abun$XO, abun$LO, sep="+")

# Cast into wide-format data
wp.abun <- abun %>% 
  pivot_wider(names_from = Stage, values_from = FTU)

wp.abun <- wp.abun[, c("Taxon_corrected","Stage_earliest","Stage_latest","Current_status",
                       "Danian","Selandian","Thanetian",
                       "Ypresian","Lutetian","Bartonian","Priabonian",
                       "Rupelian","Chattian","Aquitanian","Burdigalian",
                       "Langhian","Serravallian","Tortonian","Messinian",
                       "Zanclean","Piacenzian","Gelasian","Calabrian",
                       "Chibanian","Late/Upper","Recent","nonmatch")]


# Mark 0 and 1 for absence and presence data in all epochs
wp.abun$Danian[is.na(wp.abun$Danian)] <- "0" 
wp.abun$Selandian[is.na(wp.abun$Selandian)] <- "0"
wp.abun$Thanetian[is.na(wp.abun$Thanetian)] <- "0"
wp.abun$Ypresian[is.na(wp.abun$Ypresian)] <- "0"
wp.abun$Lutetian[is.na(wp.abun$Lutetian)] <- "0"
wp.abun$Bartonian[is.na(wp.abun$Bartonian)] <- "0"
wp.abun$Priabonian[is.na(wp.abun$Priabonian)] <- "0" 
wp.abun$Rupelian[is.na(wp.abun$Rupelian)] <- "0"
wp.abun$Chattian[is.na(wp.abun$Chattian)] <- "0"
wp.abun$Aquitanian[is.na(wp.abun$Aquitanian)] <- "0"
wp.abun$Burdigalian[is.na(wp.abun$Burdigalian)] <- "0"
wp.abun$Langhian[is.na(wp.abun$Langhian)] <- "0"
wp.abun$Serravallian[is.na(wp.abun$Serravallian)] <- "0"
wp.abun$Tortonian[is.na(wp.abun$Tortonian)] <- "0"
wp.abun$Messinian[is.na(wp.abun$Messinian)] <- "0"
wp.abun$Zanclean[is.na(wp.abun$Zanclean)] <- "0" 
wp.abun$Piacenzian[is.na(wp.abun$Piacenzian)] <- "0"
wp.abun$Gelasian[is.na(wp.abun$Gelasian)] <- "0"
wp.abun$Calabrian[is.na(wp.abun$Calabrian)] <- "0"
wp.abun$Chibanian[is.na(wp.abun$Chibanian)] <- "0"
wp.abun$`Late/Upper`[is.na(wp.abun$`Late/Upper`)] <- "0"
wp.abun$Recent[is.na(wp.abun$Recent)] <- "0"
wp.abun$nonmatch[is.na(wp.abun$nonmatch)] <- "0"

wp.abun <- wp.abun %>% 
  mutate(Danian = replace(Danian,Danian!="0", "1")) %>%
  mutate(Selandian = replace(Selandian,Selandian!="0", "1")) %>% 
  mutate(Thanetian = replace(Thanetian,Thanetian!="0", "1")) %>% 
  mutate(Ypresian = replace(Ypresian,Ypresian!="0", "1")) %>%  
  mutate(Lutetian = replace(Lutetian,Lutetian!="0", "1")) %>%
  mutate(Bartonian = replace(Bartonian,Bartonian!="0", "1")) %>% 
  mutate(Priabonian = replace(Priabonian,Priabonian!="0", "1")) %>%
  mutate(Rupelian = replace(Rupelian,Rupelian!="0", "1")) %>% 
  mutate(Chattian = replace(Chattian,Chattian!="0", "1")) %>% 
  mutate(Aquitanian = replace(Aquitanian,Aquitanian!="0", "1")) %>%  
  mutate(Burdigalian = replace(Burdigalian,Burdigalian!="0", "1")) %>%
  mutate(Langhian = replace(Langhian,Langhian!="0", "1")) %>%
  mutate(Serravallian = replace(Serravallian,Serravallian!="0", "1")) %>%
  mutate(Tortonian = replace(Tortonian,Tortonian!="0", "1")) %>% 
  mutate(Messinian = replace(Messinian,Messinian!="0", "1")) %>% 
  mutate(Zanclean = replace(Zanclean,Zanclean!="0", "1")) %>%  
  mutate(Piacenzian = replace(Piacenzian,Piacenzian!="0", "1")) %>%
  mutate(Gelasian = replace(Gelasian,Gelasian!="0", "1")) %>%
  mutate(Calabrian = replace(Calabrian,Calabrian!="0", "1")) %>%  
  mutate(Chibanian = replace(Chibanian,Chibanian!="0", "1")) %>%
  mutate(`Late/Upper` = replace(`Late/Upper`,`Late/Upper`!="0", "1")) %>%
  mutate(Recent = replace(Recent,Recent!="0", "1")) %>%
  mutate(nonmatch = replace(nonmatch,nonmatch!="0", "1"))

## Fill in non-matches
wp.abun <- wp.abun %>% 
  rowwise() %>% 
  mutate(
    Danian = if_else(Stage_earliest == "Danian" & nonmatch == "1", "1", Danian),
    Danian = if_else(Stage_latest == "Danian" & nonmatch == "1", "1", Danian),
    Selandian = if_else(Stage_earliest == "Selandian" & nonmatch == "1", "1", Selandian),
    Selandian = if_else(Stage_latest == "Selandian" & nonmatch == "1", "1", Selandian),
    Thanetian = if_else(Stage_earliest == "Thanetian" & nonmatch == "1", "1", Thanetian),
    Thanetian = if_else(Stage_latest == "Thanetian" & nonmatch == "1", "1", Thanetian),
    Ypresian = if_else(Stage_earliest == "Ypresian" & nonmatch == "1", "1", Ypresian),
    Ypresian = if_else(Stage_latest == "Ypresian" & nonmatch == "1", "1", Ypresian),
    Lutetian = if_else(Stage_earliest == "Lutetian" & nonmatch == "1", "1", Lutetian),
    Lutetian = if_else(Stage_latest == "Lutetian" & nonmatch == "1", "1", Lutetian),
    Bartonian = if_else(Stage_earliest == "Bartonian" & nonmatch == "1", "1", Bartonian),
    Bartonian = if_else(Stage_latest == "Bartonian" & nonmatch == "1", "1", Bartonian),
    Priabonian = if_else(Stage_earliest == "Priabonian" & nonmatch == "1", "1", Priabonian),
    Priabonian = if_else(Stage_latest == "Priabonian" & nonmatch == "1", "1", Priabonian),
    Rupelian = if_else(Stage_earliest == "Rupelian" & nonmatch == "1", "1", Rupelian),
    Rupelian = if_else(Stage_latest == "Rupelian" & nonmatch == "1", "1", Rupelian),
    Chattian = if_else(Stage_earliest == "Chattian" & nonmatch == "1", "1", Chattian),
    Chattian = if_else(Stage_latest == "Chattian" & nonmatch == "1", "1", Chattian),
    Aquitanian = if_else(Stage_earliest == "Aquitanian" & nonmatch == "1", "1", Aquitanian),
    Aquitanian = if_else(Stage_latest == "Aquitanian" & nonmatch == "1", "1", Aquitanian),
    Burdigalian = if_else(Stage_earliest == "Burdigalian" & nonmatch == "1", "1", Burdigalian),
    Burdigalian = if_else(Stage_latest == "Burdigalian" & nonmatch == "1", "1", Burdigalian),
    Langhian = if_else(Stage_earliest == "Langhian" & nonmatch == "1", "1", Langhian),
    Langhian = if_else(Stage_latest == "Langhian" & nonmatch == "1", "1", Langhian),
    Serravallian = if_else(Stage_earliest == "Serravallian" & nonmatch == "1", "1", Serravallian),
    Serravallian = if_else(Stage_latest == "Serravallian" & nonmatch == "1", "1", Serravallian),
    Tortonian = if_else(Stage_earliest == "Tortonian" & nonmatch == "1", "1", Tortonian),
    Tortonian = if_else(Stage_latest == "Tortonian" & nonmatch == "1", "1", Tortonian),
    Messinian = if_else(Stage_earliest == "Messinian" & nonmatch == "1", "1", Messinian),
    Messinian = if_else(Stage_latest == "Messinian" & nonmatch == "1", "1", Messinian),
    Zanclean = if_else(Stage_earliest == "Zanclean" & nonmatch == "1", "1", Zanclean),
    Zanclean = if_else(Stage_latest == "Zanclean" & nonmatch == "1", "1", Zanclean),
    Piacenzian = if_else(Stage_earliest == "Piacenzian" & nonmatch == "1", "1", Piacenzian),
    Piacenzian = if_else(Stage_latest == "Piacenzian" & nonmatch == "1", "1", Piacenzian),
    Gelasian = if_else(Stage_earliest == "Gelasian" & nonmatch == "1", "1", Gelasian),
    Gelasian = if_else(Stage_latest == "Gelasian" & nonmatch == "1", "1", Gelasian),
    Calabrian = if_else(Stage_earliest == "Calabrian" & nonmatch == "1", "1", Calabrian),
    Calabrian = if_else(Stage_latest == "Calabrian" & nonmatch == "1", "1", Calabrian),
    Chibanian = if_else(Stage_earliest == "Chibanian" & nonmatch == "1", "1", Chibanian),
    Chibanian = if_else(Stage_latest == "Chibanian" & nonmatch == "1", "1", Chibanian),
    `Late/Upper` = if_else(Stage_earliest == "Late/Upper" & nonmatch == "1", "1", `Late/Upper`),
    `Late/Upper` = if_else(Stage_latest == "Late/Upper" & nonmatch == "1", "1", `Late/Upper`),
    Recent = if_else(Stage_earliest == "Recent" & nonmatch == "1", "1", Recent),
    Recent = if_else(Stage_latest == "Recent" & nonmatch == "1", "1", Recent))

# Range through via "extant" status
# Automatically range-through extant taxa
wp.abun <- wp.abun %>% 
  rowwise() %>% 
  mutate(
    Selandian = case_when(Danian == "1" && Current_status == "Extant" ~ "1", TRUE ~ Selandian),
    Thanetian = case_when(Selandian == "1" && Current_status == "Extant" ~ "1", TRUE ~ Thanetian),
    Ypresian = case_when(Thanetian == "1" && Current_status == "Extant" ~ "1", TRUE ~ Ypresian),
    Lutetian = case_when(Ypresian == "1" && Current_status == "Extant" ~ "1", TRUE ~ Lutetian),
    Bartonian = case_when(Lutetian == "1" && Current_status == "Extant" ~ "1", TRUE ~ Bartonian),
    Priabonian = case_when(Bartonian == "1" && Current_status == "Extant" ~ "1", TRUE ~ Priabonian),
    Rupelian = case_when(Priabonian == "1" && Current_status == "Extant" ~ "1", TRUE ~ Rupelian),
    Chattian = case_when(Rupelian == "1" && Current_status == "Extant" ~ "1", TRUE ~ Chattian),
    Aquitanian = case_when(Chattian == "1" && Current_status == "Extant" ~ "1", TRUE ~ Aquitanian),
    Burdigalian = case_when(Aquitanian == "1" && Current_status == "Extant" ~ "1", TRUE ~ Burdigalian),
    Langhian = case_when(Burdigalian == "1" && Current_status == "Extant" ~ "1", TRUE ~ Langhian),
    Serravallian = case_when(Langhian == "1" && Current_status == "Extant" ~ "1", TRUE ~ Serravallian),
    Tortonian = case_when(Serravallian == "1" && Current_status == "Extant" ~ "1", TRUE ~ Tortonian),
    Messinian = case_when(Tortonian == "1" && Current_status == "Extant" ~ "1", TRUE ~ Messinian),
    Zanclean = case_when(Messinian == "1" && Current_status == "Extant" ~ "1", TRUE ~ Zanclean),
    Piacenzian = case_when(Zanclean == "1" && Current_status == "Extant" ~ "1", TRUE ~ Piacenzian),
    Gelasian = case_when(Piacenzian == "1" && Current_status == "Extant" ~ "1", TRUE ~ Gelasian),
    Calabrian = case_when(Gelasian == "1" && Current_status == "Extant" ~ "1", TRUE ~ Calabrian),
    Chibanian = case_when(Calabrian == "1" && Current_status == "Extant" ~ "1", TRUE ~ Chibanian),
    `Late/Upper` = case_when(Chibanian == "1" && Current_status == "Extant" ~ "1", TRUE ~ `Late/Upper`),
    Recent = case_when(`Late/Upper` == "1" && Current_status == "Extant" ~ "1", TRUE ~ Recent)
  )

baskets.range <- wp.abun %>% 
  select(-c(Stage_earliest,Stage_latest,nonmatch))

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

baskets.range$Danian[is.na(baskets.range$Danian)] <- "0" 
baskets.range$Selandian[is.na(baskets.range$Selandian)] <- "0"
baskets.range$Thanetian[is.na(baskets.range$Thanetian)] <- "0"
baskets.range$Ypresian[is.na(baskets.range$Ypresian)] <- "0"
baskets.range$Lutetian[is.na(baskets.range$Lutetian)] <- "0"
baskets.range$Bartonian[is.na(baskets.range$Bartonian)] <- "0"
baskets.range$Priabonian[is.na(baskets.range$Priabonian)] <- "0" 
baskets.range$Rupelian[is.na(baskets.range$Rupelian)] <- "0"
baskets.range$Chattian[is.na(baskets.range$Chattian)] <- "0"
baskets.range$Aquitanian[is.na(baskets.range$Aquitanian)] <- "0"
baskets.range$Burdigalian[is.na(baskets.range$Burdigalian)] <- "0"
baskets.range$Langhian[is.na(baskets.range$Langhian)] <- "0"
baskets.range$Serravallian[is.na(baskets.range$Serravallian)] <- "0"
baskets.range$Tortonian[is.na(baskets.range$Tortonian)] <- "0"
baskets.range$Messinian[is.na(baskets.range$Messinian)] <- "0"
baskets.range$Zanclean[is.na(baskets.range$Zanclean)] <- "0" 
baskets.range$Piacenzian[is.na(baskets.range$Piacenzian)] <- "0"
baskets.range$Gelasian[is.na(baskets.range$Gelasian)] <- "0"
baskets.range$Calabrian[is.na(baskets.range$Calabrian)] <- "0"
baskets.range$Chibanian[is.na(baskets.range$Chibanian)] <- "0"
baskets.range$`Late/Upper`[is.na(baskets.range$`Late/Upper`)] <- "0"
baskets.range$Recent[is.na(baskets.range$Recent)] <- "0"

# Fill in any gaps in range automatically
# Function to count non-consecutive 1s
count_non_consecutive_ones <- function(row) {
  consecutive_ones <- rle(row)$lengths
  num_non_consecutive_ones <- sum(consecutive_ones > 1)
  return(num_non_consecutive_ones)
}

# Apply the function row-wise
num_non_consecutive <- apply(baskets.range[, -1], 1, count_non_consecutive_ones)

# Function to fill in non-consecutive 1s
fill_non_consecutive_ones <- function(row) {
  ones_index <- which(row == "1")
  if (length(ones_index) <= 1) {
    return(row)
  }
  for (i in 1:(length(ones_index)-1)) {
    if (!is.na(ones_index[i]) && !is.na(ones_index[i+1])) {
      if ((ones_index[i+1] - ones_index[i]) > 1) {
        row[(ones_index[i]+1):(ones_index[i+1]-1)] <- "1"
      }
    }
  }
  return(row)
}

# Apply the function row-wise
filled_baskets <- t(apply(baskets.range[, -1], 1, function(row) {
  if (sum(row == "1", na.rm = TRUE) > 0 && !all(is.na(row))) {
    fill_non_consecutive_ones(row)
  } else {
    row
  }
}))
filled_baskets <- cbind(baskets.range[, 1, drop = FALSE], filled_baskets)

# Mark and remove species that only occur in Recent
selected_taxa <- filled_baskets %>%
  filter(
    Recent == "1" &
      Danian == "0" &
      Selandian == "0" &
      Thanetian == "0" &
      Ypresian == "0" &
      Lutetian == "0" &
      Bartonian == "0" &
      Priabonian == "0" &
      Rupelian == "0" &
      Chattian == "0" &
      Aquitanian == "0" &
      Burdigalian == "0" &
      Langhian == "0" &
      Serravallian == "0" &
      Tortonian == "0" &
      Messinian == "0" &
      Zanclean == "0" &
      Piacenzian == "0" &
      Gelasian == "0" &
      Calabrian == "0" &
      Chibanian == "0" &
      `Late/Upper` == "0"
  ) %>%
  select(Taxon_corrected)
# Filter taxa that are NOT present in "selected_taxa"
filtered_baskets.range <- filled_baskets %>%
  anti_join(selected_taxa, by = "Taxon_corrected")

# Make taxa row names & tidy dataframe to form final matrix
baskets <- filtered_baskets.range %>% 
  remove_rownames %>% 
  column_to_rownames(var="Taxon_corrected") %>% 
  select(-Current_status)

# Form functional taxonomic units (FTUs), grouping by species to account for intraspecific variability
FTU <- data %>% 
  group_by(Taxon_corrected) %>% 
  distinct(CH,CW,CE,LC,XO,LO)

# Order FTUs alphabetically to help ensure matches between species-trait matrix and occurrence matrix later
FTU <- FTU[order(FTU$Taxon_corrected),]

# Filter taxa that are NOT present in "selected_taxa"
FTU_filtered <- FTU %>%
  anti_join(selected_taxa, by = "Taxon_corrected")

## Run iterations & calculate mean by taxon
# Set up model to run on n iterations
registerDoParallel(cores = 10)
getDoParWorkers()
res.StageVar.null<-NULL
res.StageVar.null<- lapply(1:1000, function(x){
  # Form sliced species-trait matrix
  sharks_traits.nullvar <- FTU_filtered %>%
    group_by(Taxon_corrected) %>%
    slice_sample(n=1, replace = FALSE) %>% 
    remove_rownames
  
  # Form species-trait matrix
  sharks_traits_null <- sharks_traits.nullvar %>% 
    ungroup() %>%
    mutate(ID = c(1:537))
  sharks_traits_null$ID <- numbers_to_words(sharks_traits_null$ID)
  
  sharks_traits.null <- sharks_traits_null %>% 
    column_to_rownames(var = "ID") %>% 
    select(-c(Taxon_corrected))
  
  # Form randomised number column in baskets
  baskets.null <- filtered_baskets.range %>%
    ungroup() %>% 
    mutate(vector = sample(1:537, 537, replace= FALSE))
  
  # Order 'vector' and set to column names
  baskets.null <- baskets.null[order(baskets.null$vector),]
  baskets.null$vector <- numbers_to_words(baskets.null$vector)
  
  baskets.null <- baskets.null %>%
    remove_rownames() %>% 
    column_to_rownames(var = "vector") %>% 
    select(-c(Taxon_corrected,Current_status))
  
  # Transpose to form final occurrence matrix
  baskets_sharks_null <- t(baskets.null)
  
  # Make matrix and numeric
  baskets_sharks_weights.null <- data.matrix(baskets_sharks_null, rownames.force = NA)
  class(baskets_sharks_weights.null) <- "numeric"
  
  # Summarise & extract occurrences
  asb_sp_sharks_summ.null <- mFD::asb.sp.summary(asb_sp_w = baskets_sharks_weights.null)
  
  asb_sp_sharks_occ.null <- asb_sp_sharks_summ.null$asb_sp_occ
  
  ############################################################################################
  # Make trait_cat matrix
  sharks_traits_cat.null <- tibble(trait_name = c("CH", "CW", "CE", "LC", "XO", "LO"),
                                  trait_type = c("O", "O", "N", "N", "N", "N"),
                                  trait_weight = c(0.5,0.5,1,1,0.33,0.67))
  
  ############################################################################################
  # Form the functional entities
  sp_to_fe_sharks.null <- mFD::sp.to.fe(
    sp_tr       = sharks_traits.null, 
    tr_cat      = sharks_traits_cat.null, 
    fe_nm_type  = "fe_rank", 
    check_input = TRUE)
  
  # Calculate and display metrics
  alpha_fd_fe_sharks.null <- alpha.fd.fe(asb_sp_sharks_occ.null,
                                        sp_to_fe_sharks.null,
                                        ind_nm = c("fred", "fored", "fvuln"),
                                        check_input = TRUE,
                                        details_returned = TRUE)
  
  # Form dataframe
  FEmetrics.null <- as.data.frame(alpha_fd_fe_sharks.null$asb_fdfe) %>% 
    tibble::rownames_to_column("Stage")
  FEmetrics.null$Stage <- ordered(FEmetrics.null$Stage, 
                                 levels=c("Danian","Selandian","Thanetian",
                                          "Ypresian","Lutetian","Bartonian","Priabonian",
                                          "Rupelian","Chattian","Aquitanian","Burdigalian",
                                          "Langhian","Serravallian","Tortonian","Messinian",
                                          "Zanclean","Piacenzian","Gelasian","Calabrian",
                                          "Chibanian","Late/Upper","Recent"))
  
  # Construct trait distance matrix
  sp_dist_sharks.null <- mFD::funct.dist(
    sp_tr         = sharks_traits.null,
    tr_cat        = sharks_traits_cat.null,
    metric        = "gower",
    scale_euclid  = "scale_center",
    ordinal_var   = "classic",
    weight_type   = "equal",
    stop_if_NA    = TRUE)
  
  # Assess quality
  fspaces_quality_sharks.null <- mFD::quality.fspaces(
    sp_dist             = sp_dist_sharks.null,
    maxdim_pcoa         = 10,
    deviation_weighting = c("absolute","squared"),
    fdist_scaling       = FALSE,
    fdendro             = "average")
  
  # Return coordinates of each axis
  sp_faxes_coord_sharks.null <- fspaces_quality_sharks.null$"details_fspaces"$"sp_pc_coord"
  
  # Form alpha diversity indices
  alpha_fd_indices_sharks.null <- mFD::alpha.fd.multidim(
    sp_faxes_coord   = sp_faxes_coord_sharks.null[ , c("PC1", "PC2", "PC3")],
    asb_sp_w         = baskets_sharks_weights.null,
    ind_vect         = c("fric", "fori", "fspe"),
    scaling          = TRUE,
    check_input      = TRUE,
    details_returned = TRUE)
  
  fd_ind_values_sharks.null <- alpha_fd_indices_sharks.null$"functional_diversity_indices"
  
  ## Form dataframe
  FD.null <- as.data.frame(fd_ind_values_sharks.null) %>% 
    tibble::rownames_to_column("Stage")
  FD.null$Stage <- ordered(FD.null$Stage, levels=c("Danian","Selandian","Thanetian",
                                                 "Ypresian","Lutetian","Bartonian","Priabonian",
                                                 "Rupelian","Chattian","Aquitanian","Burdigalian",
                                                 "Langhian","Serravallian","Tortonian","Messinian",
                                                 "Zanclean","Piacenzian","Gelasian","Calabrian",
                                                 "Chibanian","Late/Upper","Recent"))
  
  # Form list to merge datasets
  FDind.null = list(FEmetrics.null,FD.null)
  
  # Output
  FDind.null %>% 
    reduce(inner_join, by = "Stage")
  
  #close loop
})

stopImplicitCluster()

# Merge lists into 1 dataframe

res_df_StageNull <- res.StageVar.null %>% 
  bind_rows()

# Remove sp_richn, which is the same as nb_sp - redundant
res_df_StageNull$sp_richn<-NULL

# Format dataframe to be loaded for comparison plots
Res_FDmetrics_StageNull<- res_df_StageNull %>% 
  select(Stage:fspe)
save(Res_FDmetrics_StageNull, file = "~/Stage_null_variation_metrics.RData")

# Melt data
FDmetrics_long_StageNull<- melt(Res_FDmetrics_StageNull, id.vars= "Stage")
save(FDmetrics_long_StageNull, file = "~/Stage_null_variation_long_metrics.RData")

# Form dataframe of mean, median and standard deviation of all FD metrics
Stage_null <- Res_FDmetrics_StageNull %>% 
  group_by(Stage) %>%
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

Stage_null$Stage <- ordered(Stage_null$Stage, levels=c("Danian","Selandian","Thanetian",
                                                     "Ypresian","Lutetian","Bartonian","Priabonian",
                                                     "Rupelian","Chattian","Aquitanian","Burdigalian",
                                                     "Langhian","Serravallian","Tortonian","Messinian",
                                                     "Zanclean","Piacenzian","Gelasian","Calabrian",
                                                     "Chibanian","Late/Upper","Recent"))

# Save iteration data
save(Stage_null, file = "~/Mean_Stage_Null_metrics.RData")
