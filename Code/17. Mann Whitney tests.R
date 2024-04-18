########################################################################################################################################
# 17. Mann-Whitney tests
## This R code calculates pairwise Mann-Whitney u-tests between epochs
## it produces supplementary table S8
#######################################################################################################################################

## Import packages
library(tidyverse)

# Load empirical data
load(file = "~/Taxon_variation_metrics.RData")
load(file = "~/Taxon_variation_long_metrics.RData")

# Filter by FD metric
FDmetrics_FE <- FDmetrics_long_TaxonVar %>% 
  filter(variable == "nb_fe")
FDmetrics_Fred <- FDmetrics_long_TaxonVar %>% 
  filter(variable == "fred")
FDmetrics_Fored <- FDmetrics_long_TaxonVar %>% 
  filter(variable == "fored")
FDmetrics_Fvuln <- FDmetrics_long_TaxonVar %>% 
  filter(variable == "fvuln")
FDmetrics_FRic <- FDmetrics_long_TaxonVar %>% 
  filter(variable == "fric")
FDmetrics_Fori <- FDmetrics_long_TaxonVar %>% 
  filter(variable == "fori")
FDmetrics_Fspe <- FDmetrics_long_TaxonVar %>% 
  filter(variable == "fspe")

# Filter by epoch
FE_Taxonvar_pal <- FDmetrics_FE %>% 
  filter(Epoch == "Palaeocene")
FRed_Taxonvar_pal <- FDmetrics_Fred %>% 
  filter(Epoch == "Palaeocene")
FOred_Taxonvar_pal <- FDmetrics_Fored %>% 
  filter(Epoch == "Palaeocene")
FRic_Taxonvar_pal <- FDmetrics_FRic %>% 
  filter(Epoch == "Palaeocene")
FOri_Taxonvar_pal <- FDmetrics_Fori %>% 
  filter(Epoch == "Palaeocene")
FSpe_Taxonvar_pal <- FDmetrics_Fspe %>% 
  filter(Epoch == "Palaeocene")

FE_Taxonvar_eo <- FDmetrics_FE %>% 
  filter(Epoch == "Eocene")
FRed_Taxonvar_eo <- FDmetrics_Fred %>% 
  filter(Epoch == "Eocene")
FOred_Taxonvar_eo <- FDmetrics_Fored %>% 
  filter(Epoch == "Eocene")
FRic_Taxonvar_eo <- FDmetrics_FRic %>% 
  filter(Epoch == "Eocene")
FOri_Taxonvar_eo <- FDmetrics_Fori %>% 
  filter(Epoch == "Eocene")
FSpe_Taxonvar_eo <- FDmetrics_Fspe %>% 
  filter(Epoch == "Eocene")

FE_Taxonvar_oli <- FDmetrics_FE %>% 
  filter(Epoch == "Oligocene")
FRed_Taxonvar_oli <- FDmetrics_Fred %>% 
  filter(Epoch == "Oligocene")
FOred_Taxonvar_oli <- FDmetrics_Fored %>% 
  filter(Epoch == "Oligocene")
FRic_Taxonvar_oli <- FDmetrics_FRic %>% 
  filter(Epoch == "Oligocene")
FOri_Taxonvar_oli <- FDmetrics_Fori %>% 
  filter(Epoch == "Oligocene")
FSpe_Taxonvar_oli <- FDmetrics_Fspe %>% 
  filter(Epoch == "Oligocene")

FE_Taxonvar_mio <- FDmetrics_FE %>% 
  filter(Epoch == "Miocene")
FRed_Taxonvar_mio <- FDmetrics_Fred %>% 
  filter(Epoch == "Miocene")
FOred_Taxonvar_mio <- FDmetrics_Fored %>% 
  filter(Epoch == "Miocene")
FRic_Taxonvar_mio <- FDmetrics_FRic %>% 
  filter(Epoch == "Miocene")
FOri_Taxonvar_mio <- FDmetrics_Fori %>% 
  filter(Epoch == "Miocene")
FSpe_Taxonvar_mio <- FDmetrics_Fspe %>% 
  filter(Epoch == "Miocene")

FE_Taxonvar_plio <- FDmetrics_FE %>% 
  filter(Epoch == "Pliocene")
FRed_Taxonvar_plio <- FDmetrics_Fred %>% 
  filter(Epoch == "Pliocene")
FOred_Taxonvar_plio <- FDmetrics_Fored %>% 
  filter(Epoch == "Pliocene")
FRic_Taxonvar_plio <- FDmetrics_FRic %>% 
  filter(Epoch == "Pliocene")
FOri_Taxonvar_plio <- FDmetrics_Fori %>% 
  filter(Epoch == "Pliocene")
FSpe_Taxonvar_plio <- FDmetrics_Fspe %>% 
  filter(Epoch == "Pliocene")

FE_Taxonvar_ple <- FDmetrics_FE %>% 
  filter(Epoch == "Pleistocene")
FRed_Taxonvar_ple <- FDmetrics_Fred %>% 
  filter(Epoch == "Pleistocene")
FOred_Taxonvar_ple <- FDmetrics_Fored %>% 
  filter(Epoch == "Pleistocene")
FRic_Taxonvar_ple <- FDmetrics_FRic %>% 
  filter(Epoch == "Pleistocene")
FOri_Taxonvar_ple <- FDmetrics_Fori %>% 
  filter(Epoch == "Pleistocene")
FSpe_Taxonvar_ple <- FDmetrics_Fspe %>% 
  filter(Epoch == "Pleistocene")

FE_Taxonvar_rec <- FDmetrics_FE %>% 
  filter(Epoch == "Recent")
FRed_Taxonvar_rec <- FDmetrics_Fred %>% 
  filter(Epoch == "Recent")
FOred_Taxonvar_rec <- FDmetrics_Fored %>% 
  filter(Epoch == "Recent")
FRic_Taxonvar_rec <- FDmetrics_FRic %>% 
  filter(Epoch == "Recent")
FOri_Taxonvar_rec <- FDmetrics_Fori %>% 
  filter(Epoch == "Recent")
FSpe_Taxonvar_rec <- FDmetrics_Fspe %>% 
  filter(Epoch == "Recent")

## Mann-Whitney tests
# FE
wilcox.test(FE_Taxonvar_pal$value,FE_Taxonvar_eo$value)
wilcox.test(FE_Taxonvar_eo$value,FE_Taxonvar_oli$value)
wilcox.test(FE_Taxonvar_oli$value,FE_Taxonvar_mio$value)
wilcox.test(FE_Taxonvar_mio$value,FE_Taxonvar_plio$value)
wilcox.test(FE_Taxonvar_plio$value,FE_Taxonvar_ple$value)
wilcox.test(FE_Taxonvar_ple$value,FE_Taxonvar_rec$value)

# FRed
wilcox.test(FRed_Taxonvar_pal$value,FRed_Taxonvar_eo$value)
wilcox.test(FRed_Taxonvar_eo$value,FRed_Taxonvar_oli$value)
wilcox.test(FRed_Taxonvar_oli$value,FRed_Taxonvar_mio$value)
wilcox.test(FRed_Taxonvar_mio$value,FRed_Taxonvar_plio$value)
wilcox.test(FRed_Taxonvar_plio$value,FRed_Taxonvar_ple$value)
wilcox.test(FRed_Taxonvar_ple$value,FRed_Taxonvar_rec$value)

# FOred
wilcox.test(FOred_Taxonvar_pal$value,FOred_Taxonvar_eo$value)
wilcox.test(FOred_Taxonvar_eo$value,FOred_Taxonvar_oli$value)
wilcox.test(FOred_Taxonvar_oli$value,FOred_Taxonvar_mio$value)
wilcox.test(FOred_Taxonvar_mio$value,FOred_Taxonvar_plio$value)
wilcox.test(FOred_Taxonvar_plio$value,FOred_Taxonvar_ple$value)
wilcox.test(FOred_Taxonvar_ple$value,FOred_Taxonvar_rec$value)

# FRic
wilcox.test(FRic_Taxonvar_pal$value,FRic_Taxonvar_eo$value)
wilcox.test(FRic_Taxonvar_eo$value,FRic_Taxonvar_oli$value)
wilcox.test(FRic_Taxonvar_oli$value,FRic_Taxonvar_mio$value)
wilcox.test(FRic_Taxonvar_mio$value,FRic_Taxonvar_plio$value)
wilcox.test(FRic_Taxonvar_plio$value,FRic_Taxonvar_ple$value)
wilcox.test(FRic_Taxonvar_ple$value,FRic_Taxonvar_rec$value)

# FOri
wilcox.test(FOri_Taxonvar_pal$value,FOri_Taxonvar_eo$value)
wilcox.test(FOri_Taxonvar_eo$value,FOri_Taxonvar_oli$value)
wilcox.test(FOri_Taxonvar_oli$value,FOri_Taxonvar_mio$value)
wilcox.test(FOri_Taxonvar_mio$value,FOri_Taxonvar_plio$value)
wilcox.test(FOri_Taxonvar_plio$value,FOri_Taxonvar_ple$value)
wilcox.test(FOri_Taxonvar_ple$value,FOri_Taxonvar_rec$value)

# FSpe
wilcox.test(FSpe_Taxonvar_pal$value,FSpe_Taxonvar_eo$value)
wilcox.test(FSpe_Taxonvar_eo$value,FSpe_Taxonvar_oli$value)
wilcox.test(FSpe_Taxonvar_oli$value,FSpe_Taxonvar_mio$value)
wilcox.test(FSpe_Taxonvar_mio$value,FSpe_Taxonvar_plio$value)
wilcox.test(FSpe_Taxonvar_plio$value,FSpe_Taxonvar_ple$value)
wilcox.test(FSpe_Taxonvar_ple$value,FSpe_Taxonvar_rec$value)
