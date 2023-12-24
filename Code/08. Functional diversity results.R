###################################################################################################################
# 08. Functional diversity results 
## This R code provides the plots and results comparing the empirical analyses and the null model
## It further compares net changes in taxonomic richness between successive epochs in the empirical and resampled samples
## it produces Figure 2, Table 2, S4 and S5
###################################################################################################################

## Import packages
library(tidyverse)
library(deeptime)
library(writexl)
library(cowplot)
library(xfun)
library(Hmisc)
library(boot)

# Load data
load(file="~/Cleaned data.RData")

# Load empirical and null model results
## Empirical analyses
load(file = "~/Taxon_variation_metrics.RData")
load(file = "~/Taxon_variation_long_metrics.RData")
load(file = "~/Median_Taxon_metrics.RData")
## Null model
load(file = "~/Taxon_variation_null.RData")
load(file = "~/Taxon_variation_null_long.RData")

## Null model set up for violin plots
FEmetrics_null_Taxonvar <- FDmetrics_null_long_taxonvar %>% 
  filter(variable == "nb_fe")
FRedmetrics_null_Taxonvar <- FDmetrics_null_long_taxonvar %>% 
  filter(variable == "fred")
FOredmetrics_null_Taxonvar <- FDmetrics_null_long_taxonvar %>% 
  filter(variable == "fored")
FVulnmetrics_null_Taxonvar <- FDmetrics_null_long_taxonvar %>% 
  filter(variable == "fvuln")
FRicmetrics_null_Taxonvar <- FDmetrics_null_long_taxonvar %>% 
  filter(variable == "fric")
FOrimetrics_null_Taxonvar <- FDmetrics_null_long_taxonvar %>% 
  filter(variable == "fori")
FSpemetrics_null_Taxonvar <- FDmetrics_null_long_taxonvar %>% 
  filter(variable == "fspe")

# Plot Figure 2
## Produce geological time scale (for FSpe plots)
epochs_custom <- data.frame(
  name = c("Palaeocene","Eocene","Oligocene","Miocene","Pliocene","Pleistocene","Recent"), 
  max_age = c(1.5,2.5,3.5,4.5,5.5,6.5,7.6),
  min_age = c(0.4,1.5,2.5,3.5,4.5,5.5,6.5), 
  color = c("#FBA75F","#FDB46C","#FDC07A","#FFFF90","#FFFF99","#FFF2AE","#FEF2E0")
)

## Functional entities
FDmetrics_FE <- FDmetrics_long_TaxonVar %>% 
  filter(variable == "nb_fe")
FE_null_variation <- ggplot(data=FEmetrics_null_Taxonvar, aes(x= Epoch, y= value))+
  geom_violin(trim=TRUE, fill='grey', color="grey", alpha = 0.25)+
  scale_x_discrete(limits=c("Palaeocene","Eocene","Oligocene","Miocene", "Pliocene", "Pleistocene", "Recent"))+
  geom_boxplot(data = FDmetrics_FE, aes(x = Epoch, y = value),
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
FDmetrics_Fred <- FDmetrics_long_TaxonVar %>% 
  filter(variable == "fred")
FRed_null_variation <- ggplot(data=FRedmetrics_null_Taxonvar, aes(x= Epoch, y= value))+
  geom_violin(trim=TRUE, fill='grey', color="grey", alpha = 0.25)+
  scale_x_discrete(limits=c("Palaeocene","Eocene","Oligocene","Miocene", "Pliocene","Pleistocene","Recent"))+
  geom_boxplot(data = FDmetrics_Fred, aes(x = Epoch, y = value), 
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
FDmetrics_Fored <- FDmetrics_long_TaxonVar %>% 
  filter(variable == "fored")
FOred_null_variation <- ggplot(data=FOredmetrics_null_Taxonvar, aes(x= Epoch, y= value))+
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
  theme(axis.text = element_text(size= 6.5, color= "black"),
        axis.title= element_text(size= 8), panel.background= element_rect(fill= "white"))+
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent"))+
  theme(axis.text.x = element_blank())

## Functional vulnerability
FDmetrics_Fvuln <- FDmetrics_long_TaxonVar %>% 
  filter(variable == "fvuln")
FVuln_null_variation <- ggplot(data=FVulnmetrics_null_Taxonvar, aes(x= Epoch, y= value))+
  geom_violin(trim=TRUE, fill='grey', color="grey", alpha = 0.25)+
  scale_x_discrete(limits=c("Palaeocene","Eocene","Oligocene","Miocene", "Pliocene", "Pleistocene","Recent"))+
  geom_boxplot(data = FDmetrics_Fvuln, aes(x = Epoch, y = value), 
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
FDmetrics_FRic <- FDmetrics_long_TaxonVar %>% 
  filter(variable == "fric")
FRic_null_variation <- ggplot(data=FRicmetrics_null_Taxonvar, aes(x= Epoch, y= value))+
  geom_violin(trim=TRUE, fill='grey', color="grey", alpha = 0.25)+
  scale_x_discrete(limits=c("Palaeocene","Eocene","Oligocene","Miocene", "Pliocene", "Pleistocene","Recent"))+
  geom_boxplot(data = FDmetrics_FRic, aes(x = Epoch, y = value), 
               fill=c("#FBA75F","#FDB46C","#FDC07A","#FFFF90","#FFFF99","#FFF2AE","#FEF2E0"))+
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
  theme(axis.text.x = element_blank())

## Functional originality
FDmetrics_Fori <- FDmetrics_long_TaxonVar %>% 
  filter(variable == "fori")
FOri_null_variation <- ggplot(data=FOrimetrics_null_Taxonvar, aes(x= Epoch, y= value))+
  geom_violin(trim=TRUE, fill='grey', color="grey", alpha = 0.25)+
  scale_x_discrete(limits=c("Palaeocene","Eocene","Oligocene","Miocene", "Pliocene", "Pleistocene","Recent"))+
  geom_boxplot(data = FDmetrics_Fori, aes(x = Epoch, y = value), 
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
FDmetrics_Fspe <- FDmetrics_long_TaxonVar %>% 
  filter(variable == "fspe")
FSpe_null_variation <- ggplot(data=FSpemetrics_null_Taxonvar, aes(x= Epoch, y= value))+
  geom_violin(trim=TRUE, fill='grey', color="grey", alpha = 0.25)+
  scale_x_discrete(limits=c("Palaeocene","Eocene","Oligocene","Miocene", "Pliocene", "Pleistocene","Recent"))+
  geom_boxplot(data = FDmetrics_Fspe, aes(x = Epoch, y = value), 
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


# Plot everything together  - produces Figure 2
Fig_2 <- plot_grid(FE_null_variation,
          FRed_null_variation,FOred_null_variation,
          FRic_null_variation,
          FOri_null_variation,FSpe_null_variation,
          labels= c("(a)","(b)","(c)","(d)","(e)","(f)"), 
          label_size = 10,align = "hv", label_fontface = "bold",  nrow=6)


# Calculate Z scores between empirical analyses and null model per epoch
# Form empty dataframe
emp.null<- as.data.frame(matrix(data= NA,nrow= 7, ncol= 16, dimnames= list(c("Palaeocene", "Eocene", "Oligocene", "Miocene", "Pliocene", "Pleistocene","Recent"),
                                                                           c("nb_sp","nb_fe", "fred", "fored", "fvuln", "fric", "fori","fspe",
                                                                             "nb_sp_Z","nb_fe_Z", "fred_Z", "fored_Z", "fvuln_Z", "fric_Z", "fori_Z","fspe_Z"))))
# Separate by epoch for each analysis
## Empirical 
FDmetrics_taxonvar_Pal<- Res_FDmetrics_TaxonVar[grep("Palaeocene", Res_FDmetrics_TaxonVar$Epoch), (2:9)]
FDmetrics_taxonvar_Eo<- Res_FDmetrics_TaxonVar[grep("Eocene", Res_FDmetrics_TaxonVar$Epoch), (2:9)]
FDmetrics_taxonvar_Oli<- Res_FDmetrics_TaxonVar[grep("Oligocene", Res_FDmetrics_TaxonVar$Epoch), (2:9)]
FDmetrics_taxonvar_Mio<- Res_FDmetrics_TaxonVar[grep("Miocene", Res_FDmetrics_TaxonVar$Epoch), (2:9)]
FDmetrics_taxonvar_Plio<- Res_FDmetrics_TaxonVar[grep("Pliocene", Res_FDmetrics_TaxonVar$Epoch), (2:9)]
FDmetrics_taxonvar_Ple<- Res_FDmetrics_TaxonVar[grep("Pleistocene", Res_FDmetrics_TaxonVar$Epoch), (2:9)]
FDmetrics_taxonvar_Rec<- Res_FDmetrics_TaxonVar[grep("Recent", Res_FDmetrics_TaxonVar$Epoch), (2:9)]
## Null model
Null_FDmetrics_taxonvar_Pal<- Null_FDmetrics_taxonvar[grep("Palaeocene", Null_FDmetrics_taxonvar$Epoch), (2:9)]
Null_FDmetrics_taxonvar_Eo<- Null_FDmetrics_taxonvar[grep("Eocene", Null_FDmetrics_taxonvar$Epoch), (2:9)]
Null_FDmetrics_taxonvar_Oli<- Null_FDmetrics_taxonvar[grep("Oligocene", Null_FDmetrics_taxonvar$Epoch), (2:9)]
Null_FDmetrics_taxonvar_Mio<- Null_FDmetrics_taxonvar[grep("Miocene", Null_FDmetrics_taxonvar$Epoch), (2:9)]
Null_FDmetrics_taxonvar_Plio<- Null_FDmetrics_taxonvar[grep("Pliocene", Null_FDmetrics_taxonvar$Epoch), (2:9)]
Null_FDmetrics_taxonvar_Ple<- Null_FDmetrics_taxonvar[grep("Pleistocene", Null_FDmetrics_taxonvar$Epoch), (2:9)]
Null_FDmetrics_taxonvar_Rec<- Null_FDmetrics_taxonvar[grep("Recent", Null_FDmetrics_taxonvar$Epoch), (2:9)]

# Make empirical dataframes - means or medians of empirical & null models
FDemp <- Res_FDmetrics_TaxonVar %>% 
  group_by(Epoch) %>%
  summarise(Sp_med = median(nb_sp),
            FE_med = median(nb_fe),
            Red_med = median(fred),
            Ored_med = median(fored),
            Vul_med = median(fvuln),
            FRic_med = median(fric),
            Fori_med = median(fori),
            Fspe_med = median(fspe))
FDemp <- FDemp %>% 
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
    emp.null[g,e]<- as.numeric(FDemp[g,e])-as.numeric(FDnull[g,e])
    
  }
}


#write function to calculate differences between empirical and null model values
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
  emp.null[1,i+8]<- (emp.null[1,i]- median(N_Pal_null_Slopes[,i], na.rm = TRUE))/sd(N_Pal_null_Slopes[,i], na.rm= TRUE)}

for (i in 1:8){
  emp.null[2,i+8]<- (emp.null[2,i]- median(N_Eo_null_Slopes[,i],  na.rm = TRUE))/sd(N_Eo_null_Slopes[,i],  na.rm = TRUE)}

for (i in 1:8){
  emp.null[3,i+8]<- (emp.null[3,i]- median(N_Oli_null_Slopes[,i], na.rm = TRUE))/sd(N_Oli_null_Slopes[,i],  na.rm = TRUE)}

for (i in 1:8){
  emp.null[4,i+8]<- (emp.null[4,i]- median(N_Mio_null_Slopes[,i],  na.rm = TRUE))/sd(N_Mio_null_Slopes[,i],  na.rm = TRUE)}

for (i in 1:8){
  emp.null[5,i+8]<- (emp.null[5,i]- median(N_Plio_null_Slopes[,i],  na.rm = TRUE))/sd(N_Plio_null_Slopes[,i],  na.rm = TRUE)}

for (i in 1:8){
  emp.null[6,i+8]<- (emp.null[6,i]- median(N_Ple_null_Slopes[,i],  na.rm = TRUE))/sd(N_Ple_null_Slopes[,i],  na.rm = TRUE)}

for (i in 1:8){
  emp.null[7,i+8]<- (emp.null[7,i]- median(N_Rec_null_Slopes[,i],  na.rm = TRUE))/sd(N_Rec_null_Slopes[,i],  na.rm = TRUE)}


emp.null

# Extract results as excel worksheet - produces Tables 2 and S4
save(emp.null, file = "~/Empirical_null_FD metrics.Rdata")
write_xlsx(emp.null, "~/Empirical_null_FD differences.xlsx")

##################################################################################
## Calculate CI and significance of empirical taxonomic richness vs resampled taxonomic richness
# Load long empirical data
load(file = "~/Taxon_variation_long_metrics.RData")

# Load long resampled data
load(file = "~/Taxon_variation_resampling_long.RData")

## Filter by epochs
# Empirical
PalEo_empirical <- FDmetrics_long_TaxonVar %>% 
  filter(variable == "nb_sp") %>% 
  filter(Epoch == c("Palaeocene", "Eocene"))
EoOli_empirical <- FDmetrics_long_TaxonVar %>% 
  filter(variable == "nb_sp") %>% 
  filter(Epoch == c("Eocene", "Oligocene"))
OliMio_empirical <- FDmetrics_long_TaxonVar %>% 
  filter(variable == "nb_sp") %>% 
  filter(Epoch == c("Oligocene", "Miocene"))
MioPlio_empirical <- FDmetrics_long_TaxonVar %>% 
  filter(variable == "nb_sp") %>% 
  filter(Epoch == c("Miocene", "Pliocene"))
PlioPle_empirical <- FDmetrics_long_TaxonVar %>% 
  filter(variable == "nb_sp") %>% 
  filter(Epoch == c("Pliocene", "Pleistocene"))
PleRec_empirical <- FDmetrics_long_TaxonVar %>% 
  filter(variable == "nb_sp") %>% 
  filter(Epoch == c("Pleistocene", "Recent"))
# Resampled
PalEo_resampled <- FDmetrics_resamp_long_taxonvar %>% 
  filter(variable == "nb_sp") %>% 
  filter(Epoch == c("Palaeocene", "Eocene"))
EoOli_resampled <- FDmetrics_resamp_long_taxonvar %>% 
  filter(variable == "nb_sp") %>% 
  filter(Epoch == c("Eocene", "Oligocene"))
OliMio_resampled <- FDmetrics_resamp_long_taxonvar %>% 
  filter(variable == "nb_sp") %>% 
  filter(Epoch == c("Oligocene", "Miocene"))
MioPlio_resampled <- FDmetrics_resamp_long_taxonvar %>% 
  filter(variable == "nb_sp") %>% 
  filter(Epoch == c("Miocene", "Pliocene"))
PlioPle_resampled <- FDmetrics_resamp_long_taxonvar %>% 
  filter(variable == "nb_sp") %>% 
  filter(Epoch == c("Pliocene", "Pleistocene"))
PleRec_resampled <- FDmetrics_resamp_long_taxonvar %>% 
  filter(variable == "nb_sp") %>% 
  filter(Epoch == c("Pleistocene", "Recent"))

## Calculate net changes in each group
# Empirical
PalEo_emp <- PalEo_empirical %>%
  group_by(grp = ceiling(row_number() / 2)) %>%
  pivot_wider(names_from = Epoch, values_from = value) %>%
  mutate(diff_value = Eocene - Palaeocene) %>%
  select(grp, diff_value) %>%
  ungroup()
EoOli_emp <- EoOli_empirical %>%
  group_by(grp = ceiling(row_number() / 2)) %>%
  pivot_wider(names_from = Epoch, values_from = value) %>%
  mutate(diff_value = Oligocene - Eocene) %>%
  select(grp, diff_value) %>%
  ungroup()
OliMio_emp <- OliMio_empirical %>%
  group_by(grp = ceiling(row_number() / 2)) %>%
  pivot_wider(names_from = Epoch, values_from = value) %>%
  mutate(diff_value = Miocene - Oligocene) %>%
  select(grp, diff_value) %>%
  ungroup()
MioPlio_emp <- MioPlio_empirical %>%
  group_by(grp = ceiling(row_number() / 2)) %>%
  pivot_wider(names_from = Epoch, values_from = value) %>%
  mutate(diff_value = Pliocene - Miocene) %>%
  select(grp, diff_value) %>%
  ungroup()
PlioPle_emp <- PlioPle_empirical %>%
  group_by(grp = ceiling(row_number() / 2)) %>%
  pivot_wider(names_from = Epoch, values_from = value) %>%
  mutate(diff_value = Pleistocene - Pliocene) %>%
  select(grp, diff_value) %>%
  ungroup()
PleRec_emp <- PleRec_empirical %>%
  group_by(grp = ceiling(row_number() / 2)) %>%
  pivot_wider(names_from = Epoch, values_from = value) %>%
  mutate(diff_value = Recent - Pleistocene) %>%
  select(grp, diff_value) %>%
  ungroup()
# Resampling
PalEo_resamp <- PalEo_resampled %>%
  group_by(grp = ceiling(row_number() / 2)) %>%
  pivot_wider(names_from = Epoch, values_from = value) %>%
  mutate(diff_value = Eocene - Palaeocene) %>%
  select(grp, diff_value) %>%
  ungroup()
EoOli_resamp <- EoOli_resampled %>%
  group_by(grp = ceiling(row_number() / 2)) %>%
  pivot_wider(names_from = Epoch, values_from = value) %>%
  mutate(diff_value = Oligocene - Eocene) %>%
  select(grp, diff_value) %>%
  ungroup()
OliMio_resamp <- OliMio_resampled %>%
  group_by(grp = ceiling(row_number() / 2)) %>%
  pivot_wider(names_from = Epoch, values_from = value) %>%
  mutate(diff_value = Miocene - Oligocene) %>%
  select(grp, diff_value) %>%
  ungroup()
MioPlio_resamp <- MioPlio_resampled %>%
  group_by(grp = ceiling(row_number() / 2)) %>%
  pivot_wider(names_from = Epoch, values_from = value) %>%
  mutate(diff_value = Pliocene - Miocene) %>%
  select(grp, diff_value) %>%
  ungroup()
PlioPle_resamp <- PlioPle_resampled %>%
  group_by(grp = ceiling(row_number() / 2)) %>%
  pivot_wider(names_from = Epoch, values_from = value) %>%
  mutate(diff_value = Pleistocene - Pliocene) %>%
  select(grp, diff_value) %>%
  ungroup()
PleRec_resamp <- PleRec_resampled %>%
  group_by(grp = ceiling(row_number() / 2)) %>%
  pivot_wider(names_from = Epoch, values_from = value) %>%
  mutate(diff_value = Recent - Pleistocene) %>%
  select(grp, diff_value) %>%
  ungroup()

## Calculate differences
PalEo_med <- median(PalEo_emp$diff_value)
EoOli_med <- median(EoOli_emp$diff_value)
OliMio_med <- median(OliMio_emp$diff_value)
MioPlio_med <- median(MioPlio_emp$diff_value)
PlioPle_med <- median(PlioPle_emp$diff_value)
PleRec_med <- median(PleRec_emp$diff_value)

PalEo_net <- PalEo_resamp$diff_value-PalEo_med
EoOli_net <- EoOli_resamp$diff_value-EoOli_med
OliMio_net <- OliMio_resamp$diff_value-OliMio_med
MioPlio_net <- MioPlio_resamp$diff_value-MioPlio_med
PlioPle_net <- PlioPle_resamp$diff_value-PlioPle_med
PleRec_net <- PleRec_resamp$diff_value-PleRec_med

# Confidence interval tests
mean_cl_normal(PalEo_net)
mean_cl_normal(EoOli_net)
mean_cl_normal(OliMio_net)
mean_cl_normal(MioPlio_net)
mean_cl_normal(PlioPle_net)
mean_cl_normal(PleRec_net)

# Test for significance
## Function to calculate the median statistic
median_stat <- function(data, indices) {
  median(data[indices])
}

# Perform bootstrap for the median for each dataset
bootstrap_results_PalEo <- boot(PalEo_net, median_stat, R = 1000)
bootstrap_results_EoOli <- boot(EoOli_net, median_stat, R = 1000)
bootstrap_results_OliMio <- boot(OliMio_net, median_stat, R = 1000)
bootstrap_results_MioPlio <- boot(MioPlio_net, median_stat, R = 1000)
bootstrap_results_PlioPle <- boot(PlioPle_net, median_stat, R = 1000)
bootstrap_results_PleRec <- boot(PleRec_net, median_stat, R = 1000)

# Calculate p-values
observed_median_PalEo <- median(PalEo_net)
mean(bootstrap_results_PalEo$t >= observed_median_PalEo)

observed_median_EoOli <- median(EoOli_net)
mean(bootstrap_results_EoOli$t >= observed_median_EoOli)

observed_median_OliMio <- median(OliMio_net)
mean(bootstrap_results_OliMio$t >= observed_median_OliMio)

observed_median_MioPlio <- median(MioPlio_net)
mean(bootstrap_results_MioPlio$t >= observed_median_MioPlio)

observed_median_PlioPle <- median(PlioPle_net)
mean(bootstrap_results_PlioPle$t >= observed_median_PlioPle)

observed_median_PleRec <- median(PleRec_net)
mean(bootstrap_results_PleRec$t >= observed_median_PleRec)
