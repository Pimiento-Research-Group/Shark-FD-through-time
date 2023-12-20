# R code for the analysis of: The rise and fall of shark functional diversity over the last 66 million years

## Authors
Jack A. Cooper and Catalina Pimiento

## Introduction
This folder contains R scripts for the functional diversity analysis of sharks through the Cenozoic era. 

The following packages are required for:
- Data manipulation & handling: ```readxl```, ```tidyverse```, ```doBy```, ```reshape2```, ```tibble```, ```data.table```, ```writexl``` and ```xfun```.
- Analyses: ```DescTools```, ```corrr```, ```mFD```, ```foreach```, ```parallel```, ```doParallel```, ```Hmisc``` and ```boot```.
- Visualisation: ```ggsci```, ```scales```, ```cowplot```, ```RColorBrewer```, ```deeptime```, ```grid```, ```gridExtra``` and ```patchwork```.

These packages can be downloaded and installed using the following commands:
``` {r}
install.packages(readxl)
install.packages(tidyverse)
install.packages(doBy)
install.packages(reshape2)
install.packages(DescTools)
install.packages(corrr)
install.packages(ggsci)
install.packages(tibble)
install.packages(scales)
install.packages(cowplot)
install.packages(RColorBrewer)
install.packages(mFD)
install.packages(data.table)
install.packages(deeptime)
install.packages(grid)
install.packages(gridExtra)
install.packages(foreach)
install.packages(parallel)
install.packages(doParallel)
install.packages(writexl)
install.packages(xfun)
install.packages(Hmisc)
install.packages(boot)
install.packages(patchwork)
```

## Codes
A total of 13 R scripts were produced containing all analyses. The scripts are as follows:

01. **Data treatment**: This code prepares the shark tooth data (Data S1) for analyses in all other codes
- It produces 1 Rdata file: Cleaned data.RData 
02. **Data exploration**: This code explores taxonomic data and geographic distribution of the data; as well as testing correlations between dental characters and tooth position
- It produces Figure S1 and Table S2
03. **Functional space**: This code produces the functional space of Cenozoic sharks based on their dental characters
- It produces Figure 1, Figure S2-S4 and Table S3
- It further produces 2 Rdata files: one for excluding Recent sharks without a fossil record - Filtered taxa.Rdata, and one of the functional space's species-trait matrix - Species-trait matrix.Rdata
04. **Empirical functional diversity analyses**: This code produces functional diversity analyses based on our recorded tooth data
- It produces 3 Rdata files - the functional diversity metric results in wide and long form: Taxon_variation_metrics.RData, Taxon_variation_long_metrics.RData; as well as the median and standard deviation results: Mean_Taxon_metrics.RData
05. **Null model**: This code produces the null model used to assess expected functional diversity metrics based on taxonomic richness
- It produces 2 Rdata files - the null model results in wide and long form: Taxon_variation_null.RData, Taxon_variation_null_long.RData
06. **Extant shark functional space**: This code produces the functional space and functional richness/specialisation results of the total diversity of Recent shark species recorded in Data S2
- It produces Figure S8
07. **Randomisation results**: This code produces the random simulations on shark functional diversity data; one where time bin data are resampled based on the smallest sample size, and one simulating randomised taxonomic losses to assess the relationship between taxonomic richness and functional richness
- It produces Figure S5 and S6
- It further produces 3 Rdata files: the results of the resampling procedure in wide and long format - Taxon_variation_resampling_metrics.RData, Taxon_variation_resampling_long.RData - and the taxonomic loss sequence - Species richness sequence.RData
08. **Functional diversity results**: This code plots shark functional diversity through time, assesses if this deviates from null expectations via Z-scores, and examines if empirical changes in taxonomic richness are significantly different from when this data are resampled
- It produces Figure 2, Table 2, S4 and S5 
- It further produces an Rdata and xlsx file used to produce Table 2 - Empirical_null_FD_metrics.Rdata and Empirical_null_FD differences.xlsx
09. **Recent resampling analyses**: This code assesses if our Recent results are reliable by re-analysing Recent functional diversity with all Recent teeth and random resampling
- It produces 4 Rdata files - 2 providing the results of the "Recent-plus" sampling in wide and long format - Full_Recent.Rdata, Full_Recent_long.Rdata - and 2 providing the results of the "Recent-plus resampled" sampling in wide and long format - Recent_resampling.Rdata, Recent_resampling_long.Rdata
10. **Recent resampling results**: This code plots the Recent resampling results and analyses their deviation from null expectations
- It produces Figure S7
- It further produces 2 Rdata files of a null model of the "Recent-plus" sample in wide and long format - Full_Recent_null.Rdata, Full_Recent_null_long.Rdata
11. **Sensitivity tests**: This code produces sensitivity tests in which all functional diversity analyses are re-done with individual traits removed.
- It produces Figure S9
- It further produces 18 Rdata files - the functional diversity results in wide, long and median/standard deviation format for each of the six dental characters examined
12. **Contribution analyses**: This code analyses whether extinct or extant sharks contribute more to functional diversity, and assesses the most functionally original and specialised sharks of the Cenozoic
- It produces Figure 3 and S10
- It further produces 1 Rdata file: FUn_Fsp metrics.RData
13. **FD contributors through time**: This code assesses functional diversity contributions of sharks through time
- It reveals which sharks were the most functionally original and specialised for each time bin

All code and data are placed in the following folders:

- **Code**: containing all 13 R scripts
- **Data**: containing all Rdata and xlsx outputs from the code. Rdata files can be loaded directly into replications of our analyses using the load() command
