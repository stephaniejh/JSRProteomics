# MSstatsTMT_2025_09_v3.1.3_sp_w_iso_rerunNov2024.R

# This pulls in the dataframe from Nov 2024 analysis sp w iso (9855 seq).
# compared to v3.1.1.2 (sp w iso - 9878 seq)

# t1 - Nov 2024 PSMs data, with same settings as (08NOV2024 v1) & same as v3.1.2 t7
# Same results as 08NOV2024 ANS! 
# PD is the issue?

# 1 - LOAD PACKAGES ~  * ~  * ~ * ~ 
library(MSstatsTMT) 
library(MSstats) 
library(readr) 
library(dplyr) 


# Read in PSMs (same file used for 08NOV2024 - for ANS poster)
raw.pd <- read_tsv("input/Steph_01NOV2024_all_EE-SH_M-F_MSstatsTMT/NOV2024_PSMs__Steph_TMTpro-16plex_all.txt")

# Annotation File
annotation <- read_csv(
  "input/anno/Steph_TMT_all_EE&SH_F&M_final_annotation_file_v2_2025_08.csv")


# CHECKS
# total number of unique protein names
length(unique(raw.pd$"Protein Accessions")) #5170

# total number of unique protein names
length(unique(raw.pd$"Master Protein Accessions")) #4549

# total number of unique peptide names
length(unique(raw.pd$"Annotated Sequence")) #50040

# unique Spectrum.File, which is TMT run.
unique(raw.pd$"Spectrum File") # should print out 96 files

# no. of observations with contaminants TRUE vs FALSE
raw.pd %>%
  count(Contaminant) # 763459 FALSE & 4767 TRUE

# Remove contaminants!
# raw PSM pd file with contaminants removed. 
raw.pd.nc <- raw.pd %>% 
  filter(Contaminant == FALSE)

# no. of observations with contaminants TRUE vs FALSE
raw.pd.nc %>%
  count(Contaminant) # should be all FALSE

# check after removal
# total number of unique protein names
length(unique(raw.pd.nc$"Protein Accessions")) #5125

# total number of unique protein names
length(unique(raw.pd.nc$"Master Protein Accessions")) #4524

# total number of unique peptide names
length(unique(raw.pd.nc$"Annotated Sequence")) #49900





# 3 - PDtoMSstatsTMTFormat() ~  * ~  * ~ * ~ 
input.pd <- PDtoMSstatsTMTFormat(input = raw.pd.nc, 
                                 annotation = annotation,
                                 which.proteinid = "Master Protein Accessions",
                                 useNumProteinsColumn = TRUE, 
                                 useUniquePeptide = TRUE,
                                 rmPSM_withfewMea_withinRun = TRUE,
                                 rmProtein_with1Feature = FALSE,
                                 summaryforMultipleRows = sum,
                                 use_log_file = TRUE,
                                 append = FALSE)


# 4 - proteinSummarization() ~  * ~  * ~ * ~
quant.msstats <- proteinSummarization(input.pd,
                                      method="msstats",
                                      global_norm=TRUE,
                                      reference_norm=TRUE,
                                      remove_norm_channel = TRUE,
                                      remove_empty_channel = TRUE,
                                      MBimpute = TRUE, #only 4 MSstats method?
                                      maxQuantileforCensored = NULL,
                                      use_log_file = TRUE,
                                      append = FALSE) 

# # save psum
# save(quant.msstats,
#      file = 'output/MSstatsTMT/psum/2025_09_10_sp_w_iso_rrN24_Protein_Summarization_v3.1.3_t1.rda')
# 
# write_csv(quant.msstats$FeatureLevelData,
#           file='output/MSstatsTMT/psum/2025_09_10_sp_w_iso_rrN24_ProteinAbnd_featurelvl_v3.1.3_t1.csv')
# 
# write_csv(quant.msstats$ProteinLevelData,
#           file='output/MSstatsTMT/psum/2025_09_10_sp_w_iso_rrN24_ProteinAbnd_proteinlvl_v3.1.3_t1.csv')


# note w/ protein accession - model fitting for 3847 proteins (v3.1.3 t1)
# 5 - groupComparisonTMT() ~  * ~  * ~ * ~
test.pairwise.msstats <- groupComparisonTMT(data = quant.msstats, 
                                            contrast.matrix = "pairwise", # EE vs SH
                                            moderated = TRUE, # default F
                                            adj.method = "BH", # multiple compa adj
                                            save_fitted_models = TRUE) #otherwise not saved


# # save comp
# save(test.pairwise.msstats,
#      file = 'output/MSstatsTMT/comp/2025_09_10_sp_w_iso_rrN24_test_pairwise_msstats_v3.1.3_t1.rda')
# 
# #write_csv gives weird ' apostrophe by value but write.csv doesn't?
# write.csv(test.pairwise.msstats$ComparisonResult,
#           file='output/MSstatsTMT/comp/2025_09_10_sp_w_iso_rrN24_test_pairwise_msstats_v3.1.3_t1.csv')


# volcano plot
groupComparisonPlots(data = test.pairwise.msstats$ComparisonResult,
                     type="VolcanoPlot",
                     dot.size = 2,
                     text.size = 3,
                     width = 900,
                     height = 800,
                     address = "")