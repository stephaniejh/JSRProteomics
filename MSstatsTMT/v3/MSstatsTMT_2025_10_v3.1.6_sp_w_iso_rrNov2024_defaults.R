# MSstatsTMT_2025_10_v3.1.6_sp_w_iso_rrNov2024_defaults.R
# 2025_10_02
# stephaniejh

# ~ * ` ~ * `
# This pulls in the dataframe from Nov 2024 analysis sp w iso (9855 seq).


# t11 - runs Nov 2024 sp w iso data with all MSstatsTMT defaults (EE vs SH)
# to compare with 2025 e.g. t8 1pspg etc, 


# 1 - LOAD PACKAGES ~  * ~  * ~ * ~ 
library(MSstatsTMT) 
library(MSstats) 
library(readr) 
library(dplyr) 



# 2 - LOAD FILES ~  * ~  * ~ * ~ 
# Set WD to source file location

# Read in PSMs (same file used for 08NOV2024 - for ANS poster) - 768226 x 54
raw.pd <- read_tsv("input/Steph_01NOV2024_all_EE-SH_M-F_MSstatsTMT/NOV2024_PSMs__Steph_TMTpro-16plex_all.txt")


# Annotation File
annotation <- read_csv(
  "input/anno/Steph_TMT_all_EE&SH_F&M_final_annotation_file_v2_2025_08.csv")



# 3 - PDtoMSstatsTMTFormat() ~  * ~  * ~ * ~ 
# all defaults - 5746608 x 11
input.pd.t11 <- PDtoMSstatsTMTFormat(input = raw.pd, 
                                 annotation = annotation,
                                 which.proteinid = "Protein Accessions",
                                 useNumProteinsColumn = TRUE, 
                                 useUniquePeptide = TRUE,
                                 rmPSM_withfewMea_withinRun = TRUE,
                                 rmProtein_with1Feature = FALSE,
                                 summaryforMultipleRows = sum,
                                 use_log_file = TRUE,
                                 append = FALSE)


# 4 - proteinSummarization() ~  * ~  * ~ * ~
# all defaults | Protein level data - 472208 x 8
quant.msstats.t11 <- proteinSummarization(input.pd.t11,
                                      method="msstats",
                                      global_norm=TRUE,
                                      reference_norm=TRUE,
                                      remove_norm_channel = TRUE,
                                      remove_empty_channel = TRUE,
                                      MBimpute = TRUE, #only 4 MSstats method?
                                      maxQuantileforCensored = NULL,
                                      use_log_file = TRUE,
                                      append = FALSE) 

# # # save psum
# save(quant.msstats.t11,
#      file = 'output/MSstatsTMT/psum/2025_10_02_sp_w_iso_rrN24_def_Protein_Summarization_v3.1.6_t11.rda')
# 
# write_csv(quant.msstats.t11$FeatureLevelData,
#           file='output/MSstatsTMT/psum/2025_10_02_sp_w_iso_rrN24_def_ProteinAbnd_featurelvl_v3.1.6_t11.csv')
# 
# write_csv(quant.msstats.t11$ProteinLevelData,
#           file='output/MSstatsTMT/psum/2025_10_02_sp_w_iso_rrN24_def_ProteinAbnd_proteinlvl_v3.1.6_t11.csv')


# note w/ protein accession - 
# 5 - groupComparisonTMT() ~  * ~  * ~ * ~
# all defaults
# testing for 3508 proteins
test.pairwise.msstats.t11 <- groupComparisonTMT(data = quant.msstats.t11, 
                                            contrast.matrix = "pairwise", # EE vs SH
                                            moderated = FALSE, # default F
                                            adj.method = "BH", # multiple compa adj
                                            save_fitted_models = TRUE) #otherwise not saved


# # save comp
# save(test.pairwise.msstats.t11,
#      file = 'output/MSstatsTMT/comp/2025_10_02_sp_w_iso_rrN24_def_test_pairwise_msstats_v3.1.6_t11.rda')
# 
# #write_csv gives weird ' apostrophe by value but write.csv doesn't?
# write.csv(test.pairwise.msstats.t11$ComparisonResult,
#           file='output/MSstatsTMT/comp/2025_10_02_sp_w_iso_rrN24_def_test_pairwise_msstats_v3.1.6_t11.csv')


# volcano plot
groupComparisonPlots(data = test.pairwise.msstats.t11$ComparisonResult,
                     type="VolcanoPlot",
                     dot.size = 2,
                     text.size = 3,
                     width = 950,
                     height = 900,
                     address = "")



