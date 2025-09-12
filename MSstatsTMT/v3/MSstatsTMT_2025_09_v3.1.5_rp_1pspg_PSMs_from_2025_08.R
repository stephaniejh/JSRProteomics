# MSstatsTMT_2025_09_v3.1.5_rp_1pspg_PSMs_from_2025_08
# 2025_09_12
# stephaniejh

# ~ * ` ~ * `
# This script uses the PSMs from 2025_08_Steph_TMTpro_PSMs_MSstatsTMT_EEvsSH_for_PhD which searched the runs in PD2.4
# using the SwissProt db with isoforms (>9878) sequences. 
# This is just a quick redo of the v2 data, but with just basic comparisons between EE vs SH. 

# This MSstatsTMT analysis is used to compare results with the normal SwissProt (v3.1.1)
# to 2025_08 rat proteome (1pspg >20000) (v2) + recheck again Nov 2024 ANS (sp+iso) (v1)

# t8 - all MSstatsTMT defaults.  EE vs SH 
# t9 - similar to Nov 2024 MSstats ANS. EE vs SH (raw.pd.nc, mas pro, mod = TRUE)

# This script is similar to MSstatsTMT_2025_09_v3.1.1.1_sp.Rmd but as a less detailed .R
# Just for a quick check! 



# 1 - LOAD PACKAGES ~  * ~  * ~ * ~ 
library(MSstatsTMT) 
library(MSstats) 
library(readr) 
library(dplyr) 



# 2 - LOAD FILES ~  * ~  * ~ * ~ 
# Set WD to source file location

# PSMS
mydir = 'input/2025_08_Steph_TMTpro_PSMs_MSstatsTMT_EEvsSH_for_PhD/txt_exports'

myfiles_PSMs = list.files(path = mydir,
                          pattern = ".*PSM.*\\.txt$",
                          full.names = TRUE)
#batch read
df_list_PSMs <- setNames(lapply(myfiles_PSMs, read_tsv),
                         myfiles_PSMs)

# bind to 1
raw.pd <- bind_rows(df_list_PSMs,
                    .id = 'pdresults_file_id') 

# # write out combined PSMs
write_tsv(raw.pd, "output/PSMs/combined_PSMs_2025_09_12_rp_1pspg.txt")
write_csv(raw.pd, "output/PSMs/combined_PSMs_2025_09_12_rp_1pspg.csv")


# Annotation File
annotation <- read_csv(
  "input/anno/Steph_TMT_all_EE&SH_F&M_final_annotation_file_v2_2025_08.csv")



# 3 - PDtoMSstatsTMTFormat() ~  * ~  * ~ * ~ 
# all defaults 
input.pd.t8 <- PDtoMSstatsTMTFormat(input = raw.pd, 
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
# all defaults
quant.msstats.t8 <- proteinSummarization(input.pd.t8,
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
# save(quant.msstats.t8,
#      file = 'output/MSstatsTMT/psum/2025_09_12_rp_1pspg_Protein_Summarization_v3.1.5_t8.rda')
# 
# write_csv(quant.msstats.t8$FeatureLevelData,
#           file='output/MSstatsTMT/psum/2025_09_12_rp_1pspg_ProteinAbnd_featurelvl_v3.1.5_t8.csv')
# 
# write_csv(quant.msstats.t8$ProteinLevelData,
#           file='output/MSstatsTMT/psum/2025_09_12_rp_1pspg_ProteinAbnd_proteinlvl_v3.1.5_t8.csv')


# note w/ protein accession - 
# 5 - groupComparisonTMT() ~  * ~  * ~ * ~
# all defaults
# testing for 7746 proteins
test.pairwise.msstats.t8 <- groupComparisonTMT(data = quant.msstats.t8, 
                                            contrast.matrix = "pairwise", # EE vs SH
                                            moderated = FALSE, # default F
                                            adj.method = "BH", # multiple compa adj
                                            save_fitted_models = TRUE) #otherwise not saved


# # save comp
# save(test.pairwise.msstats.t8,
#      file = 'output/MSstatsTMT/comp/2025_09_12_rp_1pspg_test_pairwise_msstats_v3.1.5_t8.rda')
# 
# #write_csv gives weird ' apostrophe by value but write.csv doesn't?
# write.csv(test.pairwise.msstats.t8$ComparisonResult,
#           file='output/MSstatsTMT/comp/2025_09_12_rp_1pspg_test_pairwise_msstats_v3.1.5_t8.csv')


# volcano plot
groupComparisonPlots(data = test.pairwise.msstats.t8$ComparisonResult,
                     type="VolcanoPlot",
                     dot.size = 2,
                     text.size = 3,
                     width = 950,
                     height = 900,
                     address = "")



##  * ~ ` *  ~ * ` * ~ ` *  ~ * ` * ~ ` *  ~ * ` * ~ ` *  ~ * `
## ~ T9  - Similar to Nov2024 ANS MSstatsTMT analysis
# - NC & master protein & Mod = TRUE

# no. of observations with contaminants TRUE vs FALSE
raw.pd %>%
  count(Contaminant) # 

# raw PSM pd file with contaminants removed. 
raw.pd.nc <- raw.pd %>% 
  filter(Contaminant == FALSE)

# no. of observations with contaminants TRUE vs FALSE
raw.pd.nc %>%
  count(Contaminant) # should be all FALSE

input.pd.t9 <- PDtoMSstatsTMTFormat(input = raw.pd.nc, # t9 no contaminants
                                 annotation = annotation,
                                 which.proteinid = "Master Protein Accessions", # t9 master
                                 useNumProteinsColumn = TRUE, 
                                 useUniquePeptide = TRUE,
                                 rmPSM_withfewMea_withinRun = TRUE,
                                 rmProtein_with1Feature = FALSE,
                                 summaryforMultipleRows = sum,
                                 use_log_file = TRUE,
                                 append = FALSE)


quant.msstats.t9 <- proteinSummarization(input.pd.t9,
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
# save(quant.msstats.t9,
#      file = 'output/MSstatsTMT/psum/2025_09_12_rp_1pspg_Protein_Summarization_v3.1.5_t9.rda')
# 
# write_csv(quant.msstats.t9$FeatureLevelData,
#           file='output/MSstatsTMT/psum/2025_09_12_rp_1pspg_ProteinAbnd_featurelvl_v3.1.5_t9.csv')
# 
# write_csv(quant.msstats.t9$ProteinLevelData,
#           file='output/MSstatsTMT/psum/2025_09_12_rp_1pspg_ProteinAbnd_proteinlvl_v3.1.5_t9.csv')


# t9 - model fits for 8048 proteins (why is it more for master protein?)
test.pairwise.msstats.t9 <- groupComparisonTMT(data = quant.msstats.t9, 
                                            contrast.matrix = "pairwise", #EE vs SH
                                            moderated = TRUE, # t9
                                            adj.method = "BH", # multiple compa adj
                                            save_fitted_models = TRUE) #otherwise not saved


# # save comp
# save(test.pairwise.msstats.t9,
#      file = 'output/MSstatsTMT/comp/2025_09_09_12_rp_1pspg_test_pairwise_msstats_v3.1.5_t9.rda')
# 
# #write_csv gives weird ' apostrophe by value but write.csv doesn't?
# write.csv(test.pairwise.msstats.t9$ComparisonResult,
#           file='output/MSstatsTMT/comp/2025_09_09_12_rp_1pspg_test_pairwise_msstats_v3.1.5_t9.csv')


# volcano plot
groupComparisonPlots(data = test.pairwise.msstats.t9$ComparisonResult,
                     type="VolcanoPlot",
                     dot.size = 2,
                     text.size = 3,
                     width = 950,
                     height = 850,
                     address = "")
