#MSstatsTMT_2025_09_v3.1.2_sp_w_iso
# 2025_09_09
# stephaniejh

# ~ * ` ~ * `
# This script uses the PSMs from 2025_09_sp_w_iso which searched the runs in PD2.4
# using the SwissProt db with isoforms (>9878) sequences. 
# This MSstatsTMT analysis is used to compare results with the normal SwissProt (v3.1.1)
# to 2025_08 rat proteome (1pspg >20000) (v2) + recheck again Nov 2024 ANS (sp+iso) (v1)

# t3 - sp w iso - MSstatsTMT defaults - just basic check - housing ee vs sh
# t4 - same except pairwise comp does moderated

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
mydir = 'input/2025_09_sp_w_iso_Steph_TMTpro_PSMs_MSstatsTMT_EEvsSH_for_PhD/txt_exports'

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
# write_tsv(raw.pd, "output/PSMs/combined_PSMs_2025_09_05_sp_w_iso.txt")
# write_csv(raw.pd, "output/PSMs/combined_PSMs_2025_09_05_sp_w_iso.csv")


# Annotation File
annotation <- read_csv(
  "input/anno/Steph_TMT_all_EE&SH_F&M_final_annotation_file_v2_2025_08.csv")



# 3 - PDtoMSstatsTMTFormat() ~  * ~  * ~ * ~ 
# all defaults 
input.pd <- PDtoMSstatsTMTFormat(input = raw.pd, 
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
#      file = 'output/MSstatsTMT/psum/2025_09_09_sp_w_iso_Protein_Summarization_v3.1.2_t3.rda')
# 
# write_csv(quant.msstats$FeatureLevelData,
#           file='output/MSstatsTMT/psum/2025_09_09_sp_w_iso_ProteinAbnd_featurelvl_v3.1.2_t3.csv')
# 
# write_csv(quant.msstats$ProteinLevelData,
#           file='output/MSstatsTMT/psum/2025_09_09_sp_w_iso_ProteinAbnd_proteinlvl_v3.1.2_t3.csv')


# note w/ protein accession - model fitting for 4230 proteins
# 5 - groupComparisonTMT() ~  * ~  * ~ * ~
# all defaults
test.pairwise.msstats <- groupComparisonTMT(data = quant.msstats, 
                                            contrast.matrix = "pairwise", # 3 - EE vs SH
                                            moderated = FALSE, # default F
                                            adj.method = "BH", # multiple compa adj
                                            save_fitted_models = TRUE) #otherwise not saved


# # save comp
# save(test.pairwise.msstats,
#      file = 'output/MSstatsTMT/comp/2025_09_09_sp_w_iso_test_pairwise_msstats_v3.1.2_t3.rda')
# 
# #write_csv gives weird ' apostrophe by value but write.csv doesn't?
# write.csv(test.pairwise.msstats$ComparisonResult,
#           file='output/MSstatsTMT/comp/2025_09_09_sp_w_iso_test_pairwise_msstats_v3.1.2_t3.csv')


# volcano plot
groupComparisonPlots(data = test.pairwise.msstats$ComparisonResult,
                     type="VolcanoPlot",
                     dot.size = 2,
                     text.size = 3,
                     width = 800,
                     height = 800,
                     address = "")



## ~ T4  -  all same except try moderated = TRUE ~ ~ ~ ~ ~ 
# load in protein sum from above
# groupComparisonTMT() ~  * ~  * ~ * ~
# all defaults except 
test.pairwise.msstats.t4 <- groupComparisonTMT(data = quant.msstats, 
                                            contrast.matrix = "pairwise", # 3 - EE vs SH
                                            # moderated = FALSE, # default F
                                            moderated = TRUE, # t4
                                            adj.method = "BH", # multiple compa adj
                                            save_fitted_models = TRUE) #otherwise not saved


# # save comp
# save(test.pairwise.msstats.t4,
#      file = 'output/MSstatsTMT/comp/2025_09_09_sp_w_iso_test_pairwise_msstats_v3.1.2_t4.rda')
# 
# #write_csv gives weird ' apostrophe by value but write.csv doesn't?
# write.csv(test.pairwise.msstats.t4$ComparisonResult,
#           file='output/MSstatsTMT/comp/2025_09_09_sp_w_iso_test_pairwise_msstats_v3.1.2_t4.csv')


# volcano plot
groupComparisonPlots(data = test.pairwise.msstats.t4$ComparisonResult,
                     type="VolcanoPlot",
                     dot.size = 2,
                     text.size = 3,
                     width = 800,
                     height = 800,
                     address = "")



