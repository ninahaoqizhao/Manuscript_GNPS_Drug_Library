library(tidyverse)
library(readxl)
library(httr)

##################################################################
# STEP ONE: SEARCH FOR DRUG ANALOGS ON METABOLOMICS REPOSITORIES #
##################################################################

# Run repository-scale search of the GNPS Drug Library with fastMASST in batch: https://github.com/robinschmid/microbe_masst/tree/v1.2.0
# Enable analog mode in the fastMASST search
# Output for each spectrum is a tsv files containing all the USI matched to the spectrum
# The resulting files can be downloaded from: 10.5281/zenodo.16423040

# Read in the analog fastMASST results
folder_path <- "/Users/ninazhao/PycharmProjects/masst_only/output/20250310_DrugLib_Updated_Analog_Search"
file_list <- list.files(folder_path, pattern = "_analog_matches.tsv$", full.names = TRUE)
process_file <- function(file) {
  df <- read_tsv(file, col_types = cols(.default = "c"))
  df <- df %>%
    mutate(`Delta Mass` = as.numeric(`Delta Mass`))
  df_filtered <- df %>%
    filter(abs(`Delta Mass`) > 0.9)
  ccmslib_value <- str_extract(basename(file), "CCMSLIB\\d+")
  df_filtered <- df_filtered %>%
    mutate(CCMSLIB = ccmslib_value)
  return(df_filtered)
} # filter for matches with delta mass > 0.9 (so that it is an analog)
combined_df <- bind_rows(lapply(file_list, process_file))

# read in druglib metadata
druglib_metadata <- read.csv("druglib_metadata/GNPS_Drug_Library_Metadata_Drugs.csv")
druglib_metadata_filtered <- druglib_metadata[!duplicated(druglib_metadata$gnps_libid),]
druglib_metadata_filtered <- druglib_metadata_filtered %>%
  mutate(split_inchikey = sub("-.*", "", inchikey))
druglib_metadata_exogenous <- druglib_metadata_filtered %>% dplyr::filter(!grepl("Background|confidence|Endogenous|Food", chemical_source))
# how many has analog
length(unique(combined_df$CCMSLIB))/length(unique(all_exogenous_drug$gnps_libid))


##############################################
# STEP TWO: FILTER BY EXPLAINABLE DELTA MASS #
##############################################

explained_delta_mass <- read_xlsx("AnalogLib_Development/Delta_Mass_filters/Delta_mass_inclusion.xlsx") # delta mass relevant to drug metabolism or mass spectrometry adducts 
all_detected_delta <- data.frame(detected_delta = sort(unique(combined_df$Delta.Mass)))
all_detected_delta$include <- rep(FALSE, nrow(all_detected_delta))
for (i in 1:nrow(all_detected_delta)){
  for (j in 1:nrow(explained_delta_mass)){
    if (abs(all_detected_delta$detected_delta[i] - explained_delta_mass$`mz delta (=drug - analog)`[j]) < 0.005){
      all_detected_delta$include[i] <- TRUE
    }
  }
} # find detect delta mass within 0.01 Da of explainable delta mass

# create a dataframe for all filtering information
include_delta <- all_detected_delta$detected_delta[all_detected_delta$include == TRUE]
new_analog_all <- combined_df %>% left_join(druglib_metadata_filtered, by = c("CCMSLIB" = "gnps_libid"))

# exclude delta mass not in the inclusion list
new_analog_all_delta_mass_filtered <- new_analog_all %>% dplyr::filter(new_analog_all$Delta.Mass %in% include_delta)
USI_to_download <- data.frame(USI = unique(new_analog_all_delta_mass_filtered$USI))
#write.csv(USI_to_download, "analogMASST/20250313_analog_reprocess_William_code/20250313_USI_to_download.csv")

#####################################################
# STEP THREE: DOWNLOAD SPECTRA OF THE REMAINING USI #
#####################################################
# achieved with matchms; see "AnalogLib_Development/USI_download/DrugAnalog_usi_download_.ipynb"
# results: "AnalogLib_Development/USI_download/Drug_analog_USI_spectra_RAW.mgf"

# double-check the modified cosine using matchms.similarity.ModifiedCosine
analog_usi_all <- readLines("AnalogLib_Development/USI_download/Drug_analog_USI_spectra_RAW.mgf") 
analog_usi_info <- data.frame(
  analog_scans = as.numeric(sub("SCANS=", "", analog_usi_all[grepl("SCANS=", analog_usi_all)])),
  analog_pepmass = as.numeric(sub("PEPMASS=", "", analog_usi_all[grepl("PEPMASS=", analog_usi_all)])),
  analog_usi = sub("TITLE=", "", analog_usi_all[grepl("TITLE=", analog_usi_all)]),
  analog_charge = rep(1, 433274)
)
drug_spectra <- readLines("final_mgf/GNPS_Drug_Library_Spectra_Drugs_and_Metabolites.mgf")
drug_spectra_info <- data.frame(
  drug_pepmass = as.numeric(sub("PEPMASS=", "", drug_spectra[grepl("PEPMASS=", drug_spectra)])),
  drug_CCMSLIBID = sub("SPECTRUMID=", "", drug_spectra[grepl("SPECTRUMID=", drug_spectra)]),
  drug_charge = sub("CHARGE=", "", drug_spectra[grepl("CHARGE=", drug_spectra)])
)
drug_spectra_info$drug_charge[drug_spectra_info$drug_charge == "-1"] <- 1
drug_spectra_info$drug_charge[drug_spectra_info$drug_charge == "0"] <- 1
drug_spectra_info$drug_charge[drug_spectra_info$drug_charge == "0+"] <- 1
drug_spectra_info$drug_charge[drug_spectra_info$drug_charge == "1"] <- 1
drug_spectra_info$drug_charge[drug_spectra_info$drug_charge == "1+"] <- 1
drug_spectra_info$drug_charge[drug_spectra_info$drug_charge == "2"] <- 2
drug_spectra_info$drug_charge[drug_spectra_info$drug_charge == "2+"] <- 2
drug_spectra_info$drug_charge[drug_spectra_info$drug_charge == "3"] <- 3
drug_spectra_info$drug_charge <- as.numeric(drug_spectra_info$drug_charge)
drug_analog_spectra <- c(drug_spectra, analog_usi_all)
# writeLines(drug_analog_spectra, "AnalogLib_Development/USI_download/ModifiedCosine_check/input_all_spectra_drug_analog_usi.mgf") 

analog_drug_matches <- new_analog_all_delta_mass_filtered %>% dplyr::select(c("Delta.Mass", "USI", "CCMSLIB"))
analog_drug_matches_all <- analog_drug_matches %>%
  left_join(analog_usi_info, by = c("USI" = "analog_usi")) %>%
  left_join(drug_spectra_info, by = c("CCMSLIB" = "drug_CCMSLIBID"))
analog_drug_matches_all$delta_diff <- ((analog_drug_matches_all$drug_pepmass * analog_drug_matches_all$drug_charge - 1.0078 * analog_drug_matches_all$drug_charge) - (analog_drug_matches_all$analog_pepmass - 1.0078) - analog_drug_matches_all$Delta.Mass)            
analog_drug_matches_all_filtered <- analog_drug_matches_all %>% dplyr::filter(!is.na(analog_pepmass))
match_check<- data.frame(
  analog = analog_drug_matches_all_filtered$USI,
  drug = analog_drug_matches_all_filtered$CCMSLIB
)
# write.csv(match_check, "AnalogLib_Development/USI_download/ModifiedCosine_check/drug_analog_usi_pairs.csv", row.names = F)
# Modified Cosine was checked with matchms package; see “AnalogLib_Development/USI_download/ModifiedCosine_check/mod_cosine_check.py”

# filter drug-analogUSI pairs based on ModifiedCosine results
mod_cosine_check_result <- read.csv("AnalogLib_Development/USI_download/ModifiedCosine_check/mod_cosine_check_output.csv")
pairs_passing_mod_cosine <- mod_cosine_check_result %>% dplyr::filter(score > 0.65) 
colnames(pairs_passing_mod_cosine)[3:4] <- c("matchms_score", "matchms_peak_numbers")

pairs_passing_mod_cosine_all_info <- pairs_passing_mod_cosine %>% left_join(analog_drug_matches_all, by = c("analog" = "USI", "drug" = "CCMSLIB"))
pairs_passing_mod_cosine_all_info$analog_charge[abs(pairs_passing_mod_cosine_all_info$delta_diff) > 0.1] <- 2 # large delta mass difference caused by double charge for analogs
pairs_passing_mod_cosine_all_info <- pairs_passing_mod_cosine_all_info %>%
  mutate(delta_diff_modified = 
           (drug_pepmass * drug_charge - 1.0078 * drug_charge) -
           (analog_pepmass * analog_charge - 1.0078 * analog_charge) -
           Delta.Mass)

#############################################################
# STEP FOUR: REMOVE DRUG ANALOGS THAT MATCH TO GNPS LIBRARY #
#############################################################
# library matches were retrived by comparing the above analog mgf with GNPS2 default library
# job link: https://gnps2.org/status?task=dabe0537e415416c9b2efadd9280b945
gnps_lib_match <- read.csv("AnalogLib_Development/Library_Match_filters/drug_analog_usi_with_gnps_lib_top1.tsv", sep = "\t")
pairs_passing_gnps_lib_check <- pairs_passing_mod_cosine_all_info %>% 
  dplyr::filter(!(analog_scans %in% gnps_lib_match$X.Scan.))

######################################
# STEP FIVE: DEDUPLICATE WITH FALCON #
######################################
usi_to_include <- unique(pairs_passing_gnps_lib_check$analog)
usi_charge_two <- unique(pairs_passing_gnps_lib_check$analog[pairs_passing_gnps_lib_check$analog_charge == 2])

# function to filter mgf files for remaining USI
clean_mgf_usi <- function(input_file, usi_charge_two, usi_to_include) {
  blocks_to_keep <- vector("list", length(input_file) / 10) # Preallocate for speed
  current_block <- character()  # Store the current block as a vector of lines
  block_count <- 0  # Counter for processed blocks
  start_time <- Sys.time()  # Start time tracking
  keep_index <- 1  # Track index for blocks_to_keep
  
  for (line in input_file) {
    current_block <- c(current_block, line)  # Append line to current block
    
    if (startsWith(line, "END IONS")) {
      usi_line <- current_block[grepl("TITLE=", current_block)]  # Find TITLE= line
      usi <- sub("TITLE=", "", usi_line)  # Extract USI
      
      if (usi %in% usi_charge_two) {
        current_block <- sub("CHARGE=1", "CHARGE=2", current_block, fixed = TRUE)  # Modify charge
      }
      
      if (usi %in% usi_to_include) {
        blocks_to_keep[[keep_index]] <- paste0(current_block, collapse = "\n")  # Store processed block
        keep_index <- keep_index + 1
      }
      
      block_count <- block_count + 1  # Increment block count
      
      # Print elapsed time every 100 blocks
      if (block_count %% 1000 == 0) {
        cat(sprintf("Processed %d blocks. Elapsed time: %.2f seconds\n", block_count, as.numeric(difftime(Sys.time(), start_time, units = "secs"))))
      }
      
      current_block <- character()  # Reset block
    }
  }
  
  cleaned_content <- paste0(blocks_to_keep[1:(keep_index - 1)], collapse = "\n\n")  # Combine only non-empty elements
  return(cleaned_content)
}
usi_raw <- readLines("AnalogLib_Development/USI_download/Drug_analog_USI_spectra_RAW.mgf")
usi_modified <- clean_mgf_usi(usi_raw, usi_charge_two, usi_to_include)
# writeLines(usi_modified, "AnalogLib_Development/Falcon/DrugAnalog_USI_passing_delta_modcosine_gnpslib_filters_for_falcon.mgf")
# falcon was run in command line with the output mgf file
# parameters: --precursor_tol 50 ppm --fragment_tol 0.05 --eps 0.1 --min_peaks 3 --min_mz_range 1 --min_mz 50 --max_mz 1500

# filter after falcon
# for clustered USI, choose the representative one as the USI matching to most drug structures
pairs_passing_gnps_lib_check_full <- pairs_passing_gnps_lib_check %>% left_join(druglib_metadata_filtered, by = c("drug" = "gnps_libid"))
usi_to_drug_inchikey_count <- pairs_passing_gnps_lib_check_full %>%
  group_by(analog) %>%
  dplyr::summarise(inchikey_count = n_distinct(split_inchikey))

falcon_results <- read.csv("AnalogLib_Development/Falcon/falcon_cleaned.csv")
falcon_results$analog_scans <- as.numeric(sub("mzspec:USI000000:20250319_analog_USI_passing_delta_modcosine_gnpslib_filters_for_falcon:scan:", "", falcon_results$identifier))
falcon_results <- falcon_results %>%
  left_join(analog_usi_info %>% dplyr::select(c("analog_scans", "analog_usi")), by = "analog_scans") %>%
  left_join(usi_to_drug_inchikey_count, by = c("analog_usi" = "analog"))
falcon_results_usi_keepresult <- falcon_results %>%
  group_by(cluster) %>%
  slice_max(inchikey_count, n = 1, with_ties = FALSE) %>%
  ungroup() # selecting the representative USI for each cluster

pairs_passing_falcon_full <- pairs_passing_gnps_lib_check_full %>% 
  dplyr::filter(analog %in% falcon_results_usi_keepresult$analog_usi) # only keep the representative USI

# creating final names for the analogs
pairs_passing_falcon_full$delta_mass_final <- sprintf("%.2f", 0-pairs_passing_falcon_full$Delta.Mass)
pairs_passing_falcon_full$name_delta_mass <- paste0("'", pairs_passing_falcon_full$name_parent_compound, 
                                                    " (Delta Mass:", 
                                                    pairs_passing_falcon_full$delta_mass_final,
                                                    ")'")

# function to merge all pharmacological information for analogs matching to multiple drugs
combine_unique <- function(x) {
  if (is.numeric(x)) {
    return(as.character(paste(unique(x), collapse = "|")))  # Convert numeric values to character
  } else {
    unique_values <- unique(unlist(strsplit(as.character(x), "\\|")))  # Split, flatten, and remove duplicates
    return(paste(unique_values, collapse = "|"))  # Recombine with "|"
  }
}

# Summarizing Data
unique_analog_metadata <- pairs_passing_falcon_full %>%
  group_by(analog) %>%
  dplyr::summarise(across(everything(), combine_unique), .groups = "drop") 

unique_analog_metadata <- unique_analog_metadata %>%
  mutate(temp = analog) %>% 
  separate(temp, into = c("FILENAME", "EXACTSCAN"), sep = ":scan:", extra = "merge")

##############################################################
# STEP SIX: UPLOAD LIBRARY TO GNPS TO ENABLE DATASET TESTING #
##############################################################

# function to filter mgf files for remaining USI
clean_mgf_usi_2 <- function(input_file, usi_to_include) {
  blocks_to_keep <- vector("list", length(input_file) / 10) # Preallocate for speed
  current_block <- character()  # Store the current block as a vector of lines
  block_count <- 0  # Counter for processed blocks
  start_time <- Sys.time()  # Start time tracking
  keep_index <- 1  # Track index for blocks_to_keep
  
  for (line in input_file) {
    current_block <- c(current_block, line)  # Append line to current block
    
    if (startsWith(line, "END IONS")) {
      usi_line <- current_block[grepl("TITLE=", current_block)]  # Find TITLE= line
      usi <- sub("TITLE=", "", usi_line)  # Extract USI
      
      if (usi %in% usi_to_include) {
        blocks_to_keep[[keep_index]] <- paste0(current_block, collapse = "\n")  # Store processed block
        keep_index <- keep_index + 1
      }
      
      block_count <- block_count + 1  # Increment block count
      
      # Print elapsed time every 100 blocks
      if (block_count %% 1000 == 0) {
        cat(sprintf("Processed %d blocks. Elapsed time: %.2f seconds\n", block_count, as.numeric(difftime(Sys.time(), start_time, units = "secs"))))
      }
      
      current_block <- character()  # Reset block
    }
  }
  
  cleaned_content <- paste0(blocks_to_keep[1:(keep_index - 1)], collapse = "\n\n")  # Combine only non-empty elements
  return(cleaned_content)
}
usi_raw <- readLines("AnalogLib_Development/Falcon/DrugAnalog_USI_passing_delta_modcosine_gnpslib_filters_for_falcon.mgf")
usi_modified <- clean_mgf_usi_2(usi_raw, unique_analog_metadata$analog)
# writeLines(usi_modified, "AnalogLib_Development/Falcon/DrugAnalog_usi_Final_GNPS_upload.mgf")
# gnps library generation job link: https://gnps.ucsd.edu/ProteoSAFe/status.jsp?task=0e2d4e5137c44ab0b2d43d24e71fe4e4

###############################
# STEP SEVEN: DATASET TESTING #
###############################
#removing analog detected in >10% of samples where no drugs were expected

# extract gnps library id of the new analogs
gnps_format_library <- readLines("AnalogLib_Development/Dataset_testing/20230326_updated_drug_analog_before_filtering_GNPS_FORMAT.mgf")
analog_gnps_info <- data.frame(
  analog_scans = as.numeric(sub("SCANS=", "", gnps_format_library[grepl("SCANS=", gnps_format_library)])),
  analog_ID = sub("SPECTRUMID=", "", gnps_format_library[grepl("SPECTRUMID=", gnps_format_library)])
)
analog_gnps_info$analog_scans <- as.character(analog_gnps_info$analog_scans)
analog_metadata <- left_join(analog_gnps_info, unique_analog_metadata, by = c("analog_scans" = "analog_scans"))
analog_metadata$chemical_source <- "Drug_analog"
analog_metadata_selected <- analog_metadata %>%
  dplyr::select(c("analog_ID", "name_delta_mass", "pharmacologic_class", "therapeutic_area", "therapeutic_indication", "mechanism_of_action"))

# function to download library annotation results
download_annotation <- function(task_id) {
  url <- paste0("https://gnps2.org/resultfile?task=", task_id, "&file=nf_output/merged_results_with_gnps.tsv")
  directory_path <- "AnalogLib_Development/Dataset_testing/"
  output_file <- file.path(directory_path, paste0(task_id, ".tsv"))
  response <- GET(url)
  writeBin(content(response, "raw"), output_file)
  ftable <- read_tsv(output_file)
}

# function to calculate number of file occurrences of the analogs in the dataset
cal_occurrence <- function(df_libmatch){
  df_libmatch_summary <- df_libmatch %>% 
    dplyr::filter(MQScore > 0.9 | SharedPeaks > 4) %>%
    group_by(SpectrumID, Compound_Name) %>%
    dplyr::summarise(occurrence = n_distinct(SpectrumFile), .groups = 'drop') %>%
    left_join(analog_metadata_selected, by = c("SpectrumID" = "analog_ID"))
  return(df_libmatch_summary)
}
# function to add analog with high occurrence into exclusion list
add_exclusion <- function(df_libmatch_summary, occurrence_threshold, keep_name, exclude_name){
  df_exclusion <- df_libmatch_summary %>%
    dplyr::filter(occurrence > occurrence_threshold & #exclude analogs with >10% detection 
                    !grepl(keep_name, name_delta_mass)) #keep drugs known to be used in this cohort
  df_exclusion <- analog_metadata_selected %>% 
    dplyr::filter(analog_ID %in% df_exclusion$SpectrumID |
                    grepl(exclude_name, name_delta_mass))
  return(df_exclusion)
}

#GATES_milk_MSV000090877
GATES_milk_MSV000090877 <- download_annotation("7fdeff0051f641e482ac44a7a41e1fd3")
GATES_milk_MSV000090877_summary <- cal_occurrence(GATES_milk_MSV000090877)
exclusion_list <- add_exclusion(df_libmatch_summary = GATES_milk_MSV000090877_summary,
                                occurrence_threshold = 49,
                                keep_name = "nortriptyline|gamithromycin|sulfadimethoxine|framycetin|sulfamerazine",
                                exclude_name = "sclareolide|miltefosine|medroxyprogesterone")
analog_metadata$chemical_source[grepl("nortriptyline|sulfadimethoxine", analog_metadata$name_delta_mass)] <- "Drug_analog|Background/QAQC"

#AGP_MSV000080673
AGP_MSV000080673 <- download_annotation("fff37bcf157d4517a50b850f55dee949")
AGP_MSV000080673_summary <- cal_occurrence(AGP_MSV000080673)
exclusion_list <- add_row(exclusion_list,
                          add_exclusion(df_libmatch_summary = AGP_MSV000080673_summary,
                                occurrence_threshold = 280,
                                keep_name = "sulfathiazole",
                                exclude_name = "icosapent|testosterone|rosiptor"))
analog_metadata$chemical_source[grepl("sulfathiazole", analog_metadata$name_delta_mass)] <- "Drug_analog|Background/QAQC"

# AD mice_MSV000093168 receiving abx
AD_mice_MSV000093168 <- download_annotation("3e5cc81f050d478ba0039faebd0f4bae")
AD_mice_MSV000093168_summary <- cal_occurrence(AD_mice_MSV000093168)
exclusion_list <- add_row(exclusion_list,
                          add_exclusion(df_libmatch_summary = AD_mice_MSV000093168_summary,
                                        occurrence_threshold = 32,
                                        keep_name = "sulfadimethoxine|framycetin",
                                        exclude_name = "2-acetamido-2-deoxy-d-mannopyranose"))
analog_metadata$chemical_source[grepl("sulfadimethoxine", analog_metadata$name_delta_mass)] <- "Drug_analog|Background/QAQC"

# AD mice_MSV000093230 receiving abx
AD_mice_MSV000093230 <- download_annotation("2b6f72bc05c14de8895dee21d2076b8e")
AD_mice_MSV000093230_summary <- cal_occurrence(AD_mice_MSV000093230)
exclusion_list <- add_row(exclusion_list,
                          add_exclusion(df_libmatch_summary = AD_mice_MSV000093230_summary,
                                        occurrence_threshold = 54,
                                        keep_name = "sulfadimethoxine|framycetin",
                                        exclude_name = "mirogabalin"))
analog_metadata$chemical_source[grepl("sulfadimethoxine", analog_metadata$name_delta_mass)] <- "Drug_analog|Background/QAQC"

# AD mice_MSV000096039 receiving abx
AD_mice_MSV000096039 <- download_annotation("a27c4e45c5f241a98b94abfe63cdd12a")
AD_mice_MSV000096039_summary <- cal_occurrence(AD_mice_MSV000096039)
exclusion_list <- add_row(exclusion_list,
                          add_exclusion(df_libmatch_summary = AD_mice_MSV000096039_summary,
                                        occurrence_threshold = 48,
                                        keep_name = "sulfadimethoxine|vancomycin|cefalexin|benzylpenicillin|penciclovir|fusafungine",
                                        exclude_name = "ca-074|oglufanide"))
analog_metadata$chemical_source[grepl("sulfadimethoxine", analog_metadata$name_delta_mass)] <- "Drug_analog|Background/QAQC"

# AD mice_MSV000096038 receiving abx
AD_mice_MSV000096038 <- download_annotation("cdaf98f5c09343d4aae05d3afc00e599")
AD_mice_MSV000096038_summary <- cal_occurrence(AD_mice_MSV000096038)
exclusion_list <- add_row(exclusion_list,
                          add_exclusion(df_libmatch_summary = AD_mice_MSV000096038_summary,
                                        occurrence_threshold = 43,
                                        keep_name = "sulfadimethoxine|vancomycin|cefalexin",
                                        exclude_name = "fsldkjflsdjfoiwejf"))
analog_metadata$chemical_source[grepl("sulfadimethoxine", analog_metadata$name_delta_mass)] <- "Drug_analog|Background/QAQC"


# AD mice_MSV000096037 receiving abx
AD_mice_MSV000096037 <- download_annotation("009a5833cd5145cd8d49e93c33de74ff")
AD_mice_MSV000096037_summary <- cal_occurrence(AD_mice_MSV000096037)
exclusion_list <- add_row(exclusion_list,
                          add_exclusion(df_libmatch_summary = AD_mice_MSV000096037_summary,
                                        occurrence_threshold = 45,
                                        keep_name = "sulfadimethoxine|cefalexin",
                                        exclude_name = "dyphylline|su9516|hydrindantin|ethylestrenol|brazilin"))
analog_metadata$chemical_source[grepl("sulfadimethoxine", analog_metadata$name_delta_mass)] <- "Drug_analog|Background/QAQC"


# AD mice_MSV000096036 receiving abx
AD_mice_MSV000096036 <- download_annotation("d47484bd29af4fad9e80a09d23a055d3")
AD_mice_MSV000096036_summary <- cal_occurrence(AD_mice_MSV000096036)
exclusion_list <- add_row(exclusion_list,
                          add_exclusion(df_libmatch_summary = AD_mice_MSV000096036_summary,
                                        occurrence_threshold = 41,
                                        keep_name = "sulfadimethoxine|sulfamerazine|sulfameter|cefalexin|nortriptyline|sulfamoxole",
                                        exclude_name = "ethyl 3-aminobenzoate"))

# AD mice_MSV000096035 receiving abx
AD_mice_MSV000096035 <- download_annotation("8ee38c757025492daaadbc7341c27dba")
AD_mice_MSV000096035_summary <- cal_occurrence(AD_mice_MSV000096035)
exclusion_list <- add_row(exclusion_list,
                          add_exclusion(df_libmatch_summary = AD_mice_MSV000096035_summary,
                                        occurrence_threshold = 39,
                                        keep_name = "sulfamerazine|nortriptyline",
                                        exclude_name = "jgkgkgjgjgygjd"))

# 6ppd mice_MSV000091363
PPD_mice_MSV000091363 <- download_annotation("294a63bbdeba49189c1c01033e5ac6ca")
PPD_mice_MSV000091363_summary <- cal_occurrence(PPD_mice_MSV000091363)
exclusion_list <- add_row(exclusion_list,
                          add_exclusion(df_libmatch_summary = PPD_mice_MSV000091363_summary,
                                        occurrence_threshold = 28,
                                        keep_name = "nortriptyline",
                                        exclude_name = "rimantadine|hymecromone|progesterone"))

# MPRINT ADMIN milk_MSV000091520
MPRINT_milk_MSV000091520 <- download_annotation("5dbe1ad7216a49c29560e6f118168c58")
MPRINT_milk_MSV000091520_summary <- cal_occurrence(MPRINT_milk_MSV000091520)
exclusion_list <- add_row(exclusion_list,
                          add_exclusion(df_libmatch_summary = MPRINT_milk_MSV000091520_summary,
                                        occurrence_threshold = 300,
                                        keep_name = "sulfadimethoxine|sulfamoxole",
                                        exclude_name = "acarbose"))

# filter the analog library based on the exclusion list
clean_mgf_ccmslib <- function(input_file, libid_to_include) {
  blocks_to_keep <- vector("list", length(input_file) / 10) # Preallocate for speed
  current_block <- character()  # Store the current block as a vector of lines
  block_count <- 0  # Counter for processed blocks
  start_time <- Sys.time()  # Start time tracking
  keep_index <- 1  # Track index for blocks_to_keep
  
  for (line in input_file) {
    current_block <- c(current_block, line)  # Append line to current block
    
    if (startsWith(line, "END IONS")) {
      libid_line <- current_block[grepl("SPECTRUMID=", current_block)]  
      libid <- sub("SPECTRUMID=", "", libid_line)  # Extract USI
      
      if (libid %in% libid_to_include) {
        blocks_to_keep[[keep_index]] <- paste0(current_block, collapse = "\n")  # Store processed block
        keep_index <- keep_index + 1
      }
      
      block_count <- block_count + 1  # Increment block count
      
      # Print elapsed time every 100 blocks
      if (block_count %% 1000 == 0) {
        cat(sprintf("Processed %d blocks. Elapsed time: %.2f seconds\n", block_count, as.numeric(difftime(Sys.time(), start_time, units = "secs"))))
      }
      
      current_block <- character()  # Reset block
    }
  }
  
  cleaned_content <- paste0(blocks_to_keep[1:(keep_index - 1)], collapse = "\n\n")  # Combine only non-empty elements
  return(cleaned_content)
}

# modify the mgf
include_libid <- analog_metadata_selected$analog_ID[!(analog_metadata_selected$analog_ID %in% exclusion_list$analog_ID)]
mgf_raw <- readLines("AnalogLib_Development/Dataset_testing/20230326_updated_drug_analog_before_filtering_GNPS_FORMAT.mgf")
usi_modified <- clean_mgf_ccmslib(mgf_raw,  include_libid)
# writeLines(usi_modified, "20230326_updated_drug_analog_GNPS_FORMAT_first_round_filtered.mgf")

# produce analog metadata
analog_metadata_filtered <- analog_metadata %>% dplyr::filter(!(analog_ID %in% exclusion_list$analog_ID))
analog_metadata_selected <- analog_metadata %>%
  dplyr::select(c("analog_ID", "name_delta_mass", "chemical_source", "pharmacologic_class", "therapeutic_area", "therapeutic_indication", "mechanism_of_action"))
# write.csv(analog_metadata_filtered, "AnalogLib_Development/Dataset_testing/20250327_analog_metadata_filtered.csv")

############################################
# STEP EIGHT: FINAL UPLOAD TO GNPS LIBRARY #
############################################
# combine drug analogs from repository-scale molecular networking
# mgf files: Analog_from_Repo_MN.mgf; metadata: Analog_from_Repo_MN_metadata.csv
# gnps job link to create the final drug analog library: https://gnps.ucsd.edu/ProteoSAFe/status.jsp?task=071202b98c134e3fb88a1f4cfddb5bae
