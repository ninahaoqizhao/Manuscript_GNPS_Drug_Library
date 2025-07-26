# Packages
library(tidyverse)
library(networkD3)
library(plyr)
library(htmlwidgets)
library(webshot)
library(viridis)
library(ComplexHeatmap)
library(circlize)
library(ggplotify)
library(ggpubr)
library(readxl)
library(ggplot2)
library(patchwork)
library(Spectra)
library(msdata)
library(MsBackendMgf)

####################################################################
# Sankey plot for drug library metadata based on number of spectra #
####################################################################

drug_metadata <- read.csv("druglib_metadata/GNPS_Drug_Library_Metadata_Drugs.csv")
drug_metadata[is.na(drug_metadata)] <- "NA"
drug_metadata$therapeutic_area[drug_metadata$therapeutic_area == "NA"] <- "unspecified_area"

# simplify combinatorial exposure sources into single category for the figure
cleanup_source <- drug_metadata$chemical_source
cleanup_source[grepl("Medical", drug_metadata$chemical_source)] <- "Medical"
cleanup_source[grepl("Personal Care", drug_metadata$chemical_source) | grepl("Industrial", drug_metadata$chemical_source)] <- "PersonalCare/Industrial"
cleanup_source[grepl("Food", drug_metadata$chemical_source)] <- "Food"
cleanup_source[grepl("Endogenous", drug_metadata$chemical_source)] <- "Endogenous"
cleanup_source[grepl("metabolite", drug_metadata$chemical_source)] <- "Drug metabolite"
cleanup_source[grepl("QAQC|Low", drug_metadata$chemical_source)] <- "Medical"
drug_metadata <- cbind(drug_metadata, cleanup_source)

# for drugs with multiple therapeutic areas, count them in all areas
metadata_separate <- tidyr::separate_rows(drug_metadata, therapeutic_area, sep = "\\|")
# build count table for the linkages
counts_first <- ddply(metadata_separate, .(metadata_separate$cleanup_source, metadata_separate$therapeutic_area), nrow)
metadata_neurology <- metadata_separate[metadata_separate$therapeutic_area == "neurology/psychiatry",]
metadata_fda_separate <- tidyr::separate_rows(metadata_neurology, pharmacologic_class, sep = "\\|")
counts_second <- ddply(metadata_fda_separate, .(metadata_fda_separate$therapeutic_area, metadata_fda_separate$pharmacologic_class), nrow)
counts_second <- counts_second[counts_second$V1 > 250 &
                                 counts_second$`metadata_fda_separate$pharmacologic_class` != "no_match",] # visualize pharmacologic class with more than 250 spectra
counts_third <- ddply(metadata_fda_separate, .(metadata_fda_separate$pharmacologic_class, metadata_fda_separate$name_compound), nrow)
counts_third <- counts_third[counts_third$V1 > 70 &
                               counts_third$`metadata_fda_separate$pharmacologic_class` %in% counts_second$`metadata_fda_separate$pharmacologic_class` &
                               counts_third$`metadata_fda_separate$pharmacologic_class` != "no_match",] # visualize compounds with more than 70 spectra

names(counts_first) <- c("source", "target", "value")
names(counts_second) <- c("source", "target", "value")
names(counts_third) <- c("source", "target", "value")
counts <- rbind(counts_first, counts_second, counts_third)

# build a connection data frame with a list of flows 
links <- data.frame(source=counts$source, target=counts$target, value=counts$value)
# create a node data frame
nodes <- data.frame(name=c(as.character(links$source), as.character(links$target)) %>% unique())
# reformat the data frame to be compatible with networkD3
links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1
# Sankey plot based on spectral numbers
sankeyNetwork(Links = links, Nodes = nodes,
              Source = "IDsource", Target = "IDtarget",
              Value = "value", NodeID = "name", 
              sinksRight = F, fontSize = 16, fontFamily = "arial", nodePadding = 8)

#############################################################################
# Sankey plot for drug library metadata based on numbers of unique compound #
#############################################################################
drug_metadata <- read.csv("druglib_metadata/GNPS_Drug_Library_Metadata_Drugs.csv")
drug_metadata[is.na(drug_metadata)] <- "NA"
drug_metadata$therapeutic_area[drug_metadata$therapeutic_area == "NA"] <- "unspecified_area"
metadata_w_structure <- drug_metadata[(drug_metadata$inchikey != "NA"),]
custom_aggregate <- function(x) {
  if (is.numeric(x)) {
    return(paste(unique(x), collapse = "|"))
  } else {
    return(paste(unique(x), collapse = "|"))
  }
}
# Group by 'inchikey' and aggregate all columns
metadata_structure_combine <- metadata_w_structure %>%
  group_by(inchikey) %>%
  summarise_all(custom_aggregate)

cleanup_source <- metadata_structure_combine$chemical_source
cleanup_source[grepl("Medical", metadata_structure_combine$chemical_source)] <- "Medical"
cleanup_source[grepl("Personal Care", metadata_structure_combine$chemical_source) | grepl("Industrial", metadata_structure_combine$chemical_source)] <- "PersonalCare/Industrial"
cleanup_source[grepl("Food", metadata_structure_combine$chemical_source)] <- "Food"
cleanup_source[grepl("Endogenous", metadata_structure_combine$chemical_source)] <- "Endogenous"
cleanup_source[grepl("metabolite", metadata_structure_combine$chemical_source)] <- "Drug metabolite"
cleanup_source[grepl("QAQC|Low", metadata_structure_combine$chemical_source)] <- "Medical"
metadata_structure_combine <- cbind(metadata_structure_combine, cleanup_source)

metadata_separate <- tidyr::separate_rows(metadata_structure_combine, therapeutic_area, sep = "\\|")

#choose top 15 therapeutic areas to plot
top_14 <- names(sort(table(metadata_separate$therapeutic_area), decreasing = TRUE)[1:14])
counts_first <- ddply(metadata_separate, .(metadata_separate$cleanup_source, metadata_separate$therapeutic_area), nrow)
names(counts_first) <- c("source", "target", "value")
counts_first$target[!(counts_first$target %in% top_15)] <- "other"

metadata_neurology <- metadata_separate[metadata_separate$therapeutic_area == "neurology/psychiatry",]
metadata_neurology_fda_separate <- tidyr::separate_rows(metadata_neurology, pharmacologic_class, sep = "\\|")
metadata_neurology_fda_separate$pharmacologic_class[metadata_neurology_fda_separate$pharmacologic_class == "antiepileptic drug (AED)"] <- "antiepileptic drug" # simplify names for figure
metadata_neurology_fda_separate$pharmacologic_class[metadata_neurology_fda_separate$pharmacologic_class == "selective serotonin reuptake inhibitor (SSRI)"] <- "SSRI"
metadata_neurology_fda_separate$pharmacologic_class[metadata_neurology_fda_separate$pharmacologic_class == "HMG-CoA reductase inhibitor (statin)"] <- "statins"
metadata_neurology_fda_separate$pharmacologic_class[metadata_neurology_fda_separate$pharmacologic_class == "serotonin and norepinephrine reuptake inhibitor (SNRI)"] <- "SNRI"
metadata_neurology_fda_separate$pharmacologic_class[metadata_neurology_fda_separate$pharmacologic_class == "tricyclic antidepressant (TCA)"] <- "tricyclic antidepressant"

counts_second <- ddply(metadata_neurology_fda_separate, .(metadata_neurology_fda_separate$therapeutic_area, metadata_neurology_fda_separate$pharmacologic_class), nrow)
counts_second <- counts_second[counts_second$V1 > 10 &
                                 counts_second$`metadata_neurology_fda_separate$pharmacologic_class` != "no_match",] # visualize pharmacologic class with more than ten drugs
metadata_ssri <- metadata_neurology_fda_separate[metadata_neurology_fda_separate$pharmacologic_class == "SSRI",]
metadata_ssri$name_parent_compound[grepl("paroxetine", metadata_ssri$name_parent_compound)] <- "paroxetine"
metadata_ssri$name_parent_compound[grepl("fluvoxamine", metadata_ssri$name_parent_compound)] <- "fluvoxamine"
metadata_ssri$name_parent_compound[grepl("citalopram", metadata_ssri$name_parent_compound)] <- "citalopram"
metadata_ssri$name_parent_compound[grepl("sertraline", metadata_ssri$name_parent_compound)] <- "sertraline"
metadata_ssri$name_parent_compound[grepl("trazodone", metadata_ssri$name_parent_compound)] <- "trazodone"
metadata_ssri$name_parent_compound[grepl("fluoxetine", metadata_ssri$name_parent_compound)] <- "fluoxetine"

counts_third <- ddply(metadata_ssri, .(metadata_ssri$pharmacologic_class, metadata_ssri$name_parent_compound), nrow)

names(counts_second) <- c("source", "target", "value")
names(counts_third) <- c("source", "target", "value")
counts <- rbind(counts_first, counts_second, counts_third)

# A connection data frame is a list of flows with intensity for each flow
links <- data.frame(source=counts$source, target=counts$target, value=counts$value)
# create a node data frame
nodes <- data.frame(name=c(as.character(links$source), as.character(links$target)) %>% unique())
# reformat the data frame to be compatible with networkD3
links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1
# Sankey plot based on spectral numbers
sankeyNetwork(Links = links, Nodes = nodes,
              Source = "IDsource", Target = "IDtarget",
              Value = "value", NodeID = "name", 
              sinksRight = F, fontSize = 16, fontFamily = "arial", nodePadding = 8)

###########################################################
# Modifiner to predict modification sites of drug analogs #
###########################################################
analoglib_metadata <- read.csv("druglib_metadata/GNPS_Drug_Library_Metadata_Drug_Analogs_Updated.csv")
drug_metadata <- read.csv("druglib_metadata/GNPS_Drug_Library_Metadata_Drugs.csv")
drug_metadata <- drug_metadata[!duplicated(drug_metadata$gnps_libid), ]
analog_drug_pairs <- analoglib_metadata %>%
  separate_rows(parent_drug_libid, sep = "\\|")
modifinder <- analog_drug_pairs %>% dplyr::select(c("analog_libid", "parent_drug_libid")) %>%
  left_join(drug_metadata %>% dplyr::select(c("gnps_libid", "smiles")), by = c("parent_drug_libid" = "gnps_libid"))
modifinder$analog_libid <- paste0("mzspec:GNPS:GNPS-LIBRARY:accession:", modifinder$analog_libid)
modifinder$parent_drug_libid <- paste0("mzspec:GNPS:GNPS-LIBRARY:accession:", modifinder$parent_drug_libid)
colnames(modifinder) <- c("analog_usi", "drug_usi", "drug_smiles")
# write.csv(modifinder, "20250409_updated_drug_analoglib_modifinder_pairs.csv")
# The .csv file was used for batch ModiFinder processing in gnps2
# job link: https://gnps2.org/status?task=3f23c8ae79c34a6e86e7d64fc6919b76

###########################################################
# Identify potential instrument artifacts in drug analogs #
###########################################################

# STEP ONE: PERFORM REPOSITORY SCALE SEARCH OF DRUG AND ANALOG TO FIND CO-OCCURRENCE IN THE SAME FILE

# Run repository-scale search of the GNPS Drug Library with fastMASST in batch: https://github.com/robinschmid/microbe_masst/tree/v1.2.0
# Output for each spectrum is a tsv files containing all the USI matched to the spectrum
# The resulting files can be downloaded from: 10.5281/zenodo.16423040

# Function to read in each file and filter for zero delta mass (exact match)
process_file <- function(file) {
  df <- read_tsv(file, col_types = cols(.default = "c"))
  if (nrow(df) == 0) {return(NULL)}
  
  df_filtered <- df %>%
    mutate(`Delta Mass` = as.numeric(`Delta Mass`)) %>%
    filter(abs(`Delta Mass`) == 0)
  
  if (nrow(df_filtered) == 0) {return(NULL)}
  
  ccmslib_value <- stringr::str_extract(basename(file), "CCMSLIB\\d+")
  df_filtered <- df_filtered %>% mutate(CCMSLIB = ccmslib_value)
  
  return(df_filtered)
}
# Process all files and combine into one dataframe
drug_masst_file_list <- list.files("/Users/ninazhao/PycharmProjects/masst_only_William/output/20250329_Drug_Exogenous_analog_off", pattern = "_matches.tsv$", full.names = TRUE)
drug_masst <- bind_rows(lapply(drug_masst_file_list, process_file))
drug_masst_filtered <- drug_masst %>% 
  dplyr::filter(!grepl("MSV000084314|MSV000089210", drug_masst$USI) 
                & drug_masst$Cosine >= 0.9) # remove two datasets for repository scale molecular networking, use only cosine > 0.9 match
drug_masst_filtered$short_usi <- sub(":scan.*", "", drug_masst_filtered$USI)
colnames(drug_masst_filtered) <- paste("drug_", colnames(drug_masst_filtered), sep = "")

# masst results for all drug analogs
analog_masst_file_list <- list.files("/Users/ninazhao/PycharmProjects/masst_only_William/output/20250408_Updated_Drug_Analogs_analog_search_off", pattern = "_matches.tsv$", full.names = TRUE)
analog_masst <- bind_rows(lapply(analog_masst_file_list, process_file))
analog_masst_filtered <- analog_masst %>% 
  dplyr::filter(!grepl("MSV000084314|MSV000089210", analog_masst$USI) 
                & analog_masst$Cosine >= 0.9) 
analog_masst_filtered$short_usi <- sub(":scan.*", "", analog_masst_filtered$USI)
colnames(analog_masst_filtered) <- paste("analog_", colnames(analog_masst_filtered), sep = "")

# drug and drug analog metadata
drug_metadata <- read.csv("druglib_metadata/GNPS_Drug_Library_Metadata_Drugs.csv")
analog_metadata <- read.csv("druglib_metadata/GNPS_Drug_Library_Metadata_Drug_Analogs_Updated.csv")
analog_metadata_filtered <- analog_metadata %>% dplyr::select(c("analog_libid", "parent_drug_libid", "chemical_source", "delta_mass")) %>%
  dplyr::filter(!is.na(analog_metadata$parent_drug_libid) &
                  !grepl("Low", analog_metadata$chemical_source))

# collect analog co-occurred with drug
extract_all_libid_of_matched_drugs <- function(row) {
  ccmslib_ids <- unlist(strsplit(row, "\\|"))
  matched_smiles <- drug_metadata %>%
    filter(gnps_libid %in% ccmslib_ids, !is.na(smiles)) %>%
    pull(smiles)
  matched_ccmslibs <- drug_metadata %>%
    filter(smiles %in% matched_smiles | gnps_libid %in% ccmslib_ids) %>%
    pull(gnps_libid)
  unique_ccmslibs <- unique(matched_ccmslibs)
  return(paste(unique_ccmslibs, collapse = "|"))
} # function to find all drugs with the same smile matched to the analog

analog_metadata_filtered <- analog_metadata_filtered %>%
  mutate(matched_drug_all_libid = sapply(parent_drug_libid, extract_all_libid_of_matched_drugs))

# function to find all files where paired drug and drug analogs co-occurred
find_df <- function(row) {
  analog <- row[["analog_libid"]]
  drugs <- str_split(row[["matched_drug_all_libid"]], "\\|")[[1]]
  
  analog_datafile <- analog_masst_filtered %>%
    filter(analog_CCMSLIB == analog) %>%
    pull(analog_short_usi) %>%
    unique()
  
  drug_datafile <- c()
  for (d in drugs) {
    drug_datafile_tmp <- drug_masst_filtered %>%
      filter(drug_CCMSLIB == d) %>%
      pull(drug_short_usi)
    drug_datafile <- union(drug_datafile, drug_datafile_tmp)
  }
  
  common_datafile <- intersect(analog_datafile, drug_datafile)
  
  df2_tmp <- analog_masst_filtered %>%
    filter(analog_CCMSLIB == analog, analog_short_usi %in% common_datafile)
  
  df3_tmp <- drug_masst_filtered %>%
    filter(drug_CCMSLIB %in% drugs, drug_short_usi %in% common_datafile)
  
  df_tmp <- inner_join(df2_tmp, df3_tmp, by = c("analog_short_usi" = "drug_short_usi"), keep = T) %>%
    distinct(analog_USI, drug_USI, .keep_all = TRUE) %>%
    select(analog_CCMSLIB, analog_USI, analog_short_usi,
           drug_CCMSLIB, drug_USI, drug_short_usi)
  
  return(df_tmp)
}

# Apply to all analogs
df_final <- tibble(analog_CCMSLIB = character(), analog_USI = character(), analog_short_usi = character(),
  drug_CCMSLIB = character(), drug_USI = character(), drug_short_usi = character())

for (i in 1:nrow(analog_metadata_filtered)) {
  temp <- find_df(analog_metadata_filtered[i,])
  df_final <- bind_rows(df_final, temp)
  print(c(i, dim(df_final)))
}
colnames(df_final)[4] <- "drug_GNPSID"
length(unique(df_final$analog_GNPSID))/13932
# write.csv(df_final, "ISF_analysis/Drug_analog_co_occurrence_input_for_RT_correlation.csv", row.names = F)

# STEP TWO: PEAK SHAPE CORRELATION FOR DRUG AND ANALOG THAT CO-OCCURRED IN THE SAME FILE
# all raw data were downloaded to local server to run peak shape correlation; see scripts in "ISF_analysis/isf_main.py

correlation_output <- read.csv("ISF_analysis/Drug_analog_RT_correlation_output.csv")
correlation_output_filtered <- correlation_output %>% dplyr::filter(!is.na(ppc)) # raw data download failed for some datasets 
max_R_each_analog <- correlation_output_filtered %>%
  group_by(analog_GNPSID) %>%
  dplyr::summarise(max_R = max(ppc, na.rm = TRUE)) 
max_R_each_analog$ISF <- ifelse(max_R_each_analog$max_R > 0.9, "yes", "no")
sum(max_R_each_analog$ISF == "yes")
# calculate the maximum correlation coefficient observed for each analog-drug pair, in any dataset
# if it is higher than 0.9, it is a potential artifact

# STEP THREE: FOR ANALOGS WITHOUT FILE CO-OCCURRENCE WITH DRUGS, ESTIMATE ARTIFACT POSSIBILITY BASED ON MS/MS FRAGMENTS
extract_all_libid_of_matched_drugs <- function(row) {
  ccmslib_ids <- unlist(strsplit(row, "\\|"))
  matched_smiles <- drug_metadata %>%
    filter(gnps_libid %in% ccmslib_ids, smiles != "N/A") %>%
    pull(smiles)
  matched_ccmslibs <- drug_metadata %>%
    filter(smiles %in% matched_smiles | gnps_libid %in% ccmslib_ids) %>%
    pull(gnps_libid)
  unique_ccmslibs <- unique(matched_ccmslibs)
  return(paste(unique_ccmslibs, collapse = "|"))
}

# negative delta masses
drug_analog <- Spectra("final_mgf/GNPS_Drug_Library_Spectra_Drug_Analogs_Updated.mgf", source = MsBackendMgf())
analog_precusor <- data.frame(analog_id = drug_analog@backend@spectraData@listData[["SPECTRUMID"]], precursor_mz = drug_analog@backend@spectraData@listData[["precursorMz"]])
analog_metadata_filtered_neg <- analog_metadata_filtered %>% left_join(analog_precusor, by = c("analog_libid" = "analog_id")) %>%
  dplyr::filter(grepl("-", analog_metadata_filtered$delta_mass))

GNPS_drug <- Spectra("final_mgf/GNPS_Drug_Library_Spectra_Drugs_and_Metabolites.mgf", source = MsBackendMgf())

# check if the precursor mass of the analog matches any fragment of the drug
analog_intensity_ratio_all <- list()
for (i in 1:nrow(analog_metadata_filtered_neg)){
  analog_id <- analog_metadata_filtered_neg$analog_libid[i]
  analog_mass <- analog_metadata_filtered_neg$precursor_mz[i]
  drug_ids <- unlist(strsplit(analog_metadata_filtered_neg$matched_drug_all_libid[i], "\\|"), use.names = F)
  analog_intensity_ratio <- numeric()
  for (j in 1:length(drug_ids)){
    drug_tag <- which(GNPS_drug@backend@spectraData@listData[["SPECTRUMID"]] == drug_ids[j])
    drug_spectra <- data.frame(mz = GNPS_drug@backend@spectraData@listData[["mz"]]@listData[[drug_tag]],
                               intensity = GNPS_drug@backend@spectraData@listData[["intensity"]]@listData[[drug_tag]])
    drug_spectra_analog <- drug_spectra %>% 
      dplyr::filter((analog_mass - 0.005) < drug_spectra$mz & (analog_mass + 0.005) > drug_spectra$mz)
    ifelse(nrow(drug_spectra_analog) == 0, analog_intensity_ratio[j] <- 0, analog_intensity_ratio[j] <- max(drug_spectra_analog$intensity)/max(drug_spectra$intensity))
  }
  names(analog_intensity_ratio) <- drug_ids
  analog_intensity_ratio_all[[analog_id]] <- analog_intensity_ratio
}

summary_analog_intensity_neg <- data.frame(
  median = sapply(analog_intensity_ratio_all, median),
  max = sapply(analog_intensity_ratio_all, max)
)

# positive delta masses
analog_metadata_filtered_pos <- analog_metadata_filtered %>% left_join(analog_precusor, by = c("analog_libid" = "analog_id")) %>%
  dplyr::filter(!grepl("-", analog_metadata_filtered$delta_mass))

# check if the precursor mass of the drug matches any fragment of the analog
drug_intensity_ratio_all <- list()
for (i in 1:nrow(analog_metadata_filtered_pos)){
  analog_id <- analog_metadata_filtered_pos$analog_libid[i]
  analog_tag <- which(drug_analog@backend@spectraData@listData[["SPECTRUMID"]] == analog_id)
  analog_spectra <- data.frame(mz = drug_analog@backend@spectraData@listData[["mz"]]@listData[[analog_tag]],
                             intensity = drug_analog@backend@spectraData@listData[["intensity"]]@listData[[analog_tag]])
  # analog_mass <- analog_metadata_filtered_pos$precursor_mz[i]
  drug_ids <- unlist(strsplit(analog_metadata_filtered_pos$matched_drug_all_libid[i], "\\|"), use.names = F)
  drug_intensity_ratio <- numeric()
  
  for (j in 1:length(drug_ids)){
    drug_tag <- which(GNPS_drug@backend@spectraData@listData[["SPECTRUMID"]] == drug_ids[j])
    drug_mass <- GNPS_drug@backend@spectraData@listData[["precursorMz"]][drug_tag]
    analog_spectra_drug <- analog_spectra %>% 
      dplyr::filter((drug_mass - 0.005) < analog_spectra$mz & (drug_mass + 0.005) > analog_spectra$mz)
    ifelse(nrow(analog_spectra_drug) == 0, drug_intensity_ratio[j] <- 0, drug_intensity_ratio[j] <- max(analog_spectra_drug$intensity)/max(analog_spectra$intensity))
  }
  names(drug_intensity_ratio) <- drug_ids
  drug_intensity_ratio_all[[analog_id]] <- drug_intensity_ratio
}

summary_analog_intensity_pos <- data.frame(
  median = sapply(drug_intensity_ratio_all, median),
  max = sapply(drug_intensity_ratio_all, max)
)

# if the fragment matches, consider as a potential instrument artifact
summary_analog_intensity_all <- rbind(summary_analog_intensity_neg, summary_analog_intensity_pos)
summary_analog_intensity_all$ISF <- ifelse(summary_analog_intensity_all$max > 0, "yes", "no")
summary_analog_intensity_all <- summary_analog_intensity_all %>% rownames_to_column("analog_GNPSID")
# sum(summary_analog_intensity_all$ISF == "yes")
MS_based_ISF <- summary_analog_intensity_all %>% dplyr::filter(!summary_analog_intensity_all$analog_GNPSID %in% max_R_each_analog$analog_GNPSID)
MS_based_ISF <- MS_based_ISF %>% mutate(source = rep("MS/MS", nrow(MS_based_ISF))) %>% dplyr::select(-c("median", "max"))
ISF_summary <- max_R_each_analog %>% dplyr::select(c("analog_GNPSID", "ISF")) %>% mutate(source = rep("RT", nrow(max_R_each_analog)))
ISF_summary <- rbind(ISF_summary, MS_based_ISF)

# STEP THREE: FOR DELTA MASS THAT CAN ONLY BE EXPLAINED AS ARTIFACTS, MANUALLY ASSIGN AS ARTIFACTS
analog_metadata_filtered_ISF <- analog_metadata_filtered %>% left_join(ISF_summary, by = c("analog_libid" = "analog_GNPSID"))
analog_metadata_filtered_ISF$modified_delta_mass <- paste0("|", analog_metadata_filtered_ISF$delta_mass, "|")
# +1 as C isotope
analog_metadata_filtered_ISF$ISF[grepl("\\|1\\|", analog_metadata_filtered_ISF$modified_delta_mass)] <- "yes"
analog_metadata_filtered_ISF$source[grepl("\\|1\\|", analog_metadata_filtered_ISF$modified_delta_mass)] <- "isotope"
# +2 as S/Cl/Br/O isotope
analog_metadata_filtered_ISF$ISF[grepl("\\|2\\|", analog_metadata_filtered_ISF$modified_delta_mass)] <- "yes"
analog_metadata_filtered_ISF$source[grepl("\\|2\\|", analog_metadata_filtered_ISF$modified_delta_mass)] <- "isotope"
# +21.98 as Na adduct
analog_metadata_filtered_ISF$ISF[grepl("\\|21.98\\|", analog_metadata_filtered_ISF$modified_delta_mass)] <- "yes"
analog_metadata_filtered_ISF$source[grepl("\\|21.98\\|", analog_metadata_filtered_ISF$modified_delta_mass)] <- "adduct"
# +37.95 as Ca adduct
analog_metadata_filtered_ISF$ISF[grepl("\\|37.95\\|", analog_metadata_filtered_ISF$modified_delta_mass)] <- "yes"
analog_metadata_filtered_ISF$source[grepl("\\|37.95\\|", analog_metadata_filtered_ISF$modified_delta_mass)] <- "adduct"
# +37.96 as K adduct
analog_metadata_filtered_ISF$ISF[grepl("\\|37.96\\|", analog_metadata_filtered_ISF$modified_delta_mass)] <- "yes"
analog_metadata_filtered_ISF$source[grepl("\\|37.96\\|", analog_metadata_filtered_ISF$modified_delta_mass)] <- "adduct"

# STEP FOUR: CONSTRUCT SANKEY PLOT
data <- analog_metadata_filtered_ISF %>% 
  dplyr::select(c("analog_libid", "delta_mass", "ISF", "source"))
data$adduct <- data$source
data$adduct[data$source != "isotope" & data$source != "adduct"] <- "others"
data$co_occur <- data$source
data$judgement <- data$ISF
data$judgement[data$source == "RT" & data$ISF == "yes"] <- "R-square > 0.9"
data$judgement[data$source == "RT" & data$ISF == "no"] <- "R-square < 0.9"
data$judgement[data$source == "MS/MS" & data$ISF == "yes"] <- "analog/drug precursor as drug/analog fragments"
data$judgement[data$source == "MS/MS" & data$ISF == "no"] <- "no overlaps in precursor/fragments of drugs/analogs"
data$judgement[data$source == "isotope"] <- "isotope"
data$judgement[data$source == "adduct"] <- "adduct"
data$final_class <- data$ISF
data$final_class[data$ISF == "yes" & grepl("-", data$delta_mass)] <- "Possible In-source Fragments"
data$final_class[data$ISF == "yes" & !grepl("-", data$delta_mass)] <- "Possible Adducts"
data$final_class[data$ISF == "no"] <- "Non-instrumental artifacts"
data$final_class[data$source == "isotope"] <- "isotope"
data_final <- data %>% dplyr::select(-c("delta_mass", "source"))
data_final$analog <- "analog"

# Aggregate the data to calculate frequencies
aggregated_data <- data_final %>%
  group_by(analog, adduct, co_occur, judgement, final_class, ISF) %>%
  dplyr::summarise(freq = n(), .groups = "drop")
head(aggregated_data)
aggregated_data$judgement <- sub("isotope", "isotope_2", aggregated_data$judgement)
aggregated_data$judgement <- sub("adduct", "adduct_2", aggregated_data$judgement)
aggregated_data$co_occur <- sub("isotope", "isotope_1", aggregated_data$co_occur)
aggregated_data$co_occur <- sub("adduct", "adduct_1", aggregated_data$co_occur)
aggregated_data$final_class <- sub("isotope", "isotope_3", aggregated_data$final_class)

# Prepare nodes
nodes <- data.frame(name = unique(c(aggregated_data$analog, aggregated_data$adduct, aggregated_data$co_occur,
  aggregated_data$judgement, aggregated_data$final_class)))

# Prepare links
links <- data.frame(
  source = match(aggregated_data$analog, nodes$name) - 1,  # Analog to Adduct
  target = match(aggregated_data$adduct, nodes$name) - 1,
  value = aggregated_data$freq,
  group = aggregated_data$ISF  # Color by ISF
)
links <- rbind(
  links,
  data.frame(
    source = match(aggregated_data$adduct, nodes$name) - 1,  # Adduct to Co-occur
    target = match(aggregated_data$co_occur, nodes$name) - 1,
    value = aggregated_data$freq,
    group = aggregated_data$ISF
  ),
  data.frame(
    source = match(aggregated_data$co_occur, nodes$name) - 1,  # Co-occur to Judgement
    target = match(aggregated_data$judgement, nodes$name) - 1,
    value = aggregated_data$freq,
    group = aggregated_data$ISF
  ),
  data.frame(
    source = match(aggregated_data$judgement, nodes$name) - 1,  # Judgement to Final Class
    target = match(aggregated_data$final_class, nodes$name) - 1,
    value = aggregated_data$freq,
    group = aggregated_data$ISF
  )
)
link_color <- 'd3.scaleOrdinal()
                .domain(["yes", "no"])
                .range(["rgba(249, 200, 34)", "rgba(135, 37, 89)"])' 
node_color <- 'd3.scaleOrdinal()
                .domain(["white"])
                .range(["#ffffff"])' 

# Set NodeGroup to fill nodes with white
nodes$group <- "white"  # Add a dummy group for white color fill
sankey <- sankeyNetwork(
  Links = links, Nodes = nodes, Source = "source", Target = "target", Value = "value", NodeID = "name",
  LinkGroup = "group", NodeGroup = "group", fontSize = 12, nodeWidth = 30, nodePadding = 40, colourScale = link_color)
# Save the Sankey plot as an HTML file
# saveWidget(sankey, "sankey_plot.html", selfcontained = TRUE)

######################################
# delta mass frequency visualization #
######################################
analog_metadata <- read.csv("druglib_metadata/GNPS_Drug_Library_Metadata_Drug_Analogs_Updated.csv")
analog_metadata_separate <- analog_metadata %>% dplyr::select(c("analog_libid", "delta_mass")) %>%
  separate_rows(delta_mass, sep = "\\|") 
delta_count <- as.data.frame(table(analog_metadata_separate$delta_mass))
delta_count$atom <- rep("others",nrow(delta_count))
colnames(delta_count) <- c("delta_mass", "number_of_drug_analog", "atom")
delta_count$delta_mass <- as.numeric(as.character(delta_count$delta_mass))

# add atomic composition for frequent delta masses
delta_count$atom[delta_count$delta_mass == 14.02 | delta_count$delta_mass == -14.02] <- "1C, 2H" 
delta_count$atom[delta_count$delta_mass == 17.03 | delta_count$delta_mass == -17.03] <- "1N, 3H" 
delta_count$atom[delta_count$delta_mass == -15.99 | delta_count$delta_mass == 15.99] <- "1O" 
delta_count$atom[delta_count$delta_mass == 18.01 | delta_count$delta_mass == -18.01] <- "1O, 2H" 
delta_count$atom[delta_count$delta_mass == 28.03 | delta_count$delta_mass == -28.03] <- "2C, 4H" 
delta_count$atom[delta_count$delta_mass == 30.01 | delta_count$delta_mass == -30.01] <- "1C, 2H, 1O" 

delta_count$atom <- factor(delta_count$atom, levels = c("1C, 2H", "1N, 3H", "1O", "1O, 2H", "2C, 4H", "1C, 2H, 1O", "others"))

barplot <- ggplot(delta_count, aes(x=delta_mass, weight=number_of_drug_analog, fill = atom)) +
  geom_histogram(binwidth = 1, boundary = 0, color = "white", size = 0) +
  scale_fill_manual(values = c("darkred", "#339933", "#0066CC", "#66cccc", "orange", "pink", "black")) + 
  xlab("Delta mass") +
  ylab("Spectra number") +
  scale_x_continuous(limits = c(-100, 100), breaks = c(-100, -50, 0, 50, 100))+
  scale_y_continuous(limits = c(0, 2300), breaks = c(0, 500, 1000, 1500, 2000), expand = c(0, 0))+
  theme_classic() +
  theme(axis.ticks.length=unit(.25, "cm"),
        axis.line.x.bottom=element_line(size=1.2), axis.line.y.left = element_line(size=1.2),
        axis.ticks.x.bottom = element_line(size=1.2), axis.ticks.y.left = element_line(size=1.2),
        axis.text.y.left = element_text(size = 20), axis.text.x.bottom = element_text(size = 20),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.position = c(1, 0.7),  # Align legend to the right edge of the plot area
        legend.justification = c(1, 0.5),
        legend.title = element_blank(),  # Remove legend title
        plot.margin = margin(5.5, 5.5, 5.5, 5.5))
# ggsave(plot = barplot, file = "analog_delta_mass_barplot.svg", width = 10, height = 6)

#######################################################
# drug detection in different disease classes in ReDU #
#######################################################

# read in GNPS Drug Library masst matches
drug_masst <- read.csv("masst_results_entire_library/druglib_parent_drug_fasst_matches.csv")
drug_masst <- drug_masst[,-1]
analog_masst <- read.csv("masst_results_entire_library/druglib_drug_analogs_fasst_matches.csv")
all_masst <- rbind(drug_masst, analog_masst)

all_masst_filtered <- all_masst %>% 
  dplyr::filter(!grepl("MSV000084314|MSV000089210", all_masst$USI)
                & all_masst$Cosine >= 0.8
                & all_masst$Matching.Peaks >= 6
                & all_masst$Delta.Mass == 0.00) %>%
  dplyr::select(c("USI", "GNPSID")) # filter detections for cosine higher than 0.8 and match peaks more than 6
all_masst_filtered$short_usi <- sub(":scan.*", "", all_masst_filtered$USI)
all_masst_filtered <- all_masst_filtered[,c("short_usi", "GNPSID")]
colnames(all_masst_filtered) <- c("short_usi", "drug_GNPSID")

# read in ReDU metadata on disease and tissue types
redu <- read.csv("masst_results_entire_library/all_redu_sample_information.tsv", sep = "\t")
redu_gnps <- redu %>% dplyr::filter(redu$DataSource == "GNPS")
extract_usi_redu <- function(x) {
  part1 <- str_extract(x, ".*:")
  part2 <- str_extract(x, "(?<=/)[^/]+(?=\\.)")
  result <- paste0(part1, part2)
  return(result)
}
redu_gnps_short_usi <- redu_gnps %>% mutate(short_usi = sapply(USI, extract_usi_redu))
redu_human_disease <- redu_gnps_short_usi %>% dplyr::filter(grepl("Homo", redu_gnps_short_usi$NCBITaxonomy) & redu_gnps_short_usi$DOIDCommonName != "missing value")
redu_human_disease_filtered <- redu_human_disease %>% dplyr::select(c("short_usi", "DOIDCommonName", "UBERONBodyPartName"))
redu_human_disease_filtered$disease_bodypart <- paste0(redu_human_disease_filtered$DOIDCommonName, " | ", redu_human_disease_filtered$UBERONBodyPartName) 

# join drug masst detection with redu metadata
all_masst_disease_redu <- all_masst_filtered %>% inner_join(redu_human_disease_filtered, by = "short_usi")

# read in drug and drug analog metadata
drug_metadata <- read.csv("druglib_metadata/GNPS_Drug_Library_drug_metadata.csv")
analog_metadata <- read.csv("druglib_metadata/GNPS_Drug_Library_drug_analog_metadata.csv")
analog_id_match <- read.csv("masst_results_entire_library/drug_analoglib_ID_match.csv") 
analog_metadata <- analog_metadata %>% left_join(analog_id_match, by = c("analog_libid" = "new_id"))
drug_metadata_nodup <- drug_metadata[!duplicated(drug_metadata$gnps_libid),]
drug_metadata_selected <- drug_metadata_nodup %>% dplyr::select(c("gnps_libid", "name_compound", "chemical_source", "pharmacologic_class",
                                                                  "therapeutic_area", "therapeutic_indication", "mechanism_of_action"))
analog_metadata_selected <- analog_metadata %>% dplyr::select(c("analog_libid", "name_analog", "chemical_source", "pharmacologic_class",
                                                                "therapeutic_area", "therapeutic_indication", "mechanism_of_action"))
colnames(analog_metadata_selected) <- colnames(drug_metadata_selected)
drug_metadata_selected <- rbind(drug_metadata_selected, analog_metadata_selected)

# join drug masst detection with drug metadata
drug_human_disease_metadata <- all_masst_disease_redu %>% left_join(drug_metadata_selected, by = c("drug_GNPSID" = "gnps_libid"))
drug_human_disease_metadata_filtered <- drug_human_disease_metadata %>%
  dplyr::filter(!grepl("QAQC|Low|Endoge|Food|Personal|Industrial", drug_human_disease_metadata$chemical_source)
                & !grepl("NA|no_match", drug_human_disease_metadata$therapeutic_area)
                & !is.na(drug_human_disease_metadata$therapeutic_area))

# separate antibiotics, antifungal, antivirals, simplify the disease area for visulization
drug_human_disease_metadata_filtered$therapeutic_area[grepl("HIV|influenza|nucleoside", drug_human_disease_metadata_filtered$pharmacologic_class)] <- "antivirals"
drug_human_disease_metadata_filtered$therapeutic_area[grepl("antifungal", drug_human_disease_metadata_filtered$pharmacologic_class)] <- "antifungals"
drug_human_disease_metadata_filtered$therapeutic_area[grepl("cardiology", drug_human_disease_metadata_filtered$therapeutic_area)] <- "cardiology"
drug_human_disease_metadata_filtered$therapeutic_area[grepl("antimalarial", drug_human_disease_metadata_filtered$pharmacologic_class)] <- "anti-parasites|cardiology" # this is quinine and quinidine
drug_human_disease_metadata_filtered$therapeutic_area[drug_human_disease_metadata_filtered$therapeutic_area == "dental|infectious disease"] <- "antiseptic and disinfectant"
drug_human_disease_metadata_filtered$therapeutic_area[grepl("sulfone|antimicrobial|antibac|tetracycline", drug_human_disease_metadata_filtered$pharmacologic_class)] <- "antibiotics"
drug_human_disease_metadata_filtered$therapeutic_area[grepl("parasite|leishmaniasis", drug_human_disease_metadata_filtered$therapeutic_indication)] <- "anti-parasites"
drug_human_disease_metadata_filtered$therapeutic_area[grepl("corticosteroid", drug_human_disease_metadata_filtered$pharmacologic_class)] <- "corticosteroids"

# get unique therapeutic area to count detections
all_area <- unique(drug_human_disease_metadata_filtered$therapeutic_area)
unique_disease_area <- strsplit(all_area, "\\|")
unique_disease_area <- unlist(unique_disease_area)
unique_disease_area <- unique(unique_disease_area)

# aggregate all therapeutic areas detected in each sample 
drug_area_separate <- tidyr::separate_rows(drug_human_disease_metadata_filtered, therapeutic_area, sep = "\\|")
custom_aggregate <- function(x) {
  return(paste(unique(x), collapse = "|"))
}
drug_usi_aggregate <- drug_area_separate %>%
  group_by(short_usi) %>%
  summarise_all(custom_aggregate)

# summarize drug detection based on theraputic area in each sample
result <- data.frame(matrix(FALSE, nrow = nrow(drug_usi_aggregate), ncol = length(unique_disease_area)))
colnames(result) <- unique_disease_area
rownames(result) <- drug_usi_aggregate$short_usi

drug_usi_aggregate$therapeutic_area_fixed <- paste0("|", drug_usi_aggregate$therapeutic_area, "|")
for(i in 1:nrow(result)){
  for (j in 1:ncol(result)){
    if(grepl(paste0("|", unique_disease_area[j], "|"), drug_usi_aggregate$therapeutic_area_fixed[drug_usi_aggregate$short_usi == rownames(result)[i]], fixed = TRUE)){
      result[i,j] <- TRUE
    }
  }
}

# calculate detection frequency of drugs in each disease type
result_filtered <- result %>% rownames_to_column("short_usi") %>% left_join(drug_usi_aggregate %>% dplyr::select(c("short_usi", "disease_bodypart")), by = "short_usi")
result_freq_per_area <- result_filtered %>%
  dplyr::select(-"short_usi") %>%
  group_by(disease_bodypart) %>%
  dplyr::summarise(across(everything(), ~ sum(. == TRUE) / n(), .names = "ratio_{col}")) %>%
  column_to_rownames("disease_bodypart")

# create heatmap visualization
heatmap_data <- as.matrix(result_freq_per_area[c("inflammatory bowel disease | feces", "Kawasaki disease | blood plasma", "dental caries | saliva", "COVID-19 | blood plasma",
                                                 "psoriasis | arm skin", "psoriasis | head or neck skin", "psoriasis | skin of trunk",
                                                 "acquired immunodeficiency syndrome | blood plasma", "human immunodeficiency virus infectious disease | blood plasma",
                                                 "Alzheimer's disease | blood plasma", "Alzheimer's disease | cerebrospinal fluid", "Crohn's disease | feces", 
                                                 "diabetes mellitus | feces", "diabetes mellitus | urine", "hypertension | urine", "obesity | saliva", "primary bacterial infectious disease | feces",
                                                 "primary bacterial infectious disease | nasal cavity", "primary bacterial infectious disease | saliva", "primary bacterial infectious disease | skin of trunk",
                                                 "ulcerative colitis | feces"),
                                               c("ratio_antibiotics", "ratio_antifungals", "ratio_antivirals", "ratio_cardiology", 
                                                 "ratio_neurology/psychiatry", "ratio_gastroenterology", "ratio_corticosteroids",
                                                 "ratio_rheumatology", "ratio_dermatology", "ratio_pulmonary", "ratio_antiseptic and disinfectant"
                                                 )]) # order the rows and columns
col_fun <- colorRamp2(c(0, 0.5, 1), c("#336699", "#FFFF99", "#FF3333"))
masst_disease_heatmap <- 
  Heatmap(heatmap_data, name = "Detection_Frequency", col = col_fun, 
          cluster_rows = F, cluster_columns = F, rect_gp = gpar(col = "lightgray", lwd = 1), row_names_side = "right")
# heatmap_ggplot <- as.ggplot(grid.grabExpr(draw(masst_disease_heatmap)))
# ggsave(file="drug_disease_area_masst_heatmap.svg", plot=heatmap_ggplot, width=10, height=8)

############################################################################
# HNRC cohort analysis - drug and drug metabolite/analog detection heatmaps#
############################################################################
# Due to human subject protection constraint, metadata for the HNRC cohort will be provided upon request to HNRC: https://hnrp.hivresearch.ucsd.edu.

# read feature table and annotations
feature <- read.csv("HNRC_analysis/HDRC_quant.csv")
feature_filter <- feature %>% column_to_rownames("row.ID") %>% dplyr::select(contains("Sample"))
annotation <- read.csv("HNRC_analysis/HNRC_DrugLib_annotations.tsv", sep = "\t") # https://gnps2.org/status?task=0cf8b2ea91ce4cc4bdd10b47cbcd5cd1
annotation_filtered <- annotation %>% dplyr::filter(MQScore > 0.9 | SharedPeaks > 4) %>% 
  dplyr::select(c("X.Scan.", "SpectrumID", "Compound_Name")) 
colnames(annotation_filtered) <- c("FeatureID", "SpectrumID", "Compound_Name")

# read in Drug Library metadata
druglib_metadata <- read.csv("druglib_metadata/GNPS_Drug_Library_Metadata_Drugs.csv")
druglib_metadata_nodup <- druglib_metadata[!duplicated(druglib_metadata$gnps_libid),]
analoglib_metadata <- read.csv("druglib_metadata/GNPS_Drug_Library_Metadata_Drug_Analogs_Updated.csv")
druglib_metadata_filtered <- druglib_metadata_nodup %>% 
  dplyr::select(c("gnps_libid", "name_compound", "name_parent_compound", "chemical_source", "pharmacologic_class", "therapeutic_area", "therapeutic_indication"))
analoglib_metadata_filtered <- analoglib_metadata %>% 
  dplyr::select(c("analog_libid","name_analog","name_connected_compound", "chemical_source", "pharmacologic_class", "therapeutic_area", "therapeutic_indication"))
colnames(analoglib_metadata_filtered) <- colnames(druglib_metadata_filtered)
all_metadata <- rbind(druglib_metadata_filtered, analoglib_metadata_filtered)
annotation_metadata <-  merge(annotation_filtered, all_metadata, by.x = "SpectrumID", by.y = "gnps_libid", all.x = T, all.y=F)
exo_drug <- annotation_metadata %>% dplyr::filter(!grepl("Background|confidence|Endogenous|Food", annotation_metadata$chemical_source))
# write.csv(exo_drug,"HNRC_analysis/HNRC_exogenous_drug_for_manual_name_cleaning.csv", row.names = F)
# The exported table was manually cleaned up by:
# (a) Clean up the "name_parent_compound" column to consistent names for parent drugs;
# (b) Assign annotation type for each unique drug as "drug only", "drug with metabolites/analogs", "metabolites only", or "analog only"
# (c) Reject annotation for feature 13703 - this is nevirapine annotation that were rejected by standard check
# (d) For drugs with metabolites/analogs, their pharmacologic classes were manually cleaned up for better visualization
# (e) Results is provided in "HNRC_exogenous_drug_manually_cleaned.csv"

# read in the annotation table after manual cleanup
exo_drug_cleaned <- read.csv("HNRC_analysis/HNRC_exogenous_drug_manually_cleaned.csv")
unique_combos <- exo_drug_cleaned %>%
  distinct(cleanup_name_drug_only, Group)
table(unique_combos$Group) # 66 out of 175 drugs are detected with metabolites/analogs

exo_drug_with_analog <- exo_drug_cleaned %>% dplyr::filter(Group == "both") %>% arrange(Final_class)
summary_analog_number <- exo_drug_with_analog %>% dplyr::select(c("cleanup_name_drug_only", "chemical_source")) %>%
  group_by(cleanup_name_drug_only, chemical_source) %>%
  dplyr::summarise(count = n(), .groups = 'drop') %>%
  pivot_wider(names_from = chemical_source, values_from = count, values_fill = 0)
drugs_to_plot <- summary_analog_number$cleanup_name_drug_only[summary_analog_number$Medical > 0 & summary_analog_number$Drug_analog > 1] #plot drug with at least two analogs

exo_drug_with_analog$FeatureID <- as.character(exo_drug_with_analog$FeatureID)
feature_annotation <- exo_drug_with_analog %>% inner_join((feature_filter %>% rownames_to_column("FeatureID")), by = "FeatureID")
feature_annotation_filtered <- feature_annotation %>% column_to_rownames("FeatureID") %>% dplyr::select(contains(c(".mzML", "cleanup_")))
feature_annotation_intensity_filtered <- feature_annotation_filtered %>% dplyr::select(-cleanup_name_drug_only)

# filter for drugs detected in >10% of samples for better visualization
feature_annotation_intensity_filtered[feature_annotation_intensity_filtered < 10000] <- 0
feature_annotation_intensity_filtered <- feature_annotation_intensity_filtered %>% rownames_to_column("FeatureID") %>% 
  left_join(feature_annotation_filtered %>% rownames_to_column("FeatureID") %>% dplyr::select(c("FeatureID", "cleanup_name_drug_only")), by = "FeatureID") %>%
  column_to_rownames("FeatureID")
exo_drug_group <- feature_annotation_intensity_filtered %>% aggregate(. ~ cleanup_name_drug_only, sum) %>% column_to_rownames("cleanup_name_drug_only")
drug_keep <- intersect(drugs_to_plot, rownames(exo_drug_group)[rowSums(exo_drug_group>0)>32]) # plot drugs detected in >10% of samples and with at least two analogs
drug_keep <- c(drug_keep, "terbinafine") # terbinafine was added in microbial cultures, so included in Figure 2a

# produce the heatmap for drugs and analogs
exo_drug_heatmap <- feature_annotation %>% dplyr::filter(feature_annotation$cleanup_name_drug_only %in% drug_keep) 
exo_drug_heatmap_clean <- exo_drug_heatmap %>% dplyr::select(contains(".mzML"))
exo_drug_heatmap_filtered <- exo_drug_heatmap_clean[,colSums(exo_drug_heatmap_clean)>0]
exo_drug_heatmap_log <- log10(exo_drug_heatmap_filtered+1)
exo_drug_heatmap_log <- as.matrix(exo_drug_heatmap_log)

col_fun <- colorRamp2(c(0, 4, 5, 8.5), c("#336699", "#336699", "#FFFF99", "#FF3333"))
exo_drug_heatmap <- Heatmap(exo_drug_heatmap_log, name = "Log_Peak_Area", col = col_fun, cluster_rows = F, show_column_dend = F, show_column_names = F, show_row_names = F, 
                            row_split = factor(exo_drug_heatmap$cleanup_name_drug_only, levels = unique(exo_drug_heatmap$cleanup_name_drug_only)), 
                            row_gap = unit(0.3, "mm"), cluster_row_slices = F,
                            row_title_rot = 0, row_title_gp = gpar(fontsize = 24),
                            row_title_side = "right",
                            heatmap_legend_param = list(title_gp = gpar(fontsize = 24), labels_gp = gpar(fontsize = 24))
) +
  Heatmap(exo_drug_heatmap$chemical_source, name = "annotation", width = unit(3, "mm"),
          col = c("Medical" = "#c25759", "Drug metabolite" = "#edb8b0", "Drug_analog" = "#F5DFDB"),
          heatmap_legend_param = list(title_gp = gpar(fontsize = 24), labels_gp = gpar(fontsize = 24))
  ) +
  Heatmap(exo_drug_heatmap$Final_class, name = "class", width = unit(3, "mm"),
          col = c("antibiotics" = "#990033", "antidepressant" = "#CC6633", "antifungal" = "#CD9600", "antihistamine" = "#66CC99",
                  "antimalaria" = "#336666", "antiplatelet drug" = "#333333", "antiseptic" = "#339933", "anxiolytics" = "#99CC99", "atypical antipsychotic" = "#B5E2FF", 
                  "cardiology medication" = "#3399CC", "HIV medication" = "#330066", "insomnia medication" = "#9999CC", 
                  "muscle relaxant" = "#999999", "proton pump inhibitor" = "#666600", "statin" = "#CCCC00"),
          heatmap_legend_param = list(title_gp = gpar(fontsize = 24), labels_gp = gpar(fontsize = 24))
  )
heatmap_ggplot <- as.ggplot(grid.grabExpr(draw(exo_drug_heatmap)))
# ggsave(file="exo_drug_heatmap.png", plot=heatmap_ggplot, width=15.2, height=11, dpi = 800)

# heatmap darunavir
darunavir <- feature_annotation %>% dplyr::filter(cleanup_name_drug_only == "darunavir")
darunavir$final_name <- paste(darunavir$cleanup_name, darunavir$FeatureID, sep = " | ")
darunavir_filtered <- darunavir %>% column_to_rownames("final_name") %>% dplyr::select(contains(c(".mzML")))
darunavir_filtered[darunavir_filtered < 10000] <- 0
darunavir_heatmap <- darunavir_filtered[,colSums(darunavir_filtered)>0]
darunavir_heatmap_log <- log10(darunavir_heatmap+1)
darunavir_heatmap_log <- as.matrix(darunavir_heatmap_log)

col_fun <- colorRamp2(c(0, 4, 5, 8.5), c("#336699", "#336699", "#FFFF99", "#FF3333"))
darunavir_heatmap <- Heatmap(darunavir_heatmap_log, name = "Log_Peak_Area", col = col_fun, 
                             cluster_rows = F, show_column_dend = F, show_column_names = F,
                             row_names_gp = gpar(fontsize = 14),
                             row_gap = unit(5, "mm"), cluster_row_slices = F,
                             show_heatmap_legend = F
)
heatmap_ggplot <- as.ggplot(grid.grabExpr(draw(darunavir_heatmap)))
# ggsave(file="darunavir_heatmap.svg", plot=heatmap_ggplot, width=10, height=4)

# check co-occurrence between darunavir analog and the parent drug
darunavir_heatmap_t <- as.data.frame(t(darunavir_heatmap))
darunavir_heatmap_t$parent_detection <- darunavir_heatmap_t$`darunavir | 20602` > 0
detection_w_parent <- colSums((darunavir_heatmap_t[darunavir_heatmap_t$parent_detection == TRUE,] > 0)) / colSums(darunavir_heatmap_t >0)
min(detection_w_parent)
max(detection_w_parent)
median(detection_w_parent)

##############################################
# people with HIV drug cluster analysis #
##############################################

# HIV drug annotations combined with peak area
exo_drug_cleaned <- read.csv("HNRC_analysis/HNRC_exogenous_drug_manually_cleaned.csv")
feature <- read.csv("HNRC_analysis/HDRC_quant.csv")
feature_filter <- feature %>% column_to_rownames("row.ID") %>% dplyr::select(contains("Sample"))
exo_drug_cleaned$FeatureID <- as.character(exo_drug_cleaned$FeatureID)
hiv_drug_detection <- exo_drug_cleaned %>% inner_join((feature_filter %>% rownames_to_column("FeatureID")), by = "FeatureID") %>%
  dplyr::filter((grepl("HIV", therapeutic_indication) | grepl("atazanavir", therapeutic_indication)) & !grepl("analog_only", Group))
hiv_drug_detection_filtered <- hiv_drug_detection %>% dplyr::select(contains(c("FeatureID", ".mzML"))) %>% column_to_rownames("FeatureID")

hiv_drug_detection_noise_filtered <- hiv_drug_detection_filtered
hiv_drug_detection_noise_filtered[hiv_drug_detection_noise_filtered < 50000] <- 0
hiv_drug_detection$FeatureID <- as.character(hiv_drug_detection$FeatureID)
hiv_drug_detection_noise_filtered_final <- hiv_drug_detection_noise_filtered %>% rownames_to_column("FeatureID") %>%
  left_join(hiv_drug_detection %>% dplyr::select(contains(c("Feature", "cleanup"))), by = "FeatureID") %>%
  dplyr::select(-FeatureID)

# sum peak areas of drug, metabolites, and analogs
hiv_drug_group <- hiv_drug_detection_noise_filtered_final %>%
  aggregate(. ~ cleanup_name_drug_only, sum) %>% column_to_rownames("cleanup_name_drug_only")

# join HIV drug peak area with participant metadata
process_string <- function(x) {
  # Extract the substring starting from "Sample_" and before ".mzML"
  extracted <- sub(".*(Sample_.*)\\.mzML.*", "\\1", x)
  # Replace "_" with "." after "Sample_"
  modified <- sub("Sample_", "Sample.", gsub("_", ".", extracted))
  return(modified)
}
colnames(hiv_drug_group) <- sapply(colnames(hiv_drug_group), process_string)
hiv_drug_log <- log10(hiv_drug_group+1)
metadata <- read_xlsx("HNRC_analysis/hnrp_mibi_master_metadata_072424_updated.xlsx")  #metadata for the HNRC cohort will be provided upon request to HNRC: https://hnrp.hivresearch.ucsd.edu
metadata$SampleID <- sub("11135.", "", metadata$SampleID)
metadata_filtered <- metadata %>% dplyr::select(c("SampleID", "hiv_status_clean", "arv_status"))
hiv_drug_log_filtered <- hiv_drug_log %>% dplyr::select(any_of(metadata_filtered$SampleID))
metadata_filtered <- metadata_filtered %>% dplyr::filter(metadata_filtered$SampleID %in% colnames(hiv_drug_log_filtered)) %>%
  arrange(match(SampleID, colnames(hiv_drug_log_filtered)))
all(metadata_filtered$SampleID == colnames(hiv_drug_log_filtered))

# heatmap for all people (HIV+ and HIV-) and all metadata
metadata_filtered_hm <- metadata_filtered %>% mutate_if(is.character, ~ ifelse(is.na(.), "unknown", .))
col_fun <- colorRamp2(c(0, 4.5, 5.5, 8.5), c("#336699","#336699", "#FFFF99", "#FF3333"))
ha = HeatmapAnnotation(
  hiv_status = metadata_filtered_hm$hiv_status_clean,
  arv_status = metadata_filtered_hm$arv_status,
  col = list(
    hiv_status = c("HIV+" = "darkblue", "HIV-" = "gold"),
    arv_status = c("HAART" = "#922091", "non-HAART" = "#ED708D", "No ARVs" = "#00BFC4", "ARV Naive" = "#068E38", "unknown" = "darkgray")
  )
)

hiv_drug_log_filtered <- as.matrix(hiv_drug_log_filtered)
hiv_all_heatmap <- Heatmap(hiv_drug_log_filtered, clustering_distance_columns = "euclidean", clustering_method_columns = "ward.D2",
                           name = "Log_Peak_Area", col = col_fun, cluster_rows = T, show_column_dend = T, 
                           show_column_names = F, show_row_names = T, 
                           row_gap = unit(0.3, "mm"), cluster_row_slices = F,
                           row_title_rot = 0, row_title_gp = gpar(fontsize = 8),
                           row_title_side = "right",
                           heatmap_legend_param = list(title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 8)),
                           bottom_annotation = ha)
heatmap_ggplot <- as.ggplot(grid.grabExpr(draw(hiv_all_heatmap)))
# ggsave(file="hiv_all_heatmap.svg", plot=heatmap_ggplot, width=13, height=4)

#####################################################
# feature table for HCA analysis in people with HIV #
#####################################################
all(metadata_filtered$SampleID == colnames(hiv_drug_log_filtered))
hiv_pos_drug_log <- hiv_drug_log_filtered %>% as.data.frame() %>%
  dplyr::select(any_of(metadata_filtered$SampleID[metadata_filtered$hiv_status_clean == "HIV+"])) %>%
  dplyr::filter(rowSums(. > 0) > 22) # only include drugs detected in >10% people in the clustering analysis
metadata_hiv_pos_only <- metadata_filtered_hm %>% 
  dplyr::filter(metadata_filtered_hm$hiv_status_clean == "HIV+") %>%
  arrange(match(SampleID, colnames(hiv_pos_drug_log)))
all(metadata_hiv_pos_only$SampleID == colnames(hiv_pos_drug_log))

# Clustering Step 1: Produce the dendrogram
col_fun <- colorRamp2(c(0, 4.5, 5.5, 8.5), c("#336699","#336699", "#FFFF99", "#FF3333"))
hiv_pos_drug_log <- as.matrix(hiv_pos_drug_log)
HIV_pos_heatmap <- Heatmap(hiv_pos_drug_log, clustering_distance_columns = "euclidean", clustering_method_columns = "ward.D2",
                           name = "Log_Peak_Area", col = col_fun, cluster_rows = T, show_column_dend = T, 
                           show_column_names = F, show_row_names = T, 
                           row_gap = unit(0.3, "mm"), cluster_row_slices = F,
                           row_title_rot = 0, row_title_gp = gpar(fontsize = 8),
                           row_title_side = "right",
                           heatmap_legend_param = list(title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 8)))
column_dendrogram <- column_dend(HIV_pos_heatmap) 
# Step 2: Cut the dendrogram into 4 groups
column_clusters <- cutree(as.hclust(column_dendrogram), k = 4)
# Step 3: Create a dataframe with sample names (columns) and their assigned clusters
sample_cluster_df <- data.frame(Sample = colnames(hiv_pos_drug_log), Cluster = column_clusters)
sample_cluster_df <- sample_cluster_df %>% arrange(match(Sample, colnames(hiv_pos_drug_log)))
sample_cluster_df$Cluster <- gsub(3,5, sample_cluster_df$Cluster)
sample_cluster_df$Cluster <- gsub(2,3, sample_cluster_df$Cluster)
sample_cluster_df$Cluster <- gsub(5,2, sample_cluster_df$Cluster)
all(sample_cluster_df$Sample == colnames(hiv_pos_drug_log))
# Step 4: Add the cluster information to the heatmap for splitting
split_heatmap <- Heatmap(hiv_pos_drug_log, 
                         clustering_distance_columns = "euclidean",  
                         clustering_method_columns = "ward.D2",  
                         name = "Log_Peak_Area", col = col_fun, 
                         cluster_rows = T, show_column_dend = T, show_column_names = F, show_row_dend = F,
                         column_dend_height = unit(2, "cm"),
                         show_row_names = T, row_gap = unit(0.3, "mm"), cluster_row_slices = F,
                         row_title_rot = 0, row_title_gp = gpar(fontsize = 8),
                         row_title_side = "right",
                         heatmap_legend_param = list(title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 8)),
                         column_split = factor(sample_cluster_df$Cluster, levels = c(2,1,3,4)))  
# heatmap_ggplot <- as.ggplot(grid.grabExpr(draw(split_heatmap)))
# ggsave(file="hiv_drug_clustered.svg", plot=heatmap_ggplot, width=10, height=4)

table(sample_cluster_df$Cluster)

# test N-acyl lipids by patient clusters 
annotated_lipids <- read.delim("HNRC_analysis/HNRC_N_acyl_lipids_annotation.tsv")
annotated_lipids$final_name <- paste(annotated_lipids$X.Scan., "_", annotated_lipids$Compound_Name, sep = "")

process_string <- function(x) {
  # Extract the substring starting from "Sample_" and before ".mzML"
  extracted <- sub(".*(Sample_.*)\\.mzML.*", "\\1", x)
  # Replace "_" with "." after "Sample_"
  modified <- sub("Sample_", "Sample.", gsub("_", ".", extracted))
  return(modified)
}
colnames(feature_filter) <- sapply(colnames(feature_filter), process_string)
annotated_lipids$X.Scan. <- as.character(annotated_lipids$X.Scan.)
feature_acyl_lipid <- feature_filter %>% rownames_to_column("FeatureID") %>% 
  right_join(annotated_lipids %>% dplyr::select(c("X.Scan.", "final_name")), by = c("FeatureID" = "X.Scan.")) %>%
  column_to_rownames("final_name") %>% dplyr::select(-FeatureID)

# filter for samples from people with HIV
feature_HIV_pos <- feature_acyl_lipid %>% dplyr::select(any_of(sample_cluster_df$Sample)) 
feature_normalized <- as.data.frame(t(apply(feature_HIV_pos, 1, function(x) x / max(x))))
feature_normalized_final <- feature_normalized %>% dplyr::filter((rowSums(.>0)) > 22) %>% t() %>% as.data.frame() %>%
  rownames_to_column("Sample") %>% left_join(sample_cluster_df, by = "Sample") %>% dplyr::select(-Sample) 
# only perform analysis on N-acyl lipids detected in >10% samples

# Reshape the data from wide to long format
data_box_long <- tidyr::pivot_longer(feature_normalized_final, cols = -Cluster, names_to = "Feature", values_to = "Value")

# Kruskal-Wallis followed by Wilcoxon with Benjamini-Hochberg adjustment
test_feature <- unique(data_box_long$Feature)
for(i in 1:length(test_feature)){
  temp <- data_box_long[data_box_long$Feature == test_feature[i],]
  print(test_feature[i])
  print(kruskal.test(Value ~ Cluster,
                     data = temp))
}
for(i in 1:length(test_feature)){
  temp <- data_box_long[data_box_long$Feature == test_feature[i],]
  print(test_feature[i])
  print(pairwise.wilcox.test(temp$Value, temp$Cluster, p.adjust.method = 'BH'))
}

# Box plot for features that are significantly different among the drug exposure groups
feature_final <- feature_normalized_final %>% dplyr::select(
  c("Cluster",
    "717_Candidate Histamine-C2:0 (delta mass:42.0104)",
    "894_Candidate Histamine-C2:0 (delta mass:42.0106)",
    "2118_Candidate Histamine-C4:0 (delta mass:70.0418)",
    "2112_Candidate Histamine-C4:0 (delta mass:70.0418)",
    "5981_Candidate Dopamine-C2:0 (delta mass:42.0105)",
    "5106_Candidate Spermidine-C5:0 (delta mass:84.0581)",
    "5534_Candidate Tyrosine-C2:0 (delta mass:42.0107)",
    "32710_Candidate Leucine-C18:1 (delta mass:264.2453)",
    "4169_Candidate N-acetylputrescine-C2:0 (delta mass:42.0103)"
  ))

data_box_long_sig <- tidyr::pivot_longer(feature_final, cols = -Cluster, names_to = "Feature", values_to = "Value")
data_box_long_sig$Cluster <- factor(data_box_long_sig$Cluster, levels = c("1", "2", "3", "4"))

# simplify N-acyl lipid names
data_box_long_sig$Feature <- gsub("717_Candidate Histamine-C2:0 (delta mass:42.0104)", "Histamine-C2:0", data_box_long_sig$Feature, fixed = T)
data_box_long_sig$Feature <- gsub("894_Candidate Histamine-C2:0 (delta mass:42.0106)", "Histamine-C2:0-isomer", data_box_long_sig$Feature, fixed = T)
data_box_long_sig$Feature <- gsub("2118_Candidate Histamine-C4:0 (delta mass:70.0418)", "Histamine-C4:0", data_box_long_sig$Feature, fixed = T)
data_box_long_sig$Feature <- gsub("2112_Candidate Histamine-C4:0 (delta mass:70.0418)", "Histamine-C4:0-isomer", data_box_long_sig$Feature, fixed = T)
data_box_long_sig$Feature <- gsub("5981_Candidate Dopamine-C2:0 (delta mass:42.0105)", "Dopamine-C2:0", data_box_long_sig$Feature, fixed = T)
data_box_long_sig$Feature <- gsub("5106_Candidate Spermidine-C5:0 (delta mass:84.0581)", "Spermidine-C5:0", data_box_long_sig$Feature, fixed = T)
data_box_long_sig$Feature <- gsub("5534_Candidate Tyrosine-C2:0 (delta mass:42.0107)", "Tyrosine-C2:0", data_box_long_sig$Feature, fixed = T)
data_box_long_sig$Feature <- gsub("32710_Candidate Leucine-C18:1 (delta mass:264.2453)", "Leucine-C18:1", data_box_long_sig$Feature, fixed = T)
data_box_long_sig$Feature <- gsub("4169_Candidate N-acetylputrescine-C2:0 (delta mass:42.0103)", "N-acetylputrescine-C2:0", data_box_long_sig$Feature, fixed = T)

# Create box plots
my_comparisons <- list( c("1", "4"),  c("2", "4"), c("3", "4"))
data_box_long_sig$Feature <- factor(data_box_long_sig$Feature,
                                    levels = c("Histamine-C2:0",
                                               "Histamine-C2:0-isomer",
                                               "Histamine-C4:0",
                                               "Histamine-C4:0-isomer",
                                               "Dopamine-C2:0",
                                               "Spermidine-C5:0",
                                               "Tyrosine-C2:0",
                                               "Leucine-C18:1",
                                               "N-acetylputrescine-C2:0"))
acyl_lipids_boxplot <- ggplot(data_box_long_sig, aes(x = Cluster, y = Value, fill = Cluster)) +
  geom_boxplot(outlier.shape=NA, alpha = 0.6, size = 1) +
  geom_point(position = position_jitterdodge(jitter.width = 0.8), size = 4, shape = 20, aes(color = Cluster), alpha = 0.4) +
  facet_wrap(~ Feature, scales = "free_y", ncol = 3) + 
  scale_color_manual(values = c("#922091","darkred", "#ED708D","#336699")) +
  scale_fill_manual(values = c("#922091","darkred", "#ED708D","#336699")) +
  theme_classic() +
  theme(axis.ticks.length=unit(.15, "cm"),
        axis.line.x.bottom=element_line(linewidth=1.0), axis.line.y.left = element_line(linewidth=1.0),
        axis.ticks.x.bottom = element_blank(), axis.ticks.y.left = element_line(linewidth=1.0),
        axis.text.y.left = element_blank(), axis.text.x.bottom = element_blank(),
        axis.title.y.left = element_blank(), axis.title.x.bottom = element_blank(),
        strip.background = element_blank(), strip.text = element_text(size = 12), face = "bold") + 
  stat_compare_means(comparisons = my_comparisons, hide.ns = F, symnum.args = list(cutpoints = c(0, 0.01, 0.05, Inf), symbols = c("**", "*", "ns")))
# ggsave(file="acyl_lipid_boxplot.svg", plot=acyl_lipids_boxplot, width=10, height=7)

# test N-acyl lipids by ARV_status in clinical metadata
feature_normalized_arv_status <- feature_normalized %>% dplyr::filter((rowSums(.>0)) > 22) %>% t() %>% as.data.frame() %>%
  rownames_to_column("sample") %>% left_join(metadata_hiv_pos_only[,c("SampleID", "arv_status")], by = c("sample" = "SampleID")) %>% dplyr::select(-sample) %>%
  dplyr::filter(arv_status != "unknown")

# Reshape the data from wide to long format
data_box_long_arv_status <- tidyr::pivot_longer(feature_normalized_arv_status, cols = -arv_status, names_to = "Feature", values_to = "Value")

# Kruskal-Wallis followed by Wilcoxon with Benjamini-Hochberg adjustment
feature_test_arv_status <- unique(data_box_long_arv_status$Feature)
for(i in 1:length(feature_test_arv_status)){
  temp <- data_box_long_arv_status[data_box_long_arv_status$Feature == feature_test_arv_status[i],]
  print(feature_test_arv_status[i])
  print(kruskal.test(Value ~ arv_status,
                     data = temp))
}

feature_final_arv_status <- feature_normalized_arv_status %>% dplyr::select(
  c("arv_status",
    "717_Candidate Histamine-C2:0 (delta mass:42.0104)", 
    "894_Candidate Histamine-C2:0 (delta mass:42.0106)", 
    "2118_Candidate Histamine-C4:0 (delta mass:70.0418)",
    "2112_Candidate Histamine-C4:0 (delta mass:70.0418)", 
    "5981_Candidate Dopamine-C2:0 (delta mass:42.0105)", 
    "5106_Candidate Spermidine-C5:0 (delta mass:84.0581)", 
    "5534_Candidate Tyrosine-C2:0 (delta mass:42.0107)", 
    "32710_Candidate Leucine-C18:1 (delta mass:264.2453)", 
    "4169_Candidate N-acetylputrescine-C2:0 (delta mass:42.0103)"
  ))
data_box_long_sig_arv_status <- tidyr::pivot_longer(feature_final_arv_status, cols = -arv_status, names_to = "Feature", values_to = "Value")
data_box_long_sig_arv_status$arv_status <- factor(data_box_long_sig_arv_status$arv_status, levels = c("ARV Naive", "No ARVs", "non-HAART", "HAART"))

data_box_long_sig_arv_status$Feature <- gsub("717_Candidate Histamine-C2:0 (delta mass:42.0104)", "Histamine-C2:0", data_box_long_sig_arv_status$Feature, fixed = T)
data_box_long_sig_arv_status$Feature <- gsub("894_Candidate Histamine-C2:0 (delta mass:42.0106)", "Histamine-C2:0-isomer", data_box_long_sig_arv_status$Feature, fixed = T)
data_box_long_sig_arv_status$Feature <- gsub("2118_Candidate Histamine-C4:0 (delta mass:70.0418)", "Histamine-C4:0", data_box_long_sig_arv_status$Feature, fixed = T)
data_box_long_sig_arv_status$Feature <- gsub("2112_Candidate Histamine-C4:0 (delta mass:70.0418)", "Histamine-C4:0-isomer", data_box_long_sig_arv_status$Feature, fixed = T)
data_box_long_sig_arv_status$Feature <- gsub("5981_Candidate Dopamine-C2:0 (delta mass:42.0105)", "Dopamine-C2:0", data_box_long_sig_arv_status$Feature, fixed = T)
data_box_long_sig_arv_status$Feature <- gsub("5106_Candidate Spermidine-C5:0 (delta mass:84.0581)", "Spermidine-C5:0", data_box_long_sig_arv_status$Feature, fixed = T)
data_box_long_sig_arv_status$Feature <- gsub("5534_Candidate Tyrosine-C2:0 (delta mass:42.0107)", "Tyrosine-C2:0", data_box_long_sig_arv_status$Feature, fixed = T)
data_box_long_sig_arv_status$Feature <- gsub("32710_Candidate Leucine-C18:1 (delta mass:264.2453)", "Leucine-C18:1", data_box_long_sig_arv_status$Feature, fixed = T)
data_box_long_sig_arv_status$Feature <- gsub("4169_Candidate N-acetylputrescine-C2:0 (delta mass:42.0103)", "N-acetylputrescine-C2:0", data_box_long_sig_arv_status$Feature, fixed = T)

# Create box plots
data_box_long_sig_arv_status$Feature <- factor(data_box_long_sig_arv_status$Feature,
                                               levels = c("Histamine-C2:0",
                                                          "Histamine-C2:0-isomer",
                                                          "Histamine-C4:0",
                                                          "Histamine-C4:0-isomer",
                                                          "Dopamine-C2:0",
                                                          "Spermidine-C5:0",
                                                          "Tyrosine-C2:0",
                                                          "Leucine-C18:1",
                                                          "N-acetylputrescine-C2:0"))
acyl_lipids_boxplot_arv_status <- ggplot(data_box_long_sig_arv_status, aes(x = arv_status, y = Value, fill = arv_status)) +
  geom_boxplot(outlier.shape=NA, alpha = 0.6, size = 1) +
  geom_point(position = position_jitterdodge(jitter.width = 0.8), size = 2, shape = 20, aes(color = arv_status), alpha = 0.4) +
  facet_wrap(~ Feature, scales = "free_y", ncol = 5) + 
  scale_color_manual(values = c("#922091","darkred", "#ED708D","#336699")) +
  scale_fill_manual(values = c("#922091","darkred", "#ED708D","#336699")) +
  theme_classic() +
  theme(axis.ticks.length=unit(.15, "cm"),
        axis.line.x.bottom=element_line(linewidth=1.0), axis.line.y.left = element_line(linewidth=1.0),
        axis.ticks.x.bottom = element_blank(), axis.ticks.y.left = element_line(linewidth=1.0),
        axis.text.y.left = element_blank(), axis.text.x.bottom = element_blank(),
        axis.title.y.left = element_blank(), axis.title.x.bottom = element_blank(),
        strip.background = element_blank(), strip.text = element_text(size = 12), face = "bold")
# ggsave(file="acyl_lipid_boxplot_arv.svg", plot=acyl_lipids_boxplot_arv_status, width=9, height=4)

#############################################
# microbial culture boxplot of drug analogs #
#############################################

feature_microbe <- read.csv("HNRC_analysis/SynCom_quant.csv")
feature_filter <- feature_microbe %>% column_to_rownames("row.ID") %>% dplyr::select(contains("mzML"))
colnames(feature_filter) <- sub(".mzML.Peak.area","", colnames(feature_filter))
metadata <- read.csv("HNRC_analysis/SynCom_metadata.csv")
feature_filter_analogs <- feature_filter %>% rownames_to_column("FeatureID") %>%
  dplyr::filter(FeatureID %in% c("6222", "7256", "13875", "13809", "13613", "11007", "12995", "11849", "11178")) %>%
  column_to_rownames("FeatureID") %>% t() %>% as.data.frame() %>% rownames_to_column("sample") %>%
  left_join(metadata, by = "sample") # analog features shared between HNRC and microbial cultures that are not ISF/adducts

colnames(feature_filter_analogs)[2:10] <- c("F6222", "F11849", "F11007", "F7256", "F11178", "F13809", "F13613", "F13875", "F12995")
all_analogs <- data.frame(analog = c("F6222", "F7256", "F13875", "F13809", "F13613", "F11007", "F12995", "F11849", "F11178"), 
                          drug_group = c(4, 4, 3, 3, 3, 2, 2, 1, 3))
boxplot <- data.frame(final_peak_area = NULL, final_group = NULL, feature = NULL)

# blank subtraction with BHI media only features
for (i in 1:9){
  if (i<9){
    temp <- feature_filter_analogs %>% dplyr::filter(drug_group == all_analogs$drug_group[i])
    mean_BHI_0 <- mean(temp %>% dplyr::filter(group == "BHI" & time == "T0") %>% dplyr::pull(all_analogs$analog[i])) 
    mean_BHI_72 <- mean(temp %>% dplyr::filter(group == "BHI" & time == "T72") %>% dplyr::pull(all_analogs$analog[i])) 
    blk_sub_peak_area <- rep(NA, 42)
    raw_peak_area <- temp %>% dplyr::pull(all_analogs$analog[i])
    blk_sub_peak_area[temp$time == "T0"] <- raw_peak_area[temp$time == "T0"] - mean_BHI_0
    blk_sub_peak_area[temp$time == "T72"] <- raw_peak_area[temp$time == "T72"] - mean_BHI_72
    blk_sub_peak_area[blk_sub_peak_area<0] <- 0
    result <- data.frame(final_peak_area = blk_sub_peak_area, final_group = temp$final_group, feature = rep(all_analogs$analog[i], 42))
    boxplot <- rbind(boxplot, result)
  }
  else{
    temp <- feature_filter_analogs %>% dplyr::filter(drug_group == all_analogs$drug_group[i])
    result <- data.frame(final_peak_area = temp$F11178, final_group = temp$final_group, feature = rep(all_analogs$analog[i], 42))
    boxplot <- rbind(boxplot, result)
  }
} 

boxplot_final <- boxplot %>% dplyr::filter(!grepl("BHI",boxplot$final_group))
boxplot_final$final_group <- factor(boxplot_final$final_group, levels = c(
  "SynCom_T0", "SynCom_T72", "SynCom reduced_T0", "SynCom reduced_T72",
  "Fast growers_T0", "Fast growers_T72", "Medium fast growers_T0", "Medium fast growers_T72",
  "Medium slow growers_T0", "Medium slow growers_T72", "Slow growers_T0", "Slow growers_T72"
))

y_axis_limits <- list(
  "F11178" = c(0, 3000000), "F12995" = c(0, 30000), "F13875" = c(0, 40000), "F6222" = c(0, 300000), "F7256" = c(0, 300000),
  "F13809" = c(0, 60000), "F13613" = c(0, 80000), "F11007" = c(0, 120000), "F11849" = c(0, 50000))
plots <- lapply(names(y_axis_limits), function(feat) {
  ggplot(subset(boxplot_final, feature == feat), aes(x = final_group, y = final_peak_area, fill = final_group)) +
    geom_boxplot(outlier.shape=NA, alpha = 0.7, size = 1, col = "black", linewidth = 0.7) +
    geom_point(position = position_dodge(width = 0.75), size = 4, shape = 21, aes(color = final_group), col = "black", stroke = 0.7) +
    scale_color_manual(values = c("#e3716e","#e3716e","#7ac7e2","#7ac7e2","#54beaa","#54beaa","#eca680","#eca680","#f7df87","#f7df87","#af8fd0","#af8fd0")) +
    scale_fill_manual(values = c("#e3716e","#e3716e","#7ac7e2","#7ac7e2","#54beaa","#54beaa","#eca680","#eca680","#f7df87","#f7df87","#af8fd0","#af8fd0")) +
    coord_cartesian(ylim = y_axis_limits[[feat]]) +  # Set y-axis limits
    theme_classic() +
    theme(axis.ticks.length=unit(.15, "cm"),
          axis.line.x.bottom=element_line(linewidth=0.7), axis.line.y.left = element_line(linewidth=0.7),
          axis.ticks.x.bottom = element_line(linewidth = 0.7), axis.ticks.y.left = element_line(linewidth=0.7),
          axis.text.y.left = element_text(size = 12, hjust = -0.2), axis.text.x.bottom = element_text(angle = 45, hjust = 1),
          axis.title.y.left = element_blank(), axis.title.x.bottom = element_blank(),
          legend.position = "none")
})
combined_plot_SynCom <- wrap_plots(plots, ncol = 3)
# ggsave("combined_plot_SynCom.svg", plot = combined_plot_SynCom, width = 12, height = 20)

# Wilcox test to compare T0 and T72
boxplot_final <- boxplot_final %>% mutate(microbe_group = str_split(final_group, "_", simplify = TRUE)[, 1])
results_list <- list()
test_feature <- unique(boxplot_final$feature)
test_group <- unique(boxplot_final$microbe_group)
for (i in 1:length(test_feature)) {
  for (j in 1:length(test_group)) {
    temp <- boxplot_final %>%
      filter(feature == test_feature[i], microbe_group == test_group[j])
    
    # Only run test if we have at least 2 groups to compare
    if (length(unique(temp$final_group)) > 1 && nrow(temp) > 1) {
      test <- wilcox.test(final_peak_area ~ final_group, data = temp, exact = TRUE)
      p_val <- test$p.value
    } else {
      p_val <- NA  # insufficient data
    }
    
    # Store results
    results_list[[length(results_list) + 1]] <- data.frame(
      feature = test_feature[i],
      microbe_group = test_group[j],
      p_value = p_val
    )
  }
}
pval_df <- do.call(rbind, results_list) # feature 13613 excluded due to no RT match; feature 11849 excluded due to no significant peak area increase from T0 to T72

#######################
# wastewater analysis #
#######################
feature <- read.csv("wastewater_analysis/FBMN_Barcelona_WWTP_pos_quant.csv")
feature_filter <- feature %>% column_to_rownames("row.ID") %>% dplyr::select(contains(".mzML"))
colnames(feature_filter) <- sub(".mzML.Peak.area", "", colnames(feature_filter))
metadata <- data.frame(filename = colnames(feature_filter))
metadata$EDTA <- ifelse(grepl("EDTA", metadata$filename), "yes", "no")
metadata <- metadata %>%  mutate(class = str_extract(filename, "^[^_]+_[^_]+"),  
                                 dates = str_match(filename, "^[^_]+_[^_]+_([0-9]{4})")[, 2],
                                 month = case_when(
                                   str_starts(dates, "03") ~ "March",
                                   str_starts(dates, "04") ~ "April",
                                   str_starts(dates, "05") ~ "May-June",
                                   str_starts(dates, "06") ~ "May-June"))
# drug library annotation
annotation <- read.csv("wastewater_analysis/Barcelona_WWTP_annotation_pos.tsv", sep = "\t") #https://gnps2.org/status?task=4929dd75212743f48fa0383c5563456b
annotation_filtered <- annotation %>% dplyr::filter(MQScore > 0.9 | SharedPeaks > 4) %>% 
  dplyr::select(c("X.Scan.", "SpectrumID", "Compound_Name")) 
colnames(annotation_filtered) <- c("FeatureID", "SpectrumID", "Compound_Name")
annotation_filtered$FeatureID <- as.character(annotation_filtered$FeatureID)

druglib_metadata <- read.csv("wastewater_analysis/GNPS_Drug_Library_drug_metadata.csv")
druglib_metadata_filtered <- druglib_metadata %>% 
  dplyr::select(c("gnps_libid", "name_parent_compound", "chemical_source", "pharmacologic_class", "therapeutic_area", "therapeutic_indication"))
analoglib_metadata <- read.csv("wastewater_analysis/GNPS_Drug_Library_drug_analog_metadata.csv")
analoglib_metadata_filtered <- analoglib_metadata %>% 
  dplyr::select(c("analog_libid", "name_connected_compound", "chemical_source", "pharmacologic_class", "therapeutic_area", "therapeutic_indication"))
colnames(analoglib_metadata_filtered)[1:2] <- c("gnps_libid", "name_parent_compound")
all_metadata <- rbind(druglib_metadata_filtered, analoglib_metadata_filtered)
annotation_metadata <-  merge(annotation_filtered, all_metadata, by.x = "SpectrumID", by.y = "gnps_libid", all.x = T, all.y=F)
annotation_metadata <- annotation_metadata[!duplicated(annotation_metadata$FeatureID),]

annotation_metadata$FeatureID <- as.character(annotation_metadata$FeatureID)
feature_annotation <- annotation_metadata %>% inner_join((feature_filter %>% rownames_to_column("FeatureID")), by = "FeatureID")
feature_annotation$chemical_source[feature_annotation$FeatureID == 8313] <- "Medical"
feature_annotation$chemical_source[feature_annotation$FeatureID == 21396] <- "Medical"
feature_annotation$chemical_source[feature_annotation$FeatureID == 7794] <- "Medical"
feature_annotation$chemical_source[feature_annotation$FeatureID == 6895] <- "Medical" # these are internal standards used in Dorrestein Lab. They should not be rejected when analyzing data from other labs

exo_drug <- feature_annotation %>% dplyr::filter(!grepl("Background|confidence|Endogenous|Food", feature_annotation$chemical_source))
#write.csv(exo_drug, "wastewater_analysis/exo_drug_manual_check.csv")

exo_drug <- read.csv("wastewater_analysis/exo_drug_manual_checked.csv") # cleanup name of parent drugs
exo_drug$group_plot <- "others"
exo_drug$group_plot[grepl("hypertension",exo_drug$therapeutic_indication)] <- "antihypertensives"
exo_drug$group_plot[grepl("HIV",exo_drug$therapeutic_indication)] <- "HIV medications"
exo_drug$group_plot[grepl("depression",exo_drug$therapeutic_indication)] <- "antidepressants"
exo_drug$group_plot[grepl("seizures",exo_drug$therapeutic_indication)] <- "antiepileptics"
exo_drug$group_plot[grepl("antimicrobial|antibacterial",exo_drug$pharmacologic_class) & !grepl("veterinary",exo_drug$pharmacologic_class)] <- "antibiotics"
exo_drug$group_plot[grepl("cough suppressant",exo_drug$therapeutic_indication)] <- "cough suppressant"

exo_drug_class_grouped <- exo_drug %>% dplyr::select(any_of(c("group_plot", metadata$filename[!grepl("BL", metadata$filename)]))) %>%
  dplyr::filter(group_plot != "others") %>%
  group_by(group_plot) %>%
  dplyr::summarise(across(everything(), sum)) %>%
  column_to_rownames("group_plot")

exo_drug_class_grouped_normalized <- as.data.frame(t(apply(exo_drug_class_grouped, 1, function(x) x / max(x))))
exo_drug_class_grouped_normalized <- exo_drug_class_grouped_normalized %>% rownames_to_column("groups_plot")

# Create boxplot
exo_drug_class_grouped_long <- exo_drug_class_grouped_normalized %>%
  pivot_longer(cols = -groups_plot, names_to = "sample", values_to = "abundance")
combined_data <- exo_drug_class_grouped_long %>%
  left_join(metadata, by = c("sample" = "filename")) 
combined_data$month <- factor(combined_data$month, levels = c("March", "April", "May-June"))
WWTP_boxplot <- ggplot(combined_data, aes(x = class, y = abundance, fill = month)) +
  geom_boxplot(outlier.shape=NA, alpha = 0.8, size = 1) +
  geom_point(position = position_jitterdodge(jitter.width = 0), size = 3, shape = 20, aes(color = month), alpha = 0.6) +
  facet_wrap(~ groups_plot, scales = "free_y", ncol = 3) + 
  scale_color_manual(values = c("#EB716B","#EE8227", "#79ADD6")) +
  scale_fill_manual(values = c("#EB716B","#EE8227", "#79ADD6")) +
  labs(x = "Class") +
  theme_classic() +
  theme(axis.ticks.length=unit(.15, "cm"),
        axis.line.x.bottom=element_line(linewidth=1.0), axis.line.y.left = element_line(linewidth=1.0),
        axis.ticks.x.bottom = element_line(linewidth=1.0), axis.ticks.y.left = element_line(linewidth=1.0),
        axis.text.y.left = element_blank(), axis.text.x.bottom = element_text(size = 12),
        axis.title.y.left = element_blank(), axis.title.x.bottom = element_blank(),
        strip.background = element_blank(), strip.text = element_text(size = 12)) 
# ggsave(file="WWTP_boxplot.svg", plot=WWTP_boxplot, width=10, height=6)

# Kruskal-Wallis followed by Wilcoxon with Benjamini-Hochberg adjustment
group_test <- unique(combined_data$groups_plot)
plant <- unique(combined_data$class)
for(i in 1:length(group_test)){
  for (j in 1:3){
    temp <- combined_data[combined_data$groups_plot == group_test[i] & combined_data$class == plant[j],]
    print(paste0(group_test[i], "___", plant[j]))
    print(kruskal.test(abundance ~ month,
                       data = temp))
  }}
for(i in 1:length(group_test)){
  for (j in 1:3){
    temp <- combined_data[combined_data$groups_plot == group_test[i] & combined_data$class == plant[j],]
    print(paste0(group_test[i], "___", plant[j]))
    print(pairwise.wilcox.test(temp$abundance, temp$month, p.adjust.method = 'BH'))
  }
}
