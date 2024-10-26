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

####################################################################
# Sankey plot for drug library metadata based on number of spectra #
####################################################################

drug_metadata <- read.csv("druglib_metadata/GNPS_Drug_Library_drug_metadata.csv")
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
drug_metadata <- read.csv("druglib_metadata/GNPS_Drug_Library_drug_metadata.csv")
drug_metadata[is.na(drug_metadata)] <- "NA"
drug_metadata$therapeutic_area[drug_metadata$therapeutic_area == "NA"] <- "unspecified_area"
metadata_w_structure <- drug_metadata[(drug_metadata$inchikey != "N/A"),]
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
counts_first <- ddply(metadata_separate, .(metadata_separate$cleanup_source, metadata_separate$therapeutic_area), nrow)
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

names(counts_first) <- c("source", "target", "value")
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

##########################################################################
# co-occurrence frequency of drug and drug analogs in fastMASST searches #
##########################################################################

# Code to run repository-scale search of the GNPS Drug Library with fastMASST in batch: https://github.com/robinschmid/microbe_masst/tree/v1.2.0
# Output for each spectrum is a tsv files containing all the USI matched to the spectrum
# All the tsv file has been pasted together into a single csv provided with this code
# The resulting files are too large to upload to GitHub, request access from haz072@health.ucsd.edu

# masst results for all parent drugs
drug_masst <- read.csv("masst_results_entire_library/druglib_parent_drug_fasst_matches.csv")
drug_masst_filtered <- drug_masst %>% dplyr::select(colnames(drug_masst)[-1]) %>%
  dplyr::filter(!grepl("MSV000084314|MSV000089210", drug_masst$USI)) # the excluded datasets are re-analysis repositories from the lab - not real data
drug_masst_filtered$short_usi <- sub(":scan.*", "", drug_masst_filtered$USI)
colnames(drug_masst_filtered) <- paste("drug_", colnames(drug_masst_filtered), sep = "")
drug_masst_filtered <- drug_masst_filtered %>% dplyr::select(c("drug_GNPSID", "drug_short_usi"))

# masst results for all drug analogs
analog_masst <- read.csv("masst_results_entire_library/druglib_drug_analogs_fasst_matches.csv")
analog_masst_filtered <- analog_masst %>% dplyr::select(colnames(analog_masst)[-1]) %>%
  dplyr::filter(!grepl("MSV000084314|MSV000089210", analog_masst$USI))
analog_masst_filtered$short_usi <- sub(":scan.*", "", analog_masst_filtered$USI)
colnames(analog_masst_filtered) <- paste("analog_", colnames(analog_masst_filtered), sep = "")
analog_masst_filtered <- analog_masst_filtered %>% dplyr::select(c("analog_GNPSID", "analog_short_usi"))

# drug and drug analog metadata
drug_metadata <- read.csv("druglib_metadata/GNPS_Drug_Library_drug_metadata.csv")
analog_metadata <- read.csv("druglib_metadata/GNPS_Drug_Library_drug_analog_metadata.csv")
analog_id_match <- read.csv("masst_results_entire_library/drug_analoglib_ID_match.csv") 
# fasst search was performed with older version of the analog library, which was re-uploaded to GNPS after filtering to generate the final ccmslibid 
analog_metadata <- analog_metadata %>% left_join(analog_id_match, by = c("analog_libid" = "new_id"))
analog_metadata_filtered <- analog_metadata %>% dplyr::select(c("old_id", "parent_drug_libid", "chemical_source", "delta_mass")) %>%
  dplyr::filter(!is.na(analog_metadata$parent_drug_libid))

analog_masst_metadata <- inner_join(analog_masst_filtered, analog_metadata_filtered, by = c("analog_GNPSID" = "old_id"))

# group all USI matched to a specific drug spectrum
custom_aggregate <- function(x) {
  return(paste(unique(x), collapse = "|"))
}
drug_masst_filename_cat <- drug_masst_filtered %>% 
  group_by(drug_short_usi) %>% 
  summarise_all(custom_aggregate)

# check if drug analog showed up in the same data file as the parent drug
analog_masst_metadata$parent_check <- rep(FALSE, nrow(analog_masst_metadata))
for (i in 1:nrow(analog_masst_metadata)){
  if (analog_masst_metadata$analog_short_usi[i] %in% drug_masst_filename_cat$drug_short_usi){
    temp_df <- drug_masst_filename_cat %>% 
      dplyr::filter(drug_masst_filename_cat$drug_short_usi == analog_masst_metadata$analog_short_usi[i])
    if (grepl(analog_masst_metadata$parent_drug_libid[i], temp_df$drug_GNPSID)){
      analog_masst_metadata$parent_check[i] <- TRUE
    }
  }
}

# check if drug analog appear at least in one file with parent drug
cooccurrence_cal <- analog_masst_metadata %>% dplyr::select(c("analog_GNPSID", "parent_check"))
cooccurrence_check <- cooccurrence_cal %>%
  group_by(analog_GNPSID) %>%
  dplyr::summarize(cooccurrence = any(parent_check == TRUE))

# group cooccurrence check based on delta mass
analog_masst_final <- cooccurrence_check %>% 
  left_join(analog_metadata %>% dplyr::select(c("old_id", "delta_mass")), 
            by = c("analog_GNPSID" = "old_id"))
delta_separate <- tidyr::separate_rows(analog_masst_final, delta_mass, sep = "\\|")
cooccurrence_by_delta <- as.data.frame(table(delta_separate$delta_mass))
colnames(cooccurrence_by_delta) <- c("delta_mass", "number_of_drug_analog")
cooccurrence_by_delta$masst_detection <- rep(NA, nrow(cooccurrence_by_delta))
for(i in cooccurrence_by_delta$delta_mass){
  temp <- delta_separate %>% dplyr::filter(delta_separate$delta_mass == i)
  cooccurrence_by_delta$masst_detection[cooccurrence_by_delta$delta_mass == i] <- sum(temp$cooccurrence == TRUE)
}
write.csv(cooccurrence_by_delta, "masst_results_entire_library/drug_analog_cooccurrence_with_drug_by_delta_mass.csv", row.names = F)


######################################
# delta mass frequency visualization #
######################################

cooccurrence_by_delta <- read.csv("masst_results_entire_library/drug_analog_cooccurrence_with_drug_by_delta_mass.csv")
cooccurrence_by_delta$atom <- rep("unknown",nrow(cooccurrence_by_delta))

# add atomic composition for frequent delta masses
cooccurrence_by_delta$atom[cooccurrence_by_delta$delta_mass == -15.99 | cooccurrence_by_delta$delta_mass == 15.99] <- "1O" 
cooccurrence_by_delta$atom[cooccurrence_by_delta$delta_mass == 14.02 | cooccurrence_by_delta$delta_mass == -14.02] <- "1C, 2H" 
cooccurrence_by_delta$atom[cooccurrence_by_delta$delta_mass == 1] <- "C isotope" 
cooccurrence_by_delta$atom[cooccurrence_by_delta$delta_mass == 2] <- "O/S/Cl/Br isotope" 
cooccurrence_by_delta$atom[cooccurrence_by_delta$delta_mass == 28.03 | cooccurrence_by_delta$delta_mass == -28.03] <- "2C, 4H" 
cooccurrence_by_delta$atom[cooccurrence_by_delta$delta_mass == 18.01 | cooccurrence_by_delta$delta_mass == -18.01] <- "1O, 2H" 
cooccurrence_by_delta$atom[cooccurrence_by_delta$delta_mass == 27.99 | cooccurrence_by_delta$delta_mass == -27.99] <- "1C, 1O" 
cooccurrence_by_delta$atom[cooccurrence_by_delta$delta_mass == 2.02 | cooccurrence_by_delta$delta_mass == -2.02] <- "2H" 

cooccurrence_by_delta$atom <- factor(cooccurrence_by_delta$atom, levels = c("1O", "1C, 2H", "C isotope", "O/S/Cl/Br isotope", "2C, 4H", "1O, 2H","1C, 1O","2H","unknown"))

barplot <- ggplot(cooccurrence_by_delta, aes(x=delta_mass, weight=number_of_drug_analog, fill = atom)) +
  geom_histogram(binwidth = 1, boundary = 0.5, color = "white", size = 0.5) +
  scale_fill_manual(values = c("darkred", "#339933", "#0066CC", "#66cccc", "orange", "pink",
                               "#330066", "#9999CC", "black")) + 
  xlab("Delta mass") +
  ylab("Spectra number") +
  scale_x_continuous(limits = c(-100, 100), breaks = c(-100, -50, 0, 50, 100))+
  scale_y_continuous(limits = c(0, 420), breaks = c(0, 100, 200, 300, 400), expand = c(0, 0))+
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
# ggsave(plot = plot, file = "analog_delta_mass_barplot.svg", width = 10, height = 6)

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
# fasst search was performed with older version of the analog library, which was re-uploaded to GNPS after filtering to generate the final ccmslibid 
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
annotation <- read.csv("HNRC_analysis/HNRC_DrugLib_annotations.tsv", sep = "\t") # https://gnps.ucsd.edu/ProteoSAFe/status.jsp?task=90b6edcfb8fc43dab54d2d7951a3291a
annotation_filtered <- annotation %>% dplyr::filter(MQScore > 0.9 | SharedPeaks > 4) %>% 
  dplyr::select(c("X.Scan.", "SpectrumID", "Compound_Name")) 
colnames(annotation_filtered) <- c("FeatureID", "SpectrumID", "Compound_Name")

# read in Drug Library metadata
druglib_metadata <- read.csv("druglib_metadata/GNPS_Drug_Library_drug_metadata.csv")
druglib_metadata_nodup <- druglib_metadata[!duplicated(druglib_metadata$gnps_libid),]
analoglib_metadata <- read.csv("druglib_metadata/GNPS_Drug_Library_drug_analog_metadata.csv")
druglib_metadata_filtered <- druglib_metadata_nodup %>% 
  dplyr::select(c("gnps_libid", "name_compound", "name_parent_compound", "chemical_source", "pharmacologic_class", "therapeutic_area", "therapeutic_indication"))
analoglib_metadata_filtered <- analoglib_metadata %>% 
  dplyr::select(c("analog_libid","name_analog","name_connected_compound", "chemical_source", "pharmacologic_class", "therapeutic_area", "therapeutic_indication"))
colnames(analoglib_metadata_filtered) <- colnames(druglib_metadata_filtered)
all_metadata <- rbind(druglib_metadata_filtered, analoglib_metadata_filtered)
annotation_metadata <-  merge(annotation_filtered, all_metadata, by.x = "SpectrumID", by.y = "gnps_libid", all.x = T, all.y=F)
exo_drug <- annotation_metadata %>% dplyr::filter(!grepl("Background|confidence|Endogenous|Food", annotation_metadata$chemical_source))
write.csv(exo_drug,"HNRC_analysis/HNRC_exogenous_drug_for_manual_name_cleaning.csv")
# The exported table was manually cleaned up by:
# (a) Clean up the "name_parent_compound" column to consistent names for parent drugs;
# (b) Assign annotation type for each unique drug as "drug only", "drug with metabolites/analogs", "metabolites only", or "analog only"
# (c) Reject annotation for feature 9562, 9198, 9388, 8088 - they are quinine or quindine related annotations that were rejected by standard check
# (d) Reject annotation for feature 13703 - this is nevirapine annotation that were rejected by standard check
# (e) For drugs with metabolites/analogs, their pharmacologic classes were manually cleaned up for better visualization
# (f) Results is provided in "HNRC_exogenous_drug_manually_cleaned.csv"

# read in the annotation table after manual cleanup
exo_drug_cleaned <- read.csv("HNRC_analysis/HNRC_exogenous_drug_manually_cleaned.csv")
exo_drug_with_analog <- exo_drug_cleaned %>% dplyr::filter(Group == "both") %>% arrange(Final_class)

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
drug_keep <- rownames(exo_drug_group)[rowSums(exo_drug_group>0)>32]

# produce the heatmap for drugs and analogs
exo_drug_heatmap <- feature_annotation %>% dplyr::filter(feature_annotation$cleanup_name_drug_only %in% drug_keep) 
exo_drug_heatmap_clean <- exo_drug_heatmap %>% dplyr::select(contains(".mzML"))
exo_drug_heatmap_filtered <- exo_drug_heatmap_clean[,colSums(exo_drug_heatmap_clean)>0]
exo_drug_heatmap_log <- log10(exo_drug_heatmap_filtered+1)
exo_drug_heatmap_log <- as.matrix(exo_drug_heatmap_log)

col_fun <- colorRamp2(c(0, 4, 5, 8), c("#336699", "#336699", "#FFFF99", "#FF3333"))
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
  Heatmap(exo_drug_heatmap$Final_class, name = "class", width = unit(5, "mm"),
          col = c("alpha adrenergic blocker" = "#669900", "Alzheimer's disease medication" = "#999966", "antiarrhythmic" = "#000000",
                  "antibiotics" = "#990033", "antidepressant" = "#CC6633", "antifungal" = "#003366", "antihistamine" = "#66CC99",
                  "antimalaria" = "#336666", "antiseptic" = "#333333", "anxiolytics" = "#339933", "atypical antipsychotic" = "#99CC99", 
                  "cardiology medication" = "#3399CC", "histamine-2 receptor antagonist" = "#66CCCC",
                  "HIV medication" = "#330066", "insomnia medication" = "#9999CC", 
                  "muscle relaxant" = "#999999", "proton pump inhibitor" = "#666600", "statin" = "#CCCCCC"),
          heatmap_legend_param = list(title_gp = gpar(fontsize = 24), labels_gp = gpar(fontsize = 24))
  )
# heatmap_ggplot <- as.ggplot(grid.grabExpr(draw(exo_drug_heatmap)))
# ggsave(file="exo_drug_heatmap.png", plot=heatmap_ggplot, width=8, height=11.5, dpi = 800)

# heatmap darunavir
darunavir <- feature_annotation %>% dplyr::filter(cleanup_name_drug_only == "darunavir")
darunavir$final_name <- paste(darunavir$cleanup_name, darunavir$FeatureID, sep = " | ")
darunavir_filtered <- darunavir %>% column_to_rownames("final_name") %>% dplyr::select(contains(c(".mzML")))
darunavir_filtered[darunavir_filtered < 10000] <- 0
darunavir_heatmap <- darunavir_filtered[,colSums(darunavir_filtered)>0]
darunavir_heatmap_log <- log10(darunavir_heatmap+1)
darunavir_heatmap_log <- as.matrix(darunavir_heatmap_log)

col_fun <- colorRamp2(c(0, 4, 5, 8), c("#336699", "#336699", "#FFFF99", "#FF3333"))
darunavir_heatmpa <- Heatmap(darunavir_heatmap_log, name = "Log_Peak_Area", col = col_fun, 
                             cluster_rows = F, show_column_dend = F, show_column_names = F,
                             row_names_gp = gpar(fontsize = 14),
                             row_gap = unit(5, "mm"), cluster_row_slices = F,
                             show_heatmap_legend = F
)
# heatmap_ggplot <- as.ggplot(grid.grabExpr(draw(darunavir_heatmpa)))
# ggsave(file="darunavir_heatmap.svg", plot=heatmap_ggplot, width=10, height=4)

##############################################
# HNRC people with HIV drug cluster analysis #
##############################################

# HIV drug annotations combined with peak area
exo_drug_cleaned <- read.csv("HNRC_analysis/HNRC_exogenous_drug_manually_cleaned.csv")
feature <- read.csv("HNRC_analysis/HDRC_quant.csv")
feature_filter <- feature %>% column_to_rownames("row.ID") %>% dplyr::select(contains("Sample"))
exo_drug_cleaned$FeatureID <- as.character(exo_drug_cleaned$FeatureID)
hiv_drug_detection <- exo_drug_cleaned %>% inner_join((feature_filter %>% rownames_to_column("FeatureID")), by = "FeatureID") %>%
  dplyr::filter(grepl("HIV", therapeutic_indication) | grepl("atazanavir", therapeutic_indication))
hiv_drug_detection_filtered <- hiv_drug_detection %>% dplyr::select(contains(c("FeatureID", ".mzML"))) %>% column_to_rownames("FeatureID")

hiv_drug_detection_noise_filtered <- hiv_drug_detection_filtered
hiv_drug_detection_noise_filtered[hiv_drug_detection_noise_filtered < 10000] <- 0
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
metadata <- read_xlsx("HNRC_analysis/hnrp_mibi_master_metadata_072424_updated.xlsx")
metadata$SampleID <- sub("11135.", "", metadata$SampleID)
metadata_filtered <- metadata %>% dplyr::select(c("SampleID", "hiv_status_clean", "arv_status"))
hiv_drug_log_filtered <- hiv_drug_log %>% dplyr::select(any_of(metadata_filtered$SampleID))
metadata_filtered <- metadata_filtered %>% dplyr::filter(metadata_filtered$SampleID %in% colnames(hiv_drug_log_filtered)) %>%
  arrange(match(SampleID, colnames(hiv_drug_log_filtered)))
all(metadata_filtered$SampleID == colnames(hiv_drug_log_filtered))

# heatmap for all people (HIV+ and HIV-) and all metadata
metadata_filtered_hm <- metadata_filtered %>% mutate_if(is.character, ~ ifelse(is.na(.), "unknown", .))
col_fun <- colorRamp2(c(0, 4, 5, 8), c("#336699","#336699", "#FFFF99", "#FF3333"))
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
# heatmap_ggplot <- as.ggplot(grid.grabExpr(draw(hiv_all_heatmap)))
# ggsave(file="hiv_all_heatmap.svg", plot=heatmap_ggplot, width=10, height=4)

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
col_fun <- colorRamp2(c(0, 4, 5, 8), c("#336699","#336699", "#FFFF99", "#FF3333"))
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
                         column_split = factor(sample_cluster_df$Cluster, levels = c(1,2,3,4)))  
# heatmap_ggplot <- as.ggplot(grid.grabExpr(draw(split_heatmap)))
# ggsave(file="hiv_drug_clustered.svg", plot=heatmap_ggplot, width=10, height=2.7)

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
  temp <- data_box_long_sig[data_box_long_sig$Feature == a[i],]
  print(a[i])
  print(pairwise.wilcox.test(temp$Value, temp$Cluster, p.adjust.method = 'BH'))
}

# Box plot for features that are significantly different among the drug exposure groups
feature_final <- feature_normalized_final %>% dplyr::select(
  c("Cluster",
    "717_Candidate Histamine-C2:0 (delta mass:42.0104)",
    "894_Candidate Histamine-C2:0 (delta mass:42.0106)",
    "2118_Candidate Histamine-C4:0 (delta mass:70.0418)",
    "4689_Candidate Histidine-C5:0 (delta mass:84.0572)",
    "4813_Candidate Putrescine-C5:0 (delta mass:84.0576)",
    "5534_Candidate Tyrosine-C2:0 (delta mass:42.0107)",
    "18173_Candidate Leucine-C9:4 (delta mass:132.0572)",
    "31692_Candidate Leucine-C14:0 (delta mass:210.1977)",
    "32710_Candidate Leucine-C18:1 (delta mass:264.2453)",
    "5981_Candidate Dopamine-C2:0 (delta mass:42.0105)",
    "4169_Candidate N-acetylputrescine-C2:0 (delta mass:42.0103)",
    "4361_Candidate Lysine-C5:0 (delta mass:84.0573)"
  ))

data_box_long_sig <- tidyr::pivot_longer(feature_final, cols = -Cluster, names_to = "Feature", values_to = "Value")
data_box_long_sig$Cluster <- factor(data_box_long_sig$Cluster, levels = c("1", "2", "3", "4"))

# simplify N-acyl lipid names
data_box_long_sig$Feature <- gsub("717_Candidate Histamine-C2:0 (delta mass:42.0104)", "Histamine-C2:0", data_box_long_sig$Feature, fixed = T)
data_box_long_sig$Feature <- gsub("894_Candidate Histamine-C2:0 (delta mass:42.0106)", "Histamine-C2:0-isomer", data_box_long_sig$Feature, fixed = T)
data_box_long_sig$Feature <- gsub("2118_Candidate Histamine-C4:0 (delta mass:70.0418)", "Histamine-C4:0", data_box_long_sig$Feature, fixed = T)
data_box_long_sig$Feature <- gsub("4689_Candidate Histidine-C5:0 (delta mass:84.0572)", "Histidine-C5:0", data_box_long_sig$Feature, fixed = T)
data_box_long_sig$Feature <- gsub("4813_Candidate Putrescine-C5:0 (delta mass:84.0576)", "Putrescine-C5:0", data_box_long_sig$Feature, fixed = T)
data_box_long_sig$Feature <- gsub("5534_Candidate Tyrosine-C2:0 (delta mass:42.0107)", "Tyrosine-C2:0", data_box_long_sig$Feature, fixed = T)
data_box_long_sig$Feature <- gsub("18173_Candidate Leucine-C9:4 (delta mass:132.0572)", "Leucine-C9:4", data_box_long_sig$Feature, fixed = T)
data_box_long_sig$Feature <- gsub("31692_Candidate Leucine-C14:0 (delta mass:210.1977)", "Leucine-C14:0", data_box_long_sig$Feature, fixed = T)
data_box_long_sig$Feature <- gsub("32710_Candidate Leucine-C18:1 (delta mass:264.2453)", "Leucine-C18:1", data_box_long_sig$Feature, fixed = T)
data_box_long_sig$Feature <- gsub("5981_Candidate Dopamine-C2:0 (delta mass:42.0105)", "Dopamine-C2:0", data_box_long_sig$Feature, fixed = T)
data_box_long_sig$Feature <- gsub("4169_Candidate N-acetylputrescine-C2:0 (delta mass:42.0103)", "N-acetylputrescine-C2:0", data_box_long_sig$Feature, fixed = T)
data_box_long_sig$Feature <- gsub("4361_Candidate Lysine-C5:0 (delta mass:84.0573)", "Lysine-C5:0", data_box_long_sig$Feature, fixed = T)

# Create box plots
my_comparisons <- list( c("1", "4"),  c("2", "4"), c("3", "4"))
data_box_long_sig$Feature <- factor(data_box_long_sig$Feature,
                                    levels = c("Histamine-C2:0",
                                               "Histamine-C2:0-isomer",
                                               "Histamine-C4:0",
                                               "Histidine-C5:0",
                                               "Putrescine-C5:0",
                                               "Tyrosine-C2:0",
                                               "Leucine-C9:4",
                                               "Leucine-C14:0",
                                               "Leucine-C18:1",
                                               "Dopamine-C2:0",
                                               "N-acetylputrescine-C2:0",
                                               "Lysine-C5:0"))
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
# ggsave(file="acyl_lipid_boxplot.svg", plot=acyl_lipids_boxplot, width=10, height=8)

# test N-acyl lipids by ARV_status in clinical metadata
feature_normalized_arv_status <- feature_normalized %>% dplyr::filter((rowSums(.>0)) > 22) %>% t() %>% as.data.frame() %>%
  rownames_to_column("sample") %>% left_join(metadata_hiv_pos_only[,c("SampleID", "arv_status")], by = c("sample" = "SampleID")) %>% dplyr::select(-sample) %>%
  dplyr::filter(!is.na(arv_status))

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
    "4689_Candidate Histidine-C5:0 (delta mass:84.0572)",
    "4813_Candidate Putrescine-C5:0 (delta mass:84.0576)",
    "5534_Candidate Tyrosine-C2:0 (delta mass:42.0107)",
    "18173_Candidate Leucine-C9:4 (delta mass:132.0572)",
    "31692_Candidate Leucine-C14:0 (delta mass:210.1977)",
    "32710_Candidate Leucine-C18:1 (delta mass:264.2453)",
    "5981_Candidate Dopamine-C2:0 (delta mass:42.0105)",
    "4169_Candidate N-acetylputrescine-C2:0 (delta mass:42.0103)",
    "4361_Candidate Lysine-C5:0 (delta mass:84.0573)"
  ))
data_box_long_sig_arv_status <- tidyr::pivot_longer(feature_final_arv_status, cols = -arv_status, names_to = "Feature", values_to = "Value")
data_box_long_sig_arv_status$arv_status <- factor(data_box_long_sig_arv_status$arv_status, levels = c("ARV Naive", "No ARVs", "non-HAART", "HAART"))

data_box_long_sig_arv_status$Feature <- gsub("717_Candidate Histamine-C2:0 (delta mass:42.0104)", "Histamine-C2:0", data_box_long_sig_arv_status$Feature, fixed = T)
data_box_long_sig_arv_status$Feature <- gsub("894_Candidate Histamine-C2:0 (delta mass:42.0106)", "Histamine-C2:0-isomer", data_box_long_sig_arv_status$Feature, fixed = T)
data_box_long_sig_arv_status$Feature <- gsub("2118_Candidate Histamine-C4:0 (delta mass:70.0418)", "Histamine-C4:0", data_box_long_sig_arv_status$Feature, fixed = T)
data_box_long_sig_arv_status$Feature <- gsub("4689_Candidate Histidine-C5:0 (delta mass:84.0572)", "Histidine-C5:0", data_box_long_sig_arv_status$Feature, fixed = T)
data_box_long_sig_arv_status$Feature <- gsub("4813_Candidate Putrescine-C5:0 (delta mass:84.0576)", "Putrescine-C5:0", data_box_long_sig_arv_status$Feature, fixed = T)
data_box_long_sig_arv_status$Feature <- gsub("5534_Candidate Tyrosine-C2:0 (delta mass:42.0107)", "Tyrosine-C2:0", data_box_long_sig_arv_status$Feature, fixed = T)
data_box_long_sig_arv_status$Feature <- gsub("18173_Candidate Leucine-C9:4 (delta mass:132.0572)", "Leucine-C9:4", data_box_long_sig_arv_status$Feature, fixed = T)
data_box_long_sig_arv_status$Feature <- gsub("31692_Candidate Leucine-C14:0 (delta mass:210.1977)", "Leucine-C14:0", data_box_long_sig_arv_status$Feature, fixed = T)
data_box_long_sig_arv_status$Feature <- gsub("32710_Candidate Leucine-C18:1 (delta mass:264.2453)", "Leucine-C18:1", data_box_long_sig_arv_status$Feature, fixed = T)
data_box_long_sig_arv_status$Feature <- gsub("5981_Candidate Dopamine-C2:0 (delta mass:42.0105)", "Dopamine-C2:0", data_box_long_sig_arv_status$Feature, fixed = T)
data_box_long_sig_arv_status$Feature <- gsub("4169_Candidate N-acetylputrescine-C2:0 (delta mass:42.0103)", "N-acetylputrescine-C2:0", data_box_long_sig_arv_status$Feature, fixed = T)
data_box_long_sig_arv_status$Feature <- gsub("4361_Candidate Lysine-C5:0 (delta mass:84.0573)", "Lysine-C5:0", data_box_long_sig_arv_status$Feature, fixed = T)

# Create box plots
data_box_long_sig_arv_status$Feature <- factor(data_box_long_sig_arv_status$Feature,
                                               levels = c("Histamine-C2:0",
                                                          "Histamine-C2:0-isomer",
                                                          "Histamine-C4:0",
                                                          "Histidine-C5:0",
                                                          "Putrescine-C5:0",
                                                          "Tyrosine-C2:0",
                                                          "Leucine-C9:4",
                                                          "Leucine-C14:0",
                                                          "Leucine-C18:1",
                                                          "Dopamine-C2:0",
                                                          "N-acetylputrescine-C2:0",
                                                          "Lysine-C5:0"))
acyl_lipids_boxplot_arv_status <- ggplot(data_box_long_sig_arv_status, aes(x = arv_status, y = Value, fill = arv_status)) +
  geom_boxplot(outlier.shape=NA, alpha = 0.6, size = 1) +
  geom_point(position = position_jitterdodge(jitter.width = 0.8), size = 2, shape = 20, aes(color = arv_status), alpha = 0.4) +
  facet_wrap(~ Feature, scales = "free_y", ncol = 6) + 
  scale_color_manual(values = c("#922091","darkred", "#ED708D","#336699")) +
  scale_fill_manual(values = c("#922091","darkred", "#ED708D","#336699")) +
  theme_classic() +
  theme(axis.ticks.length=unit(.15, "cm"),
        axis.line.x.bottom=element_line(linewidth=1.0), axis.line.y.left = element_line(linewidth=1.0),
        axis.ticks.x.bottom = element_blank(), axis.ticks.y.left = element_line(linewidth=1.0),
        axis.text.y.left = element_blank(), axis.text.x.bottom = element_blank(),
        axis.title.y.left = element_blank(), axis.title.x.bottom = element_blank(),
        strip.background = element_blank(), strip.text = element_text(size = 12), face = "bold")
# ggsave(file="acyl_lipid_boxplot_arv.svg", plot=acyl_lipids_boxplot_arv_status, width=10, height=4)

geom_point(position = position_jitterdodge(jitter.width = 0.8), size = 4, shape = 20, aes(color = Cluster), alpha = 0.4) +
  
#############################################
# microbial culture boxplot of drug analogs #
#############################################
feature_microbe <- read.csv("HNRC_analysis/SynCom_quant.csv")
feature_filter <- feature_microbe %>% column_to_rownames("row.ID") %>% dplyr::select(contains("mzML"))
colnames(feature_filter) <- sub(".mzML.Peak.area","", colnames(feature_filter))
metadata <- read.csv("HNRC_analysis/SynCom_metadata.csv")
feature_filter_analogs <- feature_filter %>% rownames_to_column("FeatureID") %>%
  dplyr::filter(FeatureID %in% c("12995", "13875", "6222", "7256", "11178")) %>%
  column_to_rownames("FeatureID") %>% t() %>% as.data.frame() %>% rownames_to_column("sample") %>%
  left_join(metadata, by = "sample")
colnames(feature_filter_analogs)[2:6] <- c("F6222", "F7256", "F11178", "F13875", "F12995")

all_analogs <- data.frame(analog = c("F12995", "F13875", "F6222", "F7256", "F11178"), drug_group = c(2, 3, 4, 4, 3))
boxplot <- data.frame(final_peak_area = NULL, final_group = NULL, feature = NULL)

# blank subtraction with BHI media only features
for (i in 1:5){
  if (i<5){
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

# List of y-axis limits for each feature
y_axis_limits <- list(
  "F11178" = c(0, 3000000),
  "F12995" = c(0, 30000),
  "F13875" = c(0, 40000),
  "F6222" = c(0, 300000),
  "F7256" = c(0, 300000)
  # Add more features and their corresponding limits here
)

# Create a list to store plots
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

# Combine the plots
combined_plot_SynCom <- wrap_plots(plots, ncol = 2)
ggsave("combined_plot_SynCom.svg", plot = combined_plot_SynCom, width = 12, height = 12)

# Wilcox test to compare T0 and T72
boxplot_final <- boxplot_final %>% mutate(microbe_group = str_split(final_group, "_", simplify = TRUE)[, 1])
test_feature <- unique(boxplot_final$feature)
test_group <- unique(boxplot_final$microbe_group)
for(i in 1:5){
  for (j in 1:6){
    temp <- boxplot_final[boxplot_final$feature == test_feature[i] & boxplot_final$microbe_group == test_group[j],]
    print(test_feature[i])
    print(test_group[j])
    print(wilcox.test(final_peak_area ~ final_group,
                      data = temp, exact = TRUE))
  }
}