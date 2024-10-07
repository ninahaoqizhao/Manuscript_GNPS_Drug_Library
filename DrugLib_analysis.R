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


############################################
# Bar plot for analog library mass offsets #
############################################

analog_metadata <- read.csv("druglib_metadata/GNPS_Drug_Library_drug_analog_metadata.csv")
extract_unique_delta_mass <- function(input_string) {
  delta_mass_values <- str_extract_all(input_string, "(?<=delta Mass:)-?\\d+\\.\\d+")
  all_values <- c(unlist(delta_mass_values))
  all_values <- as.numeric(all_values)
  all_values <- formatC(round(all_values, 2), format='f', digits=2)
  unique_values <- unique(all_values)
  cooccurrence_by_delta <- paste(unique_values, collapse = "|")
  return(cooccurrence_by_delta)
}
negate_values <- function(x) {
  values <- strsplit(x, "\\|")[[1]]
  negated_values <- as.numeric(values) * -1
  negated_string <- paste(negated_values, collapse = "|")
  return(negated_string)
}
cleaned_delta$cleaned_delta_mass <- sapply(cleaned_delta$cleaned_delta_mass, negate_values)
cleaned_delta_separate <- tidyr::separate_rows(cleaned_delta, cleaned_delta_mass, sep = "\\|")
cooccurrence_by_delta <- as.data.frame(table(cleaned_delta_separate$cleaned_delta_mass))
colnames(cooccurrence_by_delta) <- c("delta_mass", "freq")
cooccurrence_by_delta$masst_detection <- rep(NA, nrow(cooccurrence_by_delta))
for(i in cooccurrence_by_delta$delta_mass){
  temp <- cleaned_delta_separate %>% dplyr::filter(cleaned_delta_separate$cleaned_delta_mass == i)
  cooccurrence_by_delta$masst_detection[cooccurrence_by_delta$delta_mass == i] <- sum(temp$masst_detection_ratio > 0)
}
cooccurrence_by_delta$detection_freq <- cooccurrence_by_delta$masst_detection/cooccurrence_by_delta$freq
write.csv(cooccurrence_by_delta, "application/drug_in_masst/20240711_analog_masst_detection_freq_by_delta_mass.csv")

# delta mass frequency visualization
cooccurrence_by_delta$atom <- rep("unknown",nrow(cooccurrence_by_delta))

#explain: freq vs delta mass plot
cooccurrence_by_delta$atom[cooccurrence_by_delta$delta_mass == -15.99 | cooccurrence_by_delta$delta_mass == 15.99] <- "1O" 
cooccurrence_by_delta$atom[cooccurrence_by_delta$delta_mass == 14.02 | cooccurrence_by_delta$delta_mass == -14.02] <- "1C, 2H" 
cooccurrence_by_delta$atom[cooccurrence_by_delta$delta_mass == 1] <- "C isotope" 
cooccurrence_by_delta$atom[cooccurrence_by_delta$delta_mass == 2] <- "O/S/Cl/Br isotope" 
cooccurrence_by_delta$atom[cooccurrence_by_delta$delta_mass == 28.03 | cooccurrence_by_delta$delta_mass == -28.03] <- "2C, 4H" 
cooccurrence_by_delta$atom[cooccurrence_by_delta$delta_mass == 18.01 | cooccurrence_by_delta$delta_mass == -18.01] <- "1O, 2H" 
cooccurrence_by_delta$atom[cooccurrence_by_delta$delta_mass == 27.99 | cooccurrence_by_delta$delta_mass == -27.99] <- "1C, 1O" 
cooccurrence_by_delta$atom[cooccurrence_by_delta$delta_mass == 2.02 | cooccurrence_by_delta$delta_mass == -2.02] <- "2H" 



cooccurrence_by_delta$atom <- factor(cooccurrence_by_delta$atom, levels = c("1O", "1C, 2H", "C isotope", "O/S/Cl/Br isotope", "2C, 4H", "1O, 2H","1C, 1O","2H","unknown"))
cooccurrence_by_delta$delta_mass <- as.numeric(as.character(cooccurrence_by_delta$delta_mass))
class(cooccurrence_by_delta$delta_mass)

plot <- ggplot(cooccurrence_by_delta, aes(x=delta_mass, weight=freq, fill = atom)) +
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

ggsave(plot = plot, file = "suspect_delta_mass_5.svg", width = 10, height = 6)

# test unique delta in different drug groups
cleaned_delta_filtered <- cleaned_delta %>% dplyr::filter(cleaned_delta$Kine_FDA_drug_class != "no_match" &
                                                            !grepl("QAQC|Low", cleaned_delta$chemical_source))
cooccurrence_by_delta <- cleaned_delta_filtered %>%
  group_by(Kine_FDA_drug_class, cleaned_delta_mass) %>%
  summarise(count = n()) %>%
  ungroup()



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

# Code to run the GNPS Drug Library with fastMASST in batch: https://github.com/robinschmid/microbe_masst/tree/v1.2.0
# Output for each spectrum is a tsv files containing all the USI matched to the spectrum
# All the tsv file has been pasted together into a single csv provided with this code

# masst results for all parent drugs
drug_masst <- read.csv("masst_results_entire_library/druglib_parent_drug_fasst_matches.csv")
drug_masst_filtered <- drug_masst %>% dplyr::select(colnames(drug_masst)[-1]) %>%
  dplyr::filter(!grepl("MSV000084314|MSV000089210", drug_masst$USI))
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
# fasst search was performed with older version of the analog library, which was re-uploaded to GNPS aftre filtering to generate the final ccmslibid 
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
