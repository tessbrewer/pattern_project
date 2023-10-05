library(ggtree)
library(ape)
library(stringr)
library(stringi)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(phytools)
library(scales)
library(ggnewscale)
library(plyr)
library(ggpattern)
library(ggseqlogo)
library(seqinr)
library(cowplot)
library(ggrepel)
library(gggenes)

#Pattern project functions =================================
setwd('~/Documents/pattern_project/data/') #wherever you saved the data folder

#used to load output data from python scripts
load_data <- function(path, file_pattern) { 
  files <- dir(path, file_pattern, full.names=TRUE)
  slashes <- str_count(files[1], '/')
  sample <- sapply(strsplit(as.character(files), '/'), '[', slashes+1)
  file_extension <- sapply(strsplit(as.character(sample), '\\.'), '[', str_count(sample[1], '\\.') + 1)
  sample <- gsub(gsub('*', '', file_pattern), '', sample)
  if (file_extension[1] == 'csv'){
    tables <- lapply(files, read.csv, header=TRUE)
  } else {
    tables <- lapply(files, read.delim, header=TRUE)
  }
  tables <- Map(cbind, tables, family = sample)
  df <- do.call(rbind, tables)
  df$family <- gsub('_PP', '', df$family)
  return(df)
}

#used to load interpro data 
load_interpro <- function(path, file_pattern) { 
  files <- dir(path, file_pattern, full.names=TRUE)
  slashes <- str_count(files[1], '/')
  sample <- sapply(strsplit(as.character(files), '/'), '[', slashes+1)
  sample <- sapply(strsplit(as.character(sample), '_'), '[', 1)
  file_extension <- sapply(strsplit(as.character(sample), '\\.'), '[', str_count(sample[1], '\\.') + 1)
  tables <- lapply(files, read.delim, header=FALSE)
  tables <- Map(cbind, tables, protein = sample)
  df <- do.call(rbind, tables)
  colnames(df) <- c("protein_id", "protein_md5", "protein_length", "analysis",
                    "signature_id", "signature_desc", 'pfam_desc', "start", "end", "score",
                    "status", "date", "interpro_id", "interpro_desc", "go_terms", 'cyc_terms', 'protein')
  df$cyc_terms <- NULL
  return(df)
}

#formats data how phylANOVA wants it
phylogenetic_anova <- function(tree, groups_to_test, valuez){
  test_values <- as.vector(valuez)
  names(test_values) <- names(groups_to_test) #MUST BE IN SAME ORDER
  anova_results <- phylANOVA(tree, groups_to_test, test_values, posthoc = FALSE)
  return(anova_results)
}

#creates data frame of stalls / gene presence of every silix family
produce_nonidential_stalls <- function(data, file_pattern, sequence_pattern){
  conserved_stalls <- load_data(data, file_pattern)
  if (TRUE %in% names(table(grepl('Planctopirus', conserved_stalls$gene)))){
    plancto = TRUE
    colnames(conserved_stalls) <- c('x', 'gene', 'genome', 'conserved_stall', 'position', 'family')
  } else {
    plancto = FALSE
    colnames(conserved_stalls) <- c('x', 'genome', 'gene', 'conserved_stall', 'position', 'family')
  } #files from planctos used different column names :(
  conserved_stalls$stalled <- ifelse(grepl(sequence_pattern, conserved_stalls$conserved_stall), 2, 1)
  conserved_stalls <- conserved_stalls[order(conserved_stalls$stalled), ]
  # conserved_stalls <- conserved_stalls[order(conserved_stalls$stalled, decreasing = TRUE), ]
  #to be counted the motif must not be present in any copy! if they have multiple copies and one has the motif, one 
  #doesn't, the one with the motif is used
  conserved_table <- dcast(conserved_stalls[!duplicated(paste(conserved_stalls$genome, conserved_stalls$family)), ], genome ~ family, value.var = 'stalled')
  #only one family per genome is kept, with the order deciding -> first entry is kept
  row.names(conserved_table) <- conserved_table$genome
  conserved_table$genome <- NULL
  conserved_table[is.na(conserved_table)] <- 0 #does not contain a gene of the family
  nonidentical_stalls <- Filter(function(x) (length(unique(x)) > 1), conserved_table)
  #filters out columns (families) that have no variation (all stalled, unstalled, or gene absent)
  return(nonidentical_stalls)
}

#checks for significant co-occurrence of sequence patterns with specified phylogenetic group
run_tests <- function(data, file_pattern, sequence_pattern, group, group_name){
  #load data
  print(paste('Running program with', sequence_pattern, 'and', group_name))
  nonidentical_stalls <- produce_nonidential_stalls(data, file_pattern, sequence_pattern)
  #check for plancto and get groups to test for 
  plotted_map$testing_groups <- ifelse(plotted_map$taxon_oid %in% group, 'yes', 'no')
  if (TRUE %in% names(table(grepl('Planctopirus', row.names(nonidentical_stalls))))){
    plancto = TRUE
    nonidentical_stalls$hgt <- plotted_map$testing_groups[match(row.names(nonidentical_stalls), plotted_map$species_code)]
  } else {
    plancto = FALSE
    nonidentical_stalls$hgt <- plotted_map$testing_groups[match(row.names(nonidentical_stalls), plotted_map$taxon_oid)]
  }
  #sub-setting to criteria in methods (keeps p-value correction from being too stringent)
  check_havegene <- aggregate(. ~ hgt, nonidentical_stalls, function(x){length(which(x != 0))})
  check_havepattern <- aggregate(. ~ hgt, nonidentical_stalls, function(x){length(which(x == 2))})
  check <- rbind(check_havegene, check_havepattern)
  check$hgt <- c('no_have_gene', 'yes_have_gene', 'no_have_pp', 'yes_have_pp')
  check <- dcast(melt(check, id.var = 'hgt'), variable ~ hgt, value.var = 'value', id.var = 'hgt')
  check$no <- as.numeric(check$no_have_pp)/as.numeric(check$no_have_gene) # non-HGT groups ratio w/ pp motif
  check$yes <- as.numeric(check$yes_have_pp)/as.numeric(check$yes_have_gene)  # HGT groups ratio w/ pp motif
  check$no_withgene <- check$no_have_gene / (nrow(subset(plotted_map, testing_groups == 'no'))) # non-HGT groups ratio w/ GENE
  check$yes_withgene <- check$yes_have_gene / (nrow(subset(plotted_map, testing_groups == 'yes'))) # HGT groups ratio w/ GENE
  check$ko <- annot$Var2[match(check$variable, annot$Var1)]
  check$description <- kegg_key$description[match(check$ko, kegg_key$ko_number)]
  to_check <- as.character(subset(check, no > 0.85 & no_withgene > 0.70 & yes <= 0.20 | #motif is poorly conserved in HGT group
                                         no > 0.85 & no_withgene > 0.70 & is.nan(yes) | #protein is absent from HGT group
                                         no > 0.85 & no_withgene > 0.70 & yes_withgene <= 0.20)$variable) #protein is rare in HGT group
  #okay now to testing
  testing_table <- nonidentical_stalls[, colnames(nonidentical_stalls) %in% to_check]
  #groupings
  if (plancto == FALSE) {
    groups_to_test <- as.vector(subset(plotted_map, taxon_oid %in% row.names(testing_table))$testing_groups)
    names(groups_to_test) <- subset(plotted_map, taxon_oid %in% row.names(testing_table))$taxon_oid
  } else {
    groups_to_test <- as.vector(subset(plotted_map, species_code %in% row.names(testing_table))$testing_groups)
    names(groups_to_test) <- subset(plotted_map, species_code %in% row.names(testing_table))$species_code
  }
  #order between groups and values must be the same, function assumes it!
  groups_to_test <- groups_to_test[order(names(groups_to_test))]
  testing_table <- testing_table[order(row.names(testing_table)),]
  if (plancto == FALSE) {
    testing_table <- subset(testing_table, row.names(testing_table) %in% plotted_map$taxon_oid)
  } else {
    testing_table <- subset(testing_table, row.names(testing_table) %in% plotted_map$species_code)
  }
  table(row.names(testing_table) == names(groups_to_test))
  if (plancto == FALSE) {
    names(groups_to_test) <- plotted_map$species_code[match(names(groups_to_test), plotted_map$taxon_oid)]
    row.names(testing_table) <- names(groups_to_test)
    #tree has species_code as ID, it's a bit of a mess, yeah
  }
  table(row.names(testing_table) == names(groups_to_test))
  test_tree <- drop.tip(phy = plotted_tree,
                        tip = plotted_tree$tip.label[!plotted_tree$tip.label %in%
                                                       names(groups_to_test)])
  #run test and collect results
  raw_results <- apply(testing_table, 2, function(x){phylogenetic_anova(test_tree, groups_to_test, x)})
  total_results <- data.frame(matrix(unlist(raw_results), nrow = length(raw_results), byrow = TRUE))[, 1:2]
  colnames(total_results) <- c('F', 'Pf') #first two are important
  row.names(total_results) <- names(raw_results)
  total_results$p_adjust <- round(p.adjust(total_results$Pf, method = 'fdr'), 2)
  #add annotations
  total_results$ko <- annot$Var2[match(row.names(total_results), annot$Var1)]
  total_results$c_level <- kegg_key$c_level[match(total_results$ko, kegg_key$ko_number)]
  total_results$d_level <- kegg_key$description[match(total_results$ko, kegg_key$ko_number)]
  total_results$gene <- ifelse(is.na(total_results$d_level),
                               clusters$family[match(row.names(total_results), clusters$family)],
                               kegg_key$short_description[match(total_results$ko, kegg_key$ko_number)])
  total_results <- total_results[order(total_results$p_adjust),]
  total_results$sequence <- sequence_pattern
  total_results$group <- group_name
  print(paste('Total significant hits', nrow(total_results[total_results$p_adjust <= 0.1,])))
  cat(c(''), sep="\n\n")
  return(total_results)
}

#specific notation for ggrepel, making multicolored plot labels
phantom <- function(x) {
  paste0('phantom("', x, '")')
}
#specific notation for ggrepel, making multicolored plot labels
quotes <- function(x) {
  paste0('"', x, '"')
}

#removes redundant pairs left behind after cophenetic()
remove_pairs <- function(practice, first, second){
  practice$first_pair <- paste(practice[[first]], practice[[second]], sep = '_')
  practice$second_pair <- paste(practice[[second]], practice[[first]], sep = '_')
  for (i in 1:nrow(practice)){
    location <- grep(practice$first_pair[i], practice$second_pair)
    if (is.na(location[1])){
      next
    } else {
      practice <- practice[-location, ]
    }
  }
  practice$second_pair <- NULL
  return(practice)
}

#drawing a polygon around points of interest
find_hull <- function(x) {
  x[chull(x$phylogenetic_distance, x$efp_distance), ]
}


#Pattern project commonly used files =================================

#key to kegg annotation descriptions
kegg_key <- read.delim("kegg_key.txt", quote="", fill=FALSE)
kegg_key$short_description <- sapply(strsplit(as.character(kegg_key$description), ';'), '[', 1)
kegg_key$short_description<- trimws(ifelse(grepl(',', kegg_key$short_description), 
                                           sapply(strsplit(as.character(kegg_key$short_description), ','), '[', 1), 
                                           kegg_key$short_description))

#kegg doesn't annotate earP, so using the pFAM annotation
earP <- read.csv("earP_pfam10093.csv", header = FALSE,
                 col.names = c('genome', 'pfam', 'gene', 'description', 'evalue', 'bitscore')) 
earP <- subset(earP, gene != "none")

#mapping file for all genomes
mapping <- read.delim("final_mapping.txt")

#efp mapping file for all genomes
efp_mapping <- read.delim("efp_mapping.txt")


