## Planctomycetes specific analyses for pattern project paper
##
## Purpose: Do analyses for the plancto genomes and make figures 4, S3, S4, and S6, 
##
## Author: Dr. Tess Brewer ;)

# Load functions and color palettes =================================
source('scripts/pattern_project_functions.R') #make sure to change working directory in this script first

#palette for the figures
colorz <- c("#002954", "#2373A2", "#4DCCFF", 
            "#118038", "#01ff70",
            "#fcff4e", "#FDA257", "#FF0066")

# File set-up =================================

plotted_tree <- read.tree("planctomycetes/planctomycetes_species.tre")
plotted_map <- subset(mapping, species_code %in% plotted_tree$tip.label & phylum_gtdbtk == 'Planctomycetota')
plotted_tree <- drop.tip(plotted_tree, tip = plotted_tree$tip.label[!plotted_tree$tip.label %in% plotted_map$species_code])
row.names(plotted_map) <- plotted_map$species_code
plotted_map <- plotted_map[order(row.names(plotted_map)),]

#nodes to place the green plus and red x on
HGT = subset(plotted_map, grepl('EpmAB', efp_modification))$species_code
HGT_transfer2 <- getMRCA(plotted_tree, subset(plotted_map, grepl('EarP', efp_modification))$species_code)
HGT_transfer1 <- getMRCA(plotted_tree, subset(plotted_map, grepl('EpmAB', efp_modification))$species_code)
ancestral_loss1 <- getMRCA(plotted_tree, subset(plotted_map, grepl('EpmAB', efp_modification))$species_code)
ancestral_loss2 <- getMRCA(plotted_tree, subset(plotted_map, grepl('EarP', efp_modification))$species_code)

plotted_map$hgt_efp <- ifelse(plotted_map$species_code %in% HGT, 'yes', 'no')
plotted_map$color <- ifelse(plotted_map$class_gtdbtk == 'Planctomycetes', plotted_map$family_gtdbtk, plotted_map$class_gtdbtk)

#base tree
node_interest <- findMRCA(plotted_tree, tips = as.character(subset(plotted_map, hgt_efp == 'yes')$species_code))
#class 55
full_tree <- ggtree(plotted_tree, open.angle=12, color = "#666666") %<+% plotted_map +
  ggtitle('Species tree') +
  geom_tiplab(aes(label = species), size = 4, align=F, family="Avenir", offset = 0.03) +
  geom_nodepoint(aes(subset=(node == ancestral_loss1)), 
                 shape=4, size=5, color='#D74840', stroke = 2) +
  geom_nodepoint(aes(subset=(node == ancestral_loss2)), 
                 shape=4, size=5, color='#D74840', stroke = 2, position = position_nudge(x = -0.06)) +
  geom_nodepoint(aes(subset=(node == HGT_transfer1)), 
                 shape=3, size=5, color='#008B52', stroke = 2, position = position_nudge(x = -0.06)) +
  geom_nodepoint(aes(subset=(node == HGT_transfer2)), 
                 shape=3, size=5, color='#008B52', stroke = 2, position = position_nudge(x = -0.12)) +
  scale_fill_manual(values = rev(colorz), name = 'Family (GTDB-TK)', 
                    breaks = c("Isosphaeraceae", "Pirellulaceae", "Planctomycetaceae", "Thermoguttaceae", 
                               "Unclassified", "Brocadiae", "Phycisphaerae", "Verrucomicrobiae")) +
  scale_size_continuous(name = 'Genome size (MB)', range = c(3, 7)) +
  geom_tippoint(aes(fill = color, size = genome_size/1000000), pch = 21, color = "#666666") +
  xlim_tree(1.75) +
  geom_treescale() +
  theme(text = element_text(family="Avenir", size=16, color = "black")) +
  guides(fill = guide_legend(override.aes = list(size=5)))


# Making figure s6 =================================

efp_mapping <- subset(efp_mapping, taxon_oid %in% plotted_map$taxon_oid)

#position 2 and 11 are "-" in all planctomycetes sequences, toss them (because they were aligned with non-planctos)
#just doing it for aesthetic reasons really
first_part <- paste(substring(efp_mapping$spaced_loop, 1, 1), substring(efp_mapping$spaced_loop, 3, 6), sep = '')
efp_mapping$final_loop <- paste(first_part, substring(efp_mapping$spaced_loop, 8, 10), sep = '')
efp_mapping$final_loop <- gsub('TPSA-RGA', 'T-PSARGA', efp_mapping$final_loop)
table(str_count(efp_mapping$final_loop)) #check they're all the same length

#for this figure need to differentiate between the loop sequence and modifications of 'efp 1' and 'efp 2'
plotted_map$actino_loop_prime <- subset(efp_mapping, modification != 'Unmodified YeiP')$final_loop[match(plotted_map$taxon_oid, subset(efp_mapping, modification != 'Unmodified YeiP')$taxon_oid)]
plotted_map$actino_loop_second <- ifelse(plotted_map$taxon_oid %in% subset(efp_mapping, modification == 'Unmodified YeiP')$taxon_oid,
                                         subset(efp_mapping, modification == 'Unmodified YeiP')$final_loop[match(plotted_map$taxon_oid, subset(efp_mapping, modification == 'Unmodified YeiP')$taxon_oid)]
                                         , '')
plotted_map$modification_prime <- subset(efp_mapping, modification != 'Unmodified YeiP')$modification[match(plotted_map$taxon_oid, subset(efp_mapping, modification != 'Unmodified YeiP')$taxon_oid)]
plotted_map$modification_second <- ifelse(plotted_map$efp_count > 1, 'YeiP', '')

#labelling stuff -- this is all specific to ggrepel (adding underline of conserved residue)
#complicated becasue there are two efp for many planctos
plotted_map$label_most1 <- paste(quotes(substr(plotted_map$actino_loop_prime, 0, 5)), 
                                 phantom('A'), 
                                 quotes(substr(plotted_map$actino_loop_prime, 7, 8)), sep = " * ")
plotted_map$label_residue1 <- paste(phantom(substr(plotted_map$actino_loop_prime, 0, 5)), 
                                    quotes(substr(plotted_map$actino_loop_prime, 6, 6)), 
                                    phantom(substr(plotted_map$actino_loop_prime, 7, 8)), sep = " * ")

#secondary 
plotted_map$label_most2 <- paste(quotes(substr(plotted_map$actino_loop_second, 0, 5)), 
                                 phantom('A'), 
                                 quotes(substr(plotted_map$actino_loop_second, 7, 8)), sep = " * ")
plotted_map$label_residue2 <- paste(phantom(substr(plotted_map$actino_loop_second, 0, 5)), 
                                    quotes(substr(plotted_map$actino_loop_second, 6, 6)), 
                                    phantom(substr(plotted_map$actino_loop_second, 7, 8)), sep = " * ")

#need to reattach plotted_map
#add the pluses and minuses next time
full_tree <- ggtree(plotted_tree, open.angle=12, color = "#666666") %<+% plotted_map +
  ggtitle('Species tree') +
  geom_hilight(node = ancestral_loss1, extend = 1.5, alpha = 0.2) +
  geom_tiplab(aes(label = species), size = 4, align=F, family="Avenir", offset = 0.03) +
  geom_nodepoint(aes(subset=(node == ancestral_loss1)), 
                 shape=4, size=5, color='#D74840', stroke = 2) +
  geom_nodepoint(aes(subset=(node == ancestral_loss2)), 
                 shape=4, size=5, color='#D74840', stroke = 2, position = position_nudge(x = -0.06)) +
  geom_nodepoint(aes(subset=(node == HGT_transfer1)), 
                 shape=3, size=5, color='#008B52', stroke = 2, position = position_nudge(x = -0.06)) +
  geom_nodepoint(aes(subset=(node == HGT_transfer2)), 
                 shape=3, size=5, color='#008B52', stroke = 2, position = position_nudge(x = -0.12)) +
  scale_fill_manual(values = rev(colorz), name = 'Family (GTDB-Tk)', 
                    breaks = c("Isosphaeraceae", "Pirellulaceae", "Planctomycetaceae", "Thermoguttaceae", 
                               "Unclassified", "Brocadiae", "Phycisphaerae", "Verrucomicrobiae")) +
  scale_size_continuous(name = 'Genome size (MB)', range = c(3, 7)) +
  geom_tippoint(aes(fill = color, size = genome_size/1000000), pch = 21, color = "#666666") +
  xlim_tree(1.75) +
  geom_treescale() +
  theme(text = element_text(family="Avenir", size=16, color = "black")) +
  guides(fill = guide_legend(override.aes = list(size=5)))

#efp detail tree
trial <- plotted_map[, c("modification_prime", "modification_second")]
new_plot <- full_tree + new_scale_fill()

figure_s6 <- gheatmap(new_plot, trial, width=0.18, offset = 0.6, hjust = 0, color=NA,
                           colnames_position = 'top', colnames_offset_x = 0, colnames = FALSE, font.size = 4, 
                           colnames_offset_y = -0.1, colnames_angle = 45, family = 'Avenir') +
  xlim_tree(1.1) +
  geom_text_repel(aes(x = 1.95, label = label_residue1), color = "red", parse = TRUE, seed = 1, hjust = 'left', family = "mono", size = 4) +
  geom_text_repel(aes(x = 1.95, label = label_most1), color = "black", parse = TRUE, seed = 1, hjust = 'left', family = "mono", size = 4) +
  geom_text_repel(aes(x = 2.15, label = label_residue2), color = "red", parse = TRUE, seed = 1, hjust = 'left', family = "mono", size = 4) +
  geom_text_repel(aes(x = 2.15, label = label_most2), color = "black", parse = TRUE, seed = 1, hjust = 'left', family = "mono", size = 4) +
  scale_fill_manual(values = c('#BC80BD', '#8dde9f', '#ededed', '#FB8072'), name = 'EFP modification', na.value="white") +
  theme(text = element_text(family="Avenir", size=18, color = "black"))

# ggsave("figure_s6.svg", plot = figure_s6, width = 14, height = 12)

# Load annotations =================================

#kegg annotations of genes in planctomycetes genomes
ko_annotations <- read.csv("planctomycetes/Planctomycetota_gathered_kos.csv")
ko_annotations <- ko_annotations[order(ko_annotations$bit_score, decreasing = TRUE), ]
ko_annotations <- ko_annotations[!duplicated(ko_annotations$gene_oid), ]
ko_annotations$ko_id <- gsub('KO:K', 'K', ko_annotations$ko_id)

#silix clusters
clusters <- read.delim('planctomycetes/planctomycetes_silix.fnodes',header = FALSE, col.names = c('family', 'genome'))
clusters$gene <- sapply(strsplit(as.character(clusters$genome), '_'), '[', 2)
clusters$taxon_oid <- sapply(strsplit(as.character(clusters$genome), '\\.'), '[', 3)
clusters$taxon_oid <- sapply(strsplit(as.character(clusters$taxon_oid), '_'), '[', 1)
clusters$ko <- ko_annotations$ko_id[match(paste(clusters$taxon_oid, clusters$gene), paste(ko_annotations$genome, ko_annotations$gene_oid))]
clusters$ko[is.na(clusters$ko)] <- 'unassigned'
annot <- as.data.frame(table(clusters$family, clusters$ko))
annot <- annot[order(annot$Freq, decreasing = TRUE), ]
annot <- annot[!duplicated(annot$Var1), ]
annot$description <- ifelse(annot$Var2 == 'unassigned', as.character(annot$Var1), 
                            kegg_key$short_description[match(annot$Var2, kegg_key$ko_number)])


# Running tests =================================
run_data <- list(data = c(rep('planctomycetes/pp_search/', 2)),
                 file_pattern = c(rep('_conserved_stall.csv', 2)),
                 sequence_pattern = c(rep('..PP.', 2)),
                 group_name = c('HGT', 'Planctomycetaceae_only'),
                 group = c(list(subset(plotted_map, hgt_efp == 'yes')$taxon_oid),
                           list(subset(plotted_map, family_gtdbtk == 'Planctomycetaceae')$taxon_oid)))

#just run the function in a loop
#some element of randomness to this test, sometimes comes up with different levels of significance
all_results <- data.frame()
for (i in 1:length(run_data[[1]])){
  results <- run_tests(run_data$data[i], run_data$file_pattern[i], 
                       run_data$sequence_pattern[i], run_data$group[[i]], run_data$group_name[i])
  results$family <- row.names(results)
  row.names(results) <- NULL #forcing duplicate row.names otherwise and R will rename adding 1,2... to end
  all_results <- rbind(all_results, results)
  rm(results)
}

results_final <- subset(all_results, p_adjust <= 0.1)


#to look directly at the table for particular family of proteins, for example lon
stall_data <- load_data(run_data$data, run_data$file_pattern)
kegg_key[grep('^lon', kegg_key$description),] # get kegg ko number for protein of interest
annot[annot$Var2 == 'K01870',] # get silix family for that ko, some kos have multiple families
look <- subset(stall_data, family == 'FAM003801') # #extract from stall data


# Take a look at hits / make figure 4 =================================

#get pp motif presence / absence / gene absence info
nonidentical_stalls <- produce_nonidential_stalls(run_data$data[1], run_data$file_pattern[1], run_data$sequence_pattern[1])
fam_to_plot <- nonidentical_stalls[, c(colnames(nonidentical_stalls) %in% results_final$family |
                                       colnames(nonidentical_stalls) == 'FAM002091')] #ileS is interesting to include for figure
fam_to_plot$hgt <- plotted_map$hgt_efp[match(row.names(fam_to_plot), plotted_map$species_code)]

#automate ordering so families with similar patterns are grouped together and everything looks nice
orderingz <- subset(melt(aggregate(. ~ hgt, fam_to_plot, mean)), hgt == 'yes')
orderingz <- droplevels(orderingz$variable[order(orderingz$value, decreasing = TRUE)])
fam_to_plot$hgt <- NULL #need to get rid of it for future dataframe

#convert pp motif presence to categorical and label correctly
plotting_data <- as.data.frame(sapply(fam_to_plot, function(x){str_replace_all(x, c('0' = 'gene absent', '1' = 'absent', '2' = 'present'))}))
row.names(plotting_data) <- row.names(fam_to_plot)
plotting_data <- plotting_data[, orderingz]
annot$description[annot$Var1 == 'FAM002091'] <- 'IleS2'
annot$description[annot$Var1 == 'FAM003949'] <- 'IleS1'
annot$description[annot$Var1 == 'FAM003801'] <- 'ValS'
annot$description[annot$Var1 == 'FAM004963'] <- 'Lon'


colnames(plotting_data) <- paste(annot$description[match(colnames(plotting_data), annot$Var1)])

new_plot <- full_tree + new_scale_fill()
figure_4 <- gheatmap(new_plot, plotting_data, width=0.35, offset = 0.5, color = "white", hjust = 0,
                     colnames_position = 'top', colnames_offset_x = 0, font.size = 5,
                     colnames_offset_y = -0.1, colnames_angle = 45, family = 'Avenir') +
  xlim_tree(2) +
  # ggplot2::ylim(0, 25) +
  ggplot2::coord_cartesian(clip = "off") +
  scale_fill_manual(values = c("#80B1D3", "#BC80BD", "#D9D9D9"), 
                    breaks = c('present', 'absent', 'gene absent'), 
                    name = 'Conserved PP motif', na.value = 'white') +
  theme(text = element_text(family="Avenir", size=18, color = "black"))

# ggsave("figure_4.svg", plot=figure_4, width = 14, height = 10)
# I did two aesthetic edits afterwards for this figure in an image editor
# (enlarged taxonomic circles in legend and shortened span of highlighted clades)


# Make EFP tree =================================
efp_mapping <- read.delim("efp_mapping.txt")

efp_tree <- read.tree("planctomycetes/planctomycetes_efp.tre")
efp_tree <- root(efp_tree, 'Roseimicrobium.gellanilyticum.2770939505_2771563938')

plotted_efp <- subset(efp_mapping, species_code %in% plotted_map$species_code)
plotted_efp$color <- plotted_map$color[match(plotted_efp$species_code, plotted_map$species_code)]

plotted_efp <- plotted_efp[, c('gene', colnames(plotted_efp)[-grep('^gene$', colnames(plotted_efp))])] #needs to be first column to link to tree
plotted_tree_efp <- drop.tip(phy = efp_tree, tip = efp_tree$tip.label[!efp_tree$tip.label %in% plotted_efp$gene])

earp_type_HGT <- getMRCA(plotted_tree_efp, subset(plotted_efp, modification == 'EarP')$gene)
epmab_type_HGT <- getMRCA(plotted_tree_efp, subset(plotted_efp, modification == 'EpmAB')$gene)
yeip_type_HGT <- getMRCA(plotted_tree_efp, subset(plotted_efp, modification == 'Unmodified YeiP')$gene)

efp_plot <- ggtree(plotted_tree_efp, open.angle=12, color = "#666666") %<+% plotted_efp +
  ggtitle('EFP tree') +
  geom_hilight(node = earp_type_HGT, extend = 1.5, alpha = 0.2) +
  geom_hilight(node = epmab_type_HGT, extend = 1.5, alpha = 0.2) +
  geom_hilight(node = yeip_type_HGT, extend = 1.5, alpha = 0.2) +
  geom_tiplab(aes(label = species), geom = 'text', size = 4, family = 'Avenir', align = T) +
  geom_tippoint(aes(fill = color), pch = 21, size = 5) +
  scale_fill_manual(values = rev(colorz), name = 'Family (GTDB-Tk)',
                    breaks = c("Isosphaeraceae", "Pirellulaceae", "Planctomycetaceae", "Thermoguttaceae",
                               "Unclassified", "Brocadiae", "Phycisphaerae", "Verrucomicrobiae")) +
  xlim_tree(3) +
  geom_treescale() +
  theme(legend.justification = c(0.025, 0.975), legend.position = c(0.025, 0.975), legend.title=element_blank()) +
  theme(text = element_text(family="Avenir", size=18, color = "black")) 


# Add EFP gene synteny and make figure S3 =================================
synteny <- read.csv('planctomycetes/Planctomycetota_gathered_refseq.csv')
kegg_data <- read.csv('planctomycetes/Planctomycetota_gathered_kos.csv')
synteny$ko_id <- kegg_data$ko_id[match(synteny$gene_oid, kegg_data$gene_oid)]

plotted_synteny <- setNames(data.frame(matrix(ncol = ncol(synteny) + 1, nrow = 0)), c(colnames(synteny), 'efp_geneoid'))
#extracts genes on either side of EFP -- can edit it to look for other genes ;)
for (efp in subset(plotted_efp, phylum == 'Planctomycetota')$gene_oid){
  print(efp)
  range_get = 4
  working <- subset(synteny, gene_oid < (efp + range_get) & gene_oid > (efp - range_get))
  working$efp_geneoid <- efp
  plotted_synteny <- rbind(plotted_synteny, working)
}

#clean-up some stuff
plotted_synteny$efp_id <- efp_mapping$gene[match(plotted_synteny$efp_geneoid, efp_mapping$gene_oid)]
plotted_synteny <- subset(plotted_synteny, efp_id %in% plotted_tree_efp$tip.label)
plotted_synteny$family <- clusters$family[match(plotted_synteny$gene_oid, clusters$gene)]
plotted_synteny$gene_name <- ifelse(is.na(plotted_synteny$ko_id) | plotted_synteny$ko_id == 'K06955', 
                                    plotted_synteny$family, 
                                    kegg_key$short_description[match(plotted_synteny$ko, kegg_key$ko_number)])
plotted_synteny$gene_name <- ifelse(plotted_synteny$type == 'tRNA', 'tRNA', plotted_synteny$gene_name)
plotted_synteny$gene_name[plotted_synteny$gene_name == 'FAM001846'] <- 'epmB'

#annotate earP specifically, unannotated by kegg
plotted_synteny$gene_name[plotted_synteny$gene_name %in% earP$gene] <- 'earP'

#for whatever reason these families have heterogeneous kegg annotations
multi_annot <- subset(as.data.frame(table(plotted_synteny$family, plotted_synteny$gene_name)), Freq != 0)
multi_annot <- names(table(multi_annot$Var1))[table(multi_annot$Var1) > 1]
plotted_synteny$gene_name <- ifelse(plotted_synteny$family %in% multi_annot, plotted_synteny$family, plotted_synteny$gene_name)

#composite plot, so need to reorder with efp plot order, this part gets that order 
plotted_synteny$efp_order <- make.unique(plotted_synteny$efp_id)
efp_order <- make.unique(rep(get_taxa_name(efp_plot), 1, each = (range_get * 2) - 1))
plotted_synteny$efp_order <- factor(plotted_synteny$efp_order, levels = efp_order)
plotted_synteny <- plotted_synteny[order(plotted_synteny$efp_order),]
#janky but it works to order the species as they are in the tree -- just gotta pop it off at the end
new_ids <- paste(rep(letters[1:26], each = 182), rep(letters[1:26], each = (range_get * 2) - 1, 2),
                               plotted_synteny$efp_id, sep = '_')
new_ids <- new_ids[1:nrow(plotted_synteny)]
plotted_synteny$efp_id <- new_ids
plotted_synteny <- plotted_synteny[, c('efp_id', colnames(plotted_synteny)[-grep('efp_id', colnames(plotted_synteny))])]
#should be first column to match

#to get both strands orientated in the same direction
backwards = subset(plotted_synteny, gene_name == 'efp' & strand == '+')$efp_id
plotted_synteny$start <- ifelse(! plotted_synteny$efp_id %in% backwards, -plotted_synteny$start, plotted_synteny$start)
plotted_synteny$end <- ifelse(! plotted_synteny$efp_id %in% backwards, -plotted_synteny$end, plotted_synteny$end)
plotted_synteny$orientation <- ifelse(plotted_synteny$strand == '+', 1, 0)

#these genes are huge and mess up the spacing, can leave in just looks nicer without them 
plotted_synteny <- subset(plotted_synteny, ! gene_oid %in% c(2691571115, 2638254724, 2730647146, 2730647145, 2887510974, 2758449130))

#color palette
countz <- table(plotted_synteny$gene_name)
plotted_synteny$color <- ifelse(plotted_synteny$gene_name %in% names(tail(sort(countz), n = 12)) | 
                                  plotted_synteny$gene_name %in% c('earP', 'epmB', 'epmA'), plotted_synteny$gene_name,
                                ifelse(plotted_synteny$gene_name %in% c('rfbB', 'rfbD', 'rfbC'), 'rfb', 'Other'))
gene_palette <- c('#ffb3ba', "#BC80BD", '#91DDF2', '#91DDF2', 
                  '#ffdfba', "white", '#FFFFB3', "#CCEBC5",
                  '#BDB2FF', '#DFFDFF', '#FEC868', '#FD8a8a', "#D9D9D9", '#A0C3D2')

#gggenes specific step -- have to add other variables called w/ ggplot cause it's just putting invisible new genes
#so efp is aligned in the center of each synteny
dummies <- make_alignment_dummies(plotted_synteny,
                                  aes(xmin = start, xmax = end, y = efp_id, id = gene_name), on = "efp")
dummies$orientation <- 1
dummies$color <- 'Other'

correct_order <- ggplot(plotted_synteny, aes(xmin = start, xmax = end, y = efp_id, 
                                             fill = color, forward = orientation, label = gene_name)) +
  geom_gene_arrow(arrow_body_height = grid::unit(3, "mm"), arrowhead_height = grid::unit(5, "mm"),
                  arrowhead_width = grid::unit(5, "mm")) +
  geom_gene_label(min.size = 7, grow = FALSE, reflow = FALSE, height = grid::unit(7, "mm"),
                  padding.y = grid::unit(0, "mm")) +
  geom_blank(data = dummies) +
  ylab('') +
  facet_wrap(~ efp_id, scales = "free", ncol = 1) +
  scale_fill_manual(values = gene_palette, name = 'Gene name') +
  theme(panel.border = element_rect(fill = NA, colour = NA)) +
  theme(panel.spacing = unit(-5, "mm")) +
  theme(legend.position = 'NONE') +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x =element_blank(), 
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank(), strip.background = element_blank(),
        strip.text.x = element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank())

figure_s3 <- cowplot::plot_grid(efp_plot, correct_order, rel_widths = c(1, 1), rel_heights = c(1, 0.9), rows = 1)
#ggsave("figure_s3.svg", plot=figure_2, width = 14, height = 14)
# in a figure editor I: manually adjusted clade highlights to look nicer & labelled HGT EFP types


# Species tree vs. phylogenetic tree and making figure s4 =================================
efp_distances <- subset(melt(cophenetic(plotted_tree_efp)), Var1 != Var2)
efp_distances$species_code1 <- sapply(strsplit(as.character(efp_distances$Var1), '_'), '[', 1)
efp_distances$species_code2 <- sapply(strsplit(as.character(efp_distances$Var2), '_'), '[', 1)
efp_distances <- subset(efp_distances, species_code1 != species_code2 & species_code1 %in% plotted_efp$species_code)
efp_distances$pair <- paste(efp_distances$species_code1, efp_distances$species_code2)
efp_distances$family_1 <- plotted_map$class_gtdbtk[match(efp_distances$species_code1, plotted_map$species_code)]
efp_distances$family_1 <- factor(efp_distances$family_1,
                                 levels = c('Planctomycetes', names(table(subset(plotted_map, class_gtdbtk != 'Planctomycetes')$class_gtdbtk))))
efp_distances <- efp_distances[order(efp_distances$family_1), ]
phylogeny <- subset(melt(cophenetic(plotted_tree)), Var1 != Var2)
efp_distances$phylogenetic_distance <- phylogeny$value[match(efp_distances$pair, paste(phylogeny$Var1, phylogeny$Var2))]

# pruned_phy <- remove_pairs(efp_distances, 'Var1', 'Var2')
pruned_phy <- efp_distances
colnames(pruned_phy) <- c('efp_1', 'efp_2', 'efp_distance', 'species_code1', 'species_code2', 'pair_code', 'family_1', "phylogenetic_distance")
pruned_phy$family_2 <- plotted_map$class_gtdbtk[match(pruned_phy$species_code2, plotted_map$species_code)]
pruned_phy$fam_color <- paste(pruned_phy$family_1, pruned_phy$family_2, sep = ' x ')
pruned_phy <- subset(pruned_phy, grepl('Verruco', fam_color) == FALSE)

#group EFP
pruned_phy$efp_fam1 <- plotted_efp$modification[match(pruned_phy$efp_1, plotted_efp$gene)]
pruned_phy$efp_fam2 <- plotted_efp$modification[match(pruned_phy$efp_2, plotted_efp$gene)]
pruned_phy$efp_residue1 <- plotted_efp$conserved_residue[match(pruned_phy$efp_1, plotted_efp$gene)]
pruned_phy$efp_residue2 <- plotted_efp$conserved_residue[match(pruned_phy$efp_2, plotted_efp$gene)]

#combine same ids just different order
pruned_phy$comparision <- ifelse(pruned_phy$efp_fam1 == pruned_phy$efp_fam2, 'Self', 
                                 paste(pruned_phy$efp_fam1, pruned_phy$efp_fam2, sep = ' x '))
pruned_phy$comparision <- stri_replace_all_regex(pruned_phy$comparision,
                       pattern=c('Unmodified YeiP x EarP', 'Unmodified YeiP x EpmAB', 'Unknown x EarP', 'EpmAB x EarP', 'Unmodified YeiP x Unknown', 'Unknown x EpmAB'),
                       replacement=c('EarP x Unmodified YeiP', 'EpmAB x Unmodified YeiP', 'EarP x Unknown', 'EarP x EpmAB', 'Unknown x Unmodified YeiP', 'EpmAB x Unknown'),
                       vectorize=FALSE)

# get maximum values to plot to keep the plot area a perfect square, so 1:1 line splits on the diagonal 
max = ifelse(max(pruned_phy$efp_distance) > max(pruned_phy$phylogenetic_distance), 
             max(pruned_phy$efp_distance), 
             max(pruned_phy$phylogenetic_distance))

pruned_phy$efp_fam1 <- factor(pruned_phy$efp_fam1, levels = c('EarP', 'EpmAB', 'Unknown', 'Unmodified YeiP'))

figure_s4 <- ggplot(pruned_phy, aes(y = phylogenetic_distance, x = efp_distance)) +
  facet_wrap(~paste(efp_fam1, 'EFP type', sep = ' ')) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = 'purple', linewidth= 1.01) +
  geom_point(aes(fill = comparision), pch = 21, color = '#666666', size = 5, alpha = 1) +
  xlim(0, max + 0.05*max) + #keep it a square
  ylim(0, max + 0.05*max) + #keep it a square
  xlab("EFP distance") +
  ylab("Phylogenetic distance") +
  scale_fill_brewer(palette = 'Spectral', name = 'EFP type comparision') +
  theme_linedraw() +
  theme(strip.background= element_rect(colour = "black", fill = "grey")) +
  theme(strip.text = element_text(family="Avenir", size=20, color = "black")) +
  theme(axis.text = element_text(family="Avenir", size=20, color = "black")) + 
  theme(axis.title = element_text(family="Avenir", size=20, color = "black")) +
  theme(text = element_text(family="Avenir", size=20, color = "black")) +
  theme(legend.text = element_text(family="Avenir", size=16, color = "black")) 

#ggsave("figure_s4.svg", plot=figure_s4, width = 11, height = 8)


