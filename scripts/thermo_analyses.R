## Thermotogota specific analyses for pattern project paper
##
## Purpose: Do analyses for the thermotogota genomes and make figures 2, 3, S1, S5
##
## Author: Dr. Tess Brewer ;)

# Load functions and color palettes =================================
source('scripts/pattern_project_functions.R') #make sure to change working directory in this script first

#palette for the figures
colorz <- c("#2373A2", "#01ff70",
            "#fcff4e", "#FDA257", "#FF0066")

# File set-up =================================

#subset for thermotogota
root = "Acetomicrobium.thermoterrenum.2599185146"
plotted_tree <- read.tree('thermotogota/thermotoga_species.tre')
plotted_tree <- root(plotted_tree, root)
plotted_map <- subset(mapping, species_code %in% plotted_tree$tip.label & species_code != root)
plotted_map <- plotted_map[, c('species_code', colnames(plotted_map)[-grep('species_code', colnames(plotted_map))])]
#those ids have to be first column to link with the tree (also tree tip labels)

plotted_map <- plotted_map[order(row.names(plotted_map)),]
plotted_tree <- drop.tip(plotted_tree, plotted_tree$tip.label[ ! plotted_tree$tip.label %in% plotted_map$species_code])
# plotted_tree$tip.label <- gsub('\\_', '\\.', plotted_tree$tip.label)
plotted_map$hgt_efp <- ifelse(plotted_map$family_gtdbtk == 'Petrotogaceae' & plotted_map$genus_gtdbtk != 'Marinitoga', 'yes', 'no')
plotted_map$hgt_group <- ifelse(plotted_map$species_code %in% subset(plotted_map, efp_modification == 'EarP')$species_code, 'HGT Event 2',
                                ifelse(plotted_map$species_code %in% subset(plotted_map, family_gtdbtk == 'Petrotogaceae' & 
                                                genus_gtdbtk != 'Marinitoga')$species_code, 'HGT Event 1', ''))
#nodes to place the green plus and red x on
HGT_transfer2 <- getMRCA(plotted_tree, subset(plotted_map, efp_modification == 'EarP')$species_code)
HGT_transfer1 <- getMRCA(plotted_tree, subset(plotted_map, family_gtdbtk == 'Petrotogaceae' & genus_gtdbtk != 'Marinitoga' & efp_modification != 'EarP')$species_code)
ancestral_loss <- getMRCA(plotted_tree, subset(plotted_map, hgt_efp == 'yes')$species_code)

#base tree
full_tree <- ggtree(plotted_tree, open.angle=12, color = "#666666") %<+% plotted_map +
  ggtitle('Species tree') +
  geom_tiplab(aes(label = species), size = 6, align=F, family="Avenir", offset = 0.025) +
  geom_tippoint(aes(fill = family_gtdbtk, size = genome_size/1000000), pch = 21, color = "#666666") +
  geom_nodepoint(aes(subset=(node == ancestral_loss)), shape=4, size=5, color='#D74840', stroke = 2) +
  geom_nodepoint(aes(subset=(node %in% c(HGT_transfer1, HGT_transfer2))), shape=3, size=5, color='#008B52', stroke = 2) +
  geom_treescale() +
  scale_fill_manual(values = colorz, name = 'Family (GTDB-Tk)') +
  scale_size_continuous(name = 'Genome size (MB)', range = c(3, 6)) +
  theme(text = element_text(family="Avenir", size=20, color = "black")) +
  theme(legend.text = element_text(family="Avenir", size=17, color = "black")) +
  guides(fill = guide_legend(override.aes = list(size=5)))

# Making figure s5 =================================

#add label of EFP conserved residue region and make conserved residue red
plotted_map$label_red <- paste(quotes(substr(plotted_map$actino_loop, 0, 3)), phantom('A'), quotes(substr(plotted_map$actino_loop, 5, 7)), sep = " * ")
plotted_map$label_black <- paste(phantom(substr(plotted_map$actino_loop, 0, 3)), quotes(substr(plotted_map$actino_loop, 4, 4)), phantom(substr(plotted_map$actino_loop, 5, 7)), sep = " * ")

#need to reattach edited plotted_map
full_tree <- ggtree(plotted_tree, open.angle=12, color = "#666666") %<+% plotted_map +
  ggtitle('Species tree') +
  geom_hilight(node = ancestral_loss, extend = 1.5, alpha = 0.2) +
  geom_tiplab(aes(label = species), size = 6, align=F, family="Avenir", offset = 0.025) +
  geom_tippoint(aes(fill = family_gtdbtk, size = genome_size/1000000), pch = 21, color = "#666666") +
  geom_nodepoint(aes(subset=(node == ancestral_loss)), shape=4, size=5, color='#D74840', stroke = 2) +
  geom_nodepoint(aes(subset=(node %in% c(HGT_transfer1, HGT_transfer2))), shape=3, size=5, color='#008B52', stroke = 2) +
  geom_treescale() +
  scale_fill_manual(values = colorz, name = 'Family (GTDB-Tk)') +
  scale_size_continuous(name = 'Genome size (MB)', range = c(3, 6)) +
  theme(text = element_text(family="Avenir", size=20, color = "black")) +
  theme(legend.text = element_text(family="Avenir", size=17, color = "black")) +
  guides(fill = guide_legend(override.aes = list(size=5)))


heatmap_layer <- as.data.frame(plotted_map[, c("efp_modification")])
row.names(heatmap_layer) <- plotted_map$species_code
heat_plot <- full_tree + new_scale_fill()

figure_s5 <- gheatmap(heat_plot, heatmap_layer, width=0.08, offset = 0.3, hjust = 0, color=NA,
         colnames_position = 'top', colnames_offset_x = 0, colnames = FALSE, font.size = 4, 
         colnames_offset_y = -0.1, colnames_angle = 45, family = 'Avenir') +
  xlim_tree(1.1) +
  geom_text_repel(aes(x = 1.08, label = label_black), color = "red", parse = TRUE, seed = 1, hjust = 'left', family = "mono", size = 5) +
  geom_text_repel(aes(x = 1.08, label = label_red), color = "black", parse = TRUE, seed = 1, hjust = 'left', family = "mono", size = 5) +
  scale_fill_manual(values = c('#BC80BD', '#ededed'), name = 'EFP modification') +
  theme(text = element_text(family="Avenir", size=18, color = "black"))

# ggsave("figure_s5.svg", plot=figure_s5, width = 14, height = 8)


# Load annotations =================================

#kegg annotations of genes in thermotogota genomes
ko_annotations <- load_data("thermotogota/ko_annotations/", "*.ko.tab.txt")
ko_annotations <- ko_annotations[order(ko_annotations$bit_score, decreasing = TRUE), ]
ko_annotations <- ko_annotations[!duplicated(ko_annotations$gene_oid), ]
ko_annotations$ko_id <- gsub('KO:K', 'K', ko_annotations$ko_id)

#silix clusters
clusters <- read.delim("thermotogota/thermotogota_silix.fnodes", header = FALSE, col.names = c('family', 'gene'))
clusters$genome <- sapply(strsplit(as.character(clusters$gene), '_'), '[', 2)
clusters$gene <- sapply(strsplit(as.character(clusters$gene), '_'), '[', 1)
clusters$ko <- ko_annotations$ko_id[match(paste(clusters$genome, clusters$gene), paste(ko_annotations$family, ko_annotations$gene_oid))]
clusters$ko[is.na(clusters$ko)] <- 'unassigned'
annot <- as.data.frame(table(clusters$family, clusters$ko))
annot <- annot[order(annot$Freq, decreasing = TRUE), ]
annot <- annot[!duplicated(annot$Var1), ]
annot$description <- ifelse(annot$Var2 == 'unassigned', as.character(annot$Var1), 
                            kegg_key$short_description[match(annot$Var2, kegg_key$ko_number)])


# Running tests =================================
run_data <- list(data = 'thermotogota/pp_search/',
                 file_pattern = '_conserved_stall.csv',
                 sequence_pattern = '..PP.',
                 group_name = 'HGT',
                 group = c(list(subset(plotted_map, hgt_efp == 'yes')$taxon_oid)))

#just run the function in a loop
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

#to look directly at the table for particular family of proteins, for example ycaJ
stall_data <- load_data(run_data$data[[1]], run_data$file_pattern[[1]])
kegg_key[grep('ycaJ', kegg_key$short_description),] # get kegg ko number for protein of interest
annot[annot$Var2 == 'K07478',] # get silix family for that ko, some kos have multiple families
look <- subset(stall_data, family == 'FAM000117') # #extract from stall data
table(look$conserved_stall)

# Take a look at hits / make figure 3 =================================

#get pp motif presence / absence / gene absence info
nonidentical_stalls <- produce_nonidential_stalls(run_data$data[1], run_data$file_pattern[1], run_data$sequence_pattern[1])
fam_to_plot <- nonidentical_stalls[, colnames(nonidentical_stalls) %in% results_final$family]
# fam_to_plot <- nonidentical_stalls[, colnames(nonidentical_stalls) %in% c('FAM000117', 'FAM001406')] #to look at particular families of interest
fam_to_plot$hgt <- plotted_map$hgt_efp[match(row.names(fam_to_plot), plotted_map$taxon_oid)]

#automate ordering so families with similar patterns are grouped together and everything looks nice
orderingz <- subset(melt(aggregate(. ~ hgt, fam_to_plot, mean)), hgt == 'yes')
orderingz <- droplevels(orderingz$variable[order(orderingz$value, decreasing = TRUE)])
fam_to_plot$hgt <- NULL #need to get rid of it for future dataframe

#convert pp motif presence to categorical and label correctly
plotting_data <- as.data.frame(sapply(fam_to_plot, function(x){str_replace_all(x, c('0' = 'gene absent', '1' = 'absent', '2' = 'present'))}))
row.names(plotting_data) <- row.names(fam_to_plot)
row.names(plotting_data) <- plotted_map$species_code[match(row.names(plotting_data), plotted_map$taxon_oid)]
plotting_data <- plotting_data[, orderingz]
colnames(plotting_data) <- paste(annot$description[match(colnames(plotting_data), annot$Var1)])
colnames(plotting_data) <- gsub('FAM000', 'FAM', colnames(plotting_data)) #all those zeroes are a bummer

new_plot <- full_tree + new_scale_fill()
figure_3 <- gheatmap(new_plot, plotting_data, width=0.7, offset = 0.5, color = "white", hjust = 0,
                  colnames_position = 'top', colnames_offset_x = 0, font.size = 5,
                  colnames_offset_y = -0.1, colnames_angle = 45, family = 'Avenir') +
  xlim_tree(1.5) +
  ggplot2::ylim(0, 25) +
  ggplot2::coord_cartesian(clip = "off") +
  scale_fill_manual(values = c("#80B1D3", "#BC80BD", "#D9D9D9"), 
                    breaks = c('present', 'absent', 'gene absent'), 
                    name = 'Conserved PP motif', na.value = 'white') +
  theme(text = element_text(family="Avenir", size=18, color = "black"))

#ggsave("figure_3.svg", plot=figure_3, width = 14, height = 8)

# I did two aesthetic edits afterwards for this figure in an image editor
# (enlarged taxonomic circles in legend and shortened span of highlighted clades)

# Make EFP tree =================================
efp_mapping <- read.delim("efp_mapping.txt")
efp_tree <- read.tree("thermotogota/thermotogota_efp.tre")
efp_tree <- drop.tip(efp_tree, efp_tree$tip.label[!efp_tree$tip.label %in% efp_mapping$gene])

plotted_efp <- subset(efp_mapping, species_code %in% plotted_map$species_code)
plotted_efp <- plotted_efp[, c('gene', colnames(plotted_efp)[-grep('^gene$', colnames(plotted_efp))])] #needs to be first column to link to tree
plotted_tree_efp <- drop.tip(phy = efp_tree, tip = efp_tree$tip.label[!efp_tree$tip.label %in% plotted_efp$gene])

unknown_type_HGT <- getMRCA(plotted_tree_efp, subset(plotted_efp, grepl('Petrotoga', species) | grepl('Defluviitoga', species))$gene)

efp_plot <- ggtree(plotted_tree_efp, open.angle=12, color = "#666666", size = 0.75) %<+% plotted_efp +
  ggtitle('EFP tree') +
  geom_treescale(x = 0.9, y = 8) +
  geom_hilight(node = c(2), extend = 4, alpha = 0.2) +
  geom_hilight(node = c(1), extend = 4, alpha = 0.2) +
  geom_hilight(node = unknown_type_HGT, extend = 1.5, alpha = 0.2) +
  geom_tiplab(aes(label = species), geom = 'text', size = 6, hjust = -0.075, family = 'Avenir', align = T) +
  geom_tippoint(aes(fill = family), pch = 21, size = 5) +
  scale_fill_manual(values = colorz[1:4], name = 'Family (GTDB-Tk)') +
  theme(legend.justification = c(0.025, 0.925), legend.position = c(0.025, 0.925)) +
  theme(text = element_text(family="Avenir", size=20, color = "black")) +
  theme(plot.title = element_text(margin=margin(0,0,-20,0))) +
  theme(legend.text = element_text(family="Avenir", size=17, color = "black")) +
  xlim_tree(2.5) + 
  ggplot2::coord_cartesian(clip = "off")
  
  
# Add EFP gene synteny and make figure 2 =================================
synteny <- read.csv('thermotogota/gathered_refseq.csv')
kegg_data <- read.csv('thermotogota/gathered_kos.csv')
synteny$ko_id <- kegg_data$ko_id[match(synteny$gene_oid, kegg_data$gene_oid)]

plotted_synteny <- setNames(data.frame(matrix(ncol = ncol(synteny), nrow = 0)), colnames(synteny))
#extracts genes on either side of EFP
for (efp in efp_mapping$gene_oid){
  range_get = 4
  working <- subset(synteny, gene_oid < (efp + range_get) & gene_oid > (efp - range_get))
  plotted_synteny <- rbind(plotted_synteny, working)
}

plotted_synteny$full_name <- kegg_key$description[match(plotted_synteny$ko, kegg_key$ko_number)]
plotted_synteny$gene_name <- kegg_key$short_description[match(plotted_synteny$ko, kegg_key$ko_number)]
#annotate earP specifically, unannotated by kegg
plotted_synteny$gene_name[plotted_synteny$gene_oid %in% earP$gene] <- 'earP'

#annotate some poorly named things
plotted_synteny$gene_name[grepl('YloU', plotted_synteny$attributes)] <- 'yloU'
plotted_synteny$gene_name[is.na(plotted_synteny$gene_name)] <- 'Other'
plotted_synteny$gene_name[grepl('K07133', plotted_synteny$gene_name)] <- 'Other'
plotted_synteny$gene_name[grepl('E2.2.1.1', plotted_synteny$gene_name)] <- 'tktA'
plotted_synteny$gene_name[grepl('NTH', plotted_synteny$gene_name)] <- tolower(plotted_synteny$gene_name[grepl('NTH', plotted_synteny$gene_name)])


plotted_synteny$efp_id <- efp_mapping$gene[match(plotted_synteny$genome, efp_mapping$taxon_oid)]
plotted_synteny$efp_order <- make.unique(plotted_synteny$efp_id)

#composite plot, so need to reorder with efp plot order, this part gets that order 
efp_order <- make.unique(rep(get_taxa_name(efp_plot), 1, each = (range_get * 2) - 1))
plotted_synteny$efp_order <- factor(plotted_synteny$efp_order, levels = efp_order)
plotted_synteny <- plotted_synteny[order(plotted_synteny$efp_order),]
#janky but it works to order the species as they are in the tree -- just gotta pop it off at the end
plotted_synteny$efp_id <- paste(rep(letters[1:length(table(plotted_synteny$genome))], each = (range_get * 2) - 1),
                                plotted_synteny$efp_id, sep = '_')

plotted_synteny <- plotted_synteny[, c(ncol(plotted_synteny), 1:(ncol(plotted_synteny)-1))]
plotted_synteny$species <- plotted_map$species[match(plotted_synteny$genome, plotted_map$taxon_oid)]
countz <- table(plotted_synteny$gene_name)
plotted_synteny$color <- ifelse(plotted_synteny$gene_name %in% names(tail(sort(countz), n = 6)) | 
                                  plotted_synteny$gene_name %in% c('deoC', 'cheY', 'earP', 'yloU'), plotted_synteny$gene_name,
                                ifelse(plotted_synteny$gene_name %in% c('rfbB', 'rfbD', 'rfbC'), 'rfb', 'Other'))
plotted_synteny$species_code <- plotted_map$species_code[match(plotted_synteny$genome, plotted_map$taxon_oid)]
plotted_synteny <- plotted_synteny[, c('efp_id', colnames(plotted_synteny)[-grep('efp_id', colnames(plotted_synteny))])]
names(table(plotted_synteny$color))
gene_palette <- c("#FCCDE5", '#ffdfba', '#ffb3ba',  "#BC80BD", "white", '#91DDF2', "#D9D9D9", '#FFFFB3', "#BEBADA", "#CCEBC5")

#to get both strands orientated in the same direction
backwards = subset(plotted_synteny, gene_name == 'efp' & strand == '+')$genome
plotted_synteny$start <- ifelse(! plotted_synteny$genome %in% backwards, -plotted_synteny$start, plotted_synteny$start)
plotted_synteny$end <- ifelse(! plotted_synteny$genome %in% backwards, -plotted_synteny$end, plotted_synteny$end)

#last clean up step
plotted_synteny <- subset(plotted_synteny, !gene_oid %in% c(2651318000, 2651318001)) #janky numbering
plotted_synteny$orientation <- ifelse(plotted_synteny$strand == '+', 1, 0)

#gggenes specific step -- have to add other variables called w/ ggplot cause it's just putting invisible new genes to flush out efp
dummies <- make_alignment_dummies(plotted_synteny,
                                  aes(xmin = start, xmax = end, y = efp_id, id = gene_name), on = "efp")
dummies$orientation <- 1
dummies$color <- 'Other'

correct_order <- ggplot(plotted_synteny, aes(xmin = start, xmax = end, y = efp_id, 
                            fill = color, forward = orientation, label = gene_name)) +
  geom_gene_arrow(arrow_body_height = grid::unit(4, "mm"), arrowhead_height = grid::unit(6, "mm"),
                  arrowhead_width = grid::unit(6, "mm")) +
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
        strip.text.x = element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())

figure_2 <- cowplot::plot_grid(efp_plot, NULL, correct_order, rel_widths = c(1.5, 0.10, 1), rows = 1)
# ggsave("figure_2.svg", plot=figure_2, width = 20, height = 8)
# in a figure editor I: manually made scale break in one tree branch, extended highlight of clades & labelled HGT EFP types


# Species tree vs. phylogenetic tree and making figure s1 =================================

phylogeny <- subset(melt(cophenetic(plotted_tree)), Var1 != Var2)
# melt species tree and remove self distances (species 1 to species 1 is distance of 0)
phylogeny$family_1 <- plotted_map$family_gtdbtk[match(phylogeny$Var1, plotted_map$species_code)]
phylogeny$family_1 <- factor(phylogeny$family_1, levels = c('Petrotogaceae', 'Fervidobacteriaceae', 'Kosmotogaceae', 'Thermotogaceae'))
phylogeny <- phylogeny[order(phylogeny$family_1), ]
pruned_phy <- remove_pairs(phylogeny, 'Var1', 'Var2')
# remove doubles (ex: cophenetic gives both species 1 vs species 2 and species 2 vs species 1)
efp_distances <- subset(melt(cophenetic(plotted_tree_efp)), Var1 != Var2)
# melt efp tree and remove self distances (species 1 to species 1 is distance of 0)
efp_distances$pair <- paste(sapply(strsplit(as.character(efp_distances$Var1), '_'), '[', 1), 
                            sapply(strsplit(as.character(efp_distances$Var2), '_'), '[', 1),
                            sep = '_')
colnames(pruned_phy) <- c('taxa_1', 'taxa_2', 'phylogenetic_distance', 'family_1', 'pair_code')
pruned_phy$efp_distance <- efp_distances$value[match(pruned_phy$pair_code, efp_distances$pair)]
pruned_phy$family_2 <- plotted_map$family_gtdbtk[match(pruned_phy$taxa_2, plotted_map$species_code)]
pruned_phy$fam_color <- paste(pruned_phy$family_1, pruned_phy$family_2, sep = ' x ')

pruned_phy <- subset(pruned_phy, taxa_1 != strsplit(root, '_')[[1]][1]
                     & taxa_2 != strsplit(root, '_')[[1]][1])
#remove root from this figure

max = ifelse(max(pruned_phy$efp_distance) > max(pruned_phy$phylogenetic_distance), 
             max(pruned_phy$efp_distance), 
             max(pruned_phy$phylogenetic_distance))
# get maximum values to plot to keep the plot area a perfect square, so 1:1 line splits on the diagonal 

pruned_phy$highlight_species <- plotted_map$hgt_group[match(pruned_phy$taxa_1, plotted_map$species_code)]
pruned_phy$highlight_species_secondary <- plotted_map$hgt_group[match(pruned_phy$taxa_2, plotted_map$species_code)]
#highlight HGT stuff but only the comparision of HGT groups versus non HGT groups

micro.hulls <- ddply(subset(pruned_phy, highlight_species != '' &
                                        highlight_species != highlight_species_secondary), 
                                        "highlight_species", find_hull)
#make the polygon outlines in the figure

figure_s2 <- ggplot(pruned_phy, aes(y = phylogenetic_distance, x = efp_distance)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = 'purple', linewidth = 1.01) +
  geom_polygon(data = micro.hulls, aes(color = highlight_species), alpha = 0) +
  geom_point(aes(fill = fam_color), pch = 21, color = '#666666', size = 3, alpha = 1) +
  geom_text(data = subset(pruned_phy, ! duplicated(highlight_species)), aes(label = highlight_species), 
            check_overlap = FALSE, vjust = -6, family="Avenir", size=6, color = "black") +
  xlim(0, max + 0.05*max) +
  ylim(0, max + 0.05*max) +
  #keep it a square
  xlab("EFP distance") +
  ylab("Phylogenetic distance") +
  scale_color_manual(values = c('#727272', '#727272'), guide = 'none') +
  scale_fill_brewer(palette = 'Spectral', name = '') +
  theme_linedraw() +
  theme(legend.justification = c(0.025, 0.985), legend.position = c(0.025, 0.985), legend.title=element_blank()) +
  theme(axis.text = element_text(family="Avenir", size=20, color = "black")) + 
  theme(axis.title = element_text(family="Avenir", size=20, color = "black")) +
  theme(text = element_text(family="Avenir", size=20, color = "black")) +
  theme(legend.text = element_text(family="Avenir", size=15, color = "black")) 

# ggsave("figure_s2.svg", plot=figure_s2, width = 8.5, height = 8.5)
