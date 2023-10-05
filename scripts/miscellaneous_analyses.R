## Miscellaneous analyses for pattern project paper
##
## Purpose: Make figures 1, S1, 5, and Table S1
##
## Author: Dr. Tess Brewer ;)

# Load functions  =================================
source('scripts/pattern_project_functions.R') #make sure to change working directory in this script first

# Make figure 1 =================================

mapping <- mapping[, c('taxon_oid', colnames(mapping)[-grep('taxon_oid', colnames(mapping))])]
#phylogenetic tree for full 3000 dataset has taxon_oid as tip labels, has to be first column to link tree to mapping file
mapping$modification_tree <- stri_replace_all_regex(mapping$efp_modification,
                                                        pattern=names(table(mapping$efp_modification)),
                                                        replacement=c('EarP', 'EarP', 'EpmAB', 'EpmAB + EarP', 'EpmAB', 'EpmAB',
                                                                      'EpmAB', 'EpmABC', 'EpmABC + EarP', 'EpmABC + EarP', 'EpmABC', 'EpmABC',
                                                                      'Unknown', 'Unknown', 'Unknown', 'Unmodified-proline-loop', 'YmfI'),
                                                        vectorize=FALSE)
mapping$modification_tree[grepl('Unmodified YeiP', mapping$modification_tree)] <- 'EpmABC + EarP'

phylogeny_tree <- read.tree("full_species.tre")
phylogeny_tree <- root(phylogeny_tree, "GCF_000009185.1_ASM918v1_genomic") # Haloquadratum walsbyi (archaea)
phylogeny_tree <- drop.tip(phy = phylogeny_tree, 
                           tip = phylogeny_tree$tip.label[!phylogeny_tree$tip.label %in% 
                                                         mapping$taxon_oid])

# unifying species in mapping file and tree
plotted_tree <- drop.tip(phy = phylogeny_tree, 
                         tip = phylogeny_tree$tip.label[!phylogeny_tree$tip.label %in% 
                                                          mapping$taxon_oid])
plotted_map <- subset(mapping, taxon_oid %in% plotted_tree$tip.label)

# picking which phyla should be indicated on legend
countz <- table(plotted_map$phylum_gtdbtk)
plotted_map$color <- ifelse(plotted_map$phylum_gtdbtk %in% names(tail(sort(countz), n = 8)) | 
                            plotted_map$phylum_gtdbtk %in% c('Chloroflexota', 'Planctomycetota', 'Thermotogota'), 
                            plotted_map$phylum_gtdbtk, 'Other')
palettez <- brewer.pal(n = (length(table(plotted_map$color)) - 1), name = 'Spectral')
palettez <- append(palettez, 'white', after = (grep('Other', names(table(plotted_map$color))) - 1)) #other should be white

# getting nodes for the many many MANY efp HGT events
a_event <- getMRCA(plotted_tree, tip = as.character(subset(plotted_map, phylum_gtdbtk == 'Spirochaetota' & grepl('EpmAB|EarP', efp_modification))$taxon_oid))
b_event <- getMRCA(plotted_tree, tip = as.character(subset(plotted_map, phylum_gtdbtk == 'Planctomycetota' & grepl('EpmAB', efp_modification))$taxon_oid))
b_eventb <- getMRCA(plotted_tree, tip = as.character(subset(plotted_map, phylum_gtdbtk == 'Verrucomicrobiota' & grepl('EpmAB', efp_modification))$taxon_oid))
c_event <- getMRCA(plotted_tree, tip = as.character(subset(plotted_map, phylum_gtdbtk == 'Fibrobacterota' & grepl('EpmAB', efp_modification))$taxon_oid))
d_event <- getMRCA(plotted_tree, tip = as.character(subset(plotted_map, phylum_gtdbtk == 'Fusobacteriota' & grepl('EarP', efp_modification))$taxon_oid))
e_event <- getMRCA(plotted_tree, tip = as.character(subset(plotted_map, phylum_gtdbtk == 'Thermotogota' & grepl('EarP', efp_modification))$taxon_oid))
f_event <- getMRCA(plotted_tree, tip = as.character(subset(plotted_map, phylum_gtdbtk == 'Chloroflexota' & grepl('EpmAB', efp_modification))$taxon_oid))
g_event <- getMRCA(plotted_tree, tip = as.character(subset(plotted_map, phylum_gtdbtk %in% c('Aquificota', 'Campylobacterota') & grepl('EpmAB', efp_modification))$taxon_oid))
h_event <- getMRCA(plotted_tree, tip = as.character(subset(plotted_map, phylum_gtdbtk %in% c('Desulfobacterota', 'Desulfuromonadota', 'Myxococcota')
                                                           & grepl('EpmAB', efp_modification))$taxon_oid))
i_event <- getMRCA(plotted_tree, tip = as.character(subset(plotted_map, phylum_gtdbtk == 'Alphaproteobacteria' & grepl('EpmAB', efp_modification))$taxon_oid))
j_event <- getMRCA(plotted_tree, tip = as.character(subset(plotted_map, order_gtdbtk == 'Pseudomonadales' & grepl('EarP', efp_modification))$taxon_oid))
j_eventb <- getMRCA(plotted_tree, tip = as.character(subset(plotted_map, family_gtdbtk %in% c('Alteromonadaceae', 'Shewanellaceae') & grepl('EarP', efp_modification))$taxon_oid))
j_eventc <- getMRCA(plotted_tree, tip = as.character(subset(plotted_map, family_gtdbtk == 'Succinivibrionaceae' & grepl('EarP', efp_modification))$taxon_oid))
j_eventd <- getMRCA(plotted_tree, tip = as.character(subset(plotted_map, family_gtdbtk == 'Aeromonadaceae' & grepl('EarP', efp_modification))$taxon_oid))

tree_plot <- ggtree(plotted_tree, open.angle=20, color = "#666666", layout = "fan") %<+% plotted_map +
  theme(text = element_text(family="Avenir", size=16, color = "black")) +
  guides(color = guide_legend(override.aes = list(size=5)))

# remove long branch of intracellular genomes for readability
long_branch <- getMRCA(plotted_tree, tip = as.character(subset(mapping, order_gtdbtk %in% c('Mycoplasmatales', 'Acholeplasmatales'))$taxon_oid))
tree_plot <- ggtree::collapse(tree_plot, long_branch)

# add clade highlights where HGT occurred (have to do after collapsing branch, otherwise throws off position)
tree_plot <- tree_plot + 
  geom_hilight(node = c(a_event, b_event, b_eventb, c_event, d_event, e_event, f_event, 
                                         g_event, h_event, i_event, 
                                         j_event, j_eventb, j_eventc, j_eventd), extendto = 2, alpha = 0.2)
# add tips last so they sit on top of clade highlight
tree_plot <- tree_plot + 
  geom_tippoint(aes(fill = color), pch = 21, size =3, color = "#666666") +
  scale_fill_manual(values = palettez, name = 'Phylum (GTDB-Tk)')
  
#add heatmap with simplified modifications (epmAB / epmABC / earP / YmfI / unmodified loop / unknown)
trial <- as.data.frame(plotted_map[, c("modification_tree")])
row.names(trial) <- plotted_map$taxon_oid

new_plot <- tree_plot + new_scale_fill()
figure_1 <- gheatmap(new_plot, trial, width=0.07, color=NA, colnames = FALSE) +
  scale_fill_manual(values = c('#BC80BD', '#8dde9f', '#3dc85d', '#8dcdde', '#2e88a0', '#ededed', '#FFFFB3', 'pink'), name = 'EFP modification') +
  theme(text = element_text(family="Avenir", size=18, color = "black"))

# ggsave("figure_1.svg", plot=figure_1, width = 14, height = 12 )
# labels added in inkscape

# Supplemental figures for clustering families =================================

#how many EFP are in each family
countzi <- as.data.frame(table(efp_mapping$fam))
countzi$Freq_final <- paste('n = ', countzi$Freq, sep = '')

efp_mapping$plotted_phylum <- ifelse(efp_mapping$phylum %in% names(tail(sort(countz), n = 5)) | 
                              efp_mapping$phylum %in% c('Chloroflexota', 'Planctomycetota', 'Thermotogota', 'Acidobacteriota', 'Verrucomicrobiota', 'Spirochaetota'), 
                            efp_mapping$phylum, 'Other')

#not the best way to do this but it's fine
plot_bar <- as.data.frame(table(efp_mapping$fam, efp_mapping$modification))
plot_bar$type <- c('Post-translational modification')
plot_bar1 <- as.data.frame(table(efp_mapping$fam, efp_mapping$conserved_residue))
plot_bar1$type <- c('Conserved amino acid')
plot_bar2 <- as.data.frame(table(efp_mapping$fam, efp_mapping$plotted_phylum))
plot_bar2$type <- c('Phylogenetic diversity')
plot_bar <- rbind(plot_bar, plot_bar1, plot_bar2)
rm(plot_bar1, plot_bar2)
plot_bar$numeric_count <- countzi$Freq[match(plot_bar$Var1, countzi$Var1)]
plot_bar$family_count <- countzi$Freq_final[match(plot_bar$Var1, countzi$Var1)]
plot_bar$y_label <- paste(plot_bar$Var1, plot_bar$family_count)
plot_bar$proportional <- plot_bar$Freq / plot_bar$numeric_count

#color palette and assigning grey to other
color_vals <- names(table(subset(plot_bar, type == 'Phylogenetic diversity')$Var2))
color_pal <- c(brewer.pal(11, 'Spectral'))
color_pal <- append(color_pal, 'grey', (grep('Other', color_vals) - 1))
plot_bar$Var1 <- gsub('fam', '', plot_bar$Var1)
plot_bar$Var1 <- gsub('^0', '', plot_bar$Var1)

#substitute three letter codes for one letter
plot_bar$Var2 <- stri_replace_all_regex(plot_bar$Var2, 
                                        pattern = c('^-$','^G$', '^K$', '^M$', '^N$', '^Q$', '^R$', '^S$'), 
                                        replacement = c('-','Gly', 'Lys', 'Met', 'Asn', 'Gln', 'Arg', 'Ser'), vectorize = FALSE)
#specific to donut chart
hsize <- 2
plot_bar <- plot_bar %>% 
  mutate(x = hsize)

amino_acid <- ggplot(subset(plot_bar, type == 'Conserved amino acid'), aes(x = hsize, y = proportional, fill = Var2)) +
  facet_wrap(~ y_label, ncol = 1) +
  geom_col(color = '#666666') +
  # ggtitle('Conserved amino acid', ) +
  geom_text(aes(label = Var1), y = 0, x = 0.2, family = 'Helvetica', fontface = 'bold', 
            size = 10, color = "#666666") +
  geom_text(aes(label = family_count), y = 0, x = 0, vjust = -0.9, angle = 90, family = 'Avenir',
            size = 6, color = "#666666") +
  coord_polar(theta = "y") +
  xlim(c(0.2, hsize + 0.5)) + 
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm")) +
  scale_fill_manual(values = c(brewer.pal(8, 'Set3')), name = 'Conserved amino acid') +
  theme(panel.background = element_rect(fill = "white"), panel.grid = element_blank(),
        axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank()) +
  theme(strip.background = element_blank(), strip.text.x = element_blank()) + 
  theme(text = element_text(family="Avenir", size=18, color = "#666666"))

modi <- ggplot(subset(plot_bar, type == 'Post-translational modification'), aes(x = hsize, y = proportional, fill = Var2)) +
  facet_wrap(~ y_label, ncol = 1) +
  geom_col(color = '#666666') +
  # ggtitle('Post-translational modification', ) +
  geom_text(aes(label = Var1), y = 0, x = 0.2, family = 'Helvetica', fontface = 'bold', 
            size = 10, color = "#666666") +
  # geom_text(aes(label = family_count), y = 0, x = -3, family = 'Avenir',
  #           size = 6, color = "#666666") +
  coord_polar(theta = "y") +
  xlim(c(0.2, hsize + 0.5)) + 
  scale_fill_manual(values = c('#BC80BD', '#8dde9f', '#8ddec7', '#ededed', '#8dcdde', '#FFFFB3', 'pink'), name = 'Post-translational modification') +
  theme(panel.background = element_rect(fill = "white"), panel.grid = element_blank(),
        axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank()) +
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm")) +
  theme(strip.background = element_blank(), strip.text.x = element_blank()) + 
  theme(text = element_text(family="Avenir", size=18, color = "#666666")) +
  theme(plot.title = element_text(family="Avenir", size=26, color = "#666666", hjust = 0.5))

phylo_plot <- ggplot(subset(plot_bar, type == 'Phylogenetic diversity'), aes(x = hsize, y = proportional, fill = Var2)) +
  facet_wrap(~ y_label, ncol = 1) +
  geom_col(color = '#666666') +
  # ggtitle('Phylogenetic diversity', ) +
  geom_text(aes(label = Var1), y = 0, x = 0.2, family = 'Helvetica', fontface = 'bold', 
            size = 10, color = "#666666") +
  # geom_text(aes(label = family_count), y = 0, x = -3, family = 'Avenir',
  #           size = 6, color = "#666666") +
  coord_polar(theta = "y") +
  xlim(c(0.2, hsize + 0.5)) + 
  scale_fill_manual(values = color_pal, name = 'Phylogenetic diversity') +
  theme(panel.background = element_rect(fill = "white"), panel.grid = element_blank(),
        axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank()) +
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm")) +
  theme(strip.background = element_blank(), strip.text.x = element_blank()) + 
  theme(text = element_text(family="Avenir", size=18, color = "#666666")) +
  theme(plot.title = element_text(family="Avenir", size=26, color = "#666666", hjust = 0.5))

# figure_s1 <- plot_grid(amino_acid, modi, phylo_plot, align = 'vh', nrow = 3)
#remove legends and reposition
legendz <- plot_grid(get_legend(amino_acid), get_legend(modi), get_legend(phylo_plot), align = 'h', ncol = 1)

figure_s1 <- plot_grid(amino_acid + theme(legend.position = 'NONE'), 
                       modi + theme(legend.position = 'NONE'), 
                       phylo_plot + theme(legend.position = 'NONE'),
                       legendz,
                       align = 'v', nrow = 1, rel_widths = c(1, 1, 1, 0.5))
# ggsave("figure_S1.svg", plot=figure_s1, width = 12, height = 16)
# I made <several> aesthetic changes in a figure editor, should probably spend the time to improve in ggplot instead

# If you want to zoom in on any group =================================

peek <- getMRCA(plotted_tree, tip = as.character(subset(mapping, phylum_gtdbtk %in% c('Aquificota', 'Campylobacterota'))$taxon_oid))

sub_tree <- drop.tip(phy = phylogeny_tree, 
                     tip = phylogeny_tree$tip.label[!phylogeny_tree$tip.label %in% 
                                                      get_taxa_name(tree_view = tree_plot, node = peek)])
sub_map <- subset(mapping, taxon_oid %in% sub_tree$tip.label)

peek_tree_plot <- ggtree(sub_tree, open.angle=12, color = "#666666") %<+% sub_map +
  geom_tiplab(aes(label = species), offset = 0.01) +
  geom_tippoint(aes(fill = family_gtdbtk), pch = 21, size =4, color = "#666666") +
  theme(text = element_text(family="Avenir", size=16, color = "black")) +
  guides(color = guide_legend(override.aes = list(size=5)))

trial <- as.data.frame(sub_map[, c("modification_tree")])
row.names(trial) <- sub_map$taxon_oid
new_plot <- peek_tree_plot + new_scale_fill()
gheatmap(new_plot, trial, width=0.05, offset = 0.4, hjust = 0, color=NA,
         colnames_position = 'top', colnames_offset_x = 0, colnames = FALSE, font.size = 4, 
         colnames_offset_y = -0.1, colnames_angle = 45, family = 'Avenir') +
  xlim_tree(1.8) +
  theme(text = element_text(family="Avenir", size=18, color = "black"))


# Check conservation of all significant prolines/proteins =================================

lon_stalls <- load_data('overall_stall_conservation/', '*_conserved_stall.csv')
lon_stalls$stalled <- ifelse(grepl('..PP.', lon_stalls$conserved_stall) == FALSE & lon_stalls$conserved_stall != '-----', 1, 2)
# 1 = no PP, 2 = PP present
lon_stalls <- subset(lon_stalls, genome %in% mapping$taxon_oid)
lon_stalls <- lon_stalls[order(lon_stalls$stalled, decreasing = TRUE), ]
lon_stalls$kegg <- sapply(strsplit(as.character(lon_stalls$family), '\\.'), '[', 1)
lon_stalls$desc <- kegg_key$short_description[match(lon_stalls$kegg, kegg_key$ko_number)]

stall_count <- dcast(lon_stalls, family ~ stalled, value.var = 'desc', fun.aggregate =  length)
stall_count$percent_stalled <- round(stall_count$`2` / (stall_count$`1` + stall_count$`2`), 3)
stall_count$kegg <- sapply(strsplit(as.character(stall_count$family), '\\.'), '[', 1)
stall_count$desc <- kegg_key$short_description[match(stall_count$kegg, kegg_key$ko_number)]

check <- subset(mapping, taxon_oid %in% subset(lon_stalls, family == 'K01873.02' & stalled == 2)$genome)
#check which genomes have whatever

# Getting consensus motif of non-PP stalls =================================
table(lon_stalls$family)
table(lon_stalls$desc)
fam = subset(lon_stalls, desc == 'lon')$family[1]
check <- subset(lon_stalls, family == fam & stalled == 1)
print(paste(c(nrow(check), 'genomes missing polyproline motif from', fam), collapse = ' '))
print(paste(c(nrow(subset(lon_stalls, family == fam & stalled == 2)), 'genomes contain polyproline motif from', fam), collapse = ' '))

alignment_name = paste(c('overall_stall_conservation/alignments/aligned_', fam, '.fa'), collapse = '')
look <- as.matrix(read.alignment(alignment_name, format = 'fasta'))
#fish out partial matches to genomes of interest
matching <- as.vector(check$genome)
matches <- grep(paste(matching, collapse="|"), 
                row.names(look), value = TRUE)
look <- look[row.names(look) %in% matches, ]
place <- check$position[[1]]
consensus(look[, (place-1):(place+3)])


# figure 5B =================================
kegg_list = c('K01338.01', 'K01870.01', 'K01870.02', 'K01873.01', 'K03696.01', 'K03798.01')
lon_stalls <- subset(lon_stalls, family %in% kegg_list)

#make color scheme, prolines are pink
cs1 = make_col_scheme(chars=c('P'), 
                      groups = c('proline'), 
                      cols=c('#FF0066'))
kegg_list <- factor(kegg_list, levels = c('K01873.01', 'K01870.01', 'K01870.02', 'K03798.01', 'K01338.01', 'K03696.01'))
kegg_list <- kegg_list[order(kegg_list)]

#for loop to collect geom_logo plots for every logo in figure 5, just found this the easiest way to do it
plot_list <- list()
for (i in 1:length(kegg_list)){
  fam = kegg_list[i]
  if (fam != kegg_list[[length(kegg_list)]]){
    plot_list[[i]] <- ggplot(subset(lon_stalls, family == fam)) +
      geom_logo(subset(lon_stalls, family == fam)$conserved_stall, method = 'prob', col_scheme = cs1, color = '#666666') + 
      ggtitle('') +
      ylab('') +
      xlab('') +
      theme_classic() +
      theme(plot.margin=grid::unit(c(0,0,0,0), "mm")) +
      theme(plot.title = element_blank()) +
      theme(legend.position = 'NONE') +
      theme(axis.line = element_blank()) +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
      theme(axis.text = element_text(family="Avenir", size=18, color = '#666666')) +
      theme(text = element_text(family="Avenir", size=20, color = '#666666'))}
  else{
    #only want the x-axis label on the last motif graph
    plot_list[[i]] <- ggplot(subset(lon_stalls, family == fam)) +
      geom_logo(subset(lon_stalls, family == fam)$conserved_stall, method = 'prob', col_scheme = cs1, color = '#666666') + 
      ggtitle('') +
      xlab('Polyproline motif') +
      ylab('') +
      theme_classic() +
      theme(plot.margin=grid::unit(c(0,0,0,0), "mm")) +
      theme(plot.title = element_blank()) +
      theme(legend.position = 'NONE') +
      theme(axis.line.y = element_blank()) +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
      theme(text = element_text(family="Avenir", size=28, color = "#666666"))
    
  }
  rm(fam)
}

figure_5b <- plot_grid(plotlist = plot_list, ncol = 1)


# figure 5A =================================
pfam_fams <- read.delim("pfam_clans.tsv", header = FALSE,
                        col.names = c('pfam', 'clan', 'clan_desc', 'pfam_desc', 'detailed_description'))
interpro <- load_interpro("interpro_data/", '*_interpro.tsv')
interpro$polyproline <- ifelse(interpro$protein == 'ClpC', 601,
                               ifelse(interpro$protein == 'Lon', 371, 
                                      ifelse(interpro$protein == 'FtsH', 206,
                                             ifelse(interpro$protein == 'IleS1', 56,
                                                    ifelse(interpro$protein == 'IleS2', 49, 41)))))
# locations of polyproline motif in specific proteins plotted 
# (all from Kosmotoga arenicorallina except for IleS2, which is from Mesotoga prima)

interpro$pfam_desc <- ifelse(interpro$signature_id %in% subset(pfam_fams, clan %in% c('CL0023', 'CL0671'))$pfam, 
                             pfam_fams$clan_desc[match(interpro$signature_id, pfam_fams$pfam)], interpro$pfam_desc)
interpro$pfam_desc <- gsub(' \\(DUF5915\\)', '', interpro$pfam_desc)

# set up position of y-axis of each protein
interpro$y_pos <- ifelse(interpro$protein == 'ClpC', 1,
                         ifelse(interpro$protein == 'Lon', 1.015, 
                                ifelse(interpro$protein == 'FtsH', 1.03,
                                       ifelse(interpro$protein == 'IleS1', 1.06,
                                              ifelse(interpro$protein == 'IleS2', 1.045, 1.075)))))
interpro <- interpro[order(interpro$start), ]

protein_domains <- ggplot(interpro) +
  geom_rect(aes(xmax = protein_length, ymin=y_pos - 0.003, ymax=y_pos + 0.003), xmin = 1, color = '#716e8c', fill = '#9a98ae') +
  geom_rect(data = subset(interpro, analysis == 'Pfam'), aes(xmin = start, xmax = end, fill = pfam_desc, ymin=y_pos - 0.004, ymax=y_pos + 0.004), color = '#716e8c') +
  geom_rect(data = subset(interpro, interpro_id == 'IPR001270' | signature_id == 'PIRSR001174-2.1'), aes(xmin = start, xmax = end, ymin=y_pos - 0.004, ymax=y_pos + 0.004),
            color = '#716e8c', fill = 'white') +
  geom_rect(data = subset(interpro, interpro_id == 'IPR001412' | signature_desc == 'Active site'), aes(xmin = start, xmax = end, ymin=y_pos - 0.004, ymax=y_pos + 0.004),
            color = '#716e8c', fill = 'white') +
  geom_point(aes(x = polyproline, y = y_pos), pch = 4, stroke = 3, size = 3, color = 'light grey') +
  geom_point(aes(x = polyproline, y = y_pos), pch = 4, stroke = 2, size = 3, color = '#FF0066') +
  geom_text(aes(label = protein, x = -90, y = y_pos), family="Avenir", size = 11, color = "#666666") +
  xlab("Amino acid number") +
  ylab("") +  
  xlim(-110, 1050) +
  scale_fill_manual(values = c(brewer.pal(11, 'Spectral'), '#40366f')) +
  theme_classic() +
  theme(legend.position = 'bottom', legend.title = element_blank()) +
  theme(text = element_text(family="Avenir", size=28, color = "#666666")) +
  theme(axis.text = element_text(family="Avenir", size=28, color = "#666666")) + 
  theme(legend.text = element_text(family="Avenir", size=24, color = "#666666")) + 
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  guides(fill = guide_legend(nrow = 4))

legendz <- get_legend(protein_domains)
first_plotz <- plot_grid(protein_domains + theme(legend.position = 'NONE'), figure_5b, align = 'h', axis = 'bt', nrow = 1, rel_widths = c(2, 0.75))
# combine figure 5A and 5B
figure_5 <- plot_grid(first_plotz, legendz, nrow = 2, rel_heights = c(1, 0.20))
# put legend on the bottom of combined figure

# ggsave("figure_5.svg", plot=figure_5, width = 18, height = 12)
