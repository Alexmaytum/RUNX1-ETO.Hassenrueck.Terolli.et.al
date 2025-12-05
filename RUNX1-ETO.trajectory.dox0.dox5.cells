# Script written by Dr Peter Keane


library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(ggplot2)
library(tidyverse)

#############################################################################################################################

setwd("")

#############################################################################################################################

# Load cells
cells<- readRDS(file = 'KT29Subset.rds') %>% subset(dox %in% c('dox0', 'dox5'))

# Convert to a Monocle3 Object
cds<- as.cell_data_set(x = cells)

# Cluster cells and calculate trajectories
cds<- cluster_cells(cds = cds)
cds<- learn_graph(cds = cds, close_loop = FALSE, use_partition = FALSE)

# Plot trajectory
plot_cells(cds = cds, 
           trajectory_graph_segment_size = 2, 
           label_branch_points = FALSE, 
           label_leaves = FALSE, 
           label_cell_groups = FALSE, 
           cell_size = 0.05, 
           cell_stroke = 0.75, 
           color_cells_by = 'seurat_clusters')

#############################################################################################################################

plot_cells(cds = cds, 
           trajectory_graph_segment_size = 2, 
           label_branch_points = FALSE, 
           label_leaves = FALSE, 
           label_cell_groups = FALSE, 
           label_principal_points = TRUE, 
           cell_size = 0.05, 
           cell_stroke = 0.75, 
           color_cells_by = 'seurat_clusters')

#############################################################################################################################

# Pseudotime
cds<- order_cells(cds = cds, root_pr_nodes = c('Y_5', 'Y_27'))

# Plot trajectory
plot_cells(cds = cds, 
           trajectory_graph_segment_size = 2, 
           label_branch_points = FALSE, 
           label_leaves = FALSE, 
           label_cell_groups = FALSE, 
           cell_size = 0.05, 
           cell_stroke = 0.75, 
           color_cells_by = 'pseudotime', 
           label_principal_points = FALSE, 
           label_roots = FALSE)

#############################################################################################################################

# Make publication ready figures

# Plot trajectory
g<- plot_cells(cds = cds, 
               trajectory_graph_segment_size = 2, 
               label_branch_points = FALSE, 
               label_leaves = FALSE, 
               label_cell_groups = FALSE, 
               color_cells_by = 'seurat_clusters', 
               label_roots = FALSE,
               cell_size = 0.5) 

g$layers[[1]][['aes_params']]$alpha<- 0

g<- g + 
  theme_classic() +
  theme(axis.text = element_text(size = 18), axis.title = element_text(size = 18)) +
  theme(axis.line = element_line(linewidth = 1)) +
  theme(axis.ticks = element_line(linewidth = 1), axis.ticks.length = unit(0.25, 'cm')) +
  scale_x_continuous(limits = c(-10.5,10.5), breaks = seq(-10,10,5)) +
  scale_y_continuous(limits = c(-10.5,10.5), breaks = seq(-10,10,5))


# Plot pseudotime
g<- plot_cells(cds = cds, 
               trajectory_graph_segment_size = 2, 
               trajectory_graph_color = 'grey50',
               label_branch_points = FALSE, 
               label_leaves = FALSE, 
               label_cell_groups = FALSE, 
               color_cells_by = 'pseudotime', 
               label_principal_points = FALSE, 
               label_roots = FALSE,
               cell_size= 0.5)

g$layers[[1]][['aes_params']]$alpha<- 0

g<- g + 
  theme_classic() +
  theme(axis.text = element_text(size = 24), axis.title = element_text(size = 24)) +
  theme(axis.line = element_line(linewidth = 1)) +
  theme(axis.ticks = element_line(linewidth = 1), axis.ticks.length = unit(0.25, 'cm')) +
  theme(legend.key.height = unit(1.25, 'cm'), legend.text = element_text(size = 14),) +
  theme(legend.title.position = 'left', legend.title = element_text(size = 24, angle = 90, hjust = 0.5)) +
  scale_x_continuous(limits = c(-10.5,10.5), breaks = seq(-10,10,5)) +
  scale_y_continuous(limits = c(-10.5,10.5), breaks = seq(-10,10,5)) 

ggsave(g, file = 'pseudotime_results/Dox0_and_Dox5_Pseudotime.png', dpi = 800, height = 8, width = 9.5)

#############################################################################################################################

# Save Monocle3 object for later
saveRDS(cds, file = 'saved_R_objects/Dox0_and_Dox5_Monocle3.rds')

#############################################################################################################################
