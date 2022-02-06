
#  1. Set working directory.
setwd("BarCluster_paper")


#2. Establish needed variables with sample names and alpha analysis values.
source("Paper/extractionScripts/Constants.R")

#3. Get clusters across range of alphas from 0 to 1 and save intermediate output in extracted data folder. Plot long Sankey as in Figure 1C and Figure S2A.
source("Paper/extractionScripts/Long_alluvia.R")

#4. Perform marker identification for cluster and reorganization markers as well as gene set overrepresentation analysis, saves to extracted data folder.
source("Paper/extractionScripts/Marker_analysis.R")

#5. Generate plots of cluster sizes as in Figure 1D and Figure S3.
source("Paper/plotScripts/Cluster_size_plots.R")

#6. Generate Cohen's kappa plots as in Figure S2B.
source("Paper/plotScripts/Confusion_stats.R")

#7. Generate short Sankeys as in Figure 1E.
source("Paper/plotScripts/Short_alluvia.R")

#8. Generate top cluster marker fidelity boxplots (Figure 2A and Figure S5B) and venn diagrams of all marker overlap (Figure S4)
source("Paper/plotScripts/Marker_box_and_venn.R")

#9. Generate heatmaps of marker fidelity for top cluster markers as in Figure S5A
source("Paper/plotScripts/Marker_turnover_hm.R")

#10. Generate alluvia of interesting markers as in Figure 2B-D.
source("Paper/plotScripts/Marker_sankey.R")

#11. Generate sample Sankeys for the reorganization analysis (Figure 3B).
source("Paper/plotScripts/Reorg_sankey.R")

#12. Generate reorganization marker overrepresentation analysis heat maps (Figure 3C).
source("Paper/plotScripts/geneset_hm.R")

#13. Generate warped UMAPs from sample data and real data (Figure 4).
source("Paper/plotScripts/Fig4.R")
source("Paper/plotScripts/Fig4_clusters.R")

#14. Generate UMAP and Sankeys as in Figure 5, showing the combined hybrid clustering and warp factor.
source("Paper/plotScripts/Labeled_umaps.R")

#15. Generate supplemental analyses showing how BarCluster works (Figure S1).
  source("Paper/plotScripts/Model_edge_weights.R") # show curves for how model influences edge weight at beta = 0.1
  source("Paper/plotScripts/grid_graph.R") # simulation of network graphs with alpha
