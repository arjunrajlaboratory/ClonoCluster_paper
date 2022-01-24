# BarCluster_paper

## Intro

Raw data and analysis scripts for the *BarCluster* paper. Includes additional figures generated and analyses that were not present in the publication.

For the actual distributed open source software, including worked examples, check out [this repository](https://github.com/leeprichman/BarCluster).

## Get the raw data

Due to size limitations, the raw data for this study cannot be directy hosted on GitHub. I have compressed it and uploaded it with the release of this package. You will need to download and install [`pixz`](https://github.com/vasi/pixz) to uncompress this package.

```
pixz -d BarCluster_raw_data.tpxz
```

Then untar the package.

```
tar -xvf BarCluster_raw_data.tar
```

Now copy the two folders into the repository.

```
cp Data_genes ./BarCluster_paper/
```



## Repo contents

  * `Data_barcodes/` - barcode assignments for each dataset

  * `Data_genes/` - normalized/scaled count matrices for each dataset

  * `Figs/` - Output figures from analysis scripts, divided by subfolder for each script.

  * `Scripts/` - R scripts to be sourced that generate intermediate processed data and figures, see below for walkthrough.

## Analysis order

  1. Establish variables with sample names and alpha analysis values

  ```
  source("Constants.R")
  ```

  2. Get clusters across range of alphas from 0 to 1 and save intermediate output. Plot long form alluvia.

  ```
source("Long_alluvia.R")
  ```

  3. Generate plots of cluster sizes.

  ```
  source("Cluster_size_plots.R")
  ```

  4. Generate Cohen's kappa plots.

  ```
  source("Confusion_stats.R")
  ```

  5. Generate short alluvia
  ```
  source("Short_alluvia.R")
  ```
  6. Perform marker identification for cluster and reorganization markers as well as gene set overrepresentation analysis
  ```
  source("Marker_analysis.R")
  ```

  7. Generate marker top cluster marker fidelity boxplots and venn diagrams of all marker overlap
  ```
  source("Marker_box_and_venn.R")
  ```

  8. Generate heatmaps of marker fidelity for top cluster markers
  ```
  source("Marker_turnover_hm.R")
  ```

  9. Generate alluvia of interesting markers.
  ```
  source("Marker_sankey.R")
  ```

  10. Generate sample Sankeys for the reorganization analysis
  ```
  source("Reorg_sankey.R")
  ```

  11. Generate reorganization marker overrepresentation analysis heat maps
  ```
  source("geneset_hm.R")
  ```

  12. Generate warped UMAPs from sample data and real data, some saved in the warped/ folder
  ```
  source("Fig4.R")
  source("Fig4_clusters.R")
  ```

  13. Generate supplemental analyses showing how BarCluster works
  ```
  source("Model_edge_weights.R") # show curves for how model influences edge weight at beta = 0.1
  source("grid_graph.R") # simulation of network graphs with alpha
  ```

## Citation

[biorxiv link](#intro)

## Contact

Open an issue on this GitHub repository, contact @leeprichman, @arjunrajlaboratory, or email [myself](mailto:leeprichman@gmail.com) or [Dr. Arjun Raj](arjunrajlaboratory@gmail.com).
