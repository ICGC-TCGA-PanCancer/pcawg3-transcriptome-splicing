# Code relating to analysis of alternative splicing events - ICGC / PCAWG-3
This public repository contains the research and analysis code that was used to generate the results
presented in the publication "Genomic basis for RNA alterations in cancer."

The code is organized across several directoris:

* **`event_calling`**
The scripts used to call the alternative splicing events from aligned RNA-Seq data
* **`exonizations`**
Code related to detecting exonization events from the splicing graph data
* **`filter_events`**
Code filtering for potentially functional events (used for gene centric analysis)
* **`gene_centric_tables`** 
Code summarizing the alternative splicing events of all samples into a splice-outlier table
* **`icgc_anno`, `icgc_colors`, `icgc_utils`**
Helper and utility code to take care of logistics and metadata
* **`junction_db`**
Scripts to collect all aligned splice junctions into a joint database
* **`pca_tsne`**
Scripts using dimensionality reduction for visualization
* **`stats_events`**
Code summarizing the event statistics
* **`utils`**
Utility and helper functions
