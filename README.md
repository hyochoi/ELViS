
<h1 align="center">
  <br>
  <br>
  ELViS
  <br>
</h1> 

<sup><b>R4.3 & Bioc3.18 (Ubuntu & Windows) : </b></sup> [![R-CMD-check-bioc](https://github.com/hyochoi/ELViS/actions/workflows/check-bioc_3.18.yml/badge.svg)](https://github.com/hyochoi/ELViS/actions/workflows/check-bioc_3.18.yml)
<br>
<sup><b>R4.4 & Bioc3.19 (Ubuntu & Windows) : </b></sup> [![R-CMD-check-bioc](https://github.com/hyochoi/ELViS/actions/workflows/check_bioc.yml/badge.svg)](https://github.com/hyochoi/ELViS/actions/workflows/check_bioc.yml)

---------------------------

An R Package for ***E***stimating Copy Number ***L***evels of ***Vi***ral Genome ***S***egments Using Base-Resolution Read Depth Profile.

<p align="center">
  <a href="#description">DESCRIPTION</a> •
  <a href="#installation">INSTALLATION</a> •
  <a href="#how-to-use">HOW TO USE</a> •
</p>

## DESCRIPTION
* Base-resolution copy number analysis of viral genome. Utilizes base-resolution read depth data over viral genome to find copy number segments with two-dimensional segmentation approach.

* Provides publish-ready figures, including
    - histograms of read depths
    - coverage line plots over viral genome annotated with copy number change events and viral genes
    - heatmaps showing multiple types of data with integrative clustering of samples.


## INSTALLATION

To install this package, start R (version "4.5") and enter:

```{r , echo=TRUE, eval=FALSE}
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ELViS")
```


## HOW TO USE

Please refer to [vignette](https://jyleebioinfo.github.io/ELViS/).
