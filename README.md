
<h1 align="center">
  <br>
  <br>
  ELViS
  <br>
</h1>

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

```
# Install BiocManager if not already
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

# Install devtools if not already
if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
}

# install release version
## Ensure repos includes both CRAN and Bioconductor repositories
if( BiocManager::version() >= "3.19" ){
   BiocManager::install("ELViS")
}else{
   options(repos = BiocManager::repositories())
   devtools::install_github("https://github.com/hyochoi/ELViS.git",ref="master_bioc_le_3.18 ")
}

# install development version
## Ensure repos includes both CRAN and Bioconductor repositories
options(repos = BiocManager::repositories())
if( BiocManager::version() >= "3.19" ){
   devtools::install_github("https://github.com/hyochoi/ELViS.git",ref="devel")
}else{
   devtools::install_github("https://github.com/hyochoi/ELViS.git",ref="devel_bioc_le_3.18")
}

```


## HOW TO USE

Please refer to vignette.
