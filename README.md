
<h1 align="center">
  <br>
  <br>
  ELViS
  <br>
</h1>

<h4 align="center">An R Package for <u><b>E</b></u>stimating Copy Number <u><b>L</b></u>evels of <u><b>Vi</b></u>ral Genome <u><b>S</b></u>egments Precisely at Base-Resolution.</h4>

<p align="center">
  <a href="#installation">INSTALLATION</a> •
  <a href="#how-to-use">How To Use</a> •
</p>

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
'
# Ensure repos includes both CRAN and Bioconductor repositories
options(repos = BiocManager::repositories())
devtools::install_github("https://github.com/hyochoi/ELViS.git")
```


## How To Use

To clone and run this application, you'll need [Git](https://git-scm.com) and [Node.js](https://nodejs.org/en/download/) (which comes with [npm](http://npmjs.com)) installed on your computer. From your command line:
