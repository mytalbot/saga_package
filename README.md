
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SAGA <img src="https://talbotsr.com/saga_package/logo.png" align="right" height="139" />

### Surrogate Assay for Genotoxicity Assessment

In our paper **“Predicting genotoxicity of integrating viral vectors for
stem cell gene therapy using gene expression-based machine learning”**
we report an improved in vitro test to determine the risk of insertional
mutagenesis of integrating vectors for gene therapy. SAGA builds on the
well accepted cell culture protocol of the in vitro immortalization
assay (IVIM) but screens for the deregulation of oncogenic gene
expression signatures. We demonstrate a new bioinformatic approach to
correctly classify the mutagenic potential of retroviral vectors used in
previous and current clinical trials.

## Installation

The package requires **R\>=3.6**.

You can download the development version from
[GitHub](https://github.com/mytalbot/saga_package/) or use:

``` r
# install.packages("devtools")
devtools::install_github("mytalbot/saga_package", build_vignettes = TRUE)
library(saga)
```

If you want the Vignette installed in R as well, set the
build\_vignettes = TRUE (access browseVignettes(“saga”)). Otherwise, the
Vignette will only be accessible via the website.

If you don’t have or want to use devtools, SAGA can be installed using
the source file (i.e. in RStudio). [It can be downloaded
here.](https://github.com/mytalbot/saga_package/tree/master/sourcefiles)

### Important Note

Please note, that the SAGA package requires dependencies. Some functions
can mask each other (depending on the local R setup). This may cause
**warnings** (not errors). These can be ignored since they do not hamper
the function of the SAGA package. Warnings may be caused mainly by the
HTSanalyzeR and mnormt dependencies.

**Since the mnormt package has not yet been updated to R4.0**, the only
solution to work with saga on R4.0 is the installation of a previous
version of mnormt. Make sure not to update mnormt when installing saga
with
devtools\!

``` r
PackageUrl <- "https://cran.r-project.org/src/contrib/Archive/mnormt/mnormt_1.5-7.tar.gz"
install.packages(PackageUrl, repos = NULL,type="source")
```

### SAGA core data

You’ll also need the SAGA core data. They are too large for the main
package. Download them here (also from GitHub):

``` r
devtools::install_github("mytalbot/sagadata")
library(sagadata)
```

## Examples

The following working examples demonstrate the general workflow of SAGA
analysis.

You can download an example data set of four arrays plus a ready-to-run
SampleInformation.txt file from the GitHub repository. Just download the
files to a folder on your machine, unzip them and specify the path (see
below) to the selected folder in the script below. The sample arrays had
to be zipped because of the GitHub size limitations for files.

Note: Make sure that the phenoTest package is sourced as
library(phenoTest) before running the script. This is required for GESEA
analysis.

### Using the wrapper function

``` r
library(phenoTest) # this is mandatory (for GESEA)!
library(saga)
library(sagadata)

path           <- "path to your samples"

### saga_gentargets; automatically generates an empty (!) sample information file
# Modify the columns: Filename, Batch, Group, Vector, TrueLabel
# targets      <- saga_gentargets(smplpath=path)

################################################################################
### Wrapper function - all in one
### to use: uncomment function & execute
################################################################################
mySAGAres      <- saga_wrapper(samplepath, showModel=0, doGESEA=1)
```

### More detailed SAGA analysis

``` r
library(phenoTest) # this is mandatory!
library(saga)
library(sagadata)
### Using the wrapper function
################################################################################
### Or use these single functions
################################################################################
### saga_import
rawdata        <- saga_import(smplpath=path, showjoint=1)
SAGA_RAW       <- rawdata$SAGA_RAW
TEST_RAW       <- rawdata$TEST_RAW
pData.joint    <- rawdata$pData.joint
pData.Test     <- rawdata$pData.Test

### normalize saga data (plotnumber=1 for normalized boxplot, 2 for tSNE plot)
normalized     <- saga_norm(SAGA_RAW, TEST_RAW, pData.joint, plotnumber=1)
qunorm.SAGA    <- normalized$qunorm.SAGA
matrix.SAGA.qn <- normalized$matrix.SAGA.qn
matrix.test.qn <- normalized$matrix.test.qn

### remove batch effects
# plotnumber = 1 (tSNE of first 2 dimensions)
# plotnumber = 2 (PCA of first 2 dimensions)
batchnorm      <- saga_batch(matrix.SAGA.qn, matrix.test.qn, rawdata, pData.joint, plotnumber=1)
index          <- batchnorm$index
matrix.SAGA    <- batchnorm$matrix.SAGA
matrix.test    <- batchnorm$matrix.test

### Sampling and model data collection
model          <- saga_sampling(matrix.SAGA, matrix.test)
matrix.train   <- model$matrix.train
labels.train   <- model$labels.train
matrix.unknown <- model$matrix.unknown

### Array predictions with optimized SVM parameters (default settings)
output         <- saga_predict(samplepath, matrix.train, labels.train, matrix.unknown,
                               pData.Test, writeFile=1, showRoc=0)
classes        <- output$predictions
classes

### GESEA
gesea_results  <- saga_gesea(smplpath=path, saveResults=0)
```
