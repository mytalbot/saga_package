######################################
# SAGA Script
# v.1.0.1
# by the SAGA team
######################################

library(saga)

samplepath     <- "C:/MHH Bleich/Rothe/SAGA/saga_samples"

### saga_gentargets; automatically generates an empty (!) sample information file
# Modify manually according to your requirements
#targets        <- saga_gentargets(samplepath)

################################################################################
### Wrapper function - all in one
### to use: uncomment function & execute
################################################################################
mySAGAres      <- saga_wrapper(samplepath, showModel=0, doGESEA=1)


################################################################################
### Or use these single functions
################################################################################
### saga_import
rawdata        <- saga_import(samplepath, showjoint=1)
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
gesea_results  <- saga_gesea(samplepath, saveResults=1)



#############################
######## HELP Files ########
############################
# to see the help overview
help(package = "saga", help_type = "html")

# to see the vignette
vignette("saga_vignette")



