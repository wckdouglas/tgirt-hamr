#!/bin/env Rscript

suppressMessages(library(getopt))
suppressMessages(library(doMC))

#=============      set variables         ===========
opt <- getopt(matrix(c(
    'model',     'm', 1,"character",
    'resultFile','o', 1, "character",
    'threads',   't', 2, "numeric",
    'dt',        'i', 1, 'character',
    'enzyme',    'e', 1, 'character',
    'hyp',       'h', 1, 'character',
    'seqErr',    's', 1, 'numeric',
    'pCutOff',   'p', 1, 'numeric',
    'devMode',   'v',2, 'numeric'),
    byrow=T,ncol=4))
model <- opt$model
resultFile <- opt$resultFile
registerDoMC(cores = opt$threads)
predictTable <- opt$dt 
hyp <- opt$hyp
enzyme <- opt$enzyme
seqErr <- opt$seqErr
pCutOff <- opt$pCutOff

if (!is.null(opt$devMode)){
    devMode = 1
}else{
    devMode=0
}

#============== run program ========================
message('Using: ',model,' for ', enzyme)
tgirthamr::tgirthamr(predictTable,model,enzyme,seqErr,pCutOff,resultFile,hyp,devMode)
message('Finished: ',model,' for ',enzyme)
