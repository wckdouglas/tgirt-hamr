#!/bin/env Rscript

suppressMessages(library(readr))
suppressMessages(library(tidyr))
suppressMessages(library(stringr))
suppressMessages(library(Rcpp))
suppressMessages(library(getopt))
suppressMessages(library(caret))
suppressMessages(library(doMC))
suppressMessages(library(kernlab))
suppressMessages(library(dplyr))
suppressMessages(library(tgirthamr))

fitControl <- trainControl(method = 'LOOCV',
    allowParallel = T)

modeling <- function(base,df, model){
    message('start modeling..')
    columns <- c('A','C','T','G','deletion','label')
    columns <- columns[!grepl(base, columns)]
    df <- df[df$ref==base,columns]
    #split data
    trainMat <- select(df, -label)
    trainClass <- factor(df$label)
    if (length(unique(trainClass)) > 1){
        message ('Start training ',model,' ',base)
        modelFit <- train(y = trainClass, x = trainMat, method = model, trControl = fitControl)
        message('Trained ',model, ' for ',base)
        return(modelFit)
    }else{
        message ('Skipped ',model,' ',base)
        return(unique(trainClass))
    }
}

predicting <- function(base,model,df){
    #subsetting data
    message ('start prediction')
    df <- subset(df,ref == base )
    columns <- c('A','C','T','G','deletion')
    columns <- columns[!grepl(base,columns)]
    dataMat <- df[df$ref==base,columns]
    df$label <- predict(model$finalModel,dataMat,type = 'class')
    return (df)
}

fitAndPredict <- function(base,dataTable,predictTable,model){
    modelFit <- modeling(base,dataTable, model)
    message('Tuned model')
    if(!is.character(modelFit)){
        return(predicting(base,modelFit,predictTable))
    }else{
        return (mutate(df,label = modelFit))
    }
}

filterSets <-function(df,hyp){
    if (hyp == 'hyp1') { df = subset(df,padj1==1)}
    else if (hyp == 'hyp2') { df = subset(df,padj2==1)}
    return(df)
}

main <- function(predictTable,model,enzyme,seqErr,
                pCutOff,resultFile,hyp,dbpath,devMode) {
    dataTable <- str_c(dbpath,'/',enzyme,'Table.tsv') %>%
        read_tsv(col_type= 'cncnnnnnnncc') %>%
        transformPredict(seqErr,pCutOff,binomTest) %>%
        mutate(label =  mergeType(as.character(abbrev)))  %>%
        filterSets(hyp) %>%
        group_by(label) %>% 
        do(data.frame(count = nrow(.),
                    A = .$A, 
                    C = .$C, 
                    T=.$T,G=.$G,
                    ref = .$ref,
                    deletion = .$deletion)) %>%
        filter(count > 2) %>%
        select(-count) %>%
        tbl_df

    predictTable <- predictTable %>%
        read_tsv %>%
        transformPredict(seqErr,pCutOff,binomTest) %>%
        filterSets(hyp) 
    message('Read Data!')
    
    bases = as.character(unique(predictTable$ref))
    tablename <- resultFile
    result <- lapply(bases,fitAndPredict,dataTable,predictTable,model)  %>%
        do.call(rbind,.) 

    if (length(result)<1){stop("No modification sites! \n")}
    else{
        if (devMode==1){
            result %>%
                select(chrom, start, end, ref, cov, strand, A, C, T, G, deletion, label) %>%
                write.table(tablename, sep='\t',quote=F,row.names=F,col.names=F)
        }else{
            result %>%
                select(chrom, start, end, ref, cov, strand, label) %>%
                write.table(tablename, sep='\t',quote=F,row.names=F,col.names=F)
        }
    }
}

#=============      set variables         ===========
opt <- getopt(matrix(c(
    'model',     'm', 1,"character",
    'resultFile','o', 1, "character",
    'threads',   't', 1, "numeric",
    'dt',        'i', 1, 'character',
    'enzyme',    'e', 1, 'character',
    'hyp',       'h', 1, 'character',
    'seqErr',    's', 1, 'numeric',
    'pCutOff',   'p', 1, 'numeric',
    'dbpath',    'd', 1, 'character',
    'dir',       'f',2,'character',
    'devMode',   'v',2, 'numeric'),
    byrow=T,ncol=4))
model <- opt$model
resultFile <- opt$resultFile
registerDoMC(cores = opt$threads)
predictTable <- opt$dt 
hyp <- opt$hyp
enzyme <- opt$enzyme
path <- opt$dir

seqErr <- opt$seqErr
pCutOff <- opt$pCutOff
dbpath <- opt$dbpath
sourceCpp(str_c(path,'function.cpp',sep='/'))
if (!is.null(opt$devMode)){
    devMode = 1
}else{
    devMode=0
}

#============== run program ========================
message('Using: ',model,' for ', enzyme)
main(predictTable,model,enzyme,seqErr,pCutOff,resultFile,hyp,dbpath,devMode)
message('Finished: ',model,' for ',enzyme)
