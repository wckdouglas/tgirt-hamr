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

#fitControl <- trainControl(method = 'repeatedcv',
#	repeats = 10,
#	number = 6,
#	allowParallel = T)
fitControl <- trainControl(method = 'LOOCV',
	allowParallel = T)

binomTest <- function(a,b,seqErr){
	return(pbinom(b-a,b, p = 1- seqErr, lower.tail=T))
}

modeling <- function(base,df, model){
	message('start modeling..')
    columns <- c('A','C','T','G','deletion','label')
    columns <- columns[!grepl(base, columns)]
	df <- df[df$ref==base,columns]
	#split data
	trainMat <- select(df, -label)
	trainClass <- factor(df$label)
	if (length(unique(trainClass)) > 1){
		message ('Start training ',model,' ',base,'\n')
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
    if (nrow(df)>2){
        columns <- c('A','C','T','G','deletion')
        columns <- columns[!grepl(base,columns)]
	    dataMat <- df[df$ref==base,columns]
	    df$label <- predict(model$finalModel,dataMat,type = 'class')
		return (df)
	}else{
		return(0)
	}
}

fitAndPredict <- function(base,dataTable,predictionTable,model){
	modelFit <- modeling(base,dataTable, model)
	message('Tuned model')
	if(!is.character(modelFit)){
		df <- predicting(base,modelFit,predictionTable)
		return(df)
	}else{
		return (df %>% mutate(label = modelFit))
	}
}

filterSets <-function(df,hyp){
	if (hyp == 'hyp1') { df = subset(df,padj1==1)}
	else if (hyp == 'hyp2') { df = subset(df,padj2==1)}
    return(df)
}

main <- function(predictTable,model,enzyme,seqErr,pCutOff,resultFile,hyp,dbpath) {
	dataTable <- str_c(dbpath,'/',enzyme,'Table.tsv') %>%
		read_tsv(col_type= 'cncnnnnnnncc') %>%
		transformDF(seqErr,pCutOff,binomTest) %>%
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
		transformDF(seqErr,pCutOff,binomTest) %>%
		filterSets(hyp) 
	message('Read Data!')
	
	bases = c('A','C','T','G')
	tablename <- resultFile
	result <- lapply(bases,fitAndPredict,dataTable,predictTable,model) 
	result <- result[sapply(result,function(x) is.data.frame(x))] %>%
		do.call(rbind,.) 

    if (length(result)<1){stop("No modification sites! \n")}
    else{
        result %>%
#	        select(chrom, start, end, ref, cov, strand, A, C, T, G, deletion, label) %>%
	        select(chrom, start, end, ref, cov, strand, label) %>%
		    write.table(tablename, sep='\t',quote=F,row.names=F,col.names=F)
    }
	return (0)
}

#=============      set variables         ===========
opt <- getopt(matrix(c(
	'model',	 'm', 1,"character",
	'resultFile','o', 1, "character",
	'threads',   't', 1, "numeric",
	'dt',        'i', 1, 'character',
	'enzyme',    'e', 1, 'character',
	'hyp',       'h', 1, 'character',
	'seqErr',    's', 1, 'numeric',
	'pCutOff',   'p', 1, 'numeric',
	'dbpath',    'd', 1, 'character',
	'dir',       'f',2,'character'),
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

#============== run program ========================
message('Using: ',model,' for ', enzyme)
main(predictTable,model,enzyme,seqErr,pCutOff,resultFile,hyp,dbpath)
message('Finished: ',model,' for ',enzyme)
