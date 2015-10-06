#!/bin/env Rscript

suppressMessages(library(data.table))
suppressMessages(library(tidyr))
suppressMessages(library(stringr))
suppressMessages(library(Rcpp))
suppressMessages(library(getopt))
suppressMessages(library(caret))
suppressMessages(library(doMC))
suppressMessages(library(kernlab))
suppressMessages(library(dplyr))

fitControl <- trainControl(method = 'repeatedcv',
	repeats = 10,
	number = 6,
	allowParallel = T)

binomTest <- function(a,b,seqErr){
	return(pbinom(b-a,b, p = 1- seqErr, lower.tail=T))
}

filterBase <- function(dataframe,base){
	#subsetting data
	dataframe <- dataframe %>% filter(ref == base )
	if (base == 'A'){
		df <- subset(dataframe,select=c(label,C,T,G,deletion))
	}else if(base == 'C'){
		df <- subset(dataframe,select=c(label,A,T,G,deletion))
	}else if(base == 'T'){
		df <- subset(dataframe,select=c(label,A,C,G,deletion))
	}else if(base == 'G'){
		df <- subset(dataframe,select=c(label,A,T,C,deletion))
	}
	return (df)
}

modeling <- function(base,df, model){
	message('start modeling..')
	df <- df %>% filterBase(base)
	#split data
	trainMat <- select(df, -label)
	trainClass <- factor(df$label)
	if (length(unique(trainClass)) > 1){
		message ('Start training',model,base,'\n')
		modelFit <- train(y = trainClass, 
				x = trainMat,
				method = model,
				trControl = fitControl)
        message('Trained ',model, ' for ',base)
		return(modelFit)
	}else{
		message ('Skipped ',model,' ',base)
	}
}


predicting <- function(base,model,df){
	#subsetting data
    message ('start prediction')
	df <- filter(df,ref == base )
	if (nrow(df) != 0){
		if (base == 'A'){
			dataMat <- subset(df,select=c(C,T,G,deletion))
		}else if(base == 'C'){
			dataMat <- subset(df,select=c(A,T,G,deletion))
		}else if(base == 'T'){
			dataMat <- subset(df,select=c(A,C,G,deletion))
		}else if(base == 'G'){
			dataMat <- subset(df,select=c(A,T,C,deletion))
		}
		df$label <- predict(model$finalModel,dataMat,type = 'class')
		return (df)
	}else{
		return(0)
	}
}

fitAndPredict <- function(base,dataTable,predictionTable,model){
	modelFit <- modeling(base,dataTable, model)
	message('Established model')
	if(!is.null(modelFit)){
		df <- predicting(base,modelFit,predictionTable)
		return(df)
	}else{
		return (0)
	}
}

filterSets <-function(df,hyp){
	if (hyp == 'hyp1') { df = df %>% filter(padj1==1)}
	else if (hyp == 'hyp2') { df = df %>% filter(padj2==1)}
    return(df)
}

main <- function(predictTable,model,enzyme,seqErr,pCutOff,resultFile,hyp,dbpath) {
	dataTable <- paste0(dbpath,
					'/',
					enzyme,
					'Table.tsv')
	dataTable <- dataTable %>%
		fread(colClasses=c('character','numeric','character',
						   rep('numeric',7),rep('character',2))) %>%
		transformDF(seqErr,pCutOff,binomTest) %>%
		mutate(label =  mergeType(as.character(abbrev)))  %>%
		filterSets(hyp) %>%
        tbl_df()
	dataTable <- dataTable %>%
		group_by(label) %>%
		summarize(count = n()) %>%
		filter(count > 6) %>%
		inner_join(dataTable)

	predictTable <- predictTable %>%
		fread() %>%
		transformDF(seqErr,pCutOff,binomTest) %>%
		filterSets(hyp) 
	message('Read Data!')
	
	bases = c('A','C','T','G')
	tablename <- resultFile
	result <- lapply(bases,fitAndPredict,dataTable,predictTable,model) 
	result <- result[sapply(result,function(x) is.data.frame(x))] %>%
				do.call(rbind,.) %>%
				select(chrom, start, end, ref, cov, strand, A, C, T, G, deletion, label) %>%
				write.table(tablename, sep='\t',quote=F,row.names=F,col.names=F)
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
sourceCpp(paste(path,'function.cpp',sep='/'))

#============== run program ========================
message('Using: ',model,' for ', enzyme)
main(predictTable,model,enzyme,seqErr,pCutOff,resultFile,hyp,dbpath)
message('Finished: ',model,' for ',enzyme)
