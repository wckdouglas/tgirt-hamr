#!/usr/bin/env Rscript

suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(stringr))
suppressMessages(library(Rcpp))
suppressMessages(library(getopt))
suppressMessages(library(caret))
suppressMessages(library(doMC))
suppressMessages(library(kernlab))

fitControl <- trainControl(method = 'repeatedcv',
	repeats = 10,
	number = 6,
	allowParallel = T)

binomTest <- function(a,b){
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
	cat('start modeling\n')
	df <- df %>% filterBase(base)
	#split data
	trainMat <- select(df, -label)
	trainClass <- factor(df$label)
	if (length(unique(trainClass)) > 1){
		cat ('Start training',model,base,'\n')
		modelFit <- train(y = trainClass, 
				x = trainMat,
				method = model,
				trControl = fitControl)
		return(modelFit)
	}else{
		cat ('Skipped',model,base,'\n')
	}
}


predicting <- function(base,model,df){
	#subsetting data
	cat('start prediction\n')
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
	cat('Established model\n')
	if(!is.null(modelFit)){
		df <- predicting(base,modelFit,predictionTable)
		return(df)
	}else{
		return (0)
	}
}

transformTable <- function(df,seqErr,pCutOff){
	df <- df %>%
		mutate(mismatch = A + C + T +G,
			A = A / cov,
			T = T / cov,
			G = G / cov,
			C = C / cov,
			deletion = deletion / cov,
			p1 = mapply(binomTest, mismatch,cov),
			padj1 = FDRcontrol(p1,pCutOff),
			het = heterozygote(A,C,T,G,cov),
			p2 = mapply(binomTest,het,cov),
			padj2 = FDRcontrol(p2,pCutOff))  %>%
			normalized()
	return(df)
}

filterSets <-function(df,hyp){
	if (hyp == 'hyp1') { df = df[padj1==1]}
	else if (hyp == 'hyp2') { df = df[padj2==1]}
}

main <- function(predictTable,model,enzyme,seqErr,pCutOff,resultFile,hyp,dbpath) {
	dataTable <- paste0(dbpath,
					'/',
					enzyme,
					'Table.tsv')
	dataTable <- dataTable %>%
		fread(colClasses=c('character','numeric','character',
						   rep('numeric',7),rep('character',2))) %>%
		transformTable(seqErr,pCutOff) %>%
		rename(label = abbrev) %>%
		mutate(label =  mergeType(label))  %>%
		filterSets('hyp1') 
	dataTable <- dataTable %>%
		group_by(label) %>%
		summarize(count = n()) %>%
		filter(count > 10) %>%
		inner_join(dataTable)

	predictTable <- predictTable %>%
		fread() %>%
		transformTable(seqErr,pCutOff) %>%
		filterSets(hyp) 
	cat('Read Data!\n')
	
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
cat('Using:',model,'for', enzyme,'\n')
main(predictTable,model,enzyme,seqErr,pCutOff,resultFile,hyp,dbpath)
cat ('Finished:',model,'for',enzyme,'\n')
