# TGIRT-HAMR

TGIRT-HAMR stands for thermostable group II intron reverse transcriptase (TGIRT) high-throughput annotation of modification using RNA-seq. This software adopted the algorithm in [HAMR](http://www.ncbi.nlm.nih.gov/pubmed/24149843) for predicting post-transcriptional modifications in TGIRT. 

INSTALL:

	git clone https://github.com/wckdouglas/tgirt-hamr.git
	cd tgirt-hamr
	make

Run:
	
	python $(tgirt-hamr-path)/script/tgirt-hamr.py	
	usage: python tgirt-hamr.py [options] -e <TGIRT> -i <input bam> -o <output bed> -r <reference fasta>
	[options]
	-e, --enzyme=<gsi>|<tei>            This can be <gsi> or <tei> for GsI-IIc and TeI4c libraries respectively
	-i, --inBam=<bamfile>               input bam file
	-o, --outBed=<bedfile>              output bed file
	-r, --refFasta=<fasta file>         reference fasta file
	-p, --cores=<int>                   number of cores to use                                         default: 1
	-y, --hyp=<hyp1>|<hyp2>             hypothesis to use, can be <hyp1> or <hyp2>                     default: hyp2
	-s, --seqErr=<float>                sequencing error probability                                   default: 0.01
	-t, --pThreshold=<float>            False discovery rate cut off                                   default: 0.01
	-m, --model=<knn>                   model for prediction, can be anything available in R::caret    default: knn
	-q, --qual=<int>                    base quality cut off for retaining bases from read             default: 33
	-c, --cov=<int>                     coverage cutoff for position to retain                         default: 10
	-h, --help                          print out usage

Dependencies:    
softwares: R >=3.1.0, python = 2.7, samtools >= 1.2, gcc > 4.4.5     
R-package: dplyr, caret, readr, stringr, tidyr, Rcpp, getopt, doMC, tgirthamr
