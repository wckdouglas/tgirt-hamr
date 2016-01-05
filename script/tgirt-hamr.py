#!/bin/env python

from __future__ import division
import os
import sys
import time
import argparse

def pileup(bamFile, refFasta, resultFile, programpath, qualThresh, covThresh):
    start = time.time()
    command = 'samtools mpileup -d 1000000 -f %s %s ' %(refFasta,bamFile)+\
        '| %s/bin/pileup2bed - %i %i > %s' %(programpath, qualThresh, covThresh, resultFile)
    print command
    os.system(command)
    usedtime= time.time() - start
    sys.stderr.write('Pileup used time: %.3f min\n' %(usedtime/60))
    return 0

def prediction(inFile,resultBed, cores, programpath, model, seqErr,pThreshold,enzyme, hyp, devMode):
    start = time.time()
    model = 'knn'
    if devMode == 1:
        dev = '--dev'
    else:
        dev = ' '
    command = 'Rscript %s/src/prediction.R -o %s -t %s ' %(programpath,resultBed,cores)+\
                '-i %s -e %s -h %s -s %.4f -p %.4f ' %(inFile,enzyme,hyp,seqErr,pThreshold)+\
                '-m %s -d %s/table -f %s/src '  %(model, programpath,programpath) +\
                dev
    print command
    os.system(command)
    usedtime= time.time() - start
    sys.stderr.write('Prediction used time: %.3f min\n' %(usedtime/60))
    return 0

def main():
    programpath = '/'.join(os.path.abspath(__file__).split('/')[:-2])
    programname = sys.argv[0]
    parser = argparse.ArgumentParser(description='Pipeline for running TGIRT-HAMR, predicting modifications in TGIRT-seq data.')
    parser.add_argument('-i','--inBam', help='input bam file', required=True)
    parser.add_argument('-o','--outBed', help='output bed file', required=True)
    parser.add_argument('-r','--refFasta', help='reference fasta file', required=True)
    parser.add_argument('-e','--enzyme', help='This can be <gsi> or <tei> for GsI-IIc and TeI4c libraries respectively', choices=['gsi','tei'], required=True)
    parser.add_argument('-p','--cores', help='number of cores to use (default: 1)', default = 1, type=int)
    parser.add_argument('-y','--hyp', help='hypothesis to use, can be <hyp1> or <hyp2> (default: hyp2)', default = 'hyp2', choices = ['hyp1','hyp2'])
    parser.add_argument('-s','--seqErr', help='sequencing error probability (default: 0.01)', default = 0.01, type=float)
    parser.add_argument('-t','--pThreshold', help='False discovery rate cut off (default: 0.01)', default = 0.01, type=float)
    parser.add_argument('-m','--model', help='model for prediction, can be anything available in R::caret (default: knn)', default = 'knn')
    parser.add_argument('-q','--qual', help='base quality cut off for retaining bases from read (default: 33)', default = 33 , type=int)
    parser.add_argument('-c','--cov', help='coverage cutoff for position to retain (default: 10)', default = 10 , type=int)
    parser.add_argument('-v','--dev', help='developer mode',  action='store_true')
    
    args = parser.parse_args()
    bamFile = args.inBam
    outBedFile = args.outBed
    refFasta = args.refFasta
    enzyme = args.enzyme
    cores = args.cores
    hyp = args.hyp
    seqErr = args.seqErr
    pThreshold = args.pThreshold
    model = args.model
    qual = args.qual
    cov = args.cov
    dev = args.dev

    devMode = 1 if dev else 0
    message = '##################################\n' +\
              '#Using seqErr:        %.3f\n' %(seqErr) +\
              '#Using cores:         %s\n' %(cores) +\
              '#Using FDR threshold: %.3f\n' %(pThreshold) + \
              '#Using reference: %s\n' %(refFasta) + \
              '#################################\n'
    sys.stderr.write(message)
	
    # ========================   start pipeline =====================================
    start = time.time()
    sample = outBedFile.split('/')[-1].rsplit('.',1)[0]
    resultpath = '/'.join(outBedFile.split('/')[:-1])
    if resultpath == '':
        resultpath = './'

    tempFile = resultpath + '/' + sample + '.tsv'
    #extract pileup
    pileup(bamFile, refFasta, tempFile, programpath, qual, cov)
    #prediction
    prediction(tempFile, outBedFile, cores, programpath, model, \
            seqErr,pThreshold,enzyme, hyp, devMode)

    #remove temp files
    if devMode== 0:
        try:
            os.remove(tempFile)
        except OSError:
            pass
	
    # print summary
    usedTime = time.time() - start
    sys.stderr.write('Finsihed pipeline in %.3f min\n' %(usedTime/60))
    return 0

if __name__ == '__main__':
	main()
