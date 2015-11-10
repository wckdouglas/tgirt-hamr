#!/bin/env python

from __future__ import division
import os
import sys
import time
import getopt

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

def usage(programname):
    sys.stderr.write('usage: python %s [options] -e <TGIRT> -i <input bam> -o <output bed> -r <reference fasta>\n' %programname)
    sys.stderr.write('[options]\n')
    sys.stderr.write('-e, --enzyme <gsi>|<tei>            This can be <gsi> or <tei> for GsI-IIc and TeI4c libraries respectively\n')
    sys.stderr.write('-i, --inBam <bamfile>               input bam file\n')
    sys.stderr.write('-o, --outBed <bedfile>              output bed file\n')
    sys.stderr.write('-r, --refFasta <fasta file>         reference fasta file\n')
    sys.stderr.write('-p, --cores <int>                   number of cores to use                                         default: 1 \n')
    sys.stderr.write('-y, --hyp <hyp1>|<hyp2>             hypothesis to use, can be <hyp1> or <hyp2>                     default: hyp2\n')
    sys.stderr.write('-s, --seqErr <float>                sequencing error probability                                   default: 0.01\n')
    sys.stderr.write('-t, --pThreshold <float>            False discovery rate cut off                                   default: 0.01\n')
    sys.stderr.write('-m, --model <knn>                   model for prediction, can be anything available in R::caret    default: knn\n')
    sys.stderr.write('-q, --qual <int>                    base quality cut off for retaining bases from read             default: 33\n')
    sys.stderr.write('-c, --cov <int>                     coverage cutoff for position to retain                         default: 10\n')
    sys.stderr.write('-d, --dev                           retain intermediate files\n')
    sys.stderr.write('-h, --help                          print out usage\n')
    sys.exit()

def main():
    programpath = '/'.join(os.path.abspath(__file__).split('/')[:-2])
    programname = sys.argv[0]
    try:
        opts, args = getopt.getopt(sys.argv[1:], \
                "e:i:o:r:p:y:s:t:m:q:c:hv", \
                ['enzyme=','inBam=', 'outBed=','refFasta=','cores=',\
                'hyp=','seqErr=','pThreshold=','model=','qual=','cov=','help','dev'])
    except getopt.GetoptError as err:
        usage(programname)
    cores = 1
    hyp = 'hyp2'
    seqErr = 0.01
    pThreshold = 0.01
    model = 'knn'
    bamFile = ''
    outBedFile = ''
    refFasta = ''
    enzyme = ''
    cov = 10
    qual = 33
    devMode = 0

    # ==================           read args =======================================
    for option, arg in opts:
        if option in ('-h','--help'):
            usage(programname)
        elif option in ('-i','--inBam'):
	        bamFile = arg
        elif option in ('-o','--outBed'):
	        outBedFile = arg
        elif option in ('-r','--refFasta'):
            refFasta = arg
        elif option in ('-e','--enzyme'):
            enzyme = arg
        elif option in ('-p','--cores'):
            cores = arg
        elif option in ('-y','--hyp'):
            hyp = arg
        elif option in ('-s','--seqErr'):
            seqErr = float(arg)
        elif option in ('-t','--pThreshold'):
            pThreshold = float(arg)
        elif option in ('-m','--model'):
            model = arg
        elif option in ('-q','--qual'):
            qual = int(arg)
        elif option in ('-c','--cov'):
            cov = int(arg)
        elif option in ('-v','--dev'):
            devMode = 1
        else:
            assert False, "unused option %s" %option
    message = '##################################\n' +\
              '#Using seqErr:        %.3f\n' %(seqErr) +\
              '#Using cores:         %s\n' %(cores) +\
              '#Using FDR threshold: %.3f\n' %(pThreshold) + \
              '#################################\n'
    sys.stderr.write(message)
	
    if bamFile == '' or outBedFile == '' or refFasta == '' or enzyme == '':
        sys.stderr.write('options: -i, -r, -e, -o are required!!!\n')
        usage(programname)

    # ========================   start pipeline =====================================
    start = time.time()
    sample = outBedFile.split('/')[-1].rsplit('.',1)[0]
    resultpath = '/'.join(outBedFile.split('/')[:-1])
    if resultpath == '':
        resultpath = './'

    tempFile = resultpath + '/' + sample + '.tsv'
    #extract pileup
    #pileup(bamFile, refFasta, tempFile, programpath, qual, cov)
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
