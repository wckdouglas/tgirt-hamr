CPP=g++
DIR=mkdir -p
RSCRIPT=`which Rscript`
CPP_FLAG=-I./include -std=c++11
LINK=https://gist.githubusercontent.com/wckdouglas/19db0cb30dfb387b1669/raw/4a9aee4b9e7f245f9088b2406be6c168b785e797/install_packages.R

all: directory pileup2bed installRpackage

directory:
	$(DIR) bin

pileup2bed:
	$(CPP) src/pileup2bed.cpp -o bin/pileup2bed $(CPP_FLAG)

installRpackage:
	curl $(LINK) | $(RSCRIPT) -
