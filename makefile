cpp=g++
mkdir=mkdir -p
Rscript=Rscript

all: directory pileup2bed installRpackage

directory:
	$(mkdir) bin

pileup2bed:
	$(cpp) src/pileup2bed.cpp -o bin/pileup2bed

installRpackage:
	$(Rscript) src/install_packages.R
