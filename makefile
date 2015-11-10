cpp=g++
mkdir=mkdir -p
Rscript=`which Rscript `
installLink=https://gist.githubusercontent.com/wckdouglas/19db0cb30dfb387b1669/raw/4a9aee4b9e7f245f9088b2406be6c168b785e797/install_packages.R

all: directory pileup2bed installRpackage

directory:
	$(mkdir) bin

pileup2bed:
	$(cpp) src/pileup2bed.cpp -o bin/pileup2bed

installRpackage:
	curl $(installLink) | $(Rscript) -
