cpp=g++
mkdir=mkdir -p

all: directory pileup2bed

directory:
	$(mkdir) bin

pileup2bed:
	$(cpp) src/pileup2bed.cpp -o bin/pileup2bed
