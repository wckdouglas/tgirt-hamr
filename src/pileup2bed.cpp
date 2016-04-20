#include <iostream>
#include <vector>
#include <sstream>
#include <stdlib.h>
#include <string.h> 
#include <ctype.h>
#include <cstdlib>
#include <fstream>
#include <cassert>
#include <algorithm>
#include "stringManipulation.h"
#include "pileupFix.h"
#include "stringVector.h"

//print usage 
int usage(char *argv[])
{
	cerr << "usage: "<< argv[0] << " <filename>|<stdin> ";
	cerr << "<quality threshold> <coverage threshold>" << endl;
	cerr << endl;
	cerr << "Takes in mpileup result format:"<<endl;
	cerr << "samtools mpileup -f <ref.fa> <bamFile> | ";
	cerr << argv[0] << " - ";
	cerr << "<quality threshold> <coverage threshold>" << endl;
	cerr << endl;
	return 0;
}

//printing all variables
void printTable(string chrom, string start, string ref, 
				int A, int C, int T, int G, 
				int a, int c, int t, int g,
				int insertion, int deletion, int refCount, int refCountrev)
{
	char strand ;
	int mispos = atoi(start.c_str()) - 1;
	int end = mispos + 1;
	int cov;
	string realRef;
	int adenosine, thymine, cytidine, guanosine;
	int referenceCount;
	int print = 0;
	transform(ref.begin(), ref.end(), ref.begin(), ::toupper);
	// clean up forward/reverse strand data
	if (A + C + T + G > 0)
	{
		strand = '+';
		realRef = ref;
		adenosine = A;
		cytidine = C;
		thymine = T;
		guanosine = G;
		referenceCount = refCount;
		print = 1;
	}
	else if (a + c + t + g > 0)
	{
		strand = '-';
		realRef = reverse_complement(ref);
		cov = a + c + t + g + refCountrev;
		adenosine = a;
		cytidine = c;
		thymine = t;
		guanosine = g;
		referenceCount = refCountrev;
		print = 1;
	}
	//print all variables 
	if (print == 1)
	{
		cov = adenosine + cytidine + thymine + guanosine + referenceCount;
		cout << chrom << '\t' << mispos << '\t' << end << '\t';
		cout << realRef << '\t' << cov << '\t' << strand << '\t';
		cout << adenosine << '\t' << cytidine << '\t' << thymine << '\t' << guanosine << '\t';
		cout << insertion << '\t' << deletion << '\n';
	}
}

// processing lines with mismatches 
void extractMismatches(string reads, string baseQuals, int cov, 
		string transcriptID, string mispos, 
		string ref, int qualThreshold, int coverageThreshold)
{
	int start = 0, end = 0, i = 0;
	int A = 0, C = 0, T = 0, G = 0, N = 0; 
	int a = 0, c = 0, t = 0, g = 0, n = 0; 
	int qual;
	int insertion = 0, deletion = 0, current = 0;
	int refCount = 0, refCountrev = 0;
	fixpileup(A, C, T, G, N,
			a, c, t, g, n,
			deletion, insertion, reads, baseQuals,
			qualThreshold, cov, refCount, refCountrev, start, end);
	cov = cov +  deletion - n - N;
	
	if (cov > coverageThreshold && (refCount + refCountrev) !=  cov)
	{
		printTable(transcriptID, mispos, ref,
					A, C, T, G, a, c, t, g, 
					insertion, deletion, refCount, refCountrev);
		assert (A + T + G + C + 
				a + c + t + g + 
				refCount + refCountrev + deletion == cov);
	}
		
}


// extract from each line different columns
// and give them to further processing
void processLine( stringList columns, int qualThreshold, int coverageThreshold) 
{
	if (columns[2] != "N" && columns[2] != "." && columns[2] != "_")
	{
	    string transcriptID, pos, ref, reads, baseQuals;
	    int cov;
	    if (columns.size() == 6) 
	    {
	        cov = atoi(columns[3].c_str());
	        if (cov > coverageThreshold)
	        { 
	            transcriptID = columns[0];
	            pos = columns[1];
	            ref = columns[2];
	            reads = columns[4];
	            baseQuals = columns[5];
	            assert ( baseQuals.length() == cov ) ;
				extractMismatches(reads, baseQuals, cov, transcriptID, 
					pos, ref, qualThreshold, coverageThreshold);
            }
        }
    }
}


// if lines are read from file,
// this function takes in and open the file and 
// parse it line by line
void readFile(const char* filename, int qualThreshold, int coverageThreshold)
{
    ifstream myfile(filename);
    for (string line; getline(myfile, line);)
    {
        stringList columns = split(line,'\t');
        processLine(columns, qualThreshold, coverageThreshold);
    }
}

// if lines are read from stdin,
// this function takes in and open the file and 
// parse it line by line
void readStream(int qualThreshold, int coverageThreshold)
{
    for (string line; getline(cin, line);)
    {
        stringList columns = split(line,'\t');
        processLine(columns, qualThreshold, coverageThreshold);
    }
}

void printHeader()
{
    cout << "chrom\tstart\tend\tref\tcov" ;
	cout << "\tstrand\tA\tC\tT\tG\t";
	cout << "insertion\tdeletion\n";
}

// main function
int main(int argc, char *argv[])
{
    ios::sync_with_stdio(false);
    // warnings
    if (argc != 4)
    {
        usage(argv);
	return 0;
    }

    // create modified RNA index
    int qualThreshold, coverageThreshold;
    qualThreshold = atoi(argv[2]);
    coverageThreshold = atoi(argv[3]);
    // read lines
    printHeader();
    if (strcmp(argv[1],"-") == 0)
    {
	cerr << "Reading from stdin" << endl;
        readStream(qualThreshold, coverageThreshold);
    }
    else
    {
        const char* filename = argv[1];
	cerr << "Reading from: " << filename << endl;
        readFile(filename, qualThreshold, coverageThreshold);
    }
    return 0;
}
