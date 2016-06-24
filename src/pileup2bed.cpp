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
#include <unordered_map>
#include "stringManipulation.h"
#include "pileupFix.h"
#include "stringVector.h"

typedef unordered_map<string , int> count_dict;

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

void printFunction(count_dict counts, string chrom, int mispos, int end,
		    string realRef, int cov, string strand, int insertion, int deletion)
{
    cout << chrom << '\t' << mispos << '\t' << end << '\t';
    cout << realRef << '\t' << cov << '\t' << strand << '\t';
    cout << counts["A"] << '\t' << counts["C"] << '\t' << counts["T"] << '\t' << counts["G"] << '\t';
    cout << insertion << '\t' << deletion << '\n';
}

//printing all variables
count_dict printTable(string chrom, string start, string ref, 
		int A, int C, int T, int G, 
		int a, int c, int t, int g,
		int insertion, int deletion, int refCount, int refCountrev)
{
	string strand ;
	int mispos = atoi(start.c_str());
	int end = mispos + 1;
	int cov;
	string realRef;
	int referenceCount;
	int print = 0;
	count_dict counts;
	transform(ref.begin(), ref.end(), ref.begin(), ::toupper);
	// clean up forward/reverse strand data
	if (A + C + T + G + refCount > 0)
	{
	    strand = "+";
	    realRef = ref;
	    referenceCount = refCount;
	    counts["A"] = A;
	    counts["C"] = C;
	    counts["G"] = G;
	    counts["T"] = T;
	    counts[realRef] = referenceCount;
	    cov =  counts["A"] + counts["C"] + counts["G"] + counts["T"];
	    printFunction(counts, chrom, mispos, end,
		    realRef, cov, strand, insertion, deletion);
	}

	if (a + c + t + g + refCountrev > 0)
	{
	    strand = "-";
	    realRef = reverse_complement(ref);
	    referenceCount = refCountrev;
	    counts["A"] = a;
	    counts["C"] = c;
	    counts["G"] = g;
	    counts["T"] = t;
	    counts[realRef] = referenceCount;
	    cov =  counts["A"] + counts["C"] + counts["G"] + counts["T"];
	    printFunction(counts, chrom, mispos, end,
		    realRef, cov, strand, insertion, deletion);
	}
	return counts;
}

// processing lines with mismatches 
void extractMismatches(string reads, string baseQuals, int cov, 
		string transcriptID, string mispos, 
		string ref, int qualThreshold, int coverageThreshold, string line)
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
	
	if (cov > coverageThreshold)
	{
	    count_dict countTable = printTable(transcriptID, mispos, ref,
		    A, C, T, G, a, c, t, g, 
		    insertion, deletion, refCount, refCountrev);
	    int infer_coverage = a + c + t + g + A + C + T + G + n + N + refCountrev + refCount;
	    bool condition( infer_coverage + deletion != cov);
	    if(condition)
	    {
		cerr << '\n';
		cerr << "#######   ASSERTION #########" << '\n';
		cerr << line << '\n';
		cerr << "A:" << A << "  C:" << C << "  G: " << G << " T: " << T << '\n';
		cerr << "a:" << a << "  c:" << c << "  g: " << g << " t: " << t << '\n';
		cerr << "N:" << N << "  n:" << n << '\n';
		cerr << "coverage:" << cov <<  "   infer coverage: " << infer_coverage <<  '\n';
		cerr << "reverse ref: " << refCountrev << "    ref: " << refCount << '\n';
		abort();
	    }
	}
		
}


// extract from each line different columns
// and give them to further processing
void processLine( stringList columns, int qualThreshold, int coverageThreshold, string line) 
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
					pos, ref, qualThreshold, coverageThreshold, line);
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
        processLine(columns, qualThreshold, coverageThreshold, line);
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
        processLine(columns, qualThreshold, coverageThreshold, line);
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
    cerr << "Using quality threshold:  "<< qualThreshold << '\n';
    cerr << "Using coverage threshold: "<< coverageThreshold << '\n';
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
