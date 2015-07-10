#include <Rcpp.h>
#include <stdlib.h>
using namespace Rcpp;
using namespace std;

// define pairs
typedef pair<double,int> pairType;
typedef vector<string> stringList;


// define method for sorting
// advised by william chow
// https://williamchowhk.wordpress.com/2015/05/06/c-stdsort-objects/ 
struct comparison{
	inline bool operator() (const pairType& a, const pairType& b)
   	{
       	return a.first < b.first;
   	}
};

// function for sorting vector
NumericVector indexSort(NumericVector p)
{
    	// sorting  a vector in ascending order and return 
    	// a vector of the index of the sorted vector 
	int i, size = p.size();
	NumericVector resultVec(size);
	vector<pairType> newVec(size);

	//sorting
	for (i = 0; i < size; i++)
	{
		newVec[i] = make_pair(p[i],i);
	}
	sort(newVec.begin(),newVec.end(),comparison());

	//pushing back index vector
	for (i = 0; i < size; i ++)
	{
		resultVec[i] = newVec[i].second;
	}
	
	return(resultVec);
}


// [[Rcpp::export]]
// using a vector of p-values and false discovery rate threshold to 
// label all the significant results
// where 1 is significant,  is insignificant
NumericVector FDRcontrol(NumericVector p, double alpha)
{
	float threshold1, threshold2;
    	int i = 0, size = p.size();
    
	NumericVector passed(size);
	NumericVector index = indexSort(p);

    	while ( i < size)
    	{
        	threshold1 = alpha * (i+1) / size;
	        threshold2 = alpha * (i+2) / size;
        	if (threshold1 > p[index[i]])
        	{
        		if (threshold2 < p[index[i+1]])
            		{
                		for (int j = 0; j <= i; j++)
                		{
                    		passed[index[j]] = 1;
                		}
            		}
        	} 
        	i++;
    	}
    	return(passed);
}



//[[Rcpp::export]]
stringList string_split(stringList x, char sep, int num)
{
	int i, j, numChar, start, length, numStrings = x.size();
    	stringList result(numStrings);

    	for (i = 0 ; i < numStrings; i++ )
    	{
        	numChar = x[i].length();
		vector< int > sepPos(1); // which character is separater
		j = 0;
        	for (j = 0 ; j < numChar; j++)
        	{
			if (x[i][j] == sep)
            		{
		        	sepPos.push_back(j); 
            		}
        	}
		sepPos.push_back(numChar); // add last character number
		
		if (num==1)
		{	
			start = sepPos[num-1] ;
			length = sepPos[num]-sepPos[num-1];
		}
		else
		{
			start = sepPos[num-1] + 1;
			length = sepPos[num]-sepPos[num-1] -1;
		}
		result[i] = x[i].substr(start,length);
    	}
    	return result;
}

//[[Rcpp::export]]
stringList mergeType(stringList realbase)
{
    int i, size = realbase.size();
    stringList out(size);

    for (i = 0 ; i < size; i++)
    {
		string base = realbase[i];
        if (base == "m2G" || base == "m2,2G")
        {
            out[i] = "m2G|m2,2G";
        }
        else if (realbase[i] == "m6A" || base == "m6,6A")
		{
			out[i] = "m6A|m6,6A";
		}
		else 
        {
            out[i] = realbase[i];
        }
    }
    return out;
}

//[[Rcpp::export]]
vector<double> heterozygote(NumericVector A, NumericVector C, 
							NumericVector T, NumericVector G,
							NumericVector cov)
{
	int i, correctCount, size = A.size();
	vector<double> result(size);
	vector <double> baseCount(5);
	for (i = 0; i < size ; i ++)
	{
		baseCount[0] = A[i];
		baseCount[1] = C[i];
		baseCount[2] = T[i];
		baseCount[3] = G[i];
		baseCount[4] = 1 - accumulate(baseCount.begin(),baseCount.end() - 1,0.0);
		sort(baseCount.begin(),baseCount.end(),greater<double>());
		result[i] =  (1 - ( baseCount[0] + baseCount[1])) * cov[i];
	}
	return result;
}

//[[Rcpp::export]]
DataFrame normalized(DataFrame df)
{
	NumericVector A = df["A"] ;
	NumericVector C = df["C"] ;
	NumericVector T = df["T"] ;
	NumericVector G = df["G"] ;
	NumericVector cov = df["cov"] ;
	NumericVector mismatch = df["mismatch"];
	int size = A.size();
	NumericVector newA(size), newC(size), newT(size), newG(size);
	for (int i = 0 ; i < size ; i ++)
	{
		newA[i] = A[i] * cov[i] / mismatch[i];
		newC[i] = C[i] * cov[i] / mismatch[i];
		newT[i] = T[i] * cov[i] / mismatch[i];
		newG[i] = G[i] * cov[i] / mismatch[i];
	}
	df["A"] = newA;
	df["T"] = newT;
	df["C"] = newC;
	df["G"] = newG;
	return df;
}
