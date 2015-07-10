#include <string>
#include "stringVector.h"

using namespace std;

int countDigits(int number) 
{
	if (number < 10) 
	{
		return 1;
	}
	int count = 0;
	while (number > 0) 
	{
		number /= 10;
		count++;
	}
	return count;
}

char reverseStrandcomplement(char base)
{
	char basepair;
	if (base == 'a') 
	{
		basepair = 't';
	}
	else if (base == 't') 
	{
		basepair = 'a';
	}
	else if (base == 'c') 
	{
		basepair = 'g';
	}
	else if (base == 'g') 
	{
		basepair = 'c';
	}
	else
	{
		basepair = base;
	}
	return basepair;
}

void baseCount(int &A, int &C, int &T,int &G, int &N,
				int &a, int &c, int &t, int &g,int &n,
				char readPos, int &refCount, int &refCountrev)
{
	if (readPos == 'A')
	{	
		A ++;
	}
	else if (readPos == 'C')
	{
		C ++;
	}
	else if(readPos == 'G')
	{
		G ++;
	}
	else if (readPos == 'T')
	{
		T ++;
	}
	else if (readPos == 'N')
	{
		N ++;
	}
	else if (readPos == 'a')
	{
		a ++;
	}
	else if (readPos == 'c')
	{
		c ++;
	}
	else if (readPos == 't')
	{
		t ++;
	}
	else if (readPos == 'g')
	{
		g ++;
	}
	else if (readPos == 'n')
	{
		n ++;
	}
	else if (readPos == '.')
	{
		refCount ++;
	}
	else if (readPos == ',')
	{
		refCountrev ++;
	}
}

void fixpileup(int &A,int &C, int &T, int &G, int &N, 
				int &a, int &c, int &t, int&g, int &n,
				int &deletion, int &insertion, 
				string reads, string baseQuals,
				int qualThreshold, int &cov, 
				int &refCount, int &refCountrev,
				int &start, int &end)
{
	int i = 0, j = 0, current = 0, qual;
	char readPos, fixedReadPos;
    while (i < reads.length())
    {
		readPos = reads.at(i);
        if (readPos == '+')
		//insertion
        {
			i ++ ; 
			current = 0;
			insertion ++;
			while (isdigit(reads.at(i)))
			{
				current += current * 10 + (reads[i]-'0');
				i++;
			}
			i += current - countDigits(current);
        }
		else if (readPos == '-')
		// deletion
		{
			i ++ ; 
			current = 0;
			deletion ++;
			while (isdigit(reads.at(i)))
			{
				current += current * 10 + (reads[i]-'0');
				i++;
			}
			i += current - countDigits(current);
		}
        else if (readPos == '^')
        {
            i ++;
			if (islower(reads.at(i+1)) || reads.at(i+1) == ',')
			{
				start ++;
			}
			else if (isupper(reads.at(i+1)) || reads.at(i+1) == '.')
			{
				end ++;
			}
        }
        else if (readPos == '$')
        {
			if (islower(reads.at(i-1)) || reads.at(i-1) == ',')
			{
				end ++;
			}
			else if (isupper(reads.at(i-1)) || reads.at(i-1) == '.')
			{
				start ++;
			}
        }
		else if (readPos == '<' || readPos == '>')
		{
			cov --;
		}
		else 
		{
            qual = baseQuals[j] - 33 ;
			j++;
			if (qual < qualThreshold || readPos == '*')
			{
				cov --;
			}
			else 
			{
				fixedReadPos = reverseStrandcomplement(readPos);
				baseCount(A,C,T,G,N,
						a,c,t,g,n,
						fixedReadPos,refCount,refCountrev);
			}
		}
		i++;
    }
}
