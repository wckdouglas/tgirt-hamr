#include <string>
#include <vector>
#include <sstream>
#include "stringVector.h"

using namespace std;

//split function to split line with desired deliminator
stringList split(const string &s, char delim) 
{
	stringList result; 
    stringstream ss(s);
    string item;
    while (getline(ss, item, delim)) 
    {
        result.push_back(item);
    }
    return result;
}


string reverse_complement(string base)
{
	string basepair;
	if (base.compare("A") == 0) 
	{
		basepair = "T";
	}
	else if (base.compare("T") == 0) 
	{
		basepair = "A";
	}
	else if (base.compare("C") == 0) 
	{
		basepair = "G";
	}
	else if (base.compare("G") == 0) 
	{
		basepair = "C";
	}
	return basepair;
}

