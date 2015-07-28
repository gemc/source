/// \file string_utilities.h
/// Set of string manipulation functions: \n
/// \author \n Maurizio Ungaro
/// \author mail: ungaro@jlab.org\n\n\n

#ifndef string_utilities_H
#define string_utilities_H

// Qt4 headers
#include <QString>
#include <QVariant>

// C++ headers
#include <vector>
#include <iostream>
#include <string>
#include <sstream>
using namespace std;


inline string stringify(double x)
{
	ostringstream o;
	o << x;
	return o.str();
}

inline string stringify(int x)
{
	ostringstream o;
	o << x;
	return o.str();
}


// from QVariant to string
// OS dependance
inline string qv_tostring(QVariant input)
{
	
	string output;
	
	// In Windows we need to initialize the string from toAscii
	// In Posix we need the toStdString
	#ifdef _MSC_VER
		output = input.toString().toLatin1();
	#else
		output = input.toString().toStdString();
	#endif
	
	return output;
}

// from QString to string
// OS dependance
inline string qs_tostring(QString input)
{
	
	string output;
	
	// In Windows we need to initialize the string from toAscii
	// In Posix we need the toStdString
	#ifdef _MSC_VER
		output = input.toLatin1();
	#else
		output = input.toStdString();
	#endif
	
	return output;
}


//replaces a character to another characters.
inline string replaceCharWithChars(string input, string x, string y)
{
	
	string output = "";
	
	for(unsigned int k=0; k<input.size(); k++)
	{
		int replace = 1;
		for(unsigned int j=0; j<x.size(); j++)
		{
			
			if(input[k] == x[j])
			{
				output.append(y);
				replace = 0;
			}
		}
		if(replace) output += input[k];
	}
	
	return output;
	
}

//replaces a string that has more than one character to another string
inline string replaceCharsWithChars(string input, string x, string y)
{
	
	string output = "";
	unsigned long k = x.length();
	unsigned long i = 0;
	
	
	
	while(i<input.size())
	{
		unsigned int j=0;
		unsigned int l=0;
		while(j<k)
		{
			if(input[i+j] == x[j])
				l++;
			j++;
		}
		
		if(l==k)
		{
			output.append(y);
			i = i+k;
		}
		else
		{
			output = output+input[i];
			i++;
		}
	}
	
	
	return output;
}

vector< vector<string> > dimensionstype(string);    ///< Returns dimensions nomenclature for different solid type
double get_number(string,int warn_no_unit=0);       ///< Returns number with dimension from string, i.e. 100*cm
string TrimSpaces(string);                          ///< Removes leading and trailing spaces
vector<string> get_strings_except(string, string);  ///< returns a vector of strings from a stringstream, space is delimiter, ignore string with second argument
vector<string> get_strings(string);                 ///< returns a vector of strings from a stringstream, space is delimiter
vector<string> get_strings(string, string);         ///< returns a vector of strings from a stringstream, second string is delimiter
void print_vstring(vector<string>);                 ///< prints each element of a string vector
vector<string> get_info(string);                    ///< get information from strings such as "5*GeV, 2*deg, 10*deg"
string get_variation(string);                       ///< parse variation name from string
bool is_main_variation(string);                     ///< returns 1 if the string "main:" is found on the input

ostream &operator<<(ostream &stream, map<string, string>);  ///< overload << for map<string, string>


// returns a double from a QVariant
inline double get_number(QVariant input)
{
	return get_number(TrimSpaces(qv_tostring(input)));
}


// returns a double from a string
// notice, atof should work but it doesn't work on some Mac
// atof may also not be thread safe
// solving all this with stringstream
inline double stringToDouble(string v)
{
	stringstream ss(TrimSpaces(v));
	double d;
	ss >> d;
	return d;
}


// returns a double from a string
// notice, atof should work but it doesn't work on some Mac
// atof may also not be thread safe
// solving all this with stringstream
// from QString to string
// OS dependance
inline double qs_toDouble(QString input)
{	
	return stringToDouble(qs_tostring(input));
}

#endif





