/// \file options.h
/// Defines the option class.\n
/// The main argv options
/// are filled into a map<string opts>.
/// \author \n Maurizio Ungaro
/// \author mail: ungaro@jlab.org\n\n\n

#ifndef OPTIONS_H
#define OPTIONS_H


// Qt4 headers
#include <QDomDocument>
#include <QString>
#include <QFile>

// C++ headers
#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <ctime>
#include <fstream>
using namespace std;

/// \class aopt
/// <b> aopt </b>\n\n
/// Option class.\n
/// - arg:  double assigned to argument.
/// - args: string assigned to argument.
/// - name: name to be displayed for the argument variable.
/// - help: help for the argument variable.
/// - type: 0 = number, 1 = string
class aopt
{
public:
	double  arg;  ///< double assigned to argument.
	string args;  ///< string assigned to argument.
	string name;  ///< name to be displayed for the argument variable.
	string help;  ///< help for the argument variable.
	int    type;  ///< 0 = number, 1 = string
	string ctgr;  ///< help category
	int    repe;  ///< if this is set to 1: then this option can be repeated
					  ///< if this is set to 0: command line will always overwrite the gcard
	
	aopt() {repe = 0;}
	
public:
	void printSetting();
};


/// \class options
/// <b> options </b>\n\n
/// This is the general options class.
/// It contains a map of opt where the key is
/// the option string given at command line or gcard\n
class goptions
{
public:
	
	goptions();
	~goptions(){;}
	
	virtual void setGoptions();                 ///< Function to fill optMap
	void setOptions();                          ///< Define option Map
	void scanGcard(string file);                ///< Scan option file for options
	int setOptMap(int argc, char **args);       ///< Sets map from command line arguments
	int setOptMap(int argc, char **args, int i) ///< Sets map from command line arguments - include ignoreNotFound
	{
		ignoreNotFound = i;
		return setOptMap(argc, args);
	}
	
	map<string, aopt> optMap;              ///< Options map
	map<string, string> getOptMap();       ///< Returns a map<string, string> with all options and values
	
	vector<aopt> getArgs(string);          ///< get a vector of arguments matching a string
	
	int ignoreNotFound;                    ///< if this is not 0, it means the options
														///  could be shared with another application. So ignore
														///  the case that options are not found in the current application
};


#endif





