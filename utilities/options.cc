// gemc headers
#include "options.h"
#include "string_utilities.h"

// C++ headers
#include <cstdio>
#include <set>
#include <cstdlib>

// initialize general options
goptions::goptions()
{
	optMap["LOG_VERBOSITY"].arg  = 0;
	optMap["LOG_VERBOSITY"].help = "Controls General Log Verbosity.";
	optMap["LOG_VERBOSITY"].name = "Log Verbosity";
	optMap["LOG_VERBOSITY"].type = 0;
	optMap["LOG_VERBOSITY"].ctgr = "verbosity";
}


void goptions::scanGcard(string file)
{
	// If found, parse the <options> section of the file.
	QDomDocument domDocument;
    
	cout << " >> Parsing " << file << " for options: \n";
    
	QFile gcard(file.c_str());
    
	if( !gcard.exists() )
	{
		cout << " >>  gcard: " << file <<" not found. Exiting." << endl;
		exit(0);
	}
	
	// opening gcard and filling domDocument
	if(!domDocument.setContent(&gcard))
	{
		gcard.close();
		cout << " >>  gcard format for file <" << file << "> is wrong - check XML syntax. Exiting." << endl;
		exit(0);
	}
	gcard.close();
	
	// setting this gcard file as the one in the option
	// this way the only gcard file with detectors in it
	// is the last one in the command line
	// however there should be only one anyway
	optMap["gcard"].args = file;
	
	/// reading gcard file
	QDomElement docElem = domDocument.documentElement();
	QDomNode n;
	
	map<string, int> count;
	for(map<string, aopt>::iterator itm = optMap.begin(); itm != optMap.end(); itm++)
		count[itm->first] = 0;
	
	///< looping over options
	n = docElem.firstChild();
	while(!n.isNull())
	{
		QDomElement e = n.toElement();                ///< converts the node to an element.
		if(!e.isNull())                               ///< the node really is an element.
			if(e.tagName().toStdString() == "option")     ///< selecting "option" nodes
			{
				int found = 0;
				for(map<string, aopt>::iterator itm = optMap.begin();itm != optMap.end(); itm++)
				{
					// looking for a valid option. If a two instances of the same option exist
					// the string __REPETITION__ will be appended
					if(e.attributeNode("name").value().toStdString() == itm->first )
					{

						found = 1;
						count[itm->first] += 1;

						// first time it finds it
						if(count[itm->first] == 1)
						{
							itm->second.args =      e.attributeNode("value").value().toStdString();
							itm->second.arg  = stringToDouble(e.attributeNode("value").value().toStdString());

							itm->second.printSetting();
						}

						else
						{
							string new_opt = itm->first + "__REPETITION__" + stringify(count[itm->first] - 1);
							optMap[new_opt].args = e.attributeNode("value").value().toStdString();
							optMap[new_opt].arg  = stringToDouble(e.attributeNode("value").value().toStdString());
							optMap[new_opt].name = itm->second.name;
							optMap[new_opt].help = itm->second.help;
							optMap[new_opt].type = itm->second.type;
							optMap[new_opt].ctgr = itm->second.ctgr;
							optMap[new_opt].printSetting();
						}
						break;
					}
				}
				if( found == 0 )
				{
					cout << "  !! Error: The option in the gcard file "
					<< e.attributeNode("name").value().toStdString()
					<< " is not known to this system. Please check your spelling." << endl;
					exit(3);
				}
			}
		
		// now looking for child arguments
		QDomNode nn= e.firstChild();
		while( !nn.isNull() && e.tagName().toStdString() == "option")
		{
			QDomElement ee = nn.toElement();
			int found=0;
			for(map<string, aopt>::iterator itm = optMap.begin(); itm != optMap.end(); itm++)
			{
				if(ee.tagName().toStdString() == itm->first )
				{
					itm->second.args= ee.attributeNode("value").value().toStdString();
					itm->second.arg = stringToDouble(ee.attributeNode("value").value().toStdString());
					itm->second.printSetting();
					found = 1;
				}
			}
			if( found == 0 )
			{
				cout << "  !! Error: The option in the gcard file "
				<< e.attributeNode("name").value().toStdString()
				<< " is not known to this system. Please check your spelling." << endl;
				exit(3);
			}
			nn = nn.nextSibling();
		}
		n = n.nextSibling();
	}
}



int goptions::setOptMap(int argc, char **argv)
{
	// Check the command line for the -gcard=file option.
	// Then call the parser for the gcard file that reads all the options.
	// This must be done BEFORE the commandline is parsed, so that command line
	// options override the options in the file.
	
	
	// Look for -gcard special option:
	size_t pos;
	for(int i=1;i<argc;i++)
	{
		string arg = argv[i];
		pos = arg.find("gcard=");
		if(pos != string::npos)
		{
			scanGcard(arg.substr(pos+6));
		}
	}
	
	// if no gcard option is passed, checking that one of the option is a gcard file
	for(int i=1;i<argc;i++)
	{
		string arg = argv[i];
		pos = arg.find(".gcard");
		if(pos != string::npos)
		{
			ifstream my_file(argv[i]);
			if(my_file)
			{
				scanGcard(argv[i]);
			}
		}
	}
	
    
	set<string> category;
	// Filling Categories
	for(map<string, aopt>::iterator itm = optMap.begin(); itm != optMap.end(); itm++)
		if(category.find(itm->second.ctgr) == category.end()) category.insert(itm->second.ctgr);
	
	
	// -help-all
	cout << endl;
	for(int i=1; i<argc; i++)
	{
		string arg = argv[i];
		string com = "-help-all";
		if(arg == com)
		{
			cout <<  "    Usage: -Option=<option>" << endl << endl;
			cout <<  "    Options:" <<  endl << endl ;
			
			for(map<string, aopt>::iterator itm = optMap.begin(); itm != optMap.end(); itm++)
				cout <<  "   > Option " <<  itm->first << ": " << itm->second.help << endl;
			
			cout << endl << endl;
			exit(0);
		}
	}
	
	
	// -help
	for(int i=1; i<argc; i++)
	{
		string arg = argv[i];
		string com1 = "-help";
		string com2 = "-h";
		if(arg == com1 || arg == com2)
		{
			cout <<  endl << endl;
			cout <<  "    Help Options:" <<  endl << endl ;
			cout <<  "   >  -help-all:  all available options. " <<  endl << endl;
			for(set<string>::iterator itcat = category.begin(); itcat != category.end(); itcat++)
			{
				cout <<  "   >  -help-" << *itcat << "     ";
				cout.width(15);
				cout << *itcat << " options." << endl;
			}
			cout << endl << endl;
			exit(0);
		}
	}
	
	
	// -help-option
	for(int i=1; i<argc; i++)
	{
		string arg = argv[i];
		for(set<string>::iterator itcat = category.begin(); itcat != category.end(); itcat++)
		{
			string com = "-help-" + *itcat;
			if(arg == com)
			{
				cout << endl << endl <<  "   ## " << *itcat << " ## " << endl << endl;
				cout << "    Usage: -Option=<option>" << endl << endl;
				for(map<string, aopt>::iterator itm = optMap.begin(); itm != optMap.end(); itm++)
					if(itm->second.ctgr == *itcat)
						cout <<  "   > " <<  itm->first << ": " << itm->second.help << endl;
				cout << endl << endl;
				exit(0);
			}
		}
	}
	
	
	// -help-html
	for(int i=1; i<argc; i++)
	{
		string arg = argv[i];
		for(set<string>::iterator itcat = category.begin(); itcat != category.end(); itcat++)
		{
			string com = "-help-html";
			if(arg == com)
			{
				ofstream hf;
				hf.open("options.html");
				hf << "<html>" << endl;
				hf << "	<STYLE TYPE=\"text/css\">" << endl;
				hf << "<!--" << endl;
				hf << ".pretty-table" << endl;
				hf << "{" << endl;
				hf << "	padding: 0;" << endl;
				hf << "	margin: 0;" << endl;
				hf << "	border-collapse: collapse;" << endl;
				hf << "	border: 1px solid #333;" << endl;
				hf << "	font-family: \"Trebuchet MS\", Verdana, Arial, Helvetica, sans-serif;" << endl;
				hf << "	font-size: 0.8em;" << endl;
				hf << "	color: #000;" << endl;
				hf << "	background: #bcd0e4;" << endl;
				hf << "}" << endl;
				hf << ".pretty-table caption" << endl;
				hf << "{" << endl;
				hf << "	caption-side: bottom;" << endl;
				hf << "	font-size: 0.9em;" << endl;
				hf << "	font-style: italic;" << endl;
				hf << "	text-align: right;" << endl;
				hf << "	padding: 0.5em 0;" << endl;
				hf << "}" << endl;
				hf << ".pretty-table th, .pretty-table td" << endl;
				hf << "{" << endl;
				hf << "	border: 1px dotted #666;" << endl;
				hf << "	padding: 0.5em;" << endl;
				hf << "	text-align: left;" << endl;
				hf << "	color: #632a39;" << endl;
				hf << "}" << endl;
				hf << ".pretty-table th[scope=col]" << endl;
				hf << "{" << endl;
				hf << "	color: #000;" << endl;
				hf << "	background-color: #8fadcc;" << endl;
				hf << "	text-transform: uppercase;" << endl;
				hf << "	font-size: 0.9em;" << endl;
				hf << "	border-bottom: 2px solid #333;" << endl;
				hf << "	border-right: 2px solid #333;" << endl;
				hf << "}" << endl;
				hf << ".pretty-table th+th[scope=col]" << endl;
				hf << "{" << endl;
				hf << "	color: #009;" << endl;
				hf << "	background-color: #7d98b3;" << endl;
				hf << "	border-right: 1px dotted #666;" << endl;
				hf << "}" << endl;
				hf << ".pretty-table th[scope=row]" << endl;
				hf << "{" << endl;
				hf << "	background-color: #b8cfe5;" << endl;
				hf << "	border-right: 2px solid #333;" << endl;
				hf << "}" << endl;
				hf << "pre{font-family:Helvetica;font-size:12pt}" << endl;
				
				hf << "--->" << endl;
				hf << "</STYLE>" << endl;
				hf << "</head>" << endl;
				hf << "<body>" << endl;
				hf << "<br><br>" << endl;
				hf << "<center>" << endl;
				hf << "<h1> GEMC options</h1>" << endl;
				hf << "</center>" << endl;
				hf << "<br><br><br>" << endl;
				hf << "<table cellsize=20>" << endl;
				hf << "<tr><td>" << endl;
				hf << "<table class=\"pretty-table\">" << endl;
				hf << "<caption>options. This table is produced with the option: -help-html </caption>" << endl;
				hf << "<tr><th scope=\"col\" >Category</th>" << endl;
				hf << "    <th scope=\"col\" >Option</th>" << endl;
				hf << "    <th scope=\"col\" >Help</th></tr>" << endl;
				for(set<string>::iterator itcat = category.begin(); itcat != category.end(); itcat++)
					for(map<string, aopt>::iterator itm = optMap.begin(); itm != optMap.end(); itm++)
						if(itm->second.ctgr == *itcat)
						{
							hf << "<tr><th scope=\"row\">";
							hf << *itcat ;
							
							hf << "</th> <td>";
							hf << itm->first;
							hf << "</td><td><pre>" << endl;
							hf << itm->second.help;
							hf << "</pre></td></tr>" << endl;
						}
				
				hf << "</table>" << endl;
				hf << "</td>" << endl;
				hf << "<td>" << endl;
				hf << "</table>" << endl;
				hf << " </body></html>";
				
				hf.close();
				exit(0);
			}
		}
	}
	
	map<string, int> count;
	for(map<string, aopt>::iterator itm = optMap.begin(); itm != optMap.end(); itm++)
		count[itm->first] = 0;
	
	for(int i=1; i<argc; i++)
	{
		string arg = argv[i];
		int found=0;
		for(map<string, aopt>::iterator itm = optMap.begin(); itm != optMap.end(); itm++)
		{
			string com = "-" + itm->first + "=";
			string comp;
			comp.assign(arg, 0, arg.find("=") + 1);
			
			// skip if argument is a file
			ifstream my_file(argv[i]);
			if(my_file)
				found = 1;

			
			
			if(comp == com)
			{
				found=1;
				count[itm->first] += 1;

				string opts;
				opts.assign(arg, com.size(), arg.size()-com.size());
				
				// first time 
				if(count[itm->first] == 1)
				{
					itm->second.args = opts;
					itm->second.arg  = stringToDouble(opts);
					itm->second.printSetting();
				}
				
				
				if(count[itm->first]>1)
				{
					string new_opt = itm->first + "__REPETITION__" + stringify(count[itm->first]-1);
					optMap[new_opt].args  = opts;
					optMap[new_opt].arg   = stringToDouble(opts);
					optMap[new_opt].name  = itm->second.name;
					optMap[new_opt].help  = itm->second.help;
					optMap[new_opt].type  = itm->second.type;
					optMap[new_opt].ctgr  = itm->second.ctgr;
					optMap[new_opt].printSetting();
				}
				break;
			}
		}
		
		// For MAC OS X, we want to ignore the -psn_# type argument. This argument is added by
		// the system when launching an application as an "app", and # contains the process id.
		if( found == 0 && strncmp(argv[i],"-psn_", 4) !=0 && ignoreNotFound == 0)
		{
			cout << " The argument " << argv[i] << " is not known to this system / or file not found. Continuing anyway.\n\n";
			// exit(2);
		}
	}
	
	cout << endl;
	
	return 1;
}

vector<aopt> goptions::getArgs(string opt)
{
	vector<aopt> options;
	map<string, int> count;
	for(map<string, aopt>::iterator itm = optMap.begin();itm != optMap.end(); itm++)
	{
		if(itm->first.find(opt) != string::npos)
		{
			options.push_back(itm->second);
		}
	}
	
	return options;
}


// Returns a map<string, string> with all options and values
map<string, string> goptions::getOptMap()
{
	map<string, string> optmap;
	
	for(map<string, aopt>::iterator it = optMap.begin(); it != optMap.end(); it++)
	{
		string key = "option " + it->first;
		if(it->second.type == 0) optmap[key] = stringify(it->second.arg);
		else                     optmap[key] = it->second.args;
	}
	
	return optmap;
}




// print the option setting
void aopt::printSetting()
{
	if(name.find("gemc card file") != string::npos) return;
	
	
	cout << "  > " << name << " set to: ";
	if(type) cout << args;
	else     cout << arg;
	cout << endl;
	
}









