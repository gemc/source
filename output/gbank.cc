/// \file gbank.cc
/// Contains:
/// - read_banks: reads the Hit Process Map and
/// builds the gemc bank class map.\n\n
/// \author \n Maurizio Ungaro
/// \author mail: ungaro@jlab.org\n\n\n

// Qt headers
#include <QtSql>

// gemc headers
#include "gbank.h"
#include "string_utilities.h"
#include "utils.h"

// Variable Type is two chars.
// The first char:
//  R for raw integrated variables
//  D for dgt integrated variables
//  S for raw step by step variables
//  M for digitized multi-hit variables
//  V for voltage(time) variables

// The second char:
// i for integers
// d for doubles
// s for strings


map<string, gBank> read_banks(goptions gemcOpt, map<string, string> allSystems)
{
	double verbosity     = gemcOpt.optMap["BANK_VERBOSITY"].arg ;
	
	// geant4 information, integrated over
	map<string, gBank> banks;
	
	gBank abank;
	
	// event header
	abank = gBank(HEADER_BANK_TAG, "header", "Data Header Bank");
	abank.load_variable("time",       1, "Ns", "Time-stamp");
	abank.load_variable("runNo",      2, "Ni", "Run Number");
	abank.load_variable("evn",        3, "Ni", "Event Number");
	abank.load_variable("evn_type",   4, "Ni", "Event Type. 1 for physics events, 10 for scaler. Negative sign for MC.");
	abank.load_variable("beamPol",    5, "Nd", "Beam Polarization");
	abank.load_variable("var1",       6, "Nd", "User defined. In LUND this was: Target Polarization");
	abank.load_variable("var2",       7, "Nd", "User defined. In LUND this was: Number of nucleons in the target");
	abank.load_variable("var3",       8, "Nd", "User defined. In LUND this was: Number of protons in the target");
	abank.load_variable("var4",       9, "Nd", "User defined. In LUND this was: Bjorken x");
	abank.load_variable("var5",      10, "Nd", "User defined. In LUND this was: Fraction of energy loss");
	abank.load_variable("var6",      11, "Nd", "User defined. In LUND this was: W square");
	abank.load_variable("var7",      12, "Nd", "User defined. In LUND this was: Q square");
	abank.load_variable("var8",      13, "Nd", "User defined. In LUND this was: Energy loss");
	abank.orderNames();
	banks["header"] = abank;
	
	// generated particle infos
	abank =  gBank(GENERATED_PARTICLES_BANK_TAG, "generated", "Generated Particles");
	abank.load_variable("pid",    1,  "Ni", "Particle ID");
	abank.load_variable("px",     2,  "Nd", "x component of momentum");
	abank.load_variable("py",     3,  "Nd", "y component of momentum");
	abank.load_variable("pz",     4,  "Nd", "z component of momentum");
	abank.load_variable("vx",     5,  "Nd", "x component of vertex");
	abank.load_variable("vy",     6,  "Nd", "y component of vertex");
	abank.load_variable("vz",     7,  "Nd", "z component of vertex");
	abank.orderNames();
	banks["generated"] = abank;
	
	// particle summary infos
	// this is a daughter bank of the generated particle infos
	abank =  gBank(GENERATED_SUMMARY_BANK_TAG, "psummary", "Generated Particles Summary");
	abank.load_variable("dname",    1,  "Ns", "Detector Name");
	abank.load_variable("stat",     2,  "Ni", "Number of Hits in the detector");
	abank.load_variable("etot",     3,  "Nd", "Total Energy Deposited");
	abank.load_variable("t",        4,  "Nd", "Fastest Time on Detector");
	abank.load_variable("nphe",     5,  "Ni", "Number of Photoelectrons");
	abank.orderNames();
	banks["psummary"] = abank;
	
	
	// all hit bank information have a variable with num 0 "hitn" = hit number;
	
	// geant4 raw integrated
	// common for all banks
	abank =  gBank(RAWINT_ID, "raws", "Geant4 true information integrated over the hit");
	abank.load_variable("pid",     1,   "Ri", "ID of the first particle entering the sensitive volume");
	abank.load_variable("mpid",    2,   "Ri", "ID of the mother of the first particle entering the sensitive volume");
	abank.load_variable("tid",     3,   "Ri", "Track ID of the first particle entering the sensitive volume");
	abank.load_variable("mtid",    4,   "Ri", "Track ID of the mother of the first particle entering the sensitive volume");
	abank.load_variable("otid",    5,   "Ri", "Track ID of the original track that generated the first particle entering the sensitive volume");
	abank.load_variable("trackE",  6,   "Rd", "Energy of the track");
	abank.load_variable("totEdep", 7,   "Rd", "Total Energy Deposited (in MeV");
	abank.load_variable("avg_x",   8,   "Rd", "Average X position in the global reference system (in mm)");
	abank.load_variable("avg_y",   9,   "Rd", "Average Y position in the global reference system (in mm)");
	abank.load_variable("avg_z",  10,   "Rd", "Average Z position in the global reference system (in mm)");
	abank.load_variable("avg_lx", 11,   "Rd", "Average X position in the local reference system (in mm)");
	abank.load_variable("avg_ly", 12,   "Rd", "Average Y position in the local reference system (in mm)");
	abank.load_variable("avg_lz", 13,   "Rd", "Average Z position in the local reference system (in mm)");
	abank.load_variable("px",     14,   "Rd", "x component of momentum of the particle entering the sensitive volume");
	abank.load_variable("py",     15,   "Rd", "y component of momentum of the particle entering the sensitive volume");
	abank.load_variable("pz",     16,   "Rd", "z component of momentum of the particle entering the sensitive volume");
	abank.load_variable("vx",     17,   "Rd", "x component of point of origin of the particle entering the sensitive volume");
	abank.load_variable("vy",     18,   "Rd", "y component of point of origin of the particle entering the sensitive volume");
	abank.load_variable("vz",     19,   "Rd", "z component of point of origin of the particle entering the sensitive volume");
	abank.load_variable("mvx",    20,   "Rd", "x component of point of origin the mother of the particle entering the sensitive volume");
	abank.load_variable("mvy",    21,   "Rd", "y component of point of origin of the mother of the particle entering the sensitive volume");
	abank.load_variable("mvz",    22,   "Rd", "z component of point of origin of the mother of the particle entering the sensitive volume");
	abank.load_variable("avg_t",  23,   "Rd", "Average time");
	abank.load_variable("nsteps", 24,   "Ri", "Number of geant4 steps");
	abank.load_variable("procID", 25,   "Ri", "Process that created the particle. It's an integer described at gemc.jlab.org");
	abank.load_variable("hitn",   99,   "Ri", "Hit Number");
	abank.orderNames();
	banks["raws"] = abank;
	
	
	// geant4 raw step by step
	// common for all banks
	abank =  gBank(RAWSTEP_ID, "allraws", "Geant4 true information step by step");
	abank.load_variable("pid",     1,   "Si", "ID of the particle in the sensitive volume");
	abank.load_variable("mpid",    2,   "Si", "ID of the mother in the sensitive volume");
	abank.load_variable("tid",     3,   "Si", "Track ID of the particle in the sensitive volume");
	abank.load_variable("mtid",    4,   "Si", "Track ID of the mother of the particle in the sensitive volume");
	abank.load_variable("otid",    5,   "Si", "Track ID of the original track that generated the particle in the sensitive volume");
	abank.load_variable("trackE",  6,   "Sd", "Energy of the track");
	abank.load_variable("edep",    7,   "Sd", "Energy Deposited");
	abank.load_variable("x",       8,   "Sd", "X position in global reference system");
	abank.load_variable("y",       9,   "Sd", "Y position in global reference system");
	abank.load_variable("z",      10,   "Sd", "Z position in global reference system");
	abank.load_variable("lx",     11,   "Sd", "X position in local reference system");
	abank.load_variable("ly",     12,   "Sd", "Y position in local reference system");
	abank.load_variable("lz",     13,   "Sd", "Z position in local reference system");
	abank.load_variable("px",     14,   "Sd", "x component of momentum of the particle in the sensitive volume");
	abank.load_variable("py",     15,   "Sd", "y component of momentum of the particle in the sensitive volume");
	abank.load_variable("pz",     16,   "Sd", "z component of momentum of the particle in the sensitive volume");
	abank.load_variable("vx",     17,   "Sd", "x component of primary vertex of the particle in the sensitive volume");
	abank.load_variable("vy",     18,   "Sd", "y component of primary vertex of the particle in the sensitive volume");
	abank.load_variable("vz",     19,   "Sd", "z component of primary vertex of the particle in the sensitive volume");
	abank.load_variable("mvx",    20,   "Sd", "x component of primary vertex of the mother of the particle in the sensitive volume");
	abank.load_variable("mvy",    21,   "Sd", "y component of primary vertex of the mother of the particle in the sensitive volume");
	abank.load_variable("mvz",    22,   "Sd", "z component of primary vertex of the mother of the particle in the sensitive volume");
	abank.load_variable("t",      23,   "Sd", "time");
	abank.load_variable("stepn",  98,   "Si", "step index");
	abank.load_variable("hitn",   99,   "Si", "Hit Number");
	abank.orderNames();
	banks["allraws"] = abank;
	
	
	
	// flux bank integrated digitized infos
	// flux digitized provide just one "digitized" variable, the detector id
	abank =  gBank(FLUX_BANK_TAG, "flux", "Geant4 flux digitized information integrated over the hit");
	abank.load_variable("hitn",   99,  "Di", "Hit Number");
	abank.load_variable("id",     1,  "Di", "ID of flux element");
	abank.orderNames();
	banks["flux"] = abank;
	
	
	// Loading all banks related to a system
	// then checking that all sensitive detectors have a bank
	for(map<string, string>::iterator sit = allSystems.begin(); sit != allSystems.end(); sit++)
	{
		string systemName    = sit->first;
		string systemFactory = sit->second;
		
		if(systemName == "flux") continue;
		
		// text factory
		if(systemFactory == "TEXT")
		{
			
			string fname = systemName + "__bank.txt";
			ifstream IN(fname.c_str());
			if(!IN)
			{
				// if file is not found, maybe it's in the GEMC_DATA_DIR directory
				if(getenv("GEMC_DATA_DIR")  != NULL)
				{
					fname = (string) getenv("GEMC_DATA_DIR") + "/" + fname;
					IN.open(fname.c_str());
				}
				
			}
			// now file should be loaded
			if(IN)
			{
				if(verbosity > 1)
					cout << "   > Loading bank TEXT definitions for <" << systemName << ">." << endl;
				
				
				// first get all banks for this system
				vector<string> banksForSystem;
				while(!IN.eof())
				{
					string dbline;
					getline(IN, dbline);
					
					if(!dbline.size()) continue;
					
					gtable gt(get_strings(dbline, "|"));
					
					if(gt.data.size())
						if(gt.data[1] == "bankid")
							banksForSystem.push_back(gt.data[0]); // 0: bank name
				}
				// rewind IN
				IN.clear();
				IN.seekg(0);
				
				// now loading bank and variables
				for(unsigned b=0; b<banksForSystem.size(); b++)
				{
					while(!IN.eof())
					{
						string dbline;
						getline(IN, dbline);
						
						if(!dbline.size()) continue;
						
						gtable gt(get_strings(dbline, "|"));
						
						// the bankid entry is always the first
						if(gt.data.size())
							if(gt.data[0] == banksForSystem[b])
							{
								string bname = gt.data[0];             // 0: bank name
								string vname = gt.data[1];             // 0: variable name
								string desc  = gt.data[2];             // 2: bank/variable description
								int num      = get_number(gt.data[3]); // 3: variable num is bank id
								string type  = gt.data[4];             // 4: variable type
								
								if(vname == "bankid")
								{
									abank =  gBank(num, bname, desc);
								}
								else
								{
									abank.load_variable(vname, num, type, desc);
									
								}
							}
					}
					abank.orderNames();
					banks[banksForSystem[b]] = abank;
					IN.clear();
					IN.seekg(0);
				}
				
				IN.close();
			}
			else
			{
				if(verbosity>2)
					cout << "  !!! Error: Failed to open system bank file " << fname
					     << ". Maybe the filename doesn't exist? Exiting." << endl;
			}
		}
		
		if(systemFactory == "MYSQL")
		{
			// connection to the DB
			QSqlDatabase db = openGdb(gemcOpt);
			string tname    = systemName + "__bank";
			// hardcoding original variation for now
			string variation = "original";
			
			if(verbosity > 1) cout << "   > Loading MYSQL definitions for <" << systemName << ">." << endl;
			
			// first getting all banks from system
			string dbexecute  = "select bankname, name from " + tname ;
			dbexecute += " where variation ='" + variation + "'";
			
			QSqlQuery q;
			if(!q.exec(dbexecute.c_str()))
			{
				cout  << " !!! Failed to execute MYSQL query " << dbexecute <<  ". This is a fatal error. Exiting." << endl;
				qDebug() << q.lastError();
				exit(0);
			}
			
			// Warning if nothing is found
			if(q.size() == 0 && verbosity)
			{
				cout << "  ** WARNING: system  \"" << systemName << "\" not found in MYSQL database "
					 << " for variation " << variation << endl << endl;
			}
			
			vector<string> banksForSystem;
			while (q.next())
			{
				// Reading variables
				
				// 0: variable name
				string name = qv_tostring(q.value(1));
				
				if(name == "bankid")
					banksForSystem.push_back(qv_tostring(q.value(0)));
	
			}
			
			
			for(unsigned b=0; b<banksForSystem.size(); b++)
			{
				// re-executing the query and loading variables
				dbexecute  = "select bankname, name, num, type, description from " + tname ;
				dbexecute += " where variation ='" + variation + "'";
				dbexecute += " and bankname = '" + banksForSystem[b]  + "'";
				q.exec(dbexecute.c_str());
				while (q.next())
				{
					string bname       = qv_tostring(q.value(0));  // 0: bank name
					string vname       = qv_tostring(q.value(1));  // 1: variable name
					int num = get_number(qv_tostring(q.value(2))); // 2: variable num is bank id
					string type =        qv_tostring(q.value(3));  // 3: variable type
					string desc =        qv_tostring(q.value(4));  // 4: variable description
					
					// need to make sure bankid is the first entry
					// that comes out of the mysql query
					if(vname == "bankid")
					{
						abank =  gBank(num, bname, desc);
					}
					else
					{
						abank.load_variable(vname, num, type, desc);
					}
				}
				abank.orderNames();
				banks[banksForSystem[b]] = abank;
			
			}
			
				
			// closing DB connection
			closeGdb(db);
		}
		
	}
	if(verbosity > 3)
	{
		for(map<string, gBank>::iterator it = banks.begin(); it != banks.end(); it++)
		cout << it->second;
	}
	
	return banks;
}

// Load a variable in the bank definition
void gBank::load_variable(string n, int i, string t, string d)
{
	// adding the variable index to order the map
	name.push_back(n);
	id.push_back(i);
	type.push_back(t);
	description.push_back(d);
}


int gBank::getVarId(string bank)
{
	for(unsigned int i=0; i<name.size(); i++)
	{
		if(name[i].find(bank) == 0) return id[i];
	}
	return -1;
}



string gBank::getVarType(string var)
{
	for(unsigned int i=0; i<name.size(); i++)
	if(name[i] == var && type[i].length() == 2)
	{
		if(type[i].find("i") == 1) return "i";
		if(type[i].find("d") == 1) return "d";
	}
	
	return "na";
}

int gBank::getVarBankType(string var)
{
	for(unsigned int i=0; i<name.size(); i++)
	{
		if(name[i] == var  && type[i].length() == 2)
		{
			if(type[i].find("R") == 0) return RAWINT_ID;
			if(type[i].find("D") == 0) return DGTINT_ID;
			if(type[i].find("S") == 0) return RAWSTEP_ID;
			if(type[i].find("M") == 0) return DGTMULTI_ID;
			if(type[i].find("V") == 0) return VOLTAGETIME_ID;
		}
	}
	return 0;
}

// order names based on their ID
void gBank::orderNames()
{
	int minId = 1000;
	int maxId = 0;
	
	// first find min, max ID
	for(unsigned i=0; i<id.size(); i++)
	{
		if(id[i] < minId) minId = id[i];
		if(id[i] > maxId) maxId = id[i];
	}
	
	int j = 0;
	for(int i=minId; i<=maxId; i++)
	{
		for(unsigned k=0; k<id.size(); k++)
		{
			if(i == id[k])
			{
				orderedNames[j++] = name[k];
			}
		}
	}
}

// get bank definitions (all)
gBank getBankFromMap(string name, map<string, gBank>* banksMap)
{
	if(banksMap->find(name) == banksMap->end())
	{
		cout << "   !!! Error: >" << name << "< bank definitions not found. Exiting." << endl;
		exit(0);
	}
	
	return (*banksMap)[name];
}



// get dgt bank definitions
gBank getDgtBankFromMap(string name, map<string, gBank>* banksMap)
{
	gBank thisBank, dgtBank;
	if(banksMap->find(name) == banksMap->end())
	{
		cout << "   !!! Error: >" << name << "< bank definitions not found. Exiting." << endl;
		exit(0);
	}
	else
	{
		// thisBank may have definitions other than DGT
		// so I'm extracting just the DGT variables from it.
		thisBank = (*banksMap)[name];
		
		dgtBank = gBank(DGTINT_ID, thisBank.bankName, thisBank.bdescription);
		for(unsigned int i=0; i<thisBank.name.size(); i++)
		{
			if(thisBank.getVarBankType(thisBank.name[i]) == DGTINT_ID)
			{
				dgtBank.load_variable(thisBank.name[i], thisBank.id[i], thisBank.type[i], thisBank.description[i]);
			}
		}
		
	}
	dgtBank.orderNames();
	return dgtBank;
}

// Overloaded "<<" for the class 'bank'
ostream &operator<<(ostream &stream, gBank bank)
{
	cout << " >> Bank " << bank.bankName << " loaded with id " << bank.idtag << " : " << bank.bdescription << endl;
	
	for(unsigned i = 0; i<bank.name.size(); i++)
	{
		cout << "  > Variable : " << bank.name[i] << "\t id: " << bank.id[i] << "\t type: " << bank.type[i] << " : " << bank.description[i] << endl;
	}
	cout << endl;
	return stream;
}

