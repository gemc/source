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
#include "gemcUtils.h"

// Variable Type is two chars.
// The first char:
//  R for raw integrated
//  D for dgt integrated
//  S for raw step by step
//  M for digitized multi-hit
//  C for charge, time
//  N for not relevant

// The second char:
// i for integers
// l for long
// d for doubles
// s for strings


map <string, gBank> read_banks(goptions gemcOpt, map <string, string> allSystems) {
    double verbosity = gemcOpt.optMap["BANK_VERBOSITY"].arg;

    // geant4 information, integrated over
    map <string, gBank> banks;

    gBank abank;

    // event header
    abank = gBank(HEADER_BANK_TAG, "header", "Data Header Bank");
    abank.load_variable("time", 1, "Ns", "Time-stamp");
    abank.load_variable("runNo", 2, "Ni", "Run Number");
    abank.load_variable("evn", 3, "Ni", "Event Number");
    abank.load_variable("evn_type", 4, "Ni", "Event Type. 1 for physics events, 10 for scaler. Negative sign for MC.");
    abank.load_variable("beamPol", 5, "Nd", "Beam Polarization");
    abank.orderNames();
    banks["header"] = abank;

    // event header
    abank = gBank(USER_HEADER_BANK_TAG, "userHeader", "User Header Bank");
    abank.orderNames();
    banks["userHeader"] = abank;

    // generated particle infos
    abank = gBank(GENERATED_PARTICLES_BANK_TAG, "generated", "Generated Particles");
    abank.load_variable("pid", 1, "Ni", "Particle ID");
    abank.load_variable("px", 2, "Nd", "x component of momentum");
    abank.load_variable("py", 3, "Nd", "y component of momentum");
    abank.load_variable("pz", 4, "Nd", "z component of momentum");
    abank.load_variable("vx", 5, "Nd", "x component of vertex");
    abank.load_variable("vy", 6, "Nd", "y component of vertex");
    abank.load_variable("vz", 7, "Nd", "z component of vertex");
    abank.load_variable("time", 8, "Nd", "time");
    abank.load_variable("multiplicity", 9, "Ni", "number of particles at this vertex");
    abank.orderNames();
    banks["generated"] = abank;

    // RF info
    abank = gBank(RF_BANK_TAG, "rf", "RF Signals");
    abank.load_variable("id", 1, "Ni", "RF ID");
    abank.load_variable("rf", 2, "ND", "RF ID");
    banks["rf"] = abank;

    // particle summary infos
    // this is a daughter bank of the generated particle infos
    abank = gBank(GENERATED_SUMMARY_BANK_TAG, "psummary", "Generated Particles Summary");
    abank.load_variable("dname", 1, "Ns", "Detector Name");
    abank.load_variable("stat", 2, "Ni", "Number of Hits in the detector");
    abank.load_variable("etot", 3, "Nd", "Total Energy Deposited");
    abank.load_variable("t", 4, "Nd", "Fastest Time on Detector");
    abank.load_variable("nphe", 5, "Ni", "Number of Photoelectrons");
    abank.load_variable("upx", 10, "Nd", "unsmeared px");
    abank.load_variable("upy", 11, "Nd", "unsmeared py");
    abank.load_variable("upz", 12, "Nd", "unsmeared pz");
    abank.load_variable("spx", 13, "Nd", "smeared px");
    abank.load_variable("spy", 14, "Nd", "smeared py");
    abank.load_variable("spz", 15, "Nd", "smeared pz");
    abank.orderNames();
    banks["psummary"] = abank;


    // all hit bank information have a variable with num 0 "hitn" = hit number;

    // geant4 raw integrated
    // common for all banks
    abank = gBank(RAWINT_ID, "raws", "Geant4 true information integrated over the hit");
    abank.load_variable("pid", 1, "Ri", "ID of the first particle entering the sensitive volume");
    abank.load_variable("mpid", 2, "Ri", "ID of the mother of the first particle entering the sensitive volume");
    abank.load_variable("tid", 3, "Ri", "Track ID of the first particle entering the sensitive volume");
    abank.load_variable("mtid", 4, "Ri", "Track ID of the mother of the first particle entering the sensitive volume");
    abank.load_variable("otid", 5, "Ri", "Track ID of the original track that generated the first particle entering the sensitive volume");
    abank.load_variable("trackE", 6, "Rd", "Energy of the track");
    abank.load_variable("totEdep", 7, "Rd", "Total Energy Deposited (in MeV");
    abank.load_variable("avg_x", 8, "Rd", "Average X position in the global reference system (in mm)");
    abank.load_variable("avg_y", 9, "Rd", "Average Y position in the global reference system (in mm)");
    abank.load_variable("avg_z", 10, "Rd", "Average Z position in the global reference system (in mm)");
    abank.load_variable("avg_lx", 11, "Rd", "Average X position in the local reference system (in mm)");
    abank.load_variable("avg_ly", 12, "Rd", "Average Y position in the local reference system (in mm)");
    abank.load_variable("avg_lz", 13, "Rd", "Average Z position in the local reference system (in mm)");
    abank.load_variable("px", 14, "Rd", "x component of momentum of the particle entering the sensitive volume");
    abank.load_variable("py", 15, "Rd", "y component of momentum of the particle entering the sensitive volume");
    abank.load_variable("pz", 16, "Rd", "z component of momentum of the particle entering the sensitive volume");
    abank.load_variable("vx", 17, "Rd", "x component of point of origin of the particle entering the sensitive volume");
    abank.load_variable("vy", 18, "Rd", "y component of point of origin of the particle entering the sensitive volume");
    abank.load_variable("vz", 19, "Rd", "z component of point of origin of the particle entering the sensitive volume");
    abank.load_variable("mvx", 20, "Rd", "x component of point of origin the mother of the particle entering the sensitive volume");
    abank.load_variable("mvy", 21, "Rd", "y component of point of origin of the mother of the particle entering the sensitive volume");
    abank.load_variable("mvz", 22, "Rd", "z component of point of origin of the mother of the particle entering the sensitive volume");
    abank.load_variable("avg_t", 23, "Rd", "Average time");
    abank.load_variable("nsteps", 24, "Ri", "Number of geant4 steps");
    abank.load_variable("procID", 25, "Ri", "Process that created the particle. It's an integer described at gemc.jlab.org");
    abank.load_variable("hitn", 99, "Ri", "Hit Number");
    abank.orderNames();
    banks["raws"] = abank;


    // geant4 raw step by step
    // common for all banks
    abank = gBank(RAWSTEP_ID, "allraws", "Geant4 true information step by step");
    abank.load_variable("pid", 1, "Si", "ID of the particle in the sensitive volume");
    abank.load_variable("mpid", 2, "Si", "ID of the mother in the sensitive volume");
    abank.load_variable("tid", 3, "Si", "Track ID of the particle in the sensitive volume");
    abank.load_variable("mtid", 4, "Si", "Track ID of the mother of the particle in the sensitive volume");
    abank.load_variable("otid", 5, "Si", "Track ID of the original track that generated the particle in the sensitive volume");
    abank.load_variable("trackE", 6, "Sd", "Energy of the track");
    abank.load_variable("edep", 7, "Sd", "Energy Deposited");
    abank.load_variable("x", 8, "Sd", "X position in global reference system");
    abank.load_variable("y", 9, "Sd", "Y position in global reference system");
    abank.load_variable("z", 10, "Sd", "Z position in global reference system");
    abank.load_variable("lx", 11, "Sd", "X position in local reference system");
    abank.load_variable("ly", 12, "Sd", "Y position in local reference system");
    abank.load_variable("lz", 13, "Sd", "Z position in local reference system");
    abank.load_variable("px", 14, "Sd", "x component of momentum of the particle in the sensitive volume");
    abank.load_variable("py", 15, "Sd", "y component of momentum of the particle in the sensitive volume");
    abank.load_variable("pz", 16, "Sd", "z component of momentum of the particle in the sensitive volume");
    abank.load_variable("vx", 17, "Sd", "x component of primary vertex of the particle in the sensitive volume");
    abank.load_variable("vy", 18, "Sd", "y component of primary vertex of the particle in the sensitive volume");
    abank.load_variable("vz", 19, "Sd", "z component of primary vertex of the particle in the sensitive volume");
    abank.load_variable("mvx", 20, "Sd", "x component of primary vertex of the mother of the particle in the sensitive volume");
    abank.load_variable("mvy", 21, "Sd", "y component of primary vertex of the mother of the particle in the sensitive volume");
    abank.load_variable("mvz", 22, "Sd", "z component of primary vertex of the mother of the particle in the sensitive volume");
    abank.load_variable("t", 23, "Sd", "time");
    abank.load_variable("stepn", 98, "Si", "step index");
    abank.load_variable("hitn", 99, "Si", "Hit Number");
    abank.orderNames();
    banks["allraws"] = abank;

    // geant4 raw step by step
    // common for all banks
    abank = gBank(CHARGE_TIME_ID, "chargeTime", "charge and time as seen by the electronics");
    abank.load_variable("id", 1, "Ci", "hit identifier");
    abank.load_variable("q", 2, "Cd", "charge as seen by electronics");
    abank.load_variable("t", 3, "Cd", "time as seen by electronics");
    abank.load_variable("stepi", 98, "Ci", "step index");
    abank.load_variable("hitn", 99, "Ci", "Hit Number");
    banks["chargeTime"] = abank;


    // flux bank digitized infos
    // flux digitized provide just one "digitized" variable, the detector id
    abank = gBank(FLUX_BANK_TAG, "flux", "Geant4 flux digitized information");
    abank.load_variable("hitn", 99, "Di", "Hit Number");
    abank.load_variable("sector", 1, "Di", "ID of flux element");
    abank.load_variable("layer", 2, "Di", "ID of flux element");
    abank.load_variable("component", 3, "Di", "ID of flux element");
    abank.load_variable("ADC_order", 4, "Di", "ID of flux element");
    abank.load_variable("ADC_ADC", 5, "Di", "ID of flux element");
    abank.load_variable("ADC_time", 6, "Di", "ID of flux element");
    abank.orderNames();
    banks["flux"] = abank;

    // mirror bank digitized infos
    // mirror digitized provide just one "digitized" variable, the detector id
    abank = gBank(MIRROR_BANK_TAG, "mirror", "Geant4 mirror digitized information");
    abank.load_variable("hitn", 99, "Di", "Hit Number");
    abank.load_variable("id", 1, "Di", "ID of flux element");
    abank.orderNames();
    banks["mirror"] = abank;

    // counter bank integrated digitized infos
    // flux digitized provide just one "digitized" variable, the detector id
    abank = gBank(COUNTER_BANK_TAG, "counter", "Geant4 counter digitized information");
    abank.load_variable("id", 1, "Di", "ID of counter element");
    abank.load_variable("hitn", 99, "Di", "Hit Number");
    abank.load_variable("ngamma", 10, "Di", "number of gamma");
    abank.load_variable("nep", 11, "Di", "number of electrons");
    abank.load_variable("nem", 12, "Di", "number of positrons");
    abank.load_variable("npip", 13, "Di", "number of pi+");
    abank.load_variable("npim", 14, "Di", "number of pi-");
    abank.load_variable("npi0", 15, "Di", "number of pi0");
    abank.load_variable("nkp", 16, "Di", "number of k+");
    abank.load_variable("nkm", 17, "Di", "number of k-");
    abank.load_variable("nk0", 18, "Di", "number of k0");
    abank.load_variable("nproton", 19, "Di", "number of protons");
    abank.load_variable("nneutron", 20, "Di", "number of neutrons");
    abank.load_variable("nopticalphoton", 21, "Di", "number of optical photons");
    abank.orderNames();
    banks["counter"] = abank;


    // ancestors bank
    // Information about ancestral trajectories
    abank = gBank(ANCESTORS_BANK_TAG, "ancestors", "Geant4 ancestors information");
    abank.load_variable("pid", 1, "Ri", "ID of the ancestor");
    abank.load_variable("tid", 2, "Ri", "Track ID of the ancestor");
    abank.load_variable("mtid", 3, "Ri", "Track ID of the mother of the ancestor");
    abank.load_variable("trackE", 4, "Rd", "Energy of the ancestor");
    abank.load_variable("px", 5, "Rd", "x component of momentum of the ancestor");
    abank.load_variable("py", 6, "Rd", "y component of momentum of the ancestor");
    abank.load_variable("pz", 7, "Rd", "z component of momentum of the ancestor");
    abank.load_variable("vx", 5, "Rd", "x component of vertex of the ancestor");
    abank.load_variable("vy", 6, "Rd", "y component of vertex of the ancestor");
    abank.load_variable("vz", 7, "Rd", "z component of vertex of the ancestor");
    abank.orderNames();
    banks["ancestors"] = abank;


    // Loading all banks related to a system
    // then checking that all sensitive detectors have a bank
    for (map<string, string>::iterator sit = allSystems.begin(); sit != allSystems.end(); sit++) {
        string systemName = sit->first;
        string systemFactory = sit->second;

        // these are already loaded
        if (systemName == "flux" || systemName == "mirror" || systemName == "counter") continue;

        // text factory
        if (systemFactory == "TEXT") {
            string fname = systemName + "__bank.txt";
            ifstream IN(fname.c_str());
            if (!IN) {
                // if file is not found, maybe it's in the GEMC_DATA_DIR directory
                if (getenv("GEMC_DATA_DIR") != nullptr) {
                    fname = (string) getenv("GEMC_DATA_DIR") + "/" + fname;
                    IN.open(fname.c_str());
                }
            }

            // now file should be loaded
            if (IN) {
                if (verbosity > 1) {
                    cout << "   > Loading bank TEXT definitions for <" << systemName << ">." << endl;
                }

                // first get all banks for this system
                vector <string> banksForSystem;
                while (!IN.eof()) {
                    string dbline;
                    getline(IN, dbline);

                    if (!dbline.size()) continue;

                    gtable gt(getStringVectorFromStringWithDelimiter(dbline, "|"));

                    if (gt.data.size())
                        if (gt.data[1] == "bankid")
                            banksForSystem.push_back(gt.data[0]); // 0: bank name
                }
                // rewind IN
                IN.clear();
                IN.seekg(0);

                // now loading bank and variables
                for (unsigned b = 0; b < banksForSystem.size(); b++) {
                    while (!IN.eof()) {
                        string dbline;
                        getline(IN, dbline);

                        if (!dbline.size()) continue;

                        gtable gt(getStringVectorFromStringWithDelimiter(dbline, "|"));

                        // the bankid entry is always the first
                        if (gt.data.size())
                            if (gt.data[0] == banksForSystem[b]) {
                                string bname = gt.data[0];            // 0: bank name
                                string vname = gt.data[1];            // 1: variable name
                                string desc = gt.data[2];             // 2: bank/variable description
                                int num = get_number(gt.data[3]);     // 3: variable num is bank id
                                string type = gt.data[4];             // 4: variable type

                                // define the bank with its general ID, name and description
                                if (vname == "bankid") {
                                    abank = gBank(num, bname, desc);
                                } else {
                                    // load each variable
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
            } else {
                if (verbosity > 2) {
                    cout << "  !!! Warning: Failed to open system bank file " << fname
                         << ". Maybe the filename doesn't exist?." << endl;
                }
            }
        }

        if (systemFactory == "MYSQL") {

        }


        if (systemFactory == "SQLITE") {

            // connection to the DB
            QSqlDatabase db = openGdb(gemcOpt);
            if (verbosity > 1) cout << "   > Loading SQLITE definitions for <" << systemName << ">." << endl;

            // first select all unique bank_name for this system
            string dbexecute = "select DISTINCT bank_name from banks where system = '" + systemName + "'";

            // executing query - will exit if not successful.
            QSqlQuery q;
            if (!q.exec(dbexecute.c_str())) {
                cout  << " !!! Failed to execute SQLITE query " << dbexecute << ". This is a fatal error. Exiting." << endl;
                qDebug() << q.lastError();
                exit(1);
            }
            vector<string> banksForSystem;
           // pushing back bank names
            while (q.next()) {
                banksForSystem.push_back(qv_tostring(q.value(0)));
            }


            // now loading bank and variables
            for (unsigned b = 0; b < banksForSystem.size(); b++) {
                string bname = banksForSystem[b];
                cout << "  > Loading bank definitions for <" << bname << ">." << endl;

                string dbexecute = "select variable_name, int_id, type, description from banks where bank_name = '" + bname + "'";

                // executing query - will exit if not successful.
                QSqlQuery q;
                if (!q.exec(dbexecute.c_str())) {
                    cout << " !!! Failed to execute SQLITE query " << dbexecute << ". This is a fatal error. Exiting." << endl;
                    qDebug() << q.lastError();
                    exit(1);
                }
                // Warning if nothing is found
                if (q.size() == 0 && verbosity) {
                    cout << "  ** WARNING: bank definitions for \"" << systemName << "\" not found." << endl << endl;
                }

                while (q.next()) {
                    //string bname = qv_tostring(q.value(0));  // 0: bank (system) name
                    string vname = qv_tostring(q.value(0));         // 0: variable name
                    int num = get_number(qv_tostring(q.value(1)));  // 1: variable num is bank id
                    string type = qv_tostring(q.value(2));          // 2: variable type
                    string desc = qv_tostring(q.value(3));          // 3: variable description



                    // need to make sure bankid is the first entry
                    // that comes out of the sqlite query
                    // define the bank with its general ID, name and description
                    if (vname == "bankid") {
                        abank = gBank(num, bname, desc);
                    } else {
                        // load each variable
                        abank.load_variable(vname, num, type, desc);
                    }
                }
                abank.orderNames();
                banks[bname] = abank;
            }
        }

    }
    if (verbosity > 3) {
        for (map<string, gBank>::iterator it = banks.begin(); it != banks.end(); it++)
            cout << it->second;
    }

    return banks;
}

// Load a variable in the bank definition
void gBank::load_variable(string n, int i, string t, string d) {
    // adding the variable index to order the map
    name.push_back(n);
    gid.push_back(i);
    type.push_back(t);
    description.push_back(d);
}


int gBank::getVarId(string bank) {
    for (unsigned int i = 0; i < name.size(); i++) {
        if (name[i].find(bank) == 0) return gid[i];
    }
    return -1;
}


string gBank::getVarType(string var) {
    for (unsigned int i = 0; i < name.size(); i++) {
        if (name[i] == var && type[i].length() == 2) {
            if (type[i].find("l") == 1) return "l";
            if (type[i].find("i") == 1) return "i";
            if (type[i].find("d") == 1) return "d";
        }
    }
    return "na";
}

int gBank::getVarBankType(string var) {
    for (unsigned int i = 0; i < name.size(); i++) {
        if (name[i] == var && type[i].length() == 2) {
            if (type[i].find("R") == 0) return RAWINT_ID;
            if (type[i].find("D") == 0) return DGTINT_ID;
            if (type[i].find("S") == 0) return RAWSTEP_ID;
            if (type[i].find("M") == 0) return DGTMULTI_ID;
            if (type[i].find("C") == 0) return CHARGE_TIME_ID;
        }
    }
    return 0;
}


// order names based on their ID
void gBank::orderNames() {
    int minId = 1000;
    int maxId = 0;

    // first find min, max ID
    for (unsigned i = 0; i < gid.size(); i++) {
        if (gid[i] < minId) minId = gid[i];
        if (gid[i] > maxId) maxId = gid[i];
    }

    int j = 0;
    for (int i = minId; i <= maxId; i++) {
        for (unsigned k = 0; k < gid.size(); k++) {
            if (i == gid[k]) {
                orderedNames[j++] = name[k];
            }
        }
    }
}

// get bank definitions (all)
gBank getBankFromMap(string name, map <string, gBank> *banksMap) {
    if (banksMap->find(name) == banksMap->end()) {
        cout << "   !!! Error: >" << name << "< bank definitions not found. Exiting." << endl;
        exit(1);
    }

    return (*banksMap)[name];
}


// get dgt bank definitions
gBank getDgtBankFromMap(string name, map <string, gBank> *banksMap) {
    gBank thisBank, dgtBank;
    if (banksMap->find(name) == banksMap->end()) {
        cout << "   !!! Error: >" << name << "< dgt bank definitions not found. Exiting." << endl;
        exit(1);
    } else {
        // thisBank may have definitions other than DGT
        // so I'm extracting just the DGT variables from it.
        thisBank = (*banksMap)[name];

        dgtBank = gBank(DGTINT_ID, thisBank.bankName, thisBank.bdescription);
        for (unsigned int i = 0; i < thisBank.name.size(); i++) {
            if (thisBank.getVarBankType(thisBank.name[i]) == DGTINT_ID) {
                dgtBank.load_variable(thisBank.name[i], thisBank.gid[i], thisBank.type[i], thisBank.description[i]);
            }
        }
    }
    dgtBank.orderNames();
    return dgtBank;
}


// Overloaded "<<" for the class 'bank'
ostream &operator<<(ostream &stream, gBank bank) {
    cout << " >> Bank " << bank.bankName << " loaded with id " << bank.idtag << " : " << bank.bdescription << endl;

    for (unsigned i = 0; i < bank.name.size(); i++) {
        cout << "  > Variable : " << bank.name[i] << "\t id: " << bank.gid[i] << "\t type: " << bank.type[i] << "\t Description: " << bank.description[i] << endl;
    }
    cout << endl;
    return stream;
}
