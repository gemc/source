// gemc headers
#include "hipo_output.h"
#include "gemcUtils.h"

// C++ headers
#include <fstream>

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"

using namespace CLHEP;

// hipo4
#include "hipo4/writer.h"

// ccdb
#include <CCDB/Calibration.h>
#include <CCDB/Model/Assignment.h>
#include <CCDB/CalibrationGenerator.h>

using namespace ccdb;


map<string, double> hipo_output::fieldScales = {};
int hipo_output::rasterInitialized = -99;
double hipo_output::rasterP0[2] = {0, 0};
double hipo_output::rasterP1[2] = {0, 0};

// record the simulation conditions
// the format is a string for each variable
// the num is 0 for the mother bank, 1 for the data
void hipo_output::recordSimConditions(outputContainer *output, map <string, string> sims) {
    vector <string> data;

    // writing both key and argument as one string
    for (map<string, string>::iterator it = sims.begin(); it != sims.end(); it++) {
        if (it->first != "JSON") {
            data.push_back(it->first + ":  " + it->second + "  ");
        }
    }

    // writing both key and argument as one string
    for (map<string, string>::iterator it = sims.begin(); it != sims.end(); it++) {
        if (it->first == "option ACTIVEFIELDS") {

            vector <string> fieldNames = getStringVectorFromString(it->second);

            for (auto &fieldName: fieldNames) {

                double scaleFactor = 1;

                vector <aopt> FIELD_SCALES_OPTION = output->gemcOpt.getArgs("SCALE_FIELD");
                for (unsigned int f = 0; f < FIELD_SCALES_OPTION.size(); f++) {
                    vector <string> scales = getStringVectorFromStringWithDelimiter(FIELD_SCALES_OPTION[f].args, ",");
                    if (scales.size() == 2) {
                        if (scales[0].find(fieldName) != string::npos) {
                            scaleFactor = get_number(scales[1]);
                            // scale to 1 unless set below
                            fieldScales[trimSpacesFromString(fieldName)] = scaleFactor;
                            data.push_back("field" + trimSpacesFromString(fieldName) + " scale:  " + to_string(scaleFactor));
                        }
                    }
                }
            }

            string hardcodedBinaryTorusOptionName = "binary_torus";
            string hardcodedBinarySolenOptionName = "binary_solenoid";
            double scaleFactor = 1;
            vector <aopt> FIELD_SCALES_OPTION = output->gemcOpt.getArgs("SCALE_FIELD");
            for (unsigned int f = 0; f < FIELD_SCALES_OPTION.size(); f++) {
                vector <string> scales = getStringVectorFromStringWithDelimiter(FIELD_SCALES_OPTION[f].args, ",");
                if (scales.size() == 2) {
                    if (scales[0].find(hardcodedBinaryTorusOptionName) != string::npos) {
                        scaleFactor = get_number(scales[1]);
                        // scale to 1 unless set below
                        fieldScales[trimSpacesFromString(hardcodedBinaryTorusOptionName)] = scaleFactor;
                        data.push_back("field" + trimSpacesFromString(hardcodedBinaryTorusOptionName) + " scale:  " + to_string(scaleFactor));
                    } else if (scales[0].find(hardcodedBinarySolenOptionName) != string::npos) {
                        scaleFactor = get_number(scales[1]);
                        // scale to 1 unless set below
                        fieldScales[trimSpacesFromString(hardcodedBinarySolenOptionName)] = scaleFactor;
                        data.push_back("field" + trimSpacesFromString(hardcodedBinaryTorusOptionName) + " scale:  " + to_string(scaleFactor));
                    }
                }
            }

        }
    }

    string bigData;
    for (auto b: data) {
        bigData += b + string("\n");
    }

    // file need to be opened after user configuration is added
    // output->hipoWriter->addUserConfig("GEMC::config",  bigData);

    output->initializeHipo(true);


}

// returns detectorID from map, given hitType
int hipo_output::getDetectorID(string hitType) {

    if (detectorID.find(hitType) != detectorID.end()) {
        return detectorID[hitType];
    } else {
        cout << " Error: " << hitType << " has no detector id. Exiting. " << endl;
        exit(601);
    }
}

// returns hipo name from true info var name
string hipo_output::getHipoVariableName(string trueInfoVar) {

    if (trueInfoNamesMap.find(trueInfoVar) != trueInfoNamesMap.end()) {
        return trueInfoNamesMap[trueInfoVar];
    } else {
        return trueInfoVar;
    }
}


// instantiates hipo event
// write run::config bank
void hipo_output::writeHeader(outputContainer *output, map<string, double> data, gBank bank) {
    int verbosity = int(output->gemcOpt.optMap["BANK_VERBOSITY"].arg);

    //	for(auto &fieldScale: fieldScales) {
    //		cout << ">" << fieldScale.first << "<" << " scaled by: " << fieldScale.second << endl;
    //	}

    // this will never enter the second condition because hipo_output is instantiated every event
    if (outEvent == nullptr) {
        outEvent = new hipo::event(1024 * 1024 * 2);
        //	cout << " Event Size before reset: " << outEvent->getSize() << endl;
    } else {
        outEvent->reset();
        //	cout << " Event Size after reset: " << outEvent->getSize() << endl;
    }

    // Create runConfigBank with 1 row based on schema
    // second argument is number of hits
    hipo::bank runConfigBank(output->hipoSchema->runConfigSchema, 1);

    runConfigBank.putInt("run", 0, data["runNo"]);
    runConfigBank.putInt("event", 0, data["evn"]);

    // time in seconds
    time_t t = std::time(0);
    int now = static_cast<int> (t);
    runConfigBank.putInt("unixtime", 0, now);

    // other infos
    runConfigBank.putInt("trigger", 0, 0);
    runConfigBank.putFloat("timestamp", 0, 0);
    runConfigBank.putInt("type", 0, 0);

    // solenoid and torus scales, if present
    if (fieldScales.find("TorusSymmetric") != fieldScales.end()) {
        runConfigBank.putFloat("torus", 0, fieldScales["TorusSymmetric"]);
        //		cout << "TorusSymmetric scaled by: " << fieldScales["TorusSymmetric"] << endl;
    } else if (fieldScales.find("binary_torus") != fieldScales.end()) {
        runConfigBank.putFloat("torus", 0, fieldScales["binary_torus"]);
    } else {
        runConfigBank.putFloat("torus", 0, 0);
    }
    if (fieldScales.find("clas12-newSolenoid") != fieldScales.end()) {
        runConfigBank.putFloat("solenoid", 0, fieldScales["clas12-newSolenoid"]);
    } else if (fieldScales.find("binary_solenoid") != fieldScales.end()) {
        runConfigBank.putFloat("solenoid", 0, fieldScales["binary_solenoid"]);
    } else {
        runConfigBank.putFloat("solenoid", 0, 0);
    }


    if (verbosity > 2) {
        runConfigBank.show();
    }

    outEvent->addStructure(runConfigBank);


    // RASTER constant initialization
    if (rasterInitialized == -99) {

        // loading raster p0 and p1 from CCDB
        string digiVariation = output->gemcOpt.optMap["DIGITIZATION_VARIATION"].args;
        int runno = output->gemcOpt.optMap["RUNNO"].arg;
        string database = "/calibration/raster/adc_to_position";

        string connection = "mysql://clas12reader@clasdb.jlab.org/clas12";

        if (getenv("CCDB_CONNECTION") != nullptr) {
            connection = (string) getenv("CCDB_CONNECTION");
        }

        vector <vector<double>> dbdata;

        unique_ptr <Calibration> calib(CalibrationGenerator::CreateCalibration(connection));
        database = database + ":" + to_string(runno) + ":" + digiVariation;
        cout << " Connecting to " << connection << database << " to retrive raster parameters" << endl;

        dbdata.clear();
        calib->GetCalib(dbdata, database);

        if (dbdata.size() == 2) {
            for (unsigned row = 0; row < dbdata.size(); row++) {
                rasterP0[row] = dbdata[row][3];
                rasterP1[row] = dbdata[row][4];
            }

            cout << " Raster Parameters: p0(x,y) = (" << rasterP0[0] << ", " << rasterP0[1] << "), p1(x, y) = (" << rasterP1[0] << ", " << rasterP1[1] << ")" << endl;
        }
        rasterInitialized = 1;
    }

}


// write user infos header
// this include the LUND file header entries (10) plus user entries if present
void hipo_output::writeUserInfoseHeader(outputContainer *output, map<string, double> data) {

    if (data.size() > 0) {

        int lundStandardSize = 10;

        int userBankSize = data.size() - lundStandardSize;

        hipo::bank mcEventHeaderBank(output->hipoSchema->mcEventHeader, 1);
        hipo::bank userLundBank(output->hipoSchema->userLund, userBankSize);

        int index = -1;
        for (auto it = data.begin(); it != data.end(); it++) {
            index += 1;
            if (index < lundStandardSize) {
                switch (index) {
                    case 0:
                        mcEventHeaderBank.putShort("npart", 0, (short) it->second);
                        break;
                    case 1:
                        mcEventHeaderBank.putShort("atarget", 0, (short) it->second);
                        break;
                    case 2:
                        mcEventHeaderBank.putShort("ztarget", 0, (short) it->second);
                        break;
                    case 3:
                        mcEventHeaderBank.putFloat("ptarget", 0, (float) it->second);
                        break;
                    case 4:
                        mcEventHeaderBank.putFloat("pbeam", 0, (float) it->second);
                        break;
                    case 5:
                        mcEventHeaderBank.putShort("btype", 0, (short) it->second);
                        break;
                    case 6:
                        mcEventHeaderBank.putFloat("ebeam", 0, (float) it->second);
                        break;
                    case 7:
                        mcEventHeaderBank.putShort("targetid", 0, (short) it->second);
                        break;
                    case 8:
                        mcEventHeaderBank.putShort("processid", 0, (short) it->second);
                        break;
                    case 9:
                        mcEventHeaderBank.putFloat("weight", 0, (float) it->second);
                        break;

                    default:
                        break;
                }

            } else {
                userLundBank.putFloat("userVar", index - lundStandardSize, (float) it->second);
            }
        }

        outEvent->addStructure(mcEventHeaderBank);
        if (userBankSize > 0) {
            outEvent->addStructure(userLundBank);
        }
    }
}


void hipo_output::writeRFSignal(outputContainer *output, FrequencySyncSignal rfsignals, gBank bank) {
    int verbosity = int(output->gemcOpt.optMap["BANK_VERBOSITY"].arg);

    // Create runRFBank with 1 row based on schema
    // second argument is number of hits

    vector <oneRFOutput> rfs = rfsignals.getOutput();

    // only the first rf output is written
    // beware: in gemc there are 2 rf outputs
    vector<int> ids = rfs.front().getIDs();
    vector<double> times = rfs.front().getValues();

    hipo::bank runRFBank(output->hipoSchema->runRFSchema, ids.size());

    for (unsigned i = 0; i < ids.size(); i++) {
        runRFBank.putShort("id", i, (short) ids[i]);
        runRFBank.putFloat("time", i, (float) times[i]);
    }

    if (verbosity > 2) {
        runRFBank.show();
    }

    outEvent->addStructure(runRFBank);
}

void hipo_output::writeGenerated(outputContainer *output, vector <generatedParticle> MGP, map <string, gBank> *banksMap, vector <userInforForParticle> userInfo) {
    double MAXP = output->gemcOpt.optMap["NGENP"].arg;
    int verbosity = int(output->gemcOpt.optMap["BANK_VERBOSITY"].arg);

    vector<int> pid;
    vector<double> px;
    vector<double> py;
    vector<double> pz;
    vector<double> vx;
    vector<double> vy;
    vector<double> vz;
    vector<double> btime;

    for (unsigned i = 0; i < MAXP && i < MGP.size(); i++) {

        int my_pid = MGP[i].PID;

        // more user friendly values
        if (my_pid == 1000010020) {
            my_pid = 45;  // deuteron
        } else if (my_pid == 1000010030) {
            my_pid = 46;  // triton
        } else if (my_pid == 1000020040) {
            my_pid = 47;  // alpha
        } else if (my_pid == 1000020030) {
            my_pid = 49;  // He3
        }

        pid.push_back(my_pid);

        px.push_back(MGP[i].momentum.getX() / MeV);
        py.push_back(MGP[i].momentum.getY() / MeV);
        pz.push_back(MGP[i].momentum.getZ() / MeV);
        vx.push_back(MGP[i].vertex.getX() / mm);
        vy.push_back(MGP[i].vertex.getY() / mm);
        vz.push_back(MGP[i].vertex.getZ() / mm);
        btime.push_back(MGP[i].time);
    }

    hipo::bank geantParticleBank(output->hipoSchema->geantParticle, pid.size());

    for (unsigned i = 0; i < pid.size(); i++) {

        geantParticleBank.putInt("pid", i, pid[i]);
        geantParticleBank.putFloat("px", i, (float) px[i] / 1000.0);  // in hipo the units are GeV
        geantParticleBank.putFloat("py", i, (float) py[i] / 1000.0);  // in hipo the units are GeV
        geantParticleBank.putFloat("pz", i, (float) pz[i] / 1000.0);  // in hipo the units are GeV
        geantParticleBank.putFloat("vx", i, (float) vx[i] / 10.0);    // in hipo the units are cm
        geantParticleBank.putFloat("vy", i, (float) vy[i] / 10.0);    // in hipo the units are cm
        geantParticleBank.putFloat("vz", i, (float) vz[i] / 10.0);    // in hipo the units are cm
        geantParticleBank.putFloat("vt", i, (float) btime[i]);

    }

    if (verbosity > 2) {
        geantParticleBank.show();
    }

    outEvent->addStructure(geantParticleBank);


    hipo::bank lundParticleBank(output->hipoSchema->lundParticle, userInfo.size());

    // p is particle index
    for (unsigned p = 0; p < userInfo.size(); p++) {
        for (unsigned u = 0; u < userInfo[p].infos.size(); u++) {
            switch (u) {
                case 0:
                    lundParticleBank.putByte("index", p, userInfo[p].infos[u]);
                    break;
                case 1:
                    lundParticleBank.putFloat("lifetime", p, (float) userInfo[p].infos[u]);
                    break;
                case 2:
                    lundParticleBank.putByte("type", p, userInfo[p].infos[u]);
                    break;
                case 3:
                    lundParticleBank.putInt("pid", p, (int) userInfo[p].infos[u]);
                    break;
                case 4:
                    lundParticleBank.putByte("parent", p, userInfo[p].infos[u]);
                    break;
                case 5:
                    lundParticleBank.putByte("daughter", p, userInfo[p].infos[u]);
                    break;
                case 6:
                    lundParticleBank.putFloat("px", p, (float) userInfo[p].infos[u]);
                    break;
                case 7:
                    lundParticleBank.putFloat("py", p, (float) userInfo[p].infos[u]);
                    break;
                case 8:
                    lundParticleBank.putFloat("pz", p, (float) userInfo[p].infos[u]);
                    break;
                case 9:
                    lundParticleBank.putFloat("energy", p, (float) userInfo[p].infos[u]);
                    break;
                case 10:
                    lundParticleBank.putFloat("mass", p, (float) userInfo[p].infos[u]);
                    break;
                case 11:
                    lundParticleBank.putFloat("vx", p, (float) userInfo[p].infos[u]);
                    break;
                case 12:
                    lundParticleBank.putFloat("vy", p, (float) userInfo[p].infos[u]);
                    break;
                case 13:
                    lundParticleBank.putFloat("vz", p, (float) userInfo[p].infos[u]);
                    break;

                default:
                    break;
            }

        }
    }
    outEvent->addStructure(lundParticleBank);



    // raster:
    // sector=0
    // layer=0
    // order=0
    // ADC=0
    // time=0
    //
    // given vx, vy of the first particle
    // component = 1=vx 2=vy
    // vx = p0(0) + p1(0)*pedestal
    // vy = p0(1) + p1(1)*pedestal
    // ped = (vx - p0) / p1
    // p0, p1 from  /calibration/raster/adc_to_position


    short components[2] = {1, 2};
    short peds[2];

    peds[0] = (short) ((vx[0] / cm - rasterP0[0]) / rasterP1[0]);
    peds[1] = (short) ((vy[0] / cm - rasterP0[1]) / rasterP1[1]);

    hipo::bank rasterBank(output->hipoSchema->rasterADCSchema, 2);

    // zero var infos
    for (int j = 0; j < 2; j++) {
        rasterBank.putByte("sector", j, 0);
        rasterBank.putByte("layer", j, 0);
        rasterBank.putShort("component", j, components[j]);
        rasterBank.putByte("order", j, 0);
        rasterBank.putInt("ADC", j, 0);
        rasterBank.putFloat("time", j, 0);
        rasterBank.putShort("ped", j, peds[j]);
    }
    outEvent->addStructure(rasterBank);

}

void hipo_output::writeAncestors(outputContainer *output, vector <ancestorInfo> ainfo, gBank bank) {
    vector<int> pid;
    vector<int> tid;
    vector<int> mtid;
    vector<double> trackE;
    vector<double> px;
    vector<double> py;
    vector<double> pz;
    vector<double> vx;
    vector<double> vy;
    vector<double> vz;

    for (unsigned i = 0; i < ainfo.size(); i++) {
        pid.push_back(ainfo[i].pid);
        tid.push_back(ainfo[i].tid);
        mtid.push_back(ainfo[i].mtid);
        trackE.push_back(ainfo[i].trackE);
        px.push_back(ainfo[i].p.getX() / MeV);
        py.push_back(ainfo[i].p.getY() / MeV);
        pz.push_back(ainfo[i].p.getZ() / MeV);
        vx.push_back(ainfo[i].vtx.getX() / MeV);
        vy.push_back(ainfo[i].vtx.getY() / MeV);
        vz.push_back(ainfo[i].vtx.getZ() / MeV);
    }


}

void hipo_output::initBank(outputContainer *output, gBank thisHitBank, int what) {

}

void hipo_output::prepareEvent(outputContainer *output, map<string, double> *configuration) {
    int verbosity = int(output->gemcOpt.optMap["BANK_VERBOSITY"].arg);
    int nBankEntries = 0;
    lastHipoTrueInfoBankIndex = 0;

    for (auto &conf: *configuration) {

        if (verbosity > 1) {
            cout << " Hipo preparing " << conf.first << "  with " << conf.second << " hits " << endl;
        }
        nBankEntries = nBankEntries + conf.second;
    }

    if (verbosity > 1) {
        cout << " Total true info bank entries: " << nBankEntries << endl;
    }

    hipo::schema trueInfoSchema = output->hipoSchema->trueInfoSchema;

    trueInfoBank = new hipo::bank(trueInfoSchema, nBankEntries);

}

void hipo_output::writeG4RawIntegrated(outputContainer *output, vector <hitOutput> HO, string hitType, map <string, gBank> *banksMap) {
    if (HO.size() == 0) return;
    int verbosity = int(output->gemcOpt.optMap["BANK_VERBOSITY"].arg);

    gBank thisHitBank = getBankFromMap(hitType, banksMap);
    gBank rawBank = getBankFromMap("raws", banksMap);

    // perform initializations if necessary
    initBank(output, thisHitBank, RAWINT_ID);

    // we only need the first hit to get the definitions
    map<string, double> raws = HO[0].getRaws();

    int detectorID = getDetectorID(hitType);

    // looping over the loaded banknames (this to make sure we only publish the ones declared
    // we actually never used this steps and it's actually cumbersome. To be removed in gemc3
    for (auto &bankName: rawBank.orderedNames) {

        string bname = bankName.second;
        int bankId = rawBank.getVarId(bname);       // bankId is num
        int bankType = rawBank.getVarBankType(bname); // bankType: 1 = raw 2 = dgt

        if (raws.find(bname) != raws.end() && bankId > 0 && bankType == RAWINT_ID) {

            // looping over the hits
            for (unsigned int nh = 0; nh < HO.size(); nh++) {

                int hipoBankIndex = lastHipoTrueInfoBankIndex + nh;

                map<string, double> theseRaws = HO[nh].getRaws();

                trueInfoBank->putByte("detector", hipoBankIndex, detectorID);

                for (auto &thisVar: theseRaws) {

                    // found data match to bank definition
                    if (thisVar.first == bname) {

                        string varType = rawBank.getVarType(thisVar.first);

                        string hipoName = getHipoVariableName(bname);

                        if (hipoName == "hitn") {
                            thisVar.second = nh + 1;
                        }

                        if (varType == "i") {
                            trueInfoBank->putInt(hipoName.c_str(), hipoBankIndex, thisVar.second);
                        } else if (varType == "d") {
                            trueInfoBank->putFloat(hipoName.c_str(), hipoBankIndex, thisVar.second);
                        }

                        if (verbosity > 2) {
                            cout << " Hit Type: " << hitType << ", detector id: " << detectorID << ", hit index " << nh << ", bank hit index " << hipoBankIndex << ", name " << bname << ", hname "
                                 << hipoName << ", value: " << thisVar.second << ", raw/dgt: " << bankType << ", type: " << varType << endl;
                        }
                    }
                }
            }

        }
    }
    lastHipoTrueInfoBankIndex = lastHipoTrueInfoBankIndex + HO.size();


    // loop over variable names
    for (map<int, string>::iterator it = rawBank.orderedNames.begin(); it != rawBank.orderedNames.end(); it++) {

        int bankId = rawBank.getVarId(it->second);
        int bankType = rawBank.getVarBankType(it->second);

        // we only need the first hit to get the definitions
        map<string, double> raws = HO[0].getRaws();

        if (raws.find(it->second) != raws.end() && bankId > 0 && bankType == RAWINT_ID) {
            vector<double> thisVar;
            for (unsigned int nh = 0; nh < HO.size(); nh++) {
                map<string, double> theseRaws = HO[nh].getRaws();
                thisVar.push_back(theseRaws[it->second]);
            }
        }
    }
}


void hipo_output::writeG4DgtIntegrated(outputContainer *output, vector <hitOutput> HO, string hitType, map <string, gBank> *banksMap) {
    if (HO.size() == 0) return;
    int verbosity = int(output->gemcOpt.optMap["BANK_VERBOSITY"].arg);

    gBank thisHitBank = getBankFromMap(hitType, banksMap);
    gBank dgtBank = getDgtBankFromMap(hitType, banksMap);

    // perform initializations if necessary
    initBank(output, thisHitBank, DGTINT_ID);

    // we only need the first hit to get the definitions
    map<string, double> dgts = HO[0].getDgtz();

    bool hasADCBank = false;
    bool hasTDCBank = false;
    bool hasWF136Bank = false;

    hipo::schema detectorADCSchema = output->hipoSchema->getSchema(hitType, 0);
    hipo::schema detectorTDCSchema = output->hipoSchema->getSchema(hitType, 1);
    hipo::schema detectorWF136Schema = output->hipoSchema->getSchema(hitType, 2);
    hipo::bank detectorADCBank(detectorADCSchema, HO.size());
    hipo::bank detectorTDCBank(detectorTDCSchema, HO.size());
    hipo::bank detectorWF136Bank(detectorWF136Schema, HO.size());

    // check if there is at least one adc or tdc var
    // and if the schema is valid
    for (auto &bankName: dgtBank.orderedNames) {

        // flag ADC content if any variable has ADC_ prefix AND detectorADCSchema exists
        if (bankName.second.find("ADC_") != string::npos) {
            if (detectorADCSchema.getEntryName(0) != "empty") {
                hasADCBank = true;
            }
        }

        // flag TDC content if any variable has TDC_ prefix detectorTDCSchema schema exists
        if (bankName.second.find("TDC_") != string::npos) {

            if (detectorTDCSchema.getEntryName(0) != "empty") {
                hasTDCBank = true;
            }
        }

        // flag WF136 content if any variable has WF136_ prefix detectorWF136Schema schema exists
        if (bankName.second.find("WF136_") != string::npos) {

            if (detectorWF136Schema.getEntryName(0) != "empty") {
                hasWF136Bank = true;
            }
        }
    }

    if (verbosity > 2) {
        if (hasADCBank) {
            cout << hitType << " has ADC bank." << endl;
        }
        if (hasTDCBank) {
            cout << hitType << " has TDC bank." << endl;
        }
        if (hasWF136Bank) {
            cout << hitType << " has WF136 bank." << endl;
        }
    }


    // looping over the loaded banknames (this to make sure we only publish the ones declared
    // we actually never used this steps and it's actually cumbersome. To be removed in gemc3
    for (auto &bankName: dgtBank.orderedNames) {

        string bname = bankName.second;
        int bankId = dgtBank.getVarId(bname);       // bankId is num
        int bankType = dgtBank.getVarBankType(bname); // bankType: 1 = raw 2 = dgt

        if (dgts.find(bname) != dgts.end() && bankId > 0 && bankType == DGTINT_ID) {

            if (hasADCBank) {

                // looping over the hits
                for (unsigned int nh = 0; nh < HO.size(); nh++) {

                    map<string, double> theseDgts = HO[nh].getDgtz();
                    for (auto &thisVar: theseDgts) {

                        // found data match to bank definition
                        if (thisVar.first == bname) {

                            string varType = dgtBank.getVarType(thisVar.first);

                            // sector, layer, component are common in adc/tdc so their names are w/o prefix
                            // sector, layers are "Bytes"
                            if (bname == "sector" || bname == "layer" || bname == "order") {
                                detectorADCBank.putByte(bname.c_str(), nh, thisVar.second);
                            } else if (bname == "component") {
                                detectorADCBank.putShort(bname.c_str(), nh, thisVar.second);
                            } else {
                                // all other ADC vars must begin with "ADC_"
                                if (bname.find("ADC_") == 0) {
                                    string adcName = bname.substr(4);
                                    if (varType == "i") {
                                        detectorADCBank.putInt(adcName.c_str(), nh, thisVar.second);
                                    } else if (varType == "d") {
                                        detectorADCBank.putFloat(adcName.c_str(), nh, thisVar.second);
                                    } else if (varType == "l") {
                                        detectorADCBank.putLong(adcName.c_str(), nh, thisVar.second);
                                    }
                                }
                            }

                            if (verbosity > 2) {
                                cout << "hit index " << nh << ", name " << bname << ", value: " << thisVar.second << ", raw/dgt: " << bankType << ", type: " << varType << endl;
                            }
                        }
                    }
                }

            }

            if (hasTDCBank) {

                // looping over the hits
                for (unsigned int nh = 0; nh < HO.size(); nh++) {
                    map<string, double> theseDgts = HO[nh].getDgtz();
                    for (auto &thisVar: theseDgts) {

                        // found data match to bank definition
                        if (thisVar.first == bname) {

                            string varType = dgtBank.getVarType(thisVar.first);

                            // sector, layer, component are common in adc/tdc so their names are w/o prefix
                            // sector, layers are "Bytes"
                            if (bname == "sector" || bname == "layer" || bname == "order") {
                                detectorTDCBank.putByte(bname.c_str(), nh, thisVar.second);
                            } else if (bname == "component") {
                                detectorTDCBank.putShort(bname.c_str(), nh, thisVar.second);
                            } else {
                                // all other TDC vars must begin with "ADC_"
                                if (bname.find("TDC_") == 0) {
                                    string adcName = bname.substr(4);
                                    if (varType == "i") {
                                        detectorTDCBank.putInt(adcName.c_str(), nh, thisVar.second);
                                    } else if (varType == "d") {
                                        detectorTDCBank.putFloat(adcName.c_str(), nh, thisVar.second);
                                    } else if (varType == "l") {
                                        detectorTDCBank.putLong(adcName.c_str(), nh, thisVar.second);
                                    }
                                }

                            }

                            // do not repeat message logged above if hasADCBank was true
                            if (verbosity > 2 && !hasADCBank) {
                                cout << "hit index " << nh << ", name " << thisVar.first << ", value: " << thisVar.second << ", raw/dgt: " << bankType << ", type: " << varType << endl;
                            }
                        }
                    }
                }
            }


            if (hasWF136Bank) {
                // looping over the hits
                for (unsigned int nh = 0; nh < HO.size(); nh++) {

                    map<string, double> theseDgts = HO[nh].getDgtz();
                    for (auto &thisVar: theseDgts) {

                        // found data match to bank definition
                        if (thisVar.first == bname) {

                            string varType = dgtBank.getVarType(thisVar.first);

                            // sector, layer, component are common in adc/tdc so their names are w/o prefix
                            // sector, layers are "Bytes"
                            if (bname == "sector" || bname == "layer" || bname == "order") {
                                detectorWF136Bank.putByte(bname.c_str(), nh, thisVar.second);
                            } else if (bname == "component") {
                                detectorWF136Bank.putShort(bname.c_str(), nh, thisVar.second);
                            } else if (bname == "WF136_timestamp") {
                                detectorWF136Bank.putLong("timestamp", nh, thisVar.second);
                            } else {
                                // all other ADC vars must begin with "ADC_"
                                if (bname.find("WF136_s") == 0) {
                                    // sample number is the string following "WF136_s" converted to int
                                    int sample_value = stoi(bname.substr(7));
                                    string wfname = "s" + to_string(sample_value);
                                    detectorWF136Bank.putShort(wfname.c_str(), nh, thisVar.second);
                                }
                            }

                            if (verbosity > 2) {
                                cout << "hit index " << nh << ", name " << bname << ", value: " << thisVar.second << ", raw/dgt: " << bankType << ", type: " << varType << endl;
                            }
                        }
                    }
                }
            }

        }
    }
    if (hasADCBank) {
        if (verbosity > 2) {
            detectorADCBank.show();
        }
        outEvent->addStructure(detectorADCBank);
    }
    if (hasTDCBank) {
        if (verbosity > 2) {
            detectorTDCBank.show();
        }
        outEvent->addStructure(detectorTDCBank);
    }
    if (hasWF136Bank) {
        if (verbosity > 2) {
            detectorWF136Bank.show();
        }
        outEvent->addStructure(detectorWF136Bank);
    }

}


// index 0: hit number
// index 1: step index
// index 2: charge at electronics
// index 3: time at electronics
// index 4: vector of identifiers - have to match the translation table
void hipo_output::writeChargeTime(outputContainer *output, vector <hitOutput> HO, string hitType, map <string, gBank> *banksMap) {
    if (HO.size() == 0) return;

    gBank thisHitBank = getBankFromMap(hitType, banksMap);
    gBank chargeTimeBank = getBankFromMap("chargeTime", banksMap);

    initBank(output, thisHitBank, CHARGE_TIME_ID);

    // collecting vectors from each hit into one big array
    vector<double> allHitN;
    vector<double> allStep;
    vector<double> allCharge;
    vector<double> allTime;
    vector<double> allID;

    for (unsigned int nh = 0; nh < HO.size(); nh++) {

        map<int, vector<double> > thisChargeTime = HO[nh].getChargeTime();

        vector<double> thisHitN = thisChargeTime[0];
        vector<double> thisStep = thisChargeTime[1];
        vector<double> thisCharge = thisChargeTime[2];
        vector<double> thisTime = thisChargeTime[3];
        vector<double> thisID = thisChargeTime[4];


        // hit number
        if (thisHitN.size() != 1) {
            cout << "  !! Error: hit number should not be a vector. Bank: " << hitType << endl;
            exit(1);
        }

        for (auto h: thisHitN) allHitN.push_back(h);
        for (auto h: thisStep) allStep.push_back(h);
        for (auto h: thisCharge) allCharge.push_back(h);
        for (auto h: thisTime) allTime.push_back(h);
        for (auto h: thisID) allID.push_back(h);
    }


}


void hipo_output::writeG4RawAll(outputContainer *output, vector <hitOutput> HO, string hitType, map <string, gBank> *banksMap) {
}

void hipo_output::writeFADCMode1(map<int, vector<hitOutput> > HO, int ev_number) {
}


void hipo_output::writeFADCMode1(outputContainer *output, vector <hitOutput> HO, int ev_number) {
}


// write fadc mode 7 (integrated mode) - jlab hybrid banks. This uses the translation table to write the crate/slot/channel
void hipo_output::writeFADCMode7(outputContainer *output, vector <hitOutput> HO, int ev_number) {
}

void hipo_output::writeEvent(outputContainer *output) {
    outEvent->addStructure(*trueInfoBank);

    output->hipoWriter->addEvent(*outEvent);

}
