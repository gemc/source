// gemc headers
#include "outputFactory.h"
#include "hipoSchemas.h"

HipoSchema::HipoSchema() {

    cout << " Defining Hipo4 schemas..." << endl;
    //---------------------------------------------------------
    // Define a Schema for particle bank and detector bank
    // Schema constructor takes:
    // name - name of the schema (bank)
    // groupid - a 16 bit number which identifies a group.
    //           this way banks can be grouped in skimming
    // itemid  - and item number (8 bit) to just have a
    //           unique number associated with bank.

    // CLAS12 schemas are defined at https://github.com/JeffersonLab/clas12-offline-software/tree/development/etc/bankdefs/hipo4
    // Detectors: https://github.com/JeffersonLab/clas12-offline-software/blob/development/etc/bankdefs/hipo4/data.json
    runConfigSchema = hipo::schema("RUN::config", 10000, 11);
    runRFSchema = hipo::schema("RUN::rf", 10000, 12);
    trueInfoSchema = hipo::schema("MC::True", 40, 4);

    // generators
    geantParticle = hipo::schema("MC::Particle", 40, 2);
    mcEventHeader = hipo::schema("MC::Event", 40, 1);
    userLund = hipo::schema("MC::User", 40, 5);
    lundParticle = hipo::schema("MC::Lund", 40, 3);

    // flux
    fluxADCSchema = hipo::schema("FLUX::adc", 22200, 20);


    // detectors
    alertAhdcTDCSchema = hipo::schema("AHDC::tdc", 22400, 12);
    alertAhdcADCSchema = hipo::schema("AHDC::adc", 22400, 11);
    alertAhdcWF10Schema = hipo::schema("AHDC::wf:10", 22400, 13);
    alertAtofADCSchema = hipo::schema("ATOF::adc", 22500, 11);
    bandADCSchema = hipo::schema("BAND::adc", 22100, 11);
    bandTDCSchema = hipo::schema("BAND::tdc", 22100, 12);
    bmtADCSchema = hipo::schema("BMT::adc", 20100, 11);
    bstADCSchema = hipo::schema("BST::adc", 20200, 11);
    cndADCSchema = hipo::schema("CND::adc", 20300, 11);
    cndTDCSchema = hipo::schema("CND::tdc", 20300, 12);
    ctofADCSchema = hipo::schema("CTOF::adc", 20400, 11);
    ctofTDCSchema = hipo::schema("CTOF::tdc", 20400, 12);
    dcTDCSchema = hipo::schema("DC::tdc", 20600, 12);
    dcDOCASchema = hipo::schema("DC::doca", 20600, 14);
    ecalADCSchema = hipo::schema("ECAL::adc", 20700, 11);
    ecalTDCSchema = hipo::schema("ECAL::tdc", 20700, 12);
    fmtADCSchema = hipo::schema("FMT::adc", 20800, 11);
    ftcalADCSchema = hipo::schema("FTCAL::adc", 21000, 11);
    fthodoADCSchema = hipo::schema("FTHODO::adc", 21100, 11);
    ftofADCSchema = hipo::schema("FTOF::adc", 21200, 11);
    ftofTDCSchema = hipo::schema("FTOF::tdc", 21200, 12);
    ftrkTDCSchema = hipo::schema("FTTRK::adc", 21300, 11);
    htccADCSchema = hipo::schema("HTCC::adc", 21500, 11);
    htccTDCSchema = hipo::schema("HTCC::tdc", 21500, 12);
    ltccADCSchema = hipo::schema("LTCC::adc", 21600, 11);
    ltccTDCSchema = hipo::schema("LTCC::tdc", 21600, 12);
    rfADCSchema = hipo::schema("RF::adc", 21700, 11);
    rfTDCSchema = hipo::schema("RF::adc", 21700, 12);
    richTDCSchema = hipo::schema("RICH::tdc", 21800, 12);
    rtpcADCSchema = hipo::schema("RTPC::adc", 21900, 11);
    rtpcPOSSchema = hipo::schema("RTPC::pos", 21900, 14);
    helADCSchema = hipo::schema("HEL::adc", 22000, 11);
    helFLIPSchema = hipo::schema("HEL::flip", 22000, 12);
    helONLINESchema = hipo::schema("HEL::online", 22000, 13);
    urwellADCSchema = hipo::schema("URWELL::adc", 22300, 11);
    rawADCSchema = hipo::schema("RAW::adc", 20000, 11);
    rawTDCSchema = hipo::schema("RAW::tdc", 20000, 12);
    rawSCALERSchema = hipo::schema("RAW::scaler", 20000, 13);
    rawVTPSchema = hipo::schema("RAW::vtp", 20000, 14);
    rawEPICSSchema = hipo::schema("RAW::epics", 20000, 15);
    rasterADCSchema = hipo::schema("RASTER::adc", 22200, 11);

    // Defining structure of the schema (bank)
    // The columns in the banks (or leafs, if you like ROOT)
    // are given comma separated with type after the name.
    // Available types : I-integer, S-short, B-byte, F-float,
    //                   D-double, L-long

    runConfigSchema.parse("run/I, event/I, unixtime/I, trigger/L, timestamp/L, type/B,mode/B, torus/F, solenoid/F");
    runRFSchema.parse("id/S, time/F");
    trueInfoSchema.parse(
            "detector/B, pid/I, mpid/I, tid/I, mtid/I, otid/I, trackE/F, totEdep/F, avgX/F, avgY/F, avgZ/F, avgLx/F, avgLy/F, avgLz/F, px/F, py/F, pz/F, vx/F, vy/F, vz/F, mvx/F, mvy/F, mvz/F, avgT/F, nsteps/I, procID/I, hitn/I");
    rasterADCSchema.parse("sector/B, layer/B, component/S, order/B, ADC/I, time/F, ped/S");

    // generators
    geantParticle.parse("pid/I, px/F, py/F, pz/F, vx/F, vy/F, vz/F, vt/F");
    mcEventHeader.parse("npart/S, atarget/S, ztarget/S, ptarget/F, pbeam/F, btype/S, ebeam/F, targetid/S, processid/S, weight/F");
    userLund.parse("userVar/F");
    lundParticle.parse("index/B, lifetime/F, type/B, pid/I, parent/B, daughter/B, px/F, py/F, pz/F, energy/F, mass/F, vx/F, vy/F, vz/F");

    // flux
    fluxADCSchema.parse("sector/B, layer/B, component/S, order/B, ADC/I, amplitude/I, time/F, ped/S");


    // detectors
    alertAhdcADCSchema.parse("sector/B, layer/B, component/S, order/B, ADC/I, time/F, ped/S, integral/I, timestamp/L");
    alertAhdcTDCSchema.parse("sector/B, layer/B, component/S, order/B, TDC/I, ped/S");
    alertAhdcWF10Schema.parse("sector/B, layer/B, component/S, order/B, timestamp/F, s1/S, s2/S, s3/S, s4/S, s5/S, s6/S, s7/S, s8/S, s9/S, s10/S");
    alertAtofADCSchema.parse("sector/B, layer/B, component/S, order/B, ADC/I, time/F, ped/S");

    bandADCSchema.parse("sector/B, layer/B, component/S, order/B, ADC/I, amplitude/I, time/F, ped/S");
    bandTDCSchema.parse("sector/B, layer/B, component/S, order/B, TDC/I");

    bmtADCSchema.parse("sector/B, layer/B, component/S, order/B, ADC/I, time/F, ped/S, integral/I, timestamp/L");
    fmtADCSchema.parse("sector/B, layer/B, component/S, order/B, ADC/I, time/F, ped/S, integral/I, timestamp/L");
    bstADCSchema.parse("sector/B, layer/B, component/S, order/B, ADC/I, time/F, ped/S, timestamp/L");
    cndADCSchema.parse("sector/B, layer/B, component/S, order/B, ADC/I, time/F, ped/S");
    cndTDCSchema.parse("sector/B, layer/B, component/S, order/B, TDC/I");
    ctofADCSchema.parse("sector/B, layer/B, component/S, order/B, ADC/I, time/F, ped/S");
    ctofTDCSchema.parse("sector/B, layer/B, component/S, order/B, TDC/I");

    dcTDCSchema.parse("sector/B, layer/B, component/S, order/B, TDC/I");

    // need to add pcal to this
    ecalADCSchema.parse("sector/B, layer/B, component/S, order/B, ADC/I, time/F, ped/S");
    ecalTDCSchema.parse("sector/B, layer/B, component/S, order/B, TDC/I");
    ftcalADCSchema.parse("sector/B, layer/B, component/S, order/B, ADC/I, time/F, ped/S");
    fthodoADCSchema.parse("sector/B, layer/B, component/S, order/B, ADC/I, time/F, ped/S");
    ftrkTDCSchema.parse("sector/B, layer/B, component/S, order/B, ADC/I, time/F, ped/S, integral/I, timestamp/L");
    ftofADCSchema.parse("sector/B, layer/B, component/S, order/B, ADC/I, time/F, ped/S");
    ftofTDCSchema.parse("sector/B, layer/B, component/S, order/B, TDC/I");
    htccADCSchema.parse("sector/B, layer/B, component/S, order/B, ADC/I, time/F, ped/S");
    htccTDCSchema.parse("sector/B, layer/B, component/S, order/B, TDC/I");
    ltccADCSchema.parse("sector/B, layer/B, component/S, order/B, ADC/I, time/F, ped/S");
    ltccTDCSchema.parse("sector/B, layer/B, component/S, order/B, TDC/I");

    rfADCSchema.parse("sector/B, layer/B, component/S, order/B, ADC/I, time/F, ped/S");
    rfTDCSchema.parse("sector/B, layer/B, component/S, order/B, TDC/I");

    richTDCSchema.parse("sector/B, layer/B, component/S, order/B, TDC/I");
    rtpcADCSchema.parse("sector/B, layer/B, component/S, order/B, ADC/I, time/F, ped/S");
    rtpcPOSSchema.parse("step/I, time/F, energy/F, posx/F, posy/F, posz/F, phi/F, tid/F");

    helADCSchema.parse("sector/B, layer/B, component/S, order/B, ADC/I, time/F, ped/S");
    helFLIPSchema.parse("run/I, event/I, timestamp/L, helicity/B, helicityRaw/B, pair/B, pattern/B, status/B");
    helONLINESchema.parse("helicity/B, helicityRaw/B");

    urwellADCSchema.parse("sector/B, layer/B, component/S, order/B, ADC/I, time/F, ped/S");

    rawADCSchema.parse("crate/B, slot/B, channel/S, order/B, ADC/I, time/F, ped/S");
    rawTDCSchema.parse("crate/B, slot/B, channel/S, order/B, TDC/I");
    rawSCALERSchema.parse("crate/B, slot/B, channel/S, helicity/B, quartet/B, value/L");
    rawVTPSchema.parse("crate/B, word/I");
    rawEPICSSchema.parse("json/B");

    emptySchema.parse("empty/B");

    schemasToLoad["RUN::config"] = runConfigSchema;
    schemasToLoad["RUN::rf"] = runRFSchema;
    schemasToLoad["MC::True"] = trueInfoSchema;

    // generators
    schemasToLoad["MC::Particle"] = geantParticle;
    schemasToLoad["MC::Event"] = mcEventHeader;
    schemasToLoad["MC::User"] = userLund;
    schemasToLoad["MC::Lund"] = lundParticle;

    // flux
    schemasToLoad["FLUX::adc"] = fluxADCSchema;

    // The names corresponds to the hit process routine names, capitalized
    schemasToLoad["AHDC::adc"] = alertAhdcADCSchema;
    schemasToLoad["AHDC::wf:10"] = alertAhdcWF10Schema;
    schemasToLoad["AHDC::tdc"] = alertAhdcTDCSchema;
    schemasToLoad["ATOF::adc"] = alertAtofADCSchema;
    schemasToLoad["BAND::adc"] = bandADCSchema;
    schemasToLoad["BAND::tdc"] = bandTDCSchema;
    schemasToLoad["BMT::adc"] = bmtADCSchema;
    schemasToLoad["BST::adc"] = bstADCSchema;
    schemasToLoad["CND::adc"] = cndADCSchema;
    schemasToLoad["CND::tdc"] = cndTDCSchema;
    schemasToLoad["CTOF::adc"] = ctofADCSchema;
    schemasToLoad["CTOF::tdc"] = ctofTDCSchema;
    schemasToLoad["DC::tdc"] = dcTDCSchema;
    schemasToLoad["ECAL::adc"] = ecalADCSchema;
    schemasToLoad["ECAL::tdc"] = ecalTDCSchema;
    schemasToLoad["FMT::adc"] = fmtADCSchema;
    schemasToLoad["FTCAL::adc"] = ftcalADCSchema;
    schemasToLoad["FTHODO::adc"] = fthodoADCSchema;
    schemasToLoad["FTTRK::adc"] = ftrkTDCSchema;
    schemasToLoad["FTOF::adc"] = ftofADCSchema;
    schemasToLoad["FTOF::tdc"] = ftofTDCSchema;
    schemasToLoad["HTCC::adc"] = htccADCSchema;
    schemasToLoad["HTCC::tdc"] = htccTDCSchema;
    schemasToLoad["LTCC::adc"] = ltccADCSchema;
    schemasToLoad["LTCC::tdc"] = ltccTDCSchema;
    schemasToLoad["RICH::tdc"] = richTDCSchema;
    schemasToLoad["RTPC::adc"] = rtpcADCSchema;
    schemasToLoad["RTPC::pos"] = rtpcPOSSchema;
    schemasToLoad["HEL::flip"] = helFLIPSchema;
    schemasToLoad["RASTER::adc"] = rasterADCSchema;
    schemasToLoad["URWELL::adc"] = urwellADCSchema;

    cout << " Done defining Hipo4 schemas." << endl;

}

#include <cstring>

// type: 0 = adc, 1 = tdc
hipo::schema HipoSchema::getSchema(string schemaName, int type) {

    string schemaType;
    if (type == 0) schemaType = "adc";
    if (type == 1) schemaType = "tdc";
    if (type == 2) schemaType = "wf:10";

    string toUpperS = schemaName;
    transform(toUpperS.begin(), toUpperS.end(), toUpperS.begin(), ::toupper);
    string thisSchema = toUpperS + "::" + schemaType;

    if (schemasToLoad.find(thisSchema) != schemasToLoad.end()) {
        return schemasToLoad[thisSchema];
    } else if (schemaName == "ft_cal") {
        return ftcalADCSchema;
    } else if (schemaName == "ft_hodo") {
        return fthodoADCSchema;
    } else if (schemaName == "ft_trk") {
        return ftrkTDCSchema;
    } else {
        if (non_registered_detectors(schemaName, type)) {
            cout << " SCHEMA " << schemaName << " " << " not found for type " << type << " = " << schemaType << endl;
        }
        return emptySchema;
    }
}


bool HipoSchema::non_registered_detectors(string schemaName, int type) {


    if (type == 0) {  // non adc detectors:
        if (schemaName == "dc" || schemaName == "rich") {
            return false;
        }
    } else if (type == 1) { // non tdc detectors
        if (schemaName == "bmt" || schemaName == "fmt" || schemaName == "rtpc" || schemaName == "bst" || schemaName == "atof" || schemaName == "urwell" || schemaName == "flux") {
            return false;
        }
    } else if (type == 2) { // non wf:10 detectors
        if (schemaName == "atof" || schemaName == "band" || schemaName == "bmt" || schemaName == "fmt" || schemaName == "ftm"
            || schemaName == "dc" || schemaName == "bst" || schemaName == "cnd" || schemaName == "ctof" || schemaName == "ecal"
            || schemaName == "ftof" || schemaName == "ft_cal" || schemaName == "ft_hodo" || schemaName == "ft_trk"
            || schemaName == "htcc" || schemaName == "ltcc" || schemaName == "rich" || schemaName == "rtpc" || schemaName == "urwell" || schemaName == "flux") {
            return false;
        }
    }

    return true;
}


void outputContainer::initializeHipo(bool openFile) {

    if (!openFile) {
        hipoSchema = new HipoSchema();

        cout << " Initializing hipoSchema" << endl;
        hipoWriter = new hipo::writer();

        // Open a writer and register schemas with the writer.
        // The schemas have to be added to the writer before openning
        // the file, since they are written into the header of the file.
        for (auto &schema: hipoSchema->schemasToLoad) {
            hipoWriter->getDictionary().addSchema(schema.second);
            // hipoWriter->addUserConfig("gemc","{\"version\": \"4.4.2\", \"beam\": \"e-,10.6 GeV\"}");
        }
    } else {
        hipoWriter->open(outFile.c_str());
    }
}
