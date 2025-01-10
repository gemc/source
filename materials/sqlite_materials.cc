// Qt headers
#include <QtSql>

// gemc headers
#include "material_factory.h"
#include "sqlite_materials.h"
#include "string_utilities.h"
#include "gemcUtils.h"

// mlibrary
#include "gstring.h"

using namespace gstring;

// G4 headers
#include "G4Element.hh"
#include "G4NistManager.hh"
#include "G4OpBoundaryProcess.hh"


map<string, G4Material *> sqlite_materials::initMaterials(runConditions rc, goptions opts) {

    string hd_msg = opts.optMap["LOG_MSG"].args + " SQLITE Materials Factory: >> ";
    double verbosity = opts.optMap["MATERIAL_VERBOSITY"].arg;

    double runno_arg = opts.optMap["RUNNO"].arg;

    map <string, material> mymats;                        // material map

    // first check if there's at least one detector with SQLITE factory
    if (!check_if_factory_is_needed(rc.detectorConditionsMap, "SQLITE"))
        return materialsFromMap(mymats);

    // connection to the DB
    QSqlDatabase db = openGdb(opts);

    // Looping over detectorConditionsMap for detector names
    // To each detector is associated a material and (optional) opt properties
    for (map<string, detectorCondition>::iterator it = rc.detectorConditionsMap.begin(); it != rc.detectorConditionsMap.end(); it++) {
        // building materials belonging to detectors that are tagged with SQLITE factory
        if (it->second.get_factory() != "SQLITE")
            continue;

        if (verbosity)
            cout << hd_msg << " Initializing " << it->second.get_factory() << " for detector " << it->first << endl;

        string dname = it->first;
        string variation = get_variation(it->second.get_variation());

        // geometry was empty. TODO: this mechanism is a bit clunky . Solution: add a detector with existence = 0 in the geometry table
        if (variation == "empty") variation = "default";

        int run = it->second.get_run_number();
        if (runno_arg != -1) run = runno_arg; // if RUNNO is set (different from -1), use it
        int run_number = get_sql_run_number(db, dname, variation, run, "materials");

        string dbexecute = "select name, description, density, ncomponents, components, photonEnergy, indexOfRefraction, ";
        dbexecute += "absorptionLength, reflectivity, efficiency, fastcomponent, slowcomponent, ";
        dbexecute += "scintillationyield, resolutionscale, fasttimeconstant, slowtimeconstant, yieldratio, rayleigh, birkconstant, "; // rayleigh and birk constant were not here?
        dbexecute += "mie, mieforward, miebackward, mieratio from materials";
        dbexecute += " where variation ='" + variation + "'";
        dbexecute += " and run = " + stringify(run_number);
        dbexecute += " and system = '" + dname + "'";

        // executing query - will exit if not successful.
        QSqlQuery q;
        if (!q.exec(dbexecute.c_str())) {
            cout << hd_msg << "  Failed to execute SQLITE query " << dbexecute << ". This is a fatal error. Exiting." << endl;
            qDebug() << q.lastError();
            exit(1);
        }
        // Warning if nothing is found
        if (q.size() == 0 && verbosity) {
            cout << "  ** WARNING: material for system \"" << dname << "\" not found with variation \"" << variation << endl << endl;
        }

        while (q.next()) {
            material thisMat(trimSpacesFromString(qv_tostring( q.value(0))));         // name
            thisMat.desc = qv_tostring(q.value(1));                                   // description
            thisMat.density = q.value(2).toDouble();                                  // density
            thisMat.ncomponents = q.value(3).toInt();                                 // number of components
            thisMat.componentsFromString(qv_tostring(q.value(4)));                    // component + quantity list
            thisMat.opticalsFromString(qv_tostring(q.value(5)), "photonEnergy");
            thisMat.opticalsFromString(qv_tostring(q.value(6)), "indexOfRefraction");
            thisMat.opticalsFromString(qv_tostring(q.value(7)), "absorptionLength");
            thisMat.opticalsFromString(qv_tostring(q.value(8)), "reflectivity");
            thisMat.opticalsFromString(qv_tostring(q.value(9)), "efficiency");

            // scintillation
            thisMat.opticalsFromString(qv_tostring(q.value(10)), "fastcomponent");
            thisMat.opticalsFromString(qv_tostring(q.value(11)), "slowcomponent");
            thisMat.scintillationyield = q.value(12).toDouble();
            thisMat.resolutionscale = q.value(13).toDouble();
            thisMat.fasttimeconstant = q.value(14).toDouble();
            thisMat.slowtimeconstant = q.value(15).toDouble();
            thisMat.yieldratio = q.value(16).toDouble();
            thisMat.opticalsFromString(qv_tostring(q.value(17)), "rayleigh");

            // Birk Constant
            thisMat.birkConstant = q.value(18).toDouble();

            // Mie scattering
            thisMat.opticalsFromString(qv_tostring(q.value(19)), "mie");
            thisMat.mieforward = q.value(20).toDouble();
            thisMat.miebackward = q.value(21).toDouble();
            thisMat.mieratio = q.value(22).toDouble();

            mymats[thisMat.name] = thisMat;

        }
    }

    cout << endl;

    map < string, G4Material * > returnMap = materialsFromMap(mymats);
    if (verbosity > 0) printMaterials(returnMap);

    return returnMap;
}
