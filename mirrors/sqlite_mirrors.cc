// Qt headers
#include <QtSql>

// gemc headers
#include "mirrors_factory.h"
#include "sqlite_mirrors.h"
#include "string_utilities.h"
#include "gemcUtils.h"

// mlibrary
#include "gstring.h"

using namespace gstring;

map<string, mirror *> sqlite_mirrors::initMirrors(runConditions rc, goptions opts) {
    string hd_msg = opts.optMap["LOG_MSG"].args + "  > SQLITE mirrors Factory: >> ";
    double verbosity = opts.optMap["MIRROR_VERBOSITY"].arg;
    double runno_arg = opts.optMap["RUNNO"].arg;

    map < string, mirror * > mymirs;  // mirror map

    // first check if there's at least one detector with SQLITE factory
    if (!check_if_factory_is_needed(rc.detectorConditionsMap, "SQLITE")) { return mymirs; }

    // connection to the DB
    QSqlDatabase db = openGdb(opts);

    // Looping over detectorConditionsMap for detector names
    // To each detector is associated a mirror table and an opt properties table
    for (map<string, detectorCondition>::iterator it = rc.detectorConditionsMap.begin(); it != rc.detectorConditionsMap.end(); it++) {
        // building mirrors belonging to detectors that are tagged with SQLITE factory
        if (it->second.get_factory() != "SQLITE") { continue; }

        string dname = it->first;
        string variation = get_variation(it->second.get_variation());
        int run = it->second.get_run_number();
        if (runno_arg != -1) run = runno_arg; // if RUNNO is set (different from -1), use it
        int run_number = get_sql_run_number(db, dname, variation, run, "mirrors");

        // if run number is -1, the detector is not in the DB. Return empty map
        if (run_number == -1) {
            return mymirs;
        }

        if (verbosity) {
            cout << hd_msg << " Importing SQLITE mirrors for detector " << dname << " with variation \""
                 << variation << "\", run number requested: " << run << ",  sqlite run number used: " << run_number << endl;
        }

        string dbexecute = "select name, description, type, finish, model, border, ";
        dbexecute += " mat_opt_props, photon_energy, refraction_index, reflectivity, efficiency, ";
        dbexecute += " specular_lobe, specular_spike, backscatter, sigma_alpha from mirrors";
        dbexecute += " where variation ='" + variation + "'";
        dbexecute += " and run = " + stringify(run_number);
        dbexecute += " and system = '" + dname + "'";

        // executing query - will exit if not successfull.
        QSqlQuery q;
        if (!q.exec(dbexecute.c_str())) {
            cout << hd_msg << "  Failed to execute SQLITE query " << dbexecute << ". This is a fatal error. Exiting." << endl;
            qDebug() << q.lastError();
            exit(1);
        }
        // Warning if nothing is found
        if (q.size() == 0 && verbosity) {
            cout << "  ** WARNING: mirror for system \"" << dname << "\" not found with variation \"" << variation << endl << endl;
        }

        while (q.next()) {
            mirror *thisMir = new mirror(trimSpacesFromString(qv_tostring(q.value(0))));         // name
            thisMir->desc = qv_tostring(q.value(1));                                             // description
            thisMir->type = qv_tostring(q.value(2));                                             // type
            thisMir->finish = qv_tostring(q.value(3));                                           // finish
            thisMir->model = qv_tostring(q.value(4));                                            // model
            thisMir->border = qv_tostring(q.value(5));                                           // border
            thisMir->maptOptProps = qv_tostring(q.value(6));                                     // maptOptProps
            thisMir->opticalsFromString(qv_tostring(q.value(7)), "photonEnergy");
            thisMir->opticalsFromString(qv_tostring(q.value(8)), "indexOfRefraction");
            thisMir->opticalsFromString(qv_tostring(q.value(9)), "reflectivity");
            thisMir->opticalsFromString(qv_tostring(q.value(10)), "efficiency");
            thisMir->opticalsFromString(qv_tostring(q.value(11)), "specularlobe");
            thisMir->opticalsFromString(qv_tostring(q.value(12)), "specularspike");
            thisMir->opticalsFromString(qv_tostring(q.value(13)), "backscatter");
            thisMir->sigmaAlpha = q.value(14).toDouble();

            mymirs[thisMir->name] = thisMir;
        }
    }

    cout << endl;

    if (verbosity > 0) printMirrors(mymirs);

    return mymirs;
}
