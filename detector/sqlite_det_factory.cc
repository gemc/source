// Qt headers
#include <QtSql>

// gemc headers
#include "sqlite_det_factory.h"
#include "gemcUtils.h"

map <string, detector> sqlite_det_factory::loadDetectors() {

    string hd_msg = gemcOpt.optMap["LOG_MSG"].args + "  > SQLITE Detector Factory: >> ";
    double verbosity = gemcOpt.optMap["GEO_VERBOSITY"].arg;

    double runno_arg = gemcOpt.optMap["RUNNO"].arg;
    string catch_v = gemcOpt.optMap["CATCH"].args;

    map <string, detector> dets;

    // first check if there's at least one detector with SQLITE factory
    if (!check_if_factory_is_needed(RC.detectorConditionsMap, factoryType)) { return dets; }

    // connection to the DB
    QSqlDatabase db = openGdb(gemcOpt);

    // connection is ok. Loading tables now
    // building detectors that are tagged with MYSQL factory
    for (map<string, detectorCondition>::iterator it = RC.detectorConditionsMap.begin(); it != RC.detectorConditionsMap.end(); it++) {

        if (it->second.get_factory() != factoryType)
            continue;

        string dname = it->first;
        string variation = get_variation(it->second.get_variation());
        int run = it->second.get_run_number();
        if (runno_arg != -1) run = runno_arg; // if RUNNO is set (different from -1), use it
        int run_number = get_sql_run_number(db, dname, variation, run, "geometry");

        // if run number is -1, the detector is not in the DB. exit with error
        if (run_number == -1) {
            cout << hd_msg << " !!! Detector \"" << dname << "\" not found in the SQLITE database for run number <" << run_number << ">. Exiting." << endl;
            exit(1);
        }

        if (verbosity) {
            cout << hd_msg << " Importing SQLITE Detector: " << dname << " with <" << factoryType
                 << "> factory, variation \"" << variation << "\", run number requested: " << run << ",  sqlite run number used: " << run_number << endl;
        }

        string dbexecute = "select name, mother, description, pos, rot, col, type, ";
        dbexecute += "dimensions, material, magfield, ncopy, pmany, exist, ";
        dbexecute += "visible, style, sensitivity, hitType, identity from geometry";
        dbexecute += " where variation ='" + variation + "'";
        dbexecute += " and run = " + stringify(run_number);
        dbexecute += " and system = '" + dname + "'";

        // executing query - will exit if not successfull.
        QSqlQuery q;
        if (!q.exec(dbexecute.c_str())) {
            cout << hd_msg << " !!! Failed to execute SQLITE query " << dbexecute << ". This is a fatal error. Exiting." << endl;
            qDebug() << q.lastError();
            exit(1);
        }
        // Warning if nothing is found
        if (q.size() == 0 && verbosity) {
            cout << "  ** WARNING: detector \"" << dname << "\" not found with variation \"" << variation << "\" for run number " << run << endl << endl;
        }

        // else loading parameters from DB
        while (q.next()) {
            gtable gt;

            // there should be 18 values in the MYSQL table
            for (int i = 0; i < 18; i++) {
                gt.add_data(q.value(i));
            }

            // adding additional info: system, factory, variation, run number
            gt.add_data(dname);
            gt.add_data((string) "SQLITE");
            gt.add_data(variation);
            gt.add_data(stringify(run_number));

            // big warning if detector already exist
            // detector is NOT loaded if already existing
            if (dets.find(gt.data[0]) != dets.end()) {
                cout << endl << " *** WARNING! A detector >" << gt.data[0] << " exists already. Keeping original, not loading this instance. " << endl << endl;
            } else {
                dets[gt.data[0]] = get_detector(gt, gemcOpt, RC);
            }

            if (verbosity > 3) { cout << get_detector(gt, gemcOpt, RC); }
        }
    }

    cout << endl;

    return dets;
}
