// gemc
#include "sqlite_cad_factory.h"
#include "gemcUtils.h"

// c++
#include <dirent.h>

map <string, detector> sqlitecad_det_factory::loadDetectors() {

    string hd_msg = " >> SQLITE CAD Factory: >> ";
    double verbosity = gemcOpt.optMap["GEO_VERBOSITY"].arg;
    string catch_v = gemcOpt.optMap["CATCH"].args;
    double runno_arg = gemcOpt.optMap["RUNNO"].arg;

    map <string, detector> dets;

    // first check if there's at least one detector with SQLITE factory
    if (!check_if_factory_is_needed(RC.detectorConditionsMap, factoryType)) { return dets; }

    // connection to the DB
    QSqlDatabase db = openGdb(gemcOpt);

    // connection is ok. Loading tables now
    // building detectors that are tagged with MYSQL factory
    for (map<string, detectorCondition>::iterator it = RC.detectorConditionsMap.begin(); it != RC.detectorConditionsMap.end(); it++) {

        if (it->second.get_factory() != factoryType) { continue; }

        string dname = it->first;
        string variation = get_variation(it->second.get_variation());

        int run = it->second.get_run_number();
        if (runno_arg != -1) run = runno_arg; // if RUNNO is set (different from -1), use it
        int run_number = get_sql_run_number(db, dname, variation, run, "cad");


        // if run number is -1, the detector is not in the DB. exit with error
        if (run_number == -1) {
            cout << hd_msg << " !!! Detector \"" << dname << "\" not found in the SQLITE CAD database for run number <" << run_number << ">. Exiting." << endl;
            exit(1);
        }

        if (verbosity) {
            cout << hd_msg << " Importing SQLITE CAD Detector: " << dname << " with <" << factoryType
                 << "> factory, variation \"" << variation << "\", run number requested: " << run << ",  sqlite run number used: " << run_number << endl;
        }


        string dbexecute = "select name, subdir, sensitivity, hit_type, identifiers, " ;
        dbexecute += "visible, style, position, rotation, mfield, mother, material, color from cad" ;
        dbexecute += " where variation ='" + variation + "'";
        dbexecute += " and run = " + stringify(run_number) ;
        dbexecute += " and system = '" + dname  + "'";

        // executing query - will exit if not successfull.
        QSqlQuery q;
        if (!q.exec(dbexecute.c_str())) {
            cout << hd_msg << " !!! Failed to execute CAD SQLITE query " << dbexecute << ". This is a fatal error. Exiting." << endl;
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

            string name = q.value(0).toString().toStdString();
            string cad_subdir = q.value(1).toString().toStdString();
            string sensitivity = q.value(2).toString().toStdString();
            string hit_type = q.value(3).toString().toStdString();
            string identifiers = q.value(4).toString().toStdString();
            string visible = q.value(5).toString().toStdString();
            string style = q.value(6).toString().toStdString();
            string position = q.value(7).toString().toStdString();
            string rotation = q.value(8).toString().toStdString();
            string mfield = q.value(9).toString().toStdString();
            string mother = q.value(10).toString().toStdString();
            string material = q.value(11).toString().toStdString();
            string color = q.value(12).toString().toStdString();

            string filename = checkFormat(cad_subdir + "/" + name);
            cout << " > Parsing 2 CAD volume from " << filename << endl;

            gt.add_data(q.value(0));                        // 1 name
            gt.add_data(mother);                            // 2 mother volume
            gt.add_data((string) " cadImported");                    // 3 description
            gt.add_data(position);                          // 4 position
            gt.add_data(rotation);                          // 5 rotation
            gt.add_data(color);                             // 6 color
            gt.add_data((string) "cadImport");              // 7 type
            gt.add_data((string) "0");                      // 8 dimensions
            gt.add_data(material);                          // 9 material is aluminum by defaul
            gt.add_data(mfield);                            // 10 magnetic field
            gt.add_data((string) "0");                      // 11 copy number
            gt.add_data((string) "0");                      // 12 pmany
            gt.add_data((string) "1");                      // 13 activation flag
            gt.add_data(visible);                           // 14 visibility
            gt.add_data(style);                             // 15 style
            gt.add_data(sensitivity);                       // 16 sensitivity
            gt.add_data(hit_type);                          // 17 hit_type
            gt.add_data(identifiers);                       // 18 identifiers
            gt.add_data(dname);                             // 19 system is dname (can be a path)
            gt.add_data((string) "SQLITECAD");              // 20 factory
            gt.add_data(filename);                          // 21 variation is the full filename
            gt.add_data(stringify(run_number));             // 22 run number

            // big warning if detector already exist
            // detector is NOT loaded if already existing
            if (dets.find(gt.data[0]) != dets.end()) {
                cout << endl << " *** WARNING! A detector >" << gt.data[0] << " exists already. Keeping original, not loading this instance. " << endl << endl;
            } else {
                dets[gt.data[0]] = get_detector(gt, gemcOpt, RC);
            }

            if (verbosity > 2) { cout << get_detector(gt, gemcOpt, RC); }
        }


//        DIR *thisDir = opendir(cad_subdir.c_str());
//        if (thisDir == nullptr) {
//            if (getenv("GEMC_DATA_DIR") != nullptr) {
//                string maybeHere = (string) getenv("GEMC_DATA_DIR") + "/" + dname + "/" + cad_subdir;
//                thisDir = opendir(maybeHere.c_str());
//            }
//        }
//
//        if (thisDir != nullptr) {
//            struct dirent *thisDirent = readdir(thisDir);
//            while (thisDirent != nullptr) {
//
//                // removing 4 char from filename. Extension must be 3 letters long + period
//                string thisName = thisDirent->d_name;
//
//                // removing extension including .
//                if (thisName.size() > 4) {
//                    string thisRootFileName(thisName.begin(), thisName.end() - 4);
//
//                    string thisFileName = checkFormat(dname + thisRootFileName);
//
//                    if (thisFileName != "na") {
//                        cadFiles.push_back(thisFileName);
//                    }
//                }
//            }
//            closedir(thisDir);
//
//        } else {
//            cout << " !! Error: CAD directory <" << cad_subdir
//                 << "> cannot be read. Did you set GEMC_DATA_DIR to point to the location containing the experiments folder? Exiting." << endl;
//            exit(1);
//        }

        // this is a directory, use name directly

//        if (cadFiles.size() == 0) {
//            cout << " !! Error: SQLITECAD system " << it->first << " subdir " << cad_subdir
//                 << " is not a directory containing cad files. Exiting." << endl;
//            exit(1);
//        }


        // looping over all cad files for this system
//        for (const auto &cf: cadFiles) {
//            // there is only one solid / cad file
//            // using gtable to create it






//
//
//
//                    string volumeName = e.attribute("name").toStdString();
//                    string mother = e.attribute("mother").toStdString();
//                    string color = e.attribute("color").toStdString();
//                    string material = e.attribute("material").toStdString();
//                    string sensitivity = e.attribute("sensitivity").toStdString();
//                    string hitType = e.attribute("hitType").toStdString();
//                    string identifiers = e.attribute("identifiers").toStdString();
//                    string visible = e.attribute("visible").toStdString();
//                    string style = e.attribute("style").toStdString();
//                    string position = e.attribute("position").toStdString();
//                    string rotation = e.attribute("rotation").toStdString();
//                    string mfield = e.attribute("mfield").toStdString();
//                    detector_in_gxml.push_back(volumeName);
//                    // assigning attributes to volume
//                    if (dets.find(volumeName) != dets.end()) {
//
//                        if (verbosity > 3)
//                            cout << " Modifying attributes for volume: " << volumeName;
//
//                        if (color != "") {
//                            if (verbosity > 3)
//                                cout << " color: " << color;
//
//                            G4Colour thisCol = gcol(color);
//                            dets[volumeName].VAtts = G4VisAttributes(thisCol);
//                        }
//
//                        visible == "no" ? dets[volumeName].VAtts.SetVisibility(false) : dets[volumeName].VAtts.SetVisibility(true);
//                        style == "wireframe" ? dets[volumeName].VAtts.SetForceWireframe(true) : dets[volumeName].VAtts.SetForceSolid(true);
//
//                        if (sensitivity != "") {
//                            if (verbosity > 3)
//                                cout << " sensitivity: " << sensitivity;
//
//                            // setting system as sensitivity, so the hit definitions can be loaded
//                            // this should be modified later
//                            dets[volumeName].system = sensitivity;
//                            dets[volumeName].sensitivity = sensitivity;
//
//                            // identifier must be defined
//                            if (identifiers != "") {
//                                if (verbosity > 3)
//                                    cout << " identifiers: " << identifiers;
//
//                                dets[volumeName].identity = get_identifiers(identifiers);
//
//                            } else {
//                                cout << " !! Error: volume " << volumeName << " has sensitivity but not identifier. " << endl;
//                            }
//                        }
//
//                        if (hitType != "") {
//                            if (verbosity > 3)
//                                cout << " hitType: " << hitType;
//                            dets[volumeName].hitType = hitType;
//                        }
//
//                        if (material != "") {
//                            if (verbosity > 3)
//                                cout << " material: " << material;
//                            dets[volumeName].material = material;
//                        }
//
//                        if (position != "") {
//                            if (verbosity > 3)
//                                cout << " position: " << position;
//                            dets[volumeName].pos = calc_position(position);
//                        }
//
//                        if (rotation != "") {
//                            if (verbosity > 3)
//                                cout << " rotation: " << rotation;
//                            dets[volumeName].rot = calc_rotation(rotation, volumeName);
//                        }
//
//                        if (mother != "") {
//                            if (verbosity > 3)
//                                cout << " mother: " << mother;
//                            dets[volumeName].mother = mother;
//                        }
//                        if (mfield != "") {
//                            if (verbosity > 3)
//                                cout << " mfield: " << mfield;
//                            dets[volumeName].magfield = mfield;
//                        }
//
//                        if (verbosity > 3)
//                            cout << endl;
//                    }
//                }
//            }
//        }
//    }


        for (const auto &dd: dets) {
            if (verbosity > 3)
                cout << dd.second;
        }

        // remove all detectors in dets that are not in the sql DB
//    for (const auto &dd: dets) {
//        if (find(detector_in_gxml.begin(), detector_in_gxml.end(), dd.first) == detector_in_gxml.end()) {
//            cout << " >>  Detector " << dd.first << " not found in gxml file. It will be removed from GEMC." << endl;
//            dets[dd.first].exist = 0;
//        }
//    }

    }
    return dets;
}

string sqlitecad_det_factory::checkFormat(string file) {
    // trying STL
    string filename = file + ".stl";
    ifstream fstl(filename.c_str());
    if (fstl.good()) {
        fstl.close();
        return filename;
        // trying .PLY
    } else {
        filename = file + ".ply";
        ifstream fply(filename.c_str());
        if (fply.good()) {
            fply.close();
            return filename;
            // trying OBJ
        } else {
            filename = file + ".obj";
            ifstream fobj(filename.c_str());
            if (fobj.good()) {
                fobj.close();
                return filename;
            } else {
                return "na";
            }
        }
    }

    return "na";
}
