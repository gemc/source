// gemc
#include "cad_det_factory.h"
#include "gemcUtils.h"

// c++
#include <dirent.h>

map <string, detector> cad_det_factory::loadDetectors() {

    string hd_msg = " >> CAD Factory: >> ";
    double verbosity = gemcOpt.optMap["GEO_VERBOSITY"].arg;
    string catch_v = gemcOpt.optMap["CATCH"].args;

    map <string, detector> dets;

    // first check if there's at least one detector with CAD factory
    if (!check_if_factory_is_needed(RC.detectorConditionsMap, factoryType))
        return dets;

    // looking for possible gxml files
    vector <string> possibleGXML;

    // there is at least one build with this factory.
    // building all detectors that are tagged with CAD factory
    for (map<string, detectorCondition>::iterator it = RC.detectorConditionsMap.begin(); it != RC.detectorConditionsMap.end(); it++) {
        if (it->second.get_factory() != factoryType)
            continue;

        string dname = it->first;
        string variation = it->second.get_variation();

        vector <string> cadFiles;

        // checking if "/" is the end of the cad name
        if (strcmp(&dname.back(), "/") == 0) {
            string dir(dname.begin(), dname.end() - 1);
            DIR *thisDir = opendir(dir.c_str());

            if (thisDir == nullptr) {
                if (getenv("GEMC_DATA_DIR") != nullptr) {
                    string maybeHere = (string) getenv("GEMC_DATA_DIR") + "/" + dir;
                    thisDir = opendir(maybeHere.c_str());
                }
            }

            if (thisDir != nullptr) {
                possibleGXML.push_back(dname + "cad_" + variation + ".gxml");
                struct dirent *thisDirent = readdir(thisDir);
                while (thisDirent != nullptr) {

                    // removing 4 char from filename. Extension must be 3 letters long + period
                    string thisName = thisDirent->d_name;

                    // removing extension including .
                    if (thisName.size() > 4) {
                        string thisRootFileName(thisName.begin(), thisName.end() - 4);

                        string thisFileName = checkFormat(dname + thisRootFileName);

                        if (thisFileName != "na") {
                            cadFiles.push_back(thisFileName);
                        } else {
                            // trying GEMC_DATA_DIR
                            if (getenv("GEMC_DATA_DIR") != nullptr) {
                                string envLoc = (string) getenv("GEMC_DATA_DIR") + "/";
                                thisFileName = checkFormat(envLoc + dname + thisRootFileName);
                                if (thisFileName != "na") {
                                    cadFiles.push_back(thisFileName);
                                }
                            }
                        }
                    }

                    thisDirent = readdir(thisDir);
                }

            } else {
                cout << " !! Error: directory " << dir << " cannot be read. Did you set GEMC_DATA_DIR to point to the location containing the experiments folder? Exiting." << endl;
                exit(1);
            }
            closedir(thisDir);

            // this is a directory, use name directly
        } else {
            string thisFileName = checkFormat(dname);
            if (thisFileName != "na")
                cadFiles.push_back(thisFileName);
        }

        if (cadFiles.size() == 0) {
            cout << " !! Error: cad system " << it->first << " is not a cad file or a directory containing cad files. Exiting." << endl;
            exit(1);
        }

        if (verbosity)
            cout << hd_msg << " Importing Cad Detector from: <" << dname << "> with " << factoryType << " factory. " << endl;

        // looping over all cad files for this system
        for (const auto &cf: cadFiles) {
            // there is only one solid / cad file
            // using gtable to create it
            gtable gt;

            // cf is at least 5 chars because it has 3 letters extension
            string detPN(cf.begin(), cf.end() - 4);

            // removing path
            string detN = detPN.substr(detPN.find_last_of("/") + 1);

            gt.add_data(detN);                            // 1 name
            gt.add_data((string) "root");                // 2 mother volume
            gt.add_data(detN + (string) " cadImported");  // 3 description
            gt.add_data((string) "0*cm 0*cm 0*cm");      // 4 position
            gt.add_data((string) "0*deg 0*deg 0*deg");   // 5 rotation
            gt.add_data((string) "2222aa");              // 6 color
            gt.add_data((string) "cadImport");           // 7 type
            gt.add_data((string) "0");                   // 8 dimensions
            gt.add_data((string) "G4_Al");               // 9 material is aluminum by defaul
            gt.add_data((string) "no");                  // 10 magnetic field
            gt.add_data((string) "0");                   // 11 copy number
            gt.add_data((string) "0");                   // 12 pmany
            gt.add_data((string) "1");                   // 13 activation flag
            gt.add_data((string) "1");                   // 14 visibility
            gt.add_data((string) "1");                   // 15 style
            gt.add_data((string) "no");                  // 16 sensitivity
            gt.add_data((string) "no");                  // 17 hit_type
            gt.add_data((string) "");                    // 18 identifiers
            gt.add_data(dname);                          // 19 system is dname (can be a path)
            gt.add_data((string) "CAD");                 // 20 factory
            gt.add_data(cf);                              // 21 variation is the full filename
            gt.add_data((string) "1");                   // 22 run number

            dets[gt.data[0]] = get_detector(gt, gemcOpt, RC);

        }
    }

    bool no_gxml_found = true;
    vector <string> detector_in_gxml;
    for (unsigned g = 0; g < possibleGXML.size(); g++) {

        QFile *gxml = new QFile(possibleGXML[g].c_str());

        if (!gxml->exists()) {
            // trying GEMC_DATA_DIR
            string envLoc = (string) getenv("GEMC_DATA_DIR") + "/" + possibleGXML[g];
            gxml = new QFile(envLoc.c_str());
        }

        // checking that file exists
        if (!gxml->exists()) {
            cout << hd_msg << " " << possibleGXML[g] << " not found. Cad volumes will be imported as original." << endl;
        } else {
            no_gxml_found = false;
            QDomDocument domDocument;
            // opening gcard and filling domDocument
            if (!domDocument.setContent(gxml)) {
                cout << " >>  xml format for file <" << possibleGXML[g] << "> is wrong - check XML syntax. Exiting." << endl;
                exit(1);
            }
            gxml->close();

            QDomNodeList volumes = domDocument.firstChildElement().elementsByTagName("volume");
            for (int i = 0; i < volumes.count(); i++) {
                QDomNode elm = volumes.at(i);
                if (elm.isElement()) {
                    QDomElement e = elm.toElement();
                    string volumeName = e.attribute("name").toStdString();
                    string mother = e.attribute("mother").toStdString();
                    string color = e.attribute("color").toStdString();
                    string material = e.attribute("material").toStdString();
                    string sensitivity = e.attribute("sensitivity").toStdString();
                    string hitType = e.attribute("hitType").toStdString();
                    string identifiers = e.attribute("identifiers").toStdString();
                    string visible = e.attribute("visible").toStdString();
                    string style = e.attribute("style").toStdString();
                    string position = e.attribute("position").toStdString();
                    string rotation = e.attribute("rotation").toStdString();
                    string mfield = e.attribute("mfield").toStdString();
                    detector_in_gxml.push_back(volumeName);
                    // assigning attributes to volume
                    if (dets.find(volumeName) != dets.end()) {

                        if (verbosity > 3)
                            cout << " Modifying attributes for volume: " << volumeName;

                        if (color != "") {
                            if (verbosity > 3)
                                cout << " color: " << color;

                            G4Colour thisCol = gcol(color);
                            dets[volumeName].VAtts = G4VisAttributes(thisCol);
                        }

                        visible == "no" ? dets[volumeName].VAtts.SetVisibility(false) : dets[volumeName].VAtts.SetVisibility(true);
                        style == "wireframe" ? dets[volumeName].VAtts.SetForceWireframe(true) : dets[volumeName].VAtts.SetForceSolid(true);

                        if (sensitivity != "") {
                            if (verbosity > 3)
                                cout << " sensitivity: " << sensitivity;

                            // setting system as sensitivity, so the hit definitions can be loaded
                            // this should be modified later
                            dets[volumeName].system = sensitivity;
                            dets[volumeName].sensitivity = sensitivity;

                            // identifier must be defined
                            if (identifiers != "") {
                                if (verbosity > 3)
                                    cout << " identifiers: " << identifiers;

                                dets[volumeName].identity = get_identifiers(identifiers);

                            } else {
                                cout << " !! Error: volume " << volumeName << " has sensitivity but not identifier. " << endl;
                            }
                        }

                        if (hitType != "") {
                            if (verbosity > 3)
                                cout << " hitType: " << hitType;
                            dets[volumeName].hitType = hitType;
                        }

                        if (material != "") {
                            if (verbosity > 3)
                                cout << " material: " << material;
                            dets[volumeName].material = material;
                        }

                        if (position != "") {
                            if (verbosity > 3)
                                cout << " position: " << position;
                            dets[volumeName].pos = calc_position(position);
                        }

                        if (rotation != "") {
                            if (verbosity > 3)
                                cout << " rotation: " << rotation;
                            dets[volumeName].rot = calc_rotation(rotation, volumeName);
                        }

                        if (mother != "") {
                            if (verbosity > 3)
                                cout << " mother: " << mother;
                            dets[volumeName].mother = mother;
                        }
                        if (mfield != "") {
                            if (verbosity > 3)
                                cout << " mfield: " << mfield;
                            dets[volumeName].magfield = mfield;
                        }

                        if (verbosity > 3)
                            cout << endl;
                    }
                }
            }
        }
    }

    // if no_gxml_found, exit with error
    if (no_gxml_found) {
        cerr << " Error: no GXML file found among all choices: " << endl;
        for (auto &dd : possibleGXML) {
            cerr << "  - " << dd << endl;
        }
        exit (1);
    }

    for (const auto &dd: dets) {
        if (verbosity > 3)
            cout << dd.second;
    }

    // remove all detectors in dets that are not in the gxml file
    for (const auto &dd: dets) {
        if (find(detector_in_gxml.begin(), detector_in_gxml.end(), dd.first) == detector_in_gxml.end()) {
            cout << " >>  Detector " << dd.first << " not found in gxml file. It will be removed from GEMC." << endl;
            dets[dd.first].exist = 0;
        }
    }

    return dets;
}

string cad_det_factory::checkFormat(string file) {
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
