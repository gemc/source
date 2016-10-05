// gemc
#include "cad_det_factory.h"
#include "utils.h"

// c++
#include <dirent.h>

map<string, detector> cad_det_factory::loadDetectors()
{

	string hd_msg     = " >> CAD Factory: >> ";
	double verbosity  = gemcOpt.optMap["GEO_VERBOSITY"].arg;
	string catch_v    = gemcOpt.optMap["CATCH"].args;

	map<string, detector> dets;
	
  	// first check if there's at least one detector with CAD factory
	if(!check_if_factory_is_needed(RC.detectorConditionsMap, factoryType))
		return dets;


	// saving directories so we can look for cad.gxml on those locations
	vector<string> possibleDirs;

  	// there is at least one build with this factory.
	// building all detectors that are tagged with CAD factory
	for(map<string, detectorCondition>::iterator it=RC.detectorConditionsMap.begin(); it != RC.detectorConditionsMap.end() ; it++) {
		if(it->second.get_factory() != factoryType )
			continue;

		string dname = it->first;

		vector<string> cadFiles;

		// checking wildcards
		if(strcmp(&dname.back(), "/") == 0) {
			string dir(dname.begin(), dname.end()-1);
			DIR *thisDir = opendir(dir.c_str());

			if(thisDir != NULL) {
				possibleDirs.push_back(dir + "/");
				struct dirent *thisDirent = readdir(thisDir);
				while (thisDirent != NULL){

					// removing 4 char from filename. Extension must be 3 letters long + period
					string thisName = thisDirent->d_name;

					if(thisName.size()>4) {
						string thisRootFileName(thisName.begin(), thisName.end()-4);

						string thisFileName = checkFormat(dname + thisRootFileName);

						if(thisFileName != "na")
							cadFiles.push_back(thisFileName);
					}

					thisDirent = readdir(thisDir);
				}

			} else {
				cout << " !! Error: directory " << dir << " cannot be read. Exiting." << endl;
				exit(0);
			}
			closedir(thisDir);

		// no wildcard, use name directly
		} else {
			string thisFileName = checkFormat(dname);
			if(thisFileName != "na")
				cadFiles.push_back(thisFileName);
		}

		if(cadFiles.size() == 0) {
			cout << " !! Error: cad system " << it->first << " is not a cad file or a directory containing cad files. Exiting." << endl;
			exit(0);
		}

		if(verbosity)
			cout <<  hd_msg << " Importing Cad Detector from: <" <<  dname << "> with " << factoryType << " factory. "  << endl;

		// looping over all cad files for this system
		for(const auto &cf : cadFiles) {
			// there is only one solid / cad file
			// using gtable to create it
			gtable gt;

			// cf is at least 5 chars because it has 3 letters extension
			string detPN(cf.begin(), cf.end()-4);

			// removing path
			string detN = detPN.substr(detPN.find_last_of("/") + 1);

			gt.add_data(detN);                            // 1 name
			gt.add_data( (string) "root");                // 2 mother volume
			gt.add_data(detN + (string) " cadImported");  // 3 description
			gt.add_data( (string) "0*cm 0*cm 0*cm");      // 4 position
			gt.add_data( (string) "0*deg 0*deg 0*deg");   // 5 rotation
			gt.add_data( (string) "2222aa");              // 6 color
			gt.add_data( (string) "cadImport");           // 7 type
			gt.add_data( (string) "0");                   // 8 dimensions
			gt.add_data( (string) "G4_Al");               // 9 material is aluminum by defaul
			gt.add_data( (string) "no");                  // 10 magnetic field
			gt.add_data( (string) "0");                   // 11 copy number
			gt.add_data( (string) "0");                   // 12 pmany
			gt.add_data( (string) "1");                   // 13 activation flag
			gt.add_data( (string) "1");                   // 14 visibility
			gt.add_data( (string) "1");                   // 15 style
			gt.add_data( (string) "no");                  // 16 sensitivity
			gt.add_data( (string) "no");                  // 17 hit_type
			gt.add_data( (string) "");                    // 18 identifiers
			gt.add_data( dname);                          // 19 system is dname (can be a path)
			gt.add_data( (string) "CAD");                 // 20 factory
			gt.add_data(cf);                              // 21 variation is the full filename
			gt.add_data( (string) "1");                   // 22 run number

			dets[gt.data[0]] = get_detector(gt, gemcOpt, RC);

		}
	}

	// parsing attribute modifications. All cad imported volumes are stored in cad.gxml
	string fname = "cad.gxml";

	QFile *gxml = new QFile(fname.c_str());

	if( !gxml->exists() ) {

		for(auto &dirs: possibleDirs) {

			gxml = NULL;

			string fullName = dirs + fname;
			gxml = new QFile(fullName.c_str());

			// first found, exit
			if( gxml->exists() ) {
				continue;
			}
		}
	}
	// checking again
	if( !gxml->exists() ) {
			cout << hd_msg << " " << fname << " not found. All cad volumes are imported as original." << endl;
	} else {

		QDomDocument domDocument;
		// opening gcard and filling domDocument
		if(!domDocument.setContent(gxml)) {
			cout << " >>  xml format for file <" << fname << "> is wrong - check XML syntax. Exiting." << endl;
			exit(0);
		}
		gxml->close();

		QDomNodeList volumes = domDocument.firstChildElement().elementsByTagName("volume");
		for(int i = 0; i < volumes.count(); i++) {
			QDomNode elm = volumes.at(i);
			if(elm.isElement()) {
				QDomElement e = elm.toElement();
				string volumeName   = e.attribute("name").toStdString();
				string color        = e.attribute("color").toStdString();
				string material     = e.attribute("material").toStdString();
				string sensitivity  = e.attribute("sensitivity").toStdString();
				string identifiers  = e.attribute("identifiers").toStdString();
				string visible      = e.attribute("visible").toStdString();
				string style        = e.attribute("style").toStdString();
				string position     = e.attribute("position").toStdString();
				string rotation     = e.attribute("rotation").toStdString();

				// assigning attributes to volume
				if(dets.find(volumeName) != dets.end()) {

					if(verbosity>3)
						cout << " Modifying attributes for volume: " << volumeName ;


					if(color != "") {
						if(verbosity>3)
							cout << " color: " << color ;

						G4Colour thisCol = gcol(color);
						dets[volumeName].VAtts = G4VisAttributes(thisCol);
					}

					visible == "no" ? dets[volumeName].VAtts.SetVisibility(false) : dets[volumeName].VAtts.SetVisibility(true);
					style == "wireframe" ?  dets[volumeName].VAtts.SetForceWireframe(true) : dets[volumeName].VAtts.SetForceSolid(true);

					if(sensitivity != "") {
						if(verbosity>3)
							cout << " sensitivity: " << sensitivity ;

						// same hitType as sensitivity
						// setting system as sensitivity, so the hit definitions can be loaded
						// this should be modified later
						dets[volumeName].system = sensitivity;
						dets[volumeName].sensitivity = sensitivity;
						dets[volumeName].hitType = sensitivity;

						// identifier must be defined
						if(identifiers != "") {
							if(verbosity>3)
								cout << " identifiers: " << identifiers ;

							dets[volumeName].identity = get_identifiers(identifiers);

						} else {
							cout << " !! Error: volume " << volumeName << " has sensitivity but not identifier. " << endl;
						}
					}

					if(material != "") {
						if(verbosity>3)
						cout << " material: " << material ;
						dets[volumeName].material = material;
					}

					if(position != "") {
						if(verbosity>3)
							cout << " position: " << position ;
						dets[volumeName].pos = calc_position(position);
					}

					if(rotation != "") {
						if(verbosity>3)
							cout << " rotation: " << rotation ;
						dets[volumeName].rot = calc_rotation(rotation, volumeName);
					}

					if(verbosity>3)
						cout << endl;
				}
			}
		}
	}

	for(const auto &dd : dets)
		if(verbosity>3)
			cout << dd.second;

 	return dets;
}

string cad_det_factory::checkFormat(string file)
{
	// trying STL
	string filename = file + ".stl";
	ifstream fstl(filename.c_str());
	if(fstl.good()) {
		fstl.close();
		return filename;
		// trying .PLY
	} else {
		filename = file + ".ply";
		ifstream fply(filename.c_str());
		if(fply.good()) {
			fply.close();
			return filename;
			// trying OBJ
		} else {
			filename = file + ".obj";
			ifstream fobj(filename.c_str());
			if(fobj.good()) {
				fobj.close();
				return filename;
				// trying OBJ
			} else {
				return "na";
			}
		}
	}

	return "na";
}





