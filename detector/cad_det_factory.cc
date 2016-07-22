// gemc headers
#include "cad_det_factory.h"
#include "utils.h"


map<string, detector> cad_det_factory::loadDetectors()
{

	string hd_msg     = " >> CAD Factory: >> ";
	double verbosity  = gemcOpt.optMap["GEO_VERBOSITY"].arg;
	string catch_v    = gemcOpt.optMap["CATCH"].args;

	map<string, detector> dets;
	
  	// first check if there's at least one detector with CAD factory
	if(!check_if_factory_is_needed(RC.detectorConditionsMap, factoryType))
		return dets;

  	// there is at least one build with this factory.
	// building all detectors that are tagged with CAD factory
	for(map<string, detectorCondition>::iterator it=RC.detectorConditionsMap.begin(); it != RC.detectorConditionsMap.end() ; it++) {
		if(it->second.get_factory() != factoryType )
			continue;
		
		string dname     = it->first;
		if(verbosity)
			cout <<  hd_msg << " Importing Detector: <" <<  dname << "> with " << factoryType << " factory. "  << endl;

		// checking for the filename
		// this will be the detector variation
		string filename = dname + ".stl";
		ifstream fstl(filename.c_str());
		if(!fstl.good()) {
			// checking .ply format
			filename = dname + ".ply";
			ifstream fply(filename.c_str());
			if(!fply.good()) {
				// checking .obj format
				filename = dname + ".obj";
				ifstream fobj(filename.c_str());

				if(!fply.good()) {
					cout << " !! Error: no " << dname << ".stl or " << dname << ".ply file found. Exiting." << endl;
					exit(0);
				}
			}
		}

		// there is only one solid / cad file
		// using gtable to create it
		gtable gt;

		gt.add_data(dname);                            // 1 name
		gt.add_data( (string) "root");                 // 2 mother volume
		gt.add_data(dname + (string) " cadImported");  // 3 description
		gt.add_data( (string) "0*cm 0*cm 0*cm");       // 4 position
		gt.add_data( (string) "0*deg 0*deg 0*deg");    // 5 rotation
		gt.add_data( (string) "2222aa");               // 6 color
		gt.add_data( (string) "cadImport");            // 7 type
		gt.add_data( (string) "0");                    // 8 dimensions
		gt.add_data( (string) "G4_Al");                // 9 material is aluminum by defaul
		gt.add_data( (string) "no");                   // 10 magnetic field
		gt.add_data( (string) "0");                    // 11 copy number
		gt.add_data( (string) "0");                    // 12 pmany
		gt.add_data( (string) "1");                    // 13 activation flag
		gt.add_data( (string) "1");                    // 14 visibility
		gt.add_data( (string) "1");                    // 15 style
		gt.add_data( (string) "no");                   // 16 sensitivity
		gt.add_data( (string) "no");                   // 17 hit_type
		gt.add_data( (string) "");                     // 18 identifiers
		gt.add_data( (string) "dname");                // 19 system
		gt.add_data( (string) "CAD");                  // 20 factory
		gt.add_data(filename);                         // 21 variation
		gt.add_data( (string) "1");                    // 22 run number

		dets[gt.data[0]] = get_detector(gt, gemcOpt, RC);
	}

	// parsing attribute modifications. All cad imported volumes are stored in cad.gxml
	string fname = "cad.gxml";

	QFile gxml(fname.c_str());

	if( !gxml.exists() ) {
		cout << hd_msg << " " << fname << " not found. All cad volumes are imported as original." << endl;
	} else {

		QDomDocument domDocument;
		// opening gcard and filling domDocument
		if(!domDocument.setContent(&gxml)) {
			cout << " >>  xml format for file <" << fname << "> is wrong - check XML syntax. Exiting." << endl;
			exit(0);
		}
		gxml.close();

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

				// assigning attributes to volume
				if(dets.find(volumeName) != dets.end()) {

					//if(verbosity>3)
					cout << " Modifying attributes for volume: " << volumeName ;


					if(color != "") {
						//if(verbosity>3)
						cout << " color: " << color ;

						G4Colour thisCol = gcol(color);
						dets[volumeName].VAtts = G4VisAttributes(thisCol);
					}

					visible == "no" ? dets[volumeName].VAtts.SetVisibility(false) : dets[volumeName].VAtts.SetVisibility(true);
					style == "wireframe" ?  dets[volumeName].VAtts.SetForceWireframe(true) : dets[volumeName].VAtts.SetForceSolid(true);

					if(sensitivity != "") {
						cout << " sensitivity: " << sensitivity ;

						// same hitType as sensitivity
						dets[volumeName].sensitivity = sensitivity;
						dets[volumeName].hitType = sensitivity;

						// identifier must be defined
						if(identifiers != "") {
							cout << " identifiers: " << identifiers ;
							dets[volumeName].identity = get_identifiers(identifiers);

						} else {
							cout << " !! Error: volume " << volumeName << " has sensitivity but not identifier. " << endl;
						}
					}

					if(material != "") {
						//if(verbosity>3)
						cout << " material: " << material ;
						dets[volumeName].material = material;
					}

					//if(verbosity>3)
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

