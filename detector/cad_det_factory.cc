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
	for(map<string, detectorCondition>::iterator it=RC.detectorConditionsMap.begin(); it != RC.detectorConditionsMap.end() ; it++)
	{
		if(it->second.get_factory() != factoryType )
			continue;
		
		string dname     = it->first;
		if(verbosity)
			cout <<  hd_msg << " Importing Detector: <" <<  dname << "> with " << factoryType << " factory. "  << endl;





		// parsing attribute modifications
		string fname = dname + ".gxml";

		QFile gxml(fname.c_str());

		if( !gxml.exists() ) {
			cout << " > " << fname << " not found. All volumes in " << fname << " are imported as original." << endl;
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
	}

	for(const auto &dd : dets)
		if(verbosity>3)
			cout << dd.second;


 	return dets;
}

