
// gemc headers
#include "gdml_det_factory.h"
#include "utils.h"

// geant4
#include "G4GDMLParser.hh"
#include "G4LogicalVolume.hh"

map<string, detector> gdml_det_factory::loadDetectors()
{

	string hd_msg     = " >> GDML Factory: >> ";
	double verbosity  = gemcOpt.optMap["GEO_VERBOSITY"].arg;
	string catch_v    = gemcOpt.optMap["CATCH"].args;

	map<string, detector> dets;
	
  	// first check if there's at least one detector with GDML factory
	if(!check_if_factory_is_needed(RC.detectorConditionsMap, factoryType))
		return dets;
	
  	// there is at least one build with this factory.
	// building all detectors that are tagged with GDML factory
	for(map<string, detectorCondition>::iterator it=RC.detectorConditionsMap.begin(); it != RC.detectorConditionsMap.end() ; it++)
	{
		if(it->second.get_factory() != factoryType )
			continue;
		
		string dname     = it->first;
		if(verbosity)
			cout <<  hd_msg << " Importing Detector: <" <<  dname << "> with " << factoryType << " factory. "  << endl;



		// parsing GDML file to build the map.
		// This will store in the detector class all the logical and physical volumes
		// At Construct() time, only the physical volumes to be G4PVPlaced will be
		// the ones with mom = World

		G4GDMLParser *parser = new G4GDMLParser();
		string gname     = dname + ".gdml";

		// parsing G4 volumes
		parser->Read(gname, 0);

		// the volume name has to be "World"
		// its origin are "root" coordinate
		G4LogicalVolume* gdmlWorld = parser->GetVolume("World");

		// first daughters: these volumes will be the ones with mother = "root"
		for(int d=0; d<gdmlWorld->GetNoDaughters (); d++) {


			string thisDetName = gdmlWorld->GetDaughter(d)->GetLogicalVolume()->GetName();
			string momName =  gdmlWorld->GetDaughter(d)->GetMotherLogical()->GetName() ;

			detector thisDet = get_detector(gdmlWorld->GetDaughter(d),  gemcOpt, RC);
			thisDet.mother = "root";
			thisDet.variation = gname;

			if(verbosity>2) cout << thisDet ;

			// now browsing for daughters
			G4LogicalVolume* firstDaughter = gdmlWorld->GetDaughter(d)->GetLogicalVolume();
			for(int gd=0; gd<firstDaughter->GetNoDaughters (); gd++) {
				cout <<  " gd: " << firstDaughter->GetDaughter(d)->GetLogicalVolume()->GetName() << endl;
			}

			dets[thisDet.name] = thisDet;
		}



		// parsing attribute modifications
		string fname = dname + ".gxml";

		QFile gxml(fname.c_str());

		if( !gxml.exists() ) {
			cout << " > " << fname << " not found. All volumes in " << fname << " are imported as original." << endl;
		} else {

			QDomDocument domDocument;
			// opening gcard and filling domDocument
			if(!domDocument.setContent(&gxml))
			{
				cout << " >>  xml format for file <" << fname << "> is wrong - check XML syntax. Exiting." << endl;
				exit(0);
			}
			gxml.close();

			QDomNodeList volumes = domDocument.firstChildElement().elementsByTagName("volume");
			for(int i = 0; i < volumes.count(); i++)
			{
				QDomNode elm = volumes.at(i);
				if(elm.isElement())
				{
					QDomElement e = elm.toElement();
					string volumeName   = e.attribute("name").toStdString();
					string color        = e.attribute("color").toStdString();
					string material     = e.attribute("material").toStdString();
					string sensitivity  = e.attribute("sensitivity").toStdString();
					string identifiers  = e.attribute("identifiers").toStdString();

					if(verbosity>3) {
						cout << " Volume: " << volumeName
						<< " color: " << color
						<< " material: " << material
						<< " sensitivity: " << sensitivity
						<< " identifiers: " << identifiers << endl;
					}
				}
				
			}
		}

 	}






 	return dets;
}

