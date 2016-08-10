// Qt headers
#include <QtSql>

// gemc headers
#include "material_factory.h"
#include "mysql_materials.h"
#include "string_utilities.h"
#include "utils.h"

// mlibrary
#include "gstring.h"
using namespace gstring;

// G4 headers
#include "G4Element.hh"
#include "G4NistManager.hh"
#include "G4OpBoundaryProcess.hh"


map<string, G4Material*> mysql_materials::initMaterials(runConditions rc, goptions opts)
{
	
	string hd_msg    = opts.optMap["LOG_MSG"].args + " MYSQL Materials Factory: >> ";
	double verbosity = opts.optMap["MATERIAL_VERBOSITY"].arg;
	
	map<string, material> mymats;                        // material map
	
	// first check if there's at least one detector with MYSQL factory
	if(!check_if_factory_is_needed(rc.detectorConditionsMap, "MYSQL"))
		return materialsFromMap(mymats);

	// connection to the DB
	QSqlDatabase db = openGdb(opts);
	
	// Looping over detectorConditionsMap for detector names
	// To each detector is associated a material and (optional) opt properties
	for(map<string, detectorCondition>::iterator it=rc.detectorConditionsMap.begin(); it != rc.detectorConditionsMap.end(); it++)
	{
		// building materials belonging to detectors that are tagged with MYSQL factory
		if(it->second.get_factory() != "MYSQL")
			continue;

		if(verbosity)
			cout << hd_msg << " Initializing " << it->second.get_factory() << " for detector " << it->first << endl;

		// only add "main" if it's the main variation
		string dname     = it->first ;
		string tname     = dname + "__materials";
		string variation = get_variation(it->second.get_variation());
		if(is_main_variation(it->second.get_variation()))
			tname += "_main";

		string dbexecute  = "select name, description, density, ncomponents, components, photonEnergy, indexOfRefraction, ";
				 dbexecute += "absorptionLength, reflectivity, efficiency, fastcomponent, slowcomponent, ";
				 dbexecute += "scintillationyield, resolutionscale, fasttimeconstant, slowtimeconstant, yieldratio from " + tname;
				 dbexecute += " where variation ='" + variation + "'";
				
		// executing query - will exit if not successfull.
		QSqlQuery q;
		if(!q.exec(dbexecute.c_str()))
		{
			cout << hd_msg << "  Failed to execute MYSQL query " << dbexecute <<  ". This is a fatal error. Exiting." << endl;
     		qDebug() << q.lastError();
			exit(0);
		}
		// Warning if nothing is found
		if(q.size() == 0 && verbosity)
		{
			cout << "  ** WARNING: material for system \"" << dname << "\" not found with variation \"" << variation << endl << endl;
		}
		
		while (q.next())
		{
			material thisMat(trimSpacesFromString(qv_tostring( q.value(0))));         // name
			thisMat.desc               = qv_tostring(q.value(1));           // description
			thisMat.density            =             q.value(2).toDouble(); // density
			thisMat.ncomponents        =             q.value(3).toInt();    // number of components
			thisMat.componentsFromString(qv_tostring(q.value(4)));          // component + quantity list
			thisMat.opticalsFromString(qv_tostring(  q.value(5)), "photonEnergy");
			thisMat.opticalsFromString(qv_tostring(  q.value(6)), "indexOfRefraction");
			thisMat.opticalsFromString(qv_tostring(  q.value(7)), "absorptionLength");
			thisMat.opticalsFromString(qv_tostring(  q.value(8)), "reflectivity");
			thisMat.opticalsFromString(qv_tostring(  q.value(9)), "efficiency");
			
			// scintillation
			thisMat.opticalsFromString(qv_tostring(  q.value(10)), "fastcomponent");
			thisMat.opticalsFromString(qv_tostring(  q.value(11)), "slowcomponent");
			thisMat.scintillationyield =             q.value(12).toDouble(); // scintillationyield
			thisMat.resolutionscale    =             q.value(13).toDouble(); // resolutionscale
			thisMat.fasttimeconstant   =             q.value(14).toDouble(); // fasttimeconstant
			thisMat.slowtimeconstant   =             q.value(15).toDouble(); // slowtimeconstant
			thisMat.yieldratio         =             q.value(16).toDouble(); // yieldratio
			thisMat.opticalsFromString(qv_tostring(  q.value(17)), "rayleigh");
			
			mymats[thisMat.name] = thisMat;
			
		}
	}
	
	// closing DB connection
	closeGdb(db);
	cout << endl;
	
	map<string, G4Material*> returnMap = materialsFromMap(mymats);
	if(verbosity>0) printMaterials(returnMap);
	
	return returnMap;
}




