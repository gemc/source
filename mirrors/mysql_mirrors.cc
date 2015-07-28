// Qt headers
#include <QtSql>

// gemc headers
#include "mirrors_factory.h"
#include "mysql_mirrors.h"
#include "string_utilities.h"
#include "utils.h"


map<string, mirror*> mysql_mirrors::initMirrors(runConditions rc, goptions opts)
{
	
	string hd_msg    = opts.optMap["LOG_MSG"].args + " MYSQL mirrors Factory: >> ";
	double verbosity = opts.optMap["MIRROR_VERBOSITY"].arg;
	
	map<string, mirror*> mymirs;  // mirror map
	
	// first check if there's at least one detector with MYSQL factory
	if(!check_if_factory_is_needed(rc.detectorConditionsMap, "MYSQL"))
		return mymirs;

	// connection to the DB
	QSqlDatabase db = openGdb(opts);
	
	// Looping over detectorConditionsMap for detector names
	// To each detector is associated a mirror table and an opt properties table
	for(map<string, detectorCondition>::iterator it=rc.detectorConditionsMap.begin(); it != rc.detectorConditionsMap.end(); it++)
	{
		// building mirrors belonging to detectors that are tagged with MYSQL factory
		if(it->second.get_factory() != "MYSQL")
			continue;

		if(verbosity)
			cout << hd_msg << " Initializing " << it->second.get_factory() << " for detector " << it->first << endl;

		// only add "main" if it's the main variation
		string dname     = it->first ;
		string tname     = dname + "__mirrors";
		string variation = get_variation(it->second.get_variation());
		if(is_main_variation(it->second.get_variation()))
			tname += "_main";

		string dbexecute = "select name, desc, type, finish, model, border, maptOptProps, photonEnergy, indexOfRefraction, reflectivity, efficiency, specularlobe, specularspike, backscatter from " + tname;
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
			cout << "  ** WARNING: mirror for system \"" << dname << "\" not found with variation \"" << variation << endl << endl;
		}
		
		while (q.next())
		{
			mirror *thisMir = new mirror(TrimSpaces(qv_tostring( q.value(0))));         // name
			thisMir->desc          =                 qv_tostring(q.value(1));           // description
			thisMir->type          =                 qv_tostring(q.value(2));           // type
			thisMir->finish        =                 qv_tostring(q.value(3));           // finish
			thisMir->model         =                 qv_tostring(q.value(4));           // model
			thisMir->border        =                 qv_tostring(q.value(5));           // border
			thisMir->maptOptProps  =                 qv_tostring(q.value(6));           // maptOptProps
			thisMir->opticalsFromString(qv_tostring(             q.value(7)),  "photonEnergy");
			thisMir->opticalsFromString(qv_tostring(             q.value(8)),  "indexOfRefraction");
			thisMir->opticalsFromString(qv_tostring(             q.value(9)),  "reflectivity");
			thisMir->opticalsFromString(qv_tostring(             q.value(10)), "efficiency");
			thisMir->opticalsFromString(qv_tostring(             q.value(11)), "specularlobe");
			thisMir->opticalsFromString(qv_tostring(             q.value(12)), "specularspike");
			thisMir->opticalsFromString(qv_tostring(             q.value(13)), "backscatter");
			
			mymirs[thisMir->name] = thisMir;
			
		}
	}
	
	// closing DB connection
	closeGdb(db);
	cout << endl;
	
 	if(verbosity>0) printMirrors(mymirs);

	return mymirs;
}




