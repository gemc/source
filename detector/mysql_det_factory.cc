// Qt headers
#include <QtSql>

// gemc headers
#include "mysql_det_factory.h"
#include "utils.h"

map<string, detector> mysql_det_factory::loadDetectors()
{
	string hd_msg      = "  > MYSQL Detector Factory: >> ";
	double verbosity   = gemcOpt.optMap["GEO_VERBOSITY"].arg;
	string catch_v     = gemcOpt.optMap["CATCH"].args;

	map<string, detector> dets;
	
	// first check if there's at least one detector with MYSQL factory
	if(!check_if_factory_is_needed(RC.detectorConditionsMap, factoryType))
		return dets;

	// connection to the DB
	QSqlDatabase db = openGdb(gemcOpt);
			
	// connection is ok. Loading tables now
	// building detectors that are tagged with MYSQL factory
	for(map<string, detectorCondition>::iterator it=RC.detectorConditionsMap.begin(); it != RC.detectorConditionsMap.end(); it++)
	{
		if(it->second.get_factory() != factoryType)
			continue;
		
		string dname     = it->first;
		string tname     = dname + "__geometry";
		string variation = get_variation(it->second.get_variation());
		if(is_main_variation(it->second.get_variation())) 
			tname += "_main";
		int    run       = it->second.get_run_number();
			
		if(verbosity)
			cout <<  hd_msg << " Importing Detector: " <<  dname << " with <" << factoryType 
			     << "> factory, variation \"" << variation << "\", run number: " << run << endl;
				
		// getting the last id entry that matches variation and run number range
		int last_id = getLastId(db, tname, variation, run);
		
		string dbexecute  = "select name, mother, description, pos, rot, col, type, ";
		       dbexecute += "dimensions, material, magfield, ncopy, pmany, exist, ";
		       dbexecute += "visible, style, sensitivity, hitType, identity from " + tname;
		       dbexecute += " where variation ='" + variation + "'";
		       dbexecute += " and rmin <= " + stringify(run)  + " and rmax >= " + stringify(run)  ;
		       dbexecute += " and id = " + stringify(last_id);
						
		// executing query - will exit if not successfull.
		QSqlQuery q;
		if(!q.exec(dbexecute.c_str()))
		{
			cout << hd_msg << " !!! Failed to execute MYSQL query " << dbexecute <<  ". This is a fatal error. Exiting." << endl;
     		qDebug() << q.lastError();
			exit(0);
		}
		// Warning if nothing is found
		if(q.size() == 0 && verbosity)
		{
			cout << "  ** WARNING: detector \"" << dname << "\" not found with variation \"" << variation << "\" for run number " << run << endl << endl;
		}
		
						
		// else loading parameters from DB
		while (q.next())
		{
			gtable gt;
			
			// there should be 18 values in the MYSQL table
			for(int i=0; i<18; i++)
		     gt.add_data(q.value(i));

			// adding additional info: system, factory, variation, run number  
			gt.add_data(dname);
			gt.add_data((string) "MYSQL");
			gt.add_data(variation);
			gt.add_data(stringify(run));
			
			
			// big warning if detector already exist
			// detector is NOT loaded if already existing
			if(dets.find(gt.data[0]) != dets.end())
			{
				cout << endl <<  " *** WARNING! A detector >" << gt.data[0] << " exists already. Keeping original, not loading this instance. " << endl << endl;
			}
			else
			{
				dets[gt.data[0]] = get_detector(gt, gemcOpt, RC);
			}

			if(verbosity > 3)
				cout << get_detector(gt, gemcOpt, RC);
		}
	}

	
	// closing DB connection
	closeGdb(db);
	cout << endl;

	return dets;
}








