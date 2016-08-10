// Qt headers
#include <QtSql>

// gemc headers
#include "parameter_factory.h"
#include "mysql_parameters.h"
#include "string_utilities.h"
#include "utils.h"

// mlibrary
#include "gstring.h"
using namespace gstring;

map<string, double> mysql_parameters::loadParameters(goptions opts, runConditions RC)
{
	string hd_msg    = "  > MYSQL Parameters: >> ";
	double verbosity = opts.optMap["PARAMETER_VERBOSITY"].arg;
	
	map<string, double> GParameters;   // parameters maps
	// first check if there's at least one detector with MYSQL factory
	if(!check_if_factory_is_needed(RC.detectorConditionsMap, factoryType))
		return GParameters;
	
	// connection to the DB
	openGdb(opts);
	QSqlDatabase db = QSqlDatabase::database("gemc");
	
	// Looping over detectorConditionsMap for detector names
	// To each detector is associated a material table and an opt properties table
	for(map<string, detectorCondition>::iterator it=RC.detectorConditionsMap.begin(); it != RC.detectorConditionsMap.end(); it++)
	{
		if(it->second.get_factory() != factoryType)
			continue;
		
		string dname     = it->first ;
		string tname     = dname + "__parameters";
		string variation = get_variation(it->second.get_variation());
		if(is_main_variation(it->second.get_variation()))
			tname += "_main";
		int    run       = it->second.get_run_number();
		
		if(verbosity)
			cout <<  hd_msg << " Importing Parameters for Detector: " <<  dname
			<< " with " << factoryType << " factory, variation " << variation << endl;
		
		string dbexecute  = "select name, value, units, description from " + tname;
		dbexecute += " where variation ='" + variation + "'";
		dbexecute += " and rmin <= " + stringify(run)  + " and rmax >= " + stringify(run)  ;
		dbexecute += " order by id desc limit 1 " ;
		
		// executing query - will exit if not successfull.
		QSqlQuery q;
		if(!q.exec(dbexecute.c_str()))
		{
			cout << hd_msg << "  Errror! Failed to execute MYSQL query " << dbexecute << endl;
			exit(0);
		}
		// Warning if nothing is found
		if(q.size() == 0 && verbosity)
		{
			cout << "  ** WARNING: detector \"" << dname << "\" not found with variation \"" << variation << "\" for run number " << run << endl << endl;
		}
		
		while (q.next())
		{
			
			gtable gt;
			gt.add_data(dname + "/" + trimSpacesFromString(qv_tostring(q.value(0))));
			gt.add_data(q.value(1));
			gt.add_data(q.value(2));
			gt.add_data(q.value(3));
			
			GParameters[gt.data[0]] = get_par_value(gt);
			
			if(verbosity > 1)
				log_value(gt, factoryType);
		}
	}
	
  	// closing DB connection
	db.close();
	QSqlDatabase::removeDatabase("QMYSQL");
	
	return GParameters;
}
