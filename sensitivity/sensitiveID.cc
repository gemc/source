// Qt headers
#include <QtSql>

// gemc headers
#include "sensitiveID.h"
#include "utils.h"

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;


sensitiveID::sensitiveID(string SD, goptions gemcOpt, string factory, string variation, string s)
{
	double verbosity = gemcOpt.optMap["HIT_VERBOSITY"].arg;
	name             = SD;
	thisFactory      = factory + " " + variation;
	system           = s;
	
	// iF SD is FLUX, returns special sensitiveID
	if(SD == "flux")
	{
		description = "generic flux detector";
		identifiers.push_back("id");
		signalThreshold = 0;
		timeWindow      = 0;
		prodThreshold   = 1*mm;  // standard 1 mm production threshold
		maxStep         = 1*mm;
		return;
	}
	
	
	// MYSQL sensitivity infos
	if(factory == "TEXT")
	{
		string fname = system + "__hit_" + variation + ".txt";
		if(verbosity > 1) cout << "   > Loading TEXT definitions for <" << SD << ">..." << endl;

		ifstream IN(fname.c_str());
		if(!IN)
		{
			// if file is not found, maybe it's in the GEMC_DATA_DIR directory
			if(getenv("GEMC_DATA_DIR")  != NULL)
			{
				string maybeHere = (string) getenv("GEMC_DATA_DIR") + "/" + fname;
				
				IN.open(maybeHere.c_str());
				if(!IN && verbosity > 2)
				{
					cout << "  !!! Error: Failed to open hit file " << fname << " for sensitive detector: >"
					     << SD << "<. Maybe the filename doesn't exist?" << endl;
				}
			}
			if(!IN && verbosity > 2)
			{
				cout << "  !!! Error: Failed to open hit file " << fname << " for sensitive detector: >"
					 << SD << "<. Maybe the filename doesn't exist?" << endl;
			}
		}

		if(IN)
		{
			while(!IN.eof())
			{
				string dbline;
				getline(IN, dbline);

				if(!dbline.size())
					continue;

				gtable gt(get_strings(dbline, "|"));

				if(gt.data.size())
					if(gt.data[0] == SD)
					{
						// Reading variables
						// 0 is system name, by construction is SD

						// 1: description
						description = gt.data[1];

						// 2: Identifiers
						vector<string> ids = get_strings(gt.data[2]);
						for(unsigned i=0; i<ids.size(); i++)
							identifiers.push_back(ids[i]);

						// 3: Minimum Energy Cut for processing the hit
						signalThreshold = get_number(gt.data[3], 1);

						// 4: Time Window
						timeWindow = get_number(gt.data[4], 1);

						// 5: Production Threshold in the detector
						prodThreshold = get_number(gt.data[5], 1);

						// 6: Maximum Acceptable Step in the detector
						maxStep = get_number(gt.data[6], 1);

						// 7: rise time of the PMT signal
						riseTime = get_number(gt.data[7], 1);

						// 8: fall time of the PMT signal
						fallTime = get_number(gt.data[8], 1);

						// 9: from MeV to mV constant
						mvToMeV = get_number(gt.data[9]);

						// 10: pedestal
						pedestal = get_number(gt.data[10]);

						// 11: time from PMT face to signal
						delay = get_number(gt.data[11]);

					}
			}
			IN.close();
		}
		// default values if file is not present
		else
		{
			// 1: description
			description     = "unknown";
			signalThreshold = 1;
			timeWindow      = 100;
			prodThreshold   = 1;
			maxStep         = 10;
			riseTime        = 10;
			fallTime        = 20;
			mvToMeV         = 100;
			pedestal        = 100;
			delay           = 100;
		}
		if(verbosity > 3)
			// this will print the sensitive detector properties
			cout << *this << endl;;
		
		return;
	}
	
	
	
	// MYSQL sensitivity infos
	if(factory == "MYSQL")
	{
		// connection to the DB
		QSqlDatabase db = openGdb(gemcOpt);
		string tname    = system + "__hit";
		
		if(verbosity > 1) cout << "   > Loading MYSQL definitions for <" << SD << ">..." ;
		
		string dbexecute  = "select name, description, identifiers, signalThreshold, timeWindow, prodThreshold, maxStep, riseTime, fallTime, mvToMeV, pedestal, delay from " + tname ;
		dbexecute += " where variation ='" + variation + "'";
		dbexecute += " and name = '" + SD  + "'";
		
		QSqlQuery q;
		if(!q.exec(dbexecute.c_str()))
		{
			cout << " !!! Failed to execute MYSQL query " << dbexecute <<  ". This is a fatal error. Exiting." << endl;
			qDebug() << q.lastError();
			exit(0);
		}
		// Warning if nothing is found
		if(q.size() == 0 && verbosity)
		{
			cout << "  ** WARNING: sensitive detector \"" << SD << "\" not found in factory " << factory
			 << " for variation " << variation << endl << endl;
		}
		
		// else loading parameters from DB
		while (q.next())
		{
			// Reading variables
			// 0 is system name, by construction is SD
			
			// 1: description
			description = qv_tostring(q.value(1));
			
			// 2: Identifiers
			vector<string> ids = get_strings(qv_tostring(q.value(2)));
			for(unsigned i=0; i<ids.size(); i++)
				identifiers.push_back(ids[i]);
			
			// 3: Minimum Energy Cut for processing the hit
			signalThreshold = get_number(qv_tostring(q.value(3)));
			
			// 4: Time Window
			timeWindow = get_number(qv_tostring(q.value(4)));
			
			// 5: Production Threshold in the detector
			prodThreshold = get_number(qv_tostring(q.value(5)));
			
			// 6: Maximum Acceptable Step in the detector
			maxStep = get_number(qv_tostring(q.value(6)));
			
			// 7: rise time of the PMT signal
			riseTime = get_number(qv_tostring(q.value(7)));
			
			// 8: fall time of the PMT signal
			fallTime = get_number(qv_tostring(q.value(8)));
			
			// 9: from MeV to mV constant
			mvToMeV = get_number(qv_tostring(q.value(9)));
			
			// 10: pedestal
			pedestal = get_number(qv_tostring(q.value(10)));
			
			// 11: time from PMT face to signal
			delay = get_number(qv_tostring(q.value(11)));
						
		}
		
		// closing DB connection
		closeGdb(db);
		if(verbosity > 3)
			cout << *this << endl;
		
		return;
	}
	
}


ostream &operator<<(ostream &stream, sensitiveID SD)
{
	cout << "  > Sensitive detector " << SD.name << ": " << endl << endl ;
	for(unsigned int i=0; i<SD.identifiers.size(); i++)
	{
		cout << "  identifier element name:  " << SD.identifiers[i] << endl;
	}
	cout << endl << "  >  Signal Threshold:        " << SD.signalThreshold/MeV << " MeV." << endl ;
	cout         << "  >  Production Threshold:    " << SD.prodThreshold/mm << " mm." << endl ;
	cout         << "  >  Time Window for:         " << SD.timeWindow/ns << " ns." << endl ;
	cout         << "  >  Maximum Acceptable Step: " << SD.maxStep/mm  << " mm." << endl ;
	cout         << "  >  Signal Rise Time:        " << SD.riseTime/ns << " ns." << endl ;
	cout         << "  >  Signal Fall Time:        " << SD.fallTime/ns << " ns." << endl ;
	cout         << "  >  Signal MeV to mV:        " << SD.mvToMeV     << "  mV/MeV" << endl ;
	cout         << "  >  Signal Pedestal:         " << SD.pedestal    << " mV" << endl ;
	cout         << "  >  Signal Delay             " << SD.delay/ns    << " ns." << endl ;
	cout << endl;
	
	return stream;
}












