// gemc headers
#include "run_conditions.h"
#include "string_utilities.h"
#include "utils.h"

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;

// geant4 
#include "Randomize.hh"


runConditions::runConditions(goptions gemcOpt)
{
	string hd_msg    = gemcOpt.optMap["LOG_MSG"].args + " gcard: >> " ;
	int run_number   = gemcOpt.optMap["RUNNO"].arg;
	
	
	// adding detector defined with the DF directive Detector/Factory
	vector<string> dfopt = get_info(gemcOpt.optMap["DF"].args);
	if(dfopt[0] != "no" && dfopt.size() > 1)
	{
		detectorConditionsMap[dfopt[0]] = detectorCondition(dfopt[1]);
	}
	
	string file = gemcOpt.optMap["gcard"].args;
	if(file=="no") return;
	
	QFile gcard(file.c_str());
	// reading gcard and checking the consistency
	if (!domDocument.setContent(&gcard))
	{
		gcard.close();
		cout << hd_msg << " gcard format is wrong - check XML syntax. Exiting." << endl;
		exit(0);
	}
	gcard.close();
	
	// start reading gcard file
	QDomElement docElem = domDocument.documentElement();
	
	QDomNode n = docElem.firstChild();
	// looping over detector  entries
	while(!n.isNull())
	{
		QDomElement e = n.toElement();                    // converts the node to an element.
		if(!e.isNull())                                   // the node really is an element.
			if(qv_tostring(e.tagName()) == "detector")    // selecting "detector" nodes
			{
				detectorCondition dCon;
				dCon.set_system(assignAttribute(e,     "name", "na"));
				dCon.set_factory(assignAttribute(e,    "factory", "na"));
				dCon.set_variation(assignAttribute(e,  "variation", "main"));  // default variation is "main"
				dCon.set_run_number(assignAttribute(e, "run_number", run_number));
				
				// Initializing Position and rotation
				dCon.set_position("0*cm",  "0*cm",  "0*cm");
				dCon.set_rotation("0*deg", "0*deg", "0*deg");
				
				
				QDomNode nn = e.firstChild();               ///< checking detector attributes
				while(!nn.isNull())
				{
					QDomElement ee = nn.toElement();
					if(!ee.isNull())
					{
						if(qv_tostring(ee.tagName()) == "position")  ///< Initializing Position
						{
							dCon.set_position(qv_tostring(ee.attributeNode("x").value()),
											  qv_tostring(ee.attributeNode("y").value()),
											  qv_tostring(ee.attributeNode("z").value()) );
						}
						if(qv_tostring(ee.tagName()) == "rotation")  ///< Initializing Rotation
						{
							dCon.set_rotation(qv_tostring(ee.attributeNode("x").value()),
											  qv_tostring(ee.attributeNode("y").value()),
											  qv_tostring(ee.attributeNode("z").value()) );
						}
						if(qv_tostring(ee.tagName()) == "existence")  ///< Initializing Existance
						{
							dCon.set_existance(qv_tostring(ee.attributeNode("exist").value()));
						}
					}
					nn = nn.nextSibling();
				}
				
				// inserting this detector in the map
				detectorConditionsMap[assignAttribute(e, "name", "na")] = dCon;
				
			}
		
		n = n.nextSibling();
	}
	
	cout << endl;
}

runConditions::~runConditions(){}


void detectorCondition::set_position(string X, string Y, string Z)
{
	pos.setX(get_number(X,1));
	pos.setY(get_number(Y,1));
	pos.setZ(get_number(Z,1));
}

void detectorCondition::set_rotation(string X, string Y, string Z)
{
	rot = G4RotationMatrix(G4ThreeVector(1, 0, 0),
						   G4ThreeVector(0, 1, 0),
						   G4ThreeVector(0, 0, 1));
	
	rot.rotateX(get_number(X,1));
	rot.rotateY(get_number(Y,1));
	rot.rotateZ(get_number(Z,1));
	
	vrot.setX(get_number(X,1));
	vrot.setY(get_number(Y,1));
	vrot.setZ(get_number(Z,1));
	
}

void detectorCondition::set_existance(string exist)
{
	presentFlag = true;
	if(exist == "no" || exist == "NO" || exist == "No")
		is_present = 0;
}


map<string, string> runConditions::getDetectorConditionsMap()
{
	map<string, string> detmap;
	
	// filling name, rotation, position modifications from gcard
	for(map<string, detectorCondition>::iterator it = detectorConditionsMap.begin(); it != detectorConditionsMap.end(); it++)
	{
		
		detmap[it->first] = " is loaded with factory " +  it->second.get_factory()
		+ ", variation "     + it->second.get_variation()
		+ " and run number " + stringify(it->second.get_run_number());
		
		if(it->second.get_position().mag2() != 0)
		{
			string key = "local shift for " + it->first;
			detmap[key] = "(" + stringify(it->second.get_position().x()/mm) + ", "
			+ stringify(it->second.get_position().y()/mm) + ", "
			+ stringify(it->second.get_position().z()/mm) + ")mm";
		}
		
		if(it->second.get_vrotation().mag2() != 0)
		{
			string key = "local rotation for " + it->first;
			detmap[key] = "(" + stringify(it->second.get_vrotation().x()/degree) + ", "
			+ stringify(it->second.get_vrotation().y()/degree) + ", "
			+ stringify(it->second.get_vrotation().z()/degree) + ")deg";
		}
	}
	
	return detmap;
}

int check_if_factory_is_needed(map<string, detectorCondition> dcon, string factory)
{
	int isneeded = 0;
	
	for(map<string, detectorCondition>::iterator it=dcon.begin(); it != dcon.end(); it++)
	{
		if(it->second.get_factory() == factory)
			isneeded++;
	}
	return isneeded;
}

int runConditions::get_run_number(string detector)
{
	if(detectorConditionsMap.find(detector) != detectorConditionsMap.end())
		return detectorConditionsMap[detector].get_run_number();
	
	else return 0;
}


string runConditions::get_variation(string detector)
{
	if(detectorConditionsMap.find(detector) != detectorConditionsMap.end())
		return detectorConditionsMap[detector].get_variation();
	
	else return "na";
}

string runConditions::get_system(string detector)
{
	if(detectorConditionsMap.find(detector) != detectorConditionsMap.end())
		return detectorConditionsMap[detector].get_system();
	
	else return "na";
		
}

// get_systems return a map<string, string>
map<string, string> runConditions::get_systems()
{
	map<string, string> allSystems;
	
	for(map<string, detectorCondition>::iterator it=detectorConditionsMap.begin(); it != detectorConditionsMap.end(); it++)
		allSystems[it->first] = it->second.get_factory();
		
	return allSystems;
}



runWeights::runWeights(goptions opts)
{
	string fname     = opts.optMap["RUN_WEIGHTS"].args;
	int nevts        = opts.optMap["N"].arg;
	startEvent       = opts.optMap["EVN"].arg;
	defaultRunNumber = opts.optMap["RUNNO"].arg;
	
	isNewRun = FALSE;
	
	if(nevts==0) return;
	
	// by default there is only 1 weight, the run number is defaultRunNumber for all events
	if(fname == "no")
	{
		w[defaultRunNumber]  = 1;
		n[defaultRunNumber]  = nevts;
		return;
	}
	
	ifstream in(fname.c_str());
	
	if(!in)
	{
		cerr << " !!! Can't open input file " << fname.c_str() << ". Exiting. " << endl;
		exit(1);
	}
	
	// filling weight map
	while (!in.eof())
	{
		int rr;
		double ww;
		in >> rr >> ww;
		w[rr] = ww;
		n[rr] = 0;
	}
	
	// now randomizing the run / event map based on number of events
	for(int i=0; i<nevts; i++)
	{
		double random = G4UniformRand();
		double ww = 0;
		for(map<int, double>::iterator it = w.begin(); it != w.end(); it++)
		{
			ww += it->second;
			if(random <= ww)
			{
				n[it->first]++;
				break;
			}
		}
	}
		
	cout << " > Run weights table loaded: " << endl;
	for(map<int, double>::iterator it = w.begin(); it != w.end(); it++)
		cout << "    - run: " << it->first << "\t weight: " << w[it->first] << "\t  n. events: " << n[it->first] << endl;
	
}


int runWeights::getRunNumber(int evn)
{
	int dn = evn - startEvent ;
	
	double nn = 0;
	
	for(map<int, int>::iterator it = n.begin(); it != n.end(); it++)
	{
		nn += it->second;
		if(dn < nn)
		{
			if(it->first != runNo)
			{
				isNewRun = TRUE;
				runNo = it->first;
			}
			else
			{
				isNewRun = FALSE;
			}
			return it->first;
		}
	}
	
	// default comes from the option map
	runNo = defaultRunNumber;
	
	return runNo;
}












