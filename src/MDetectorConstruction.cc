// G4 headers
#include "G4Box.hh"
#include "G4GeometryManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4ProductionCuts.hh"
#include "G4PVPlacement.hh"
#include "G4RunManager.hh"
#include "G4SolidStore.hh"
#include "G4VisAttributes.hh"
#include "G4SDManager.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4RegionStore.hh"

// gemc headers
#include "MDetectorConstruction.h"
#include "MDetectorMessenger.h"

// C++ headers
#include <sstream>
using namespace std;

MDetectorConstruction::MDetectorConstruction(goptions Opts)
{
	gemcOpt = Opts;
}

MDetectorConstruction::~MDetectorConstruction()
{
	;
}


G4VPhysicalVolume* MDetectorConstruction::Construct()
{
	string hd_msg     = gemcOpt.optMap["LOG_MSG"].args + " >> Construction:";
	double VERB       = gemcOpt.optMap["G4P_VERBOSITY"].arg ;
	string catch_v    = gemcOpt.optMap["CATCH"].args;
	string hall_mat   = gemcOpt.optMap["HALL_MATERIAL"].args;
	string hall_field = gemcOpt.optMap["HALL_FIELD"].args;
	
	// Clean old geometry, if any
	G4GeometryManager::GetInstance()->OpenGeometry();
	G4PhysicalVolumeStore::GetInstance()->Clean();
	G4LogicalVolumeStore::GetInstance()->Clean();
	G4SolidStore::GetInstance()->Clean();
	
	// Experimental hall is a 20mx2 = 40 meters box
	// dimensions coming from HALL_DIMENSIONS options
	(*hallMap)["root"].create_solid(gemcOpt, hallMap);
	(*hallMap)["root"].create_logical_volume(mats, gemcOpt);
	(*hallMap)["root"].create_physical_volumes(gemcOpt, NULL);
	hasMagfield((*hallMap)["root"]);
	
	if(VERB > 3 || catch_v == "root") cout << hd_msg << "    " << (*hallMap)["root"] ;
	
	
	cout << hd_msg <<  " Building Detector from Geometry STL Map..." << endl;
	
	
	// ########################################################################
	// Resetting Detector Map "scanned". Propagating "exist" to all generations
	// ########################################################################
	if(VERB > 2) cout << hd_msg << " Mapping Physical Detector..." << endl << endl;
	
	for(map<string, detector>::iterator i =  hallMap->begin(); i!=hallMap->end(); i++)
	{
		if(VERB > 3) cout << hd_msg << " Scanning Detector " << i->first << " - existance: " << i->second.exist << endl;
		
		// Find the mother up to "root" - disable kid if ancestor does not exist
		detector mother = findDetector(i->second.mother);
		
		while(mother.name != "akasha" && mother.name != "notfound")
		{
			if(mother.exist == 0)
			{
				if(VERB > 2) cout << hd_msg <<   "\t" << i->second.mother  << " is not activated. Its child " << i->second.name
					<< " will be disactivated as well." << endl;
				i->second.exist = 0;
			}
			mother = findDetector(mother.mother);
		}
		if(i->first != "root") i->second.scanned = 0;

		// scanning for replica, replicants physical volumes are not built
		if(i->second.type.find("ReplicaOf:") != string::npos)
		{
			stringstream ops;
			string operands(i->second.type, 10, i->second.type.size());
			ops << operands;
			string original;
			ops >> original;

			replicants.insert(original);
		}
		
	}
	
	
	// ########################################################################
	// Building Solids, Logical Volumes, Physical Volumes from the detector Map
	// ########################################################################
	vector<string> relatives;
	string mom, kid;
	
	vector<string> regions;  // all volumes for which mom is "root"

	for( map<string, detector>::iterator i =  hallMap->begin(); i!=hallMap->end(); i++)
	{
		// don't build anything if the exist flag is not set
		if(i->second.exist == 0) continue;
		
		// put the first volume in relatives
		// typically it's the first in alphabetical order
		if(i->first != "root") relatives.push_back(i->second.name);
		
		while(relatives.size() > 0)
		{
			detector kid = findDetector(relatives.back());
			detector mom = findDetector(kid.mother);
			
			// if the mother system is different than the kid system
			// then the kid will define a new region
			if(mom.system != kid.system && kid.material != "Component")
				regions.push_back(kid.name);

			
			// Mom doesn't exists in the hallMap. Stopping everything.
			if(mom.name != "akasha"  && mom.name == "notfound")
			{
				cout << hd_msg << "  Mom was not found for <" << relatives.back() << ">. "
				     << " We have a No Child Left Behind policy. Exiting. " << endl << endl;
				exit(0);
			}
			
			// output the Geneaology
			if(VERB > 3 || kid.name.find(catch_v) != string::npos)
			{
				for(unsigned int ir=0; ir<relatives.size()-1; ir++) cout << "\t";
				cout << hd_msg << " Checking " << kid.name << ", child of " << mom.name
				               << ", for a living ancestor. "
				               << " This Geneaology Depth is " << relatives.size() << "." << endl;
			}
			
			// Mom is built, kid not built yet.
			if(kid.scanned == 0 && mom.scanned == 1)
			{
				if(VERB > 3 || kid.name.find(catch_v) != string::npos)
				{
					for(unsigned int ir=0; ir<relatives.size()-1; ir++) cout << "\t";
					cout << hd_msg << "  Found:  " << kid.name
					     << " is not built yet but its mommie " << mom.name << " is."
					     << " Building " << kid.name << " of type: "  << kid.type << "..." << endl;
				}
				
				// Check kid dependency on copies
				if(kid.type.find("CopyOf") == 0)
				{
					
					stringstream ops;
					string operands(kid.type, 6, kid.type.size());
					ops << operands;
					string original;
					ops >> original;
					
					detector dorig = findDetector(original);
					// if dependency is not built yet, then
					// add it to the relative list
					if(dorig.scanned == 0)
					{
						relatives.push_back(original);
						if(VERB > 3  || kid.name.find(catch_v) != string::npos)
						{
							for(unsigned int ir=0; ir<relatives.size()-1; ir++) cout << "\t";
							cout << hd_msg << kid.name << " is copied volume. "
							<< " Must build: " << original << " first " << endl;
						}
					}
					// otherwise can build the kid
					else
					{
						buildDetector(kid.name);
						kid.scanned = 1;
					}
				}
				
				// Check kid dependency on replicas
				if(kid.type.find("ReplicaOf:") == 0)
				{
					
					stringstream ops;
					string operands(kid.type, 10, kid.type.size());
					ops << operands;
					string original;
					ops >> original;
					
					detector dorig = findDetector(original);
					// if dependency is not built yet, then
					// add it to the relative list
					if(dorig.scanned == 0)
					{
						relatives.push_back(TrimSpaces(original));
						
						if(VERB > 3  || kid.name.find(catch_v) != string::npos)
						{
							for(unsigned int ir=0; ir<relatives.size()-1; ir++) cout << "\t";
							cout << hd_msg << kid.name << " is a replica volume of  " << original
							<< ". Must build: " << original << " first " << endl;
						}
					}
					// otherwise can build the kid
					else
					{
						buildDetector(kid.name);
						kid.scanned = 1;
					}
				}
				
				
				
				// Check kid dependency on operations
				else if(kid.type.find("Operation:") == 0)
				{
					int sstart = 10;
					if(kid.type.find("Operation:~") == 0 || kid.type.find("Operation:@") == 0 ) sstart = 11;
					
					string operation(kid.type, sstart, kid.type.size());
					vector<string> operands = get_strings(replaceCharWithChars(operation, "-+*/", " "));
					string firstop  = operands[0];
					string secondop = operands[1];
					
					// if dependency is not built yet, then
					// add it to the relative list
					detector dsecondop = findDetector(secondop);
					if(dsecondop.scanned == 0)
					{
						relatives.push_back(secondop);
						if(VERB > 3  || kid.name.find(catch_v) != string::npos)
						{
							for(unsigned int ir=0; ir<relatives.size()-1; ir++) cout << "\t";
							cout << hd_msg << kid.name << " is the result of an operation. "
							<< " Must build: " << secondop << " first " << endl;
						}
					}
					detector dfirstop = findDetector(firstop);
					if(dfirstop.scanned == 0)
					{
						relatives.push_back(firstop);
						if(VERB > 3  || kid.name.find(catch_v) != string::npos)
						{
							for(unsigned int ir=0; ir<relatives.size()-1; ir++) cout << "\t";
							cout << hd_msg << kid.name << " is the result of an operation. "
							<< " Must build: " << firstop << " first " << endl;
						}
					}
					// otherwise can build the kid
					else if(dsecondop.scanned == 1 && dfirstop.scanned == 1)
					{
						buildDetector(kid.name);
						kid.scanned = 1;
					}
				}
				// no dependencies found, build the kid
				else
				{
					buildDetector(kid.name);
					kid.scanned = 1;
				}
			}
			
			// if the kid still doesn't exists and its mom doesn't exist.
			// adding mom to the relatives list
			if(kid.scanned == 0 && mom.scanned == 0)
			{
				relatives.push_back(kid.mother);
			}
			
			// the kid has been built. Can go down one step in geneaology
			if(kid.scanned == 1 && relatives.size())
			{
				if(VERB > 3 || kid.name.find(catch_v) != string::npos)
					cout << hd_msg  << " " <<  kid << " is built." <<  endl << endl;
				
				relatives.pop_back();
			}
		}
	}
	
	// build mirrors
	buildMirrors();
	
	
	// assign regions
	// includes root
	regions.push_back("root");
	assignRegions(regions);
	
	
	return (*hallMap)["root"].GetPhysical();
}

#include "G4UserLimits.hh"

void MDetectorConstruction::isSensitive(detector detect)
{
	string hd_msg  = gemcOpt.optMap["LOG_MSG"].args + " Sensitivity: >> ";
	double VERB    = gemcOpt.optMap["HIT_VERBOSITY"].arg ;
	string catch_v = gemcOpt.optMap["CATCH"].args;
	
	string sensi   = detect.sensitivity;
	
	if(sensi != "no" && sensi.find("mirror:") == string::npos)
	{
		G4SDManager* SDman = G4SDManager::GetSDMpointer();
		if (!SDman) throw "Can't get G4SDManager";
		
		if(VERB > 5 || detect.name.find(catch_v) != string::npos)
			cout << hd_msg  << " " <<  detect.name << " has assigned sensitivity: "  << sensi << endl;
		
		// The key for the SD, Region, PC maps is the same so we can only check SD
		map<string, sensitiveDetector*>::iterator itr = SeDe_Map.find(sensi);
		// if not found, add new sensitive detector
		if(itr == SeDe_Map.end())
		{
			if(VERB > 3 || detect.name.find(catch_v) != string::npos)
				cout << endl << hd_msg  << " Sensitive detector <" << sensi
				     << "> doesn't exist yet. Adding <" << sensi << ">. " << endl;
			
			// passing detector infos to access factory, runMin, runMax and variation
			SeDe_Map[sensi] = new sensitiveDetector(sensi, gemcOpt, detect.factory, detect.run, detect.variation, detect.system);
			
			// Pass Detector Map Pointer to Sensitive Detector
			SeDe_Map[sensi]->hallMap        = hallMap;
			
			SDman->AddNewDetector( SeDe_Map[sensi]);
		}
		detect.setSensitivity(SeDe_Map[sensi]);
		
		// Setting Max Acceptable Step for this SD
		detect.SetUserLimits(new G4UserLimits(SeDe_Map[sensi]->SDID.maxStep, SeDe_Map[sensi]->SDID.maxStep));
	}
}


void MDetectorConstruction::hasMagfield(detector detect)
{
	string hd_msg    = gemcOpt.optMap["LOG_MSG"].args + " Magnetic Field: >> ";
	double verbosity = gemcOpt.optMap["FIELD_VERBOSITY"].arg ;
	string catch_v   = gemcOpt.optMap["CATCH"].args;
	string field     = gemcOpt.optMap["NO_FIELD"].args;
	
	if(field == "all" || detect.name.find(field) != string::npos)
		return;
	
	string magf   = detect.magfield;
	
	if(magf != "no")
	{
		map<string, gfield>::iterator itr = fieldsMap->find(magf);
		
		if(itr == fieldsMap->end())
		{
			cout << hd_msg << " Electro-Magnetic Field <" << magf
			<< "> is not defined. Exiting." << endl;
			exit(0);
		}
		
		activeFields.insert(magf);
		detect.AssignMFM(itr->second.get_MFM());
		
		if((verbosity > 1 && verbosity != 99) || detect.name.find(catch_v) != string::npos )
			cout << hd_msg  << " Field <" <<  magf << "> is built and assigned to " << detect.name << "." << endl;
	}
}


void MDetectorConstruction::updateGeometry()
{
	cout << "Updating geometry... " << endl;
	G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

detector MDetectorConstruction::findDetector(string name)
{
	map<string, detector>::iterator it = hallMap->find(name);
	if(it != hallMap->end())
		return it->second;
	
	detector notfound;
	notfound.name   = "notfound";
	
	return notfound;
}

void MDetectorConstruction::buildDetector(string name)
{
	detector kid = findDetector(name);
	detector mom = findDetector(kid.mother);

	if(kid.name != "notfound" || mom.name != "notfound")
	{
		
		// handling replicas
		if(kid.type.find("ReplicaOf:") == 0)
		{
			// get the replicant detector
			stringstream ops;
			string operands(kid.type, 10, kid.type.size());
			ops << operands;
			string repName;
			ops >> repName;
			
			detector rep = findDetector(repName);
			
			if(rep.name != "notfound" )
			{
				// creating the replicas volume
				(*hallMap)[kid.name].create_replicas(gemcOpt, mom.GetLogical(), rep);
			}
			else
			{
				cout << "   Attention: " << kid.name << " not found. " << endl;
			}
		}
		else
		{
			(*hallMap)[kid.name].create_solid(gemcOpt, hallMap);
			
			// creating logical volume
			if((*hallMap)[kid.name].create_logical_volume(mats, gemcOpt))
			{
				
				// creating physical volume unless it is in the replicant set
				if(replicants.find(kid.name) == replicants.end())
					(*hallMap)[kid.name].create_physical_volumes(gemcOpt, mom.GetLogical());
				
				
				isSensitive((*hallMap)[kid.name]);  // if sensitive, associate sensitive detector
				hasMagfield((*hallMap)[kid.name]);  // if it has magnetic field, create or use MFM
			}
			(*hallMap)[kid.name].scanned = 1;
		}
	}
	
	if(kid.name == "notfound")
		cout << "   Attention: " << kid.name << " not found. " << endl;
	if(mom.name == "notfound" )
		cout << "   Attention: " << mom.name << " not found. " << endl;
}


void MDetectorConstruction::buildMirrors()
{
	string hd_msg  = gemcOpt.optMap["LOG_MSG"].args + " Mirrors: >> ";
	double VERB    = gemcOpt.optMap["MIRROR_VERBOSITY"].arg ;
	string catch_v = gemcOpt.optMap["CATCH"].args;
	
	vector<G4OpticalSurface*>          mirrorSurfaces;
	vector<G4MaterialPropertiesTable*> mirrorsMPT;
	vector<G4LogicalBorderSurface*>    mirrorLBorderSurf;
	vector<G4LogicalSkinSurface*>      mirrorLSkinSurf;

	for( map<string, detector>::iterator it =  hallMap->begin(); it != hallMap->end(); it++)
	{
		string name = it->first;

		string mirrorString = "no";
		if(it->second.sensitivity.find("mirror:") != string::npos)
		{
			vector<string> mmirror = get_strings(it->second.sensitivity);
			if(mmirror.size() == 2)
				mirrorString = mmirror[1];
		}
		
		if(mirrorString != "no")
		{
			if(VERB > 5 || name.find(catch_v) != string::npos)
				cout << hd_msg  << " " <<  name << " has assigned mirror: "  << mirrorString << endl;
			
			map<string, mirror*>::iterator itr = mirs->find(mirrorString);
			
			// must find the mirror in the mirror map
			if(itr != mirs->end())
			{
				if(VERB > 3 || name.find(catch_v) != string::npos)
					cout << endl << hd_msg  << " Mirror <" << mirrorString << "> found for " << name << "." << endl;
				
				// checking that the border volume exist
				// for non SkinSurface types
				string borderv     = itr->second->border;
				string mirrorSurf = "mirrorBorderSurface_for_" + name + "_and_" + borderv;
				
				if(!(*hallMap)[borderv].GetPhysical() && borderv != "SkinSurface")
				{
					cout << hd_msg << " !! Error: border volume " << borderv << " is not found for volume " << name << ". Exiting." << endl;
					exit(0);
				}
				
				else if(borderv == "SkinSurface")
				{
					mirrorSurf = "mirrorSkinSurface_for_" + name;
				}
			
				mirrorSurfaces.push_back(new G4OpticalSurface(mirrorSurf));
				
				// If a material was given for this interface,
				// AND a material properties table has been defined for said mirror,
				// use this instead of a new table!

				string maptOptProps = itr->second->maptOptProps;
				if( maptOptProps != "notDefined" && ( (*mats)[maptOptProps] )->GetMaterialPropertiesTable() )
				{
					mirrorsMPT.push_back( ( (*mats)[maptOptProps] )->GetMaterialPropertiesTable() );
				}
				else
				{
					mirrorsMPT.push_back(new G4MaterialPropertiesTable());
					
					vector<double> photonEnergy      = itr->second->photonEnergy;
					if(photonEnergy.size())
					{
						vector<double> indexOfRefraction = itr->second->indexOfRefraction;
						vector<double> reflectivity      = itr->second->reflectivity;
						vector<double> efficiency        = itr->second->efficiency;
						vector<double> specularlobe      = itr->second->specularlobe;
						vector<double> specularspike     = itr->second->specularspike;
						vector<double> backscatter       = itr->second->backscatter;
						
						// array size must be the same for all
						// properties by construction, assured by the API
						unsigned peneSize = photonEnergy.size();
												
						G4double pene[peneSize];
						G4double var[peneSize];
						
						for(unsigned i=0; i<peneSize; i++)
							pene[i] = photonEnergy[i];
						
						if(indexOfRefraction.size())
						{
							for(unsigned i=0; i<peneSize; i++) var[i] = indexOfRefraction[i];
							mirrorsMPT.back()->AddProperty("RINDEX", pene, var, peneSize);
						}
						if(reflectivity.size())
						{
							for(unsigned i=0; i<peneSize; i++) var[i] = reflectivity[i];
							mirrorsMPT.back()->AddProperty("REFLECTIVITY", pene, var, peneSize);
						}
						if(efficiency.size())
						{
							for(unsigned i=0; i<peneSize; i++) var[i] = efficiency[i];
							mirrorsMPT.back()->AddProperty("EFFICIENCY", pene, var, peneSize);
						}
						if(specularlobe.size())
						{
							for(unsigned i=0; i<peneSize; i++)var[i] = specularlobe[i];
							mirrorsMPT.back()->AddProperty("SPECULARLOBECONSTANT", pene, var, peneSize);
						}
						if(specularspike.size())
						{
							for(unsigned i=0; i<peneSize; i++)var[i] = specularspike[i];
							mirrorsMPT.back()->AddProperty("SPECULARSPIKECONSTANT", pene, var, peneSize);
						}
						if(backscatter.size())
						{
							for(unsigned i=0; i<peneSize; i++)var[i] = backscatter[i];
							mirrorsMPT.back()->AddProperty("BACKSCATTERCONSTANT", pene, var, peneSize);
						}
					}
					else
					{
						cout << " !! Fatal error: no optical property material, and no optical properties for mirror "
						     << itr->second->name << endl;
						exit(0);
					}
				}
				mirrorSurfaces.back()->SetMaterialPropertiesTable(mirrorsMPT.back());
								
				// surface type
				string surfaceType = itr->second->type;
				if(surfaceType == "dielectric_metal")      mirrorSurfaces.back()->SetType(dielectric_metal);
				if(surfaceType == "dielectric_dielectric") mirrorSurfaces.back()->SetType(dielectric_dielectric);
				if(surfaceType == "dielectric_LUT")        mirrorSurfaces.back()->SetType(dielectric_LUT);
			
				// surface finish
				string surfaceFinish = itr->second->finish;
				if(surfaceFinish == "polished")             mirrorSurfaces.back()->SetFinish(polished);              // smooth perfectly polished surface
				if(surfaceFinish == "polishedfrontpainted") mirrorSurfaces.back()->SetFinish(polishedfrontpainted);  // smooth top-layer (front) paint
				if(surfaceFinish == "polishedbackpainted")  mirrorSurfaces.back()->SetFinish(polishedbackpainted);   // same is 'polished' but with a back-paint
				
				if(surfaceFinish == "ground")               mirrorSurfaces.back()->SetFinish(ground);                // rough surface
				if(surfaceFinish == "groundfrontpainted")   mirrorSurfaces.back()->SetFinish(groundfrontpainted);    // rough top-layer (front) paint
				if(surfaceFinish == "groundbackpainted")    mirrorSurfaces.back()->SetFinish(groundbackpainted);     // same as 'ground' but with a back-paint
				
				
				if(surfaceFinish == "polishedlumirrorair")  mirrorSurfaces.back()->SetFinish(polishedlumirrorair);   // mechanically polished surface, with lumirror
				if(surfaceFinish == "polishedlumirrorglue") mirrorSurfaces.back()->SetFinish(polishedlumirrorglue);  // mechanically polished surface, with lumirror & meltmount

				
				// surface model
				string surfaceModel = itr->second->model;
				if(surfaceModel == "glisur")  mirrorSurfaces.back()->SetModel(glisur);   // original GEANT3 model
				if(surfaceModel == "unified") mirrorSurfaces.back()->SetModel(unified);  // UNIFIED model
				if(surfaceModel == "LUT")     mirrorSurfaces.back()->SetModel(LUT);      // Look-Up-Table model

				
				if(borderv=="SkinSurface")
				{
					mirrorLSkinSurf.push_back(new G4LogicalSkinSurface(name,
											  (*hallMap)[name].GetLogical(), mirrorSurfaces.back()));
				}
				else
				{
					
					mirrorLBorderSurf.push_back(new G4LogicalBorderSurface(name,
																 (*hallMap)[borderv].GetPhysical(),
																 (*hallMap)[name].GetPhysical(), mirrorSurfaces.back()));
				}
				
				
				if(VERB > 3 || name.find(catch_v) != string::npos)
				{
					cout << hd_msg  << " " <<  name << " is a mirror:" << endl;
					cout << "                             > Border Volume: "      << borderv << endl;
					cout << "                             > Surface Type: "       << surfaceType << endl;
					cout << "                             > Finish: "             << surfaceFinish << endl;
					cout << "                             > Model: "              << surfaceModel << endl;
					
					// why it's not dumping all properties?
					mirrorsMPT.back()->DumpTable();
				}
			}
			else
			{
				cout << " !! Fatal error: mirror <" << mirrorString << "> not found for " << it->second.name << "." << endl;
				exit(0);

			}
		}
	}
}

void MDetectorConstruction::assignRegions(vector<string> volumes)
{
	double VERB    = gemcOpt.optMap["HIT_VERBOSITY"].arg ;

	for(unsigned int i=0; i<volumes.size(); i++)
	{
		// looking in the sensitive detector map for the SD with matching system
		for(map<string, sensitiveDetector*>::iterator itr = SeDe_Map.begin(); itr != SeDe_Map.end(); itr++)
		{
			detector regionDet = findDetector(volumes[i]);

			if(regionDet.system == itr->second->SDID.system)
			{

				// region name is volume + system
				string regionName = volumes[i] + "_" + get_info(regionDet.system, "/").back();
				
				
				map<string, G4Region*>::iterator itrr = SeRe_Map.find(regionName);
				
				// if not found, add new region. Otherwise it's already there.
				if(itrr == SeRe_Map.end())
				{
					// Creating G4 Region, assigning Production Cut to it.
					// assigning the logical volume to the region (this will apply the region to all daughters)
					SeRe_Map[regionName] = new G4Region(regionName);
					SeRe_Map[regionName]->AddRootLogicalVolume(regionDet.GetLogical());

					SePC_Map[regionName] = new G4ProductionCuts;
					SePC_Map[regionName] ->SetProductionCut(itr->second->SDID.prodThreshold);
					
					if(VERB > 3)
						cout << "  Region " << regionName << " activated for volume " << regionDet.name << " with range: " << itr->second->SDID.prodThreshold << endl;
					
					SeRe_Map[regionName]->SetProductionCuts(SePC_Map[regionName]);
					
				}
			}
		}		
	}
}










