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
#include "G4GDMLParser.hh"
#include "G4NistManager.hh"

// cadmesh
#include "CADMesh.hh"

// gemc headers
#include "MDetectorConstruction.h"

// mlibrary
#include "gstring.h"
using namespace gstring;

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
	double geo_verb   = gemcOpt.optMap["GEO_VERBOSITY"].arg ;
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
	(*hallMap)["root"].scanned = 1;

	if(VERB > 3 || catch_v == "root") cout << hd_msg << "    " << (*hallMap)["root"] ;
	
	
	cout << hd_msg <<  " Building Detector from Geometry STL Map..." << endl;
	
	
	// ########################################################################
	// Resetting Detector Map "scanned". Propagating "exist" to all generations
	// ########################################################################
	if(VERB > 2) cout << hd_msg << " Mapping Physical Detector..." << endl << endl;
	
	for(map<string, detector>::iterator i =  hallMap->begin(); i!=hallMap->end(); i++)
	{
		if(VERB > 3) cout << hd_msg << " Native Scanning Detector " << i->first << " - existance: " << i->second.exist << endl;
		
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
	string mom, kid;


	// CAD imports
	scanCadDetectors(VERB, catch_v);
	scanDetectors(VERB, catch_v);

	while (nativeRelativesOfCad.size() > 0 || cadRelativesOfNative.size() > 0 || remainingCad.size() > 0 || remainingNative.size() > 0) {


		if(VERB > 10) {

			cout << " native size " << nativeRelativesOfCad.size() << " cad size  " <<  cadRelativesOfNative.size()
			<<  " remaining cad size "  << remainingCad.size() << " remaining native size "  << remainingNative.size() << endl;

			for(auto &c : nativeRelativesOfCad) {
				cout << " native " << c << endl;
			}

			for(auto &c : cadRelativesOfNative) {
				cout << " cad " << c << endl;
			}

			for(auto &c : remainingCad) {
				cout << " remaining cad " << c << endl;
			}
			for(auto &c : remainingNative) {
				cout << " remaining native  "  << c << endl;
			}
		}

		scanCadDetectors(VERB, catch_v);
		scanDetectors(VERB, catch_v);
	}

	// now build GDML volumes.
	set<string> gdmlAlreadyProcessed;
	for(auto &dd : *hallMap) {

		string filename = dd.second.variation;

		if(dd.second.factory == "gdml" && filename != "gdml") {

			// checking that we didn't already processed this file
			if(gdmlAlreadyProcessed.find(filename) == gdmlAlreadyProcessed.end()) {

				if(VERB > 1)
				 cout << "  > Parsing GDML Physical volumes from " << filename << endl;

				// parsing G4 volumes
				G4GDMLParser *parser = new G4GDMLParser();
				parser->Read(filename, 0);

				G4PhysicalVolumeStore::DeRegister(parser->GetWorldVolume());

				// the volume name has to be "World"
				// its oririn are "root" coordinate
				G4LogicalVolume* gdmlWorld = parser->GetVolume("World");

				// only daughters of World will be a new G4PVPlacement in root
				for(int d=0; d<gdmlWorld->GetNoDaughters (); d++) {

					string thisDetName = gdmlWorld->GetDaughter(d)->GetLogicalVolume()->GetName();
					G4LogicalVolume* thisLogical = parser->GetVolume(thisDetName.c_str());

					thisLogical->SetVisAttributes((*hallMap)[thisDetName].VAtts);

					string materialName = trimSpacesFromString((*hallMap)[thisDetName].material);

					if(mats->find(materialName) == mats->end() ) {
						G4NistManager* matman = G4NistManager::Instance();
						if(matman->FindOrBuildMaterial(materialName)) (*mats)[materialName] = matman->FindOrBuildMaterial(materialName);
					}
					thisLogical->SetMaterial((*mats)[materialName]);

					// assigning the logical volume to the detector
					(*hallMap)[thisDetName].SetLogical(thisLogical);
					isSensitive((*hallMap)[thisDetName]);  // if sensitive, associate sensitive detector


					// has to be in the same scope, otherwise parser loses all the pointers
					(*hallMap)[thisDetName].SetPhysical(new G4PVPlacement(&(*hallMap)[thisDetName].rot,
																		  (*hallMap)[thisDetName].pos,
																		  thisLogical,
																		  thisDetName.c_str(),
																		  (*hallMap)["root"].GetLogical(),
																		  false,
																		  0)
														);
					(*hallMap)[thisDetName].scanned = 1;
					// cout << " why are these different " << thisDetName << " " << (*hallMap)[thisDetName].GetPhysical() << "  difference  " << dd.second.GetPhysical()  << endl;
				}
				delete gdmlWorld;
				gdmlAlreadyProcessed.insert(filename);
			}
		}
	}



	// build mirrors
	buildMirrors();


	// assign regions
	// includes root
	regions.push_back("root");
	assignRegions(regions);
	

	// now output det information if verbosity or catch is given
	for(auto &dd : *hallMap) {
		if(geo_verb > 3 || dd.second.name.find(catch_v) != string::npos)
			cout << dd.second ;

	}

	return (*hallMap)["root"].GetPhysical();
}

#include "G4UserLimits.hh"

void MDetectorConstruction::isSensitive(detector detect)
{
	string hd_msg  = gemcOpt.optMap["LOG_MSG"].args + " Sensitivity: >> ";
	double VERB    = gemcOpt.optMap["HIT_VERBOSITY"].arg ;
	string catch_v = gemcOpt.optMap["CATCH"].args;
	
	string sensi   = detect.sensitivity;

	// mirrors are recorded as flux
	if(sensi.find("mirror:") != string::npos) sensi = "flux";

	if(sensi != "no" )
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
			SeDe_Map[sensi]->hallMap = hallMap;
			
			SDman->AddNewDetector( SeDe_Map[sensi]);
		}
		detect.setSensitivity(SeDe_Map[sensi]);
		
		// Setting Max Acceptable Step for this SD
		detect.SetUserLimits(new G4UserLimits(SeDe_Map[sensi]->SDID.maxStep, SeDe_Map[sensi]->SDID.maxStep));
	}
}

void MDetectorConstruction::buildCADDetector(string dname, string filename, int VERB)
{
	// filename has been already verified to exist?
	if(VERB > 1)
		cout << "  > Parsing CAD volume from " << filename << endl;

	CADMesh * mesh = new CADMesh((char *) filename.c_str());
	mesh->SetScale(mm);
	mesh->SetReverse(false);

	// solid
	G4VSolid *cad_solid = mesh->TessellatedMesh();

	// material
	string materialName = trimSpacesFromString((*hallMap)[dname].material);

	if(mats->find(materialName) == mats->end() ) {
		G4NistManager* matman = G4NistManager::Instance();
		if(matman->FindOrBuildMaterial(materialName)) (*mats)[materialName] = matman->FindOrBuildMaterial(materialName);
	}

	// logical
	G4LogicalVolume *cad_logical = new G4LogicalVolume(cad_solid, (*mats)[materialName],dname, 0, 0, 0);
	cad_logical->SetVisAttributes((*hallMap)[dname].VAtts);

	// assigning the logical volume to the detector
	(*hallMap)[dname].SetLogical(cad_logical);
	isSensitive((*hallMap)[dname]);  // if sensitive, associate sensitive detector


	string mom = (*hallMap)[dname].mother;
	if(hallMap->find(mom) == hallMap->end()) {
		cout << " Error: mom " << mom << " not found for " << dname << endl;
		exit(0);
	}

	detector thisMom =  (*hallMap)[mom];


	// has to be in the same scope, otherwise parser loses all the pointers
	(*hallMap)[dname].SetPhysical(new G4PVPlacement(&(*hallMap)[dname].rot,
														  (*hallMap)[dname].pos,
														  cad_logical,
														  dname.c_str(),
														  thisMom.GetLogical(),
														  false,
														  0)
										);
	(*hallMap)[dname].scanned = 1;

}

void MDetectorConstruction::hasMagfield(detector detect)
{
	string hd_msg    = gemcOpt.optMap["LOG_MSG"].args + " Magnetic Field: >> ";
	double verbosity = gemcOpt.optMap["FIELD_VERBOSITY"].arg ;
	string catch_v   = gemcOpt.optMap["CATCH"].args;
	string field     = gemcOpt.optMap["NO_FIELD"].args;


	//cout << field << " " <<

	if(field == "all" || detect.magfield.find(field) != string::npos)
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
			vector<string> mmirror = getStringVectorFromString(it->second.sensitivity);
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
					cout << hd_msg << " !! Error: border volume >" << borderv << "< is not found for volume " << name << ". Exiting." << endl;
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

				// sigma alpha
				// Looking at the code processes/optical/src/G4OpBoundaryProcess.cc it
				// looks like sigmaAlpha only works with opticalSurfaces, not skin surfaces.
				double sigmaAlpha = itr->second->sigmaAlpha;
				if(sigmaAlpha != -1) {
					mirrorSurfaces.back()->SetSigmaAlpha(sigmaAlpha);
				}

				// documentation on sigma alhpa:
				// http://geant4.slac.stanford.edu/UsersWorkshop/PDF/Peter/OpticalPhoton.pdf
				// http://hypernews.slac.stanford.edu/HyperNews/geant4/get/opticalphotons/397.html

				if(borderv=="SkinSurface") {
					mirrorLSkinSurf.push_back(new G4LogicalSkinSurface(name,
											  (*hallMap)[name].GetLogical(), mirrorSurfaces.back()));
				} else {
					mirrorLBorderSurf.push_back(new G4LogicalBorderSurface(name,
																		   (*hallMap)[name].GetPhysical(),
																		   (*hallMap)[borderv].GetPhysical(), mirrorSurfaces.back()));
				}
				

				if(VERB > 3 || name.find(catch_v) != string::npos)
				{
					cout << hd_msg  << " " <<  name << " is a mirror:" << endl;
					cout << "                             > Border Volume: "      << borderv << endl;
					cout << "                             > Surface Type: "       << surfaceType << endl;
					cout << "                             > Finish: "             << surfaceFinish << endl;
					cout << "                             > Model: "              << surfaceModel << endl;
					cout << "                             > Sigma Alpha: "        << sigmaAlpha << endl;

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
	int fastmcMode = gemcOpt.optMap["FASTMCMODE"].arg;

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

					// protecting against a too low production threshold
					double productionThreshold = itr->second->SDID.prodThreshold;
					if(fastmcMode > 0) productionThreshold = 5000;

					if(productionThreshold < 0.00001) {
						cout << " !! Warning: production threshold for " << regionName << " is  " << productionThreshold << "mm."
						<< " That is too low. Overwriting it with 1mm." << endl;
						productionThreshold = 1;
						itr->second->SDID.prodThreshold = productionThreshold;
					}

					SePC_Map[regionName] ->SetProductionCut(productionThreshold);


					if(VERB > 3)
						cout << "  Region " << regionName << " activated for volume " << regionDet.name << " with range: " << itr->second->SDID.prodThreshold << endl;
					
					SeRe_Map[regionName]->SetProductionCuts(SePC_Map[regionName]);
					
				}
			}
		}		
	}
}




void MDetectorConstruction::scanDetectors(int VERB, string catch_v)
{
	string hd_msg = "  Native Scanning: ";
	vector<string> relatives;

	for( map<string, detector>::iterator i =  hallMap->begin(); i!=hallMap->end(); i++)
	{
		// don't build anything if the exist flag is not set
		if(i->second.exist == 0 || i->second.scanned == 1 || i->second.factory != "TEXT") continue;

		// put the volume in relatives to fill it
		// if everything is good, it will be built right away
		if(i->first != "root") relatives.push_back(i->second.name);

		while(relatives.size() > 0)
		{
			detector kid = findDetector(relatives.back());
			detector mom = findDetector(kid.mother);
			// cout << " ASD " << kid.name << " " << kid.mother <<  " " << kid.scanned << " " << mom.scanned << " " << mom.factory << endl;

			// production cut affects all volumes in a system rather than just the sensitive volumes
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
			if(kid.scanned == 0 && mom.scanned == 1) {
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
						if(dorig.factory == "TEXT") {
							relatives.push_back(original);
						} else {
							// skipping this until later
							cadRelativesOfNative.push_back(original);
							remainingNative.push_back(kid.name);
							relatives.pop_back();
						}
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
				} else if(kid.type.find("ReplicaOf:") == 0) // Check kid dependency on replicas
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
						relatives.push_back(trimSpacesFromString(original));

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
				} else if(kid.type.find("Operation:") == 0)  // Check kid dependency on operations
				{
					int sstart = 10;
					if(kid.type.find("Operation:~") == 0 || kid.type.find("Operation:@") == 0 ) sstart = 11;

					string operation(kid.type, sstart, kid.type.size());
					vector<string> operands = getStringVectorFromString(replaceCharInStringWithChars(operation, "-+*/", " "));
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
				else {
					buildDetector(kid.name);
					kid.scanned = 1;
				}
			}

			// if the kid still doesn't exists and its mom doesn't exist.
			// adding mom to the relatives list
			else if(kid.scanned == 0 && mom.scanned == 0)
			{
				if(mom.factory == "TEXT") {
					// we can still build this unless the mother is inside remainingNative
					if(find(remainingNative.begin(), remainingNative.end(), kid.mother) == remainingNative.end()) {
						relatives.push_back(kid.mother);
					} else {
						relatives.pop_back();
						remainingNative.push_back(kid.name);
					}
				} else {
					// need to build this later, popping this from the local list
					cadRelativesOfNative.push_back(kid.mother);
					remainingNative.push_back(kid.name);
					relatives.pop_back();

				}
			}

			// the kid has been built. Can go down one step in geneaology
			if(kid.scanned == 1) {
				if(VERB > 3 || kid.name.find(catch_v) != string::npos)
					cout << hd_msg  << " " <<  kid.name << " is built." <<  endl << endl;

				if(relatives.size()) {
					relatives.pop_back();
				}

				// if this volume was a parent of a CAD volume, removing it from the list
				auto nativeRelativesOfCadIT = find(nativeRelativesOfCad.begin(), nativeRelativesOfCad.end(), kid.name);
				while(nativeRelativesOfCadIT != nativeRelativesOfCad.end()) {
					nativeRelativesOfCad.erase(nativeRelativesOfCadIT);
					nativeRelativesOfCadIT = find(nativeRelativesOfCad.begin(), nativeRelativesOfCad.end(), kid.name);
				}

				// if this volume was on hold but now it is build, removing it from the remainingNative
				auto remainingNativeIT = find(remainingNative.begin(), remainingNative.end(), kid.name);
				while(remainingNativeIT != remainingNative.end()) {
					remainingNative.erase(remainingNativeIT);
					remainingNativeIT = find(remainingNative.begin(), remainingNative.end(), kid.name);
				}

			}


		}
	}
}




void MDetectorConstruction::scanCadDetectors(int VERB, string catch_v)
{
	string hd_msg = "  CAD Scanning: ";
	vector<string> relatives;
	// building these first in case we want to make copies of these
	for(auto &dd : *hallMap) {

		// don't build anything if the exist flag is not set
		if(dd.second.exist == 0 || dd.second.scanned == 1 || dd.second.factory != "CAD" || dd.first == "root" ) continue;


		string thisDetName = dd.first;


		// put the volume in relatives to fill it
		// if everything is good, it will be built right away
		if(thisDetName != "root") relatives.push_back(dd.second.name);

		while(relatives.size() > 0)
		{
			detector kid = findDetector(relatives.back());
			detector mom = findDetector(kid.mother);
			// cout << " ASD " << kid.name << " " << kid.mother <<  " " << kid.scanned << " " << mom.scanned << " " << mom.factory << endl;

			// Mom doesn't exists in the hallMap. Stopping everything.
			if(mom.name != "akasha"  && mom.name == "notfound") {
				cout << hd_msg << "  Mom was not found for <" << relatives.back() << ">. "
				<< " We have a No Child Left Behind policy. Exiting. " << endl << endl;
				exit(0);
			}

			// Mom is built, kid not built yet. Build kid
			if(kid.scanned == 0 && mom.scanned == 1) {
				string filename = dd.second.variation;
				buildCADDetector(thisDetName, filename, VERB);
				kid.scanned = 1;

			} else if(kid.scanned == 0 && mom.scanned == 0 ) {
				if(mom.factory == "CAD") {
					// we can still build this unless the mother is inside remainingCad
					if(find(remainingCad.begin(), remainingCad.end(), kid.mother) == remainingCad.end()) {
						relatives.push_back(kid.mother);
					} else {
						relatives.pop_back();
						remainingCad.push_back(kid.name);
				}

				} else {
					// need to build this later, popping this from the local list
					nativeRelativesOfCad.push_back(kid.mother);
					remainingCad.push_back(kid.name);
					relatives.pop_back();
				}
				// the kid has been built. Can go down one step in geneaology
			} else if(kid.scanned == 1 && relatives.size()) {

				if(VERB > 3 || kid.name.find(catch_v) != string::npos)
					cout << hd_msg  << " " <<  kid.name << " is built." <<  endl << endl;

				if(relatives.back() == kid.name) {
					relatives.pop_back();
				}


				// if this volume was a parent of a native volume, removing it from the list
				auto cadRelativesOfNativeIT = find(cadRelativesOfNative.begin(), cadRelativesOfNative.end(), kid.name);
				while(cadRelativesOfNativeIT != cadRelativesOfNative.end()) {
					cadRelativesOfNative.erase(cadRelativesOfNativeIT);
					cadRelativesOfNativeIT = find(cadRelativesOfNative.begin(), cadRelativesOfNative.end(), kid.name);
				}

				// if this volume was on hold but now it is build, removing it from the remainingCad
				auto remainingCadIT = find(remainingCad.begin(), remainingCad.end(), kid.name);
				while(remainingCadIT != remainingCad.end()) {
					remainingCad.erase(remainingCadIT);
					remainingCadIT = find(remainingCad.begin(), remainingCad.end(), kid.name);
				}

			}
		}
	}
	
	
	
	
	
}


