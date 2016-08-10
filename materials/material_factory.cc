// c++ headers
#include <set>

// gemc headers
#include "material_factory.h"
#include "cpp_materials.h"
#include "mysql_materials.h"
#include "text_materials.h"
#include "string_utilities.h"

// mlibrary
#include "gstring.h"
using namespace gstring;

// G4 headers
#include "G4NistManager.hh"

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;

materials *getMaterialFactory(map<string, materialFactory> *factory, string materialsMethod)
{
	
	if(factory->find(materialsMethod) == factory->end())
	{
		cout << endl << endl << "  >>> WARNING: " << materialsMethod << " NOT FOUND IN Material Factory Map." << endl;
		return NULL;
	}
	
	return (*factory)[materialsMethod]();
}

map<string, materialFactory> registerMaterialFactories()
{
	
	map<string, materialFactory> materialMethodMap;
	
	// CPP initialization
	materialMethodMap["CPP"]   = &cpp_materials::createMaterials;
	
	// MYSQL initialization
	materialMethodMap["MYSQL"] = &mysql_materials::createMaterials;
	
	// TEXT initialization
	materialMethodMap["TEXT"] = &text_materials::createMaterials;
		
	return materialMethodMap;
}

void printMaterials(map<string, G4Material*> matMap)
{
	for(map<string, G4Material*>::iterator it = matMap.begin(); it != matMap.end(); it++)
	{
		cout << "    - material: " << it->first << " " << it->second << endl;
		if(it->second->GetMaterialPropertiesTable())
		{
			cout << "    - Optical Properties for " << it->first << ":" << endl;
			it->second->GetMaterialPropertiesTable()->DumpTable();
			cout << endl << endl;
		}
	}
}

void material::componentsFromString(string s)
{
	stringstream comps(s);
	
	for(int e=0; e<ncomponents; e++)
	{
		string thisComp;
		double thisFrac;
		comps >> thisComp >> thisFrac;
		components.push_back(thisComp);
		fracs.push_back(thisFrac);
	}
}

void material::opticalsFromString(string property, string what)
{
	
	string checkProperty = trimSpacesFromString(property);
	if(checkProperty == "none") return;
	else
	{
		stringstream comps(property);
	
		while(!comps.eof())
		{
			string c;
			comps >> c ;
			string trimmedC = trimSpacesFromString(c);
			
			if(what == "photonEnergy")
				photonEnergy.push_back(get_number(trimmedC));
			
			if(what == "indexOfRefraction")
				indexOfRefraction.push_back(get_number(trimmedC));

			if(what == "absorptionLength")
				absorptionLength.push_back(get_number(trimmedC));
			
			if(what == "reflectivity")
				reflectivity.push_back(get_number(trimmedC));
			
			if(what == "efficiency")
				efficiency.push_back(get_number(trimmedC));

			// scintillation
			if(what == "fastcomponent")
				fastcomponent.push_back(get_number(trimmedC));
			
			if(what == "slowcomponent")
				slowcomponent.push_back(get_number(trimmedC));
			
			if(what == "scintillationyield")
				scintillationyield = get_number(trimmedC);
			
			if(what == "resolutionscale")
				resolutionscale = get_number(trimmedC);
			
			if(what == "fasttimeconstant")
				fasttimeconstant = get_number(trimmedC);
			
			if(what == "slowtimeconstant")
				slowtimeconstant = get_number(trimmedC);
			
			if(what == "yieldratio")
				yieldratio = get_number(trimmedC);
			
			if(what == "rayleigh")
				rayleigh.push_back(get_number(trimmedC));
		}
	
		if(what == "rayleigh")
		{
			// yieldratio is the last quantity to be loaded
			// now we can check the vector sizes for comparison
			// if no match, resetting quantities
			if(indexOfRefraction.size() != photonEnergy.size()) indexOfRefraction.clear();
			if(absorptionLength.size()  != photonEnergy.size()) absorptionLength.clear();
			if(reflectivity.size()      != photonEnergy.size()) reflectivity.clear();
			if(efficiency.size()        != photonEnergy.size()) efficiency.clear();
			if(fastcomponent.size()     != photonEnergy.size()) fastcomponent.clear();
			if(slowcomponent.size()     != photonEnergy.size()) slowcomponent.clear();
			if(rayleigh.size()          != photonEnergy.size()) rayleigh.clear();
		}
	
	}
}




// We start with the CPP definitions
// Then load the MYSQL, TEXT and GDML
// When all the detectors materials are moved to TEXT/MYSQL/GDML
// the CPP factory should be empty
map<string, G4Material*> buildMaterials(map<string, materialFactory> materialFactoryMap, goptions go, runConditions rc)
{
	// Loading CPP def
	materials *materialSelectedFactory = getMaterialFactory(&materialFactoryMap, "CPP");
	map<string, G4Material*> mats = materialSelectedFactory->initMaterials(rc, go);
		
	// adding MYSQL
	materials *mysqlFactory = getMaterialFactory(&materialFactoryMap, "MYSQL");
	map<string, G4Material*> mysqlMats = mysqlFactory->initMaterials(rc, go);
	for(map<string, G4Material*>::iterator it = mysqlMats.begin(); it != mysqlMats.end(); it++)
		mats[it->first] = it->second;
	
	// adding TEXT
	materials *textFactory = getMaterialFactory(&materialFactoryMap, "TEXT");
	map<string, G4Material*> textMats = textFactory->initMaterials(rc, go);
	for(map<string, G4Material*>::iterator it = textMats.begin(); it != textMats.end(); it++)
		mats[it->first] = it->second;

	return mats;
}

map<string, G4Material*>  materials::materialsFromMap(map<string, material> mmap)
{
	G4NistManager* matman = G4NistManager::Instance();   // material G4 Manager
	set<string> nistMap;
	
	// building STL map of existing material names
	// so I can use "find" later
	vector<G4String> allMats = matman->GetNistMaterialNames();
	for(unsigned int j=0; j<allMats.size(); j++)
		nistMap.insert(allMats[j]);
	
	vector<G4String> allEls = matman->GetNistElementNames();
	for(unsigned int j=0; j<allEls.size(); j++)
	{
		nistMap.insert(allEls[j]);
		//	cout << allEls[j] << " element " << endl;
	}
	map<string, G4Material*> mats = materialsWithIsotopes();
	
	// vector of Optical Properties
	vector<G4MaterialPropertiesTable*> optTable;
	
	
	//// If the materials is molecular constructed each quantity is
	//// an integer >= 1
	//// The G4 elements symbols *do not* have G4 in them
	//// Example:
	//// H2O->AddElement(H, natoms=2);
	//// H2O->AddElement(O, natoms=1);
	//// These elements are pre-lodaded by geant4.
	//// In g4 users can also define their own element - not supported by gemc yet.
	////
	//// If the materials is material-constructed each quantity is
	//// a percentage
	//// The materials could have G4 in their name if they come from the G4 database
	//// Example:
	//// Gas->AddMaterial(matman->FindOrBuildMaterial("G4_H"), 0.2);
	//// all percentages must add to 1.

	// if all the materials / elements exist, adding it to the material list.
	// otherwise continue until the map is empty
	

	// doing maxMapSize max iterations.
	int maxMapSize  = 10000;
	int maxMapIndex = 0;
	
	// removing materials already built from the map
	// can't remove materials while in the map loop, so bookkeeping them
	// in a list to be deleted outside the loop
	vector<string> toDelete;

	while (mmap.size() && maxMapIndex < maxMapSize)
	{
		maxMapIndex++;
		
		// this is the outside loop - deleting all materials added in the inside loop
		// then clearing the map again
		for(unsigned i=0; i<toDelete.size(); i++)
			mmap.erase(toDelete[i]);
		toDelete.clear();
		
		
		for(map<string, material>::iterator it = mmap.begin(); it != mmap.end() && mmap.size(); it++)
		{
			// first check if the components are available
			// if not continue, they may be build later
			bool allExist = 1;
			for(unsigned int i=0; i<it->second.components.size(); i++)
			{
				string compName = it->second.components[i];
				// first check if it's mats
				if(mats.find(compName) == mats.end())
				{
					// not in mats, check if it's the Nist Manager
					if(nistMap.find(compName) == nistMap.end())
					{
						// cout << " Material " << compName << " does not exist yet. Continuing..." <<  endl;
						allExist = 0;
					}
				}
			}
			
			// material elements do not exist yet - continue
			if(!allExist) continue;
			
			else
			// elements exist, build material and remove it from from mmap
			{
				mats[it->first] = new G4Material(it->first, it->second.density*g/cm3, it->second.ncomponents);
				
				// now check if the components are elements or material
				// They will either sum to exactly 1 or > 1
				double totComps = 0;
				for(unsigned int i=0; i<it->second.fracs.size(); i++)
					totComps += it->second.fracs[i];
				
				
				// fractional components - must add to one exactly
				if(fabs(totComps-1) < 0.00001)
				{
					for(unsigned int i=0; i<it->second.fracs.size(); i++)
					{
						string compName = it->second.components[i];
												
						// existing material
						if(mats.find(compName) != mats.end())
							mats[it->first]->AddMaterial(mats[compName], it->second.fracs[i]);
						else
							// G4 Material DB
							mats[it->first]->AddMaterial(matman->FindOrBuildMaterial(compName), it->second.fracs[i]);
					}
				}
				// molecular composition
				else if(totComps > 1)
				{
					for(unsigned int i=0; i<it->second.fracs.size(); i++)
					{
						
						string compName = it->second.components[i];
						if(it->second.fracs[i] < 1)
						{
							cout << " The number of atoms of " << compName << " is " << it->second.fracs[i] << " but it should be an integer. Exiting" << endl;
							exit(0);
						}
						mats[it->first]->AddElement(matman->FindOrBuildElement(compName), (int) it->second.fracs[i]);
					}
				}
				else
				{
					cout << " Warning: the sum of all components for >"  << it->first << "< does not add to one." << endl;
				}
				
				// check for optical properties
				unsigned nopts = it->second.photonEnergy.size();
				if(nopts)
				{
					optTable.push_back(new G4MaterialPropertiesTable());
					
					double penergy[nopts];
					for(unsigned i=0; i<nopts; i++)
						penergy[i] = it->second.photonEnergy[i];

					// index of refraction
					if(it->second.indexOfRefraction.size() == nopts)
					{
						double ior[nopts];
						for(unsigned i=0; i<nopts; i++)
						{
							ior[i] = it->second.indexOfRefraction[i];
						}
						optTable.back()->AddProperty("RINDEX", penergy, ior, nopts );
					}
					
					// absorption length
					if(it->second.absorptionLength.size() == nopts)
					{
						double abs[nopts];
						for(unsigned i=0; i<nopts; i++)
							abs[i] = it->second.absorptionLength[i];
						optTable.back()->AddProperty("ABSLENGTH", penergy, abs, nopts );
					}
					
					// reflectivity
					if(it->second.reflectivity.size() == nopts)
					{
						double ref[nopts];
						for(unsigned i=0; i<nopts; i++)
							ref[i] = it->second.reflectivity[i];
						
						optTable.back()->AddProperty("REFLECTIVITY", penergy, ref, nopts );
					}
					
					// efficiency
					if(it->second.efficiency.size() == nopts)
					{
						double eff[nopts];
						for(unsigned i=0; i<nopts; i++)
							eff[i] = it->second.efficiency[i];
						
						optTable.back()->AddProperty("EFFICIENCY", penergy, eff, nopts );
					}
					
					// fastcomponent
					if(it->second.fastcomponent.size() == nopts)
					{
						double fastc[nopts];
						for(unsigned i=0; i<nopts; i++)
							fastc[i] = it->second.fastcomponent[i];
						
						optTable.back()->AddProperty("FASTCOMPONENT", penergy, fastc, nopts );
					}
					
					// slowcomponent
					if(it->second.slowcomponent.size() == nopts)
					{
						double slowc[nopts];
						for(unsigned i=0; i<nopts; i++)
							slowc[i] = it->second.slowcomponent[i];
						
						optTable.back()->AddProperty("SLOWCOMPONENT", penergy, slowc, nopts );
					}
					
					// rayleigh scattering
					if(it->second.rayleigh.size() == nopts)
					{
						double ray[nopts];
						for(unsigned i=0; i<nopts; i++)
							ray[i] = it->second.rayleigh[i];
						
						optTable.back()->AddProperty("RAYLEIGH", penergy, ray, nopts);
					}
					

					
					// scintillationyield
					if(it->second.scintillationyield != -1)
						optTable.back()->AddConstProperty("SCINTILLATIONYIELD", it->second.scintillationyield);
					
					// resolutionscale
					if(it->second.resolutionscale != -1)
						optTable.back()->AddConstProperty("RESOLUTIONSCALE", it->second.resolutionscale);
					
					// fasttimeconstant
					if(it->second.fasttimeconstant != -1)
						optTable.back()->AddConstProperty("FASTTIMECONSTANT", it->second.fasttimeconstant*ns);
					
					// slowtimeconstant
					if(it->second.slowtimeconstant != -1)
						optTable.back()->AddConstProperty("SLOWTIMECONSTANT", it->second.slowtimeconstant*ns);
					
					// yieldratio
					if(it->second.yieldratio != -1)
						optTable.back()->AddConstProperty("YIELDRATIO", it->second.yieldratio);

					mats[it->first]->SetMaterialPropertiesTable(optTable.back());
				}
			
				// tagging material to be removed from map
				toDelete.push_back(it->first);
			}
		}
	}
	
	return mats;
}


// build materials with standard isotopes
map<string, G4Material*>  materialsWithIsotopes()
{
	map<string, G4Material*> mats;

	G4MaterialTable* matTable = (G4MaterialTable*) G4Material::GetMaterialTable();
	
	bool already_defined = 0;
	for(unsigned i=0; i<matTable->size(); ++i)
	{
		G4Material* thisMat = (*(matTable))[i];
		
		if(thisMat->GetName() == "LD2")
		{
			// LD2 already defined, so we can add the isotopes in our map here
			mats["LD2"] = thisMat;
			already_defined = 1;
		}
		if(thisMat->GetName() == "helium3Gas")
		{
			// LD2 already defined, so we can add the isotopes in our map here
			mats["helium3Gas"] = thisMat;
		}
		if(thisMat->GetName() == "deuteriumGas")
		{
			// LD2 already defined, so we can add the isotopes in our map here
			mats["deuteriumGas"] = thisMat;
		}
		if(thisMat->GetName() == "ND3")
		{
			// LD2 already defined, so we can add the isotopes in our map here
			mats["ND3"] = thisMat;
		}
	
	}
	
	if(already_defined)
		return mats;
	
	// isotopes not yet defined, defining them for the first time
	
	int Z, N;
	double a;
	
	
	// ----  gas and liquid deuterium

	// Deuteron isotope
	G4Isotope* deuteron  = new G4Isotope("deuteron", Z=1, N=2, a=2.0141018*g/mole);
	
	// Deuterium element
	G4Element* deuterium = new G4Element("deuterium", "deuterium", 1);
	deuterium->AddIsotope(deuteron, 1);
	
	// Deuterium gas
	mats["deuteriumGas"] = new G4Material("deuteriumGas", 0.000452*g/cm3, 1, kStateGas, 294.25*kelvin);
	mats["deuteriumGas"]->AddElement(deuterium, 1);
	
	// Liquid Deuterium
	mats["LD2"] = new G4Material("LD2", 0.169*g/cm3, 1, kStateLiquid, 22.0*kelvin);
	mats["LD2"]->AddElement(deuterium, 2);


	// ---- helium 3 gas
	
	// helion isotope
	G4Isotope* helion  = new G4Isotope("helion", Z=2, N=3, a=3.0160293*g/mole);
	
	// helium 3 element
	G4Element* helium3 = new G4Element("helium3", "helium3", 1);
	helium3->AddIsotope(helion, 1);

	// helium 3 material gas
	// Density at 21.1°C (70°F) : 0.1650 kg/m3
	mats["helium3Gas"] = new G4Material("helium3Gas",  0.1650*mg/cm3, 1, kStateGas, 294.25*kelvin);
	mats["helium3Gas"]->AddElement(helium3, 1);
	
	
	
	// ammonia
	G4Element* Nitro   = new G4Element("Nitro",  "N",  Z=7,  a=14.01*g/mole);
	mats["ND3"] = new G4Material("ND3", 1.007*g/cm3, 2, kStateLiquid, 1.0*kelvin);
	mats["ND3"]->AddElement(Nitro, 1);
	mats["ND3"]->AddElement(deuterium, 3);
	

	
	return mats;
}


