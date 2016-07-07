// gemc headers
#include "material_factory.h"
#include "cpp_materials.h"

// G4 headers
#include "G4Element.hh"
#include "G4NistManager.hh"
#include "G4OpBoundaryProcess.hh"


// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;


map<string, G4Material*> cpp_materials::initMaterials(runConditions rc, goptions opts)
{
	
	// Many of these elements are already in the G4 Material Database as materials:
	// http://geant4.cern.ch/UserDocumentation/UsersGuides/ForApplicationDeveloper/html/apas08.html
	// When building materials using fractional masses, one can use the G4 Materials Database
	//
	// Elements are still needed when materials are formed using
	// molecular components
	// Example for Water:
	// H2O->AddElement(elH, natoms=2);
	// H2O->AddElement(elO, natoms=1);
	
	
	map<string, G4Material*> MMats;
	
	
	G4double a, z, density;
	G4int nel;
	G4String symbol;

	
	G4Element* Ar  = new G4Element("Argon",     "Ar", z=18, a=    39.95*g/mole);
	G4Element* Al  = new G4Element("Aluminum",  "Al", z=13, a=   26.982*g/mole);
	G4Element* C   = new G4Element("Carbon",    "C",  z=6,  a=    12.01*g/mole);
	G4Element* F   = new G4Element("Fluorine",  "F",  z=9,  a=  18.9984*g/mole);
	G4Element* H   = new G4Element("Hydrogen",  "H",  z=1,  a=     1.01*g/mole);
	G4Element* O   = new G4Element("Oxygen",    "O",  z=8,  a=    16.00*g/mole);
	G4Element* Si  = new G4Element("Silicon",   "Si", z=14, a=    28.09*g/mole);

	G4NistManager* matman = G4NistManager::Instance();
	

	G4Material *Air_Opt = new G4Material("Air_Opt",   density= 1.29*mg/cm3, nel=2);
	Air_Opt->AddMaterial(matman->FindOrBuildMaterial("G4_N"), 70.*perCent);
	Air_Opt->AddMaterial(matman->FindOrBuildMaterial("G4_O"), 30.*perCent);
	
	// this material will kill every tracks that touch it
	G4Material *Kryptonite = new G4Material("Kryptonite", density= 0.00000001*mg/cm3, nel=1);
	Kryptonite->AddElement(Ar, 100.*perCent);

	// aluminum honeycomb core  (it is actually made of alloy Alu-Alloy 3003 (AlMnCu) - still use by HPS ECAL
	G4Material *AlHoneycomb = new G4Material("AlHoneycomb", z=13, a=  26.982*g/mole, density =  0.13*g/cm3);
	



















	
	



	


	


	
	









	// polyethylene
	
		// Nema G10:
	G4Material* NEMAG10 = new G4Material("NEMAG10", 1.70*g/cm3, nel=4);
	NEMAG10 -> AddElement(Si, 1);
	NEMAG10 -> AddElement(O , 2);
	NEMAG10 -> AddElement(C , 3);
	NEMAG10 -> AddElement(H , 3);
	
	// Ar (70) CO2 (30) STP
	G4double density_Ar = 1.7823*mg/cm3 ;
	G4Material* Argon = new G4Material("Argon"  , density_Ar, nel=1);
	Argon->AddElement(Ar, 1);
	
	G4double density_CO2 = 1.977*mg/cm3;
	G4Material* CO2 = new G4Material("CO2", density_CO2, nel=2);
	CO2->AddElement(C, 1);
	CO2->AddElement(O, 2);
	
	G4double density_ArCO2 = .7*density_Ar + .3*density_CO2;
	G4Material *ArCO2 = new G4Material("GEMgas", density_ArCO2, nel=2);
	ArCO2->AddMaterial(Argon, 0.7*density_Ar/density_ArCO2) ;
	ArCO2->AddMaterial(CO2, 0.3*density_CO2/density_ArCO2) ;
	
	
	//RICH MATERIAL
	G4Material *C6F14 = new G4Material("C6F14", density=1.680*g/cm3, nel=2);
	C6F14->AddElement(C,  6);
	C6F14->AddElement(F, 14);
	
	G4Material *C5F12 = new G4Material("C5F12", density=1.680*g/cm3, nel=2);
	C5F12->AddElement(C ,   5);
	C5F12->AddElement(F ,  12);
	
	G4Material* H2O = new G4Material("H20", density=1.000*g/cm3, nel=2);
	H2O->AddElement(H, 2);
	H2O->AddElement(O, 1);

	G4Material *Quartz = new G4Material("Quartz", density= 4.400*g/cm3, nel=2);
	Quartz->AddMaterial(matman->FindOrBuildMaterial("G4_Si"), 1);
	Quartz->AddElement(O,  2);
	
	G4Material *Methane = new G4Material("Methane", density= 0.667*kg/m3, nel=2);    //gas at 20 degree, 1atm
	Methane->AddElement(C, 1);
	Methane->AddElement(H, 4);
	
	G4Material *Alumi = new G4Material("Alumi", density=2.7*g/cm3, nel=1);
	Alumi->AddElement(Al, 1);
	
	G4Material *Glass = new G4Material("Glass", density=1.032*g/cm3, nel=2);
	Glass->AddMaterial(matman->FindOrBuildMaterial("G4_C"), 91.533*perCent);
	Glass->AddMaterial(matman->FindOrBuildMaterial("G4_H"),  8.467*perCent);
	
	
	// END OF RICH MATERIAL
	
	
	
	// Various Vacuums
	// 1 torr is 1/760 atmospheric pressure, that is 1.29*mg/cm3
	G4Material *vacuum_m9 = new G4Material("vacuum_m9", density= 1.68e-12*mg/cm3, nel=2);
	vacuum_m9->AddMaterial(matman->FindOrBuildMaterial("G4_N"), 70.*perCent);
	vacuum_m9->AddMaterial(matman->FindOrBuildMaterial("G4_O"), 30.*perCent);
	
	G4Material *vacuum_m3 = new G4Material("vacuum_m3", density= 1.68e-6*mg/cm3, nel=2);
	vacuum_m3->AddMaterial(matman->FindOrBuildMaterial("G4_N"), 70.*perCent);
	vacuum_m3->AddMaterial(matman->FindOrBuildMaterial("G4_O"), 30.*perCent);
	
	
	
	
	
	MMats["Air_Opt"]        = Air_Opt;
	MMats["Alumi"]          = Alumi;
	MMats["Aluminum"]       = matman->FindOrBuildMaterial("G4_Al");
	MMats["C5F12"]          = C5F12;
	MMats["C6F14"]          = C6F14;
	MMats["Concrete"]       = matman->FindOrBuildMaterial("G4_CONCRETE");       // concrete
	MMats["Copper"]         = matman->FindOrBuildMaterial("G4_Cu");
	MMats["Glass"]          = Glass;
	MMats["Gold"]           = matman->FindOrBuildMaterial("G4_Au");
	MMats["H2O"]            = H2O;
	MMats["Iron"]           = matman->FindOrBuildMaterial("G4_Fe");
	MMats["Kryptonite"]     = matman->FindOrBuildMaterial("Kryptonite");
	MMats["Lead"]           = matman->FindOrBuildMaterial("G4_Pb");
	MMats["LH2"]            = matman->FindOrBuildMaterial("G4_lH2");;
	MMats["Methane"]        = Methane;
	MMats["Nickel"]         =  matman->FindOrBuildMaterial("G4_Ni");
	MMats["Quartz"]         = Quartz;
	MMats["Silicium"]       = matman->FindOrBuildMaterial("G4_Si");
	MMats["Silicon"]        = matman->FindOrBuildMaterial("G4_Si");
	MMats["Teflon"]         = matman->FindOrBuildMaterial("G4_TEFLON");
	MMats["Tungsten"]       = matman->FindOrBuildMaterial("G4_W");
	MMats["Vacuum"]         = matman->FindOrBuildMaterial("G4_Galactic");
	MMats["vacuum_m9"]      = vacuum_m9;
	MMats["vacuum_m3"]      = vacuum_m3;
	MMats["AlHoneycomb"]    = AlHoneycomb;
	MMats["NEMAG10"]        = NEMAG10;
	MMats["Argon"]          = Argon;
	MMats["CO2"]            = CO2;
	MMats["ArCO2"]          = ArCO2;

	
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// Materials Optical Properties
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%
	



	
	// Air Reflection
	const G4int nEntries_Air = 2;
	G4double PhotonEnergy_Air[nEntries_Air]    = { 2.034*eV , 4.136*eV };
	G4double RefractiveIndex_Air[nEntries_Air] = { 1.00, 1.00 };
	
	G4MaterialPropertiesTable* Air_MPT = new G4MaterialPropertiesTable();
	Air_MPT->AddProperty("RINDEX", PhotonEnergy_Air, RefractiveIndex_Air, nEntries_Air);
	MMats["Air_Opt"]->SetMaterialPropertiesTable(Air_MPT);


	
	

	
	// Aerogel properties
	const G4int nEntries_Aerogel = 2;
	G4double PhotonEnergy_Aerogel[nEntries_Aerogel]    = { 2.034*eV, 4.136*eV };
	G4double RefractiveIndex_Aerogel[nEntries_Aerogel] = { 1.03  , 1.03   };
	G4double Absorption_Aerogel[nEntries_Aerogel]      = { 10*m ,  3*cm };
	
	G4MaterialPropertiesTable* Aerogel_MPT = new G4MaterialPropertiesTable();
	Aerogel_MPT->AddProperty("RINDEX",        PhotonEnergy_Aerogel, RefractiveIndex_Aerogel, nEntries_Aerogel);
	Aerogel_MPT->AddProperty("ABSLENGTH",     PhotonEnergy_Aerogel, Absorption_Aerogel,      nEntries_Aerogel);
	Aerogel_MPT->AddConstProperty("SCINTILLATIONYIELD", 10./MeV);
	Aerogel_MPT->AddConstProperty("RESOLUTIONSCALE",    1.0);
	Aerogel_MPT->AddConstProperty("YIELDRATIO",         0.8);
	MMats["Aerogel"]->SetMaterialPropertiesTable(Aerogel_MPT);
	
	
	// RICH
	const G4int nEntries_Rich = 14;
	
	G4double PhotonEnergy_Rich[nEntries_Rich] =
	{ 2.21*eV, 2.30*eV, 2.38*eV, 2.48*eV,
		2.58*eV, 2.70*eV, 2.82*eV, 2.95*eV,
		3.10*eV, 3.26*eV, 3.44*eV, 3.65*eV,
		3.88*eV, 4.13*eV};
	
	
	G4double C5F12_Rind[nEntries_Rich] =
	{ 1.23862, 1.23884, 1.23906, 1.23933,
		1.23962, 1.23998, 1.24035, 1.24078,
		1.24130, 1.24189, 1.24259, 1.24346,
		1.24448, 1.24567};
	
	G4double C5F12_Abs[nEntries_Rich] =
	{ 2000.*mm, 2000.*mm, 2000.*mm, 2000.*mm,
		2000.*mm, 2000.*mm, 2000.*mm, 2000.*mm,
		2000.*mm, 2000.*mm, 2000.*mm, 2000.*mm,
		2000.*mm, 2000.*mm};
	
	
	G4MaterialPropertiesTable* C5F12_MPT = new G4MaterialPropertiesTable();
	C5F12_MPT->AddProperty("RINDEX",     PhotonEnergy_Rich, C5F12_Rind, nEntries_Rich);
	C5F12_MPT->AddProperty("ABSLENGTH",  PhotonEnergy_Rich, C5F12_Abs,  nEntries_Rich);
	MMats["C5F12"]->SetMaterialPropertiesTable(C5F12_MPT);
	

	G4double C6F14_Rind[nEntries_Rich] =
	{ 1.21501, 1.21656, 1.21794, 1.21966,
		1.22138, 1.22344, 1.22550, 1.22774,
		1.23032, 1.23307, 1.23617, 1.23978,
		1.24374, 1.24804};
	
	G4double C6F14_Abs[nEntries_Rich] =
	{ 2000.*mm, 2000.*mm, 2000.*mm, 2000.*mm,
		2000.*mm, 2000.*mm, 2000.*mm, 2000.*mm,
		2000.*mm, 2000.*mm, 2000.*mm, 2000.*mm,
		2000.*mm, 2000.*mm};
	
	G4MaterialPropertiesTable* C6F14_MPT = new G4MaterialPropertiesTable();
	C6F14_MPT->AddProperty("RINDEX",     PhotonEnergy_Rich, C6F14_Rind, nEntries_Rich);
	C6F14_MPT->AddProperty("ABSLENGTH",  PhotonEnergy_Rich, C6F14_Abs,  nEntries_Rich);
	MMats["C6F14"]->SetMaterialPropertiesTable(C6F14_MPT);
	
	
	G4double Quartz_Rind[nEntries_Rich] =
	{ 1.505, 1.509, 1.511, 1.515,
		1.520, 1.525, 1.528, 1.527,
		1.522, 1.512, 1.505, 1.492,
		1.471, 1.503};
	
	G4double Quartz_Abs[nEntries_Rich] =
	{ 550.7*mm, 530.7*mm, 590.1*mm, 490.7*mm,
		470.7*mm, 520.3*mm, 500.0*mm, 470.7*mm,
		450.5*mm, 270.5*mm, 190.1*mm,  60.9*mm,
		10.6*mm,   4.0*mm};
	
	G4MaterialPropertiesTable* Quartz_MPT = new G4MaterialPropertiesTable();
	Quartz_MPT->AddProperty("RINDEX",     PhotonEnergy_Rich, Quartz_Rind, nEntries_Rich);
	Quartz_MPT->AddProperty("ABSLENGTH",  PhotonEnergy_Rich, Quartz_Abs,  nEntries_Rich);
	MMats["Quartz"]->SetMaterialPropertiesTable(Quartz_MPT);
	
	G4double Methane_Rind[nEntries_Rich] =
	{ 1., 1., 1., 1., 1., 1., 1.,
		1., 1., 1., 1., 1., 1., 1.};
	
	G4double Methane_Abs[nEntries_Rich] =
	{4000.*cm,4000.*cm,4000.*cm, 4000.*cm,
		4000.*cm,4000.*cm,4000.*cm, 4000.*cm,
		4000.*cm,4000.*cm};
	
	
	G4MaterialPropertiesTable* Methane_MPT = new G4MaterialPropertiesTable();
	Methane_MPT->AddProperty("RINDEX",     PhotonEnergy_Rich, Methane_Rind, nEntries_Rich);
	Methane_MPT->AddProperty("ABSLENGTH",  PhotonEnergy_Rich, Methane_Abs,  nEntries_Rich);
	MMats["Methane"]->SetMaterialPropertiesTable(Methane_MPT);
	
	
	G4double Alumi_Rind[nEntries_Rich] =
	{ 1., 1., 1., 1., 1., 1., 1.,
		1., 1., 1., 1., 1., 1., 1.};
	
	G4double Alumi_Abs[nEntries_Rich] =
	{0., 0., 0., 0., 0., 0., 0.,
		0., 0., 0., 0., 0., 0., 0.};
	
	G4double Alumi_Effi[nEntries_Rich] =
	{0., 1., 1., 1., 1., 1., 1.,
		1., 1., 1., 1., 1., 1., 0.};
	
	G4double Alumi_Refl[nEntries_Rich] =
	{0., 0., 0., 0., 0., 0., 0.,
		0., 0., 0., 0., 0., 0., 0.};
	
	G4MaterialPropertiesTable* Alumi_MPT = new G4MaterialPropertiesTable();
	Alumi_MPT->AddProperty("RINDEX",       PhotonEnergy_Rich, Alumi_Rind, nEntries_Rich);
	Alumi_MPT->AddProperty("ABSLENGTH",    PhotonEnergy_Rich, Alumi_Abs,  nEntries_Rich);
	Alumi_MPT->AddProperty("EFFICIENCY",   PhotonEnergy_Rich, Alumi_Effi, nEntries_Rich);
	Alumi_MPT->AddProperty("REFLECTIVITY", PhotonEnergy_Rich, Alumi_Refl, nEntries_Rich);
	MMats["Alumi"]->SetMaterialPropertiesTable(Alumi_MPT);
	
	G4double Glass_Rind[nEntries_Rich] =
	{1.49, 1.49, 1.49, 1.49,
		1.49, 1.49, 1.49, 1.49,
		1.49, 1.49, 1.49, 1.49,
		1.49, 1.49};
	
	G4double Glass_Abs[nEntries_Rich] =
	{4200.0*mm, 4200.0*mm, 4200.0*mm, 4200.0*mm,
		4200.0*mm, 4200.0*mm, 4200.0*mm, 4200.0*mm,
		4200.0*mm, 4200.0*mm, 4200.0*mm, 4200.0*mm,
		4200.0*mm, 4200.0*mm};
	
	
	G4MaterialPropertiesTable* Glass_MPT = new G4MaterialPropertiesTable();
	Glass_MPT->AddProperty("RINDEX",     PhotonEnergy_Rich, Glass_Rind, nEntries_Rich);
	Glass_MPT->AddProperty("ABSLENGTH",  PhotonEnergy_Rich, Glass_Abs,  nEntries_Rich);
	MMats["Glass"]->SetMaterialPropertiesTable(Glass_MPT);
	
	
	return MMats;
}
