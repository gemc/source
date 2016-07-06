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
	G4Element* Ba  = new G4Element("Barium",    "Ba", z=56, a=  137.327*g/mole);
	G4Element* C   = new G4Element("Carbon",    "C",  z=6,  a=    12.01*g/mole);
	G4Element* Cr  = new G4Element("Chromium",  "Cr", z=24, a=    52.00*g/mole);
	G4Element* F   = new G4Element("Fluorine",  "F",  z=9,  a=  18.9984*g/mole);
	G4Element* H   = new G4Element("Hydrogen",  "H",  z=1,  a=     1.01*g/mole);
	G4Element* He  = new G4Element("Helium",    "He", z=2,  a=   4.0026*g/mole);
	G4Element* He3 = new G4Element("Helium3",  "He3", z=2,  a=3.0160293*g/mole);
	G4Element* Mn  = new G4Element("Manganese", "Mn", z=25, a=    54.94*g/mole);
	G4Element* N   = new G4Element("Nitrogen",  "N",  z=7,  a=    14.01*g/mole);
	G4Element* O   = new G4Element("Oxygen",    "O",  z=8,  a=    16.00*g/mole);
	G4Element* Ox  = new G4Element("Ox",        "Ox", z=8,  a=     16.0*g/mole);
	G4Element* Ni  = new G4Element("Nickel",    "Ni", z=28, a=    58.70*g/mole);
	G4Element* Pb  = new G4Element("Lead",      "Pb", z=82, a=   207.19*g/mole);
	G4Element* Si  = new G4Element("Silicon",   "Si", z=14, a=    28.09*g/mole);
	G4Element* Sr  = new G4Element("Strontium", "Sr", z=38, a=    87.62*g/mole);
	G4Element* Wf  = new G4Element("Wf",        "Wf", z=74, a=   183.85*g/mole);
	
	G4NistManager* matman = G4NistManager::Instance();
	

	G4Material *Air_Opt = new G4Material("Air_Opt",   density= 1.29*mg/cm3, nel=2);
	Air_Opt->AddMaterial(matman->FindOrBuildMaterial("G4_N"), 70.*perCent);
	Air_Opt->AddMaterial(matman->FindOrBuildMaterial("G4_O"), 30.*perCent);
	
	// this material will kill every tracks that touch it
	G4Material *Kryptonite = new G4Material("Kryptonite", density= 0.00000001*mg/cm3, nel=1);
	Kryptonite->AddElement(Ar, 100.*perCent);

	// aluminum honeycomb core  (it is actually made of alloy Alu-Alloy 3003 (AlMnCu) - still use by HPS ECAL
	G4Material *AlHoneycomb = new G4Material("AlHoneycomb", z=13, a=  26.982*g/mole, density =  0.13*g/cm3);
	



















	
	


	G4Material *DCgas    = new G4Material("DCgas",              density = 1.8*mg/cm3, nel=3);
	DCgas->AddElement(Ar, 90*perCent);
	DCgas->AddMaterial(matman->FindOrBuildMaterial("G4_O"),  6.6*perCent);
	DCgas->AddMaterial(matman->FindOrBuildMaterial("G4_C"),  3.4*perCent);
	
	
	
	G4Material *FTinsfoam    = new G4Material("FTinsfoam",     density = 34*kg/m3, nel=4);
	FTinsfoam->AddMaterial(matman->FindOrBuildMaterial("G4_C"), 60.*perCent);
	FTinsfoam->AddMaterial(matman->FindOrBuildMaterial("G4_H"), 10.*perCent);
	FTinsfoam->AddMaterial(matman->FindOrBuildMaterial("G4_N"), 10.*perCent);
	FTinsfoam->AddMaterial(matman->FindOrBuildMaterial("G4_O"), 20.*perCent);
	
	G4Material *Noryl	= new G4Material("Noryl",	    density = 1.06*g/cm3, nel=3);
	Noryl->AddMaterial(matman->FindOrBuildMaterial("G4_C"), 47.06*perCent );
	Noryl->AddMaterial(matman->FindOrBuildMaterial("G4_H"), 47.06*perCent );
	Noryl->AddMaterial(matman->FindOrBuildMaterial("G4_O"),  5.88*perCent );
	
	
	G4Material *svtwirebond    = new G4Material("svtwirebond",  density = 2.69*g/cm3, nel=2);
	svtwirebond->AddElement(Al, 99*perCent);
	svtwirebond->AddMaterial(matman->FindOrBuildMaterial("G4_Si"), 1*perCent);
	
	G4Material *MMGas = new G4Material("MMGas",   density = (1.662*0.95+2.489*0.05)*mg/cm3, nel=3);
	MMGas->AddElement(Ar, 95.0*perCent);
	MMGas->AddMaterial(matman->FindOrBuildMaterial("G4_H"), 0.173414*5.0*perCent);
	MMGas->AddMaterial(matman->FindOrBuildMaterial("G4_C"), 0.826586*5.0*perCent);
	
	
	
	
	/*
	 G4Material *MMGas = new G4Material("MMGas",   density = 1.17*mg/cm3, nel=4);
	 MMGas->AddElement(Ne, 79.0*perCent);
	 MMGas->AddMaterial(matman->FindOrBuildMaterial("G4_H"), 0.2011*11.0*perCent);
	 MMGas->AddElement(C, (0.7989*11.0+0.1365*10.0)*perCent);
	 MMGas->AddElement(F, 0.8635*10.0*perCent);
	 */
	/*
	 G4Material *MMGas = new G4Material("MMGas",   density = 1.87*mg/cm3, nel=4);
	 MMGas->AddElement(Ar, 95.0*perCent);
	 MMGas->AddMaterial(matman->FindOrBuildMaterial("G4_H"), 0.1734*2.0*perCent);
	 MMGas->AddElement(C, (0.8266*2.0+0.1365*3.0)*perCent);
	 MMGas->AddElement(F, 0.8635*3.0*perCent);
	 */
	
	
	//G4double MMMeshTransparency = (19./50.)*(19./50.);
	G4double MMMeshTransparency = 1.0;
	
	G4Material *MMMesh = new G4Material("MMMesh",   density = 8.02*MMMeshTransparency*g/cm3, nel=5);
	MMMesh->AddElement(Mn, 0.02);
	MMMesh->AddMaterial(matman->FindOrBuildMaterial("G4_Si"), 0.01);
	MMMesh->AddElement(Cr, 0.19);
	MMMesh->AddElement(Ni, 0.10);
	MMMesh->AddMaterial(matman->FindOrBuildMaterial("G4_Fe"), 0.68);
	
	/*
	 G4Material *MMMesh       = new G4Material("MMMesh",  z=28, a=   58.70*g/mole,  density =  8.902*MMMeshTransparency*g/cm3);
	 */
	
	G4Material *MMMylar = new G4Material("MMMylar",   density = 1.40*g/cm3, nel=3);
	MMMylar->AddMaterial(matman->FindOrBuildMaterial("G4_H"), 0.041958);
	MMMylar->AddMaterial(matman->FindOrBuildMaterial("G4_C"), 0.625017);
	MMMylar->AddMaterial(matman->FindOrBuildMaterial("G4_O"), 0.333025);
	
	/*
	 G4Material *MMMylar = new G4Material("MMMylar",   density = 8.02*MMMeshTransparency*g/cm3, nel=5);
	 MMMylar->AddElement(Mn, 0.02);
	 MMMylar->AddMaterial(matman->FindOrBuildMaterial("G4_Si"), 0.01);
	 MMMylar->AddElement(Cr, 0.19);
	 MMMylar->AddElement(Ni, 0.10);
	 MMMylar->AddMaterial(matman->FindOrBuildMaterial("G4_Fe"), 0.68);
	 */
	//G4Material *MMMylar       = new G4Material("MMMylar",  z=28, a=   58.70*g/mole,  density =  8.902*g/cm3);
	G4Material* He4_1atm = new G4Material( "He4_1atm",  density = 1.*0.1786*mg/cm3, nel=1 );
	He4_1atm->AddElement( He, 100.0*perCent );
	G4Material* He4_2atm = new G4Material( "He4_2atm",  density = 2.*0.1786*mg/cm3, nel=1 );
	He4_2atm->AddElement( He, 100.0*perCent );
	G4Material* He4_3atm = new G4Material( "He4_3atm",  density = 3.*0.1786*mg/cm3, nel=1 );
	He4_3atm->AddElement( He, 100.0*perCent );
	G4Material* He4_7atm = new G4Material( "He4_7atm",  density = 7.*0.1786*mg/cm3, nel=1 );
	He4_7atm->AddElement( He, 100.0*perCent );
	
	G4Material* PbWO4    = new G4Material( "PbWO4",  density = 8.28*g/cm3, nel=3 );
	PbWO4->AddElement( Pb,   1./6.*100.*perCent );
	PbWO4->AddElement( Wf,   1./6.*100.*perCent );
	PbWO4->AddElement( Ox,   4./6.*100.*perCent );
	
	G4Material* SemiMirror    = new G4Material( "SemiMirror",  density = 8.28*g/cm3, nel=3 );
	SemiMirror->AddElement( Pb,   1./6.*100.*perCent );
	SemiMirror->AddElement( Wf,   1./6.*100.*perCent );
	SemiMirror->AddElement( Ox,   4./6.*100.*perCent );
	
	//for pol He3 target
	//He3 gas target
	G4Material* He3_10amg = new G4Material( "He3_10amg",  density = 10.*0.1345*mg/cm3, nel=1 );
	//0.1345=44.6(amg=mol/m3)*3.016(g/mol)
	He3_10amg->AddElement( He3, 100.0*perCent );
	
	//He3 cell glass
	G4Material *BariumOxide = new G4Material("BariumOxide", density=5.72*g/cm3, nel=2);
	BariumOxide->AddElement(Ba, 1);
	BariumOxide->AddElement(O,  1);
	
	G4Material *StrontiumOxide = new G4Material("StrontiumOxide", density=4.7*g/cm3, nel=2);
	StrontiumOxide->AddElement(Sr, 1);
	StrontiumOxide->AddElement(O,  1);
	

	//proton pol target NH3
	//solid NH3
	G4double density_NH3_solid = 0.817*g/cm3;
	G4Material *NH3_solid = new G4Material("NH3_solid", density_NH3_solid, nel=2);
	NH3_solid->AddElement(H, 3);
	NH3_solid->AddElement(N, 1);
	
	G4Material *He4_liquid = new G4Material( "He4_liquid",  density = 0.145*g/cm3, nel=1 );
	He4_liquid->AddElement( He, 100.0*perCent );
	
	//SolidNH3(55%)+LiquidHe(45%) in volumn?
	//density = mLiquidHeD*(1.0-mNH3VolumnRatio)+mSolidNH3D*mNH3VolumnRatio;
	//density_NH3He = 1.0/((1.0-mNH3WeightRatio)/mLiquidHeD+mNH3WeightRatio/mSolidNH3D);
	G4double density_NH3He = (0.817*0.55+0.145*0.45)*g/cm3;
	G4Material *NH3He = new G4Material("NH3He", density_NH3He, nel=2);
	NH3He->AddMaterial(NH3_solid, 0.8732);
	NH3He->AddMaterial(He4_liquid, 1-0.8732);
	
	//Beryllium oxide
	G4Material *BerylliumOxide =  matman->FindOrBuildMaterial("G4_BERYLLIUM_OXIDE");  // BeO 3.02g/cm3
	
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
	
	//G4Material *Methane = new G4Material("Methane", density= 0.422*g/cm3, nel=2);  //liquid
	//G4Material *Methane = new G4Material("Methane", density= 0.717*kg/m3, nel=2);  //gas at 0 degree, 1atm
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
	MMats["DCgas"]          = DCgas;
	MMats["FTinsfoam"]      = FTinsfoam;
	MMats["Glass"]          = Glass;
	MMats["Gold"]           = matman->FindOrBuildMaterial("G4_Au");
	MMats["H2O"]            = H2O;
	MMats["He4_1atm"]       = He4_1atm ;
	MMats["He4_2atm"]       = He4_2atm ;
	MMats["He4_3atm"]       = He4_3atm ;
	MMats["He4_7atm"]       = He4_7atm ;
	MMats["He3_10amg"]      = He3_10amg ;
	MMats["Iron"]           = matman->FindOrBuildMaterial("G4_Fe");
	MMats["Kryptonite"]     = matman->FindOrBuildMaterial("Kryptonite");
	MMats["Lead"]           = matman->FindOrBuildMaterial("G4_Pb");
	MMats["LH2"]            = matman->FindOrBuildMaterial("G4_lH2");;
	MMats["Methane"]        = Methane;
	MMats["MMGas"]          = MMGas;
	MMats["MMMesh"]         = MMMesh;
	MMats["MMMylar"]        = MMMylar;
	MMats["Nickel"]         =  matman->FindOrBuildMaterial("G4_Ni");
	MMats["Noryl"]          = Noryl ;
	MMats["PbWO4"]          = PbWO4 ;
	MMats["Quartz"]         = Quartz;
	MMats["SemiMirror"]     = SemiMirror ;
	MMats["Silicium"]       = matman->FindOrBuildMaterial("G4_Si");
	MMats["Silicon"]        = matman->FindOrBuildMaterial("G4_Si");
	MMats["svtwirebond"]    = svtwirebond;
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
	MMats["NH3_solid"]      = NH3_solid;
	MMats["He4_liquid"]     = He4_liquid;
	MMats["NH3He"]          = NH3He;
	MMats["BerylliumOxide"] = BerylliumOxide;
	
	
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

	G4double PhotonEnergy_SemiMirror[nEntries_Air]    = { 2.034*eV , 4.136*eV };
	G4double RefractiveIndex_SemiMirror[nEntries_Air] = { 5.00, 5.00 };
	G4double Absorption_SemiMirror[nEntries_Air]      = { 100.0*m    , 100.0*m  };
	
	G4MaterialPropertiesTable* SemiMirrorMPT = new G4MaterialPropertiesTable();
	SemiMirrorMPT->AddProperty("RINDEX",     PhotonEnergy_SemiMirror, RefractiveIndex_SemiMirror, 2);
	SemiMirrorMPT->AddProperty("ABSLENGTH",  PhotonEnergy_SemiMirror, Absorption_SemiMirror,      2);
	MMats["SemiMirror"]->SetMaterialPropertiesTable(SemiMirrorMPT);
	
	
	
	

	
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
