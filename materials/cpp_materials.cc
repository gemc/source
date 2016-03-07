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
	G4int zz;
	G4int nn;
	G4double abundance;
	G4String symbol;
	G4int nIso;
	
	
	G4Element* As  = new G4Element("Arsen"   ,  "As", z=33, a=    74.92*g/mole);
	G4Element* Ar  = new G4Element("Argon",     "Ar", z=18, a=    39.95*g/mole);
	G4Element* Al  = new G4Element("Aluminum",  "Al", z=13, a=   26.982*g/mole);
	G4Element* Ba  = new G4Element("Barium",    "Ba", z=56, a=  137.327*g/mole);
	G4Element* C   = new G4Element("Carbon",    "C",  z=6,  a=    12.01*g/mole);
	G4Element* Cr  = new G4Element("Chromium",  "Cr", z=24, a=    52.00*g/mole);
	G4Element* Cl  = new G4Element("Chlorine",  "Cl", z=17, a=   35.453*g/mole);
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
	G4Element* W   = new G4Element("Tungsten",  "W",  z=74, a=   183.85*g/mole);
	G4Element* Wf  = new G4Element("Wf",        "Wf", z=74, a=   183.85*g/mole);
	
	G4NistManager* matman = G4NistManager::Instance();
	
	
	// common materials
	G4Material *Air = new G4Material("Air", density= 1.29*mg/cm3, nel=2);
	Air->AddMaterial(matman->FindOrBuildMaterial("G4_N"), 70.*perCent);
	Air->AddMaterial(matman->FindOrBuildMaterial("G4_O"), 30.*perCent);
    
	G4Material *Air_Opt = new G4Material("Air_Opt",   density= 1.29*mg/cm3, nel=2);
	Air_Opt->AddMaterial(matman->FindOrBuildMaterial("G4_N"), 70.*perCent);
	Air_Opt->AddMaterial(matman->FindOrBuildMaterial("G4_O"), 30.*perCent);
	
	G4Material *Air_Opt2 = new G4Material("Air_Opt2", density= 1.29*mg/cm3, nel=2);
	Air_Opt2->AddMaterial(matman->FindOrBuildMaterial("G4_N"), 70.*perCent);
	Air_Opt2->AddMaterial(matman->FindOrBuildMaterial("G4_O"), 30.*perCent);
	
	// this material will kill every tracks that touch it
	G4Material *Kryptonite = new G4Material("Kryptonite", density= 0.00000001*mg/cm3, nel=1);
	Kryptonite->AddElement(Ar, 100.*perCent);
	
	
	
	
	// Materials Below are organized by detector
	// Each detector have names and email associated with it
	
	
	
	// BST Support Material: Copper + HFCB
	// All dimensions in mm
	//          Size                 %
	// Gold:    0.0600              1.6
	// Copper:  0.1397 + 2.5       68.9
	// Polymide 0.2154              5.6    (Kapton)
	// FR4      0.9144             23.9
	// Total:   1.3295 + 2.5 = 3.8925
	//	G4Material *BSTSupport = new G4Material("BSTSupport",   density = 9.0*g/cm3, nel=3);
	//	BSTSupport->AddElement(G4_Au,        1.6*perCent);
	//	BSTSupport->AddElement(G4_Cu,       68.9*perCent);
	//	BSTSupport->AddElement(Kapton,       5.6*perCent);
	//	BSTSupport->AddElement(CarbonFiber, 23.9*perCent);
	
	
	
	
	
	
	
	
	
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// Ring Imaging CHerenkov (RICH)
	// By Marco Contalbrigo,
	//    Amhed
	// mcontalb@fe.infn.it
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// Forward Time-of-Flight (FTOF)
	// By Dan Carman,
	//    Maurizio Ungaro
	// carman@jlab.org
	// ungaro@jlab.org
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	

	
	
	
	
	G4Material *Aerogel = new G4Material("Aerogel",   density= 0.14*g/cm3, nel=2);
	Aerogel->AddMaterial(matman->FindOrBuildMaterial("G4_Si"), 30.*perCent);
	Aerogel->AddMaterial(matman->FindOrBuildMaterial("G4_O"),  70.*perCent);
	

	
	
	G4Isotope* D  = new G4Isotope("Deuteron", zz=1, nn=2, a= 2.0141018*g/mole);
	G4Element* elD = new G4Element("Deuterium",symbol="elD", nIso = 1);
	elD->AddIsotope(D, abundance=100.*perCent);
	
	G4Material *Carbon       = new G4Material("Carbon",       z=6,  a=  12.011*g/mole,  density =  2.265*g/cm3);
	G4Material *LHe          = new G4Material("LHe",          z=2,  a=  4.0026*g/mole,  density =  0.125*g/cm3);
	
	// aluminum honeycomb core  (it is actually made of alloy Alu-Alloy 3003 (AlMnCu)
	G4Material *AlHoneycomb = new G4Material("AlHoneycomb", z=13, a=  26.982*g/mole, density =  0.13*g/cm3);
	
	// a denser aluminum honeycomb core for comparison
	G4Material *AlHoneycomb2 = new G4Material("AlHoneycomb2", z=13, a=  26.982*g/mole, density =  1.00*g/cm3);
	
	// Densimet-180 alloy (roughly)
	G4Material *Densimet = new G4Material("Densimet",   density = 18.0*g/cm3, nel=3);
	Densimet->AddElement(W, 90*perCent);
	Densimet->AddElement(Ni, 5*perCent);
	Densimet->AddMaterial(matman->FindOrBuildMaterial("G4_Fe"), 5*perCent);
	
	
	
	
	G4Material *LeadTungsten = new G4Material("LeadTungsten",   density = 8.28*g/cm3, nel=3);
	LeadTungsten->AddElement(Pb,1);
	LeadTungsten->AddElement(W, 1);
	LeadTungsten->AddElement(O, 4);
	
	
	G4Material* LeadOxide   = matman->FindOrBuildMaterial("G4_LEAD_OXIDE");       // PbO
	G4Material* SilicOxide  = matman->FindOrBuildMaterial("G4_SILICON_DIOXIDE");  // SiO2	2.32g/cm3
	G4Material* PotasOxide  = matman->FindOrBuildMaterial("G4_POTASSIUM_OXIDE");  // K2O
	G4Material* SodMonOxide = matman->FindOrBuildMaterial("G4_SODIUM_MONOXIDE");  // Na2O
	
	G4Material* As2O3 = new G4Material("As2O3", density= 3.738*g/cm3, nel=2); // As2O3
	As2O3->AddElement( As, 2 );
	As2O3->AddElement( O,  3 );
	
	G4Material* LgTF1 = new G4Material("LgTF1", density= 3.86*g/cm3 , nel=5 );   // Lead-Glass
	LgTF1->AddMaterial( LeadOxide   , 0.5120 );
	LgTF1->AddMaterial( SilicOxide  , 0.4130 );
	LgTF1->AddMaterial( PotasOxide  , 0.0422 );
	LgTF1->AddMaterial( SodMonOxide , 0.0278 );
	LgTF1->AddMaterial( As2O3       , 0.0050 );
	
	
	
	
	
//	G4Material *Scintillator = new G4Material("Scintillator",   density = 1.032*g/cm3, nel=2);
//	Scintillator->AddElement(C, 9);
//	Scintillator->AddElement(H, 10);
	
	
	
	
	
	G4Material *ScintillatorB = new G4Material("ScintillatorB",   density = 1.032*g/cm3, nel=2);
	ScintillatorB->AddElement(C, 9);
	ScintillatorB->AddElement(H, 10);
	ScintillatorB->GetIonisation()->SetBirksConstant(0.126*mm/MeV);
	
	
	
	
	
	
	G4Material *StainlessSteel = new G4Material("StainlessSteel",   density = 8.02*g/cm3, nel=5);
	StainlessSteel->AddElement(Mn, 0.02);
	StainlessSteel->AddMaterial(matman->FindOrBuildMaterial("G4_Si"), 0.01);
	StainlessSteel->AddElement(Cr, 0.19);
	StainlessSteel->AddElement(Ni, 0.10);
	StainlessSteel->AddMaterial(matman->FindOrBuildMaterial("G4_Fe"), 0.68);
	
	G4Material *Kapton = new G4Material("Kapton",   density = 1.42*g/cm3, nel=4);
	Kapton->AddMaterial(matman->FindOrBuildMaterial("G4_H"), 0.026362);
	Kapton->AddMaterial(matman->FindOrBuildMaterial("G4_C"), 0.691133);
	Kapton->AddMaterial(matman->FindOrBuildMaterial("G4_N"), 0.073270);
	Kapton->AddMaterial(matman->FindOrBuildMaterial("G4_O"), 0.209235);
	
	
	G4Material *OptScint     = new G4Material("OptScint",       density = 1.032*g/cm3, nel=2);
	OptScint->AddElement(C, 9);
	OptScint->AddElement(H, 10);
	
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
	
	G4Material *AluminiumOxide =  matman->FindOrBuildMaterial("G4_ALUMINUM_OXIDE");  // Al2O3	3.97g/cm3
	G4Material *CalciumOxide =  matman->FindOrBuildMaterial("G4_CALCIUM_OXIDE");  // CaO	3.3g/cm3
	
	G4Material *Glass_GE180 = new G4Material("Glass_GE180", density=2.76*g/cm3, nel=5);
	Glass_GE180->AddMaterial(SilicOxide, 60.8*perCent);   //add 0.5 to make total as 1
	Glass_GE180->AddMaterial(BariumOxide, 18.2*perCent);
	Glass_GE180->AddMaterial(AluminiumOxide, 14.3*perCent);
	Glass_GE180->AddMaterial(CalciumOxide, 6.5*perCent);
	Glass_GE180->AddMaterial(StrontiumOxide, 0.25*perCent);
	//http://galileo.phys.virginia.edu/research/groups/spinphysics/glass_properties.html
	
	// 	G4Material *Glass_GE180 = new G4Material("Glass_GE180", density=2.76*g/cm3, nel=6);
	// 	Glass_GE180->AddMaterial(matman->FindOrBuildMaterial("G4_Si"), 1./3.*60.3*perCent);
	// 	Glass_GE180->AddElement(Ba, 1./2.*18.2*perCent);
	// 	Glass_GE180->AddElement(Al, 2./5.*14.3*perCent);
	// 	Glass_GE180->AddElement(Ca, 1./2.*6.5*perCent);
	// 	Glass_GE180->AddElement(Sr, 1./2.*0.25*perCent);
	// 	Glass_GE180->AddMaterial(matman->FindOrBuildMaterial("G4_O"), (2./3.*60.3+1./2.*18.2+3./5.*14.3+1./2.*6.5+1./2.*0.25)*perCent);
	
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
	
	// for GEM
	// NOMEX:
	G4Material* NOMEX_pure = new G4Material("NOMEX_pure", density = 1.38*g/cm3, 5);
	NOMEX_pure -> AddMaterial(matman->FindOrBuildMaterial("G4_H"),0.04);
	NOMEX_pure -> AddMaterial(matman->FindOrBuildMaterial("G4_C"),0.54);
	NOMEX_pure -> AddMaterial(matman->FindOrBuildMaterial("G4_N"),0.09);
	NOMEX_pure -> AddMaterial(matman->FindOrBuildMaterial("G4_O"),0.10);
	NOMEX_pure -> AddElement(Cl,0.23);
	
	G4Material* NOMEX = new G4Material("NOMEX",density = 0.04*g/cm3, 2);
	NOMEX -> AddMaterial(NOMEX_pure,0.45);
	NOMEX -> AddMaterial(Air,0.55);
	
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
	
	G4Material *AerogelQuartz = new G4Material("AerogelQuartz", density=2.200*g/cm3, nel=2);
	AerogelQuartz->AddMaterial(matman->FindOrBuildMaterial("G4_Si"), 1);
	AerogelQuartz->AddElement(O,  2);
	
	G4Material *RichAerogel3 = new G4Material("RichAerogel3", density=0.158*g/cm3, nel=2);
	RichAerogel3->AddMaterial(AerogelQuartz, 97.0*perCent);
	RichAerogel3->AddMaterial(H2O,          3.0*perCent);
	
	G4Material *RichAerogel5 = new G4Material("RichAerogel5", density=0.25*g/cm3, nel=2);
	RichAerogel5->AddMaterial(AerogelQuartz, 97.0*perCent);
	RichAerogel5->AddMaterial(H2O,          3.0*perCent);
	
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
	
	
	
	
	
	MMats["AerogelQuartz"]  = AerogelQuartz;
	MMats["Air"]            = Air;
	MMats["Air_Opt"]        = Air_Opt;
	MMats["Air_Opt2"]       = Air_Opt2;
	MMats["Alumi"]          = Alumi;
	MMats["Aluminum"]       = matman->FindOrBuildMaterial("G4_Al");
	MMats["Aerogel"]        = Aerogel;
	MMats["C5F12"]          = C5F12;
	MMats["C6F14"]          = C6F14;
	MMats["Carbon"]         = Carbon;
	MMats["Concrete"]       = matman->FindOrBuildMaterial("G4_CONCRETE");       // concrete
	MMats["Copper"]         = matman->FindOrBuildMaterial("G4_Cu");
	MMats["DCgas"]          = DCgas;
	MMats["FTinsfoam"]      = FTinsfoam;
	MMats["Glass"]          = Glass;
	MMats["Glass_GE180"]    = Glass_GE180;
	MMats["Gold"]           = matman->FindOrBuildMaterial("G4_Au");
	MMats["H2O"]            = H2O;
	MMats["He4_1atm"]       = He4_1atm ;
	MMats["He4_2atm"]       = He4_2atm ;
	MMats["He4_3atm"]       = He4_3atm ;
	MMats["He4_7atm"]       = He4_7atm ;
	MMats["He3_10amg"]      = He3_10amg ;
	MMats["Iron"]           = matman->FindOrBuildMaterial("G4_Fe");
	MMats["Kapton"]         = Kapton;
	MMats["Kryptonite"]     = matman->FindOrBuildMaterial("Kryptonite");
	MMats["Lead"]           = matman->FindOrBuildMaterial("G4_Pb");
	MMats["LeadTungsten"]   = LeadTungsten;
	MMats["LgTF1"]          = LgTF1;                                     // Lead-Glass
	MMats["LH2"]            = matman->FindOrBuildMaterial("G4_lH2");;
	MMats["LHe"]            = LHe;
	MMats["Methane"]        = Methane;
	MMats["MMGas"]          = MMGas;
	MMats["MMMesh"]         = MMMesh;
	MMats["MMMylar"]        = MMMylar;
	MMats["Nickel"]         =  matman->FindOrBuildMaterial("G4_Ni");
	MMats["Noryl"]          = Noryl ;
	MMats["PbWO4"]          = PbWO4 ;
	MMats["Quartz"]         = Quartz;
	MMats["OptScint"]       = OptScint;
	MMats["RichAerogel3"]   = RichAerogel3;
	MMats["RichAerogel5"]   = RichAerogel5;
//	MMats["Scintillator"]   = Scintillator;
	MMats["ScintillatorB"]  = ScintillatorB;
	MMats["SemiMirror"]     = SemiMirror ;
	MMats["Silicium"]       = matman->FindOrBuildMaterial("G4_Si");
	MMats["Silicon"]        = matman->FindOrBuildMaterial("G4_Si");
	MMats["StainlessSteel"] = StainlessSteel;
	MMats["svtwirebond"]    = svtwirebond;
	MMats["Teflon"]         = matman->FindOrBuildMaterial("G4_TEFLON");
	MMats["Tungsten"]       = matman->FindOrBuildMaterial("G4_W");
	MMats["Vacuum"]         = matman->FindOrBuildMaterial("G4_Galactic");
	MMats["vacuum_m9"]      = vacuum_m9;
	MMats["vacuum_m3"]      = vacuum_m3;
	MMats["Water"]          = matman->FindOrBuildMaterial("G4_WATER");
	MMats["AlHoneycomb"]    = AlHoneycomb;
	MMats["AlHoneycomb2"]   = AlHoneycomb2;
	MMats["Densimet"]       = Densimet;
	MMats["NOMEX_pure"]     = NOMEX_pure;
	MMats["NOMEX"]          = NOMEX;
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
	
	// Optical Scintillator
	const G4int nEntries_OptScint = 2;
	
	G4double PhotonEnergy_OptScint[nEntries_OptScint]    = { 2.034*eV , 4.136*eV };
	G4double RefractiveIndex_OptScint[nEntries_OptScint] = { 1.49     , 1.49   };
	G4double Absorption_OptScint[nEntries_OptScint]      = { 4.5*m    , 4.5*m  };
	G4double ScintilFast_OptScint[nEntries_OptScint]     = { 1.00     , 1.00   };
	G4double ScintilSlow_OptScint[nEntries_OptScint]     = { 10.00    , 10.00  };
	
	G4MaterialPropertiesTable* OptScint_MPT = new G4MaterialPropertiesTable();
	OptScint_MPT->AddProperty("RINDEX",        PhotonEnergy_OptScint, RefractiveIndex_OptScint, nEntries_OptScint);
	OptScint_MPT->AddProperty("ABSLENGTH",     PhotonEnergy_OptScint, Absorption_OptScint,      nEntries_OptScint);
	OptScint_MPT->AddProperty("FASTCOMPONENT", PhotonEnergy_OptScint, ScintilFast_OptScint,     nEntries_OptScint);
	OptScint_MPT->AddProperty("SLOWCOMPONENT", PhotonEnergy_OptScint, ScintilSlow_OptScint,     nEntries_OptScint);
	
	OptScint_MPT->AddConstProperty("SCINTILLATIONYIELD", 1./MeV);
	OptScint_MPT->AddConstProperty("RESOLUTIONSCALE",    1.0);
	OptScint_MPT->AddConstProperty("FASTTIMECONSTANT",   2.1*ns);
	OptScint_MPT->AddConstProperty("SLOWTIMECONSTANT",   10.*ns);
	OptScint_MPT->AddConstProperty("YIELDRATIO",         0.1);
	MMats["OptScint"]->SetMaterialPropertiesTable(OptScint_MPT);
	
	// Lead Glass
	const G4int nEntries = 15;
	G4double PhotonEnergyTF1[nEntries] =
	{ 1.99988723*eV, 2.06655002*eV, 2.13781047*eV,
		2.21416092*eV, 2.29616666*eV, 2.38448095*eV,
		2.47986007*eV, 2.58318758*eV, 2.69550014*eV,
		2.81802297*eV, 2.95221448*eV, 3.09982514*eV,
		3.26297379*eV, 3.44425011*eV, 3.64685321*eV };
	
	G4double RefractiveIndex1[nEntries] =
	{ 1.6436, 1.6454, 1.6474, 1.6496, 1.6521,
		1.6548, 1.6579, 1.6614, 1.6653, 1.6697,
		1.6748, 1.6806, 1.6874, 1.6952, 1.7045 };
	
	G4double Absorption1[nEntries] =
	{ 7.0*m, 7.0*m, 7.0*m, 7.0 *m, 7.0 *m,
		7.0*m, 6.0*m, 5.0*m, 3.5 *m, 2.75*m,
		1.8*m, 1.0*m, 0.5*m, 0.15*m, 0.08*m };
	
	G4MaterialPropertiesTable* TF1glass = new G4MaterialPropertiesTable();
	TF1glass->AddProperty("RINDEX",       PhotonEnergyTF1, RefractiveIndex1, nEntries);
	TF1glass->AddProperty("ABSLENGTH",    PhotonEnergyTF1, Absorption1,      nEntries);
	MMats["LgTF1"]->SetMaterialPropertiesTable(TF1glass);
	
	
	// Air Reflection
	const G4int nEntries_Air = 2;
	G4double PhotonEnergy_Air[nEntries_Air]    = { 2.034*eV , 4.136*eV };
	G4double RefractiveIndex_Air[nEntries_Air] = { 1.00, 1.00 };
	
	G4MaterialPropertiesTable* Air_MPT = new G4MaterialPropertiesTable();
	Air_MPT->AddProperty("RINDEX", PhotonEnergy_Air, RefractiveIndex_Air, nEntries_Air);
	MMats["Air_Opt"]->SetMaterialPropertiesTable(Air_MPT);
	// 	MMats["Air"]->SetMaterialPropertiesTable(Air_MPT);
	
	
	// Air2 Reflection
	G4double RefractiveIndex_Air2[nEntries_Air] = { 1.00, 1.00 };
	G4double Absorption_Air2[nEntries_Air]      = { 0.0*mm    , 0.0*mm  };
	
	G4MaterialPropertiesTable* Air_MPT2 = new G4MaterialPropertiesTable();
	Air_MPT2->AddProperty("RINDEX",     PhotonEnergy_Air, RefractiveIndex_Air2, nEntries_Air);
	Air_MPT2->AddProperty("ABSLENGTH",  PhotonEnergy_Air, Absorption_Air2,      nEntries_Air);
	MMats["Air_Opt2"]->SetMaterialPropertiesTable(Air_MPT2);
	
	G4double PhotonEnergy_SemiMirror[nEntries_Air]    = { 2.034*eV , 4.136*eV };
	G4double RefractiveIndex_SemiMirror[nEntries_Air] = { 5.00, 5.00 };
	G4double Absorption_SemiMirror[nEntries_Air]      = { 100.0*m    , 100.0*m  };
	
	G4MaterialPropertiesTable* SemiMirrorMPT = new G4MaterialPropertiesTable();
	SemiMirrorMPT->AddProperty("RINDEX",     PhotonEnergy_SemiMirror, RefractiveIndex_SemiMirror, 2);
	SemiMirrorMPT->AddProperty("ABSLENGTH",  PhotonEnergy_SemiMirror, Absorption_SemiMirror,      2);
	MMats["SemiMirror"]->SetMaterialPropertiesTable(SemiMirrorMPT);
	
	
	
	
	// Water Properties
	const G4int nEntries_Water = 32;
	G4double PhotonEnergy_Water[nEntries_Water] =
	{ 2.034*eV, 2.068*eV, 2.103*eV, 2.139*eV,
		2.177*eV, 2.216*eV, 2.256*eV, 2.298*eV,
		2.341*eV, 2.386*eV, 2.433*eV, 2.481*eV,
		2.532*eV, 2.585*eV, 2.640*eV, 2.697*eV,
		2.757*eV, 2.820*eV, 2.885*eV, 2.954*eV,
		3.026*eV, 3.102*eV, 3.181*eV, 3.265*eV,
		3.353*eV, 3.446*eV, 3.545*eV, 3.649*eV,
		3.760*eV, 3.877*eV, 4.002*eV, 4.136*eV };
	
	G4double RefractiveIndex_Water[nEntries_Water] =
	{ 1.3435  , 1.344   , 1.3445  , 1.345   ,
		1.3455  , 1.346   , 1.3465  , 1.347   ,
		1.3475  , 1.348   , 1.3485  , 1.3492  ,
		1.35    , 1.3505  , 1.351   , 1.3518  ,
		1.3522  , 1.3530  , 1.3535  , 1.354   ,
		1.3545  , 1.355   , 1.3555  , 1.356   ,
		1.3568  , 1.3572  , 1.358   , 1.3585  ,
		1.359   , 1.3595  , 1.36    , 1.3608   };
	
	G4double Absorption_Water[nEntries_Water] =
	{  3.448*m ,  4.082*m,  6.329*m,  9.174*m,
		12.346*m , 13.889*m, 15.152*m, 17.241*m,
		18.868*m , 20.000*m, 26.316*m, 35.714*m,
		45.455*m , 47.619*m, 52.632*m, 52.632*m,
		55.556*m , 52.632*m, 52.632*m, 47.619*m,
		45.455*m , 41.667*m, 37.037*m, 33.333*m,
		30.000*m , 28.500*m, 27.000*m, 24.500*m,
		22.000*m , 19.500*m, 17.500*m, 14.500*m };
	
	G4double ScintilFast_Water[nEntries_Water] =
	{ 1.00    ,  1.00   ,  1.00   ,  1.00   ,
		1.00    ,  1.00   ,  1.00   ,  1.00   ,
		1.00    ,  1.00   ,  1.00   ,  1.00   ,
		1.00    ,  1.00   ,  1.00   ,  1.00   ,
		1.00    ,  1.00   ,  1.00   ,  1.00   ,
		1.00    ,  1.00   ,  1.00   ,  1.00   ,
		1.00    ,  1.00   ,  1.00   ,  1.00   ,
		1.00    ,  1.00   ,  1.00   ,  1.00    };
	
	G4double ScintilSlow_Water[nEntries_Water] =
	{ 0.01    ,  1.00   ,  2.00   ,  3.00   ,
		4.00    ,  5.00   ,  6.00   ,  7.00   ,
		8.00    ,  9.00   ,  8.00   ,  7.00   ,
		6.00    ,  4.00   ,  3.00   ,  2.00   ,
		1.00    ,  0.01   ,  1.00   ,  2.00   ,
		3.00    ,  4.00   ,  5.00   ,  6.00   ,
		7.00    ,  8.00   ,  9.00   ,  8.00   ,
		7.00    ,  6.00   ,  5.00   ,  4.00    };
	
	
	G4MaterialPropertiesTable* Water_MPT = new G4MaterialPropertiesTable();
	Water_MPT->AddProperty("RINDEX",        PhotonEnergy_Water, RefractiveIndex_Water, nEntries_Water);
	Water_MPT->AddProperty("ABSLENGTH",     PhotonEnergy_Water, Absorption_Water,      nEntries_Water);
	Water_MPT->AddProperty("FASTCOMPONENT", PhotonEnergy_Water, ScintilFast_Water,     nEntries_Water);
	Water_MPT->AddProperty("SLOWCOMPONENT", PhotonEnergy_Water, ScintilSlow_Water,     nEntries_Water);
	
	Water_MPT->AddConstProperty("SCINTILLATIONYIELD", 10./MeV);
	Water_MPT->AddConstProperty("RESOLUTIONSCALE",    1.0);
	Water_MPT->AddConstProperty("FASTTIMECONSTANT",   1.*ns);
	Water_MPT->AddConstProperty("SLOWTIMECONSTANT",  10.*ns);
	Water_MPT->AddConstProperty("YIELDRATIO",         0.8);
	
	MMats["Water"]->SetMaterialPropertiesTable(Water_MPT);
	MMats["H2O"]->SetMaterialPropertiesTable(Water_MPT);
	
	
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
	
	const G4int RichAerogel_Len=50;
	G4double RichAerogel_PhoE[RichAerogel_Len]=
	{   1.87855*eV,  1.96673*eV,  2.05490*eV,  2.14308*eV,  2.23126*eV,
		2.31943*eV,  2.40761*eV,  2.49579*eV,  2.58396*eV,  2.67214*eV,
		2.76032*eV,  2.84849*eV,  2.93667*eV,  3.02485*eV,  3.11302*eV,
		3.20120*eV,  3.28938*eV,  3.37755*eV,  3.46573*eV,  3.55391*eV,
		3.64208*eV,  3.73026*eV,  3.81844*eV,  3.90661*eV,  3.99479*eV,
		4.08297*eV,  4.17114*eV,  4.25932*eV,  4.34750*eV,  4.43567*eV,
		4.52385*eV,  4.61203*eV,  4.70020*eV,  4.78838*eV,  4.87656*eV,
		4.96473*eV,  5.05291*eV,  5.14109*eV,  5.22927*eV,  5.31744*eV,
		5.40562*eV,  5.49380*eV,  5.58197*eV,  5.67015*eV,  5.75833*eV,
		5.84650*eV,  5.93468*eV,  6.02286*eV,  6.11103*eV,  6.19921*eV };
	G4double RichAerogel_Rind3[RichAerogel_Len]=
	{   1.02825,  1.02829,  1.02834,  1.02839,  1.02844,
		1.02849,  1.02854,  1.02860,  1.02866,  1.02872,
		1.02878,  1.02885,  1.02892,  1.02899,  1.02906,
		1.02914,  1.02921,  1.02929,  1.02938,  1.02946,
		1.02955,  1.02964,  1.02974,  1.02983,  1.02993,
		1.03003,  1.03014,  1.03025,  1.03036,  1.03047,
		1.03059,  1.03071,  1.03084,  1.03096,  1.03109,
		1.03123,  1.03137,  1.03151,  1.03166,  1.03181,
		1.03196,  1.03212,  1.03228,  1.03244,  1.03261,
		1.03279,  1.03297,  1.03315,  1.03334,  1.03354 };
	G4double RichAerogel_Abs[RichAerogel_Len]=
	{   17.5000*cm,  17.7466*cm,  17.9720*cm,  18.1789*cm,  18.3694*cm,
		18.5455*cm,  18.7086*cm,  18.8602*cm,  19.0015*cm,  19.1334*cm,
		19.2569*cm,  19.3728*cm,  19.4817*cm,  19.5843*cm,  19.6810*cm,
		19.7725*cm,  19.8590*cm,  19.9410*cm,  20.0188*cm,  20.0928*cm,
		18.4895*cm,  16.0174*cm,  13.9223*cm,  12.1401*cm,  10.6185*cm,
		9.3147*cm,   8.1940*cm,   7.2274*cm,   6.3913*cm,   5.6659*cm,
		5.0347*cm,   4.4841*cm,   4.0024*cm,   3.5801*cm,   3.2088*cm,
		2.8817*cm,   2.5928*cm,   2.3372*cm,   2.1105*cm,   1.9090*cm,
		1.7296*cm,   1.5696*cm,   1.4266*cm,   1.2986*cm,   1.1837*cm,
		1.0806*cm,   0.9877*cm,   0.9041*cm,   0.8286*cm,   0.7603*cm };
	G4double RichAerogel_Scat[RichAerogel_Len]=
	{   23.4256*cm,  19.4987*cm,  16.3612*cm,  13.8302*cm,  11.7702*cm,
		10.0798*cm,  8.6823*cm,   7.5188*cm,   6.5439*cm,   5.7219*cm,
		5.0251*cm,   4.4312*cm,   3.9225*cm,   3.4847*cm,   3.1064*cm,
		2.7780*cm,   2.4919*cm,   2.2417*cm,   2.0221*cm,   1.8288*cm,
		1.6580*cm,   1.5067*cm,   1.3723*cm,   1.2525*cm,   1.1455*cm,
		1.0497*cm,   0.9637*cm,   0.8864*cm,   0.8166*cm,   0.7536*cm,
		0.6965*cm,   0.6448*cm,   0.5977*cm,   0.5549*cm,   0.5159*cm,
		0.4802*cm,   0.4475*cm,   0.4176*cm,   0.3901*cm,   0.3649*cm,
		0.3417*cm,   0.3203*cm,   0.3005*cm,   0.2822*cm,   0.2653*cm,
		0.2497*cm,   0.2352*cm,   0.2217*cm,   0.2092*cm,   0.1975*cm };
	
	G4double RichAerogel_Rind5[RichAerogel_Len]; // shift rindex from 1.03 to 1.05
	for (int iii=0; iii<RichAerogel_Len; iii++)
		RichAerogel_Rind5[iii]=RichAerogel_Rind3[iii]+(1.05-1.03);
	
	G4MaterialPropertiesTable* RichAerogel3_MPT = new G4MaterialPropertiesTable();
	RichAerogel3_MPT->AddProperty("RINDEX",    RichAerogel_PhoE, RichAerogel_Rind3, RichAerogel_Len);
	RichAerogel3_MPT->AddProperty("ABSLENGTH", RichAerogel_PhoE, RichAerogel_Abs,   RichAerogel_Len);
	RichAerogel3_MPT->AddProperty("RAYLEIGH",  RichAerogel_PhoE, RichAerogel_Scat,  RichAerogel_Len);
	MMats["RichAerogel3"]->SetMaterialPropertiesTable(RichAerogel3_MPT);
	
	G4MaterialPropertiesTable* RichAerogel5_MPT = new G4MaterialPropertiesTable();
	RichAerogel5_MPT->AddProperty("RINDEX",    RichAerogel_PhoE, RichAerogel_Rind5, RichAerogel_Len);
	RichAerogel5_MPT->AddProperty("ABSLENGTH", RichAerogel_PhoE, RichAerogel_Abs,   RichAerogel_Len);
	RichAerogel5_MPT->AddProperty("RAYLEIGH",  RichAerogel_PhoE, RichAerogel_Scat,  RichAerogel_Len);
	MMats["RichAerogel5"]->SetMaterialPropertiesTable(RichAerogel5_MPT);
	
	
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
