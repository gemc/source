// Qt headers
#include <QtSql>

// gemc headers
#include "material_factory.h"
#include "gdml_materials.h"
#include "string_utilities.h"
#include "utils.h"

// G4 headers
#include "G4Element.hh"
#include "G4NistManager.hh"
#include "G4OpBoundaryProcess.hh"

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;

map<string, G4Material*> gdml_materials::initMaterials(runConditions rc, goptions opts)
{

	string hd_msg    = opts.optMap["LOG_MSG"].args + " GDML Materials: >> ";
	double verbosity = opts.optMap["MATERIAL_VERBOSITY"].arg;


	map<string, G4Element*> elementsMap;
	map<string, G4Material*> materialsMap;

	// first check if there's at least one detector with MYSQL factory
	if(!check_if_factory_is_needed(rc.detectorConditionsMap, "GDML"))
		return materialsMap;

	// building materials from GDML file
	for(map<string, detectorCondition>::iterator it=rc.detectorConditionsMap.begin(); it != rc.detectorConditionsMap.end() ; it++)
	{
		if(it->second.get_factory() != "GDML" )
			continue;
		
		if(verbosity)
			cout << hd_msg << " Initializing " << it->second.get_factory() << " for detector " << it->first << endl;

		string dname     = it->first;
		string fname = dname + ".gdml";

		// If found, parse the gdml file
		QFile gdet(fname.c_str());
	    
		if( !gdet.exists() )
		{
			cout << hd_msg << "  Failed to open geometry file " << fname << " for system: " << dname << ". Maybe the filename doesn't exist? Exiting." << endl;
			exit(0);
		}
	  
		QDomDocument domDocument;
		// opening gcard and filling domDocument
		if(!domDocument.setContent(&gdet))
		{
			cout << hd_msg << "  Failed to open geometry file " << fname << " for system: " << dname << ". Wrong XML syntax. Exiting." << endl;
			exit(0);
		}
		gdet.close();

		QDomElement docElem = domDocument.documentElement(); // reading gcard file
		QDomNode n = docElem.firstChild();                   // looping over first tags
		while(!n.isNull())
		{
			///< converts the node to an element.
			QDomElement e = n.toElement();                		   
			// if the node is an element with tag materials
			if(!e.isNull() && e.tagName().toStdString() == "materials")     
			{
				QDomNode nn= e.firstChild();

				while( !nn.isNull() )
				{
					QDomElement ee = nn.toElement();
					//parses and adds elements to elementsMap
					if(ee.tagName().toStdString() == "element")
					{
						
						string ename = assignAttribute(ee, "name", "noname");
						int Z = assignAttribute(ee, "Z", 0);
						double molar_mass = 0;	
						
						QDomNode nnn= ee.firstChild();
						
						while( !nnn.isNull() )
						{
							QDomElement eee = nnn.toElement();

							if(eee.tagName().toStdString() == "atom")
								molar_mass = assignAttribute(eee, "value", 0.0);
							
							nnn=nnn.nextSibling();
							cout << ename << " Z: "<< Z << " molar_mass: " << molar_mass << endl; 
						}
	
						elementsMap[ename] = new G4Element(ename, ename, Z, molar_mass*g/mole);
					}

					if(ee.tagName().toStdString() == "material")
					{			
						vector<double> fractions;  
						vector<int> composite; 
 						vector<string> refs;

						string ename      = assignAttribute(ee, "name", "noname");
						int Z             = assignAttribute(ee, "Z", 0);
						double molar_mass = 0;	
						double density    = 0;

						// filling molar mass, density and fractions 
						QDomNode nnn= ee.firstChild();
						while( !nnn.isNull() )
						{
							QDomElement eee = nnn.toElement();

							if(eee.tagName().toStdString() == "D")
								density = assignAttribute(eee, "value", 0.0);	
						
							if(eee.tagName().toStdString() == "atom")
							{
								molar_mass = assignAttribute(eee, "value", 0.0);
							}
							
							if(eee.tagName().toStdString() == "fraction")
							{
								fractions.push_back(assignAttribute(eee, "n", 0.0));
								refs.push_back(assignAttribute(eee, "ref", "noref"));
							}
							if(eee.tagName().toStdString() == "composite")
							{
								fractions.push_back(assignAttribute(eee, "n", 0));
								refs.push_back(assignAttribute(eee, "ref", "noref"));
							}
							nnn=nnn.nextSibling();
						}
						//If the material doesn't have any fraction node, it is an element
						if(fractions.size() == 0)
							elementsMap[ename] = new G4Element(ename, ename, Z, molar_mass*g/mole);
						// defines a simple material
						if (Z>0 && molar_mass >0)
							materialsMap[ename]= new G4Material(ename, Z, molar_mass*g/mole, density*g/cm3);
						
						// defines molecules and mixtures
						if(fractions.size() != 0)
						{	
							cout << "fraction.size() " << fractions.size() << endl;					
							materialsMap[ename] =  new G4Material(ename, density*g/cm3, (int) fractions.size());
						
							for(unsigned int i=0; i<fractions.size(); i++)
							{
								cout << refs[i] << "     elementsMap :    " << elementsMap[refs[i]] << endl;
                                                                if(fractions[i] >= 1) materialsMap[ename]->AddElement(elementsMap[refs[i]], (int)fractions[i]);
								else                  materialsMap[ename]->AddElement(elementsMap[refs[i]], fractions[i]);
							}
						}

						
					}			
			 		
					nn=nn.nextSibling();
				}
			}
			n=n.nextSibling();
		}
	
	}


  
 	if(verbosity>0) printMaterials(materialsMap); 
	return materialsMap;
}
