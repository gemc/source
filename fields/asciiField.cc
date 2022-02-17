// gemc headers
#include "fieldFactory.h"

extern "C" {
#include "asciiField.h"
}

#include "string_utilities.h"
#include "gemcUtils.h"
#include "magfieldio.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "magfieldio.h"

// mlibrary
#include "gstring.h"
using namespace gstring;


// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;

bool asciiField::isEligible(string file)
{
  ifstream IN(file.c_str());
  
  if(!IN.is_open())
    return 0;
  
  string first;
  IN >> first;
  IN.close();
  
  if(strcmp (first.c_str(), "<mfield>") != 0)
    return 0;

  return 1;
}

// load field definitions
gfield asciiField::loadField(string file, goptions opts)
{
  gfield gf(opts);

  ifstream IN(file.c_str());
  string content = "";
  string stop    = "";
  
  while(IN.good() && stop != "</mfield>")
    {
      IN >> stop;
      content += stop + " " ;
    }
  IN.close();

  // sucking up all that string into a domdocument
  QDomDocument domDocument;
  domDocument.setContent(QString(content.c_str()));
  
  
  QDomElement docElem = domDocument.documentElement();
  QDomNode n = docElem.firstChild();
  while(!n.isNull())
    {
      QDomElement e = n.toElement();       // converts the node to an element.
      if(!e.isNull())                      // check that the node really is an element. 
	{            
	  if(e.tagName().toStdString() == "description")  ///< selecting "description" nodes
	    {
	      gf.name        = assignAttribute(e, "name", "na");
	      gf.factory     = assignAttribute(e, "factory", "na");
	      gf.description = assignAttribute(e, "comment", "no comment");
	    }
	  
	  if(e.tagName().toStdString() == "symmetry")     ///< selecting "symmetry" nodes
	    {
	      gf.symmetry    = assignAttribute(e, "type",   "na");
	      gf.format      = assignAttribute(e, "format", "na");
	    }
	  
	  // simple symmetry, looking for uniform field definition
	  if(gf.format == "simple" && gf.symmetry == "uniform")
	    {
	      if(e.tagName().toStdString() == "dimension")     ///< selecting "dimension" nodes
		{
		  string units = "*" + assignAttribute(e, "units", "gauss");
		  gf.dimensions  = assignAttribute(e, "bx", "0") + units + " ";
		  gf.dimensions += assignAttribute(e, "by", "0") + units + " " ; 
		  gf.dimensions += assignAttribute(e, "bz", "0") + units ; 
		}
	    }


	  // simple symmetry, looking for multipole field definition
	  if(gf.format == "simple" && gf.symmetry == "multipole")
	    {
	      if(e.tagName().toStdString() == "dimension")     ///< selecting "dimension" nodes
		{
		  gf.dimensions  = assignAttribute(e, "Npole", "0") + " ";
		  gf.dimensions += assignAttribute(e, "scale", "0") + "*" + assignAttribute(e, "Bunit",   "gauss") + " ";
		  gf.dimensions += assignAttribute(e, "x", "0")     + "*" + assignAttribute(e, "XYZunit", "cm")    + " ";
		  gf.dimensions += assignAttribute(e, "y", "0")     + "*" + assignAttribute(e, "XYZunit", "cm")    + " ";
		  gf.dimensions += assignAttribute(e, "z", "0")     + "*" + assignAttribute(e, "XYZunit", "cm")    + " ";
		  gf.dimensions += assignAttribute(e, "rot", "0")   + "*" + assignAttribute(e, "ROTunit", "deg")   + " ";
		  gf.dimensions += assignAttribute(e, "ROTaxis", "Y");
		}
	    }
	  
	  // map symmetry, looking for map field definition
	  if(gf.format == "map")
	    {
	      if(!gf.map) gf.map = new gMappedField(file, gf.symmetry);
	      
	      // selecting "map" nodes
	      // selecting "coordinate" nodes
	      if(e.tagName().toStdString() == "map")    
		{
		  QDomNode nn= e.firstChild();
		  while(!nn.isNull())
		    {
		      QDomElement ee = nn.toElement();
		      if(ee.tagName().toStdString() == "coordinate")    
			{
			  QDomNode nnn= ee.firstChild();
			  while(!nnn.isNull())
			    {
			      QDomElement eee = nnn.toElement();
			      
			      if(eee.tagName().toStdString() == "first" || eee.tagName().toStdString() == "second" || eee.tagName().toStdString() == "third")
				{
				  string name = assignAttribute(eee, "name", "na");
				  int np      = (int) assignAttribute(eee, "npoints", 0);
				  string unit = assignAttribute(eee, "units", "mm");
				  // loading min and max with their unit
				  double min  = assignAttribute(eee, "min", 0.0)*get_number("1*" + unit);
				  double max  = assignAttribute(eee, "max", 0.0)*get_number("1*" + unit);
				  int speed = 0;
				  if(eee.tagName().toStdString() == "second") speed = 1;
				  if(eee.tagName().toStdString() == "third")  speed = 2;
				  gf.map->coordinates.push_back(gcoord(name, np, min, max, unit, speed));
				}
			      nnn=nnn.nextSibling();
			    }
			}
		      
		      /// selecting "field" nodes. Default unit is gauss
		      if(ee.tagName().toStdString() == "field")   
			{
			  gf.map->unit = assignAttribute(ee, "unit", "gauss");
			}
		      nn = nn.nextSibling();
		    }
		}
	    }
	  
	}
      n = n.nextSibling();
    }
  // initialize field and field map
  gf.initialize(opts);
  
  // rescaling dimensions for uniform field
  if(gf.scaleFactor != 1 && gf.format == "simple" && gf.symmetry == "uniform")
    {
      vector<string> olddim = getStringVectorFromString(gf.dimensions);
      string newdim;
      for(unsigned int d=0; d<olddim.size(); d++)
	{
	  // default unit is "megatesla" since for volts is megavolt
	  newdim += stringify(get_number(olddim[d])*1000*gf.scaleFactor) + "*T " ;
	}
      gf.dimensions = trimSpacesFromString(newdim);
    }
  
  // rescaling dimensions for multipole field
  if(gf.scaleFactor != 1 && gf.format == "simple" && gf.symmetry == "multipole")
    {
      vector<string> olddim = getStringVectorFromString(gf.dimensions);
      string newdim;
      for(unsigned int d=0; d<olddim.size(); d++)
	{
	  if (d==1) newdim += stringify(gf.scaleFactor*atof(olddim[d].substr(0, olddim[d].find("*")).c_str()))
		      + "*" + trimSpacesFromString(olddim[d].substr(olddim[d].find("*")+1, olddim[d].find("*") + 20)) + " ";
	  else newdim += olddim[d]+" ";
	}
      gf.dimensions = trimSpacesFromString(newdim);
    }
  
  return gf;
}



// load field map
void asciiField::loadFieldMap(gMappedField* map, double v)
{
  cout << "  > Loading field map from " << map->identifier << " with symmetry: " << map->symmetry << endl;

  // dipole field
  if(map->symmetry == "dipole-x" || map->symmetry == "dipole-y" || map->symmetry == "dipole-z")
    loadFieldMap_Dipole(map, v);
  
  // cylindrical field
  else if(map->symmetry == "cylindrical-x" || map->symmetry == "cylindrical-y" || map->symmetry == "cylindrical-z")
    loadFieldMap_Cylindrical(map, v);

  // phi-segmented field

  else if(map->symmetry == "phi-segmented")
    loadFieldMap_phiSegmented(map, v);
  
  else if(map->symmetry == "cartesian_3D" || map->symmetry == "cartesian_3D_quadrant")
    loadFieldMap_cartesian3d(map, v);
  
  else {cout << "can't recognize the field symmetry "<< map->symmetry << endl; exit(0);}
}


void asciiField::loadFieldMap(gclas12BinaryMappedField* map, double v) {


  char *solenoidPath = (char*) malloc(255);
  char *torusSymmetricPath = (char*) malloc(255);
  char *torusFullPath = (char*) malloc(255);
  
  const char *dataDir;

  dataDir = "/w/hallb_scshelf2102/clas12/jnewton/binary/data/fieldmaps";

  sprintf(solenoidPath, "%s/Symm_solenoid_r601_phi1_z1201_13June2018.dat", dataDir);
  
  if(isSymmetric==true)  {
    sprintf(torusSymmetricPath, "%s/Symm_torus_r2501_phi16_z251_24Apr2018.dat",dataDir);//Absolute Path To Symmetric Torus
    symmetricTorus = initializeTorus("/w/hallb_scshelf2102/clas12/jnewton/binary/data/fieldmaps/Symm_torus_r2501_phi16_z251_24Apr2018.dat");
  }
  
  else  {
    sprintf(torusFullPath, "%s/Full_torus_r251_phi181_z251_03March2020.dat",dataDir);//Absolute Path To Full Torus
    fullTorus = initializeTorus("/w/hallb_scshelf2102/clas12/jnewton/binary/data/fieldmaps/Full_torus_r251_phi181_z251_03March2020.dat");
  }
  

}



// Examples

/*
 
clas12 solenoid:
<mfield>
<description name="clas12-solenoid" factory="ASCII" comment="clas12 superconducting solenoid"/>
<symmetry type="cylindrical-z" format="map"/>
<map>
<coordinate>
<first  name="transverse"    npoints="601"   min="0"  max="3" units="m"/>
<second name="longitudinal"  npoints="1201"  min="-3" max="3" units="m"/>
</coordinate>
<field unit="T"/>
</map>
</mfield>
clas12 torus:
<mfield>
<description name="clas12-torus" factory="ASCII" comment="clas12 superconducting torus"/>
<symmetry type="phi-segmented" format="map""/>
<map>
<coordinate>
<first  name="azimuthal"     npoints="61"   min="0"   max="30"  units="deg"/>
<second name="transverse"    npoints="126"  min="0"   max="500" units="cm"/>
<third  name="longitudinal"  npoints="126"  min="100" max="600" units="cm"/>
</coordinate>
<field unit="kilogauss"/>
</map>
</mfield>
Example of uniform field:
 <mfield>
 <description name="uniform" factory="ASCII" comment="Uniform 10 T Magnetic Field along x-axis"/>
 <symmetry type="uniform" format="simple"/>
 <dimension bx="10" by="0" bz="0" units="T"/>
</mfield>
*/













