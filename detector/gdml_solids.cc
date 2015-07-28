
// c++ headers
#include <string>
#include "utils.h"
#include "gdml_det_factory.h"

using namespace std;

// Qt4 headers
#include <QDomDocument>


string rad_to_deg (string radValue) //converts a radian value to a degree value
{
	string newValue = radValue;

	if (radValue == "TWOPI")
	{
		newValue = "360";
		
	}
	if (radValue == "PI")
	{
		newValue = "180";
		
	}
	if (radValue == "HALFPI")
	{
		newValue = "90";
		
	}

	return newValue;
	

}

string rad_to_deg_u (string radValue, string radUnit) /// converts the unit radian to degree if there is radian value.
{
	string newUnit = radUnit;
	
	if (radValue == "360" || radValue == "180" || radValue == "90" )
	{
		newUnit = "deg";
		
	}


	return newUnit;
	

}

//This function parses solid nodes and get values using gconstantMap and resolve_string.
//Return value is a string with dimension values and units.
string get_dimensions(QDomNode solidNode, map<string, double> gconstantMap)
{
	string dimension;
				
	QDomElement ee = solidNode.toElement();
				
				
	if(ee.tagName().toStdString() == "box"){
				
		string x = assignAttribute(ee, "x", "0");
		string y = assignAttribute(ee, "y", "0");
		string z = assignAttribute(ee, "z", "0");
				
		string c_x = stringify(resolve_string(x, gconstantMap));
		string c_y = stringify(resolve_string(y, gconstantMap));
		string c_z = stringify(resolve_string(z, gconstantMap));

		string lunits = assignAttribute(ee, "lunit", "mm");

		dimension =  c_x + "*" + lunits + " " + c_y + "*" + lunits + " " + c_z + "*" + lunits;
					
	}
	if(ee.tagName().toStdString() == "cone"){

		string rmin1 = assignAttribute(ee, "rmin1", "0");
    		string rmax1 = assignAttribute(ee, "rmax1", "0");
		string rmin2 = assignAttribute(ee, "rmin2", "0");
             				string rmax2 = assignAttribute(ee, "rmax2", "0");
		string z = assignAttribute(ee, "z", "0");
	
		string startphi = rad_to_deg (assignAttribute(ee, "startphi", "0"));
		string deltaphi = rad_to_deg (assignAttribute(ee, "deltaphi", "0"));

		string c_rmin1 = stringify(resolve_string(rmin1, gconstantMap));
		string c_rmax1 = stringify(resolve_string(rmax1, gconstantMap));
		string c_rmin2 = stringify(resolve_string(rmin2, gconstantMap));
		string c_rmax2 = stringify(resolve_string(rmax2, gconstantMap));
		string c_z = stringify(resolve_string(z, gconstantMap));
		string c_startphi = stringify(resolve_string(startphi, gconstantMap));
		string c_deltaphi = stringify(resolve_string(deltaphi, gconstantMap));

		string lunits = assignAttribute(ee, "lunit", "mm");
    		string aunits = assignAttribute(ee, "aunit", "deg");

		if(rad_to_deg_u (c_startphi, aunits) == "deg" || rad_to_deg_u (c_deltaphi, aunits) == "deg")
			aunits = "deg";
					
		dimension = c_rmin1 + "*" + lunits + " " + c_rmax1 + "*" + lunits + " " + c_rmin2 + "*" + lunits + " "+  c_rmax2 + "*" + lunits + " " + c_z + "*" + lunits + " " + 
						    c_startphi + "*" + aunits + " " + c_deltaphi + "*" + aunits;
	}

	if(ee.tagName().toStdString() == "para"){
			
		string x = assignAttribute(ee, "x", "0");
    		string y = assignAttribute(ee, "y", "0");
		string z = assignAttribute(ee, "z", "0");

		string alpha = rad_to_deg (assignAttribute(ee, "alpha", "0"));
		string theta = rad_to_deg (assignAttribute(ee, "theta", "0"));
		string phi = rad_to_deg (assignAttribute(ee, "phi", "0"));

		string c_x = stringify(resolve_string(x, gconstantMap));
		string c_y = stringify(resolve_string(y, gconstantMap));
		string c_z = stringify(resolve_string(z, gconstantMap));
		string c_alpha = stringify(resolve_string(alpha, gconstantMap));
		string c_theta = stringify(resolve_string(theta, gconstantMap));
		string c_phi= stringify(resolve_string(phi, gconstantMap));

		string lunits = assignAttribute(ee, "lunit", "mm");
    		string aunits = assignAttribute(ee, "aunit", "deg");

		if(rad_to_deg_u (c_alpha, aunits) == "deg" || rad_to_deg_u (c_theta, aunits) == "deg" || rad_to_deg_u (c_phi, aunits) == "deg")
			aunits = "deg";

		dimension =  c_x + "*" + lunits + " " + c_y +"*" + lunits + " " + c_z + "*" + lunits + " " + c_alpha + "*" + aunits + " " + c_theta + "*" + aunits + " " + c_phi + "*" + aunits;
	}

	if(ee.tagName().toStdString() == "sphere"){

		string rmin = assignAttribute(ee, "rmin", "0");
    		string rmax = assignAttribute(ee, "rmax", "0");
				
		string startphi = rad_to_deg (assignAttribute(ee, "startphi", "0"));
		string deltaphi =rad_to_deg ( assignAttribute(ee, "deltaphi", "0"));
		string starttheta = rad_to_deg (assignAttribute(ee, "starttheta", "0"));
		string deltatheta = rad_to_deg (assignAttribute(ee, "deltatheta", "0"));

		string c_rmin = stringify(resolve_string(rmin, gconstantMap));
		string c_rmax = stringify(resolve_string(rmax, gconstantMap));
		string c_startphi = stringify(resolve_string(startphi, gconstantMap));
		string c_deltaphi = stringify(resolve_string(deltaphi, gconstantMap));
		string c_starttheta = stringify(resolve_string(starttheta, gconstantMap));
		string c_deltatheta = stringify(resolve_string(deltatheta, gconstantMap));

				
		string lunits = assignAttribute(ee, "lunit", "mm");
    		string aunits = assignAttribute(ee, "aunit", "deg");

					
		if(rad_to_deg_u (c_startphi, aunits) == "deg" || rad_to_deg_u (c_deltaphi, aunits) == "deg" || rad_to_deg_u (c_starttheta, aunits) == "deg" ||rad_to_deg_u (c_deltatheta, aunits)=="deg" )
			aunits = "deg";

		dimension = c_rmin + "*" + lunits + " " + c_rmax + "*" + lunits + " " 
							+ c_startphi + "*" + aunits + " " + c_deltaphi + "*" + aunits + " " + c_starttheta + "*" + aunits + " " + c_deltatheta + "*" + aunits;
					
	}
				
	//can't find arb8 and trap.

	if(ee.tagName().toStdString() == "trd"){
					
		string x1 = assignAttribute(ee, "x1", "0");
    		string y1 = assignAttribute(ee, "y1", "0");
		string x2 = assignAttribute(ee, "x2", "0");
    		string y2 = assignAttribute(ee, "y2", "0");
		string z = assignAttribute(ee, "z", "0");

		string c_x1 = stringify(resolve_string(x1, gconstantMap));
		string c_y1 = stringify(resolve_string(y1, gconstantMap));
		string c_x2 = stringify(resolve_string(x2, gconstantMap));
		string c_y2 = stringify(resolve_string(y2, gconstantMap));
		string c_z = stringify(resolve_string(z, gconstantMap));
			
		string lunits = assignAttribute(ee, "lunit", "mm");

		dimension = c_x1 + "*" + lunits + " " + c_y1 + "*" + lunits + " " + c_x2 + "*" + lunits + " " + c_y2 + "*" + lunits + " " + c_z + "*" + lunits;
	}         
 
	if(ee.tagName().toStdString() == "tube"){
					
		string rmin = assignAttribute(ee, "rmin", "0");
    		string rmax = assignAttribute(ee, "rmax", "0");
		string z = assignAttribute(ee, "z", "0");

		string startphi = rad_to_deg (assignAttribute(ee, "startphi", "0"));
		string deltaphi = rad_to_deg (assignAttribute(ee, "deltaphi", "0"));

		string c_rmin = stringify(resolve_string(rmin, gconstantMap));
		string c_rmax = stringify(resolve_string(rmax, gconstantMap));
		string c_z = stringify(resolve_string(z, gconstantMap));
		string c_startphi = stringify(resolve_string(startphi, gconstantMap));
		string c_deltaphi = stringify(resolve_string(deltaphi, gconstantMap));

		string lunits = assignAttribute(ee, "lunit", "mm");
    		string aunits = assignAttribute(ee, "aunit", "deg");

					
		if(rad_to_deg_u (c_startphi, aunits) == "deg" || rad_to_deg_u (c_deltaphi, aunits) == "deg")
			aunits = "deg";
				

		dimension =c_rmin + "*" + lunits + " " + c_rmax + "*" + lunits + " " +c_z + "*" + lunits + " " + c_startphi + "*" + aunits + " " + c_deltaphi + "*" + aunits;
	}           
				
	//union, subtraction, polycone, telleselated
				   
  if(ee.tagName().toStdString() == "torus"){

		string rmin = assignAttribute(ee, "rmin", "0");
    		string rmax = assignAttribute(ee, "rmax", "0");
		string rtor = assignAttribute(ee, "rtor", "0");
					
		string startphi = rad_to_deg (assignAttribute(ee, "startphi", "0"));
		string deltaphi = rad_to_deg (assignAttribute(ee, "deltaphi", "0"));

		string c_rmin = stringify(resolve_string(rmin, gconstantMap));
		string c_rmax = stringify(resolve_string(rmax, gconstantMap));
		string c_rtor = stringify(resolve_string(rtor, gconstantMap));
		string c_startphi = stringify(resolve_string(startphi, gconstantMap));
		string c_deltaphi = stringify(resolve_string(deltaphi, gconstantMap));

		string lunits = assignAttribute(ee, "lunit", "mm");
    		string aunits = assignAttribute(ee, "aunit", "deg");

		if(rad_to_deg_u (c_startphi, aunits) == "deg" || rad_to_deg_u (c_deltaphi, aunits) == "deg")
			aunits = "deg";

				
		dimension = c_rmin + "*" + lunits + " " + c_rmax + "*" + lunits + " " + c_rtor + "*" + lunits + " " + c_startphi + "*" + aunits + " " + c_deltaphi + "*" + aunits;
	}           	
				
	if(ee.tagName().toStdString() == "orb"){
					
		string rmax = assignAttribute(ee, "rmax", "0");
		string c_rmax = stringify(resolve_string(rmax, gconstantMap));

		string lunits = assignAttribute(ee, "lunit", "mm");

		dimension = c_rmax +  "*" + lunits;
	}           	
	//polyhedra

	if(ee.tagName().toStdString() == "hype"){
			
		string rmin = assignAttribute(ee, "rmin", "0");
  		string rmax = assignAttribute(ee, "rmax", "0");
				
		string inst = rad_to_deg (assignAttribute(ee, "inst", "0"));
		string outst = rad_to_deg (assignAttribute(ee, "outst", "0"));
					
		string z = assignAttribute(ee, "z", "0");

		string c_rmin = stringify(resolve_string(rmin, gconstantMap));
		string c_rmax = stringify(resolve_string(rmax, gconstantMap));
		string c_inst = stringify(resolve_string(inst, gconstantMap));
		string c_outst = stringify(resolve_string(outst, gconstantMap));
		string c_z = stringify(resolve_string(z, gconstantMap));

		string lunits = assignAttribute(ee, "lunit", "mm");
    		string aunits = assignAttribute(ee, "aunit", "deg");
	
		if(rad_to_deg_u (c_inst, aunits) == "deg" || rad_to_deg_u (c_outst, aunits) == "deg")
			aunits = "deg";


		dimension = c_rmin + "*" + lunits + " " + c_rmax + "*" + lunits + " "  + c_inst + "*" + aunits + " " + c_outst + "*" + aunits + " " +  c_z + "*" + lunits;
	}      

	if(ee.tagName().toStdString() == "eltube"){
			
		string x = assignAttribute(ee, "dx", "0");
    		string y = assignAttribute(ee, "dy", "0");
		string z = assignAttribute(ee, "dz", "0");

		string c_x = stringify(resolve_string(x, gconstantMap));
		string c_y = stringify(resolve_string(y, gconstantMap));
		string c_z = stringify(resolve_string(z, gconstantMap));
				
		string lunits = assignAttribute(ee, "lunit", "mm");

		dimension =  c_x + "*" + lunits + " " + c_y + "*" + lunits + " " + c_z + "*" + lunits;
				
	}  

	if(ee.tagName().toStdString() == "ellipsoid"){
			
		string ax = assignAttribute(ee, "ax", "0");
    		string by = assignAttribute(ee, "by", "0");
		string cz = assignAttribute(ee, "cz", "0");
		string zcut1 = assignAttribute(ee, "zcut1", "0");
		string zcut2 = assignAttribute(ee, "zcut2", "0");

		string c_ax = stringify(resolve_string(ax, gconstantMap));
		string c_by = stringify(resolve_string(by, gconstantMap));
		string c_cz = stringify(resolve_string(cz, gconstantMap));
		string c_zcut1 = stringify(resolve_string(zcut1, gconstantMap));
		string c_zcut2 = stringify(resolve_string(zcut2, gconstantMap));

		string lunits = assignAttribute(ee, "lunit", "mm");

		dimension = c_ax +"*" + lunits + " "  + c_by + "*" + lunits + " " + c_cz + "*" + lunits+ " " + c_zcut1 + "*" + lunits + " " + c_zcut2 + "*" + lunits;
					
	}  
	//elcone
				
//------------------------------followed parameters in solids.gdml

	if(ee.tagName().toStdString() == "paraboloid"){
			
		string rlo = assignAttribute(ee, "rlo", "0");
    		string rhi = assignAttribute(ee, "rhi", "0");
		string z = assignAttribute(ee, "z", "0");
			
		string c_rlo = stringify(resolve_string(rlo, gconstantMap));
		string c_rhi = stringify(resolve_string(rhi, gconstantMap));
		string c_z = stringify(resolve_string(z, gconstantMap));

		string lunits = assignAttribute(ee, "lunit", "mm");

		dimension = c_rlo + "*" + lunits + " " + c_rhi + "*" + lunits + " " + c_z + "*" + lunits;
					
	}  

	if(ee.tagName().toStdString() == "tet"){
					
		string vertex1 = assignAttribute(ee, "vertex1", "0");
		string vertex2 = assignAttribute(ee, "vertex2", "0");
		string vertex3 = assignAttribute(ee, "vertex3", "0");
		string vertex4 = assignAttribute(ee, "vertex4", "0");

		string c_vertex1 = stringify(resolve_string(vertex1, gconstantMap));
		string c_vertex2 = stringify(resolve_string(vertex2, gconstantMap));
		string c_vertex3 = stringify(resolve_string(vertex3, gconstantMap));
		string c_vertex4 = stringify(resolve_string(vertex4, gconstantMap));

		dimension = c_vertex1 + ", " + c_vertex2 + ", "  + c_vertex3 + ", " + c_vertex4;
	}      

	if(ee.tagName().toStdString() == "twistedbox"){
			
		string PhiTwist = rad_to_deg (assignAttribute(ee, "PhiTwist", "0"));
				
		string x = assignAttribute(ee, "x", "0");
    		string y = assignAttribute(ee, "y", "0");
		string z = assignAttribute(ee, "z", "0");

		string c_PhiTwist = stringify(resolve_string(PhiTwist, gconstantMap));
		string c_x = stringify(resolve_string(x, gconstantMap));
		string c_y = stringify(resolve_string(y, gconstantMap));
		string c_z = stringify(resolve_string(z, gconstantMap));
				
		string lunits = assignAttribute(ee, "lunit", "mm");
		string aunits = assignAttribute(ee, "aunit", "deg");
	
		rad_to_deg_u (c_PhiTwist, aunits);

		dimension = c_PhiTwist + "*" + aunits + " " + c_x + "*" + lunits + " " + c_y + "*" + lunits + " " + c_z + "*" + lunits;
					
	}      
				
	if(ee.tagName().toStdString() == "twistedtrd"){
			
		string PhiTwist = rad_to_deg (assignAttribute(ee, "PhiTwist", "0"));
		
		string x1 = assignAttribute(ee, "x1", "0");
    		string y1 = assignAttribute(ee, "y1", "0");
		string x2 = assignAttribute(ee, "x2", "0");
    		string y2 = assignAttribute(ee, "y2", "0");
		string z = assignAttribute(ee, "z", "0");

		string c_PhiTwist = stringify(resolve_string(PhiTwist, gconstantMap));
		string c_x1 = stringify(resolve_string(x1, gconstantMap));
		string c_y1 = stringify(resolve_string(y1, gconstantMap));
		string c_x2 = stringify(resolve_string(x2, gconstantMap));
		string c_y2 = stringify(resolve_string(y2, gconstantMap));
		string c_z = stringify(resolve_string(z, gconstantMap));

		string lunits = assignAttribute(ee, "lunit", "mm");
		string aunits = assignAttribute(ee, "aunit", "deg");
	
		rad_to_deg_u (c_PhiTwist, aunits);

		dimension = c_PhiTwist  + "*" + aunits + " " + c_x1 + "*" + lunits + " "  + c_x2 + "*" + lunits + " " + c_y1 + "*" + lunits + " " +  c_y2 + "*" + lunits + " " + c_z + "*" + lunits;

	} 

	if(ee.tagName().toStdString() == "twistedtrap"){
			
		string PhiTwist = rad_to_deg (assignAttribute(ee, "PhiTwist", "0"));
					
		string z = assignAttribute(ee, "z", "0");
    		string theta = rad_to_deg (assignAttribute(ee, "theta", "0"));
		string phi = rad_to_deg (assignAttribute(ee, "phi", "0"));

    		string y1 = assignAttribute(ee, "y1", "0");
		string y2= assignAttribute(ee, "y2", "0");
		string x1 = assignAttribute(ee, "x1", "0");
   		string x2 = assignAttribute(ee, "x2", "0");
		string x3 = assignAttribute(ee, "x3", "0");
    		string x4 = assignAttribute(ee, "x4", "0");

		string alpha = rad_to_deg (assignAttribute(ee, "alpha", "0"));

		string c_PhiTwist = stringify(resolve_string(PhiTwist, gconstantMap));
		string c_z = stringify(resolve_string(z, gconstantMap));
		string c_theta = stringify(resolve_string(theta, gconstantMap));
		string c_phi = stringify(resolve_string(phi, gconstantMap));
		string c_y1 = stringify(resolve_string(y1, gconstantMap));
		string c_y2 = stringify(resolve_string(y2, gconstantMap));
		string c_x1 = stringify(resolve_string(x1, gconstantMap));
		string c_x2 = stringify(resolve_string(x2, gconstantMap));
		string c_x3 = stringify(resolve_string(x3, gconstantMap));
		string c_x4 = stringify(resolve_string(x4, gconstantMap));
				
		string lunits = assignAttribute(ee, "lunit", "mm");
		string aunits = assignAttribute(ee, "aunit", "deg");

		if(rad_to_deg_u (c_PhiTwist, aunits) == "deg" || rad_to_deg_u (c_theta, aunits) == "deg" || rad_to_deg_u (c_phi, aunits) == "deg" || rad_to_deg_u (c_PhiTwist, aunits) == "deg")
			aunits = "deg";

		dimension = c_PhiTwist  + "*" + aunits + " " + c_z + "*" + lunits + " "  + c_theta + "*" + aunits + " " + c_phi +  "*" + aunits + " "
							+  c_y1 + "*" + lunits + " " + c_y2 + "*" + lunits + " " + c_x1 + "*" + lunits + " "  + c_x2 + "*" + lunits + " " 
							+ c_x3 + "*" + lunits + " " +  c_x4 + "*" + lunits;

	}

	if(ee.tagName().toStdString() == "twistedtubs"){
			
		string endinnerrad = assignAttribute(ee, "endinnerrad", "0");
    		string endouterrad = assignAttribute(ee, "endouterrad", "0");
		string zlen = assignAttribute(ee, "zlen", "0");

    		string twistedangle = rad_to_deg (assignAttribute(ee, "twistedangle", "0"));
		string phi = rad_to_deg (assignAttribute(ee, "phi", "0"));
			
		string c_endinnerrad = stringify(resolve_string(endinnerrad, gconstantMap));
		string c_endouterrad = stringify(resolve_string(endouterrad, gconstantMap));
		string c_zlen = stringify(resolve_string(zlen, gconstantMap));
		string c_twistedangle = stringify(resolve_string(twistedangle, gconstantMap));
		string c_phi = stringify(resolve_string(phi, gconstantMap));

		string lunits = assignAttribute(ee, "lunit", "mm");
   		string aunits = assignAttribute(ee, "aunit", "deg");


		if(rad_to_deg_u (c_phi, aunits) == "deg" || rad_to_deg_u (c_twistedangle, aunits) == "deg" ) aunits = "deg";

		dimension = c_endinnerrad  + "*mm " + c_endouterrad + "*mm "  + c_zlen + "*mm " +c_phi  + "*deg " + c_twistedangle  + "*deg ";

	}
 

	return dimension;

	
}


