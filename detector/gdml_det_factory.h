#ifndef GDML_DET_FACTORY_H
#define GDML_DET_FACTORY_H 1

// gemc headers
#include "detector_factory.h"

// Qt4 headers
#include <QDomDocument>
#include <QtWidgets>
#include <QString>

#include <map>
#include <string>
#include <iostream>
using namespace std;

class gphysV
{
	public:
		gphysV()
		{
			positionRef = "";
			rotationRef = "";
		}
		gphysV(string m, string s)
		{
			positionRef = m;
			rotationRef = s;		
		}

	public:
		string positionRef;
		string rotationRef;

};
//ostream &operator<<(ostream &stream,map<string, gphysV> gphysVMap );


class glogicV
{
	public:
		glogicV()
		{
			materialRef = "";
			solidRef    = "";
		}
		glogicV(string m, string s)
		{
			materialRef = m;
			solidRef    = s;		
		}

	public:
		string materialRef;
		string solidRef;

};

//ostream &operator<<(ostream &stream, map<string, glogicV> glogicVMap);

class gposition
{
	public:
		gposition(){;}
		gposition(double a, double b, double c, string l)
		{
			x    = a;
			y    = b;
			z    = c;
			unit = l;
		}
	
	public:
		double x, y, z;
		string unit;

	inline string get_dimensions()
		{
		
		return stringify(x) + "*" + unit + " " + stringify(y) + "*" + unit + " " + stringify(z) + "*" + unit; 
		}
	
	

};

//ostream &operator<<(ostream &stream, map<string, gposition> gpositionsMap);


class grotation
{
	public:
		grotation(){;}
		grotation(string a, string b, string c, string u)
		{
			x    = a;
			y    = b;
			z    = c;
			unit = u;
			
		}
	public:
		string x, y, z;
		string unit;
	
	inline string get_dimensions()
		{
		
			return x + "*" +unit + " " + y + "*" + unit + " " + z + "*" + unit; 
		}
	
};

//ostream &operator<<(ostream &stream,map<string, grotation> grotationsMap );




class gsolid
{
	public:
		gsolid()
		{
			type      = "na";
			dimension = "0";
		}
		gsolid(string t, string d)
		{

			type      = t;
			dimension = d;
		}
	
	public: 
		// dimension follow the same order 
		// as in the genat4 documentation
		// each number is multiplied by the unit
		// default units are length=mm angle=rad
		// example for box w/o units:
		// name="b100" x="10.0" y="10.0" z="10.0"
		// dimension = "10.0*mm  10.0*mm 10.0*mm"
		// example for tube 
		// z="1000.0" rmax="100.0" deltaphi="TWOPI" aunit="rad"
		// rmin is not defined so it's zero
		// lunit is not defined so it's the default (mm)
		// startphi is not defined so it's 0
		// deltaphi is defined as 2pi so I put 360 degrees
		// dimension = "0*mm 100.0*mm 1000.0*mm 0*deg 360*deg"
		
		
		string type;
		string dimension;

};

class goperation
{
	public:
		
		goperation()
		{
			first         = "";
			second        = "";
			renamed_second = "";
			pos           = "";
			rot           = "";

		}
		
		goperation(string f, string s,string u, string p, string r)
		{
			first = f;
			second = s;
			renamed_second = u;
			pos = p;
			rot = r;
		}

	public:
		string first;
		string second;
		string renamed_second;
		string pos;
		string rot;
};


//ostream &operator<<(ostream &stream,map<string, gsolid> gsolidsMap );


// you want to create a mpa<string, glogicV> and a map<string, gposition>
// where the key is the name of the objects
class gdml_det_factory : public detectorFactory
{
	public:
		gdml_det_factory(){;}
	
		// load all detectors that matches factorytype
		map<string, detector> loadDetectors();      
		
		// initialize factorytype, option and runcondition classes
		void initFactory(goptions, runConditions, string);		
		
		static detectorFactory *createFactory() 
		{
			return new gdml_det_factory; 
		}
};


string get_dimensions(QDomNode solidNode, map<string, double> gconstantMap);

string togType(string);
string rad_to_deg (string radValue);
string rad_to_deg_u (string radValue, string radUnit);
double resolve_strings(string constant_name, string v, map<string, double> gconstantMap);
//map<string, double> resolve_strings(map<string, string> stringConstantMap);
double resolve_string(string constString, map<string, double> gconstantMap);
string p_calculation(vector<string> paren_calculation);

#endif
