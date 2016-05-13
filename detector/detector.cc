// G4 headers
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Ellipsoid.hh"
#include "G4IntersectionSolid.hh"
#include "G4NistManager.hh"
#include "G4Para.hh"
#include "G4Polycone.hh"
#include "G4Polyhedra.hh"
#include "G4Torus.hh"
#include "G4Trd.hh"
#include "G4Trap.hh"
#include "G4Tubs.hh"
#include "G4CutTubs.hh"
#include "G4EllipticalTube.hh"
#include "G4Paraboloid.hh"
#include "G4Hype.hh"
#include "G4Sphere.hh"
#include "G4GenericTrap.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4PVReplica.hh"
#include "G4UnitsTable.hh"
#include "G4RotationMatrix.hh"
#include "G4TwoVector.hh"

// gemc headers
#include "detector.h"

// C++ headers
#include <string>
#include <vector>
using namespace std;

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;

detector::detector()
{
	SolidV    = NULL;
	LogicV    = NULL;
	PhysicalV = NULL;
}

int detector::create_solid(goptions gemcOpt, map<string, detector> *Map)
{
	int built = 0;
	if(SolidV) delete SolidV;
	
	string hd_msg  = gemcOpt.optMap["LOG_MSG"].args + " Solid: >> ";
	double VERB    = gemcOpt.optMap["G4P_VERBOSITY"].arg ;
	string catch_v = gemcOpt.optMap["CATCH"].args;
	
	if(type.find("ReplicaOf") != string::npos)
	{
		if(VERB>4 || name.find(catch_v) != string::npos)
			cout << hd_msg << " " << name << " is a Replica. Solid Volume will not be built." << endl;
		return 0;
	}
	
	// ####
	// Box
	// ####
	if(type == "Box")
	{
		if(dimensions.size() != 3)
		{
			cout << hd_msg << " Fatal Error: the number of dimensions for " << name
			<< " is " << dimensions.size() <<  ":" << endl;
			for(unsigned int i=0; i<dimensions.size(); i++) cout << "      dimension " << i + 1 << ": " <<  dimensions[i] << endl;
			cout << "      This does not match a G4Box. Exiting" << endl << endl;
			exit(0);
		}
		
		// checking BOX dimensions
		for(unsigned int i=0; i<dimensions.size(); i++)
			if(dimensions[i] == 0)
				cout <<   "   !!! Warning: BOX has one side null!" << endl;
		
		SolidV = new G4Box(name,             ///< name
						   dimensions[0],    ///< half length in X
						   dimensions[1],    ///< half length in Y
						   dimensions[2]);   ///< half length in Z
		
		built = 1;
	}
	
	// ##############
	// Parallelepiped
	// ##############
	if(type == "Parallelepiped")
	{
		if(dimensions.size() != 6)
		{
			cout << hd_msg << " Fatal Error: the number of dimensions for " << name
			<< " is " << dimensions.size() <<  ":" << endl;
			for(unsigned int i=0; i<dimensions.size(); i++) cout << "      dimension " << i + 1 << ": " <<  dimensions[i] << endl;
			cout << "      This does not match a G4Para. Exiting" << endl << endl;
			exit(0);
		}
		
		SolidV = new G4Para(name,            ///< name
							dimensions[0],   ///< half length in X
							dimensions[1],   ///< half length in Y
							dimensions[2],   ///< half length in Z
							dimensions[3],   ///< Angle with y axis
							dimensions[4],   ///< Polar angle of the line joining the centres of the faces at -dz and +dz in z
							dimensions[5]);  ///< Azimuthal angle of the line joining the centres of the faces at -dz and +dz in z
		
		built = 1;
	}
	
	
	// ######
	// Sphere
	// ######
	if(type == "Sphere")
	{
		if(dimensions.size() != 6)
		{
			cout << hd_msg << " Fatal Error: the number of dimensions for " << name
			<< " is " << dimensions.size() <<  ":" << endl;
			for(unsigned int i=0; i<dimensions.size(); i++) cout << "      dimension " << i + 1 << ": " <<  dimensions[i] << endl;
			cout << "      This does not match a G4Sphere. Exiting." << endl << endl;
			exit(0);
		}
		
		SolidV = new G4Sphere(name,            ///< name
							  dimensions[0],   ///< Inner radius
							  dimensions[1],   ///< Outer radius
							  dimensions[2],   ///< Starting Phi angle of the segment
							  dimensions[3],   ///< Delta Phi angle of the segment
							  dimensions[4],   ///< Starting Theta angle of the segment
							  dimensions[5]);  ///< Delta Theta angle of the segment
		
		built = 1;
	}
	
	
	// #########
	// Ellipsoid
	// #########
	if(type == "Ellipsoid")
	{
		if(dimensions.size() != 5)
		{
			cout << hd_msg << " Fatal Error: the number of dimensions for " << name
			<< " is " << dimensions.size() <<  ":" << endl;
			for(unsigned int i=0; i<dimensions.size(); i++) cout << "      dimension " << i + 1 << ": " <<  dimensions[i] << endl;
			cout << "      This does not match a G4Ellipsoid. Exiting." << endl << endl;
			exit(0);
		}
		
		SolidV = new G4Ellipsoid(name,            ///< name
								 dimensions[0],   ///< Semiaxis in X
								 dimensions[1],   ///< Semiaxis in Y
								 dimensions[2],   ///< Semiaxis in Z
								 dimensions[3],   ///< lower cut plane level, z
								 dimensions[4]);  ///< upper cut plane level, z
		
		built = 1;
	}
	
	// ##########
	// Paraboloid
	// ##########
	if(type == "Paraboloid")
	{
		if(dimensions.size() != 3)
		{
			cout << hd_msg << " Fatal Error: the number of dimensions for " << name
			<< " is " << dimensions.size() <<  ":" << endl;
			for(unsigned int i=0; i<dimensions.size(); i++) cout << "      dimension " << i + 1 << ": " <<  dimensions[i] << endl;
			cout << "      This does not match a G4Paraboloid. Exiting." << endl << endl;
			exit(0);
		}
		
		SolidV = new G4Paraboloid(name,            ///< name
								  dimensions[0],   ///< Half length Z
								  dimensions[1],   ///< Radius at -Dz
								  dimensions[2]);  ///< Radius at +Dz greater than R1
		
		built = 1;
	}
	
	
	// ##################
	// Hyperbolic Profile
	// ##################
	if(type == "Hype")
	{
		if(dimensions.size() != 5)
		{
			cout << hd_msg << " Fatal Error: the number of dimensions for " << name
			<< " is " << dimensions.size() <<  ":" << endl;
			for(unsigned int i=0; i<dimensions.size(); i++) cout << "      dimension " << i + 1 << ": " <<  dimensions[i] << endl;
			cout << "      This does not match a G4Hype. Exiting." << endl << endl;
			exit(0);
		}
		
		SolidV = new G4Hype(name,            ///< name
							dimensions[0],   ///< Inner radius
							dimensions[1],   ///< Outer radius
							dimensions[2],   ///< Inner stereo angle in radians
							dimensions[3],   ///< Outer stereo angle in radians
							dimensions[4]);  ///< Half length in Z
		
		built = 1;
	}
	
	
	
	// ####
	// Tube
	// ####
	if(type == "Tube")
	{
		if(dimensions.size() != 5)
		{
			cout << hd_msg << " Fatal Error: the number of dimensions for " << name
			<< " is " << dimensions.size() <<  ":" << endl;
			for(unsigned int i=0; i<dimensions.size(); i++) cout << "      dimension " << i + 1 << ": " <<  dimensions[i] << endl;
			cout << "      This does not match a G4Tubs. Exiting" << endl << endl;
			exit(0);
		}
		
		SolidV = new G4Tubs(name,        ///< name
							dimensions[0],   ///< Inner radius
							dimensions[1],   ///< Outer radius
							dimensions[2],   ///< Half length in z
							dimensions[3],   ///< The starting phi angle
							dimensions[4]);  ///< Delta Phi angle of the segment
		
		built = 1;
	}
	
	// ####
	// Cut Tube
	// ####
	if(type == "CTube")
	{
		if(dimensions.size() != 11)
		{
			cout << hd_msg << " Fatal Error: the number of dimensions for " << name
			<< " is " << dimensions.size() <<  ":" << endl;
			for(unsigned int i=0; i<dimensions.size(); i++) cout << "      dimension " << i + 1 << ": " <<  dimensions[i] << endl;
			cout << "      This does not match a G4CutTubs. Exiting" << endl << endl;
			exit(0);
		}

		SolidV = new G4CutTubs(name,  ///< name
				       dimensions[0],   ///< Inner radius
				       dimensions[1],   ///< Outer radius
				       dimensions[2],   ///< Half length in z
				       dimensions[3],   ///< The starting phi angle
				       dimensions[4],   ///< Delta Phi angle of the segment
						 G4ThreeVector(dimensions[5], dimensions[6], dimensions[7]),        ///< Outside Normal at -z
				       G4ThreeVector(dimensions[8], dimensions[9], dimensions[10]));      ///< Outside Normal at +z
		
		built = 1;
	}

	// ###############
	// G4ElipticalTube
	// ###############
	if(type == "EllipticalTube" || type == "Eltu")
	{
		if(dimensions.size() != 3)
		{
			cout << hd_msg << " Fatal Error: the number of dimensions for " << name
			<< " is " << dimensions.size() <<  ":" << endl;
			for(unsigned int i=0; i<dimensions.size(); i++) cout << "      dimension " << i + 1 << ": " <<  dimensions[i] << endl;
			cout << "      This does not match a G4ElipticalTube. Exiting" << endl << endl;
			exit(0);
		}
		
		SolidV = new G4EllipticalTube(name,
									  dimensions[0],  // dx: The equation of the surface in x/y is 1.0 = (x/dx)**2 +(y/dy)**2
									  dimensions[1],  // dy
									  dimensions[2]); // dz = half length in z
		built = 1;
		
	}
	
	
	
	// ####
	// Cone
	// ####
	if(type == "Cons")
	{
		if(dimensions.size() != 7)
		{
			cout << hd_msg << " Fatal Error: the number of dimensions for " << name
			<< " is " << dimensions.size() <<  ":" << endl;
			for(unsigned int i=0; i<dimensions.size(); i++) cout << "      dimension " << i + 1 << ": " <<  dimensions[i] << endl;
			cout << "      This does not match a G4Cons. Exiting" << endl << endl;
			exit(0);
		}
		
		SolidV = new G4Cons(name,            ///< name
							dimensions[0],   ///< Inner radius 1
							dimensions[1],   ///< Outer radius 1
							dimensions[2],   ///< Inner radius 2
							dimensions[3],   ///< Outer radius 2
							dimensions[4],   ///< Half length in z
							dimensions[5],   ///< The starting phi angle
							dimensions[6]);  ///< Delta Phi angle of the segment
		
		built = 1;
	}
	
	// #####
	// Torus
	// #####
	if(type == "Torus")
	{
		if(dimensions.size() != 5)
		{
			cout << hd_msg << " Fatal Error: the number of dimensions for " << name
			<< " is " << dimensions.size() <<  ":" << endl;
			for(unsigned int i=0; i<dimensions.size(); i++) cout << "      dimension " << i + 1 << ": " <<  dimensions[i] << endl;
			cout << "      This does not match a G4Torus. Exiting" << endl << endl;
			exit(0);
		}
		SolidV = new G4Torus(name,           ///< name
							 dimensions[0],   ///< Inner radius
							 dimensions[1],   ///< Outer radius
							 dimensions[2],   ///< Swept radius of torus
							 dimensions[3],   ///< pSPhi: Starting Phi angle in radians (fSPhi+fDPhi<=2PI, fSPhi>-2PI)
							 dimensions[4]);  ///< pDPhi: Delta angle of the segment in radians
		
		built = 1;
	}
	
	
	
	// #########
	// Trapezoid
	// #########
	if(type == "Trd")
	{
		if(dimensions.size() != 5)
		{
			cout << hd_msg << " Fatal Error: the number of dimensions for " << name
			<< " is " << dimensions.size() <<  ":" << endl;
			for(unsigned int i=0; i<dimensions.size(); i++) cout << "      dimension " << i + 1 << ": " <<  dimensions[i] << endl;
			cout << "      This does not match a G4Trd. Exiting" << endl << endl;
			exit(0);
		}
		SolidV = new G4Trd(name,             ///< name
						   dimensions[0],    ///< Half-length along x at the surface positioned at -dz
						   dimensions[1],    ///< Half-length along x at the surface positioned at +dz
						   dimensions[2],    ///< Half-length along y at the surface positioned at -dz
						   dimensions[3],    ///< Half-length along y at the surface positioned at +dz
						   dimensions[4]);   ///< Half-length along z axis
		built = 1;
	}
	
	
	// ##################
	// Inclined Trapezoid
	// ##################
	if(type == "ITrd")
	{
		if(dimensions.size() != 7)
		{
			cout << hd_msg << " Fatal Error: the number of dimensions for " << name
			<< " is " << dimensions.size() <<  ":" << endl;
			for(unsigned int i=0; i<dimensions.size(); i++) cout << "      dimension " << i + 1 << ": " <<  dimensions[i] << endl;
			cout << "      This does not match a ITrd. Exiting" << endl << endl;
			exit(0);
		}
		
		double alph_xz = dimensions[5];
		double alph_yz = dimensions[6];
		double x = tan(alph_xz);
		double y = tan(alph_yz);
		double r = sqrt(x*x + y*y);
		
		// Calculating the 11 G4Trap parameters
		double pDz    = dimensions[4];
		double pTheta = atan2(r, 1);
		double pPhi   = atan2(y, x);
		double pDy1   = dimensions[2];
		double pDx1   = dimensions[0];
		double pDx2   = pDx1;
		double pAlp1  = 0;
		double pDy2   = dimensions[3];
		double pDx3   = dimensions[1];
		double pDx4   = pDx3;
		double pAlp2  = 0;
		
		SolidV = new G4Trap(name,     ///< name
							pDz,      ///< Half z length
							pTheta,   ///< Polar angle of the line joining the centres of the faces at -/+pDz
							pPhi,     ///< pPhi
							pDy1,     ///< Half y length at -pDz
							pDx1,     ///< Half x length at -pDz, y=-pDy1
							pDx2,     ///< Half x length at -pDz, y=+pDy
							pAlp1,    ///< Angle with respect to the y axis from the centre of the side (lower endcap)
							pDy2,     ///< Half y length at +pDz
							pDx3,     ///< Half x length at +pDz, y=-pDy2
							pDx4,     ///< Half x length at +pDz, y=+pDy2
							pAlp2);   ///< Angle with respect to the y axis from the centre of the side (upper endcap)
		
		built = 1;
	}
	
	
	// ########################
	// Geant4 Generic Trapezoid
	// ########################
	if(type == "G4Trap")
	{
		if(dimensions.size() != 11)
		{
			cout << hd_msg << " Fatal Error: the number of dimensions for " << name
			<< " is " << dimensions.size() <<  ":" << endl;
			for(unsigned int i=0; i<dimensions.size(); i++) cout << "      dimension " << i + 1 << ": " <<  dimensions[i] << endl;
			cout << "      This does not match a G4Trd. Exiting" << endl << endl;
			exit(0);
		}
		SolidV = new G4Trap(name,            ///< name
							dimensions[0],    ///< Half z length
							dimensions[1],    ///< Polar angle of the line joining the centres of the faces at -/+pDz
							dimensions[2],    ///< pPhi
							dimensions[3],    ///< Half y length at -pDz
							dimensions[4],    ///< Half x length at -pDz, y=-pDy1
							dimensions[5],    ///< Half x length at -pDz, y=+pDy
							dimensions[6],    ///< Angle with respect to the y axis from the centre of the side (lower endcap)
							dimensions[7],    ///< Half y length at +pDz
							dimensions[8],    ///< Half x length at +pDz, y=-pDy2
							dimensions[9],    ///< Half x length at +pDz, y=+pDy2
							dimensions[10]);  ///< Angle with respect to the y axis from the centre of the side (upper endcap)
		built = 1;
	}
	
	// ##########################
	// Geant4 Arbitrary Trapezoid
	// ##########################
	if(type == "G4GenericTrap")
	{
		
		int nz = (dimensions.size() -1 ) / 2;
		double dz = dimensions[0];
		vector<G4TwoVector> vertices;
		
		for(int v=0; v<nz; v++)
			vertices.push_back(G4TwoVector(dimensions[2*v+1], dimensions[2*v+2]));
		
		SolidV = new G4GenericTrap(name,    ///< name
							         dz,    ///< Half z length
					           vertices);   ///< The (x,y) coordinates of vertices
		built = 1;
	}
	
	// #######################################
	// G4Trap Constructor with 4+4 coordinates
	// #######################################
	if(type == "G4TrapC")
	{
		if(dimensions.size() != 24)
		{
			cout << hd_msg << " Fatal Error: the number of dimensions for " << name
			<< " is " << dimensions.size() <<  ":" << endl;
			for(unsigned int i=0; i<dimensions.size(); i++) cout << "      dimension " << i + 1 << ": " <<  dimensions[i] << endl;
			cout << "      This does not match a G4Trd. Exiting" << endl << endl;
			exit(0);
		}
		
		G4ThreeVector points[8];
		for(int v=0; v<8; v++)
		{
			points[v].setX(dimensions[v*3+0]/cm * cm);
			points[v].setY(dimensions[v*3+1]/cm * cm);
			points[v].setZ(dimensions[v*3+2]/cm * cm);
		}
		
		SolidV = new G4Trap(name,        // name
							points);      // coordinates
		
		built = 1;
	}
	
	
	
	// ####################
	// Polyhedra (PGON)
	// First G4 constructor
	// WARNING: The order is NOT the same
	// as the geant4 constructor.
	// Let's revert to something that is compatible?
	// ####################
	if(type == "Pgon")
	{
		if(dimensions.size() < 8)
		{
			cout << hd_msg << " Fatal Error: the number of dimensions for " << name
			<< " is " << dimensions.size() <<  ":" << endl;
			for(unsigned int i=0; i<dimensions.size(); i++) cout << "      dimension " << i + 1 << ": " <<  dimensions[i] << endl;
			cout << "      This does not match a G4Polyhedra. Exiting" << endl << endl;
			exit(0);
		}
		int nsides  = (int) dimensions[2];  ///< number of sides
		int zplanes = (int) dimensions[3];  ///< number of planes in z directions
		if(nsides < 1)
		{
			cout << hd_msg << " Fatal Error: no sides for " << name
			<< "... should be a G4Polyhedra. Exiting" << endl << endl;
			exit(0);
		}
		double* zPlane = new double[zplanes];
		double* rInner = new double[zplanes];
		double* rOuter = new double[zplanes];
		
		for(int zpl=0; zpl<zplanes; zpl++)
		{
			rInner[zpl] = dimensions[4 + 0*zplanes + zpl] ;
			rOuter[zpl] = dimensions[4 + 1*zplanes + zpl] ;
			zPlane[zpl] = dimensions[4 + 2*zplanes + zpl] ;
		}
		SolidV = new G4Polyhedra(name,           ///< name
								 dimensions[0],  ///< Initial Phi starting angle
								 dimensions[1],  ///< Total Phi angle
								 nsides,         ///< Number of sides
								 zplanes,        ///< Number of z planes
								 zPlane,         ///< z coordinate of corners
								 rInner,         ///< Tangent distance to inner surface
								 rOuter);        ///< Tangent distance to outer surface
		built = 1;
	}
	
	
	// ####################
	// Polyhedra (PCON)
	// First G4 constructor
	// ####################
	if(type == "Polycone")
	{
		if(dimensions.size() < 7)
		{
			cout << hd_msg << " Fatal Error: the number of dimensions for " << name
			<< " is " << dimensions.size() <<  ":" << endl;
			for(unsigned int i=0; i<dimensions.size(); i++) cout << "      dimension " << i + 1 << ": " <<  dimensions[i] << endl;
			cout << "      This does not match a G4Polycone. Exiting" << endl << endl;
			exit(0);
		}
		int zplanes = (int) dimensions[2];  ///< number of planes in z directions
		
		double* zPlane = new double[zplanes];
		double* rInner = new double[zplanes];
		double* rOuter = new double[zplanes];
		
		for(int zpl=0; zpl<zplanes; zpl++)
		{
			// notice that the order should be different
			// why didn't I put z first in the constructor as in geant4?
			// might have to change it later
			rInner[zpl] = dimensions[3 + 0*zplanes + zpl] ;
			rOuter[zpl] = dimensions[3 + 1*zplanes + zpl] ;
			zPlane[zpl] = dimensions[3 + 2*zplanes + zpl] ;
		}
		SolidV = new G4Polycone(name,            ///< name
								dimensions[0],  ///< Initial Phi starting angle
								dimensions[1],  ///< Total Phi angle
								zplanes,        ///< Number of z planes
								zPlane,         ///< z coordinate of corners
								rInner,         ///< Tangent distance to inner surface
								rOuter);        ///< Tangent distance to outer surface
		built = 1;
	}
	
	
	// ############################
	// CopyPlacement:
	// Point LogicV to the original
	// ############################
	if(type.find("CopyOf") != string::npos && type.find("CopyOf") == 0)
	{
		hd_msg  = gemcOpt.optMap["LOG_MSG"].args + " Copy: >> ";
		string original(type, 6, 190);
		
		// Look for original
		map<string, detector>::iterator it = (*Map).find(TrimSpaces(original));
		if(it == (*Map).end())
		{
			cout <<  hd_msg << " <" << original << "> not found. Exiting." << endl << endl;
			exit(0);
		}
		else
		{
			if(VERB>4 || name.find(catch_v) != string::npos)
			{
				cout << hd_msg << " " << name << " is a copy of <" << TrimSpaces(original) << ">. Pointing to its logical volume." << endl;
			}
			SetLogical(it->second.GetLogical());
		}
		built = 1;
	}
	
	
	// ################
	// Solid Operations
	// ################
	if(type.find("Operation:") != string::npos && type.find("Operation:") == 0)
	{
		hd_msg  = gemcOpt.optMap["LOG_MSG"].args + " Operation: >> ";
		
		bool translationFirst = false;
		
		// If Operation:~ it will perform the translation first
		size_t posTld = type.find("~");
		if( posTld != string::npos )
		{
			translationFirst = true;
			type.replace( posTld, 1, " " );
		}
		
		//////////////////////////////////
		bool absolutecoordinates = false;
		
		// If Operation:@ it will assume that position of second object is given in the common mother volume of both objects
		
		size_t pos_at = type.find("@");
		if( pos_at != string::npos )
		{
		    absolutecoordinates = true;
		    type.replace( pos_at, 1, " " );
		}
		///////////////////////////////////////////////
		
		
		size_t posp = 0;
		size_t posm = 0;
		size_t post = 0;
		size_t pos  = 0;
		// harcoded max size of 200 here
		string operation(type, 10, 190);
		posp = operation.find("+");
		posm = operation.find("-");
		post = operation.find("*");
		if     (posp != string::npos) pos = posp;
		else if(posm != string::npos) pos = posm;
		else if(post != string::npos) pos = post;
		if(!posp && !posm && !post)
		{
			cout << hd_msg << " Operation " << operation << " for " << name << " not recognized. Exiting." << endl;
			exit(0);
		}
		
		// Locating solids
		string solid1, solid2;
		string tsolid1, tsolid2;
		solid1.assign(operation, 0,     pos);
		solid2.assign(operation, pos+1, operation.size());
		tsolid1 = TrimSpaces(solid1);
		tsolid2 = TrimSpaces(solid2);
		
		// Locating second solid transformation
		map<string, detector>::iterator it1 = (*Map).find(tsolid1);
		map<string, detector>::iterator it2 = (*Map).find(tsolid2);
		if(it1 == (*Map).end())
		{
			cout <<  hd_msg << " " << tsolid1 << " Not found. Exiting." << endl << endl;
			exit(0);
		}
		if(it2 == (*Map).end())
		{
			cout <<  hd_msg << " " << tsolid2 << " Not found. Exiting." << endl << endl;
			exit(0);
		}
		
		// Define rotational and translational transformations then combine them
		G4RotationMatrix rotate    = it2->second.rot ;
		G4ThreeVector    translate = it2->second.pos;
		G4RotationMatrix invRot    =  rotate.invert() ;
		G4Transform3D    transf1( invRot, G4ThreeVector( 0, 0, 0 ) );
		G4Transform3D    transf2( G4RotationMatrix(), translate );
		G4Transform3D    transform = transf2 * transf1 ;
		
		if( absolutecoordinates && TrimSpaces(it1->second.mother) == TrimSpaces(it2->second.mother) )
		{
			//assume that second object position and rotation are given in absolute (mother) coordinates:
			
			G4RotationMatrix invrot1 = (it1->second.rot).inverse();
			G4RotationMatrix rotate2 = it2->second.rot;
			
			G4ThreeVector net_translation = it2->second.pos - it1->second.pos;
			
			// The net rotation should be the INVERSE of the rotation of object 2 relative to object 1.
			// rotate1/2 is the rotation of object 1/2 relative to their common mother.
			// If R is the rotation of 2 relative to 1, then R2 = R * R1 --> R = R2 R1^-1
			G4RotationMatrix invnet_rotation = (rotate2 * invrot1).invert();
			
			// In order to express the relative position of object 2 in the coordinate system of object one, we must rotate it
			// applying the same rotation as that used to position object 1, according to the GEANT4 framework:
			// I do not quite understand WHY this works, but through trial and error, I have
			// discovered that the combination of operations below is what works:
			net_translation *= it1->second.rot;
			transform = G4Transform3D( invnet_rotation, net_translation );
			//We don't want there to be any possibility to overwrite "transform" in this special case, so we force translationFirst to false here:
			translationFirst = false;
		}
		
		
		// If there was tilda in the operation string then the rotation and translation are switched
		// with respect to the default behaviour of G4UnionSolid with separate rotational matrix
		// and translatin vector
		if( translationFirst )
		{
			transform = transf1 * transf2 ;
		}
		
		if(posp != string::npos)
		{
			SolidV = new G4UnionSolid(       name, it1->second.GetSolid(), it2->second.GetSolid(), transform );
		}
		if(posm != string::npos)
		{
			SolidV = new G4SubtractionSolid( name, it1->second.GetSolid(), it2->second.GetSolid(), transform );
		}
		if(post != string::npos)
		{
			SolidV = new G4IntersectionSolid(name, it1->second.GetSolid(), it2->second.GetSolid(), transform );
		}
		
		if(VERB>4 || name.find(catch_v) != string::npos)
		{
			cout << hd_msg << " " << name << " is the  " << (pos==posp ? " sum " : " difference ") << " of " << tsolid1 << " and " << tsolid2 << endl;;
		}
		built = 1;
	}
	
	
	
	
	if(VERB>4 || name.find(catch_v) != string::npos)
	{
		cout << hd_msg << " " << name << " solid " << type << " built." << endl;
	}
	
	
	if(built==0)
	{
		cout << hd_msg << " " << name << " solid >" << type << "< not recognized. Exiting." << endl;
		exit(0);
	}
	return 1;
}


int detector::create_logical_volume(map<string, G4Material*> *MMats, goptions gemcOpt)
{
	string hd_msg  = gemcOpt.optMap["LOG_MSG"].args + " Logical: >> ";
	double VERB    = gemcOpt.optMap["G4P_VERBOSITY"].arg ;
	string catch_v = gemcOpt.optMap["CATCH"].args;
	string defmat  = gemcOpt.optMap["DEFAULT_MATERIAL"].args;

	vector<aopt> changeMatOptions = gemcOpt.getArgs("SWITCH_MATERIALTO");
	for (unsigned int f = 0; f < changeMatOptions.size(); f++)
	{
		vector < string > oldNewMats = get_strings(changeMatOptions[f].args, ",");
		if(oldNewMats.size() == 2)
		{
			// oldNewMats[0] = old
			// oldNewMats[1] = new
			if(material == TrimSpaces(oldNewMats[0]))
				material = TrimSpaces(oldNewMats[1]);
		}
	}

	// don't build the logical volumes for components or replicas
	if(material == "Component" || material == "OfReplica")
	{
		if(VERB>4 || name.find(catch_v) != string::npos)
			cout << hd_msg << " " << name << " is a Solid Component or a Replicant. Logical Volume will not be built." << endl;
		return 0;
	}
	
	// Check if Material Exists
	map<string, G4Material*>::iterator i = MMats->find(material);
	
	// if material is not defined, look in G4 table.
	if(i == MMats->end() && LogicV == 0)
	{
		G4NistManager* matman = G4NistManager::Instance();
		if(matman->FindOrBuildMaterial(material)) (*MMats)[material] = matman->FindOrBuildMaterial(material);
	}
	i = MMats->find(material);
	
	// if material is still not defined, use air
	if(i == MMats->end() && LogicV == 0)
	{
		if(defmat == "none")
		{
			cout << hd_msg << " Warning: material >" << material << "< is not defined. Exiting" << endl;
			cout << hd_msg << " You can set the DEFAULT_MATERIAL flag to replace an undefined material. " << endl;
			exit(0);
		}
		else
		{
			material = defmat;
			if(MMats->find(material)== MMats->end())
			{
				cout << hd_msg << " Warning: " << defmat << " set with DEFAULT_MATERIAL is not found. Exiting" << endl;
				exit(0);
				
			}
		}
	}
	
	// Logical Volume Basic Constructor
	// If LogicV exists already, this is a copy
	if(LogicV == 0)
		LogicV = new G4LogicalVolume(SolidV, (*MMats)[material], name, 0, 0, 0, true);
	
	if(name == "root") LogicV->SetVisAttributes(G4VisAttributes::GetInvisible());
	else LogicV->SetVisAttributes(VAtts);
	
	if(VERB>4 || name.find(catch_v) != string::npos)
	{
		cout << hd_msg << " " << name << " Logical Volume built." << endl;
	}
	
	return 1;
}


int detector::create_physical_volumes(goptions gemcOpt, G4LogicalVolume *mamma)
{
	string hd_msg  = gemcOpt.optMap["LOG_MSG"].args + " Physical: >> ";
	double VERB    = gemcOpt.optMap["G4P_VERBOSITY"].arg ;
	bool   OVERL   = gemcOpt.optMap["CHECK_OVERLAPS"].arg > 0 ;
	string catch_v = gemcOpt.optMap["CATCH"].args;
	if(PhysicalV) delete PhysicalV;
	
	// don't build physical volumes for components or replicas.
	// Replicas are built in the dedicated routine
	if(material == "Component" || material == "OfReplica")
	{
		if(VERB>4 || name.find(catch_v) != string::npos)
			cout << hd_msg << " " << name << " is a Solid Component or a Replicant. Physical Volume will not be built, or replicas will be built instead." << endl;
		return 1;
	}
	
	
	if(name == "root")
		PhysicalV = new G4PVPlacement(0,       ///< rotation
									  G4ThreeVector(),   ///< translation
									  LogicV,            ///< logical volume
									  name.c_str(),      ///< name
									  0,                 ///< Mother Logical Volume
									  false,             ///< no boolean operation
									  0);                ///< copy number
	
	else
		PhysicalV = new G4PVPlacement(&rot,       ///< rotation
									  pos,                  ///< translation
									  LogicV,               ///< logical volume
									  name.c_str(),         ///< name
									  mamma,                ///< Mother Logical Volume
									  false,                ///< pMany (for future use)
									  ncopy,                ///< ncopy
									  OVERL);               ///< Checks Volume Overlapping at Placement time
	
	if(VERB>4 || name.find(catch_v) != string::npos)
	{
		if(mamma)
			cout << hd_msg << " " << name << " Physical Volume(s) built inside " << mamma->GetName() << "." << endl;
	}
	
	return 1;
}

int detector::create_replicas(goptions gemcOpt, G4LogicalVolume *mamma, detector replicant)
{
	string hd_msg  = gemcOpt.optMap["LOG_MSG"].args + " Physical: >> ";
	double VERB    = gemcOpt.optMap["G4P_VERBOSITY"].arg ;
	string catch_v = gemcOpt.optMap["CATCH"].args;
	
	// deleting the replicant detector physicsal volume as well
	if(replicant.PhysicalV) delete replicant.PhysicalV;
	
	// reset this physical volume
	if(PhysicalV) delete PhysicalV;
	
	
	if(type.find("ReplicaOf:") == 0)
	{
		if(dimensions.size() != 4)
		{
			cout << hd_msg << " Fatal Error: the number of parameters for " << name
				  << " is " << dimensions.size() <<  ":" << endl;
			for(unsigned int i=0; i<dimensions.size(); i++)
				cout << "      parameter " << i + 1 << ": " <<  dimensions[i] << endl;
			cout << "      This does not match a G4Replicas (4). Exiting" << endl << endl;
			exit(0);
		}
		
		EAxis pAxis;
		if(dimensions[0] == 1) pAxis = kXAxis;
		if(dimensions[0] == 2) pAxis = kYAxis;
		if(dimensions[0] == 3) pAxis = kZAxis;
		int nreps     = (int) get_number(dimensions[1]);
		double width  = get_number(dimensions[2]);
		double offset = get_number(dimensions[3]);
		
		// the logical volume is built
		// the mother is built as well
		PhysicalV = new G4PVReplica(name,   ///< name
											 replicant.LogicV, ///< Logical Volume to be replicated
											 mamma,            ///< Mother Logical Volume
											 pAxis,            ///< EAxis of copy - can be kXAxis,kYAxis,kZAxis,kRho,kRadial3D,kPhi
											 nreps,            ///< Number of repetitions
											 width,            ///< Width of repetitions
											 offset);          ///< Offset of repetitions
		
	}
	
	if(VERB>4 || name.find(catch_v) != string::npos)
	{
		if(mamma)
			cout << hd_msg << " " << name << " Physical Volume(s) built inside " << mamma->GetName() << "." << endl;
	}
	
	return 1;
}



// Returns dimensions nomenclature for different solid type
vector< vector<string> > dimensionstype(string solidtype)
{
	vector< vector<string> > dtypes;
	vector<string> dt;
	dt.resize(2);
	
	if(solidtype == "Box")
	{
		dt[0] = "half length in X";
		dt[1] = "Length";
		dtypes.push_back(dt);
		dt[0] = "half length in Y";
		dtypes.push_back(dt);
		dt[0] = "half length in Z";
		dtypes.push_back(dt);
	}
	
	if(solidtype == "Sphere")
	{
		dt[1] = "Length";
		dt[0] = "Inner radius";
		dtypes.push_back(dt);
		dt[0] = "Outer radius";
		dtypes.push_back(dt);
		dt[1] = "Angle";
		dt[0] = "Starting Phi angle of the segment";
		dtypes.push_back(dt);
		dt[0] = "Delta Phi angle of the segment";
		dtypes.push_back(dt);
		dt[0] = "Starting Theta angle of the segment";
		dtypes.push_back(dt);
		dt[0] = "Delta Theta angle of the segment";
		dtypes.push_back(dt);
	}
	if(solidtype == "Tube")
	{
		dt[1] = "Length";
		dt[0] = "Inner radius";
		dtypes.push_back(dt);
		dt[0] = "Outer radius";
		dtypes.push_back(dt);
		dt[0] = "Half length in z";
		dtypes.push_back(dt);
		dt[1] = "Angle";
		dt[0] = "Starting Phi angle";
		dtypes.push_back(dt);
		dt[0] = "Delta Phi angle";
		dtypes.push_back(dt);
	}
	if(solidtype == "Trd")
	{
		dt[1] = "Length";
		dt[0] = "Half-length along x at the surface at -dz";
		dtypes.push_back(dt);
		dt[0] = "Half-length along x at the surface at +dz";
		dtypes.push_back(dt);
		dt[0] = "Half-length along y at the surface at -dz";
		dtypes.push_back(dt);
		dt[0] = "Half-length along y at the surface at +dz";
		dtypes.push_back(dt);
		dt[0] = "dz: Half-length along z axis";
		dtypes.push_back(dt);
	}
	if(solidtype == "Cons")
	{
		dt[1] = "Length";
		dt[0] = "Inner radius at start";
		dtypes.push_back(dt);
		dt[0] = "Outer radius at start";
		dtypes.push_back(dt);
		dt[0] = "Inner radius at end";
		dtypes.push_back(dt);
		dt[0] = "Outer radius at end";
		dtypes.push_back(dt);
		dt[0] = "Half length in z";
		dtypes.push_back(dt);
		dt[1] = "Angle";
		dt[0] = "Starting Phi angle";
		dtypes.push_back(dt);
		dt[0] = "Delta Phi angle";
		dtypes.push_back(dt);
	}
	if(solidtype == "G4Trap")
	{
		dt[1] = "Length";
		dt[0] = "Half z length ";
		dtypes.push_back(dt);
		dt[1] = "Angle";
		dt[0] = "Polar angle of the line joining the centres of the faces";
		dtypes.push_back(dt);
		dt[0] = "Azimuthal angle of the line joining the centre of the face";
		dtypes.push_back(dt);
		dt[1] = "Length";
		dt[0] = "Half y length at -pDz";
		dtypes.push_back(dt);
		dt[0] = "Half x length of the side at y=-pDy1, -pDz";
		dtypes.push_back(dt);
		dt[0] = "Half x length of the side at y=+pDy1, -pDz";
		dtypes.push_back(dt);
		dt[1] = "Angle";
		dt[0] = "Angle to the y axis from the centre (lower endcap) ";
		dtypes.push_back(dt);
		dt[1] = "Length";
		dt[0] = "Half y length at +pDz";
		dtypes.push_back(dt);
		dt[0] = "Half x length of the side at y=-pDy2, +pDz";
		dtypes.push_back(dt);
		dt[0] = "Half x length of the side at y=+pDy2, +pDz";
		dtypes.push_back(dt);
		dt[1] = "Angle";
		dt[0] = "Angle to the y axis from the centre (upper endcap) ";
		dtypes.push_back(dt);
	}
	
	if(solidtype == "G4EllipticalTube")
	{
		dt[1]="Lenght";
		dt[0]="Half length x";
		dtypes.push_back(dt);
		dt[1]="Lenght";
		dt[0]="Half length y";
		dtypes.push_back(dt);
		dt[1]="Lenght";
		dt[0]="Half length z";
		dtypes.push_back(dt);
	}
	
	if(solidtype == "Hype")
	{
		dt[1]="Length";
		dt[0]="Inner radius";
		dtypes.push_back(dt);
		dt[0] = "Outer radius";
		dtypes.push_back(dt);
		dt[1] = "Angle";
		dt[0] = "Inner stereo angle";
		dtypes.push_back(dt);
		dt[0] = "Outer stereo angle";
		dtypes.push_back(dt);
		dt[1]="Length";
		dt[0] = "Half length in z";
		dtypes.push_back(dt);
	}
	
	if(solidtype == "Parallelepiped")
	{
		dt[1]="Length";
		dt[0]="Half length in x";
		dtypes.push_back(dt);
		dt[0] = "Half length in y";
		dtypes.push_back(dt);
		dt[0] = "Half length in z";
		dtypes.push_back(dt);
		dt[1] = "Angle";
		dt[0] = "Angle formed by the y axis and the plane joining the centre of the faces parallel to the z-x plane at -dy and +dy" ;
		dtypes.push_back(dt);
		dt[0] = "Polar angle of the line joining the centres of the faces at -dz and +dz in z ";
		dtypes.push_back(dt);
		dt[0] = "Azimuthal angle of the line joining the centres of the faces at -dz and +dz in z ";
		dtypes.push_back(dt);
	}
	
	if(solidtype == "Torus")
	{
		dt[1]="Length";
		dt[0]="Inner radius";
		dtypes.push_back(dt);
		dt[0] = "Outer radius";
		dtypes.push_back(dt);
		dt[0] = "Swept radius of torus";
		dtypes.push_back(dt);
		dt[1] = "Angle";
		dt[0] = "Starting Phi angle";
		dtypes.push_back(dt);
		dt[0] = "Outer stereo angle";
		dtypes.push_back(dt);
	}
	
	if(solidtype == "Ellipsoid")
	{
		dt[1]="Length";
		dt[0]="Semiaxis in X";
		dtypes.push_back(dt);
		dt[0] = "Semiaxis in Y ";
		dtypes.push_back(dt);
		dt[0] = "Semiaxis in Z ";
		dtypes.push_back(dt);
		dt[0] = "lower cut plane level, z";
		dtypes.push_back(dt);
		dt[0] = "upper cut plane level, z";
		dtypes.push_back(dt);
	}
	
	return dtypes;
}


ostream &operator<<(ostream &stream, detector Detector)
{
	cout  << endl;
	cout << "   Detector name:  "  << Detector.name        << "  -  " <<  Detector.description << endl;
	cout << "   Mother:  "         << Detector.mother                                 << endl;
	cout << "   Position (cm):  "  << Detector.pos/cm                                 << endl;
	cout << "   Rotation:       "  << Detector.rot                                    << endl;
	cout << "   Color:  "          << Detector.VAtts.GetColour()                      << endl;
	cout << "   Type:  "           << Detector.type                                   << endl;
	vector< vector<string> > dtypes = dimensionstype( Detector.type.c_str() );
	
	if(dtypes.size() != Detector.dimensions.size() && Detector.type.find("CopyOf") != 0)
	{
		for(unsigned int i=0; i<Detector.dimensions.size(); i++)
			cout << "   Size " << i + 1 << ":  "  << Detector.dimensions[i]  << endl;
	}
	if(dtypes.size() == Detector.dimensions.size() && Detector.type.find("CopyOf") != 0)
		for(unsigned int i=0; i<dtypes.size(); i++)
			cout << "   " << dtypes[i][0] << ":  "  << G4BestUnit(Detector.dimensions[i], dtypes[i][1])  << endl;
	
	cout << "   Material:  "       << Detector.material                               << endl;
	cout << "   Magnetic Field:  " << Detector.magfield                               << endl;
	cout << "   Copy Number:  "    << Detector.ncopy                                  << endl;
	cout << "   Activated: "       << ( Detector.exist==1 ?   "yes"   : "no" )        << endl;
	cout << "   Visible: "         << ( Detector.visible==1 ? "yes"   : "no" )        << endl;
	cout << "   Style: "           << ( Detector.style==1 ?   "solid" : "wireframe" ) << endl;
	cout << "   Sensitivity: "     << Detector.sensitivity                            << endl;
	if(Detector.sensitivity != "no")
		cout << "   hitType: "       << Detector.hitType                               << endl;
	
	if(Detector.identity.size())
		cout << Detector.identity ;
	
	cout << endl;
	
	return stream;
}

bool detector::operator == (const detector& D) const
{
	// Name uniquely identifies a volume
	if(D.name == this->name)
		return true;
	else return false;
}


























