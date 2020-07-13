// gemc include
#include "mappedField.h"



void gMappedField::GetFieldValue(const double point[3], double *bField) const
{
	static int FIRST_ONLY;
	
	// displacement point
	double dpoint[3] = {point[0] - mapOrigin[0], point[1] - mapOrigin[1], point[2] - mapOrigin[2]};
	
	double rpoint[3] = {dpoint[0], dpoint[1], dpoint[2]};

//	// this rotated only the coordinates but we need to rotate the fields not the coordinate
//	if(mapRotation[0] != 0) {
//		double yPrime = yRotX(rpoint);
//		double zPrime = zRotX(rpoint);
//		rpoint[1] = yPrime;
//		rpoint[2] = zPrime;
//	}
//
//	if(mapRotation[1] != 0) {
//		double xPrime = xRotY(rpoint);
//		double zPrime = zRotY(rpoint);
//		rpoint[0] = xPrime;
//		rpoint[2] = zPrime;
//	}
//
//	if(mapRotation[2] != 0) {
//		double xPrime = xRotZ(rpoint);
//		double yPrime = yRotZ(rpoint);
//		rpoint[0] = xPrime;
//		rpoint[1] = yPrime;
//	}


	bField[0] = bField[1] = bField[2] = 0;

	// dipole field
	if(symmetry == "dipole-z" || symmetry == "dipole-y" || symmetry == "dipole-x") {
		GetFieldValue_Dipole(rpoint, bField, FIRST_ONLY);
	}
	// phi-symmetric cylindrical field
	else if(symmetry == "cylindrical-z" || symmetry == "cylindrical-y" || symmetry == "cylindrical-x") {
		GetFieldValue_Cylindrical(rpoint, bField, FIRST_ONLY);
	// phi-segmented
	} else 	if(symmetry == "phi-segmented") {
		GetFieldValue_phiSegmented(rpoint, bField, FIRST_ONLY);
	// cartesian_3D
	} else 	if(symmetry == "cartesian_3D" || symmetry == "cartesian_3D_quadrant"){
		GetFieldValue_cartesian3d(rpoint, bField, FIRST_ONLY);
	}
	




	if(verbosity == 99) FIRST_ONLY = 99;
	
}

void gMappedField::RotateField( double *Bfield) const  {

	// rotating the fields 
	if(mapRotation[0] != 0) {
		double yPrime = yRotX(Bfield);
		double zPrime = zRotX(Bfield);
		Bfield[1] = yPrime;
		Bfield[2] = zPrime;
	}

	if(mapRotation[1] != 0) {
		double xPrime = xRotY(Bfield);
		double zPrime = zRotY(Bfield);
		Bfield[0] = xPrime;
		Bfield[2] = zPrime;
	}

	if(mapRotation[2] != 0) {
		double xPrime = xRotZ(Bfield);
		double yPrime = yRotZ(Bfield);
		Bfield[0] = xPrime;
		Bfield[1] = yPrime;
	}
}




gcoord gMappedField::getCoordinateWithSpeed(int speed)
{
	gcoord dummy("na", 0, 0, 0, "na", 0);
	
	for(unsigned int i=0; i<coordinates.size(); i++)
		if(coordinates[i].speed == speed) return coordinates[i];

	return dummy;
}


gcoord gMappedField::getCoordinateWithName(string name)
{
	gcoord dummy("na", 0, 0, 0, "na", 0);
	
	for(unsigned int i=0; i<coordinates.size(); i++)
		if(coordinates[i].name == name) return coordinates[i];

	return dummy;
}



void gMappedField::initializeMap()
{
	// startMap and cellSize and np class member variables
	// are ordered and stored as pointers in order to speed up GetFieldValue

	// dipole field
	// first index is longitudinal
	if(symmetry == "dipole-x" || symmetry == "dipole-y" || symmetry == "dipole-z")
	{
		startMap = new double[2];
		endMap = new double[2];		
		cellSize = new double[2];
		np       = new unsigned int[2];

		np[0]       = getCoordinateWithName("longitudinal").np;
		np[1]       = getCoordinateWithName("transverse").np;
		startMap[0] = getCoordinateWithName("longitudinal").min;
		startMap[1] = getCoordinateWithName("transverse").min;
		endMap[0] = getCoordinateWithName("longitudinal").max;
		endMap[1] = getCoordinateWithName("transverse").max;		
		cellSize[0] = (getCoordinateWithName("longitudinal").max - startMap[0]) / (np[0] - 1);
		cellSize[1] = (getCoordinateWithName("transverse").max   - startMap[1]) / (np[1] - 1);
	}
	
	// phi-symmetric cylindrical field
	// first index is transverse
	if(symmetry == "cylindrical-x" || symmetry == "cylindrical-y" || symmetry == "cylindrical-z")
	{
		startMap = new double[2];
		endMap = new double[2];				
		cellSize = new double[2];
		np       = new unsigned int[2];
		
		np[0]       = getCoordinateWithName("transverse").np;
		np[1]       = getCoordinateWithName("longitudinal").np;
		startMap[0] = getCoordinateWithName("transverse").min;
		startMap[1] = getCoordinateWithName("longitudinal").min;
		endMap[0] = getCoordinateWithName("longitudinal").max;
		endMap[1] = getCoordinateWithName("transverse").max;				
		cellSize[0] = (getCoordinateWithName("transverse").max   - startMap[0]) / (np[0] - 1);
		cellSize[1] = (getCoordinateWithName("longitudinal").max - startMap[1]) / (np[1] - 1);
	}
	
	// phi-segmented cylindrical field
	// first index is transverse
	if(symmetry == "phi-segmented")
	{
		startMap = new double[3];
		endMap = new double[3];				
		cellSize = new double[3];
		np       = new unsigned int[3];
		
		np[0]       = getCoordinateWithName("azimuthal").np;
		np[1]       = getCoordinateWithName("transverse").np;
		np[2]       = getCoordinateWithName("longitudinal").np;
		startMap[0] = getCoordinateWithName("azimuthal").min;
		startMap[1] = getCoordinateWithName("transverse").min;
		startMap[2] = getCoordinateWithName("longitudinal").min;
		endMap[0] = getCoordinateWithName("azimuthal").max;
		endMap[1] = getCoordinateWithName("transverse").max;				
		endMap[2] = getCoordinateWithName("longitudinal").max;		
		cellSize[0] = (getCoordinateWithName("azimuthal").max    - startMap[0]) / (np[0] - 1);
		cellSize[1] = (getCoordinateWithName("transverse").max   - startMap[1]) / (np[1] - 1);
		cellSize[2] = (getCoordinateWithName("longitudinal").max - startMap[2]) / (np[2] - 1);
	}
	
	//cartesian_3D
	if(symmetry == "cartesian_3D" || symmetry == "cartesian_3D_quadrant")
	{
		startMap = new double[3];
		endMap = new double[3];				
		cellSize = new double[3];
		np       = new unsigned int[3];
		
		np[0]       = getCoordinateWithName("X").np;
		np[1]       = getCoordinateWithName("Y").np;
		np[2]       = getCoordinateWithName("Z").np;
		startMap[0] = getCoordinateWithName("X").min;
		startMap[1] = getCoordinateWithName("Y").min;
		startMap[2] = getCoordinateWithName("Z").min;
		endMap[0] = getCoordinateWithName("X").max;
		endMap[1] = getCoordinateWithName("Y").max;				
		endMap[2] = getCoordinateWithName("Z").max;				
		cellSize[0] = (getCoordinateWithName("X").max    - startMap[0]) / (np[0] - 1);
		cellSize[1] = (getCoordinateWithName("Y").max   - startMap[1]) / (np[1] - 1);
		cellSize[2] = (getCoordinateWithName("Z").max - startMap[2]) / (np[2] - 1);
	}	
	
	// setting rotation sin and cosines
	sinAlpha = sin(mapRotation[0]);
	cosAlhpa = cos(mapRotation[0]);
	sinBeta = sin(mapRotation[1]);
	cosBeta = cos(mapRotation[1]);
	sinGamma = sin(mapRotation[2]);
	cosGamma = cos(mapRotation[2]);
	
}






