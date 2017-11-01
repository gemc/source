// gemc include
#include "mappedField.h"



void gMappedField::GetFieldValue( const double point[3], double *Bfield) const
{
	static int FIRST_ONLY;
	
	double Point[3];
	Point[0] = point[0] - mapOrigin[0];
	Point[1] = point[1] - mapOrigin[1];
	Point[2] = point[2] - mapOrigin[2];
	
	Bfield[0] = Bfield[1] = Bfield[2] = 0;

	// dipole field
	if(symmetry == "dipole-z" || symmetry == "dipole-y" || symmetry == "dipole-x") {
		GetFieldValue_Dipole(Point, Bfield, FIRST_ONLY);
	}
	else if(symmetry == "cylindrical-z" || symmetry == "cylindrical-y" || symmetry == "cylindrical-x") {
	// phi-symmetric cylindrical field
		GetFieldValue_Cylindrical(Point, Bfield, FIRST_ONLY);
	} else 	if(symmetry == "phi-segmented") {
	// phi-segmented
		GetFieldValue_phiSegmented(Point, Bfield, FIRST_ONLY);
	}

	if(verbosity == 99)
		FIRST_ONLY = 99;


	// cout << " CHECK FIELD " << " (" << Bfield[0]/gauss << ",  " << Bfield[1]/gauss << ",  " << Bfield[2]/gauss << ") gauss " << endl;
	
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
		cellSize = new double[2];
		np       = new unsigned int[2];

		np[0]       = getCoordinateWithName("longitudinal").np;
		np[1]       = getCoordinateWithName("transverse").np;
		startMap[0] = getCoordinateWithName("longitudinal").min;
		startMap[1] = getCoordinateWithName("transverse").min;
		cellSize[0] = (getCoordinateWithName("longitudinal").max - startMap[0]) / (np[0] - 1);
		cellSize[1] = (getCoordinateWithName("transverse").max   - startMap[1]) / (np[1] - 1);
	}
	
	// phi-symmetric cylindrical field
	// first index is transverse
	if(symmetry == "cylindrical-x" || symmetry == "cylindrical-y" || symmetry == "cylindrical-z")
	{
		startMap = new double[2];
		cellSize = new double[2];
		np       = new unsigned int[2];
		
		np[0]       = getCoordinateWithName("transverse").np;
		np[1]       = getCoordinateWithName("longitudinal").np;
		startMap[0] = getCoordinateWithName("transverse").min;
		startMap[1] = getCoordinateWithName("longitudinal").min;
		cellSize[0] = (getCoordinateWithName("transverse").max   - startMap[0]) / (np[0] - 1);
		cellSize[1] = (getCoordinateWithName("longitudinal").max - startMap[1]) / (np[1] - 1);
	}
	
	// phi-segmented cylindrical field
	// first index is transverse
	if(symmetry == "phi-segmented")
	{
		startMap = new double[3];
		cellSize = new double[3];
		np       = new unsigned int[3];
		
		np[0]       = getCoordinateWithName("azimuthal").np;
		np[1]       = getCoordinateWithName("transverse").np;
		np[2]       = getCoordinateWithName("longitudinal").np;
		startMap[0] = getCoordinateWithName("azimuthal").min;
		startMap[1] = getCoordinateWithName("transverse").min;
		startMap[2] = getCoordinateWithName("longitudinal").min;
		cellSize[0] = (getCoordinateWithName("azimuthal").max    - startMap[0]) / (np[0] - 1);
		cellSize[1] = (getCoordinateWithName("transverse").max   - startMap[1]) / (np[1] - 1);
		cellSize[2] = (getCoordinateWithName("longitudinal").max - startMap[2]) / (np[2] - 1);
	}
	
}






