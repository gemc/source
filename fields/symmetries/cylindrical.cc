// gemc headers
#include "asciiField.h"
#include "utils.h"

// 2D phi-symmetric cylindrical field
// Dependent on 2 cartesian coordinates (transverse, longitudinal)
// Expressed in Cylindrical coordinate
// The values can be loaded from the map in any order
// as long as their speed is ordered.
// The values are indexed as B1_2D[transverse][longi]
// The field is two dimensional, ordered in the class as B1=BT, B2=BL

// from fieldFactory:  load field map
void asciiField::loadFieldMap_Cylindrical(gMappedField* map, double verbosity)
{
	setlocale(LC_NUMERIC, "en_US");

	int np_1     =  map->getCoordinateWithSpeed(0).np;
	int np_2     =  map->getCoordinateWithSpeed(1).np;
	double min1  =  map->getCoordinateWithSpeed(0).min;
	double min2  =  map->getCoordinateWithSpeed(1).min;
	double cell1 = (map->getCoordinateWithSpeed(0).max - min1)/(np_1 - 1);
	double cell2 = (map->getCoordinateWithSpeed(1).max - min2)/(np_2 - 1);

	// Allocate memory. [LONGI][TRANSVERSE]
	// as initialized in the map
	map->B1_2D = new double*[map->np[0]];
	map->B2_2D = new double*[map->np[0]];
	for (unsigned int i = 0; i < map->np[0]; ++i)
	{
		map->B1_2D[i] = new double[map->np[1]];
		map->B2_2D[i] = new double[map->np[1]];
	}

	double unit1 = get_number("1*" + map->getCoordinateWithSpeed(0).unit);
	double unit2 = get_number("1*" + map->getCoordinateWithSpeed(1).unit);

	double scale = map->scaleFactor*get_number("1*" + map->unit);

	double d1, d2, b1, b2;

	// progress bar
	int barWidth = 50;

	// using fscanf instead of c++ to read file, this is a lot faster
	char ctmp[100];
	string tmp;
	FILE *fp = fopen (map->identifier.c_str(), "r");

	// ignoring header
	while(tmp != "</mfield>")
	{
		fscanf(fp, "%s", ctmp);
		tmp = string(ctmp);
	}

	// now reading map values
	for(int i1 = 0; i1<np_1 ; i1++)
	{
		for(int i2 = 0; i2<np_2 ; i2++)
		{
			fscanf(fp, "%lg %lg %lg %lg", &d1, &d2, &b1, &b2);

			d1 *= unit1;
			d2 *= unit2;
			b1 *= scale;
			b2 *= scale;

			if(verbosity>4 && verbosity != 99)
				cout << "  Loading Map: coordinates (" << d1 << ", " << d2 << ")   values: (" << b1 << ", " << b2 << ")." << endl;

			// checking map consistency for first coordinate
			if( (min1  + i1*cell1 - d1)/d1 > 0.001)
			{
				cout << "   !! Error: first coordinate index wrong. Map point should be " <<  min1  + i1*cell1
				<< " but it's  " << d1 << " instead."  << endl;
				cout << " min1: " << min1 << " cell1: " << cell1 << " unit1: " << unit1 << " i1: " << i1 << endl;
			}
			// checking map consistency for second coordinate
			if( (min2  + i2*cell2 - d2)/d2 > 0.001)
			{
				cout << "   !! Error: second coordinate index wrong. Map point should be " <<  min2  + i2*cell2
				<< " but it's  " << d2 << " instead." << endl;
				cout << " min2: " << min2 << " cell2: " << cell2 << " unit2: " << unit2 << " i2: " << i2 << endl;
			}

			// calculating index
			unsigned t1 = (unsigned) floor( ( d1 - min1 + cell1/2 ) / ( cell1 ) ) ;
			unsigned t2 = (unsigned) floor( ( d2 - min2 + cell2/2 ) / ( cell2 ) ) ;

			// The values are indexed as B1_2D[transverse][longi]
			if(   map->getCoordinateWithSpeed(0).name == "transverse"
				&& map->getCoordinateWithSpeed(1).name == "longitudinal")
			{
				map->B1_2D[t1][t2] = b1;
				map->B2_2D[t1][t2] = b2;
			}
			if(   map->getCoordinateWithSpeed(0).name == "longitudinal"
				&& map->getCoordinateWithSpeed(1).name == "transverse")
			{
				map->B1_2D[t2][t1] = b1;
				map->B2_2D[t2][t1] = b2;
			}
		}

		cout << "    [";
		double progress = (double)i1/(double)np_1;
		int pos = progress*barWidth;
		for (int i = 0; i < barWidth; ++i)
		{
			if      (i < pos)  cout << "=";
			else if (i == pos) cout << ">";
			else cout << " ";
		}
		cout << "] " << int(progress * 100.0) << " %\r";
		cout.flush();
	}
	fclose(fp);

	cout << endl;
}

// from mappedField:  GetFieldValue
void gMappedField::GetFieldValue_Cylindrical( const double x[3], double *Bfield, int FIRST_ONLY) const
{
	double LC  = 0;    // longitudinal
	double TC  = 0;    // transverse
	double phi = 0;    // phi angle

	// map plane is in ZX, phi on X axis
	if(symmetry == "cylindrical-z")
	{
		LC  = x[2];
		TC  = sqrt(x[0]*x[0] + x[1]*x[1]);
		phi = G4ThreeVector(x[0], x[1], x[2]).phi();
	}
	// map plane is in XY, phi on Y axis
	else if(symmetry == "cylindrical-x")
	{
		LC  = x[0];
		TC  = sqrt(x[1]*x[1] + x[2]*x[2]);
		phi = G4ThreeVector(x[2], x[0], x[1]).phi();
	}
	// map plane is in XZ, phi on Z axis
	else if(symmetry == "cylindrical-y")
	{
		LC  = x[1];
		TC  = sqrt(x[0]*x[0] + x[2]*x[2]);
		phi = G4ThreeVector(x[1], x[2], x[0]).phi();
	}

	// map indexes, bottom of the cell
	// 0 to np-1
	unsigned int IT = floor( ( TC - startMap[0] ) / cellSize[0] );
	unsigned int IL = floor( ( LC - startMap[1] ) / cellSize[1] );

	// outside map, returning no field
	if (TC < startMap[0] || LC < startMap[1]) return;

	// outside map, returning no field
	if(IT>=np[0]-1 || IL>=np[1]-1) return;


	// no interpolation
	if(interpolation == "none")
	{
		// checking if the point is closer to the top of the cell
		if( fabs( startMap[0] + IT*cellSize[0] - TC) > fabs( startMap[0] + (IT+1)*cellSize[0] - TC)  ) IT++;
		if( fabs( startMap[1] + IL*cellSize[1] - LC) > fabs( startMap[1] + (IL+1)*cellSize[1] - LC)  ) IL++;

		if(symmetry == "cylindrical-z")
		{
			Bfield[0] = B1_2D[IT][IL] * cos(phi);
			Bfield[1] = B1_2D[IT][IL] * sin(phi);
			Bfield[2] = B2_2D[IT][IL];
		}
		else if(symmetry == "cylindrical-x")
		{
			Bfield[0] = B2_2D[IT][IL];
			Bfield[1] = B1_2D[IT][IL] * cos(phi);
			Bfield[2] = B1_2D[IT][IL] * sin(phi);
		}
		else if(symmetry == "cylindrical-y")
		{
			Bfield[1] = B2_2D[IT][IL];
			Bfield[0] = B1_2D[IT][IL] * sin(phi);
			Bfield[2] = B1_2D[IT][IL] * cos(phi);
		}
	}
	else if (interpolation == "linear")
	{
		// relative positions within cell
		double xtr = (TC - (startMap[0] + IT*cellSize[0])) / cellSize[0];
		double xlr = (LC - (startMap[1] + IL*cellSize[1])) / cellSize[1];

		// first field component interpolation
		double b10 = B1_2D[IT][IL]   * (1.0 - xtr) + B1_2D[IT+1][IL]   * xtr;
		double b11 = B1_2D[IT][IL+1] * (1.0 - xtr) + B1_2D[IT+1][IL+1] * xtr;
		double b1  = b10 * (1.0 - xlr) + b11 * xlr;

		// second field component interpolation
		double b20 = B2_2D[IT][IL]   * (1.0 - xtr) + B2_2D[IT+1][IL]   * xtr;
		double b21 = B2_2D[IT][IL+1] * (1.0 - xtr) + B2_2D[IT+1][IL+1] * xtr;
		double b2 = b20 * (1.0 - xlr) + b21 * xlr;

		if(symmetry == "cylindrical-z")
		{
			Bfield[0] = b1 * cos(phi);
			Bfield[1] = b1 * sin(phi);
			Bfield[2] = b2;
		}
		else if(symmetry == "cylindrical-x")
		{
			Bfield[0] = b2;
			Bfield[1] = b1 * cos(phi);
			Bfield[2] = b1 * sin(phi);
		}
		else if(symmetry == "cylindrical-y")
		{
			Bfield[1] = b2;
			Bfield[0] = b1 * sin(phi);
			Bfield[2] = b1 * cos(phi);
		}
	}
	else
	{
		cout << "  !! Unkown field interpolation method >" << interpolation << "<" << endl;
		return;
	}

	// we don't worry about computer speed
	// if verbosity is set this high
	// so we can output units as well
	if(verbosity>3 && FIRST_ONLY != 99)
	{
		cout << "  > Track position in magnetic field: "
		<< "("  << (x[0] + mapOrigin[0])/cm << ", "
		<< (x[1] + mapOrigin[1])/cm << ", "
		<< (x[2] + mapOrigin[2])/cm << ") cm,  " << endl;
		cout << "    Cylindrical: ";
		cout << "loc. pos. = ("    << x[0]/cm << ", " << x[1]/cm << ", " << x[2]/cm << ") cm,  ";
		cout << "tr="    << TC/cm << "cm,   long=" << LC/cm << "cm, phi=" << phi/deg << ", ";
		cout << "IT="   << IT << "   ";
		cout << "IL="   << IL << ",  ";
		cout << "B = ("   << Bfield[0]/gauss << ",  " << Bfield[1]/gauss << ",  " << Bfield[2]/gauss << ") gauss " << endl;
	}
}






