// gemc headers
#include "asciiField.h"
#include "utils.h"

// 1D-dipole field
// Dependent on 2 cartesian coordinates (longitudinal, transverse)
// uniform in the other coordinate
// The values can be loaded from the map in any order
// as long as their speed is ordered.
// The values are indexed as B1_2D[longi][transverse]

// from fieldFactory:  load field map
void asciiField::loadFieldMap_Dipole(gMappedField* map, double verbosity)
{
	setlocale(LC_NUMERIC, "en_US");

	int np_1     =  map->getCoordinateWithSpeed(0).np;
	int np_2     =  map->getCoordinateWithSpeed(1).np;
	double min1  =  map->getCoordinateWithSpeed(0).min;
	double min2  =  map->getCoordinateWithSpeed(1).min;
	double max1  =  map->getCoordinateWithSpeed(0).max;
	double max2  =  map->getCoordinateWithSpeed(1).max;
	double cell1 = (max1 - min1)/(np_1 - 1);
	double cell2 = (max2 - min2)/(np_2 - 1);

	// Allocate memory. [LONGI][TRANSVERSE]
	// as initialized in the map
	map->B1_2D = new double*[map->np[0]];
	for (unsigned int i = 0; i < map->np[0]; ++i)
		map->B1_2D[i] = new double[map->np[1]];


	double unit1 = get_number("1*" + map->getCoordinateWithSpeed(0).unit);
	double unit2 = get_number("1*" + map->getCoordinateWithSpeed(1).unit);

	double scale = map->scaleFactor*get_number("1*" + map->unit);

	double d1, d2, b;

	// progress bar
	int barWidth = 50;

	// using fscanf instead of c++ to read file, this is a lot faster
	char ctmp[100];
	string tmp;
	FILE *fp = fopen (map->identifier.c_str(), "r");

	// ignoring header
	while(tmp != "</mfield>")
	{
		if(fscanf(fp, "%s", ctmp) != 0)
			tmp = string(ctmp);
	}

	// now reading map
	// values as read from map
	for(int i1 = 0; i1<np_1 ; i1++)
	{
		for(int i2 = 0; i2<np_2 ; i2++)
		{
			if(fscanf(fp, "%lg %lg %lg", &d1, &d2, &b) != 0) {

				d1 *= unit1;
				d2 *= unit2;

				b  *= scale;

				if(verbosity>4 && verbosity != 99)
					cout << "  Loading Map: coordinates (" << d1 << ", " << d2 << ")   value: " << b << endl;


				// checking map consistency for first coordinate
				if( (min1  + i1*cell1 - d1)/d1 > 0.001)
				{
					cout << "   !! Error:  coordinate index wrong. Map first point should be " <<  min1  + i1*cell1
					<< " but it's  " << d1 << " instead. Cell size: " << cell1 << "  min: " << min1 << " max : " << max1 << "  index: " << i1 << endl;
				}
				// checking map consistency for second coordinate
				if( (min2  + i2*cell2 - d2)/d2 > 0.001)
				{
					cout << "   !! Error:  coordinate index wrong. Map second point should be " <<  min2  + i2*cell2
					<< " but it's  " << d2 << " instead. Cell size: " << cell2 << "  min: " << min2 << " max : " << max2 << "  index: " << i2 << endl;
				}

				// calculating index
				unsigned t1 = (unsigned) floor( ( d1 - min1 + cell1/2 ) / ( cell1 ) ) ;
				unsigned t2 = (unsigned) floor( ( d2 - min2 + cell2/2 ) / ( cell2 ) ) ;

				// The values are indexed as B1_2D[longi][transverse]
				if(   map->getCoordinateWithSpeed(0).name == "longitudinal"
					&& map->getCoordinateWithSpeed(1).name == "transverse")
					map->B1_2D[t1][t2] = b;

				if(   map->getCoordinateWithSpeed(0).name == "transverse"
					&& map->getCoordinateWithSpeed(1).name == "longitudinal")
					map->B1_2D[t2][t1] = b;
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
void gMappedField::GetFieldValue_Dipole( const double x[3], double *Bfield, int FIRST_ONLY) const
{
	double LC = 0;     	// longitudinal
	double TC = 0;     	// transverse

	if(symmetry == "dipole-x")
	{
		TC  = x[1];
		LC  = x[2];
	}
	if(symmetry == "dipole-y")
	{
		TC  = x[0];
		LC  = x[2];
	}
	if(symmetry == "dipole-z")
	{
		TC  = x[0];
		LC  = x[1];
	}

	// map indexes, bottom of the cell
	unsigned int IL = floor( ( LC - startMap[0] ) / cellSize[0] );
	unsigned int IT = floor( ( TC - startMap[1] ) / cellSize[1] );

	// checking if the point is closer to the top of the cell
	if( fabs( startMap[0] + IL*cellSize[0] - LC) > fabs( startMap[0] + (IL+1)*cellSize[0] - LC)  ) IL++;
	if( fabs( startMap[1] + IT*cellSize[1] - TC) > fabs( startMap[1] + (IT+1)*cellSize[1] - TC)  ) IT++;

	// outside map, returning no field
	if(IL>=np[0] || IT>=np[1]) return;

	// no interpolation
	if(interpolation == "none")
	{
		if(symmetry == "dipole-x") Bfield[0] = B1_2D[IL][IT];
		if(symmetry == "dipole-y") Bfield[1] = B1_2D[IL][IT];
		if(symmetry == "dipole-z") Bfield[2] = B1_2D[IL][IT];
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
		cout << "    Dipole: ";
		cout << "loc. pos. = ("    << x[0]/cm << ", " << x[1]/cm << ", " << x[2]/cm << ") cm,  ";
		cout << "IT="   << IT << "   ";
		cout << "IL="   << IL << ",  ";
		cout << "B = ("   << Bfield[0]/gauss << ",  " << Bfield[1]/gauss << ",  " << Bfield[2]/gauss << ") gauss " << endl;

	}
}






