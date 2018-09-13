// gemc headers
#include "asciiField.h"
#include "utils.h"


// 3D field in cartesian coordinates. Field itself is in cartesian coordinates.
// Dependent on 3 cartesian coordinates (Y, X, Z) 
// The values can be loaded from the map in any order as long as their speed is ordered.
// The values are indexed as B1_3D[X][Y][Z],B2_3D[X][Y][Z],B3_3D[X][Y][Z]
// The field is three dimensional, ordered in the class as B1=Bx, B2=By, B3=Bz

//symmetry "cartesian_3D" is for full 3D map
//symmetry "cartesian_3D_quadrant" is for 3D map covering only the 1st quadrant where x>=0 && y>=0

// from fieldFactory:  load field map
void asciiField::loadFieldMap_cartesian3d(gMappedField* map, double verbosity)
{
	setlocale(LC_NUMERIC, "en_US");

	int np_1     =  map->getCoordinateWithSpeed(0).np;
	int np_2     =  map->getCoordinateWithSpeed(1).np;
	int np_3     =  map->getCoordinateWithSpeed(2).np;
	double min1  =  map->getCoordinateWithSpeed(0).min;
	double min2  =  map->getCoordinateWithSpeed(1).min;
	double min3  =  map->getCoordinateWithSpeed(2).min;
	double cell1 = (map->getCoordinateWithSpeed(0).max - min1)/(np_1 - 1);
	double cell2 = (map->getCoordinateWithSpeed(1).max - min2)/(np_2 - 1);
	double cell3 = (map->getCoordinateWithSpeed(2).max - min3)/(np_3 - 1);

	// Allocate memory. [AZI][TRANSVERSE][LONGI]
	// as initialized in the map
	map->B1_3D = new double**[map->np[0]];
	map->B2_3D = new double**[map->np[0]];
	map->B3_3D = new double**[map->np[0]];
	for (unsigned i = 0; i < map->np[0]; ++i)
	{
		map->B1_3D[i] = new double*[map->np[1]];
		map->B2_3D[i] = new double*[map->np[1]];
		map->B3_3D[i] = new double*[map->np[1]];
		for (unsigned j = 0; j < map->np[1]; ++j)
		{
			map->B1_3D[i][j] = new double[map->np[2]];
			map->B2_3D[i][j] = new double[map->np[2]];
			map->B3_3D[i][j] = new double[map->np[2]];

		}
	}

	double unit1 = get_number("1*" + map->getCoordinateWithSpeed(0).unit);
	double unit2 = get_number("1*" + map->getCoordinateWithSpeed(1).unit);
	double unit3 = get_number("1*" + map->getCoordinateWithSpeed(2).unit);

	double scale = map->scaleFactor*get_number("1*" + map->unit);

	double d1, d2, d3, b1, b2, b3;

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
			for(int i3 = 0; i3<np_3 ; i3++)
			{
				if(fscanf(fp, "%lg %lg %lg %lg %lg %lg", &d1, &d2, &d3, &b1, &b2, &b3) !=0 ) {

					d1 *= unit1;
					d2 *= unit2;
					d3 *= unit3;
					b1 *= scale;
					b2 *= scale;
					b3 *= scale;

					// checking map consistency for first coordinate
					if( (min1  + i1*cell1 - d1)/d1 > 0.001) {
						cout << "   !! Error:  coordinate index wrong. Map point should be " <<  min1  + i1*cell1
						<< " but it's  " << d1 << " instead." << endl;
					}
					// checking map consistency for second coordinate
					if( (min2  + i2*cell2 - d2)/d2 > 0.001) {
						cout << "   !! Error:  coordinate index wrong. Map point should be " <<  min2  + i2*cell2
						<< " but it's  " << d2 << " instead." << endl;
					}

					// checking map consistency for third coordinate
					if( (min3  + i3*cell3 - d3)/d3 > 0.001) {
						cout << "   !! Error:  coordinate index wrong. Map point should be " <<  min2  + i2*cell2
						<< " but it's  " << d2 << " instead." << endl;
					}


					// calculating index
					// this
					unsigned t1 = (unsigned) floor( ( d1 - min1 + cell1/2 ) / ( cell1 ) ) ;
					unsigned t2 = (unsigned) floor( ( d2 - min2 + cell2/2 ) / ( cell2 ) ) ;
					unsigned t3 = (unsigned) floor( ( d3 - min3 + cell3/2 ) / ( cell3 ) ) ;


					// The values are indexed as B1_3D[X][Y][Z]
					if(   map->getCoordinateWithSpeed(0).name == "X"
					   && map->getCoordinateWithSpeed(1).name == "Y"
					   && map->getCoordinateWithSpeed(2).name == "Z" ) {
						map->B1_3D[t1][t2][t3] = b1;
						map->B2_3D[t1][t2][t3] = b2;
						map->B3_3D[t1][t2][t3] = b3;
					} else if(   map->getCoordinateWithSpeed(0).name == "X"
							  && map->getCoordinateWithSpeed(1).name == "Z"
							  && map->getCoordinateWithSpeed(2).name == "Y" ) {
						t1 = (unsigned) floor( ( d1 - min1 + cell1/2 ) / ( cell1 ) ) ;
						t2 = (unsigned) floor( ( d3 - min3 + cell3/2 ) / ( cell3 ) ) ;
						t3 = (unsigned) floor( ( d2 - min2 + cell2/2 ) / ( cell2 ) ) ;

						map->B1_3D[t1][t3][t2] = b1;
						map->B2_3D[t1][t3][t2] = b2;
						map->B3_3D[t1][t3][t2] = b3;
					} else if(   map->getCoordinateWithSpeed(0).name == "Z"
							  && map->getCoordinateWithSpeed(1).name == "X"
							  && map->getCoordinateWithSpeed(2).name == "Y" ) {
						t1 = (unsigned) floor( ( d3 - min3 + cell3/2 ) / ( cell3 ) ) ;
						t2 = (unsigned) floor( ( d1 - min1 + cell1/2 ) / ( cell1 ) ) ;
						t3 = (unsigned) floor( ( d2 - min2 + cell2/2 ) / ( cell2 ) ) ;

						map->B1_3D[t3][t1][t2] = b1;
						map->B2_3D[t3][t1][t2] = b2;
						map->B3_3D[t3][t1][t2] = b3;
					} else if(   map->getCoordinateWithSpeed(0).name == "Z"
							  && map->getCoordinateWithSpeed(1).name == "Y"
							  && map->getCoordinateWithSpeed(2).name == "X" ) {
						t1 = (unsigned) floor( ( d3 - min3 + cell3/2 ) / ( cell3 ) ) ;
						t2 = (unsigned) floor( ( d2 - min2 + cell2/2 ) / ( cell2 ) ) ;
						t3 = (unsigned) floor( ( d1 - min1 + cell1/2 ) / ( cell1 ) ) ;

						map->B1_3D[t3][t2][t1] = b1;
						map->B2_3D[t3][t2][t1] = b2;
						map->B3_3D[t3][t2][t1] = b3;
					} else if(   map->getCoordinateWithSpeed(0).name == "Y"
							  && map->getCoordinateWithSpeed(1).name == "Z"
							  && map->getCoordinateWithSpeed(2).name == "X" ) {
						t1 = (unsigned) floor( ( d2 - min2 + cell2/2 ) / ( cell2 ) ) ;
						t2 = (unsigned) floor( ( d3 - min3 + cell3/2 ) / ( cell3 ) ) ;
						t3 = (unsigned) floor( ( d1 - min1 + cell1/2 ) / ( cell1 ) ) ;

						map->B1_3D[t2][t3][t1] = b1;
						map->B2_3D[t2][t3][t1] = b2;
						map->B3_3D[t2][t3][t1] = b3;
					} else if(   map->getCoordinateWithSpeed(0).name == "Y"
							  && map->getCoordinateWithSpeed(1).name == "X"
							  && map->getCoordinateWithSpeed(2).name == "Z" ) {
						t1 = (unsigned) floor( ( d2 - min2 + cell2/2 ) / ( cell2 ) ) ;
						t2 = (unsigned) floor( ( d1 - min1 + cell1/2 ) / ( cell1 ) ) ;
						t3 = (unsigned) floor( ( d3 - min3 + cell3/2 ) / ( cell3 ) ) ;
						
						map->B1_3D[t2][t1][t3] = b1;
						map->B2_3D[t2][t1][t3] = b2;
						map->B3_3D[t2][t1][t3] = b3;
					}

					if(verbosity>4 && verbosity != 99) {
						cout << "  Loading Map: coordinates (" << d1 << ", " << d2 << ", " << d3 << ")   values: (" << b1 << ", " << b2 << ", " << b3 << ")";
						cout << ",  indexes (t1, t2, t3) = ("  << t1 << ", " << t2 << ", " << t3 << ") " ;
						cout << ",  array sizes = ("  << np_1 << ", " << np_2 << ", " << np_3 << ") " << endl;
					}


				}
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
void gMappedField::GetFieldValue_cartesian3d( const double x[3], double *Bfield, int FIRST_ONLY) const
{

	double xx = x[0];
	double yy = x[1];
	double zz = x[2];
  
	double XX = 0;
	double YY = 0;
	double ZZ = 0;
	
	if(symmetry == "cartesian_3D"){
		XX = xx; 	  YY = yy; 	  ZZ = zz;
	}else if (symmetry == "cartesian_3D_quadrant"){
	  if (xx>=0 && yy>=0)	{ XX = xx; 	  YY = yy; 	  ZZ = zz;}
	  if (xx>=0 && yy<0)	{ XX = -yy; 	  YY = xx; 	  ZZ = zz;}
	  if (xx<0 && yy<0)	{ XX = -xx; 	  YY = -yy; 	  ZZ = zz;}
	  if (xx<0 && yy>=0)	{ XX = yy; 	  YY = -xx; 	  ZZ = zz;}	  
	  
	  if (XX<0 || YY <0) {cout << "xx " << xx << " yy " << yy << " zz " << zz << " XX " << XX << " YY " << YY << " ZZ " << ZZ << endl; return;}	  
	}	
	  
	// map indexes, bottom of the cell
	// 0 to np-1
	unsigned int IXX = floor( ( XX - startMap[0] ) / cellSize[0] );
	unsigned int IYY = floor( ( YY - startMap[1] ) / cellSize[1] );
	unsigned int IZZ = floor( ( ZZ - startMap[2] ) / cellSize[2] );	

	// outside map, returning no field
	if (XX < startMap[0] || YY < startMap[1] || ZZ < startMap[2]) return;
	if (XX >= endMap[0] || YY >= endMap[1] || ZZ >= endMap[2]) return;
	
	// no interpolation
	if(interpolation == "none")
	{
		// checking if the point is closer to the top of the cell
		if( fabs( startMap[0] + IXX*cellSize[0] - XX) > fabs( startMap[0] + (IXX+1)*cellSize[0] - XX)  ) IXX++;
		if( fabs( startMap[0] + IYY*cellSize[0] - YY) > fabs( startMap[0] + (IYY+1)*cellSize[0] - YY)  ) IYY++;
		if( fabs( startMap[0] + IZZ*cellSize[0] - ZZ) > fabs( startMap[0] + (IZZ+1)*cellSize[0] - ZZ)  ) IZZ++;
		
		Bfield[0] = B1_3D[IXX][IYY][IZZ];
		Bfield[1] = B2_3D[IXX][IYY][IZZ];
		Bfield[2] = B3_3D[IXX][IYY][IZZ];
	}
	else if (interpolation == "linear")
	{
		// relative positions within cell
		double XXR = (XX - (startMap[0] + IXX*cellSize[0])) / cellSize[0];
		double YYR = (YY - (startMap[1] + IYY*cellSize[1])) / cellSize[1];
		double ZZR = (ZZ - (startMap[2] + IZZ*cellSize[2])) / cellSize[2];		

		// field component interpolation
		double B1 = B1_3D[IXX][IYY][IZZ] * (1.0 - XXR) + B1_3D[IXX][IYY][IZZ] * XXR;
		double B2 = B2_3D[IXX][IYY][IZZ] * (1.0 - YYR) + B2_3D[IXX][IYY][IZZ] * YYR;
		double B3 = B3_3D[IXX][IYY][IZZ] * (1.0 - ZZR) + B3_3D[IXX][IYY][IZZ] * ZZR;

		Bfield[0] = B1;
		Bfield[1] = B2;
		Bfield[2] = B3;
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
		cout << "  > Track position in magnetic field map, with displacement and rotations (x,y,z)/cm:"
		<< "("  << x[0]/cm << ", "
		<< x[1]/cm << ", "
		<< x[2]/cm << ") cm,  " << endl;
		cout << "    3D XYZ: ";
		cout << "loc. pos. = ("    << x[0]/cm << ", " << x[1]/cm << ", " << x[2]/cm << ") cm,  ";
		cout << "XX="    << XX/cm << "cm,   YY=" << YY/cm << "cm, ZZ=" << ZZ/cm << ", ";
		cout << "IXX="   << IXX << "   ";
		cout << "IYY="   << IYY << "   ";
		cout << "IZZ="   << IZZ << ",  ";		
		cout << "B = ("   << Bfield[0]/gauss << ",  " << Bfield[1]/gauss << ",  " << Bfield[2]/gauss << ") gauss " << endl;
	}
}
