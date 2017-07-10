// gemc headers
#include "asciiField.h"
#include "utils.h"


// phi-segmented 3D field in cylindrical coordinates. Field itself is in cartesian coordinates.
// Dependent on 3 cartesian coordinates (transverse, azimuthal, longitudinal) expressed in cylindrical coordinate
// The values can be loaded from the map in any order as long as their speed is ordered.
// The values are indexed as B1_3D[transverse][longi]
// The field is two dimensional, ordered in the class as B1=BT, B2=BL

// from fieldFactory:  load field map
void asciiField::loadFieldMap_phiSegmented(gMappedField* map, double verbosity)
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


					// The values are indexed as B1_3D[AZI][TRANSVERSE][LONGI]
					if(   map->getCoordinateWithSpeed(0).name == "azimuthal"
					   && map->getCoordinateWithSpeed(1).name == "transverse"
					   && map->getCoordinateWithSpeed(2).name == "longitudinal" ) {
						map->B1_3D[t1][t2][t3] = b1;
						map->B2_3D[t1][t2][t3] = b2;
						map->B3_3D[t1][t2][t3] = b3;
					} else if(   map->getCoordinateWithSpeed(0).name == "azimuthal"
							  && map->getCoordinateWithSpeed(1).name == "longitudinal"
							  && map->getCoordinateWithSpeed(2).name == "transverse" ) {
						t1 = (unsigned) floor( ( d1 - min1 + cell1/2 ) / ( cell1 ) ) ;
						t2 = (unsigned) floor( ( d3 - min3 + cell3/2 ) / ( cell3 ) ) ;
						t3 = (unsigned) floor( ( d2 - min2 + cell2/2 ) / ( cell2 ) ) ;

						map->B1_3D[t1][t3][t2] = b1;
						map->B2_3D[t1][t3][t2] = b2;
						map->B3_3D[t1][t3][t2] = b3;
					} else if(   map->getCoordinateWithSpeed(0).name == "longitudinal"
							  && map->getCoordinateWithSpeed(1).name == "azimuthal"
							  && map->getCoordinateWithSpeed(2).name == "transverse" ) {
						t1 = (unsigned) floor( ( d3 - min3 + cell3/2 ) / ( cell3 ) ) ;
						t2 = (unsigned) floor( ( d1 - min1 + cell1/2 ) / ( cell1 ) ) ;
						t3 = (unsigned) floor( ( d2 - min2 + cell2/2 ) / ( cell2 ) ) ;

						map->B1_3D[t3][t1][t2] = b1;
						map->B2_3D[t3][t1][t2] = b2;
						map->B3_3D[t3][t1][t2] = b3;
					} else if(   map->getCoordinateWithSpeed(0).name == "longitudinal"
							  && map->getCoordinateWithSpeed(1).name == "transverse"
							  && map->getCoordinateWithSpeed(2).name == "azimuthal" ) {
						t1 = (unsigned) floor( ( d3 - min3 + cell3/2 ) / ( cell3 ) ) ;
						t2 = (unsigned) floor( ( d2 - min2 + cell2/2 ) / ( cell2 ) ) ;
						t3 = (unsigned) floor( ( d1 - min1 + cell1/2 ) / ( cell1 ) ) ;

						map->B1_3D[t3][t2][t1] = b1;
						map->B2_3D[t3][t2][t1] = b2;
						map->B3_3D[t3][t2][t1] = b3;
					} else if(   map->getCoordinateWithSpeed(0).name == "transverse"
							  && map->getCoordinateWithSpeed(1).name == "longitudinal"
							  && map->getCoordinateWithSpeed(2).name == "azimuthal" ) {
						t1 = (unsigned) floor( ( d2 - min2 + cell2/2 ) / ( cell2 ) ) ;
						t2 = (unsigned) floor( ( d3 - min3 + cell3/2 ) / ( cell3 ) ) ;
						t3 = (unsigned) floor( ( d1 - min1 + cell1/2 ) / ( cell1 ) ) ;

						map->B1_3D[t2][t3][t1] = b1;
						map->B2_3D[t2][t3][t1] = b2;
						map->B3_3D[t2][t3][t1] = b3;
					} else if(   map->getCoordinateWithSpeed(0).name == "transverse"
							  && map->getCoordinateWithSpeed(1).name == "azimuthal"
							  && map->getCoordinateWithSpeed(2).name == "longitudinal" ) {
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
void gMappedField::GetFieldValue_phiSegmented( const double x[3], double *Bfield, int FIRST_ONLY) const
{

	// this routine will rotate the coordinate to segment local
	// coordinates to find the map indexes.
	// Then it will rotate the field back to the original position

	double mfield[3] = {0,0,0};

	double aC = 0;  // azimuthal
	double tC = 0;  // transverse
	double lC = 0;  // longitudinal

	double aLC;            	///< phi local  to the first segment
	unsigned aI, tI, lI;  	///< indexes of the vector map


	aC = atan2( x[1], x[0] )*rad;       ///< phi in radians
	if( aC < 0 ) aC += 360*deg;

	tC = sqrt(x[0]*x[0] + x[1]*x[1]); 	///< R
	lC = x[2];                          ///< Z

	// Rotating the point to within the map limit
	int segment = 1;
	aLC = aC;
	while (aLC/deg > 30)
	{
		aLC -= 60*deg;
		segment++;
	}

	// need the delta phi from the local coordinate, so we can rotate the field back to the original point
	double dphi = aC - aLC;

	// for the map index we need the absolute value
	double aaLC = fabs(aLC);

	// map indexes, bottom of the cell
	// 0 to np-1
	aI = floor( ( aaLC - startMap[0] ) / cellSize[0] );
	tI = floor( ( tC   - startMap[1] ) / cellSize[1] );
	lI = floor( ( lC   - startMap[2] ) / cellSize[2] );

	// checking if the point is closer to the top of the cell
	if( fabs( startMap[0] + aI*cellSize[0] - aaLC) > fabs( startMap[0] + (aI + 1)*cellSize[0] - aaLC)  ) aI++;
	if( fabs( startMap[1] + tI*cellSize[1] - tC)   > fabs( startMap[1] + (tI + 1)*cellSize[1] - tC)    ) tI++;
	if( fabs( startMap[2] + lI*cellSize[2] - lC)   > fabs( startMap[2] + (lI + 1)*cellSize[2] - lC)    ) lI++;


	// outside map, returning no field
	if (aC < startMap[0] || tC < startMap[1] || lC < startMap[2]) return;


	// outside map, returning no field
	if(aI >= np[0]-1 || tI >= np[1]-1 || lI >= np[2]-1) return;

	// positive on the right side of the segment
	int sign = (aLC >= 0 ? 1 : -1);


	// no interpolation
	if(interpolation == "none")
	{
		// checking if the point is closer to the top of the cell
		if( fabs( startMap[0] + aI*cellSize[0] - aaLC) > fabs( startMap[0] + (aI+1)*cellSize[0] - aaLC)  ) aI++;
		if( fabs( startMap[1] + tI*cellSize[1] - tC)   > fabs( startMap[1] + (tI+1)*cellSize[1] - tC)    ) tI++;
		if( fabs( startMap[1] + lI*cellSize[1] - lC)   > fabs( startMap[1] + (lI+1)*cellSize[1] - lC)    ) lI++;

		// Field at local point
		mfield[0] = B1_3D[aI][tI][lI];
		mfield[1] = B2_3D[aI][tI][lI];
		mfield[2] = B3_3D[aI][tI][lI];
	}

	else if (interpolation == "linear")
	{
		// relative positions within cell
		double xaz = (aaLC - (startMap[0] + aI*cellSize[0])) / cellSize[0];
		double xtr = (tC   - (startMap[1] + tI*cellSize[1])) / cellSize[1];
		double xlr = (lC   - (startMap[2] + lI*cellSize[2])) / cellSize[2];

		// first field component interpolation

		// first, interpolate along azimutal
		double b100 = B1_3D[aI][tI][lI]     * (1-xaz) + B1_3D[aI+1][tI][lI]     * xaz;
		double b101 = B1_3D[aI][tI][lI+1]   * (1-xaz) + B1_3D[aI+1][tI][lI+1]   * xaz;
		double b110 = B1_3D[aI][tI+1][lI]   * (1-xaz) + B1_3D[aI+1][tI+1][lI]   * xaz;
		double b111 = B1_3D[aI][tI+1][lI+1] * (1-xaz) + B1_3D[aI+1][tI+1][lI+1] * xaz;

		// second, interpolate along transverse
		double b10 = b100 * (1 - xtr) + b110*xtr;
		double b11 = b101 * (1 - xtr) + b111*xtr;

		// finally interpolate along longitudinal
		mfield[0] = b10 * (1 - xlr) + b11 * xlr;


		// second field component interpolation

		// first, interpolate along azimutal
		double b200 = B2_3D[aI][tI][lI]     * (1-xaz) + B2_3D[aI+1][tI][lI]     * xaz;
		double b201 = B2_3D[aI][tI][lI+1]   * (1-xaz) + B2_3D[aI+1][tI][lI+1]   * xaz;
		double b210 = B2_3D[aI][tI+1][lI]   * (1-xaz) + B2_3D[aI+1][tI+1][lI]   * xaz;
		double b211 = B2_3D[aI][tI+1][lI+1] * (1-xaz) + B2_3D[aI+1][tI+1][lI+1] * xaz;

		// second, interpolate along transverse
		double b20 = b200 * (1 - xtr) + b210*xtr;
		double b21 = b201 * (1 - xtr) + b211*xtr;

		// finally interpolate along longitudinal
		mfield[1] = b20 * (1 - xlr) + b21 * xlr;

		// third field component interpolation

		// first, interpolate along azimutal
		double b300 = B3_3D[aI][tI][lI]     * (1-xaz) + B3_3D[aI+1][tI][lI]     * xaz;
		double b301 = B3_3D[aI][tI][lI+1]   * (1-xaz) + B3_3D[aI+1][tI][lI+1]   * xaz;
		double b310 = B3_3D[aI][tI+1][lI]   * (1-xaz) + B3_3D[aI+1][tI+1][lI]   * xaz;
		double b311 = B3_3D[aI][tI+1][lI+1] * (1-xaz) + B3_3D[aI+1][tI+1][lI+1] * xaz;

		// second, interpolate along transverse
		double b30 = b300 * (1 - xtr) + b310*xtr;
		double b31 = b301 * (1 - xtr) + b311*xtr;

		// finally interpolate along longitudinal
		mfield[2] = b30 * (1 - xlr) + b31 * xlr;
	}
	else
	{
		cout << "  !! Unkown field interpolation method >" << interpolation << "<" << endl;
		return;
	}

	// Rotating the field back to original point
	Bfield[0] =  sign*mfield[0] * cos(dphi/rad) - mfield[1] * sin(dphi/rad);
	Bfield[1] =  sign*mfield[0] * sin(dphi/rad) + mfield[1] * cos(dphi/rad);
	Bfield[2] =  sign*mfield[2];

	if(verbosity>3 && FIRST_ONLY != 99)
	{

		cout << "  > Track position in magnetic field: "
		<< "("  << (x[0] + mapOrigin[0])/cm << ", "
		<< (x[1] + mapOrigin[1])/cm << ", "
		<< (x[2] + mapOrigin[2])/cm << ") cm,  " << endl;

		cout << "  > Track position in magnetic field (phi, r, z): "
		<< "("  << aC/deg << " deg, "
		<< tC/cm  << " cm, "
		<< lC/cm << " cm)  " << endl;

		cout << "  > Local coordinates (phi, segment, delta phi) "
		<< "("  << aaLC/deg << " deg, "
		<< segment  << ", "
		<< dphi/deg << " deg)  " << endl;

		cout << "  > Map Indexes (phi, r, z) "
		<< "("  << aI << ", "
		<< tI << ", "
		<< lI << ")  " << endl;

		cout << "  > Field Values (local) (Bx, By, Bz) "
		<< "("  << mfield[0]/tesla << ", "
		<< mfield[1]/tesla << ", "
		<< mfield[2]/tesla << ") tesla " << endl;

		cout << "  > Field Values (absolute) (Bx, By, Bz) "
		<< "("  << Bfield[0]/tesla << ", "
		<< Bfield[1]/tesla << ", "
		<< Bfield[2]/tesla << ") tesla " << endl;

		cout << endl;
	}
	
}






