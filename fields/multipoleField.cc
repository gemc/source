// gemc headers
#include "multipoleField.h"

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;

multipoleField::multipoleField(int Npole, G4double scale, G4double x,
		G4double y, G4double z, G4double rot, string ROTaxis)
{
	polenumber = Npole;
	strength   = scale;
	origin[0]  = x;
	origin[1]  = y;
	origin[2]  = z;
	rotation   = rot;
	rotaxis    = ROTaxis;

	if (rotaxis != "X" && rotaxis != "Y" && rotaxis != "Z")
	{
		cout	<< "!!! Error: multipole field has rot axis along X or Y or Z, while you have axis "
				<< rotaxis << endl; exit(-1);
	}

}

multipoleField::~multipoleField() {}

// ------------------------------------------------------------------------

void multipoleField::GetFieldValue(const G4double pos[4], G4double *B) const
{

	G4ThreeVector x0(pos[0], pos[1], pos[2]);
	G4ThreeVector x1(origin[0], origin[1], origin[2]);
	G4ThreeVector x2;
	G4ThreeVector x0_local;
	if (rotaxis=="X")
	{
		x2 = G4ThreeVector(0*cm,0*cm,1*cm).rotateX(rotation)+x1;
		x0_local = (x0 - x1).rotateX(-rotation);
	}
	else if (rotaxis=="Y")
	{
		x2 = G4ThreeVector(0*cm,0*cm,1*cm).rotateY(rotation)+x1;
		x0_local = (x0 - x1).rotateY(-rotation);
	}
	else if (rotaxis=="Z")
	{
		x2 = G4ThreeVector(0*cm,0*cm,1*cm).rotateZ(rotation)+x1;
		x0_local = (x0 - x1).rotateZ(-rotation);
	}
	else
	{
		cout	<< "!!! Error: multipole field has rot axis along X or Y or Z, while you have axis: "
				<< rotaxis << endl;
		exit(-1);
	}

	G4double r = (x2 - x1).cross(x1 - x0).mag() / (x2 - x1).mag(); //distance from x0 to line x1-x2
	G4double phi = atan2(x0_local.y(), x0_local.x());

	G4ThreeVector B_local;
	if (polenumber == 2)
	{
		B_local.setX(0);
		B_local.setY(strength);
		B_local.setZ(0);
	}
	else
	{
		int a = polenumber / 2 - 1;
		B_local.setX(strength * pow(r/m, a) * sin(a * phi));
		B_local.setY(strength * pow(r/m, a) * cos(a * phi));
		B_local.setZ(0);
	}

	G4ThreeVector B_lab=B_local;
	if (rotaxis=="X")		{B_lab.rotateX(rotation);}
	else if (rotaxis=="Y")	{B_lab.rotateY(rotation);}
	else if (rotaxis=="Z")	{B_lab.rotateZ(rotation);}
	else {cout	<< "!!! Error: multipole field has rot axis along X or Y or Z, while you have axis "
				<< rotaxis << endl; exit(-1);
	}

	B[0]=B_lab.x();
	B[1]=B_lab.y();
	B[2]=B_lab.z();

}
