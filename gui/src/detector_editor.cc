// Qt headers
#include <QtWidgets>

// G4 headers
#include "G4UnitsTable.hh"

// gemc headers
#include "detector_editor.h"
#include "images/g4Box_xpm.h"
#include "images/g4Cons_xpm.h"
#include "images/g4BREPSolidPCone_xpm.h"
#include "images/g4Sphere_xpm.h"
#include "images/g4Trd_xpm.h"
#include "images/g4Trap_xpm.h"
#include "images/g4Tubs_xpm.h"
#include "images/g4Torus_xpm.h"
#include "images/g4EllipticalTube_xpm.h"
#include "images/g4Copy_xpm.h"
#include "images/g4Operation_xpm.h"
#include "string_utilities.h"
#include "utils.h"

// C++ headers
#include <iostream>
#include <sstream>
using namespace std;

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;


void descriptionTab::update_detector(detector *Det)
{
	det = Det;

	if(det)
	{
		// detector name, description
		string d_name = "Name: '" + det->name + "', daughter of '" + det->mother + "'";
		nameLabel->setText(d_name.c_str());
		
		string d_desc = "Description: " + det->description ;
		descLabel->setText(d_desc.c_str());
		
		// position and rotation
		string xyzpos = "( " +
						stringify(det->pos.getX()/mm) + ", " +
		                stringify(det->pos.getY()/mm) + ", " +
		                stringify(det->pos.getZ()/mm) + ") mm";

		placeEdit->setText(xyzpos.c_str());
		
		
		string xyzrot = "( " +
		                stringify(det->rot.getPhi()/deg) + ", " +
		                stringify(det->rot.getTheta()/deg) + ", " +
						stringify(det->rot.getPsi()/deg) + ") deg ";
		
		rotEdit->setText(xyzrot.c_str());
		
		// solid type and all dimensions
		string typeDims = det->type.c_str();
		vector< vector<string> > dimTypes = dimensionstype(det->type);
		for(unsigned int i=0; i<dimTypes.size(); i++)
		{
			stringstream Dim;
			string d1, d2, dtot;
			Dim << G4BestUnit(det->dimensions[i], dimTypes[i][1]);
			Dim >> d1 >> d2;
			dtot = d1 + "*" + d2 ;
			typeDims = typeDims + "\n" + dimTypes[i][0] + ": " + dtot;
		}
		solidType->setText(typeDims.c_str());
		solidType->setStyleSheet("font: 12pt;");

		
		if(det->type == "Box")        solidpic = QPixmap(aBox_xpm);
		if(det->type == "Tube")       solidpic = QPixmap(aTubs_xpm);
		if(det->type == "Cons")       solidpic = QPixmap(aCons_xpm);
		if(det->type == "G4Trap")     solidpic = QPixmap(aTrap_xpm);
		if(det->type == "ITrd")       solidpic = QPixmap(aTrap_xpm);
		if(det->type == "Trd")        solidpic = QPixmap(aTrd_xpm);
		if(det->type == "Sphere")     solidpic = QPixmap(aSphere_xpm);
		if(det->type == "Torus")      solidpic = QPixmap(aTorus_xpm);
		if(det->type == "Polycone")   solidpic = QPixmap(aBREPSolidPCone_xpm);
		if(det->type == "Ellipsoid")  solidpic = QPixmap(aEllipticalTube_xpm);
		if(det->type.find("CopyOf")     != string::npos) solidpic = QPixmap(copy_xpm);
		if(det->type.find("Operation:") != string::npos) solidpic = QPixmap(operation_xpm);
		// the images are already scaled to 150 so no need to do this
		// solidPicL->setPixmap(solidpic.scaledToWidth(150.0));
		solidPicL->setPixmap(solidpic);

		
		
		string what;
		// material
		what = "Material: " + det->material;
		matLabel->setText(what.c_str());

				
		if(det->magfield == "no")
			what = "Magnetic Field: inherited" ;
		else
			what = "Magnetic Field: " + det->magfield;
		mgnLabel->setText(what.c_str());
		
		what = "Sensitivity: " + det->sensitivity + "    Hit Process: " + det->hitType;
		sensHitLabel->setText(what.c_str());

		if(det->identity.size())
		{
			what = "id: " ;
			for(unsigned int i=0; i<det->identity.size(); i++)
				what = what + det->identity[i].name + " " + stringify(det->identity[i].id) + "   ";
			
			idLabel->setText(what.c_str());
		}

		
	}
}


descriptionTab::descriptionTab(QWidget *parent):QWidget(parent)
{
	//	Detector = detect;
	nameLabel = new QLabel(tr("Placeholder for Volume Name"));
	descLabel = new QLabel(tr("Placeholder for Volume description"));
	nameLabel->setWordWrap(true);
	descLabel->setWordWrap(true);
	
	// type and dimensions
	solidType = new QLabel("placeholder for type");
	solidType->setFrameStyle(QFrame::Panel | QFrame::Sunken);
	solidType->setWordWrap(true);
	solidType->setMinimumWidth (200);
	
	solidpic = QPixmap(aBox_xpm);
	solidPicL = new QLabel();
	solidPicL->setPixmap(solidpic);

	QSplitter *splitter     = new QSplitter(Qt::Horizontal);
	
	QWidget *leftWidget = new QWidget(splitter);
	QVBoxLayout *leftLayout = new QVBoxLayout(leftWidget);
	leftLayout->addWidget(solidType);
	
	
	QWidget *rightWidget = new QWidget(splitter);
	QVBoxLayout *rightLayout = new QVBoxLayout(rightWidget);
	rightLayout->addWidget(solidPicL);

	
	// position
	QLabel *placeLabel = new QLabel("(X,Y,Z) Position");
	placeEdit = new QLineEdit(tr("placeholder for (X,Y,Z) pos"));

	// rotation
	QLabel *rotLabel = new QLabel("(phi, theta, psi) Euler rotation");
	rotEdit = new QLineEdit("placeholder for (X,Y,Z) rot");
		

	matLabel     = new QLabel(tr("Placeholder for material "));
	mgnLabel     = new QLabel(tr("Placeholder for magnetic field"));
	sensHitLabel = new QLabel(tr("Placeholder for sensitivity, Hit process"));
	idLabel      = new QLabel(tr("Placeholder for identifier"));
	idLabel->setWordWrap(true);
	idLabel->setMaximumWidth (300);
	
	
	
	// all layouts
	QVBoxLayout *mainLayout = new QVBoxLayout;
	
	mainLayout->addWidget(nameLabel);
	mainLayout->addWidget(descLabel);
	mainLayout->addWidget(splitter);
	
	mainLayout->addWidget(placeLabel);
	mainLayout->addWidget(placeEdit);
	
	mainLayout->addWidget(rotLabel);
	mainLayout->addWidget(rotEdit);

	mainLayout->addWidget(matLabel);
	mainLayout->addWidget(mgnLabel);
	mainLayout->addWidget(sensHitLabel);
	mainLayout->addWidget(idLabel);
	
	mainLayout->addStretch(10);
	
	setLayout(mainLayout);
	
	// dont connect these just yet, coordinate not set yet?
//	connect ( placeXEdit , SIGNAL( textChanged (QString) ),  this, SLOT( change_placement() ) );
//	connect ( placeYEdit , SIGNAL( textChanged (QString) ),  this, SLOT( change_placement() ) );
//	connect ( placeZEdit , SIGNAL( textChanged (QString) ),  this, SLOT( change_placement() ) );	
}



void descriptionTab::change_dimension()
{
//	vector<double> newdimensions;
//	vector< vector<string> > dimTypes = dimensionstype(det->type);
//	for(unsigned int i=0; i<dimTypes.size(); i++)
//		newdimensions.push_back(get_number(qs_tostring(dimTypesEdit[i]->text())));
//	det->dimensions = newdimensions;
	
}


void descriptionTab::change_placement()
{
	
//	double x = get_number(qs_tostring(placeXEdit->text()));
//	double y = get_number(qs_tostring(placeYEdit->text()));
//	double z = get_number(qs_tostring(placeZEdit->text()));
		
//	Detector->pos = G4ThreeVector(x, y, z);
}












