// gemc headers
#include "vtxControls.h"
#include "string_utilities.h"

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;


vtxControls::vtxControls(goptions *Opts, string type)
{	
	// vector of string - filled from the various option
	vector<string> values;
	string units;

	
	// x
	vxLE = new QLineEdit(tr("0"));
	vxLE->setMaximumWidth(50);
	
	// Δr
	rvrLE = new QLineEdit(tr("0"));
	rvrLE->setMaximumWidth(50);

	
	QHBoxLayout *vxLayout = new QHBoxLayout;
	vxLayout->addWidget(new QLabel(tr("vX: ")));
	vxLayout->addWidget(vxLE);
	vxLayout->addWidget(new QLabel(tr("Δr: ")));
	vxLayout->addWidget(rvrLE);
	vxLayout->addStretch(1);

	
	// y
	vyLE = new QLineEdit(tr("0"));
	vyLE->setMaximumWidth(50);

	// Δz
	rvzLE = new QLineEdit(tr("0"));
	rvzLE->setMaximumWidth(50);
	
	
	QHBoxLayout *vyLayout = new QHBoxLayout;
	vyLayout->addWidget(new QLabel(tr("vY: ")));
	vyLayout->addWidget(vyLE);
	vyLayout->addWidget(new QLabel(tr("Δz: ")));
	vyLayout->addWidget(rvzLE);
	vyLayout->addStretch(1);

	// z
	vzLE = new QLineEdit(tr("0"));
	vzLE->setMaximumWidth(50);
	
	//  units
	vtxUnits  = new QComboBox;
	vtxUnits->addItem(tr("cm"));
	vtxUnits->addItem(tr("mm"));
	vtxUnits->addItem(tr("m"));
	vtxUnits->setMaximumWidth(55);

	QHBoxLayout *vzLayout = new QHBoxLayout;
	vzLayout->addWidget(new QLabel(tr("vZ: ")));
	vzLayout->addWidget(vzLE);
	vzLayout->addWidget(new QLabel(tr("Units: ")));
	vzLayout->addWidget(vtxUnits);
	vzLayout->addStretch(1);

		


	QHBoxLayout *rvzLayout = new QHBoxLayout;
	rvzLayout->addStretch(1);

	
	
	// kinematic layout
	QVBoxLayout *vLayout = new QVBoxLayout;
	vLayout->addLayout(vxLayout);
	vLayout->addLayout(vyLayout);
	vLayout->addLayout(vzLayout);

	vertexGroup = new QGroupBox(tr("Vertex"));
	vertexGroup->setLayout(vLayout);
	
	
//	if(type == "Primary")
//		vertexGroup->setStyleSheet(" * { background-color: rgb(220, 230, 240);} QLabel {background-color: transparent; }");
//	
//	if(type == "Lumi1")
//		vertexGroup->setStyleSheet(" * { background-color: rgb(220, 245, 235);} QLabel {background-color: transparent; }");
//	
//	if(type == "Lumi2")
//		vertexGroup->setStyleSheet(" * { background-color: rgb(195, 240, 195);} QLabel {background-color: transparent; }");
	
}





string vtxControls::get_vertex(double verbosity)
{
	string vx = qs_tostring(vxLE->text());
	string vy = qs_tostring(vyLE->text());
	string vz = qs_tostring(vzLE->text());
	
	string vU = qs_tostring(vtxUnits->currentText());
	
	string beam_v = "( " + vx + ", " + vy + ", " + vz + ") " + vU;
	
	if(verbosity)
		cout << "  - GUI BEAM_V: " << beam_v << endl;
	
	return beam_v;
}


void vtxControls::set_vertex(string vtxOption)
{
	vector<string> values  = get_info(vtxOption);
	
	if(values.size() == 4)
	{
		string 		units = TrimSpaces(values[3]);

		// default values of units are from the GUI
		double vx = get_number(values[0] + "*" + units)/cm;
		double vy = get_number(values[1] + "*" + units)/cm;
		double vz = get_number(values[2] + "*" + units)/cm;
		
		vxLE->setText(QString(stringify(vx).c_str()));
		vyLE->setText(QString(stringify(vy).c_str()));
		vzLE->setText(QString(stringify(vz).c_str()));
	}
}


string vtxControls::get_rvertex(double verbosity)
{
	string dr = qs_tostring(rvrLE->text());
	string dz = qs_tostring(rvzLE->text());
	
	string vU = qs_tostring(vtxUnits->currentText());
	
	string spread_v = "( " + dr  + ", " + dz + ") " + vU;
	
	if(verbosity)
		cout << "  - GUI SPREAD_V: " << spread_v << endl;
	
	return spread_v;
}

void vtxControls::set_rvertex(string vtxOption)
{
	vector<string> values  = get_info(vtxOption);
	
	if(values.size() == 3)
	{
		string 		units = TrimSpaces(values[2]);
		
		// default values of units are from the GUI
		double dvr = get_number(values[0] + "*" + units)/cm;
		double dvz = get_number(values[1] + "*" + units)/cm;
		
		rvrLE->setText(QString(stringify(dvr).c_str()));
		rvzLE->setText(QString(stringify(dvz).c_str()));
	}
}






