// gemc headers
#include "momControls.h"
#include "string_utilities.h"

// mlibrary
#include "gstring.h"
using namespace gstring;

// G4 headers
#include "G4ParticleTable.hh"

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;

momControls::momControls(goptions *Opts, string type)
{
	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();   ///< Geant4 Particle Table
	
	// vector of string - filled from the various option
	vector<string> values;
	
	// Getting momentum from option value
	if(type == "Primary") values    = get_info(Opts->optMap["BEAM_P"].args);
	if(type == "Lumi1")   values    = get_info(Opts->optMap["LUMI_P"].args);
	if(type == "Lumi2")   values    = get_info(Opts->optMap["LUMI2_P"].args);
	string pname     = trimSpacesFromString(values[0]);

	// check that particle is in particle table
	if(!particleTable->FindParticle(pname))
		cout << "  !! GUI Control: Particle " << pname << " not found in G4 table." << endl << endl;
	
	// filling drop down menu with all particles availabe
	beam_particle = new QComboBox;
	// most important particles on top
	beam_particle->addItem(pname.c_str());
	beam_particle->addItem("e-");
	beam_particle->addItem("proton");
	beam_particle->addItem("pi+");
	beam_particle->addItem("pi-");
	beam_particle->addItem("kaon+");
	beam_particle->addItem("kaon-");
	beam_particle->addItem("mu+");
	beam_particle->addItem("mu-");
	for(int i=0; i<particleTable->entries(); i++)
		beam_particle->addItem(particleTable->GetParticleName(i).c_str());
	
	// Particle Type Layout
	QHBoxLayout *beam_particleLayout = new QHBoxLayout;
	beam_particleLayout->addSpacing(40);
	beam_particleLayout->addWidget(new QLabel(tr("Particle Type:")));
	beam_particleLayout->addWidget(beam_particle);
	beam_particleLayout->addStretch(1);
	
	
	
	// momentum
	momLE = new QLineEdit(tr("10"));
	momLE->setMaximumWidth(100);

	rmomLE = new QLineEdit(tr("0"));
	rmomLE->setMaximumWidth(100);
	
	// momentum units
	momUnits  = new QComboBox;
	momUnits->addItem(tr("GeV"));
	momUnits->addItem(tr("MeV"));
	momUnits->addItem(tr("eV"));
	momUnits->addItem(tr("KeV"));
	momUnits->addItem(tr("TeV"));

	QHBoxLayout *momLayout = new QHBoxLayout;
	momLayout->addWidget(new QLabel(tr("p: ")));
	momLayout->addWidget(momLE);
	momLayout->addWidget(new QLabel(tr(" ± ")));
	momLayout->addWidget(rmomLE);
	momLayout->addWidget(momUnits);
	momLayout->addStretch(1);

	
	// theta
	theLE = new QLineEdit(tr("0"));
	theLE->setMaximumWidth(100);
	
	rtheLE = new QLineEdit(tr("0"));
	rtheLE->setMaximumWidth(100);
	
	// theta units
	theUnits  = new QComboBox;
	theUnits->addItem(tr("deg"));
	theUnits->addItem(tr("rad"));
	
	QHBoxLayout *theLayout = new QHBoxLayout;
	theLayout->addWidget(new QLabel(tr("θ: ")));
	theLayout->addWidget(theLE);
	theLayout->addWidget(new QLabel(tr(" ± ")));
	theLayout->addWidget(rtheLE);
	theLayout->addWidget(theUnits);
	theLayout->addStretch(1);
	
	
	
	// phi
	phiLE = new QLineEdit(tr("0"));
	phiLE->setMaximumWidth(100);
	
	rphiLE = new QLineEdit(tr("0"));
	rphiLE->setMaximumWidth(100);
	
	// phi units
	phiUnits  = new QComboBox;
	phiUnits->addItem(tr("deg"));
	phiUnits->addItem(tr("rad"));
	
	QHBoxLayout *phiLayout = new QHBoxLayout;
	phiLayout->addWidget(new QLabel(tr("φ: ")));
	phiLayout->addWidget(phiLE);
	phiLayout->addWidget(new QLabel(tr(" ± ")));
	phiLayout->addWidget(rphiLE);
	phiLayout->addWidget(phiUnits);
	phiLayout->addStretch(1);

	
	// kinematic layout
	QVBoxLayout *mLayout = new QVBoxLayout;
	mLayout->addLayout(beam_particleLayout);
	mLayout->addSpacing(10);
	mLayout->addLayout(momLayout);
	mLayout->addLayout(theLayout);
	mLayout->addLayout(phiLayout);

	momentumGroup = new QGroupBox(tr("Momentum:"));
	momentumGroup->setLayout(mLayout);
	
	
//	if(type == "Primary")
//		momentumGroup->setStyleSheet(" * { background-color: rgb(220, 230, 240);} QLabel {background-color: transparent; }");
//	
//	if(type == "Lumi1")
//		momentumGroup->setStyleSheet(" * { background-color: rgb(220, 245, 235);} QLabel {background-color: transparent; }");
//	
//	if(type == "Lumi2")
//		momentumGroup->setStyleSheet(" * { background-color: rgb(195, 240, 195);} QLabel {background-color: transparent; }");
	
}





string momControls::get_momentum(double verbosity)
{
	string particle = qs_tostring(beam_particle->currentText());
	string mom  = qs_tostring(momLE->text());
	string momU = qs_tostring(momUnits->currentText());
	string the  = qs_tostring(theLE->text());
	string theU = qs_tostring(theUnits->currentText());
	string phi  = qs_tostring(phiLE->text());
	string phiU = qs_tostring(phiUnits->currentText());
	
	string beam_p = particle + ", " +
	mom + "*" + momU + ", " +
	the + "*" + theU + ", " +
	phi + "*" + phiU;
	
	if(verbosity)
		cout << "  - GUI BEAM_P: " << beam_p << endl;
	
	return beam_p;
}


void momControls::set_momentum(string momOption)
{
	vector<string> values  = get_info(momOption);
	
	if(values.size() == 4)
	{
		string pname = trimSpacesFromString(values[0]);
		
		// default values of units are from the GUI
		double mom   = get_number(values[1])/GeV;
		double theta = get_number(values[2])/deg;
		double phi   = get_number(values[3])/deg;

		beam_particle->setCurrentIndex(beam_particle->findText(QString(pname.c_str())));
		momLE->setText(QString(stringify(mom).c_str()));
		theLE->setText(QString(stringify(theta).c_str()));
		phiLE->setText(QString(stringify(phi).c_str()));
	}
}




string momControls::get_rmomentum(double verbosity)
{

	string rmom = qs_tostring(rmomLE->text());
	string momU = qs_tostring(momUnits->currentText());
	string rthe = qs_tostring(rtheLE->text());
	string theU = qs_tostring(theUnits->currentText());
	string rphi = qs_tostring(rphiLE->text());
	string phiU = qs_tostring(phiUnits->currentText());
	
	string spread_p = rmom + "*" + momU + ", " +
	rthe + "*" + theU + ", " +
	rphi + "*" + phiU;
	
	if(verbosity)
		cout << "  - GUI SPREAD_P: " << spread_p << endl;
	return spread_p;
}

void momControls::set_rmomentum(string momOption)
{
	vector<string> values  = get_info(momOption);
	
	if(values.size() == 3)
	{
		// default values of units are from the GUI
		double dmom   = get_number(values[0])/GeV;
		double dtheta = get_number(values[1])/deg;
		double dphi   = get_number(values[2])/deg;
				
		rmomLE->setText(QString(stringify(dmom).c_str()));
		rtheLE->setText(QString(stringify(dtheta).c_str()));
		rphiLE->setText(QString(stringify(dphi).c_str()));
	}
}









