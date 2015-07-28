// gemc headers
#include "primaryTab.h"
#include "momControls.h"
#include "string_utilities.h"

// G4 headers
#include "G4ParticleTable.hh"

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;


lumiEvent::lumiEvent(goptions *Opts, QWidget *parent, string type)
{
	// Getting Luminosity Event infos from command line
	vector<string> values  = get_info(Opts->optMap["LUMI_EVENT"].args);
	int lnevents = (int) get_number(values[0]);
	double twindow = get_number(values[1]);
	double tbunch  = get_number(values[2]);
		
	if(type == "Lumi2")
	{
		values  = get_info(Opts->optMap["LUMI2_EVENT"].args);
		lnevents = (int) get_number(values[0]);
		tbunch  = get_number(values[1]);
	}
	
	QLabel *neventLabel      = new QLabel(tr("N. Particles/Event"));
	QLabel *timewindowLabel = NULL;
	if( type == "Primary") timewindowLabel = new QLabel(tr("Dummy, not displayed\n"));
	if( type == "Lumi1")   timewindowLabel = new QLabel(tr("Time Window"));
	if( type == "Lumi2")   timewindowLabel = new QLabel(tr(""));
	QLabel *bunchLabel       = new QLabel(tr("Bunch Time"));
	neventLabel->setFixedHeight(10);
	timewindowLabel->setFixedHeight(10);
	bunchLabel->setFixedHeight(10);
	QHBoxLayout *lumiHLayout = new QHBoxLayout;
	lumiHLayout->addSpacing(45);
	lumiHLayout->addWidget(neventLabel);
	lumiHLayout->addSpacing(25);
	lumiHLayout->addWidget(timewindowLabel);
	lumiHLayout->addSpacing(20);
	lumiHLayout->addWidget(bunchLabel);
	
	nevents     = new QLineEdit();
	timewindow  = new QLineEdit();
	time_bunch  = new QLineEdit();
	nevents->setMaximumWidth(100);
	timewindow->setMaximumWidth(100);
	time_bunch->setMaximumWidth(100);
	nevents->setMinimumHeight(18);
	timewindow->setMinimumHeight(18);
	time_bunch->setMinimumHeight(18);
	QHBoxLayout *lumiHLayout2 = new QHBoxLayout;
	lumiHLayout2->addWidget(nevents);
	if( type == "Lumi1") lumiHLayout2->addWidget(timewindow);
	if( type == "Lumi2") lumiHLayout2->addSpacing(135);
	lumiHLayout2->addWidget(time_bunch);
		
	QString L_ne = stringify(lnevents).c_str();
	nevents->setText(L_ne);
	
	QString L_tw = stringify(twindow/ns).c_str() + QString("*ns");
	timewindow->setText(L_tw);
	
	QString L_tb = stringify(tbunch/ns).c_str() +  QString("*ns");
	time_bunch->setText(L_tb);
	
	connect ( nevents     , SIGNAL( textChanged(const QString & ) ), parent, SLOT( changePars() ) );
	connect ( timewindow  , SIGNAL( textChanged(const QString & ) ), parent, SLOT( changePars() ) );
	connect ( time_bunch  , SIGNAL( textChanged(const QString & ) ), parent, SLOT( changePars() ) );
	
	LumiGroup = new QGroupBox(tr(""));
	QVBoxLayout *LLayout = new QVBoxLayout;
	LLayout->addLayout(lumiHLayout);
	LLayout->addLayout(lumiHLayout2);
	LLayout->addStretch(1);
	LumiGroup->setMaximumHeight(52);

	LumiGroup->setLayout(LLayout);
	
}



particleTab::particleTab(QWidget *parent, goptions *Opts, string t) : QWidget(parent)
{
	
	type = t;
	verbosity = (int) Opts->optMap["GEN_VERBOSITY"].arg;
	
	// luminosity event
	lEvent = new lumiEvent(Opts, parent, type);
	
	// momentum
	momControl = new momControls(Opts, type);
	
	connect ( momControl->beam_particle, SIGNAL( currentIndexChanged(int) )  , parent, SLOT( changePars() ) );
	connect ( momControl->momLE        , SIGNAL( textEdited(const QString &)), parent, SLOT( changePars() ) );
	connect ( momControl->rmomLE       , SIGNAL( textEdited(const QString &)), parent, SLOT( changePars() ) );
	connect ( momControl->momUnits     , SIGNAL( currentIndexChanged(int) )  , parent, SLOT( changePars() ) );
	connect ( momControl->theLE        , SIGNAL( textEdited(const QString &)), parent, SLOT( changePars() ) );
	connect ( momControl->rtheLE       , SIGNAL( textEdited(const QString &)), parent, SLOT( changePars() ) );
	connect ( momControl->theUnits     , SIGNAL( currentIndexChanged(int) )  , parent, SLOT( changePars() ) );
	connect ( momControl->phiLE        , SIGNAL( textEdited(const QString &)), parent, SLOT( changePars() ) );
	connect ( momControl->rphiLE       , SIGNAL( textEdited(const QString &)), parent, SLOT( changePars() ) );
	connect ( momControl->theUnits     , SIGNAL( currentIndexChanged(int) )  , parent, SLOT( changePars() ) );

	
	// vertex
	vtxControl = new vtxControls(Opts, type);
	
	connect ( vtxControl->vxLE     , SIGNAL( textEdited(const QString &)), parent, SLOT( changePars() ) );
	connect ( vtxControl->vyLE     , SIGNAL( textEdited(const QString &)), parent, SLOT( changePars() ) );
	connect ( vtxControl->vzLE     , SIGNAL( textEdited(const QString &)), parent, SLOT( changePars() ) );
	connect ( vtxControl->vtxUnits , SIGNAL( currentIndexChanged(int) )  , parent, SLOT( changePars() ) );
	connect ( vtxControl->rvrLE    , SIGNAL( textEdited(const QString &)), parent, SLOT( changePars() ) );
	connect ( vtxControl->rvzLE    , SIGNAL( textEdited(const QString &)), parent, SLOT( changePars() ) );

	
		
	// All together
	QVBoxLayout *mainLayout = new QVBoxLayout;
	if(type == "Lumi1" || type == "Lumi2")
		mainLayout->addWidget(lEvent->LumiGroup);
		
	mainLayout->addWidget(momControl->momentumGroup);
	mainLayout->addWidget(vtxControl->vertexGroup);
	mainLayout->addStretch(1);
	setLayout(mainLayout);
	
}



string particleTab::get_lumi_event(double t)
{
	string L_ne = qs_tostring(lEvent->nevents->text());
	string L_tw = qs_tostring(lEvent->timewindow->text());
	string L_tb = qs_tostring(lEvent->time_bunch->text());
	
	string lumi_event = L_ne + ", " + L_tw + ", " + L_tb;
	if(t != 0) lumi_event =  L_ne + ", " + stringify(t) + "*ns, " + L_tb;
	
	if(verbosity)
		cout << "  - GUI Generator Settings: " << type <<  " lumi event " << lumi_event << endl;

	return lumi_event;
}

double particleTab::get_time_window()
{
	return get_number(qs_tostring(lEvent->timewindow->text()));
}












