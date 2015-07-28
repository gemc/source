// G4 headers
#include "G4ParticleTable.hh"

// Qt headers
#include <QtWidgets>

// gemc headers
#include "detector.h"
#include "run_control.h"
#include "string_utilities.h"

// C++ headers
#include <ctime>

run_control::run_control(QWidget *parent, goptions *Opts) : QWidget(parent)
{
	gemcOpt   = Opts;
	verbosity  = gemcOpt->optMap["LOG_VERBOSITY"].arg ;


	pbeamtab  = new particleTab(this, gemcOpt, "Primary");
	lbeamtab  = new particleTab(this, gemcOpt, "Lumi1");
	lbeamtab2 = new particleTab(this, gemcOpt, "Lumi2");
	
	
	QTabWidget* BeamType = new QTabWidget;
	BeamType->addTab(pbeamtab,  tr("Generator"));
	BeamType->addTab(lbeamtab,  tr("Beam 1"));
	BeamType->addTab(lbeamtab2, tr("Beam 2"));
		
	BeamType->setStyleSheet(QString("QTabBar::tab:selected { background-color:  rgb(124,  164, 210); } "));

	QVBoxLayout *mainLayout = new QVBoxLayout;
	mainLayout->addWidget(BeamType);
	setLayout(mainLayout);

	// Momentum
	pbeamtab->momControl->set_momentum(gemcOpt->optMap["BEAM_P"].args);
	pbeamtab->momControl->set_rmomentum(gemcOpt->optMap["SPREAD_P"].args);
		
	lbeamtab->momControl->set_momentum(gemcOpt->optMap["LUMI_P"].args);
	lbeamtab->momControl->set_rmomentum(gemcOpt->optMap["LUMI_SPREAD_P"].args);
	
	lbeamtab2->momControl->set_momentum(gemcOpt->optMap["LUMI2_P"].args);
	lbeamtab2->momControl->set_rmomentum(gemcOpt->optMap["LUMI2_SPREAD_P"].args);

	
	// vertex
	pbeamtab->vtxControl->set_vertex(gemcOpt->optMap["BEAM_V"].args);
	pbeamtab->vtxControl->set_rvertex(gemcOpt->optMap["SPREAD_V"].args);
	
	lbeamtab->vtxControl->set_vertex(gemcOpt->optMap["LUMI_V"].args);
	lbeamtab->vtxControl->set_rvertex(gemcOpt->optMap["LUMI_SPREAD_V"].args);
	
	lbeamtab2->vtxControl->set_vertex(gemcOpt->optMap["LUMI2_V"].args);
	lbeamtab2->vtxControl->set_rvertex(gemcOpt->optMap["LUMI2_SPREAD_V"].args);
	
}

run_control::~run_control()
{
	string hd_msg = gemcOpt->optMap["LOG_MSG"].args ;
	
	if(verbosity>2)
		cout << hd_msg << " Run Control Widget Deleted." << endl;
}



void run_control::changePars()
{
	// Momentum
	gemcOpt->optMap["BEAM_P"].args         = pbeamtab->momControl->get_momentum(verbosity);
	gemcOpt->optMap["SPREAD_P"].args       = pbeamtab->momControl->get_rmomentum(verbosity);
	
	gemcOpt->optMap["LUMI_P"].args         = lbeamtab->momControl->get_momentum(verbosity);
	gemcOpt->optMap["LUMI_SPREAD_P"].args  = lbeamtab->momControl->get_rmomentum(verbosity);
	
	gemcOpt->optMap["LUMI2_P"].args        = lbeamtab2->momControl->get_momentum(verbosity);
	gemcOpt->optMap["LUMI2_SPREAD_P"].args = lbeamtab2->momControl->get_rmomentum(verbosity);
	

	// Vertex
	gemcOpt->optMap["BEAM_V"].args         = pbeamtab->vtxControl->get_vertex(verbosity);
	gemcOpt->optMap["SPREAD_V"].args       = pbeamtab->vtxControl->get_rvertex(verbosity);
	
	gemcOpt->optMap["LUMI_V"].args         = lbeamtab->vtxControl->get_vertex(verbosity);
	gemcOpt->optMap["LUMI_SPREAD_V"].args  = lbeamtab->vtxControl->get_rvertex(verbosity);
	
	gemcOpt->optMap["LUMI2_V"].args        = lbeamtab2->vtxControl->get_vertex(verbosity);
	gemcOpt->optMap["LUMI2_SPREAD_V"].args = lbeamtab2->vtxControl->get_rvertex(verbosity);
	

	// Luminosity Event
	gemcOpt->optMap["LUMI_EVENT"].args     = lbeamtab->get_lumi_event(0);
	gemcOpt->optMap["LUMI2_EVENT"].args    = lbeamtab2->get_lumi_event(lbeamtab->get_time_window());
	
}















