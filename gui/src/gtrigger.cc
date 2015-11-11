// Qt headers
#include <QtWidgets>

// gemc headers
#include "gtrigger.h"

gtrigger::gtrigger(QWidget *parent, goptions *Opts, map<string, sensitiveDetector*> SD_Map) : QWidget(parent)
{
	gemcOpt  = Opts;
	SeDe_Map = SD_Map;

	
//	//  Layout: just a big scroll area with all the signals in
//	QScrollArea *scrollArea = new QScrollArea;
//	scrollArea->setBackgroundRole(QPalette::Dark);

	
	// Vertical Splitter - Top and Bottom layouts
	QSplitter *vMainsplitter     = new QSplitter(this);
	
	
	QComboBox *whatToShow = new QComboBox;
	whatToShow->addItem("Voltages and Triggers");
	whatToShow->addItem("Voltages only");
	whatToShow->addItem("Triggers only");

	
	vMainsplitter->addWidget(whatToShow);
	
}


gtrigger::~gtrigger()
{
	string hd_msg = gemcOpt->optMap["LOG_MSG"].args ;
	double VERB   = gemcOpt->optMap["GEO_VERBOSITY"].arg ;
	if(VERB>2)
		cout << hd_msg << " Signal Widget Deleted." << endl;
	
}



