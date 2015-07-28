// Qt headers
#include <QtWidgets>

// gemc headers
#include "physicsListGui.h"

// G4 headers
#include "G4UIcommandTree.hh"

physicsList::physicsList(QWidget *parent, goptions *Opts) : QWidget(parent)
{
	gemcOpt = Opts;
	UImanager  = G4UImanager::GetUIpointer();
	
	//  Layout:
	//
	//  + +----------------------------+ +
	//  | |            |               | |
	//  | |  Phys List |  Description  | |
	//  | |            |               | |
	//  + +----------------------------+ +
 	
	
	// Vertical Splitter - Top and Bottom layouts
	QSplitter *splitter     = new QSplitter(Qt::Vertical);
	
	// top layout
	QWidget   *topWidget    = new QWidget(splitter);
	QSplitter *treesplitter = new QSplitter(Qt::Horizontal);
	
	// treesplitter size
	QList<int> tlist;
	tlist.append( 400 );
	tlist.append( 400 );
	treesplitter->setSizes(tlist);
	
	QVBoxLayout *layoutTop = new QVBoxLayout(topWidget);
	layoutTop->addWidget(treesplitter);
	
	
	// all layouts
	QVBoxLayout *mainLayout = new QVBoxLayout;
	mainLayout->addWidget(splitter);
	setLayout(mainLayout);
	
}




physicsList::~physicsList()
{
	string hd_msg = gemcOpt->optMap["LOG_MSG"].args ;
	double VERB   = gemcOpt->optMap["PHYS_VERBOSITY"].arg ;
	if(VERB>2)
		cout << hd_msg << " g4 Physics List Widget Deleted." << endl;	
}









