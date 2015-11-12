#ifndef gsignal_H
#define gsignal_H 1

// Qt headers
#include <QtWidgets>

// gemc headers
#include "options.h"
#include "utils.h"
#include "sensitiveDetector.h"
#include "graph.h"

// G4 headers
//#include "G4UImanager.hh"

// C++ headers
#include <string>
#include <map>
using namespace std;



// Class definition
class gsignal : public QWidget
{
	// metaobject required for non-qt slots
	Q_OBJECT
		
	public:
		gsignal(QWidget *parent, goptions*, map<string, sensitiveDetector*>);
		~gsignal();
	
		// we need to double check that
		// some variables are available to display
		string WRITE_INTRAW;
		int SAVE_ALL_MOTHERS;
	
		map<string, sensitiveDetector*> SeDe_Map;
		goptions *gemcOpt;

		string signalChoice;    // what to plot
		vector<string> availableSignals;
		QComboBox *qcsignal;    // combo box
	
	
	private:
		QTreeWidget *hitList;
		QTreeWidget *hitData;
		QTreeWidget *var_choice;

		QLinearGradient HitGrad;
		QBrush          HitBrush;
		
		graph *graphView;
		
		
	public slots:
		
		void chooseVariable(int);
		void createHitListTree();         // creates Sensitive Detector / Hits Tree
		void createSignalsTree();         // creates signals tree
		
};

#endif








