#ifndef gsignal_H
#define gsignal_H 1

// Qt headers
#include <QtWidgets>

// gemc headers
#include "options.h"
#include "utils.h"
#include "sensitiveDetector.h"

// G4 headers
#include "G4UImanager.hh"

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
		
		int xorig, yorig;       // origin of the axis
		int xaxil, yaxil;       // axis length
		double xmin, ymin;      // graph minima
		double xmax, ymax;      // graph maxima
		double inside;          // how much inside the graph will be
		double dx, dy, DX, DY;  // lowercase: graph limits; uppercase: scene inside limits
		int nticksx, nticksy;
		
		string signalChoice;    // what to plot
		vector<string> availableSignals;
		QComboBox *qcsignal;    // combo box
	
		void plots_bg(string xtit, string ytit, vector<double> x, vector<double> y, string title);  // draw axis, ticks and labels
		void osci_bg( string xtit, string ytit, vector<double> x, vector<double> y, string title);  // draw axis, ticks, labels and oscilloscope background
		void plot_graph(vector<double> x, vector<double> y, vector<int> pid);
		
		goptions *gemcOpt;
		// we need to double check that
		// some variables are available to display
		string WRITE_INTRAW;
		int SAVE_ALL_MOTHERS;
	
		G4UImanager  *UImanager;
		map<string, sensitiveDetector*> SeDe_Map;
		
		map<int, QPen> pcolors;
		
	private:
		QTreeWidget *gsignals;
		QTreeWidget *s_detectors;
		QTreeWidget *var_choice;

		QLinearGradient HitGrad;
		QBrush          HitBrush;
		
		QGraphicsView *view;
		QGraphicsScene *scene;
		
		
	public slots:
		
		void chooseVariable(int);
		QTreeWidget* CreateSDetsTree();           // creates Sensitive Detector / Hits Tree
		QTreeWidget* CreateSignalsTree();         // creates signals tree
		
};

#endif








