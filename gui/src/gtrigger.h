#ifndef gtrigger_H
#define gtrigger_H 1

// Qt headers
#include <QtWidgets>

// gemc headers
#include "options.h"
#include "sensitiveDetector.h"
#include "graph.h"


// Class definition
class gtrigger : public QWidget
{
	// metaobject required for non-qt slots
	Q_OBJECT
		
	public:
		gtrigger(QWidget *parent, goptions*, map<string, sensitiveDetector*>);
		~gtrigger();
		
		map<string, sensitiveDetector*> SeDe_Map;
		goptions *gemcOpt;

		void createGraphs();

	private:
		QScrollArea *scrollArea;
		QSplitter *vGraphsplitter;
		string VOLTAGES;
		int plotChoice;

		void createSummary();
		int determineNhits();
	
	
	public slots:
		void choosePlots(int);


};

#endif








