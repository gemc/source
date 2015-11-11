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

	
	private:
	
		
	public slots:
		
	
};

#endif








