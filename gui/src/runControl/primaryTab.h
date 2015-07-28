#ifndef particleTab_H
#define particleTab_H 1

// Qt headers
#include <QtWidgets>

// gemc headers
#include "options.h"
#include "momControls.h"
#include "vtxControls.h"




// Luminosity Event class with particle type and
// time window info
class lumiEvent : public QWidget
{
	Q_OBJECT
		
	public:
		lumiEvent(goptions *Opts, QWidget*, string);
	
		QLineEdit *nevents;
		QLineEdit *timewindow;
		QLineEdit *time_bunch;
		
		QGroupBox *LumiGroup;
};




// Class definition
class particleTab : public QWidget
{
	Q_OBJECT
	
	public:
		particleTab(QWidget*, goptions* , string);
		
		// self identification: Primary, Lumi1 or Lumi2
		string type;
		int verbosity;
	
		QComboBox *beam_particle;
		
		lumiEvent   *lEvent;
	
		momControls *momControl;
		vtxControls *vtxControl;	
	
		string get_lumi_event(double);
		double get_time_window();
};



#endif






