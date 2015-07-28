#ifndef momControls_H
#define momControls_H 1

// Qt headers
#include <QtWidgets>

// gemc headers
#include "options.h"


class momControls : public QWidget
{
	Q_OBJECT
	
	public:
		
		momControls(goptions *Opts, string type);
		
		string get_momentum(double verbosity);
		string get_rmomentum(double verbosity);

		void set_momentum(string momOption);
		void set_rmomentum(string momOption);
	
	
	public:
		QComboBox *beam_particle;
		
		// user input for momentum
		QLineEdit *momLE, *rmomLE;
		QLineEdit *theLE, *rtheLE;
		QLineEdit *phiLE, *rphiLE;

		// user input for units
		QComboBox *momUnits, *theUnits, *phiUnits;
	
		// accessed in primary Tab
		QGroupBox *momentumGroup;
	
};


#endif






