#ifndef vtxControls_H
#define vtxControls_H 1

// Qt headers
#include <QtWidgets>

// gemc headers
#include "options.h"


class vtxControls : public QWidget
{
	Q_OBJECT
	
	public:
		
		vtxControls(goptions *Opts, string type);

		string get_vertex(double verbosity);
		string get_rvertex(double verbosity);

		void set_vertex(string vtxOption);
		void set_rvertex(string vtxOption);
	
	

	public:
	
		// user input for vertex
		QLineEdit *vxLE;
		QLineEdit *vyLE;
		QLineEdit *vzLE;
		
		QLineEdit *rvrLE;
		QLineEdit *rvzLE;

	
		// user input for units
		QComboBox *vtxUnits;
	
		// accessed in primary Tab
		QGroupBox *vertexGroup;
	
};

#endif






