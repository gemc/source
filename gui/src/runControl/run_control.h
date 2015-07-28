#ifndef run_control_H
#define run_control_H 1

// Qt headers
#include <QWidget>
#include <QSlider>
#include <QLabel>
#include <QComboBox>
#include <QtWidgets>

// gemc headers
#include "options.h"
#include "primaryTab.h"

// G4 headers
#include "G4UImanager.hh"

// C++ headers
#include <string>
#include <map>
using namespace std;


class run_control : public QWidget
{
	// metaobject required for non-qt slots
	Q_OBJECT
	
	public:
		run_control(QWidget *parent, goptions*);
	   ~run_control();
		
		goptions *gemcOpt;

	public slots:
		void changePars();

	private:
		particleTab *pbeamtab;
		particleTab *lbeamtab;
		particleTab *lbeamtab2;
			
		double verbosity;
	
};

#endif
