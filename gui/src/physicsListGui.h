#ifndef physics_list_H
#define physics_list_H 1

// Qt headers
#include <QWidget>
#include <QSlider>
#include <QComboBox>
#include <QTextEdit>
#include <QtWidgets>

// gemc headers
#include "options.h"

// G4 headers
#include "G4UImanager.hh"

// C++ headers
#include <string>
#include <map>
using namespace std;



// Class definition

class physicsList : public QWidget
{
	// metaobject required for non-qt slots
	Q_OBJECT
	
	public:
		physicsList(QWidget *parent, goptions*);
	   ~physicsList();
		
		goptions *gemcOpt;
		G4UImanager  *UImanager;

	private:  

	private slots:
  


};

#endif








