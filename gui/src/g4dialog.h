#ifndef g4dialog_H
#define g4dialog_H 1

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

class g4dialog : public QWidget
{
	// metaobject required for non-qt slots
	Q_OBJECT
	
	public:
		g4dialog(QWidget *parent, goptions*);
	   ~g4dialog();
		
		goptions *gemcOpt;
		G4UImanager  *UImanager;

	private:  
    	QListWidget *fCommandHistoryArea;
    	QLineEdit   *fCommandArea;
    	QTextEdit   *fHelpArea;
    	QTreeWidget *fHelpTreeWidget;

	private slots:
    	void CommandHistoryCallback();        // calls back a command from the history box
    	void CommandEnteredCallback();        // command entered manually
  
    	QTreeWidget* CreateHelpTree();                                                  // creates help Tree
    	void CreateChildTree(QTreeWidgetItem *aParent,G4UIcommandTree *aCommandTree);   // creates branch
    	void HelpTreeClicCallback();                                                    // displays help on a command
    	QString GetCommandList (const G4UIcommand *aCommand);                           // gets the G4 help on a command
    	void HelpTreeDoubleClicCallback();                                              // sets a double clicked command to command line


};

#endif








