#ifndef gemc_MAINGUI_H
#define gemc_MAINGUI_H

// Qt headers
#include <QtWidgets>

// G4 headers
#include "G4RunManager.hh"

// gemc headers
#include "options.h"
#include "runControl/run_control.h"
#include "camera_control.h"
#include "detector_tree.h"
#include "infos.h"
#include "g4dialog.h"
#include "gsignal.h"
#include "gtrigger.h"


// Class definition

class gemcMainWidget : public QWidget
{
  // metaobject required for non-qt slots
	Q_OBJECT
	
	public:
		gemcMainWidget(goptions*, map<string, sensitiveDetector*>, map<string, detector>*, map<string, G4Material*>);
		~gemcMainWidget();
		  
		goptions *gemcOpt;

	public slots:
		void changePage(QListWidgetItem *current, QListWidgetItem *previous);
		void change_background(QListWidgetItem*);

	private:
    	gsignal   *gsig;
		gtrigger  *gtrig;
		QLineEdit *nEvents;
		QTimer *gtimer;
		int playing;        // controls cycling event
	
	private slots:
		
		// we need to delete runManager before quitting the qApp
		// looks like no need to delete the visManager
		void  gemc_quit(){
			delete G4RunManager::GetRunManager();
			cout << endl << " Arrivederci! " << endl << endl;
			qApp->quit();
		}
		void stopBeam();
		
		void beamOn();
		void beamOnCycle();
	
	private:
		void createIcons();
		
		QListWidget *contentsWidget;
		QStackedWidget *pagesWidget;
	
		vector<QListWidgetItem*> buttons;
		QListWidgetItem* addItem(string text, QPixmap pixmap);
	
		QColor i_default;
		QColor i_hovered;
};

#endif



