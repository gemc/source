// Qt headers
#include <QApplication>
#include <QLineEdit>

// gemc headers
#include "gemc_MainGui.h"
#include "string_utilities.h"

// All icon images must be 256x256
#include "images/aboutControl_xpm.h"
#include "images/cameraControl_xpm.h"
#include "images/detectorControl_xpm.h"
#include "images/dialogControl_xpm.h"
#include "images/generatorControl_xpm.h"
#include "images/physicsControl_xpm.h"
#include "images/signalsControl_xpm.h"
#include "images/trigger_xpm.h"

// C++ headers
#include <iostream>
#include <string>
using namespace std;

gemcMainWidget::gemcMainWidget(goptions *Opts, G4RunManager *rm,  map<string, sensitiveDetector*> SDM, map<string, detector> *hMap, map<string, G4Material*> matMap)
{
	gemcOpt    = Opts;
	RM = rm;

	gtimer = new QTimer(this);
    connect(gtimer, SIGNAL(timeout()), this, SLOT(beamOnCycle()));
	
	// revisit this, why we need it?
	double qt_mode = gemcOpt->optMap["USE_GUI"].arg ;
	
	contentsWidget = new QListWidget;
	contentsWidget->setViewMode(QListView::IconMode);
	
	// icon size
	contentsWidget->setIconSize(QSize(55, 55));
	contentsWidget->setMovement(QListView::Static);

	// icon container sizes
	contentsWidget->setMinimumWidth(76);
	contentsWidget->setMaximumWidth(76);
	contentsWidget->setMinimumHeight(600);
	contentsWidget->setMaximumHeight(600);
	
	// makes all icon the same size
	contentsWidget->setUniformItemSizes(1);
	
	// activate tracking (for mouse over effect)
	contentsWidget->setMouseTracking(1);
	contentsWidget->setSpacing(2);

	// Content:
	// Generator Control
	// Camera Control
	// Infos Control
	// G4 Dialog Control
	// Signal Control
	// Physics List Control
  
	pagesWidget = new QStackedWidget;
	pagesWidget->addWidget(new run_control     (this, gemcOpt));  // for some reason run_control is very slow
	pagesWidget->addWidget(new camera_control  (this, gemcOpt));
	pagesWidget->addWidget(new detector_tree   (this, *gemcOpt, RM, hMap,  matMap));
	pagesWidget->addWidget(new infos           (this, gemcOpt));
	pagesWidget->addWidget(new g4dialog        (this, gemcOpt));
	pagesWidget->addWidget(gsig  = new gsignal (this, gemcOpt, SDM));
	pagesWidget->addWidget(gtrig = new gtrigger(this, gemcOpt, SDM));
	pagesWidget->setMinimumWidth(550);
	pagesWidget->setMaximumWidth(550);
	pagesWidget->setMinimumHeight(600);
	pagesWidget->setMaximumHeight(600);
	
	createIcons();
	contentsWidget->setCurrentRow(0);
	
	// revisit this, why we need it?
	if(qt_mode > 1)
		contentsWidget->setCurrentRow(1);
	
	
	QPushButton *runButton = new QPushButton(tr("Run"));
	runButton->setToolTip("Run 1 event");
	runButton->setIcon(style()->standardIcon(QStyle::SP_MediaPlay));
	
	QPushButton *cycleButton = new QPushButton(tr("Cycle"));
	cycleButton->setToolTip("Run 1 event every 2 seconds");
	cycleButton->setIcon(style()->standardIcon(QStyle::SP_BrowserReload));
	
	QPushButton *stopButton = new QPushButton(tr("Stop"));
	stopButton->setToolTip("Stops running events");
	stopButton->setIcon(style()->standardIcon(QStyle::SP_MediaStop));

	QPushButton *closeButton = new QPushButton(tr("Exit"));
	closeButton->setToolTip("Quit GEMC");
	closeButton->setIcon(style()->standardIcon(QStyle::SP_TitleBarCloseButton));
	
	
	QLabel *nEventsLabel = new QLabel(tr("N. Events:"));
	nEvents = new QLineEdit(tr("1"));
	nEvents->setMaximumWidth(50);

	
	connect(closeButton, SIGNAL(clicked()), this, SLOT(gemc_quit()));
	connect(runButton,   SIGNAL(clicked()), this, SLOT(beamOn()));
	connect(cycleButton, SIGNAL(clicked()), this, SLOT(beamOnCycle()));
	connect(stopButton,  SIGNAL(clicked()), this, SLOT(stopBeam()));
	
	
	QHBoxLayout *horizontalLayout = new QHBoxLayout;
	horizontalLayout->addWidget(contentsWidget);
	horizontalLayout->addWidget(pagesWidget, 1);
	
	QHBoxLayout *buttonsLayout = new QHBoxLayout;
	buttonsLayout->addWidget(nEventsLabel);
	buttonsLayout->addWidget(nEvents);
	buttonsLayout->addWidget(runButton);
	buttonsLayout->addWidget(cycleButton);
	buttonsLayout->addWidget(stopButton);
	buttonsLayout->addStretch(1);
	buttonsLayout->addWidget(closeButton);
	
	
	QVBoxLayout *mainLayout = new QVBoxLayout;
	mainLayout->addLayout(buttonsLayout);
	mainLayout->addLayout(horizontalLayout);
	setLayout(mainLayout);
	
	setWindowTitle(tr("Config Dialog"));

	return;
}



void gemcMainWidget::changePage(QListWidgetItem *current, QListWidgetItem *previous)
{
	if (!current)
		current = previous;
	
	int thisIndex = contentsWidget->row(current);
	pagesWidget->setCurrentIndex(thisIndex);
	
	if(thisIndex == 5)
	{
		gsig->createHitListTree();
	}
	if(thisIndex == 6)
	{
		gtrig->createGraphs();
	}

}

void gemcMainWidget::createIcons()
{
	i_default = QColor( 255,  255, 255);
	i_hovered = QColor( 200,  200, 240);
	
	buttons.push_back(addItem("Generator", QPixmap(generatorControl_xpm)));
	buttons.push_back(addItem("Camera",    QPixmap(cameraControl_xpm)));
	buttons.push_back(addItem("Detector",  QPixmap(detectorControl_xpm)));
	buttons.push_back(addItem("Infos",     QPixmap(aboutControl_xpm)));
	buttons.push_back(addItem("G4Dialog",  QPixmap(dialogControl_xpm)));
	buttons.push_back(addItem("Signals",   QPixmap(signalsControl_xpm)));
	buttons.push_back(addItem("Trigger",   QPixmap(trigger_xpm)));
	buttons.push_back(addItem("Physics",   QPixmap(physicsControl_xpm)));

	connect(contentsWidget,
					SIGNAL(currentItemChanged(QListWidgetItem *, QListWidgetItem*)),
          this, SLOT(changePage(QListWidgetItem *, QListWidgetItem*)));
	
	connect(contentsWidget,	SIGNAL(itemEntered(QListWidgetItem *)), this, SLOT(change_background(QListWidgetItem *)) );
	
}

void gemcMainWidget::change_background(QListWidgetItem* item)
{
	for(int i=0; i<contentsWidget->count(); i++)
		contentsWidget->item(i)->setBackground(QBrush(i_default));
	
	item->setBackground(QBrush(i_hovered));

}



QListWidgetItem* gemcMainWidget::addItem(string text, QPixmap pixmap)
{
	QListWidgetItem *newItem = new QListWidgetItem(contentsWidget);
	newItem->setFont(QFont("Helvetica", 9, QFont::Normal));
	newItem->setIcon(pixmap);
	newItem->setText(tr(text.c_str()));
	newItem->setBackground(QBrush(i_default));
	newItem->setTextAlignment(Qt::AlignHCenter);
	newItem->setFlags( Qt::ItemIsSelectable | Qt::ItemIsEnabled);
	
	return newItem;
}

void gemcMainWidget::beamOn()
{
	G4UImanager *uim = G4UImanager::GetUIpointer();
	
	int nevs = atoi(qs_tostring(nEvents->text()).c_str());

	char command[100];
	sprintf(command, "/run/beamOn %d", nevs);
	uim->ApplyCommand(command);

	if(pagesWidget->currentIndex() == 5)
	{
		gsig->createHitListTree();
	}
	if(pagesWidget->currentIndex() == 6)
	{
		gtrig->createGraphs();
	}

}

void gemcMainWidget::beamOnCycle()
{
	
	playing = 1;
	gtimer->start(2000);
	
	G4UImanager *uim = G4UImanager::GetUIpointer();
	
	char command[100];
	sprintf(command, "/vis/viewer/refresh");
	uim->ApplyCommand(command);
	
	sprintf(command, "/run/beamOn 1");
	uim->ApplyCommand(command);
	
}

void gemcMainWidget::stopBeam()
{
	playing = 0;
	gtimer->stop();
}


gemcMainWidget::~gemcMainWidget()
{
	string hd_msg = gemcOpt->optMap["LOG_MSG"].args;


	cout << endl;
}





