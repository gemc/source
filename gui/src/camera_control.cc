// Qt headers
#include <QtWidgets>

// windows.h and gl.h in windows
#ifdef _MSC_VER
	#include <windows.h>	
  #include <gl.h>
#endif


// gemc headers
#include "camera_control.h"
#include "string_utilities.h"

// C++ headers

camera_control::camera_control(QWidget *parent, goptions *Opts) : QWidget(parent)
{
	gemcOpt = Opts;
	UImanager  = G4UImanager::GetUIpointer();
	
	
	QGroupBox *anglesGroup = new QGroupBox(tr("Camera Control"));
	anglesGroup->setMinimumWidth(450);

	// move camera or detector
	QLabel *moveLabel = new QLabel(tr("Move:"));
	moveCombo = new QComboBox;
	moveCombo->addItem(tr("Light"));
	moveCombo->addItem(tr("Detector"));
	
	
	// projection: Orthogonal / Perspective
	QLabel *projLabel = new QLabel(tr("Projection:"));
	projCombo = new QComboBox;
	projCombo->addItem(tr("Orthogonal"));
	projCombo->addItem(tr("Perspective 30"));
	projCombo->addItem(tr("Perspective 45"));
	projCombo->addItem(tr("Perspective 60"));
	connect ( projCombo   , SIGNAL( currentIndexChanged (int) ), this, SLOT( set_perspective(int)    ) );
	
	QHBoxLayout *mpLayout = new QHBoxLayout;
	mpLayout->addWidget(moveLabel);
	mpLayout->addWidget(moveCombo);
	mpLayout->addSpacing(100);
	mpLayout->addWidget(projLabel);
	mpLayout->addWidget(projCombo);
	
	theta_hall = phi_hall = 0;  // initial angles

	// Theta Controls
	QLabel *thetaLabel = new QLabel(tr("theta"));
	theta_slider = new QSlider(Qt::Horizontal);
	theta_slider->setTickInterval(1);
	theta_slider->setRange(0, 180);
	
	QLCDNumber *Theta_LCD = new QLCDNumber(this);
	Theta_LCD->setFont(QFont("Helvetica", 32, QFont::Bold));
	Theta_LCD->setMaximumSize(45, 45);
	Theta_LCD->setSegmentStyle(QLCDNumber::Flat);
	for(int t=0; t<=180; t+=30)
	{
		char tmp[100];
		sprintf(tmp, "%d", t);
		ThetaSet.push_back(tmp);
	}
	
	QComboBox *ThetaCombo = new QComboBox;
	for(unsigned int i=0; i<ThetaSet.size(); i++)  ThetaCombo->addItem(tr(ThetaSet[i].c_str()));
	ThetaCombo->setMaximumSize(60, 45);
	QHBoxLayout *thetaLayout = new QHBoxLayout;
	thetaLayout->addWidget(thetaLabel);
	thetaLayout->addWidget(theta_slider);
	thetaLayout->addWidget(Theta_LCD);
	thetaLayout->addWidget(ThetaCombo);

	connect ( theta_slider , SIGNAL( valueChanged        (int) ),     this,   SLOT( change_theta(int)   ) );
	connect ( theta_slider , SIGNAL( valueChanged        (int) ),  Theta_LCD, SLOT( display(int)        ) );
	connect ( ThetaCombo   , SIGNAL( currentIndexChanged (int) ),     this,   SLOT( set_theta(int)      ) );
	connect ( ThetaCombo   , SIGNAL( currentIndexChanged (int) ),     this,   SLOT( change_theta_s(int) ) );
	
	// Phi Controls
	QLabel *phiLabel = new QLabel(tr("phi"));
	phi_slider = new QSlider(Qt::Horizontal);
	phi_slider->setTickInterval(1);
	phi_slider->setRange(0, 360);
	QLCDNumber *Phi_LCD = new QLCDNumber(this);
	Phi_LCD->setFont(QFont("Helvetica", 32, QFont::Bold));
	Phi_LCD->setMaximumSize(45, 45);
	Phi_LCD->setSegmentStyle(QLCDNumber::Flat);
	for(int t=0; t<=360; t+=30)
	{
		char tmp[100];
		sprintf(tmp, "%d", t);
		PhiSet.push_back(tmp);
	}
	QComboBox *PhiCombo = new QComboBox;
	for(unsigned int i=0; i<PhiSet.size(); i++)  PhiCombo->addItem(tr(PhiSet[i].c_str()));
	PhiCombo->setMaximumSize(60, 45);
	QHBoxLayout *phiLayout = new QHBoxLayout;
	phiLayout->addWidget(phiLabel);
	phiLayout->addWidget(phi_slider);
	phiLayout->addWidget(Phi_LCD);
	phiLayout->addWidget(PhiCombo);
	
	connect ( phi_slider , SIGNAL( valueChanged        (int) ),     this, SLOT( change_phi(int)   ) );
	connect ( phi_slider , SIGNAL( valueChanged        (int) ),  Phi_LCD, SLOT( display(int)      ) );
	connect ( PhiCombo   , SIGNAL( currentIndexChanged (int) ),     this, SLOT( set_phi(int)      ) );
	connect ( PhiCombo   , SIGNAL( currentIndexChanged (int) ),     this, SLOT( change_phi_s(int) ) );

	QVBoxLayout *anglesLayout = new QVBoxLayout;
	anglesLayout->addLayout(mpLayout);
	anglesLayout->addSpacing(12);
	anglesLayout->addLayout(thetaLayout);
	anglesLayout->addSpacing(2);
	anglesLayout->addLayout(phiLayout);
	anglesGroup->setLayout(anglesLayout);
	


	
	
	// x slice
	sliceXEdit = new QLineEdit(tr("0"));
	sliceXEdit->setMaximumWidth(100);

	sliceXActi = new QCheckBox(tr("&Active: "));
	sliceXActi->setChecked(false);

	sliceXInve = new QCheckBox(tr("&Invert: "));
	sliceXInve->setChecked(false);

	QHBoxLayout *sliceXLayout = new QHBoxLayout;
	sliceXLayout->addWidget(new QLabel(tr("X: ")));
	sliceXLayout->addWidget(sliceXEdit);
	sliceXLayout->addStretch(1);
	sliceXLayout->addWidget(sliceXActi);
 	sliceXLayout->addWidget(sliceXInve);
	sliceXLayout->addStretch(1);

 	
	// y slice
	sliceYEdit = new QLineEdit(tr("0"));
	sliceYEdit->setMaximumWidth(100);
	
	sliceYActi = new QCheckBox(tr("&Active: "));
	sliceYActi->setChecked(false);
	
	sliceYInve = new QCheckBox(tr("&Invert: "));
	sliceYInve->setChecked(false);
	
	QHBoxLayout *sliceYLayout = new QHBoxLayout;
	sliceYLayout->addWidget(new QLabel(tr("Y: ")));
	sliceYLayout->addWidget(sliceYEdit);
	sliceYLayout->addStretch(1);
	sliceYLayout->addWidget(sliceYActi);
 	sliceYLayout->addWidget(sliceYInve);
	sliceYLayout->addStretch(1);


	
	// z slice
	sliceZEdit = new QLineEdit(tr("0"));
	sliceZEdit->setMaximumWidth(100);
	
	sliceZActi = new QCheckBox(tr("&Active: "));
	sliceZActi->setChecked(false);

	sliceZInve = new QCheckBox(tr("&Invert: "));
	sliceZInve->setChecked(false);

	QHBoxLayout *sliceZLayout = new QHBoxLayout;
	sliceZLayout->addWidget(new QLabel(tr("Z: ")));
	sliceZLayout->addWidget(sliceZEdit);
	sliceZLayout->addStretch(1);
 	sliceZLayout->addWidget(sliceZActi);
 	sliceZLayout->addWidget(sliceZInve);
	sliceZLayout->addStretch(1);

	
	connect ( sliceXEdit , SIGNAL( textChanged(const QString &) ), this, SLOT( slice() ) );
	connect ( sliceYEdit , SIGNAL( textChanged(const QString &) ), this, SLOT( slice() ) );
	connect ( sliceZEdit , SIGNAL( textChanged(const QString &) ), this, SLOT( slice() ) );

	connect ( sliceXActi , SIGNAL( stateChanged(int) ), this, SLOT( slice() ) );
	connect ( sliceYActi , SIGNAL( stateChanged(int) ), this, SLOT( slice() ) );
	connect ( sliceZActi , SIGNAL( stateChanged(int) ), this, SLOT( slice() ) );
	

	QPushButton *clearSliceButton = new QPushButton(tr("Clear Slices"));
	clearSliceButton->setToolTip("Clear Slice Planes");
	clearSliceButton->setIcon(style()->standardIcon(QStyle::SP_DialogResetButton));
	connect ( clearSliceButton , SIGNAL(clicked()), this, SLOT( clearSlice() ) );

	
	// slices layout
	QVBoxLayout *sliceLayout = new QVBoxLayout;
	sliceLayout->addLayout(sliceXLayout);
	sliceLayout->addLayout(sliceYLayout);
	sliceLayout->addLayout(sliceZLayout);
	sliceLayout->addWidget(clearSliceButton);
	
	// slices group
	QGroupBox *sliceGroup = new QGroupBox(tr("Slices   [mm]"));
	sliceGroup->setLayout(sliceLayout);


	
	
	QLabel *antialiasingLabel = new QLabel(tr("Anti-Aliasing"));
	aliasing = new QComboBox;
	aliasing->addItem(tr("OFF"));
	aliasing->addItem(tr("ON"));
	QHBoxLayout *aliasingLayout = new QHBoxLayout;
	aliasingLayout->addWidget(antialiasingLabel);
	aliasingLayout->addWidget(aliasing);
	connect ( aliasing   , SIGNAL( currentIndexChanged (int) ), this, SLOT( switch_antialiasing(int)    ) );
	
	QLabel *sides_per_circlesLabel = new QLabel(tr("Sides per circle"));
	sides_per_circle = new QComboBox;
	sides_per_circle->addItem(tr("50"));
	sides_per_circle->addItem(tr("100"));
	sides_per_circle->addItem(tr("200"));
	sides_per_circle->addItem(tr("500"));
	sides_per_circle->setCurrentIndex(1);
	QHBoxLayout *sides_per_circleLayout = new QHBoxLayout;
	sides_per_circleLayout->addWidget(sides_per_circlesLabel);
	sides_per_circleLayout->addWidget(sides_per_circle);
	connect ( sides_per_circle   , SIGNAL( currentIndexChanged (int) ), this, SLOT( switch_sides_per_circle(int)    ) );
	
	QLabel *auxiliaryEdgesLabel = new QLabel(tr("Auxiliary Edges"));
	auxiliary = new QComboBox;
	auxiliary->addItem(tr("OFF"));
	auxiliary->addItem(tr("ON"));
	QHBoxLayout *auxedgesLayout = new QHBoxLayout;
	auxedgesLayout->addWidget(auxiliaryEdgesLabel);
	auxedgesLayout->addWidget(auxiliary);
	connect ( auxiliary   , SIGNAL( currentIndexChanged (int) ), this, SLOT( switch_auxiliary_edges(int)    ) );
	
	QVBoxLayout *voptionsLayout = new QVBoxLayout;
	voptionsLayout->addLayout(aliasingLayout);
	voptionsLayout->addLayout(sides_per_circleLayout);
	voptionsLayout->addLayout(auxedgesLayout);
	
	// options group
	QGroupBox *vOptionsGroup = new QGroupBox(tr("Visualization Options"));
	vOptionsGroup->setLayout(voptionsLayout);
	
	
	QHBoxLayout *sliceOptionsLayout = new QHBoxLayout;
	sliceOptionsLayout->addWidget(sliceGroup);
	sliceOptionsLayout->addWidget(vOptionsGroup);
	
	
	// Explode Controls
	// slider is from 1 to 2.5 in intervals of 0.05 (30 total)
	explodeSlider = new QSlider(Qt::Horizontal);
	explodeSlider->setTickInterval(1);
	explodeSlider->setRange(0, 30);
	connect ( explodeSlider , SIGNAL( valueChanged(int) ), this, SLOT( explode(int) ) );
	
	QHBoxLayout *explodeLayout = new QHBoxLayout;
	explodeLayout->addWidget(explodeSlider);

	QGroupBox *explodeGroup = new QGroupBox(tr("Explode"));
	explodeGroup->setLayout(explodeLayout);


	
	// All together
	QVBoxLayout *mainLayout = new QVBoxLayout;
	mainLayout->addWidget(anglesGroup);
	mainLayout->addLayout(sliceOptionsLayout);
	mainLayout->addSpacing(6);
	mainLayout->addWidget(explodeGroup);
	mainLayout->addStretch(1);
	setLayout(mainLayout);
	
	
	
}

void camera_control::slice()
{
	UImanager->ApplyCommand("/vis/viewer/clearCutawayPlanes");
	UImanager->ApplyCommand("/vis/viewer/clearCutawayPlanes");
	UImanager->ApplyCommand("/vis/viewer/clearCutawayPlanes");
	UImanager->ApplyCommand("/vis/viewer/clearCutawayPlanes");

	UImanager->ApplyCommand("/vis/viewer/set/cutawayMode intersection");
	
	char command[200] = "";
		
	if(sliceXActi->isChecked() )
	{
		sprintf(command, "/vis/viewer/addCutawayPlane  %d   0  0 mm %d 0 0 ", (int) qs_toDouble(sliceXEdit->text()), sliceXInve->isChecked() ? -1 : 1 );
		UImanager->ApplyCommand(command);
	}
	if(sliceYActi->isChecked() )
	{
		sprintf(command, "/vis/viewer/addCutawayPlane   0  %d  0 mm 0 %d 0 ", (int) qs_toDouble(sliceYEdit->text()), sliceYInve->isChecked() ? -1 : 1);
		UImanager->ApplyCommand(command);
	}
	if(sliceZActi->isChecked() )
	{
		sprintf(command, "/vis/viewer/addCutawayPlane   0   0 %d mm 0 0 %d ", (int) qs_toDouble(sliceZEdit->text()), sliceZInve->isChecked() ? -1 : 1);
		UImanager->ApplyCommand(command);
	}



}




void camera_control::clearSlice()
{
	// making sure it catches all planes, I thought
	// one command would be enough for all
	UImanager->ApplyCommand("/vis/viewer/clearCutawayPlanes");
	UImanager->ApplyCommand("/vis/viewer/clearCutawayPlanes");
	UImanager->ApplyCommand("/vis/viewer/clearCutawayPlanes");
	UImanager->ApplyCommand("/vis/viewer/clearCutawayPlanes");
	
	
	
	sliceXActi->setChecked(false);
	sliceYActi->setChecked(false);
	sliceZActi->setChecked(false);

	
}


void camera_control::change_theta(int theta)
{
	char command[100];
	theta_hall = theta - 1;

	if(qs_tostring(moveCombo->currentText()) == "Detector")
		sprintf(command,"/vis/viewer/set/viewpointThetaPhi %d %d.", theta_hall, phi_hall);

	if(qs_tostring(moveCombo->currentText()) == "Light")
		sprintf(command,"/vis/viewer/set/lightsThetaPhi %d %d.", theta_hall, phi_hall);

	UImanager->ApplyCommand(command);
}

void camera_control::set_theta(int index)
{
	char command[100];
	theta_hall = atoi(ThetaSet[index].c_str());
	sprintf(command,"/vis/viewer/set/viewpointThetaPhi %d %d.", theta_hall, phi_hall);
	UImanager->ApplyCommand(command);
}

void camera_control::change_theta_s(int theta)
{
	theta_slider->setSliderPosition(atoi(ThetaSet[theta].c_str()));
}

void camera_control::change_phi(int phi)
{
	char command[100];
	phi_hall = phi - 1;


	if(qs_tostring(moveCombo->currentText()) == "Detector")
		sprintf(command,"/vis/viewer/set/viewpointThetaPhi %d %d", theta_hall, phi_hall);

	if(qs_tostring(moveCombo->currentText()) == "Light")
		sprintf(command,"/vis/viewer/set/lightsThetaPhi %d %d.", theta_hall, phi_hall);
		
	UImanager->ApplyCommand(command);
}

void camera_control::change_phi_s(int phi)
{
	phi_slider->setSliderPosition(atoi(PhiSet[phi].c_str()));
}

void camera_control::set_phi(int index)
{
	char command[100];
	phi_hall = atoi(PhiSet[index].c_str());
	sprintf(command,"/vis/viewer/set/viewpointThetaPhi %d %d.", theta_hall, phi_hall);
	UImanager->ApplyCommand(command);
}


void camera_control::explode(int boom)
{
	char command[100];

	double xf = 1 +  ((double) boom)/15.0;
	
	sprintf(command,"/vis/viewer/set/explodeFactor %3.2f", xf);
	
	UImanager->ApplyCommand(command);
}

void camera_control::set_perspective(int index)
{
	char command[100];
	double angles[4] = { 0,   30,  45, 60};
	string which[4]  = {"o", "p", "p", "p"};

	sprintf(command,"/vis/viewer/set/projection %s %4.2f ",which[index].c_str(),  angles[index]);
	UImanager->ApplyCommand(command);
}

void camera_control::switch_antialiasing(int index)
{
	if(index == 0)
	{
		glDisable (GL_LINE_SMOOTH);
		glDisable (GL_POLYGON_SMOOTH);
	}
	else
	{
		glEnable (GL_LINE_SMOOTH);
		glHint (GL_LINE_SMOOTH_HINT, GL_NICEST);
		glEnable (GL_POLYGON_SMOOTH);
		glHint (GL_POLYGON_SMOOTH_HINT, GL_NICEST);
	}
	UImanager->ApplyCommand("/vis/viewer/flush");
}

void camera_control::switch_auxiliary_edges(int index)
{
	if(index == 0)
	{
		UImanager->ApplyCommand("/vis/viewer/set/auxiliaryEdge 0");
		UImanager->ApplyCommand("/vis/viewer/set/hiddenEdge 0");
	}
	else
	{
		UImanager->ApplyCommand("/vis/viewer/set/auxiliaryEdge 1");
		UImanager->ApplyCommand("/vis/viewer/set/hiddenEdge 1");
	}
	UImanager->ApplyCommand("/vis/viewer/flush");
}

void camera_control::switch_sides_per_circle(int index)
{
	char command[100];
	int sides[4] = { 50,   100,  200, 500};
	
	sprintf(command,"/vis/viewer/set/lineSegmentsPerCircle %d ", sides[index]);
	UImanager->ApplyCommand(command);
	UImanager->ApplyCommand("/vis/viewer/flush");
}


void camera_control::update_angles()
{
// 	  G4VViewer* currentViewer = visManager->GetCurrentViewer();
// 	  theta_hall = (int) floor(currentViewer->GetViewParameters().GetViewpointDirection().getTheta()/deg);
// 	  phi_hall   = (int) floor(currentViewer->GetViewParameters().GetViewpointDirection().getPhi()/deg );
// 	  TransfSli[0]->setSliderPosition ( theta_hall + 1);
// 	  TransfSli[1]->setSliderPosition ( phi_hall + 1);
}
	
	
camera_control::~camera_control()
{
	string hd_msg = gemcOpt->optMap["LOG_MSG"].args ;
	double VERB   = gemcOpt->optMap["GEO_VERBOSITY"].arg ;
	if(VERB>2)
		cout << hd_msg << " Camera Control Widget Deleted." << endl;
	
}


