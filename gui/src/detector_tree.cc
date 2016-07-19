// Qt headers
#include <QActionGroup>
#include <QColorDialog>
#include <QMenu>

// gemc headers
#include "detector_tree.h"
#include "MDetectorConstruction.h"
#include "MPrimaryGeneratorAction.h"
#include "options.h"
#include "string_utilities.h"
#include "detector_editor.h"

// C++ headers
#include <sys/stat.h>
#include <vector>
#include <string>

detector_tree::detector_tree(QWidget *parent, goptions Opts, G4RunManager* RM, map<string, detector> *Map,
                             map<string, G4Material*> MMap) : QWidget(parent)
{
	gemcOpt    = Opts;
	runManager = RM;
	UImanager  = G4UImanager::GetUIpointer();
	Hall_Map   = Map;
	MMats      = MMap;
	createActions();
	
	//  Layout:
	//
	//  + +-------------------+ +
	//  | |         |         | |
	//  | |  Tree   |  infos  | |
	//  | |         |         | |
	//  | +-------------------+ |
	//  | |         |         | |
	//  | |         |         | |
	//  | |         |         | |
	//  | |         |         | |
	//  | +-------------------+ |
	//  | |         |         | |
	//  | |         |         | |
	//  | |         |         | |
	//  | +-------------------+ |
	
	
	
	
	// Detector Tree
	treeWidget = new QTreeWidget();
	treeWidget->setColumnCount(1);
	QStringList headers;
	headers << tr("Volumes")   ;
	treeWidget->setHeaderLabels(headers);
	
	// in case the geometry changes need to reload the tree
	// in case it doesn't, no need to load the tree
	tree_map = read_geometry(treeWidget);
	treeWidget->setContextMenuPolicy( Qt::ActionsContextMenu );
	treeWidget-> addAction( Switch_visibility );
	treeWidget-> addAction( Switch_visibility_daughters );
	treeWidget-> addAction( Switch_wiresolid );
	treeWidget-> addAction( Switch_color );
	treeWidget-> addAction( Write_GDML_File );
	treeWidget-> addAction( Write_WRL_File );
	connect(treeWidget, SIGNAL(itemSelectionChanged ()), this, SLOT(show_detector()));

	// detector infos
	dTab = new descriptionTab();

	
	QSplitter *splitter     = new QSplitter(Qt::Horizontal);

	QWidget *leftWidget = new QWidget(splitter);
	leftWidget->setMinimumWidth(300);
	QVBoxLayout *leftLayout = new QVBoxLayout(leftWidget);
	leftLayout->addWidget(treeWidget);
	

	// writing single volumes in GDML or WRL
	QHBoxLayout *writeSingleVolumeButtons = new QHBoxLayout();

	writeToGDML = new QPushButton(tr("detector to gdml"));
	writeToGDML->setEnabled(false);
	writeToGDML->setToolTip("Creates a gdml file for the highlighted detector");
	writeToGDML->setIcon(style()->standardIcon(QStyle::SP_DialogSaveButton));
	connect ( writeToGDML , SIGNAL(clicked()), this, SLOT( set_gdml_name() ));
	writeSingleVolumeButtons->addWidget(writeToGDML);

	writeToWRL = new QPushButton(tr("detector to wrl"));
	writeToWRL->setEnabled(false);
	writeToWRL->setToolTip("Creates a gdml file for the highlighted detector");
	writeToWRL->setIcon(style()->standardIcon(QStyle::SP_DialogSaveButton));
	connect ( writeToWRL , SIGNAL(clicked()), this, SLOT( set_wrl_name() ));
	writeSingleVolumeButtons->addWidget(writeToWRL);

	
	// writing all volumes in GDML or WRL
	QHBoxLayout *writeAllVolumeButtons = new QHBoxLayout();
	
	QPushButton *writeToGDMLA = new QPushButton(tr("all to gdml"));
	writeToGDMLA->setToolTip("Creates a gdml file for all detectors");
	writeToGDMLA->setIcon(style()->standardIcon(QStyle::SP_DialogSaveButton));
	connect ( writeToGDMLA , SIGNAL(clicked()), this, SLOT( set_gdml_nameAll() ));
	writeAllVolumeButtons->addWidget(writeToGDMLA);

	QPushButton *writeToWRLA = new QPushButton(tr("all to wrl"));
	writeToWRLA->setToolTip("Creates a wrl file for all detectors");
	writeToWRLA->setIcon(style()->standardIcon(QStyle::SP_DialogSaveButton));
	connect ( writeToWRLA , SIGNAL(clicked()), this, SLOT( set_wrl_nameAll() ));
	writeAllVolumeButtons->addWidget(writeToWRLA);
	
	// inspect detectors will open a new window
	showDetInNewWindow= new QPushButton(tr("inspect detector"));
	showDetInNewWindow->setEnabled(false);
	showDetInNewWindow->setToolTip("Creates a gdml file for all detectors");
	showDetInNewWindow->setIcon(style()->standardIcon(QStyle::SP_ComputerIcon));
	connect ( showDetInNewWindow , SIGNAL(clicked()), this, SLOT( inspectDetector() ) );

	
	// adding all together on the left layout
	leftLayout->addLayout(writeSingleVolumeButtons);
	leftLayout->addLayout(writeAllVolumeButtons);
	leftLayout->addWidget(showDetInNewWindow);
	
	
	QWidget *rightWidget = new QWidget(splitter);
	QVBoxLayout *rightLayout = new QVBoxLayout(rightWidget);
	rightLayout->addWidget(dTab);

		
	// all layouts
	QVBoxLayout *mainLayout = new QVBoxLayout;
	mainLayout->addWidget(splitter);
	setLayout(mainLayout);
	
}

detector_tree::~detector_tree()
{
	string hd_msg = gemcOpt.optMap["LOG_MSG"].args ;
	double VERB   = gemcOpt.optMap["GEO_VERBOSITY"].arg ;
	if(VERB>2)
		cout << hd_msg << " Detector Widget Tree Deleted." << endl;
	
	delete Switch_visibility;
	delete Switch_wiresolid;
	delete treeWidget;
}

map<string, tree_item> detector_tree::read_geometry(QTreeWidget *motherWidget)
{
	string hd_msg = gemcOpt.optMap["LOG_MSG"].args + " Detector Tree >> " ;
	double VERB   = gemcOpt.optMap["GEO_VERBOSITY"].arg ;
	
	map<string, tree_item> tree_map;
	tree_item item;
	
	
	if(VERB > 2)
	{
		cout << endl;
		cout << hd_msg <<  " Building Detector Widget Tree from Geometry STL Map..." << endl << endl;
		cout << hd_msg << " Experimenta Hall has " << Hall_Map->size() << " components." << endl << endl;
	}
	
	
	// fills the volumes tree map
	// initialize the Hall_Map member QTreeWidgetItem pointer and built variable to zero
	for(map<string, detector>::iterator i = Hall_Map->begin(); i != Hall_Map->end(); i++)
	{
		item.volume    = i->second.name;
		item.mother    = i->second.mother;
		item.treeItem  = NULL;
		item.scanned   = 0;
		item.sensitive = 0;
		if(i->second.sensitivity != "no")
			item.sensitive = 1;
		
		item.exist     = i->second.exist;
		item.visible   = i->second.visible;
		item.wiresolid = i->second.style;

//		if(i->second.GetPhysical() ) cout << i->second.name << " has physical " << i->second.GetPhysical() << endl;
//		else cout << i->second.name << " has NOT physical " << endl;

		if(i->second.name != "root" && i->second.GetPhysical())
		{
			tree_map.insert(map<string, tree_item>::value_type(i->second.name, item));
			if(VERB > 3) cout << hd_msg << " Mapping Detector " << i->first << endl;
		}
	}
	
	if(VERB > 3) cout << endl << endl;
	
	ActiveGrad = QLinearGradient(QPointF(1, 100), QPointF(180, 70));
	ActiveGrad.setColorAt(0, QColor( 80,  80, 255));
	ActiveGrad.setColorAt(1, QColor(245, 245, 245));
	
	SensitiveGrad = QLinearGradient(QPointF(1, 100), QPointF(180, 20));
	SensitiveGrad.setColorAt(0, QColor(255, 80,  80));
	SensitiveGrad.setColorAt(1, QColor(245, 245, 245));
	
	NonActiveGrad = QLinearGradient(QPointF(1, 100), QPointF(10, 30));
	NonActiveGrad.setColorAt(0, QColor(255, 20, 20));
	NonActiveGrad.setColorAt(1, QColor(245, 245, 245));
	
	NonVisibleGrad = QLinearGradient(QPointF(1, 100), QPointF(60, 30));
	NonVisibleGrad.setColorAt(0, QColor(5, 5, 5));
	NonVisibleGrad.setColorAt(1, QColor(240, 240, 240));
	
	ActiveBrush     = QBrush(ActiveGrad);
	SensitiveBrush  = QBrush(SensitiveGrad);
	NonActiveBrush  = QBrush (NonActiveGrad);
	NonVisibleBrush = QBrush (NonVisibleGrad);
	
	
	// now build the actual QTreeWidgetItem
	// if a kid has no mom, its mom will be built first
	vector<string> relatives;
	string mom, kid;
	if(VERB > 2) cout << endl << hd_msg << " Reordering Tree Elements..." << endl << endl;
	for(map<string, tree_item>::iterator i = tree_map.begin(); i != tree_map.end(); i++)
	{
		if(i->first != "" && i->first != "root")     relatives.push_back(i->second.volume);
		
		while(relatives.size() > 0 && relatives.size() < 100)
		{
			mom = tree_map[relatives.back()].mother;
			kid = tree_map[relatives.back()].volume;
			
			
			if(VERB > 3)
			{
				for(unsigned int i=0; i<relatives.size()-1; i++) cout << "\t";
				cout << hd_msg << " Checking " << kid << ", child of " << mom
				<< ", for a living ancestor. This Geneaology Depth is " << relatives.size() << "." << endl;
			}
			
			// Mom is <root>. Can build the kid.
			if(tree_map[kid].scanned == 0 && mom.find("root", 0) != string::npos)
			{
				if(VERB > 3)
				{
					for(unsigned int i=0; i<relatives.size()-1; i++) cout << "\t";
					cout << hd_msg << "  Found:  " << relatives.back()
					<< " has <root> as mom and it's not built yet. Building " << kid << "..."  << endl;
				}
				tree_map[kid].treeItem = new QTreeWidgetItem(motherWidget);
				tree_map[kid].treeItem->setText(0,  kid.c_str());
				tree_map[kid].scanned = 1;
				
			}
			
			// Mom (not root) is built, kid not built yet. Build the kid.
			if(tree_map[kid].scanned == 0 && tree_map[mom].scanned == 1)
			{
				if(VERB > 3)
				{
					for(unsigned int i=0; i<relatives.size()-1; i++) cout << "\t";
					cout << hd_msg << "  Found:  " << kid
					<< " is not built yet but its mommie " << mom << " is. Building " << kid << "..."  << endl;
				}
				
				tree_map[kid].treeItem = new QTreeWidgetItem(tree_map[mom].treeItem);
				tree_map[kid].treeItem->setText(0,  kid.c_str());
				tree_map[kid].scanned = 1;
			}
			
			// if the kid still doesn't exists it means its mom doesn't exist. Need to go up one level.
			if(tree_map[kid].scanned == 0)  relatives.push_back(mom);
			
			// the kid has been built. Can go down one step in geneaology
			if(tree_map[kid].scanned == 1 && relatives.size())
			{
				if(VERB > 3)
					cout << hd_msg  << " " <<  kid << " is built." <<  endl << endl;
				
				if(tree_map[kid].exist == 1)     tree_map[kid].treeItem->setBackground(0, ActiveBrush );
				if(tree_map[kid].exist == 0)     tree_map[kid].treeItem->setBackground(0, NonActiveBrush );
				if(tree_map[kid].visible == 1)   tree_map[kid].treeItem->setBackground(0, ActiveBrush );
				if(tree_map[kid].sensitive == 1) tree_map[kid].treeItem->setBackground(0, SensitiveBrush );
				if(tree_map[kid].visible == 0)   tree_map[kid].treeItem->setBackground(0, NonVisibleBrush );
				relatives.pop_back();
			}
		}
	}
	
	return tree_map;
}


void detector_tree::switch_visibility()
{
	QTreeWidgetItem *CurrentItem = treeWidget->currentItem();
	string name  = qs_tostring(CurrentItem->text(0));
	tree_map[name].visible == 1 ? tree_map[name].visible = 0 : tree_map[name].visible = 1;
	tree_map[name].visible == 1 ? CurrentItem->setBackground(0, ActiveBrush ) : CurrentItem->setBackground(0, NonVisibleBrush);
	char command[100];
	sprintf(command, "/vis/geometry/set/visibility %s 0 %d", name.c_str(), tree_map[name].visible);
	UImanager->ApplyCommand(command);
}

void detector_tree::switch_visibility_daughters()
{
	string name  = qs_tostring(treeWidget->currentItem()->text(0));
	int vis = 0;
	for(map<string, tree_item>::iterator it = tree_map.begin(); it != tree_map.end(); it++)
		if(it->second.mother == name)
		{
			it->second.visible == 1 ? it->second.visible = 0 : it->second.visible = 1;
			it->second.visible == 1 ? it->second.treeItem->setBackground(0, ActiveBrush )
			: it->second.treeItem->setBackground(0, NonVisibleBrush);
			vis = it->second.visible;
		}
	
	char command[100];
	sprintf(command, "/vis/geometry/set/daughtersInvisible %s 0 %d", name.c_str(), !vis);
	UImanager->ApplyCommand(command);	
}

void detector_tree::inspectDetector()
{
	string name  = qs_tostring(treeWidget->currentItem()->text(0));
	
	char command[100];
	sprintf(command,"/vis/viewer/set/lineSegmentsPerCircle 100 ");
	UImanager->ApplyCommand(command);
	sprintf(command, "/vis/open OGL");
	UImanager->ApplyCommand(command);
	sprintf(command, "/vis/specify %s", name.c_str());
	UImanager->ApplyCommand(command);
	sprintf(command, "/vis/viewer/set/background .85 .95 .98 1");
	UImanager->ApplyCommand(command);
}





void detector_tree::switch_wiresolid()
{
	string name  = qs_tostring(treeWidget->currentItem()->text(0));
	tree_map[name].wiresolid == 1 ? tree_map[name].wiresolid = 0 : tree_map[name].wiresolid = 1;
	
	char command[100];
	if( tree_map[name].wiresolid == 1) sprintf(command, "/vis/geometry/set/forceSolid     %s 0 1", name.c_str());
	if( tree_map[name].wiresolid == 0) sprintf(command, "/vis/geometry/set/forceWireframe %s 0 1", name.c_str());
	UImanager->ApplyCommand(command);
}


void detector_tree::switch_color()
{
	int r, g, b;
	string name  = qs_tostring(treeWidget->currentItem()->text(0));
	QColor color = QColorDialog::getColor(Qt::green, this);
	color.getRgb(&r, &g, &b);
	char command[100];
	sprintf(command, "/vis/geometry/set/colour %s 0 %3.2f %3.2f %3.2f", name.c_str(), r/255.0, g/255.0, b/255.0);
	UImanager->ApplyCommand(command);
}



// maybe add a apply button to apply the changes?
void detector_tree::show_detector()
{
	QTreeWidgetItem *CurrentItem = treeWidget->currentItem();
	string name  = qs_tostring(CurrentItem->text(0));
	
	detector detect = (*Hall_Map)[name];
	dTab->update_detector(&detect);

	writeToGDML->setEnabled(true);
	writeToWRL->setEnabled(true);
	showDetInNewWindow->setEnabled(true);

	
		// have to overload "=" to make this happen
//		 if((*Hall_Map)[name] != detect)
//		cout << hd_msg << " Old detector: "  <<  (*Hall_Map)[name];
		(*Hall_Map)[name] = detect;
//		cout << hd_msg << " New detector: "  <<  (*Hall_Map)[name];
	
		// need mother to remove it
	//	string mother = (*Hall_Map)[name].mother;
	//	//(*Hall_Map)[name].create_solid(gemcOpt, Hall_Map);
	//	//(*Hall_Map)[name].create_logical_volume(MMats, gemcOpt);
	//	//(*Hall_Map)[mother].RemoveDaughter( (*Hall_Map)[name].GetPhysical());
	//	//(*Hall_Map)[name].create_physical_volumes(gemcOpt, (*Hall_Map)[mother].GetLogical());
		(*Hall_Map)[name].SetTranslation(detect.pos);
	//
	//// runManager->DefineWorldVolume(  (*Hall_Map)["root"].GetPhysical() );
	//// runManager->GeometryHasBeenModified();
	//
	//	char command[100];
	//	sprintf(command, "/vis/geometry/set/visibility root 0 0");
	//// sprintf(command, "/vis/geometry/set/visibility %s 0 %d", name.c_str(), 1);
	//	UImanager->ApplyCommand(command);
	
}


void detector_tree::set_gdml_name()
{
	write_gdml_file(qs_tostring(treeWidget->currentItem()->text(0)));
}

void detector_tree::set_gdml_nameAll()
{
	write_gdml_file("root");
}

void detector_tree::set_wrl_name()
{
	write_wrl_file(qs_tostring(treeWidget->currentItem()->text(0)));
}

void detector_tree::set_wrl_nameAll()
{
	write_wrl_file("root");
}



void detector_tree::write_gdml_file(string name)
{
	string fileout = name + ".gdml";
	
	detector detect = (*Hall_Map)[name];
	
	struct stat stFileInfo;
	if(stat(fileout.c_str(),&stFileInfo)==0)   // Check if file exists already.
	{
		if(remove(fileout.c_str())){ // Remove file, otherwise parser throws trap.
			cout << "ERROR -- Could not remove file " << fileout << " cannot write new one.\n";
			return;
		}
	}
	
	G4GDMLParser parser;
	parser.Write(fileout.c_str(),detect.GetPhysical(),false);	
}



void detector_tree::write_wrl_file(string name)
{
	char command[100];

	if(name != "root")
	{
		sprintf(command, "/vis/specify %s", name.c_str());
		UImanager->ApplyCommand(command);
	}

	sprintf(command, "/vis/open VRML2FILE");
	UImanager->ApplyCommand(command);
	
	sprintf(command, "/vis/viewer/flush");
	UImanager->ApplyCommand(command);
}


				
void detector_tree::change_placement()
{
	cout << " YAY! " << endl;
}


void detector_tree::createActions()
{
	// Toggle visibility
	Switch_visibility = new QAction(tr("&Switch Visibility"), this);
	Switch_visibility->setShortcut(tr("Ctrl+V"));
	Switch_visibility->setStatusTip(tr("Switch Visibility"));
	connect(Switch_visibility, SIGNAL(triggered()), this, SLOT(switch_visibility()));
	
	// Toggle daughters' visibility
	Switch_visibility_daughters = new QAction(tr("&Switch Visibility of Daughters"), this);
	Switch_visibility_daughters->setShortcut(tr("Ctrl+D"));
	Switch_visibility_daughters->setStatusTip(tr("Switch Visibility for the daughters"));
	connect(Switch_visibility_daughters, SIGNAL(triggered()), this, SLOT(switch_visibility_daughters()));
	
	// Toggle wireframe/solid
	Switch_wiresolid = new QAction(tr("&Switch Wireframe / Solid"), this);
	Switch_wiresolid->setShortcut(tr("Ctrl+W"));
	Switch_wiresolid->setStatusTip(tr("Switch Wireframe / Solid"));
	connect(Switch_wiresolid, SIGNAL(triggered()), this, SLOT(switch_wiresolid()));
	
	// Pick Volume Color
	Switch_color = new QAction(tr("&Pick Volume Color"), this);
	Switch_color->setShortcut(tr("Ctrl+L"));
	Switch_color->setStatusTip(tr("Pick Volume Color"));
	connect(Switch_color, SIGNAL(triggered()), this, SLOT(switch_color()));
    
	// Write GDML File
	Write_GDML_File = new QAction(tr("Write GDML File"), this);
	Write_GDML_File->setShortcut(tr("Ctrl+G"));
	Write_GDML_File->setStatusTip(tr("Write file name.gdml to current directory"));
	connect(Write_GDML_File, SIGNAL(triggered()), this, SLOT(set_gdml_name()));

	// Write WRL File
	Write_WRL_File = new QAction(tr("Write WRL File"), this);
	Write_WRL_File->setShortcut(tr("Ctrl+R"));
	Write_WRL_File->setStatusTip(tr("Write file g4_#.wrl to current directory"));
	connect(Write_WRL_File, SIGNAL(triggered()), this, SLOT(set_wrl_name()));

}



