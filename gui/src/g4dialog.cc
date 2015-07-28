// Qt headers
#include <QtWidgets>

// gemc headers
#include "g4dialog.h"

// G4 headers
#include "G4UIcommandTree.hh"

g4dialog::g4dialog(QWidget *parent, goptions *Opts) : QWidget(parent)
{
	gemcOpt = Opts;
	UImanager  = G4UImanager::GetUIpointer();
	
	//  Layout:
	//
	//  + +-------------------+ +
	//  | |         |         | |
	//  | |  Tree   |  Help   | |
	//  | |         |         | |
	//  | +-------------------+ |
	//  | +-------------------+ |
	//  | |                   | |
	//  | |   promt history   | |
	//  | +-------------------+ |
	//  | +-------------------+ |
	//  | |   > promt area    | |
	//  | +-------------------+ |
	//  +-----------------------+
	
	
	// Vertical Splitter - Top and Bottom layouts
	QSplitter *splitter     = new QSplitter(Qt::Vertical);
	
	
	// top layout
	QWidget   *topWidget    = new QWidget(splitter);
	QSplitter *treesplitter = new QSplitter(Qt::Horizontal);
	
	
	// Tree + Help VBOX assembly
	QVBoxLayout *helpLayout = new QVBoxLayout(treesplitter);
	
	// Left: The help tree
	fHelpTreeWidget = new QTreeWidget();
	fHelpTreeWidget = CreateHelpTree();
	if(fHelpTreeWidget)
    helpLayout->addWidget(fHelpTreeWidget);
	
	
	// Right: the help on individual commands
	fHelpArea = new QTextEdit(treesplitter);
	fHelpArea->setReadOnly(true);
	
	
	// treesplitter size
	QList<int> tlist;
	tlist.append( 400 );
	tlist.append( 400 );
	treesplitter->setSizes(tlist);
	
	
	QVBoxLayout *layoutTop = new QVBoxLayout(topWidget);
	layoutTop->addWidget(treesplitter);
	
	// bottom layout
	QWidget     *bottomWidget = new QWidget(splitter);
	
	// history area
	fCommandHistoryArea = new QListWidget();
	fCommandHistoryArea->setSelectionMode(QAbstractItemView::SingleSelection);
	fCommandHistoryArea->installEventFilter(this);
	connect(fCommandHistoryArea, SIGNAL(itemSelectionChanged()), SLOT(CommandHistoryCallback()));
	
	
	// manual command
	// QLabel *fCommandLabel = new QLabel("Press Enter to Execute Command:");
	fCommandArea = new QLineEdit();
	fCommandArea->installEventFilter(this);
	fCommandArea->activateWindow();
	fCommandArea->setFocusPolicy ( Qt::StrongFocus );
	fCommandArea->setFocus(Qt::TabFocusReason);
	connect(fCommandArea, SIGNAL(returnPressed()), SLOT(CommandEnteredCallback()));
	
	// putting all together
	QVBoxLayout *layoutBottom = new QVBoxLayout(bottomWidget);
	layoutBottom->addWidget(new QLabel("History:"));
	layoutBottom->addWidget(fCommandHistoryArea);
	layoutBottom->addWidget(new QLabel("Press Enter to Execute Command:"));
	layoutBottom->addWidget(fCommandArea);
	
	
	// splitter size
	QList<int> list;
	list.append( 500 );
	list.append( 300 );
	splitter->setSizes(list);
	
	
	
	// all layouts
	QVBoxLayout *mainLayout = new QVBoxLayout;
	mainLayout->addWidget(splitter);
	setLayout(mainLayout);	
}




g4dialog::~g4dialog()
{
	string hd_msg = gemcOpt->optMap["LOG_MSG"].args ;
	double VERB   = gemcOpt->optMap["GEO_VERBOSITY"].arg ;
	if(VERB>2)
		cout << hd_msg << " g4 Dialog Widget Deleted." << endl;
	
}


// execute history item
void g4dialog::CommandHistoryCallback()
{
	QListWidgetItem* item =  NULL;
	if (!fCommandHistoryArea)
    return ;
	
	
	
	QList<QListWidgetItem *> list = fCommandHistoryArea->selectedItems();
	if(list.isEmpty())
    return;
	
	item = list.first();
	if(!item)
    return;
	
	fCommandArea->setText(item->text());
	
}



// execute G4 manual command
void g4dialog::CommandEnteredCallback()
{
	
	if(fCommandArea->text().trimmed() != "")
	{
		fCommandHistoryArea->addItem(fCommandArea->text());
		
		UImanager->ApplyCommand(fCommandArea->text().toStdString().c_str());
		
		fCommandHistoryArea->clearSelection();
		fCommandHistoryArea->setCurrentItem(NULL);
		fCommandArea->setText("");
	}
}




QTreeWidget* g4dialog::CreateHelpTree()
{
	if(UImanager==NULL) return NULL;
	G4UIcommandTree *treeTop = UImanager->GetTree();
	
	
	// build widget
	fHelpTreeWidget->setSelectionMode(QAbstractItemView::SingleSelection);
	
	QStringList labels;
	labels << QString("Command");
	fHelpTreeWidget->setHeaderLabels(labels);
	
	G4int treeSize = treeTop->GetTreeEntry();
	QTreeWidgetItem * newItem;
	
	for (int a=0;a<treeSize;a++)
	{
		// Creating new item
		newItem = new QTreeWidgetItem(fHelpTreeWidget);
		newItem->setText(0, QString((char*)(treeTop->GetTree(a+1)->GetPathName()).data()).trimmed());
		
		// look for childs
		CreateChildTree(newItem, treeTop->GetTree(a+1));
	}
	
	
	connect(fHelpTreeWidget, SIGNAL(itemSelectionChanged() ),                    this, SLOT(HelpTreeClicCallback() ) );
	connect(fHelpTreeWidget, SIGNAL(itemDoubleClicked (QTreeWidgetItem*, int) ), this, SLOT(HelpTreeDoubleClicCallback()));
	
	return fHelpTreeWidget;
}




void g4dialog::CreateChildTree(QTreeWidgetItem *aParent,G4UIcommandTree *aCommandTree)
{
	if (aParent == NULL) return;
	if (aCommandTree == NULL) return;
	
	
	// Creating new item
	QTreeWidgetItem * newItem;
	
	// Get the Sub directories
	for (int a=0;a<aCommandTree->GetTreeEntry();a++)
	{
		
		newItem = new QTreeWidgetItem(aParent);
		newItem->setText(0,QString((char*)(aCommandTree->GetTree(a+1)->GetPathName()).data()).trimmed());
		
		CreateChildTree(newItem,aCommandTree->GetTree(a+1));
	}
	
	// Get the Commands
	for (int a=0;a<aCommandTree->GetCommandEntry();a++)
	{
		QStringList stringList;
		newItem = new QTreeWidgetItem(aParent);
		newItem->setText(0, QString((char*)(aCommandTree->GetCommand(a+1)->GetCommandPath()).data()).trimmed());
		newItem->setExpanded(false);
		
	}
}


// displays help on the right help area
void g4dialog::HelpTreeClicCallback()
{
	QTreeWidgetItem* item =  NULL;
	if(!fHelpTreeWidget || !fHelpArea)
    return;
	
	QList<QTreeWidgetItem *> list = fHelpTreeWidget->selectedItems();
	if(list.isEmpty())
    return;
	
	item = list.first();
	if(!item)
    return;
	
	
	if(UImanager==NULL) return;
	G4UIcommandTree *treeTop = UImanager->GetTree();
	
	string itemText = item->text(0).toStdString();
	G4UIcommand* command = treeTop->FindPath(itemText.c_str());
	
	// if it's a valid command, display the help
	if(command)
	{
		fHelpArea->setText(GetCommandList(command));
	}
	else
	{
		// this is a command
		G4UIcommandTree* path = treeTop->FindCommandTree(itemText.c_str());
		if(path)
		{
			// this is not a command, this is a sub directory
			// We display the Title
			fHelpArea->setText(path->GetTitle().data());
		}
	}
}


QString g4dialog::GetCommandList (const G4UIcommand *aCommand)
{
	
	QString txt ="";
	if (aCommand == NULL)
    return txt;
	
	G4String commandPath   = aCommand->GetCommandPath();
	G4String rangeString   = aCommand->GetRange();
	G4int n_guidanceEntry  = aCommand->GetGuidanceEntries();
	G4int n_parameterEntry = aCommand->GetParameterEntries();
	
	if ((commandPath == "") &&
		(rangeString == "") &&
		(n_guidanceEntry == 0) &&
		(n_parameterEntry == 0)) {
		return txt;
	}
	
	if((commandPath.length()-1)!='/')
	{
		txt += "Command " + QString((char*)(commandPath).data()) + "\n";
	}
	txt += "Guidance :\n";
	
	for( G4int i_thGuidance=0; i_thGuidance < n_guidanceEntry; i_thGuidance++ )
	{
		txt += QString((char*)(aCommand->GetGuidanceLine(i_thGuidance)).data()) + "\n";
	}
	if( ! rangeString.isNull() )
	{
		txt += " Range of parameters : " + QString((char*)(rangeString).data()) + "\n";
	}
	if( n_parameterEntry > 0 )
	{
		G4UIparameter *param;
		
		// Re-implementation of G4UIparameter.cc
		for( G4int i_thParameter=0; i_thParameter<n_parameterEntry; i_thParameter++ )
		{
			param = aCommand->GetParameter(i_thParameter);
			txt += "\nParameter : " + QString((char*)(param->GetParameterName()).data()) + "\n";
			if( ! param->GetParameterGuidance().isNull() )
			txt += QString((char*)(param->GetParameterGuidance()).data())+ "\n" ;
			txt += " Parameter type  : " + QString(QChar(param->GetParameterType())) + "\n";
			if(param->IsOmittable())
			{
				txt += " Omittable       : True\n";
			}
			else
			{
				txt += " Omittable       : False\n";
			}
			if( param->GetCurrentAsDefault() )
			{
				txt += " Default value   : taken from the current value\n";
			}
			else if( ! param->GetDefaultValue().isNull() )
			{
				txt += " Default value   : " + QString((char*)(param->GetDefaultValue()).data())+ "\n";
			}
			if( ! param->GetParameterRange().isNull() )
			{
				txt += " Parameter range : " + QString((char*)(param->GetParameterRange()).data())+ "\n";
			}
			if( ! param->GetParameterCandidates().isNull() )
			{
				txt += " Candidates      : " + QString((char*)(param->GetParameterCandidates()).data())+ "\n";
			}
		}
	}
	return txt;
}



void g4dialog::HelpTreeDoubleClicCallback()
{
	HelpTreeClicCallback();
	
	QTreeWidgetItem* item =  NULL;
	
	if(!fHelpTreeWidget || !fHelpArea)
    return ;
	
	QList<QTreeWidgetItem *> list = fHelpTreeWidget->selectedItems();
	if(list.isEmpty())
    return;
	
	item = list.first();
	if(!item)
    return;
	
	
	fCommandArea->clear();
	fCommandArea->setText(item->text(0));
}









