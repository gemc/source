// Qt headers
#include <QtWidgets>

// gemc headers
#include "gsignal.h"
#include "Hit.h"
#include "string_utilities.h"

// G4 headers
#include "G4UIcommandTree.hh"
#include "Randomize.hh"

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;

gsignal::gsignal(QWidget *parent, goptions *Opts, map<string, sensitiveDetector*> SD_Map) : QWidget(parent)
{
	gemcOpt = Opts;
	WRITE_INTRAW     = replaceCharWithChars(gemcOpt->optMap["INTEGRATEDRAW"].args, ",", "  ");
	SAVE_ALL_MOTHERS = gemcOpt->optMap["SAVE_ALL_MOTHERS"].arg;
	
	UImanager  = G4UImanager::GetUIpointer();
	SeDe_Map   = SD_Map;
	
	//  Layout:
	//
	//  + +-----------+------+--------+
	//  | |           |      |        |
	//  | |  Element  | Var  | Var    |
	//  | |   Choice  | List | Choice |
	//  | |           |      |        |
	//  | +-----------+------+--------+
	//  | +---------------------------+
	//  | |                           |
	//  | |           Graph           |
	//  | +---------------------------|
	//  +-----------------------------+
	
	// Vertical Splitter - Top and Bottom layouts
	QSplitter *splitter     = new QSplitter(Qt::Vertical);
	
	// top layout
	QWidget   *topWidget    = new QWidget(splitter);
	QSplitter *treesplitter = new QSplitter(Qt::Horizontal);
	
	QVBoxLayout *topLayout = new QVBoxLayout(treesplitter);
	
	// Left: The sensitive detectors hits tree
	s_detectors = new QTreeWidget();
	s_detectors = CreateSDetsTree();
	if(s_detectors)
		topLayout->addWidget(s_detectors);
	
	
	// Center: The individual signals choice
	gsignals = new QTreeWidget();
	gsignals = CreateSignalsTree();
	if(gsignals)
		topLayout->addWidget(gsignals);
	
	
	availableSignals.push_back("E Dep.");
	availableSignals.push_back("Trk ID");
	availableSignals.push_back("Orig. Trk");
	availableSignals.push_back("Mom");
	availableSignals.push_back("<x>");
	availableSignals.push_back("<y>");
	availableSignals.push_back("<z>");
	availableSignals.push_back("<lx>");
	availableSignals.push_back("<ly>");
	availableSignals.push_back("<lz>");
	availableSignals.push_back("<vx>");
	availableSignals.push_back("<vy>");
	availableSignals.push_back("<vz>");
	availableSignals.push_back("<mvx>");
	availableSignals.push_back("<mvy>");
	availableSignals.push_back("<mvz>");
	availableSignals.push_back("Voltage");
	
	qcsignal = new QComboBox;
	for(unsigned v=0; v<availableSignals.size(); v++)
		qcsignal->addItem(tr(availableSignals[v].c_str()));
	
	qcsignal->setMaximumWidth(300);
	connect ( qcsignal   , SIGNAL( currentIndexChanged (int) ), this, SLOT( chooseVariable(int)    ) );
	
	// treesplitter size
	QList<int> tlist;
	tlist.append( 300 );
	tlist.append( 300 );
	treesplitter->setSizes(tlist);
	
	QVBoxLayout *layoutTop = new QVBoxLayout(topWidget);
	layoutTop->addWidget(qcsignal);
	layoutTop->addWidget(treesplitter);
	
	
	// bottom layout
	QWidget     *bottomWidget = new QWidget(splitter);
	bottomWidget->setMinimumSize(600, 450);
	bottomWidget->setMaximumSize(600, 450);
	
	graphView = new graph(this);
	
	// Graph axis origins and legth, number of ticks each axis
	//                 xorig  yorig  xlength ylength nticksx nticksy
	graphView->setAxis(15,   350,    445,    310,     5,     5);
	// inside shift of the axis ticks, and factor that multiplies the ticks size
	graphView->setInside(20, 2, 2);
	// particle colors follow gemc settings
	// red:   positive
	// gray:  neutral
	// green: negative
	//
	// second argument of QPen is thickness of pencil
	graphView->pcolors[2112] = QPen(Qt::black,             3);   // neutrons: black
	graphView->pcolors[22]   = QPen(Qt::blue,              3);   // photons: blue
	graphView->pcolors[11]   = QPen(Qt::cyan,              3);   // electrons: cyan
	graphView->pcolors[2212] = QPen(QColor(240, 80, 80),   3);   // protons: orange
	graphView->pcolors[211]  = QPen(Qt::magenta,           3);   // pi+: magenta
	graphView->pcolors[-211] = QPen(Qt::yellow,            3);   // pi-: yellow
	graphView->pcolors[-11]  = QPen(Qt::red,               3);   // positrons: positive - red
	graphView->pcolors[0]    = QPen(Qt::blue,              3);   // optical photons: blue
	graphView->pcolors[13  ] = QPen(QColor(0,125, 0   ),   3);   // Muon+ - dark green
	graphView->pcolors[-13 ] = QPen(QColor(0,250, 0   ),   3);   // Muon- - light green
	graphView->pcolors[1000] = QPen(Qt::black,             3);   // neutrons: black
	
	
	// putting all together
	QVBoxLayout *layoutBottom = new QVBoxLayout(bottomWidget);
	layoutBottom->addWidget(new QLabel("Signal:"));
	layoutBottom->addWidget(graphView);
	
	
	// splitter size
	QList<int> list;
	list.append( 420 );
	list.append( 450 );
	splitter->setSizes(list);
	
	
	// all layouts
	QVBoxLayout *mainLayout = new QVBoxLayout;
	mainLayout->addWidget(splitter);
	setLayout(mainLayout);
	
	
	// Has Hit Gradient Color
	HitGrad = QLinearGradient(QPointF(1, 100), QPointF(180, 20));
	HitGrad.setColorAt(0, QColor(255, 80,  80));
	HitGrad.setColorAt(1, QColor(245, 245, 245));
	HitBrush  = QBrush(HitGrad);
	
	// initialize signalChoice to Energy deposited
	signalChoice = "E Dep.";
}

QTreeWidget* gsignal::CreateSDetsTree()
{
	s_detectors->clear();
	s_detectors->setSelectionMode(QAbstractItemView::SingleSelection);
	QStringList labels;
	labels << QString("Hits List");
	s_detectors->setHeaderLabels(labels);
	
	
	for(map<string, sensitiveDetector*>::iterator it = SeDe_Map.begin(); it!= SeDe_Map.end(); it++)
	{
		MHitCollection *MHC = it->second->GetMHitCollection();
		int nhits = 0;
		if(MHC) nhits = MHC->GetSize();
		
		// Creating sensitive detectors name tree if it's different than "mirrors"
		if(it->first != "mirror")
		{
			QTreeWidgetItem *newItem = new QTreeWidgetItem(s_detectors);
			newItem->setText(0, QString(it->first.c_str()));
			
			if(nhits)
			{
				newItem->setBackground(0, HitBrush );
				string snhits = it->first + "   " + stringify(nhits)  + " hit";
				if(nhits>1) snhits += "s";
				newItem->setText(0, QString(snhits.c_str()));
				
				// if the last sensitive identifier is nphe then
				// visualization screen is different:
				// need to visualize number of photoelectrons only
				for(int h=0; h<nhits; h++)
				{
					MHit *aHit = (*MHC)[h];
					
					int nsteps = aHit->GetPos().size();
					QTreeWidgetItem  *newHit = new QTreeWidgetItem(newItem);
					
					string hitindex = "Hit n. " + stringify(h+1) + "  nsteps: " +  stringify(nsteps) ;
					
					// if last sensitive element is nphe_pmt then writing number of photo-electrons
					if(!it->second->SDID.identifiers.size())
						cout << "   !!! Error: no identifiers found for SD >" << it->second->SDID.name  << "<" << endl;
					
					if(it->second->SDID.identifiers.back().find("nphe") != string::npos)
						hitindex = "Hit n. " + stringify(h+1) + "  nphe: " +  stringify(nsteps) ;
					
					newHit->setText(0, QString(hitindex.c_str()));
				}
			}
		}
	}
	connect(s_detectors, SIGNAL(itemSelectionChanged() ),   this, SLOT(CreateSignalsTree() ) );
	return s_detectors;
}



void gsignal::chooseVariable(int index)
{
	signalChoice = availableSignals[index];
	CreateSignalsTree();
	
}



QTreeWidget* gsignal::CreateSignalsTree()
{
	gsignals->clear();
	gsignals->setSelectionMode(QAbstractItemView::SingleSelection);
	QStringList labels;
	labels << QString("Signal");
	gsignals->setHeaderLabels(labels);
	
	QTreeWidgetItem* item =  NULL;
	if(s_detectors)
	{
		QList<QTreeWidgetItem *> list = s_detectors->selectedItems();
		if(!list.isEmpty())
		{
			item = list.first();
			
			// trick to look for the index of this item
			// will return zero if it's the sensitive detector name
			// this index is also the index of the hit
			
			// itext looks like this:
			// Hit n. 1  nsteps: 6
			string itext = item->text(0).toStdString();
			
			stringstream sitext(itext);
			string a, b, c;
			
			// index is the hit index
			unsigned int index;
			sitext >> a >> b >> c;
			index = atoi(c.c_str());
			
			string SD;
			string sDetector;
			sensitiveDetector *MSD = NULL;
			if(index==0)
			{
				SD = a;
				sDetector = SD;
				MSD = SeDe_Map[SD];
			}
			else
			{
				stringstream mtext(item->parent()->text(0).toStdString());
				mtext >> SD;
				sDetector = SD;
				MSD = SeDe_Map[SD];
				SD += "  Hit n. " + c;
			}
			
			QTreeWidgetItem * newHit = new QTreeWidgetItem(gsignals);
			newHit->setText(0, QString(SD.c_str()));
			newHit->setExpanded(1);
			
			MHitCollection *MHC = MSD->GetMHitCollection();
			
			if(MHC)
				if(index>0 && index<=MHC->GetSize())
				{
					MHit* aHit = (*MHC)[index-1];
					int nsteps = aHit->GetPos().size();
					
					// photo electron tubes are treated differently than
					// normal sensitive detectors
					if(MSD->SDID.identifiers.back().find("nphe") != string::npos)
					{
						SD += "  nphe: " + stringify(nsteps);
						
						newHit->setText(0, QString(SD.c_str()));
						
						vector<double>        ene = aHit->GetEs();
						vector<double>       time = aHit->GetTime();
						vector<int>           pid = aHit->GetPIDs();
						
						vector<double> lambda;
						// 1 nm = 1240 / eV
						for(unsigned int i=0; i<ene.size(); i++)
							lambda.push_back(1240.0/(ene[i]/eV));
						
						
						vector<identifier> identi = aHit->GetId();
						string title = "";
						for(unsigned int i=0; i<identi.size(); i++)
							title += identi[i].name + " " + stringify(identi[i].id) + "   " ;
						
						QTreeWidgetItem * EneI = new QTreeWidgetItem(newHit);
						EneI->setText(0, QString("   Wavelenght[nm]    pid      Time[ns] "));
						EneI->setExpanded(1);
						
						QTreeWidgetItem * EneItems;
						for(int i=0; i<nsteps; i++)
						{
							EneItems = new QTreeWidgetItem(EneI);
							char etext[200];
							sprintf(etext, "      %4.1f             %d        %5.4f     ", lambda[i], pid[i], time[i]);
							EneItems->setText(0, QString(etext));
							EneItems->setTextAlignment(1, Qt::AlignJustify);
							EneItems->setTextAlignment(2, Qt::AlignJustify);
						}
						graphView->plots_bg("time [ns]", "Wavelenght [nm]", time, lambda, title);
						graphView->plot_graph(time, lambda, pid);
					}
					else
					{
						vector<double> signal;
						vector<double> time = aHit->GetTime();
						vector<int>     pid = aHit->GetPIDs();
						
						if(signalChoice == "E Dep.")
							signal = aHit->GetEdep();
						
						if(signalChoice == "Trk ID")
							signal = convertVintVdouble(aHit->GetTIds());
						
						
						// this info is available only if WRITE_INTRAW is enable
						if(signalChoice == "Orig. Trk")
						{
							if(WRITE_INTRAW.find(sDetector) != string::npos)
							{
								signal = convertVintVdouble(aHit->GetoTrackIds());
							}
							else
							{
								vector<int> trks = aHit->GetTIds();
								for(unsigned int i=0; i<trks.size(); i++)
									signal.push_back(-1.0);
							}
						}
						
						if(signalChoice == "<x>")
						{
							vector<G4ThreeVector> pos = aHit->GetPos();
							for(unsigned int i=0; i<pos.size(); i++)
								signal.push_back(pos[i].x());
						}
						if(signalChoice == "<y>")
						{
							vector<G4ThreeVector> pos = aHit->GetPos();
							for(unsigned int i=0; i<pos.size(); i++)
								signal.push_back(pos[i].y());
						}
						if(signalChoice == "<z>")
						{
							vector<G4ThreeVector> pos = aHit->GetPos();
							for(unsigned int i=0; i<pos.size(); i++)
								signal.push_back(pos[i].z());
						}
						if(signalChoice == "<lx>")
						{
							vector<G4ThreeVector> pos = aHit->GetLPos();
							for(unsigned int i=0; i<pos.size(); i++)
								signal.push_back(pos[i].x());
						}
						if(signalChoice == "<ly>")
						{
							vector<G4ThreeVector> pos = aHit->GetLPos();
							for(unsigned int i=0; i<pos.size(); i++)
								signal.push_back(pos[i].y());
						}
						if(signalChoice == "<lz>")
						{
							vector<G4ThreeVector> pos = aHit->GetLPos();
							for(unsigned int i=0; i<pos.size(); i++)
								signal.push_back(pos[i].z());
						}
						if(signalChoice == "<vx>")
						{
							vector<G4ThreeVector> pos = aHit->GetVerts();
							for(unsigned int i=0; i<pos.size(); i++)
								signal.push_back(pos[i].x());
						}
						if(signalChoice == "<vy>")
						{
							vector<G4ThreeVector> pos = aHit->GetVerts();
							for(unsigned int i=0; i<pos.size(); i++)
								signal.push_back(pos[i].y());
						}
						if(signalChoice == "<vz>")
						{
							vector<G4ThreeVector> pos = aHit->GetVerts();
							for(unsigned int i=0; i<pos.size(); i++)
								signal.push_back(pos[i].z());
						}
						if(signalChoice == "<mvx>")
						{
							if(SAVE_ALL_MOTHERS)
							{
								vector<G4ThreeVector> pos = aHit->GetmVerts();
								for(unsigned int i=0; i<pos.size(); i++)
									signal.push_back(pos[i].x());
							}
							else
							{
								vector<int> trks = aHit->GetTIds();
								for(unsigned int i=0; i<trks.size(); i++)
									signal.push_back(0.0);
							}
						}
						if(signalChoice == "<mvy>")
						{
							if(SAVE_ALL_MOTHERS)
							{
								vector<G4ThreeVector> pos = aHit->GetmVerts();
								for(unsigned int i=0; i<pos.size(); i++)
									signal.push_back(pos[i].y());
							}
							else
							{
								vector<int> trks = aHit->GetTIds();
								for(unsigned int i=0; i<trks.size(); i++)
									signal.push_back(0.0);
							}
						}
						if(signalChoice == "<mvz>")
						{
							if(SAVE_ALL_MOTHERS)
							{
								vector<G4ThreeVector> pos = aHit->GetmVerts();
								for(unsigned int i=0; i<pos.size(); i++)
									signal.push_back(pos[i].z());
							}
							else
							{
								vector<int> trks = aHit->GetTIds();
								for(unsigned int i=0; i<trks.size(); i++)
									signal.push_back(0.0);
							}
						}
						
						if(signalChoice == "Mom")
						{
							vector<G4ThreeVector> mom = aHit->GetMoms();
							for(unsigned int i=0; i<mom.size(); i++)
								signal.push_back(mom[i].mag());
						}
						if(signalChoice == "Voltage")
						{
							time   = aHit->getSignalT();
							signal = aHit->getSignalV();
							nsteps = time.size();
							
							pid.clear();
							// black pen if it's a signal
							for(unsigned int i=0; i<time.size(); i++)
								pid.push_back(1000);
						}
						
						
						SD += "  nsteps: " + stringify(nsteps);
						newHit->setText(0, QString(SD.c_str()));
						
						vector<identifier> identi = aHit->GetId();
						string title = "";
						
						for(unsigned int i=0; i<identi.size(); i++)
							title += identi[i].name + " " + stringify(identi[i].id) + "   " ;
						
						QTreeWidgetItem * signalI = new QTreeWidgetItem(newHit);
						char ttext[200];
						sprintf(ttext, "%s    pid   Time[ns] ", signalChoice.c_str());
						signalI->setText(0, QString(ttext));
						signalI->setExpanded(1);
						
						QTreeWidgetItem * signalItems;
						for(int i=0; i<nsteps; i++)
						{
							signalItems = new QTreeWidgetItem(signalI);
							char etext[200];
							sprintf(etext, "%6.5f      %d       %5.4f", signal[i], pid[i], time[i]);
							signalItems->setText(0, QString(etext));
							signalItems->setTextAlignment(1, Qt::AlignJustify);
							signalItems->setTextAlignment(2, Qt::AlignJustify);
						}
						graphView->plots_bg("time [ns]", signalChoice.c_str(), time, signal, title);
						graphView->plot_graph(time, signal, pid);
						
					}
				}
		}
	}
	
	return gsignals;
}



gsignal::~gsignal()
{
	string hd_msg = gemcOpt->optMap["LOG_MSG"].args ;
	double VERB   = gemcOpt->optMap["GEO_VERBOSITY"].arg ;
	if(VERB>2)
		cout << hd_msg << " Signal Widget Deleted." << endl;
	
}


















