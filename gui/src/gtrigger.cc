// Qt headers
#include <QtWidgets>

// gemc headers
#include "gtrigger.h"

gtrigger::gtrigger(QWidget *parent, goptions *Opts, map<string, sensitiveDetector*> SD_Map) : QWidget(parent)
{
	gemcOpt  = Opts;
	SeDe_Map = SD_Map;
	VOLTAGES = replaceCharWithChars(gemcOpt->optMap["SIGNALVT"].args, ",", "  ");
	plotChoice = 0;  // start with both plots
	
	// top: choose what to show
	QComboBox *whatToShow = new QComboBox;
	whatToShow->addItem("Voltages only");
	whatToShow->addItem("Triggers only");
	whatToShow->addItem("Voltages and Triggers");
	whatToShow->addItem("Summary");
	connect ( whatToShow   , SIGNAL( currentIndexChanged (int) ), this, SLOT( choosePlots(int)    ) );

	
	// need an overal scroll area with all the signals
	scrollArea = new QScrollArea;
	scrollArea->setMinimumSize(545, 565);

	vGraphsplitter  = new QSplitter(this);
	vGraphsplitter->setOrientation(Qt::Vertical);
	vGraphsplitter->setHandleWidth(0);
	scrollArea->setWidget(vGraphsplitter);

	
	// Vertical Splitter - Top and Bottom layouts
	QSplitter *vMainsplitter     = new QSplitter(this);
	vMainsplitter->setOrientation(Qt::Vertical);
	vMainsplitter->addWidget(whatToShow);
	vMainsplitter->addWidget(scrollArea);
	vMainsplitter->setMinimumSize(545, 595);
}


int gtrigger::determineNhits()
{
	// clear all widgets
	for(int i=0; i<vGraphsplitter->count(); i++)
		vGraphsplitter->widget(i)->deleteLater() ;
	

	// determine how many total hits
	unsigned nHits = 0;
	for(map<string, sensitiveDetector*>::iterator it = SeDe_Map.begin(); it!= SeDe_Map.end(); it++)
	{
		MHitCollection *MHC = it->second->GetMHitCollection();
		if(MHC)
			if(VOLTAGES.find(it->first) != string::npos)
			{
				int nhits = MHC->GetSize();
				for(int h=0; h<nhits; h++)
				{
					MHit *aHit = (*MHC)[h];
					if(aHit->diditpassTrigger())
						nHits++;
				}
			}
	}
	return nHits;
}

void gtrigger::createGraphs()
{
	int graphVsize   = 150;
	int triggerValue = 3500;
	
	unsigned ntotHits = determineNhits();
	
	// now creating all the graphs
	vGraphsplitter->setMinimumSize(540, ntotHits*graphVsize);
	vector<graph*> allGraphs(ntotHits);

	for(unsigned h=0; h<ntotHits; h++)
	{
		allGraphs[h] = new graph();
		
		// Graph axis origins and legth, number of ticks each axis
		//                 xorig  yorig  xlength ylength nticksx nticksy
		allGraphs[h]->setAxis(25,   graphVsize-20,    435,   graphVsize-30,     0,    3);
		// inside shift of the axis ticks, and factor that multiplies the ticks size
		allGraphs[h]->setInside(5, 2, 2);
		allGraphs[h]->setPensWidth(2, 1);
		allGraphs[h]->setMinimumSize(545, graphVsize);
		allGraphs[h]->setMaximumSize(545, graphVsize);

		vGraphsplitter->addWidget(allGraphs[h]);
		vGraphsplitter->setCollapsible(h,0);
	}
	

	// now plotting all hits
	// notice detector index
	int dIndex = 0;
	for(map<string, sensitiveDetector*>::iterator it = SeDe_Map.begin(); it!= SeDe_Map.end(); it++)
	{
		MHitCollection *MHC = it->second->GetMHitCollection();
		int nhits = 0;
		int ntriggers = 0;
		if(MHC) nhits = MHC->GetSize();
		if(VOLTAGES.find(it->first) != string::npos)
		{
			for(int h=0; h<nhits; h++)
			{
				MHit *aHit = (*MHC)[h];
				if(aHit->diditpassTrigger() == 0)
					continue;
				
				ntriggers++;
				
				string title = it->first + " ";
				
				vector<identifier> identi = aHit->GetId();
				for(unsigned int i=0; i<identi.size(); i++)
					title += identi[i].name + " " + stringify(identi[i].id) + "   " ;
				
				vector<int> time = aHit->getQuantumT();
				vector<int> qadc = aHit->getQuantumQ();
				vector<int> trig = aHit->getQuantumTR();
				
				// pid needed to differentiate colors
				vector<int> pid;
				vector<int> tpid;
				
				vector<double> dtime;
				vector<double> dqadc;
				vector<double> dtrig;
				vector<double> dttime;
				
				int goingUp   = 0;
				int goingDown = triggerValue;
				
				for(unsigned int i=0; i<time.size(); i++)
				{
					dtime.push_back( (double) time[i]);
					dttime.push_back((double) time[i]);
					dqadc.push_back((double) qadc[i]);
					dtrig.push_back((double) trig[i]);
					pid.push_back(2212);
					tpid.push_back(1000);
					
					// trick to add vertical transition lines between 0 and trigger
					if(trig[i] > 0)
					{
						while(goingUp < triggerValue)
						{
							goingUp += triggerValue / 100;
							dttime.push_back((double) time[i] + 0.001);
							dtrig.push_back(goingUp);
							tpid.push_back(1000);
						}
						goingDown = triggerValue;
					}
					else
					{
						if(goingUp >= triggerValue)
						{
							while(goingDown > 0)
							{
								goingDown -= triggerValue / 100;
								dttime.push_back((double) time[i] + 0.001);
								dtrig.push_back(goingDown);
								tpid.push_back(1000);
								goingUp = 0;
							}
						}
					}
				}
				
				int gIndex = dIndex + ntriggers - 1;
				allGraphs[gIndex]->setFixedAxis(0, 500, 0, 4000);
				allGraphs[gIndex]->plots_bg("time [ns]", "", dtime, dqadc, title);
				if(plotChoice == 2 || plotChoice == 0) allGraphs[gIndex]->plot_graph(dtime, dqadc, pid);
				if(plotChoice == 2 || plotChoice == 1) allGraphs[gIndex]->plot_graph(dttime, dtrig, tpid);

			}
		
			dIndex += ntriggers ;
		}
	}
}

void gtrigger::createSummary()
{
	int triggerValue = 40;
	
	unsigned ntotHits = determineNhits();
	
	
	graph* summaryGraph = new graph();
	// Graph axis origins and legth, number of ticks each axis
	//                 xorig  yorig  xlength ylength nticksx nticksy
	summaryGraph->setAxis(25,   475,    435,   450,     5,    5);
	// inside shift of the axis ticks, and factor that multiplies the ticks size
	summaryGraph->setInside(5, 2, 2);
	summaryGraph->setPensWidth(2, 2);
	summaryGraph->setMinimumSize(545, 560);
	summaryGraph->setMaximumSize(545, 560);
	
	
	vGraphsplitter->setMinimumSize(540, 570);
	vGraphsplitter->addWidget(summaryGraph);
	vGraphsplitter->addWidget(summaryGraph);


	int ntriggers = 0;
	// graph colors, index.
	int gcolors[8] = {22, 11, 211, 2212, -11, 13, -13, -211};
	int igcolor = 0;
	double xTitlePos = 350;  // position of detector labels
	double yTitlePos = 500;
	double dyTitlePos = 40;
	for(map<string, sensitiveDetector*>::iterator it = SeDe_Map.begin(); it!= SeDe_Map.end(); it++)
	{
		MHitCollection *MHC = it->second->GetMHitCollection();
		int nhits = 0;
		if(MHC) nhits = MHC->GetSize();
		if(VOLTAGES.find(it->first) != string::npos)
		{
			for(int h=0; h<nhits; h++)
			{
				MHit *aHit = (*MHC)[h];
				if(aHit->diditpassTrigger() == 0)
					continue;
				
				ntriggers++;
				
				
				// title set for debugging purposes only
				string title = it->first + " ";
				
				vector<identifier> identi = aHit->GetId();
				for(unsigned int i=0; i<identi.size(); i++)
					title += identi[i].name + " " + stringify(identi[i].id) + "   " ;
				
				vector<int> time = aHit->getQuantumT();
				vector<int> trig = aHit->getQuantumTR();
				
				// pid needed to differentiate colors
				vector<int> pid;
				
				vector<double> dtime;
				vector<double> dtrig;
				
			
				for(unsigned int i=0; i<time.size(); i++)
				{
					dtime.push_back((double) time[i]);
					
					if(trig[i])
						dtrig.push_back(triggerValue*ntriggers);
					else
						dtrig.push_back(-10000);
					
					pid.push_back(gcolors[igcolor]);
					
				}
				
				
				summaryGraph->setFixedAxis(0, 500, 0, ntotHits*triggerValue);
				if(ntriggers == 1)
				{
					summaryGraph->plots_bg("time [ns]", "", dtime, dtrig, "Triggers SUmmary");
				}
				summaryGraph->plot_graph(dtime, dtrig, pid);
				
			}
			// labels at the top
			yTitlePos = 400 - igcolor*dyTitlePos;
			summaryGraph->plotLabel(it->first, gcolors[igcolor], xTitlePos, yTitlePos, 32);

			igcolor++;
			if(igcolor == 8) igcolor = 0;

			
		}
	}

	
}

gtrigger::~gtrigger()
{
	string hd_msg = gemcOpt->optMap["LOG_MSG"].args ;
	double VERB   = gemcOpt->optMap["GEO_VERBOSITY"].arg ;
	if(VERB>2)
		cout << hd_msg << " Signal Widget Deleted." << endl;
	
}

void gtrigger::choosePlots(int index)
{
	plotChoice = index;
	if(plotChoice != 3) createGraphs();
	else createSummary();
}


