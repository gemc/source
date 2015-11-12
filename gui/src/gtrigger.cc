// Qt headers
#include <QtWidgets>

// gemc headers
#include "gtrigger.h"

gtrigger::gtrigger(QWidget *parent, goptions *Opts, map<string, sensitiveDetector*> SD_Map) : QWidget(parent)
{
	gemcOpt  = Opts;
	SeDe_Map = SD_Map;
	VOLTAGES = replaceCharWithChars(gemcOpt->optMap["SIGNALVT"].args, ",", "  ");

	// top: choose what to show
	QComboBox *whatToShow = new QComboBox;
	whatToShow->addItem("Voltages and Triggers");
	whatToShow->addItem("Voltages only");
	whatToShow->addItem("Triggers only");

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



void gtrigger::createGraphs()
{
	int graphVsize = 150;
	
	// clear all widgets first
	for(unsigned i=0; i<vGraphsplitter->count(); i++)
	{
		vGraphsplitter->widget(i)->hide();
		//		delete vGraphsplitter->widget(i);
		vGraphsplitter->widget(i)->deleteLater() ;
	}
	
	// first determine how many total hits
	unsigned ntotHits = 0;
	for(map<string, sensitiveDetector*>::iterator it = SeDe_Map.begin(); it!= SeDe_Map.end(); it++)
	{
		MHitCollection *MHC = it->second->GetMHitCollection();
		if(MHC)
			if(VOLTAGES.find(it->first) != string::npos)
				ntotHits += MHC->GetSize();
	}
	
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
		if(MHC) nhits = MHC->GetSize();
		if(VOLTAGES.find(it->first) != string::npos)
		{
			for(int h=0; h<nhits; h++)
			{
				MHit *aHit = (*MHC)[h];
				
				string title = it->first + " ";
				
				vector<identifier> identi = aHit->GetId();
				for(unsigned int i=0; i<identi.size(); i++)
					title += identi[i].name + " " + stringify(identi[i].id) + "   " ;
				
				vector<int> time = aHit->getQuantumT();
				vector<int> qadc = aHit->getQuantumQ();
				// pid needed to differentiate colors
				vector<int> pid;
				
				vector<double> dtime;
				vector<double> dqadc;
				for(unsigned int i=0; i<time.size(); i++)
				{
					dtime.push_back((double) time[i]);
					dqadc.push_back((double) qadc[i]);
					pid.push_back(2212);
					
				}
				
				
				allGraphs[dIndex + h]->plots_bg("bunch [4*ns]", "Voltage", dtime, dqadc, title);
				allGraphs[dIndex + h]->plot_graph(dtime, dqadc, pid);
				
			}
		
		
			dIndex += nhits ;
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



