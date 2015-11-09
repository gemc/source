#include "graph.h"

// C++ headers
#include <iostream>



graph::graph(QWidget *parent)
{
	scene = new QGraphicsScene(parent);
}



void graph::plots_bg(string xtit, string ytit, vector<double> x, vector<double> y, string title)
{
	
	scene->clear();
	// title
	QGraphicsSimpleTextItem *Title = new QGraphicsSimpleTextItem(QString(title.c_str()));
	scene->addItem(Title);
	QFont sansFont("Times-Roman", 18);
	Title->setFont(sansFont);
	Title->moveBy( -(double) title.length()*18.0, 10); // this should more or less center it
	
	if(x.size()<1) return;
	
	// axis
	QGraphicsLineItem *xaxis  = new QGraphicsLineItem( xorig-60      ,  yorig, xorig+xaxil-60     ,  yorig);
	//	QGraphicsLineItem *xaxisa = new QGraphicsLineItem( xorig+xaxil,  yorig, xorig+xaxil - 15,  yorig + 4);
	//	QGraphicsLineItem *xaxisb = new QGraphicsLineItem( xorig+xaxil,  yorig, xorig+xaxil - 15,  yorig - 4);
	xaxis->setPen( QPen(Qt::black, 3));
	//	xaxisa->setPen(QPen(Qt::black, 2));
	//	xaxisb->setPen(QPen(Qt::black, 2));
	scene->addItem(xaxis);
	//	scene->addItem(xaxisa);
	//	scene->addItem(xaxisb);
	
	
	QGraphicsLineItem *yaxis  = new QGraphicsLineItem(xorig-60, yorig+4, xorig-60,  yorig-yaxil);
	//	QGraphicsLineItem *yaxisa = new QGraphicsLineItem(xorig,  yorig-yaxil, xorig - 4,  yorig-yaxil + 15);
	//	QGraphicsLineItem *yaxisb = new QGraphicsLineItem(xorig,  yorig-yaxil, xorig + 4,  yorig-yaxil + 15);
	yaxis->setPen( QPen(Qt::black, 3));
	//	yaxisa->setPen(QPen(Qt::black, 3));
	//	yaxisb->setPen(QPen(Qt::black, 3));
	scene->addItem(yaxis);
	//	scene->addItem(yaxisa);
	//	scene->addItem(yaxisb);
	
	
	// labels
	QGraphicsSimpleTextItem *xlab = new QGraphicsSimpleTextItem(QString(xtit.c_str()));
	QGraphicsSimpleTextItem *ylab = new QGraphicsSimpleTextItem(QString(ytit.c_str()));
	scene->addItem(xlab);
	scene->addItem(ylab);
	xlab->moveBy(xorig+xaxil-60, yorig + 10);
	ylab->moveBy(xorig-30, yorig-yaxil+50);
	ylab->setRotation(-90);
	
	// calculating minima and maxima
	xmin = 100000;
	ymin = 100000;
	for(unsigned int i=0; i<x.size(); i++) if(x[i] < xmin) xmin = x[i];
	for(unsigned int i=0; i<y.size(); i++) if(y[i] < ymin) ymin = y[i];
	xmax = -100000;
	ymax = -100000;
	for(unsigned int i=0; i<x.size(); i++) if(x[i] > xmax) xmax = x[i];
	for(unsigned int i=0; i<y.size(); i++) if(y[i] > ymax) ymax = y[i];
	
	if(ymin == ymax)
	{
		if(ymin < 0) ymin *= 1.0001;
		else ymin *= 0.9999;
		
		if(ymax > 0) ymax *= 1.0001;
		else ymax *= 0.9999;
	}
	if(ymin == 0 && ymax == 0)
	{
		ymin = -1;
		ymax = 1;
	}
	
	// putting 5 vertical axis - dashed lines - between xmin and xmax
	dx = (xmax - xmin);
	dy = (ymax - ymin);
	
	// lab should be at least size 10
	char lab[12];
	QGraphicsSimpleTextItem *alab;
	
	QGraphicsLineItem *xtick;
	for(int a=0; a<nticksx+1; a++)
	{
		xtick = new QGraphicsLineItem(xorig + inside + a*DX/nticksx, yorig + 4, xorig + inside + a*DX/nticksx,  yorig-yaxil+inside);
		xtick->setPen( QPen(Qt::blue, 1, Qt::DashDotLine));
		scene->addItem(xtick);
	}
	
	for(int a=0; a<nticksx; a++)
	{
		sprintf(lab, "%4.3f", xmin+ a*dx/nticksx);
		alab = new QGraphicsSimpleTextItem(QString(lab));
		QFont sansFont("Helvetica", 12);
		alab->setFont(sansFont);
		alab->moveBy(xorig + inside + a*DX/nticksx - 12, yorig + 8);
		scene->addItem(alab);
	}
	
	QGraphicsLineItem *ytick;
	for(int a=0; a<nticksy+1; a++)
	{
		ytick = new QGraphicsLineItem(xorig - 4, yorig - inside - a*DY/nticksy, xorig + xaxil - inside,  yorig - inside - a*DY/nticksy);
		ytick->setPen( QPen(Qt::blue, 1, Qt::DashDotLine));
		scene->addItem(ytick);
	}
	
	
	for(int a=0; a<nticksy; a++)
	{
		if     (fabs(ymin + a*dy/nticksy) < 0.01)
			sprintf(lab, "%5.4f", ymin+ a*dy/nticksy);
		else if(fabs(ymin + a*dy/nticksy) < 1)
			sprintf(lab, "%5.3f", ymin+ a*dy/nticksy);
		else if(fabs(ymin + a*dy/nticksy) < 10)
			sprintf(lab, "%5.2f", ymin+ a*dy/nticksy);
		else if(fabs(ymin + a*dy/nticksy) < 100)
			sprintf(lab, "%5.1f", ymin+ a*dy/nticksy);
		else if(fabs(ymin + a*dy/nticksy) > 1000)
			sprintf(lab, "%2.1e", ymin+ a*dy/nticksy);
		else
			sprintf(lab, "%5.0f", ymin+ a*dy/nticksy);
		
		alab = new QGraphicsSimpleTextItem(QString(lab));
		QFont sansFont("Helvetica", 12);
		alab->setFont(sansFont);
		alab->moveBy(xorig - 80, yorig - inside - a*DY/nticksy - 6);
		scene->addItem(alab);
	}
	
	
}

void graph::plot_graph(vector<double> x, vector<double> y, vector<int> pid)
{
	if(x.size()<2) return;
	
	QGraphicsRectItem    *rect;
	
	for(unsigned int i=0; i<x.size(); i++)
	{
		if(pcolors.find(pid[i]) == pcolors.end())
			cout << " Attention: color not found for: " << pid[i] << endl;
		else
		{
			rect = scene->addRect(0, 0, 4, 4, pcolors[pid[i]]);
			rect->moveBy(xorig + inside + DX*(x[i]-xmin)/dx, yorig - inside - DY*(y[i]-ymin)/dy);
		}
	}
}
