#ifndef graph_H
#define graph_H 1

// Qt headers
#include <QtWidgets>

// C++ headers
#include <string>
using namespace std;


// Class definition
class graph : public QGraphicsScene
{
	// metaobject required for non-qt slots
	Q_OBJECT
	
	public:
	
		int xorig, yorig;       // origin of the axis
		int xaxil, yaxil;       // axis length
		int nticksx, nticksy;
	
		double xmin, ymin;      // graph minima
		double xmax, ymax;      // graph maxima
		double dx, dy;          // lowercase: graph limits
	
		double inside;          // how much inside the graph will be
		double DX, DY;          // uppercase: scene inside limits

		void setAxis(int a, int b, int c, int d, int e, int f){xorig = a; yorig = b; xaxil = c; yaxil = d; nticksx = e; nticksy=f;}
		void setAxisLimits(double a, double b, double c, double d, double e, double f){xmin = a; ymin = b; xmax = c; ymax = d; dx = e; dy = f;}
		void setInside(double a, double b, double c){inside = a; DX = b; DY = c;}
	
		map<int, QPen> pcolors;

		graph(QWidget *parent);
	   ~graph(){;}

		QGraphicsScene *scene;
	
		void plots_bg(string xtit, string ytit, vector<double> x, vector<double> y, string title);  // draw axis, ticks and labels
		void plot_graph(vector<double> x, vector<double> y, vector<int> pid);

	
};


#endif
