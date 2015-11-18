#ifndef graph_H
#define graph_H 1

// Qt headers
#include <QtWidgets>

// C++ headers
#include <string>
using namespace std;


// Class definition
class graph : public QGraphicsView
{
	// metaobject required for non-qt slots
	Q_OBJECT
	
	public:
		graph(QWidget *parent = 0);
	  ~graph(){;}
	
		QGraphicsScene *scene;
	
		int xorig, yorig;       // origin of the axis
		int xaxil, yaxil;       // axis length
		int nticksx, nticksy;
	
		double xmin, ymin;      // graph minima
		double xmax, ymax;      // graph maxima
		double dx, dy;          // data deltas
		int fixedAxis;          // if 1, limits are provided by user
	
		double inside;          // how much inside the ticks line will be
		double DX, DY;          // graph deltas

		double axisPenWidth, dataPenWidth;
	
		void setAxis(int a, int b, int c, int d, int e, int f){xorig = a; yorig = b; xaxil = c; yaxil = d; nticksx = e; nticksy=f;}
		void setInside(double a, double b, double c){inside = a; DX = xaxil-a*b; DY = yaxil-a*c;}
		void setDataAxisLimits(vector<double>, vector<double>);
		void setPensWidth(int a, int b){axisPenWidth = a; dataPenWidth = b;}
		void setFixedAxis(int xmi, int xma, int ymi, int yma);
	
		map<int, QPen> pcolors;
	
		void plots_bg(string xtit, string ytit, vector<double> x, vector<double> y, string title);  // draw axis, ticks and labels
		void plot_graph(vector<double> x, vector<double> y, vector<int> pid);
		void plotLabel(string, int, double, double, int);
	
};


#endif
