#ifndef camera_control_H
#define camera_control_H 1

// Qt headers
#include <QtWidgets>

// gemc headers
#include "options.h"

// G4 headers
#include "G4UImanager.hh"

// C++ headers
#include <string>
#include <map>
using namespace std;


// Class definition
class camera_control : public QWidget
{
	// metaobject required for non-qt slots
	Q_OBJECT
	
	public:
		camera_control(QWidget *parent, goptions*);
	 ~camera_control();
		
		goptions *gemcOpt;
		G4UImanager  *UImanager;

	private:
		int theta_hall, phi_hall;    ///< theta, phi of the camera
		
		QSlider *theta_slider;       ///< Theta Slider
		QSlider *phi_slider;         ///< Phi Slider
		vector<string> ThetaSet;     ///< strings to ThetaCombo LCD display
		vector<string> PhiSet;       ///< strings to PhiCombo LCD display
		
		QSlider  *explodeSlider;     ///< Pan Factor Slider

		QComboBox *moveCombo;
		QComboBox *detvCombo;
		QComboBox *projCombo;
	
		QComboBox *aliasing;
		QComboBox *sides_per_circle;
		QComboBox *auxiliary;

	
		// slices
		QLineEdit *sliceXEdit, *sliceYEdit, *sliceZEdit;
		QCheckBox *sliceXActi, *sliceYActi, *sliceZActi;
		QCheckBox *sliceXInve, *sliceYInve, *sliceZInve;
	
	private slots:
		void change_theta(int);
		void change_theta_s(int);
		void set_theta(int);
		void change_phi(int);
		void change_phi_s(int);
		void set_phi(int);
		void update_angles();
		
		void slice();
		void clearSlice();

		void explode(int);
		
		void printPNG();
		void printEPS();
		void printPDF();
	
		void set_perspective(int);
		void switch_antialiasing(int);
		void switch_auxiliary_edges(int);
		void switch_sides_per_circle(int);
		
};

#endif








