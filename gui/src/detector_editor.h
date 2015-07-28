/// \file detector_editor.h
/// Defines the detector_editor class.\n
/// detector_editor is a QDialog formed by
/// Placement, Dimensions, Sensitivity, Magnetic Field Tabs
/// \author \n Maurizio Ungaro
/// \author mail: ungaro@jlab.org\n\n\n

#ifndef DETECTOR_EDITOR_H
#define DETECTOR_EDITOR_H


// Qt headers
//#include <QDialog>
//#include <QDialogButtonBox>
#include <QLabel>
#include <QLineEdit>
//#include <QTabWidget>
#include <QRadioButton>

// gemc headers
#include "detector.h"


/// \class descriptionTab
/// <b> descriptionTab </b>\n\n
/// This tab contains the description
/// of the volume, including name, mother,
/// position, rotation, type, dimensions,
/// materials, etc
class descriptionTab : public QWidget
{
	Q_OBJECT
	
	public:
		
		descriptionTab(QWidget *parent = 0);  ///< Constructor
		
		detector *det;           ///< pointer to detector
		
		QLabel*    nameLabel;    ///< label for name
		QLabel*    descLabel;    ///< label for description

		QLabel *solidType;       ///< label for solid type
		QPixmap solidpic;        ///< solid pic from G4
		QLabel *solidPicL;       ///< container for solid pic

		QLabel*    placeLabel;   ///< label for position
		QLineEdit* placeEdit;    ///< line editor for position
		
		QLabel*    rotLabel;     ///< label for rotation
		QLineEdit* rotEdit;      ///< line editor for rotation
	
		QLabel*    dimLabel;     ///< label for dimensions
		QLineEdit* dimEdit;      ///< line editor for dimensions
	
		QLabel*    matLabel;     ///< label for material
		QLabel*    mgnLabel;     ///< label for magnetic field
		QLabel*    sensHitLabel; ///< label for sensitivity, Hit process
		QLabel*    idLabel;      ///< label for identifier
	
		QLabel*    twLabel;      ///< label for time window
		QLabel*    pcutLabel;    ///< label for production cut
		
	
	private slots:
		void change_placement();                ///< changes coordinates/rotation of the detector
		void change_dimension();              ///< changes the detector dimensions

	public slots:
		void update_detector(detector *Det);    ///< change tab according to detector
	
	
};



#endif
