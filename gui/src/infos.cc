// Qt headers
#include <QtWidgets>

// gemc headers
#include "infos.h"
#include "docs/about.h"
#include "docs/particles_color.h"
#include "docs/lund.h"

// Notice:
// For images, use png format.

infos::infos(QWidget *parent, goptions *Opts) : QWidget(parent)
{
	gemcOpt = Opts;
	
	QTabWidget* InfoType = new QTabWidget;
	abouttab   = new AboutTab(this,   gemcOpt);
	pcolorstab = new PColorsTab(this, gemcOpt);
	lundtab    = new LUNDTab(this,    gemcOpt);
	
	InfoType->addTab(abouttab,   tr("About"));
	InfoType->addTab(pcolorstab, tr("Particle Colors"));
	InfoType->addTab(lundtab,    tr("LUND format"));
	
	QVBoxLayout *mainLayout = new QVBoxLayout;
	mainLayout->addWidget(InfoType);
	setLayout(mainLayout);
}


	
infos::~infos()
{
	string hd_msg = gemcOpt->optMap["LOG_MSG"].args ;
	double VERB   = gemcOpt->optMap["GEO_VERBOSITY"].arg ;
	if(VERB>2)
		cout << hd_msg << " Infos Widget Deleted." << endl;
	
}


AboutTab::AboutTab(QWidget *parent, goptions *Opts) : QWidget(parent)
{
	QTextEdit *view = new QTextEdit(parent);
    view->setHtml( load_doc() );
	
	QVBoxLayout *mLayout = new QVBoxLayout;
	mLayout->addWidget(view);
	setLayout(mLayout);
}

PColorsTab::PColorsTab(QWidget *parent, goptions *Opts) : QWidget(parent)
{
	QTextEdit *view = new QTextEdit(parent);
	view->setHtml(load_doc());
	
	QVBoxLayout *mLayout = new QVBoxLayout;
	mLayout->addWidget(view);
	setLayout(mLayout);
}

LUNDTab::LUNDTab(QWidget *parent, goptions *Opts) : QWidget(parent)
{
	QTextEdit *view = new QTextEdit(parent);
	view->setHtml(load_doc());
	
	QVBoxLayout *mLayout = new QVBoxLayout;
	mLayout->addWidget(view);
	setLayout(mLayout);
}


















