#ifndef about_H
#define about_H 1

/// \file about.h
/// about.h is a documentation file that reports the version of GEMC and
/// the libraries: CLHEP, GEANT4, Qt4, EVIO, XERCESC.
/// \author \n Maurizio Ungaro
/// \author mail: ungaro@jlab.org\n\n\n

// Notice:
// For images, use png format.

/// \class AboutTab
/// <b> AboutTab </b>\n\n
/// This is the a documentation class that provides library versions:\n
/// - GEMC
/// - CLHEP
/// - Qt4
/// - GEANT4
/// - EVIO
/// - XERCESC

extern char* GEMC_VERSION;

QString AboutTab::load_doc()
{
	string doc = "<html><body><center>";
	doc.append("<br>");
	doc.append(GEMC_VERSION);
	doc.append("<br>");
	doc.append("<br>");
	doc.append("<h3> Author, Support: </h3>");
	doc.append("Maurizio Ungaro ");
	doc.append("<a target='_blank' href='mailto:ungaro@jlab.org'>ungaro@jlab.org</a>");
	doc.append("<br><br><h3> Website: </h3>");
	doc.append("<a target='_blank' href=\"http://gemc.jlab.org\">gemc.jlab.org</a><br>");
	doc.append("</body></center></html>");
	
	return doc.c_str();
}
#endif
