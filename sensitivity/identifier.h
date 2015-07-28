/// \file identifier.h
/// Defines the identifier class.\n
/// \author \n Maurizio Ungaro
/// \author mail: ungaro@jlab.org\n\n\n

#ifndef identifier_H
#define identifier_H 1

// G4 headers
#include "G4VTouchable.hh"

// C++ headers
#include <string>
#include <vector>
using namespace std;


/// \class identifier
/// <b> Identifier</b>\n\n
/// - name: identifer name. Must match volume name if <i> rule = ncopy </i>.
/// - rule:
///    -# ncopy: the copy number of the volume will be used for identification.
///    -# manual: id be used for identification.
/// - id: manually assigned id.
class identifier
{
	public:
		identifier()
		{
			time       = 0;
			TimeWindow = 0;
			id_sharing = 1;
		}
		
		~identifier(){;}
		
		string       name;   ///< Name of the detector
		string       rule;   ///< "manual" or "ncopy"
		int            id;   ///< manually assing ID. 0 if "ncopy" (will be set at hit processing time)
		double       time;   ///< Time of the first step
		double TimeWindow;   ///< Time Window. If abs(steptime - time) is smaller than TimeWindow, it's the same hit
		int       TrackId;   ///< If Time Window is 0, it's a "flux" detector: if it's the same track, it's the same hit. Different track, different hit.
		double id_sharing;   ///< A single step can generate multiple identifiers. This variable represent the percentage sharing of the current identifier. Sum must be normalized to 1
		
	public:
		friend ostream &operator<<(ostream &stream, vector<identifier>);       ///< Overloaded "<<" for the class 'identifier'
		bool operator== (const identifier& I) const;                           ///< Overloaded "==" operator for the class 'identifier'
		bool operator<  (const identifier& I) const;                           ///< Overloaded "<"  operator for the class 'identifier'
		bool operator>  (const identifier& I) const;                           ///< Overloaded ">"  operator for the class 'identifier'
		bool operator<= (const identifier& I) const {return !(*this >  I);}    ///< Overloaded "<=" operator for the class 'identifier'
		bool operator>= (const identifier& I) const {return !(*this <  I);}    ///< Overloaded ">=" operator for the class 'identifier'
		bool operator!= (const identifier& I) const {return !(*this == I);}    ///< Overloaded "!=" operator for the class 'identifier'
};


// move this somewhere?
vector<identifier> SetId(vector<identifier>, G4VTouchable*, double, double, int);  ///< Sets the ncopy ID accordingly to Geant4 Volumes copy number. Sets time, TimeWindow, TrackId

// returns vector of identifier from stringstream
vector<identifier> get_identifiers(string var);

#endif

