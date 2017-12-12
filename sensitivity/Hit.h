/// \file Hit.h
/// Defines the MHit class.\n
/// \author \n Maurizio Ungaro
/// \author mail: ungaro@jlab.org\n\n\n

#ifndef MHit_H
#define MHit_H 1

// G4 headers
#include "G4ThreeVector.hh"
#include "G4VHit.hh"

// gemc headers
#include "detector.h"
#include "sensitiveID.h"

// C++ headers
#include <iostream>
using namespace std;


// TODO: Many of these quantities could be calculated as needed if the Touchable history is given.


// Class definition
class MHit : public G4VHit
{

public:
	MHit();
	
	// electronic noise hit
	MHit(double edep, double time, vector<identifier> identity, int pid);
	// background hit
	MHit(double edep, double time, int nphe, vector<identifier> identity);

	virtual ~MHit();
	const MHit& operator=(const MHit&){return *this;}

	void Draw();

	G4Colour colour_touch, colour_hit, colour_passby;

private:
	// all these infos are recorded
	// in each step of the hit
	vector<G4ThreeVector>   pos;    ///< Hit Position (Global)
	vector<G4ThreeVector>  Lpos;    ///< Hit Positions (Local to the Volume)
	vector<G4ThreeVector>  vert;    ///< Primary Vertex of track in each step
	vector<double>         edep;    ///< Energy Deposited
	vector<double>           dx;    ///< Length of the step
	vector<double>         time;    ///< Time from the start of event
	vector<G4ThreeVector>   mom;    ///< Momentum of the Track
	vector<double>            E;    ///< Energy of the track
	vector<int>               q;    ///< Charge of the particle in each step
	vector<int>             PID;    ///< particle ID in each step
	vector<int>            mPID;    ///< Mother particle ID in each step
	vector<int>         trackID;    ///< G4Track ID in each step
	vector<int>        mtrackID;    ///< Mother G4Track ID in each step
	vector<int>        otrackID;    ///< Original G4Track ID in each step
	vector<G4ThreeVector> mvert;    ///< Primary Vertex of the track's mother
	vector<string> materialName;    ///< Material name
	vector<int>       processID;    ///< Process that originated this step
    vector<double>         mgnf;    ///< magnetic field

	vector<detector>  Detectors;    ///< Detectors Hit. It might be a vector if multiple detectors have the same identifier

	vector<identifier> identity;    ///< Identity
	sensitiveID SID;                ///< Sensitive ID has detector information like  signalThreshold, timeWindow, prodThreshold, maxStep, riseTime, fallTime, mvToMeV

	vector<double> signalT;         ///< Time vector
	vector<double> signalV;         ///< Voltage Vector

	vector<int> quantumT;           ///< quantized Time vector
	vector<int> quantumQ;           ///< quantized charge (ADC) Vector
	vector<int> quantumTR;          ///< Trigger based on QuantumQ

	int hasTrigger;                 ///< is 1 if this hit produces a signal above threshold

public:
	// infos filled in Sensitive Detector
	inline void SetPos(G4ThreeVector xyz)       { pos.push_back(xyz); }
	inline vector<G4ThreeVector> GetPos()       { return pos; }
	inline G4ThreeVector GetLastPos()           { if(pos.size()) return pos[pos.size()-1]; else return G4ThreeVector(0,0,0); }

	inline void SetLPos(G4ThreeVector xyz)      { Lpos.push_back(xyz); }
	inline vector<G4ThreeVector> GetLPos()      { return Lpos; }

	inline void SetVert(G4ThreeVector ver)      { vert.push_back(ver); }
	inline G4ThreeVector GetVert()              { return  vert[0]; }
	inline vector<G4ThreeVector> GetVerts()     { return  vert; }

	inline void SetEdep(double depe)            { edep.push_back(depe); }
	inline vector<double> GetEdep()             { return edep; }

	inline void SetDx(double Dx)                { dx.push_back(Dx); }
	inline vector<double> GetDx()               { return dx; }

	inline void SetMgnf(double m)               { mgnf.push_back(m); }
	inline vector<double> GetMgnf()             { return  mgnf; }

    inline void SetTime(double ctime)           { time.push_back(ctime); }
    inline vector<double> GetTime()             { return  time; }

	inline void SetMom(G4ThreeVector pxyz)      { mom.push_back(pxyz); }
	inline G4ThreeVector GetMom()               { return mom[0]; }
	inline vector<G4ThreeVector> GetMoms()      { return mom; }

	inline void SetE(double ene)                { E.push_back(ene); }
	inline double GetE()                        { return E[0]; }
	inline vector<double> GetEs()               { return E; }

	inline void SetTrackId(int tid)             { trackID.push_back(tid); }
	inline int GetTId()                         { return trackID[0]; }
	inline vector<int> GetTIds()                { return trackID; }

	inline vector<identifier> GetId()           { return identity; }
	inline void SetId(vector<identifier> iden)  { identity = iden; }

	inline void SetDetector(detector det)       {Detectors.push_back(det);}
	inline vector<detector> GetDetectors()      {return Detectors;}
	inline detector GetDetector()               {return Detectors[0];}

	inline void SetPID(int pid)                 { PID.push_back(pid); }
	inline int GetPID()                         { return PID[0]; }
	inline vector<int> GetPIDs()                { return PID; }

	inline void SetCharge(int Q)                { q.push_back(Q); }
	inline int GetCharge()                      { return q[0]; }
	inline vector<int> GetCharges()             { return q; }

	// infos filled in MEvent Action
	inline void SetmTrackId(int tid)            { mtrackID.push_back(tid); }
	inline void SetmTrackIds(vector<int> tid)   { mtrackID = tid; }
	inline int GetmTrackId()                    { return mtrackID[0]; }
	inline vector<int> GetmTrackIds()           { return mtrackID; }

	inline void SetoTrackId(int tid)            { otrackID.push_back(tid); }
	inline void SetoTrackIds(vector<int> tid)   { otrackID = tid; }
	inline int GetoTrackId()                    { return otrackID[0]; }
	inline vector<int> GetoTrackIds()           { return otrackID; }

	inline void SetmPID(int mpid)               { mPID.push_back(mpid); }
	inline void SetmPIDs(vector<int> mpid)      { mPID = mpid; }
	inline int GetmPID()                        { return mPID[0]; }
	inline vector<int> GetmPIDs()               { return mPID; }

	inline void SetmVert(G4ThreeVector ver)          { mvert.push_back(ver); }
	inline void SetmVerts(vector<G4ThreeVector> ver) { mvert = ver; }
	inline G4ThreeVector GetmVert()                  { return  mvert[0]; }
	inline vector<G4ThreeVector> GetmVerts()         { return  mvert; }

	inline void SetMatName(string mname)           { materialName.push_back(mname); }
	inline void SetMatNames(vector<string> mnames) { materialName = mnames; }
	inline string GetMatName()                     { return  materialName[0]; }
	inline vector<string> GetMatNames()            { return  materialName; }

	inline void SetProcID(int procID)          { processID.push_back(procID); }
	inline void SetProcID(vector<int> procIDs) { processID = procIDs; }
	inline int GetProcID()                     { return  processID[0]; }
	inline vector<int> GetProcIDs()            { return  processID; }

	inline void SetSDID(sensitiveID s)   { SID = s; }
	inline sensitiveID GetSDID()         { return SID; }


	inline void setSignal(map< double, double > VT)
	{
		signalT.clear();
		signalV.clear();

		for(map< double, double >::iterator it = VT.begin(); it!=VT.end(); it++)
		{
			signalT.push_back(it->first);
			signalV.push_back(it->second);
		}
	}

	inline vector<double> getSignalT(){return signalT;}
	inline vector<double> getSignalV(){return signalV;}

	inline void setQuantum(map< int, int > QS)
	{
		quantumT.clear();
		quantumQ.clear();

		for(map< int, int >::iterator it = QS.begin(); it!=QS.end(); it++)
		{
			quantumT.push_back(it->first);
			quantumQ.push_back(it->second);
		}
	}

	inline vector<int> getQuantumT() {return quantumT;}
	inline vector<int> getQuantumQ() {return quantumQ;}
	inline vector<int> getQuantumTR(){return quantumTR;}
	inline void setQuantumTR(vector<int> t)   { quantumTR = t; }

	// trigger
	inline void passedTrigger(){hasTrigger = 1;}
	inline int diditpassTrigger(){return hasTrigger;}

	int isElectronicNoise;          ///< 1 if this is an electronic noise hit
	int isBackgroundHit;            ///< 1 if this hit a background hit

};


#include "G4THitsCollection.hh"
typedef G4THitsCollection<MHit> MHitCollection;

#endif












