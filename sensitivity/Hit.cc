// G4 headers
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4VisAttributes.hh"
#include "G4ParticleTable.hh"

// gemc headers
#include "Hit.h"

MHit::MHit() : G4VHit()
{
	// If the energy is above threshold, red circle.
	// If the energy is below threshold, blue smaller circle.
	// If no energy deposited, green smaller circle.

	colour_touch  = G4Colour(0.0, 0.0, 1.0);
	colour_hit    = G4Colour(1.0, 0.0, 0.0);
	colour_passby = G4Colour(0.0, 1.0, 0.0);

	hasTrigger = 0;
	isElectronicNoise = 0;
	isBackgroundHit = 0;

}

MHit::~MHit() {;}

void MHit::Draw()
{
	G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
	if(pVVisManager)
	{
		G4Circle circle(pos[0]);
		circle.SetFillStyle(G4Circle::filled);


		// summing all energies
		double Etot = 0;
		for(unsigned int i=0; i<edep.size(); i++)
			Etot += edep[i];

		if(Etot > SID.signalThreshold)
		{
			circle.SetVisAttributes(G4VisAttributes(colour_hit));
			circle.SetScreenSize(5);
		}
		else if(Etot > 0 && Etot <= SID.signalThreshold)
		{
			circle.SetVisAttributes(G4VisAttributes(colour_touch));
			circle.SetScreenSize(4);
		}
		else if(Etot == 0)
		{
			circle.SetVisAttributes(G4VisAttributes(colour_passby));
			circle.SetScreenSize(3);
		}

		if(PID[0] != 0) {
			pVVisManager->Draw(circle);
		}
	}
}

MHit::MHit(double energy, double tim, vector<identifier> vid, int pid)
{
	isElectronicNoise = 1;

	pos.push_back(G4ThreeVector(0,0,0));
	Lpos.push_back(G4ThreeVector(0,0,0));
	vert.push_back(G4ThreeVector(0,0,0));
	edep.push_back(energy);
	dx.push_back(0);
	time.push_back(tim);
	mom.push_back(G4ThreeVector(0,0,0));
	E.push_back(0);
	q.push_back(0);
	PID.push_back(pid);
	mPID.push_back(0);
	trackID.push_back(-1);
	mtrackID.push_back(-1);
	otrackID.push_back(-1);
	mvert.push_back(G4ThreeVector(0,0,0));
	materialName.push_back("noise");
	processID.push_back(999);
	mgnf.push_back(0);

	identity = vid;

}


MHit::MHit(double energy, double tim, int nphe, vector<identifier> vid)
{
	isBackgroundHit = 1;

	pos.push_back(G4ThreeVector(0,0,0));
	Lpos.push_back(G4ThreeVector(0,0,0));
	vert.push_back(G4ThreeVector(0,0,0));
	edep.push_back(energy);
	dx.push_back(0);
	time.push_back(tim);
	mom.push_back(G4ThreeVector(0,0,0));
	E.push_back(0);
	q.push_back(nphe);
	PID.push_back(0);
	mPID.push_back(0);
	trackID.push_back(-1);
	mtrackID.push_back(-1);
	otrackID.push_back(-1);
	mvert.push_back(G4ThreeVector(0,0,0));
	materialName.push_back("backgroundHit");
	processID.push_back(999);
	mgnf.push_back(0);

	identity = vid;

}








