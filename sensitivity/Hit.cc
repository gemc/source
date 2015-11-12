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
		
		pVVisManager->Draw(circle);
	}
}

