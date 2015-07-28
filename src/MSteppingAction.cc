// G4 headers
#include "G4ParticleTypes.hh"

// gemc headers
#include "MSteppingAction.h"

MSteppingAction::MSteppingAction(goptions Opt)
{
	gemcOpt   = Opt;
	energyCut = gemcOpt.optMap["ENERGY_CUT"].arg;
	max_x_pos = gemcOpt.optMap["MAX_X_POS"].arg;
	max_y_pos = gemcOpt.optMap["MAX_Y_POS"].arg;
	max_z_pos = gemcOpt.optMap["MAX_Z_POS"].arg;
	
//	oldpos = G4ThreeVector(0,0,0);
//	nsame  = 0;
}

MSteppingAction::~MSteppingAction(){ cout << " > Closing Stepping Action." << endl;}


void MSteppingAction::UserSteppingAction(const G4Step* aStep)
{
	G4ThreeVector   pos   = aStep->GetPostStepPoint()->GetPosition();      ///< Global Coordinates of interaction
	G4Track*        track = aStep->GetTrack();
	string          volname(aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName());
	
	if(fabs(pos.x()) > max_x_pos ||
	   fabs(pos.y()) > max_y_pos ||
	   fabs(pos.z()) > max_z_pos ) track->SetTrackStatus(fStopAndKill);   ///< Killing track if outside of interest region
	
	if(track->GetKineticEnergy() < energyCut )
		track->SetTrackStatus(fStopAndKill);
	
	// Anything passing material "Kryptonite" is killed
	if(track->GetMaterial()->GetName() == "Kryptonite")
	{
		track->SetTrackStatus(fStopAndKill);
	}
	
	
	if(track->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition())
	{
		if(track->GetLogicalVolumeAtVertex()->GetMaterial()->GetName() == "SemiMirror")
		{
			track->SetTrackStatus(fStopAndKill);
		}
	}
	
	
//	// checking if a step is stuck in the same position
//	// for more than 10 steps
//    // this should be revisited
//	if(sqrt( (pos - oldpos).x() *  (pos - oldpos).x() +
//			 (pos - oldpos).y() *  (pos - oldpos).y() +
//			 (pos - oldpos).z() *  (pos - oldpos).z() ) < 0.0000001*mm)
//		
//	{
//		nsame++;
//		if(nsame > 100)
//		{
//			cout << " Track is stuck. PID: " <<  track->GetDefinition()->GetPDGEncoding() << " Volume: " << volname ;
//			cout << " Last step at :  " << pos << " old step was at " << oldpos <<  "    track id: " << track->GetTrackID() << endl;
//			
//			cout << ". Killing this track. " << endl;
//			
//			track->SetTrackStatus(fStopAndKill);
//			
//		}
//	}
//	else
//		nsame = 0;
//	
//	oldpos = pos;

	// limiting steps in one volume to 1000
	// int nsteps = aStep->GetTrack()->GetCurrentStepNumber();
	//if(nsteps >= 1000) aStep->GetTrack()->SetTrackStatus(fStopAndKill);

}

