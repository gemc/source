// gemc headers
#include "particleFactory.h"

// geant4 headers
#include "G4ParticleTable.hh"
#include "Randomize.hh"


gparticle::gparticle(G4ThreeVector vrt, G4ThreeVector mom, int id, G4ThreeVector polar)
{
	v   = vrt;
	p   = mom;
	pid = id;
	
	G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
	G4ParticleDefinition *pdef = particleTable->FindParticle(id);
	
	name = pdef->GetParticleName();
	
	// setting polarization
	double polDeg   = polar.x();
	double polTheta = polar.y();
	double polPhi   = polar.z();
	
	double partPol = 0.0;
	double polCast = 100.0 * G4UniformRand();
	if( polCast <= polDeg ) partPol = 1;
	double polX = partPol * sin( polTheta/rad ) * cos( polPhi/rad );
	double polY = partPol * sin( polTheta/rad ) * sin( polPhi/rad );
	double polZ = partPol * cos( polTheta/rad );
	polarization = G4ThreeVector( polX, polY, polZ );
	
}




///< Overloaded "<<" for gparticle class. Dumps infos on screen.
ostream &operator<<(ostream &stream, gparticle p)
{
	cout << "  - Particle <" << p.name << ">  (id " << p.pid << ")"
	<< " momentum: "     << p.p
	<< " vertex: "       << p.v
	<< " polarization: " << p.polarization
	<< " start time: "   << p.vtime << endl;
	
	return stream;
}



