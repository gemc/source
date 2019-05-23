#ifndef GammaNuclearPhysics_h
#define GammaNuclearPhysics_h 1

#include "globals.hh"
#include "G4VPhysicsConstructor.hh"


class GammaNuclearPhysics : public G4VPhysicsConstructor
{
  public:
    GammaNuclearPhysics(const G4String& name="GammaNuclearPhysics");
   ~GammaNuclearPhysics();

  public:
    virtual void ConstructParticle() { };
    virtual void ConstructProcess();
};

#endif

