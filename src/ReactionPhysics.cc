#include "ReactionPhysics.hh"

#include "ReactionProcess.hh"
#include "G4GenericIon.hh"
#include "G4Alpha.hh"

#include "globals.hh"
#include "G4PhysicsListHelper.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
namespace MdmPpacSim
{
    G4_DECLARE_PHYSCONSTR_FACTORY(ReactionPhysics);

    ReactionPhysics::ReactionPhysics(G4int)
        : G4VPhysicsConstructor("ReactionPhysics")
    {
    }

    ReactionPhysics::ReactionPhysics(const G4String &name)
        : G4VPhysicsConstructor(name)
    {
    }

    ReactionPhysics::~ReactionPhysics()
    {
    }

    void ReactionPhysics::ConstructParticle()
    {
        G4GenericIon::GenericIon();
        G4Alpha::Alpha();
    }

    void ReactionPhysics::ConstructProcess()
    {
        ReactionProcess *reactionProcess = new ReactionProcess();
        reactionProcess->SetParams(fReactionParams);
        reactionProcess->SetWaveFunction(fWaveFunction);
        G4PhysicsListHelper::GetPhysicsListHelper()->RegisterProcess(reactionProcess, G4GenericIon::GenericIon());
        G4PhysicsListHelper::GetPhysicsListHelper()->RegisterProcess(reactionProcess, G4Alpha::Alpha());
    }

}