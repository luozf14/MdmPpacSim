#ifndef ReactionPhysics_h
#define ReactionPhysics_h 1

#include "G4VPhysicsConstructor.hh"
#include <string>
#include <map>
namespace MdmPpacSim
{
    class ReactionPhysics : public G4VPhysicsConstructor
    {
    public:
        ReactionPhysics(G4int verbose = 1);
        ReactionPhysics(const G4String &name);
        virtual ~ReactionPhysics();

        virtual void ConstructParticle();
        virtual void ConstructProcess();

        void SetReactionParams(std::map<std::string, G4int> params)
        {
            fReactionParams = params;
        };

    private:
        std::map<std::string, G4int> fReactionParams;
    };
}
#endif
