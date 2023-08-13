#ifndef ReactionPhysics_h
#define ReactionPhysics_h 1

#include "G4VPhysicsConstructor.hh"
#include <string>
#include <map>
#include "TGraph.h"
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

        void SetWaveFunction(TGraph *waveFunc)
        {
            fWaveFunction = waveFunc;
        };

    private:
        std::map<std::string, G4int> fReactionParams;
        TGraph *fWaveFunction;
    };
}
#endif
