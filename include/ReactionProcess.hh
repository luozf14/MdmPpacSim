#ifndef ReactionProcess_h
#define ReactionProcess_h

#include "G4VDiscreteProcess.hh"
#include "TGraph.h"
namespace MdmPpacSim
{
    class ReactionProcess : public G4VDiscreteProcess
    {
    public:
        ReactionProcess(const G4String &name = "Reaction");
        ~ReactionProcess();

        G4double GetMeanFreePath(const G4Track &, G4double,
                                 G4ForceCondition *);
        G4VParticleChange *PostStepDoIt(const G4Track &, const G4Step &);

        void StartTracking(G4Track *);

        G4double GetLightProductMass()
        {
            return fLightProductMass;
        }
        G4double GetLightProductCharge()
        {
            return fLightProductCharge;
        }

        void SetParams(std::map<std::string, G4int>);
        void SetWaveFunction(TGraph*);

        G4VParticleChange *TwoBody(const G4Track &aTrack, const G4Step &aStep, int Zt, int At, int Z1, int A1, int Z2, int A2, double Ex1, double Ex2);

        G4VParticleChange *TwoBodyDecay(const G4Track &aTrack, const G4Step &aStep,
                                        int Zt, int At,             // target
                                        int Z1, int A1, double Ex1, // light recoil
                                        int Z2, int A2, double Ex2,
                                        int flag); // heavy recoil

        G4VParticleChange *TwoBodyAlphaDecay(const G4Track &aTrack, const G4Step &aStep,
                                             int Zt, int At,              // target
                                             int Z1, int A1, double Ex1,  // light recoil
                                             int Z2, int A2, double Ex2); // heavy recoil

        G4VParticleChange *TwoBodyNeutronDecay(const G4Track &aTrack, const G4Step &aStep,
                                               int Zt, int At,              // target
                                               int Z1, int A1, double Ex1,  // light recoil
                                               int Z2, int A2, double Ex2); // heavy recoil

        G4VParticleChange *TrojanHorse(const G4Track &aTrack, const G4Step &aStep,
                                       int Zt, int At,              // target
                                       int Z1, int A1, double Ex1,  // light recoil
                                       int Z2, int A2, double Ex2,  // heavy recoil
                                       int Zs, int As, double Exs); // spectator

        G4VParticleChange *TrojanHorseInverse(const G4Track &aTrack, const G4Step &aStep,
                                       int Zt, int At,              // target
                                       int Z1, int A1, double Ex1,  // light recoil
                                       int Z2, int A2, double Ex2,  // heavy recoil
                                       int Zs, int As, double Exs); // spectator

    private:
        G4double fTargetMass;
        G4double fTargetCharge;
        G4double fLightProductMass;
        G4double fLightProductCharge;
        G4double fHeavyProductMass;
        G4double fHeavyProductCharge;
        G4double fReactionPosition;
        TGraph *fWaveFunction;
    };
}
#endif
