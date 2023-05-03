#ifndef MdmPpacSimBeamEmittance_h
#define MdmPpacSimBeamEmittance_h 1
#include "globals.hh"
#include "G4ThreeVector.hh"
namespace MdmPpacSim
{
    class BeamEmittance
    {
    public:
        BeamEmittance(G4double emittance);
        virtual ~BeamEmittance();

        void InitBeamEmittance();

        void SetBeamEmittance(G4double emittance) { fBeamEmittance = emittance; }
        G4double GetBeamEmittance() { return fBeamEmittance; }
        void SetBeamSigmaX(G4double sigmaX) { fBeamSigmaX = sigmaX; }
        G4double GetBeamSigmaX() { return fBeamSigmaX; }
        void SetBeamSigmaY(G4double sigmaY) { fBeamSigmaY = sigmaY; }
        G4double GetBeamSigmaY() { return fBeamSigmaY; }
        void SetBeamPhiX(G4double phiX) { fBeamPhiX = phiX; }
        G4double GetBeamPhiX() { return fBeamPhiX; }
        void SetBeamPhiY(G4double phiY) { fBeamPhiY = phiY; }
        G4double GetBeamPhiY() { return fBeamPhiY; }

        std::array<G4double, 4> GetBeam();

    private:
        G4double fBeamEmittance;
        G4double fBeamSigmaX;
        G4double fBeamSigmaY;
        G4double fBeamPhiX;
        G4double fBeamPhiY;

        G4double fBeta_x;
        G4double fTan2_x;
        G4double fGamma_x;
        G4double fAlpha_x;

        G4double fBeta_y;
        G4double fTan2_y;
        G4double fGamma_y;
        G4double fAlpha_y;
    };
}
#endif
