#include "BeamEmittance.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "Randomize.hh"

namespace MdmPpacSim
{
    BeamEmittance::BeamEmittance(G4double emittance)
    {
        fBeamEmittance = emittance;
        fBeamSigmaX = std::sqrt(emittance);
        fBeamSigmaY = std::sqrt(emittance);
        fBeamPhiX = 0. * deg;
        fBeamPhiY = 0. * deg;
    }

    BeamEmittance::~BeamEmittance()
    {
    }

    void BeamEmittance::InitBeamEmittance()
    {
        // x
        fBeta_x = fBeamSigmaX * fBeamSigmaX / fBeamEmittance;
        fTan2_x = std::pow(std::tan(2. * fBeamPhiX), 2.);
        if (0. * deg < fBeamPhiX && fBeamPhiX < 45. * deg)
        {
            fGamma_x = (fBeta_x * (fTan2_x + 2.) - 2. * std::sqrt((fTan2_x + 1.) * fBeta_x * fBeta_x - fTan2_x)) / fTan2_x;
            fAlpha_x = -std::sqrt(fBeta_x * fGamma_x - 1.);
        }
        else if (45. * deg < fBeamPhiX && fBeamPhiX < 90. * deg)
        {
            fGamma_x = (fBeta_x * (fTan2_x + 2.) + 2. * std::sqrt((fTan2_x + 1.) * fBeta_x * fBeta_x - fTan2_x)) / fTan2_x;
            fAlpha_x = -std::sqrt(fBeta_x * fGamma_x - 1.);
        }
        else if (90. * deg < fBeamPhiX && fBeamPhiX < 135. * deg)
        {
            fGamma_x = (fBeta_x * (fTan2_x + 2.) + 2. * std::sqrt((fTan2_x + 1.) * fBeta_x * fBeta_x - fTan2_x)) / fTan2_x;
            fAlpha_x = std::sqrt(fBeta_x * fGamma_x - 1.);
        }
        else if (135. * deg < fBeamPhiX && fBeamPhiX < 180. * deg)
        {
            fGamma_x = (fBeta_x * (fTan2_x + 2.) - 2. * std::sqrt((fTan2_x + 1.) * fBeta_x * fBeta_x - fTan2_x)) / fTan2_x;
            fAlpha_x = std::sqrt(fBeta_x * fGamma_x - 1.);
        }
        else if (fBeamPhiX == (0. * deg) || fBeamPhiX == (90. * deg))
        {
            fGamma_x = 1. / fBeta_x;
            fAlpha_x = 0;
        }
        if ((std::tan(2. * fBeamPhiX) - (2. * fAlpha_x / (fGamma_x - fBeta_x))) / std::tan(2. * fBeamPhiX) > 1e-3)
        {
            std::cerr << "Error: tan2phiX_theory!=tan2phiX_real" << std::endl;
            exit(1);
        }

        // y
        fBeta_y = fBeamSigmaY * fBeamSigmaY / fBeamEmittance;
        fTan2_y = std::pow(std::tan(2 * fBeamPhiY), 2.);
        if (0. * deg < fBeamPhiY && fBeamPhiY < 45. * deg)
        {
            fGamma_y = (fBeta_y * (fTan2_y + 2.) - 2. * std::sqrt((fTan2_y + 1.) * fBeta_y * fBeta_y - fTan2_y)) / fTan2_y;
            fAlpha_y = -std::sqrt(fBeta_y * fGamma_y - 1.);
        }
        else if (45. * deg < fBeamPhiY && fBeamPhiY < 90. * deg)
        {
            fGamma_y = (fBeta_y * (fTan2_y + 2.) + 2. * std::sqrt((fTan2_y + 1.) * fBeta_y * fBeta_y - fTan2_y)) / fTan2_y;
            fAlpha_y = -std::sqrt(fBeta_y * fGamma_y - 1.);
        }
        else if (90. * deg < fBeamPhiY && fBeamPhiY < 135. * deg)
        {
            fGamma_y = (fBeta_y * (fTan2_y + 2.) + 2. * std::sqrt((fTan2_y + 1.) * fBeta_y * fBeta_y - fTan2_y)) / fTan2_y;
            fAlpha_y = std::sqrt(fBeta_y * fGamma_y - 1.);
        }
        else if (135. * deg < fBeamPhiY && fBeamPhiY < 180. * deg)
        {
            fGamma_y = (fBeta_y * (fTan2_y + 2.) - 2. * std::sqrt((fTan2_y + 1.) * fBeta_y * fBeta_y - fTan2_y)) / fTan2_y;
            fAlpha_y = std::sqrt(fBeta_y * fGamma_y - 1.);
        }
        else if (fBeamPhiY == (0. * deg) || fBeamPhiY == (90. * deg))
        {
            fGamma_y = 1. / fBeta_y;
            fAlpha_y = 0;
        }
        if ((std::tan(2. * fBeamPhiY) - (2. * fAlpha_y / (fGamma_y - fBeta_y))) / std::tan(2. * fBeamPhiY) > 1e-3)
        {
            std::cerr << "Error: tan2phiY_theory!=tan2phiY_real" << std::endl;
            exit(1);
        }
    }

    std::array<G4double, 4> BeamEmittance::GetBeam()
    {
        G4double x0 = G4RandGauss::shoot(0, fBeamSigmaX / std::sqrt(6.));
        G4double p_x = G4RandGauss::shoot(0, fBeamSigmaX / std::sqrt(6.));
        G4double x0_prime = (p_x - fAlpha_x * x0) / fBeta_x;
        while (fGamma_x * x0 * x0 + 2. * fAlpha_x * x0 * x0_prime + fBeta_x * x0_prime * x0_prime > fBeamEmittance)
        {
            x0 = G4RandGauss::shoot(0, fBeamSigmaX / std::sqrt(6.));
            p_x = G4RandGauss::shoot(0, fBeamSigmaX / std::sqrt(6.));
            x0_prime = (p_x - fAlpha_x * x0) / fBeta_x;
        }

        G4double y0 = G4RandGauss::shoot(0, fBeamSigmaY / std::sqrt(6.));
        G4double p_y = G4RandGauss::shoot(0, fBeamSigmaY / std::sqrt(6.));
        G4double y0_prime = (p_y - fAlpha_y * x0) / fBeta_y;
        while (fGamma_y * y0 * y0 + 2. * fAlpha_y * y0 * y0_prime + fBeta_y * y0_prime * y0_prime > fBeamEmittance)
        {
            y0 = G4RandGauss::shoot(0, fBeamSigmaY / std::sqrt(6.));
            p_y = G4RandGauss::shoot(0, fBeamSigmaY / std::sqrt(6.));
            y0_prime = (p_y - fAlpha_y * y0) / fBeta_y;
        }

        std::array<G4double, 4> beam = {x0, x0_prime, y0, y0_prime};
        return beam;
    }
}