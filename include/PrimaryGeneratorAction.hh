//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file PrimaryGeneratorAction.hh
/// \brief Definition of the MdmPpacSim::PrimaryGeneratorAction class

#ifndef MdmPpacSimPrimaryGeneratorAction_h
#define MdmPpacSimPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "globals.hh"

#include "TF1.h"

class G4ParticleGun;
class G4Event;
class G4Box;
/// The primary generator action class with particle gun.
///
/// The default kinematic is a 6 MeV gamma, randomly distribued
/// in front of the phantom across 80% of the (X,Y) phantom size.

namespace MdmPpacSim
{
    class BeamEmittance;
    class RunAction;

    class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
    {
    public:
        PrimaryGeneratorAction(RunAction *, G4bool, G4int);
        ~PrimaryGeneratorAction() override;

        // method from the base class
        void GeneratePrimaries(G4Event *) override;

        // method to access particle gun
        const G4ParticleGun *GetParticleGun() const { return fParticleGun; }
        void SetBeamEnergy(G4double e) { fBeamEnergy = e; }
        G4double GetBeamEnergy() { return fBeamEnergy; }

        void SetIsTargetChamber(G4bool isTargetChamber) { fIsTargetChamber = isTargetChamber; }
        void SetBeamEnergyDistribution();
        

    private:
        RunAction *fRunAction = nullptr;
        G4bool fIsTargetChamber;
        G4int fProcessNum;
        G4ParticleGun *fParticleGun = nullptr; // pointer a to G4 gun class
        G4double fBeamEnergy;
        BeamEmittance *fBeamEmittance = nullptr;
        TF1 *fBeamEnergyDistri = nullptr;

        G4int fNumGenerated;
        G4int fTransmitted;
        G4double fEnergy;
        G4double fPositionX;
        G4double fPositionY;
        G4double fAngleX;
        G4double fAngleY;
        G4int fMass;
        G4int fCharge;
    };

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
