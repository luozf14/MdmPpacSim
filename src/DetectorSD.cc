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
/// \file DetectorSD.cc
/// \brief Implementation of the MdmPpacSim::DetectorSD class

#include "DetectorSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

namespace MdmPpacSim
{

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    DetectorSD::DetectorSD(const G4String &name,
                           const G4String &hitsCollectionName)
        : G4VSensitiveDetector(name)
    {
        collectionName.insert(hitsCollectionName);
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    DetectorSD::~DetectorSD()
    {
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    void DetectorSD::Initialize(G4HCofThisEvent *hce)
    {
        // Create hits collection

        fHitsCollection = new DetectorHitsCollection(SensitiveDetectorName, collectionName[0]);

        // Add this collection in hce

        G4int hcID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
        hce->AddHitsCollection(hcID, fHitsCollection);
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    G4bool DetectorSD::ProcessHits(G4Step *aStep,
                                   G4TouchableHistory *)
    {
        // energy deposit
        G4double edep = aStep->GetTotalEnergyDeposit();

        if (edep == 0.)
            return false;

        DetectorHit *newHit = new DetectorHit();

        newHit->SetTrackID(aStep->GetTrack()->GetTrackID());
        newHit->SetParticleName(aStep->GetTrack()->GetParticleDefinition()->GetParticleName());
        newHit->SetParticleCharge(aStep->GetTrack()->GetParticleDefinition()->GetAtomicNumber());
        newHit->SetParticleMass(aStep->GetTrack()->GetParticleDefinition()->GetAtomicMass());
        newHit->SetEnergyDeposit(edep);
        newHit->SetHitTime(aStep->GetTrack()->GetGlobalTime());
        newHit->SetKineticEnergy(aStep->GetTrack()->GetKineticEnergy());
        newHit->SetPosition(aStep->GetPostStepPoint()->GetPosition());
        G4ThreeVector translatedPosition = aStep->GetPreStepPoint()->GetPosition() - aStep->GetPreStepPoint()->GetPhysicalVolume()->GetObjectTranslation();
        G4ThreeVector rotatedPosition = aStep->GetPreStepPoint()->GetPhysicalVolume()->GetObjectRotation()->inverse().operator()(translatedPosition);
        newHit->SetLocalPosition(rotatedPosition);
        newHit->SetMomentum(aStep->GetTrack()->GetMomentum());
        G4ThreeVector localMomentum = aStep->GetPreStepPoint()->GetPhysicalVolume()->GetObjectRotation()->inverse().operator()(aStep->GetTrack()->GetMomentum());
        newHit->SetLocalMomentum(localMomentum);
        G4ThreeVector momentumDirection = aStep->GetTrack()->GetDynamicParticle()->GetMomentumDirection();
        G4double angleX = std::atan(momentumDirection[0] / momentumDirection[2]) * 180. / M_PI;
        G4ThreeVector hitVolumePosition = aStep->GetPreStepPoint()->GetPhysicalVolume()->GetObjectTranslation();
        G4double angleY = std::asin(rotatedPosition[1] / aStep->GetPostStepPoint()->GetPosition().mag()) * 180. / M_PI;
        // G4double angleY = std::atan(rotatedPosition[1] / hitVolumePosition.mag()) * 180 / M_PI;
        newHit->SetAngleX(angleX);
        newHit->SetAngleY(angleY);
        fHitsCollection->insert(newHit);

        // newHit->Print();

        return true;
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    void DetectorSD::EndOfEvent(G4HCofThisEvent *)
    {
        return;
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
