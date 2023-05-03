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
/// \file DetectorHit.hh
/// \brief Definition of the MdmPpacSim::DetectorHit class

#ifndef MdmPpacSimDetectorHit_h
#define MdmPpacSimDetectorHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "tls.hh"

namespace MdmPpacSim
{

    /// Detector hit class
    ///
    /// It defines data members to store the trackID, energy deposit, kinetic energy, momentum,
    /// and position of charged particles in a selected volume:
    /// - fTrackID, fEnergyDeposit, fKineticEnergy, fMomentum, fPosition

    class DetectorHit : public G4VHit
    {
    public:
        DetectorHit();
        DetectorHit(const DetectorHit &) = default;
        ~DetectorHit() override;

        // operators
        DetectorHit &operator=(const DetectorHit &) = default;
        G4bool operator==(const DetectorHit &) const;

        inline void *operator new(size_t);
        inline void operator delete(void *);

        // methods from base class
        void Print() override;

        // Set methods
        void SetTrackID(G4int track) { fTrackID = track; };
        void SetParticleName(G4String name) { fParticleName = name; };
        void SetParticleCharge(G4int charge) { fParticleCharge = charge; };
        void SetParticleMass(G4int mass) { fParticleMass = mass; };
        void SetEnergyDeposit(G4double de) { fEnergyDeposit = de; };
        void SetKineticEnergy(G4double ke) { fKineticEnergy = ke; };
        void SetHitTime(G4double time) { fHitTIme = time; };
        void SetMomentum(G4ThreeVector xyz) { fMomentum = xyz; };
        void SetLocalMomentum(G4ThreeVector xyz) { fLocalMomentum = xyz; };
        void SetPosition(G4ThreeVector xyz) { fPosition = xyz; };
        void SetLocalPosition(G4ThreeVector xyz) { fLocalPosition = xyz; };
        void SetAngleX(G4double angle) { fAngleX = angle; };
        void SetAngleY(G4double angle) { fAngleY = angle; };

        // Get methods
        G4int GetTrackID() const { return fTrackID; };
        G4String GetParticleName() const { return fParticleName; };
        G4int GetParticleCharge() const { return fParticleCharge; };
        G4int GetParticleMass() const { return fParticleMass; };
        G4double GetEnergyDeposit() const { return fEnergyDeposit; };
        G4double GetKineticEnergy() const { return fKineticEnergy; };
        G4double GetHitTime() const { return fHitTIme; };
        G4ThreeVector GetMomentum() const { return fMomentum; };
        G4ThreeVector GetLocalMomentum() const { return fLocalMomentum; };
        G4ThreeVector GetPosition() const { return fPosition; };
        G4ThreeVector GetLocalPosition() const { return fLocalPosition; };
        G4double GetAngleX() const { return fAngleX; };
        G4double GetAngleY() const { return fAngleY; };

    private:
        G4int fTrackID = -1;
        G4String fParticleName;
        G4int fParticleCharge = -1;
        G4int fParticleMass = -1;
        G4double fEnergyDeposit = 0.;
        G4double fKineticEnergy = 0.;
        G4double fHitTIme = 0.;
        G4ThreeVector fMomentum;
        G4ThreeVector fLocalMomentum;
        G4ThreeVector fPosition;
        G4ThreeVector fLocalPosition;
        G4double fAngleX;
        G4double fAngleY;
    };

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    typedef G4THitsCollection<DetectorHit> DetectorHitsCollection;

    extern G4ThreadLocal G4Allocator<DetectorHit> *DetectorHitAllocator;

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    inline void *DetectorHit::operator new(size_t)
    {
        if (!DetectorHitAllocator)
            DetectorHitAllocator = new G4Allocator<DetectorHit>;
        return (void *)DetectorHitAllocator->MallocSingle();
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    inline void DetectorHit::operator delete(void *hit)
    {
        DetectorHitAllocator->FreeSingle((DetectorHit *)hit);
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}

#endif
