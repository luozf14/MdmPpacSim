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
/// \file PrimaryGeneratorAction.cc
/// \brief Implementation of the MdmPpacSim::PrimaryGeneratorAction class

#include "PrimaryGeneratorAction.hh"
#include "BeamEmittance.hh"
#include "HistoManager.hh"
#include "RunAction.hh"

#include "G4Box.hh"
#include "G4ChargedGeantino.hh"
#include "G4IonTable.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4RootAnalysisReader.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
using G4AnalysisReader = G4RootAnalysisReader;
namespace MdmPpacSim
{

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  PrimaryGeneratorAction::PrimaryGeneratorAction(RunAction *runAction,
                                                 G4bool isTargetChamber,
                                                 G4int procNum)
      : fRunAction(runAction), fIsTargetChamber(isTargetChamber),
        fProcessNum(procNum)
  {
    G4int n_particle = 1;
    fParticleGun = new G4ParticleGun(n_particle);

    // default particle kinematic
    G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition *particle =
        particleTable->FindParticle("chargedgeantino");
    fParticleGun->SetParticleDefinition(particle);
    G4ThreeVector particleDirection(0., 0., 1);
    fParticleGun->SetParticleMomentumDirection(particleDirection);
    fNumGenerated = 0;
    if (fIsTargetChamber)
    {
      fBeamEmittance = new BeamEmittance(24.0e-6); //[m*rad]
      fBeamEmittance->SetBeamSigmaX(5.0e-3);       //[m]
      fBeamEmittance->SetBeamSigmaY(5.0e-3);       //[m]
      fBeamEmittance->SetBeamPhiX(44. * deg);
      fBeamEmittance->SetBeamPhiY(0. * deg);
      fBeamEmittance->InitBeamEmittance(); // always call this function after
                                           // setting beam parameters
    }
    else
    {
      auto analysisReader = G4AnalysisReader::Instance();
      analysisReader->SetFileName("Stage1_" + std::to_string(fProcessNum) +
                                  ".root");
      G4int ntupleId = analysisReader->GetNtuple("Stage1Data");
      analysisReader->SetNtupleIColumn(ntupleId, "SlitBoxTransmitted",
                                       fTransmitted);
      analysisReader->SetNtupleDColumn(ntupleId, "SlitBoxEnergy", fEnergy);
      analysisReader->SetNtupleDColumn(ntupleId, "MDMPositionX", fPositionX);
      analysisReader->SetNtupleDColumn(ntupleId, "MDMPositionY", fPositionY);
      analysisReader->SetNtupleDColumn(ntupleId, "MDMAngleX", fAngleX);
      analysisReader->SetNtupleDColumn(ntupleId, "MDMAngleY", fAngleY);
      analysisReader->SetNtupleIColumn(ntupleId, "SlitBoxMass", fMass);
      analysisReader->SetNtupleIColumn(ntupleId, "SlitBoxCharge", fCharge);
    }
  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  PrimaryGeneratorAction::~PrimaryGeneratorAction()
  {
    delete fParticleGun;
    if (fIsTargetChamber)
      delete fBeamEmittance;
  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  void PrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent)
  {
    if (fIsTargetChamber)
    {
      G4ParticleDefinition *particle = fParticleGun->GetParticleDefinition();
      if (particle == G4ChargedGeantino::ChargedGeantino())
      {
        G4ParticleDefinition *ion;
        G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
        G4int Z = 8, A = 18;
        if (particleTable->GetIonTable()->FindIon(Z, A, 0.0))
          ion = particleTable->GetIonTable()->FindIon(Z, A, 0.0);
        else
          ion = particleTable->GetIonTable()->GetIon(Z, A, 0.0);
        // ion = particleTable->FindParticle("alpha");
        fParticleGun->SetParticleDefinition(ion);
      }

      G4double energy = fBeamEnergy * G4UniformRand();
      fParticleGun->SetParticleEnergy(energy * MeV);
      // fParticleGun->SetParticleEnergy(G4RandGauss::shoot(
      // fBeamEnergy * MeV, 0.005 * fBeamEnergy * MeV / 2.355));

      // pencil beam
      // fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1));
      // fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,-3.*mm));

      // beam emittance
      std::array<G4double, 4> beam = fBeamEmittance->GetBeam();
      G4double z0 = -2. * mm;
      fParticleGun->SetParticlePosition(
          G4ThreeVector(beam[0] * m, beam[2] * m, z0));
      G4ThreeVector momentumDirection =
          G4ThreeVector(beam[1], beam[3], 1.0).unit();
      fParticleGun->SetParticleMomentumDirection(momentumDirection);

      fRunAction->GetHistoManager()->SetBeam(
          fParticleGun->GetParticleEnergy(), fParticleGun->GetParticlePosition(),
          fParticleGun->GetParticleMomentumDirection());
    }
    else
    {
      auto analysisReader = G4AnalysisReader::Instance();
      analysisReader->GetNtupleRow();
      fNumGenerated++;

      fEnergy = (fTransmitted) ? fEnergy : 0;
      fMass = (fTransmitted) ? fMass : 4;
      fCharge = (fTransmitted) ? fCharge : 2;

      G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
      G4ParticleDefinition *ion;
      if (particleTable->GetIonTable()->FindIon(fCharge, fMass, 0.0))
        ion = particleTable->GetIonTable()->FindIon(fCharge, fMass, 0.0);
      else
        ion = particleTable->GetIonTable()->GetIon(fCharge, fMass, 0.0);
      fParticleGun->SetParticleDefinition(ion);

      fParticleGun->SetParticlePosition(
          G4ThreeVector(fPositionX * cm, fPositionY * cm, 0.));
      G4double dircX = std::tan(fAngleX / 180. * M_PI);
      G4double dircY =
          std::tan(fAngleY / 180. * M_PI) * std::sqrt(dircX * dircX + 1.);
      fParticleGun->SetParticleMomentumDirection(
          G4ThreeVector(dircX, dircY, 1.).unit());

      fParticleGun->SetParticleEnergy(fEnergy * MeV);
    }
    fParticleGun->GeneratePrimaryVertex(anEvent);

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  }
} // namespace MdmPpacSim
