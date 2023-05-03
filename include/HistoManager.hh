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
/// \file include/HistoManager.hh
/// \brief Definition of the HistoManager class
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef HistoManager_h
#define HistoManager_h 1

#include "g4analysis.hh"
#include "globals.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
namespace MdmPpacSim
{

    class HistoManager
    {
    public:
        HistoManager();
        ~HistoManager();

        void Book(G4int);
        void Save();

        void SetIsTargetChamber(G4bool isTargetChamber) { fIsTargetChamber = isTargetChamber; }

        void SetDeltaE(G4double edep, G4double time)
        {
            DeltaEEdep = edep;
            DeltETime = time;
        }

        void SetSi(G4double edep, G4double time, G4ThreeVector hitPosition, G4ThreeVector hitLocalPosition,
                   G4int frontStripNo, G4int backStripNo)
        {
            SiEdep = edep;
            SiTime = time;
            SiHitPosition = hitPosition;
            SiHitLocalPosition = hitLocalPosition;
            SiFrontStripNo = frontStripNo;
            SiBackStripNo = backStripNo;
        }

        void SetSlitBox(G4bool accepted, G4bool transmitted, G4ThreeVector hitPosition, G4ThreeVector hitLocalPosition,
                        G4ThreeVector hitMomentum, G4ThreeVector hitLocalMomentum,
                        G4int charge, G4int mass, G4double energy, G4double time, G4double scatteredAngleX, G4double scatteredAngleY)
        {
            SlitBoxAccepted = accepted;
            SlitBoxTransmitted = transmitted;
            SlitBoxHitPosition = hitPosition;
            SlitBoxHitLocalPosition = hitLocalPosition;
            SlitBoxHitMomentum = hitMomentum;
            SlitBoxHitLocalMomentum = hitLocalMomentum;
            SlitBoxCharge = charge;
            SlitBoxMass = mass;
            SlitBoxEnergy = energy;
            SlitBoxTime = time;
            ScatteredAngleX = scatteredAngleX;
            ScatteredAngleY = scatteredAngleY;
        }

        void SetMDMTraceResult(G4double positionX, G4double positionY,
                               G4double angleX, G4double angleY)
        {
            MDMPositionX = positionX;
            MDMPositionY = positionY;
            MDMAngleX = angleX;
            MDMAngleY = angleY;
        }

        void SetBeam(G4double energy, G4ThreeVector position, G4ThreeVector momentumDirection)
        {
            BeamEnergy = energy;
            BeamPosition = position;
            BeamMomentumDirection = momentumDirection;
        }

        void SetPpac1(G4double energy, G4double time, G4double positionX, G4double positionY)
        {
            PPAC1Energy = energy;
            PPAC1Time = time;
            PPAC1PositionX = positionX;
            PPAC1PositionY = positionY;
        }

        void SetPpac2(G4double energy, G4double time, G4double positionX, G4double positionY, G4bool completed, G4double tof)
        {
            PPAC2Energy = energy;
            PPAC2Time = time;
            PPAC2PositionX = positionX;
            PPAC2PositionY = positionY;
            Completed = completed;
            PPACTof = tof;
        }

        void FillNtuple();

    private:
        G4bool fIsTargetChamber;
        G4bool fFactoryOn;

        // DeltaE detector
        G4double DeltaEEdep;
        G4double DeltETime;

        // Si detector
        G4double SiEdep;
        G4double SiTime;
        G4ThreeVector SiHitPosition;
        G4ThreeVector SiHitLocalPosition;
        G4int SiFrontStripNo;
        G4int SiBackStripNo;

        // Slit box
        G4bool SlitBoxAccepted;
        G4bool SlitBoxTransmitted;
        G4ThreeVector SlitBoxHitPosition;
        G4ThreeVector SlitBoxHitLocalPosition;
        G4ThreeVector SlitBoxHitMomentum;
        G4ThreeVector SlitBoxHitLocalMomentum;
        G4int SlitBoxCharge;
        G4int SlitBoxMass;
        G4double SlitBoxEnergy;
        G4double SlitBoxTime;

        // MDMTrace results
        G4double MDMPositionX;
        G4double MDMPositionY;
        G4double MDMAngleX;
        G4double MDMAngleY;
        G4double ScatteredAngleX;
        G4double ScatteredAngleY;

        // Beam
        G4double BeamEnergy;
        G4ThreeVector BeamPosition;
        G4ThreeVector BeamMomentumDirection;

        // PPAC1
        G4double PPAC1Energy;
        G4double PPAC1Time;
        G4double PPAC1PositionX;
        G4double PPAC1PositionY;

        // PPAC2
        G4double PPAC2Energy;
        G4double PPAC2Time;
        G4double PPAC2PositionX;
        G4double PPAC2PositionY;

        // PPAC
        G4bool Completed;
        G4double PPACTof;

    };

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
}
#endif
