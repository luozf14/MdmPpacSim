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
/// \file EventAction.cc
/// \brief Implementation of the MdmPpacSim::EventAction class

#include "EventAction.hh"
#include "RunAction.hh"
#include "HistoManager.hh"
#include "DetectorHit.hh"
#include "MDMTrace.h"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

namespace MdmPpacSim
{
    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    EventAction::EventAction(RunAction *runAction)
        : G4UserEventAction(), fRunAction(runAction)
    {
        fMDMTrace = new MDMTrace();
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    EventAction::~EventAction()
    {
        delete fMDMTrace;
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    void EventAction::BeginOfEventAction(const G4Event *)
    {
        G4SDManager *sdManager = G4SDManager::GetSDMpointer();
        if (fIsTargetChamber)
        {
            fHCID_DeltaE = sdManager->GetCollectionID("DeltaEHitsCollection");
            fHCID_Si = sdManager->GetCollectionID("SiHitsCollection");
            fHCID_SlitBox = sdManager->GetCollectionID("SlitBoxHitsCollection");
            fMDMTrace->SetMDMAngle(fMdmAngleInDeg);
            fMDMTrace->SetMDMProbe(fMdmDipoleProbe, fMdmMultipoleProbe);
        }
        else
        {
            fHCID_Ppac1 = sdManager->GetCollectionID("Ppac1HitsCollection");
            fHCID_Ppac2 = sdManager->GetCollectionID("Ppac2HitsCollection");
        }
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    void EventAction::EndOfEventAction(const G4Event *aEvent)
    {
        fHistoManager = fRunAction->GetHistoManager();

        G4HCofThisEvent *hce = aEvent->GetHCofThisEvent();
        if (!hce)
        {
            G4ExceptionDescription msg;
            msg << "No hits collection of this event found.\n";
            return;
        }

        if (fIsTargetChamber)
        {
            // DeltaE
            DetectorHitsCollection *hcDeltaE = (DetectorHitsCollection *)aEvent->GetHCofThisEvent()->GetHC(fHCID_DeltaE);
            G4int nofHitsDeltaE = hcDeltaE->GetSize();
            G4double deltaEEdep = 0.;
            G4double deltaETime = 0.;
            for (int i = 0; i < nofHitsDeltaE; i++)
            {
                deltaEEdep += (*hcDeltaE)[i]->GetEnergyDeposit();
                deltaETime += (*hcDeltaE)[i]->GetHitTime();
            }
            deltaEEdep = (nofHitsDeltaE > 0) ? G4RandGauss::shoot(deltaEEdep, fSiDetectorEnergyResolution * deltaEEdep / 2.355) : 0.;
            deltaETime = (nofHitsDeltaE > 0) ? (deltaETime / nofHitsDeltaE + 0.25 * ns * G4RandFlat::shoot()) : 0.;
            fHistoManager->SetDeltaE(deltaEEdep, deltaETime);

            // Si detector
            DetectorHitsCollection *hcSi = (DetectorHitsCollection *)aEvent->GetHCofThisEvent()->GetHC(fHCID_Si);
            G4int nofHitsSi = hcSi->GetSize();
            G4double siEdep = 0.;
            G4double siTime = 0.;
            G4ThreeVector siHitPos(0., 0., 0.);
            G4ThreeVector siHitLocalPos(0., 0., 0.);
            G4int siMass = 0;
            G4int siCharge = 0;
            G4int siTrackID = 0;
            for (int i = 0; i < nofHitsSi; i++)
            {
                siEdep += (*hcSi)[i]->GetEnergyDeposit();
                siTime += (*hcSi)[i]->GetHitTime();
                siHitPos += (*hcSi)[i]->GetPosition();
                siHitLocalPos += (*hcSi)[i]->GetLocalPosition();
                siMass += (*hcSi)[i]->GetParticleMass();
                siCharge += (*hcSi)[i]->GetParticleCharge();
                siTrackID += (*hcSi)[i]->GetTrackID();
            }
            siEdep = (nofHitsSi > 0) ? G4RandGauss::shoot(siEdep, fSiDetectorEnergyResolution * siEdep / 2.355) : 0.;
            siTime = (nofHitsSi > 0) ? (siTime / nofHitsDeltaE + 0.25 * ns * G4RandFlat::shoot()) : 0.;
            siHitPos = (nofHitsSi > 0) ? siHitPos / nofHitsSi : G4ThreeVector(0., 0., 0.);
            siHitLocalPos = (nofHitsSi > 0) ? siHitLocalPos / nofHitsSi : G4ThreeVector(0., 0., 0.);
            G4int siFrontStripNo = std::floor((siHitLocalPos.x() / cm + 2.5) / (5. / 16.));
            G4int siBackStripNo = std::floor((siHitLocalPos.y() / cm + 2.5) / (5. / 16.));
            siMass = (nofHitsSi > 0) ? siMass / nofHitsSi : 0;
            siCharge = (nofHitsSi > 0) ? siCharge / nofHitsSi : 0;
            siTrackID = (nofHitsSi > 0) ? siTrackID / nofHitsSi : 0;
            fHistoManager->SetSi(siEdep, siTime, siHitPos, siHitLocalPos, siFrontStripNo, siBackStripNo, siMass, siCharge, siTrackID);

            // Slit box
            DetectorHitsCollection *hcSlitBox = (DetectorHitsCollection *)aEvent->GetHCofThisEvent()->GetHC(fHCID_SlitBox);
            G4int nofHitsSlitBox = hcSlitBox->GetSize();
            G4bool slitBoxAccepted = false;
            G4bool slitBoxTransmitted = false;
            G4ThreeVector slitBoxHitPosition(0., 0., 0.);
            G4ThreeVector slitBoxHitLocalPosition(0., 0., 0.);
            G4ThreeVector slitBoxHitMomentum(0., 0., 0.);
            G4ThreeVector slitBoxHitLocalMomentum(0., 0., 0.);
            G4int slitBoxCharge = 0;
            G4int slitBoxMass = 0;
            G4double slitBoxEnergy = 0.;
            G4double slitBoxTime = 0.;
            G4double scatteredAngleX = 0.;
            G4double scatteredAngleY = 0.;
            G4double mdmPositionX = 0.;
            G4double mdmPositionY = 0.;
            G4double mdmAngleX = 0.;
            G4double mdmAngleY = 0.;
            G4int slitBoxTrackID = 0;
            if (nofHitsSlitBox > 0)
            {
                slitBoxAccepted = true;
                slitBoxHitPosition = (*hcSlitBox)[nofHitsSlitBox - 1]->GetPosition();
                slitBoxHitLocalPosition = (*hcSlitBox)[nofHitsSlitBox - 1]->GetLocalPosition();
                slitBoxHitMomentum = (*hcSlitBox)[nofHitsSlitBox - 1]->GetMomentum();
                slitBoxHitLocalMomentum = (*hcSlitBox)[nofHitsSlitBox - 1]->GetLocalMomentum();
                slitBoxCharge = (*hcSlitBox)[nofHitsSlitBox - 1]->GetParticleCharge();
                slitBoxMass = (*hcSlitBox)[nofHitsSlitBox - 1]->GetParticleMass();
                slitBoxEnergy = (*hcSlitBox)[nofHitsSlitBox - 1]->GetKineticEnergy();
                slitBoxTime = (*hcSlitBox)[nofHitsSlitBox - 1]->GetHitTime();
                slitBoxTrackID = (*hcSlitBox)[nofHitsSlitBox - 1]->GetTrackID();

                fMDMTrace->SetScatteredMass(slitBoxMass);
                if (slitBoxCharge == 6 && slitBoxMass == 12)
                {
                    fMDMTrace->SetScatteredCharge(slitBoxCharge - 1);
                }
                else if (slitBoxCharge == 6 && slitBoxMass == 13)
                {
                    fMDMTrace->SetScatteredCharge(slitBoxCharge - 1);
                }
                else if (slitBoxCharge == 8 && slitBoxMass == 16)
                {
                    fMDMTrace->SetScatteredCharge(slitBoxCharge - 1);
                }
                else
                {
                    fMDMTrace->SetScatteredCharge(slitBoxCharge);
                }
                fMDMTrace->SetScatteredEnergy(slitBoxEnergy);
                scatteredAngleX = (*hcSlitBox)[nofHitsSlitBox - 1]->GetAngleX(); // [deg]
                scatteredAngleY = (*hcSlitBox)[nofHitsSlitBox - 1]->GetAngleY(); // [deg]
                fMDMTrace->SetScatteredAngle(scatteredAngleX, scatteredAngleY);
                // G4cout << "Set MDM Angle: " << fMDMTrace->GetMDMAngle() << G4endl;
                // G4cout << "Set Scattered Angle: " << fMDMTrace->GetScatteredAngle() << ", "
                //        << "Energy: " << fMDMTrace->GetScatteredEnergy() << G4endl;
                // G4cout << "Set Scattered Mass: " << fMDMTrace->GetScatteredMass() << ", "
                //        << "Charge: " << fMDMTrace->GetScatteredCharge() << G4endl;
                fMDMTrace->SendRay();
                mdmPositionX = fMDMTrace->GetFirstWireX();   // [cm]
                mdmPositionY = fMDMTrace->GetFirstWireY();   // [cm]
                mdmAngleX = fMDMTrace->GetFirstWireXAngle(); // [deg]
                mdmAngleY = fMDMTrace->GetFirstWireYAngle(); // [deg]
                if (std::abs(mdmPositionX) <= 20. && std::abs(mdmPositionY) <= 5.)
                {
                    slitBoxTransmitted = true;
                }
                else
                {
                    slitBoxTransmitted = false;
                }
            }
            else
            {
                slitBoxAccepted = false;
            }
            fHistoManager->SetSlitBox(slitBoxAccepted, slitBoxTransmitted, slitBoxHitPosition, slitBoxHitLocalPosition, slitBoxHitMomentum, slitBoxHitLocalMomentum, slitBoxCharge, slitBoxMass, slitBoxEnergy, slitBoxTime, scatteredAngleX, scatteredAngleY, slitBoxTrackID);
            fHistoManager->SetMDMTraceResult(mdmPositionX, mdmPositionY, mdmAngleX, mdmAngleY);
        }
        else
        {
            // PPAC1
            DetectorHitsCollection *hcPpac1 = (DetectorHitsCollection *)aEvent->GetHCofThisEvent()->GetHC(fHCID_Ppac1);
            G4int nofHitsPpac1 = hcPpac1->GetSize();
            G4ThreeVector ppac1HitLocalPosition(0., 0., 0.);
            G4double ppac1HitEnergy = 0.;
            G4double ppac1HitTime = 0.;
            for (int i = 0; i < nofHitsPpac1; i++)
            {
                ppac1HitEnergy += (*hcPpac1)[i]->GetEnergyDeposit();
                ppac1HitTime += (*hcPpac1)[i]->GetHitTime();
                ppac1HitLocalPosition += (*hcPpac1)[i]->GetLocalPosition();
            }
            ppac1HitLocalPosition = (nofHitsPpac1 > 0) ? ppac1HitLocalPosition / nofHitsPpac1 : G4ThreeVector(0., 0., 0.);
            // ppac1HitTime = (nofHitsPpac1 > 0) ? (ppac1HitTime / nofHitsPpac1 + fTdcResolutionInNs * ns * G4RandFlat::shoot()) : 0.;
            ppac1HitTime = (nofHitsPpac1 > 0) ? G4RandGauss::shoot(ppac1HitTime / nofHitsPpac1, fTdcResolutionInNs * ns / 2.355) : 0.;
            ppac1HitEnergy = (nofHitsPpac1 > 0) ? ppac1HitEnergy / nofHitsPpac1 : 0.;
            G4double ppac1PosX = (nofHitsPpac1 > 0) ? G4RandGauss::shoot(ppac1HitLocalPosition.x(), fPpacPositionResolutionInMm * mm / 2.355) : 0.;
            G4double ppac1PosY = (nofHitsPpac1 > 0) ? G4RandGauss::shoot(ppac1HitLocalPosition.y(), fPpacPositionResolutionInMm * mm / 2.355) : 0.;
            fHistoManager->SetPpac1(ppac1HitEnergy, ppac1HitTime, ppac1PosX / cm, ppac1PosY / cm);

            // PPAC2
            DetectorHitsCollection *hcPpac2 = (DetectorHitsCollection *)aEvent->GetHCofThisEvent()->GetHC(fHCID_Ppac2);
            G4int nofHitsPpac2 = hcPpac2->GetSize();
            G4ThreeVector ppac2HitLocalPosition(0., 0., 0.);
            G4double ppac2HitEnergy = 0.;
            G4double ppac2HitTime = 0.;
            G4bool completed = false;
            G4double ppacTof = 0.;
            for (int i = 0; i < nofHitsPpac2; i++)
            {
                ppac2HitEnergy += (*hcPpac2)[i]->GetEnergyDeposit();
                ppac2HitTime += (*hcPpac2)[i]->GetHitTime();
                ppac2HitLocalPosition += (*hcPpac2)[i]->GetLocalPosition();
            }
            ppac2HitLocalPosition = (nofHitsPpac2 > 0) ? ppac2HitLocalPosition / nofHitsPpac2 : G4ThreeVector(0., 0., 0.);
            // ppac2HitTime = (nofHitsPpac2 > 0) ? (ppac2HitTime / nofHitsPpac2 + fTdcResolutionInNs * ns * G4RandFlat::shoot()) : 0.;
            ppac2HitTime = (nofHitsPpac1 > 0) ? G4RandGauss::shoot(ppac2HitTime / nofHitsPpac2, fTdcResolutionInNs * ns / 2.355) : 0.;
            ppac2HitEnergy = (nofHitsPpac2 > 0) ? ppac2HitEnergy / nofHitsPpac2 : 0.;
            G4double ppac2PosX = (nofHitsPpac2 > 0) ? G4RandGauss::shoot(ppac2HitLocalPosition.x(), fPpacPositionResolutionInMm * mm / 2.355) : 0.;
            G4double ppac2PosY = (nofHitsPpac2 > 0) ? G4RandGauss::shoot(ppac2HitLocalPosition.y(), fPpacPositionResolutionInMm * mm / 2.355) : 0.;
            completed = (nofHitsPpac1 > 0 && nofHitsPpac2 > 0) ? true : false;
            ppacTof = (completed) ? ppac2HitTime - ppac1HitTime : 0.;
            fHistoManager->SetPpac2(ppac2HitEnergy, ppac2HitTime, ppac2PosX / cm, ppac2PosY / cm, completed, ppacTof);
        }
        fHistoManager->FillNtuple();
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
    void EventAction::SetEventActionParams(std::map<std::string, G4double> eventParams)
    {
        for (auto it : eventParams)
        {
            if (it.first == "SiDetectorEnergyResolution")
            {
                fSiDetectorEnergyResolution = it.second;
            }
            else if (it.first == "TdcResolutionInNs")
            {
                fTdcResolutionInNs = it.second;
            }
            else if (it.first == "PpacPositionResolutionInMm")
            {
                fPpacPositionResolutionInMm = it.second;
            }
            else if (it.first == "MdmMultipoleProbe")
            {
                fMdmMultipoleProbe = it.second;
            }
            else if (it.first == "MdmDipoleProbe")
            {
                fMdmDipoleProbe = it.second;
            }
            else if (it.first == "MdmAngleInDeg")
            {
                fMdmAngleInDeg = it.second;
            }
        }
    }
}
