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
/// \file EventAction.hh
/// \brief Definition of the MdmPpacSim::EventAction class

#ifndef MdmPpacSimEventAction_h
#define MdmPpacSimEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

/// Event action class
///

namespace MdmPpacSim
{

    class RunAction;
    class HistoManager;
    class MDMTrace;

    class EventAction : public G4UserEventAction
    {
    public:
        EventAction(RunAction *);
        ~EventAction() override;

        void BeginOfEventAction(const G4Event *event) override;
        void EndOfEventAction(const G4Event *event) override;

        void SetEventActionParams(std::map<std::string, G4double> eventActionParams);
        void SetIsTargetChamber(G4bool isTargetChamber) { fIsTargetChamber = isTargetChamber; }

    private:
        RunAction *fRunAction = nullptr;
        HistoManager *fHistoManager = nullptr;
        MDMTrace *fMDMTrace=nullptr;
        std::map<std::string, G4double> fEventActionParams;
        G4bool fIsTargetChamber;
        G4double fSiDetectorEnergyResolution;
        G4double fTdcResolutionInNs;
        G4double fPpacPositionResolutionInMm;
        G4double fMdmMultipoleProbe;
        G4double fMdmDipoleProbe;
        G4double fMdmAngleInDeg;
        G4int fHCID_DeltaE;
        G4int fHCID_Si;
        G4int fHCID_SlitBox;
        G4int fHCID_Ppac1;
        G4int fHCID_Ppac2;
    };

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
