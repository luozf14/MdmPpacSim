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
/// \file ActionInitialization.cc
/// \brief Implementation of the MdmPpacSim::ActionInitialization class

#include "ActionInitialization.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"
#include "HistoManager.hh"
#include "MDMTrace.h"

#include "G4Threading.hh"
namespace MdmPpacSim
{

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    ActionInitialization::ActionInitialization()
    {
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    ActionInitialization::~ActionInitialization()
    {
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    void ActionInitialization::BuildForMaster() const
    {
        HistoManager *histoMan = new HistoManager();
        // G4int threadId = G4Threading::G4GetThreadId();
        // G4cout << "--->In ActionInitialization::BuildForMaster(), threadId=" << threadId << G4endl;
        std::string commandShell = "cp rayin.dat rayin_"+std::to_string(fProcessNum)+".dat";
        std::system(commandShell.c_str());
        MDMTrace *mdmTrace = new MDMTrace(fProcessNum);
        RunAction *runAction = new RunAction(histoMan,mdmTrace);
        runAction->SetIsTargetChamber(fIsTargetChamber);
        runAction->SetProcessNum(fProcessNum);
        SetUserAction(runAction);
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    void ActionInitialization::Build() const
    {
        HistoManager *histoMan = new HistoManager();
        // G4int threadId = G4Threading::G4GetThreadId();
        // G4cout << "--->In ActionInitialization::Build(), threadId=" << threadId << G4endl;
        std::string commandShell = "cp rayin.dat rayin_"+std::to_string(fProcessNum)+".dat";
        std::system(commandShell.c_str());
        MDMTrace *mdmTrace = new MDMTrace(fProcessNum);
        RunAction *runAction = new RunAction(histoMan,mdmTrace);
        runAction->SetIsTargetChamber(fIsTargetChamber);
        runAction->SetProcessNum(fProcessNum);
        SetUserAction(runAction);

        PrimaryGeneratorAction *pgAction = new PrimaryGeneratorAction(runAction, fIsTargetChamber, fProcessNum);
        pgAction->SetBeamEnergy(fBeamEnergy);
        SetUserAction(pgAction);

        EventAction *eventAction = new EventAction(runAction);
        eventAction->SetIsTargetChamber(fIsTargetChamber);
        eventAction->SetEventActionParams(fEventActionParams);
        SetUserAction(eventAction);

        SetUserAction(new SteppingAction());
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
