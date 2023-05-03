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
/// \file RunAction.hh
/// \brief Definition of the MdmPpacSim::RunAction class

#ifndef MdmPpacSimRunAction_h
#define MdmPpacSimRunAction_h 1

#include "G4UserRunAction.hh"
#include "G4Accumulable.hh"
#include "globals.hh"

class G4Run;

/// Run action class
///
/// In EndOfRunAction(), it calculates the dose in the selected volume
/// from the energy deposit accumulated via stepping and event actions.
/// The computed dose is then printed on the screen.

namespace MdmPpacSim
{
    class HistoManager;

    class RunAction : public G4UserRunAction
    {
    public:
        RunAction(HistoManager *);
        ~RunAction() override;

        void BeginOfRunAction(const G4Run *) override;
        void EndOfRunAction(const G4Run *) override;

        HistoManager *GetHistoManager() { return fHistoManager; }
        void SetIsTargetChamber(G4bool isTargetChamber) { fIsTargetChamber = isTargetChamber; }
        void SetProcessNum(G4int num) { fProcessNum = num; }

    private:
        // G4Accumulable<G4double> fEdep;
        G4bool fIsTargetChamber;
        HistoManager *fHistoManager = nullptr;
        G4int fProcessNum;
    };

}

#endif
