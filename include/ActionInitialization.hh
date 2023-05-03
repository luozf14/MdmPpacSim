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
/// \file ActionInitialization.hh
/// \brief Definition of the MdmPpacSim::ActionInitialization class

#ifndef MdmPpacSimActionInitialization_h
#define MdmPpacSimActionInitialization_h 1

#include "G4VUserActionInitialization.hh"
#include "globals.hh"

/// Action initialization class.

namespace MdmPpacSim
{

    class ActionInitialization : public G4VUserActionInitialization
    {
    public:
        ActionInitialization();
        ~ActionInitialization() override;

        void BuildForMaster() const override;
        void Build() const override;
        void SetEventActionParams(std::map<std::string, G4double> params) { fEventActionParams = params; }
        void SetIsTargetChamber(G4bool isTargetChamber) { fIsTargetChamber = isTargetChamber; }
        void SetBeamEnergy(G4double e) { fBeamEnergy = e; }
        G4double GetBeamEnergy() { return fBeamEnergy; }
        void SetProcessNum(G4int proc) { fProcessNum = proc; }
        G4int GetProcessNum() { return fProcessNum; }

    private:
        G4int fProcessNum;
        G4bool fIsTargetChamber;
        G4double fBeamEnergy;
        std::map<std::string, G4double> fEventActionParams;
    };

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
