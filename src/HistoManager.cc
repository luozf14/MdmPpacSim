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
/// \file src/HistoManager.cc
/// \brief Implementation of the HistoManager class
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "HistoManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
namespace MdmPpacSim
{

    HistoManager::HistoManager()
        : fFactoryOn(false)
    {
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    HistoManager::~HistoManager()
    {
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    void HistoManager::Book(G4int processNum)
    {
        // Create or get analysis manager
        // The choice of analysis technology is done via selection of a namespace
        // in HistoManager.hh
        G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();

        if (!fFactoryOn)
        {
            //
            analysisManager->SetDefaultFileType("root");
            analysisManager->SetVerboseLevel(1);
            // Only merge in MT mode to avoid warning when running in Sequential mode
            if (G4Threading::IsMultithreadedApplication())
            {
                analysisManager->SetNtupleMerging(true);
            }
            analysisManager->SetCompressionLevel(6);
            // Create directories
            //   analysisManager->SetHistoDirectoryName("histo");
            //   analysisManager->SetNtupleDirectoryName("ntuple");
        }

        // Open an output file
        //
        if (fIsTargetChamber)
        {
            G4bool fileOpen = analysisManager->OpenFile("Stage1_" + std::to_string(processNum));
            if (!fileOpen)
            {
                G4cerr << "\n---> HistoManager::Book(): cannot open "
                       << analysisManager->GetFileName() << G4endl;
                return;
            }

            if (!fFactoryOn)
            {
                // Create ntuples.
                // Ntuples ids are generated automatically starting from 0.
                // The start value can be changed by:
                // analysisManager->SetFirstMtupleId(1);

                // Create 1st ntuple (id = 0)
                analysisManager->CreateNtuple("Stage1Data", "Simulation data of stage 1");

                analysisManager->CreateNtupleDColumn("DeltaEEdep");  // column Id = 0
                analysisManager->CreateNtupleDColumn("DeltETime");   // column Id = 1
                analysisManager->CreateNtupleIColumn("DeltECharge"); // column Id = 2
                analysisManager->CreateNtupleIColumn("DeltEMass");   // column Id = 3

                analysisManager->CreateNtupleDColumn("SiEdep");               // column Id = 4
                analysisManager->CreateNtupleDColumn("SiTime");               // column Id = 5
                analysisManager->CreateNtupleIColumn("SiCharge");             // column Id = 6
                analysisManager->CreateNtupleIColumn("SiMass");               // column Id = 7
                analysisManager->CreateNtupleDColumn("SiHitPosition.x");      // column Id = 8
                analysisManager->CreateNtupleDColumn("SiHitPosition.y");      // column Id = 9
                analysisManager->CreateNtupleDColumn("SiHitPosition.z");      // column Id = 10
                analysisManager->CreateNtupleDColumn("SiHitLocalPosition.x"); // column Id = 11
                analysisManager->CreateNtupleDColumn("SiHitLocalPosition.y"); // column Id = 12
                analysisManager->CreateNtupleDColumn("SiHitLocalPosition.z"); // column Id = 13
                analysisManager->CreateNtupleIColumn("SiFrontStripNo");       // column Id = 14
                analysisManager->CreateNtupleIColumn("SiBackStripNo");        // column Id = 15

                analysisManager->CreateNtupleIColumn("SlitBoxAccepted");           // column Id = 16
                analysisManager->CreateNtupleIColumn("SlitBoxTransmitted");        // column Id = 17
                analysisManager->CreateNtupleDColumn("SlitBoxHitPosition.x");      // column Id = 18
                analysisManager->CreateNtupleDColumn("SlitBoxHitPosition.y");      // column Id = 19
                analysisManager->CreateNtupleDColumn("SlitBoxHitPosition.z");      // column Id = 20
                analysisManager->CreateNtupleDColumn("SlitBoxHitLocalPosition.x"); // column Id = 21
                analysisManager->CreateNtupleDColumn("SlitBoxHitLocalPosition.y"); // column Id = 22
                analysisManager->CreateNtupleDColumn("SlitBoxHitLocalPosition.z"); // column Id = 23
                analysisManager->CreateNtupleDColumn("SlitBoxHitMomentum.x");      // column Id = 24
                analysisManager->CreateNtupleDColumn("SlitBoxHitMomentum.y");      // column Id = 25
                analysisManager->CreateNtupleDColumn("SlitBoxHitMomentum.z");      // column Id = 26
                analysisManager->CreateNtupleDColumn("SlitBoxHitLocalMomentum.x"); // column Id = 27
                analysisManager->CreateNtupleDColumn("SlitBoxHitLocalMomentum.y"); // column Id = 28
                analysisManager->CreateNtupleDColumn("SlitBoxHitLocalMomentum.z"); // column Id = 29
                analysisManager->CreateNtupleIColumn("SlitBoxCharge");             // column Id = 30
                analysisManager->CreateNtupleIColumn("SlitBoxMass");               // column Id = 31
                analysisManager->CreateNtupleDColumn("SlitBoxEnergy");             // column Id = 32
                analysisManager->CreateNtupleDColumn("SlitBoxTime");               // column Id = 33
                analysisManager->CreateNtupleDColumn("ScatteredAngleX");           // column Id = 34
                analysisManager->CreateNtupleDColumn("ScatteredAngleY");           // column Id = 35

                analysisManager->CreateNtupleDColumn("MDMPositionX"); // column Id = 36
                analysisManager->CreateNtupleDColumn("MDMPositionY"); // column Id = 37
                analysisManager->CreateNtupleDColumn("MDMAngleX");    // column Id = 38
                analysisManager->CreateNtupleDColumn("MDMAngleY");    // column Id = 39

                analysisManager->CreateNtupleDColumn("BeamEnergy");              // column Id = 40
                analysisManager->CreateNtupleDColumn("BeamPosition.x");          // column Id = 41
                analysisManager->CreateNtupleDColumn("BeamPosition.y");          // column Id = 42
                analysisManager->CreateNtupleDColumn("BeamPosition.z");          // column Id = 43
                analysisManager->CreateNtupleDColumn("BeamMomentumDirection.x"); // column Id = 44
                analysisManager->CreateNtupleDColumn("BeamMomentumDirection.y"); // column Id = 45
                analysisManager->CreateNtupleDColumn("BeamMomentumDirection.z"); // column Id = 46

                analysisManager->FinishNtuple();

                fFactoryOn = true;
            }
        }
        else
        {
            G4bool fileOpen = analysisManager->OpenFile("Stage2_" + std::to_string(processNum));
            if (!fileOpen)
            {
                G4cerr << "\n---> HistoManager::Book(): cannot open "
                       << analysisManager->GetFileName() << G4endl;
                return;
            }

            // Create 1st ntuple (id = 0)
            analysisManager->CreateNtuple("Stage2Data", "Simulation data of stage 2");

            analysisManager->CreateNtupleDColumn("PPAC1Energy");    // column Id = 0
            analysisManager->CreateNtupleDColumn("PPAC1Time");      // column Id = 1
            analysisManager->CreateNtupleDColumn("PPAC1PositionX"); // column Id = 2
            analysisManager->CreateNtupleDColumn("PPAC1PositionY"); // column Id = 3

            analysisManager->CreateNtupleDColumn("PPAC2Energy");    // column Id = 4
            analysisManager->CreateNtupleDColumn("PPAC2Time");      // column Id = 5
            analysisManager->CreateNtupleDColumn("PPAC2PositionX"); // column Id = 6
            analysisManager->CreateNtupleDColumn("PPAC2PositionY"); // column Id = 7

            analysisManager->CreateNtupleIColumn("Completed"); // column Id = 8
            analysisManager->CreateNtupleDColumn("PPACTof");   // column Id = 9

            analysisManager->FinishNtuple();

            fFactoryOn = true;
        }
        G4cout << "\n----> Output file is open in "
               << analysisManager->GetFileName() << "."
               << analysisManager->GetFileType() << G4endl;
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    void HistoManager::Save()
    {
        if (!fFactoryOn)
            return;

        G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();
        analysisManager->Write();
        analysisManager->CloseFile();

        G4cout << "\n----> Histograms and ntuples are saved\n"
               << G4endl;
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    void HistoManager::FillNtuple()
    {
        G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();

        if (fIsTargetChamber)
        {
            // Fill 1st ntuple ( ntupleId = 0)
            // DeltaE detector
            analysisManager->FillNtupleDColumn(0, DeltaEEdep);
            analysisManager->FillNtupleDColumn(1, DeltaETime);
            analysisManager->FillNtupleIColumn(2, DeltaECharge);
            analysisManager->FillNtupleIColumn(3, DeltaEMass);

            // Si detector
            analysisManager->FillNtupleDColumn(4, SiEdep);
            analysisManager->FillNtupleDColumn(5, SiTime);
            analysisManager->FillNtupleIColumn(6, SiCharge);
            analysisManager->FillNtupleIColumn(7, SiMass);
            analysisManager->FillNtupleDColumn(8, SiHitPosition.x());
            analysisManager->FillNtupleDColumn(9, SiHitPosition.y());
            analysisManager->FillNtupleDColumn(10, SiHitPosition.z());
            analysisManager->FillNtupleDColumn(11, SiHitLocalPosition.x());
            analysisManager->FillNtupleDColumn(12, SiHitLocalPosition.y());
            analysisManager->FillNtupleDColumn(13, SiHitLocalPosition.z());
            analysisManager->FillNtupleIColumn(14, SiFrontStripNo);
            analysisManager->FillNtupleIColumn(15, SiBackStripNo);

            // Slit box
            analysisManager->FillNtupleIColumn(16, SlitBoxAccepted);
            analysisManager->FillNtupleIColumn(17, SlitBoxTransmitted);
            analysisManager->FillNtupleDColumn(18, SlitBoxHitPosition.x());
            analysisManager->FillNtupleDColumn(19, SlitBoxHitPosition.y());
            analysisManager->FillNtupleDColumn(20, SlitBoxHitPosition.z());
            analysisManager->FillNtupleDColumn(21, SlitBoxHitLocalPosition.x());
            analysisManager->FillNtupleDColumn(22, SlitBoxHitLocalPosition.y());
            analysisManager->FillNtupleDColumn(23, SlitBoxHitLocalPosition.z());
            analysisManager->FillNtupleDColumn(24, SlitBoxHitMomentum.x());
            analysisManager->FillNtupleDColumn(25, SlitBoxHitMomentum.y());
            analysisManager->FillNtupleDColumn(26, SlitBoxHitMomentum.z());
            analysisManager->FillNtupleDColumn(27, SlitBoxHitLocalMomentum.x());
            analysisManager->FillNtupleDColumn(28, SlitBoxHitLocalMomentum.y());
            analysisManager->FillNtupleDColumn(29, SlitBoxHitLocalMomentum.z());
            analysisManager->FillNtupleIColumn(30, SlitBoxCharge);
            analysisManager->FillNtupleIColumn(31, SlitBoxMass);
            analysisManager->FillNtupleDColumn(32, SlitBoxEnergy);
            analysisManager->FillNtupleDColumn(33, SlitBoxTime);
            analysisManager->FillNtupleDColumn(34, ScatteredAngleX);
            analysisManager->FillNtupleDColumn(35, ScatteredAngleY);

            // MDMTrace results
            analysisManager->FillNtupleDColumn(36, MDMPositionX);
            analysisManager->FillNtupleDColumn(37, MDMPositionY);
            analysisManager->FillNtupleDColumn(38, MDMAngleX);
            analysisManager->FillNtupleDColumn(39, MDMAngleY);

            // Beam
            analysisManager->FillNtupleDColumn(40, BeamEnergy);
            analysisManager->FillNtupleDColumn(41, BeamPosition.x());
            analysisManager->FillNtupleDColumn(42, BeamPosition.y());
            analysisManager->FillNtupleDColumn(43, BeamPosition.z());
            analysisManager->FillNtupleDColumn(44, BeamMomentumDirection.x());
            analysisManager->FillNtupleDColumn(45, BeamMomentumDirection.y());
            analysisManager->FillNtupleDColumn(46, BeamMomentumDirection.z());

            analysisManager->AddNtupleRow(0);
        }
        else
        {
            // PPAC1
            analysisManager->FillNtupleDColumn(0, PPAC1Energy);
            analysisManager->FillNtupleDColumn(1, PPAC1Time);
            analysisManager->FillNtupleDColumn(2, PPAC1PositionX);
            analysisManager->FillNtupleDColumn(3, PPAC1PositionY);

            // PPAC2
            analysisManager->FillNtupleDColumn(4, PPAC2Energy);
            analysisManager->FillNtupleDColumn(5, PPAC2Time);
            analysisManager->FillNtupleDColumn(6, PPAC2PositionX);
            analysisManager->FillNtupleDColumn(7, PPAC2PositionY);

            // PPAC
            analysisManager->FillNtupleIColumn(8, Completed);
            analysisManager->FillNtupleDColumn(9, PPACTof);

            analysisManager->AddNtupleRow(0);
        }
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......