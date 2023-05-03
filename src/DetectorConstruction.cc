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
/// \file DetectorConstruction.cc
/// \brief Implementation of the MdmPpacSim::DetectorConstruction class

#include "DetectorConstruction.hh"
#include "DetectorSD.hh"

#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4UserLimits.hh"

namespace MdmPpacSim
{

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    DetectorConstruction::DetectorConstruction()
    {
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    DetectorConstruction::~DetectorConstruction()
    {
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    G4VPhysicalVolume *DetectorConstruction::Construct()
    {
        // Get nist material manager
        G4NistManager *nist = G4NistManager::Instance();

        // Option to switch on/off checking of volumes overlaps
        //
        G4bool checkOverlaps = false;
        // G4double inch = 2.54 * cm;
        // G4double atm = 101325. * pascal;

        //
        // World
        //
        G4double worldSizeXY = 2 * m;
        G4double worldSizeZ = 2 * m;
        G4Material *worldMat = (fIsTargetChamber) ? nist->FindOrBuildMaterial("G4_Galactic") : new G4Material("Pentane", (0.00385 * fPpacGasPressureInTorr + 1.666e-5) * kg / m3, nist->FindOrBuildMaterial("G4_N-PENTANE"), kStateGas, 300. * kelvin, fPpacGasPressureInTorr * 133.322 * pascal);
        G4Box *solidWorld =
            new G4Box("SolidWorld",                                            // its name
                      0.5 * worldSizeXY, 0.5 * worldSizeXY, 0.5 * worldSizeZ); // its size

        G4LogicalVolume *logicWorld =
            new G4LogicalVolume(solidWorld,    // its solid
                                worldMat,      // its material
                                "LogicWorld"); // its name

        G4VPhysicalVolume *physWorld =
            new G4PVPlacement(0,               // no rotation
                              G4ThreeVector(), // at (0,0,0)
                              logicWorld,      // its logical volume
                              "PhysWorld",     // its name
                              0,               // its mother  volume
                              false,           // no boolean operation
                              0,               // copy number
                              checkOverlaps);  // overlaps checking

        if (fIsTargetChamber)
        {
            //
            // Target
            //
            G4Material *targetMat = new G4Material("Carbon", 2.253 * g / cm3, nist->FindOrBuildMaterial("G4_C"));
            G4RotationMatrix *targetRot = new G4RotationMatrix;
            targetRot->rotateY(-fTargetRotationAngleInDeg * deg);
            G4ThreeVector targetPos = G4ThreeVector(0, 0 * cm, 0 * cm);
            G4double targetRMin = 0;
            G4double targetRMax = 3.5 * mm;
            G4double targetDz = fTargetThicknessInUm * um;

            G4Tubs *solidTarget =
                new G4Tubs("Target",                                                      // its name
                           targetRMin, targetRMax, 0.5 * targetDz, 2. * M_PI, 2. * M_PI); // its size

            G4LogicalVolume *logicTarget =
                new G4LogicalVolume(solidTarget,    // its solid
                                    targetMat,      // its material
                                    "LogicTarget"); // its name

            new G4PVPlacement(G4Transform3D(*targetRot, targetPos), // its position
                              logicTarget,                          // its logical volume
                              "PhysTarget",                         // its name
                              logicWorld,                           // its mother  volume
                              false,                                // no boolean operation
                              0,                                    // copy number
                              checkOverlaps);                       // overlaps checking

            G4double maxStep = 0.1 * targetDz;
            logicTarget->SetUserLimits(new G4UserLimits(maxStep));

            //
            // DeltaE detector
            //
            G4Material *deltaEMat = nist->FindOrBuildMaterial("G4_Si");
            G4ThreeVector deltaEPos = G4ThreeVector(0, 0, fSiDetectorDistanceInCm * cm - 0.5 * cm);
            deltaEPos.setTheta(fSiDetectorAngleInDeg * deg);
            deltaEPos.setPhi(180. * deg);
            G4RotationMatrix *deltaERot = new G4RotationMatrix;
            deltaERot->rotateY(-fSiDetectorAngleInDeg * deg);
            // Box shape
            G4double deltaEX = 5. * cm;
            G4double deltaEY = 5. * cm;
            G4double deltaEZ = 32. * um;

            G4Box *solidDeltaE =
                new G4Box("SolidDeltaE",
                          0.5 * deltaEX, 0.5 * deltaEY, 0.5 * deltaEZ);

            G4LogicalVolume *logicDeltaE =
                new G4LogicalVolume(solidDeltaE,    // its solid
                                    deltaEMat,      // its material
                                    "LogicDeltaE"); // its name

            new G4PVPlacement(G4Transform3D(*deltaERot, deltaEPos), // its position
                              logicDeltaE,                          // its logical volume
                              "PhysDeltaE",                         // its name
                              logicWorld,                           // its mother  volume
                              false,                                // no boolean operation
                              0,                                    // copy number
                              checkOverlaps);                       // overlaps checking

            //
            // E detector
            //
            G4Material *siMat = nist->FindOrBuildMaterial("G4_Si");
            G4ThreeVector siPos = G4ThreeVector(0, 0, fSiDetectorDistanceInCm * cm);
            siPos.setTheta(fSiDetectorAngleInDeg * deg);
            siPos.setPhi(180. * deg);
            G4RotationMatrix *siRot = new G4RotationMatrix;
            siRot->rotateY(-fSiDetectorAngleInDeg * deg);
            // Box shape
            G4double siX = 5. * cm;
            G4double siY = 5. * cm;
            G4double siZ = 500. * um;

            G4Box *solidSi =
                new G4Box("SolidSi",
                          0.5 * siX, 0.5 * siY, 0.5 * siZ);

            G4LogicalVolume *logicSi =
                new G4LogicalVolume(solidSi,    // its solid
                                    siMat,      // its material
                                    "LogicSi"); // its name

            new G4PVPlacement(G4Transform3D(*siRot, siPos), // its position
                              logicSi,                      // its logical volume
                              "PhysSi",                     // its name
                              logicWorld,                   // its mother  volume
                              false,                        // no boolean operation
                              0,                            // copy number
                              checkOverlaps);               // overlaps checking

            //
            // Slit Box
            //
            G4double slitBoxDistance = 63.5 * cm;
            G4Material *slitBoxMat = nist->FindOrBuildMaterial("G4_Galactic");
            // G4VSolid* solidSlitBox = new G4Box("SlitBox",slitBoxDistance*std::tan(50e-3),slitBoxDistance*std::tan(40e-3),1.*mm); //8msr
            // G4VSolid* solidSlitBox = new G4Box("SlitBox",slitBoxDistance*std::tan(0.035),slitBoxDistance*std::tan(0.035),1.*mm); //4deg by 4deg
            G4VSolid *solidSlitBox = new G4Box("SolidSlitBox", 2.27965 * cm, 2.27965 * cm, 1. * mm); // measured on 1/26/2022
            G4LogicalVolume *logicSlitBox = new G4LogicalVolume(solidSlitBox, slitBoxMat, "LogicSlitBox");
            G4RotationMatrix *slitBoxRot = new G4RotationMatrix; // rotation of daughter frame
            slitBoxRot->rotateY(fMdmAngleInDeg * deg);
            G4ThreeVector slitBoxPos = G4ThreeVector(slitBoxDistance + 0.5 * mm, 0., 0.); // position in mother frame
            slitBoxPos.setTheta(fMdmAngleInDeg * deg);
            slitBoxPos.setPhi(0. * deg);
            new G4PVPlacement(G4Transform3D(*slitBoxRot, slitBoxPos), logicSlitBox, "PhysSlitBox", logicWorld, false, 0, checkOverlaps);
        }
        else
        {
            //
            // PPAC constants
            //
            G4double ppacWidth = 40. * cm;
            G4double ppacHeight = 10. * cm;
            G4double ppacEntranceWindowThickness = 2.5 * um;
            G4double ppacCathodeAlThickness = 80. * 1e-6 * g / cm2;
            G4double ppacCathodeMylarThickness = 220. * 1e-6 * g / cm2;
            G4double ppacSpacingWindowCathode = 1. * cm;
            G4double ppacPillarWidth = 2. * mm;
            G4double ppacPillarThickness = 1. * mm;
            //
            // PPAC entrance window
            //
            G4VSolid *solidEntranceWindow = new G4Box("EntranceWindowSolid", ppacWidth / 2., ppacHeight / 2., ppacEntranceWindowThickness / 2.);
            G4LogicalVolume *logicEntranceWindow = new G4LogicalVolume(solidEntranceWindow, nist->FindOrBuildMaterial("G4_MYLAR"), "LogicPpacEntranceWindow");
            G4PVPlacement *physEntranceWindow = new G4PVPlacement(nullptr, G4ThreeVector(0., 0., 0.5 * ppacEntranceWindowThickness), logicEntranceWindow,
                                                                  "PhysPpacEntranceWindow", logicWorld, false, 0, checkOverlaps);

            //
            // PPAC cathode
            //
            // Al layer
            G4Material *cathodeAlMaterial = new G4Material("Al", 2.70 * g / cm3, nist->FindOrBuildMaterial("G4_Al"));
            G4double cathodeAlDz = ppacCathodeAlThickness / (2.70 * g / cm3);
            G4VSolid *solidCathodeAl = new G4Box("SolidCathodeAl", ppacWidth / 2., ppacHeight / 2., cathodeAlDz / 2.);
            G4LogicalVolume *logicCathodeAl = new G4LogicalVolume(solidCathodeAl, cathodeAlMaterial, "LogicCathodeAl");

            // cathode mylar layer
            G4Material *cathodeMylarMaterial = nist->FindOrBuildMaterial("G4_MYLAR");
            G4double cathodeMylarDz = ppacCathodeMylarThickness / cathodeMylarMaterial->GetDensity();
            G4cout << "cathodeMylarDz=" << G4BestUnit(cathodeMylarDz, "Length") << G4endl;
            G4VSolid *solidCathode0Mylar = new G4Box("SolidCathode0Mylar", ppacWidth / 2., ppacHeight / 2., cathodeMylarDz / 2.);
            G4LogicalVolume *logicCathode0Mylar = new G4LogicalVolume(solidCathode0Mylar, cathodeMylarMaterial, "LogicCathode0Mylar");
            G4VSolid *solidCathode1Mylar = new G4Box("SolidCathode1Mylar", ppacWidth / 2., ppacHeight / 2., cathodeMylarDz / 2.);
            G4LogicalVolume *logicCathode1Mylar = new G4LogicalVolume(solidCathode1Mylar, cathodeMylarMaterial, "LogicCathode1Mylar");

            //
            //  place PPAC cathode 1 (Al_0 + Mylar_0 + Al_1)
            //
            G4ThreeVector cathode0MylarPos(0., 0., physEntranceWindow->GetObjectTranslation().z() + ppacSpacingWindowCathode);
            new G4PVPlacement(nullptr, cathode0MylarPos, logicCathode0Mylar,
                              "PhysCathode0Mylar",
                              logicWorld, false, 0, checkOverlaps);
            new G4PVPlacement(nullptr, G4ThreeVector(0., 0., cathode0MylarPos.z() - 0.5 * (cathodeAlDz + cathodeMylarDz)), logicCathodeAl,
                              "PhysCathodeAl",
                              logicWorld, false, 0, checkOverlaps);
            new G4PVPlacement(nullptr, G4ThreeVector(0., 0., cathode0MylarPos.z() + 0.5 * (cathodeAlDz + cathodeMylarDz)), logicCathodeAl,
                              "PhysCathodeAl",
                              logicWorld, false, 1, checkOverlaps);

            //
            //  place PPAC cathode 2 (Al_2 + Mylar_1 + Al_3)
            //
            G4ThreeVector cathode1MylarPos(0., 0., cathode0MylarPos.z() + fPpacLengthInCm * cm);
            new G4PVPlacement(nullptr, cathode1MylarPos, logicCathode1Mylar,
                              "PhysCathode1Mylar",
                              logicWorld, false, 0, checkOverlaps);
            new G4PVPlacement(nullptr, G4ThreeVector(0., 0., cathode1MylarPos.z() - 0.5 * (cathodeAlDz + cathodeMylarDz)), logicCathodeAl,
                              "PhysCathodeAl",
                              logicWorld, false, 2, checkOverlaps);
            new G4PVPlacement(nullptr, G4ThreeVector(0., 0., cathode1MylarPos.z() + 0.5 * (cathodeAlDz + cathodeMylarDz)), logicCathodeAl,
                              "PhysCathodeAl",
                              logicWorld, false, 3, checkOverlaps);

            //
            // Pillars
            //
            G4Material *pillarMaterial = nist->FindOrBuildMaterial("G4_POLYCARBONATE");
            G4VSolid *solidPillar = new G4Box("SolidPillar", ppacPillarWidth / 2., ppacHeight / 2., ppacPillarThickness / 2.);
            G4LogicalVolume *logicPillar = new G4LogicalVolume(solidPillar, pillarMaterial, "LogicPillar");

            // place pillars of PPAC 1
            G4ThreeVector pillar0Pos(0., 0., cathode0MylarPos.z() - 0.5 * (cathodeAlDz + cathodeMylarDz));
            new G4PVPlacement(nullptr, G4ThreeVector(pillar0Pos.x() + ppacWidth / 4., pillar0Pos.y(), pillar0Pos.z()), logicPillar,
                              "PhysPillar0", logicWorld, false, 0, checkOverlaps);
            new G4PVPlacement(nullptr, pillar0Pos, logicPillar,
                              "PhysPillar0", logicWorld, false, 1, checkOverlaps);
            new G4PVPlacement(nullptr, G4ThreeVector(pillar0Pos.x() - ppacWidth / 4., pillar0Pos.y(), pillar0Pos.z()), logicPillar,
                              "PhysPillar0", logicWorld, false, 2, checkOverlaps);

            // place pillars of PPAC 2
            G4ThreeVector pillar1Pos(0., 0., cathode1MylarPos.z() - 0.5 * (cathodeAlDz + cathodeMylarDz));
            new G4PVPlacement(nullptr, G4ThreeVector(pillar1Pos.x() + ppacWidth / 4., pillar1Pos.y(), pillar1Pos.z()), logicPillar,
                              "PhysPillar1", logicWorld, false, 0, checkOverlaps);
            new G4PVPlacement(nullptr, pillar1Pos, logicPillar,
                              "PhysPillar1", logicWorld, false, 1, checkOverlaps);
            new G4PVPlacement(nullptr, G4ThreeVector(pillar1Pos.x() - ppacWidth / 4., pillar1Pos.y(), pillar1Pos.z()), logicPillar,
                              "PhysPillar1", logicWorld, false, 2, checkOverlaps);
        }

        //
        // always return the physical World
        //
        return physWorld;
    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    void DetectorConstruction::ConstructSDandField()
    {
        // Sensitive detector
        auto *sdManager = G4SDManager::GetSDMpointer();
        if (fIsTargetChamber)
        {
            DetectorSD *aDeltaESD = new DetectorSD("/DeltaESD", "DeltaEHitsCollection");
            sdManager->AddNewDetector(aDeltaESD);
            SetSensitiveDetector("LogicDeltaE", aDeltaESD, true);
            DetectorSD *aSiSD = new DetectorSD("/SiSD", "SiHitsCollection");
            sdManager->AddNewDetector(aSiSD);
            SetSensitiveDetector("LogicSi", aSiSD, true);
            DetectorSD *aSlitBoxSD = new DetectorSD("/SlitBoxSD", "SlitBoxHitsCollection");
            sdManager->AddNewDetector(aSlitBoxSD);
            SetSensitiveDetector("LogicSlitBox", aSlitBoxSD, true);
        }
        else
        {
            DetectorSD *aPpac1SD = new DetectorSD("/Ppac1SD", "Ppac1HitsCollection");
            sdManager->AddNewDetector(aPpac1SD);
            SetSensitiveDetector("LogicCathode0Mylar", aPpac1SD, true);
            DetectorSD *aPpac2SD = new DetectorSD("/Ppac2SD", "Ppac2HitsCollection");
            sdManager->AddNewDetector(aPpac2SD);
            SetSensitiveDetector("LogicCathode1Mylar", aPpac2SD, true);
        }
    }

    void DetectorConstruction::SetDetectorParams(std::map<std::string, G4double> detectorParams)
    {
        for (auto it : detectorParams)
        {
            if (it.first == "TargetThicknessInUm")
            {
                fTargetThicknessInUm = it.second;
            }
            else if (it.first == "TargetRotationAngleInDeg")
            {
                fTargetRotationAngleInDeg = it.second;
            }
            else if (it.first == "SiDetectorAngleInDeg")
            {
                fSiDetectorAngleInDeg = it.second;
            }
            else if (it.first == "SiDetectorDistanceInCm")
            {
                fSiDetectorDistanceInCm = it.second;
            }
            else if (it.first == "MdmAngleInDeg")
            {
                fMdmAngleInDeg = it.second;
            }
            else if (it.first == "PpacGasPressureInTorr")
            {
                fPpacGasPressureInTorr = it.second;
            }
            else if (it.first == "PpacLengthInCm")
            {
                fPpacLengthInCm = it.second;
            }
        }
    }

}
