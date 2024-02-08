#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"
#include "ReactionPhysics.hh"
#include "nlohmann/json.hpp"

// #ifdef G4MULTITHREADED
// #include "G4MTRunManager.hh"
// #else
// #include "G4RunManager.hh"
// #endif

#include "G4RunManagerFactory.hh"

#include "G4SteppingVerbose.hh"
#include "G4UImanager.hh"
#include "QBBC.hh"
#include "G4EmStandardPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4DecayPhysics.hh"
#include "G4StoppingPhysics.hh"
#include "G4StepLimiterPhysics.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "Randomize.hh"

#include "TGraph.h"

#include <fstream>
#include <iostream>

using namespace MdmPpacSim;
using json = nlohmann::json;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc, char **argv)
{
    // Detect interactive mode (if no arguments) and define UI session
    //
    if (argc != 2)
    {
        G4cerr << "--->Error: wrong input parameters!"
               << "\n--->Usage: ./MdmPpacSim <config.json>" << G4endl;
        return 0;
    }

    //
    // Parse the configure JSON file
    //
    std::string configFile = argv[1];
    G4cout << configFile << G4endl;
    std::ifstream configStream(configFile.c_str());
    json config = json::parse(configStream);
    // main level
    G4bool useGUI = config["GUI"].get<G4bool>();
    G4String runMac = config["RunMacro"].get<std::string>();
    G4int processNum = config["ProcessNumber"].get<G4int>();
    // G4int processNum = std::stoi(argv[2]);
    G4bool isTargetChamber = config["IsTargetChamber"].get<G4bool>();
    // Beam
    G4double beamEnergy = config["BeamEnergy"].get<G4double>();
    // Detector
    G4double targetThicknessInUm = config["TargetThicknessInUm"].get<G4double>();
    G4double targetRotationAngleInDeg = config["TargetRotationAngleInDeg"].get<G4double>();
    G4double siDetectorAngleInDeg = config["SiDetectorAngleInDeg"].get<G4double>();
    G4double siDetectorDistanceInCm = config["SiDetectorDistanceInCm"].get<G4double>();
    G4double mdmAngleInDeg = config["MdmAngleInDeg"].get<G4double>();
    G4double ppacGasPressureInTorr = config["PpacGasPressureInTorr"].get<G4double>();
    G4double ppacLengthInCm = config["PpacLengthInCm"].get<G4double>();
    // EventAction
    G4double siDetectorEnergyResolution = config["SiDetectorEnergyResolution"].get<G4double>();
    G4double tdcResolutionInNs = config["TdcResolutionInNs"].get<G4double>();
    G4double ppacPositionResolutionInMm = config["PpacPositionResolutionInMm"].get<G4double>();
    G4double mdmMultipoleProbe = config["MdmMultipoleProbe"].get<G4double>();
    G4double mdmDipoleProbe = config["MdmDipoleProbe"].get<G4double>();
    // Reaction
    G4int targetCharge = config["Target"][0].get<G4int>();
    G4int targetMass = config["Target"][1].get<G4int>();
    G4int lightProductCharge = config["LightProduct"][0].get<G4int>();
    G4int lightProductMass = config["LightProduct"][1].get<G4int>();
    G4int heavyProductCharge = config["HeavyProduct"][0].get<G4int>();
    G4int heavyProductMass = config["HeavyProduct"][1].get<G4int>();
    // Wave function
    G4String waveFuncFile = config["WaveFunction"].get<std::string>();

    // Read bound state wave function psi_p^2
    std::ifstream dataStream(waveFuncFile.c_str());
    G4double p, psi_p2;
    G4int nPoints = 0;
    TGraph *grPsi_p2 = new TGraph();
    grPsi_p2->SetName("waveFunction");
    while (!dataStream.eof())
    {
        dataStream >> p >> psi_p2;
        grPsi_p2->SetPoint(nPoints, p, psi_p2);
        nPoints++;
    }

    // Parameters for classes
    std::map<std::string, G4double> detectorParams; // parameters for DectectorConstruction
    detectorParams["TargetThicknessInUm"] = targetThicknessInUm;
    detectorParams["TargetRotationAngleInDeg"] = targetRotationAngleInDeg;
    detectorParams["SiDetectorAngleInDeg"] = siDetectorAngleInDeg;
    detectorParams["SiDetectorDistanceInCm"] = siDetectorDistanceInCm;
    detectorParams["MdmAngleInDeg"] = mdmAngleInDeg;
    detectorParams["PpacGasPressureInTorr"] = ppacGasPressureInTorr;
    detectorParams["PpacLengthInCm"] = ppacLengthInCm;
    std::map<std::string, G4double> eventActionParams; // parameters for EventAction
    eventActionParams["SiDetectorEnergyResolution"] = siDetectorEnergyResolution;
    eventActionParams["TdcResolutionInNs"] = tdcResolutionInNs;
    eventActionParams["PpacPositionResolutionInMm"] = ppacPositionResolutionInMm;
    eventActionParams["MdmMultipoleProbe"] = mdmMultipoleProbe;
    eventActionParams["MdmDipoleProbe"] = mdmDipoleProbe;
    eventActionParams["MdmAngleInDeg"] = mdmAngleInDeg;
    std::map<std::string, G4int> reactionParams; // parameters for ReactionPhysics
    reactionParams["TargetCharge"] = targetCharge;
    reactionParams["TargetMass"] = targetMass;
    reactionParams["LightProductCharge"] = lightProductCharge;
    reactionParams["LightProductMass"] = lightProductMass;
    reactionParams["HeavyProductCharge"] = heavyProductCharge;
    reactionParams["HeavyProductMass"] = heavyProductMass;

    // Optionally: choose a different Random engine...
    G4Random::setTheEngine(new CLHEP::MTwistEngine);
    // G4Random::setTheEngine(new CLHEP::RanecuEngine);
    // set random seed with system time
    G4long seed = time(NULL);
    CLHEP::HepRandom::setTheSeed(seed);

    // Construct the default run manager
    //
    // #ifdef G4MULTITHREADED
    //     G4MTRunManager *runManager = new G4MTRunManager;
    // #else
    // G4RunManager *runManager = new G4RunManager;
    // #endif
    G4RunManagerType runManagerType = isTargetChamber ? G4RunManagerType::Serial : G4RunManagerType::Serial;
    auto *runManager = G4RunManagerFactory::CreateRunManager(runManagerType, 2);

    // Set mandatory initialization classes
    //
    // Detector construction
    DetectorConstruction *detector = new DetectorConstruction();
    detector->SetIsTargetChamber(isTargetChamber);
    detector->SetDetectorParams(detectorParams);
    runManager->SetUserInitialization(detector);

    // Physics list
    // G4VModularPhysicsList *physicsList = new QBBC;
    G4VModularPhysicsList *physicsList = new G4VModularPhysicsList(); // build user defined physcis list
    physicsList->RegisterPhysics(new G4EmStandardPhysics());
    physicsList->RegisterPhysics(new G4RadioactiveDecayPhysics());
    physicsList->RegisterPhysics(new G4DecayPhysics());
    physicsList->RegisterPhysics(new G4StoppingPhysics());
    physicsList->RegisterPhysics(new G4StepLimiterPhysics());
    ReactionPhysics *reactionPhysics = new ReactionPhysics;
    reactionPhysics->SetReactionParams(reactionParams);
    reactionPhysics->SetWaveFunction(grPsi_p2);
    physicsList->RegisterPhysics(reactionPhysics);
    physicsList->SetVerboseLevel(0);
    runManager->SetUserInitialization(physicsList);

    // User action initialization
    ActionInitialization *actionInit = new ActionInitialization();
    actionInit->SetIsTargetChamber(isTargetChamber);
    actionInit->SetBeamEnergy(beamEnergy);
    actionInit->SetProcessNum(processNum);
    actionInit->SetEventActionParams(eventActionParams);
    runManager->SetUserInitialization(actionInit);

    // Initialize visualization
    //
    G4VisManager *visManager = new G4VisExecutive;
    // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
    // G4VisManager* visManager = new G4VisExecutive("Quiet");
    visManager->Initialize();

    // Get the pointer to the User Interface manager
    G4UImanager *UImanager = G4UImanager::GetUIpointer();

    // Process macro or start UI session
    //
    if (!useGUI)
    {
        // batch mode
        G4String command = "/control/execute ";
        G4String fileName = runMac;
        UImanager->ApplyCommand(command + fileName);
    }
    else
    {
        // interactive mode
        G4UIExecutive *ui = new G4UIExecutive(argc, argv);
        UImanager->ApplyCommand("/control/execute init_vis.mac");
        ui->SessionStart();
        delete ui;
    }

    // Job termination
    // Free the store: user actions, physics_list and detector_description are
    // owned and deleted by the run manager, so they should not be deleted
    // in the main() program !

    delete visManager;
    delete runManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
