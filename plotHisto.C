#include <iostream>
#include <vector>
#include <memory>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <TLegend.h>

// constants
Double_t Ea = 26.;
Double_t Vcb = 8.6;                      // coulomb barrier of 13C+12C
Double_t ma = 13. * 931.49 + 3.12500933; // mass of 13C
Double_t mA = 12. * 931.49;              // mass of 12C
Double_t mx = 12. * 931.49;              // mass of 12C
Double_t ms = 931.49 + 8.07131806;       // mass of neutron
Double_t mu_sx = ms * mx / (ms + mx);
// Double_t mb = 931.49 + 7.288971064;      // mass of proton
// Double_t mB = 23. * 931.49 - 9.52985352; // mass of 23Na
Double_t mb = 4. * 931.49 + 2.42491587;  // mass of 4He
Double_t mB = 20. * 931.49 - 7.04193217; // mass of 20Ne
Double_t mF = mx + mA;                   // mass of 24Mg
Double_t mu_sF = ms * mF / (ms + mF);
Double_t e_sx = mx + ms - ma;
Double_t ksx = std::sqrt(2. * mu_sx * e_sx);
Double_t siAngle = 65;
Double_t rE = 14.224; // distance between center of E detector and target
Double_t L = 5.;      // length of E detector

void plotHisto()
{
    TFile *file1 = new TFile("Stage1~0.root", "READ");
    TFile *file2 = new TFile("Stage2~0.root", "READ");

    TTree *simData1 = (TTree *)file1->Get("Stage1Data");
    TTree *simData2 = (TTree *)file2->Get("Stage2Data");

    Double_t energydE, timedE;
    Double_t energyE, timeE;
    Double_t siHitLocalPositionX, siHitLocalPositionY, siHitLocalPositionZ;
    Int_t frontNo, backNo;

    Int_t accepted, transmitted;
    Double_t energySlitBox, timeSlitBox;
    Double_t slitHitLocalPositionX, slitHitLocalPositionY, slitHitLocalPositionZ;
    Int_t charge, mass;

    Int_t completed;
    Double_t X1, Y1, X2, Y2;
    Double_t tof;

    simData1->SetBranchAddress("DeltaEEdep", &energydE);
    simData1->SetBranchAddress("DeltETime", &timedE);
    simData1->SetBranchAddress("SiEdep", &energyE);
    simData1->SetBranchAddress("SiTime", &timeE);
    simData1->SetBranchAddress("SiHitLocalPosition.x", &siHitLocalPositionX);
    simData1->SetBranchAddress("SiHitLocalPosition.y", &siHitLocalPositionY);
    simData1->SetBranchAddress("SiHitLocalPosition.z", &siHitLocalPositionZ);
    simData1->SetBranchAddress("SiFrontStripNo", &frontNo);
    simData1->SetBranchAddress("SiBackStripNo", &backNo);

    simData1->SetBranchAddress("SlitBoxAccepted", &accepted);
    simData1->SetBranchAddress("SlitBoxTransmitted", &transmitted);
    simData1->SetBranchAddress("SlitBoxEnergy", &energySlitBox);
    simData1->SetBranchAddress("SlitBoxTime", &timeSlitBox);
    simData1->SetBranchAddress("SlitBoxCharge", &charge);
    simData1->SetBranchAddress("SlitBoxMass", &mass);
    simData1->SetBranchAddress("SlitBoxHitLocalPosition.x", &slitHitLocalPositionX);
    simData1->SetBranchAddress("SlitBoxHitLocalPosition.y", &slitHitLocalPositionY);
    simData1->SetBranchAddress("SlitBoxHitLocalPosition.z", &slitHitLocalPositionZ);

    simData2->SetBranchAddress("Completed", &completed);
    simData2->SetBranchAddress("PPAC1PositionX", &X1);
    simData2->SetBranchAddress("PPAC1PositionY", &Y1);
    simData2->SetBranchAddress("PPAC2PositionX", &X2);
    simData2->SetBranchAddress("PPAC2PositionY", &Y2);
    simData2->SetBranchAddress("PPACTof", &tof);

    TH2D *h2DeltaEVsE = new TH2D("h2DeltaEVsE", "#DeltaE vs E", 100, 0, 20, 100, 0, 10);
    TH1D *h1E = new TH1D("h1E", "Si energy", 100, 0, 10);
    h1E->SetXTitle("Energy [MeV]");
    TH1D *h1SlitBoxE_light = new TH1D("h1SlitBoxE_light", "Slit box energy (light recoil)", 100, 0, 10);
    h1SlitBoxE_light->SetXTitle("Energy [MeV/u]");
    TH1D *h1SlitBoxE_heavy = new TH1D("h1SlitBoxE_heavy", "Slit box energy (heavy recoil)", 100, 0, 2);
    h1SlitBoxE_heavy->SetXTitle("Energy [MeV/u]");
    TH1D *h1SlitBoxE_heavy_trans = new TH1D("h1SlitBoxE_heavy_trans", "Slit box energy (heavy recoil, transmitted)", 100, 0, 2);
    h1SlitBoxE_heavy_trans->SetXTitle("Energy [MeV/u]");

    TH2D *h2TofVsX1 = new TH2D("h2TofVsX1", "Tof vs X1", 100, -20, 20, 100, 0, 30);
    h2TofVsX1->SetXTitle("X1 (cm)");
    h2TofVsX1->SetYTitle("Tof (ns)");

    Int_t nEntries = (Int_t)simData1->GetEntries();
    Int_t tra = 0;
    Int_t com = 0;
    for (int i = 0; i < nEntries; i++)
    {
        simData1->GetEntry(i);
        simData2->GetEntry(i);

        // if (!(frontNo >= 3 && frontNo <= 9 && backNo >= 6 && backNo <= 9))
        //     continue;

        if (accepted)
        {
            if (timeE > 1.)
            {
                h1E->Fill(energyE);
            }
            if (mass == 4 || mass == 1)
                h1SlitBoxE_light->Fill(energySlitBox / (double)mass);
            if (mass == 20 || mass == 23)
                h1SlitBoxE_heavy->Fill(energySlitBox / (double)mass);
            if (transmitted)
            {
                tra++;
                h1SlitBoxE_heavy_trans->Fill(energySlitBox / (double)mass);
                if (completed)
                {
                    com++;
                    h2TofVsX1->Fill(X1,tof);
                }
            }
        }
    }
    std::cout << "tra: " << tra << ", com: " << com << std::endl;

    TCanvas *c1 = new TCanvas("c1", "c1", 1024, 768);
    c1->cd();
    h2TofVsX1->Draw("colz");

    // TCanvas *c2 = new TCanvas("c2", "c2", 1024, 768);
    // c2->Divide(3, 1);
    // c2->cd(1);
    // h1SlitBoxE_light->Draw();
    // c2->cd(2);
    // h1SlitBoxE_heavy->Draw();
    // c2->cd(3);
    // h1SlitBoxE_heavy_trans->Draw();
}