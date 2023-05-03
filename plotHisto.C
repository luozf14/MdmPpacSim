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

Double_t ma = 4. * 931.49 + 2.42491587;
Double_t mA = 12. * 931.49;
Double_t mb = 4. * 931.49 + 2.42491587;
Double_t mB = 12 * 931.49;
Double_t mHe = 4. * 931.49 + 2.42491587;
Double_t mC12 = 12. * 931.49;
Double_t Ea = 39.972; // 40MeV alpha after 0.00035mm carbon
TVector3 beamMomentum(0., 0., std::sqrt(2. * ma * Ea));
Double_t siAngle = 81.15;
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

    TH1D *h1Q = new TH1D("h1Q", "Excitation energy (singles)", 512, -10, 15.);
    h1Q->SetXTitle("Energy [MeV]");

    TH1D *h1Q_C12 = new TH1D("h1Q_C12", "Excitation energy (C12)", 512, -10, 15.);
    h1Q_C12->SetXTitle("Energy [MeV]");

    TH1D *h1Q_Alpha = new TH1D("h1Q_Alpha", "Excitation energy (Alpha)", 512, -10, 15.);
    h1Q_Alpha->SetXTitle("Energy [MeV]");

    TH1D *h1ToF = new TH1D("h1ToF", "timeMDM - timeE", 512, 0., 100.);
    h1ToF->SetXTitle("ToF [ns]");

    TH2D *h2ToFvsQ = new TH2D("h2ToFvsQ", "ToF vs Excitation energy", 512, -10, 15., 512, 0., 100.);
    h2ToFvsQ->SetYTitle("Energy [MeV]");
    h2ToFvsQ->SetXTitle("ToF [ns]");

    TH2D *h2X2vsX1 = new TH2D("h2X2vsX1", "X2 vs X1", 512, -40., 40., 512, -40., 40.);
    h2X2vsX1->SetYTitle("X2 [cm]");
    h2X2vsX1->SetXTitle("X1 [cm]");

    // TH2D *h2PrimaryParticleTraceSpaceX = new TH2D("h2PrimaryParticleTraceSpaceX", "Primary particle X trace space", 512, -10., 10., 512, -10., 10.);
    // h2PrimaryParticleTraceSpaceX->SetXTitle("X [mm]");
    // h2PrimaryParticleTraceSpaceX->SetYTitle("X' [mrad]");

    // TH2D *h2PrimaryParticleTraceSpaceY = new TH2D("h2PrimaryParticleTraceSpaceY", "Primary particle Y trace space", 512, -10., 10., 512, -10., 10.);
    // h2PrimaryParticleTraceSpaceY->SetXTitle("Y [mm]");
    // h2PrimaryParticleTraceSpaceY->SetYTitle("Y' [mrad]");

    Int_t nEntries = (Int_t)simData1->GetEntries();
    Int_t acc = 0;
    Int_t tra = 0;

    for (int i = 0; i < nEntries; i++)
    {
        simData1->GetEntry(i);
        simData2->GetEntry(i);

        // h2PrimaryParticleTraceSpaceX->Fill(primaryParticlePositionX, primaryParticleMomentumDirectionX / primaryParticleMomentumDirectionZ * 1000.);
        // h2PrimaryParticleTraceSpaceY->Fill(primaryParticlePositionY, primaryParticleMomentumDirectionY / primaryParticleMomentumDirectionZ * 1000.);

        if (timeE < 1.)
            continue;

        if (!(frontNo >= 3 && frontNo <= 9 && backNo >= 6 && backNo <= 9))
            continue;

        // Si momentum
        Double_t xOnDetector = -(((double)frontNo - 8.) * (L / 16.) + (L / 32.));
        Double_t yOnDetector = ((double)backNo - 8.) * (L / 16.) + (L / 32.);
        Double_t xInLab = rE * TMath::Sin(siAngle / 180. * TMath::Pi()) + xOnDetector * TMath::Cos(siAngle / 180. * TMath::Pi());
        Double_t yInLab = yOnDetector;
        Double_t zInLab = rE * TMath::Cos(siAngle / 180. * TMath::Pi()) - xOnDetector * TMath::Sin(siAngle / 180. * TMath::Pi());
        TVector3 siMomentum(xInLab, yInLab, zInLab);
        Double_t siAngleX = std::atan(xInLab / zInLab);
        Double_t siAngleY = std::atan(yInLab / zInLab);
        // Double_t tarZ = 0.00035 / std::cos(siAngleX - 22.5 / 180. * 3.14159) * std::cos(siAngleX);
        // Double_t tarX = 0.00035 / std::cos(siAngleX - 22.5 / 180. * 3.14159) * std::sin(siAngleX);
        // Double_t tarY = tarZ * std::tan(siAngleY / 180. * 3.14159);
        // Double_t rangeAlphaInTarget = std::sqrt(tarX * tarX + tarY * tarY + tarZ * tarZ);
        // Double_t siEnergy = energyLossAlphaInCarbon->AddBack(Eb, rangeAlphaInTarget);
        Double_t siEnergy = energyE + energydE + 2.49E-01 * std::pow(energyE + energydE, -7.91E-01);
        siMomentum.SetMag(std::sqrt(2 * mHe * siEnergy));
        Double_t cosThetab = std::cos(siMomentum.Theta());

        //
        // Q value
        Double_t Q = (ma / mB - 1.) * Ea + (mb / mB + 1.) * siEnergy - 2. * TMath::Sqrt(ma * mb * Ea * siEnergy) * cosThetab / mB;
        // Double_t Q = (ma / mB - 1.) * Ea + (mb / mB + 1.) * siEnergy - 2. * TMath::Sqrt(ma * mb * Ea * siEnergy) * std::cos(siTheta/180.*TMath::Pi()) / mB;
        h1Q->Fill(-Q);

        // h2X2vsX1->Fill(X1, X2);
        if (accepted)
        {
            acc++;
            h1ToF->Fill(timeSlitBox - timeE);
            if (transmitted)
            {
                tra++;
                h1Q_C12->Fill(-Q);
                h2ToFvsQ->Fill(-Q, timeSlitBox - timeE);
                h1Q_Alpha->Fill(-Q);
            }
        }
    }
    cout << "transmission efficiency: " << (float)tra / (float)acc << endl;

    TCanvas *c1 = new TCanvas("c1", "c1", 1024, 768);
    c1->Divide(3, 1);
    c1->cd(1);
    h1Q->Draw();
    c1->cd(2);
    h1Q_C12->Draw();
    c1->cd(3);
    h1Q_Alpha->Draw();
    c1->Update();

    TCanvas *c2 = new TCanvas("c2", "c2", 1024, 768);
    c2->cd(1);
    h2ToFvsQ->Draw("colz");
    c2->Update();
    /*
        TCanvas *c5 = new TCanvas("c5", "c5", 1600, 800);
        c5->Divide(2, 1);
        c5->cd(1);
        h2PrimaryParticleTraceSpaceX->Draw("colz");
        c5->cd(2);
        h2PrimaryParticleTraceSpaceY->Draw("colz");
        c5->Update();
        // c5->SaveAs("c5.png");
    */
}