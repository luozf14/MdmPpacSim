#include "ReactionProcess.hh"
#include "DetectorConstruction.hh"

#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4RunManager.hh"

#include "TLorentzVector.h"
#include "TVector3.h"
namespace MdmPpacSim
{
    ReactionProcess::ReactionProcess(const G4String &processName) : G4VDiscreteProcess(processName, fHadronic)
    {
        SetProcessSubType(111);
    }

    ReactionProcess::~ReactionProcess()
    {
    }

    G4double ReactionProcess::GetMeanFreePath(const G4Track &aTrack, G4double /*previousStepSize*/, G4ForceCondition *condition)
    {

        G4String name = aTrack.GetParticleDefinition()->GetParticleName();
        //   G4cout<<"name: "<<name<<G4endl;
        G4bool test1 = (name == "alpha") ? true : false;
        G4bool test2 = (aTrack.GetNextVolume()->GetLogicalVolume()->GetName() == "LogicTarget") ? true : false;
        G4double mfp = DBL_MAX;
        if (test1 && test2 && aTrack.GetPosition().z() < fReactionPosition && aTrack.GetTrackID() == 1)
        { // for target
            mfp = 0;
        }
        *condition = NotForced;
        return mfp;
    }

    G4VParticleChange *ReactionProcess::PostStepDoIt(const G4Track &aTrack,
                                                     const G4Step &aStep)
    {

        G4StepPoint *postStepPoint = aStep.GetPostStepPoint();
        if (postStepPoint->GetStepStatus() == fGeomBoundary)
        {
            return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
        }

        aParticleChange.Initialize(aTrack);

        //   TwoBody(aTrack,aStep,fTargetCharge, fTargetMass,fLightProductCharge, fLightProductMass, fHeavyProductCharge, fHeavyProductMass, 0.0, 4.44 );

        // G4double ran = G4UniformRand();
        // if (ran < 0.4)
        // {
            TwoBodyDecay(aTrack, aStep, fTargetCharge, fTargetMass, fLightProductCharge, fLightProductMass, 0.0, fHeavyProductCharge, fHeavyProductMass, 9.87, 1);
        // }
        // else
        // {
            // TwoBodyDecay(aTrack, aStep, fTargetCharge, fTargetMass, fLightProductCharge, fLightProductMass, 0.0, fHeavyProductCharge, fHeavyProductMass, 7.65, 2);
        // }
        // else if (ran < 0.6)
        // {
        //     TwoBodyAlphaDecay(aTrack, aStep, 8, 16, 2, 4, 0.0, 8, 16, 10.36);
        // }
        // else if (ran < 0.8)
        // {
        //     TwoBodyNeutronDecay(aTrack, aStep, 6, 13, 2, 4, 0.0, 6, 13, 6.864);
        // }
        // else
        // {
        //     TwoBodyNeutronDecay(aTrack, aStep, 6, 13, 2, 4, 0.0, 6, 13, 7.55);
        // }
        return &aParticleChange;
    }

    void ReactionProcess::StartTracking(G4Track *track)
    {
        G4VProcess::StartTracking(track); // Apply base class actions
        const DetectorConstruction *detectorConstruction = static_cast<const DetectorConstruction *>(G4RunManager::GetRunManager()->GetUserDetectorConstruction());
        fReactionPosition = detectorConstruction->GetTargetThicknessInUm() * um / std::cos(detectorConstruction->GetTargetRotationAngleInDeg() / 180. * M_PI) * G4UniformRand();
    }

    void ReactionProcess::SetParams(std::map<std::string, G4int> params)
    {
        for (auto it : params)
        {
            if (it.first == "TargetCharge")
            {
                fTargetCharge = it.second;
            }
            else if (it.first == "TargetMass")
            {
                fTargetMass = it.second;
            }
            else if (it.first == "LightProductCharge")
            {
                fLightProductCharge = it.second;
            }
            else if (it.first == "LightProductMass")
            {
                fLightProductMass = it.second;
            }
            else if (it.first == "HeavyProductCharge")
            {
                fHeavyProductCharge = it.second;
            }
            else if (it.first == "HeavyProductMass")
            {
                fHeavyProductMass = it.second;
            }
        }
    }

    G4VParticleChange *ReactionProcess::TwoBody(const G4Track &aTrack, const G4Step &aStep, int Zt, int At, int Z1, int A1, int Z2, int A2, double Ex1, double Ex2)
    {

        G4double Mt, M1, M2;

        G4double energy = aTrack.GetDynamicParticle()->GetKineticEnergy() / MeV;
        G4DynamicParticle *target = new G4DynamicParticle;
        G4ParticleDefinition *targetdef;
        targetdef = G4IonTable::GetIonTable()->GetIon(Zt, At, 0.);
        target->SetDefinition(targetdef);
        Mt = targetdef->GetPDGMass() / CLHEP::amu_c2;
        G4DynamicParticle *part1 = new G4DynamicParticle;
        G4ParticleDefinition *part1def;
        if (Z1 == 0)
        {
            G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
            G4String particleName;
            part1def = particleTable->FindParticle(particleName = "neutron");
        }
        if (Z1 != 0)
        {
            part1def = G4IonTable::GetIonTable()->GetIon(Z1, A1, Ex1 * MeV);
        }
        part1->SetDefinition(part1def);
        M1 = part1def->GetPDGMass() / CLHEP::amu_c2;
        G4DynamicParticle *part2 = new G4DynamicParticle;
        G4ParticleDefinition *part2def;
        part2def = G4IonTable::GetIonTable()->GetIon(Z2, A2, Ex2 * MeV);
        part2->SetDefinition(part2def);
        M2 = part2def->GetPDGMass() / CLHEP::amu_c2;
        ////Do some kinematics
        G4double CM_theta = acos(-1. + 2. * G4UniformRand());   // 0->pi relative to the particle direction
        G4double CM_psi = G4UniformRand() * 4. * 2. * atan(1.); // 0->pi
        /// Rotate so theta/psi are relative to the original beam direction
        G4ThreeVector momentumDirection = aTrack.GetMomentumDirection(); // using rotation method from reactionphysics
        G4ThreeVector v = G4ThreeVector(0., 0., 1.).cross(momentumDirection);
        G4double rotAngle = acos(momentumDirection.z());

        G4ThreeVector dir = G4ThreeVector(sin(CM_theta) * sin(CM_psi), sin(CM_theta) * cos(CM_psi), cos(CM_theta));
        if (v.getR() > 0)
            dir.rotate(v, rotAngle); // rotate the direction to be relative to the beam axis

        // Get Q-value
        //	G4cout<<CM_theta<<G4endl;
        G4double Q_value = 0;
        Q_value = aTrack.GetDynamicParticle()->GetDefinition()->GetPDGMass() + targetdef->GetPDGMass() - (part1def->GetPDGMass() + part2def->GetPDGMass());
        // G4cout<<"Effective Q-value= "<<Q_value/MeV<<G4endl;
        // G4cout<<"E_CM = "<<Mt*energy/(M1+M2)<<G4endl;
        G4double E_cm = (Mt * energy / (M1 + M2)) + Q_value; // new CM energy MeV
        if (E_cm < 0.)
            return &aParticleChange;                                               // sub-threshold
        G4double p_1 = sqrt(2. * part1->GetMass() * E_cm * (1. * M2 / (M1 + M2))); // E_1 = m2/(m1+m2) * E_t
        G4double p_2 = sqrt(2. * part2->GetMass() * E_cm * (1. * M1 / (M1 + M2))); // E_2 = m1/(m1+m2) * E_t
        //	G4cout<<"Free momentum: "<<p_1<<"\t"<<p_2<<"\tE_CM\t"<<E_cm<<G4endl;
        G4ThreeVector p_new_1 = p_1 * dir; // new momentum of scattered Be in COM
        G4ThreeVector p_new_2 = -p_new_1;
        G4ThreeVector p_n = aTrack.GetMomentum();
        p_new_1 += p_n * (1. * M1 / (M1 + M2));
        p_new_2 += p_n * (1. * M2 / (M1 + M2));
        //	G4cout<<"Part 1:\t"<<p_new_1<<G4endl;
        //	G4cout<<"Part 2:\t"<<p_new_2<<G4endl;
        //	G4cout<<"Orig:\t"<<p_n<<G4endl;
        part1->SetMomentum(p_new_1);
        part2->SetMomentum(p_new_2);
        G4double total_mom_1 = p_new_1.getR();
        G4double total_mom_2 = p_new_2.getR();
        part1->SetKineticEnergy((total_mom_1 * total_mom_1) / (2. * part1->GetMass()));
        part2->SetKineticEnergy((total_mom_2 * total_mom_2) / (2. * part2->GetMass()));
        G4double lab_theta = p_new_1.theta();
        //	aTrack.GetDynamicParticle()->DumpInfo();
        //	G4cout<<"Alpha E: "<<(total_mom_1*total_mom_1)/(2.*part1->GetMass())<<"\tMom. vector: "<<part1->GetMomentumDirection()<<G4endl;
        //	part1->DumpInfo();
        //	part2->DumpInfo();
        G4Track *sec1 = new G4Track(part1,
                                    aTrack.GetGlobalTime(),
                                    aTrack.GetPosition());

        G4Track *sec2 = new G4Track(part2,
                                    aTrack.GetGlobalTime(),
                                    aTrack.GetPosition());
        // part1->DumpInfo();
        // part2->DumpInfo();

        aParticleChange.AddSecondary(sec1);
        aParticleChange.AddSecondary(sec2);
        aParticleChange.ProposeEnergy(0.);
        aParticleChange.ProposeTrackStatus(fStopAndKill);
        return &aParticleChange;
    }

    G4VParticleChange *ReactionProcess::TwoBodyDecay(const G4Track &aTrack, const G4Step &aStep,
                                                     int Zt, int At,             // target
                                                     int Z1, int A1, double Ex1, // light recoil
                                                     int Z2, int A2, double Ex2, // heavy recoil
                                                     int flag)
    {

        //
        // 3 steps reaction (1  + 2 decays)
        // Step1: a+A->b+B*, Step2: B*->D1+D2, Step3: D2->D21+D22
        //

        // define target
        G4DynamicParticle *target = new G4DynamicParticle;
        G4ParticleDefinition *targetDef;
        targetDef = G4IonTable::GetIonTable()->GetIon(Zt, At, 0.);
        target->SetDefinition(targetDef);
        G4double mA = targetDef->GetPDGMass();

        // define light recoil
        G4DynamicParticle *light = new G4DynamicParticle;
        G4ParticleDefinition *lightDef;
        if (Z1 == 0)
        {
            G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
            lightDef = particleTable->FindParticle("neutron");
        }
        if (Z1 != 0)
        {
            lightDef = G4IonTable::GetIonTable()->GetIon(Z1, A1, Ex1 * MeV);
        }
        light->SetDefinition(lightDef);
        G4double mb = lightDef->GetPDGMass();

        // define heavy recoil
        G4DynamicParticle *heavy = new G4DynamicParticle;
        G4ParticleDefinition *heavyDef;
        heavyDef = G4IonTable::GetIonTable()->GetIon(Z2, A2, Ex2 * MeV);
        heavy->SetDefinition(heavyDef);
        G4double mB = heavyDef->GetPDGMass();

        // define light decay product of heacy recoil
        G4DynamicParticle *d1 = new G4DynamicParticle;
        G4ParticleDefinition *d1Def = G4IonTable::GetIonTable()->GetIon(2, 4, 0.0);
        d1->SetDefinition(d1Def);
        G4double mD1 = d1Def->GetPDGMass();

        // define heavy decay product of heacy recoil
        G4DynamicParticle *d2 = new G4DynamicParticle;
        G4ParticleDefinition *d2Def = G4IonTable::GetIonTable()->GetIon(4, 8, 0.0);
        d2->SetDefinition(d2Def);
        G4double mD2 = d2Def->GetPDGMass();

        // define d21
        G4DynamicParticle *d21 = new G4DynamicParticle;
        G4ParticleDefinition *d21Def = G4IonTable::GetIonTable()->GetIon(2, 4, 0.0);
        d21->SetDefinition(d21Def);
        G4double mD21 = d21Def->GetPDGMass();

        // define d22
        G4DynamicParticle *d22 = new G4DynamicParticle;
        G4ParticleDefinition *d22Def = G4IonTable::GetIonTable()->GetIon(2, 4, 0.0);
        d22->SetDefinition(d22Def);
        G4double mD22 = d22Def->GetPDGMass();

        // get beam properties
        G4LorentzVector tempPa = aTrack.GetDynamicParticle()->Get4Momentum();
        TLorentzVector Pa(tempPa.px(), tempPa.py(), tempPa.pz(), tempPa.e());

        //
        // a+A->b+B*
        //
        // boost to CM of projectile+target (CM1)
        TLorentzVector Pt(0, 0, 0, mA);
        TLorentzVector Pcm1 = Pa + Pt;
        G4double Mcm1 = Pcm1.M();
        // calculate in CM1
        G4double mag_Pcm1 = std::sqrt((std::pow(Mcm1, 2) - std::pow(mB + mb, 2)) * (std::pow(Mcm1, 2) - std::pow(mB - mb, 2))) / 2. / Mcm1;
        G4double E_bcm1 = (Mcm1 * Mcm1 + mb * mb - mB * mB) / 2. / Mcm1;
        G4double E_Bcm1 = (Mcm1 * Mcm1 + mB * mB - mb * mb) / 2. / Mcm1;
        TVector3 vector3Pb(mag_Pcm1, 0, 0);
        vector3Pb.SetTheta(std::acos(-1. + 2. * G4UniformRand()));
        vector3Pb.SetPhi(G4UniformRand() * 2 * M_PI);
        TVector3 vector3PB = (-1.) * vector3Pb;
        TLorentzVector Pb(vector3Pb, E_bcm1);
        TLorentzVector PB(vector3PB, E_Bcm1);
        // boost b and B back to lab
        Pb.Boost(Pcm1.BoostVector());
        PB.Boost(Pcm1.BoostVector());

        //
        // B*->D1+D2
        //
        // calculate in CM of B* (CM2)
        G4double Mcm2 = PB.M();
        G4double mag_Pcm2 = std::sqrt((std::pow(Mcm2, 2) - std::pow(mD1 + mD2, 2)) * (std::pow(Mcm2, 2) - std::pow(mD1 - mD2, 2))) / 2. / Mcm2;
        G4double E_D1cm2 = (Mcm2 * Mcm2 + mD1 * mD1 - mD2 * mD2) / 2. / Mcm2;
        G4double E_D2cm2 = (Mcm2 * Mcm2 + mD2 * mD2 - mD1 * mD1) / 2. / Mcm2;
        TVector3 vector3PD1(mag_Pcm2, 0, 0);
        vector3PD1.SetTheta(std::acos(-1. + 2. * G4UniformRand()));
        vector3PD1.SetPhi(G4UniformRand() * 2 * M_PI);
        TVector3 vector3PD2 = (-1.) * vector3PD1;
        TLorentzVector PD1(vector3PD1, E_D1cm2);
        TLorentzVector PD2(vector3PD2, E_D2cm2);
        // boost PD1 and PD2 back to lab
        PD1.Boost(PB.BoostVector());
        PD2.Boost(PB.BoostVector());

        //
        // D2->D21+D22
        //
        // calculate in CM of D21 (CM3)
        G4double Mcm3 = PD2.M();
        G4double mag_Pcm3 = std::sqrt((std::pow(Mcm3, 2) - std::pow(mD21 + mD22, 2)) * (std::pow(Mcm3, 2) - std::pow(mD21 - mD22, 2))) / 2. / Mcm3;
        G4double E_D21cm3 = (Mcm3 * Mcm3 + mD21 * mD21 - mD22 * mD22) / 2. / Mcm3;
        G4double E_D22cm3 = (Mcm3 * Mcm3 + mD22 * mD22 - mD21 * mD21) / 2. / Mcm3;
        TVector3 vector3PD21(mag_Pcm3, 0, 0);
        vector3PD21.SetTheta(std::acos(-1. + 2. * G4UniformRand()));
        vector3PD21.SetPhi(G4UniformRand() * 2 * M_PI);
        TVector3 vector3PD22 = (-1.) * vector3PD21;
        TLorentzVector PD21(vector3PD21, E_D21cm3);
        TLorentzVector PD22(vector3PD22, E_D22cm3);
        // boost D21 and D22 back to lab
        PD21.Boost(PD2.BoostVector());
        PD22.Boost(PD2.BoostVector());
        // set secondary particles
        light->Set4Momentum(G4LorentzVector(G4ThreeVector(Pb.Px(), Pb.Py(), Pb.Pz()), Pb.E()));
        heavy->Set4Momentum(G4LorentzVector(G4ThreeVector(PB.Px(), PB.Py(), PB.Pz()), PB.E()));
        d1->Set4Momentum(G4LorentzVector(G4ThreeVector(PD1.Px(), PD1.Py(), PD1.Pz()), PD1.E()));
        d21->Set4Momentum(G4LorentzVector(G4ThreeVector(PD21.Px(), PD21.Py(), PD21.Pz()), PD21.E()));
        d22->Set4Momentum(G4LorentzVector(G4ThreeVector(PD22.Px(), PD22.Py(), PD22.Pz()), PD22.E()));

        // set secondary tracks
        G4Track *sec_light = new G4Track(light,
                                         aTrack.GetGlobalTime(),
                                         aTrack.GetPosition());
        G4Track *sec_heavy = new G4Track(heavy,
                                         aTrack.GetGlobalTime(),
                                         aTrack.GetPosition());
        G4Track *sec_d1 = new G4Track(d1,
                                      aTrack.GetGlobalTime(),
                                      aTrack.GetPosition());
        G4Track *sec_d21 = new G4Track(d21,
                                       aTrack.GetGlobalTime(),
                                       aTrack.GetPosition());
        G4Track *sec_d22 = new G4Track(d22,
                                       aTrack.GetGlobalTime(),
                                       aTrack.GetPosition());

        if (flag == 2)
        {
            // add secondary tracks
            aParticleChange.SetNumberOfSecondaries(4);
            aParticleChange.AddSecondary(sec_light);
            // aParticleChange.AddSecondary(sec_heavy);
            aParticleChange.AddSecondary(sec_d1);
            aParticleChange.AddSecondary(sec_d21);
            aParticleChange.AddSecondary(sec_d22);
            // G4cout << "sec_light kinetic energy: " << sec_light->GetKineticEnergy() << G4endl;
            // G4cout << "sec_d1 kinetic energy: " << sec_d1->GetKineticEnergy() << G4endl;
            // G4cout << "sec_d21 kinetic energy: " << sec_d21->GetKineticEnergy() << G4endl;
            // G4cout << "sec_d22 kinetic energy: " << sec_d22->GetKineticEnergy() << G4endl;
        }
        else if (flag == 1)
        {
            // add secondary tracks
            aParticleChange.SetNumberOfSecondaries(2);
            aParticleChange.AddSecondary(sec_light);
            aParticleChange.AddSecondary(sec_heavy);
            // aParticleChange.AddSecondary(sec_d1);
            // aParticleChange.AddSecondary(sec_d21);
            // aParticleChange.AddSecondary(sec_d22);
        }

        aParticleChange.ProposeEnergy(0.);
        aParticleChange.ProposeTrackStatus(fStopAndKill);
        return &aParticleChange;
    }

    G4VParticleChange *ReactionProcess::TwoBodyAlphaDecay(const G4Track &aTrack, const G4Step &aStep,
                                                          int Zt, int At,             // target
                                                          int Z1, int A1, double Ex1, // light recoil
                                                          int Z2, int A2, double Ex2  // heavy recoil
    )
    {

        //
        // 2 steps reaction (1  + 1 decays)
        // Step1: a+A->b+B*, Step2: B*->D1+D2 (D1 is alpha)
        //

        // define target
        G4DynamicParticle *target = new G4DynamicParticle;
        G4ParticleDefinition *targetDef;
        targetDef = G4IonTable::GetIonTable()->GetIon(Zt, At, 0.);
        target->SetDefinition(targetDef);
        G4double mA = targetDef->GetPDGMass();

        // define light recoil
        G4DynamicParticle *light = new G4DynamicParticle;
        G4ParticleDefinition *lightDef;
        if (Z1 == 0)
        {
            G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
            lightDef = particleTable->FindParticle("neutron");
        }
        if (Z1 != 0)
        {
            lightDef = G4IonTable::GetIonTable()->GetIon(Z1, A1, Ex1 * MeV);
        }
        light->SetDefinition(lightDef);
        G4double mb = lightDef->GetPDGMass();

        // define heavy recoil
        G4DynamicParticle *heavy = new G4DynamicParticle;
        G4ParticleDefinition *heavyDef;
        heavyDef = G4IonTable::GetIonTable()->GetIon(Z2, A2, Ex2 * MeV);
        heavy->SetDefinition(heavyDef);
        G4double mB = heavyDef->GetPDGMass();

        // define light decay product of heacy recoil, which is alpha
        G4DynamicParticle *d1 = new G4DynamicParticle;
        G4ParticleDefinition *d1Def = G4IonTable::GetIonTable()->GetIon(2, 4, 0.0);
        d1->SetDefinition(d1Def);
        G4double mD1 = d1Def->GetPDGMass();

        // define heavy decay product of heacy recoil, which is 12C(g.s.)
        G4DynamicParticle *d2 = new G4DynamicParticle;
        G4ParticleDefinition *d2Def = G4IonTable::GetIonTable()->GetIon(6, 12, 0.0);
        d2->SetDefinition(d2Def);
        G4double mD2 = d2Def->GetPDGMass();

        // get beam properties
        G4LorentzVector tempPa = aTrack.GetDynamicParticle()->Get4Momentum();
        TLorentzVector Pa(tempPa.px(), tempPa.py(), tempPa.pz(), tempPa.e());

        //
        // a+A->b+B*
        //
        // boost to CM of projectile+target (CM1)
        TLorentzVector Pt(0, 0, 0, mA);
        TLorentzVector Pcm1 = Pa + Pt;
        G4double Mcm1 = Pcm1.M();
        // calculate in CM1
        G4double mag_Pcm1 = std::sqrt((std::pow(Mcm1, 2) - std::pow(mB + mb, 2)) * (std::pow(Mcm1, 2) - std::pow(mB - mb, 2))) / 2. / Mcm1;
        G4double E_bcm1 = (Mcm1 * Mcm1 + mb * mb - mB * mB) / 2. / Mcm1;
        G4double E_Bcm1 = (Mcm1 * Mcm1 + mB * mB - mb * mb) / 2. / Mcm1;
        TVector3 vector3Pb(mag_Pcm1, 0, 0);
        vector3Pb.SetTheta(std::acos(-1. + 2. * G4UniformRand()));
        vector3Pb.SetPhi(G4UniformRand() * 2 * M_PI);
        TVector3 vector3PB = (-1.) * vector3Pb;
        TLorentzVector Pb(vector3Pb, E_bcm1);
        TLorentzVector PB(vector3PB, E_Bcm1);
        // boost b and B back to lab
        Pb.Boost(Pcm1.BoostVector());
        PB.Boost(Pcm1.BoostVector());

        //
        // B*->D1+D2
        //
        // calculate in CM of B* (CM2)
        G4double Mcm2 = PB.M();
        G4double mag_Pcm2 = std::sqrt((std::pow(Mcm2, 2) - std::pow(mD1 + mD2, 2)) * (std::pow(Mcm2, 2) - std::pow(mD1 - mD2, 2))) / 2. / Mcm2;
        G4double E_D1cm2 = (Mcm2 * Mcm2 + mD1 * mD1 - mD2 * mD2) / 2. / Mcm2;
        G4double E_D2cm2 = (Mcm2 * Mcm2 + mD2 * mD2 - mD1 * mD1) / 2. / Mcm2;
        TVector3 vector3PD1(mag_Pcm2, 0, 0);
        G4double ThetaCM2 = std::acos(-1. + 2. * G4UniformRand());
        G4double PhiCM2 = G4UniformRand() * 2 * M_PI;
        vector3PD1.SetTheta(ThetaCM2);
        vector3PD1.SetPhi(PhiCM2);
        TVector3 vector3PD2 = (-1.) * vector3PD1;
        TLorentzVector PD1(vector3PD1, E_D1cm2);
        TLorentzVector PD2(vector3PD2, E_D2cm2);
        // boost b and B back to lab
        PD1.Boost(PB.BoostVector());
        PD2.Boost(PB.BoostVector());

        // set secondary particles
        light->Set4Momentum(G4LorentzVector(G4ThreeVector(Pb.Px(), Pb.Py(), Pb.Pz()), Pb.E()));
        heavy->Set4Momentum(G4LorentzVector(G4ThreeVector(PB.Px(), PB.Py(), PB.Pz()), PB.E()));
        d1->Set4Momentum(G4LorentzVector(G4ThreeVector(PD1.Px(), PD1.Py(), PD1.Pz()), PD1.E()));
        d2->Set4Momentum(G4LorentzVector(G4ThreeVector(PD2.Px(), PD2.Py(), PD2.Pz()), PD2.E()));

        // set secondary tracks
        G4Track *sec_light = new G4Track(light,
                                         aTrack.GetGlobalTime(),
                                         aTrack.GetPosition());
        G4Track *sec_heavy = new G4Track(heavy,
                                         aTrack.GetGlobalTime(),
                                         aTrack.GetPosition());
        G4Track *sec_d1 = new G4Track(d1,
                                      aTrack.GetGlobalTime(),
                                      aTrack.GetPosition());
        G4Track *sec_d2 = new G4Track(d2,
                                      aTrack.GetGlobalTime(),
                                      aTrack.GetPosition());

        // add secondary tracks
        aParticleChange.SetNumberOfSecondaries(3);
        aParticleChange.AddSecondary(sec_light);
        aParticleChange.AddSecondary(sec_d1);
        aParticleChange.AddSecondary(sec_d2);

        aParticleChange.ProposeEnergy(0.);
        aParticleChange.ProposeTrackStatus(fStopAndKill);
        return &aParticleChange;
    }

    G4VParticleChange *ReactionProcess::TwoBodyNeutronDecay(const G4Track &aTrack, const G4Step &aStep,
                                                            int Zt, int At,             // target
                                                            int Z1, int A1, double Ex1, // light recoil
                                                            int Z2, int A2, double Ex2  // heavy recoil
    )
    {

        //
        // 2 steps reaction (1  + 1 decays)
        // Step1: a+A->b+B*, Step2: B*->D1+D2 (D1 is neutron)
        //

        // define target
        G4DynamicParticle *target = new G4DynamicParticle;
        G4ParticleDefinition *targetDef;
        targetDef = G4IonTable::GetIonTable()->GetIon(Zt, At, 0.);
        target->SetDefinition(targetDef);
        G4double mA = targetDef->GetPDGMass();

        // define light recoil
        G4DynamicParticle *light = new G4DynamicParticle;
        G4ParticleDefinition *lightDef;
        if (Z1 == 0)
        {
            G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
            lightDef = particleTable->FindParticle("neutron");
        }
        if (Z1 != 0)
        {
            lightDef = G4IonTable::GetIonTable()->GetIon(Z1, A1, Ex1 * MeV);
        }
        light->SetDefinition(lightDef);
        G4double mb = lightDef->GetPDGMass();

        // define heavy recoil
        G4DynamicParticle *heavy = new G4DynamicParticle;
        G4ParticleDefinition *heavyDef;
        heavyDef = G4IonTable::GetIonTable()->GetIon(Z2, A2, Ex2 * MeV);
        heavy->SetDefinition(heavyDef);
        G4double mB = heavyDef->GetPDGMass();

        // define light decay product of heacy recoil, which is neutron
        G4DynamicParticle *d1 = new G4DynamicParticle;
        G4ParticleDefinition *d1Def = G4ParticleTable::GetParticleTable()->FindParticle("neutron");
        d1->SetDefinition(d1Def);
        G4double mD1 = d1Def->GetPDGMass();

        // define heavy decay product of heacy recoil, which is 12C(g.s.)
        G4DynamicParticle *d2 = new G4DynamicParticle;
        G4ParticleDefinition *d2Def = G4IonTable::GetIonTable()->GetIon(6, 12, 0.0);
        d2->SetDefinition(d2Def);
        G4double mD2 = d2Def->GetPDGMass();

        // get beam properties
        G4LorentzVector tempPa = aTrack.GetDynamicParticle()->Get4Momentum();
        TLorentzVector Pa(tempPa.px(), tempPa.py(), tempPa.pz(), tempPa.e());

        //
        // a+A->b+B*
        //
        // boost to CM of projectile+target (CM1)
        TLorentzVector Pt(0, 0, 0, mA);
        TLorentzVector Pcm1 = Pa + Pt;
        G4double Mcm1 = Pcm1.M();
        // calculate in CM1
        G4double mag_Pcm1 = std::sqrt((std::pow(Mcm1, 2) - std::pow(mB + mb, 2)) * (std::pow(Mcm1, 2) - std::pow(mB - mb, 2))) / 2. / Mcm1;
        G4double E_bcm1 = (Mcm1 * Mcm1 + mb * mb - mB * mB) / 2. / Mcm1;
        G4double E_Bcm1 = (Mcm1 * Mcm1 + mB * mB - mb * mb) / 2. / Mcm1;
        TVector3 vector3Pb(mag_Pcm1, 0, 0);
        vector3Pb.SetTheta(std::acos(-1. + 2. * G4UniformRand()));
        vector3Pb.SetPhi(G4UniformRand() * 2 * M_PI);
        TVector3 vector3PB = (-1.) * vector3Pb;
        TLorentzVector Pb(vector3Pb, E_bcm1);
        TLorentzVector PB(vector3PB, E_Bcm1);
        // boost b and B back to lab
        Pb.Boost(Pcm1.BoostVector());
        PB.Boost(Pcm1.BoostVector());

        //
        // B*->D1+D2
        //
        // calculate in CM of B* (CM2)
        G4double Mcm2 = PB.M();
        G4double mag_Pcm2 = std::sqrt((std::pow(Mcm2, 2) - std::pow(mD1 + mD2, 2)) * (std::pow(Mcm2, 2) - std::pow(mD1 - mD2, 2))) / 2. / Mcm2;
        G4double E_D1cm2 = (Mcm2 * Mcm2 + mD1 * mD1 - mD2 * mD2) / 2. / Mcm2;
        G4double E_D2cm2 = (Mcm2 * Mcm2 + mD2 * mD2 - mD1 * mD1) / 2. / Mcm2;
        TVector3 vector3PD1(mag_Pcm2, 0, 0);
        G4double ThetaCM2 = std::acos(-1. + 2. * G4UniformRand());
        G4double PhiCM2 = G4UniformRand() * 2 * M_PI;
        vector3PD1.SetTheta(ThetaCM2);
        vector3PD1.SetPhi(PhiCM2);
        TVector3 vector3PD2 = (-1.) * vector3PD1;
        TLorentzVector PD1(vector3PD1, E_D1cm2);
        TLorentzVector PD2(vector3PD2, E_D2cm2);
        // boost b and B back to lab
        PD1.Boost(PB.BoostVector());
        PD2.Boost(PB.BoostVector());

        // set secondary particles
        light->Set4Momentum(G4LorentzVector(G4ThreeVector(Pb.Px(), Pb.Py(), Pb.Pz()), Pb.E()));
        heavy->Set4Momentum(G4LorentzVector(G4ThreeVector(PB.Px(), PB.Py(), PB.Pz()), PB.E()));
        d1->Set4Momentum(G4LorentzVector(G4ThreeVector(PD1.Px(), PD1.Py(), PD1.Pz()), PD1.E()));
        d2->Set4Momentum(G4LorentzVector(G4ThreeVector(PD2.Px(), PD2.Py(), PD2.Pz()), PD2.E()));

        // set secondary tracks
        G4Track *sec_light = new G4Track(light,
                                         aTrack.GetGlobalTime(),
                                         aTrack.GetPosition());
        G4Track *sec_heavy = new G4Track(heavy,
                                         aTrack.GetGlobalTime(),
                                         aTrack.GetPosition());
        G4Track *sec_d1 = new G4Track(d1,
                                      aTrack.GetGlobalTime(),
                                      aTrack.GetPosition());
        G4Track *sec_d2 = new G4Track(d2,
                                      aTrack.GetGlobalTime(),
                                      aTrack.GetPosition());

        // add secondary tracks
        aParticleChange.SetNumberOfSecondaries(3);
        aParticleChange.AddSecondary(sec_light);
        aParticleChange.AddSecondary(sec_d1);
        aParticleChange.AddSecondary(sec_d2);

        aParticleChange.ProposeEnergy(0.);
        aParticleChange.ProposeTrackStatus(fStopAndKill);
        return &aParticleChange;
    }
}