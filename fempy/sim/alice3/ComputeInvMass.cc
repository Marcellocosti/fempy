#include <stdio.h>
#include <stdlib.h>

#include "Riostream.h"
#include "TDatabasePDG.h"
#include "TFile.h"
#include "TH2F.h"
#include "TList.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TTree.h"

#include "ActsFatras/EventData/Barcode.hpp"

struct particle {
    int pdg;                        // pdg code
    uint64_t id;                    // index of the particle in the event
    double t_x;                     // production vertex x
    double t_y;                     // production vertex y
    double t_z;                     // production vertex z
    ROOT::Math::PxPyPzMVector p;    // reconstructed quadrimomentum
    ROOT::Math::PxPyPzMVector t_p;  // true quadrimomentum
    int mother1idx;                 // index of mother 1
    int mother2idx;                 // index of mother 2
    int mother1pdg;                 // pdg code of mother 1
};

bool AreSiblings(particle p1, particle p2) { return std::abs(static_cast<int>(p1.id - p2.id)) == 1; }

particle BuildMother2Daus(int pdg, particle p1, particle p2) {
    // Check if particles have the same mother. Only check for weak decays. In this case m1>0 and and m2=0
    // see: https://pythia.org/latest-manual/ParticleProperties.html
    int abspdg = std::abs(pdg);
    double targetMass = TDatabasePDG::Instance()->GetParticle(abspdg)->Mass();
    particle mother({});

    if (abspdg == 421) {
        if (std::abs(p1.mother1pdg) != pdg || std::abs(p2.mother1pdg) != pdg) return particle({});

        if (p1.pdg * p2.pdg >= 0) return particle({});

        if (p1.mother1idx != p2.mother1idx || p1.mother1idx <= 0 || p2.mother1idx <= 0 || p1.mother2idx != 0 ||
            p2.mother2idx != 0)
            return particle({});

        mother = particle({p1.mother1pdg, (u_int64_t)p1.mother1idx, 0, 0, 0, p1.p + p2.p, p1.t_p + p2.t_p, 0, 0, 0});

        if (std::abs(mother.t_p.M() - targetMass) / targetMass > 1.e-4) return particle({});
        if (!AreSiblings(p1, p2)) return particle({});
    } else if (abspdg == 413) {
        particle Dzero, piSoft;
        if (std::abs(p1.pdg) == 421 && std::abs(p2.pdg) == 211) {
            Dzero = p1;
            piSoft = p2;
        } else if (std::abs(p1.pdg) == 211 && std::abs(p2.pdg) == 421) {
            piSoft = p1;
            Dzero = p2;
        } else
            return particle({});

        if (std::abs(piSoft.mother1pdg) != 413 || piSoft.mother1idx <= 0 || piSoft.mother2idx != 0) return particle({});

        mother =
            particle({413, (u_int64_t)piSoft.mother1idx, 0, 0, 0, Dzero.p + piSoft.p, Dzero.t_p + piSoft.t_p, 0, 0, 0});
        if (std::abs(mother.t_p.M() - targetMass) / targetMass > 1.e-4) return particle({});
    }

    return mother;
}

particle BuildMother3Prong(int pdg, particle p1, particle p2, particle p3) {
    if (pdg == 413) {
        particle K, pi1, pi2;
        if (std::abs(p1.pdg) == 321 && std::abs(p2.pdg) == 211 && std::abs(p3.pdg)) {
            K = p1;
            pi1 = p2;
            pi2 = p3;
        } else if (std::abs(p2.pdg) == 321 && std::abs(p3.pdg) == 211 && std::abs(p1.pdg)) {
            K = p2;
            pi1 = p1;
            pi2 = p3;
        } else if (std::abs(p3.pdg) == 321 && std::abs(p1.pdg) == 211 && std::abs(p2.pdg)) {
            K = p3;
            pi1 = p1;
            pi2 = p2;
        } else
            return particle({});

        particle Dzero, piSoft;
        if (K.pdg * pi1.pdg < 0) {
            Dzero = BuildMother2Daus(421, K, pi1);
            piSoft = pi2;
        } else if (K.pdg * pi2.pdg < 0) {
            Dzero = BuildMother2Daus(421, K, pi2);
            piSoft = pi1;
        } else
            return particle({});

        if (Dzero.pdg != 421) return particle({});

        // Check if particles have the same mother. Only check for weak decays. In this case m1>0 and and m2=0
        // see: https://pythia.org/latest-manual/ParticleProperties.html
        if (piSoft.mother1idx <= 0 || piSoft.mother2idx != 0 || piSoft.mother1pdg != 413) return particle({});

        // id from acts and the motherId from pythia don't match => don't check if the D0 and piSoft are siblings
        particle mother(
            {pdg, (u_int64_t)piSoft.mother1idx, 0, 0, 0, Dzero.p + piSoft.p, Dzero.t_p + piSoft.t_p, 0, 0, 0});

        double targetMass = TDatabasePDG::Instance()->GetParticle(std::abs(413))->Mass();
        if (std::abs(mother.t_p.M() - targetMass) / targetMass > 1.e-4) return particle({});

        return mother;
    } else {
        printf("Error: pdg code %d not implemented. Exit\n", pdg);
        exit(1);
    }
    return particle({});
}

void ComputeInvMass(const char *inFileName, const char *oFileName, int pdg1, int pdg2) {
    // consts
    const double Pi = TMath::Pi();
    int minNHits = 7;

    TFile *inFile = new TFile(inFileName);

    TList *keys = (TList *)inFile->GetListOfKeys();
    if (keys->GetEntries() > 1) {
        printf("Warning: more than 1 tree found. Check\n");
    }

    // init to 0 otherwise it breaks
    std::vector<int> *particle_type = 0;
    std::vector<int> *mother1_particle_id = 0;
    std::vector<int> *mother2_particle_id = 0;
    std::vector<int> *mother1_pdg = 0;
    std::vector<int> *mother2_pdg = 0;
    std::vector<int> *nMeasurements = 0;
    std::vector<double> *px = 0;
    std::vector<double> *py = 0;
    std::vector<double> *pz = 0;

    std::vector<double> *qop = 0;  // charge over momentum
    std::vector<double> *phi = 0;
    std::vector<double> *theta = 0;

    std::vector<double> *vx = 0;
    std::vector<double> *vy = 0;
    std::vector<double> *vz = 0;

    // MC truth variables
    std::vector<double> *t_vx = 0;
    std::vector<double> *t_vy = 0;
    std::vector<double> *t_vz = 0;
    std::vector<double> *t_px = 0;
    std::vector<double> *t_py = 0;
    std::vector<double> *t_pz = 0;
    std::vector<int> *t_charge = 0;
    std::vector<uint64_t> *t_majorityParticleId = 0;

    auto treeName = std::string(keys->At(0)->GetName());
    TTree *tree = (TTree *)inFile->Get(treeName.data());
    tree->SetBranchAddress("eQOP_fit", &qop);
    tree->SetBranchAddress("ePHI_fit", &phi);
    tree->SetBranchAddress("eTHETA_fit", &theta);

    // MC truth vaariables
    tree->SetBranchAddress("particle_type", &particle_type);
    tree->SetBranchAddress("mother1_particle_id", &mother1_particle_id);
    tree->SetBranchAddress("mother2_particle_id", &mother2_particle_id);
    tree->SetBranchAddress("mother1_pdg", &mother1_pdg);
    tree->SetBranchAddress("mother2_pdg", &mother2_pdg);
    tree->SetBranchAddress("nMeasurements", &nMeasurements);
    tree->SetBranchAddress("t_vx", &t_vx);
    tree->SetBranchAddress("t_vy", &t_vy);
    tree->SetBranchAddress("t_vz", &t_vz);
    tree->SetBranchAddress("t_px", &t_px);
    tree->SetBranchAddress("t_py", &t_py);
    tree->SetBranchAddress("t_pz", &t_pz);
    tree->SetBranchAddress("t_charge", &t_charge);
    tree->SetBranchAddress("majorityParticleId", &t_majorityParticleId);

    TFile *oFile = new TFile(oFileName, "recreate");
    double massMin = 0.4;
    double massMax = 2.5;
    TString titleMassDzero = "#it{M}(K,#pi) (GeV/#it{c}^{2})";
    TString titleMassDstar = "#it{M}(K,#pi,#pi) (GeV/#it{c}^{2})";
    TString titlePt = "#it{p}_{T} (GeV/#it{c})";

    std::vector<int> lightPdg = {211, 321};
    std::vector<int> charmPdg = {421, 413};
    std::vector<int> allPdg = lightPdg;
    allPdg.insert(allPdg.end(), charmPdg.begin(), charmPdg.end());

    // names
    std::map<int, const char *> pdg2name = {
        {321, "K"},
        {211, "Pi"},
        {421, "Dzero"},
        {413, "Dstar"},
    };

    std::map<int, std::vector<particle>> particles = {
        {321, std::vector<particle>()},
        {211, std::vector<particle>()},
        {421, std::vector<particle>()},
        {413, std::vector<particle>()},
    };

    // this map also includes TH2*
    std::map<int, std::map<std::string, TH1 *>> hPartProp = {
        {321, {}},
        {211, {}},
        {421, {}},
        {413, {}},
    };

    TString name;
    TString title;
    double massBinWidth = 0.1;  // MeV/c2

    for (auto abspdg : {211, 321, 421, 413}) {
        // p
        name = Form("hP_%s", pdg2name[abspdg]);
        title = ";#it{p} (GeV/#it{c});Counts";
        hPartProp[abspdg].insert({"hP", new TH1D(name, title, 500, 0, 50)});

        // pt
        name = Form("hPt_%s", pdg2name[abspdg]);
        title = ";#it{p}_{T} (GeV/#it{c});Counts";
        hPartProp[abspdg].insert({"hPt", new TH1D(name, title, 500, 0, 50)});

        // eta
        name = Form("hEta_%s", pdg2name[abspdg]);
        title = ";#eta;Counts";
        hPartProp[abspdg].insert({"hEta", new TH1D(name, title, 500, -5, 5)});

        // phi
        name = Form("hPhi_%s", pdg2name[abspdg]);
        title = ";#varphi (rad);Counts";
        hPartProp[abspdg].insert({"hPhi", new TH1D(name, title, 500, -Pi, Pi)});

        // phi vs eta
        name = Form("hPhiVsEta_%s", pdg2name[abspdg]);
        title = ";#eta;#varphi (rad);Counts";
        hPartProp[abspdg].insert({"hPhiVsEta", new TH2F(name, title, 200, -5, 5, 200, -Pi, Pi)});

        // phi vs rProd
        name = Form("hPhiVsRprod_%s", pdg2name[abspdg]);
        title = ";#it{r}_{prod} (mm);#varphi (rad);Counts";
        hPartProp[abspdg].insert({"hPhiVsRprod", new TH2F(name, title, 200, 0, 5, 200, -Pi, Pi)});

        // P resolution
        name = Form("hResolutionP_%s", pdg2name[abspdg]);
        title = ";#it{p}^{true} (GeV/#it{c});#it{p}^{reco} (GeV/#it{c});Counts";
        hPartProp[abspdg].insert({"hResolutionP", new TH2F(name, title, 100, 0, 10, 100, 0, 10)});

        // P resolution delta
        name = Form("hResolutionDeltaP_%s", pdg2name[abspdg]);
        title = ";#it{p}^{true} (GeV/#it{c});#it{p}^{reco} - #it{p}^{true} (GeV/#it{c});Counts";
        hPartProp[abspdg].insert({"hResolutionDeltaP", new TH2F(name, title, 100, 0, 10, 100, -0.5, 0.5)});

        // P resolution percentage
        name = Form("hResolutionPercP_%s", pdg2name[abspdg]);
        title = ";#it{p}^{true} (GeV/#it{c});(#it{p}^{reco} - #it{p}^{true})/#it{p}^{true}) (%);Counts";
        hPartProp[abspdg].insert({"hResolutionPercP", new TH2F(name, title, 100, 0, 10, 100, -5, 5)});

        // Pt resolution
        name = Form("hResolutionPt_%s", pdg2name[abspdg]);
        title = ";#it{p}_{T}^{true} (GeV/#it{c});#it{p}_{T}^{reco} (GeV/#it{c});Counts";
        hPartProp[abspdg].insert({"hResolutionPt", new TH2F(name, title, 100, 0, 10, 100, 0, 10)});

        // Pt resolution delta
        name = Form("hResolutionDeltaPt_%s", pdg2name[abspdg]);
        title = ";#it{p}_{T}^{true} (GeV/#it{c});#it{p}_{T}^{reco} - #it{p}_{T}^{true} (GeV/#it{c});Counts";
        hPartProp[abspdg].insert({"hResolutionDeltaPt", new TH2F(name, title, 100, 0, 10, 100, -0.5, 0.5)});

        // Pt resolution percentage
        name = Form("hResolutionPercPt_%s", pdg2name[abspdg]);
        title = ";#it{p}_{T}^{true} (GeV/#it{c});(#it{p}_{T}^{reco} - #it{p}_{T}^{true})/#it{p}_{T}^{true} (%);Counts";
        hPartProp[abspdg].insert({"hResolutionPercPt", new TH2F(name, title, 100, 0, 10, 100, -5, 5)});

        // Eta resolution
        name = Form("hResolutionEta_%s", pdg2name[abspdg]);
        title = ";#eta^{true};#eta^{reco};Counts";
        hPartProp[abspdg].insert({"hResolutionEta", new TH2F(name, title, 200, -5, 5, 200, -5, 5)});

        // Eta resolution delta
        name = Form("hResolutionDeltaEta_%s", pdg2name[abspdg]);
        title = ";#eta^{true};#eta^{reco} - #eta^{true};Counts";
        hPartProp[abspdg].insert({"hResolutionDeltaEta", new TH2F(name, title, 200, -5, 5, 200, -0.1, 0.1)});

        // Eta resolution percentage
        name = Form("hResolutionPercEta_%s", pdg2name[abspdg]);
        title = ";#eta^{true};(#eta^{reco} - #eta^{true})/#eta^{true} (%);Counts";
        hPartProp[abspdg].insert({"hResolutionPercEta", new TH2F(name, title, 200, -5, 5, 200, -5, 5)});

        // Phi resolution
        name = Form("hResolutionPhi_%s", pdg2name[abspdg]);
        title = ";#phi^{true};#phi^{reco};Counts";
        hPartProp[abspdg].insert({"hResolutionPhi", new TH2F(name, title, 200, -Pi, Pi, 200, -Pi, Pi)});

        // Phi resolution delta
        name = Form("hResolutionDeltaPhi_%s", pdg2name[abspdg]);
        title = ";#phi^{true};#phi^{reco} - #phi^{true};Counts";
        hPartProp[abspdg].insert({"hResolutionDeltaPhi", new TH2F(name, title, 200, -Pi, Pi, 200, -Pi, Pi)});

        // Phi resolution percentage
        name = Form("hResolutionPercPhi_%s", pdg2name[abspdg]);
        title = ";#phi^{true};(#phi^{reco} - #phi^{true})/#phi^{true} (%);Counts";
        hPartProp[abspdg].insert({"hResolutionPercPhi", new TH2F(name, title, 200, -Pi, Pi, 200, -Pi, Pi)});

        if (std::find(charmPdg.begin(), charmPdg.end(), abspdg) != charmPdg.end()) {
            TString titleMass;
            int nMassBins;
            if (abspdg == 421) {
                titleMass = titleMassDzero;
                massMin = 1.5;
                massMax = 2.2;
                nMassBins = (int)std::round((massMax - massMin) * 1000 / massBinWidth);
            } else if (abspdg == 413) {
                titleMass = titleMassDstar;
                massMin = 1.6;
                massMax = 2.4;
                nMassBins = (int)std::round((massMax - massMin) * 1000 / massBinWidth);
            } else {
                exit(1);
            }

            // invariant mass
            name = Form("hInvMass_%s", pdg2name[abspdg]);
            title = ";" + titleMass + ";Counts";
            hPartProp[abspdg].insert({"hInvMass", new TH1D(name, title, nMassBins, massMin, massMax)});

            // invariant mass vs pt
            name = Form("hInvMassVsPt_%s", pdg2name[abspdg]);
            title = ";" + titlePt + ";" + titleMass + ";Counts";
            hPartProp[abspdg].insert({"hInvMassVsPt", new TH2F(name, title, 100, 0, 10, nMassBins, massMin, massMax)});
        }
    }

    for (int iEvent = 0; iEvent < tree->GetEntries(); iEvent++) {
        tree->GetEntry(iEvent);

        particles[321].clear();
        particles[211].clear();
        particles[421].clear();
        particles[413].clear();

        size_t nPart = particle_type->size();
        for (size_t iPart = 0; iPart < nPart; iPart++) {
            int pdg = (*particle_type)[iPart];
            int abspdg = std::abs(pdg);
            if (abspdg != 211 && abspdg != 321) continue;

            int nHits = (*nMeasurements)[iPart];
            if (nHits < minNHits) continue;

            double t_ppx = (*t_px)[iPart];
            double t_ppy = (*t_py)[iPart];
            double t_ppz = (*t_pz)[iPart];
            uint64_t t_pmajorityParticleId = (*t_majorityParticleId)[iPart];

            double pqop = (*qop)[iPart];
            double pphi = (*phi)[iPart];
            double ptheta = (*theta)[iPart];

            double ppx = 1. / std::abs(pqop) * sin(ptheta) * cos(pphi);
            double ppy = 1. / std::abs(pqop) * sin(ptheta) * sin(pphi);
            double ppz = 1. / std::abs(pqop) * cos(ptheta);
            double pm = TDatabasePDG::Instance()->GetParticle(abspdg)->Mass();

            uint64_t id = ActsFatras::Barcode(t_pmajorityParticleId).particle();

            auto p = ROOT::Math::PxPyPzMVector(ppx, ppy, ppz, pm);
            auto t_p = ROOT::Math::PxPyPzMVector(t_ppx, t_ppy, t_ppz, pm);

            double x = (*t_vx)[iPart]; // u. m. = mm (should be)
            double y = (*t_vy)[iPart]; // u. m. = mm (should be)
            double z = (*t_vz)[iPart]; // u. m. = mm (should be)

            // reject decay products of strange decays. This selection also removes the peak in eta-phi for pions
            double rProd = pow(x*x + y*y + z*z, 0.5); 
           
            int m1idx = (*mother1_particle_id)[iPart];
            int m2idx = (*mother2_particle_id)[iPart];

            int m1pdg = (*mother1_pdg)[iPart];
            particle part = particle({pdg, id, x, y, z, p, t_p, m1idx, m2idx, m1pdg});
            particles[abspdg].push_back(part);

            // Fill hists
            hPartProp[abspdg]["hP"]->Fill(p.P());
            hPartProp[abspdg]["hPt"]->Fill(p.Pt());
            hPartProp[abspdg]["hEta"]->Fill(p.Eta());
            hPartProp[abspdg]["hPhi"]->Fill(p.Phi());
            hPartProp[abspdg]["hPhiVsEta"]->Fill(p.Eta(), p.Phi());
            hPartProp[abspdg]["hPhiVsRprod"]->Fill(rProd, p.Phi());

            hPartProp[abspdg]["hResolutionEta"]->Fill(t_p.Eta(), p.Eta());
            hPartProp[abspdg]["hResolutionDeltaEta"]->Fill(t_p.Eta(), p.Eta() - t_p.Eta());
            hPartProp[abspdg]["hResolutionPercEta"]->Fill(t_p.Eta(), (p.Eta() - t_p.Eta()) / t_p.Eta());

            hPartProp[abspdg]["hResolutionPhi"]->Fill(t_p.Phi(), p.Phi());
            hPartProp[abspdg]["hResolutionDeltaPhi"]->Fill(t_p.Phi(), p.Phi() - t_p.Phi());
            hPartProp[abspdg]["hResolutionPercPhi"]->Fill(t_p.Phi(), (p.Phi() - t_p.Phi()) / t_p.Phi());

            hPartProp[abspdg]["hResolutionP"]->Fill(t_p.P(), p.P());
            hPartProp[abspdg]["hResolutionDeltaP"]->Fill(t_p.P(), p.P() - t_p.P());
            hPartProp[abspdg]["hResolutionPercP"]->Fill(t_p.P(), (p.P() - t_p.P()) / t_p.P() * 100);

            hPartProp[abspdg]["hResolutionPt"]->Fill(t_p.Pt(), p.Pt());
            hPartProp[abspdg]["hResolutionDeltaPt"]->Fill(t_p.Pt(), p.Pt() - t_p.Pt());
            hPartProp[abspdg]["hResolutionPercPt"]->Fill(t_p.Pt(), (p.Pt() - t_p.Pt()) / t_p.Pt() * 100);
        }

        // reconstruct Dzero
        auto kaons = particles[321];
        auto pions = particles[211];

        for (size_t iK = 0; iK < kaons.size(); iK++) {
            auto K = kaons[iK];

            for (size_t iPi = 0; iPi < pions.size(); iPi++) {
                auto Pi = pions[iPi];
                auto Dzero = BuildMother2Daus(421, K, Pi);
                if (Dzero.pdg != 421) continue;
                particles[421].push_back(Dzero);

                double m = Dzero.p.M();
                double eta = Dzero.p.Eta();
                double phi = Dzero.p.Phi();
                double p = Dzero.p.P();
                double pt = Dzero.p.Pt();

                double t_eta = Dzero.t_p.Eta();
                double t_phi = Dzero.t_p.Phi();
                double t_p = Dzero.t_p.P();
                double t_pt = Dzero.t_p.Pt();

                // kinematics
                hPartProp[421]["hP"]->Fill(p);
                hPartProp[421]["hPt"]->Fill(pt);
                hPartProp[421]["hEta"]->Fill(eta);
                hPartProp[421]["hPhi"]->Fill(phi);
                hPartProp[421]["hPhiVsEta"]->Fill(eta, phi);

                // Invariant mass
                hPartProp[421]["hInvMass"]->Fill(m);
                hPartProp[421]["hInvMassVsPt"]->Fill(pt, m);

                // Resolution
                hPartProp[421]["hResolutionEta"]->Fill(t_eta, eta);
                hPartProp[421]["hResolutionDeltaEta"]->Fill(t_eta, eta - t_eta);
                hPartProp[421]["hResolutionPercEta"]->Fill(t_eta, (eta - t_eta) / t_eta * 100);

                hPartProp[421]["hResolutionPhi"]->Fill(t_phi, phi);
                hPartProp[421]["hResolutionDeltaPhi"]->Fill(t_phi, phi - t_phi);
                hPartProp[421]["hResolutionPercPhi"]->Fill(t_phi, (phi - t_phi) / t_phi * 100);

                hPartProp[421]["hResolutionP"]->Fill(t_p, p);
                hPartProp[421]["hResolutionDeltaP"]->Fill(t_p, p - t_p);
                hPartProp[421]["hResolutionPercP"]->Fill(t_p, (p - t_p) / t_p * 100);

                hPartProp[421]["hResolutionPt"]->Fill(t_pt, pt);
                hPartProp[421]["hResolutionDeltaPt"]->Fill(t_pt, pt - t_pt);
                hPartProp[421]["hResolutionPercPt"]->Fill(t_pt, (pt - t_pt) / t_pt * 100);
            }
        }

        // for (size_t iK = 0; iK < kaons.size(); iK++) {
        //     auto K = kaons[iK];

        //     for (size_t iPi1 = 0; iPi1 < pions.size(); iPi1++) {
        //         auto Pi1 = pions[iPi1];

        //         for (size_t iPi2 = iPi1 + 1; iPi2 < pions.size(); iPi2++) {
        //             auto Pi2 = pions[iPi2];

        //             auto Dstar = BuildMother3Prong(413, K, Pi1, Pi2);
        //             if (Dstar.pdg != 413) continue;
        //             particles[413].push_back(Dstar);

        //             double m = Dstar.p.M();
        //             double eta = Dstar.p.Eta();
        //             double phi = Dstar.p.Phi();
        //             double p = Dstar.p.P();
        //             double pt = Dstar.p.Pt();

        //             double t_eta = Dstar.t_p.Eta();
        //             double t_phi = Dstar.t_p.Phi();
        //             double t_p = Dstar.t_p.P();
        //             double t_pt = Dstar.t_p.Pt();

        //             // kinematics
        //             hPartProp[413]["hP"]->Fill(p);
        //             hPartProp[413]["hPt"]->Fill(pt);
        //             hPartProp[413]["hEta"]->Fill(eta);
        //             hPartProp[413]["hPhi"]->Fill(phi);
        //             hPartProp[413]["hPhiVsEta"]->Fill(eta, phi);

        //             // Invariant mass
        //             hPartProp[413]["hInvMass"]->Fill(m);
        //             hPartProp[413]["hInvMassVsPt"]->Fill(pt, m);

        //             // Resolution
        //             hPartProp[413]["hResolutionEta"]->Fill(t_eta, eta);
        //             hPartProp[413]["hResolutionDeltaEta"]->Fill(t_eta, eta - t_eta);
        //             hPartProp[413]["hResolutionPercEta"]->Fill(t_eta, (eta - t_eta) / t_eta * 100);

        //             hPartProp[413]["hResolutionPhi"]->Fill(t_phi, phi);
        //             hPartProp[413]["hResolutionDeltaPhi"]->Fill(t_phi, phi - t_phi);
        //             hPartProp[413]["hResolutionPercPhi"]->Fill(t_phi, (phi - t_phi) / t_phi * 100);

        //             hPartProp[413]["hResolutionP"]->Fill(t_p, p);
        //             hPartProp[413]["hResolutionDeltaP"]->Fill(t_p, p - t_p);
        //             hPartProp[413]["hResolutionPercP"]->Fill(t_p, (p - t_p) / t_p * 100);

        //             hPartProp[413]["hResolutionPt"]->Fill(t_pt, pt);
        //             hPartProp[413]["hResolutionDeltaPt"]->Fill(t_pt, pt - t_pt);
        //             hPartProp[413]["hResolutionPercPt"]->Fill(t_pt, (pt - t_pt) / t_pt * 100);
        //         }
        //     }
        // }
    }  // event loop

    // Write histograms
    for (auto pdg : allPdg) {
        oFile->mkdir(pdg2name[pdg]);
        oFile->cd(pdg2name[pdg]);

        // change the name of the qa histograms and save them to file
        for (auto const &[name, hist] : hPartProp[pdg]) {
            hist->SetName(name.data());
            hist->Write();
        }
    }
    oFile->Close();
}
