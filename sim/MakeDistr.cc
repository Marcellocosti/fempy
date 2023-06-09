#include <TClonesArray.h>
#include <TDatabasePDG.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TParticle.h>
#include <TROOT.h>
#include <TRandom3.h>
#include <TStopwatch.h>
#include <TSystem.h>
#include <TTree.h>

#include <array>
#include <deque>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "Math/GenVector/Boost.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "functions.hxx"


void MakeDistr(std::string inFileName = "/home/ktas/ge86rim/phsw/fempy/sim/AnalysisResults_merged.root",
               std::string oFileName = "test.root", const int cpdg = 411, const int lpdg = 211,
               bool selDauKinem = false, bool align = false, bool cleanPairs=true) {
    int md = 15;

    std::map<std::string, TH1D *> hSE = {
        {"sc", new TH1D("hSE_sc", ";#it{k}* (GeV/#it{c});pairs", 2000, 0., 2.)},
        {"oc", new TH1D("hSE_oc", ";#it{k}* (GeV/#it{c});pairs", 2000, 0., 2.)},
    };

    std::map<std::string, TH1D *> hME = {
        {"sc", new TH1D("hME_sc", ";#it{k}* (GeV/#it{c});pairs", 2000, 0., 2.)},
        {"oc", new TH1D("hME_oc", ";#it{k}* (GeV/#it{c});pairs", 2000, 0., 2.)},
    };

    // Load dataset file
    TFile *inFile = TFile::Open(inFileName.data());
    if (!inFile) {
        printf("The file %s does not exist. Skip!", inFileName.data());
        return;
    }

    TClonesArray *particles = new TClonesArray("TParticle", 1000);

    // Load the events
    TTree *tEvents = (TTree *)inFile->Get("tEvents");
    tEvents->SetBranchStatus("particles", 1);
    tEvents->SetBranchAddress("particles", &particles);

    std::vector<FemtoParticle> partCharm{};
    std::vector<FemtoParticle> partLight{};
    std::deque<std::vector<FemtoParticle>> partBufferLight{};

    for (int iEvent = 0; iEvent < tEvents->GetEntries(); iEvent++) {
        tEvents->GetEntry(iEvent);

        partCharm.clear();
        partLight.clear();

        for (int iPart = 2; iPart < particles->GetEntriesFast(); iPart++) {
            TParticle *particle = dynamic_cast<TParticle *>(particles->At(iPart));

            if (std::abs(particle->Y()) > 2) continue;      // keep only midrapidity
            if (ProductionRadius(particle) > 1.) continue;  // to remove products of strange decays

            int pdg = particle->GetPdgCode();
            int absPdg = std::abs(pdg);

            if (absPdg == cpdg) {
                ROOT::Math::PxPyPzMVector part(particle->Px(), particle->Py(), particle->Pz(),
                                               TDatabasePDG::Instance()->GetParticle(absPdg)->Mass());
                partCharm.push_back({part, pdg, iPart, {particle->GetFirstDaughter(), particle->GetLastDaughter()}});
            } else if (absPdg == lpdg) {
                ROOT::Math::PxPyPzMVector part(particle->Px(), particle->Py(), particle->Pz(),
                                               TDatabasePDG::Instance()->GetParticle(absPdg)->Mass());
                partLight.push_back({part, pdg, iPart, {0,0}});
            }
        }

        if (align) {
            // Perform alignment of the events
            std::array<std::array<double, 3>, 3> rot = GetRotationSpheri3DFull(particles);
            partLight = Rotate(rot, partLight);
            partCharm = Rotate(rot, partCharm);
        }

        // same event
        for (const auto &charm : partCharm) {
            for (const auto &light : partLight) {
                if (cleanPairs && !IsPairClean(charm, light)) continue;
                double kStar = ComputeKstar(charm.p, light.p);
                std::string pair = charm.pdg * light.pdg > 0 ? "sc" : "oc";
                hSE[pair]->Fill(kStar);
            }
        }

        partBufferLight.push_back(partLight);

        // mixed event
        if (partBufferLight.size() < 2)  // to avoid repetitions
            continue;

        for (const auto & charm : partCharm) {
            for (size_t iME = 0; iME < partBufferLight.size() - 1; iME++) { // from 0 to last-1
                for (const auto & light : partBufferLight[iME]) {
                    double kStar = ComputeKstar(charm.p, light.p);
                    std::string pair = charm.pdg * light.pdg > 0 ? "sc" : "oc";

                    hME[pair]->Fill(kStar);
                }
            }
        }
        if (partBufferLight.size() > md) partBufferLight.pop_front();
    }

    TFile *oFile = new TFile(oFileName.data(), "recreate");

    for (auto pair : {"sc", "oc"}) {
        oFile->mkdir(pair);
        oFile->cd(pair);
        hSE[pair]->Write("hSE");
        hME[pair]->Write("hME");
    }
}