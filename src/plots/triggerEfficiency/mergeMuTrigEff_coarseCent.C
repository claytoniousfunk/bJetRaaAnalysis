// Sum the per-muon HLT_HIL3Mu12 efficiency histograms (h_muTrigEff_*, filled by
// PbPb_scan.C since 2026-09-18) over several MinBias scans and merge the
// ultra-fine 5% centrality slices into the nominal classes, so the result can be
// read by plotMuonTriggerEfficiency.C like any nominal-binned scan.
//
// Ultra-fine scheme (centrality_ultraFineCentBins.h): C0 = 0-90% inclusive,
// C1..C18 = 5% slices of 10 hiBin units. Nominal classes:
//   C1 0-10%  = slices  1- 2
//   C2 10-30% = slices  3- 6
//   C3 30-50% = slices  7-10
//   C4 50-80% = slices 11-16
//   C0 0-80%  = slices  1-16   (NOT the scan's own C0, which runs to 90%)
// Counts are summed; with prescale-1, unit-weight fills, pass/all of a merged
// class is the binomial efficiency of that class.
//
// Usage: root -l -b -q mergeMuTrigEff_coarseCent.C

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TSystem.h"
#include <vector>
#include <string>
#include <cstdio>

void mergeMuTrigEff_coarseCent()
{
  const std::string dir = "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/";
  const std::vector<std::string> files = {
    "PbPb_MinBias_Part1_caloJets_manualJEC_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_2026-9-19_ultraFineCentBins.root",
    "PbPb_MinBias_Part2_manualJEC_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_2026-9-20_ultraFineCentBins.root",
    "PbPb_MinBias_Part3_manualJEC_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_2026-9-21_ultraFineCentBins.root",
    "PbPb_MinBias_Part4_manualJEC_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_2026-9-21_ultraFineCentBins.root"};
  const std::string outDir  = "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/triggerEfficiency/";
  const std::string outFile = outDir + "PbPb_MinBias_Parts1-4_muTrigEff_nominalCentBins_2026-9-21.root";

  const int nClass = 5;
  const int sliceLo[nClass] = {1, 1, 3, 7, 11};
  const int sliceHi[nClass] = {16, 2, 6, 10, 16};
  const char *label[nClass] = {"0-80%", "0-10%", "10-30%", "30-50%", "50-80%"};
  const std::vector<std::string> names2D = {"h_muTrigEff_all", "h_muTrigEff_pass",
                                            "h_muTrigEff_all_noVeto", "h_muTrigEff_pass_noVeto"};

  std::vector<TFile*> in;
  for (auto &f : files) {
    TFile *F = TFile::Open((dir + f).c_str());
    if (!F || F->IsZombie()) { printf("cannot open %s\n", f.c_str()); return; }
    in.push_back(F);
  }

  gSystem->mkdir(outDir.c_str(), true);
  TFile out(outFile.c_str(), "RECREATE");

  for (int c = 0; c < nClass; c++) {
    for (auto &n : names2D) {
      TH2D *sum = nullptr;
      for (TFile *F : in)
        for (int s = sliceLo[c]; s <= sliceHi[c]; s++) {
          TH2D *h = (TH2D*) F->Get(Form("%s_C%d", n.c_str(), s));
          if (!h) { printf("missing %s_C%d in %s\n", n.c_str(), s, F->GetName()); return; }
          if (!sum) {
            sum = (TH2D*) h->Clone(Form("%s_C%d", n.c_str(), c));
            sum->SetDirectory(nullptr);
            sum->SetTitle(Form("%s, %s, MinBias Parts 1-4", n.c_str(), label[c]));
          } else sum->Add(h);
        }
      out.cd(); sum->Write();
      if (n == "h_muTrigEff_all")
        printf("  %-7s slices %2d-%2d  tight probes pT>15 |eta|<2: %.0f\n", label[c], sliceLo[c], sliceHi[c],
               sum->Integral(sum->GetXaxis()->FindBin(15 + 1e-6), sum->GetNbinsX() + 1,
                             sum->GetYaxis()->FindBin(-2 + 1e-6), sum->GetYaxis()->FindBin(2 - 1e-6)));
      delete sum;
    }
    TH1D *ps = nullptr;
    for (TFile *F : in)
      for (int s = sliceLo[c]; s <= sliceHi[c]; s++) {
        TH1D *h = (TH1D*) F->Get(Form("h_muTrigEff_mu12Prescale_C%d", s));
        if (!ps) { ps = (TH1D*) h->Clone(Form("h_muTrigEff_mu12Prescale_C%d", c)); ps->SetDirectory(nullptr); }
        else ps->Add(h);
      }
    out.cd(); ps->Write(); delete ps;
  }
  out.Close();
  printf("wrote %s\n", outFile.c_str());
}
