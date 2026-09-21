// Effect of a 2 GeV PF-candidate pT cut on same-event PbPb MinBias spectra, in
// coarse centrality classes.
//
// The cut (pseudoJetCandPt_min in pseudoJets.h) drops PF candidates below
// 2 GeV from the random cones and from the FastJet clustering inputs. Compared
// here, per coarse class:
//   h_fastJetPt_PF      FastJet anti-kT R=0.4 PF jets, raw pT, |eta| < 1.6,
//                       after the jetTrkMax filter
//   h_fastJetPt_PF_JEC  the same jets, JEC applied to the raw pT
//   h_pseudoJetPt       random-cone pT (R = 0.4, 100 throws per event) -- the
//                       quantity the UE maps are built from
//   h_inclRecoJetPt     forest akCs4PF reco jets. Built without the PF
//                       candidates, so the cut cannot touch it: a sanity check
//                       that the two scans ran over the same events and
//                       selection. It must come out identical.
//
// SAME EVENTS. Both scans ran over the same MinBias file list (8,291,975 events
// each, identical reco-jet entries), so the ratio isolates the cut. For the
// same reason the two spectra are strongly correlated, and the ratio panels
// carry the with-cut histogram's error only.
//
// INPUTS: the 2 GeV scan is from 2026-09-15 (git ea1f91eb, dirty tree). The
// no-cut baseline is the 2026-08-18 scan. No commit between the two changes how
// the same-event h_fastJetPt_PF, h_fastJetPt_PF_JEC or h_pseudoJetPt are
// filled (checked against git log). The 2026-08-17 scan is kept as a
// cross-check: it should agree with 2026-08-18 bin for bin in these histograms.
//
// NOT compared: the RC- and dPT-subtracted FastJet spectra. The 2 GeV scan
// subtracted maps made WITHOUT the cut (the 2026-08-17 file), so those spectra
// mix a 2 GeV jet with a 0 GeV background and mean nothing.
//
// NORMALISATION: FastJet and reco-jet spectra per event (h_vz, ungated), per
// GeV. Random cones per throw (events x 100), per GeV.
//
// Usage: root -l -b -q 'plotPFCandPtCutEffect_sameEvent_coarseCent.C'
// Run from: src/plots/pfCand/

#include <cmath>
#include "TMath.h"
#include "TH1D.h"
#include "TFile.h"
#include "TSystem.h"
#include "../../../headers/functions/divideByBinwidth.h"
#include "../../../headers/plotting/coarseCent.h"
#include "../../../headers/plotting/plotStyle.h"
#include "../../../headers/plotting/ratioPanel.h"
#include "../../../headers/plotting/shapeComparison.h"

const char *dir = "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/";
const char *fileCut   = "PbPb_MinBias_Part1_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_sameEventPFClustering_pseudoJetCandPtMin-2.0_2026-9-15_ultraFineCentBins.root";
const char *fileNoCut = "PbPb_MinBias_Part1_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_sameEventPFClustering_pseudoJetCandPtMin-0.0_2026-8-18_ultraFineCentBins.root";
const char *fileCheck = "PbPb_MinBias_Part1_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_sameEventPFClustering_pseudoJetCandPtMin-0.0_2026-8-17_ultraFineCentBins.root";

const char *outDir = "../../../figures/pfCandPtCut/";

const int    NThrowsPerEvent = 100;   // N_mixedEventsInPool random cones per event

// jet pT: fine where the cut matters most, wide in the tail
const double edgeJetPt[] = {0, 5, 10, 15, 20, 25, 30, 40, 50, 60, 80, 100, 120, 150, 200, 300, 500};
const int    nEdgeJetPt  = sizeof(edgeJetPt)/sizeof(double) - 1;

struct Spectrum {
  const char *hist;
  const char *tag;
  const char *xTitle;
  const char *yTitle;
  bool        perThrow;   // random cones: normalize per throw, not per event
  bool        logY;
  double      xMin, xMax;
  bool        variableBins;
};

const Spectrum spectra[] = {
  {"h_fastJetPt_PF",     "fastJetPtRaw", "raw #it{p}_{T}^{FastJet} [GeV]", "1/N_{evt} dN/d#it{p}_{T} [GeV^{-1}]", false, true,  0., 500., true},
  {"h_fastJetPt_PF_JEC", "fastJetPtJEC", "#it{p}_{T}^{FastJet} (JEC) [GeV]", "1/N_{evt} dN/d#it{p}_{T} [GeV^{-1}]", false, true,  0., 500., true},
  // log y: with the cut the distribution collapses into the first few bins, and
  // on a linear axis the peripheral classes show nothing
  {"h_pseudoJetPt",      "randomConePt", "random-cone #it{p}_{T} [GeV]", "1/N_{cone} dN/d#it{p}_{T} [GeV^{-1}]", true,  true,  0., 200., false},
};
const int NSpectra = sizeof(spectra)/sizeof(Spectrum);

// coarse-class spectrum, per event (or per throw), per GeV
static TH1D* spectrum(TFile *f, const Spectrum &s, int ci, const char *tag, double &mean)
{
  TH1D *h = coarseSum1(f, s.hist, ci, tag);
  if(!h){ printf("ERROR: %s missing for %s in %s\n", s.hist, coarseLabel[ci], tag); return nullptr; }
  double nEvt = coarseEvents(f, ci, "h_vz");
  if(nEvt <= 0.){ printf("ERROR: no events for %s in %s\n", coarseLabel[ci], tag); delete h; return nullptr; }
  mean = h->GetMean();
  TH1D *r = h;
  if(s.variableBins){ r = rebinTo(h, nEdgeJetPt, edgeJetPt, Form("%s_rb", h->GetName())); delete h; }
  r->Scale(1./(nEvt * (s.perThrow ? NThrowsPerEvent : 1)));
  divideByBinwidth(r);
  return r;
}

// largest |a/b - 1| over populated bins
static double maxRelDiff(TH1D *a, TH1D *b)
{
  double m = 0.;
  for(int i = 1; i <= a->GetNbinsX(); i++)
    if(b->GetBinContent(i) > 0.) m = TMath::Max(m, fabs(a->GetBinContent(i)/b->GetBinContent(i) - 1.));
  return m;
}

void plotPFCandPtCutEffect_sameEvent_coarseCent()
{
  initPlotStyle();
  gSystem->mkdir(outDir, kTRUE);

  TFile *fC = TFile::Open(Form("%s%s", dir, fileCut));
  TFile *fN = TFile::Open(Form("%s%s", dir, fileNoCut));
  TFile *fX = TFile::Open(Form("%s%s", dir, fileCheck));
  if(!fC || fC->IsZombie() || !fN || fN->IsZombie()){ printf("ERROR: cannot open the input scans\n"); return; }

  // ---- sanity checks ------------------------------------------------------
  printf("\nsanity checks (largest per-bin relative difference, coarse classes)\n");
  for(int ci = 0; ci < NCoarse; ci++){
    double nC = coarseEvents(fC, ci), nN = coarseEvents(fN, ci);
    TH1D *rC = coarseSum1(fC, "h_inclRecoJetPt", ci, "sanC"), *rN = coarseSum1(fN, "h_inclRecoJetPt", ci, "sanN");
    printf("  %-7s events cut/noCut %.0f / %.0f   reco jets (must be 0): %.2e",
           coarseLabel[ci], nC, nN, (rC && rN) ? maxRelDiff(rC, rN) : -1.);
    if(fX && !fX->IsZombie()){
      TH1D *a = coarseSum1(fN, "h_fastJetPt_PF", ci, "chkN"), *b = coarseSum1(fX, "h_fastJetPt_PF", ci, "chkX");
      printf("   FastJet PF 08-18 vs 08-17: %.2e", (a && b) ? maxRelDiff(a, b) : -1.);
      delete a; delete b;
    }
    printf("\n");
    delete rC; delete rN;
  }

  // ---- comparisons ----------------------------------------------------------
  for(int is = 0; is < NSpectra; is++){
    const Spectrum &s = spectra[is];
    printf("\n%s\n", s.hist);
    if(s.variableBins){
      printf("  %-7s %9s %9s | 2 GeV / no cut, per-event yield in raw or JEC pT windows\n", "class", "<pT> cut", "<pT> 0");
      printf("  %-7s %9s %9s |", "", "", "");
      for(double lo : {20., 50., 100.}) printf("   >%3.0f GeV", lo);
      printf("\n");
    }
    else printf("  %-7s %9s %9s %9s\n", "class", "<pT> cut", "<pT> 0", "shift");

    for(int ci = 0; ci < NCoarse; ci++){
      double mC = 0., mN = 0.;
      TH1D *hC = spectrum(fC, s, ci, "cut", mC);
      TH1D *hN = spectrum(fN, s, ci, "noCut", mN);
      if(!hC || !hN) continue;

      if(s.variableBins){
        printf("  %-7s %9.2f %9.2f |", coarseLabel[ci], mC, mN);
        for(double lo : {20., 50., 100.}){
          int b1 = hC->FindBin(lo + 1e-6);
          // integrals are per GeV now, so weight by width
          double yC = hC->Integral(b1, hC->GetNbinsX(), "width"), yN = hN->Integral(b1, hN->GetNbinsX(), "width");
          printf(" %10.3f", yN > 0. ? yC/yN : -1.);
        }
        printf("\n");
      }
      else printf("  %-7s %9.2f %9.2f %+9.2f\n", coarseLabel[ci], mC, mN, mC - mN);

      // Random cones: with the cut the two distributions barely overlap in
      // central events, so the ratio mostly shows where they cross; the mean
      // shift in the printed table is the number to quote.
      double rMaxUse = s.variableBins ? 1.6 : 5.;
      drawShapeComparison(hC, hN, s.xTitle, s.yTitle, s.logY, s.xMin, s.xMax,
                          {"PbPb 5.02 TeV MinBias, same-event PF candidates",
                           Form("%s, |#it{#eta}| < 1.6", coarseLabel[ci]),
                           s.perThrow ? "random cones, #it{R} = 0.4"
                                      : "anti-#it{k}_{T} #it{R} = 0.4 PF jets"},
                          Form("%s%s_%s_pfCandPtCut.pdf", outDir, s.tag, coarseTag[ci]),
                          "#it{p}_{T}^{cand} > 2 GeV", "no cut",
                          "2 GeV / no cut", 0., rMaxUse, RatioErr::kNumerator);
    }
  }

  printf("\nfigures in %s\n", outDir);
}
