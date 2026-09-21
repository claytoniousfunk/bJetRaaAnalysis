// Muon-tagged calo jet kinematics in pp: data vs PYTHIA.
//
// Flavour-independent check of the calo scans while the PYTHIA flavour labels
// are being fixed. Everything comes from the tagged + trigger-on 2D maps both
// scans book identically, jet pT vs eta and jet pT vs phi, so data and MC go
// through the same projection:
//   data  h_inclRecoJetPt_inclRecoJetEta_inclRecoMuonTag_triggerOn   (pp_scan)
//   MC    same name                                                   (PYTHIA_scan)
// None of these uses the jet flavour, so the broken labels in the current
// PYTHIA calo file do not affect them.
//
// SHAPES ONLY. The MC carries pThat cross-section weights with no counterpart
// in the data's counts, so every distribution is normalized to unit area:
//   pT       over [normPtMin, 500] GeV. The PYTHIA sample is pThat > 15, which
//            sculpts the MC spectrum at low jet pT; normalizing there would let
//            that region set the scale for the whole comparison.
//   eta/phi  for jets with pT >= etaPhiPtMin, the region the b-fraction fits
//            use, for the same reason.
//
// The pT spectrum is drawn on a log axis -- it falls by orders of magnitude
// over the range -- with a linear ratio panel. eta and phi are linear.
//
// The figure itself is drawShapeComparison() in headers/plotting/shapeComparison.h,
// shared with the tag-muon macro (plotTagMuonKinematics_pp_caloJets_dataVsMC.C).
//
// Usage: root -l -b -q 'plotMuTaggedJetKinematics_pp_caloJets_dataVsMC.C'
// Run from: src/plots/jetPt/

#include <cmath>
#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TSystem.h"
#include "../../../headers/functions/divideByBinwidth.h"
#include "../../../headers/plotting/plotStyle.h"
#include "../../../headers/plotting/ratioPanel.h"
#include "../../../headers/plotting/shapeComparison.h"

const char *dataFile =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/pp/"
  "pp_SingleMuon_caloJets_mu12_pTmu-15to999_tight_deltaR-40_jetTrkMaxFilter_WDecayFilter_2026-9-15.root";
const char *mcFile =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PYTHIA/"
  "PYTHIA_DiJet_caloJets_pThat-15_mu12_pTmu-15to999_tight_vzReweight_jetTrkMaxFilter_removeHYDJETjet0p35CutOnGen_WDecayFilter_2026-9-15.root";

const char *hPtEta = "h_inclRecoJetPt_inclRecoJetEta_inclRecoMuonTag_triggerOn";
const char *hPtPhi = "h_inclRecoJetPt_inclRecoJetPhi_inclRecoMuonTag_triggerOn";

const char *outDir = "../../../figures/jetKinematics/pp_caloJets/";

const double normPtMin   = 50.;   // pT normalization range starts here
const double etaPhiPtMin = 50.;   // eta and phi use jets above this

// Every edge sits on the scans' 5 GeV grid. The spectrum starts at 50 GeV
// because the pp calo data has no jets below it at all (first populated bin
// 50-55 GeV, tagged and inclusive alike), while PYTHIA reaches 30 GeV; bins
// below 50 would show a data/MC ratio of zero that says nothing about either.
const double edgePt[] = {50, 60, 70, 80, 90, 100, 120, 150, 200, 300, 500};
const int    nEdgePt  = sizeof(edgePt)/sizeof(double) - 1;
const int    etaRebin = 4;   // 160 bins of 0.02 -> 0.08
const int    phiRebin = 2;   // 100 bins -> 50

// integral of h over [lo, hi], bin-boundary safe
static double integralRange(TH1D *h, double lo, double hi)
{
  return h->Integral(h->FindBin(lo + 1e-6), h->FindBin(hi - 1e-6));
}

// pT spectrum: all eta, so every y bin including under/overflow
static TH1D* ptSpectrum(TFile *f, const char *tag)
{
  TH2D *H = nullptr; f->GetObject(hPtEta, H);
  if(!H){ printf("ERROR: %s missing in %s file\n", hPtEta, tag); return nullptr; }
  TH1D *p = H->ProjectionX(Form("pt_raw_%s", tag), 0, -1);
  p->SetDirectory(nullptr);
  TH1D *r = rebinTo(p, nEdgePt, edgePt, Form("pt_%s", tag));
  delete p;
  double n = integralRange(r, normPtMin, edgePt[nEdgePt]);
  if(n <= 0.){ printf("ERROR: %s has no jets above %.0f GeV\n", tag, normPtMin); delete r; return nullptr; }
  r->Scale(1./n);
  divideByBinwidth(r);
  return r;
}

// eta or phi for jets with pT >= etaPhiPtMin
static TH1D* angular(TFile *f, const char *hname, int rebin, const char *tag)
{
  TH2D *H = nullptr; f->GetObject(hname, H);
  if(!H){ printf("ERROR: %s missing in %s file\n", hname, tag); return nullptr; }
  int b1 = H->GetXaxis()->FindBin(etaPhiPtMin + 1e-6);
  TH1D *p = H->ProjectionY(Form("ang_%s_%s", hname, tag), b1, H->GetNbinsX() + 1);
  p->SetDirectory(nullptr);
  p->Rebin(rebin);
  if(p->Integral() <= 0.){ printf("ERROR: %s empty above %.0f GeV in %s\n", hname, etaPhiPtMin, tag); delete p; return nullptr; }
  p->Scale(1./p->Integral());
  divideByBinwidth(p);
  return p;
}

void plotMuTaggedJetKinematics_pp_caloJets_dataVsMC()
{
  initPlotStyle();
  gSystem->mkdir(outDir, kTRUE);

  TFile *fD = TFile::Open(dataFile), *fM = TFile::Open(mcFile);
  if(!fD || fD->IsZombie()){ printf("ERROR: cannot open %s\n", dataFile); return; }
  if(!fM || fM->IsZombie()){ printf("ERROR: cannot open %s\n", mcFile); return; }

  const TString line1 = "pp 5.02 TeV, anti-#it{k}_{T} #it{R} = 0.4 calo jets";
  const TString line2 = "#mu-tagged (#it{p}_{T}^{#mu} > 15 GeV), HLT_HIL3Mu12";

  // ---- pT ------------------------------------------------------------------
  TH1D *ptD = ptSpectrum(fD, "data"), *ptM = ptSpectrum(fM, "mc");
  if(ptD && ptM){
    double chi2 = drawShapeComparison(ptD, ptM, "#it{p}_{T}^{jet} [GeV]", "1/N dN/d#it{p}_{T} [GeV^{-1}]",
                                      true, edgePt[0], edgePt[nEdgePt],
                                      {line1, line2, Form("unit area over %.0f < #it{p}_{T}^{jet} < %.0f GeV", normPtMin, edgePt[nEdgePt])},
                                      Form("%smuTaggedJetPt_dataVsMC.pdf", outDir));
    printf("\npT: shapes normalized over %.0f-%.0f GeV, chi2/ndf %.2f\n", normPtMin, edgePt[nEdgePt], chi2);
    printf("  %-9s %10s\n", "jet pT", "data/MC");
    for(int b = 1; b <= ptD->GetNbinsX(); b++)
      printf("  %3.0f-%-5.0f %10.3f\n", edgePt[b-1], edgePt[b],
             ptM->GetBinContent(b) > 0. ? ptD->GetBinContent(b)/ptM->GetBinContent(b) : -1.);
  }

  // ---- eta, phi --------------------------------------------------------------
  TH1D *etaD = angular(fD, hPtEta, etaRebin, "data"), *etaM = angular(fM, hPtEta, etaRebin, "mc");
  if(etaD && etaM){
    double chi2 = drawShapeComparison(etaD, etaM, "#it{#eta}^{jet}", "1/N dN/d#it{#eta}",
                                      false, -1.6, 1.6, {line1, line2, Form("#it{p}_{T}^{jet} > %.0f GeV", etaPhiPtMin)},
                                      Form("%smuTaggedJetEta_dataVsMC.pdf", outDir));
    printf("\neta (pT > %.0f GeV): mean data %.4f  MC %.4f   RMS data %.4f  MC %.4f   chi2/ndf %.2f\n",
           etaPhiPtMin, etaD->GetMean(), etaM->GetMean(), etaD->GetRMS(), etaM->GetRMS(), chi2);
  }
  TH1D *phiD = angular(fD, hPtPhi, phiRebin, "data"), *phiM = angular(fM, hPtPhi, phiRebin, "mc");
  if(phiD && phiM){
    double chi2 = drawShapeComparison(phiD, phiM, "#it{#phi}^{jet}", "1/N dN/d#it{#phi}",
                                      false, -TMath::Pi(), TMath::Pi(), {line1, line2, Form("#it{p}_{T}^{jet} > %.0f GeV", etaPhiPtMin)},
                                      Form("%smuTaggedJetPhi_dataVsMC.pdf", outDir));
    printf("phi (pT > %.0f GeV): mean data %.4f  MC %.4f   chi2/ndf %.2f\n",
           etaPhiPtMin, phiD->GetMean(), phiM->GetMean(), chi2);
  }

  printf("\nfigures in %s\n", outDir);
}
