// Tagging-muon kinematics in pp calo jets: data vs PYTHIA.
//
// Companion to plotMuTaggedJetKinematics_pp_caloJets_dataVsMC.C: the muon that
// tags each jet, for the same jets. Five observables, all booked identically by
// both scans as 2D maps against the tagged jet's pT, filled only when the event
// fired the muon trigger:
//   muon pT      data h_mupt_recoJetPt_inclRecoMuonTag_triggerOn
//   muon eta          h_mueta_recoJetPt_inclRecoMuonTag_triggerOn
//   muon phi          h_muphi_recoJetPt_inclRecoMuonTag_triggerOn
//   ptRel             h_muptrel_recoJetPt_inclRecoMuonTag_triggerOn
//   dR(mu, jet)       h_muJetDr_recoJetPt
// PYTHIA books the same names with a flavour and template index appended;
// allJets_T0 is every tagged jet whatever its flavour at the nominal jet pT, so
// the broken flavour labels in the current PYTHIA calo file do not enter.
//
// JET SELECTION: every observable is drawn for pT >= jetPtMin (50 GeV), the
// region the b-fraction fits use. The pp calo data has no jets below 50 GeV at
// all, and the PYTHIA sample (pThat > 15) is sculpted there. ptRel and dR --
// the two that depend on the jet axis, and the ones that disagree inclusively
// -- are also drawn in each analysis jet pT bin (edgeJetPt), where the
// inclusive selection includes overflow above 500 GeV and the sliced ones stop
// at 500, as the b-fraction fits do.
//
// SHAPES ONLY, each normalized to unit area over the range it is drawn in:
// the MC's pThat weights have no counterpart in the data's counts. Overflow is
// negligible (muon pT > 100 GeV: 0.02% of data; ptRel > 10 GeV: 0.2%), but
// ptRel is drawn over the 0-5 GeV fit window, so its normalization excludes the
// tail beyond 5 GeV on purpose.
//
// Means and RMS are taken from the unrebinned projection restricted to the
// drawn range, before the variable-width rebin and bin-width division, which
// would otherwise shift them.
//
// Usage: root -l -b -q 'plotTagMuonKinematics_pp_caloJets_dataVsMC.C'
// Run from: src/plots/muonPt/

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

const char *outDir = "../../../figures/muonKinematics/pp_caloJets/";

const double jetPtMin = 50.;

// standard analysis jet pT bins, for the sliced ptRel and dR comparisons
const double edgeJetPt[] = {50, 60, 80, 100, 120, 150, 200, 300, 500};
const int    nEdgeJetPt  = sizeof(edgeJetPt)/sizeof(double) - 1;

// muon pT: 1 GeV input bins, widening where the spectrum thins
const double edgeMuPt[] = {15, 17, 20, 25, 30, 35, 40, 50, 60, 80, 100};
const int    nEdgeMuPt  = sizeof(edgeMuPt)/sizeof(double) - 1;

struct Observable {
  const char *tag;        // output name stem
  const char *dataHist;
  const char *mcHist;
  const char *xTitle;
  const char *yTitle;
  double      xMin, xMax; // drawn and normalized range
  int         rebin;      // fixed rebin factor; 0 = use edgeMuPt
  bool        logY;
  bool        sliceJetPt; // also draw in each edgeJetPt bin
};

const Observable obs[] = {
  {"tagMuonPt",  "h_mupt_recoJetPt_inclRecoMuonTag_triggerOn",
                 "h_mupt_recoJetPt_inclRecoMuonTag_triggerOn_allJets_T0",
                 "#it{p}_{T}^{#mu} [GeV]", "1/N dN/d#it{p}_{T}^{#mu} [GeV^{-1}]", 15., 100., 0, true, false},
  {"tagMuonEta", "h_mueta_recoJetPt_inclRecoMuonTag_triggerOn",
                 "h_mueta_recoJetPt_inclRecoMuonTag_triggerOn_allJets_T0",
                 "#it{#eta}^{#mu}", "1/N dN/d#it{#eta}^{#mu}", -2.0, 2.0, 5, false, false},           // 0.02 -> 0.1
  {"tagMuonPhi", "h_muphi_recoJetPt_inclRecoMuonTag_triggerOn",
                 "h_muphi_recoJetPt_inclRecoMuonTag_triggerOn_allJets_T0",
                 "#it{#phi}^{#mu}", "1/N dN/d#it{#phi}^{#mu}", -TMath::Pi(), TMath::Pi(), 2, false, false},
  {"tagMuonPtRel", "h_muptrel_recoJetPt_inclRecoMuonTag_triggerOn",
                 "h_muptrel_recoJetPt_inclRecoMuonTag_triggerOn_allJets_T0",
                 "#it{p}_{T}^{rel} [GeV]", "1/N dN/d#it{p}_{T}^{rel} [GeV^{-1}]", 0., 5., 2, false, true},  // 0.1 -> 0.2 GeV
  {"tagMuonJetDr", "h_muJetDr_recoJetPt",
                 "h_muJetDr_recoJetPt_allJets_T0",
                 "#Delta#it{R}(#mu, jet)", "1/N dN/d#Delta#it{R}", 0., 0.4, 8, false, true},             // 0.0025 -> 0.02
};
const int NObs = sizeof(obs)/sizeof(Observable);

struct Summary { double mean = 0., rms = 0., n = 0.; };

// Projection onto the observable for tagged jets with ptLo <= jet pT < ptHi
// (ptHi <= 0: no upper limit, overflow included). The summary -- mean and RMS
// over [xMin, xMax], and the weighted jet count -- is taken before any
// rebinning. Fails if a jet pT edge is not a bin boundary of the map, since the
// window would then silently widen or narrow.
static TH1D* project(TFile *f, const char *hname, const Observable &o, const char *tag,
                     double ptLo, double ptHi, Summary &s)
{
  TH2D *H = nullptr; f->GetObject(hname, H);
  if(!H){ printf("ERROR: %s missing in %s file\n", hname, tag); return nullptr; }
  TAxis *ay = H->GetYaxis();
  int by1 = ay->FindBin(ptLo + 1e-6);
  int by2 = (ptHi > 0.) ? ay->FindBin(ptHi - 1e-6) : H->GetNbinsY() + 1;
  if(fabs(ay->GetBinLowEdge(by1) - ptLo) > 1e-6 || (ptHi > 0. && fabs(ay->GetBinUpEdge(by2) - ptHi) > 1e-6)){
    printf("ERROR: %s: jet pT %.0f-%.0f GeV is not on its bin boundaries\n", hname, ptLo, ptHi); return nullptr;
  }
  TString sfx = Form("%s_%s_%.0f_%.0f", o.tag, tag, ptLo, ptHi);
  TH1D *p = H->ProjectionX(Form("p_%s", sfx.Data()), by1, by2);
  p->SetDirectory(nullptr);

  p->GetXaxis()->SetRangeUser(o.xMin + 1e-6, o.xMax - 1e-6);
  s.mean = p->GetMean(); s.rms = p->GetRMS();
  p->GetXaxis()->SetRange(0, 0);

  TH1D *r = nullptr;
  if(o.rebin > 0){ p->Rebin(o.rebin); r = p; }
  else{ r = rebinTo(p, nEdgeMuPt, edgeMuPt, Form("pr_%s", sfx.Data())); delete p; }

  s.n = r->Integral(r->FindBin(o.xMin + 1e-6), r->FindBin(o.xMax - 1e-6));
  if(s.n <= 0.){ printf("ERROR: %s empty in range for %s, jet pT %.0f-%.0f\n", hname, tag, ptLo, ptHi); delete r; return nullptr; }
  r->Scale(1./s.n);
  divideByBinwidth(r);
  return r;
}

void plotTagMuonKinematics_pp_caloJets_dataVsMC()
{
  initPlotStyle();
  gSystem->mkdir(outDir, kTRUE);

  TFile *fD = TFile::Open(dataFile), *fM = TFile::Open(mcFile);
  if(!fD || fD->IsZombie()){ printf("ERROR: cannot open %s\n", dataFile); return; }
  if(!fM || fM->IsZombie()){ printf("ERROR: cannot open %s\n", mcFile); return; }

  const TString line1 = "pp 5.02 TeV, anti-#it{k}_{T} #it{R} = 0.4 calo jets";
  const TString line2 = "#mu-tagged (#it{p}_{T}^{#mu} > 15 GeV), HLT_HIL3Mu12";

  // ---- inclusive, jet pT > jetPtMin ----------------------------------------
  printf("\ntag-muon kinematics, jets with pT > %.0f GeV (means/RMS over the drawn range)\n", jetPtMin);
  printf("%-13s %10s %10s %10s %10s %9s\n", "observable", "mean data", "mean MC", "RMS data", "RMS MC", "chi2/ndf");

  for(int i = 0; i < NObs; i++){
    const Observable &o = obs[i];
    Summary sD, sM;
    TH1D *hD = project(fD, o.dataHist, o, "data", jetPtMin, -1., sD);
    TH1D *hM = project(fM, o.mcHist,   o, "mc",   jetPtMin, -1., sM);
    if(!hD || !hM) continue;

    double chi2 = drawShapeComparison(hD, hM, o.xTitle, o.yTitle, o.logY, o.xMin, o.xMax,
                                      {line1, line2, Form("#it{p}_{T}^{jet} > %.0f GeV", jetPtMin)},
                                      Form("%s%s_dataVsMC.pdf", outDir, o.tag));
    printf("%-13s %10.4f %10.4f %10.4f %10.4f %9.2f\n", o.tag, sD.mean, sM.mean, sD.rms, sM.rms, chi2);
  }

  // ---- ptRel and dR in each analysis jet pT bin ----------------------------
  for(int i = 0; i < NObs; i++){
    const Observable &o = obs[i];
    if(!o.sliceJetPt) continue;

    printf("\n%s in jet pT bins\n", o.tag);
    printf("  %-9s %9s %10s %10s %9s %9s\n", "jet pT", "N data", "mean data", "mean MC", "data/MC", "chi2/ndf");

    for(int k = 0; k < nEdgeJetPt; k++){
      double lo = edgeJetPt[k], hi = edgeJetPt[k+1];
      Summary sD, sM;
      TH1D *hD = project(fD, o.dataHist, o, "data", lo, hi, sD);
      TH1D *hM = project(fM, o.mcHist,   o, "mc",   lo, hi, sM);
      if(!hD || !hM) continue;

      double chi2 = drawShapeComparison(hD, hM, o.xTitle, o.yTitle, o.logY, o.xMin, o.xMax,
                                        {line1, line2, Form("%.0f < #it{p}_{T}^{jet} < %.0f GeV", lo, hi)},
                                        Form("%s%s_jetPt%.0f-%.0f_dataVsMC.pdf", outDir, o.tag, lo, hi));
      printf("  %3.0f-%-5.0f %9.0f %10.4f %10.4f %9.3f %9.2f\n", lo, hi, sD.n, sD.mean, sM.mean,
             sM.mean > 0. ? sD.mean/sM.mean : -1., chi2);
    }
  }

  printf("\nfigures in %s\n", outDir);
}
