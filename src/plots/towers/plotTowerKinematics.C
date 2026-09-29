// Calorimeter tower pT (et), eta and phi from rechitanalyzerpp/tower.
//
// The branches are per-event arrays, et[n] / eta[n] / phi[n], so every tower in
// every event is filled. Distributions are per event (divided by the number of
// events) and per unit of the x variable, so files with different statistics
// compare directly. et is the transverse energy of the tower, used here as its
// pT.
//
// pT is drawn on a log y axis: the spectrum falls by ~4 decades over 0-20 GeV, so
// a linear axis would show only the first few bins. eta and phi are linear.
// eta is not flat: the tower grid is 0.087 wide in the barrel and coarser and
// irregular forward, so counts per eta bin follow the granularity, and the
// 0.5-wide bins are needed to smooth the forward
// aliasing (0.25 bins still showed it).
//
// The E_T spectrum is drawn straight from the et branch, E_T = E/cosh(eta).
// phi uses 18 bins (0.35 rad): the calorimeter has 72 phi segments barrel-side and
// 36/18 forward, so 18 bins hold a whole number of towers each; 16 aliased.
//
// pT is constructed from the tower variables rather than read from et: treating
// the tower as massless, p_x = E_T cos(phi), p_y = E_T sin(phi), and
// pT = sqrt(p_x^2 + p_y^2), which reduces to E_T. The largest difference from
// E/cosh(eta) is printed as a consistency check.
//
// Each of pT, E_T and phi is also drawn with |eta| < 1.6, the jet acceptance
// (etaMax in headers/AnalysisSetup/common.h), as *_etaTrunc.pdf. This is the
// same tower selection a |eta| < 1.6 jet can actually draw on, so it is the
// more relevant comparison to jet-level distributions than the full-eta ones
// above, which are dominated by the endcap and forward region (see the
// eta-slice fraction-stack plots in this directory).
//
// Usage: root -l -b -q 'plotTowerKinematics.C("/path/to/HiForestAOD.root")'
// Run from: src/plots/towers/

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TSystem.h"
#include "TMath.h"
#include "../../../headers/plotting/plotStyle.h"

namespace {

const char *outDir = "../../../figures/towers/";

void drawOne(TH1D *h, const char *xTitle, const char *yTitle, bool logY,
             const TString &tag, const char *outName)
{
  TCanvas *c = new TCanvas(Form("c_%s", outName), "", 800, 700);
  c->SetLeftMargin(0.15); c->SetBottomMargin(0.13); c->SetRightMargin(0.05);
  c->SetLogy(logY);

  styleH(h, hexData, markFilledCircle, 0.9);
  h->GetXaxis()->SetTitle(xTitle);
  h->GetYaxis()->SetTitle(yTitle);
  h->GetXaxis()->SetTitleSize(0.05); h->GetYaxis()->SetTitleSize(0.05);
  h->GetXaxis()->SetLabelSize(0.042); h->GetYaxis()->SetLabelSize(0.042);
  h->GetYaxis()->SetTitleOffset(1.4);
  if(logY) h->SetMinimum(h->GetMinimum(0.) * 0.5);
  else     h->SetMinimum(0.);
  h->SetMaximum(h->GetMaximum() * (logY ? 5. : 1.4));
  h->Draw("E");

  TLatex l; l.SetNDC(); l.SetTextSize(0.04);
  l.DrawLatex(0.19, 0.86, "Calorimeter towers (rechitanalyzerpp)");
  l.DrawLatex(0.19, 0.80, tag);

  savePdfTight(c, Form("%s%s.pdf", outDir, outName));
  delete c;
}

} // namespace

void plotTowerKinematics(const char *path = "/home/clayton/Downloads/HiForestAOD.root")
{
  initPlotStyle();
  gSystem->mkdir(outDir, kTRUE);

  TFile *f = TFile::Open(path);
  if(!f || f->IsZombie()){ printf("ERROR: %s did not open\n", path); return; }
  TTree *t = nullptr; f->GetObject("rechitanalyzerpp/tower", t);
  if(!t){ printf("ERROR: rechitanalyzerpp/tower missing\n"); return; }

  const Long64_t nEvt = t->GetEntries();
  const int maxN = 20000;
  Int_t n = 0;
  static Float_t e[maxN], et[maxN], eta[maxN], phi[maxN];
  t->SetBranchAddress("n", &n);
  t->SetBranchAddress("e", e);
  t->SetBranchAddress("et", et);
  t->SetBranchAddress("eta", eta);
  t->SetBranchAddress("phi", phi);

  const double etaMaxCut = 1.6;   // jet acceptance, headers/AnalysisSetup/common.h

  // 0.25 GeV bins to 45 GeV: the observed maximum is 43.5 GeV
  TH1D *hPt  = new TH1D("hTowerEt",  "", 180, 0., 45.);
  TH1D *hEta = new TH1D("hTowerEta", "", 20, -5.0, 5.0);
  TH1D *hEt  = new TH1D("hTowerEtRaw", "", 180, 0., 45.);
  TH1D *hPhi = new TH1D("hTowerPhi", "", 18, -TMath::Pi(), TMath::Pi());
  TH1D *hPtT = new TH1D("hTowerEt_trunc",  "", 180, 0., 45.);
  TH1D *hEtT = new TH1D("hTowerEtRaw_trunc", "", 180, 0., 45.);
  TH1D *hPhiT = new TH1D("hTowerPhi_trunc", "", 18, -TMath::Pi(), TMath::Pi());
  for(auto *h : {hPt, hEt, hEta, hPhi, hPtT, hEtT, hPhiT}){ h->Sumw2(); h->SetDirectory(nullptr); }

  Long64_t nTowers = 0, nTowersTrunc = 0;
  double maxDiff = 0.;
  for(Long64_t i = 0; i < nEvt; i++){
    t->GetEntry(i);
    if(n > maxN){ printf("ERROR: event %lld has %d towers, raise maxN\n", i, n); return; }
    for(int j = 0; j < n; j++){
      const double px = et[j]*TMath::Cos(phi[j]), py = et[j]*TMath::Sin(phi[j]);
      const double pt = TMath::Sqrt(px*px + py*py);
      const double d = TMath::Abs(pt - e[j]/TMath::CosH(eta[j]));
      if(d > maxDiff) maxDiff = d;
      hPt->Fill(pt); hEt->Fill(et[j]); hEta->Fill(eta[j]); hPhi->Fill(phi[j]);
      if(TMath::Abs(eta[j]) < etaMaxCut){
        hPtT->Fill(pt); hEtT->Fill(et[j]); hPhiT->Fill(phi[j]);
        nTowersTrunc++;
      }
    }
    nTowers += n;
  }
  printf("  %lld events, %lld towers (%.1f per event), %lld with |eta| < %.1f (%.1f per event)\n",
         nEvt, nTowers, (double)nTowers/nEvt, nTowersTrunc, etaMaxCut, (double)nTowersTrunc/nEvt);
  printf("  max |pT - E/cosh(eta)| = %.2e GeV\n", maxDiff);
  printf("  mean pT %.4f GeV, mean eta %+.4f, mean phi %+.4f\n",
         hPt->GetMean(), hEta->GetMean(), hPhi->GetMean());

  // per event, per unit x
  for(auto *h : {hPt, hEt, hEta, hPhi, hPtT, hEtT, hPhiT}){ h->Scale(1./nEvt); h->Scale(1., "width"); }

  const TString tag = Form("%lld events", nEvt);
  drawOne(hEt,  "tower #it{E}_{T} [GeV]", "1/N_{evt} dN/d#it{E}_{T} [GeV^{-1}]", true, tag, "towerEt");
  drawOne(hPt,  "tower #it{p}_{T} [GeV]", "1/N_{evt} dN/d#it{p}_{T} [GeV^{-1}]", true,  tag, "towerPt");
  drawOne(hEta, "tower #it{#eta}",        "1/N_{evt} dN/d#it{#eta}",             false, tag, "towerEta");
  drawOne(hPhi, "tower #it{#phi}",        "1/N_{evt} dN/d#it{#phi}",             false, tag, "towerPhi");

  const TString tagT = Form("%lld events, |#it{#eta}| < %.1f", nEvt, etaMaxCut);
  drawOne(hEtT,  "tower #it{E}_{T} [GeV]", "1/N_{evt} dN/d#it{E}_{T} [GeV^{-1}]", true, tagT, "towerEt_etaTrunc");
  drawOne(hPtT,  "tower #it{p}_{T} [GeV]", "1/N_{evt} dN/d#it{p}_{T} [GeV^{-1}]", true,  tagT, "towerPt_etaTrunc");
  drawOne(hPhiT, "tower #it{#phi}",        "1/N_{evt} dN/d#it{#phi}",             false, tagT, "towerPhi_etaTrunc");
  printf("  figures in %s\n", outDir);
}
