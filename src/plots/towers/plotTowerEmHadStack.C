// Stacked EM and hadronic tower E_T versus eta, phi and tower pT, from
// rechitanalyzerpp/tower.
//
// et = emEt + hadEt tower by tower (checked: the largest residual is float
// precision, ~1e-6 GeV), so stacking the two summed over towers reproduces the
// total exactly. That is why this is drawn as E_T summed per bin -- <sum E_T>
// per event per unit eta, phi or GeV of tower pT -- and not as the per-tower
// E_T spectrum: a spectrum counts towers per E_T value, and the EM and
// hadronic spectra are distributions of different quantities that do not add
// to the spectrum of the sum.
//
// towerEtStack_pt bins BY the tower's own pT (et), so it is an energy-flow
// spectrum (d<sum ET>/dpT), not a tower-count spectrum (dN/dpT, that is
// plotTowerKinematics.C's towerPt/towerEt) -- each tower enters at its own pT
// weighted by its own ET, so a bin's height is (roughly) its pT times the
// count there. It is the un-normalized numerator/denominator of the fraction
// in plotTowerPtFractionStack.C; this shows the same EM/hadronic split as an
// absolute energy scale instead of a fraction.
//
// Bins are the same as plotTowerKinematics.C (0.5 in eta, 18 in phi, 0.25 GeV
// to 20 in pT) so the figures line up. The black line is the total et branch,
// filled independently, and lies on top of the stack top.
//
// Usage: root -l -b -q 'plotTowerEmHadStack.C("/path/to/HiForestAOD.root")'
// Run from: src/plots/towers/

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TSystem.h"
#include "TMath.h"
#include "../../../headers/plotting/plotStyle.h"

namespace {

const char *outDir = "../../../figures/towers/";

void drawStack(TH1D *hEm, TH1D *hHad, TH1D *hTot, const char *xTitle, const char *yTitle,
               const TString &tag, const char *outName, bool logY = false)
{
  TCanvas *c = new TCanvas(Form("c_%s", outName), "", 800, 700);
  c->SetLeftMargin(0.15); c->SetBottomMargin(0.13); c->SetRightMargin(0.05);
  c->SetLogy(logY);

  const int cEm = TColor::GetColor("#E69F00"), cHad = TColor::GetColor("#0072B2");
  hEm->SetFillColor(cEm);   hEm->SetLineColor(cEm);
  hHad->SetFillColor(cHad); hHad->SetLineColor(cHad);
  styleLine(hTot, hexData, 3);

  THStack *st = new THStack(Form("st_%s", outName), "");
  st->Add(hHad); st->Add(hEm);   // hadronic at the bottom
  st->Draw("HIST");
  st->GetXaxis()->SetTitle(xTitle); st->GetYaxis()->SetTitle(yTitle);
  st->GetXaxis()->SetTitleSize(0.05); st->GetYaxis()->SetTitleSize(0.05);
  st->GetXaxis()->SetLabelSize(0.042); st->GetYaxis()->SetLabelSize(0.042);
  st->GetYaxis()->SetTitleOffset(1.4);
  st->SetMinimum(logY ? hTot->GetMinimum(0.) * 0.3 : 0.);
  st->SetMaximum(hTot->GetMaximum() * (logY ? 6. : 1.75));
  hTot->Draw("HIST SAME");

  TLegend *leg = makeLegend(0.60, 0.64, 0.93, 0.78);
  leg->AddEntry(hTot, "total #it{E}_{T} (et)", "l");
  leg->AddEntry(hEm,  "EM (emEt)", "f");
  leg->AddEntry(hHad, "hadronic (hadEt)", "f");
  leg->Draw();

  TLatex l; l.SetNDC(); l.SetTextSize(0.04);
  l.DrawLatex(0.19, 0.86, "Calorimeter towers (rechitanalyzerpp)");
  l.DrawLatex(0.19, 0.80, tag);
  l.SetTextSize(0.036);
  l.DrawLatex(0.19, 0.73, "#it{E}_{T} = #it{E}_{T}^{EM} + #it{E}_{T}^{had}");

  savePdfTight(c, Form("%s%s.pdf", outDir, outName));
  delete c;
}

} // namespace

void plotTowerEmHadStack(const char *path = "/home/clayton/Downloads/HiForestAOD.root")
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
  static Float_t et[maxN], emEt[maxN], hadEt[maxN], eta[maxN], phi[maxN];
  t->SetBranchAddress("n", &n);
  t->SetBranchAddress("et", et);
  t->SetBranchAddress("emEt", emEt);
  t->SetBranchAddress("hadEt", hadEt);
  t->SetBranchAddress("eta", eta);
  t->SetBranchAddress("phi", phi);

  const int nEta = 20, nPhi = 18, nPt = 180; const double ptMax = 45.;
  TH1D *eEm  = new TH1D("eEm",  "", nEta, -5., 5.), *eHad = new TH1D("eHad", "", nEta, -5., 5.),
       *eTot = new TH1D("eTot", "", nEta, -5., 5.);
  TH1D *pEm  = new TH1D("pEm",  "", nPhi, -TMath::Pi(), TMath::Pi()),
       *pHad = new TH1D("pHad", "", nPhi, -TMath::Pi(), TMath::Pi()),
       *pTot = new TH1D("pTot", "", nPhi, -TMath::Pi(), TMath::Pi());
  TH1D *tEm  = new TH1D("tEm",  "", nPt, 0., ptMax), *tHad = new TH1D("tHad", "", nPt, 0., ptMax),
       *tTot = new TH1D("tTot", "", nPt, 0., ptMax);
  for(auto *h : {eEm, eHad, eTot, pEm, pHad, pTot, tEm, tHad, tTot}) h->SetDirectory(nullptr);

  double sumEm = 0., sumHad = 0., sumTot = 0.;
  for(Long64_t i = 0; i < nEvt; i++){
    t->GetEntry(i);
    if(n > maxN){ printf("ERROR: event %lld has %d towers, raise maxN\n", i, n); return; }
    for(int j = 0; j < n; j++){
      eEm->Fill(eta[j], emEt[j]); eHad->Fill(eta[j], hadEt[j]); eTot->Fill(eta[j], et[j]);
      pEm->Fill(phi[j], emEt[j]); pHad->Fill(phi[j], hadEt[j]); pTot->Fill(phi[j], et[j]);
      tEm->Fill(et[j], emEt[j]);  tHad->Fill(et[j], hadEt[j]);  tTot->Fill(et[j], et[j]);
      sumEm += emEt[j]; sumHad += hadEt[j]; sumTot += et[j];
    }
  }
  printf("  %lld events; per event sum E_T: EM %.1f, had %.1f, total %.1f GeV (EM %.1f%%)\n",
         nEvt, sumEm/nEvt, sumHad/nEvt, sumTot/nEvt, 100.*sumEm/sumTot);

  for(auto *h : {eEm, eHad, eTot, pEm, pHad, pTot, tEm, tHad, tTot}){ h->Scale(1./nEvt); h->Scale(1., "width"); }

  const TString tag = Form("%lld events", nEvt);
  drawStack(eEm, eHad, eTot, "tower #it{#eta}", "#LT#Sigma #it{E}_{T}#GT / event / d#it{#eta} [GeV]",
            tag, "towerEtStack_eta");
  drawStack(pEm, pHad, pTot, "tower #it{#phi}", "#LT#Sigma #it{E}_{T}#GT / event / d#it{#phi} [GeV]",
            tag, "towerEtStack_phi");
  drawStack(tEm, tHad, tTot, "tower #it{p}_{T} [GeV]", "#LT#Sigma #it{E}_{T}#GT / event / d#it{p}_{T} [GeV^{-1}]",
            tag, "towerEtStack_pt", true);
  printf("  figures in %s\n", outDir);
}
