// Tower eta distributions in slices of tower pT (= E_T for a massless tower),
// from rechitanalyzerpp/tower.
//
// Two figures:
//   towerEta_ptSlices_shape.pdf  unit area per slice, linear axis -- compares the
//                                SHAPE of each slice
//   towerEta_ptSlices_yield.pdf  towers per event per unit eta, log y -- the slices
//                                span ~2 decades in yield, so on a linear axis the
//                                highest slices would be flat lines on zero. Log is
//                                safe here: every bin is positive.
//
// 0.5-wide eta bins as in plotTowerKinematics.C, which smooth out the tower-grid
// aliasing. The highest slice has only ~1250 towers in 100 events, so its shape
// carries visible statistical error bars.
//
// Usage: root -l -b -q 'plotTowerEta_ptSlices.C("/path/to/HiForestAOD.root")'
// Run from: src/plots/towers/

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TSystem.h"
#include "../../../headers/plotting/plotStyle.h"

namespace {

const char *outDir = "../../../figures/towers/";
const int nS = 5;
const double sLo[nS] = {0., 0.5, 1., 2., 5.};
const double sHi[nS] = {0.5, 1., 2., 5., 1e9};
const char  *sLab[nS] = {"0 - 0.5 GeV", "0.5 - 1 GeV", "1 - 2 GeV", "2 - 5 GeV", "> 5 GeV"};
// ordered so adjacent slices stay distinguishable; yellow is left out
const char  *sHex[nS] = {"#000000", "#0072B2", "#009E73", "#E69F00", "#D55E00"};
const int    sMark[nS] = {markFilledCircle, markFilledSquare, markFilledDiamond, markOpenCircle, markOpenSquare};

void draw(TH1D **h, bool logY, const char *yTitle, const TString &tag, const char *outName)
{
  // Headroom for a legend drawn INSIDE the frame, upper right: the highest
  // slice's curve dips there (see the two figures' shapes), so that corner is
  // free once the y-axis is stretched -- 2.2 decades log, 2x linear -- above the
  // tallest point actually plotted.
  TCanvas *c = new TCanvas(Form("c_%s", outName), "", 900, 700);
  c->SetLeftMargin(0.15); c->SetBottomMargin(0.13); c->SetRightMargin(0.05); c->SetTopMargin(0.06);
  c->SetLogy(logY);

  double mx = 0., mn = 1e30;
  for(int s = 0; s < nS; s++){
    mx = TMath::Max(mx, h[s]->GetMaximum());
    mn = TMath::Min(mn, h[s]->GetMinimum(0.));
  }
  for(int s = 0; s < nS; s++){
    styleH(h[s], sHex[s], sMark[s], 1.0);
    h[s]->GetXaxis()->SetTitle("tower #it{#eta}");
    h[s]->GetYaxis()->SetTitle(yTitle);
    h[s]->GetXaxis()->SetTitleSize(0.05); h[s]->GetYaxis()->SetTitleSize(0.05);
    h[s]->GetXaxis()->SetLabelSize(0.042); h[s]->GetYaxis()->SetLabelSize(0.042);
    h[s]->GetYaxis()->SetTitleOffset(1.4);
    h[s]->SetMaximum(mx * (logY ? 160. : 2.0));
    h[s]->SetMinimum(logY ? mn * 0.5 : 0.);
    h[s]->Draw(s == 0 ? "E" : "E SAME");
  }
  TLegend *leg = makeLegend(0.50, 0.68, 0.93, 0.86, 0.032);
  leg->SetHeader("tower #it{p}_{T}");
  leg->SetNColumns(2);
  for(int s = 0; s < nS; s++) leg->AddEntry(h[s], sLab[s], "pl");
  leg->Draw();

  TLatex l; l.SetNDC(); l.SetTextSize(0.038);
  l.DrawLatex(0.19, 0.90, "Calorimeter towers (rechitanalyzerpp)");
  l.DrawLatex(0.19, 0.85, tag);
  savePdfTight(c, Form("%s%s.pdf", outDir, outName));
  delete c;
}

} // namespace

void plotTowerEta_ptSlices(const char *path = "/home/clayton/Downloads/HiForestAOD.root")
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
  static Float_t et[maxN], eta[maxN];
  t->SetBranchAddress("n", &n);
  t->SetBranchAddress("et", et);
  t->SetBranchAddress("eta", eta);

  TH1D *hY[nS], *hS[nS];
  for(int s = 0; s < nS; s++){
    hY[s] = new TH1D(Form("hY_%d", s), "", 20, -5., 5.);
    hY[s]->Sumw2(); hY[s]->SetDirectory(nullptr);
  }
  for(Long64_t i = 0; i < nEvt; i++){
    t->GetEntry(i);
    if(n > maxN){ printf("ERROR: event %lld has %d towers, raise maxN\n", i, n); return; }
    for(int j = 0; j < n; j++)
      for(int s = 0; s < nS; s++)
        if(et[j] >= sLo[s] && et[j] < sHi[s]){ hY[s]->Fill(eta[j]); break; }
  }

  printf("  %-14s %10s %12s\n", "slice", "towers", "per event");
  for(int s = 0; s < nS; s++){
    const double N = hY[s]->Integral();
    printf("  %-14s %10.0f %12.1f\n", sLab[s], N, N/nEvt);
    hS[s] = (TH1D*) hY[s]->Clone(Form("hS_%d", s));
    hS[s]->SetDirectory(nullptr);
    hS[s]->Scale(1./N);                      // unit area, before the width division
    hS[s]->Scale(1., "width");
    hY[s]->Scale(1./nEvt); hY[s]->Scale(1., "width");
  }

  const TString tag = Form("%lld events", nEvt);
  draw(hS, false, "1/N dN/d#it{#eta}", tag, "towerEta_ptSlices_shape");
  draw(hY, true,  "1/N_{evt} dN/d#it{#eta}", tag, "towerEta_ptSlices_yield");
  printf("  figures in %s\n", outDir);
}
