// EM and hadronic share of tower E_T (= pT for a massless tower) versus tower pT,
// stacked to 1, from rechitanalyzerpp/tower.
//
// Layout follows bJetMuonTaggingAnalysis/src/plots/flavorFraction/projectFlavor.C:
// a stacked "Fraction" plot, 1000x700 canvas, right margin left open for the
// legend and labels, italic (font 52) variable names. Colors and PDF output
// follow this repo's conventions instead of that file's kRed/kBlue and PNG.
//
// The fraction in a tower-pT bin is (sum of emEt) / (sum of et) over the towers
// in the bin, and likewise for hadEt. Since et = emEt + hadEt tower by tower the
// two fractions add to exactly 1 in every bin. The bin is set by the tower's total
// pT, so the low-pT bins are towers with little total energy and the split there
// is set by which of the two sub-detectors fired.
//
// Usage: root -l -b -q 'plotTowerPtFractionStack.C("/path/to/HiForestAOD.root")'
// Run from: src/plots/towers/

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLatex.h"
#include "TSystem.h"
#include "../../../headers/plotting/plotStyle.h"


void plotTowerPtFractionStack(const char *path = "/home/clayton/Downloads/HiForestAOD.root")
{
  initPlotStyle();
  const char *outDir = "../../../figures/towers/";
  gSystem->mkdir(outDir, kTRUE);

  TFile *f = TFile::Open(path);
  if(!f || f->IsZombie()){ printf("ERROR: %s did not open\n", path); return; }
  TTree *t = nullptr; f->GetObject("rechitanalyzerpp/tower", t);
  if(!t){ printf("ERROR: rechitanalyzerpp/tower missing\n"); return; }

  const Long64_t nEvt = t->GetEntries();
  const int maxN = 20000;
  Int_t n = 0;
  static Float_t et[maxN], emEt[maxN], hadEt[maxN];
  t->SetBranchAddress("n", &n);
  t->SetBranchAddress("et", et);
  t->SetBranchAddress("emEt", emEt);
  t->SetBranchAddress("hadEt", hadEt);

  // 0.5 GeV bins, 0-20 GeV; above 20 GeV there are only a handful of towers
  const int nBins = 40; const double xMax = 20.;
  TH1D *hEm  = new TH1D("hEmSum",  "", nBins, 0., xMax);
  TH1D *hHad = new TH1D("hHadSum", "", nBins, 0., xMax);
  TH1D *hTot = new TH1D("hTotSum", "", nBins, 0., xMax);
  for(Long64_t i = 0; i < nEvt; i++){
    t->GetEntry(i);
    if(n > maxN){ printf("ERROR: event %lld has %d towers, raise maxN\n", i, n); return; }
    for(int j = 0; j < n; j++){
      hEm->Fill(et[j], emEt[j]); hHad->Fill(et[j], hadEt[j]); hTot->Fill(et[j], et[j]);
    }
  }

  TH1D *hsEm  = (TH1D*) hEm->Clone("hsEm");
  TH1D *hsHad = (TH1D*) hHad->Clone("hsHad");
  hsEm->Divide(hTot); hsHad->Divide(hTot);
  for(auto *h : {hsEm, hsHad}) h->SetDirectory(nullptr);

  const int cEm = TColor::GetColor("#E69F00"), cHad = TColor::GetColor("#0072B2");
  hsHad->SetFillColor(cHad); hsHad->SetLineColor(cHad); hsHad->SetLineWidth(1);
  hsEm ->SetFillColor(cEm);  hsEm ->SetLineColor(cEm);  hsEm ->SetLineWidth(1);

  THStack *hs = new THStack("hs", "");
  hs->Add(hsHad);
  hs->Add(hsEm);

  TLegend *ls = new TLegend(0.71, 0.6, 0.98, 0.9);
  ls->SetBorderSize(0);
  ls->SetTextSize(0.055);
  ls->AddEntry(hsEm,  "EM #font[52]{p}_{T}", "f");
  ls->AddEntry(hsHad, "hadronic #font[52]{p}_{T}", "f");

  TCanvas *cs = new TCanvas("cs", "cs", 1000, 700);
  cs->cd();
  TPad *ps = new TPad("ps", "ps", 0, 0, 1, 1);
  ps->SetLeftMargin(0.2);
  ps->SetBottomMargin(0.2);
  ps->SetRightMargin(0.3);
  ps->Draw();
  ps->cd();

  hs->Draw("hist");
  hs->GetYaxis()->SetTitleSize(0.065);
  hs->GetXaxis()->SetTitleSize(0.065);
  hs->GetYaxis()->SetLabelSize(0.04);
  hs->GetXaxis()->SetLabelSize(0.04);
  hs->GetYaxis()->SetTitle("Fraction");
  hs->GetXaxis()->SetTitle("#font[52]{p}_{T}^{tower} [GeV]");
  hs->SetMinimum(0.);
  hs->SetMaximum(1.0);

  ls->Draw();

  TLatex *la = new TLatex();
  la->SetTextFont(42);
  la->SetTextSize(0.036);
  la->DrawLatexNDC(0.72, 0.52, "Calorimeter towers");
  la->DrawLatexNDC(0.72, 0.47, Form("%lld events", nEvt));

  const TString out = Form("%stowerPtFractionStack.pdf", outDir);
  savePdfTight(cs, out);

  printf("  overall EM share of sum E_T: %.3f\n", hEm->Integral()/hTot->Integral());
  for(int b : {1, 2, 4, 10, 20, 40})
    printf("  pT %4.1f-%4.1f GeV: EM %.3f  had %.3f\n", hTot->GetBinLowEdge(b),
           hTot->GetBinLowEdge(b) + hTot->GetBinWidth(b), hsEm->GetBinContent(b), hsHad->GetBinContent(b));
}
