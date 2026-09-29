// plotTowerPtFractionStack.C split into |eta| regions: barrel, endcap, forward.
//
// Same stacked EM / hadronic share of E_T versus tower pT, one figure per region.
// The inclusive figure is a mixture of these, weighted by how much E_T each region
// carries in each pT bin, which is what produces its shape. Region edges are
// 1.3 and 3.0 (roughly the barrel/endcap and endcap/forward boundaries).
//
// Bins are wider than in the inclusive figure and variable, since the barrel has
// very few towers above 5 GeV in 100 events (259 at 5-10, 22 above 10). The
// fraction is a ratio of sums so bin width does not enter. Each figure prints its
// tower count per bin; treat bins with few towers as noise.
//
// Usage: root -l -b -q 'plotTowerPtFractionStack_etaRegions.C("/path/to/HiForestAOD.root")'
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


void plotTowerPtFractionStack_etaRegions(const char *path = "/home/clayton/Downloads/HiForestAOD.root")
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
  static Float_t et[maxN], emEt[maxN], hadEt[maxN], eta[maxN];
  t->SetBranchAddress("n", &n);
  t->SetBranchAddress("et", et);
  t->SetBranchAddress("emEt", emEt);
  t->SetBranchAddress("hadEt", hadEt);
  t->SetBranchAddress("eta", eta);

  const double edge[] = {0, 0.5, 1, 1.5, 2, 3, 4, 5, 7, 10, 15, 20};
  const int nB = sizeof(edge)/sizeof(double) - 1;

  const int nReg = 3;
  const char *regTag[nReg] = {"barrel", "endcap", "forward"};
  const char *regLab[nReg] = {"|#font[52]{#eta}| < 1.3", "1.3 < |#font[52]{#eta}| < 3", "|#font[52]{#eta}| > 3"};
  auto regionOf = [](double e){ double a = TMath::Abs(e); return a < 1.3 ? 0 : (a < 3. ? 1 : 2); };

  TH1D *hEm[nReg], *hHad[nReg], *hTot[nReg], *hN[nReg];
  for(int r = 0; r < nReg; r++){
    hEm[r]  = new TH1D(Form("hEm_%d",  r), "", nB, edge);
    hHad[r] = new TH1D(Form("hHad_%d", r), "", nB, edge);
    hTot[r] = new TH1D(Form("hTot_%d", r), "", nB, edge);
    hN[r]   = new TH1D(Form("hN_%d",   r), "", nB, edge);
  }
  for(Long64_t i = 0; i < nEvt; i++){
    t->GetEntry(i);
    if(n > maxN){ printf("ERROR: event %lld has %d towers, raise maxN\n", i, n); return; }
    for(int j = 0; j < n; j++){
      const int r = regionOf(eta[j]);
      hEm[r]->Fill(et[j], emEt[j]); hHad[r]->Fill(et[j], hadEt[j]);
      hTot[r]->Fill(et[j], et[j]);  hN[r]->Fill(et[j]);
    }
  }

  const int cEm = TColor::GetColor("#E69F00"), cHad = TColor::GetColor("#0072B2");
  for(int r = 0; r < nReg; r++){
    TH1D *sEm  = (TH1D*) hEm[r]->Clone(Form("sEm_%d", r));
    TH1D *sHad = (TH1D*) hHad[r]->Clone(Form("sHad_%d", r));
    sEm->Divide(hTot[r]); sHad->Divide(hTot[r]);
    sHad->SetFillColor(cHad); sHad->SetLineColor(cHad); sHad->SetLineWidth(1);
    sEm ->SetFillColor(cEm);  sEm ->SetLineColor(cEm);  sEm ->SetLineWidth(1);

    THStack *hs = new THStack(Form("hs_%d", r), "");
    hs->Add(sHad);
    hs->Add(sEm);

    TLegend *ls = new TLegend(0.71, 0.6, 0.98, 0.9);
    ls->SetBorderSize(0);
    ls->SetTextSize(0.055);
    ls->AddEntry(sEm,  "EM #font[52]{p}_{T}", "f");
    ls->AddEntry(sHad, "hadronic #font[52]{p}_{T}", "f");

    TCanvas *cs = new TCanvas(Form("cs_%d", r), "cs", 1000, 700);
    cs->cd();
    TPad *ps = new TPad(Form("ps_%d", r), "ps", 0, 0, 1, 1);
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
    la->DrawLatexNDC(0.72, 0.47, regLab[r]);
    la->DrawLatexNDC(0.72, 0.42, Form("%lld events", nEvt));
    la->DrawLatexNDC(0.72, 0.37, "#font[52]{E}_{T} weighted");

    const TString out = Form("%stowerPtFractionStack_%s.pdf", outDir, regTag[r]);
    savePdfTight(cs, out);

    printf("\n  %s: overall EM share %.3f\n  %-10s %8s %8s\n", regTag[r],
           hEm[r]->Integral()/hTot[r]->Integral(), "pT [GeV]", "towers", "EM frac");
    for(int b = 1; b <= nB; b++)
      printf("  %4.1f-%-5.1f %8.0f %8.3f\n", edge[b-1], edge[b], hN[r]->GetBinContent(b), sEm->GetBinContent(b));
  }
}
