// Inclusive-jet JES and JER, calo vs PF, in pp and per PbPb centrality class.
//
// Inputs are the per-class outputs of
// src/jetEnergyResolutionCalculator/jetEnergyResolutionCalculator_pt.C:
//   rootFiles/JES/JES_{caloJets,PFJets}_{pp,C1..C4}.root : JES_result_i
//   rootFiles/JER/JER_{caloJets,PFJets}_{pp,C1..C4}.root : JER_result_i
// "_i" is all jets (inclusive); "_b" is b jets and is not used here.
//
// WHAT mu AND sigma ARE. Both come from the distribution of
// pT^reco / pT^gen for reco jets matched to gen jets, in bins of gen jet pT:
//   mu    = the MEAN of that distribution   (calculator: mu1 = h->GetMean())
//   sigma = its STANDARD DEVIATION          (calculator: sig1 = h->GetStdDev())
// The calculator fits a Gaussian but does not report its parameters, so these
// are moments, not fit values. sigma is absolute, not divided by mu.
//
// Samples: the manual-JEC response scans of 2026-09-21, all four with jet pT
// from the JEC text files on rawpt (never the forest jtpt, which corrects calo
// jets with the AK4PF payload):
//   pp    PYTHIA dijet: ak4Calo with Spring18_ppRef5TeV_V6 L2Relative AK4Calo,
//         ak4PF with Spring18 MC L2Relative AK4PF
//   PbPb  PYTHIA+HYDJET dijet: akPu4Calo with Autumn18_HI_V8 MC L2Relative
//         AK4Calo, akCs4PF with Autumn18_HI_V8 MC L2Relative AK4PF
// C1..C4 are hiBin 0-20, 20-60, 60-100, 100-160 = 0-10, 10-30, 30-50, 50-80%
// (AnalysisSetupV2p3.h), before the scan's hiBinShift of -10.
//
// Figures, in figures/jetEnergyResolution/:
//   JESJER_caloVsPF_<class>.pdf   one per system/class (pp, 0-10 ... 50-80%):
//                                 mu on top, sigma below, calo and PF overlaid
//   JESJER_caloVsPF_allCent.pdf   2x2 overview: mu and sigma for each
//                                 collection, pp and all four classes overlaid
//
// Usage: root -l -b -q 'plotJESJER_caloVsPF.C'
// Run from: src/plots/jetEnergyResolution/

#include "../../../headers/plotting/plotStyle.h"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TSystem.h"
#include "TMath.h"
#include <cstdio>

const char *repo   = "/home/clayton/Analysis/code/bJetRaaAnalysis";
const char *outDir = "/home/clayton/Analysis/code/bJetRaaAnalysis/figures/jetEnergyResolution/";

// index 0 is pp; 1..4 are the PbPb classes
const int   NCls = 5;
const char *clsLabel[NCls] = {"pp", "0-10%", "10-30%", "30-50%", "50-80%"};
const char *clsTag[NCls]   = {"pp", "0to10pct", "10to30pct", "30to50pct", "50to80pct"};
const char *clsFile[NCls]  = {"pp", "C1", "C2", "C3", "C4"};   // file-name suffix
const char *clsSample[NCls]= {"PYTHIA pp", "PYTHIA+HYDJET 0-10%", "PYTHIA+HYDJET 10-30%",
                              "PYTHIA+HYDJET 30-50%", "PYTHIA+HYDJET 50-80%"};
// one color per class in the overview; Okabe-Ito, skipping yellow, with pp in
// purple so it reads as the reference rather than as another PbPb class
const int   clsColor[NCls] = {7, 0, 1, 2, 6};

const char *collTag[2]   = {"caloJets", "PFJets"};
// generic: the algorithm differs between pp (ak4) and PbPb (akPu4 / akCs4)
const char *collLabel[2] = {"calo jets", "PF jets"};
const char *collHex[2]   = {hexCorrected, hexMC};              // vermilion, blue
const int   collMark[2]  = {markFilledSquare, markOpenCircle};

const char *xTitle   = "#it{p}_{T}^{gen} [GeV]";

// x axis starts here; bins below it are neither drawn nor used for the y range
const double xMin = 70.;
const char *muTitle  = "#mu = #LT#it{p}_{T}^{reco}/#it{p}_{T}^{gen}#GT";
const char *sigTitle = "#sigma(#it{p}_{T}^{reco}/#it{p}_{T}^{gen})";

// JES or JER result for one collection and class; nullptr if missing
static TH1D* getResult(const char *what, int coll, int cls)
{
  TString path = Form("%s/rootFiles/%s/%s_%s_%s.root", repo, what, what, collTag[coll], clsFile[cls]);
  TFile *f = TFile::Open(path);
  if(!f || f->IsZombie()){ printf("ERROR: cannot open %s\n", path.Data()); return nullptr; }
  TH1D *h = nullptr;
  f->GetObject(Form("%s_result_i", what), h);
  if(!h){ printf("ERROR: %s_result_i missing in %s\n", what, path.Data()); f->Close(); return nullptr; }
  h = (TH1D*) h->Clone(Form("%s_%s_%s", what, collTag[coll], clsFile[cls]));
  h->SetDirectory(nullptr);
  f->Close();
  return h;
}

// y range covering every drawn histogram from xMin up, with fractional padding
static void yRange(TH1D **hs, int n, double padLo, double padHi, double &lo, double &hi)
{
  lo = 1e9; hi = -1e9;
  for(int i = 0; i < n; i++){
    if(!hs[i]) continue;
    for(int b = 1; b <= hs[i]->GetNbinsX(); b++){
      if(hs[i]->GetXaxis()->GetBinLowEdge(b) < xMin - 1e-6) continue;
      double v = hs[i]->GetBinContent(b), e = hs[i]->GetBinError(b);
      if(v == 0.) continue;   // the calculator writes 0 for skipped bins
      lo = TMath::Min(lo, v - e); hi = TMath::Max(hi, v + e);
    }
  }
  double span = hi - lo;
  lo -= padLo * span; hi += padHi * span;
}

static void styleFrame(TH1D *h, const char *yT, double lo, double hi,
                       bool showX, double scale)
{
  h->SetTitle(""); h->SetStats(0);
  h->SetMinimum(lo); h->SetMaximum(hi);
  h->GetXaxis()->SetRangeUser(xMin, h->GetXaxis()->GetXmax());
  h->GetXaxis()->SetTitle(showX ? xTitle : "");
  h->GetXaxis()->SetTitleSize(0.050 * scale); h->GetXaxis()->SetLabelSize(showX ? 0.045 * scale : 0.);
  h->GetXaxis()->SetTitleOffset(1.05);
  h->GetYaxis()->SetTitle(yT);
  h->GetYaxis()->SetTitleSize(0.050 * scale); h->GetYaxis()->SetLabelSize(0.045 * scale);
  h->GetYaxis()->SetTitleOffset(1.45 / scale);
  h->GetYaxis()->SetNdivisions(505);
}

void plotJESJER_caloVsPF()
{
  initPlotStyle();
  gSystem->mkdir(outDir, kTRUE);

  TH1D *mu[2][NCls], *sig[2][NCls];
  for(int k = 0; k < 2; k++)
    for(int c = 0; c < NCls; c++){
      mu[k][c]  = getResult("JES", k, c);
      sig[k][c] = getResult("JER", k, c);
      if(!mu[k][c] || !sig[k][c]) return;
    }

  // ---- table ----------------------------------------------------------------
  const double probe[4] = {75., 105., 155., 280.};
  printf("\n  inclusive jets, PYTHIA (pp) and PYTHIA+HYDJET: mu = mean, sigma = std. dev. of pT^reco/pT^gen\n");
  printf("  %-7s %6s | %8s %8s | %8s %8s %11s\n",
         "class", "pTgen", "mu calo", "mu PF", "sig calo", "sig PF", "sig calo/PF");
  for(int c = 0; c < NCls; c++)
    for(double p : probe){
      int b = mu[0][c]->FindBin(p);
      double mc = mu[0][c]->GetBinContent(b), mp = mu[1][c]->GetBinContent(b);
      double sc = sig[0][c]->GetBinContent(b), sp = sig[1][c]->GetBinContent(b);
      printf("  %-7s %6.0f | %8.4f %8.4f | %8.4f %8.4f %11.3f\n",
             p == probe[0] ? clsLabel[c] : "", p, mc, mp, sc, sp, sp > 0. ? sc/sp : -1.);
    }

  // ---- common y ranges ------------------------------------------------------
  // One mu range and one sigma range for EVERY figure and panel, from all ten
  // histograms (both collections, pp and the four classes), so any two plots
  // can be compared by eye. The top padding leaves room for the legends.
  double muLo, muHi, sLo, sHi;
  {
    TH1D *allMu[2*NCls], *allSig[2*NCls];
    for(int k = 0; k < 2; k++) for(int c = 0; c < NCls; c++){ allMu[k*NCls+c] = mu[k][c]; allSig[k*NCls+c] = sig[k][c]; }
    yRange(allMu,  2*NCls, 0.10, 0.60, muLo, muHi);
    muLo = TMath::Min(muLo, 0.985); muHi = TMath::Max(muHi, 1.015);   // keep unity on the axis
    yRange(allSig, 2*NCls, 0.10, 0.30, sLo, sHi);
    sLo = TMath::Max(sLo, 0.);
    printf("\n  common y ranges from %.0f GeV: mu [%.3f, %.3f]   sigma [%.3f, %.3f]\n", xMin, muLo, muHi, sLo, sHi);
  }

  // ---- per class: mu on top, sigma below ------------------------------------
  for(int c = 0; c < NCls; c++){
    TCanvas *cv = new TCanvas(Form("c_%s", clsTag[c]), "", 700, 800);
    TPad *pTop = new TPad(Form("pTop_%d", c), "", 0, 0.50, 1, 1);
    TPad *pBot = new TPad(Form("pBot_%d", c), "", 0, 0,    1, 0.50);
    pTop->SetLeftMargin(0.17); pTop->SetRightMargin(0.05); pTop->SetTopMargin(0.17); pTop->SetBottomMargin(0.02);
    pBot->SetLeftMargin(0.17); pBot->SetRightMargin(0.05); pBot->SetTopMargin(0.02); pBot->SetBottomMargin(0.19);
    pTop->Draw(); pBot->Draw();

    TH1D *hm[2] = {mu[0][c], mu[1][c]}, *hs[2] = {sig[0][c], sig[1][c]};

    pTop->cd();
    for(int k = 0; k < 2; k++){
      TH1D *h = (TH1D*) hm[k]->Clone(Form("m_%d_%d", k, c));
      styleH(h, collHex[k], collMark[k]);
      if(k == 0){ styleFrame(h, muTitle, muLo, muHi, false, 1.25); h->Draw("E1"); }
      else h->Draw("E1 same");
    }
    TLine one; one.SetLineStyle(7); one.SetLineColor(kGray + 2);
    one.DrawLine(xMin, 1., hm[0]->GetXaxis()->GetXmax(), 1.);
    for(int k = 0; k < 2; k++) hm[k]->Draw("E1 same");   // markers over the line

    TLegend *leg = makeLegend(0.52, 0.62, 0.94, 0.80, 0.055);
    for(int k = 0; k < 2; k++){ styleH(hm[k], collHex[k], collMark[k]); leg->AddEntry(hm[k], collLabel[k], "lp"); }
    leg->Draw();

    TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.060);
    la.DrawLatex(0.17, 0.90, Form("%s, inclusive jets", clsSample[c]));

    pBot->cd();
    for(int k = 0; k < 2; k++){
      TH1D *h = (TH1D*) hs[k]->Clone(Form("s_%d_%d", k, c));
      styleH(h, collHex[k], collMark[k]);
      if(k == 0){ styleFrame(h, sigTitle, sLo, sHi, true, 1.25); h->Draw("E1"); }
      else h->Draw("E1 same");
    }

    cv->SaveAs(Form("%sJESJER_caloVsPF_%s.pdf", outDir, clsTag[c]));
    delete cv;
  }

  // ---- overview: 2x2, all classes -------------------------------------------
  {
    TCanvas *cv = new TCanvas("c_all", "", 700, 800);
    cv->Divide(2, 2, 0.001, 0.001);
    // same common ranges as the per-class figures

    for(int row = 0; row < 2; row++)
      for(int k = 0; k < 2; k++){
        cv->cd(1 + row*2 + k);
        gPad->SetLeftMargin(0.20); gPad->SetRightMargin(0.04);
        gPad->SetTopMargin(0.10);  gPad->SetBottomMargin(0.15);
        for(int c = 0; c < NCls; c++){
          TH1D *src = row ? sig[k][c] : mu[k][c];
          TH1D *h = (TH1D*) src->Clone(Form("o_%d_%d_%d", row, k, c));
          styleH(h, okabeHex[clsColor[c]], collMark[k], 0.7);
          if(c == 0){
            styleFrame(h, row ? sigTitle : muTitle, row ? sLo : muLo, row ? sHi : muHi, true, 1.0);
            h->GetYaxis()->SetTitleOffset(1.75);
            h->Draw("E1");
            if(!row){
              TLine one; one.SetLineStyle(7); one.SetLineColor(kGray + 2);
              one.DrawLine(xMin, 1., h->GetXaxis()->GetXmax(), 1.);
            }
          }
          else h->Draw("E1 same");
        }
        TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.055);
        la.DrawLatex(0.20, 0.925, collLabel[k]);
      }

    // one legend for the classes, in the top-right (PF mu) pad
    cv->cd(2);
    TLegend *leg = makeLegend(0.55, 0.50, 0.95, 0.88, 0.050);
    for(int c = 0; c < NCls; c++){
      TH1D *h = (TH1D*) mu[1][c]->Clone(Form("legc_%d", c));
      styleH(h, okabeHex[clsColor[c]], markFilledCircle, 0.8);
      leg->AddEntry(h, clsLabel[c], "lp");
    }
    leg->Draw();

    cv->SaveAs(Form("%sJESJER_caloVsPF_allCent.pdf", outDir));
    delete cv;
  }

  printf("\n  figures in %s\n", outDir);
}
