// Calo vs PF jets per Z, PbPb/pp, four nominal centrality classes, with the
// CORRECT jet energy calibration on each side.
//
//   calo : jetsPerZ_caloJets_10pctBins.C(true) -- manual AK4Calo JEC on every
//          spectrum, built from the ultra-fine scans and merged to the nominal
//          classes (0-10 = slices 1-2, 10-30 = 3-6, 30-50 = 7-10,
//          50-80 = 11-16), trigger-only stitch.
//   PF   : calculateJetsPerZ_caloJets.C(false,...,false) -- trigger-only,
//          forest JEC.
//
// WHY THE TWO USE DIFFERENT JEC SETTINGS, AND WHY THAT IS CORRECT. The forest
// bug is that it applies the AK4PF payload to CALO jets. For PF jets AK4PF is
// the right payload, so the forest calibration is already correct there and a
// manual-JEC PF rescan would reproduce it. Only the calo side needed fixing.
//
// This replaces the earlier calo-vs-PF comparison, which took its calo curve
// from calculateJetsPerZ_caloJets.C's forest-JEC input list. That curve was
// mis-calibrated and sat ~25% high in peripheral classes: its 50-80% read 0.784
// at 200-300 GeV where the Z-weighted average of the 10% classes gives 0.63.
//
// BINNING. The calo side is now produced on the PF chain's fixed axis
// (useFixedAxis in jetsPerZ_caloJets_10pctBins.C), so above the calo Jet80
// threshold both share {120, 150, 200, 300, 500} and a genuine bin-by-bin
// calo/PF ratio is formed in the lower pad. PF bins below 120 have no calo
// partner and are drawn in the upper pad only.
//
// Usage, from src/newFractionCalculation/:
//   root -l -b -q plotJetsPerZ_caloVsPF_manualJEC_coarse.C

#include "../../headers/plotting/plotStyle.h"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include <cstdio>

namespace {

const char *f_calo = "./rootFiles/JetsPerZ/jetsPerZ_caloJets_coarseBins.root";
const char *f_PF   = "./rootFiles/JetsPerZ/histograms_JetsPerZ_PFJets_trigOnly_noUnfold.root";
const char *outPath = "/home/clayton/Analysis/code/bJetRaaAnalysis/figures/"
                      "JetsPerZ/JetsPerZ_caloVsPF_manualJEC_coarse.pdf";

const int   NCls = 4;
const char *caloTag[NCls] = {"0to10","10to30","30to50","50to80"};
const char *label  [NCls] = {"PbPb 0-10%","PbPb 10-30%","PbPb 30-50%","PbPb 50-80%"};

// Nothing below 150 GeV is shown, on either curve.
//
// The PF chain takes pp from MinBias below 150 (its handover is hard-coded
// there) and that normalisation undershoots by a factor of ~1.9: converted back
// to jet counts, the PF chain gives 9.71e6 in 120-150 where the raw pp Jet80
// sample holds 1.87e7. Above 150, where PF switches to the triggered sample,
// both chains reproduce their raw counts to better than 1% (calo 7.969e6 vs
// 7.969e6 raw; PF 7.474e6 vs 7.463e6) and calo/PF = 1.066 against the raw
// 1.068. So 120-150 was not a measurement and the earlier ~0.82 point there was
// an artifact of that normalisation, not a calo/PF difference.
const double quoteFromPt = 150.;

const double ptLo = 150., ptHi = 500.;

TGraphErrors *toGraph(TH1D *h)
{
  TGraphErrors *g = new TGraphErrors();
  if(!h) return g;
  int n = 0;
  for(int i = 1; i <= h->GetNbinsX(); i++){
    if(h->GetBinContent(i) <= 0.) continue;       // unfilled bins are not points
    const double lo = h->GetBinLowEdge(i), hi = h->GetBinLowEdge(i+1);
    if(lo < quoteFromPt - 1e-6) continue;         // see quoteFromPt above
    g->SetPoint(n, 0.5*(lo+hi), h->GetBinContent(i));
    g->SetPointError(n, 0.5*(hi-lo), h->GetBinError(i));
    n++;
  }
  return g;
}

} // namespace

void plotJetsPerZ_caloVsPF_manualJEC_coarse()
{
  initPlotStyle();

  TFile *fc = TFile::Open(f_calo);
  TFile *fp = TFile::Open(f_PF);
  if(!fc || fc->IsZombie()){
    printf("ERROR: run jetsPerZ_caloJets_10pctBins.C(true) first\n"); return; }
  if(!fp || fp->IsZombie()){
    printf("ERROR: run calculateJetsPerZ_caloJets.C(false,false,false,false) first\n"); return; }

  printf("\n  calo (manual JEC) vs PF (forest JEC), trigger-only\n");

  TCanvas *c = new TCanvas("cCmp", "", 1000, 950);

  for(int k = 0; k < NCls; k++){
    const int col = k % 2, row = k / 2;
    const double x0 = col/2., x1 = (col+1)/2.;
    const double y1 = 1. - row/2., y0 = 1. - (row+1)/2.;

    const double split = 0.33;
    const double ys = y0 + split*(y1-y0);

    c->cd();
    TPad *pB = new TPad(Form("pB%d", k), "", x0, y0, x1, ys);
    pB->SetLeftMargin(0.16); pB->SetRightMargin(0.03);
    pB->SetTopMargin(0.02);  pB->SetBottomMargin(0.34);
    pB->Draw();

    TPad *p = new TPad(Form("p%d", k), "", x0, ys, x1, y1);
    p->SetLeftMargin(0.16); p->SetRightMargin(0.03);
    p->SetTopMargin(0.07);  p->SetBottomMargin(0.02);
    p->Draw(); p->cd();

    TH1D *fr = new TH1D(Form("fr%d", k), "", 1, ptLo, ptHi);
    fr->SetStats(0);
    fr->GetYaxis()->SetRangeUser(0., 1.2);
    fr->GetXaxis()->SetLabelSize(0.);
    fr->GetYaxis()->SetTitle("(jets per Z)_{PbPb} / (jets per Z)_{pp}");
    fr->GetYaxis()->SetTitleSize(0.062);
    fr->GetYaxis()->SetLabelSize(0.058);
    fr->GetYaxis()->SetTitleOffset(1.20);
    fr->GetYaxis()->SetNdivisions(505);
    fr->Draw("axis");

    TLine *one = new TLine(ptLo, 1., ptHi, 1.);
    one->SetLineStyle(2); one->SetLineColor(kGray+2); one->Draw();

    TH1D *hc = (TH1D*) fc->Get(Form("ratio_%s", caloTag[k]));
    TH1D *hp = (TH1D*) fp->Get(Form("r_C%d_r", k+1));

    TGraphErrors *gc = toGraph(hc);
    TGraphErrors *gp = toGraph(hp);

    const int cCalo = TColor::GetColor(okabeHex[0]);   // black
    const int cPF   = TColor::GetColor(okabeHex[6]);   // vermillion
    gc->SetLineColor(cCalo); gc->SetMarkerColor(cCalo);
    gc->SetMarkerStyle(markFilledCircle); gc->SetMarkerSize(1.3); gc->SetLineWidth(2);
    gp->SetLineColor(cPF);   gp->SetMarkerColor(cPF);
    gp->SetMarkerStyle(markOpenSquare);   gp->SetMarkerSize(1.3); gp->SetLineWidth(2);

    gp->Draw("pz same");
    gc->Draw("pz same");

    TLatex t; t.SetNDC(); t.SetTextFont(42);
    t.SetTextSize(0.055); t.DrawLatex(0.21, 0.87, label[k]);

    if(k == 0){
      TLegend *leg = new TLegend(0.40, 0.16, 0.95, 0.33);
      leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.042);
      leg->AddEntry(gc, "calo jets, manual AK4Calo JEC", "lp");
      leg->AddEntry(gp, "PF jets, forest JEC", "lp");
      leg->Draw();
    }

    // ---- ratio pad: only bins the two share -------------------------------
    pB->cd();
    TH1D *fb = new TH1D(Form("fb%d", k), "", 1, ptLo, ptHi);
    fb->SetStats(0);
    fb->GetYaxis()->SetRangeUser(0.75, 1.08);
    fb->GetXaxis()->SetTitle("jet p_{T} [GeV]");
    fb->GetYaxis()->SetTitle("calo / PF");
    fb->GetXaxis()->SetTitleSize(0.125); fb->GetXaxis()->SetLabelSize(0.115);
    fb->GetYaxis()->SetTitleSize(0.110); fb->GetYaxis()->SetLabelSize(0.105);
    fb->GetXaxis()->SetTitleOffset(1.15);
    fb->GetYaxis()->SetTitleOffset(0.62);
    fb->GetYaxis()->SetNdivisions(505);
    fb->Draw("axis");
    TLine *ob = new TLine(ptLo, 1., ptHi, 1.);
    ob->SetLineStyle(2); ob->SetLineColor(kGray+2); ob->Draw();

    printf("\n  %s\n     pT bin        calo        PF     calo/PF\n", label[k]);
    TGraphErrors *gr = new TGraphErrors();
    int nr = 0;
    if(hc && hp){
      for(int i = 1; i <= hc->GetNbinsX(); i++){
        const double lo = hc->GetBinLowEdge(i), hi = hc->GetBinLowEdge(i+1);
        const double a = hc->GetBinContent(i), ea = hc->GetBinError(i);
        if(a <= 0.) continue;
        if(lo < quoteFromPt - 1e-6) continue;
        // pair by EDGE, not index: the two axes share edges above 120 but the
        // PF axis starts lower, so the indices are offset
        const int j = hp->FindBin(0.5*(lo+hi));
        if(j < 1 || j > hp->GetNbinsX()) continue;
        if(fabs(hp->GetBinLowEdge(j) - lo) > 1e-6 ||
           fabs(hp->GetBinLowEdge(j+1) - hi) > 1e-6) continue;
        const double b = hp->GetBinContent(j), eb = hp->GetBinError(j);
        if(b <= 0.) continue;
        const double r = a/b;
        const double er = r*sqrt(pow(ea/a,2) + pow(eb/b,2));
        gr->SetPoint(nr, 0.5*(lo+hi), r);
        gr->SetPointError(nr, 0.5*(hi-lo), er);
        nr++;
        printf("   %4.0f-%4.0f  %10.4f %10.4f %9.3f +- %.3f\n", lo, hi, a, b, r, er);
      }
    }
    gr->SetLineColor(TColor::GetColor(okabeHex[0]));
    gr->SetMarkerColor(TColor::GetColor(okabeHex[0]));
    gr->SetMarkerStyle(markFilledCircle); gr->SetMarkerSize(1.2); gr->SetLineWidth(2);
    gr->Draw("pz same");
  }

  c->SaveAs(outPath);
  printf("\n  figure: %s\n", outPath);
}
