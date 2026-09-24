// Inclusive b-jet fraction: (corrected b-jet spectrum) / (corrected inclusive
// jet spectrum), formed in each system SEPARATELY -- pp and each PbPb class.
// This is not a PbPb/pp ratio; nothing cancels between systems here.
//
// NUMERATOR   bJetSpectrum_* from src/calculateBJetsPerZ/calculateBJetsPerZ.cc:
//             muon-tagged jets x b purity (with the fake-ptRel T3 subtraction),
//             Bayes-unfolded on the variable-width b-jet response with 60-80
//             acting as buffer bins, then divided by corrFactor_1.
// DENOMINATOR inclSpectrum_* from src/newFractionCalculation/calculateRAA.C:
//             stitched, fake-jet subtracted, Bayes-unfolded on the buffered
//             coarse response, 1 iteration.
//
// Both are per Z and per GeV. Since the fraction is taken within one system,
// N_Z cancels exactly and the result is a pure yield fraction -- it does not
// depend on the per-Z normalisation being right, only on it being the same in
// numerator and denominator, which it is (same dimuon window, same efficiency).
//
// ETA CONVENTION, and why it had to be handled. calculateRAA.C divides its
// spectra by 3.2 (the |eta| < 1.6 range) and the b-jet chain does not. Both
// cover the same acceptance, so it is a units choice, but a fraction formed
// across the two conventions would come out 3.2x too small. calculateRAA.C now
// writes inclSpectrum_* with that factor undone, so both sides are per Z per
// GeV with no eta division.
//
// BINNING. The b-jet axis is {60,70,80,90,100,120,150,200,300,500}; the
// inclusive axis is that with 20 and 40 prepended (the unfolding buffer). Bins
// are paired by EDGE, not index -- the indices are offset by two.
//
// QUOTED FROM 80 GeV. The b-jet 60-70 and 70-80 bins are unfolding buffers, not
// results.
//
// Usage, from src/calculateBJetsPerZ/:
//   root -l -b -q plotInclusiveBJetFraction.C

#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TLatex.h"
#include <cstdio>

#include "../../headers/plotting/plotStyle.h"

namespace {

const char *fB = "/home/clayton/Analysis/code/bJetRaaAnalysis/src/calculateBJetsPerZ/"
                 "rootFiles/BJetPbPbToPP/Data/histograms_BJetPbPbToPP.root";
const char *fI = "/home/clayton/Analysis/code/bJetRaaAnalysis/src/newFractionCalculation/"
                 "rootFiles/JetsPerZ/histograms_JetsPerZ_lightJets_rebinned.root";
const char *outPath = "/home/clayton/Analysis/code/bJetRaaAnalysis/figures/"
                      "flavorFraction/inclusiveBJetFraction.pdf";

const int   NSys = 5;
const char *sys  [NSys] = {"pp", "C1", "C2", "C3", "C4"};
const char *label[NSys] = {"pp", "PbPb 0-10%", "PbPb 10-30%", "PbPb 30-50%", "PbPb 50-80%"};
// Okabe-Ito; index 4 (yellow) skipped
const int   colIdx[NSys] = {0, 5, 2, 3, 6};
const int   mrkIdx[NSys] = {markFilledCircle, markFilledSquare, markFilledDiamond,
                            markOpenCircle, markOpenSquare};

const double quoteFromPt = 80.;

} // namespace

void plotInclusiveBJetFraction()
{
  initPlotStyle();

  TFile *fb = TFile::Open(fB);
  TFile *fi = TFile::Open(fI);
  if(!fb || fb->IsZombie()){ printf("cannot open %s\n", fB); return; }
  if(!fi || fi->IsZombie()){ printf("cannot open %s\n", fI); return; }

  TCanvas *c = new TCanvas("cFrac", "", 800, 700);
  c->SetLeftMargin(0.15); c->SetRightMargin(0.04);
  c->SetTopMargin(0.07);  c->SetBottomMargin(0.13);

  TH1D *fr = new TH1D("frFrac", "", 1, quoteFromPt, 500.);
  fr->SetStats(0);
  // 0-10% reaches 0.168 at 80-90 GeV. An earlier 0.14 maximum clipped that point
  // and the 0.157 at 90-100 straight off the figure without any indication --
  // the two largest values in the measurement were invisible.
  fr->GetYaxis()->SetRangeUser(0., 0.19);
  fr->GetXaxis()->SetTitle("jet p_{T} [GeV]");
  fr->GetYaxis()->SetTitle("inclusive b-jet fraction");
  fr->GetXaxis()->SetTitleOffset(1.15);
  fr->GetYaxis()->SetTitleOffset(1.35);
  fr->Draw("axis");

  // Upper right: every curve falls with pT, so that corner is the only region
  // free of points. At upper left the legend sat on top of the 0-10% curve.
  TLegend *leg = new TLegend(0.58, 0.62, 0.95, 0.88);
  leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.036);

  for(int s = 0; s < NSys; s++){
    TH1D *hb = (TH1D*) fb->Get(Form("bJetSpectrum_%s", sys[s]));
    TH1D *hi = (TH1D*) fi->Get(Form("inclSpectrum_%s", sys[s]));
    if(!hb || !hi){ printf("  missing histograms for %s\n", label[s]); continue; }

    printf("\n  %s\n     pT bin        b/Z/GeV     incl/Z/GeV     fraction\n", label[s]);

    TGraphErrors *g = new TGraphErrors();
    int n = 0;
    for(int i = 1; i <= hb->GetNbinsX(); i++){
      const double lo = hb->GetBinLowEdge(i), hi_ = hb->GetBinLowEdge(i+1);
      if(lo < quoteFromPt - 1e-6) continue;          // unfolding buffer bin
      const double b = hb->GetBinContent(i);
      if(b <= 0.) continue;

      // pair by edge; the inclusive axis has two extra bins at the bottom
      const int j = hi->FindBin(0.5*(lo + hi_));
      if(j < 1 || j > hi->GetNbinsX()) continue;
      if(fabs(hi->GetBinLowEdge(j) - lo) > 1e-6 ||
         fabs(hi->GetBinLowEdge(j+1) - hi_) > 1e-6) continue;
      const double d = hi->GetBinContent(j);
      if(d <= 0.) continue;

      const double f = b/d;
      const double rb = hb->GetBinError(i)/b;
      const double rd = hi->GetBinError(j)/d;
      const double ef = f*sqrt(rb*rb + rd*rd);      // treated as uncorrelated

      g->SetPoint(n, 0.5*(lo + hi_), f);
      g->SetPointError(n, 0.5*(hi_ - lo), ef);
      n++;
      printf("   %4.0f-%4.0f  %12.5g %14.5g %9.4f +- %.4f\n", lo, hi_, b, d, f, ef);
    }

    const int col = TColor::GetColor(okabeHex[colIdx[s]]);
    g->SetLineColor(col); g->SetMarkerColor(col);
    g->SetMarkerStyle(mrkIdx[s]); g->SetMarkerSize(1.3); g->SetLineWidth(2);
    g->Draw("pz same");
    leg->AddEntry(g, label[s], "lp");
  }
  leg->Draw();

  TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.034);
  t.DrawLatex(0.15, 0.945, "corrected b jets / corrected inclusive jets, each system");

  c->SaveAs(outPath);
  printf("\n  figure: %s\n", outPath);
}
