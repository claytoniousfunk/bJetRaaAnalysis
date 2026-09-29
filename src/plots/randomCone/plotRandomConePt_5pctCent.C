// Random-cone pT distributions from the calo-tower scan, one panel per 5%
// centrality slice, plus the mean and RMS of the cone pT against centrality.
//
// INPUT. h_pseudoJetPt_C{1..18} of PbPb_caloTowerAnalyzer.C. The scan books 18
// slices of 10 hiBin units; hiBin is in 0.5% units, so slice i is centrality
// (i-1)*5 - i*5 %, covering the full ultraFine range 0-90%. C0 is the scan's
// own inclusive bin, not the sum of the slices, and is not used.
//
// WHAT THE CONE IS. Sum of tower E_T (above the scan's tower threshold, inside
// its clustering acceptance) within dR_max_pfcand of a random (eta, phi) with
// |eta| < 1.6. With doEventMixing on, the towers are drawn from a pool of
// N_mixedEventsInPool events of the SAME centrality slice, not from the event
// itself -- read the file's `provenance` TNamed for the flags that produced it.
// Each event throws N_mixedEventsInPool cones, all filled at the event weight,
// so dividing by the histogram's own integral (as here) is the per-cone density.
//
// STATISTICS. This is the full MinBias Part1 scan (2026-09-28), ~515-545k
// events per slice and 51-55M cones per slice -- not the earlier 51-event
// bootstrap test. Two things still worth remembering when reading it:
//  * Cones from one event share the same mixed-event donor pool, so the 100
//    cones per event are NOT independent draws. The error bars here are
//    sqrt(N) of the cone counts and so understate the true uncertainty --
//    though with 500k+ events per slice this is a much smaller effect than it
//    was on the bootstrap test.
//  * The histogram has 5 GeV bins, 0-500 GeV. In the most peripheral slices
//    the mean cone pT is under 1 GeV, so the shape sits in the first bin or
//    two even with full statistics; the mean and RMS remain exact (the
//    histogram carries its own sums) and are the reliable peripheral content.
//
// Usage: root -l -b -q 'plotRandomConePt_5pctCent.C("/path/to/file.root")'
// Run from: src/plots/randomCone/

#include <cmath>
#include "TFile.h"
#include "TH1D.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLatex.h"
#include "TSystem.h"
#include "TNamed.h"
#include "../../../headers/plotting/plotStyle.h"

void plotRandomConePt_5pctCent(const char *path =
    "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/"
    "PbPb_MinBias_Part1_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_"
    "mixedEventPFClustering_fastJetResamples-100_pseudoJetCandPtMin-0.0_towers_"
    "etMin-0.30_etaCluster-3.0_bkgMapMaking_2026-9-28_ultraFineCentBins.root")
{
  initPlotStyle();
  const char *outDir = "../../../figures/randomCone/";
  gSystem->mkdir(outDir, kTRUE);

  TFile *f = TFile::Open(path);
  if(!f || f->IsZombie()){ printf("ERROR: %s did not open\n", path); return; }

  const int nSl = 18;
  TH1D *h[nSl+1] = {nullptr};
  double nEvt[nSl+1] = {0};
  for(int i = 1; i <= nSl; i++){
    f->GetObject(Form("h_pseudoJetPt_C%d", i), h[i]);
    TH1D *nt = nullptr; f->GetObject(Form("h_nTower_C%d", i), nt);   // one entry per event
    if(!h[i] || !nt){ printf("ERROR: C%d histograms missing\n", i); return; }
    nEvt[i] = nt->GetEntries();
  }
  if(TNamed *p = (TNamed*) f->Get("provenance"))
    printf("  provenance: %.60s...\n", p->GetTitle());
  else printf("  no provenance TNamed in this file (pre-2026-09-09 scan)\n");

  // ------------------------- 6x3 grid of distributions ----------------------
  TCanvas *cg = new TCanvas("cg", "", 1500, 900);
  cg->Divide(6, 3, 0.0, 0.0);

  printf("\n  %-8s %8s %10s %9s %9s\n", "cent", "events", "cones", "mean", "RMS");
  for(int i = 1; i <= nSl; i++){
    cg->cd(i);
    gPad->SetLeftMargin(0.20); gPad->SetBottomMargin(0.18);
    gPad->SetRightMargin(0.04); gPad->SetTopMargin(0.04);

    TH1D *d = (TH1D*) h[i]->Clone(Form("d_%d", i));
    d->SetDirectory(nullptr);
    const double N = d->Integral();
    const double mean = h[i]->GetMean(), rms = h[i]->GetRMS();   // exact: from the stored sums
    printf("  %2d-%-5d %8.0f %10.0f %9.2f %9.2f\n", (i-1)*5, i*5, nEvt[i], N, mean, rms);
    d->Scale(1./N);
    d->Scale(1., "width");

    // x range: through the last filled bin plus one, at least 4 bins
    int last = d->GetNbinsX();
    while(last > 1 && h[i]->GetBinContent(last) <= 0.) last--;
    const int hiBin = TMath::Max(4, TMath::Min(d->GetNbinsX(), last + 1));
    d->GetXaxis()->SetRange(1, hiBin);

    d->SetTitle("");
    styleH(d, hexData, markFilledCircle, 0.7);
    d->GetXaxis()->SetTitle("cone #it{p}_{T} [GeV]");
    d->GetYaxis()->SetTitle("1/N dN/d#it{p}_{T} [GeV^{-1}]");
    d->GetXaxis()->SetTitleSize(0.065); d->GetYaxis()->SetTitleSize(0.065);
    d->GetXaxis()->SetLabelSize(0.058); d->GetYaxis()->SetLabelSize(0.058);
    d->GetYaxis()->SetTitleOffset(1.35); d->GetXaxis()->SetTitleOffset(1.15);
    d->GetXaxis()->SetNdivisions(505); d->GetYaxis()->SetNdivisions(505);
    d->SetMinimum(0.); d->SetMaximum(d->GetMaximum() * 2.1);
    d->Draw("E");

    TLatex l; l.SetNDC(); l.SetTextFont(42);
    l.SetTextSize(0.095); l.DrawLatex(0.28, 0.87, Form("%d - %d%%", (i-1)*5, i*5));
    l.SetTextSize(0.062);
    l.DrawLatex(0.28, 0.77, Form("mean %.1f GeV", mean));
    l.DrawLatex(0.28, 0.68, Form("RMS %.1f GeV", rms));
  }
  cg->cd(0);
  TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.012);
  t.DrawLatex(0.005, 0.003, "PbPb MinBias Part1, full scan; error bars are sqrt(N) of correlated (mixed-pool) cones, 5 GeV bins");
  savePdfTight(cg, Form("%srandomConePt_5pctCent_grid.pdf", outDir));
  delete cg;

  // ------------------------- mean and RMS vs centrality --------------------
  TGraphErrors *gM = new TGraphErrors(nSl), *gR = new TGraphErrors(nSl);
  for(int i = 1; i <= nSl; i++){
    const double x = (i-0.5)*5.;
    gM->SetPoint(i-1, x, h[i]->GetMean()); gM->SetPointError(i-1, 2.5, 0.);
    gR->SetPoint(i-1, x, h[i]->GetRMS());  gR->SetPointError(i-1, 2.5, 0.);
  }
  TCanvas *cm = new TCanvas("cm", "", 900, 700);
  cm->SetLeftMargin(0.15); cm->SetBottomMargin(0.13); cm->SetRightMargin(0.05);
  TH1D *fr = new TH1D("fr", "", 18, 0., 90.);
  fr->SetStats(0); fr->SetMinimum(0.); fr->SetMaximum(60.);
  fr->GetXaxis()->SetTitle("centrality [%]"); fr->GetYaxis()->SetTitle("cone #it{p}_{T} [GeV]");
  fr->GetXaxis()->SetTitleSize(0.05); fr->GetYaxis()->SetTitleSize(0.05);
  fr->GetXaxis()->SetLabelSize(0.042); fr->GetYaxis()->SetLabelSize(0.042);
  fr->GetYaxis()->SetTitleOffset(1.3);
  fr->Draw();
  const int cM = TColor::GetColor(hexData), cR = TColor::GetColor("#D55E00");
  gM->SetMarkerStyle(markFilledCircle); gM->SetMarkerColor(cM); gM->SetLineColor(cM);
  gR->SetMarkerStyle(markFilledSquare); gR->SetMarkerColor(cR); gR->SetLineColor(cR);
  gM->SetMarkerSize(1.1); gR->SetMarkerSize(1.0); gM->SetLineWidth(2); gR->SetLineWidth(2);
  gM->Draw("P SAME"); gR->Draw("P SAME");
  TLegend *lg = makeLegend(0.60, 0.72, 0.93, 0.85);
  lg->AddEntry(gM, "mean", "pl"); lg->AddEntry(gR, "RMS", "pl");
  lg->Draw();
  TLatex tl; tl.SetNDC(); tl.SetTextSize(0.038);
  tl.DrawLatex(0.19, 0.86, "PbPb random cones, calo towers");
  tl.DrawLatex(0.19, 0.80, "MinBias Part1, full scan");
  savePdfTight(cm, Form("%srandomConePt_meanRMS_vsCent.pdf", outDir));
  printf("\n  figures in %s\n", outDir);
}
