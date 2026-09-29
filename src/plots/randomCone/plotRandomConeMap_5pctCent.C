// Random-cone pT maps in (eta, phi), one panel per 5% centrality slice.
//
// INPUT. h_randConeEtaPhi_C{1..18} of PbPb_caloTowerAnalyzer.C: a TProfile2D of
// the mean cone pT at the (eta, phi) where the cone was thrown, 32 eta bins
// (0.1) x 64 phi bins over |eta| < 1.6, |phi| < pi. Slice i is centrality
// (i-1)*5 - i*5 %, covering the full ultraFine range 0-90%; C0 (the scan's
// inclusive bin) is not used. This is the same object PbPb_caloTowerAnalyzer.C
// loads as its underlying-event background map (bkgMapFileOverride in
// headers/AnalysisSetup/pseudoJets.h) once merged/pointed at.
//
// STATISTICS. This is the full MinBias Part1 scan (2026-09-28), 51-55M cones
// per slice over 2048 native cells -- ~25000 cones/cell on average, so unlike
// the 51-event bootstrap test this does NOT need rebinning to fill the map;
// (rebinEta, rebinPhi) default to 1 (no rebin) and are left as arguments only
// in case a future low-statistics file needs it again.
//
// Cells with no cones are left blank (white), which is different from a cell
// whose cones all had zero pT (drawn as the lowest colour). Each panel has its
// own colour scale, since the central-to-peripheral range is >100x and a
// shared scale would show only the first few panels; the scale maximum is
// printed in each panel.
//
// Colour map: kViridis (perceptually uniform, colour-blind safe).
//
// Usage: root -l -b -q 'plotRandomConeMap_5pctCent.C("/path/to/file.root")'
// Run from: src/plots/randomCone/

#include <cmath>
#include "TFile.h"
#include "TProfile2D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
#include "../../../headers/plotting/plotStyle.h"

void plotRandomConeMap_5pctCent(const char *path =
    "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/"
    "PbPb_MinBias_Part1_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_"
    "mixedEventPFClustering_fastJetResamples-100_pseudoJetCandPtMin-0.0_towers_"
    "etMin-0.30_etaCluster-3.0_bkgMapMaking_2026-9-28_ultraFineCentBins.root",
    int rebinEta = 1, int rebinPhi = 1)
{
  initPlotStyle();
  gStyle->SetPalette(kViridis);
  const char *outDir = "../../../figures/randomCone/";
  gSystem->mkdir(outDir, kTRUE);

  TFile *f = TFile::Open(path);
  if(!f || f->IsZombie()){ printf("ERROR: %s did not open\n", path); return; }

  const int nSl = 18;
  TCanvas *c = new TCanvas("cmap", "", 1500, 900);
  c->Divide(6, 3, 0.0, 0.0);

  printf("  rebin (%d, %d)\n  %-8s %6s %8s %9s %9s\n", rebinEta, rebinPhi,
         "cent", "cones", "filled", "min", "max");
  for(int i = 1; i <= nSl; i++){
    TProfile2D *p = nullptr; f->GetObject(Form("h_randConeEtaPhi_C%d", i), p);
    if(!p){ printf("ERROR: h_randConeEtaPhi_C%d missing\n", i); return; }
    TProfile2D *r = p->Rebin2D(rebinEta, rebinPhi, Form("rMap_%d", i));

    // mean per cell, with cells that hold no cones marked below the axis minimum
    TH2D *m = r->ProjectionXY(Form("mMap_%d", i), "B");   // "B": binomial-free, plain means
    m->SetDirectory(nullptr);
    int nFilled = 0; double lo = 1e30, hi = -1e30;
    for(int bx = 1; bx <= r->GetNbinsX(); bx++)
      for(int by = 1; by <= r->GetNbinsY(); by++){
        if(r->GetBinEntries(r->GetBin(bx, by)) > 0.){
          const double v = r->GetBinContent(bx, by);
          m->SetBinContent(bx, by, v); nFilled++;
          lo = TMath::Min(lo, v); hi = TMath::Max(hi, v);
        }
        else m->SetBinContent(bx, by, -1.);
      }
    printf("  %2d-%-5d %6.0f %5d/%-3d %9.2f %9.2f\n", (i-1)*5, i*5, p->GetEntries(),
           nFilled, r->GetNbinsX()*r->GetNbinsY(), nFilled ? lo : 0., nFilled ? hi : 0.);

    c->cd(i);
    gPad->SetLeftMargin(0.19); gPad->SetBottomMargin(0.18);
    gPad->SetRightMargin(0.22); gPad->SetTopMargin(0.05);
    m->SetTitle("");
    m->SetStats(0);
    m->SetMinimum(0.); m->SetMaximum(nFilled ? hi * 1.02 : 1.);
    m->GetXaxis()->SetTitle("#it{#eta}"); m->GetYaxis()->SetTitle("#it{#phi}");
    m->GetXaxis()->SetTitleSize(0.07); m->GetYaxis()->SetTitleSize(0.07);
    m->GetXaxis()->SetLabelSize(0.06); m->GetYaxis()->SetLabelSize(0.06);
    m->GetXaxis()->SetTitleOffset(1.0); m->GetYaxis()->SetTitleOffset(1.2);
    m->GetZaxis()->SetLabelSize(0.055);
    m->GetXaxis()->SetNdivisions(505); m->GetYaxis()->SetNdivisions(505);
    m->Draw("COLZ");

    TLatex l; l.SetNDC(); l.SetTextFont(42); l.SetTextColor(kBlack);
    l.SetTextSize(0.09); l.DrawLatex(0.24, 0.89, Form("%d - %d%%", (i-1)*5, i*5));
  }
  c->cd(0);
  TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.012);
  t.DrawLatex(0.005, 0.003, "PbPb MinBias Part1, full scan; mean cone #it{p}_{T} [GeV], own scale per panel, blank = no cones");
  savePdfTight(c, Form("%srandomConeMap_5pctCent_grid.pdf", outDir));
  printf("\n  figure in %s\n", outDir);
}
