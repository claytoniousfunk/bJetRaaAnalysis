// Random-cone energy maps with and without the 2 GeV PF-candidate cut.
//
// h_randConeEtaPhi is a TProfile2D: in each (eta, phi) cell, the mean pT of
// R = 0.4 random cones thrown there. These are the maps the RC-subtracted
// FastJet spectra subtract, so they are what bkgMapFile() in pseudoJets.h
// selects by PF-candidate cut.
//
// Per coarse centrality class, three maps:
//   no cut   pseudoJetCandPt_min = 0, same-event MinBias scan of 2026-08-17
//            (the map used for 0 GeV scans)
//   2 GeV    pseudoJetCandPt_min = 2, same-event MinBias scan of 2026-09-15
//   ratio    2 GeV / no cut, cell by cell -- flat if the cut only rescales the
//            background, structured if it changes the eta-phi shape
// Coarse classes are built by adding the fine-slice profiles, which averages
// cones correctly (TProfile2D::Add weights by entries).
//
// Each map has its own color scale: the two cuts differ by a factor of ~5, so
// a shared scale would leave the 2 GeV map a single color.
//
// Usage: root -l -b -q 'plotRCMaps_pfCandPtCut.C'
// Run from: src/plots/pfCand/

#include <cmath>
#include "TMath.h"
#include "TFile.h"
#include "TH2D.h"
#include "TProfile2D.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
#include "../../../headers/plotting/coarseCent.h"
#include "../../../headers/plotting/plotStyle.h"

const char *dir = "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/";
const char *fileNoCut = "PbPb_MinBias_Part1_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_sameEventPFClustering_pseudoJetCandPtMin-0.0_2026-8-17_ultraFineCentBins.root";
const char *fileCut   = "PbPb_MinBias_Part1_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_sameEventPFClustering_pseudoJetCandPtMin-2.0_2026-9-15_ultraFineCentBins.root";

const char *mapName = "h_randConeEtaPhi";
const char *outDir  = "../../../figures/pfCandPtCut/";

// mean-cone-pT map for a coarse class, as a plain TH2D of cell means
static TH2D* classMap(TFile *f, int ci, const char *tag)
{
  TH2D *sum = coarseSum2(f, mapName, ci, tag);   // a TProfile2D underneath
  if(!sum){ printf("ERROR: %s missing for %s in %s\n", mapName, coarseLabel[ci], tag); return nullptr; }
  TProfile2D *prof = dynamic_cast<TProfile2D*>(sum);
  if(!prof){ printf("ERROR: %s is not a TProfile2D\n", mapName); delete sum; return nullptr; }
  TH2D *m = prof->ProjectionXY(Form("map_%s_%d", tag, ci));
  m->SetDirectory(nullptr);
  delete sum;
  return m;
}

// mean and relative RMS of the cell means (how non-uniform the map is)
static void cellStats(TH2D *h, double &mean, double &relRms)
{
  double s = 0., s2 = 0.; int n = 0;
  for(int ix = 1; ix <= h->GetNbinsX(); ix++)
    for(int iy = 1; iy <= h->GetNbinsY(); iy++){
      double v = h->GetBinContent(ix, iy); s += v; s2 += v*v; n++;
    }
  mean = n ? s/n : 0.;
  double var = n ? s2/n - mean*mean : 0.;
  relRms = (mean > 0. && var > 0.) ? sqrt(var)/mean : 0.;
}

static void drawMap(TH2D *h, const char *zTitle, const char *line1, const char *line2,
                    int palette, const char *outPath)
{
  gStyle->SetPalette(palette);
  TCanvas *c = new TCanvas(Form("c_%s", h->GetName()), "", 700, 800);
  c->SetLeftMargin(0.13); c->SetRightMargin(0.22);
  c->SetTopMargin(0.16);  c->SetBottomMargin(0.12);
  h->SetTitle("");
  h->SetStats(0);
  h->GetXaxis()->SetTitle("#it{#eta}");  h->GetXaxis()->SetTitleSize(0.045); h->GetXaxis()->SetLabelSize(0.040);
  h->GetYaxis()->SetTitle("#it{#phi}");  h->GetYaxis()->SetTitleSize(0.045); h->GetYaxis()->SetLabelSize(0.040);
  h->GetYaxis()->SetTitleOffset(1.1);
  h->GetZaxis()->SetTitle(zTitle); h->GetZaxis()->SetTitleSize(0.040); h->GetZaxis()->SetLabelSize(0.035);
  h->GetZaxis()->SetTitleOffset(1.75);
  h->GetZaxis()->SetMaxDigits(3);
  h->Draw("COLZ");

  TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.036);
  la.DrawLatex(0.13, 0.94, "PbPb 5.02 TeV MinBias, same-event");
  la.DrawLatex(0.13, 0.895, "random cones (#it{R} = 0.4)");
  la.DrawLatex(0.13, 0.85, Form("%s, %s", line1, line2));

  c->SaveAs(outPath);
  delete c;
}

void plotRCMaps_pfCandPtCut()
{
  initPlotStyle();
  gSystem->mkdir(outDir, kTRUE);

  TFile *fN = TFile::Open(Form("%s%s", dir, fileNoCut));
  TFile *fC = TFile::Open(Form("%s%s", dir, fileCut));
  if(!fN || fN->IsZombie() || !fC || fC->IsZombie()){ printf("ERROR: cannot open the map files\n"); return; }

  printf("\nrandom-cone energy maps: mean cone pT over the eta-phi cells, and the cell-to-cell spread\n");
  printf("  %-7s | %9s %8s | %9s %8s | %11s %8s\n",
         "class", "no cut", "rel.RMS", "2 GeV", "rel.RMS", "ratio mean", "rel.RMS");

  for(int ci = 0; ci < NCoarse; ci++){
    TH2D *mN = classMap(fN, ci, "noCut");
    TH2D *mC = classMap(fC, ci, "cut2");
    if(!mN || !mC) continue;

    TH2D *mR = (TH2D*) mC->Clone(Form("ratio_%d", ci));
    mR->SetDirectory(nullptr);
    mR->Divide(mN);

    double aN, sN, aC, sC, aR, sR;
    cellStats(mN, aN, sN); cellStats(mC, aC, sC); cellStats(mR, aR, sR);
    printf("  %-7s | %9.2f %8.3f | %9.2f %8.3f | %11.4f %8.3f\n",
           coarseLabel[ci], aN, sN, aC, sC, aR, sR);

    drawMap(mN, "mean random-cone #it{p}_{T} [GeV]", coarseLabel[ci], "no PF-candidate cut",
            kViridis, Form("%srcMap_noCut_%s.pdf", outDir, coarseTag[ci]));
    drawMap(mC, "mean random-cone #it{p}_{T} [GeV]", coarseLabel[ci], "#it{p}_{T}^{cand} > 2 GeV",
            kViridis, Form("%srcMap_cut2_%s.pdf", outDir, coarseTag[ci]));
    drawMap(mR, "2 GeV / no cut", coarseLabel[ci], "ratio of mean cone #it{p}_{T}",
            kCividis, Form("%srcMap_ratio_%s.pdf", outDir, coarseTag[ci]));

    delete mN; delete mC; delete mR;
  }

  printf("\nfigures in %s\n", outDir);
}
