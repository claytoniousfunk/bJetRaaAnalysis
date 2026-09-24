// Calo jets per Z, PbPb/pp, one pad per 10% centrality class.
//
// The 3x3 counterpart of plotJetsPerZ_caloVsPF_noUnfold.C's per-class panels,
// with the PF curve absent: there are no PbPb HardProbes (Jet80/Jet100) PF
// scans in ultra-fine centrality bins, only calo ones, so a like-for-like
// calo-vs-PF comparison cannot be built at this binning. The two apparent PF
// matches on disk are PbPb_MinBias_Part1_hiBinReweightToHardProbesJet80_*,
// which are MinBias files reweighted to the HardProbes hiBin distribution, not
// triggered samples. Adding the PF side needs a PbPb HardProbes Jet80 + Jet100
// scan with centrality_ultraFineCentBins.h.
//
// Input is jetsPerZ_caloJets_10pctBins.C's output: trigger-only stitch (MinBias
// dropped from both systems), manual AK4Calo JEC on every spectrum, N_Z from
// the 2026-09-24 ultra-fine SingleMuon scan.
//
// 60-90% carry 218 / 80 / 22 Z and are labelled statistics-limited rather than
// hidden. Their y range is the same as everyone else's here, so the 80-90%
// top bin runs off the pad -- that is deliberate, it is a 1.29 +- 0.25 point on
// 22 Z and should not set the scale for eight good classes.
//
// Usage, from src/newFractionCalculation/:
//   root -l -b -q plotJetsPerZ_caloJets_10pctPanels.C

#include "../../headers/plotting/plotStyle.h"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TGraphErrors.h"
#include "TLatex.h"
#include "TLine.h"
#include <cstdio>

namespace {

const char *inPath = "./rootFiles/JetsPerZ/jetsPerZ_caloJets_10pctBins.root";
const char *outPath = "/home/clayton/Analysis/code/bJetRaaAnalysis/figures/"
                      "JetsPerZ/JetsPerZ_caloJets_10pctPanels.pdf";

const int   NCls = 9;
const char *tag  [NCls] = {"0to10","10to20","20to30","30to40","40to50",
                           "50to60","60to70","70to80","80to90"};
const char *label[NCls] = {"0-10%","10-20%","20-30%","30-40%","40-50%",
                           "50-60%","60-70%","70-80%","80-90%"};
const double nZ  [NCls] = {7275, 5042, 3267, 1939, 1046, 502, 218, 80, 22};
const long  minZ = 400;

const double ptLo = 80., ptHi = 500.;

} // namespace

void plotJetsPerZ_caloJets_10pctPanels()
{
  initPlotStyle();

  TFile *f = TFile::Open(inPath);
  if(!f || f->IsZombie()){
    printf("ERROR: run jetsPerZ_caloJets_10pctBins.C first (%s)\n", inPath);
    return;
  }

  TCanvas *c = new TCanvas("cPanels", "", 1200, 1150);

  for(int k = 0; k < NCls; k++){
    const int col = k % 3, row = k / 3;
    const double x0 = col/3., x1 = (col+1)/3.;
    const double y1 = 1. - row/3., y0 = 1. - (row+1)/3.;

    c->cd();
    TPad *p = new TPad(Form("p%d", k), "", x0, y0, x1, y1);
    p->SetLeftMargin(0.17); p->SetRightMargin(0.03);
    p->SetTopMargin(0.06);  p->SetBottomMargin(0.15);
    p->Draw(); p->cd();

    TH1D *fr = new TH1D(Form("fr%d", k), "", 1, ptLo, ptHi);
    fr->SetStats(0);
    fr->GetYaxis()->SetRangeUser(0., 1.25);
    fr->GetXaxis()->SetTitle("calo jet p_{T} [GeV]");
    fr->GetYaxis()->SetTitle("(jets per Z)_{PbPb} / (jets per Z)_{pp}");
    fr->GetXaxis()->SetTitleSize(0.055); fr->GetYaxis()->SetTitleSize(0.050);
    fr->GetXaxis()->SetLabelSize(0.050); fr->GetYaxis()->SetLabelSize(0.050);
    fr->GetXaxis()->SetTitleOffset(1.15);
    fr->GetYaxis()->SetTitleOffset(1.55);
    fr->GetYaxis()->SetNdivisions(505);
    fr->Draw("axis");

    TLine *one = new TLine(ptLo, 1., ptHi, 1.);
    one->SetLineStyle(2); one->SetLineColor(kGray+2); one->Draw();

    TH1D *h = (TH1D*) f->Get(Form("ratio_%s", tag[k]));
    if(h){
      // points at bin centres, empty bins dropped rather than drawn at zero
      TGraphErrors *g = new TGraphErrors();
      int n = 0;
      for(int i = 1; i <= h->GetNbinsX(); i++){
        if(h->GetBinContent(i) <= 0.) continue;
        const double lo = h->GetBinLowEdge(i), hi = h->GetBinLowEdge(i+1);
        g->SetPoint(n, 0.5*(lo+hi), h->GetBinContent(i));
        g->SetPointError(n, 0.5*(hi-lo), h->GetBinError(i));
        n++;
      }
      const bool thin = nZ[k] < minZ;
      const int col0 = TColor::GetColor(okabeHex[thin ? 6 : 0]);
      g->SetLineColor(col0); g->SetMarkerColor(col0);
      g->SetMarkerStyle(thin ? markOpenCircle : markFilledCircle);
      g->SetMarkerSize(1.3); g->SetLineWidth(2);
      g->Draw("pz same");
    }

    TLatex t; t.SetNDC(); t.SetTextFont(42);
    t.SetTextSize(0.070); t.DrawLatex(0.23, 0.86, label[k]);
    t.SetTextSize(0.045);
    t.DrawLatex(0.23, 0.79, Form("N_{Z} = %.0f", nZ[k]));
    if(nZ[k] < minZ){
      t.SetTextSize(0.040);
      t.DrawLatex(0.23, 0.73, "statistics-limited");
    }
  }

  c->SaveAs(outPath);
  printf("\n  figure: %s\n", outPath);
}
