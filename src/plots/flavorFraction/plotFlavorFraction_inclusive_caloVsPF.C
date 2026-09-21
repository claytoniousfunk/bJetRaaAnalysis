// Flavor make-up of INCLUSIVE jets vs jet pT, PYTHIA calo jets vs PF jets.
//
// Inclusive counterpart of plotFlavorFraction_muTagged_caloVsPF.C. Same samples,
// same categories, same binning and the same stacked presentation, so the two
// figures can be read side by side; the only change is the source histogram and
// therefore the jet population.
//
//   muon-tagged  h_muptrel_recoJetPt_inclRecoMuonTag_triggerOn_<flavor>_T0,
//                one TH2 per flavor, projected over the whole ptRel axis
//   inclusive    h_inclRecoJetPt_flavor, a single TH2 of jet pT against the
//                flavor INTEGER, sliced on the flavor axis
//
// The flavor axis is 27 bins over [-5, 22], one per integer, so a category is a
// sum over the bins holding its values. Bins are found with FindBin(v + 0.5) to
// land in the middle of the bin rather than on an edge.
//
// No muon requirement and no trigger gating here, so this is the flavor make-up
// of the jet sample itself. Comparing it against the muon-tagged version is
// what separates "the two flavor definitions disagree" from "muon tagging
// selects different jets in the two collections".
//
// The two samples use DIFFERENT flavor definitions, which is the point of the
// comparison:
//   calo  refparton_flavor -- matched parton; no b-hadron information, so no
//         gluon-splitting (bGS) label, and jets with no matched parton are x
//   PF    jtPartonFlavor + bHadronNumber -- bGS separated out
// They are also different PYTHIA samples (calo pThat > 15; PF the 2026-02
// pThat > 30 file with trigger-efficiency weights) and different jet
// collections, so not every difference is the flavor definition.
//
// The x category is not a small residual on the calo side: refparton_flavor
// leaves 27% of inclusive calo jets unassigned, against 0.1% for
// jtPartonFlavor. See src/plots/jetPt/jetCollection/plotCaloPFFlavorMatch_PYTHIA.C,
// which matches the two collections jet by jet in a forest that has both and
// finds that 30% of PF b jets are labeled gluons by refparton_flavor.
//
// Categories: b (b and bbar), bGS, c, uds, g, x (no flavor match). Fractions
// are of all inclusive jets in the bin, so they sum to 1 including x; the table
// prints that sum as a closure check.
//
// Usage: root -l -b -q 'plotFlavorFraction_inclusive_caloVsPF.C'
// Run from: src/plots/flavorFraction/

#include <cmath>
#include <vector>
#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TSystem.h"
#include "../../../headers/plotting/plotStyle.h"

const char *caloFile =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PYTHIA/"
  "PYTHIA_DiJet_caloJets_pThat-15_mu12_pTmu-15to999_tight_vzReweight_jetTrkMaxFilter_removeHYDJETjet0p35CutOnGen_WDecayFilter_2026-9-15.root";
const char *pfFile =
  "/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/scanningOutput/PYTHIA/latest/"
  "PYTHIA_DiJet_pThat-30_mu12_pTmu-15_tight_mu12TriggerEfficiencyCorrection_vzReweight_jetTrkMaxFilter_removeHYDJETjet0p35CutOnGen_2026-2-11.root";

const char *histName = "h_inclRecoJetPt_flavor";
const char *outDir   = "../../../figures/flavorFraction/";

const double edgeJetPt[] = {50, 60, 80, 100, 120, 150, 200, 300, 500};
const int    nEdgeJetPt  = sizeof(edgeJetPt)/sizeof(double) - 1;

// drawn bottom to top
const int NCat = 6;
const char *catLabel[NCat] = {"#it{b}", "#it{b} (g#rightarrow#it{b}#bar{#it{b}})", "#it{c}", "#it{uds}", "#it{g}", "unmatched (x)"};
const char *catTag[NCat]   = {"b", "bGS", "c", "uds", "g", "x"};
const char *catHex[NCat]   = {"#D55E00", "#E69F00", "#009E73", "#56B4E9", "#0072B2", "#CC79A7"};
// flavor integers summed into each category, as PYTHIA_scan.C assigns them:
// +-1..5 quarks, 17 = b from gluon splitting, 21 = gluon, 0 = no match
const std::vector<std::vector<int>> catFlavors = {
  {5, -5}, {17}, {4, -4}, {1, -1, 2, -2, 3, -3}, {21}, {0}
};

// jet counts per jet pT bin, summed over the flavor integers of one category
static TH1D* jetPtCounts(TH2D *H, const std::vector<int> &vals, const char *name)
{
  TH1D *sum = nullptr;
  for(int v : vals){
    int iy = H->GetYaxis()->FindBin(v + 0.5);   // mid-bin, never an edge
    if(iy < 1 || iy > H->GetNbinsY()) continue;
    TH1D *p = H->ProjectionX(Form("%s_v%d", name, v), iy, iy);
    p->SetDirectory(nullptr);
    if(!sum) sum = p;
    else { sum->Add(p); delete p; }
  }
  if(!sum) return nullptr;
  TH1D *r = rebinTo(sum, nEdgeJetPt, edgeJetPt, name);
  delete sum;
  return r;
}

// fractions per category; frac[c] is a TH1D over edgeJetPt
static bool fractions(TFile *f, const char *tag, TH1D **frac)
{
  TH2D *H = nullptr; f->GetObject(histName, H);
  if(!H){ printf("ERROR: %s missing in the %s file\n", histName, tag); return false; }

  // denominator: the whole flavor axis, underflow and overflow included, so a
  // flavor code outside the categories would show up as a sum below 1
  TH1D *pAll = H->ProjectionX(Form("all_%s", tag), 0, -1);
  pAll->SetDirectory(nullptr);
  TH1D *all = rebinTo(pAll, nEdgeJetPt, edgeJetPt, Form("allR_%s", tag));
  delete pAll;

  for(int c = 0; c < NCat; c++){
    frac[c] = new TH1D(Form("frac_%s_%s", catTag[c], tag), "", nEdgeJetPt, edgeJetPt);
    frac[c]->SetDirectory(nullptr);
    TH1D *h = jetPtCounts(H, catFlavors[c], Form("cnt_%s_%s", catTag[c], tag));
    if(h){ frac[c]->Add(h); delete h; }
    for(int b = 1; b <= nEdgeJetPt; b++){
      double n = all->GetBinContent(b);
      frac[c]->SetBinContent(b, n > 0. ? frac[c]->GetBinContent(b)/n : 0.);
      frac[c]->SetBinError(b, 0.);   // a stacked fraction plot draws no errors
    }
  }
  delete all;
  return true;
}

static void drawStack(TPad *pad, TH1D **frac, const char *title)
{
  pad->cd();
  // right margin holds the legend
  pad->SetLeftMargin(0.15); pad->SetRightMargin(0.30);
  pad->SetTopMargin(0.10);  pad->SetBottomMargin(0.12);

  // cumulative from the bottom; draw the tallest first so each lower layer
  // paints over the one above it
  TH1D *cum[NCat];
  for(int c = 0; c < NCat; c++){
    cum[c] = (TH1D*) frac[c]->Clone(Form("%s_cum", frac[c]->GetName()));
    cum[c]->SetDirectory(nullptr);
    if(c > 0) cum[c]->Add(cum[c-1]);
    int col = TColor::GetColor(catHex[c]);
    cum[c]->SetFillColor(col); cum[c]->SetLineColor(col); cum[c]->SetStats(0);
  }
  TH1D *frame = cum[NCat-1];
  frame->SetTitle("");
  frame->SetMinimum(0.); frame->SetMaximum(1.);
  frame->GetXaxis()->SetTitle("#it{p}_{T}^{jet} [GeV]");
  frame->GetXaxis()->SetTitleSize(0.050); frame->GetXaxis()->SetLabelSize(0.042);
  frame->GetXaxis()->SetNdivisions(505);   // the default crowds the labels in a narrow frame
  frame->GetYaxis()->SetTitle("fraction of inclusive jets");
  frame->GetYaxis()->SetTitleSize(0.050); frame->GetYaxis()->SetTitleOffset(1.35);
  frame->GetYaxis()->SetLabelSize(0.042);
  frame->Draw("hist");
  for(int c = NCat-2; c >= 0; c--) cum[c]->Draw("hist same");
  frame->Draw("AXIS same");

  TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.048);
  la.DrawLatex(0.15, 0.925, title);
}

void plotFlavorFraction_inclusive_caloVsPF()
{
  initPlotStyle();
  gSystem->mkdir(outDir, kTRUE);

  TFile *fC = TFile::Open(caloFile), *fP = TFile::Open(pfFile);
  if(!fC || fC->IsZombie()){ printf("ERROR: cannot open %s\n", caloFile); return; }
  if(!fP || fP->IsZombie()){ printf("ERROR: cannot open %s\n", pfFile); return; }

  TH1D *fracC[NCat], *fracP[NCat];
  if(!fractions(fC, "calo", fracC) || !fractions(fP, "pf", fracP)) return;

  // ---- table ----------------------------------------------------------------
  for(int s = 0; s < 2; s++){
    TH1D **fr = s ? fracP : fracC;
    printf("\n%s\n  %-9s", s ? "PF jets, jtPartonFlavor (PYTHIA pThat > 30, 2026-02)"
                             : "calo jets, refparton_flavor (PYTHIA pThat > 15, 2026-09-15)", "jet pT");
    for(int c = 0; c < NCat; c++) printf(" %7s", catTag[c]);
    printf(" %8s %6s\n", "b+bGS", "sum");
    for(int b = 1; b <= nEdgeJetPt; b++){
      printf("  %3.0f-%-5.0f", edgeJetPt[b-1], edgeJetPt[b]);
      double tot = 0.;
      for(int c = 0; c < NCat; c++){ printf(" %7.3f", fr[c]->GetBinContent(b)); tot += fr[c]->GetBinContent(b); }
      printf(" %8.3f %6.3f\n", fr[0]->GetBinContent(b) + fr[1]->GetBinContent(b), tot);
    }
  }
  printf("\n  'sum' is a closure check: it is 1 unless a flavor code fell outside the categories.\n");

  // ---- figures: one portrait canvas per sample --------------------------------
  // Two 700x800 canvases rather than one wide one: a landscape canvas came out
  // squashed into the top of the PDF page, where 700x800 renders cleanly.
  const char *sampleTag[2]   = {"calo", "PF"};
  const char *sampleTitle[2] = {"PYTHIA calo jets, #it{refparton_flavor}", "PYTHIA PF jets, #it{jtPartonFlavor}"};
  for(int s = 0; s < 2; s++){
    TCanvas *c = new TCanvas(Form("c_flav_%s", sampleTag[s]), "", 700, 800);
    TPad *pad = new TPad(Form("pad_%s", sampleTag[s]), "", 0, 0, 1, 1);
    pad->Draw();
    drawStack(pad, s ? fracP : fracC, sampleTitle[s]);

    pad->cd();
    TLegend *leg = makeLegend(0.71, 0.45, 0.99, 0.85, 0.036);
    for(int k = NCat-1; k >= 0; k--){   // top of the stack first
      TH1D *h = (TH1D*) fracC[k]->Clone(Form("leg_%s_%d", sampleTag[s], k));
      h->SetFillColor(TColor::GetColor(catHex[k])); h->SetLineColor(TColor::GetColor(catHex[k]));
      leg->AddEntry(h, catLabel[k], "f");
    }
    leg->Draw();
    TLatex lg; lg.SetNDC(); lg.SetTextFont(42); lg.SetTextSize(0.034);
    lg.DrawLatex(0.715, 0.87, "inclusive jets");

    c->SaveAs(Form("%sflavorFraction_inclusive_%s.pdf", outDir, sampleTag[s]));
    delete c;
  }
  printf("\nfigures in %s\n", outDir);
}
