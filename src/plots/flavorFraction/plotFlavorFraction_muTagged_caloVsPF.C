// Flavour make-up of muon-tagged jets vs jet pT, PYTHIA calo jets vs PF jets.
//
// Stacked fractions per analysis jet pT bin, from the per-flavour ptRel
// templates each PYTHIA scan writes (tagged + trigger on, nominal jet pT):
//   h_muptrel_recoJetPt_inclRecoMuonTag_triggerOn_<flavour>_T0
// projected onto jet pT over the full ptRel axis. These are exactly the jets
// the b-fraction template fits are built from.
//
// The two samples use DIFFERENT flavour definitions, which is the point of the
// comparison:
//   calo  refparton_flavor -- matched parton; no b-hadron information, so no
//         gluon-splitting (bGS) label, and jets with no matched parton are x
//   PF    jtPartonFlavor + bHadronNumber -- bGS separated out
// They are also different PYTHIA samples (calo pThat > 15; PF the 2026-02
// pThat > 30 file with trigger-efficiency weights) and different jet
// collections, so not every difference is the flavour definition.
//
// Categories: b (b and bbar), bGS, c, uds, g, x (no flavour match). Fractions
// are of all tagged jets in the bin, so they sum to 1 including x.
//
// Usage: root -l -b -q 'plotFlavorFraction_muTagged_caloVsPF.C'
// Run from: src/plots/flavorFraction/

#include <cmath>
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

const char *histFmt = "h_muptrel_recoJetPt_inclRecoMuonTag_triggerOn_%s_T0";
const char *outDir  = "../../../figures/flavorFraction/";

const double edgeJetPt[] = {50, 60, 80, 100, 120, 150, 200, 300, 500};
const int    nEdgeJetPt  = sizeof(edgeJetPt)/sizeof(double) - 1;

// drawn bottom to top
const int NCat = 6;
const char *catLabel[NCat] = {"#it{b}", "#it{b} (g#rightarrow#it{b}#bar{#it{b}})", "#it{c}", "#it{uds}", "#it{g}", "unmatched (x)"};
const char *catTag[NCat]   = {"b", "bGS", "c", "uds", "g", "x"};
const char *catHex[NCat]   = {"#D55E00", "#E69F00", "#009E73", "#56B4E9", "#0072B2", "#CC79A7"};
// template flavour names summed into each category
const std::vector<std::vector<const char*>> catFlavours = {
  {"bJets"}, {"bGSJets"}, {"cJets"}, {"uJets","dJets","sJets"}, {"gJets"}, {"xJets"}
};

// weighted tagged-jet count per jet pT bin for one template flavour
static TH1D* jetPtCounts(TFile *f, const char *flavour, const char *tag)
{
  TH2D *H = nullptr; f->GetObject(Form(histFmt, flavour), H);
  if(!H){ printf("WARNING: %s missing in %s file, treated as empty\n", Form(histFmt, flavour), tag); return nullptr; }
  TH1D *p = H->ProjectionY(Form("jpt_%s_%s", flavour, tag), 0, -1);
  p->SetDirectory(nullptr);
  TH1D *r = rebinTo(p, nEdgeJetPt, edgeJetPt, Form("jptR_%s_%s", flavour, tag));
  delete p;
  return r;
}

// fractions per category; frac[c] is a TH1D over edgeJetPt
static bool fractions(TFile *f, const char *tag, TH1D **frac)
{
  TH1D *all = jetPtCounts(f, "allJets", tag);
  if(!all){ printf("ERROR: no allJets histogram in %s file\n", tag); return false; }
  for(int c = 0; c < NCat; c++){
    frac[c] = new TH1D(Form("frac_%s_%s", catTag[c], tag), "", nEdgeJetPt, edgeJetPt);
    frac[c]->SetDirectory(nullptr);
    for(const char *fl : catFlavours[c]){
      TH1D *h = jetPtCounts(f, fl, tag);
      if(h){ frac[c]->Add(h); delete h; }
    }
    for(int b = 1; b <= nEdgeJetPt; b++){
      double n = all->GetBinContent(b);
      frac[c]->SetBinContent(b, n > 0. ? frac[c]->GetBinContent(b)/n : 0.);
      frac[c]->SetBinError(b, 0.);
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
  frame->GetYaxis()->SetTitle("fraction of #mu-tagged jets");
  frame->GetYaxis()->SetTitleSize(0.050); frame->GetYaxis()->SetTitleOffset(1.35);
  frame->GetYaxis()->SetLabelSize(0.042);
  frame->Draw("hist");
  for(int c = NCat-2; c >= 0; c--) cum[c]->Draw("hist same");
  frame->Draw("AXIS same");

  TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.048);
  la.DrawLatex(0.15, 0.925, title);
}

void plotFlavorFraction_muTagged_caloVsPF()
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
    printf(" %8s\n", "b+bGS");
    for(int b = 1; b <= nEdgeJetPt; b++){
      printf("  %3.0f-%-5.0f", edgeJetPt[b-1], edgeJetPt[b]);
      for(int c = 0; c < NCat; c++) printf(" %7.3f", fr[c]->GetBinContent(b));
      printf(" %8.3f\n", fr[0]->GetBinContent(b) + fr[1]->GetBinContent(b));
    }
  }

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
    lg.DrawLatex(0.715, 0.87, "#mu-tagged, trig. on");

    c->SaveAs(Form("%sflavorFraction_muTagged_%s.pdf", outDir, sampleTag[s]));
    delete c;
  }
  printf("\nfigures in %s\n", outDir);
}
