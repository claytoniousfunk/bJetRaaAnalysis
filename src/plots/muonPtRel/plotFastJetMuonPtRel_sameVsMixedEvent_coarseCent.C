// Same-event vs. mixed-event fastJet muon ptRel (h_fastJetMuonPtRel_fastJetPt_PF_bkgSub_RC_C*),
// merged into the standard 4 coarse centrality classes and split into three
// jet pT windows, with a mixed/same ratio panel added.
//
// Adapted from plotFastJetMuonPtRel_sameVsMixedEvent.cc, which does this for a
// single hardcoded centBin/jetPt combination with no ratio. Same normalisation
// choice as that macro: NORMALISED BY N_EVENTS (h_vz integral), not by the
// histogram's own integral -- normalizeByIntegral() forces both curves to the
// same area and would make a mixed/same ratio meaningless (it measures shape
// only). Per-event normalisation is what makes the ratio panel a real yield
// comparison: same-event is genuine muon+jet production, mixed-event is the
// combinatorial background estimate built from randomly paired events, so the
// ratio panel (mixed/same) is the combinatorial-background fraction of the
// measured yield per pT bin -- how much of what you see is fake.
//
// SCAN-VERSION REQUIREMENT (added 2026-09-08). Both files must come from the
// same generation of PbPb_pfCandAnalyzer.C. Three commits that day changed what
// h_fastJetMuonPtRel_fastJetPt_PF_bkgSub_RC contains:
//
//   e52a4682  jet pT on the Y axis moved from raw rcSub to JEC-corrected, so a
//             "50 < pT < 60" window selects a different jet population
//             (JEC ~ 1.16 here). ptRel itself is unaffected -- getPtRel divides
//             by |jet|^2, so the magnitude cancels -- but the binning is not.
//   3490cf89  the tagging muon became the LEADING muon constituent instead of
//             whichever came last in constituent order.
//   3490cf89  T2/T3 gained a trigger gate (does not touch this histogram).
//
// Mixing generations folds those into the same/mixed difference this plot is
// meant to isolate: comparing the raw-scale 2026-9-1 mixed file against the
// JEC-scale 2026-9-8 one, same dataset and same histogram, the ptRel mean in a
// fixed window moves by up to 2.5x. checkScanGeneration() below refuses that
// combination rather than drawing it.
//
// The marker is whether the file carries h_mixedMuonPtRel_recoJetPt: that key
// is booked and written unconditionally by the new code, and is simply empty in
// a same-event scan, so its PRESENCE identifies the generation regardless of
// doEventMixing.
//
// INPUT FILES. Same-event: 2026-8-28 (matches the reference macro). Mixed-
// event: 2026-9-1, NOT the 2026-8-27 file the reference macro uses. Checked
// before choosing it: the 2026-8-27 mixed-event file (despite carrying no
// "_partial" suffix) has only 29539 events in slice C1 against the same-event
// file's 731408 -- a ~4% subsample, almost certainly the run that crashed
// partway through (see the "evt frac: 4%" crash log from that investigation).
// The 2026-9-1 file has 732005 events in C1, matching same-event to within
// 0.1%, and is the file produced after the null-pointer/segfault fixes. Do
// not point this at the 2026-8-27 file without re-checking its event count.
//
// STATISTICS. Same-event has thousands of entries per class/pT combination;
// mixed-event has as few as 1-5 (see the printed table). The ratio panel is
// only as good as the SMALLER of the two, so treat the 30-50% and 50-80%
// panels, and anything above 80 GeV, as indicative shapes rather than
// measurements. Bins where the same-event (denominator) projection has zero
// entries are left blank in the ratio rather than drawn as a spurious 0.
//
// ptRel rebinning is coarser than the reference macro's 11-bin scheme,
// because that scheme was tuned to same-event statistics; with mixed-event
// entries in the single digits it would put multiple points at exactly zero.
//
// Usage: root -l -b -q 'plotFastJetMuonPtRel_sameVsMixedEvent_coarseCent.C'
// Run from: src/plots/muonPtRel/

#include "../../../headers/functions/divideByBinwidth.h"

const char *sameFile =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/PbPb_SingleMuon_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_sameEventPFClustering_pseudoJetCandPtMin-0.0_2026-8-28_ultraFineCentBins.root";
const char *mixedFile =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/PbPb_SingleMuon_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_mixedEventPFClustering_pseudoJetCandPtMin-0.0_2026-9-1_ultraFineCentBins.root";

const char *histBase = "h_fastJetMuonPtRel_fastJetPt_PF_bkgSub_RC";

const char *outDir = "../../../figures/muonPtRel/";

// coarse class -> [first, last] ultra-fine slice index, standard convention
const int NClass = 4;
const int sliceLo[NClass] = { 1,  3,  7, 11};
const int sliceHi[NClass] = { 2,  6, 10, 16};
const char *classLabel[NClass] = { "0-10%", "10-30%", "30-50%", "50-80%" };
const char *classTag[NClass]   = { "0to10pct", "10to30pct", "30to50pct", "50to80pct" };

const int    NPt = 3;
const double ptLo[NPt] = {50, 60,  80};
const double ptHi[NPt] = {60, 80, 120};

const char *sameHex  = "#0072B2";   // Okabe-Ito blue
const char *mixedHex = "#D55E00";   // Okabe-Ito vermillion

const int    NEdge = 6;
double       edgeAxis[NEdge] = {0, 0.5, 1.0, 1.5, 2.5, 5.0};

// Same-event and mixed-event yields differ by 1-2 orders of magnitude
// depending on class/pT (mixed is the rare combinatorial estimate), so the
// ratio range is computed per panel below rather than fixed -- a single
// global range could not fit both a 0-10% panel (mixed/same ~0.02-0.13) and
// a depleted peripheral one in the same axis without clipping one of them.

static TH2D* mergeCoarse(TFile *f, int ci, double &Nevents){
  TH2D *sum = nullptr;
  Nevents = 0.;
  for(int si = sliceLo[ci]; si <= sliceHi[ci]; si++){
    TH2D *h = nullptr; TH1D *v = nullptr;
    f->GetObject(Form("%s_C%d", histBase, si), h);
    f->GetObject(Form("h_vz_C%d", si), v);
    if(!h || !v){ printf("WARNING: slice C%d missing from this file\n", si); continue; }
    if(!sum){ sum = (TH2D*) h->Clone(Form("sum_%d_%p", ci, (void*)f)); sum->SetDirectory(nullptr); }
    else sum->Add(h);
    Nevents += v->Integral();
  }
  return sum;
}

// true if the file was written by the post-2026-09-08 scan
bool isNewGeneration(TFile *f)
{
  return f->GetListOfKeys()->FindObject("h_mixedMuonPtRel_recoJetPt_C1") != nullptr;
}

// sameOverride / mixedOverride let a fresh scan be dropped in without editing
// the constants above.
void plotFastJetMuonPtRel_sameVsMixedEvent_coarseCent(const char *sameOverride  = nullptr,
                                                      const char *mixedOverride = nullptr,
                                                      const char *outSuffix     = "")
{
  gStyle->SetOptStat(0);
  gSystem->mkdir(outDir, kTRUE);

  const char *usedSame  = (sameOverride  && sameOverride[0])  ? sameOverride  : sameFile;
  const char *usedMixed = (mixedOverride && mixedOverride[0]) ? mixedOverride : mixedFile;

  TFile *fS = TFile::Open(gSystem->ExpandPathName(usedSame));
  TFile *fM = TFile::Open(gSystem->ExpandPathName(usedMixed));
  if(!fS || fS->IsZombie()){ printf("ERROR: cannot open %s\n", usedSame); return; }
  if(!fM || fM->IsZombie()){ printf("ERROR: cannot open %s\n", usedMixed); return; }

  // Refuse, do not warn: a plot drawn from mismatched generations looks
  // perfectly reasonable and there is nothing on it to tell the reader that
  // most of the same/mixed difference is a jet pT scale and a muon-choice
  // change rather than physics.
  bool newS = isNewGeneration(fS), newM = isNewGeneration(fM);
  printf("same-event : %s\n  generation: %s\n", usedSame,  newS ? "post-2026-09-08" : "pre-2026-09-08");
  printf("mixed-event: %s\n  generation: %s\n", usedMixed, newM ? "post-2026-09-08" : "pre-2026-09-08");
  if(newS != newM){
    printf("\nREFUSING TO PLOT: the two scans are from different generations of\n"
           "PbPb_pfCandAnalyzer.C, so h_fastJetMuonPtRel_fastJetPt_PF_bkgSub_RC\n"
           "does not mean the same thing in each (jet pT axis raw vs JEC, and a\n"
           "different choice of tagging muon). The same/mixed difference this\n"
           "plot is for would be confounded with both. See the header.\n\n"
           "Fix: rerun the %s scan with the current code.\n",
           newS ? "MIXED-event" : "SAME-event");
    return;
  }

  int colS = TColor::GetColor(sameHex);
  int colM = TColor::GetColor(mixedHex);

  printf("%-9s %8s %8s", "class", "N(same)", "N(mix)");
  for(int p = 0; p < NPt; p++) printf("  %14s", Form("%.0f-%.0f (S/M)", ptLo[p], ptHi[p]));
  printf("\n");

  for(int ci = 0; ci < NClass; ci++){

    double NevS = 0., NevM = 0.;
    TH2D *hS = mergeCoarse(fS, ci, NevS);
    TH2D *hM = mergeCoarse(fM, ci, NevM);
    if(!hS || !hM){ printf("WARNING: class %s missing, skipping\n", classLabel[ci]); continue; }

    printf("%-9s %7.0fk %7.0fk", classLabel[ci], NevS/1000., NevM/1000.);

    for(int p = 0; p < NPt; p++){
      int b1 = hS->GetYaxis()->FindBin(ptLo[p] + 1e-6);
      int b2 = hS->GetYaxis()->FindBin(ptHi[p] - 1e-6);

      TH1D *rawS = hS->ProjectionX(Form("rawS_%d_%d", ci, p), b1, b2);
      TH1D *rawM = hM->ProjectionX(Form("rawM_%d_%d", ci, p), b1, b2);
      int nS = (int) rawS->GetEntries();
      int nM = (int) rawM->GetEntries();
      printf("  %6d/%-6d", nS, nM);

      TH1D *pS = (TH1D*) rawS->Rebin(NEdge-1, Form("pS_%d_%d", ci, p), edgeAxis);
      TH1D *pM = (TH1D*) rawM->Rebin(NEdge-1, Form("pM_%d_%d", ci, p), edgeAxis);
      pS->SetDirectory(nullptr); pM->SetDirectory(nullptr);

      // per-event yield density, NOT shape-normalised -- see header
      pS->Scale(1./NevS); divideByBinwidth(pS);
      pM->Scale(1./NevM); divideByBinwidth(pM);

      TH1D *ratio = (TH1D*) pM->Clone(Form("ratio_%d_%d", ci, p));
      ratio->SetDirectory(nullptr);
      ratio->Divide(pS);
      // blank bins where the same-event (denominator) had no entries --
      // Divide() would otherwise silently leave content 0, indistinguishable
      // from a genuine null result
      for(int b = 1; b <= pS->GetNbinsX(); b++){
        if(pS->GetBinContent(b) <= 0.){ ratio->SetBinContent(b, 0.); ratio->SetBinError(b, 0.); }
      }

      // per-panel ratio range: fit every drawn point (value + error), floor
      // at 0. Deliberately NOT forced to include 1 -- mixed/same (the
      // combinatorial fraction) sits at 0.02-0.15 in the classes with good
      // statistics, and forcing headroom out to 1 would squash every real
      // point into the bottom of the pad. The y=1 line is still drawn below
      // and simply falls off the top of the sparser panels, which is itself
      // informative (mixed/same has not been shown to be small there).
      double ratioMax = 0.05;
      for(int b = 1; b <= ratio->GetNbinsX(); b++){
        double v = ratio->GetBinContent(b), e = ratio->GetBinError(b);
        if(v <= 0.) continue;
        if(v + e > ratioMax) ratioMax = v + e;
      }
      ratioMax *= 1.2;

      const double lm = 0.15, rm = 0.05, tm = 0.10, bm = 0.28, split = 0.35;
      TCanvas *c = new TCanvas(Form("c_%d_%d", ci, p), "", 700, 800);

      TPad *pUp = new TPad("pUp", "", 0, split, 1, 1);
      pUp->SetLeftMargin(lm); pUp->SetRightMargin(rm);
      pUp->SetTopMargin(tm);  pUp->SetBottomMargin(0.);
      pUp->Draw(); pUp->cd();

      double ymax = TMath::Max(pS->GetMaximum(), pM->GetMaximum());

      pS->SetLineColor(colS); pS->SetMarkerColor(colS);
      pS->SetMarkerStyle(20); pS->SetMarkerSize(1.0); pS->SetLineWidth(2);
      pS->SetTitle("");
      pS->GetXaxis()->SetLabelSize(0.);
      pS->GetYaxis()->SetTitle("(1/N_{evt}) dN^{#mu+jet}/dp_{T}^{rel}  [GeV^{-1}]");
      pS->GetYaxis()->SetTitleSize(0.052); pS->GetYaxis()->SetLabelSize(0.045);
      pS->GetYaxis()->SetTitleOffset(1.30);
      // plain decimals, not "xE-3" -- that label sits in the same top-left
      // corner as the title text drawn below and collides with it
      pS->GetYaxis()->SetNoExponent(kTRUE);
      pS->SetMinimum(0.); pS->SetMaximum(ymax * 1.5);
      pS->Draw("ep");

      pM->SetLineColor(colM); pM->SetMarkerColor(colM);
      pM->SetMarkerStyle(21); pM->SetMarkerSize(1.0); pM->SetLineWidth(2);
      pM->Draw("ep same");

      TLegend *leg = new TLegend(0.42, 0.66, 0.90, 0.86);
      leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.036);
      leg->AddEntry(pS, Form("same event  (n=%d)", nS), "lp");
      leg->AddEntry(pM, Form("mixed event  (n=%d)", nM), "lp");
      leg->Draw();

      TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.034);
      la.DrawLatex(0.19, 0.90, Form("PbPb SingleMuon (5.02 TeV), %s", classLabel[ci]));
      la.DrawLatex(0.19, 0.855, Form("%.0f < p_{T}^{jet} < %.0f GeV, |#eta^{jet}| < 1.6", ptLo[p], ptHi[p]));
      la.DrawLatex(0.19, 0.81, "p_{T}^{#mu} > 15 GeV, |#eta^{#mu}| < 2");

      c->cd();
      TPad *pDn = new TPad("pDn", "", 0, 0, 1, split);
      pDn->SetLeftMargin(lm); pDn->SetRightMargin(rm);
      pDn->SetTopMargin(0.);  pDn->SetBottomMargin(bm);
      pDn->Draw(); pDn->cd();

      const double sc = (1.-split)/split;
      ratio->SetLineColor(kBlack); ratio->SetMarkerColor(kBlack);
      ratio->SetMarkerStyle(20); ratio->SetMarkerSize(0.9); ratio->SetLineWidth(2);
      ratio->SetTitle("");
      ratio->GetXaxis()->SetTitle("muon p_{T}^{rel} [GeV]");
      ratio->GetYaxis()->SetTitle("mixed / same");
      ratio->GetXaxis()->SetTitleSize(0.052*sc); ratio->GetXaxis()->SetLabelSize(0.045*sc);
      ratio->GetYaxis()->SetTitleSize(0.052*sc); ratio->GetYaxis()->SetLabelSize(0.045*sc);
      ratio->GetYaxis()->SetTitleOffset(1.30/sc);
      ratio->GetYaxis()->SetNdivisions(505);
      ratio->SetMinimum(0.); ratio->SetMaximum(ratioMax);
      ratio->Draw("ep");

      TLine *one = new TLine(0, 1., 5., 1.);
      one->SetLineStyle(2); one->SetLineColor(kGray+2); one->Draw();

      TString out = TString(outDir) +
        Form("muonPtRel_sameVsMixed_coarseCent_%s_pt%.0fto%.0f%s.pdf", classTag[ci], ptLo[p], ptHi[p], outSuffix);
      c->SaveAs(out);
    }
    printf("\n");
  }
  printf("\nWrote plots to %s\n", outDir);
}
