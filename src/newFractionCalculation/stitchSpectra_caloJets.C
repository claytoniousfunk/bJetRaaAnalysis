// Per-centrality optimal binning, stitched jet spectra, and the PbPb/pp
// jets-per-Z ratio, for calo jets. MinBias + Jet80 + Jet100 on BOTH sides, so
// the ratio reaches the scans' 50 GeV jet cut.
//
// ------------------------------------------------ pp ENERGY SCALE (READ) ---
// The pp MinBias scan is the manualJEC one (2026-09-17): jet pT from the
// Spring18 ppRef AK4Calo files applied to rawpt. Every other calo sample here
// still takes the forest's jtpt, and the forests correct calo jets by a
// PF-sized factor. Corrected/raw, quantile-matched on the same jets:
//      pp MinBias, manualJEC   1.79 @60  1.58 @100  1.43 @150  1.32 @200
//      pp MinBias, forest      1.87      1.69       1.53       1.38
//      pp Jet80/100, SM        1.15      1.13       1.12       1.12
//      PbPb (50-80%)           1.11      1.09       1.09       1.09
// So the manual correction did NOT close the gap -- it confirms the proper calo
// scale and leaves the triggers, which have not been rescanned, ~40% low. In
// raw pT MinBias/Jet80 is flat; in corrected pT it falls 45% from 75 to 175 GeV.
// Stitching across that is knowingly inconsistent (the user chose to stitch
// as-is, 2026-09-17). pp MinBias is normalized to Jet80 just above the Jet80
// threshold, which keeps the splice continuous; the run prints what a high-pT
// window would give. This is fixed by rescanning the pp trigger samples (and
// PbPb) with useManualJEC, after which the thresholds must be re-derived.
//
// Run deriveTriggerThresholds_caloJets.C first: this macro consumes the
// thresholds it writes and uses ONLY the samples that file validates, so a
// trigger whose turn-on could not be measured is never silently trusted.
//
// ------------------------------------------------------------- STITCHING ---
// Each sample is used from its own threshold up to the threshold of the next
// one. The highest-threshold sample sets the scale (it is the least prescaled),
// and each lower sample is normalized to the one above it over the overlap
// window where BOTH are known flat: [T_above, T_above + overlapWidth]. That is
// the same chain the existing measurement uses, with the windows now measured
// rather than assumed.
//
// -------------------------------------------------------------- BINNING ----
// Derived per centrality class, independently, as the user asked: starting at
// that class's floor, 5 GeV bins are merged until the PbPb yield in the bin
// BOTH reaches the target relative statistical error AND is no narrower than
// the jet energy resolution. Edges stay on the 5 GeV grid of the scans.
//
// The resolution floor matters more than the statistical one here. These are
// trigger samples in trigger units, so a single 5 GeV bin at 110 GeV already
// holds ~1e5 counts: on statistics alone the binning came out 5 GeV wide
// almost everywhere, giving 0.2-2% errors per bin and implying a pT resolution
// the detector does not have. minWidthFrac is that floor, as a fraction of pT.
// It is a placeholder at the typical calo-jet resolution until the response
// matrices exist and it can be set from them -- it is deliberately the one
// number in this macro that is asserted rather than measured.
//
// The pp spectrum is then rebinned onto EACH class's edges before the ratio, so
// numerator and denominator always share binning even though the classes do
// not share it with each other.
//
// ------------------------------------------------------------ WHAT IS NOT --
// No unfolding (no calo response matrices yet) and no fake-jet subtraction (the
// estimate is FastJet on PF candidates; there is no calo equivalent). Both are
// stated on the figures.
//
// N_Z is jet-collection independent, so it comes from the same files as before:
// the calo SingleMuon scan in pp, and the PF SingleMuon scan in PbPb, which is
// the only complete Z sample there.
//
// Usage, from src/newFractionCalculation/:
//   root -l -b -q 'stitchSpectra_caloJets.C'

#include "../../headers/functions/divideByBinwidth.h"
#include "../../headers/functions/triggerStitch.h"
#include "../../headers/plotting/plotStyle.h"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TBox.h"
#include "TColor.h"
#include "TMath.h"
#include "TSystem.h"
#include <cstdio>
#include <vector>
#include <map>

const char *repo = "/home/clayton/Analysis/code/bJetRaaAnalysis";
const char *sib  = "/home/clayton/Analysis/code/bJetMuonTaggingAnalysis";  // read-only

const char *clsLabel[5] = {"", "0-10%", "10-30%", "30-50%", "50-80%"};
const char *clsTag[5]   = {"", "0to10pct", "10to30pct", "30to50pct", "50to80pct"};

// binning target: merge 5 GeV bins until the PbPb yield carries this relative
// statistical error, then cut a bin
const double targetRelErr = 0.05;
const double minWidthFrac = 0.10;   // bin never narrower than this fraction of pT
const double maxBinWidth  = 150.;
const double ptCeiling    = 500.;
const double overlapWidth = 40.;   // normalization window above a threshold

// Use each centrality class's PbPb thresholds for the pp samples too, so the
// numerator and denominator of every bin come from the same kind of sample
// (MinBias/MinBias, Jet80/Jet80, Jet100/Jet100). Legal because every class
// threshold sits above the pp validated one -- pp Jet80 is flat from 70 GeV and
// Jet100 from 90, against class thresholds of 85-110 and 115-140 -- so each pp
// sample is still only used where it was validated; it is just handed over to
// the next sample later than it needs to be. Set false to let pp switch at its
// own thresholds instead.
const bool matchPPToClass = true;

const double Z_lo = 75., Z_hi = 105.;
const double muEff_pp = 0.9708, muEff[5] = {1., 0.8627, 0.9069, 0.9856, 0.9778};

static bool ok = true;

static TFile* openOrFail(TString p)
{
  TFile *f = TFile::Open(p);
  if(!f || f->IsZombie()){ printf("ERROR: cannot open %s\n", p.Data()); ok = false; return nullptr; }
  return f;
}

static TH1D* getOrFail(TFile *f, const char *name, const char *newName)
{
  if(!f) return nullptr;
  TH1D *h = nullptr;
  f->GetObject(name, h);
  if(!h){ printf("ERROR: %s missing from %s\n", name, f->GetName()); ok = false; return nullptr; }
  TH1D *c = (TH1D*) h->Clone(newName);
  c->SetDirectory(nullptr);
  return c;
}

static TH1D* toSpectrum(TH1D *h, const char *name, const std::vector<double> &edges)
{
  if(!h || edges.size() < 2) return nullptr;
  TH1D *r = (TH1D*) h->Rebin((int)edges.size() - 1, name, &edges[0]);
  r->SetDirectory(nullptr);
  divideByBinwidth(r);
  r->Scale(1./3.2);   // eta range
  return r;
}

void stitchSpectra_caloJets(bool useMinBias = true)
{
  initPlotStyle();
  ok = true;

  // useMinBias = false: trigger samples only, so each side starts at its own
  // Jet80 threshold (pp 70 GeV, PbPb 85-110 by class) instead of 50. Drops the
  // reach, but also drops the one sample whose energy scale disagrees -- pp
  // MinBias is the only manualJEC scan, everything else still uses forest jtpt.
  const char *tag = useMinBias ? "" : "_trigOnly";

  TString figDir = Form("%s/figures/JetsPerZ", repo);
  gSystem->mkdir(figDir, kTRUE);
  TString outDir = "./rootFiles/JetsPerZ";
  gSystem->mkdir(outDir, kTRUE);

  std::map<TString, double> thr = stitchReadThresholds(Form("%s/triggerThresholds_caloJets.txt", outDir.Data()));
  if(thr.empty()) return;

  // ------------------------------------------------------------- inputs ---
  TFile *fpp[3], *fPb[3];
  const char *tn[3] = {"Jet60", "Jet80", "Jet100"};
  for(int i = 0; i < 3; i++){
    fpp[i] = openOrFail(Form("%s/rootFiles/scanningOuput/pp/pp_HighEGJet_caloJets_%sHLT_mu12_pTmu-15to999_tight_deltaR-40_jetTrkMaxFilter_WDecayFilter_2026-9-15.root", repo, tn[i]));
    fPb[i] = openOrFail(Form("%s/rootFiles/scanningOuput/PbPb/PbPb_HardProbes_caloJets_%sHLT_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_2026-9-15.root", repo, tn[i]));
  }
  TFile *f_pp_mu = openOrFail(Form("%s/rootFiles/scanningOuput/pp/pp_SingleMuon_caloJets_mu12_pTmu-15to999_tight_deltaR-40_jetTrkMaxFilter_WDecayFilter_2026-9-15.root", repo));
  TFile *f_Pb_mu = openOrFail(Form("%s/rootFiles/scanningOutput/PbPb/latest/PbPb_SingleMuon_mu12_pTmu-15to999_tight_mu12TriggerEfficiencyCorrection_jetTrkMaxFilter_WDecayFilter_2026-5-4.root", sib));
  // MinBias: ungated, so it needs no threshold and supplies everything below
  // Jet80. The pp file is the 2026-09-17 scan of the fixed ZeroBias forest.
  TFile *f_pp_MB = nullptr, *f_Pb_MB = nullptr;
  if(useMinBias){
    f_pp_MB = openOrFail(Form("%s/rootFiles/scanningOuput/pp/pp_MinBias_caloJets_manualJEC_mu12_pTmu-15to999_tight_deltaR-40_jetTrkMaxFilter_WDecayFilter_2026-9-17.root", repo));
    f_Pb_MB = openOrFail(Form("%s/rootFiles/scanningOuput/PbPb/PbPb_MinBias_Part1_caloJets_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_2026-9-14.root", repo));
    if(!ok) return;
    for(TFile *f : {f_pp_MB, f_Pb_MB}){
      TH1D *g = nullptr; f->GetObject("h_nEventsNoJetTrigSel", g);
      if(g && g->GetMean() > 0.5){ printf("ERROR: %s has a jet trigger applied\n", f->GetName()); return; }
    }
  }
  if(!ok) return;
  const double minBiasFloor = 50.;

  // -------------------------------------------------------- pp samples --
  // prepared here, stitched per class below: with matchPPToClass the pp
  // hand-over points are the class's, so the stitch differs class by class
  printf("\n=== pp: %svalidated trigger samples ===\n", useMinBias ? "MinBias + " : "");
  std::vector<StitchSample> ppBase;
  if(useMinBias){
    TH1D *h = getOrFail(f_pp_MB, "h_inclRecoJetPt", "ppSpec_MinBias");
    if(!h) return;
    ppBase.push_back({"MinBias", minBiasFloor, h, 0.});
    printf("  %-7s validated from %.0f GeV (ungated)\n", "MinBias", minBiasFloor);
  }
  for(int i = 0; i < 3; i++){
    TString key = Form("pp/inclusive/%s", tn[i]);
    if(!thr.count(key)){ printf("  %-7s not validated, not used\n", tn[i]); continue; }
    TH1D *h = getOrFail(fpp[i], "h_inclRecoJetPt", Form("ppSpec_%s", tn[i]));
    if(!h) return;
    ppBase.push_back({tn[i], thr[key], h, 0.});
    printf("  %-7s validated from %.0f GeV\n", tn[i], thr[key]);
  }
  std::sort(ppBase.begin(), ppBase.end(), [](const StitchSample &a, const StitchSample &b){ return a.thr < b.thr; });
  const double ppFloor = ppBase.front().thr;

  // PbPb scales, measured once in a common window (a dataset property)
  double kPbMB  = useMinBias
    ? stitchCommonScale(f_Pb_MB, fPb[1], "h_inclRecoJetPt_C%d", 4, 150., 200., "PbPb Jet80 / MinBias", true)
    : 1.;
  double kPbJ80 = stitchCommonScale(fPb[1],  fPb[2], "h_inclRecoJetPt_C%d", 4, 150., 200., "PbPb Jet100 / Jet80", true);
  if(kPbMB <= 0. || kPbJ80 <= 0.){ printf("ERROR: could not measure the PbPb scales\n"); return; }

  // ---------------------------------------------------------- PbPb stitch --
  TH1D *pbStitch[5] = {nullptr, nullptr, nullptr, nullptr, nullptr};
  TH1D *ppStitch[5] = {nullptr, nullptr, nullptr, nullptr, nullptr};
  double pbFloor[5] = {0., 0., 0., 0., 0.};
  std::vector<double> edges[5];
  for(int c = 1; c <= 4; c++){
    printf("\n=== PbPb %s ===\n", clsLabel[c]);
    std::vector<StitchSample> s;
    if(useMinBias){
      TH1D *h = getOrFail(f_Pb_MB, Form("h_inclRecoJetPt_C%d", c), Form("pbSpec_MinBias_C%d", c));
      if(!h) return;
      s.push_back({"MinBias", minBiasFloor, h, kPbMB});
      printf("  %-7s used from %.0f GeV (ungated)\n", "MinBias", minBiasFloor);
    }
    for(int i = 0; i < 3; i++){
      TString key = Form("PbPb/%s/%s", clsTag[c], tn[i]);
      if(!thr.count(key)){ printf("  %-7s not validated, not used\n", tn[i]); continue; }
      TH1D *h = getOrFail(fPb[i], Form("h_inclRecoJetPt_C%d", c), Form("pbSpec_%s_C%d", tn[i], c));
      if(!h) return;
      s.push_back({tn[i], thr[key], h, strcmp(tn[i], "Jet80") == 0 ? kPbJ80 : 0.});
      printf("  %-7s used from %.0f GeV\n", tn[i], thr[key]);
    }
    std::sort(s.begin(), s.end(), [](const StitchSample &a, const StitchSample &b){ return a.thr < b.thr; });
    pbStitch[c] = stitchTriggerSamples(s, Form("pbStitch_C%d", c), overlapWidth, ptCeiling, true);
    if(!pbStitch[c]){ ok = false; continue; }
    pbFloor[c]  = s.front().thr;

    // pp for this class: same hand-over points as the PbPb chain above when
    // matchPPToClass, otherwise pp's own validated thresholds
    std::vector<StitchSample> ppS = ppBase;
    if(matchPPToClass){
      for(size_t j = 0; j < ppS.size(); j++){
        TString key = Form("PbPb/%s/%s", clsTag[c], ppS[j].trig);
        if(thr.count(key)){
          if(thr[key] < ppS[j].thr){
            printf("  WARNING: pp %s would start at %.0f GeV, below its validated %.0f; keeping the pp value\n",
                   ppS[j].trig, thr[key], ppS[j].thr);
          }
          else ppS[j].thr = thr[key];
        }
      }
      std::sort(ppS.begin(), ppS.end(), [](const StitchSample &a, const StitchSample &b){ return a.thr < b.thr; });
    }
    printf("  pp:");
    for(size_t j = 0; j < ppS.size(); j++) printf("  %s from %.0f", ppS[j].trig, ppS[j].thr);
    printf("\n");
    ppStitch[c] = stitchTriggerSamples(ppS, Form("ppStitch_C%d", c), overlapWidth, ptCeiling, true);
    if(!ppStitch[c]){ ok = false; continue; }

    // the ratio can only start where BOTH sides are valid, and the binning must
    // hold the precision target on both: below the pp Jet80 threshold the pp
    // side is MinBias and is the sparser of the two
    double floorPt = TMath::Max(pbFloor[c], ppFloor);
    edges[c] = stitchDeriveEdges({pbStitch[c], ppStitch[c]}, floorPt, ptCeiling,
                                 targetRelErr, minWidthFrac, maxBinWidth);
    if(edges[c].size() < 2){ printf("  ERROR: could not derive binning\n"); ok = false; continue; }
    printf("  binning from %.0f GeV (%zu bins):", floorPt, edges[c].size() - 1);
    for(size_t e = 0; e < edges[c].size(); e++) printf(" %.0f", edges[c][e]);
    printf("\n");
  }
  if(!ok) return;

  // ------------------------------------------------------------- Z yields --
  TH1D *dimuPP = getOrFail(f_pp_mu, "h_dimuonMass", "dimuPP");
  if(!dimuPP) return;
  double NZ_pp = dimuPP->Integral(dimuPP->GetXaxis()->FindBin(Z_lo), dimuPP->GetXaxis()->FindBin(Z_hi)) / (muEff_pp*muEff_pp);
  double NZ[5] = {0.,0.,0.,0.,0.};
  for(int c = 1; c <= 4; c++){
    TH1D *d = getOrFail(f_Pb_mu, Form("h_dimuonMass_C%d", c), Form("dimuPb_%d", c));
    if(!d) return;
    NZ[c] = d->Integral(d->GetXaxis()->FindBin(Z_lo), d->GetXaxis()->FindBin(Z_hi)) / (muEff[c]*muEff[c]);
  }
  printf("\n  N_Z: pp %.1f", NZ_pp);
  for(int c = 1; c <= 4; c++) printf(", %s %.1f", clsLabel[c], NZ[c]);
  printf("\n");

  // --------------------------------------------------------------- ratio ---
  TH1D *ratio[5] = {nullptr,nullptr,nullptr,nullptr,nullptr};
  TH1D *spPb[5]  = {nullptr,nullptr,nullptr,nullptr,nullptr};
  TH1D *spPP[5]  = {nullptr,nullptr,nullptr,nullptr,nullptr};
  for(int c = 1; c <= 4; c++){
    spPb[c] = toSpectrum(pbStitch[c], Form("spPb_C%d", c), edges[c]);
    spPP[c] = toSpectrum(ppStitch[c], Form("spPP_C%d", c), edges[c]);   // pp matched to this class
    if(!spPb[c] || !spPP[c]){ printf("ERROR: rebinning failed in %s\n", clsLabel[c]); return; }
    spPb[c]->Scale(1./NZ[c]);
    spPP[c]->Scale(1./NZ_pp);
    ratio[c] = (TH1D*) spPb[c]->Clone(Form("ratio_C%d", c));
    ratio[c]->SetDirectory(nullptr);
    ratio[c]->Divide(spPb[c], spPP[c], 1, 1, "");
  }

  printf("\n  jets per Z, PbPb / pp, each class on its own binning\n");
  for(int c = 1; c <= 4; c++){
    printf("    %s\n", clsLabel[c]);
    for(int b = 1; b <= ratio[c]->GetNbinsX(); b++)
      printf("      %4.0f-%-5.0f  %7.3f +- %.3f\n",
             ratio[c]->GetXaxis()->GetBinLowEdge(b), ratio[c]->GetXaxis()->GetBinUpEdge(b),
             ratio[c]->GetBinContent(b), ratio[c]->GetBinError(b));
  }

  // ------------------------------------------------------------- figures ---
  const char *hexC[5] = {"", okabeHex[5], okabeHex[1], okabeHex[2], okabeHex[0]};
  const int   markC[5] = {0, markFilledCircle, markFilledSquare, markFilledDiamond, markCross};

  // 1. stitched spectra, one pad per class, with the sample boundaries marked
  {
    TCanvas *cv = new TCanvas("c_spectra", "", 700, 800);
    cv->Divide(2, 2, 0.001, 0.001);
    for(int c = 1; c <= 4; c++){
      cv->cd(c); gPad->SetLogy();
      gPad->SetLeftMargin(0.17); gPad->SetBottomMargin(0.13); gPad->SetTopMargin(0.08); gPad->SetRightMargin(0.04);
      TH1D *h = (TH1D*) spPb[c]->Clone(Form("draw_C%d", c));
      styleH(h, hexC[c], markC[c]);
      h->SetTitle("");
      h->GetXaxis()->SetTitle("jet #it{p}_{T} [GeV]");
      h->GetXaxis()->SetTitleSize(0.055); h->GetXaxis()->SetLabelSize(0.050);
      h->GetYaxis()->SetTitle("#frac{1}{#it{N}_{Z}} #frac{d#it{N}}{d#it{p}_{T}d#eta}");
      h->GetYaxis()->SetTitleSize(0.055); h->GetYaxis()->SetLabelSize(0.050);
      h->GetYaxis()->SetTitleOffset(1.35);
      h->Draw("E1");
      TH1D *hp = (TH1D*) spPP[c]->Clone(Form("drawPP_C%d", c));
      // neutral grey, not hexData: the 50-80% class is drawn in black and the
      // two series were indistinguishable in that panel
      styleH(hp, "#7F7F7F", markOpenCircle);
      hp->Draw("E1 same");
      TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.052);
      la.DrawLatex(0.30, 0.90, Form("PbPb %s", clsLabel[c]));
      if(c == 1){
        TLegend *lg = makeLegend(0.45, 0.68, 0.96, 0.86, 0.048);
        lg->AddEntry(h,  "PbPb", "lp");
        lg->AddEntry(hp, "pp (matched bins)", "lp");
        lg->Draw();
        // bottom left: the spectrum falls to the lower right and the legend
        // holds the upper right, so this corner is the free one
        la.SetTextSize(0.040);
        la.DrawLatex(0.22, 0.20, "no unfolding, no fakes subtracted");
      }
    }
    cv->SaveAs(Form("%s/stitchedSpectra_caloJets%s.pdf", figDir.Data(), tag));
    delete cv;
  }

  // 2. the ratio, one pad per class since the binning differs between them
  {
    TCanvas *cv = new TCanvas("c_ratio", "", 700, 800);
    cv->Divide(2, 2, 0.001, 0.001);
    for(int c = 1; c <= 4; c++){
      cv->cd(c);
      gPad->SetLeftMargin(0.17); gPad->SetBottomMargin(0.13); gPad->SetTopMargin(0.08); gPad->SetRightMargin(0.04);
      styleH(ratio[c], hexC[c], markC[c]);
      ratio[c]->SetTitle("");
      ratio[c]->SetMinimum(0.); ratio[c]->SetMaximum(1.4);
      ratio[c]->GetXaxis()->SetTitle("jet #it{p}_{T} [GeV]");
      ratio[c]->GetXaxis()->SetTitleSize(0.055); ratio[c]->GetXaxis()->SetLabelSize(0.050);
      ratio[c]->GetYaxis()->SetTitle("jets per #it{Z}, PbPb / pp");
      ratio[c]->GetYaxis()->SetTitleSize(0.052); ratio[c]->GetYaxis()->SetLabelSize(0.050);
      ratio[c]->GetYaxis()->SetTitleOffset(1.40);
      ratio[c]->Draw("E1");
      TLine *one = new TLine(ratio[c]->GetXaxis()->GetXmin(), 1., ratio[c]->GetXaxis()->GetXmax(), 1.);
      one->SetLineStyle(7); one->Draw();
      ratio[c]->Draw("E1 same");
      TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.052);
      la.DrawLatex(0.26, 0.86, Form("PbPb %s", clsLabel[c]));
      la.SetTextSize(0.044);
      la.DrawLatex(0.26, 0.80, Form("from %.0f GeV", ratio[c]->GetXaxis()->GetXmin()));
      // stated once, in the first pad: a Divide(2,2) canvas cannot carry a
      // header of its own, since the sub-pads are opaque and drawn in front.
      // y = 0.73 is R_CP ~ 1.1 here, above the curves and below the label.
      if(c == 1){
        la.SetTextSize(0.040);
        la.DrawLatex(0.26, 0.73, "no unfolding, no fakes subtracted");
      }
    }
    cv->SaveAs(Form("%s/jetsPerZ_caloJets_stitched%s.pdf", figDir.Data(), tag));
    delete cv;
  }

  // 2b. the same ratios overlaid on one canvas.
  // Each class carries its own binning and its own floor (the higher of that
  // class's and pp's Jet80 threshold), so the frame is drawn explicitly from the
  // lowest floor of the four -- taking the axis from whichever curve happened to
  // be drawn first would clip the classes that reach lower.
  {
    double xLo = 1e30;
    for(int c = 1; c <= 4; c++)
      if(ratio[c]) xLo = TMath::Min(xLo, ratio[c]->GetXaxis()->GetBinLowEdge(1));

    TCanvas *cv = new TCanvas("c_ratio_overlay", "", 700, 800);
    cv->SetLeftMargin(0.16); cv->SetBottomMargin(0.12);
    cv->SetTopMargin(0.13);  cv->SetRightMargin(0.05);

    TH1D *frame = new TH1D("ratioOverlayFrame", "", 1, xLo, ptCeiling);
    frame->SetStats(0);
    frame->SetMinimum(0.); frame->SetMaximum(1.4);
    frame->GetXaxis()->SetTitle("jet #it{p}_{T} [GeV]");
    frame->GetXaxis()->SetTitleSize(0.045); frame->GetXaxis()->SetLabelSize(0.040);
    frame->GetYaxis()->SetTitle("#frac{1}{#it{N}_{Z}^{PbPb}} #frac{d#it{N}_{jet}^{PbPb}}{d#it{p}_{T}} #scale[2.5]{/} #frac{1}{#it{N}_{Z}^{pp}} #frac{d#it{N}_{jet}^{pp}}{d#it{p}_{T}}");
    frame->GetYaxis()->SetTitleSize(0.038); frame->GetYaxis()->SetLabelSize(0.040);
    frame->GetYaxis()->SetTitleOffset(1.65);
    frame->Draw("axis");

    TLine *one = new TLine(xLo, 1., ptCeiling, 1.);
    one->SetLineStyle(7);
    one->Draw();

    // upper right, above the unity line: the curves all sit below ~0.85
    TLegend *leg = makeLegend(0.60, 0.68, 0.92, 0.86, 0.036);
    for(int c = 1; c <= 4; c++){
      if(!ratio[c]) continue;
      styleH(ratio[c], hexC[c], markC[c]);
      ratio[c]->Draw("E1 same");
      leg->AddEntry(ratio[c], Form("PbPb %s", clsLabel[c]), "lp");
    }
    leg->Draw();

    TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.036);
    la.DrawLatex(0.16, 0.945, "PbPb / pp 5.02 TeV, inclusive calo jets, anti-#it{k}_{T} #it{R} = 0.4");
    la.SetTextSize(0.032);
    la.DrawLatex(0.16, 0.902, "no unfolding, no fake jet subtraction");

    cv->SaveAs(Form("%s/jetsPerZ_caloJets_stitched%s_overlay.pdf", figDir.Data(), tag));
    delete cv;
  }

  // 3. binning diagnostic: bin width vs pT, per class
  {
    TCanvas *cv = new TCanvas("c_binning", "", 700, 800);
    cv->SetLeftMargin(0.15); cv->SetBottomMargin(0.12); cv->SetTopMargin(0.12); cv->SetRightMargin(0.05);
    TH1D *frame = new TH1D("binFrame", "", 100, 0., ptCeiling);
    frame->SetStats(0); frame->SetMinimum(0.); frame->SetMaximum(maxBinWidth * 1.2);
    frame->GetXaxis()->SetTitle("jet #it{p}_{T} [GeV]");
    frame->GetYaxis()->SetTitle("bin width [GeV]");
    frame->GetXaxis()->SetTitleSize(0.045); frame->GetYaxis()->SetTitleSize(0.045);
    frame->GetYaxis()->SetTitleOffset(1.45);
    frame->Draw("axis");
    TLegend *lg = makeLegend(0.20, 0.62, 0.60, 0.82, 0.036);
    for(int c = 1; c <= 4; c++){
      TH1D *w = new TH1D(Form("binW_C%d", c), "", (int)edges[c].size() - 1, &edges[c][0]);
      for(int b = 1; b <= w->GetNbinsX(); b++) w->SetBinContent(b, w->GetXaxis()->GetBinWidth(b));
      styleLine(w, hexC[c], 3);
      w->Draw("hist same");
      lg->AddEntry(w, Form("%s, from %.0f GeV", clsLabel[c], edges[c][0]), "l");
    }
    lg->Draw();
    TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.036);
    la.DrawLatex(0.15, 0.94, Form("binning derived per class for %.0f%% statistical precision", 100.*targetRelErr));
    cv->SaveAs(Form("%s/binning_caloJets%s.pdf", figDir.Data(), tag));
    delete cv;
  }

  // -------------------------------------------------------------- output ---
  TString outPath = Form("%s/jetsPerZ_caloJets_stitched%s.root", outDir.Data(), tag);
  TFile *wf = TFile::Open(outPath, "recreate");
  for(int c = 1; c <= 4; c++){
    ratio[c]->Write(Form("ratio_C%d", c));
    spPb[c]->Write(Form("spectrum_PbPb_C%d", c));
    spPP[c]->Write(Form("spectrum_pp_C%d", c));
  }
  for(int c = 1; c <= 4; c++) if(ppStitch[c]) ppStitch[c]->Write(Form("ppStitch_fine_C%d", c));
  wf->Close();

  printf("\n  figures: %s/{stitchedSpectra,jetsPerZ_caloJets_stitched,binning}_caloJets%s.pdf\n", figDir.Data(), tag);
  printf("  results: %s\n", outPath.Data());
}
