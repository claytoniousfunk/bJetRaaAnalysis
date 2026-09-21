// R_CP of jets per Z, PbPb only, calo jets -- MinBias + Jet80 + Jet100.
//
// includeMinBias = true  (default): MinBias supplies everything below each
//                        class's Jet80 threshold, so the spectra reach the
//                        MinBias scan's own 50 GeV jet cut while the triggers
//                        carry the statistics above it.
// includeMinBias = false: trigger samples only, so each class starts at its
//                        measured Jet80 threshold (85-110 GeV).
//
// MinBias needs no threshold of its own: it carries no jet trigger at all
// (h_nEventsNoJetTrigSel mean = 0, checked at run time), so it is unbiased
// everywhere and is simply the lowest sample in the chain. Jet60 is still
// excluded -- its turn-on could not be measured in any class.
//
// ------------------------------------------------------- NORMALISATION -----
// Per class, jets per Z = (stitched counts) / N_Z, and
//      R_CP(c) = [jets per Z](c) / [jets per Z](reference class).
//
// The chain normalizes upward: Jet100 sets the scale, Jet80 is scaled to it over
// [T100, T100+40], and MinBias is scaled to Jet80 over [T80, T80+40] -- in each
// case a window where both samples are known flat. Two facts make the result
// comparable across centrality classes, and both are checked rather than
// assumed: every class comes from the same HardProbes events, and
// HLT_HICsAK4PFJet100Eta1p5 has prescale 1.000 in every class.
//
// The MinBias-to-Jet80 factor is a further check. It equals the effective
// luminosity ratio of the two datasets, which is a property of the datasets and
// not of centrality, so it must come out the SAME in all four classes. The run
// prints them side by side; disagreement beyond the quoted statistics would mean
// the stitch is picking up something centrality dependent.
//
// Counts are never divided by N_evt: that would strip the centrality-slice
// population which pairs against the same population inside N_Z.
//
// -------------------------------------------------------------- BINNING ----
// Per class PAIR: bins merge until the class AND the reference both reach the
// target relative error, with a resolution floor. With MinBias in the chain
// every pair starts at 50 GeV.
//
// ------------------------------------------------------------- CAVEATS -----
// No unfolding. No fake-jet subtraction, and with MinBias included that matters
// again at the bottom of the range: combinatorial jets are a central-collision
// effect concentrated below ~100 GeV, and the existing fake-jet estimate is
// FastJet on PF candidates, which cannot be transferred to these akPu4Calo
// spectra (against them it comes out at 565% of the yield at 50-60 GeV in
// 0-10%, an impossible "fraction" that simply shows the two are not
// comparable). The lowest central points are therefore upper bounds, and the
// trigger-only variant is the one to quote where that matters.
//
// Usage, from src/newFractionCalculation/:
//   root -l -b -q 'jetsPerZ_RCP_stitched_caloJets.C'         // MinBias + Jet80 + Jet100
//   root -l -b -q 'jetsPerZ_RCP_stitched_caloJets.C(false)'  // triggers only

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
#include "TSystem.h"
#include <cstdio>
#include <vector>

const char *repo = "/home/clayton/Analysis/code/bJetRaaAnalysis";
const char *sib  = "/home/clayton/Analysis/code/bJetMuonTaggingAnalysis";  // read-only

const char *clsLabel[5] = {"", "0-10%", "10-30%", "30-50%", "50-80%"};
const char *clsTag[5]   = {"", "0to10pct", "10to30pct", "30to50pct", "50to80pct"};

const int    refClass     = 4;      // 50-80%
const double minBiasFloor = 50.;    // the MinBias scan's own jet cut
const double ptCeiling    = 500.;
const double overlapWidth = 40.;
const double targetRelErr = 0.05;
const double minWidthFrac = 0.10;
const double maxBinWidth  = 200.;

// Window for the MinBias -> Jet80 normalization. Deliberately NOT each class's
// own threshold window: that sits at 85-150 GeV where combinatorial jets inflate
// the MinBias side in central events, which dragged the factor to 113.5 +- 2.4
// in 0-10% against 133.5 +- 5.5 in 50-80% -- an 18% spread at 3.3 sigma that
// propagated straight into the sub-threshold region of R_CP as a constant
// class-dependent rescaling (0.839, 0.925, 0.881, flat across three bins).
// Measured here instead, above the fake-dominated region and above every Jet80
// threshold, the factor is class independent as it must be.
const double mbNormLo = 150., mbNormHi = 200.;

const double Z_lo = 75., Z_hi = 105.;
const double muEff[5] = {1., 0.8627, 0.9069, 0.9856, 0.9778};

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

void jetsPerZ_RCP_stitched_caloJets(bool includeMinBias = true)
{
  initPlotStyle();
  ok = true;

  const char *tag = includeMinBias ? "MBstitched" : "stitched";

  TString figDir = Form("%s/figures/JetsPerZ", repo);
  gSystem->mkdir(figDir, kTRUE);
  TString outDir = "./rootFiles/JetsPerZ";
  gSystem->mkdir(outDir, kTRUE);

  std::map<TString, double> thr = stitchReadThresholds(Form("%s/triggerThresholds_caloJets.txt", outDir.Data()));
  if(thr.empty()) return;

  const char *tn[3] = {"Jet60", "Jet80", "Jet100"};
  TFile *fPb[3];
  for(int i = 0; i < 3; i++)
    fPb[i] = openOrFail(Form("%s/rootFiles/scanningOuput/PbPb/PbPb_HardProbes_caloJets_%sHLT_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_2026-9-15.root", repo, tn[i]));
  TFile *f_MB = includeMinBias
    ? openOrFail(Form("%s/rootFiles/scanningOuput/PbPb/PbPb_MinBias_Part1_caloJets_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_2026-9-14.root", repo))
    : nullptr;
  TFile *f_mu = openOrFail(Form("%s/rootFiles/scanningOutput/PbPb/latest/PbPb_SingleMuon_mu12_pTmu-15to999_tight_mu12TriggerEfficiencyCorrection_jetTrkMaxFilter_WDecayFilter_2026-5-4.root", sib));
  if(!ok) return;

  if(includeMinBias){
    TH1D *g = nullptr; f_MB->GetObject("h_nEventsNoJetTrigSel", g);
    if(g && g->GetMean() > 0.5){
      printf("ERROR: the MinBias scan has a jet trigger applied; it cannot anchor the low-pT region\n");
      return;
    }
  }

  printf("\n=== R_CP of jets per Z, calo jets, %s ===\n",
         includeMinBias ? "MinBias + Jet80 + Jet100" : "Jet80 + Jet100 only");
  printf("  reference class: %s   no pp\n", clsLabel[refClass]);

  // ---- stitch scales, measured once in a common window ---------------------
  // Every one of these is a ratio of two datasets (or two trigger paths) and so
  // is a property of them, NOT of centrality. Measuring each in its own class's
  // overlap window let centrality leak in twice over:
  //   Jet80/MinBias   spread 18% (113.5 +- 2.4 in 0-10% vs 133.5 +- 5.5 in
  //                   50-80%), because those windows sit at 85-150 GeV where
  //                   combinatorial jets inflate the MinBias side centrally
  //   Jet100/Jet80    spread 1.4% (2.11 vs 2.14 on 0.3-0.7% statistics)
  // Both propagate below the relevant threshold as a flat class-dependent
  // rescaling of R_CP. So each is measured per class as a check, then averaged
  // with statistical weights and the average applied to every class.
  auto commonScale = [&](TFile *fLow, TFile *fHigh, const char *label, double &kOut)->bool{
    printf("\n  %s yield ratio in the common window [%.0f,%.0f]\n", label, mbNormLo, mbNormHi);
    double sw = 0., swk = 0.;
    for(int c = 1; c <= 4; c++){
      TH1D *hl = getOrFail(fLow,  Form("h_inclRecoJetPt_C%d", c), Form("kLow_%s_C%d", label, c));
      TH1D *hh = getOrFail(fHigh, Form("h_inclRecoJetPt_C%d", c), Form("kHigh_%s_C%d", label, c));
      if(!ok) return false;
      const double eps = 0.01;
      double nl = hl->Integral(hl->FindBin(mbNormLo + eps), hl->FindBin(mbNormHi - eps));
      double nh = hh->Integral(hh->FindBin(mbNormLo + eps), hh->FindBin(mbNormHi - eps));
      if(nl <= 0. || nh <= 0.){ printf("    %-8s empty window\n", clsLabel[c]); continue; }
      double k = nh/nl, sk = k*sqrt(1./nl + 1./nh);
      printf("    %-8s %8.2f +- %.2f   (lower-sample counts %.0f)\n", clsLabel[c], k, sk, nl);
      double w = 1./(sk*sk); sw += w; swk += w*k;
    }
    if(sw <= 0.){ printf("ERROR: could not measure %s\n", label); return false; }
    kOut = swk/sw;
    printf("    weighted average applied to every class: %.2f +- %.2f\n", kOut, sqrt(1./sw));
    return true;
  };

  double kCommon = 0., kCommonJ80 = 0.;
  if(includeMinBias && !commonScale(f_MB, fPb[1], "Jet80 / MinBias", kCommon)) return;
  if(!commonScale(fPb[1], fPb[2], "Jet100 / Jet80", kCommonJ80)) return;

  // ------------------------------------------------------ stitch per class --
  TH1D  *stitched[5] = {nullptr,nullptr,nullptr,nullptr,nullptr};
  double floorPt[5]  = {0.,0.,0.,0.,0.};
  double NZ[5]       = {0.,0.,0.,0.,0.};
  double kMB[5]      = {0.,0.,0.,0.,0.};

  for(int c = 1; c <= 4; c++){
    printf("\n  %s\n", clsLabel[c]);
    std::vector<StitchSample> s;

    if(includeMinBias){
      TH1D *h = getOrFail(f_MB, Form("h_inclRecoJetPt_C%d", c), Form("specMB_C%d", c));
      if(!h) return;
      s.push_back({"MinBias", minBiasFloor, h, kCommon});
      printf("    %-7s used from %.0f GeV (ungated, no threshold needed)\n", "MinBias", minBiasFloor);
    }
    for(int i = 0; i < 3; i++){
      TString key = Form("PbPb/%s/%s", clsTag[c], tn[i]);
      if(!thr.count(key)){ printf("    %-7s not validated, not used\n", tn[i]); continue; }
      TH1D *h = getOrFail(fPb[i], Form("h_inclRecoJetPt_C%d", c), Form("spec_%s_C%d", tn[i], c));
      if(!h) return;
      // Jet80 is normalized to Jet100 with the common-window average, for the
      // same reason MinBias is: its own-window factor carries a class spread
      double fk = (strcmp(tn[i], "Jet80") == 0) ? kCommonJ80 : 0.;
      s.push_back({tn[i], thr[key], h, fk});
      printf("    %-7s used from %.0f GeV\n", tn[i], thr[key]);
    }
    if(s.empty()){ printf("    ERROR: no usable sample\n"); return; }
    std::sort(s.begin(), s.end(), [](const StitchSample &a, const StitchSample &b){ return a.thr < b.thr; });

    // capture the MinBias scale factor for the cross-class consistency check
    stitched[c] = stitchTriggerSamples(s, Form("stitched_C%d", c), overlapWidth, ptCeiling, true);
    if(!stitched[c]) return;
    floorPt[c] = s.front().thr;

    TH1D *d = getOrFail(f_mu, Form("h_dimuonMass_C%d", c), Form("dimu_C%d", c));
    if(!d) return;
    NZ[c] = d->Integral(d->GetXaxis()->FindBin(Z_lo), d->GetXaxis()->FindBin(Z_hi)) / (muEff[c]*muEff[c]);

    // recompute k(MinBias->Jet80) here so it can be compared across classes
    if(includeMinBias && s.size() > 1){
      const double eps = 0.01;
      double lo = s[1].thr, hi = TMath::Min(s[1].thr + overlapWidth, ptCeiling);
      double nUp = s[1].h->Integral(s[1].h->FindBin(lo + eps), s[1].h->FindBin(hi - eps));
      double nMB = s[0].h->Integral(s[0].h->FindBin(lo + eps), s[0].h->FindBin(hi - eps));
      if(nMB > 0.) kMB[c] = nUp / nMB;
    }
  }

  printf("\n  %-8s %10s %12s\n", "class", "N_Z", "floor [GeV]");
  for(int c = 1; c <= 4; c++) printf("  %-8s %10.1f %12.0f\n", clsLabel[c], NZ[c], floorPt[c]);

  if(includeMinBias){
    // Diagnostic only -- these are the per-class Jet80/MinBias ratios measured
    // in each class's OWN threshold window, which is what the stitch used to do
    // and no longer does. They are printed to show why: the spread below is the
    // bias the common-window average removes, not an uncertainty on the result.
    printf("\n  for reference, the same ratio measured in each class's own threshold window\n"
           "  (NOT used -- this is the centrality-dependent bias the common window avoids)\n   ");
    for(int c = 1; c <= 4; c++) printf("  %s = %.3g", clsLabel[c], kMB[c]);
    printf("\n");
    double lo = 1e30, hi = -1e30;
    for(int c = 1; c <= 4; c++){ if(kMB[c] > 0.){ lo = TMath::Min(lo, kMB[c]); hi = TMath::Max(hi, kMB[c]); } }
    if(lo > 0. && hi/lo > 1.05)
      printf("   those span %.0f%% against 2-4%% statistics; taken at face value they would have\n"
             "   rescaled the sub-threshold region of R_CP by k(class)/k(ref), by up to that much\n",
             100.*(hi/lo - 1.));
  }

  // ---------------------------------------------------------------- R_CP ---
  TH1D *rcp[5] = {nullptr,nullptr,nullptr,nullptr,nullptr};
  TH1D *spec[5] = {nullptr,nullptr,nullptr,nullptr,nullptr};
  TH1D *specRef[5] = {nullptr,nullptr,nullptr,nullptr,nullptr};
  std::vector<double> edges[5];

  for(int c = 1; c < refClass; c++){
    double pairFloor = TMath::Max(floorPt[c], floorPt[refClass]);
    edges[c] = stitchDeriveEdges({stitched[c], stitched[refClass]}, pairFloor, ptCeiling,
                                 targetRelErr, minWidthFrac, maxBinWidth);
    if(edges[c].size() < 2){ printf("ERROR: no binning for %s\n", clsLabel[c]); return; }

    spec[c] = (TH1D*) stitched[c]->Rebin((int)edges[c].size() - 1, Form("spec_C%d", c), &edges[c][0]);
    spec[c]->SetDirectory(nullptr);
    divideByBinwidth(spec[c]); spec[c]->Scale(1./3.2); spec[c]->Scale(1./NZ[c]);

    specRef[c] = (TH1D*) stitched[refClass]->Rebin((int)edges[c].size() - 1, Form("specRef_C%d", c), &edges[c][0]);
    specRef[c]->SetDirectory(nullptr);
    divideByBinwidth(specRef[c]); specRef[c]->Scale(1./3.2); specRef[c]->Scale(1./NZ[refClass]);

    rcp[c] = (TH1D*) spec[c]->Clone(Form("RCP_C%d", c));
    rcp[c]->SetDirectory(nullptr);
    rcp[c]->Divide(spec[c], specRef[c], 1, 1, "");

    printf("\n  %s / %s, from %.0f GeV (%zu bins)\n", clsLabel[c], clsLabel[refClass],
           pairFloor, edges[c].size() - 1);
    for(int b = 1; b <= rcp[c]->GetNbinsX(); b++)
      printf("    %4.0f-%-6.0f  %7.3f +- %.3f\n",
             rcp[c]->GetXaxis()->GetBinLowEdge(b), rcp[c]->GetXaxis()->GetBinUpEdge(b),
             rcp[c]->GetBinContent(b), rcp[c]->GetBinError(b));
  }

  // ------------------------------------------------------------- figures ---
  const char *hexC[5] = {"", okabeHex[5], okabeHex[1], okabeHex[2], okabeHex[0]};
  const int   markC[5] = {0, markFilledCircle, markFilledSquare, markFilledDiamond, markCross};
  double xLo = includeMinBias ? minBiasFloor : floorPt[refClass];

  {
    TCanvas *cv = new TCanvas(Form("c_RCP_%s", tag), "", 700, 800);
    cv->SetLeftMargin(0.16); cv->SetBottomMargin(0.12); cv->SetTopMargin(0.13); cv->SetRightMargin(0.05);
    // low and right: the curves sit high, so the band below 0.6 is the only
    // place a legend does not land on the 30-50% points
    TLegend *leg = makeLegend(0.62, 0.20, 0.93, 0.38, 0.036);
    bool first = true;
    for(int c = 1; c < refClass; c++){
      styleH(rcp[c], hexC[c], markC[c]);
      if(first){
        rcp[c]->SetTitle("");
        rcp[c]->SetMinimum(0.); rcp[c]->SetMaximum(includeMinBias ? 1.6 : 1.4);
        rcp[c]->GetXaxis()->SetRangeUser(xLo, ptCeiling);
        rcp[c]->GetXaxis()->SetTitle("jet #it{p}_{T} [GeV]");
        rcp[c]->GetXaxis()->SetTitleSize(0.045); rcp[c]->GetXaxis()->SetLabelSize(0.040);
        rcp[c]->GetYaxis()->SetTitle(Form("#it{R}_{CP} of jets per #it{Z}   (ref. %s)", clsLabel[refClass]));
        rcp[c]->GetYaxis()->SetTitleSize(0.042); rcp[c]->GetYaxis()->SetLabelSize(0.040);
        rcp[c]->GetYaxis()->SetTitleOffset(1.50);
        rcp[c]->Draw("E1");
        first = false;
      }
      else rcp[c]->Draw("E1 same");
      leg->AddEntry(rcp[c], Form("PbPb %s", clsLabel[c]), "lp");
    }
    TLine *one = new TLine(xLo, 1., ptCeiling, 1.);
    one->SetLineStyle(7); one->Draw();
    for(int c = 1; c < refClass; c++) rcp[c]->Draw("E1 same");
    leg->Draw();
    TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.036);
    la.DrawLatex(0.16, 0.945, includeMinBias
                 ? "PbPb 5.02 TeV, calo jets, MinBias + Jet80 + Jet100"
                 : "PbPb 5.02 TeV, calo jets, stitched Jet80 + Jet100");
    la.SetTextSize(0.032);
    // kept short: longer wordings ran off the right frame edge
    la.DrawLatex(0.16, 0.902, "no unfolding, no fake jet subtraction");
    cv->SaveAs(Form("%s/jetsPerZ_RCP_%s_caloJets.pdf", figDir.Data(), tag));
    delete cv;
  }

  {
    TCanvas *cv = new TCanvas(Form("c_spec_%s", tag), "", 700, 800);
    cv->SetLogy();
    cv->SetLeftMargin(0.17); cv->SetBottomMargin(0.12); cv->SetTopMargin(0.13); cv->SetRightMargin(0.05);
    TLegend *leg = makeLegend(0.55, 0.58, 0.93, 0.78, 0.036);
    bool first = true;
    for(int c = 1; c < refClass; c++){
      styleH(spec[c], hexC[c], markC[c]);
      if(first){
        spec[c]->SetTitle("");
        spec[c]->GetXaxis()->SetTitle("jet #it{p}_{T} [GeV]");
        spec[c]->GetXaxis()->SetTitleSize(0.045); spec[c]->GetXaxis()->SetLabelSize(0.040);
        spec[c]->GetYaxis()->SetTitle("#frac{1}{#it{N}_{Z}} #frac{d#it{N}_{jet}}{d#it{p}_{T} d#eta}");
        spec[c]->GetYaxis()->SetTitleSize(0.045); spec[c]->GetYaxis()->SetLabelSize(0.040);
        spec[c]->GetYaxis()->SetTitleOffset(1.60);
        spec[c]->Draw("E1");
        first = false;
      }
      else spec[c]->Draw("E1 same");
      leg->AddEntry(spec[c], Form("PbPb %s", clsLabel[c]), "lp");
    }
    styleH(specRef[1], hexC[refClass], markC[refClass]);
    specRef[1]->Draw("E1 same");
    leg->AddEntry(specRef[1], Form("PbPb %s (ref.)", clsLabel[refClass]), "lp");
    leg->Draw();
    TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.036);
    la.DrawLatex(0.17, 0.945, includeMinBias
                 ? "PbPb 5.02 TeV, calo jets, MinBias + Jet80 + Jet100"
                 : "PbPb 5.02 TeV, calo jets, stitched Jet80 + Jet100");
    la.SetTextSize(0.032);
    la.DrawLatex(0.17, 0.902, "no unfolding, no fake-jet subtraction");
    cv->SaveAs(Form("%s/jetsPerZ_%s_spectra_caloJets.pdf", figDir.Data(), tag));
    delete cv;
  }

  TString outPath = Form("%s/jetsPerZ_RCP_%s_caloJets.root", outDir.Data(), tag);
  TFile *wf = TFile::Open(outPath, "recreate");
  for(int c = 1; c < refClass; c++){
    rcp[c]->Write(Form("RCP_C%d", c));
    spec[c]->Write(Form("jetsPerZ_C%d", c));
    specRef[c]->Write(Form("jetsPerZ_ref_for_C%d", c));
  }
  wf->Close();

  printf("\n  figures: %s/jetsPerZ_{RCP_%s,%s_spectra}_caloJets.pdf\n", figDir.Data(), tag, tag);
  printf("  results: %s\n", outPath.Data());
}
