// pp jet trigger efficiency from the MinBias sample.
//
// Reads the histograms added to pp_scan.C in commit bd910fbd and builds the
// turn-on curves two independent ways:
//
//   ABSOLUTE   eff(JetX)      = h_leadJetPt_jetX / h_leadJetPt_all
//   BOOTSTRAP  eff(JetB|JetA) = h_leadJetPt_jetA_and_jetB / h_leadJetPt_jetA
//
// The absolute method needs an unbiased sample, so it is only valid when the
// scan ran with every applyJetNNTrigger = false. That is checked here against
// h_nEventsNoJetTrigSel rather than assumed -- those flags are hard
// `continue`s, so with one enabled the denominator silently becomes the
// triggered sample and every efficiency comes out flat at 1, which looks
// entirely plausible on a plot.
//
// The bootstrap method does not need an unbiased sample and carries the
// high-pT turn-ons where MinBias runs out of statistics, but it is only valid
// on the plateau of the reference trigger JetA. Where both methods are
// available they should agree; that comparison is the point of the right
// panel.
//
// PRESCALES. These HLT paths are prescaled. A prescale draw is random and
// independent of jet pT, so it scales the efficiency by a constant 1/prescale
// without distorting the turn-on SHAPE -- the plateau sits at 1/prescale, not
// at 1. The fit therefore carries a free plateau parameter p0 and the curves
// are divided by it before plotting. p0 is reported next to the mean prescale
// from h_prescale_jetNN so the two can be checked against each other; if
// p0 * <prescale> is not close to 1, something other than prescaling is
// suppressing the plateau and the normalisation should not be trusted.
//
// ERRORS are binomial (TH1::Divide option "B"), which is correct because each
// numerator is a strict subset of its denominator. That is exact for unit
// event weights; if the scan ever runs weighted, they become approximate.
//
// Usage: root -l -b -q 'plotJetTriggerEfficiency_pp.C'
//        root -l -b -q 'plotJetTriggerEfficiency_pp.C("someOther.root")'
// Run from: src/plots/triggerEfficiency/

const char *defaultInput =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/pp/pp_MinBias_mu12_pTmu-15to999_tight_deltaR-40_jetTrkMaxFilter_WDecayFilter_2026-7-6.root";

const char *outDir  = "../../../figures/triggerEfficiency/";
const char *outName = "pp_jetTriggerEfficiency.pdf";

// Set true to build the curves against raw jet pT instead of JEC-corrected.
// Corrected is the right default: the efficiency is applied to the corrected
// analysis spectrum. Raw is a cross-check only.
const bool useRawPt = false;

// Rebin factor on the native 5 GeV binning. 2 -> 10 GeV bins, which is about
// the finest MinBias supports across the full turn-on range.
const int rebinFactor = 2;

const int   NTrig = 6;
const int   trigThresh[NTrig] = {  15,   30,   40,   60,   80,  100 };
const char *trigName[NTrig]   = {"jet15","jet30","jet40","jet60","jet80","jet100"};
const char *trigLabel[NTrig]  = {"Jet15","Jet30","Jet40","Jet60","Jet80","Jet100"};
const char *trigHex[NTrig]    = {"#999999","#0072B2","#009E73","#E69F00","#D55E00","#CC79A7"};
const int   trigMark[NTrig]   = { 24, 20, 21, 33, 34, 29 };

const double plotPtMin = 0., plotPtMax = 260.;

// Turn-on shape: plateau * 0.5 * (1 + erf((x - mu)/(sqrt(2) sigma))).
static double turnOn(double *x, double *p)
{
  return p[0] * 0.5 * (1. + TMath::Erf((x[0] - p[1]) / (TMath::Sqrt(2.) * p[2])));
}

// pT at which the NORMALISED fit first reaches target, scanned finely so the
// answer does not depend on the histogram binning.
static double crossing(TF1 *f, double target, double lo, double hi)
{
  const int N = 20000;
  double plateau = f->GetParameter(0);
  if(plateau <= 0.) return -1.;
  for(int i = 0; i <= N; i++){
    double x = lo + (hi - lo) * i / N;
    if(f->Eval(x) / plateau >= target) return x;
  }
  return -1.;
}

static TH1D* getRebinned(TFile *f, const char *name, const char *newName)
{
  TH1D *h = nullptr;
  f->GetObject(name, h);
  if(!h) return nullptr;
  TH1D *c = (TH1D*) h->Clone(newName);
  c->SetDirectory(nullptr);
  if(rebinFactor > 1) c->Rebin(rebinFactor);
  return c;
}

void plotJetTriggerEfficiency_pp(const char *inputFile = defaultInput)
{
  gStyle->SetOptStat(0);
  gSystem->mkdir(outDir, kTRUE);

  TFile *f = TFile::Open(inputFile);
  if(!f || f->IsZombie()){ printf("ERROR: cannot open %s\n", inputFile); return; }
  printf("Input: %s\n", inputFile);

  const char *ptTag = useRawPt ? "h_leadRawJetPt" : "h_leadJetPt";

  // ---- is the sample actually unbiased? ----
  bool absoluteUsable = true;
  TH1D *hBias = nullptr;
  f->GetObject("h_nEventsNoJetTrigSel", hBias);
  if(!hBias){
    printf("WARNING: h_nEventsNoJetTrigSel missing -- scan predates the trigger-efficiency\n"
           "         histograms. Cannot confirm the sample is unbiased; absolute method\n"
           "         will still be drawn but treat it with suspicion.\n");
  } else {
    double clean = hBias->GetBinContent(1), biased = hBias->GetBinContent(2);
    printf("Bias check: %.0f events with no jet-trigger selection, %.0f with one applied\n", clean, biased);
    if(biased > 0.){
      absoluteUsable = false;
      printf("*** ABSOLUTE METHOD DISABLED: the scan applied a jet-trigger selection, so\n"
             "    h_%s_all is the triggered sample, not an unbiased one. Every\n"
             "    absolute efficiency would come out ~1 by construction. Rerun with all\n"
             "    applyJetNNTrigger = false. Bootstrap method is unaffected and still shown.\n", ptTag);
    }
  }

  TH1D *hDen = getRebinned(f, Form("%s_all", ptTag), "hDen");
  if(!hDen && absoluteUsable){
    printf("ERROR: %s_all not found -- scan predates these histograms.\n", ptTag);
    absoluteUsable = false;
  }

  // ---------------- absolute efficiencies ----------------
  TH1D *effAbs[NTrig]; TF1 *fitAbs[NTrig];
  double plateau[NTrig], mu[NTrig], sigma[NTrig], p95[NTrig], p99[NTrig], meanPrescale[NTrig];
  for(int t = 0; t < NTrig; t++){
    effAbs[t] = nullptr; fitAbs[t] = nullptr;
    plateau[t] = mu[t] = sigma[t] = p95[t] = p99[t] = -1.; meanPrescale[t] = -1.;

    TH1D *hPre = nullptr;
    f->GetObject(Form("h_prescale_%s", trigName[t]), hPre);
    if(hPre && hPre->GetEntries() > 0) meanPrescale[t] = hPre->GetMean();

    if(!absoluteUsable || !hDen) continue;
    TH1D *hNum = getRebinned(f, Form("%s_%s", ptTag, trigName[t]), Form("num_%d", t));
    if(!hNum){ printf("WARNING: %s_%s missing, skipping\n", ptTag, trigName[t]); continue; }
    if(hNum->GetEntries() < 50){ printf("WARNING: %s has only %.0f entries, skipping\n", trigLabel[t], hNum->GetEntries()); continue; }

    effAbs[t] = (TH1D*) hNum->Clone(Form("effAbs_%d", t));
    effAbs[t]->SetDirectory(nullptr);
    effAbs[t]->Divide(hNum, hDen, 1., 1., "B");

    double lo = TMath::Max(5., 0.35 * trigThresh[t]);
    double hi = TMath::Min(plotPtMax, 4.0 * trigThresh[t] + 60.);
    fitAbs[t] = new TF1(Form("fitAbs_%d", t), turnOn, lo, hi, 3);
    fitAbs[t]->SetParameters(effAbs[t]->GetMaximum(), (double)trigThresh[t], 0.15 * trigThresh[t]);
    fitAbs[t]->SetParLimits(0, 1e-6, 1.5);
    fitAbs[t]->SetParLimits(1, 0.2 * trigThresh[t], 3.0 * trigThresh[t]);
    fitAbs[t]->SetParLimits(2, 0.5, 2.0 * trigThresh[t]);
    effAbs[t]->Fit(fitAbs[t], "RQ0");

    plateau[t] = fitAbs[t]->GetParameter(0);
    mu[t]      = fitAbs[t]->GetParameter(1);
    sigma[t]   = fitAbs[t]->GetParameter(2);
    p95[t]     = crossing(fitAbs[t], 0.95, lo, hi);
    p99[t]     = crossing(fitAbs[t], 0.99, lo, hi);
  }

  // ---------------- bootstrap efficiencies ----------------
  // pair t: reference = trigger t, measured = trigger t+1
  TH1D *effBoot[NTrig]; TF1 *fitBoot[NTrig];
  double p99Boot[NTrig];
  for(int t = 0; t < NTrig; t++){ effBoot[t] = nullptr; fitBoot[t] = nullptr; p99Boot[t] = -1.; }

  for(int t = 0; t + 1 < NTrig; t++){
    TH1D *hRef  = getRebinned(f, Form("%s_%s", ptTag, trigName[t]), Form("bref_%d", t));
    TH1D *hBoth = getRebinned(f, Form("%s_%s_and_%s", ptTag, trigName[t], trigName[t+1]), Form("bboth_%d", t));
    if(!hRef || !hBoth){ printf("WARNING: bootstrap pair %s&&%s missing, skipping\n", trigName[t], trigName[t+1]); continue; }
    if(hBoth->GetEntries() < 50){ printf("WARNING: bootstrap %s|%s has only %.0f entries, skipping\n", trigLabel[t+1], trigLabel[t], hBoth->GetEntries()); continue; }

    int m = t + 1;   // the trigger being measured
    effBoot[m] = (TH1D*) hBoth->Clone(Form("effBoot_%d", m));
    effBoot[m]->SetDirectory(nullptr);
    effBoot[m]->Divide(hBoth, hRef, 1., 1., "B");

    double lo = TMath::Max(5., 0.35 * trigThresh[m]);
    double hi = TMath::Min(plotPtMax, 4.0 * trigThresh[m] + 60.);
    fitBoot[m] = new TF1(Form("fitBoot_%d", m), turnOn, lo, hi, 3);
    fitBoot[m]->SetParameters(effBoot[m]->GetMaximum(), (double)trigThresh[m], 0.15 * trigThresh[m]);
    fitBoot[m]->SetParLimits(0, 1e-6, 1.5);
    fitBoot[m]->SetParLimits(1, 0.2 * trigThresh[m], 3.0 * trigThresh[m]);
    fitBoot[m]->SetParLimits(2, 0.5, 2.0 * trigThresh[m]);
    effBoot[m]->Fit(fitBoot[m], "RQ0");
    p99Boot[m] = crossing(fitBoot[m], 0.99, lo, hi);
  }

  // ---------------- report ----------------
  printf("\n=== ABSOLUTE method (%s pT) ===\n", useRawPt ? "raw" : "JEC-corrected");
  printf("%-8s %10s %10s %10s %10s %10s %12s %12s\n",
         "trigger","plateau","<presc>","plat*presc","50%% [GeV]","width","95%% [GeV]","99%% [GeV]");
  for(int t = 0; t < NTrig; t++){
    if(!fitAbs[t]){ printf("%-8s %10s\n", trigLabel[t], "--"); continue; }
    double pp = (meanPrescale[t] > 0.) ? plateau[t]*meanPrescale[t] : -1.;
    printf("%-8s %10.4f %10.2f %10.3f %10.1f %10.1f %12.1f %12.1f\n",
           trigLabel[t], plateau[t], meanPrescale[t], pp, mu[t], sigma[t], p95[t], p99[t]);
  }

  printf("\n=== BOOTSTRAP method, eff(measured | reference) ===\n");
  printf("%-18s %10s %10s %12s   %s\n","pair","plateau","50%% [GeV]","99%% [GeV]","99%% abs - 99%% boot");
  for(int m = 1; m < NTrig; m++){
    if(!fitBoot[m]) continue;
    double d = (p99[m] > 0. && p99Boot[m] > 0.) ? p99[m]-p99Boot[m] : -999.;
    printf("%-18s %10.4f %10.1f %12.1f   %s\n",
           Form("%s | %s", trigLabel[m], trigLabel[m-1]),
           fitBoot[m]->GetParameter(0), fitBoot[m]->GetParameter(1), p99Boot[m],
           d>-900. ? Form("%+.1f GeV", d) : "n/a");
  }
  printf("\nStitch guidance: use each trigger only above its 99%% point.\n");

  // ---------------- draw ----------------
  int col[NTrig];
  for(int t = 0; t < NTrig; t++) col[t] = TColor::GetColor(trigHex[t]);

  TCanvas *c = new TCanvas("cTrigEff", "", 1350, 660);
  c->Divide(2, 1, 0.001, 0.001);

  for(int pad = 0; pad < 2; pad++){
    c->cd(pad + 1);
    gPad->SetLeftMargin(0.13); gPad->SetRightMargin(0.04);
    gPad->SetTopMargin(0.13);  gPad->SetBottomMargin(0.14);

    TH1F *fr = gPad->DrawFrame(plotPtMin, 0., plotPtMax, 1.25);
    fr->GetXaxis()->SetTitle(useRawPt ? "leading jet raw p_{T} [GeV]" : "leading jet p_{T} [GeV]");
    fr->GetYaxis()->SetTitle("trigger efficiency (plateau normalized to 1)");
    fr->GetXaxis()->SetTitleSize(0.047); fr->GetXaxis()->SetLabelSize(0.041);
    fr->GetYaxis()->SetTitleSize(0.047); fr->GetYaxis()->SetLabelSize(0.041);
    fr->GetYaxis()->SetTitleOffset(1.30);

    for(double yl : {0.95, 0.99, 1.00}){
      TLine *l = new TLine(plotPtMin, yl, plotPtMax, yl);
      l->SetLineStyle(yl==1.00 ? 2 : 3); l->SetLineColor(yl==1.00 ? kGray+2 : kGray+1);
      l->Draw();
    }

    TLegend *leg = new TLegend(0.44, 0.13, 0.965, 0.45);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.034);

    bool drewAny = false;
    for(int t = 0; t < NTrig; t++){
      TH1D *e = (pad == 0) ? effAbs[t]  : effBoot[t];
      TF1  *g = (pad == 0) ? fitAbs[t]  : fitBoot[t];
      if(!e || !g) continue;
      double pl = g->GetParameter(0);
      if(pl <= 0.) continue;

      // normalize a copy so the drawn curve has its plateau at 1
      TH1D *en = (TH1D*) e->Clone(Form("norm_%d_%d", pad, t));
      en->SetDirectory(nullptr);
      en->Scale(1./pl);
      en->SetLineColor(col[t]); en->SetMarkerColor(col[t]);
      en->SetMarkerStyle(trigMark[t]); en->SetMarkerSize(0.9); en->SetLineWidth(2);
      en->GetXaxis()->SetRangeUser(plotPtMin, plotPtMax);
      en->Draw("ep same");

      TF1 *gn = new TF1(Form("gn_%d_%d", pad, t), turnOn, g->GetXmin(), g->GetXmax(), 3);
      gn->SetParameters(1.0, g->GetParameter(1), g->GetParameter(2));
      gn->SetLineColor(col[t]); gn->SetLineWidth(2); gn->SetLineStyle(1);
      gn->Draw("l same");

      double p99v = (pad == 0) ? p99[t] : p99Boot[t];
      leg->AddEntry(en, Form("%s  (99%% @ %.0f GeV)", trigLabel[t], p99v), "lp");
      drewAny = true;
    }
    if(drewAny) leg->Draw();

    TLatex lt; lt.SetNDC(); lt.SetTextSize(0.044); lt.SetTextFont(62);
    lt.DrawLatex(0.16, 0.915, pad == 0 ? "absolute (MinBias denominator)" : "bootstrap (lower-threshold reference)");

    if(pad == 0 && !absoluteUsable){
      TLatex wr; wr.SetNDC(); wr.SetTextSize(0.038); wr.SetTextColor(kRed+1);
      wr.DrawLatex(0.17, 0.60, "#splitline{DISABLED: scan applied a jet-trigger}{selection, denominator is biased}");
    }
    if(!drewAny && pad == 1){
      TLatex wr; wr.SetNDC(); wr.SetTextSize(0.038); wr.SetTextColor(kRed+1);
      wr.DrawLatex(0.20, 0.60, "no bootstrap pairs available");
    }
  }

  c->cd(0);
  TLatex lat; lat.SetNDC(); lat.SetTextSize(0.028);
  lat.DrawLatex(0.03, 0.982, "pp jet trigger efficiency, MinBias sample");

  TString out = TString(outDir) + outName;
  c->SaveAs(out);
  printf("\nSaved %s\n", out.Data());
}
