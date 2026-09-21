// HLT_HIL3Mu12 single-muon trigger efficiency vs muon pT, from the trigger
// scans (src/scanning/{pp,PbPb}/*_trigger_scan.C), with the plateau fit
//
//   f(pT) = a - (a - b) exp(-k pT),   fitted for pT >= fitMin
//
// with 0 <= b <= a enforced: the fit is done in r = b/a, i.e.
//   f(pT) = a [1 - (1 - r) exp(-k pT)],  0 <= r <= 1,
// so the curve can only rise towards its plateau. a, b = a r and k are
// reported, with the error on b propagated from the (a, r) covariance.
//
// Efficiency = muPt_trigOn / muPt_all per centrality class. The scans fill
// muPt_trigOn with weight = HLT prescale and muPt_all unweighted, so a bin with
// any prescaled entries is not a plain binomial ratio. The interval used here is
// Clopper-Pearson on the effective counts (both scaled by sumw/sumw2 of the
// pass histogram), which is exact when every weight is 1.
//
// Usage:
//   root -l -b -q 'plotMuonTriggerEfficiency.C'
// Writes PDFs to figures/triggerEfficiency/ and the fit parameters to
// figures/triggerEfficiency/HLTFitParams_<tag>.h (same variable names as
// headers/fitParameters/HLTFitParams_{pp,PbPb}.h, so it can replace them).

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TGraphAsymmErrors.h"
#include "TEfficiency.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLatex.h"
#include "TFitResult.h"
#include "TSystem.h"
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <vector>
#include <string>
#include "../../../headers/plotting/plotStyle.h"

namespace {

const std::string scanDir = "/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/scanningOutput/";
const std::string outDir  = "/home/clayton/Analysis/code/bJetRaaAnalysis/figures/triggerEfficiency/";

const double fitMin = 15.0, fitMax = 100.0;

// display / fit binning (every edge on a 1 GeV boundary of the scan histogram)
const std::vector<double> edges = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
                                   22,24,26,28,30,35,40,50,60,80,100};

// hiBin boundaries used by the scan that produced each file; C4 of the
// Nov-2024 MinBias scan runs to hiBin 180 (50-90%), later scans stop at 160.
struct Sample {
  std::string tag, file, label, trigLabel;
  int nCent;               // 1 for pp, 5 for PbPb (C0 inclusive + 4 classes)
  const char *periph;      // label for C4
  // true: read the per-muon h_muTrigEff_pass/all (pT vs eta) histograms that
  // PbPb_scan.C / pp_scan.C have filled since 2026-09-18 (prescale-1 events,
  // unweighted, other-muon veto), projected over |eta| < etaMax2D. false: the
  // old *_trigger_scan.C muPt_trigOn/muPt_all pair.
  bool from2D = false;
  const char *c0Label = "";   // label for C0 if not the default 0-90%
};

const double etaMax2D = 2.0;   // in-jet muon acceptance

const char *centLabel[5] = {"0-90%", "0-10%", "10-30%", "30-50%", ""};

// effective-count efficiency graph from a (possibly prescale-weighted) pass
// histogram and an unweighted total histogram
TGraphAsymmErrors* makeEff(TH1D *pass, TH1D *total)
{
  int n = edges.size() - 1;
  TGraphAsymmErrors *g = new TGraphAsymmErrors();
  for (int i = 0; i < n; i++) {
    int lo = total->FindBin(edges[i] + 1e-6), hi = total->FindBin(edges[i+1] - 1e-6);
    double N = 0, sw = 0, sw2 = 0;
    for (int b = lo; b <= hi; b++) {
      N   += total->GetBinContent(b);
      sw  += pass->GetBinContent(b);
      sw2 += pass->GetSumw2N() ? pass->GetSumw2()->At(b) : pass->GetBinContent(b);
    }
    if (N <= 0) continue;
    double eps = std::min(sw / N, 1.0);
    double s   = sw2 > 0 ? sw / sw2 : 1.0;      // effective-count scale
    double Ne  = N * s, ke = std::min(sw * s, Ne);
    double lo68 = TEfficiency::ClopperPearson(Ne, ke, 0.682689, false);
    double hi68 = TEfficiency::ClopperPearson(Ne, ke, 0.682689, true);
    double x = 0.5 * (edges[i] + edges[i+1]), dx = 0.5 * (edges[i+1] - edges[i]);
    int p = g->GetN();
    g->SetPoint(p, x, eps);
    g->SetPointError(p, dx, dx, std::max(eps - lo68, 0.0), std::max(hi68 - eps, 0.0));
  }
  return g;
}

struct FitOut { double p[3], e[3], chi2; int ndf; double nFit; };

FitOut fitEff(TGraphAsymmErrors *g, TF1 *f, TH1D *total)
{
  f->SetParameters(0.95, 0.5, 0.1);
  f->SetParLimits(0, 0.5, 1.0);
  f->SetParLimits(1, 0.0, 1.0);   // r = b/a
  f->SetParLimits(2, 0.0, 2.0);
  TFitResultPtr r = g->Fit(f, "RSQ0");
  FitOut o;
  for (int i = 0; i < 3; i++) { o.p[i] = f->GetParameter(i); o.e[i] = f->GetParError(i); }
  // (a, r) -> (a, b = a r)
  double a = o.p[0], rr = o.p[1];
  double cov_ar = (r.Get() && r->IsValid()) ? r->CovMatrix(0, 1) : 0.;
  o.p[1] = a * rr;
  o.e[1] = std::sqrt(std::max(0., rr * rr * o.e[0] * o.e[0] + a * a * o.e[1] * o.e[1] + 2 * a * rr * cov_ar));
  o.chi2 = f->GetChisquare(); o.ndf = f->GetNDF();
  o.nFit = total->Integral(total->FindBin(fitMin + 1e-6), total->FindBin(fitMax - 1e-6));
  return o;
}

void drawPanel(TGraphAsymmErrors *g, TF1 *f, const FitOut &o, const std::string &line1,
               const std::string &line2, const std::string &out,
               TGraphAsymmErrors *ref = nullptr, const char *refLabel = nullptr,
               const char *label = "")
{
  TCanvas c("c", "c", 700, 600);
  gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13);
  gPad->SetRightMargin(0.04); gPad->SetTopMargin(0.05);

  TH1D frame("frame", "", 100, 0, 100);
  frame.SetStats(0);
  frame.SetMinimum(0); frame.SetMaximum(1.25);
  frame.GetXaxis()->SetTitle("#it{p}_{T}^{#mu} [GeV]");
  frame.GetYaxis()->SetTitle("Trigger Efficiency");
  frame.GetXaxis()->SetTitleSize(0.05); frame.GetYaxis()->SetTitleSize(0.05);
  frame.GetXaxis()->SetLabelSize(0.045); frame.GetYaxis()->SetLabelSize(0.045);
  frame.GetYaxis()->SetTitleOffset(1.25);
  frame.Draw("axis");

  TLine l; l.SetLineStyle(7); l.SetLineColor(kGray + 2);
  l.DrawLine(fitMin, 0, fitMin, 1.05);
  l.SetLineStyle(3); l.DrawLine(0, 1, 100, 1);

  if (ref) {
    ref->SetMarkerStyle(markOpenSquare); ref->SetMarkerSize(1.0);
    int cr = TColor::GetColor(hexMC); ref->SetMarkerColor(cr); ref->SetLineColor(cr);
    ref->Draw("P same");
  }

  // fitted range solid, extrapolation below it dashed
  TF1 *fx = (TF1*)f->Clone("fx");
  fx->SetRange(0, fitMin); fx->SetLineStyle(2); fx->SetLineWidth(2);
  fx->SetLineColor(TColor::GetColor(hexCorrected));
  fx->Draw("same");
  f->SetLineColor(TColor::GetColor(hexCorrected)); f->SetLineWidth(3); f->SetLineStyle(1);
  f->Draw("same");

  g->SetMarkerStyle(markFilledCircle); g->SetMarkerSize(1.0);
  g->SetMarkerColor(kBlack); g->SetLineColor(kBlack);
  g->Draw("P same");

  TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.042);
  t.DrawLatex(0.18, 0.89, line1.c_str());
  t.DrawLatex(0.18, 0.84, line2.c_str());
  t.SetTextSize(0.036);
  t.DrawLatex(0.50, 0.40, "#it{f} = #it{a} #minus (#it{a} #minus #it{b}) e^{#minus#it{k} p_{T}}");
  t.DrawLatex(0.50, 0.34, Form("#it{a} = %.4f #pm %.4f", o.p[0], o.e[0]));
  t.DrawLatex(0.50, 0.29, Form("#it{b} = %.3f #pm %.3f", o.p[1], o.e[1]));
  t.DrawLatex(0.50, 0.24, Form("#it{k} = %.3f #pm %.3f GeV^{#minus1}", o.p[2], o.e[2]));
  t.DrawLatex(0.50, 0.19, Form("#chi^{2}/ndf = %.1f/%d = %.2f", o.chi2, o.ndf,
                               o.ndf > 0 ? o.chi2 / o.ndf : 0.0));

  if (ref) {
    TLegend *leg = makeLegend(0.49, 0.45, 0.90, 0.56, 0.036);
    leg->AddEntry(g, label, "p");
    leg->AddEntry(ref, refLabel, "p");
    leg->Draw();
  }
  c.SaveAs(out.c_str());
  delete fx;
}

} // namespace

// Measure one sample. Returns fit results per class (index = C).
std::vector<FitOut> runSample(const Sample &s, std::vector<TGraphAsymmErrors*> &graphs,
                              const std::vector<TGraphAsymmErrors*> *ref = nullptr,
                              const char *refLabel = nullptr)
{
  std::vector<FitOut> res(s.nCent);
  graphs.assign(s.nCent, nullptr);
  TFile *F = TFile::Open((s.file[0] == '/' ? s.file : scanDir + s.file).c_str());
  if (!F || F->IsZombie()) { printf("cannot open %s\n", s.file.c_str()); return res; }
  printf("\n== %s  (%s)\n", s.tag.c_str(), s.file.c_str());
  for (int c = 0; c < s.nCent; c++) {
    TH1D *pass = nullptr, *tot = nullptr;
    if(s.from2D){
      // pp has no centrality suffix
      TString sfx = s.nCent == 1 ? TString("") : TString::Format("_C%d", c);
      TH2D *P = (TH2D*)F->Get("h_muTrigEff_pass" + sfx);
      TH2D *A = (TH2D*)F->Get("h_muTrigEff_all" + sfx);
      int e1 = A->GetYaxis()->FindBin(-etaMax2D + 1e-6), e2 = A->GetYaxis()->FindBin(etaMax2D - 1e-6);
      pass = P->ProjectionX(Form("pass_%s_C%d", s.tag.c_str(), c), e1, e2);
      tot  = A->ProjectionX(Form("tot_%s_C%d",  s.tag.c_str(), c), e1, e2);
    }
    else{
      pass = (TH1D*)F->Get(Form("muPt_trigOn_C%d", c));
      tot  = (TH1D*)F->Get(Form("muPt_all_C%d", c));
    }
    TGraphAsymmErrors *g = makeEff(pass, tot);
    graphs[c] = g;
    TF1 *f = new TF1(Form("f_%s_C%d", s.tag.c_str(), c), "[0]*(1-(1-[1])*exp(-[2]*x))", fitMin, fitMax);
    res[c] = fitEff(g, f, tot);
    const FitOut &o = res[c];
    std::string cl = s.nCent == 1 ? "pp" : std::string("PbPb ") +
                     (c == 4 ? s.periph : (c == 0 && s.c0Label[0]) ? s.c0Label : centLabel[c]);
    printf("  C%d %-12s a=%.4f+-%.4f  b=%.3f+-%.3f  k=%.4f+-%.4f  chi2/ndf=%.1f/%d  N(pT>%.0f)=%.0f\n",
           c, cl.c_str(), o.p[0], o.e[0], o.p[1], o.e[1], o.p[2], o.e[2], o.chi2, o.ndf, fitMin, o.nFit);
    std::string out = outDir + "triggerEfficiency_mu12_" + s.tag + (s.nCent == 1 ? "" : Form("_C%d", c)) + ".pdf";
    drawPanel(g, f, o, cl + ", " + s.label, s.trigLabel + ", tight ID, |#eta^{#mu}| < 2.0", out,
              ref ? (*ref)[c] : nullptr, refLabel, s.label.c_str());
  }
  return res;
}

void writeHeader(const std::string &tag, const std::vector<FitOut> &r, bool isPbPb)
{
  std::string fn = outDir + "HLTFitParams_" + tag + ".h";
  FILE *h = fopen(fn.c_str(), "w");
  fprintf(h, "// written by src/plots/triggerEfficiency/plotMuonTriggerEfficiency.C (%s)\n", tag.c_str());
  fprintf(h, "// f(pT) = p0 - (p0 - p1) exp(-p2 pT), fit for %.0f < pT < %.0f GeV\n\n", fitMin, fitMax);
  for (int c = (isPbPb ? 4 : 0); c >= (isPbPb ? 1 : 0); c--) {
    std::string pre = isPbPb ? Form("_C%d", c) : "";
    for (int i = 0; i < 3; i++) fprintf(h, "double HLTFitParam%s_%d = %g;\n", pre.c_str(), i, r[c].p[i]);
    fprintf(h, "\n");
    for (int i = 0; i < 3; i++) fprintf(h, "double e_HLTFitParam%s_%d = %g;\n", pre.c_str(), i, r[c].e[i]);
    fprintf(h, "\n////////////////////////////////////////\n\n");
  }
  fclose(h);
  printf("  wrote %s\n", fn.c_str());
}

void plotMuonTriggerEfficiency()
{
  initPlotStyle();
  gSystem->mkdir(outDir.c_str(), true);

  Sample mb  {"PbPb_MinBias",   "PbPb/platinum/PbPb_MinBias_triggerEffScan_mu12_tight.root",
              "HIMinimumBias0", "HLT_HIL3Mu12", 5, "50-90%"};
  Sample hp  {"PbPb_HardProbes", "PbPb/latest/triggerEff/PbPb_HardProbes_triggerEffScan_mu12_tight_onlyOneMuonPerEvent.root",
              "HardProbes", "HLT_HIL3Mu12", 5, "50-80%"};
  Sample hp80{"PbPb_HardProbes_Jet80", "PbPb/latest/triggerEff/PbPb_HardProbes_triggerEffScan_Jet80HLT_mu12_tight_onlyOneMuonPerEvent.root",
              "HardProbes + Jet80", "HLT_HIL3Mu12", 5, "50-80%"};
  Sample zb  {"pp_ZeroBias",    "pp/obsidian/pp_ZeroBias_triggerEffScan_mu12_tight.root",
              "ZeroBias", "HLT_HIL3Mu12", 1, ""};
  Sample eg  {"pp_HighEGJet_Jet60", "pp/latest/triggerEff/pp_HighEGJet_triggerEffScan_Jet60HLT_mu12_tight_onlyOneMuonPerEvent.root",
              "HighEGJet + Jet60", "HLT_HIL3Mu12", 1, ""};

  // pp ZeroBias through pp_scan.C's h_muTrigEff histograms (2026-09-21 scan):
  // prescale-1 events only, unweighted, no second muon above 10 GeV
  Sample zb2{"pp_ZeroBias_2026-9-21",
             "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/pp/pp_MinBias_caloJets_manualJEC_mu12_pTmu-15to999_tight_deltaR-40_jetTrkMaxFilter_WDecayFilter_2026-9-21.root",
             "HIZeroBias1", "HLT_HIL3Mu12", 1, "", true};

  // PbPb MinBias Parts 1-4 (2026-09-19/20/21 scans) through PbPb_scan.C's
  // h_muTrigEff histograms, slices merged into the nominal classes by
  // mergeMuTrigEff_coarseCent.C; C0 is 0-80%
  Sample mb3{"PbPb_MinBias_Parts1-4",
             "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/triggerEfficiency/PbPb_MinBias_Parts1-4_muTrigEff_nominalCentBins_2026-9-21.root",
             "HIMinimumBias0, Parts 1-4", "HLT_HIL3Mu12", 5, "50-80%", true, "0-80%"};

  std::vector<TGraphAsymmErrors*> gHP, gHP80, gMB, gZB, gEG, gZB2, gMB3;
  auto rHP   = runSample(hp,   gHP);
  auto rHP80 = runSample(hp80, gHP80);
  auto rMB   = runSample(mb,   gMB);
  auto rZB   = runSample(zb,   gZB);
  auto rEG   = runSample(eg,   gEG);
  auto rZB2  = runSample(zb2,  gZB2);
  auto rMB3  = runSample(mb3,  gMB3);

  writeHeader("PbPb_MinBias", rMB, true);
  writeHeader("PbPb_HardProbes", rHP, true);
  writeHeader("pp_HighEGJet_Jet60", rEG, false);
  writeHeader("pp_ZeroBias_2026-9-21", rZB2, false);
  writeHeader("PbPb_MinBias_Parts1-4", rMB3, true);
}
