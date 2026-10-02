// Calo vs PF jet energy scale in PbPb and pp data, from caloVsPFJetScale_scan.C.
//
// Asks one question: is the calo jet scale (relative to PF) lower in PbPb than
// in pp, and by enough to explain the calo-jet R_AA sitting 25-35% below the
// published values in 50-80% above 200 GeV (calculateJetsPerZ.cc with
// useTAANormalization)?
//
// For each class, <calo/PF> per PF-pT bin; the bottom panel is the double ratio
//     s = <calo/PF>_PbPb / <calo/PF>_pp.
// If the PbPb calo scale is low by s relative to pp, the measured PbPb spectrum
// is (1/s) f(pT/s) ~ s^(n-1) f(pT) for f ~ pT^-n, so the calo R_AA is biased by
//     R_AA(measured) / R_AA(true) = s^(n-1).
// n is fit to the pp PF spectrum in the same window. The table at the end prints
// s and s^(n-1) per class; compare the latter with the deficit (~0.65-0.70 in
// 50-80% against ATLAS / CMS).
//
// READ BOTH CONVENTIONS. Figure 1 bins in PF pT (calo/PF), figure 2 bins in calo
// pT (PF/calo, inverted for display so both read as calo/PF). A genuine scale
// shift moves both the same way after the inversion; a PbPb-vs-pp RESOLUTION
// difference, which biases each convention toward the binning variable, moves
// them in opposite directions. Only a shift seen in both is a scale.
// Figure 3 is the same as figure 1 without JEC (raw pT): if the double ratio is
// already there, the cause is the background subtraction, not the JEC.
//
// Errors are statistical on the mean. pp and PbPb are independent, so the double
// ratio propagates both.

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TColor.h"
#include <cmath>
#include "TCanvas.h"
#include "TPad.h"
#include "TLatex.h"
#include "TNamed.h"
#include "TSystem.h"
#include "../../../../headers/plotting/plotStyle.h"
#include "../../../../headers/plotting/ratioPanel.h"

namespace cvpfPlot {

const TString dirIn  = "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/";
TString filePbPb = dirIn + "PbPb/PbPb_HardProbes_caloVsPFJetScale.root";
TString filePP   = dirIn + "pp/pp_HighEGJet_caloVsPFJetScale.root";
TString dirOut = "/home/clayton/Analysis/code/bJetRaaAnalysis/figures/jetCollection/";

const double xLo = 100., xHi = 700.;          // trigger plateau upward
const double tabWin[2][2] = {{200., 320.}, {320., 500.}};
const double fitLo = 150., fitHi = 500.;      // spectral index window

// PbPb classes drawn, as booked in the scan: C1..C4
const int   nDraw = 4;
const int   cls[nDraw]   = {1, 2, 3, 4};
const char *clsLab[nDraw] = {"0-10%", "10-30%", "30-50%", "50-80%"};
const int   clsMark[nDraw] = {markFilledSquare, markFilledDiamond, markOpenSquare, markOpenDiamond};
const char *clsHex[nDraw]  = {"#D55E00", "#E69F00", "#009E73", "#0072B2"};

// Profile of the ratio axis: mean and its error per x bin, as a graph holding
// only bins inside [xLo, xHi] with >= 10 pairs. A graph rather than a TH1 so
// that bins outside the frame or without pairs are absent, not drawn at zero
// in the pad margin. invert = true returns 1/<PF/calo> with the error
// propagated, so figure 2 reads in the same direction as figure 1.
TGraphErrors* meanProfile(TH2D *h, const char *name, bool invert = false)
{
  TGraphErrors *g = new TGraphErrors(); g->SetName(name);
  for(int i = 1; i <= h->GetNbinsX(); i++){
    const double x = h->GetXaxis()->GetBinCenter(i);
    if(x < xLo || x > xHi) continue;
    TH1D *y = h->ProjectionY(Form("%s_y%d", name, i), i, i);
    if(y->GetEntries() >= 10){
      double m = y->GetMean(), e = y->GetMeanError();
      if(invert && m > 0){ e = e/(m*m); m = 1./m; }
      const int k = g->GetN();
      g->SetPoint(k, x, m); g->SetPointError(k, 0.5*h->GetXaxis()->GetBinWidth(i), e);
    }
    delete y;
  }
  return g;
}

// PbPb / pp point by point, at x values present in both; independent samples,
// so both errors propagate.
TGraphErrors* doubleRatio(TGraphErrors *a, TGraphErrors *p, const char *name)
{
  TGraphErrors *r = new TGraphErrors(); r->SetName(name);
  for(int i = 0; i < a->GetN(); i++)
    for(int j = 0; j < p->GetN(); j++){
      if(std::fabs(a->GetX()[i] - p->GetX()[j]) > 1e-6) continue;
      const double va = a->GetY()[i], vp = p->GetY()[j];
      if(va <= 0 || vp <= 0) continue;
      const double v = va/vp, e = v*std::sqrt(std::pow(a->GetEY()[i]/va, 2) + std::pow(p->GetEY()[j]/vp, 2));
      const int k = r->GetN();
      r->SetPoint(k, a->GetX()[i], v); r->SetPointError(k, a->GetEX()[i], e);
    }
  return r;
}

void styleG(TGraphErrors *g, const char *hex, int mark, double size = 1.2)
{
  int c = TColor::GetColor(hex);
  g->SetLineColor(c); g->SetMarkerColor(c); g->SetMarkerStyle(mark);
  g->SetMarkerSize(size); g->SetLineWidth(2);
}

// mean and median of the ratio over an x window, entries pooled
void windowStats(TH2D *h, double lo, double hi, double &mean, double &err,
                 double &median, double &n, bool invert = false)
{
  TH1D *y = h->ProjectionY(Form("%s_w", h->GetName()),
                           h->GetXaxis()->FindBin(lo + 1e-6), h->GetXaxis()->FindBin(hi - 1e-6));
  n = y->GetEntries();
  mean = y->GetMean(); err = y->GetMeanError();
  double q = 0.5; median = 0.;
  if(n > 0) y->GetQuantiles(1, &median, &q);
  if(invert && mean > 0){ err = err/(mean*mean); mean = 1./mean; median = median > 0 ? 1./median : 0.; }
  delete y;
}

TH2D* get2(TFile *f, const char *name)
{
  TH2D *h = nullptr; f->GetObject(name, h);
  if(!h) printf("ERROR: %s missing from %s\n", name, f->GetName());
  return h;
}

void drawFigure(TFile *fA, TFile *fP, const char *hName, bool invert,
                const char *yTitle, const char *xTitle, const char *note, const TString &pdf,
                double yLo = 0.80, double yHi = 1.12)
{
  TH2D *hp = get2(fP, Form("%s_C0", hName));
  if(!hp) return;
  TGraphErrors *pp = meanProfile(hp, Form("pp_%s", hName), invert);
  styleG(pp, hexData, markFilledCircle);

  TGraphErrors *aa[nDraw], *dr[nDraw];
  for(int k = 0; k < nDraw; k++){
    TH2D *h = get2(fA, Form("%s_C%d", hName, cls[k]));
    if(!h) return;
    aa[k] = meanProfile(h, Form("aa%d_%s", k, hName), invert);
    styleG(aa[k], clsHex[k], clsMark[k]);
    dr[k] = doubleRatio(aa[k], pp, Form("dr%d_%s", k, hName));
    styleG(dr[k], clsHex[k], clsMark[k]);
  }

  TCanvas *c = new TCanvas(Form("c_%s", hName), "", 700, 800);
  TPad *top, *bot; splitPads(top, bot);

  top->cd(); top->SetLogx(); top->SetTicks(1, 1);
  TH1F *fr = top->DrawFrame(xLo, yLo, xHi, yHi);
  fr->GetYaxis()->SetTitle(yTitle);
  fr->GetYaxis()->SetTitleSize(0.055); fr->GetYaxis()->SetLabelSize(0.050);
  fr->GetYaxis()->SetTitleOffset(1.35);
  fr->GetXaxis()->SetLabelSize(0);
  if(pp->GetN()) pp->Draw("p same");
  for(int k = 0; k < nDraw; k++) if(aa[k]->GetN()) aa[k]->Draw("p same");
  unityLine(xLo, xHi)->Draw();
  TLegend *l = makeLegend(0.55, 0.06, 0.92, 0.42, 0.048);
  l->AddEntry(pp, "pp (ak4Calo / ak4PF)", "pe");
  for(int k = 0; k < nDraw; k++) l->AddEntry(aa[k], Form("PbPb %s", clsLab[k]), "pe");
  l->Draw();
  TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.045);
  t.DrawLatex(0.20, 0.86, "Data, matched jets, #DeltaR < 0.2, |#eta| < 1.6");
  t.DrawLatex(0.20, 0.80, note);

  bot->cd(); bot->SetLogx(); bot->SetTicks(1, 1);
  TH1F *fb = bot->DrawFrame(xLo, 0.88, xHi, 1.06);
  styleRatioAxes((TH1D*) fb, xTitle, "PbPb / pp");
  fb->GetXaxis()->SetMoreLogLabels(); fb->GetXaxis()->SetNoExponent();
  for(int k = 0; k < nDraw; k++) if(dr[k]->GetN()) dr[k]->Draw("p same");
  unityLine(xLo, xHi)->Draw();

  savePdfTight(c, pdf);
}

} // namespace cvpfPlot

void plotCaloOverPFJetScale_PbPbVsPP()
{
  using namespace cvpfPlot;
  initPlotStyle();

  TFile *fA = TFile::Open(filePbPb), *fP = TFile::Open(filePP);
  if(!fA || fA->IsZombie() || !fP || fP->IsZombie()){
    printf("ERROR: cannot open %s or %s -- hadd the condor outputs (condor_caloVsPFJetScale.py) first\n",
           filePbPb.Data(), filePP.Data());
    return;
  }
  for(TFile *f : {fA, fP}){
    TNamed *p = nullptr; f->GetObject("provenance", p);
    printf("---- %s\n%s\n", f->GetName(), p ? p->GetTitle() : "(no provenance)");
  }
  gSystem->mkdir(dirOut, kTRUE);

  drawFigure(fA, fP, "h_caloOverPF_vsPFPt", false,
             "#LT p_{T}^{calo} / p_{T}^{PF} #GT", "PF jet p_{T} [GeV]",
             "manual JEC, binned in PF p_{T}", dirOut + "caloOverPF_vsPFPt_PbPbVsPP.pdf");
  drawFigure(fA, fP, "h_PFOverCalo_vsCaloPt", true,
             "1 / #LT p_{T}^{PF} / p_{T}^{calo} #GT", "calo jet p_{T} [GeV]",
             "manual JEC, binned in calo p_{T}", dirOut + "caloOverPF_vsCaloPt_PbPbVsPP.pdf");
  drawFigure(fA, fP, "h_rawCaloOverRawPF_vsPFPt", false,
             "#LT raw p_{T}^{calo} / raw p_{T}^{PF} #GT", "PF jet p_{T} [GeV]",
             "no JEC, binned in PF p_{T}", dirOut + "rawCaloOverRawPF_vsPFPt_PbPbVsPP.pdf", 0.60, 1.05);

  // ---- spectral index from the pp PF spectrum (all PF jets, triggered sample)
  TH1D *spec = nullptr; fP->GetObject("h_pfPt_C0", spec);
  double n = 5.;
  if(spec){
    TH1D *d = (TH1D*) spec->Clone("spec_density"); d->SetDirectory(nullptr);
    d->Scale(1., "width");
    TF1 pw("pw", "[0]*pow(x/200.,-[1])", fitLo, fitHi);
    pw.SetParameters(d->GetBinContent(d->FindBin(200.)), 5.);
    d->Fit(&pw, "RQ0");
    n = pw.GetParameter(1);
    printf("\npp PF spectrum, %g-%g GeV: n = %.2f +- %.2f\n", fitLo, fitHi, n, pw.GetParError(1));
  }
  else printf("\nWARNING: h_pfPt_C0 missing from pp file, using n = 5\n");

  // ---- table ----------------------------------------------------------------
  struct Conv { const char *h; bool inv; const char *lab; };
  Conv conv[3] = {{"h_caloOverPF_vsPFPt", false, "calo/PF vs PF pT"},
                  {"h_PFOverCalo_vsCaloPt", true, "1/(PF/calo) vs calo pT"},
                  {"h_rawCaloOverRawPF_vsPFPt", false, "raw calo/raw PF vs PF pT"}};
  for(auto &cv : conv){
    printf("\n%s\n", cv.lab);
    for(auto &w : tabWin){
      double mP, eP, medP, nP;
      TH2D *hp = get2(fP, Form("%s_C0", cv.h)); if(!hp) continue;
      windowStats(hp, w[0], w[1], mP, eP, medP, nP, cv.inv);
      printf("  %3.0f-%3.0f GeV   pp: mean %.4f +- %.4f  median %.4f  (%.0f pairs)\n",
             w[0], w[1], mP, eP, medP, nP);
      for(int k = 0; k < nDraw; k++){
        double mA, eA, medA, nA;
        TH2D *h = get2(fA, Form("%s_C%d", cv.h, cls[k])); if(!h) continue;
        windowStats(h, w[0], w[1], mA, eA, medA, nA, cv.inv);
        if(nA < 10 || nP < 10 || mA <= 0 || mP <= 0){
          printf("    %-7s too few pairs (%.0f PbPb, %.0f pp)\n", clsLab[k], nA, nP); continue; }
        double s = mA/mP, es = s*std::sqrt(std::pow(eA/mA, 2) + std::pow(eP/mP, 2));
        printf("    %-7s mean %.4f  median %.4f  (%6.0f)   s = %.4f +- %.4f   s(median) = %.4f   "
               "R_AA factor s^(n-1) = %.3f\n",
               clsLab[k], mA, medA, nA, s, es, medP > 0 ? medA/medP : 0., std::pow(s, n - 1.));
      }
    }
  }
  printf("\nThe 50-80%% calo R_AA deficit to explain: ~0.65-0.70 (ours / ATLAS, 200-320 GeV).\n");
}
