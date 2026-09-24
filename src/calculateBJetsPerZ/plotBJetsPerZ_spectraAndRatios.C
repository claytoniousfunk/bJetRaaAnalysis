// Absolute b-jet spectra for the four correction configurations, with the ratio
// to the no-correction baseline beneath each one. Companion to
// plotBJetsPerZ_factorizeFakeAndUnfold.C, which shows only the ratios.
//
// pp is included as its own panel so the absolute size of each correction can
// be read off per system, not only its residue in the PbPb/pp ratio.
//
// NOTE, 2026-09-24. An earlier version of this header claimed the unfolding
// roughly halves the 80-90 GeV yield in both pp and PbPb. That was edge
// distortion, not unfolding: the b-jet axis then began at 80, exactly where the
// response matrix began, so the lowest bin had nowhere to receive migration
// from. With the axis extended to 60 and 60-80 acting as buffer bins, the
// unfolding correction above 80 is a few percent, and the fake-ptRel
// subtraction is the larger of the two corrections below ~150 GeV.
//
// LOG Y. Used here against the usual preference for linear axes, because the
// spectra fall over four and a half decades and nothing is visible otherwise.
// The ratio pads below are linear. No bin is negative in any configuration, so
// the log axis is not hiding one -- the only pathological bin is C3's 300-500,
// which is 2.8e-15 (empty) in the un-unfolded runs and is dropped throughout.

#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TPad.h"
#include <cstdio>

#include "../../headers/plotting/plotStyle.h"

namespace {

const int NCfg = 4;
const char *cfgFile[NCfg] = {"fac_none","fac_fakeOnly","fac_unfoldOnly","fac_both"};
const char *cfgLab [NCfg] = {"neither correction",
                             "fake-ptRel subtraction only",
                             "unfolding only (1 iteration)",
                             "both (nominal)"};
// Okabe-Ito; index 4 (yellow) skipped, it is unusable on white (see plotStyle.h)
const int cfgOkabe[NCfg] = {0,2,6,3};
const int cfgMark [NCfg] = {markFilledCircle,markOpenSquare,markFilledDiamond,markOpenCircle};

const int NPan = 5;
const char *pan   [NPan] = {"pp","C1","C2","C3","C4"};
const char *panLab[NPan] = {"pp","PbPb 0-10%","PbPb 10-30%","PbPb 30-50%","PbPb 50-80%"};

TFile *fIn[NCfg] = {nullptr,nullptr,nullptr,nullptr};

// The unfolding runs on the full 60-500 axis, but 60-70 and 70-80 are buffer
// bins that exist only to absorb migration from below the measurement -- they
// are not results and are never plotted. C3's top bin is empty before
// unfolding; see the header note.
const double quoteFromPt = 80.;

bool skipBin(const char *p, double lo)
{
  if(lo < quoteFromPt - 1e-6) return true;
  return TString(p) == "C3" && lo > 299.;
}

TGraphErrors *graphOf(TH1D *h, TH1D *hRef, const char *p, double dx)
{
  TGraphErrors *g = new TGraphErrors();
  int n = 0;
  for(int i = 1; i <= h->GetNbinsX(); i++){
    const double lo = h->GetBinLowEdge(i), hi = h->GetBinLowEdge(i+1);
    if(skipBin(p,lo)) continue;
    const double v = h->GetBinContent(i);
    if(v <= 0.) continue;
    double y = v, e = h->GetBinError(i);
    if(hRef){
      const double r = hRef->GetBinContent(i);
      if(r <= 0.) continue;
      y = v/r;
      e = (v != 0.) ? fabs(y)*fabs(h->GetBinError(i)/v) : 0.;  // ref held fixed
    }
    g->SetPoint(n,0.5*(lo+hi) + dx*(hi-lo),y);
    g->SetPointError(n,0.,e);
    n++;
  }
  return g;
}

void styleG(TGraphErrors *g, int cfg, double size = 1.1)
{
  const int col = TColor::GetColor(okabeHex[cfgOkabe[cfg]]);
  g->SetLineColor(col); g->SetMarkerColor(col);
  g->SetMarkerStyle(cfgMark[cfg]); g->SetMarkerSize(size); g->SetLineWidth(2);
}

} // namespace

void plotBJetsPerZ_spectraAndRatios()
{
  initPlotStyle();

  const TString dir = "/tmp/claude-1000/-home-clayton-Analysis-code-bJetRaaAnalysis/"
                      "e860c03f-c320-4045-8a31-9f3ba9c591cf/scratchpad/";
  for(int k = 0; k < NCfg; k++){
    fIn[k] = TFile::Open(dir + cfgFile[k] + ".root");
    if(!fIn[k] || fIn[k]->IsZombie()){
      printf("cannot open %s%s.root\n",dir.Data(),cfgFile[k]);
      return;
    }
  }

  // PORTRAIT, 2 columns x 3 rows. A wide 3x2 canvas came out with a third of
  // the PDF page blank: ROOT fits the canvas onto its paper, and gStyle->
  // SetPaperSize does not override that the way one would expect. A portrait
  // canvas matches the page and fills it.
  TCanvas *c = new TCanvas("cSR","",1100,1400);

  // Pads are placed by hand in canvas coordinates rather than with Divide():
  // nesting a pad pair inside a Divide() cell squeezed every panel.
  const double split = 0.34;   // fraction of each cell given to the ratio pad

  for(int j = 0; j < NPan; j++){

    const int col = j % 2, row = j / 2;
    const double x0 = col/2.,        x1 = (col+1)/2.;
    const double y1 = 1. - row/3.,   y0 = 1. - (row+1)/3.;
    const double ys = y0 + split*(y1-y0);

    c->cd();
    TPad *pTop = new TPad(Form("pTop%d",j),"",x0,ys,x1,y1);
    TPad *pBot = new TPad(Form("pBot%d",j),"",x0,y0,x1,ys);
    pTop->SetLeftMargin(0.19); pTop->SetRightMargin(0.03);
    pTop->SetTopMargin(0.07);  pTop->SetBottomMargin(0.02);
    pTop->SetLogy();
    pBot->SetLeftMargin(0.19); pBot->SetRightMargin(0.03);
    pBot->SetTopMargin(0.02);  pBot->SetBottomMargin(0.34);
    pTop->Draw(); pBot->Draw();

    // ---- absolute spectra -------------------------------------------------
    pTop->cd();
    TH1D *frT = new TH1D(Form("frT%d",j),"",1,80.,500.);
    frT->SetStats(0);
    frT->GetYaxis()->SetRangeUser(5e-5,20.);
    frT->GetYaxis()->SetTitle("b-jet spectrum  [per Z, per GeV]");
    frT->GetXaxis()->SetLabelSize(0.);
    frT->GetYaxis()->SetTitleSize(0.062); frT->GetYaxis()->SetLabelSize(0.055);
    frT->GetYaxis()->SetTitleOffset(1.25);
    frT->Draw("axis");

    for(int k = 0; k < NCfg; k++){
      TH1D *h = (TH1D*) fIn[k]->Get(Form("bJetSpectrum_%s",pan[j]));
      if(!h) continue;
      TGraphErrors *g = graphOf(h,nullptr,pan[j],0.07*(k-1.5));
      styleG(g,k);
      g->Draw("pz same");
    }

    TLatex t; t.SetNDC(); t.SetTextFont(42);
    t.SetTextSize(0.075); t.DrawLatex(0.22,0.86,panLab[j]);

    // ---- ratio to the no-correction baseline ------------------------------
    pBot->cd();
    TH1D *frB = new TH1D(Form("frB%d",j),"",1,80.,500.);
    frB->SetStats(0);
    frB->GetYaxis()->SetRangeUser(0.25,1.65);
    frB->GetXaxis()->SetTitle("jet p_{T} [GeV]");
    frB->GetYaxis()->SetTitle("/ baseline");
    frB->GetXaxis()->SetTitleSize(0.115); frB->GetXaxis()->SetLabelSize(0.100);
    frB->GetYaxis()->SetTitleSize(0.105); frB->GetYaxis()->SetLabelSize(0.100);
    frB->GetXaxis()->SetTitleOffset(1.10);
    frB->GetYaxis()->SetTitleOffset(0.72);
    frB->GetYaxis()->SetNdivisions(505);
    frB->Draw("axis");

    TLine *one = new TLine(80.,1.,500.,1.);
    one->SetLineStyle(2); one->SetLineColor(kGray+2); one->Draw();

    TH1D *hRef = (TH1D*) fIn[0]->Get(Form("bJetSpectrum_%s",pan[j]));
    for(int k = 1; k < NCfg; k++){
      TH1D *h = (TH1D*) fIn[k]->Get(Form("bJetSpectrum_%s",pan[j]));
      if(!h || !hRef) continue;
      TGraphErrors *g = graphOf(h,hRef,pan[j],0.07*(k-1.5));
      styleG(g,k,1.0);
      g->Draw("pz same");
    }
  }

  // ---- legend in the spare sixth cell ------------------------------------
  c->cd();
  TPad *pLeg = new TPad("pLeg","",0.5,0.,1.,1/3.);
  pLeg->SetFillStyle(0); pLeg->Draw(); pLeg->cd();
  TLegend *leg = new TLegend(0.08,0.52,0.98,0.86);
  leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.055);
  for(int k = 0; k < NCfg; k++){
    TGraphErrors *g = new TGraphErrors(1);
    g->SetPoint(0,-1.,-1.);
    styleG(g,k,1.3);
    leg->AddEntry(g,cfgLab[k],"lp");
  }
  leg->Draw();

  TLatex tc; tc.SetNDC(); tc.SetTextFont(42);
  tc.SetTextSize(0.062);
  tc.DrawLatex(0.08,0.92,"b-jet spectra per Z");
  tc.SetTextSize(0.046);
  tc.DrawLatex(0.08,0.42,"lower pads: ratio to the run with");
  tc.DrawLatex(0.08,0.35,"neither correction applied");
  tc.DrawLatex(0.08,0.24,"fake-ptRel changes the b purity;");
  tc.DrawLatex(0.08,0.17,"unfolding changes the spectrum");

  c->SaveAs("../../figures/JetsPerZ/bJetsPerZ_spectraAndRatios.pdf");
}
