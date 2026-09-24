// Factorize the two corrections that act on the b-jets-per-Z measurement:
//
//   fake-ptRel subtraction : removes the combinatorial muon+jet ptRel template
//                            (T3) from the measured ptRel before the template
//                            fit, so it changes the b PURITY, PbPb only.
//   unfolding              : Bayes, 1 iteration, on the variable-width b-jet
//                            response, so it changes the b-jet SPECTRUM, both
//                            systems.
//
// Four runs of calculateBJetsPerZ.cc, differing only in those two flags:
//   fac_none.root        both off  -- the baseline
//   fac_fakeOnly.root    fake sub only
//   fac_unfoldOnly.root  unfolding only
//   fac_both.root        both on   -- the nominal chain
//
// Three figures:
//   _spectra   b-jet spectrum of each configuration / baseline, per class
//   _jetsPerZ  PbPb/pp b-jets per Z, all four configurations, per class
//   _closure   whether the two corrections factorize: (both/none) against
//              (fakeOnly/none) x (unfoldOnly/none). They act on different
//              factors of the same product (purity vs spectrum), so they SHOULD
//              multiply; a departure means one correction changes what the
//              other does, which is worth knowing before either is quoted as an
//              independent systematic.
//
// C3's 300-500 bin is dropped throughout: its measured content is ~1e-16, an
// empty bin, so anything the unfolding puts there comes from the prior.

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
                             "unfolding only (1 iter)",
                             "both (nominal)"};
// Okabe-Ito: black, sky blue, vermillion, bluish green. Index 4 (yellow) is
// unusable on white and is skipped, per plotStyle.h.
const int cfgOkabe [NCfg] = {0,2,6,3};
const int cfgMark  [NCfg] = {markFilledCircle,markOpenSquare,markFilledDiamond,markOpenCircle};

const int NCls = 4;
const char *cls   [NCls] = {"C1","C2","C3","C4"};
const char *clsLab[NCls] = {"PbPb 0-10%","PbPb 10-30%","PbPb 30-50%","PbPb 50-80%"};

TFile *fIn[NCfg] = {nullptr,nullptr,nullptr,nullptr};

// The b-jet chain unfolds on the full 60-500 axis so the 60-70 and 70-80 bins
// can absorb migration from below the measurement, but those two bins exist
// only to soak up edge distortion and are not results. Nothing below 80 is
// plotted or tabulated.
const double quoteFromPt = 80.;

bool skipBin(const char *cl, double lo)
{
  if(lo < quoteFromPt - 1e-6) return true;              // buffer bin
  return TString(cl) == "C3" && lo > 299.;              // empty before unfolding
}

// Graph of h, optionally divided by hRef, with empty bins dropped.
TGraphErrors *graphOf(TH1D *h, TH1D *hRef, const char *cl, double dx)
{
  TGraphErrors *g = new TGraphErrors();
  int n = 0;
  for(int i = 1; i <= h->GetNbinsX(); i++){
    const double lo = h->GetBinLowEdge(i), hi = h->GetBinLowEdge(i+1);
    if(skipBin(cl,lo)) continue;
    const double v = h->GetBinContent(i);
    if(v == 0.) continue;
    double y = v, e = h->GetBinError(i);
    if(hRef){
      const double r = hRef->GetBinContent(i);
      if(r == 0.) continue;
      y = v/r;
      e = fabs(y)*sqrt(pow(h->GetBinError(i)/v,2));   // ref treated as fixed
    }
    // nudge x so overlapping configurations stay readable
    g->SetPoint(n,0.5*(lo+hi) + dx*(hi-lo),y);
    g->SetPointError(n,0.,e);
    n++;
  }
  return g;
}

void styleG(TGraphErrors *g, int cfg)
{
  const int col = TColor::GetColor(okabeHex[cfgOkabe[cfg]]);
  g->SetLineColor(col); g->SetMarkerColor(col);
  g->SetMarkerStyle(cfgMark[cfg]); g->SetMarkerSize(1.2); g->SetLineWidth(2);
}

// One 2x2 canvas. mode 0 = spectra ratio to baseline, 1 = jets per Z absolute.
void drawPanels(const char *cname, int mode, double yLo, double yHi,
                const char *yTitle, const char *outPdf)
{
  TCanvas *c = new TCanvas(cname,"",1000,900);
  c->Divide(2,2,0.001,0.001);

  for(int j = 0; j < NCls; j++){
    TPad *p = (TPad*) c->cd(j+1);
    p->SetLeftMargin(0.15); p->SetRightMargin(0.03);
    p->SetTopMargin(0.07);  p->SetBottomMargin(0.13);

    TH1D *fr = new TH1D(Form("fr_%s_%s",cname,cls[j]),"",1,80.,500.);
    fr->SetStats(0);
    fr->GetYaxis()->SetRangeUser(yLo,yHi);
    fr->GetXaxis()->SetTitle("jet p_{T} [GeV]");
    fr->GetYaxis()->SetTitle(yTitle);
    fr->GetXaxis()->SetTitleSize(0.048); fr->GetYaxis()->SetTitleSize(0.048);
    fr->GetXaxis()->SetLabelSize(0.043); fr->GetYaxis()->SetLabelSize(0.043);
    fr->GetXaxis()->SetTitleOffset(1.10);
    fr->GetYaxis()->SetTitleOffset(1.40);
    fr->Draw("axis");

    TLine *one = new TLine(80.,1.,500.,1.);
    one->SetLineStyle(2); one->SetLineColor(kGray+2); one->Draw();

    TLegend *leg = nullptr;
    if(j == 0){
      leg = new TLegend(0.34,0.66,0.96,0.90);
      leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.040);
    }

    const int cfg0 = (mode == 0) ? 1 : 0;   // spectra: baseline is the divisor
    for(int k = cfg0; k < NCfg; k++){
      if(!fIn[k]) continue;
      const char *stem = (mode == 0) ? "bJetSpectrum" : "bFraction";
      TH1D *h    = (TH1D*) fIn[k]->Get(Form("%s_%s",stem,cls[j]));
      TH1D *hRef = (mode == 0) ? (TH1D*) fIn[0]->Get(Form("%s_%s",stem,cls[j])) : nullptr;
      if(!h || (mode == 0 && !hRef)) continue;
      TGraphErrors *g = graphOf(h,hRef,cls[j],0.06*(k-1.5));
      styleG(g,k);
      g->Draw("pz same");
      if(leg) leg->AddEntry(g,cfgLab[k],"lp");
    }
    if(leg) leg->Draw();

    TLatex t; t.SetNDC(); t.SetTextFont(42);
    t.SetTextSize(0.052); t.DrawLatex(0.19,0.88,clsLab[j]);
  }

  c->SaveAs(outPdf);
}

} // namespace

void plotBJetsPerZ_factorizeFakeAndUnfold()
{
  initPlotStyle();

  const TString dir = "/tmp/claude-1000/-home-clayton-Analysis-code-bJetRaaAnalysis/"
                      "e860c03f-c320-4045-8a31-9f3ba9c591cf/scratchpad/";
  for(int k = 0; k < NCfg; k++){
    fIn[k] = TFile::Open(dir + cfgFile[k] + ".root");
    if(!fIn[k] || fIn[k]->IsZombie()){
      printf("cannot open %s%s.root\n",dir.Data(),cfgFile[k]);
      fIn[k] = nullptr;
    }
  }
  if(!fIn[0] || !fIn[3]){ printf("need at least the baseline and the nominal run\n"); return; }

  // ---- printed tables, so the numbers exist outside the figures ------------
  for(int j = 0; j < NCls; j++){
    printf("\n=== %s ===\n",clsLab[j]);
    printf("   pT bin     b-jet spectrum ratio to baseline        PbPb/pp b jets per Z\n");
    printf("              fakeOnly  unfoldOnly     both      none  fakeOnly unfOnly    both\n");
    TH1D *sN = (TH1D*) fIn[0]->Get(Form("bJetSpectrum_%s",cls[j]));
    TH1D *rN = (TH1D*) fIn[0]->Get(Form("bFraction_%s",cls[j]));
    if(!sN || !rN) continue;
    for(int i = 1; i <= sN->GetNbinsX(); i++){
      const double lo = sN->GetBinLowEdge(i), hi = sN->GetBinLowEdge(i+1);
      if(skipBin(cls[j],lo)) continue;
      printf("  %4.0f-%4.0f ",lo,hi);
      for(int k = 1; k < NCfg; k++){
        TH1D *s = fIn[k] ? (TH1D*) fIn[k]->Get(Form("bJetSpectrum_%s",cls[j])) : nullptr;
        const double d = sN->GetBinContent(i);
        if(s && d != 0.) printf("  %8.3f",s->GetBinContent(i)/d); else printf("      -- ");
      }
      printf("   ");
      for(int k = 0; k < NCfg; k++){
        TH1D *r = fIn[k] ? (TH1D*) fIn[k]->Get(Form("bFraction_%s",cls[j])) : nullptr;
        if(r) printf(" %7.4f",r->GetBinContent(i)); else printf("      --");
      }
      printf("\n");
    }
  }

  // ---- factorization closure ----------------------------------------------
  if(fIn[1] && fIn[2]){
    printf("\n\n=== do the two corrections factorize? ===\n");
    printf("  (both/none) vs (fakeOnly/none)x(unfoldOnly/none), on PbPb/pp b jets per Z\n");
    for(int j = 0; j < NCls; j++){
      printf("\n  %s\n     pT bin    measured   product   meas/prod\n",clsLab[j]);
      TH1D *r0 = (TH1D*) fIn[0]->Get(Form("bFraction_%s",cls[j]));
      TH1D *r1 = (TH1D*) fIn[1]->Get(Form("bFraction_%s",cls[j]));
      TH1D *r2 = (TH1D*) fIn[2]->Get(Form("bFraction_%s",cls[j]));
      TH1D *r3 = (TH1D*) fIn[3]->Get(Form("bFraction_%s",cls[j]));
      if(!r0 || !r1 || !r2 || !r3) continue;
      for(int i = 1; i <= r0->GetNbinsX(); i++){
        const double lo = r0->GetBinLowEdge(i), hi = r0->GetBinLowEdge(i+1);
        if(skipBin(cls[j],lo)) continue;
        const double b = r0->GetBinContent(i);
        if(b == 0.) continue;
        const double meas = r3->GetBinContent(i)/b;
        const double prod = (r1->GetBinContent(i)/b)*(r2->GetBinContent(i)/b);
        printf("   %4.0f-%4.0f  %9.4f %9.4f %9.4f\n",lo,hi,meas,prod,
               prod != 0. ? meas/prod : 0.);
      }
    }
  }

  drawPanels("cSpec",0,0.4,1.8,"b-jet spectrum / baseline",
             "../../figures/JetsPerZ/bJetsPerZ_factorize_spectra.pdf");
  drawPanels("cRatio",1,0.,2.6,"PbPb / pp  b jets per Z",
             "../../figures/JetsPerZ/bJetsPerZ_factorize_jetsPerZ.pdf");
}
