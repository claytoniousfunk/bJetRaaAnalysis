// b-jets per Z divided by inclusive jets per Z, both as PbPb/pp double ratios.
//
//   numerator    bFraction_C*  from calculateBJetsPerZ.cc, with BOTH b-jet
//                corrections applied: the fake-ptRel (T3) subtraction in the
//                purity fit, and Bayes unfolding, 1 iteration, on the
//                variable-width b-jet response h_..._var_bJets.
//   denominator  r_C*_r from calculateRAA.C, inclusive jets, Bayes unfolding,
//                1 iteration, on the FINE 0-500 response, rebinned onto the
//                coarse axis afterwards (useVariableBinResponse = false,
//                outputOnVarAxis = true).
//
// WHY THE TWO SIDES ARE UNFOLDED DIFFERENTLY. They should not be, and this is a
// known defect of the present figure. Unfolding the inclusive on the variable
// matrix was tried and is wrong: that matrix starts at 60 GeV in BOTH reco and
// gen, while the fine one starts at 0, so it has no truth bins to absorb jets
// of true pT below 60 that reconstruct above 80. On a steeply falling spectrum
// that content piles into the lowest truth bins -- r_C1_r at 80-90 GeV came out
// 2.267 against 0.368 from the fine matrix, a factor of 6, and the fine result
// reproduces the independent 2026-08-21 file to four decimals.
//
// The b-jet numerator's own edge distortion was fixed on 2026-09-24 by
// extending jetPtAxisEdges down to 60 GeV, so its measured and truth axes now
// both span the full response and the 60-70 / 70-80 bins act as buffers. The
// effect was large: the unfolding correction at 80-90 GeV in 0-10% went from
// x0.383 to x0.987. Those two buffer bins have no inclusive partner and are
// dropped from this ratio; the b-jet figures show them, marked.
//
// Since numerator and denominator are each already PbPb/pp, the per-Z
// normalisation N_Z cancels exactly -- the same dimuon-mass definition and the
// same efficiency factors are used on both sides.
//
// C3's 300-500 bin is dropped: the b-jet measurement there is empty before
// unfolding, so the unfolded value is prior, not data.

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

const char *fB = "/home/clayton/Analysis/code/bJetRaaAnalysis/src/calculateBJetsPerZ/"
                 "rootFiles/BJetPbPbToPP/Data/histograms_BJetPbPbToPP.root";
const char *fI = "/home/clayton/Analysis/code/bJetRaaAnalysis/src/newFractionCalculation/"
                 "rootFiles/JetsPerZ/histograms_JetsPerZ_lightJets_rebinned.root";

const int NCls = 4;
const char *cls   [NCls] = {"C1","C2","C3","C4"};
const char *clsLab[NCls] = {"PbPb 0-10%","PbPb 10-30%","PbPb 30-50%","PbPb 50-80%"};
// Okabe-Ito; index 4 (yellow) skipped, unusable on white
const int clsOkabe[NCls] = {0,5,6,3};
const int clsMark [NCls] = {markFilledCircle,markFilledSquare,markFilledDiamond,markOpenCircle};

// The b-jet chain unfolds on 60-500 so that 60-70 and 70-80 absorb migration
// from below the measurement. Those two bins are buffers, not results, and are
// dropped here -- they would in any case have no inclusive partner.
const double quoteFromPt = 80.;

bool skipBin(const char *c, double lo)
{
  if(lo < quoteFromPt - 1e-6) return true;
  return TString(c) == "C3" && lo > 299.;
}

TGraphErrors *graphOf(TH1D *h, TH1D *hDen, const char *c, double dx)
{
  TGraphErrors *g = new TGraphErrors();
  int n = 0;
  for(int i = 1; i <= h->GetNbinsX(); i++){
    const double lo = h->GetBinLowEdge(i), hi = h->GetBinLowEdge(i+1);
    if(skipBin(c,lo)) continue;
    const double v = h->GetBinContent(i);
    if(v == 0.) continue;
    double y = v, e = h->GetBinError(i);
    if(hDen){
      // The b-jet axis now starts at 60 and the inclusive one at 80, so the two
      // histograms no longer share bin INDICES. Look the denominator up by the
      // bin centre; a b-jet bin with no inclusive partner (the 60-80 buffer) is
      // skipped rather than paired with whatever index happens to line up.
      const int id = hDen->FindBin(0.5*(lo+hi));
      if(id < 1 || id > hDen->GetNbinsX()) continue;
      if(fabs(hDen->GetBinLowEdge(id) - lo) > 1e-6 ||
         fabs(hDen->GetBinLowEdge(id+1) - hi) > 1e-6) continue;   // edges differ
      const double d = hDen->GetBinContent(id);
      if(d == 0.) continue;
      y = v/d;
      // both sides carry statistical errors; treat them as uncorrelated
      const double rb = (v != 0.) ? h->GetBinError(i)/v : 0.;
      const double rd = (d != 0.) ? hDen->GetBinError(id)/d : 0.;
      e = fabs(y)*sqrt(rb*rb + rd*rd);
    }
    g->SetPoint(n,0.5*(lo+hi) + dx*(hi-lo),y);
    g->SetPointError(n,0.,e);
    n++;
  }
  return g;
}

void styleG(TGraphErrors *g, int j, int mark = -1, double size = 1.3)
{
  const int col = TColor::GetColor(okabeHex[clsOkabe[j]]);
  g->SetLineColor(col); g->SetMarkerColor(col);
  g->SetMarkerStyle(mark < 0 ? clsMark[j] : mark);
  g->SetMarkerSize(size); g->SetLineWidth(2);
}

} // namespace

void plotBOverInclusiveJetsPerZ()
{
  initPlotStyle();

  TFile *fb = TFile::Open(fB);
  TFile *fi = TFile::Open(fI);
  if(!fb || fb->IsZombie() || !fi || fi->IsZombie()){
    printf("cannot open inputs\n  b:   %s\n  inc: %s\n",fB,fI); return;
  }

  // ---------------- printed table ------------------------------------------
  printf("\nb jets per Z / inclusive jets per Z   (each PbPb/pp)\n");
  for(int j = 0; j < NCls; j++){
    TH1D *hb = (TH1D*) fb->Get(Form("bFraction_%s",cls[j]));
    TH1D *hi = (TH1D*) fi->Get(Form("r_%s_r",cls[j]));
    if(!hb || !hi){ printf("  missing histograms for %s\n",cls[j]); continue; }
    printf("\n  %s\n     pT bin       b/Z      incl/Z    b / incl\n",clsLab[j]);
    for(int i = 1; i <= hb->GetNbinsX(); i++){
      const double lo = hb->GetBinLowEdge(i), hi_ = hb->GetBinLowEdge(i+1);
      if(skipBin(cls[j],lo)) continue;
      const int id = hi->FindBin(0.5*(lo+hi_));
      const bool ok = (id >= 1 && id <= hi->GetNbinsX() &&
                       fabs(hi->GetBinLowEdge(id) - lo) < 1e-6 &&
                       fabs(hi->GetBinLowEdge(id+1) - hi_) < 1e-6);
      const double b = hb->GetBinContent(i), d = ok ? hi->GetBinContent(id) : 0.;
      printf("   %4.0f-%4.0f  %8.4f  %8.4f  ",lo,hi_,b,d);
      if(d != 0.) printf("%8.4f\n",b/d); else printf("    --\n");
    }
  }

  // ---------------- figure 1: the four ratios overlaid ---------------------
  TCanvas *c1 = new TCanvas("cBI","",800,700);
  c1->SetLeftMargin(0.17); c1->SetRightMargin(0.04);
  c1->SetTopMargin(0.07);  c1->SetBottomMargin(0.13);

  TH1D *fr = new TH1D("frBI","",1,80.,500.);
  fr->SetStats(0);
  fr->GetYaxis()->SetRangeUser(0.,5.0);
  fr->GetXaxis()->SetTitle("jet p_{T} [GeV]");
  // Spell the double ratio out. "b / inclusive, PbPb/pp" reads as though the
  // flavor ratio were taken first and then PbPb/pp; it is the other way round.
  fr->GetYaxis()->SetTitle("#frac{(b jets per Z)_{PbPb} / (b jets per Z)_{pp}}"
                           "{(incl jets per Z)_{PbPb} / (incl jets per Z)_{pp}}");
  fr->GetXaxis()->SetTitleOffset(1.15);
  fr->GetYaxis()->SetTitleOffset(1.55);
  fr->GetYaxis()->SetTitleSize(0.036);
  fr->Draw("axis");

  TLine *one = new TLine(80.,1.,500.,1.);
  one->SetLineStyle(2); one->SetLineColor(kGray+2); one->Draw();

  TLegend *leg = new TLegend(0.55,0.68,0.94,0.89);
  leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.036);

  for(int j = 0; j < NCls; j++){
    TH1D *hb = (TH1D*) fb->Get(Form("bFraction_%s",cls[j]));
    TH1D *hd = (TH1D*) fi->Get(Form("r_%s_r",cls[j]));
    if(!hb || !hd) continue;
    TGraphErrors *g = graphOf(hb,hd,cls[j],0.07*(j-1.5));
    styleG(g,j);
    g->Draw("pz same");
    leg->AddEntry(g,clsLab[j],"lp");
  }
  leg->Draw();

  TLatex t; t.SetNDC(); t.SetTextFont(42);
  t.SetTextSize(0.034);
  t.DrawLatex(0.17,0.945,"b jets vs inclusive jets, 1 unfolding iteration");
  c1->SaveAs("../../figures/JetsPerZ/bOverInclusiveJetsPerZ.pdf");

  // ---------------- figure 2: the two inputs, per class --------------------
  TCanvas *c2 = new TCanvas("cBI2","",1100,1200);
  const double split = 0.34;

  for(int j = 0; j < NCls; j++){
    const int col = j % 2, row = j / 2;
    const double x0 = col/2., x1 = (col+1)/2.;
    const double y1 = 1. - row/2., y0 = 1. - (row+1)/2.;
    const double ys = y0 + split*(y1-y0);

    c2->cd();
    TPad *pTop = new TPad(Form("tB%d",j),"",x0,ys,x1,y1);
    TPad *pBot = new TPad(Form("bB%d",j),"",x0,y0,x1,ys);
    pTop->SetLeftMargin(0.17); pTop->SetRightMargin(0.03);
    pTop->SetTopMargin(0.08);  pTop->SetBottomMargin(0.02);
    pBot->SetLeftMargin(0.17); pBot->SetRightMargin(0.03);
    pBot->SetTopMargin(0.02);  pBot->SetBottomMargin(0.32);
    pTop->Draw(); pBot->Draw();

    TH1D *hb = (TH1D*) fb->Get(Form("bFraction_%s",cls[j]));
    TH1D *hd = (TH1D*) fi->Get(Form("r_%s_r",cls[j]));
    if(!hb || !hd) continue;

    pTop->cd();
    TH1D *fT = new TH1D(Form("fT%d",j),"",1,80.,500.);
    fT->SetStats(0);
    fT->GetYaxis()->SetRangeUser(0.,1.8);
    fT->GetYaxis()->SetTitle("PbPb / pp  per Z");
    fT->GetXaxis()->SetLabelSize(0.);
    fT->GetYaxis()->SetTitleSize(0.058); fT->GetYaxis()->SetLabelSize(0.052);
    fT->GetYaxis()->SetTitleOffset(1.30);
    fT->Draw("axis");
    TLine *o1 = new TLine(80.,1.,500.,1.);
    o1->SetLineStyle(2); o1->SetLineColor(kGray+2); o1->Draw();

    TGraphErrors *gb = graphOf(hb,nullptr,cls[j],-0.06);
    TGraphErrors *gd = graphOf(hd,nullptr,cls[j],+0.06);
    styleG(gb,j,markFilledCircle,1.2);
    styleG(gd,j,markOpenSquare,1.2);
    gb->Draw("pz same"); gd->Draw("pz same");

    TLatex tp; tp.SetNDC(); tp.SetTextFont(42);
    tp.SetTextSize(0.068); tp.DrawLatex(0.22,0.86,clsLab[j]);
    if(j == 0){
      TLegend *l2 = new TLegend(0.45,0.60,0.96,0.82);
      l2->SetBorderSize(0); l2->SetFillStyle(0); l2->SetTextSize(0.055);
      l2->AddEntry(gb,"b jets","lp");
      l2->AddEntry(gd,"inclusive jets","lp");
      l2->Draw();
    }

    pBot->cd();
    TH1D *fB2 = new TH1D(Form("fB%d",j),"",1,80.,500.);
    fB2->SetStats(0);
    fB2->GetYaxis()->SetRangeUser(0.,4.6);
    fB2->GetXaxis()->SetTitle("jet p_{T} [GeV]");
    fB2->GetYaxis()->SetTitle("double ratio");
    fB2->GetXaxis()->SetTitleSize(0.110); fB2->GetXaxis()->SetLabelSize(0.098);
    fB2->GetYaxis()->SetTitleSize(0.100); fB2->GetYaxis()->SetLabelSize(0.095);
    fB2->GetXaxis()->SetTitleOffset(1.12);
    fB2->GetYaxis()->SetTitleOffset(0.74);
    fB2->GetYaxis()->SetNdivisions(505);
    fB2->Draw("axis");
    TLine *o2 = new TLine(80.,1.,500.,1.);
    o2->SetLineStyle(2); o2->SetLineColor(kGray+2); o2->Draw();

    TGraphErrors *gr = graphOf(hb,hd,cls[j],0.);
    styleG(gr,j,markFilledCircle,1.1);
    gr->Draw("pz same");
  }

  c2->SaveAs("../../figures/JetsPerZ/bOverInclusiveJetsPerZ_inputs.pdf");
}
