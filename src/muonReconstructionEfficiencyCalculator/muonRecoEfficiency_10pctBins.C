// Muon reco+ID efficiency for the nine 10%-wide centrality classes.
//
// WHY. The analysis uses four hard-coded numbers -- 0.8627 / 0.9069 / 0.9856 /
// 0.9778 for the nominal classes, 0.9708 for pp -- duplicated as literals in
// calculateRAA.C, newFractionCalculation.C, stitchSpectra_caloJets.C and
// calculateBJetsPerZ.cc. Applied to 10% classes they become a STEP function:
// 10-20% and 20-30% both take C2's value although they are 10% apart in
// centrality, and 80-90% takes C4's value although it lies outside the range
// C4 was measured in at all. N_Z must be strictly proportional to T_AA since
// every 10% class spans the same number of events, and with the step mapping it
// is not: 20-30% comes out 1.164 +- 0.035 relative to 40-50%, a 4.7 sigma
// non-closure.
//
// WHAT THIS DOES, AND WHAT IT DOES NOT. It does NOT measure the efficiency in
// fine bins -- that is impossible from the files on disk. The source is
// PYTHIAHYDJET_DiJet_scan_muonReco_mu12_tight.root, which holds FOUR centrality
// classes only (h_inclGenMuonPt_C1..C4, highest index 4), so a fine-binned
// measurement needs a PYTHIAHYDJET muonReco rescan with
// centrality_ultraFineCentBins.h. What this does is interpolate the four
// measured points to the nine class centres, which removes the step
// discontinuity and gives neighbouring 10% classes distinct values. It is an
// interim improvement over a step, not a measurement.
//
// STATISTICS. The four MC points rest on 109-201 generated muons each, so the
// input itself carries 1-3% errors -- C1 is 0.867 +- 0.026. The interpolation
// cannot be better than that, and the quoted errors propagate it.
//
// ABSCISSA. Each measured value is placed at the CENTRALITY CENTRE of the class
// it came from: C1 0-10% -> 5%, C2 10-30% -> 20%, C3 30-50% -> 40%,
// C4 50-80% -> 65%. That assumes the efficiency varies slowly enough across a
// class that its average sits at the centre, which is the same assumption the
// step mapping makes implicitly and much more crudely.
//
// EXTRAPOLATION. 70-80% and 80-90% lie beyond the last measured centre (65%).
// They are held FLAT at the C4 value rather than extrapolated: the efficiency
// is already 0.978 there and bounded above by 1, so a linear extrapolation
// would run out of physical room and invent structure the data cannot support.
// Those two classes are flagged in the output.
//
// Usage: root -l -b -q muonRecoEfficiency_10pctBins.C
// Run from: src/muonReconstructionEfficiencyCalculator/

#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TMath.h"
#include <cstdio>

#include "../../headers/plotting/plotStyle.h"

namespace {

const char *fMC = "/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/"
                  "scanningOutput/PYTHIAHYDJET/final/"
                  "PYTHIAHYDJET_DiJet_scan_muonReco_mu12_tight.root";

const double muPtMin = 15.;      // the analysis muon cut

// measured classes and the centrality centre each value is placed at
const int    NMeas          = 4;
const double measCentre[NMeas] = {5., 20., 40., 65.};
const char  *measLabel [NMeas] = {"0-10%", "10-30%", "30-50%", "50-80%"};

// the nine 10% classes
const int    NCls = 9;
const double clsCentre[NCls] = {5., 15., 25., 35., 45., 55., 65., 75., 85.};
const char  *clsLabel [NCls] = {"0-10%","10-20%","20-30%","30-40%","40-50%",
                                "50-60%","60-70%","70-80%","80-90%"};
const char  *clsTag   [NCls] = {"0to10","10to20","20to30","30to40","40to50",
                                "50to60","60to70","70to80","80to90"};

} // namespace

void muonRecoEfficiency_10pctBins()
{
  initPlotStyle();

  TFile *f = TFile::Open(fMC);
  if(!f || f->IsZombie()){ printf("cannot open %s\n", fMC); return; }

  // ---- the four measured points -------------------------------------------
  double eff[NMeas], err[NMeas];
  printf("\n  measured muon reco+ID efficiency (PYTHIAHYDJET, p_{T} > %.0f GeV)\n", muPtMin);
  printf("    class      tightReco   inclGen      eff        stat\n");
  for(int c = 1; c <= NMeas; c++){
    TH1D *n = (TH1D*) f->Get(Form("h_tightRecoMuonPt_C%d", c));
    TH1D *d = (TH1D*) f->Get(Form("h_inclGenMuonPt_C%d", c));
    if(!n || !d){ printf("    C%d histograms missing\n", c); return; }
    const int i1 = d->FindBin(muPtMin + 0.01);
    const double N = n->Integral(i1, n->GetNbinsX());
    const double D = d->Integral(i1, d->GetNbinsX());
    if(D <= 0.){ printf("    C%d empty denominator\n", c); return; }
    eff[c-1] = N/D;
    // binomial; the MC is unweighted here
    err[c-1] = sqrt(TMath::Max(eff[c-1]*(1.-eff[c-1])/D, 0.));
    printf("    %-7s %10.0f %9.0f   %.4f  +- %.4f\n",
           measLabel[c-1], N, D, eff[c-1], err[c-1]);
  }

  // ---- interpolate to the nine class centres ------------------------------
  double e10[NCls], r10[NCls];
  bool   held[NCls];
  for(int k = 0; k < NCls; k++){
    const double x = clsCentre[k];
    held[k] = false;
    if(x <= measCentre[0]){               // at or below the first centre
      e10[k] = eff[0]; r10[k] = err[0];
    }
    else if(x >= measCentre[NMeas-1]){    // beyond the last centre: hold flat
      e10[k] = eff[NMeas-1]; r10[k] = err[NMeas-1];
      held[k] = (x > measCentre[NMeas-1]);
    }
    else{
      int j = 0;
      while(j < NMeas-2 && x > measCentre[j+1]) j++;
      const double t = (x - measCentre[j])/(measCentre[j+1] - measCentre[j]);
      e10[k] = (1.-t)*eff[j] + t*eff[j+1];
      // the two inputs are independent measurements
      r10[k] = sqrt(pow((1.-t)*err[j],2) + pow(t*err[j+1],2));
    }
  }

  printf("\n  interpolated to the 10%% classes\n");
  printf("    class    centre   efficiency      stat     note\n");
  for(int k = 0; k < NCls; k++)
    printf("    %-7s %5.0f%%   %.4f  +- %.4f   %s\n", clsLabel[k], clsCentre[k],
           e10[k], r10[k], held[k] ? "held flat, beyond last measured centre" : "");

  // ---- header for the analysis macros -------------------------------------
  const char *hdr = "/home/clayton/Analysis/code/bJetRaaAnalysis/headers/"
                    "fitParameters/muonRecoEff_10pctBins.h";
  FILE *out = fopen(hdr, "w");
  if(out){
    fprintf(out, "#pragma once\n");
    fprintf(out, "// Muon reco+ID efficiency per 10%% centrality class.\n");
    fprintf(out, "// GENERATED by src/muonReconstructionEfficiencyCalculator/"
                 "muonRecoEfficiency_10pctBins.C -- do not hand-edit.\n");
    fprintf(out, "//\n");
    fprintf(out, "// Interpolated from the FOUR measured PYTHIAHYDJET classes to the\n");
    fprintf(out, "// nine class centres; it is not a fine-binned measurement. See the\n");
    fprintf(out, "// macro header for what that would require. The last two classes are\n");
    fprintf(out, "// held flat beyond the last measured centre.\n");
    fprintf(out, "const double muonRecoEff10[9] = {");
    for(int k = 0; k < NCls; k++) fprintf(out, "%s%.4f", k?", ":"", e10[k]);
    fprintf(out, "};\n");
    fprintf(out, "const double muonRecoEff10Err[9] = {");
    for(int k = 0; k < NCls; k++) fprintf(out, "%s%.4f", k?", ":"", r10[k]);
    fprintf(out, "};\n");
    fclose(out);
    printf("\n  wrote %s\n", hdr);
  }
  else printf("\n  could not write %s\n", hdr);

  // ---- figure --------------------------------------------------------------
  TCanvas *c = new TCanvas("cEff", "", 800, 650);
  c->SetLeftMargin(0.15); c->SetRightMargin(0.04);
  c->SetTopMargin(0.08);  c->SetBottomMargin(0.13);

  TH1D *fr = new TH1D("frEff", "", 1, 0., 90.);
  fr->SetStats(0);
  fr->GetYaxis()->SetRangeUser(0.78, 1.06);
  fr->GetXaxis()->SetTitle("centrality [%]");
  fr->GetYaxis()->SetTitle("muon reco + ID efficiency");
  fr->GetXaxis()->SetTitleOffset(1.15);
  fr->GetYaxis()->SetTitleOffset(1.35);
  fr->Draw("axis");

  TLine *one = new TLine(0., 1., 90., 1.);
  one->SetLineStyle(2); one->SetLineColor(kGray+2); one->Draw();

  // measured points, with the class width as the x error
  const double halfWidth[NMeas] = {5., 10., 10., 15.};
  TGraphErrors *gM = new TGraphErrors();
  for(int j = 0; j < NMeas; j++){
    gM->SetPoint(j, measCentre[j], eff[j]);
    gM->SetPointError(j, halfWidth[j], err[j]);
  }
  gM->SetLineColor(TColor::GetColor(okabeHex[0]));
  gM->SetMarkerColor(TColor::GetColor(okabeHex[0]));
  gM->SetMarkerStyle(markFilledSquare); gM->SetMarkerSize(1.5); gM->SetLineWidth(2);

  TGraphErrors *gI = new TGraphErrors();
  for(int k = 0; k < NCls; k++){
    gI->SetPoint(k, clsCentre[k], e10[k]);
    gI->SetPointError(k, 5., r10[k]);
  }
  gI->SetLineColor(TColor::GetColor(okabeHex[6]));
  gI->SetMarkerColor(TColor::GetColor(okabeHex[6]));
  gI->SetMarkerStyle(markFilledCircle); gI->SetMarkerSize(1.2); gI->SetLineWidth(2);

  gM->Draw("pz same");
  gI->Draw("pz same");

  TLegend *leg = new TLegend(0.20, 0.17, 0.72, 0.34);
  leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.034);
  leg->AddEntry(gM, "measured, 4 nominal classes (PYTHIA+HYDJET)", "lp");
  leg->AddEntry(gI, "interpolated to 10% classes", "lp");
  leg->Draw();

  TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.034);
  t.DrawLatex(0.15, 0.945, "muon reco + ID efficiency vs centrality");

  c->SaveAs("/home/clayton/Analysis/code/bJetRaaAnalysis/figures/"
            "muonReconstructionEfficiency/muonRecoEff_10pctBins.pdf");
}
