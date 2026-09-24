// HLT_HIL3Mu12 efficiency vs muon eta, pp.
//
// Motivation: the integrated pp plateau from the low-statistics ZeroBias/MinBias
// probe sample (a = 0.931 +- 0.006, ~1.9k probes) sits BELOW peripheral PbPb
// (0.970), which is backwards -- peripheral PbPb is the pp-like occupancy limit.
// The high-statistics HighEGJet+Jet60 probe sample gives 0.954 from 218k probes.
// This figure shows where the two differ in eta, so it is clear whether the
// ZeroBias deficit is localized (a real acceptance effect) or a fluctuation.
//
// ZeroBias: h_muTrigEff_pass/all (TH2 pT x eta) from pp_scan.C, prescale-1
// events only, unweighted -> genuine binomial, Clopper-Pearson intervals exact.
// Restricted here to pT > 15 GeV (the analysis muon cut, and well onto the
// plateau).
//
// HighEGJet: the legacy triggerEffScan muEta_trigOn/muEta_all (TH1, eta only).
// Two caveats, both noted on the figure: it integrates over ALL probe pT (no pT
// axis to cut on, ~31k of its 249k probes are below 15 GeV), and its numerator
// was filled prescale-WEIGHTED while the denominator was not. It is drawn as a
// statistical reference for the shape, not as the quantity to subtract.

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TMath.h"
#include "TEfficiency.h"
#include <cstdio>

#include "../../../headers/plotting/plotStyle.h"

namespace {

const char *fZB =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/pp/"
  "pp_MinBias_caloJets_manualJEC_mu12_pTmu-15to999_tight_deltaR-40_"
  "jetTrkMaxFilter_WDecayFilter_2026-9-23.root";

// read-only sibling repo; legacy scan generation
const char *fEG =
  "/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/scanningOutput/"
  "pp/latest/triggerEff/"
  "pp_HighEGJet_triggerEffScan_Jet60HLT_mu12_tight_onlyOneMuonPerEvent.root";

// Signed eta, with the muon-system boundaries as bin edges: barrel |eta|<0.9,
// barrel/endcap overlap 0.9-1.2, endcap beyond. Coarse enough that the ~1.9k
// ZeroBias probes give a meaningful error per bin.
const double etaEdge[] = {-2.4,-2.1,-1.6,-1.2,-0.9,-0.45,0.,0.45,0.9,1.2,1.6,2.1,2.4};
const int    nEta      = 12;

// Build the efficiency graph from summed pass/all counts per eta bin.
//
// TGraphAsymmErrors::Divide is NOT used here: these histograms are filled with
// Integral() sums (and, for the legacy sample, prescale weights), so its
// binomial consistency check degrades the interval to a near-meaningless one --
// it returned +0.04/-0.24 on a bin holding 140k probes. Clopper-Pearson on the
// (effective) counts is computed directly instead.
TGraphAsymmErrors *effGraph(const double *pass, const double *all)
{
  TGraphAsymmErrors *g = new TGraphAsymmErrors();
  int n = 0;
  for(int i = 0; i < nEta; i++){
    if(all[i] <= 0.) continue;
    const double lo = etaEdge[i], hi = etaEdge[i+1];
    const double eff = pass[i]/all[i];
    // round to integers: Clopper-Pearson is defined on counts
    const double N = TMath::Nint(all[i]), k = TMath::Nint(pass[i]);
    const double up = TEfficiency::ClopperPearson(N,k,0.683,true)  - eff;
    const double dn = eff - TEfficiency::ClopperPearson(N,k,0.683,false);
    g->SetPoint(n,0.5*(lo+hi),eff);
    g->SetPointError(n,0.5*(hi-lo),0.5*(hi-lo),dn,up);
    n++;
  }
  return g;
}

// Rebin a (pT x eta) pass/all pair onto etaEdge for pT > ptMin. Counts are
// unweighted here (w = 1 for this data scan), so the interval is exact.
TGraphAsymmErrors *effFrom2D(TH2D *P, TH2D *A, double ptMin, double *nOut)
{
  double pass[nEta] = {0}, all[nEta] = {0};
  const int x1 = A->GetXaxis()->FindBin(ptMin), x2 = A->GetNbinsX();
  double nTot = 0.;
  for(int iy = 1; iy <= A->GetNbinsY(); iy++){
    const double eta = A->GetYaxis()->GetBinCenter(iy);
    const int b = TMath::BinarySearch(nEta+1,etaEdge,eta);
    if(b < 0 || b >= nEta) continue;
    const double p = P->Integral(x1,x2,iy,iy);
    const double a = A->Integral(x1,x2,iy,iy);
    pass[b] += p; all[b] += a; nTot += a;
  }
  if(nOut) *nOut = nTot;
  return effGraph(pass,all);
}

// Same, from the legacy 1D eta histograms (no pT axis available).
TGraphAsymmErrors *effFrom1D(TH1D *P, TH1D *A, double *nOut)
{
  double pass[nEta] = {0}, all[nEta] = {0};
  double nTot = 0.;
  for(int i = 1; i <= A->GetNbinsX(); i++){
    const double eta = A->GetBinCenter(i);
    const int b = TMath::BinarySearch(nEta+1,etaEdge,eta);
    if(b < 0 || b >= nEta) continue;
    pass[b] += P->GetBinContent(i); all[b] += A->GetBinContent(i);
    nTot += A->GetBinContent(i);
  }
  if(nOut) *nOut = nTot;
  return effGraph(pass,all);
}

void printTable(const char *name, TGraphAsymmErrors *g)
{
  printf("\n%s\n", name);
  printf("   eta range        efficiency\n");
  for(int i = 0; i < g->GetN(); i++){
    double x, y; g->GetPoint(i,x,y);
    printf("  %+5.2f to %+5.2f   %.4f  +%.4f -%.4f\n",
           x - g->GetErrorXlow(i), x + g->GetErrorXhigh(i), y,
           g->GetErrorYhigh(i), g->GetErrorYlow(i));
  }
}

} // namespace

void plotMuonTriggerEfficiency_vsEta()
{
  initPlotStyle();

  TFile *fz = TFile::Open(fZB);
  if(!fz || fz->IsZombie()){ printf("cannot open %s\n", fZB); return; }
  TH2D *Pz = (TH2D*)fz->Get("h_muTrigEff_pass");
  TH2D *Az = (TH2D*)fz->Get("h_muTrigEff_all");
  if(!Pz || !Az){ printf("missing h_muTrigEff_pass/all in the pp scan\n"); return; }

  double nZB = 0.;
  TGraphAsymmErrors *gZB = effFrom2D(Pz,Az,15.,&nZB);

  TGraphAsymmErrors *gEG = nullptr;
  double nEG = 0.;
  TFile *fe = TFile::Open(fEG);
  if(fe && !fe->IsZombie()){
    TH1D *Pe = (TH1D*)fe->Get("muEta_trigOn_C0");
    TH1D *Ae = (TH1D*)fe->Get("muEta_all_C0");
    if(Pe && Ae) gEG = effFrom1D(Pe,Ae,&nEG);
  }

  printTable(Form("pp ZeroBias/MinBias 2026-9-23, pT > 15 GeV  (N = %.0f)", nZB), gZB);
  if(gEG) printTable(Form("pp HighEGJet + Jet60, all probe pT  (N = %.0f)", nEG), gEG);

  TCanvas *c = new TCanvas("cEta","",700,600);
  c->SetLeftMargin(0.13); c->SetRightMargin(0.04);
  c->SetTopMargin(0.07);  c->SetBottomMargin(0.13);

  TH1D *fr = new TH1D("frEta","",nEta,etaEdge);
  fr->SetStats(0);
  fr->GetYaxis()->SetRangeUser(0.60,1.10);
  fr->GetXaxis()->SetTitle("muon #eta");
  fr->GetYaxis()->SetTitle("HLT_HIL3Mu12 efficiency");
  fr->GetXaxis()->SetTitleOffset(1.15);
  fr->GetYaxis()->SetTitleOffset(1.35);
  fr->Draw("axis");

  TLine one(-2.4,1.0,2.4,1.0);
  one.SetLineStyle(2); one.SetLineColor(kGray+2); one.Draw();

  // barrel / endcap boundary, both signs
  for(int s = -1; s <= 1; s += 2){
    TLine *b = new TLine(s*1.2,0.60,s*1.2,1.10);
    b->SetLineStyle(3); b->SetLineColor(kGray+1); b->Draw();
  }

  if(gEG){
    gEG->SetLineColor(TColor::GetColor(hexMC));
    gEG->SetMarkerColor(TColor::GetColor(hexMC));
    gEG->SetMarkerStyle(markOpenSquare);
    gEG->SetMarkerSize(1.1); gEG->SetLineWidth(2);
    gEG->Draw("pz same");
  }
  gZB->SetLineColor(TColor::GetColor(hexData));
  gZB->SetMarkerColor(TColor::GetColor(hexData));
  gZB->SetMarkerStyle(markFilledCircle);
  gZB->SetMarkerSize(1.2); gZB->SetLineWidth(2);
  gZB->Draw("pz same");

  TLegend *leg = new TLegend(0.16,0.16,0.62,0.31);
  leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.033);
  leg->AddEntry(gZB,Form("ZeroBias, p_{T}^{#mu} > 15 GeV  (N = %.0f)",nZB),"lp");
  if(gEG) leg->AddEntry(gEG,Form("HighEGJet + Jet60, all p_{T}^{#mu}  (N = %.0f)",nEG),"lp");
  leg->Draw();

  TLatex t; t.SetNDC(); t.SetTextFont(42);
  t.SetTextSize(0.040); t.DrawLatex(0.13,0.945,"pp 5.02 TeV, tight probe muons");
  t.SetTextSize(0.030);
  t.DrawLatex(0.60,0.875,"dotted: barrel/endcap");

  c->SaveAs("/home/clayton/Analysis/code/bJetRaaAnalysis/figures/triggerEfficiency/"
            "triggerEfficiency_mu12_pp_vsEta.pdf");
}
