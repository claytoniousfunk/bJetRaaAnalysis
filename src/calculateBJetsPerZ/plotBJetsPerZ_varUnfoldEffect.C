// What the variable-width b-jet unfolding does to the PbPb/pp b-jets-per-Z
// ratio, one panel per centrality class.
//
// Inputs are two runs of calculateBJetsPerZ.cc differing only in
// doBJetUnfolding, saved side by side:
//   histograms_BJetPbPbToPP_noUnfold.root
//   histograms_BJetPbPbToPP_varUnfold_1iter.root   (1 Bayes iteration)
//
// The unfolded spectrum is built with h_matchedRecoJetPt_genJetPt_var_bJets on
// the response scan's variable axis, applied to the b-jet spectrum after the
// b-purity multiplication. See the block above unfoldBJetSpectrum() in
// calculateBJetsPerZ.cc for the reco/truth axis handling.
//
// C3's 300-500 bin is DROPPED: its measured content is 2.5e-16, i.e. the bin is
// empty, so the unfolded 0.226 there is the prior spilling down from the matrix
// rather than anything the data says. Plotting it would show a real
// measurement where there is none.

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

const char *fNo  = "rootFiles/BJetPbPbToPP/Data/histograms_BJetPbPbToPP_noUnfold.root";
const char *fUn  = "rootFiles/BJetPbPbToPP/Data/histograms_BJetPbPbToPP_varUnfold_1iter.root";

const char *cls[4]   = {"C1","C2","C3","C4"};
const char *clsLab[4]= {"PbPb 0-10%","PbPb 10-30%","PbPb 30-50%","PbPb 50-80%"};

// Turn a ratio histogram into a graph, skipping bins the measurement does not
// support. Points sit at the bin centre with the bin half-width as the x error.
TGraphErrors *toGraph(TH1D *h, const char *cl)
{
  TGraphErrors *g = new TGraphErrors();
  int n = 0;
  for(int i = 1; i <= h->GetNbinsX(); i++){
    if(h->GetBinContent(i) == 0.) continue;               // never filled
    const double lo = h->GetBinLowEdge(i), hi = h->GetBinLowEdge(i+1);
    if(TString(cl) == "C3" && lo > 299.) continue;        // see header note
    g->SetPoint(n,0.5*(lo+hi),h->GetBinContent(i));
    g->SetPointError(n,0.5*(hi-lo),h->GetBinError(i));
    n++;
  }
  return g;
}

} // namespace

void plotBJetsPerZ_varUnfoldEffect()
{
  initPlotStyle();

  TFile *f1 = TFile::Open(fNo);
  TFile *f2 = TFile::Open(fUn);
  if(!f1 || f1->IsZombie() || !f2 || f2->IsZombie()){
    printf("cannot open the input files -- run calculateBJetsPerZ.cc twice,\n"
           "once with doBJetUnfolding=false, and save them under the names at\n"
           "the top of this macro.\n");
    return;
  }

  TCanvas *c = new TCanvas("cUnf","",1000,900);
  c->Divide(2,2,0.001,0.001);

  for(int j = 0; j < 4; j++){

    TH1D *hNo = (TH1D*) f1->Get(Form("bFraction_%s",cls[j]));
    TH1D *hUn = (TH1D*) f2->Get(Form("bFraction_%s",cls[j]));
    if(!hNo || !hUn){ printf("missing bFraction_%s\n",cls[j]); continue; }

    TPad *p = (TPad*) c->cd(j+1);
    p->SetLeftMargin(0.14); p->SetRightMargin(0.03);
    p->SetTopMargin(0.07);  p->SetBottomMargin(0.13);

    TH1D *fr = new TH1D(Form("fr_%s",cls[j]),"",1,80.,500.);
    fr->SetStats(0);
    fr->GetYaxis()->SetRangeUser(0.,2.6);
    fr->GetXaxis()->SetTitle("jet p_{T} [GeV]");
    fr->GetYaxis()->SetTitle("PbPb / pp  b jets per Z");
    fr->GetXaxis()->SetTitleSize(0.050); fr->GetYaxis()->SetTitleSize(0.050);
    fr->GetXaxis()->SetLabelSize(0.045); fr->GetYaxis()->SetLabelSize(0.045);
    fr->GetXaxis()->SetTitleOffset(1.10);
    fr->GetYaxis()->SetTitleOffset(1.30);
    fr->Draw("axis");

    TLine *one = new TLine(80.,1.,500.,1.);
    one->SetLineStyle(2); one->SetLineColor(kGray+2); one->Draw();

    TGraphErrors *gNo = toGraph(hNo,cls[j]);
    TGraphErrors *gUn = toGraph(hUn,cls[j]);

    const int cNo = TColor::GetColor(okabeHex[6]);   // vermillion
    const int cUn = TColor::GetColor(okabeHex[0]);   // black

    gNo->SetLineColor(cNo); gNo->SetMarkerColor(cNo);
    gNo->SetMarkerStyle(markOpenSquare); gNo->SetMarkerSize(1.3); gNo->SetLineWidth(2);
    gUn->SetLineColor(cUn); gUn->SetMarkerColor(cUn);
    gUn->SetMarkerStyle(markFilledCircle); gUn->SetMarkerSize(1.3); gUn->SetLineWidth(2);

    gNo->Draw("pz same");
    gUn->Draw("pz same");

    TLatex t; t.SetNDC(); t.SetTextFont(42);
    t.SetTextSize(0.055); t.DrawLatex(0.18,0.87,clsLab[j]);

    if(j == 0){
      TLegend *leg = new TLegend(0.40,0.72,0.95,0.89);
      leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.045);
      leg->AddEntry(gNo,"not unfolded","lp");
      leg->AddEntry(gUn,"unfolded, 1 iteration","lp");
      leg->Draw();
    }
  }

  c->SaveAs("../../figures/JetsPerZ/bJetsPerZ_varUnfoldEffect.pdf");
}
