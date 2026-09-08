#include "../../../../headers/functions/divideByBinwidth.h"

TFile *f_noRhoMod, *f_withRhoMod;
TH1D *h_noRhoMod, *h_withRhoMod, *r;
TH1D *h_vz_noRhoMod, *h_vz_withRhoMod;


int centBin = 1;

TCanvas *canv;
TPad *pad_lower, *pad_upper;

TLegend *leg;

void openFiles(){

  f_noRhoMod = TFile::Open("/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/PbPb_noRhoModification_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_sameEventPFClustering_pseudoJetCandPtMin-0.0_2026-8-25_ultraFineCentBins.root");
  f_withRhoMod = TFile::Open("/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/PbPb_withRhoModification_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_sameEventPFClustering_pseudoJetCandPtMin-0.0_2026-8-25_ultraFineCentBins.root");

}

void fetchHistograms(){

  f_noRhoMod->GetObject(Form("h_inclRecoJetPt_C%i",centBin),h_noRhoMod);
  f_withRhoMod->GetObject(Form("h_inclRecoJetPt_C%i",centBin),h_withRhoMod);

  f_noRhoMod->GetObject(Form("h_vz_C%i",centBin),h_vz_noRhoMod);
  f_withRhoMod->GetObject(Form("h_vz_C%i",centBin),h_vz_withRhoMod);

}

void computeRatio(){

  r = (TH1D*) h_withRhoMod->Clone("r");
  r->Divide(h_withRhoMod, h_noRhoMod, 1, 1, "");

}

void normalizeByNEvent(){

  h_noRhoMod->Scale(1./h_vz_noRhoMod->Integral());
  h_withRhoMod->Scale(1./h_vz_withRhoMod->Integral());

}

const int NPtEdges = 8;
double newPtAxis[NPtEdges] = {20,30,40,50,60,80,100,150};

void rebin(){

  h_noRhoMod = (TH1D*) h_noRhoMod->Rebin(NPtEdges-1,"h_noRhoMod",newPtAxis);
  h_withRhoMod = (TH1D*) h_withRhoMod->Rebin(NPtEdges-1,"h_withRhoMod",newPtAxis);

  divideByBinwidth(h_noRhoMod);
  divideByBinwidth(h_withRhoMod);
  

}

void draw(){

  canv = new TCanvas("canv","canv",700,700);
  canv->cd();
  pad_upper = new TPad("pad_upper","pad_upper",0,0.4,1,1);
  pad_lower = new TPad("pad_lower","pad_lower",0,0,1,0.4);
  pad_upper->SetBottomMargin(0.);
  pad_lower->SetBottomMargin(0.2);
  pad_upper->SetTopMargin(0.1);
  pad_lower->SetTopMargin(0.);
  pad_upper->SetLeftMargin(0.15);
  pad_lower->SetLeftMargin(0.15);
  pad_upper->Draw();
  pad_lower->Draw();
  pad_upper->cd();
  h_noRhoMod->SetTitle("");
  h_noRhoMod->SetStats(0);
  h_noRhoMod->GetYaxis()->SetTitleSize(0.06);
  h_noRhoMod->GetYaxis()->SetLabelSize(0.05);
  h_noRhoMod->GetYaxis()->SetTitle("1/N^{evt} dN^{jet}/dp_{T} [GeV^{-1}]");
  h_noRhoMod->Draw();
  h_withRhoMod->Draw("same");
  pad_lower->cd();
  r->SetTitle("");
  r->SetStats(0);
  r->GetYaxis()->SetTitleSize(0.07);
  r->GetYaxis()->SetLabelSize(0.06);
  r->GetYaxis()->SetTitle("Ratio");
  r->GetYaxis()->SetRangeUser(0,1);
  r->GetXaxis()->SetTitleSize(0.07);
  r->GetXaxis()->SetLabelSize(0.06);
  r->GetXaxis()->SetTitle("Jet p_{T} [GeV]");
  r->Draw();
  
  

}

int main(){

  centBin = 2;
  
  openFiles();
  fetchHistograms();
  normalizeByNEvent();
  rebin();
  computeRatio();
  draw();

  return -1;
}
