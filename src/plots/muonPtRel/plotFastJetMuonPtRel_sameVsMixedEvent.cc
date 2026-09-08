#include "../../../headers/functions/divideByBinwidth.h"



// set f1 = sameEvent, f2 = mixedEvent
TFile *f1, *f2;
// fastJetMuonPtRel vs fastJetPt
TH2D *H1, *H2;
TH1D *h1, *h2;
// dR vs fastJetPt
TH2D *D1, *D2;
TH1D *d1, *d2;
TH1D *v1, *v2;

int centBin;

TCanvas *canv, *canv_dR;
TPad *pad, *pad_dR;
TLegend *leg;

const int Nedge = 12;
double newAxis[Nedge] = {0,0.2,0.4,0.6,0.8,1,1.4,1.8,2.2,2.8,3.8,5};

// const int Nedge = 6;
// double newAxis[Nedge] = {0,1,2,3,4,5};

const int Nedge_dR = 12;
double newAxis_dR[Nedge_dR] = {0,0.02,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5};

double jetPtLow, jetPtHigh;
TLatex *la;

void openFiles(){
  f1 = TFile::Open("/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/PbPb_SingleMuon_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_sameEventPFClustering_pseudoJetCandPtMin-0.0_2026-8-28_ultraFineCentBins.root");
  f2 = TFile::Open("/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/PbPb_SingleMuon_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_mixedEventPFClustering_pseudoJetCandPtMin-0.0_2026-8-27_ultraFineCentBins_partial.root");
}

void loadHistograms(){

  string histogramString = Form("h_fastJetMuonPtRel_fastJetPt_PF_bkgSub_RC_C%i",centBin);
  f1->GetObject(histogramString.c_str(),H1);
  f2->GetObject(histogramString.c_str(),H2);

  string histogramString_dR = Form("h_fastJetMuonDR_fastJetPt_PF_bkgSub_RC_C%i",centBin);
  f1->GetObject(histogramString_dR.c_str(),D1);
  f2->GetObject(histogramString_dR.c_str(),D2);

  string histogramString_vz = Form("h_vz_C%i",centBin);
  f1->GetObject(histogramString_vz.c_str(),v1);
  f2->GetObject(histogramString_vz.c_str(),v2);

}

void projectHistograms(double jetPtLow, double jetPtHigh){

  int binLow = H1->ProjectionY()->FindBin(jetPtLow);
  int binHigh = H1->ProjectionY()->FindBin(jetPtHigh);
  h1 = (TH1D*) H1->ProjectionX("h1",binLow,binHigh);
  h2 = (TH1D*) H2->ProjectionX("h2",binLow,binHigh);
  d1 = (TH1D*) D1->ProjectionX("d1",binLow,binHigh);
  d2 = (TH1D*) D2->ProjectionX("d2",binLow,binHigh);
  
}

void normalizeByIntegral(){

  h1->Scale(1./h1->Integral());
  h2->Scale(1./h2->Integral());
  d1->Scale(1./d1->Integral());
  d2->Scale(1./d2->Integral());

}

void normalizeByNEvent(){
  
  h1->Scale(1./v1->Integral());
  h2->Scale(1./v2->Integral());
  d1->Scale(1./v1->Integral());
  d2->Scale(1./v2->Integral());
  
}

void stylizeHistograms(){

  double lineWidth = 2;
  h1->SetLineWidth(lineWidth);
  h2->SetLineWidth(lineWidth);
  d1->SetLineWidth(lineWidth);
  d2->SetLineWidth(lineWidth);

  h1->SetLineColor(kBlue-4);
  h2->SetLineColor(kRed-4);
  d1->SetLineColor(kBlue-4);
  d2->SetLineColor(kRed-4);
  

  h1->SetStats(0);
  h2->SetStats(0);
  d1->SetStats(0);
  d2->SetStats(0);

  h1->SetTitle("");
  h2->SetTitle("");
  d1->SetTitle("");
  d2->SetTitle("");

  string xTitle = "Muon #it{p}_{T}^{rel} [GeV]";
  string xTitle_dR = "#it{#Delta}#it{R}(#it{#mu},jet)";
  string yTitle = "1/#it{N}^{#it{#mu}+jet}_{tot} d#it{N}^{#it{#mu}+jet}/d#it{p}_{T}^{rel} [GeV^{-1}]";
  string yTitle_dR = "1/#it{N}^{#it{#mu}+jet}_{tot} d#it{N}^{#it{#mu}+jet}/d(#it{#Delta}#it{R})";
  

  h1->GetXaxis()->SetTitle(xTitle.c_str());
  h2->GetXaxis()->SetTitle(xTitle.c_str());
  d1->GetXaxis()->SetTitle(xTitle_dR.c_str());
  d2->GetXaxis()->SetTitle(xTitle_dR.c_str());

  h1->GetYaxis()->SetTitle(yTitle.c_str());
  h2->GetYaxis()->SetTitle(yTitle.c_str());
  d1->GetYaxis()->SetTitle(yTitle_dR.c_str());
  d2->GetYaxis()->SetTitle(yTitle_dR.c_str());

  h1->GetXaxis()->SetTitleSize(0.04);
  h2->GetXaxis()->SetTitleSize(0.03);
  d1->GetXaxis()->SetTitleSize(0.04);
  d2->GetXaxis()->SetTitleSize(0.04);
  
  h1->GetYaxis()->SetLabelSize(0.04);
  h2->GetYaxis()->SetLabelSize(0.03);
  d1->GetYaxis()->SetLabelSize(0.04);
  d2->GetYaxis()->SetLabelSize(0.04);

  


}

void rebinHistograms(){

  h1 = (TH1D*) h1->Rebin(Nedge-1,"h1",newAxis);
  h2 = (TH1D*) h2->Rebin(Nedge-1,"h2",newAxis);
  d1 = (TH1D*) d1->Rebin(Nedge_dR-1,"d1",newAxis_dR);
  d2 = (TH1D*) d2->Rebin(Nedge_dR-1,"d2",newAxis_dR);

  divideByBinwidth(h1);
  divideByBinwidth(h2);
  divideByBinwidth(d1);
  divideByBinwidth(d2);
  
}

void drawHistograms(){

  canv = new TCanvas("canv","canv",700,700);
  canv->cd();
  pad = new TPad("pad","pad",0,0,1,1);
  pad->SetLeftMargin(0.12);
  pad->Draw();
  pad->cd();
  //h2->GetYaxis()->SetRangeUser(0,0.69);
  h1->Draw();
  h2->Draw("same");
  leg = new TLegend(0.65,0.75,0.85,0.85);
  leg->SetBorderSize(0);
  leg->AddEntry(h1,"same event");
  leg->AddEntry(h2,"mixed event");
  leg->SetTextSize(0.035);
  leg->Draw();
  la = new TLatex();
  la->SetTextFont(42);
  la->SetTextSize(0.03);
  string jetPtString = "";
  if(jetPtLow < 100.) jetPtString = Form("%2.0f < #it{p}_{T}^{jet} (GeV) < %3.0f, |#it{#eta}^{jet}| < 1.6",jetPtLow,jetPtHigh);
  else jetPtString = Form("%3.0f < #it{p}_{T}^{jet} (GeV) < %3.0f, |#it{#eta}^{jet}| < 1.6",jetPtLow,jetPtHigh);
  la->DrawLatexNDC(0.17,0.85,jetPtString.c_str());
  string muonPtString = "#it{p}_{T}^{#it{#mu}} > 15 GeV, |#it{#eta}^{#it{#mu}}| < 2";
  la->DrawLatexNDC(0.17,0.80,muonPtString.c_str());
  la->SetTextSize(0.035);
  string titleString = "PbPb SingleMuon (5.02 TeV)";
  la->DrawLatexNDC(0.13,0.92,titleString.c_str());

  ///////////////////////////////////////////////////////////////////////////////////////

  canv_dR = new TCanvas("canv_dR","canv_dR",700,700);
  canv_dR->cd();
  pad_dR = new TPad("pad_dR","pad_dR",0,0,1,1);
  pad_dR->SetLeftMargin(0.12);
  pad_dR->Draw();
  pad_dR->cd();
  d2->Draw();
  d1->Draw("same");
  

}

TFile *file_histograms = nullptr;

void saveHistograms(){

  file_histograms = (TFile*) TFile::Open(Form("../../../rootFiles/fakeMuonPtRel/fakeMuonPtRel_C%i.root",centBin),"recreate");
  h2->Write(Form("h_fakeMuonPtRel_C%i",centBin));

  
}

int main(){

  centBin = 1;
  jetPtLow = 60.;
  jetPtHigh = 80.;

  std::cout << "\n\n===== Plot FastJet Muon PtRel, Same Vs. Mixed Event =====\n\n";

  std::cout << "Opening files...";
  openFiles();
  std::cout << "done.\n";

  std::cout << "Loading histograms...";
  loadHistograms();
  std::cout << "done.\n";
  
  std::cout << "Projecting histogram from " << jetPtLow << " < pT < " << jetPtHigh << "...";
  projectHistograms(jetPtLow,jetPtHigh);
  std::cout << "done.\n";
  
  std::cout << "Normalize histograms by NEvent...";
  //normalizeByIntegral();
  normalizeByNEvent();
  std::cout << "done.\n";

  std::cout << "Rebinning histograms...";
  rebinHistograms();
  std::cout << "done.\n";

  std::cout << "Stylizing histograms...";
  stylizeHistograms();
  std::cout << "done.\n";
  
  std::cout << "Drawing histograms...";
  drawHistograms();
  std::cout << "done.\n";

  std::cout << "Saving histograms...";
  saveHistograms();
  std::cout << "done.\n";
  

  return -1;
}
