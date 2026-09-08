#include "../../../headers/functions/divideByBinwidth.h"



// set f1 = sameEvent, f2 = mixedEvent
TFile *f1, *f2;
// dR vs fastJetPt
TH2D *D1, *D2;
TH1D *d1, *d2;

int centBin;

TCanvas *canv, *canv_dR;
TPad *pad, *pad_dR;
TLegend *leg;

const int Nedge_dR = 12;
double newAxis_dR[Nedge_dR] = {0,0.02,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5};

double jetPtLow, jetPtHigh;
TLatex *la;

void openFiles(){
  f1 = TFile::Open("/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/PbPb_SingleMuon_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_mixedEventPFClustering_pseudoJetCandPtMin-0.0_2026-8-28_ultraFineCentBins_partial.root");
  f2 = TFile::Open("/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/PbPb_SingleMuon_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_mixedEventPFClustering_pseudoJetCandPtMin-0.0_2026-8-27_ultraFineCentBins_partial.root");
}

void loadHistograms(){

  string histogramString_dR_1 = Form("h_fastJetMuonDR_inclusiveClosestFastJet_C%i",centBin);
  string histogramString_dR_2 = Form("h_fastJetMuonDR_fastJetPt_PF_bkgSub_RC_C%i",centBin);
  f1->GetObject(histogramString_dR_1.c_str(),D1);
  f2->GetObject(histogramString_dR_2.c_str(),D2);

}

void projectHistograms(double jetPtLow, double jetPtHigh){

  int binLow = D1->ProjectionY()->FindBin(jetPtLow);
  int binHigh = D1->ProjectionY()->FindBin(jetPtHigh);
  d1 = (TH1D*) D1->ProjectionX("d1",binLow,binHigh);
  d2 = (TH1D*) D2->ProjectionX("d2",binLow,binHigh);

}

void normalizeByIntegral(){

  d1->Scale(1./d1->Integral());
  d2->Scale(1./d2->Integral());

}

void stylizeHistograms(){

  double lineWidth = 2;
  d1->SetLineWidth(lineWidth);
  d2->SetLineWidth(lineWidth);

  d1->SetLineColor(kBlue-4);
  d2->SetLineColor(kRed-4);
  
  d1->SetStats(0);
  d2->SetStats(0);

  d1->SetTitle("");
  d2->SetTitle("");

  string xTitle_dR = "#it{#Delta}#it{R}(#it{#mu},jet)";
  string yTitle_dR = "1/#it{N}^{#it{#mu}+jet}_{tot} d#it{N}^{#it{#mu}+jet}/d(#it{#Delta}#it{R})";
  
  d1->GetXaxis()->SetTitle(xTitle_dR.c_str());
  d2->GetXaxis()->SetTitle(xTitle_dR.c_str());

  d1->GetYaxis()->SetTitle(yTitle_dR.c_str());
  d2->GetYaxis()->SetTitle(yTitle_dR.c_str());

  d1->GetXaxis()->SetTitleSize(0.04);
  d2->GetXaxis()->SetTitleSize(0.04);
  
  d1->GetYaxis()->SetLabelSize(0.04);
  d2->GetYaxis()->SetLabelSize(0.04);

}

void rebinHistograms(){

  d1 = (TH1D*) d1->Rebin(Nedge_dR-1,"d1",newAxis_dR);
  d2 = (TH1D*) d2->Rebin(Nedge_dR-1,"d2",newAxis_dR);

  divideByBinwidth(d1);
  divideByBinwidth(d2);
  
}

void drawHistograms(){

  canv_dR = new TCanvas("canv_dR","canv_dR",700,700);
  canv_dR->cd();
  pad_dR = new TPad("pad_dR","pad_dR",0,0,1,1);
  pad_dR->SetLeftMargin(0.12);
  pad_dR->Draw();
  pad_dR->cd();
  d2->Draw();
  d1->Draw("same");
  leg = new TLegend(0.65,0.75,0.85,0.85);
  leg->SetBorderSize(0);
  leg->AddEntry(d1,"all muons");
  leg->AddEntry(d2,"in-jet muons");
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
  string candInfoString = "Mixed-event PF Candidates";
  la->DrawLatexNDC(0.17,0.75,candInfoString.c_str());
  la->SetTextSize(0.035);
  string titleString = "PbPb SingleMuon (5.02 TeV)";
  la->DrawLatexNDC(0.13,0.92,titleString.c_str());
  
  
  

}

string canvasPath = "";

void saveCanvas(){

  canvasPath = Form("../../../figures/muJetDr/fastJetMuonDR_inJetVsInclusive_C%i.pdf",centBin);
  canv_dR->SaveAs(canvasPath.c_str());

}

int main(){

  centBin = 0;
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
  
  std::cout << "Normalize histograms by integral...";
  normalizeByIntegral();
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

  std::cout << "Saving canvas...";
  saveCanvas();
  std::cout << "done.\n";
  

  

  return -1;
}
