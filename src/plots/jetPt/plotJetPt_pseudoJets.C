#include "../../../headers/functions/divideByBinwidth.h"

double integrate(TH1D *h, double lowerBound, double upperBound){

  double result = 0.0;
  double smallShift = 0.01;
  
  int a = h->FindBin(lowerBound+smallShift);
  int b = h->FindBin(upperBound-smallShift);

  for(int i = a; i <= b; i++){
    result += h->GetBinContent(i);
  }

  return result;

}

TString getCentString(bool doC16 = false,
		      bool doC15 = false,
		      bool doC14 = false,
		      bool doC13 = false,
		      bool doC12 = false,
		      bool doC11 = false,
		      bool doC10 = false,
		      bool doC9 = false,
		      bool doC8 = false,
		      bool doC7 = false,
		      bool doC6 = false,
		      bool doC5 = false,   
		      bool doC4 = true,
		      bool doC3 = false,
		      bool doC2 = false,
		      bool doC1 = false){

  TString result = "";

  if(doC16) result = "C16";
  else if(doC15) result = "C15";
  else if(doC14) result = "C14";
  else if(doC13) result = "C13";
  else if(doC12) result = "C12";
  else if(doC11) result = "C11";
  else if(doC10) result = "C10";
  else if(doC9) result = "C9";
  else if(doC8) result = "C8";
  else if(doC7) result = "C7";
  else if(doC6) result = "C6";
  else if(doC5) result = "C5";
  else if(doC4) result = "C4";
  else if(doC3) result = "C3";
  else if(doC2) result = "C2";
  else if(doC1) result = "C1";
  else{};

  return result;

}

double getCentHiBinLowerBound(bool doC16 = false,
			   bool doC15 = false,
			   bool doC14 = false,
			   bool doC13 = false,
			   bool doC12 = false,
			   bool doC11 = false,
			   bool doC10 = false,
			   bool doC9 = false,
			   bool doC8 = false,
			   bool doC7 = false,
			   bool doC6 = false,
			   bool doC5 = false,   
			   bool doC4 = true,
			   bool doC3 = false,
			   bool doC2 = false,
			   bool doC1 = false){

  double result = 0.;

  if(doC16) result = 150.;
  else if(doC15) result = 140.;
  else if(doC14) result = 130.;
  else if(doC13) result = 120.;
  else if(doC12) result = 110.;
  else if(doC11) result = 100.;
  else if(doC10) result = 90.;
  else if(doC9) result = 80.;
  else if(doC8) result = 70.;
  else if(doC7) result = 60.;
  else if(doC6) result = 50.;
  else if(doC5) result = 40.;
  else if(doC4) result = 30.;
  else if(doC3) result = 20.;
  else if(doC2) result = 10.;
  else if(doC1) result = 0.;
  else{};
  
  return result;

}

double getCentHiBinUpperBound(bool doC16 = false,
			   bool doC15 = false,
			   bool doC14 = false,
			   bool doC13 = false,
			   bool doC12 = false,
			   bool doC11 = false,
			   bool doC10 = false,
			   bool doC9 = false,
			   bool doC8 = false,
			   bool doC7 = false,
			   bool doC6 = false,
			   bool doC5 = false,   
			   bool doC4 = true,
			   bool doC3 = false,
			   bool doC2 = false,
			   bool doC1 = false){

  double result = 0.;

  if(doC16) result = 160.;
  else if(doC15) result = 150.;
  else if(doC14) result = 140.;
  else if(doC13) result = 130.;
  else if(doC12) result = 120.;
  else if(doC11) result = 110.;
  else if(doC10) result = 100.;
  else if(doC9) result = 90.;
  else if(doC8) result = 80.;
  else if(doC7) result = 70.;
  else if(doC6) result = 60.;
  else if(doC5) result = 50.;
  else if(doC4) result = 40.;
  else if(doC3) result = 30.;
  else if(doC2) result = 20.;
  else if(doC1) result = 10.;
  else{};
  
  return result;

}

TString getTitle(bool doC16 = false,
		 bool doC15 = false,
		 bool doC14 = false,
		 bool doC13 = false,
		 bool doC12 = false,
		 bool doC11 = false,
		 bool doC10 = false,
		 bool doC9 = false,
		 bool doC8 = false,
		 bool doC7 = false,
		 bool doC6 = false,
		 bool doC5 = false,   
		 bool doC4 = true,
		 bool doC3 = false,
		 bool doC2 = false,
		 bool doC1 = false){

  TString result = "";

  if(doC16) result = "PbPb 75-80%";
  else if(doC15) result = "PbPb 70-75%";
  else if(doC14) result = "PbPb 65-70%";
  else if(doC13) result = "PbPb 60-65%";
  else if(doC12) result = "PbPb 55-60%";
  else if(doC11) result = "PbPb 50-55%";
  else if(doC10) result = "PbPb 45-50%";
  else if(doC9) result = "PbPb 40-45%";
  else if(doC8) result = "PbPb 35-40%";
  else if(doC7) result = "PbPb 30-35%";
  else if(doC6) result = "PbPb 25-30%";
  else if(doC5) result = "PbPb 20-25%";
  else if(doC4) result = "PbPb 15-20%";
  else if(doC3) result = "PbPb 10-15%";
  else if(doC2) result = "PbPb 5-10%";
  else if(doC1) result = "PbPb 0-5%";
  else{};
  
  return result;

}


void stylizeHistograms(TH1D *h1,
		       TH1D *h2,
		       bool doC16,
		       bool doC15,
		       bool doC14,
		       bool doC13,
		       bool doC12,
		       bool doC11,
		       bool doC10,
		       bool doC9,
		       bool doC8,
		       bool doC7,
		       bool doC6,
		       bool doC5,
		       bool doC4,
		       bool doC3,
		       bool doC2,
		       bool doC1){

  double lineWidth = 2.0;
  h1->SetLineWidth(lineWidth);
  h2->SetLineWidth(lineWidth);

  h1->SetLineColor(kRed-4);
  h2->SetLineColor(kBlue-4);

  h1->SetStats(0);
  h2->SetStats(0);

  h1->SetTitle("");
  h2->SetTitle("");

  TString title = getTitle(doC16,
			   doC15,
			   doC14,
			   doC13,
			   doC12,
			   doC11,
			   doC10,
			   doC9,
			   doC8,
			   doC7,
			   doC6,
			   doC5,
			   doC4,
			   doC3,
			   doC2,
			   doC1);
  TString xTitle = "Jet #it{p}_{T} [GeV]";
  TString yTitle = "#frac{1}{#it{N}^{evt}} #frac{d#it{N}^{jet}}{d#it{p}_{T}} [GeV^{-1}]";

  double titleSize = 0.05;
  
  double xTitleSize = 0.046;
  double xLabelSize = 0.04;

  double yTitleSize = 0.046;
  double yLabelSize = 0.04;


  h1->SetTitleSize(titleSize);
  h2->SetTitleSize(titleSize);

  h1->SetTitle(title);
  h2->SetTitle(title);

  h1->GetXaxis()->SetLabelSize(xLabelSize);
  h2->GetXaxis()->SetLabelSize(xLabelSize);
  h1->GetXaxis()->SetTitleSize(xTitleSize);
  h2->GetXaxis()->SetTitleSize(xTitleSize);
  
  h1->GetXaxis()->SetTitle(xTitle);
  h2->GetXaxis()->SetTitle(xTitle);

  h1->GetYaxis()->SetLabelSize(yLabelSize);
  h2->GetYaxis()->SetLabelSize(yLabelSize);
  h1->GetYaxis()->SetTitleSize(yTitleSize);
  h2->GetYaxis()->SetTitleSize(yTitleSize);

  
}

void stylizeHiBinHistograms(TH1D *h1, TH1D *h2){

  h1->SetLineColor(kRed-4);
  h2->SetLineColor(kBlue-4);
  
  h1->SetLineWidth(2);
  h2->SetLineWidth(2);

  

}

void plotJetPt_pseudoJets(bool doC16 = false,
			  bool doC15 = false,
			  bool doC14 = false,
			  bool doC13 = false,
			  bool doC12 = false,
			  bool doC11 = false,
			  bool doC10 = false,
			  bool doC9 = false,
			  bool doC8 = false,
			  bool doC7 = false,
			  bool doC6 = false,
			  bool doC5 = false,
			  bool doC4 = false,
			  bool doC3 = false,
			  bool doC2 = true,
			  bool doC1 = false){


  // TFile *f_HYDJET = TFile::Open("/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/scanningOutput/HYDJET/canonical/HYDJET_pThat-unweighted_mu12_pTmu-15to999_tight_hiBinShift-10_jetTrkMaxFilter_WDecayFilter_2026-5-28_ultraFineCentBins.root");
  TFile *f_HYDJET = TFile::Open("/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/scanningOutput/HYDJET/canonical/pfCandAnalyzer/HYDJET_pThat-unweighted_mu12_pTmu-15to999_tight_hiBinShift-0_jetTrkMaxFilter_WDecayFilter_2026-6-4_ultraFineCentBins_pfCandAnalyzer_PFPT-15.root");
  
  TFile *f_MinBias = TFile::Open("/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/scanningOutput/PbPb/canonical/PbPb_MinBias_Part1_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_2026-5-29_ultraFineCentBins.root");

  TLatex *la = new TLatex();
  la->SetTextSize(0.045);

  TLine *li = new TLine();
  li->SetLineStyle(7);
  
  TH2D *H1, *H3;
  TH1D *h1, *h2, *h3; // jet histograms
  TH1D *n1, *n2; // event-counting histograms
  TH1D *h_hiBin_HYDJET, *h_hiBin_MinBias;

  f_HYDJET->GetObject("h_hiBin",h_hiBin_HYDJET);
  f_MinBias->GetObject("h_hiBin",h_hiBin_MinBias);

  stylizeHiBinHistograms(h_hiBin_HYDJET,h_hiBin_MinBias);

  TH1D *h_hiBin_HYDJET_clone = (TH1D*) h_hiBin_HYDJET->Clone("h_hiBin_HYDJET_clone");
  TH1D *h_hiBin_MinBias_clone = (TH1D*) h_hiBin_MinBias->Clone("h_hiBin_MinBias_clone");

  h_hiBin_HYDJET_clone->Scale(1./h_hiBin_HYDJET_clone->Integral());
  h_hiBin_MinBias_clone->Scale(1./h_hiBin_MinBias_clone->Integral());
  divideByBinwidth(h_hiBin_HYDJET_clone);
  divideByBinwidth(h_hiBin_MinBias_clone);
  TH1D *r_hiBin = (TH1D*) h_hiBin_MinBias_clone->Clone("r_hiBin");
  r_hiBin->Divide(h_hiBin_MinBias_clone,h_hiBin_HYDJET_clone,1,1,"");

  TCanvas *canv_hiBin = new TCanvas("canv_hiBin","hiBin distributions",700,700);
  canv_hiBin->cd();
  TPad *pad_hiBin_u = new TPad("pad_hiBin_u","pad_hiBin_u",0,0.4,1,1);
  TPad *pad_hiBin_l = new TPad("pad_hiBin_l","pad_hiBin_l",0,0,1,0.4);
  pad_hiBin_u->SetLeftMargin(0.18);
  pad_hiBin_l->SetLeftMargin(0.18);
  pad_hiBin_u->SetBottomMargin(0.);
  pad_hiBin_l->SetBottomMargin(0.15);
  pad_hiBin_u->SetTopMargin(0.1);
  pad_hiBin_l->SetTopMargin(0.);
  pad_hiBin_u->Draw();
  pad_hiBin_l->Draw();
  pad_hiBin_u->cd();
  h_hiBin_MinBias_clone->GetYaxis()->SetLabelSize(0.05);
  h_hiBin_MinBias_clone->GetYaxis()->SetTitleSize(0.06);
  h_hiBin_MinBias_clone->GetYaxis()->SetTitle("#frac{1}{#it{N}^{evt}_{tot}} #frac{d#it{N}^{evt}}{d(hiBin)}");
  h_hiBin_MinBias_clone->GetYaxis()->SetRangeUser(0.005,0.008);
  h_hiBin_MinBias_clone->SetTitle("");
  h_hiBin_MinBias_clone->SetStats(0);
  h_hiBin_MinBias_clone->Draw();
  h_hiBin_HYDJET_clone->Draw("same");
  TLegend *leg_hiBin = new TLegend(0.65,0.65,0.8,0.8);
  leg_hiBin->SetTextSize(0.055);
  leg_hiBin->SetBorderSize(0.);
  leg_hiBin->AddEntry(h_hiBin_MinBias,"PbPb MinBias");
  leg_hiBin->AddEntry(h_hiBin_HYDJET,"HYDJET");
  leg_hiBin->Draw();
  pad_hiBin_l->cd();
  r_hiBin->GetYaxis()->SetLabelSize(0.06);
  r_hiBin->GetYaxis()->SetTitleSize(0.07);
  r_hiBin->GetYaxis()->SetTitle("Data / MC");
  r_hiBin->GetYaxis()->SetRangeUser(0.75,1.25);
  r_hiBin->GetXaxis()->SetLabelSize(0.06);
  r_hiBin->GetXaxis()->SetTitleSize(0.07);
  r_hiBin->GetXaxis()->SetTitle("hiBin");
  r_hiBin->SetTitle("");
  r_hiBin->SetStats(0);
  r_hiBin->Draw();
  li->DrawLine(20,1,200,1);
    
  TString centString = getCentString(doC16,doC15,doC14,doC13,doC12,doC11,doC10,doC9,doC8,doC7,doC6,doC5,doC4,doC3,doC2,doC1);

  f_HYDJET->GetObject(Form("h_pseudoJetPt_%s",centString.Data()),h1);
  f_HYDJET->GetObject(Form("h_inclRecoJetPt_flavor_%s",centString.Data()),H3);
  h3 = (TH1D*) H3->ProjectionX("h3");
  f_HYDJET->GetObject("h_hiBin",n1);

  f_MinBias->GetObject(Form("h_inclRecoJetPt_%s",centString.Data()),h2);
  f_MinBias->GetObject("h_hiBin",n2);

  double N1 = integrate(n1,getCentHiBinLowerBound(doC16,doC15,doC14,doC13,doC12,doC11,doC10,doC9,doC8,doC7,doC6,doC5,doC4,doC3,doC2,doC1),getCentHiBinUpperBound(doC16,doC15,doC14,doC13,doC12,doC11,doC10,doC9,doC8,doC7,doC6,doC5,doC4,doC3,doC2,doC1));
  double N2 = integrate(n2,getCentHiBinLowerBound(doC16,doC15,doC14,doC13,doC12,doC11,doC10,doC9,doC8,doC7,doC6,doC5,doC4,doC3,doC2,doC1),getCentHiBinUpperBound(doC16,doC15,doC14,doC13,doC12,doC11,doC10,doC9,doC8,doC7,doC6,doC5,doC4,doC3,doC2,doC1));

  cout << "N1 = " << N1 << "\n";
  cout << "N2 = " << N2 << "\n";
  
  // colors: h1=pseudo-jets (black), h3=HYDJET reco (red), h2=PbPb MinBias (blue)
  double lineWidth = 2.0;
  h1->SetLineColor(kBlack);    h1->SetLineWidth(lineWidth); h1->SetStats(0); h1->SetTitle("");
  h3->SetLineColor(kRed+1);    h3->SetLineWidth(lineWidth); h3->SetStats(0); h3->SetTitle("");
  h2->SetLineColor(kAzure+1);  h2->SetLineWidth(lineWidth); h2->SetStats(0); h2->SetTitle("");

  TString xTitle = "Jet #it{p}_{T} [GeV]";
  for(TH1* h : {(TH1*)h1,(TH1*)h2,(TH1*)h3}){
    h->GetXaxis()->SetTitle(xTitle);
    h->GetXaxis()->SetTitleSize(0.046); h->GetXaxis()->SetLabelSize(0.04);
    h->GetYaxis()->SetTitleSize(0.046); h->GetYaxis()->SetLabelSize(0.04);
    h->SetTitle(getTitle(doC16,doC15,doC14,doC13,doC12,doC11,doC10,doC9,doC8,doC7,doC6,doC5,doC4,doC3,doC2,doC1));
  }

  const int N_edge = 11;
  double newAxis[N_edge] = {20,30,40,50,60,70,80,100,120,150,200};
  bool doRebinHistograms = true;
  if(doRebinHistograms){
    h1 = (TH1D*) h1->Rebin(N_edge-1,"h1",newAxis);
    h2 = (TH1D*) h2->Rebin(N_edge-1,"h2",newAxis);
    h3 = (TH1D*) h3->Rebin(N_edge-1,"h3",newAxis);
    divideByBinwidth(h1);
    divideByBinwidth(h2);
    divideByBinwidth(h3);
  }

  // Canvas 2 (per-event): clone before shape normalization
  TH1D *h3b = (TH1D*) h3->Clone("h3b");
  TH1D *h2b = (TH1D*) h2->Clone("h2b");
  h3b->Scale(1./N1);
  h2b->Scale(1./N2);

  // Canvas 1: per-jet shape comparison
  h1->Scale(1./h1->Integral());
  h2->Scale(1./h2->Integral());
  h3->Scale(1./h3->Integral());
  TString yTitle_perJet = "#frac{1}{#it{N}^{jet}} #frac{d#it{N}^{jet}}{d#it{p}_{T}} [GeV^{-1}]";
  h1->GetYaxis()->SetTitle(yTitle_perJet);
  h2->GetYaxis()->SetTitle(yTitle_perJet);
  h3->GetYaxis()->SetTitle(yTitle_perJet);

  // ratios relative to pseudo-jets (h1)
  TH1D *r_data = (TH1D*) h2->Clone("r_data");
  r_data->Divide(h2,h1,1,1,"B");
  r_data->SetLineColor(kAzure+1); r_data->SetLineWidth(lineWidth);

  TH1D *r_mc = (TH1D*) h3->Clone("r_mc");
  r_mc->Divide(h3,h1,1,1,"B");
  r_mc->SetLineColor(kRed+1); r_mc->SetLineWidth(lineWidth);

  TCanvas *canv = new TCanvas("canv","canv",700,700);
  canv->cd();
  TPad *pad_u = new TPad("pad_u","pad_u",0,0.4,1,1);
  TPad *pad_l = new TPad("pad_l","pad_l",0,0,1,0.4);
  pad_u->SetLeftMargin(0.15); pad_l->SetLeftMargin(0.15);
  pad_u->SetBottomMargin(0.); pad_l->SetBottomMargin(0.15);
  pad_u->SetTopMargin(0.15);  pad_l->SetTopMargin(0.);
  pad_u->SetLogy();
  pad_u->Draw(); pad_l->Draw();
  pad_u->cd();

  h1->Draw("hist");
  h3->Draw("hist same");
  h2->Draw("hist same");

  TLegend *leg = new TLegend(0.52,0.65,0.86,0.81);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.036);
  leg->AddEntry(h1,"HYDJET pseudo-jets","l");
  leg->AddEntry(h3,"HYDJET reco jets","l");
  leg->AddEntry(h2,"PbPb MinBias reco jets","l");
  leg->Draw();

  pad_l->cd();

  r_mc->GetYaxis()->SetRangeUser(0, 3.5);

  r_mc->SetTitle("");
  r_mc->GetYaxis()->SetTitleSize(0.07);
  r_mc->GetXaxis()->SetTitleSize(0.07);
  r_mc->GetYaxis()->SetLabelSize(0.06);
  r_mc->GetXaxis()->SetLabelSize(0.06);
  r_mc->GetYaxis()->SetTitle("X / pseudo-jets");
  r_mc->Draw("hist");
  r_data->Draw("hist same");
  li->DrawLine(20,1,200,1);

  TLegend *leg_r = new TLegend(0.17,0.72,0.75,0.92);
  leg_r->SetBorderSize(0);
  leg_r->SetTextSize(0.06);
  leg_r->AddEntry(r_mc,  "HYDJET reco / pseudo-jets","l");
  leg_r->AddEntry(r_data,"PbPb MinBias / pseudo-jets","l");
  leg_r->Draw();

  canv->SaveAs(Form("../../../figures/jetPt/pseudoJets/jetPt_HYDJETVsMinBias_%s.pdf",centString.Data()));

  // --- Second canvas: direct HYDJET reco vs PbPb MinBias (no pseudo-jets) ---
  h3b->SetLineColor(kRed+1);   h3b->SetLineWidth(lineWidth);
  h2b->SetLineColor(kAzure+1); h2b->SetLineWidth(lineWidth);

  TH1D *r_dataMC = (TH1D*) h2b->Clone("r_dataMC");
  r_dataMC->Divide(h2b,h3b,1,1,"B");
  r_dataMC->SetLineColor(kBlack); r_dataMC->SetLineWidth(lineWidth);

  TCanvas *canv2 = new TCanvas("canv2","canv2",700,700);
  canv2->cd();
  TPad *pad2_u = new TPad("pad2_u","pad2_u",0,0.4,1,1);
  TPad *pad2_l = new TPad("pad2_l","pad2_l",0,0,1,0.4);
  pad2_u->SetLeftMargin(0.15); pad2_l->SetLeftMargin(0.15);
  pad2_u->SetBottomMargin(0.); pad2_l->SetBottomMargin(0.15);
  pad2_u->SetTopMargin(0.15);  pad2_l->SetTopMargin(0.);
  pad2_u->SetLogy();
  pad2_u->Draw(); pad2_l->Draw();
  pad2_u->cd();

  h3b->GetYaxis()->SetTitleSize(0.046); h3b->GetYaxis()->SetLabelSize(0.04);
  h3b->GetYaxis()->SetTitle("#frac{1}{#it{N}^{evt}} #frac{d#it{N}^{jet}}{d#it{p}_{T}} [GeV^{-1}]");
  h3b->GetXaxis()->SetTitleSize(0.046); h3b->GetXaxis()->SetLabelSize(0.04);
  h3b->GetXaxis()->SetTitle(xTitle);
  h3b->SetTitle(getTitle(doC16,doC15,doC14,doC13,doC12,doC11,doC10,doC9,doC8,doC7,doC6,doC5,doC4,doC3,doC2,doC1));
  h3b->SetStats(0);
  h3b->Draw("hist");
  h2b->Draw("hist same");

  TLegend *leg2 = new TLegend(0.52,0.65,0.86,0.81);
  leg2->SetBorderSize(0);
  leg2->SetTextSize(0.036);
  leg2->AddEntry(h3b,"HYDJET reco jets","l");
  leg2->AddEntry(h2b,"PbPb MinBias reco jets","l");
  leg2->Draw();

  pad2_l->cd();
  r_dataMC->GetYaxis()->SetRangeUser(0, 3.5);
  r_dataMC->SetTitle("");
  r_dataMC->GetYaxis()->SetTitleSize(0.07);
  r_dataMC->GetXaxis()->SetTitleSize(0.07);
  r_dataMC->GetYaxis()->SetLabelSize(0.06);
  r_dataMC->GetXaxis()->SetLabelSize(0.06);
  r_dataMC->GetYaxis()->SetTitle("MinBias / HYDJET reco (per evt)");
  r_dataMC->GetXaxis()->SetTitle(xTitle);
  r_dataMC->SetStats(0);
  r_dataMC->Draw("hist");
  li->DrawLine(20,1,200,1);

  canv2->SaveAs(Form("../../../figures/jetPt/pseudoJets/jetPt_HYDJETrecoVsMinBias_%s.pdf",centString.Data()));







}
