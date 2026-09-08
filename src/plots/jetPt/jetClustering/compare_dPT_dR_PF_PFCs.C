
string get_centrality_string(int cent = 1){

  // cent from {0..16}
  string result = "";

  if(cent == 0) result = "PbPb 0-80%";
  else result = Form("PbPb %i-%i%%",5*(cent-1),5*cent);
  
  return result;

}

double getPTmin(int cent = 1){
  double pTmin = 0.;
  if(cent == 0) pTmin = 9.;
  else if(cent == 1) pTmin = 30.;
  else if(cent == 2) pTmin = 20.;
  else if(cent == 3) pTmin = 15.;
  else if(cent == 4) pTmin = 10.;
  else if(cent == 5) pTmin = 13.;
  else if(cent == 6) pTmin = 13.;
  else if(cent == 7) pTmin = 11.;
  else if(cent == 8) pTmin = 11.;
  else if(cent == 9) pTmin = 8.;
  else if(cent == 10) pTmin = 6.;
  else if(cent == 11) pTmin = 5.;
  else if(cent == 12) pTmin = 3.5;
  else if(cent == 13) pTmin = 2.5;
  else if(cent == 14) pTmin = 2.0;
  else if(cent == 15) pTmin = 1.5;
  else if(cent == 16) pTmin = 1.0;
  else{};
  return pTmin;
}

double getPTmax(int cent = 1){
  double pTmin = 0.;
  if(cent == 0) pTmin = 33.;
  else if(cent == 1) pTmin = 120.;
  else if(cent == 2) pTmin = 90.;
  else if(cent == 3) pTmin = 70.;
  else if(cent == 4) pTmin = 60.;
  else if(cent == 5) pTmin = 47.;
  else if(cent == 6) pTmin = 39.;
  else if(cent == 7) pTmin = 31.;
  else if(cent == 8) pTmin = 24.;
  else if(cent == 9) pTmin = 19.;
  else if(cent == 10) pTmin = 15.;
  else if(cent == 11) pTmin = 11.;
  else if(cent == 12) pTmin = 8.;
  else if(cent == 13) pTmin = 5.5;
  else if(cent == 14) pTmin = 3.9;
  else if(cent == 15) pTmin = 2.7;
  else if(cent == 16) pTmin = 1.9;
  else{};
  return pTmin;
}

void create_plots(int cent = 1){

  bool doCommonAxisNormalization = false;

  //TFile *f = TFile::Open("/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/PbPb_MinBias_Part1_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_sameEventPFClustering_pseudoJetCandPtMin-0.0_2026-7-28_ultraFineCentBins.root");
  //TFile *f = TFile::Open("/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/PbPb_MinBias_Part1_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_sameEventPFClustering_skipSingleConstituentJets_pseudoJetCandPtMin-0.0_2026-7-29_ultraFineCentBins_dPTMapBkgSub.root");
  TFile *f = TFile::Open("/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/PbPb_MinBias_Part1_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_sameEventPFClustering_pseudoJetCandPtMin-0.0_2026-8-4_ultraFineCentBins.root");
  
  TH1D *h_dR; // dR distribution
  TProfile2D *H_dPT; // dPT distribution

  f->GetObject(Form("h_dRmin_PF_PFCs_C%i",cent),h_dR);
  f->GetObject(Form("h_randConeEtaPhi_C%i",cent),H_dPT);
  //f->GetObject(Form("h_dPTEtaPhi_PF_PFCs_C%i",cent),H_dPT);

  string centrality_string = get_centrality_string(cent);
  
  TCanvas *canv_dR = new TCanvas("canv_dR","canv_dR",700,700);
  canv_dR->cd();
  TPad *pad_dR = new TPad("pad_dR","pad_dR",0,0,1,1);
  pad_dR->Draw();
  pad_dR->SetLeftMargin(0.15);
  pad_dR->cd();
  h_dR->SetStats(0);
  h_dR->SetTitle(centrality_string.c_str());
  h_dR->GetYaxis()->SetTitle("Entries");
  h_dR->GetXaxis()->SetTitle("#it{#Delta}R(PF,PFCs)");
  h_dR->Draw();

  TCanvas *canv_dPT = new TCanvas("canv_dPT","canv_dPT",700,700);
  canv_dPT->cd();
  TPad *pad_dPT = new TPad("pad_dPT","pad_dPT",0,0,1,1);
  pad_dPT->Draw();
  pad_dPT->SetRightMargin(0.20);
  pad_dPT->cd();
  H_dPT->SetStats(0);
  H_dPT->SetTitle(centrality_string.c_str());
  H_dPT->GetYaxis()->SetTitle("#it{#phi}");
  H_dPT->GetXaxis()->SetTitle("#it{#eta}");
  //H_dPT->GetZaxis()->SetTitle("#LT #it{p}_{T}^{PF} - #it{p}_{T}^{PFCs} #GT");
  H_dPT->GetZaxis()->SetTitle("#LT #it{p}_{T}^{RC}#GT");
  double pTmin = 0.;
  double pTmax = 120.;
  if(doCommonAxisNormalization){
    pTmin = getPTmin(cent);
    pTmax = getPTmax(cent);
    H_dPT->GetZaxis()->SetRangeUser(pTmin,pTmax);
  }
  double titleSize = 0.045;
  H_dPT->GetYaxis()->SetTitleSize(titleSize);
  H_dPT->GetXaxis()->SetTitleSize(titleSize);
  H_dPT->GetZaxis()->SetTitleSize(titleSize);
  H_dPT->Draw("colz");

  string commonAxisNormalizationStatusString = "";
  if(doCommonAxisNormalization) commonAxisNormalizationStatusString = "_commonAxisNormalization";

  canv_dR->SaveAs(Form("/home/clayton/Analysis/code/bJetRaaAnalysis/figures/jetPt/jetClustering/dR/dR_PF_PFCs_C%i.pdf",cent));
  //canv_dPT->SaveAs(Form("/home/clayton/Analysis/code/bJetRaaAnalysis/figures/jetPt/jetClustering/dPTEtaPhi/dPTEtaPhi_PF_PFCs_C%i%s.pdf",cent,commonAxisNormalizationStatusString.c_str()));
  canv_dPT->SaveAs(Form("/home/clayton/Analysis/code/bJetRaaAnalysis/figures/jetPt/jetClustering/randConeEtaPhi/randConeEtaPhi_PF_PFCs_C%i%s.pdf",cent,commonAxisNormalizationStatusString.c_str()));

  canv_dR->Close();
  canv_dPT->Close();
  f->Close();

  return;

}

void compare_dPT_dR_PF_PFCs(){
  
  for(int i = 0; i <= 16; i++){
    create_plots(i);
  }

  

}
