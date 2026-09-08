

string get_centBin_string(int centBin = 1){

  // cent from {0..16}
  string result = "";

  if(centBin == 0) result = "PbPb 0-80%";
  else result = Form("PbPb %i-%i%%",5*(centBin-1),5*centBin);
  
  return result;

}

string get_jetBin_string(int jetBin = 1){

  string result = "";
  
  if(jetBin == 0) result = "20 < #it{p}_{T}^{jet} (GeV) < 500";
  else if(jetBin == 1) result = "80 < #it{p}_{T}^{jet} (GeV) < 100";
  else if(jetBin == 2) result = "100 < #it{p}_{T}^{jet} (GeV) < 120";
  else if(jetBin == 3) result = "120 < #it{p}_{T}^{jet} (GeV) < 150";
  else if(jetBin == 4) result = "150 < #it{p}_{T}^{jet} (GeV) < 200";
  else if(jetBin == 5) result = "200 < #it{p}_{T}^{jet} (GeV) < 300";
  else if(jetBin == 6) result = "300 < #it{p}_{T}^{jet} (GeV) < 500";

  return result;

}


void plotDRMap(int centBin = 1,
	       int jetBin = 1){

  TFile *f = TFile::Open("/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/PbPb_MinBias_Part1_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_sameEventPFClustering_pseudoJetCandPtMin-0.0_2026-8-4_ultraFineCentBins.root");

  TH2D *H_dR, *H_eta;
  TH1D *h_dR, *h_eta;

  f->GetObject(Form("h_fastJetPtPF_dRmin_C%iJ%i",centBin,jetBin),H_dR);
  f->GetObject(Form("h_fastJetPtPF_etaPFCs_C%iJ%i",centBin,jetBin),H_eta);

  h_dR = (TH1D*) H_dR->ProjectionY("h_dR");
  h_eta = (TH1D*) H_eta->ProjectionY("h_eta");

  string centBin_string = get_centBin_string(centBin);
  string jetBin_string = get_jetBin_string(jetBin);

  TCanvas *canv_dR = new TCanvas("canv_dR","canv_dR",700,700);
  canv_dR->cd();
  TPad *pad_dR = new TPad("pad_dR","pad_dR",0,0,1,1);
  pad_dR->SetLeftMargin(0.15);
  pad_dR->SetRightMargin(0.15);
  pad_dR->Draw();
  pad_dR->cd();
  H_dR->SetStats(0);
  H_dR->SetTitle(Form("%s, %s",centBin_string.c_str(),jetBin_string.c_str()));
  H_dR->GetXaxis()->SetTitle("#it{p}_{T}^{PF} [GeV]");
  H_dR->GetYaxis()->SetTitle("#it{#Delta R}(PF,PFCs)");
  H_dR->GetZaxis()->SetTitle("Entries");
  H_dR->Draw("colz");

  canv_dR->SaveAs(Form("../../../../figures/jetPt/jetClustering/dRMap/fastJetPtPF_dRmin_C%iJ%i.pdf",centBin,jetBin));


  TCanvas *canv_eta = new TCanvas("canv_eta","canv_eta",700,700);
  canv_eta->cd();
  TPad *pad_eta = new TPad("pad_eta","pad_eta",0,0,1,1);
  pad_eta->SetLeftMargin(0.15);
  pad_eta->SetRightMargin(0.15);
  pad_eta->Draw();
  pad_eta->cd();
  H_eta->SetStats(0);
  H_eta->SetTitle(Form("%s, %s",centBin_string.c_str(),jetBin_string.c_str()));
  H_eta->GetXaxis()->SetTitle("#it{p}_{T}^{PF} [GeV]");
  H_eta->GetYaxis()->SetTitle("#it{#eta}^{PFCs}");
  H_eta->GetZaxis()->SetTitle("Entries");
  H_eta->Draw("colz");

  canv_eta->SaveAs(Form("../../../../figures/jetPt/jetClustering/etaMap/fastJetPtPF_etaCs_C%iJ%i.pdf",centBin,jetBin));


  TCanvas *canv_dRproj = new TCanvas("canv_dRproj","canv_dRproj",700,700);
  canv_dRproj->cd();
  TPad *pad_dRproj = new TPad("pad_dRproj","pad_dRproj",0,0,1,1);
  pad_dRproj->SetLeftMargin(0.15);
  pad_dRproj->SetBottomMargin(0.15);
  pad_dRproj->Draw();
  pad_dRproj->cd();
  h_dR->SetStats(0);
  h_dR->SetTitle(Form("%s, %s",centBin_string.c_str(),jetBin_string.c_str()));
  h_dR->GetXaxis()->SetTitle("#it{#Delta R}(PF,PFCs)");
  h_dR->GetYaxis()->SetTitle("Entries");
  h_dR->Draw();

  canv_dRproj->SaveAs(Form("../../../../figures/jetPt/jetClustering/dRproj/dRproj_C%iJ%i.pdf",centBin,jetBin));

  TCanvas *canv_etaproj = new TCanvas("canv_etaproj","canv_etaproj",700,700);
  canv_etaproj->cd();
  TPad *pad_etaproj = new TPad("pad_etaproj","pad_etaproj",0,0,1,1);
  pad_etaproj->SetLeftMargin(0.15);
  pad_etaproj->SetBottomMargin(0.15);
  pad_etaproj->Draw();
  pad_etaproj->cd();
  h_eta->SetStats(0);
  h_eta->SetTitle(Form("%s, %s",centBin_string.c_str(),jetBin_string.c_str()));
  h_eta->GetXaxis()->SetTitle("#it{#eta}^{PFCs}");
  h_eta->GetYaxis()->SetTitle("Entries");
  h_eta->Draw();

  canv_etaproj->SaveAs(Form("../../../../figures/jetPt/jetClustering/etaproj/etaproj_C%iJ%i.pdf",centBin,jetBin));




  canv_dR->Close();
  canv_eta->Close();
  canv_dRproj->Close();

  return;


}
