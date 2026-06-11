


void stylizeHistograms(TH1D *h1, TH1D *h2, TH1D *h3, TH1D *h4, TH1D *h5){

  // Colorblind-friendly colors (Paul Tol's bright palette)
  int cBlue   = TColor::GetColor("#4477AA");
  int cCyan   = TColor::GetColor("#66CCEE");
  int cGreen  = TColor::GetColor("#228833");
  int cYellow = TColor::GetColor("#CCBB44");
  int cRed    = TColor::GetColor("#EE6677");
  int cPurple = TColor::GetColor("#AA3377");

  h1->SetLineColor(cBlue);
  h2->SetLineColor(cCyan);
  h3->SetLineColor(cGreen);
  h4->SetLineColor(cYellow);
  h5->SetLineColor(cRed);


  double lineWidth = 2.;
  h1->SetLineWidth(lineWidth);
  h2->SetLineWidth(lineWidth);
  h3->SetLineWidth(lineWidth);
  h4->SetLineWidth(lineWidth);
  h5->SetLineWidth(lineWidth);


}

void jetTriggerTurnOn(){

  TFile *file_pp_jet30 = TFile::Open("/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/scanningOutput/pp/canonical/pp_HighEGJet_Jet30HLT_mu12_pTmu-15to999_tight_deltaR-40_jetTrkMaxFilter_WDecayFilter_2026-6-1.root");
  TFile *file_pp_jet30AndJet40 = TFile::Open("/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/scanningOutput/pp/canonical/pp_HighEGJet_Jet30HLT_Jet40HLT_mu12_pTmu-15to999_tight_deltaR-40_jetTrkMaxFilter_WDecayFilter_2026-6-2.root");

  TFile *file_pp_jet40 = TFile::Open("/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/scanningOutput/pp/canonical/pp_HighEGJet_Jet40HLT_mu12_pTmu-15to999_tight_deltaR-40_jetTrkMaxFilter_WDecayFilter_2026-6-1.root");
  TFile *file_pp_jet40AndJet60 = TFile::Open("/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/scanningOutput/pp/canonical/pp_HighEGJet_Jet40HLT_Jet60HLT_mu12_pTmu-15to999_tight_deltaR-40_jetTrkMaxFilter_WDecayFilter_2026-6-2.root");

  TFile *file_pp_jet60 = TFile::Open("/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/scanningOutput/pp/canonical/pp_HighEGJet_Jet60HLT_mu12_pTmu-15to999_tight_deltaR-40_jetTrkMaxFilter_WDecayFilter_2026-6-2.root");
  TFile *file_pp_jet60AndJet80 = TFile::Open("/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/scanningOutput/pp/canonical/pp_HighEGJet_Jet60HLT_Jet80HLT_mu12_pTmu-15to999_tight_deltaR-40_jetTrkMaxFilter_WDecayFilter_2026-6-2.root");
  
  TFile *file_pp_jet80 = TFile::Open("/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/scanningOutput/pp/canonical/pp_HighEGJet_Jet80HLT_mu12_pTmu-15to999_tight_deltaR-40_jetTrkMaxFilter_WDecayFilter_2026-6-2.root");
  TFile *file_pp_jet100 = TFile::Open("/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/scanningOutput/pp/canonical/pp_HighEGJet_Jet100HLT_mu12_pTmu-15to999_tight_deltaR-40_jetTrkMaxFilter_WDecayFilter_2026-6-2.root");
  TFile *file_pp_incl = TFile::Open("/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/scanningOutput/pp/canonical/pp_HighEGJet_mu12_pTmu-15to999_tight_deltaR-40_jetTrkMaxFilter_WDecayFilter_2026-6-1.root");

  TH1D *h_pp_jet15, *h_pp_jet30, *h_pp_jet40, *h_pp_jet60, *h_pp_jet80, *h_pp_jet100, *h_pp_incl;
  TH1D *h_pp_jet30AndJet40, *h_pp_jet40AndJet60, *h_pp_jet60AndJet80;

  file_pp_jet30->GetObject("h_inclRecoJetPt",h_pp_jet30);
  file_pp_jet30AndJet40->GetObject("h_inclRecoJetPt",h_pp_jet30AndJet40);
  file_pp_jet40->GetObject("h_inclRecoJetPt",h_pp_jet40);
  file_pp_jet40AndJet60->GetObject("h_inclRecoJetPt",h_pp_jet40AndJet60);
  file_pp_jet60->GetObject("h_inclRecoJetPt",h_pp_jet60);
  file_pp_jet60AndJet80->GetObject("h_inclRecoJetPt",h_pp_jet60AndJet80);
  file_pp_jet80->GetObject("h_inclRecoJetPt",h_pp_jet80);
  file_pp_jet100->GetObject("h_inclRecoJetPt",h_pp_jet100);
  file_pp_incl->GetObject("h_inclRecoJetPt",h_pp_incl);

  stylizeHistograms(h_pp_jet30, h_pp_jet40, h_pp_jet60, h_pp_jet80, h_pp_jet100);

  TH1D *r_pp_jet15, *r_pp_jet30, *r_pp_jet40, *r_pp_jet60, *r_pp_jet80, *r_pp_jet100; // ratios wrt jet100

  //r_pp_jet15 = (TH1D*) h_pp_jet15->Clone("r_pp_jet15");
  r_pp_jet30 = (TH1D*) h_pp_jet30->Clone("r_pp_jet30");
  r_pp_jet40 = (TH1D*) h_pp_jet40->Clone("r_pp_jet40");
  r_pp_jet60 = (TH1D*) h_pp_jet60->Clone("r_pp_jet60");
  r_pp_jet80 = (TH1D*) h_pp_jet80->Clone("r_pp_jet80");
  r_pp_jet100 = (TH1D*) h_pp_jet100->Clone("r_pp_jet100");

  //r_pp_jet15->Divide(h_pp_jet15,h_pp_incl,1,1,"B");
  r_pp_jet30->Divide(h_pp_jet30AndJet40,h_pp_jet40,1,1,"B");
  r_pp_jet40->Divide(h_pp_jet40AndJet60,h_pp_jet60,1,1,"B");
  r_pp_jet60->Divide(h_pp_jet60AndJet80,h_pp_jet80,1,1,"B");
  r_pp_jet80->Divide(h_pp_jet80,h_pp_incl,1,1,"B");
  r_pp_jet100->Divide(h_pp_jet100,h_pp_incl,1,1,"B");

  TLine *li = new TLine();
  li->SetLineStyle(7);

  TCanvas *canv_pp = new TCanvas("canv_pp","canv_pp",700,700);
  canv_pp->cd();
  TPad *pad_pp = new TPad("pad_pp","pad_pp",0,0,1,1);
  pad_pp->SetLeftMargin(0.15);
  pad_pp->SetBottomMargin(0.15);
  pad_pp->Draw();
  pad_pp->cd();
  r_pp_jet30->SetTitle("");
  r_pp_jet30->SetStats(0);
  r_pp_jet30->GetYaxis()->SetLabelSize(0.035);
  r_pp_jet30->GetYaxis()->SetTitleSize(0.045);
  r_pp_jet30->GetYaxis()->SetTitle("Trigger Efficiency");
  r_pp_jet30->GetXaxis()->SetLabelSize(0.035);
  r_pp_jet30->GetXaxis()->SetTitleSize(0.045);
  r_pp_jet30->GetXaxis()->SetTitle("Jet #it{p}_{T} [GeV]");
  r_pp_jet30->GetYaxis()->SetRangeUser(0,1.1);
  r_pp_jet30->Draw();
  r_pp_jet40->Draw("same");
  r_pp_jet60->Draw("same");
  r_pp_jet80->Draw("same");
  r_pp_jet100->Draw("same");

  TLegend *leg_pp = new TLegend(0.65,0.35,0.85,0.55);
  leg_pp->SetBorderSize(0);
  leg_pp->SetTextSize(0.045);
  leg_pp->AddEntry(r_pp_jet30,"Jet30+Jet40 / Jet40");
  leg_pp->AddEntry(r_pp_jet40,"Jet40+Jet60 / Jet60");
  leg_pp->AddEntry(r_pp_jet60,"Jet60+Jet80 / Jet80");
  leg_pp->AddEntry(r_pp_jet80,"Jet80 / all");
  leg_pp->AddEntry(r_pp_jet100,"Jet100 / all");
  leg_pp->Draw();

  
  li->DrawLine(20,1,500,1);
  

}
