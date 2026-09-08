



void make_plots(int cent = 1){

  TFile *f = TFile::Open("/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/PbPb_MinBias_Part1_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_sameEventPFClustering_pseudoJetCandPtMin-0.0_2026-7-28_ultraFineCentBins_dPTMapBkgSub.root");

  TH1D *h_raw, *h_JEC;

  f->GetObject(Form("h_inclRawJetPt_C%i",cent),h_raw);
  f->GetObject(Form("h_inclRecoJetPt_C%i",cent),h_JEC);

  h_raw->SetLineColor(kBlue-4);
  h_JEC->SetLineColor(kRed-4);
  h_raw->SetLineWidth(2);
  h_JEC->SetLineWidth(2);
  h_raw->SetStats(0);
  h_raw->SetTitle("");
  h_JEC->SetStats(0);
  h_JEC->SetTitle("");
  
  TH1D *r = (TH1D*) h_JEC->Clone("r");
  r->Divide(h_JEC,h_raw,1,1,"B");

  TCanvas *canv = new TCanvas("canv","canv",700,700);
  canv->cd();
  TPad *pad_l = new TPad("pad_l","pad_l",0,0,1,0.4);
  TPad *pad_u = new TPad("pad_u","pad_u",0,0.4,1,1);
  pad_l->SetTopMargin(0);
  pad_u->SetBottomMargin(0);
  pad_l->SetBottomMargin(0.15);
  pad_u->SetLogy();
  pad_l->Draw();
  pad_u->Draw();
  pad_u->cd();
  h_raw->SetStats(0);
  if(cent == 0)h_raw->SetTitle("PbPb 0-80%");
  else h_raw->SetTitle(Form("PbPb %i-%i%%",5*(cent-1),5*cent));
  h_raw->GetYaxis()->SetTitle("Entries");
  h_raw->GetYaxis()->SetTitleSize(0.05);
  h_raw->GetYaxis()->SetLabelSize(0.04);
  h_raw->Draw();
  h_JEC->Draw("same");
  TLatex *la = new TLatex();
  la->SetTextColor(kRed-4);
  la->DrawLatexNDC(0.6,0.6,"Corrected Jet #it{p}_{T}");
  la->SetTextColor(kBlue-4);
  la->DrawLatexNDC(0.6,0.52,"Raw Jet #it{p}_{T}");
  pad_l->cd();
  r->GetYaxis()->SetRangeUser(0,2.2);
  r->Draw();
  r->GetXaxis()->SetTitleSize(0.06);
  r->GetXaxis()->SetLabelSize(0.05);
  r->GetYaxis()->SetTitleSize(0.06);
  r->GetYaxis()->SetLabelSize(0.05);
  r->GetXaxis()->SetTitle("Jet #it{p}_{T} [GeV]");
  r->GetYaxis()->SetTitle("JEC / raw");
  TLine *li = new TLine();
  li->SetLineStyle(7);
  li->DrawLine(0,1,500,1);

  canv->SaveAs(Form("/home/clayton/Analysis/code/bJetRaaAnalysis/figures/jetPt/JECEffect/JECEffect_C%i.pdf",cent));

  canv->Close();
  f->Close();

  return;
}


void plotJetPt_JECEffect(){

  for(int i = 0; i <= 16; i++){
    make_plots(i);
  }

  return;

}
