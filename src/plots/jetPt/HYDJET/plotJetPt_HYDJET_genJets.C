#include "../../../../headers/functions/divideByBinwidth.h"


void stylizeHistograms(TH1D *h_C1,
		       TH1D *h_C2,
		       TH1D *h_C3,
		       TH1D *h_C4){

  double lineWidth = 2;
  h_C1->SetLineWidth(lineWidth);
  h_C2->SetLineWidth(lineWidth);
  h_C3->SetLineWidth(lineWidth);
  h_C4->SetLineWidth(lineWidth);

  h_C1->SetLineColor(kRed-4);
  h_C2->SetLineColor(kGreen+2);
  h_C3->SetLineColor(kBlue-4);
  h_C4->SetLineColor(kBlack);

  return;

}


void plotJetPt_HYDJET_genJets(){

  TFile *f = TFile::Open("/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/HYDJET/HYDJET_pThat-unweighted_mu12_pTmu-15to999_tight_hiBinShift-0_jetTrkMaxFilter_WDecayFilter_2026-8-6_ultraFineCentBins.root");

  //TFile *f = TFile::Open("/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/scanningOutput/PYTHIAHYDJET/latest/PYTHIAHYDJET_DiJet_pThat-15_mu12_pTmu-15to999_tight_vzReweight_hiBinReweightToHardProbesJet80_hiBinShift-10_jetTrkMaxFilter_doPThatCorrelationFilterTight_WDecayFilter_2026-5-12_ultraFineCentBins.root");

  TH2D *H_C1, *H_C2, *H_C3, *H_C4, *H_C5,
    *H_C6, *H_C7, *H_C8, *H_C9, *H_C10,
    *H_C11, *H_C12, *H_C13, *H_C14, *H_C15, *H_C16;

  TH2D *Hr_C1, *Hr_C2, *Hr_C3, *Hr_C4, *Hr_C5,
    *Hr_C6, *Hr_C7, *Hr_C8, *Hr_C9, *Hr_C10,
    *Hr_C11, *Hr_C12, *Hr_C13, *Hr_C14, *Hr_C15, *Hr_C16;

  TH1D *h_C1, *h_C2, *h_C3, *h_C4, *h_C5,
    *h_C6, *h_C7, *h_C8, *h_C9, *h_C10,
    *h_C11, *h_C12, *h_C13, *h_C14, *h_C15, *h_C16;

  TH1D *hr_C1, *hr_C2, *hr_C3, *hr_C4, *hr_C5,
    *hr_C6, *hr_C7, *hr_C8, *hr_C9, *hr_C10,
    *hr_C11, *hr_C12, *hr_C13, *hr_C14, *hr_C15, *hr_C16;

  TH1D *vz_C1, *vz_C2, *vz_C3, *vz_C4, *vz_C5,
    *vz_C6, *vz_C7, *vz_C8, *vz_C9, *vz_C10,
    *vz_C11, *vz_C12, *vz_C13, *vz_C14, *vz_C15, *vz_C16;

  TH1D *g_C1, *g_C2, *g_C3, *g_C4;  // grouped histograms into coarse centrality bins
  TH1D *gr_C1, *gr_C2, *gr_C3, *gr_C4;  // grouped histograms into coarse centrality bins
  TH1D *gvz_C1, *gvz_C2, *gvz_C3, *gvz_C4;

  f->GetObject("h_inclGenJetPt_flavor_C1", H_C1);
  f->GetObject("h_inclGenJetPt_flavor_C2", H_C2);
  f->GetObject("h_inclGenJetPt_flavor_C3", H_C3);
  f->GetObject("h_inclGenJetPt_flavor_C4", H_C4);
  f->GetObject("h_inclGenJetPt_flavor_C5", H_C5);
  f->GetObject("h_inclGenJetPt_flavor_C6", H_C6);
  f->GetObject("h_inclGenJetPt_flavor_C7", H_C7);
  f->GetObject("h_inclGenJetPt_flavor_C8", H_C8);
  f->GetObject("h_inclGenJetPt_flavor_C9", H_C9);
  f->GetObject("h_inclGenJetPt_flavor_C10", H_C10);
  f->GetObject("h_inclGenJetPt_flavor_C11", H_C11);
  f->GetObject("h_inclGenJetPt_flavor_C12", H_C12);
  f->GetObject("h_inclGenJetPt_flavor_C13", H_C13);
  f->GetObject("h_inclGenJetPt_flavor_C14", H_C14);
  f->GetObject("h_inclGenJetPt_flavor_C15", H_C15);
  f->GetObject("h_inclGenJetPt_flavor_C16", H_C16);

  f->GetObject("h_inclRecoJetPt_flavor_C1", Hr_C1);
  f->GetObject("h_inclRecoJetPt_flavor_C2", Hr_C2);
  f->GetObject("h_inclRecoJetPt_flavor_C3", Hr_C3);
  f->GetObject("h_inclRecoJetPt_flavor_C4", Hr_C4);
  f->GetObject("h_inclRecoJetPt_flavor_C5", Hr_C5);
  f->GetObject("h_inclRecoJetPt_flavor_C6", Hr_C6);
  f->GetObject("h_inclRecoJetPt_flavor_C7", Hr_C7);
  f->GetObject("h_inclRecoJetPt_flavor_C8", Hr_C8);
  f->GetObject("h_inclRecoJetPt_flavor_C9", Hr_C9);
  f->GetObject("h_inclRecoJetPt_flavor_C10", Hr_C10);
  f->GetObject("h_inclRecoJetPt_flavor_C11", Hr_C11);
  f->GetObject("h_inclRecoJetPt_flavor_C12", Hr_C12);
  f->GetObject("h_inclRecoJetPt_flavor_C13", Hr_C13);
  f->GetObject("h_inclRecoJetPt_flavor_C14", Hr_C14);
  f->GetObject("h_inclRecoJetPt_flavor_C15", Hr_C15);
  f->GetObject("h_inclRecoJetPt_flavor_C16", Hr_C16);

  h_C1 = (TH1D*) H_C1->ProjectionX("h_C1");
  h_C2 = (TH1D*) H_C2->ProjectionX("h_C2");
  h_C3 = (TH1D*) H_C3->ProjectionX("h_C3");
  h_C4 = (TH1D*) H_C4->ProjectionX("h_C4");
  h_C5 = (TH1D*) H_C5->ProjectionX("h_C5");
  h_C6 = (TH1D*) H_C6->ProjectionX("h_C6");
  h_C7 = (TH1D*) H_C7->ProjectionX("h_C7");
  h_C8 = (TH1D*) H_C8->ProjectionX("h_C8");
  h_C9 = (TH1D*) H_C9->ProjectionX("h_C9");
  h_C10 = (TH1D*) H_C10->ProjectionX("h_C10");
  h_C11 = (TH1D*) H_C11->ProjectionX("h_C11");
  h_C12 = (TH1D*) H_C12->ProjectionX("h_C12");
  h_C13 = (TH1D*) H_C13->ProjectionX("h_C13");
  h_C14 = (TH1D*) H_C14->ProjectionX("h_C14");
  h_C15 = (TH1D*) H_C15->ProjectionX("h_C15");
  h_C16 = (TH1D*) H_C16->ProjectionX("h_C16");

  hr_C1 = (TH1D*) Hr_C1->ProjectionX("hr_C1");
  hr_C2 = (TH1D*) Hr_C2->ProjectionX("hr_C2");
  hr_C3 = (TH1D*) Hr_C3->ProjectionX("hr_C3");
  hr_C4 = (TH1D*) Hr_C4->ProjectionX("hr_C4");
  hr_C5 = (TH1D*) Hr_C5->ProjectionX("hr_C5");
  hr_C6 = (TH1D*) Hr_C6->ProjectionX("hr_C6");
  hr_C7 = (TH1D*) Hr_C7->ProjectionX("hr_C7");
  hr_C8 = (TH1D*) Hr_C8->ProjectionX("hr_C8");
  hr_C9 = (TH1D*) Hr_C9->ProjectionX("hr_C9");
  hr_C10 = (TH1D*) Hr_C10->ProjectionX("hr_C10");
  hr_C11 = (TH1D*) Hr_C11->ProjectionX("hr_C11");
  hr_C12 = (TH1D*) Hr_C12->ProjectionX("hr_C12");
  hr_C13 = (TH1D*) Hr_C13->ProjectionX("hr_C13");
  hr_C14 = (TH1D*) Hr_C14->ProjectionX("hr_C14");
  hr_C15 = (TH1D*) Hr_C15->ProjectionX("hr_C15");
  hr_C16 = (TH1D*) Hr_C16->ProjectionX("hr_C16");

  f->GetObject("h_vz_C1", vz_C1);
  f->GetObject("h_vz_C2", vz_C2);
  f->GetObject("h_vz_C3", vz_C3);
  f->GetObject("h_vz_C4", vz_C4);
  f->GetObject("h_vz_C5", vz_C5);
  f->GetObject("h_vz_C6", vz_C6);
  f->GetObject("h_vz_C7", vz_C7);
  f->GetObject("h_vz_C8", vz_C8);
  f->GetObject("h_vz_C9", vz_C9);
  f->GetObject("h_vz_C10", vz_C10);
  f->GetObject("h_vz_C11", vz_C11);
  f->GetObject("h_vz_C12", vz_C12);
  f->GetObject("h_vz_C13", vz_C13);
  f->GetObject("h_vz_C14", vz_C14);
  f->GetObject("h_vz_C15", vz_C15);
  f->GetObject("h_vz_C16", vz_C16);
  

  g_C1 = (TH1D*) h_C1->Clone("g_C1");
  g_C1->Add(h_C2);
  gr_C1 = (TH1D*) hr_C1->Clone("gr_C1");
  gr_C1->Add(hr_C2);
  gvz_C1 = (TH1D*) vz_C1->Clone("gvz_C1");
  gvz_C1->Add(vz_C2);

  g_C2 = (TH1D*) h_C3->Clone("g_C2");
  g_C2->Add(h_C4);
  g_C2->Add(h_C5);
  g_C2->Add(h_C6);
  gr_C2 = (TH1D*) hr_C3->Clone("gr_C2");
  gr_C2->Add(hr_C4);
  gr_C2->Add(hr_C5);
  gr_C2->Add(hr_C6);
  gvz_C2 = (TH1D*) vz_C3->Clone("gvz_C2");
  gvz_C2->Add(vz_C4);
  gvz_C2->Add(vz_C5);
  gvz_C2->Add(vz_C6);

  g_C3 = (TH1D*) h_C7->Clone("g_C3");
  g_C3->Add(h_C8);
  g_C3->Add(h_C9);
  g_C3->Add(h_C10);
  gr_C3 = (TH1D*) hr_C7->Clone("gr_C3");
  gr_C3->Add(hr_C8);
  gr_C3->Add(hr_C9);
  gr_C3->Add(hr_C10);
  gvz_C3 = (TH1D*) vz_C7->Clone("gvz_C3");
  gvz_C3->Add(vz_C8);
  gvz_C3->Add(vz_C9);
  gvz_C3->Add(vz_C10);


  g_C4 = (TH1D*) h_C11->Clone("g_C4");
  g_C4->Add(h_C12);
  g_C4->Add(h_C13);
  g_C4->Add(h_C14);
  g_C4->Add(h_C15);
  g_C4->Add(h_C16);
  gr_C4 = (TH1D*) hr_C11->Clone("gr_C4");
  gr_C4->Add(hr_C12);
  gr_C4->Add(hr_C13);
  gr_C4->Add(hr_C14);
  gr_C4->Add(hr_C15);
  gr_C4->Add(hr_C16);
  gvz_C4 = (TH1D*) vz_C11->Clone("gvz_C4");
  gvz_C4->Add(vz_C12);
  gvz_C4->Add(vz_C13);
  gvz_C4->Add(vz_C14);
  gvz_C4->Add(vz_C15);
  gvz_C4->Add(vz_C16);


  

  // rebin
  const int N_edge = 7;
  double newAxis[N_edge] = {20,30,40,50,80,150,300};

  g_C1 = (TH1D*) g_C1->Rebin(N_edge-1,"g_C1",newAxis);
  g_C2 = (TH1D*) g_C2->Rebin(N_edge-1,"g_C2",newAxis);
  g_C3 = (TH1D*) g_C3->Rebin(N_edge-1,"g_C3",newAxis);
  g_C4 = (TH1D*) g_C4->Rebin(N_edge-1,"g_C4",newAxis);

  gr_C1 = (TH1D*) gr_C1->Rebin(N_edge-1,"gr_C1",newAxis);
  gr_C2 = (TH1D*) gr_C2->Rebin(N_edge-1,"gr_C2",newAxis);
  gr_C3 = (TH1D*) gr_C3->Rebin(N_edge-1,"gr_C3",newAxis);
  gr_C4 = (TH1D*) gr_C4->Rebin(N_edge-1,"gr_C4",newAxis);

  divideByBinwidth(g_C1);
  divideByBinwidth(g_C2);
  divideByBinwidth(g_C3);
  divideByBinwidth(g_C4);

  divideByBinwidth(gr_C1);
  divideByBinwidth(gr_C2);
  divideByBinwidth(gr_C3);
  divideByBinwidth(gr_C4);


  // normalize histograms
  g_C1->Scale(1./gvz_C1->Integral());
  g_C2->Scale(1./gvz_C2->Integral());
  g_C3->Scale(1./gvz_C3->Integral());
  g_C4->Scale(1./gvz_C4->Integral());

  gr_C1->Scale(1./gvz_C1->Integral());
  gr_C2->Scale(1./gvz_C2->Integral());
  gr_C3->Scale(1./gvz_C3->Integral());
  gr_C4->Scale(1./gvz_C4->Integral());


  // g_C1->Scale(1./g_C1->Integral(g_C1->FindBin(newAxis[0]),g_C1->FindBin(newAxis[N_edge-1])));
  // g_C2->Scale(1./g_C2->Integral(g_C2->FindBin(newAxis[0]),g_C2->FindBin(newAxis[N_edge-1])));
  // g_C3->Scale(1./g_C3->Integral(g_C3->FindBin(newAxis[0]),g_C3->FindBin(newAxis[N_edge-1])));
  // g_C4->Scale(1./g_C4->Integral(g_C4->FindBin(newAxis[0]),g_C4->FindBin(newAxis[N_edge-1])));
  


  stylizeHistograms(g_C1,g_C2,g_C3,g_C4);
  
  TH1D *r_C1 = (TH1D*) g_C1->Clone("r_C1");
  TH1D *r_C2 = (TH1D*) g_C2->Clone("r_C2");
  TH1D *r_C3 = (TH1D*) g_C3->Clone("r_C3");
  
  r_C1->Divide(g_C1,g_C4,1,1,"");
  r_C2->Divide(g_C2,g_C4,1,1,"");
  r_C3->Divide(g_C3,g_C4,1,1,"");
  

  TCanvas *canv = new TCanvas("canv","canv",700,700);
  canv->cd();
  TPad *pad_u = new TPad("pad_u","pad_u",0,0.4,1,1);
  TPad *pad_l = new TPad("pad_l","pad_l",0,0,1,0.4);
  pad_u->SetBottomMargin(0.0);
  pad_l->SetTopMargin(0.0);
  pad_l->SetBottomMargin(0.18);
  pad_u->SetLeftMargin(0.15);
  pad_l->SetLeftMargin(0.15);
  pad_u->SetLogy();
  pad_u->Draw();
  pad_l->Draw();
  pad_u->cd();
  g_C1->GetYaxis()->SetTitle("1/#it{N}_{evt} d#it{N} / d#it{p}_{T} [GeV^{-1}]");
  g_C1->GetXaxis()->SetRangeUser(newAxis[0],newAxis[N_edge-1]);
  g_C1->GetYaxis()->SetTitleSize(0.055);
  g_C1->GetYaxis()->SetLabelSize(0.045);
  g_C1->SetStats(0);
  g_C1->SetTitle("");
  g_C1->Draw();
  g_C2->Draw("same");
  g_C3->Draw("same");
  g_C4->Draw("same");

  pad_l->cd();
  r_C1->GetXaxis()->SetTitle("#it{p}_{T}^{genJet} [GeV]");
  r_C1->GetYaxis()->SetTitle("Cent. / Periph.");
  r_C1->GetXaxis()->SetRangeUser(newAxis[0],newAxis[N_edge-1]);
  r_C1->SetStats(0);
  r_C1->SetTitle("");
  r_C1->GetYaxis()->SetRangeUser(0,2);
  r_C1->GetYaxis()->SetTitleSize(0.075);
  r_C1->GetYaxis()->SetLabelSize(0.065);
  r_C1->GetXaxis()->SetTitleSize(0.075);
  r_C1->GetXaxis()->SetLabelSize(0.065);
  r_C1->Draw();
  r_C2->Draw("same");
  r_C3->Draw("same");

  TLine *li = new TLine();
  li->SetLineStyle(7);
  li->DrawLine(newAxis[0],1,newAxis[N_edge-1],1);

  pad_u->cd();
  TLatex *la = new TLatex();
  la->SetTextSize(0.05);
  la->DrawLatexNDC(0.55,0.8,"HYDJET 5.02 TeV, genJets");
  la->DrawLatexNDC(0.55,0.72,"50-80%");
  la->SetTextColor(kBlue-4);
  la->DrawLatexNDC(0.55,0.64,"30-50%");
  la->SetTextColor(kGreen+2);
  la->DrawLatexNDC(0.55,0.56,"10-30%");
  la->SetTextColor(kRed-4);
  la->DrawLatexNDC(0.55,0.48,"0-10%");

  canv->SaveAs("../../../../figures/jetPt/HYDJET/jetPt_HYDJET_genJets.pdf");




  TFile *f_FJ = TFile::Open("/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/fakeJets/fakeJets.root");

  TH1D *k_C1, *k_C2, *k_C3, *k_C4; // k = fakeJets
  TH1D *j_C1, *j_C2, *j_C3, *j_C4; // j = realJets (copies)

  f_FJ->GetObject("h_fakeJets_C1",k_C1);
  f_FJ->GetObject("h_fakeJets_C2",k_C2);
  f_FJ->GetObject("h_fakeJets_C3",k_C3);
  f_FJ->GetObject("h_fakeJets_C4",k_C4);

  k_C1 = (TH1D*) k_C1->Rebin(N_edge-1,"k_C1",newAxis);
  k_C2 = (TH1D*) k_C2->Rebin(N_edge-1,"k_C2",newAxis);
  k_C3 = (TH1D*) k_C3->Rebin(N_edge-1,"k_C3",newAxis);
  k_C4 = (TH1D*) k_C4->Rebin(N_edge-1,"k_C4",newAxis);

  divideByBinwidth(k_C1);
  divideByBinwidth(k_C2);
  divideByBinwidth(k_C3);
  divideByBinwidth(k_C4);

  j_C1 = (TH1D*) g_C1->Clone("j_C1");
  j_C2 = (TH1D*) g_C2->Clone("j_C2");
  j_C3 = (TH1D*) g_C3->Clone("j_C3");
  j_C4 = (TH1D*) g_C4->Clone("j_C4");


  TLegend *leg = new TLegend(0.5,0.6,0.88,0.88);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.05);
  leg->AddEntry(j_C4,"HYDJET genJets");
  leg->AddEntry(k_C4,"PbPb fakeJets");

  TCanvas *canv_fakeJets_C4 = new TCanvas("canv_fakeJets_C4","canv_fakeJets_C4",700,700);
  canv_fakeJets_C4->cd();
  TPad *pad_fakeJets_C4_u = new TPad("pad_fakeJets_C4_u","pad_fakeJets_C4_u",0,0.4,1,1);
  TPad *pad_fakeJets_C4_l = new TPad("pad_fakeJets_C4_l","pad_fakeJets_C4_l",0,0,1,0.4);
  pad_fakeJets_C4_u->SetBottomMargin(0.0);
  pad_fakeJets_C4_l->SetTopMargin(0.0);
  pad_fakeJets_C4_l->SetBottomMargin(0.18);
  pad_fakeJets_C4_u->SetLeftMargin(0.15);
  pad_fakeJets_C4_l->SetLeftMargin(0.15);
  pad_fakeJets_C4_u->SetLogy();
  pad_fakeJets_C4_u->Draw();
  pad_fakeJets_C4_l->Draw();
  pad_fakeJets_C4_u->cd();
  j_C4->SetLineColor(kBlue+2);
  k_C4->SetLineColor(kBlue-9);
  hr_C4->SetLineColor(kRed-4);
  j_C4->SetLineWidth(2);
  k_C4->SetLineWidth(2);
  k_C4->GetYaxis()->SetTitle("1/#it{N}_{evt} d#it{N} / d#it{p}_{T} [GeV^{-1}]");
  k_C4->GetXaxis()->SetRangeUser(newAxis[0],newAxis[N_edge-1]);
  k_C4->GetYaxis()->SetTitleSize(0.055);
  k_C4->GetYaxis()->SetLabelSize(0.045);
  k_C4->SetStats(0);
  k_C4->SetTitle("PbPb FakeJets vs. HYDJET genJets 50-80%");
  k_C4->Draw();
  j_C4->Draw("same");
  leg->Draw();
  TH1D *rr_C4 = (TH1D*) j_C4->Clone("rr_C4");
  rr_C4->Divide(j_C4,k_C4,1,1,"");
  pad_fakeJets_C4_l->cd();
  rr_C4->GetXaxis()->SetTitle("#it{p}_{T}^{genJet} [GeV]");
  rr_C4->GetYaxis()->SetTitle("Real / Fake");
  rr_C4->GetXaxis()->SetRangeUser(newAxis[0],newAxis[N_edge-1]);
  rr_C4->SetStats(0);
  rr_C4->SetTitle("");
  //rr_C4->GetYaxis()->SetRangeUser(0,2);
  rr_C4->GetYaxis()->SetTitleSize(0.075);
  rr_C4->GetYaxis()->SetLabelSize(0.065);
  rr_C4->GetXaxis()->SetTitleSize(0.075);
  rr_C4->GetXaxis()->SetLabelSize(0.065);
  rr_C4->Draw();
  li->DrawLine(newAxis[0],1,newAxis[N_edge-1],1);

  canv_fakeJets_C4->SaveAs("../../../../figures/jetPt/HYDJET/fakeJet_HYDJET_C4.pdf");



  TCanvas *canv_fakeJets_C3 = new TCanvas("canv_fakeJets_C3","canv_fakeJets_C3",700,700);
  canv_fakeJets_C3->cd();
  TPad *pad_fakeJets_C3_u = new TPad("pad_fakeJets_C3_u","pad_fakeJets_C3_u",0,0.4,1,1);
  TPad *pad_fakeJets_C3_l = new TPad("pad_fakeJets_C3_l","pad_fakeJets_C3_l",0,0,1,0.4);
  pad_fakeJets_C3_u->SetBottomMargin(0.0);
  pad_fakeJets_C3_l->SetTopMargin(0.0);
  pad_fakeJets_C3_l->SetBottomMargin(0.18);
  pad_fakeJets_C3_u->SetLeftMargin(0.15);
  pad_fakeJets_C3_l->SetLeftMargin(0.15);
  pad_fakeJets_C3_u->SetLogy();
  pad_fakeJets_C3_u->Draw();
  pad_fakeJets_C3_l->Draw();
  pad_fakeJets_C3_u->cd();
  j_C3->SetLineColor(kBlue+2);
  k_C3->SetLineColor(kBlue-9);
  j_C3->SetLineWidth(2);
  k_C3->SetLineWidth(2);
  k_C3->GetYaxis()->SetTitle("1/#it{N}_{evt} d#it{N} / d#it{p}_{T} [GeV^{-1}]");
  k_C3->GetXaxis()->SetRangeUser(newAxis[0],newAxis[N_edge-1]);
  k_C3->GetYaxis()->SetTitleSize(0.055);
  k_C3->GetYaxis()->SetLabelSize(0.045);
  k_C3->SetStats(0);
  k_C3->SetTitle("PbPb FakeJets vs. HYDJET genJets 30-50%");
  k_C3->Draw();
  j_C3->Draw("same");
  leg->Draw();
  TH1D *rr_C3 = (TH1D*) j_C3->Clone("rr_C3");
  rr_C3->Divide(j_C3,k_C3,1,1,"");
  pad_fakeJets_C3_l->cd();
  rr_C3->GetXaxis()->SetTitle("#it{p}_{T}^{genJet} [GeV]");
  rr_C3->GetYaxis()->SetTitle("Real / Fake");
  rr_C3->GetXaxis()->SetRangeUser(newAxis[0],newAxis[N_edge-1]);
  rr_C3->SetStats(0);
  rr_C3->SetTitle("");
  //rr_C3->GetYaxis()->SetRangeUser(0,2);
  rr_C3->GetYaxis()->SetTitleSize(0.075);
  rr_C3->GetYaxis()->SetLabelSize(0.065);
  rr_C3->GetXaxis()->SetTitleSize(0.075);
  rr_C3->GetXaxis()->SetLabelSize(0.065);
  rr_C3->Draw();
  li->DrawLine(newAxis[0],1,newAxis[N_edge-1],1);

  canv_fakeJets_C3->SaveAs("../../../../figures/jetPt/HYDJET/fakeJet_HYDJET_C3.pdf");
  

  TCanvas *canv_fakeJets_C2 = new TCanvas("canv_fakeJets_C2","canv_fakeJets_C2",700,700);
  canv_fakeJets_C2->cd();
  TPad *pad_fakeJets_C2_u = new TPad("pad_fakeJets_C2_u","pad_fakeJets_C2_u",0,0.4,1,1);
  TPad *pad_fakeJets_C2_l = new TPad("pad_fakeJets_C2_l","pad_fakeJets_C2_l",0,0,1,0.4);
  pad_fakeJets_C2_u->SetBottomMargin(0.0);
  pad_fakeJets_C2_l->SetTopMargin(0.0);
  pad_fakeJets_C2_l->SetBottomMargin(0.18);
  pad_fakeJets_C2_u->SetLeftMargin(0.15);
  pad_fakeJets_C2_l->SetLeftMargin(0.15);
  pad_fakeJets_C2_u->SetLogy();
  pad_fakeJets_C2_u->Draw();
  pad_fakeJets_C2_l->Draw();
  pad_fakeJets_C2_u->cd();
  j_C2->SetLineColor(kBlue+2);
  k_C2->SetLineColor(kBlue-9);
  j_C2->SetLineWidth(2);
  k_C2->SetLineWidth(2);
  k_C2->GetYaxis()->SetTitle("1/#it{N}_{evt} d#it{N} / d#it{p}_{T} [GeV^{-1}]");
  k_C2->GetXaxis()->SetRangeUser(newAxis[0],newAxis[N_edge-1]);
  k_C2->GetYaxis()->SetTitleSize(0.055);
  k_C2->GetYaxis()->SetLabelSize(0.045);
  k_C2->SetStats(0);
  k_C2->SetTitle("PbPb FakeJets vs. HYDJET genJets 10-30%");
  k_C2->Draw();
  j_C2->Draw("same");
  leg->Draw();
  TH1D *rr_C2 = (TH1D*) j_C2->Clone("rr_C2");
  rr_C2->Divide(j_C2,k_C2,1,1,"");
  pad_fakeJets_C2_l->cd();
  rr_C2->GetXaxis()->SetTitle("#it{p}_{T}^{genJet} [GeV]");
  rr_C2->GetYaxis()->SetTitle("Real / Fake");
  rr_C2->GetXaxis()->SetRangeUser(newAxis[0],newAxis[N_edge-1]);
  rr_C2->SetStats(0);
  rr_C2->SetTitle("");
  //rr_C2->GetYaxis()->SetRangeUser(0,2);
  rr_C2->GetYaxis()->SetTitleSize(0.075);
  rr_C2->GetYaxis()->SetLabelSize(0.065);
  rr_C2->GetXaxis()->SetTitleSize(0.075);
  rr_C2->GetXaxis()->SetLabelSize(0.065);
  rr_C2->Draw();
  li->DrawLine(newAxis[0],1,newAxis[N_edge-1],1);

  canv_fakeJets_C2->SaveAs("../../../../figures/jetPt/HYDJET/fakeJet_HYDJET_C2.pdf");


  TCanvas *canv_fakeJets_C1 = new TCanvas("canv_fakeJets_C1","canv_fakeJets_C1",700,700);
  canv_fakeJets_C1->cd();
  TPad *pad_fakeJets_C1_u = new TPad("pad_fakeJets_C1_u","pad_fakeJets_C1_u",0,0.4,1,1);
  TPad *pad_fakeJets_C1_l = new TPad("pad_fakeJets_C1_l","pad_fakeJets_C1_l",0,0,1,0.4);
  pad_fakeJets_C1_u->SetBottomMargin(0.0);
  pad_fakeJets_C1_l->SetTopMargin(0.0);
  pad_fakeJets_C1_l->SetBottomMargin(0.18);
  pad_fakeJets_C1_u->SetLeftMargin(0.15);
  pad_fakeJets_C1_l->SetLeftMargin(0.15);
  pad_fakeJets_C1_u->SetLogy();
  pad_fakeJets_C1_u->Draw();
  pad_fakeJets_C1_l->Draw();
  pad_fakeJets_C1_u->cd();
  j_C1->SetLineColor(kBlue+2);
  k_C1->SetLineColor(kBlue-9);
  j_C1->SetLineWidth(2);
  k_C1->SetLineWidth(2);
  k_C1->GetYaxis()->SetTitle("1/#it{N}_{evt} d#it{N} / d#it{p}_{T} [GeV^{-1}]");
  k_C1->GetXaxis()->SetRangeUser(newAxis[0],newAxis[N_edge-1]);
  k_C1->GetYaxis()->SetTitleSize(0.055);
  k_C1->GetYaxis()->SetLabelSize(0.045);
  k_C1->SetStats(0);
  k_C1->SetTitle("PbPb FakeJets vs. HYDJET genJets 0-10%");
  k_C1->Draw();
  j_C1->Draw("same");
  leg->Draw();
  TH1D *rr_C1 = (TH1D*) j_C1->Clone("rr_C1");
  rr_C1->Divide(j_C1,k_C1,1,1,"");
  pad_fakeJets_C1_l->cd();
  rr_C1->GetXaxis()->SetTitle("#it{p}_{T}^{genJet} [GeV]");
  rr_C1->GetYaxis()->SetTitle("Real / Fake");
  rr_C1->GetXaxis()->SetRangeUser(newAxis[0],newAxis[N_edge-1]);
  rr_C1->SetStats(0);
  rr_C1->SetTitle("");
  //rr_C1->GetYaxis()->SetRangeUser(0,2);
  rr_C1->GetYaxis()->SetTitleSize(0.075);
  rr_C1->GetYaxis()->SetLabelSize(0.065);
  rr_C1->GetXaxis()->SetTitleSize(0.075);
  rr_C1->GetXaxis()->SetLabelSize(0.065);
  rr_C1->Draw();
  li->DrawLine(newAxis[0],1,newAxis[N_edge-1],1);

  canv_fakeJets_C1->SaveAs("../../../../figures/jetPt/HYDJET/fakeJet_HYDJET_C1.pdf");
  

}
