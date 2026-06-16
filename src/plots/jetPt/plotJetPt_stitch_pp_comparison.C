#include "../../../headers/functions/divideByBinwidth.h"

// ---- new stitch (6-trigger) boundaries ---------------------------------
const double pT_15to30  = 55.;
const double pT_30to40  = 70.;
const double pT_40to60  = 100.;
const double pT_60to80  = 140.;
const double pT_80to100 = 175.;

// ---- old stitch (MinBias-anchored) boundaries --------------------------
const double pT_MB_to60  = 100.;
const double pT_60_to80  = 150.;
const double pT_80_to100 = 200.;

void setStyle(TH1D *h, int col, int style = 1){
  h->SetLineColor(col);
  h->SetMarkerColor(col);
  h->SetLineWidth(2);
  h->SetLineStyle(style);
  h->SetStats(0);
}

double scaleToMatch(TH1D *h_low, TH1D *h_high, double lo, double hi){
  double s = 0.01;
  double N_high = h_high->Integral(h_high->FindBin(lo+s), h_high->FindBin(hi-s));
  double N_low  = h_low ->Integral(h_low ->FindBin(lo+s), h_low ->FindBin(hi-s));
  if(N_low == 0) return 1.;
  h_low->Scale(N_high / N_low);
  return N_high / N_low;
}

// ---- build new stitch (Jet15→30→40→60→80→100) -------------------------
TH1D* buildNewStitch(TH1D *h15, TH1D *h30, TH1D *h40,
                     TH1D *h60, TH1D *h80, TH1D *h100){
  TH1D *h = (TH1D*) h100->Clone("h_new_stitch");
  for(int i = 1; i <= h->GetNbinsX(); i++){
    double pT = h->GetBinCenter(i);
    TH1D *src = nullptr;
    if     (pT < pT_15to30 ) src = h15;
    else if(pT < pT_30to40 ) src = h30;
    else if(pT < pT_40to60 ) src = h40;
    else if(pT < pT_60to80 ) src = h60;
    else if(pT < pT_80to100) src = h80;
    else                      src = h100;
    h->SetBinContent(i, src->GetBinContent(i));
    h->SetBinError  (i, src->GetBinError  (i));
  }
  return h;
}

// ---- build old stitch (MinBias→Jet60→Jet80→Jet100) --------------------
TH1D* buildOldStitch(TH1D *hMB, TH1D *h60, TH1D *h80, TH1D *h100){
  TH1D *h = (TH1D*) h100->Clone("h_old_stitch");
  for(int i = 1; i <= h->GetNbinsX(); i++){
    double pT = h->GetBinCenter(i);
    TH1D *src = nullptr;
    if     (pT < pT_MB_to60 ) src = hMB;
    else if(pT < pT_60_to80 ) src = h60;
    else if(pT < pT_80_to100) src = h80;
    else                       src = h100;
    h->SetBinContent(i, src->GetBinContent(i));
    h->SetBinError  (i, src->GetBinError  (i));
  }
  return h;
}


void plotJetPt_stitch_pp_comparison(){

  TString outDir  = "/home/clayton/Documents/nuclear/GroupMeeting/figures/2026-06-16/";
  const char *rbase = "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/pp/pp_HighEGJet_";
  const char *rsuf  = "HLT_mu12_pTmu-15to999_tight_deltaR-40_jetTrkMaxFilter_WDecayFilter_2026-6-16.root";

  // ---- new scan files --------------------------------------------------
  TFile *f15  = TFile::Open(Form("%sJet15%s",  rbase, rsuf));
  TFile *f30  = TFile::Open(Form("%sJet30%s",  rbase, rsuf));
  TFile *f40  = TFile::Open(Form("%sJet40%s",  rbase, rsuf));
  TFile *f60  = TFile::Open(Form("%sJet60%s",  rbase, rsuf));
  TFile *f80  = TFile::Open(Form("%sJet80%s",  rbase, rsuf));
  TFile *f100 = TFile::Open(Form("%sJet100%s", rbase, rsuf));

  // ---- MinBias file (old anchor) ---------------------------------------
  TFile *fMB = TFile::Open("/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/scanningOutput/pp/latest/pp_MinBias_mu12_pTmu-15to999_tight_deltaR-40_jetTrkMaxFilter_WDecayFilter_2026-4-7.root");

  TH1D *hR15, *hR30, *hR40, *hR60, *hR80, *hR100, *hRMB;
  f15 ->GetObject("h_inclRecoJetPt", hR15);
  f30 ->GetObject("h_inclRecoJetPt", hR30);
  f40 ->GetObject("h_inclRecoJetPt", hR40);
  f60 ->GetObject("h_inclRecoJetPt", hR60);
  f80 ->GetObject("h_inclRecoJetPt", hR80);
  f100->GetObject("h_inclRecoJetPt", hR100);
  fMB ->GetObject("h_inclRecoJetPt", hRMB);

  // ---- new stitch working copies ---------------------------------------
  TH1D *n15  = (TH1D*) hR15 ->Clone("n15");
  TH1D *n30  = (TH1D*) hR30 ->Clone("n30");
  TH1D *n40  = (TH1D*) hR40 ->Clone("n40");
  TH1D *n60  = (TH1D*) hR60 ->Clone("n60");
  TH1D *n80  = (TH1D*) hR80 ->Clone("n80");
  TH1D *n100 = (TH1D*) hR100->Clone("n100");

  double sf80  = scaleToMatch(n80,  n100, pT_80to100, 250.);
  double sf60  = scaleToMatch(n60,  n80,  pT_60to80,  pT_80to100);
  double sf40  = scaleToMatch(n40,  n60,  pT_40to60,  pT_60to80);
  double sf30  = scaleToMatch(n30,  n40,  pT_30to40,  pT_40to60);
  double sf15  = scaleToMatch(n15,  n30,  pT_15to30,  pT_30to40);

  printf("New stitch scale factors:\n");
  printf("  Jet80×%.4f  Jet60×%.4f  Jet40×%.4f  Jet30×%.4f  Jet15×%.4f\n",
         sf80,sf60,sf40,sf30,sf15);

  TH1D *h_new = buildNewStitch(n15,n30,n40,n60,n80,n100);

  // ---- old stitch working copies --------------------------------------
  TH1D *o60  = (TH1D*) hR60 ->Clone("o60");
  TH1D *o80  = (TH1D*) hR80 ->Clone("o80");
  TH1D *o100 = (TH1D*) hR100->Clone("o100");
  TH1D *oMB  = (TH1D*) hRMB ->Clone("oMB");

  double osf80 = scaleToMatch(o80,  o100, pT_80_to100, 250.);
  double osf60 = scaleToMatch(o60,  o80,  pT_60_to80,  pT_80_to100);
  double osfMB = scaleToMatch(oMB,  o60,  pT_MB_to60,  pT_60_to80);

  printf("Old stitch scale factors:\n");
  printf("  Jet80×%.4f  Jet60×%.4f  MinBias×%.4f\n", osf80,osf60,osfMB);

  TH1D *h_old = buildOldStitch(oMB,o60,o80,o100);

  // ---- divide by bin width --------------------------------------------
  TH1D *h_new_bw = (TH1D*) h_new->Clone("h_new_bw");
  TH1D *h_old_bw = (TH1D*) h_old->Clone("h_old_bw");
  divideByBinwidth(h_new_bw);
  divideByBinwidth(h_old_bw);

  setStyle(h_new_bw, kBlack,   1);
  setStyle(h_old_bw, kRed+1,   2);

  // ---- ratio ----------------------------------------------------------
  TH1D *h_ratio = (TH1D*) h_new_bw->Clone("h_ratio");
  h_ratio->Divide(h_old_bw);
  setStyle(h_ratio, kBlack, 1);

  // ---- canvas layout --------------------------------------------------
  double lm = 0.16, bm = 0.14, tm = 0.12;
  double xmin = 15., xmax = 400.;

  TLegend *leg = new TLegend(0.42,0.65,0.88,0.82);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.042);
  leg->AddEntry(h_new_bw, "New: Jet15 #oplus Jet30 #oplus Jet40 #oplus Jet60 #oplus Jet80 #oplus Jet100", "l");
  leg->AddEntry(h_old_bw, "Old: MinBias #oplus Jet60 #oplus Jet80 #oplus Jet100", "l");

  TLine *li_unity = new TLine();
  li_unity->SetLineStyle(2); li_unity->SetLineColor(kGray+1); li_unity->SetLineWidth(1);

  // ---- draw -----------------------------------------------------------
  TCanvas *c = new TCanvas("c_compare","",700,800);
  c->cd();
  TPad *pUp = new TPad("pUp","",0,0.35,1,1);
  TPad *pDn = new TPad("pDn","",0,0,1,0.35);

  pUp->SetLeftMargin(lm); pUp->SetBottomMargin(0.);
  pUp->SetTopMargin(tm);  pUp->SetRightMargin(0.05);
  pUp->SetLogy(); pUp->Draw(); pUp->cd();

  h_new_bw->SetTitle("");
  h_new_bw->GetXaxis()->SetRangeUser(xmin,xmax);
  h_new_bw->GetYaxis()->SetTitle("dN/d#it{p}_{T} [GeV^{-1}]");
  h_new_bw->GetYaxis()->SetTitleSize(0.055);
  h_new_bw->GetYaxis()->SetLabelSize(0.050);
  h_new_bw->Draw();
  h_old_bw->Draw("same");
  leg->Draw();

  TLatex *lat = new TLatex();
  lat->SetNDC(); lat->SetTextSize(0.048);
  lat->DrawLatex(lm+0.02, 1.-tm+0.03, "pp #sqrt{#it{s}} = 5.02 TeV");
  lat->DrawLatex(lm+0.02, 1.-tm-0.03, "anti-#it{k}_{T} #it{R} = 0.4");

  c->cd();
  pDn->SetLeftMargin(lm); pDn->SetTopMargin(0.);
  pDn->SetBottomMargin(0.30); pDn->SetRightMargin(0.05);
  pDn->Draw(); pDn->cd();

  h_ratio->SetTitle("");
  h_ratio->GetXaxis()->SetRangeUser(xmin,xmax);
  h_ratio->GetXaxis()->SetTitle("Jet #it{p}_{T} [GeV]");
  h_ratio->GetYaxis()->SetTitle("New / Old");
  h_ratio->GetYaxis()->SetRangeUser(0.5, 1.5);
  h_ratio->GetXaxis()->SetTitleSize(0.12);
  h_ratio->GetXaxis()->SetLabelSize(0.10);
  h_ratio->GetYaxis()->SetTitleSize(0.10);
  h_ratio->GetYaxis()->SetLabelSize(0.09);
  h_ratio->GetYaxis()->SetTitleOffset(0.6);
  h_ratio->GetYaxis()->SetNdivisions(505);
  h_ratio->Draw();
  li_unity->DrawLine(xmin,1.,xmax,1.);

  c->SaveAs(outDir + "jetPt_pp_stitch_comparison.pdf");
  c->SaveAs("/home/clayton/Analysis/code/bJetRaaAnalysis/figures/jetPt/stitch_pp/jetPt_pp_stitch_comparison.pdf");

  printf("\nDone.\n");
}
