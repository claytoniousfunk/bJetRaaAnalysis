#include "../../../headers/functions/divideByBinwidth.h"

// ---- stitch point pT boundaries ----------------------------------------
const double pT_15to30  = 55.;
const double pT_30to40  = 70.;
const double pT_40to60  = 100.;
const double pT_60to80  = 140.;
const double pT_80to100 = 175.;

// ---- colours for each trigger ------------------------------------------
const int col15  = kBlack;
const int col30  = kRed+1;
const int col40  = kOrange+1;
const int col60  = kGreen+2;
const int col80  = kAzure-3;
const int col100 = kViolet+1;

void setStyle(TH1D *h, int col){
  h->SetLineColor(col);
  h->SetMarkerColor(col);
  h->SetLineWidth(2);
  h->SetStats(0);
}

// Normalise h_low to h_high in the pT range [lo, hi]; returns scale factor
double scaleToMatch(TH1D *h_low, TH1D *h_high, double lo, double hi){
  double s = 0.01; // small shift inside bin
  double N_high = h_high->Integral(h_high->FindBin(lo+s), h_high->FindBin(hi-s));
  double N_low  = h_low ->Integral(h_low ->FindBin(lo+s), h_low ->FindBin(hi-s));
  if(N_low == 0) return 1.;
  double factor = N_high / N_low;
  h_low->Scale(factor);
  return factor;
}

// Build the stitched histogram from 6 scaled inputs
TH1D* buildStitch(TH1D *h15, TH1D *h30, TH1D *h40,
                  TH1D *h60, TH1D *h80, TH1D *h100){
  TH1D *h = (TH1D*) h100->Clone("h_stitched");
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


void plotJetPt_stitch_pp(){

  TString outDir = "/home/clayton/Documents/nuclear/GroupMeeting/figures/2026-06-16/";

  const char *rbase = "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/pp/pp_HighEGJet_";
  const char *rsuf  = "HLT_mu12_pTmu-15to999_tight_deltaR-40_jetTrkMaxFilter_WDecayFilter_2026-6-16.root";

  TFile *f15  = TFile::Open(Form("%sJet15%s",  rbase, rsuf));
  TFile *f30  = TFile::Open(Form("%sJet30%s",  rbase, rsuf));
  TFile *f40  = TFile::Open(Form("%sJet40%s",  rbase, rsuf));
  TFile *f60  = TFile::Open(Form("%sJet60%s",  rbase, rsuf));
  TFile *f80  = TFile::Open(Form("%sJet80%s",  rbase, rsuf));
  TFile *f100 = TFile::Open(Form("%sJet100%s", rbase, rsuf));

  // raw spectra
  TH1D *hR15, *hR30, *hR40, *hR60, *hR80, *hR100;
  f15 ->GetObject("h_inclRecoJetPt", hR15);
  f30 ->GetObject("h_inclRecoJetPt", hR30);
  f40 ->GetObject("h_inclRecoJetPt", hR40);
  f60 ->GetObject("h_inclRecoJetPt", hR60);
  f80 ->GetObject("h_inclRecoJetPt", hR80);
  f100->GetObject("h_inclRecoJetPt", hR100);

  // turn-on: leading jet pT for each fired trigger, from the Jet15 file
  TH1D *hTO15, *hTO30, *hTO40, *hTO60, *hTO80, *hTO100;
  f15->GetObject("h_leadJetPt_jet15",  hTO15);
  f15->GetObject("h_leadJetPt_jet30",  hTO30);
  f15->GetObject("h_leadJetPt_jet40",  hTO40);
  f15->GetObject("h_leadJetPt_jet60",  hTO60);
  f15->GetObject("h_leadJetPt_jet80",  hTO80);
  f15->GetObject("h_leadJetPt_jet100", hTO100);

  setStyle(hR15,  col15);  setStyle(hR30,  col30);
  setStyle(hR40,  col40);  setStyle(hR60,  col60);
  setStyle(hR80,  col80);  setStyle(hR100, col100);
  setStyle(hTO15, col15);  setStyle(hTO30, col30);
  setStyle(hTO40, col40);  setStyle(hTO60, col60);
  setStyle(hTO80, col80);  setStyle(hTO100,col100);

  // working copies for scaled spectra
  TH1D *h15  = (TH1D*) hR15 ->Clone("h15_sc");
  TH1D *h30  = (TH1D*) hR30 ->Clone("h30_sc");
  TH1D *h40  = (TH1D*) hR40 ->Clone("h40_sc");
  TH1D *h60  = (TH1D*) hR60 ->Clone("h60_sc");
  TH1D *h80  = (TH1D*) hR80 ->Clone("h80_sc");
  TH1D *h100 = (TH1D*) hR100->Clone("h100_sc");

  setStyle(h15, col15);  setStyle(h30, col30);
  setStyle(h40, col40);  setStyle(h60, col60);
  setStyle(h80, col80);  setStyle(h100,col100);

  // ---- stitching scale chain (Jet100 is reference) --------------------
  double sf80  = scaleToMatch(h80,  h100, pT_80to100, 250.);
  double sf60  = scaleToMatch(h60,  h80,  pT_60to80,  pT_80to100);
  double sf40  = scaleToMatch(h40,  h60,  pT_40to60,  pT_60to80);
  double sf30  = scaleToMatch(h30,  h40,  pT_30to40,  pT_40to60);
  double sf15  = scaleToMatch(h15,  h30,  pT_15to30,  pT_30to40);

  printf("Scale factors:\n");
  printf("  Jet80  × %.4f\n", sf80);
  printf("  Jet60  × %.4f\n", sf60);
  printf("  Jet40  × %.4f\n", sf40);
  printf("  Jet30  × %.4f\n", sf30);
  printf("  Jet15  × %.4f\n", sf15);

  TH1D *h_stitch = buildStitch(h15, h30, h40, h60, h80, h100);
  setStyle(h_stitch, kBlack);

  // ---- legend ---------------------------------------------------------
  auto makeLeg = [](double x1, double y1, double x2, double y2) -> TLegend* {
    TLegend *l = new TLegend(x1,y1,x2,y2);
    l->SetBorderSize(0);
    l->SetTextSize(0.040);
    return l;
  };

  // common axis limits
  double xmin = 15., xmax = 400.;
  double ymin_raw = 1., ymax_raw = 5.e7;
  double lm = 0.16, bm = 0.14, tm = 0.12;

  TLine *li = new TLine();
  li->SetLineStyle(2);
  li->SetLineWidth(2);

  // =====================================================================
  // PLOT 1: Baseline — all raw spectra overlaid
  // =====================================================================
  {
    TLegend *leg = makeLeg(0.63,0.60,0.88,0.88);
    leg->AddEntry(hR15,  "Jet15",  "l");
    leg->AddEntry(hR30,  "Jet30",  "l");
    leg->AddEntry(hR40,  "Jet40",  "l");
    leg->AddEntry(hR60,  "Jet60",  "l");
    leg->AddEntry(hR80,  "Jet80",  "l");
    leg->AddEntry(hR100, "Jet100", "l");

    TCanvas *c = new TCanvas("c_baseline","",700,700);
    TPad *p = new TPad("p_baseline","",0,0,1,1);
    p->SetLeftMargin(lm); p->SetBottomMargin(bm); p->SetTopMargin(tm);
    p->SetLogy(); p->Draw(); p->cd();

    hR15->SetTitle("");
    hR15->GetXaxis()->SetRangeUser(xmin,xmax);
    hR15->GetYaxis()->SetRangeUser(ymin_raw,ymax_raw);
    hR15->GetXaxis()->SetTitle("Jet #it{p}_{T} [GeV]");
    hR15->GetYaxis()->SetTitle("Entries");
    hR15->Draw(); hR30->Draw("same"); hR40->Draw("same");
    hR60->Draw("same"); hR80->Draw("same"); hR100->Draw("same");
    leg->Draw();

    TLatex *lat = new TLatex();
    lat->SetNDC(); lat->SetTextSize(0.038);
    lat->DrawLatex(lm+0.02, 1.-tm+0.03, "pp #sqrt{#it{s}} = 5.02 TeV");
    lat->DrawLatex(lm+0.02, 1.-tm-0.02, "anti-#it{k}_{T} #it{R} = 0.4");

    c->SaveAs(outDir + "jetPt_pp_baseline.pdf");
  }

  // =====================================================================
  // PLOT 2: Turn-on curves + ratio pad (trigger / Jet15 per bin)
  // =====================================================================
  {
    double n15 = hTO15->Integral();
    if(n15 > 0){
      hTO15 ->Scale(1./n15); hTO30 ->Scale(1./n15);
      hTO40 ->Scale(1./n15); hTO60 ->Scale(1./n15);
      hTO80 ->Scale(1./n15); hTO100->Scale(1./n15);
    }

    // ratio histograms: trigger / Jet15 bin-by-bin
    auto makeRatio = [&](TH1D *hNum, const char *name) -> TH1D* {
      TH1D *hr = (TH1D*) hNum->Clone(name);
      hr->Divide(hTO15);
      return hr;
    };
    TH1D *hR30  = makeRatio(hTO30,  "rTO30");
    TH1D *hR40  = makeRatio(hTO40,  "rTO40");
    TH1D *hR60  = makeRatio(hTO60,  "rTO60");
    TH1D *hR80  = makeRatio(hTO80,  "rTO80");
    TH1D *hR100 = makeRatio(hTO100, "rTO100");

    TLegend *leg = makeLeg(0.55,0.38,0.88,0.73);
    leg->AddEntry(hTO15,  "Jet15 (ref)", "l");
    leg->AddEntry(hTO30,  "Jet30",       "l");
    leg->AddEntry(hTO40,  "Jet40",       "l");
    leg->AddEntry(hTO60,  "Jet60",       "l");
    leg->AddEntry(hTO80,  "Jet80",       "l");
    leg->AddEntry(hTO100, "Jet100",      "l");

    TCanvas *c = new TCanvas("c_turnon","",700,900);
    c->cd();

    TPad *pUp = new TPad("pUp","",0,0.35,1,1);
    pUp->SetLeftMargin(lm); pUp->SetRightMargin(0.05);
    pUp->SetBottomMargin(0.); pUp->SetTopMargin(tm);
    pUp->SetLogy(); pUp->Draw(); pUp->cd();

    hTO15->SetTitle("");
    hTO15->GetXaxis()->SetRangeUser(xmin, xmax);
    hTO15->GetYaxis()->SetTitle("Fraction of Jet15 events");
    hTO15->GetYaxis()->SetTitleSize(0.052);
    hTO15->GetYaxis()->SetLabelSize(0.048);
    hTO15->GetYaxis()->SetTitleOffset(1.35);
    hTO15->Draw(); hTO30->Draw("same"); hTO40->Draw("same");
    hTO60->Draw("same"); hTO80->Draw("same"); hTO100->Draw("same");
    leg->Draw();

    TLatex *lat = new TLatex();
    lat->SetNDC(); lat->SetTextSize(0.050);
    lat->DrawLatex(lm+0.02, 1.-tm+0.03, "pp #sqrt{#it{s}} = 5.02 TeV");
    lat->DrawLatex(lm+0.02, 1.-tm-0.05, "Jet15-triggered events");

    li->SetLineColor(kGray+1); li->SetLineStyle(2); li->SetLineWidth(1);
    li->DrawLine(pT_15to30,  1.e-6, pT_15to30,  1.);
    li->DrawLine(pT_30to40,  1.e-6, pT_30to40,  1.);
    li->DrawLine(pT_40to60,  1.e-6, pT_40to60,  1.);
    li->DrawLine(pT_60to80,  1.e-6, pT_60to80,  1.);
    li->DrawLine(pT_80to100, 1.e-6, pT_80to100, 1.);

    c->cd();
    TPad *pDn = new TPad("pDn","",0,0,1,0.35);
    pDn->SetLeftMargin(lm); pDn->SetRightMargin(0.05);
    pDn->SetTopMargin(0.); pDn->SetBottomMargin(0.28);
    pDn->Draw(); pDn->cd();

    TLine *unity = new TLine();
    unity->SetLineStyle(2); unity->SetLineColor(kGray+1); unity->SetLineWidth(1);

    hR30->SetTitle("");
    hR30->GetXaxis()->SetRangeUser(xmin, xmax);
    hR30->GetXaxis()->SetTitle("Leading Jet #it{p}_{T} [GeV]");
    hR30->GetYaxis()->SetTitle("Trigger / Jet15");
    hR30->GetYaxis()->SetRangeUser(0., 1.4);
    hR30->GetXaxis()->SetTitleSize(0.105);
    hR30->GetXaxis()->SetLabelSize(0.090);
    hR30->GetYaxis()->SetTitleSize(0.095);
    hR30->GetYaxis()->SetLabelSize(0.085);
    hR30->GetYaxis()->SetTitleOffset(0.65);
    hR30->GetYaxis()->SetNdivisions(504);
    hR30->Draw("hist"); hR40->Draw("hist same"); hR60->Draw("hist same");
    hR80->Draw("hist same"); hR100->Draw("hist same");
    unity->DrawLine(xmin, 1., xmax, 1.);

    li->DrawLine(pT_15to30,  0., pT_15to30,  1.4);
    li->DrawLine(pT_30to40,  0., pT_30to40,  1.4);
    li->DrawLine(pT_40to60,  0., pT_40to60,  1.4);
    li->DrawLine(pT_60to80,  0., pT_60to80,  1.4);
    li->DrawLine(pT_80to100, 0., pT_80to100, 1.4);

    c->SaveAs(outDir + "jetPt_pp_turnon.pdf");
  }

  // =====================================================================
  // PLOT 3: Stitch procedure — scaled samples + stitch lines
  // =====================================================================
  {
    TLegend *leg = makeLeg(0.63,0.60,0.88,0.88);
    leg->AddEntry(h15,  "Jet15 (scaled)", "l");
    leg->AddEntry(h30,  "Jet30 (scaled)", "l");
    leg->AddEntry(h40,  "Jet40 (scaled)", "l");
    leg->AddEntry(h60,  "Jet60 (scaled)", "l");
    leg->AddEntry(h80,  "Jet80 (scaled)", "l");
    leg->AddEntry(h100, "Jet100",         "l");

    TCanvas *c = new TCanvas("c_proc","",700,700);
    TPad *p = new TPad("p_proc","",0,0,1,1);
    p->SetLeftMargin(lm); p->SetBottomMargin(bm); p->SetTopMargin(tm);
    p->SetLogy(); p->Draw(); p->cd();

    h100->SetTitle("");
    h100->GetXaxis()->SetRangeUser(xmin,xmax);
    h100->GetYaxis()->SetRangeUser(1.,ymax_raw*50.);
    h100->GetXaxis()->SetTitle("Jet #it{p}_{T} [GeV]");
    h100->GetYaxis()->SetTitle("Scaled Entries");
    h100->Draw(); h80->Draw("same"); h60->Draw("same");
    h40->Draw("same"); h30->Draw("same"); h15->Draw("same");
    leg->Draw();

    // vertical stitch lines
    li->SetLineColor(kGray+2);
    li->SetLineStyle(2); li->SetLineWidth(1);
    li->DrawLine(pT_15to30,  1., pT_15to30,  ymax_raw*50.);
    li->DrawLine(pT_30to40,  1., pT_30to40,  ymax_raw*50.);
    li->DrawLine(pT_40to60,  1., pT_40to60,  ymax_raw*50.);
    li->DrawLine(pT_60to80,  1., pT_60to80,  ymax_raw*50.);
    li->DrawLine(pT_80to100, 1., pT_80to100, ymax_raw*50.);

    TLatex *lat = new TLatex();
    lat->SetNDC(); lat->SetTextSize(0.038);
    lat->DrawLatex(lm+0.02, 1.-tm+0.03, "pp #sqrt{#it{s}} = 5.02 TeV");
    lat->DrawLatex(lm+0.02, 1.-tm-0.02, "anti-#it{k}_{T} #it{R} = 0.4");

    // label stitch pT on the line
    lat->SetNDC(false); lat->SetTextSize(14); lat->SetTextFont(43);
    lat->SetTextAngle(90); lat->SetTextColor(kGray+2);
    lat->DrawLatex(pT_15to30  + 2., 5., Form("%.0f GeV",pT_15to30));
    lat->DrawLatex(pT_30to40  + 2., 5., Form("%.0f GeV",pT_30to40));
    lat->DrawLatex(pT_40to60  + 2., 5., Form("%.0f GeV",pT_40to60));
    lat->DrawLatex(pT_60to80  + 2., 5., Form("%.0f GeV",pT_60to80));
    lat->DrawLatex(pT_80to100 + 2., 5., Form("%.0f GeV",pT_80to100));

    c->SaveAs(outDir + "jetPt_pp_stitch_procedure.pdf");
  }

  // =====================================================================
  // PLOT 4: Final stitched spectrum (÷ bin width)
  // =====================================================================
  {
    TH1D *h_final = (TH1D*) h_stitch->Clone("h_final");
    divideByBinwidth(h_final);

    TCanvas *c = new TCanvas("c_stitch","",700,700);
    TPad *p = new TPad("p_stitch","",0,0,1,1);
    p->SetLeftMargin(lm); p->SetBottomMargin(bm); p->SetTopMargin(tm);
    p->SetLogy(); p->Draw(); p->cd();

    h_final->SetTitle("");
    h_final->GetXaxis()->SetRangeUser(xmin,xmax);
    h_final->GetXaxis()->SetTitle("Jet #it{p}_{T} [GeV]");
    h_final->GetYaxis()->SetTitle("dN/d#it{p}_{T} [GeV^{-1}]");
    h_final->SetLineWidth(2);
    h_final->Draw();

    TLatex *lat = new TLatex();
    lat->SetNDC(); lat->SetTextSize(0.038);
    lat->DrawLatex(lm+0.02, 1.-tm+0.03, "pp #sqrt{#it{s}} = 5.02 TeV");
    lat->DrawLatex(lm+0.02, 1.-tm-0.02, "anti-#it{k}_{T} #it{R} = 0.4");
    lat->DrawLatex(lm+0.02, 1.-tm-0.07, "Stitched: Jet15 #oplus Jet30 #oplus Jet40 #oplus Jet60 #oplus Jet80 #oplus Jet100");

    c->SaveAs(outDir + "jetPt_pp_stitched.pdf");
    c->SaveAs("/home/clayton/Analysis/code/bJetRaaAnalysis/figures/jetPt/stitch_pp/jetPt_pp_stitched.pdf");
  }

  printf("\nDone. Output in %s\n", outDir.Data());
}
