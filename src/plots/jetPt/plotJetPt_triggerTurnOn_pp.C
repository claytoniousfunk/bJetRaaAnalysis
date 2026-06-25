// Estimate pp jet trigger turn-on curves.
// Each trigger's efficiency is measured relative to the one just below it,
// using the lower trigger's file as the reference sample (best statistics).
//   Jet30 / Jet15  from Jet15 file
//   Jet40 / Jet30  from Jet30 file
//   Jet60 / Jet40  from Jet40 file
//   Jet80 / Jet60  from Jet60 file
//   Jet100/ Jet80  from Jet80 file

void plotJetPt_triggerTurnOn_pp(){

  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);

  const char *base =
    "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/pp/"
    "pp_HighEGJet_%sHLT_mu12_pTmu-15to999_tight_deltaR-40_jetTrkMaxFilter_"
    "WDecayFilter_2026-6-16.root";

  const char *outPath =
    "/home/clayton/Analysis/code/bJetRaaAnalysis/figures/JetsPerZ/"
    "pp_triggerTurnOn.pdf";

  // {reference file trigger, numerator trigger, label, color}
  struct TurnOn { const char *refTrig; const char *numTrig; const char *label; int color; };
  TurnOn pairs[] = {
    {"Jet15", "jet30",  "Jet30 / Jet15", kAzure+1  },
    {"Jet30", "jet40",  "Jet40 / Jet30", kGreen+2  },
    {"Jet40", "jet60",  "Jet60 / Jet40", kOrange+1 },
    {"Jet60", "jet80",  "Jet80 / Jet60", kRed      },
    {"Jet80", "jet100", "Jet100 / Jet80",kViolet+1 },
  };
  const int nPairs = 5;

  TH1D *hEff[nPairs];
  for(int i = 0; i < nPairs; i++){
    TFile *f = TFile::Open(Form(base, pairs[i].refTrig));
    if(!f || f->IsZombie()){ printf("ERROR opening %s file\n", pairs[i].refTrig); return; }

    // denominator: leading jet pT for events where the reference trigger fired
    TString refLower = pairs[i].refTrig; refLower.ToLower();
    TString denName = Form("h_leadJetPt_%s", refLower.Data());
    TH1D *hD = (TH1D*) f->Get(denName);
    TH1D *hN = (TH1D*) f->Get(Form("h_leadJetPt_%s", pairs[i].numTrig));
    if(!hD || !hN){ printf("ERROR: missing histogram in %s file\n", pairs[i].refTrig); return; }
    hD->SetDirectory(nullptr); hN->SetDirectory(nullptr);

    hEff[i] = (TH1D*) hN->Clone(Form("hEff_%d", i));
    hEff[i]->SetDirectory(nullptr);
    hEff[i]->Divide(hN, hD, 1, 1, "B");
    hEff[i]->SetLineColor(pairs[i].color); hEff[i]->SetLineWidth(2);
    hEff[i]->SetMarkerColor(pairs[i].color); hEff[i]->SetMarkerStyle(20);
    hEff[i]->SetMarkerSize(0.6);
    f->Close();
  }

  TCanvas *c = new TCanvas("c","",700,600);
  c->SetLeftMargin(0.14); c->SetBottomMargin(0.14);
  c->SetRightMargin(0.05); c->SetTopMargin(0.06);
  c->SetLogx();

  TH1D *hFrame = (TH1D*) hEff[0]->Clone("hFrame");
  hFrame->Reset();
  hFrame->GetXaxis()->SetRangeUser(10., 300.);
  hFrame->GetXaxis()->SetTitle("Leading jet #it{p}_{T} [GeV]");
  hFrame->GetXaxis()->SetTitleSize(0.048); hFrame->GetXaxis()->SetTitleOffset(1.1);
  hFrame->GetXaxis()->SetLabelSize(0.040);
  hFrame->GetYaxis()->SetTitle("Trigger efficiency (relative to next lower)");
  hFrame->GetYaxis()->SetTitleSize(0.040); hFrame->GetYaxis()->SetTitleOffset(1.6);
  hFrame->GetYaxis()->SetLabelSize(0.040);
  hFrame->SetMinimum(0.); hFrame->SetMaximum(1.3);
  hFrame->SetTitle("");
  hFrame->Draw("axis");

  TLine *line1 = new TLine(10.,1.,300.,1.);
  line1->SetLineStyle(2); line1->SetLineColor(kGray+1); line1->Draw();

  for(int i = 0; i < nPairs; i++) hEff[i]->Draw("ep same");

  TLegend *leg = new TLegend(0.55,0.65,0.93,0.90);
  leg->SetBorderSize(0); leg->SetTextSize(0.038);
  leg->SetFillStyle(0); leg->SetFillColorAlpha(0, 0.);
  for(int i = 0; i < nPairs; i++) leg->AddEntry(hEff[i], pairs[i].label, "lp");
  leg->Draw();

  TLatex lat; lat.SetNDC(); lat.SetTextSize(0.038);
  lat.DrawLatex(0.16, 0.92, "pp HighEGJet, each trigger relative to next lower");

  c->SaveAs(outPath);
  printf("Saved %s\n", outPath);
}
