// Shape comparison between mixed-event and same-event FastJet anti-kT R=0.4
// jets (both on CS PF candidates), from HYDJET pfCandAnalyzer output.
// Normalizes both spectra over the visible pT range [pTlo, pThi].
// Usage: root -l -q 'plotJetPt_mixedVsSameEvent.C'

const char *inFileSame =
  "../../../rootFiles/scanningOuput/HYDJET/"
  "HYDJET_pThat-unweighted_mu12_pTmu-15to999_tight_hiBinShift-0_jetTrkMaxFilter_WDecayFilter_sameEventPseudoJets_pfCandCS_pseudoJetCandPtMin-0.0_subleadingPFCandPtMin-15_2026-6-17_ultraFineCentBins.root";

const char *inFileMixed =
  "../../../rootFiles/scanningOuput/HYDJET/"
  "HYDJET_pThat-unweighted_mu12_pTmu-15to999_tight_hiBinShift-0_jetTrkMaxFilter_WDecayFilter_mixedEventPseudoJets_pfCandCS_pseudoJetCandPtMin-0.0_subleadingPFCandPtMin-15_2026-6-17_ultraFineCentBins_partial.root";

const char *outDir = "../../../figures/jetPt/";

const double pTlo = 20.;
const double pThi = 200.;

const int NCentBins = 17;
const char *centLabel[NCentBins] = {
  "0-80%",
  "0-5%", "5-10%", "10-15%", "15-20%",
  "20-25%", "25-30%", "30-35%", "35-40%",
  "40-45%", "45-50%", "50-55%", "55-60%",
  "60-65%", "65-70%", "70-75%", "75-80%"
};

void plotOneBin(TFile *fSame, TFile *fMixed, int ci, const char *centLbl, double lm, double rm){

  TH1D *hSame = nullptr;
  fSame->GetObject(Form("h_fastJetPt_C%d", ci), hSame);
  if(!hSame){ printf("ERROR: missing h_fastJetPt_C%d in same-event file\n", ci); return; }
  hSame = (TH1D*) hSame->Clone(Form("hSame_C%d", ci));
  hSame->SetDirectory(0);

  TH1D *hMixed = nullptr;
  fMixed->GetObject(Form("h_fastJetPt_C%d", ci), hMixed);
  if(!hMixed){ printf("ERROR: missing h_fastJetPt_C%d in mixed-event file\n", ci); return; }
  hMixed = (TH1D*) hMixed->Clone(Form("hMixed_C%d", ci));
  hMixed->SetDirectory(0);

  int blo = hSame->FindBin(pTlo + 0.01);
  int bhi = hSame->FindBin(pThi - 0.01);
  hSame->Scale(1. / hSame->Integral(blo, bhi));
  hMixed->Scale(1. / hMixed->Integral(blo, bhi));

  TH1D *hRatio = (TH1D*) hMixed->Clone(Form("hRatio_C%d", ci));
  hRatio->Divide(hSame);

  TCanvas *c = new TCanvas(Form("c_C%d", ci), "", 900, 800);
  c->Divide(1, 2);

  TPad *p1 = (TPad*) c->cd(1);
  p1->SetLogy();
  p1->SetLeftMargin(lm); p1->SetRightMargin(rm);
  p1->SetBottomMargin(0.02); p1->SetTopMargin(0.08);

  hSame->SetLineColor(kBlue+1); hSame->SetLineWidth(2);
  hSame->GetXaxis()->SetRangeUser(pTlo, pThi);
  hSame->GetXaxis()->SetLabelSize(0);
  hSame->GetYaxis()->SetTitle("(1/N) dN/dp_{T} [GeV^{-1}]");
  hSame->GetYaxis()->SetTitleSize(0.06); hSame->GetYaxis()->SetLabelSize(0.055);
  hSame->GetYaxis()->SetTitleOffset(1.3);
  hSame->SetStats(0); hSame->SetTitle("");
  hSame->Draw();

  hMixed->SetLineColor(kRed+1); hMixed->SetLineWidth(2); hMixed->SetLineStyle(2);
  hMixed->Draw("same");

  TLegend *leg = new TLegend(0.48, 0.65, 0.92, 0.88);
  leg->SetBorderSize(0); leg->SetTextSize(0.057);
  leg->AddEntry(hSame,  "same-event FastJet (CS)", "l");
  leg->AddEntry(hMixed, "mixed-event FastJet (CS)", "l");
  leg->Draw();

  TLatex *lat = new TLatex(); lat->SetNDC(); lat->SetTextSize(0.055);
  lat->DrawLatex(0.22, 0.85, Form("HYDJET 5.02 TeV  |#eta^{jet}| < 1.6  cent. %s", centLbl));

  TPad *p2 = (TPad*) c->cd(2);
  p2->SetLeftMargin(lm); p2->SetRightMargin(rm);
  p2->SetTopMargin(0.02); p2->SetBottomMargin(0.18);

  hRatio->SetLineColor(kBlack); hRatio->SetLineWidth(2);
  hRatio->GetXaxis()->SetRangeUser(pTlo, pThi);
  hRatio->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
  hRatio->GetYaxis()->SetTitle("mixed / same");
  hRatio->GetYaxis()->SetRangeUser(0, 2.5);
  hRatio->GetXaxis()->SetTitleSize(0.07); hRatio->GetXaxis()->SetLabelSize(0.065);
  hRatio->GetYaxis()->SetTitleSize(0.07); hRatio->GetYaxis()->SetLabelSize(0.065);
  hRatio->GetYaxis()->SetTitleOffset(1.1);
  hRatio->GetYaxis()->SetNdivisions(505);
  hRatio->SetStats(0); hRatio->SetTitle("");
  hRatio->Draw();

  TLine *l = new TLine(pTlo, 1., pThi, 1.);
  l->SetLineStyle(2); l->SetLineColor(kGray+1); l->SetLineWidth(2); l->Draw();

  TString outFile = TString(outDir) + Form("mixedVsSame_fastJet_shape_C%d.pdf", ci);
  c->SaveAs(outFile);
  printf("Saved %s\n", outFile.Data());

  delete c;
  delete hSame;
  delete hMixed;
  delete hRatio;
}

void plotJetPt_mixedVsSameEvent(){

  TFile *fSame  = TFile::Open(inFileSame);
  if(!fSame || fSame->IsZombie()){ printf("ERROR: cannot open same-event file\n"); return; }

  TFile *fMixed = TFile::Open(inFileMixed);
  if(!fMixed || fMixed->IsZombie()){ printf("ERROR: cannot open mixed-event file\n"); return; }

  gSystem->Exec(Form("mkdir -p %s", outDir));

  double lm = 0.18, rm = 0.05;
  for(int ci = 0; ci < NCentBins; ci++){
    plotOneBin(fSame, fMixed, ci, centLabel[ci], lm, rm);
  }

  fSame->Close();
  fMixed->Close();
}
