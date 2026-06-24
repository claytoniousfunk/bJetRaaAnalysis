// Fake-jet subtraction: h_signal = h_same - h_mixed, both normalized by the
// sum of event weights (proxy: h_pseudoJetPt integral / N_generatedPseudoJets).
// Produces per-centrality plots of same-event, mixed-event, and signal spectra,
// plus the signal fraction (signal / same-event).
// Usage: root -l -q 'plotJetPt_fakeJetSubtraction.C'

const char *inFileSame =
  "../../../rootFiles/scanningOuput/HYDJET/"
  "HYDJET_pThat-unweighted_mu12_pTmu-15to999_tight_hiBinShift-0_jetTrkMaxFilter_WDecayFilter_sameEventPseudoJets_pfCandCS_pseudoJetCandPtMin-0.0_subleadingPFCandPtMin-15_2026-6-18_ultraFineCentBins.root";

const char *inFileMixed =
  "../../../rootFiles/scanningOuput/HYDJET/"
  "HYDJET_pThat-unweighted_mu12_pTmu-15to999_tight_hiBinShift-0_jetTrkMaxFilter_WDecayFilter_mixedEventPseudoJets_pfCandCS_pseudoJetCandPtMin-0.0_subleadingPFCandPtMin-15_2026-6-23_ultraFineCentBins.root";

const char *outDir = "../../../figures/jetPt/";

const double pTlo = 0.;
const double pThi = 200.;
const int N_genPJ = 100;  // N_generatedPseudoJets

const int NCentBins = 17;
const char *centLabel[NCentBins] = {
  "0-80%",
  "0-5%", "5-10%", "10-15%", "15-20%",
  "20-25%", "25-30%", "30-35%", "35-40%",
  "40-45%", "45-50%", "50-55%", "55-60%",
  "60-65%", "65-70%", "70-75%", "75-80%"
};

void plotOneBin(TFile *fSame, TFile *fMixed, int ci, const char *centLbl,
                double lm, double rm){

  // --- get fastJetPt histograms ---
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

  // --- compute sum-of-weights per centrality via pseudoJet proxy ---
  TH1D *hPJS = nullptr, *hPJM = nullptr;
  fSame->GetObject(Form("h_pseudoJetPt_C%d", ci), hPJS);
  fMixed->GetObject(Form("h_pseudoJetPt_C%d", ci), hPJM);
  if(!hPJS || !hPJM){
    printf("ERROR: missing h_pseudoJetPt_C%d\n", ci); return;
  }
  double sumW_same  = hPJS->Integral() / N_genPJ;
  double sumW_mixed = hPJM->Integral() / N_genPJ;
  if(sumW_same  <= 0 || sumW_mixed <= 0){
    printf("WARNING: zero sumW for C%d, skipping\n", ci); return;
  }

  // normalize to jets per unit event weight
  hSame ->Scale(1. / sumW_same);
  hMixed->Scale(1. / sumW_mixed);

  // --- signal = same - mixed ---
  TH1D *hSignal = (TH1D*) hSame->Clone(Form("hSignal_C%d", ci));
  hSignal->Add(hMixed, -1.);
  hSignal->SetDirectory(0);

  // --- signal fraction = signal / same ---
  TH1D *hFrac = (TH1D*) hSignal->Clone(Form("hFrac_C%d", ci));
  hFrac->Divide(hSame);
  hFrac->SetDirectory(0);

  // Okabe-Ito colorblind-safe palette
  Int_t col_blue   = TColor::GetColor("#0072B2");
  Int_t col_orange = TColor::GetColor("#E69F00");
  Int_t col_green  = TColor::GetColor("#009E73");

  TCanvas *c = new TCanvas(Form("c_C%d", ci), "", 900, 800);
  c->Divide(1, 2);

  // --- top pad: spectra ---
  TPad *p1 = (TPad*) c->cd(1);
  p1->SetLogy();
  p1->SetLeftMargin(lm); p1->SetRightMargin(rm);
  p1->SetBottomMargin(0.02); p1->SetTopMargin(0.08);

  // find y-axis range across all three histograms
  int blo = hSame->FindBin(pTlo + 0.01);
  int bhi = hSame->FindBin(pThi - 0.01);
  double ymax = 0;
  for(int b = blo; b <= bhi; b++){
    ymax = TMath::Max(ymax, hSame->GetBinContent(b));
    ymax = TMath::Max(ymax, hMixed->GetBinContent(b));
  }
  ymax *= 3.;
  double ymin = 1e-10;

  hSame->SetLineColor(col_blue);   hSame->SetLineWidth(2); hSame->SetLineStyle(1);
  hSame->GetXaxis()->SetRangeUser(pTlo, pThi);
  hSame->GetXaxis()->SetLabelSize(0);
  hSame->GetYaxis()->SetTitle("dN/dp_{T} per event [GeV^{-1}]");
  hSame->GetYaxis()->SetTitleSize(0.055); hSame->GetYaxis()->SetLabelSize(0.05);
  hSame->GetYaxis()->SetTitleOffset(1.5);
  hSame->GetYaxis()->SetRangeUser(ymin, ymax);
  hSame->SetStats(0); hSame->SetTitle("");
  hSame->Draw("E");

  hMixed->SetLineColor(col_orange); hMixed->SetLineWidth(2); hMixed->SetLineStyle(2);
  hMixed->Draw("E same");

  hSignal->SetLineColor(col_green); hSignal->SetLineWidth(2); hSignal->SetLineStyle(1);
  hSignal->Draw("E same");

  TLegend *leg = new TLegend(0.42, 0.55, 0.92, 0.88);
  leg->SetBorderSize(0); leg->SetTextSize(0.050);
  leg->AddEntry(hSame,   "same-event FastJet (CS)", "l");
  leg->AddEntry(hMixed,  "mixed-event FastJet (CS, fake jets)", "l");
  leg->AddEntry(hSignal, "signal = same #minus mixed", "l");
  leg->Draw();

  TLatex *lat = new TLatex(); lat->SetNDC(); lat->SetTextSize(0.050);
  lat->DrawLatex(0.22, 0.85, Form("HYDJET 5.02 TeV  |#eta^{jet}| < 1.6  cent. %s", centLbl));

  // --- bottom pad: signal fraction ---
  TPad *p2 = (TPad*) c->cd(2);
  p2->SetLeftMargin(lm); p2->SetRightMargin(rm);
  p2->SetTopMargin(0.02); p2->SetBottomMargin(0.18);

  hFrac->SetLineColor(col_green); hFrac->SetLineWidth(2); hFrac->SetLineStyle(1);
  hFrac->GetXaxis()->SetRangeUser(pTlo, pThi);
  hFrac->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
  hFrac->GetYaxis()->SetTitle("signal / same-event");
  hFrac->GetYaxis()->SetRangeUser(-0.5, 1.5);
  hFrac->GetXaxis()->SetTitleSize(0.07); hFrac->GetXaxis()->SetLabelSize(0.065);
  hFrac->GetYaxis()->SetTitleSize(0.065); hFrac->GetYaxis()->SetLabelSize(0.065);
  hFrac->GetYaxis()->SetTitleOffset(1.1);
  hFrac->GetYaxis()->SetNdivisions(505);
  hFrac->SetStats(0); hFrac->SetTitle("");
  hFrac->Draw("E");

  TLine *l0 = new TLine(pTlo, 0., pThi, 0.);
  l0->SetLineStyle(2); l0->SetLineColor(kGray+1); l0->SetLineWidth(2); l0->Draw();
  TLine *l1 = new TLine(pTlo, 1., pThi, 1.);
  l1->SetLineStyle(2); l1->SetLineColor(kGray+1); l1->SetLineWidth(2); l1->Draw();

  TString outFile = TString(outDir) + Form("fakeJetSubtraction_C%d.pdf", ci);
  c->SaveAs(outFile);
  printf("Saved %s\n", outFile.Data());

  delete c;
  delete hSame;
  delete hMixed;
  delete hSignal;
  delete hFrac;
}

void plotJetPt_fakeJetSubtraction(){

  TFile *fSame  = TFile::Open(inFileSame);
  if(!fSame || fSame->IsZombie()){ printf("ERROR: cannot open same-event file\n"); return; }

  TFile *fMixed = TFile::Open(inFileMixed);
  if(!fMixed || fMixed->IsZombie()){ printf("ERROR: cannot open mixed-event file\n"); return; }

  gSystem->Exec(Form("mkdir -p %s", outDir));

  double lm = 0.20, rm = 0.05;
  for(int ci = 0; ci < NCentBins; ci++){
    plotOneBin(fSame, fMixed, ci, centLabel[ci], lm, rm);
  }

  fSame->Close();
  fMixed->Close();
}
