// Compare RC-subtracted vs. unsubtracted FastJet anti-kT R=0.4 pT spectra.
// For each centrality bin, produces a two-pad canvas:
//   upper: h_fastJetPt and h_fastJetPt_bkgSub_RC overlaid (log-y)
//   lower: ratio  h_fastJetPt_bkgSub_RC / h_fastJetPt
//
// Usage (from src/plots/jetPt/pseudoJets/):
//   root -l -b -q plotFastJet_rcSubtraction.C

const char *inFile =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/PbPb_MinBias_Part1_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_mixedEventPseudoJets_pfCand_pseudoJetCandPtMin-0.0_2026-7-13_ultraFineCentBins.root";

const char *outDir =
  "../../../../figures/jetPt/pseudoJets/rcSubtraction/";

// ultraFine scheme: C0 = inclusive 0-80%, C1-C16 = 5% bins
const int NCentBins = 17;
const char *centLabel[NCentBins] = {
  "0-80% (incl.)",
  "0-5%",   "5-10%",  "10-15%", "15-20%",
  "20-25%", "25-30%", "30-35%", "35-40%",
  "40-45%", "45-50%", "50-55%", "55-60%",
  "60-65%", "65-70%", "70-75%", "75-80%"
};
const char *centSuffix[NCentBins] = {
  "C0",
  "C1",  "C2",  "C3",  "C4",
  "C5",  "C6",  "C7",  "C8",
  "C9",  "C10", "C11", "C12",
  "C13", "C14", "C15", "C16"
};

const double pTlo = 0.;
const double pThi = 300.;

const int colRaw  = kBlue  - 4;
const int colSub  = kRed   - 4;
const int colRC = kGray + 1;

static void styleHist(TH1D *h, int col, double lw = 2.)
{
  h->SetLineColor(col);
  h->SetLineWidth(lw);
  h->SetStats(0);
  h->SetTitle("");
  h->GetXaxis()->SetRangeUser(pTlo, pThi);
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetXaxis()->SetLabelSize(0.04);
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetLabelSize(0.04);
  h->GetYaxis()->SetTitleOffset(1.55);
}

void plotFastJet_rcSubtraction()
{
  gSystem->Exec(Form("mkdir -p %s", outDir));

  TFile *f = TFile::Open(inFile);
  if(!f || f->IsZombie()){
    std::cerr << "ERROR: cannot open " << inFile << "\n";
    return;
  }

  TLine *li = new TLine();
  li->SetLineStyle(7);
  li->SetLineColor(kGray + 1);

  for(int ci = 0; ci < NCentBins; ci++){

    TH1D *hRaw = nullptr;
    TH1D *hSub = nullptr;
    TH1D *hRC = nullptr;
    f->GetObject(Form("h_fastJetPt_%s",        centSuffix[ci]), hRaw);
    f->GetObject(Form("h_fastJetPt_bkgSub_RC_%s", centSuffix[ci]), hSub);
    f->GetObject(Form("h_pseudoJetPt_%s", centSuffix[ci]), hRC);

    if(!hRaw || !hSub){
      std::cerr << "WARNING: histograms not found for " << centSuffix[ci] << " — skipping\n";
      continue;
    }

    hRaw = (TH1D*) hRaw->Clone(Form("hRaw_%s", centSuffix[ci]));
    hSub = (TH1D*) hSub->Clone(Form("hSub_%s", centSuffix[ci]));
    hRC = (TH1D*) hRC->Clone(Form("hRC_%s",centSuffix[ci]));
    hRaw->SetDirectory(nullptr);
    hSub->SetDirectory(nullptr);
    hRC->SetDirectory(nullptr);

    styleHist(hRaw, colRaw);
    styleHist(hSub, colSub);
    styleHist(hRC, colRC);

    // ---- canvas with upper (spectra) and lower (ratio) pads ----------------
    TCanvas *c = new TCanvas(Form("c_%s", centSuffix[ci]), "", 700, 700);
    TPad *pUp = new TPad(Form("pUp_%s", centSuffix[ci]), "", 0, 0.35, 1, 1);
    TPad *pDn = new TPad(Form("pDn_%s", centSuffix[ci]), "", 0, 0,    1, 0.35);

    pUp->SetLeftMargin(0.15);  pUp->SetRightMargin(0.05);
    pUp->SetTopMargin(0.10);   pUp->SetBottomMargin(0.);
    pDn->SetLeftMargin(0.15);  pDn->SetRightMargin(0.05);
    pDn->SetTopMargin(0.);     pDn->SetBottomMargin(0.25);
    pUp->SetLogy();

    pUp->Draw(); pDn->Draw();

    // upper pad
    pUp->cd();

    double ymax = std::max(hRaw->GetMaximum(), hSub->GetMaximum());
    double ymin_pos = 1e30;
    for(int b = 1; b <= hRaw->GetNbinsX(); b++){
      double v = hRaw->GetBinContent(b);
      if(v > 0 && v < ymin_pos) ymin_pos = v;
    }
    if(ymin_pos > 1e29) ymin_pos = 1.;

    hRaw->GetXaxis()->SetRangeUser(pTlo, pThi);
    hRaw->GetYaxis()->SetRangeUser(ymin_pos * 0.05, ymax * 50.);
    hRaw->GetYaxis()->SetTitle("Entries");
    hRaw->GetXaxis()->SetTitle("");
    hRaw->GetXaxis()->SetLabelSize(0.);
    hRaw->Draw("hist");
    hSub->Draw("hist same");
    hRC->Draw("hist same");

    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.055);
    lat.DrawLatex(0.68, 0.82, Form("PbPb, %s", centLabel[ci]));
    lat.SetTextSize(0.045);
    lat.SetTextColor(colRaw);
    lat.DrawLatex(0.18, 0.76, "FastJet (raw)");
    lat.SetTextColor(colSub);
    lat.DrawLatex(0.18, 0.68, "FastJet (RC subtracted)");

    // lower pad: ratio
    pDn->cd();

    TH1D *hRatio = (TH1D*) hSub->Clone(Form("hRatio_%s", centSuffix[ci]));
    hRatio->Divide(hRaw);
    hRatio->SetDirectory(nullptr);

    hRatio->GetXaxis()->SetRangeUser(pTlo, pThi);
    hRatio->GetXaxis()->SetTitle("FastJet #it{p}_{T} [GeV]");
    hRatio->GetYaxis()->SetTitle("RC sub / raw");
    hRatio->GetYaxis()->SetRangeUser(0., 2.);
    hRatio->GetXaxis()->SetTitleSize(0.10);
    hRatio->GetXaxis()->SetLabelSize(0.09);
    hRatio->GetYaxis()->SetTitleSize(0.10);
    hRatio->GetYaxis()->SetLabelSize(0.09);
    hRatio->GetYaxis()->SetTitleOffset(0.65);
    hRatio->GetYaxis()->SetNdivisions(505);
    hRatio->Draw("hist");
    li->DrawLine(pTlo, 1., pThi, 1.);

    c->SaveAs(Form("%sfastJetPt_rcSub_%s.pdf", outDir, centSuffix[ci]));
    delete c;

    // ---- canvas 2: RC-subtracted FastJet vs CS reco jets (shape comparison) --
    TH1D *hCS = nullptr;
    f->GetObject(Form("h_inclRecoJetPt_%s", centSuffix[ci]), hCS);
    if(!hCS){
      std::cerr << "WARNING: h_inclRecoJetPt_" << centSuffix[ci] << " not found — skipping CS comparison\n";
      continue;
    }
    hCS = (TH1D*) hCS->Clone(Form("hCS_%s", centSuffix[ci]));
    hCS->SetDirectory(nullptr);

    // normalize to unity for shape comparison
    TH1D *hSubNorm = (TH1D*) hSub->Clone(Form("hSubNorm_%s", centSuffix[ci]));
    TH1D *hCSNorm  = (TH1D*) hCS ->Clone(Form("hCSNorm_%s",  centSuffix[ci]));
    hSubNorm->SetDirectory(nullptr);
    hCSNorm ->SetDirectory(nullptr);
    if(hSubNorm->Integral() > 0) hSubNorm->Scale(1. / hSubNorm->Integral(hSubNorm->FindBin(20),hSubNorm->FindBin(500)));
    if(hCSNorm ->Integral() > 0) hCSNorm ->Scale(1. / hCSNorm ->Integral(hCSNorm->FindBin(20),hCSNorm->FindBin(500)));

    const int colCS = kGreen + 2;
    styleHist(hSubNorm, colSub);
    styleHist(hCSNorm,  colCS);

    TCanvas *c2 = new TCanvas(Form("c2_%s", centSuffix[ci]), "", 700, 700);
    TPad *pUp2 = new TPad(Form("pUp2_%s", centSuffix[ci]), "", 0, 0.35, 1, 1);
    TPad *pDn2 = new TPad(Form("pDn2_%s", centSuffix[ci]), "", 0, 0,    1, 0.35);
    pUp2->SetLeftMargin(0.15);  pUp2->SetRightMargin(0.05);
    pUp2->SetTopMargin(0.10);   pUp2->SetBottomMargin(0.);
    pDn2->SetLeftMargin(0.15);  pDn2->SetRightMargin(0.05);
    pDn2->SetTopMargin(0.);     pDn2->SetBottomMargin(0.25);
    pUp2->SetLogy();
    pUp2->Draw(); pDn2->Draw();

    pUp2->cd();
    double ymax2 = std::max(hSubNorm->GetMaximum(), hCSNorm->GetMaximum());
    double ymin2 = 1e30;
    for(int b = 1; b <= hSubNorm->GetNbinsX(); b++){
      double v = hSubNorm->GetBinContent(b);
      if(v > 0 && v < ymin2) ymin2 = v;
    }
    if(ymin2 > 1e29) ymin2 = 1e-8;

    hSubNorm->GetXaxis()->SetRangeUser(pTlo, pThi);
    hSubNorm->GetYaxis()->SetRangeUser(ymin2 * 0.05, ymax2 * 50.);
    hSubNorm->GetYaxis()->SetTitle("Entries (norm. for #it{p}_{T} > 20 GeV)");
    hSubNorm->GetXaxis()->SetLabelSize(0.);
    hSubNorm->Draw("hist");
    hCSNorm ->Draw("hist same");

    TLatex lat2;
    lat2.SetNDC();
    lat2.SetTextSize(0.055);
    lat2.DrawLatex(0.68, 0.82, Form("PbPb, %s", centLabel[ci]));
    lat2.SetTextSize(0.045);
    lat2.SetTextColor(colSub);
    lat2.DrawLatex(0.18, 0.76, "FastJet (RC subtracted)");
    lat2.SetTextColor(colCS);
    lat2.DrawLatex(0.18, 0.68, "CS reco jets (akCs4PF)");

    pDn2->cd();
    TH1D *hRatio2 = (TH1D*) hSubNorm->Clone(Form("hRatio2_%s", centSuffix[ci]));
    hRatio2->Divide(hCSNorm);
    hRatio2->SetDirectory(nullptr);
    hRatio2->GetXaxis()->SetRangeUser(pTlo, pThi);
    hRatio2->GetXaxis()->SetTitle("Jet #it{p}_{T} [GeV]");
    hRatio2->GetYaxis()->SetTitle("RC sub / CS");
    hRatio2->GetYaxis()->SetRangeUser(0., 2.);
    hRatio2->GetXaxis()->SetTitleSize(0.10);
    hRatio2->GetXaxis()->SetLabelSize(0.09);
    hRatio2->GetYaxis()->SetTitleSize(0.10);
    hRatio2->GetYaxis()->SetLabelSize(0.09);
    hRatio2->GetYaxis()->SetTitleOffset(0.65);
    hRatio2->GetYaxis()->SetNdivisions(505);
    hRatio2->Draw("hist");
    li->DrawLine(pTlo, 1., pThi, 1.);

    c2->SaveAs(Form("%sfastJetPt_rcSub_vs_csJet_%s.pdf", outDir, centSuffix[ci]));
    delete c2;
  }

  f->Close();
  std::cout << "Figures saved to " << outDir << "\n";
}
