// Shape comparison between FastJet anti-kT R=0.4 jets (on CS PF candidates)
// and inclusive CMS reco jets, from HYDJET pfCandAnalyzer output.
// Normalizes both spectra over the visible pT range [pTlo, pThi].
// Usage: root -l -q 'plotJetPt_fastJetVsReco.C'

const char *inFile =
  "../../rootFiles/scanningOuput/HYDJET/"
  "HYDJET_pThat-unweighted_mu12_pTmu-15to999_tight_hiBinShift-0_"
  "jetTrkMaxFilter_WDecayFilter_sameEventPseudoJets_pfCandCS_"
  "pseudoJetCandPtMin-0.0_subleadingPFCandPtMin-15_2026-6-16_ultraFineCentBins.root";

const char *outDir = "../../figures/jetPt/";

const double pTlo = 20.;
const double pThi = 200.;

void plotJetPt_fastJetVsReco(){

  TFile *f = TFile::Open(inFile);
  if(!f || f->IsZombie()){ printf("ERROR: cannot open %s\n", inFile); return; }

  TH2D *h2; f->GetObject("h_inclRecoJetPt_flavor_C0", h2);
  TH1D *hReco = (TH1D*) h2->ProjectionX("h_recoJetPt");
  hReco->SetDirectory(0);

  TH1D *hFJ; f->GetObject("h_fastJetPt_C0", hFJ);
  hFJ = (TH1D*) hFJ->Clone("hFJ");
  hFJ->SetDirectory(0);
  f->Close();

  // normalize over visible range
  int blo = hReco->FindBin(pTlo + 0.01);
  int bhi = hReco->FindBin(pThi - 0.01);
  hReco->Scale(1./hReco->Integral(blo, bhi));
  hFJ->Scale(1./hFJ->Integral(blo, bhi));

  TH1D *hRatio = (TH1D*) hFJ->Clone("hRatio");
  hRatio->Divide(hReco);

  double lm = 0.18, rm = 0.05;

  TCanvas *c = new TCanvas("c","",900,800);
  c->Divide(1,2);

  TPad *p1 = (TPad*)c->cd(1);
  p1->SetLogy();
  p1->SetLeftMargin(lm); p1->SetRightMargin(rm);
  p1->SetBottomMargin(0.02); p1->SetTopMargin(0.08);

  hReco->SetLineColor(kBlue+1); hReco->SetLineWidth(2);
  hReco->GetXaxis()->SetRangeUser(pTlo, pThi);
  hReco->GetXaxis()->SetLabelSize(0);
  hReco->GetYaxis()->SetTitle("(1/N) dN/dp_{T} [GeV^{-1}]");
  hReco->GetYaxis()->SetTitleSize(0.06); hReco->GetYaxis()->SetLabelSize(0.055);
  hReco->GetYaxis()->SetTitleOffset(1.3);
  hReco->SetStats(0); hReco->SetTitle("");
  hReco->Draw("hist");

  hFJ->SetLineColor(kRed+1); hFJ->SetLineWidth(2); hFJ->SetLineStyle(2);
  hFJ->Draw("hist same");

  TLegend *leg = new TLegend(0.48, 0.65, 0.90, 0.88);
  leg->SetBorderSize(0); leg->SetTextSize(0.057);
  leg->AddEntry(hReco, "inclusive reco jets", "l");
  leg->AddEntry(hFJ,   "FastJet anti-k_{T} R = 0.4 (CS)", "l");
  leg->Draw();

  TLatex *lat = new TLatex(); lat->SetNDC(); lat->SetTextSize(0.055);
  lat->DrawLatex(0.22, 0.88, "HYDJET 5.02 TeV  |#eta^{jet}| < 1.6");

  TPad *p2 = (TPad*)c->cd(2);
  p2->SetLeftMargin(lm); p2->SetRightMargin(rm);
  p2->SetTopMargin(0.02); p2->SetBottomMargin(0.18);

  hRatio->SetLineColor(kBlack); hRatio->SetLineWidth(2);
  hRatio->GetXaxis()->SetRangeUser(pTlo, pThi);
  hRatio->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
  hRatio->GetYaxis()->SetTitle("FastJet (CS) / Reco");
  hRatio->GetYaxis()->SetRangeUser(0, 2.5);
  hRatio->GetXaxis()->SetTitleSize(0.07); hRatio->GetXaxis()->SetLabelSize(0.065);
  hRatio->GetYaxis()->SetTitleSize(0.07); hRatio->GetYaxis()->SetLabelSize(0.065);
  hRatio->GetYaxis()->SetTitleOffset(1.1);
  hRatio->GetYaxis()->SetNdivisions(505);
  hRatio->SetStats(0); hRatio->SetTitle("");
  hRatio->Draw("hist");

  TLine *l = new TLine(pTlo, 1., pThi, 1.);
  l->SetLineStyle(2); l->SetLineColor(kGray+1); l->SetLineWidth(2); l->Draw();

  gSystem->Exec(Form("mkdir -p %s", outDir));
  TString outFile = TString(outDir) + "fastJetVsReco_shape.pdf";
  c->SaveAs(outFile);
  printf("Saved %s\n", outFile.Data());
}
