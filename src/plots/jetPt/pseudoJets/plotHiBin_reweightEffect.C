// hiBin distribution, raw vs hiBin-reweighted (to match the HardProbes Jet80
// trigger), for completeness alongside
// plotFastJet_bkgSubtraction_hiBinReweightEffect_yieldVsPt.C -- same two
// files. Each curve normalized to unit integral (shape only), since the two
// files have different raw event counts.
//
// hiBin runs 0-200, lower = more central. Jet80 strongly favors central
// collisions (needs a lot of activity to produce an 80 GeV jet), so the
// reweighted distribution is expected to shift toward low hiBin relative to
// the raw (flat MinBias) one.
//
// Usage: root -l -b -q 'plotHiBin_reweightEffect.C'
// Run from: src/plots/jetPt/pseudoJets/

const char *rawFile =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/PbPb_MinBias_Part1_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_mixedEventPFClustering_pseudoJetCandPtMin-0.0_2026-8-11_ultraFineCentBins.root";
const char *reweightedFile =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/PbPb_MinBias_Part1_hiBinReweightToHardProbesJet80_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_mixedEventPFClustering_pseudoJetCandPtMin-0.0_2026-8-9_ultraFineCentBins.root";

const char *outDir  = "../../../../figures/jetPt/pseudoJets/";
const char *outName = "hiBin_reweightEffect.pdf";

const char *rawHex = "#0072B2", *rwHex = "#D55E00";
const int   rawMark = 20, rwMark = 21;

void plotHiBin_reweightEffect()
{
  gStyle->SetOptStat(0);
  gSystem->mkdir(outDir, kTRUE);

  TFile *fRaw = TFile::Open(rawFile);
  TFile *fRw  = TFile::Open(reweightedFile);
  if(!fRaw || fRaw->IsZombie()){ printf("ERROR: cannot open %s\n", rawFile); return; }
  if(!fRw  || fRw->IsZombie()){  printf("ERROR: cannot open %s\n", reweightedFile); return; }

  TH1D *hRaw0 = nullptr, *hRw0 = nullptr;
  fRaw->GetObject("h_hiBin", hRaw0);
  fRw->GetObject("h_hiBin", hRw0);
  if(!hRaw0 || !hRw0){ printf("ERROR: h_hiBin missing\n"); return; }

  TH1D *hRaw = (TH1D*) hRaw0->Clone("hRaw"); hRaw->SetDirectory(nullptr);
  TH1D *hRw  = (TH1D*) hRw0->Clone("hRw");   hRw->SetDirectory(nullptr);
  hRaw->Scale(1./hRaw->Integral());
  hRw->Scale(1./hRw->Integral());

  TH1D *hRatio = (TH1D*) hRw->Clone("hRatio");
  hRatio->Divide(hRaw);

  printf("%-10s %12s %12s\n", "hiBin", "raw (norm)", "reweighted");
  for(int b = 1; b <= hRaw->GetNbinsX(); b += 10)
    printf("%-10.0f %12.5f %12.5f\n", hRaw->GetXaxis()->GetBinLowEdge(b), hRaw->GetBinContent(b), hRw->GetBinContent(b));

  int colRaw = TColor::GetColor(rawHex), colRw = TColor::GetColor(rwHex);

  const double lm = 0.15, rm = 0.05, tm = 0.09, bm = 0.28, split = 0.35;
  TCanvas *c = new TCanvas("cHiBinReweightEffect", "", 800, 750);

  TPad *pUp = new TPad("pUp", "", 0, split, 1, 1);
  pUp->SetLeftMargin(lm); pUp->SetRightMargin(rm);
  pUp->SetTopMargin(tm);  pUp->SetBottomMargin(0.);
  pUp->SetLogy(); pUp->Draw(); pUp->cd();

  hRaw->SetLineColor(colRaw); hRaw->SetMarkerColor(colRaw);
  hRaw->SetMarkerStyle(rawMark); hRaw->SetMarkerSize(0.8); hRaw->SetLineWidth(2);
  hRaw->SetTitle("");
  hRaw->GetYaxis()->SetTitle("1/N_{evt} dN/dhiBin");
  hRaw->GetYaxis()->SetTitleSize(0.052); hRaw->GetYaxis()->SetLabelSize(0.045);
  hRaw->GetYaxis()->SetTitleOffset(1.40);
  hRaw->GetXaxis()->SetRangeUser(0., 200.);
  hRaw->SetMinimum(TMath::Max(hRaw->GetMinimum(0.), hRw->GetMinimum(0.)) * 0.3);
  hRaw->SetMaximum(TMath::Max(hRaw->GetMaximum(), hRw->GetMaximum()) * 3.);
  hRaw->Draw("ep");

  hRw->SetLineColor(colRw); hRw->SetMarkerColor(colRw);
  hRw->SetMarkerStyle(rwMark); hRw->SetMarkerSize(0.8); hRw->SetLineWidth(2);
  hRw->Draw("ep same");

  TLegend *leg = new TLegend(0.55, 0.68, 0.94, 0.86);
  leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.040);
  leg->AddEntry(hRaw, "raw hiBin",        "lp");
  leg->AddEntry(hRw,  "reweighted hiBin", "lp");
  leg->Draw();

  TLatex lat; lat.SetNDC(); lat.SetTextSize(0.036);
  lat.DrawLatex(lm, 0.955, "hiBin distribution -- reweighted to HardProbes Jet80, shape only");

  c->cd();
  TPad *pDn = new TPad("pDn", "", 0, 0, 1, split);
  pDn->SetLeftMargin(lm); pDn->SetRightMargin(rm);
  pDn->SetTopMargin(0.);  pDn->SetBottomMargin(bm);
  pDn->Draw(); pDn->cd();

  const double sc = (1.-split)/split;
  hRatio->SetLineColor(kBlack); hRatio->SetMarkerColor(kBlack);
  hRatio->SetMarkerStyle(20); hRatio->SetMarkerSize(0.7 * sc / 1.6); hRatio->SetLineWidth(2);
  hRatio->SetTitle("");
  hRatio->GetXaxis()->SetTitle("hiBin");
  hRatio->GetYaxis()->SetTitle("reweighted / raw");
  hRatio->GetXaxis()->SetTitleSize(0.052*sc); hRatio->GetXaxis()->SetLabelSize(0.045*sc);
  hRatio->GetYaxis()->SetTitleSize(0.052*sc); hRatio->GetYaxis()->SetLabelSize(0.045*sc);
  hRatio->GetYaxis()->SetTitleOffset(1.40/sc);
  hRatio->GetYaxis()->SetNdivisions(505);
  hRatio->GetXaxis()->SetRangeUser(0., 200.);
  hRatio->Draw("ep");

  TLine *one = new TLine(0., 1., 200., 1.);
  one->SetLineStyle(2); one->SetLineColor(kGray+1); one->Draw();

  TString out = TString(outDir) + outName;
  c->SaveAs(out);
  printf("\nSaved %s\n", out.Data());
}
