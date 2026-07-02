// Diagnostic plots for mixed-event PF candidate scans.
// Three figures, each showing all centrality bins overlaid:
//   1. FastJet anti-kT R=0.4 pT per event
//   2. RandomCone pT per event (h_pseudoJetPt, filled N_mixedEventsInPool times/event)
//   3. (FastJet - RandomCone) pT difference — mock UE subtraction
//
// NOTE: The current 2026-7-2 pfCand scan has h_fastJetPt empty (likely compiled
//       without -DDO_FASTJET) and h_pseudoJetPt all at zero (pool was empty).
//       The macro is written correctly for when the scan is re-run.
//
// Usage: root -l -q 'plotMixedEvent_pfCand.C'
// Run from: src/plots/jetPt/

const char *inFile =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/"
  "PbPb_MinBias_Part1_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_"
  "mixedEventPseudoJets_pfCand_pseudoJetCandPtMin-0.0_subleadingPFCandPtMin-15_2026-7-2.root";

const char *outDir = "/home/clayton/Analysis/code/bJetRaaAnalysis/figures/";

const int NCentBins     = 5;
const int N_pool        = 100;  // N_mixedEventsInPool
const double pTlo       = 0.;
const double pThi       = 150.;

const char *centLabel[NCentBins] = {
  "0-80% (incl.)", "0-10%", "10-30%", "30-50%", "50-80%"
};

static void styleHist(TH1D *h, int col, const char *xtitle, const char *ytitle,
                      double ylo, double yhi)
{
  h->SetLineColor(col);
  h->SetLineWidth(2);
  h->SetStats(0);
  h->SetTitle("");
  h->GetXaxis()->SetRangeUser(pTlo, pThi);
  h->GetXaxis()->SetTitle(xtitle);
  h->GetYaxis()->SetTitle(ytitle);
  h->GetXaxis()->SetTitleSize(0.05);  h->GetXaxis()->SetLabelSize(0.04);
  h->GetYaxis()->SetTitleSize(0.05);  h->GetYaxis()->SetLabelSize(0.04);
  h->GetYaxis()->SetTitleOffset(1.55);
  if(ylo < yhi) h->GetYaxis()->SetRangeUser(ylo, yhi);
}

static void drawCanvas(TH1D **hArr, int nBins, const char *ytitle,
                       bool logy, bool drawZeroLine, const char *label,
                       const char *outName)
{
  TCanvas *c = new TCanvas(outName, "", 800, 640);
  c->SetLeftMargin(0.15);
  c->SetBottomMargin(0.13);
  c->SetTopMargin(0.07);
  c->SetRightMargin(0.05);
  c->SetLogy(logy ? 1 : 0);

  double ymax = 0., ymin_pos = 1e30;
  for(int ci = 0; ci < nBins; ci++){
    for(int b = 1; b <= hArr[ci]->GetNbinsX(); b++){
      double v = hArr[ci]->GetBinContent(b);
      if(v > ymax) ymax = v;
      if(v > 0. && v < ymin_pos) ymin_pos = v;
    }
  }
  if(ymax == 0.) ymax = 1.;
  if(ymin_pos > 1e29) ymin_pos = 1e-6;

  double yrlo = logy ? ymin_pos * 0.05 : 0.;
  double yrhi = logy ? ymax  * 50.    : ymax * 1.45;

  // Okabe-Ito colorblind-safe palette
  int cols[5] = {
    kBlack,
    TColor::GetColor("#E69F00"),
    TColor::GetColor("#56B4E9"),
    TColor::GetColor("#009E73"),
    TColor::GetColor("#0072B2")
  };

  for(int ci = 0; ci < nBins; ci++){
    hArr[ci]->SetLineColor(cols[ci]);
    hArr[ci]->SetLineWidth(2);
    hArr[ci]->SetStats(0);
    hArr[ci]->SetTitle("");
    hArr[ci]->GetXaxis()->SetRangeUser(pTlo, pThi);
    hArr[ci]->GetXaxis()->SetTitle("p_{T} [GeV]");
    hArr[ci]->GetYaxis()->SetTitle(ytitle);
    hArr[ci]->GetXaxis()->SetTitleSize(0.05);  hArr[ci]->GetXaxis()->SetLabelSize(0.04);
    hArr[ci]->GetYaxis()->SetTitleSize(0.05);  hArr[ci]->GetYaxis()->SetLabelSize(0.04);
    hArr[ci]->GetYaxis()->SetTitleOffset(1.55);
    hArr[ci]->GetYaxis()->SetRangeUser(yrlo, yrhi);

    hArr[ci]->Draw(ci == 0 ? "hist" : "hist same");
  }

  if(drawZeroLine){
    TLine *l0 = new TLine(pTlo, 0., pThi, 0.);
    l0->SetLineStyle(2);
    l0->SetLineColor(kGray + 1);
    l0->SetLineWidth(2);
    l0->Draw();
  }

  TLegend *leg = new TLegend(0.55, 0.55, 0.93, 0.90);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.040);
  for(int ci = 0; ci < nBins; ci++)
    leg->AddEntry(hArr[ci], centLabel[ci], "l");
  leg->Draw();

  TLatex *lat = new TLatex();
  lat->SetNDC();
  lat->SetTextSize(0.040);
  lat->DrawLatex(0.16, 0.945, Form("PbPb 5.02 TeV  Mixed Event (PF cand)  %s", label));

  c->SaveAs(Form("%s%s", outDir, outName));
  printf("Saved %s%s\n", outDir, outName);

  delete c;
}

void plotMixedEvent_pfCand()
{
  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  TFile *f = TFile::Open(inFile);
  if(!f || f->IsZombie()){ printf("ERROR: cannot open %s\n", inFile); return; }

  gSystem->Exec(Form("mkdir -p %s", outDir));

  TH1D *hFJ[NCentBins], *hRC[NCentBins], *hDiff[NCentBins];

  for(int ci = 0; ci < NCentBins; ci++){
    // --- raw histograms from file ---
    TH1D *hFJ_raw = nullptr, *hRC_raw = nullptr;
    TH1D *hvz     = nullptr;
    f->GetObject(Form("h_fastJetPt_C%d",   ci), hFJ_raw);
    f->GetObject(Form("h_pseudoJetPt_C%d", ci), hRC_raw);
    f->GetObject(Form("h_vz_C%d",          ci), hvz);

    if(!hFJ_raw || !hRC_raw){
      printf("ERROR: missing histograms for C%d\n", ci);
      return;
    }

    double N_evts = hvz ? hvz->Integral() : 1.;
    if(N_evts <= 0.) N_evts = 1.;

    // --- clone and normalize per event ---
    hFJ[ci] = (TH1D*) hFJ_raw->Clone(Form("hFJ_C%d", ci));
    hRC[ci] = (TH1D*) hRC_raw->Clone(Form("hRC_C%d", ci));
    hFJ[ci]->SetDirectory(0);
    hRC[ci]->SetDirectory(0);

    hFJ[ci]->Scale(1. / N_evts);                      // jets per event per bin
    hRC[ci]->Scale(1. / (N_evts * (double)N_pool));   // cones per event per bin

    // --- difference: FastJet - RandomCone ---
    hDiff[ci] = (TH1D*) hFJ[ci]->Clone(Form("hDiff_C%d", ci));
    hDiff[ci]->SetDirectory(0);
    hDiff[ci]->Add(hRC[ci], -1.);
  }

  f->Close();

  // --- Canvas 1: FastJet pT ---
  drawCanvas(hFJ, NCentBins,
             "Anti-k_{T} R=0.4 jets / event / GeV",
             true, false,
             "anti-k_{T} R=0.4",
             "mixedEvent_pfCand_fastJetPt.pdf");

  // --- Canvas 2: RandomCone pT ---
  drawCanvas(hRC, NCentBins,
             "Random cones / event / GeV",
             true, false,
             "random cone R=0.4",
             "mixedEvent_pfCand_randomConePt.pdf");

  // --- Canvas 3: FastJet - RandomCone ---
  drawCanvas(hDiff, NCentBins,
             "(Anti-k_{T} #minus Random cone) / event / GeV",
             false, true,
             "anti-k_{T} #minus random cone",
             "mixedEvent_pfCand_difference.pdf");

  // clean up
  for(int ci = 0; ci < NCentBins; ci++){
    delete hFJ[ci]; delete hRC[ci]; delete hDiff[ci];
  }
}
