// Fine-centrality-bin twin of plotFastJet_bkgSubtraction_mixedVsSame.C:
// mixed/same for the FastJet fake-jet spectrum, one curve per ultra-fine
// centrality slice (16 x 5%) instead of summed into 4 coarse classes.
//
// Single panel, not the coarse version's spectra+ratio pair: 16 slices would
// need 32 curves (mixed and same, each slice) in a two-panel layout, which is
// unreadable. Since the coarse version already established what the absolute
// spectra look like, this one shows only the ratio -- the quantity that
// answers the motivating question, "does the mixed/same gap itself evolve
// smoothly with centrality, or does something happen specifically around
// 10-30%." Colour is a perceptually uniform (Viridis) gradient across the 16
// slices rather than 16 discrete hues, since centrality here is a continuous
// ordering, not a categorical label -- same convention as a 2D map's z-axis.
//
// See plotFastJet_bkgSubtraction_mixedVsSame.C for the full reasoning on what
// this ratio means and its same-event-contamination caveat; identical here,
// just at finer centrality resolution.
//
// Usage: root -l -b -q 'plotFastJet_bkgSubtraction_mixedVsSame_fineCent.C'
// Run from: src/plots/jetPt/pseudoJets/

const char *mixedFile =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/PbPb_MinBias_Part1_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_mixedEventPFClustering_pseudoJetCandPtMin-0.0_2026-8-11_ultraFineCentBins.root";
const char *sameFile =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/PbPb_MinBias_Part1_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_sameEventPFClustering_pseudoJetCandPtMin-0.0_2026-8-6_ultraFineCentBins.root";
const char *histBase = "h_fastJetPt_PF_bkgSub_RC";

const char *outDir  = "../../../../figures/jetPt/pseudoJets/";
const char *outName = "fastJet_bkgSub_mixedVsSame_fineCent.pdf";

const int    NSlice   = 16;
const double centStep = 5.0;

const double plotPtMin = 0., plotPtMax = 120.;
const double ratioMin = 0., ratioMax = 1.1;

// Coarser than the coarse-class version: an individual 5%-wide slice has
// roughly a quarter of the statistics of a summed class, and peripheral
// slices thin out fast above ~80 GeV (single-digit raw counts by 100 GeV in
// the most peripheral slice, checked before writing this).
const int    NEdge = 13;
double       ptEdge[NEdge] = {
  0,10,20,30,40,50,60,70,80,90,100,110,120
};

void plotFastJet_bkgSubtraction_mixedVsSame_fineCent()
{
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kViridis);
  gSystem->mkdir(outDir, kTRUE);

  TFile *fM = TFile::Open(mixedFile);
  TFile *fS = TFile::Open(sameFile);
  if(!fM || fM->IsZombie()){ printf("ERROR: cannot open %s\n", mixedFile); return; }
  if(!fS || fS->IsZombie()){ printf("ERROR: cannot open %s\n", sameFile);   return; }

  // Sample NSlice evenly-spaced colours from the current (Viridis) palette --
  // dark/central to light/peripheral, matching how a 2D map's z-axis reads.
  TColor::InitializeColors();
  TArrayI pal = TColor::GetPalette();
  int sliceCol[NSlice];
  for(int si = 0; si < NSlice; si++){
    int idx = si * (pal.GetSize() - 1) / (NSlice - 1);
    sliceCol[si] = pal[idx];
  }

  TH1D *hRatio[NSlice];
  bool useSlice[NSlice];

  printf("%-9s %10s %10s\n", "cent [%]", "N_mixed", "N_same");
  for(int si = 1; si <= NSlice; si++){
    int i = si - 1;

    TH1D *rawMixed = nullptr, *rawSame = nullptr;
    fM->GetObject(Form("%s_C%d", histBase, si), rawMixed);
    fS->GetObject(Form("%s_C%d", histBase, si), rawSame);
    TH1D *hvzM = nullptr, *hvzS = nullptr;
    fM->GetObject(Form("h_vz_C%d", si), hvzM);
    fS->GetObject(Form("h_vz_C%d", si), hvzS);
    if(!rawMixed || !rawSame || !hvzM || !hvzS){
      printf("WARNING: slice C%d missing histograms -- skipped\n", si);
      useSlice[i] = false;
      continue;
    }
    useSlice[i] = true;

    double N_mixed = hvzM->Integral(), N_same = hvzS->Integral();
    printf("%3.0f-%-5.0f %10.0f %10.0f\n", (si-1)*centStep, si*centStep, N_mixed, N_same);

    TH1D *hMixed = (TH1D*) rawMixed->Rebin(NEdge-1, Form("hMixed_%d", si), ptEdge);
    TH1D *hSame  = (TH1D*) rawSame ->Rebin(NEdge-1, Form("hSame_%d",  si), ptEdge);
    hMixed->SetDirectory(nullptr); hSame->SetDirectory(nullptr);
    hMixed->Scale(1./N_mixed);
    hSame ->Scale(1./N_same);

    hRatio[i] = (TH1D*) hMixed->Clone(Form("hRatio_%d", si));
    hRatio[i]->Divide(hSame);
  }

  // Warn rather than silently clip if a point falls outside the fixed range.
  for(int i = 0; i < NSlice; i++){
    if(!useSlice[i]) continue;
    for(int b = 1; b <= hRatio[i]->GetNbinsX(); b++){
      double lo = hRatio[i]->GetXaxis()->GetBinLowEdge(b);
      if(lo < plotPtMin || lo >= plotPtMax) continue;
      double v = hRatio[i]->GetBinContent(b);
      if(v != 0. && (v > ratioMax || v < ratioMin))
        printf("NOTE: slice C%d mixed/same at %.0f-%.0f is %.2f, outside range [%.1f,%.1f]\n",
               i+1, lo, hRatio[i]->GetXaxis()->GetBinUpEdge(b), v, ratioMin, ratioMax);
    }
  }

  const double lm = 0.13, rm = 0.14, tm = 0.09, bm = 0.13;
  TCanvas *c = new TCanvas("cFastJetMixedVsSameFine", "", 780, 620);
  c->SetLeftMargin(lm); c->SetRightMargin(rm);
  c->SetTopMargin(tm);  c->SetBottomMargin(bm);

  bool first = true;
  for(int i = 0; i < NSlice; i++){
    if(!useSlice[i]) continue;
    hRatio[i]->SetLineColor(sliceCol[i]); hRatio[i]->SetMarkerColor(sliceCol[i]);
    hRatio[i]->SetMarkerStyle(20); hRatio[i]->SetMarkerSize(0.7);
    hRatio[i]->SetLineWidth(2);
    hRatio[i]->SetTitle("");
    hRatio[i]->GetXaxis()->SetTitle("raw p_{T}^{fakeJet} [GeV]");
    hRatio[i]->GetYaxis()->SetTitle("mixed / same");
    hRatio[i]->GetXaxis()->SetTitleSize(0.045);
    hRatio[i]->GetXaxis()->SetLabelSize(0.040);
    hRatio[i]->GetYaxis()->SetTitleSize(0.045);
    hRatio[i]->GetYaxis()->SetLabelSize(0.040);
    hRatio[i]->GetYaxis()->SetTitleOffset(1.35);
    hRatio[i]->GetXaxis()->SetRangeUser(plotPtMin, plotPtMax);
    hRatio[i]->SetMinimum(ratioMin); hRatio[i]->SetMaximum(ratioMax);
    hRatio[i]->Draw(first ? "ep" : "ep same");
    first = false;
  }

  TLine *one = new TLine(plotPtMin, 1., plotPtMax, 1.);
  one->SetLineStyle(2); one->SetLineColor(kGray+1); one->Draw();

  TLatex lat; lat.SetNDC(); lat.SetTextSize(0.036);
  lat.DrawLatex(lm, 0.955, "FastJet fake-jet spectrum, raw p_{T}, RC-subtracted -- fine centrality");

  // --- hand-built colour-gradient legend (16 discrete colours, continuous read) ---
  double gx0 = 1. - rm + 0.025, gx1 = gx0 + 0.035;
  double gy0 = 1. - tm - 0.02,  gy1 = bm + 0.10;
  for(int i = 0; i < NSlice; i++){
    double y0 = gy0 - (gy0-gy1) * i     / NSlice;
    double y1 = gy0 - (gy0-gy1) * (i+1) / NSlice;
    TPave *box = new TPave(gx0, y1, gx1, y0, 0, "NDC");
    box->SetFillColor(sliceCol[i]);
    box->SetLineWidth(0);
    box->Draw();
  }
  TPave *gbord = new TPave(gx0, gy1, gx1, gy0, 0, "NDC");
  gbord->SetFillStyle(0); gbord->SetLineColor(kBlack); gbord->SetLineWidth(1);
  gbord->Draw();

  TLatex glab; glab.SetNDC(); glab.SetTextSize(0.030); glab.SetTextAlign(12);
  glab.DrawLatex(gx1 + 0.012, gy0, "0%");
  glab.DrawLatex(gx1 + 0.012, (gy0+gy1)/2., "40%");
  glab.DrawLatex(gx1 + 0.012, gy1, "80%");
  TLatex gtitle; gtitle.SetNDC(); gtitle.SetTextSize(0.030); gtitle.SetTextAngle(90); gtitle.SetTextAlign(22);
  gtitle.DrawLatex(gx1 + 0.075, (gy0+gy1)/2., "centrality");

  TString out = TString(outDir) + outName;
  c->SaveAs(out);
  printf("\nSaved %s\n", out.Data());
}
