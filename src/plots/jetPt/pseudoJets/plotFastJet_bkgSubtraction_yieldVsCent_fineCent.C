// Absolute (not ratio) fake-jet yield vs centrality, fine bins (16 x 5%).
//
// Direct follow-up to plotFastJet_bkgSubtraction_mixedVsSame_fineCent.C: that
// plot showed the mixed/same RATIO evolves smoothly with centrality, no
// discontinuity around 10-30%. This checks the companion question -- does the
// ABSOLUTE fake-jet rate itself jump discontinuously right at the 0-10%/
// 10-30% boundary, or is "much much higher in 0-10 than 10-30" just the
// smooth continuation of a steep trend (expected from underlying-event
// multiplicity growth toward central collisions)?
//
// Mean number of RC-subtracted FastJet-found jets per event, integrated above
// ptThreshold, plotted vs centrality for both the same-event and mixed-event
// estimates. ptThreshold defaults to 40 GeV raw -- matches signalJetPtCut
// used elsewhere this session -- to avoid the near-zero-pT algorithmic-
// remnant spike (FastJet clustering here has no pT floor: pseudoJetCandPtMin
// = 0.0) swamping the integral with entries that aren't real fake-jet
// candidates in the sense relevant to contaminating the reco spectrum feeding
// unfolding. Vertical dashed lines mark the coarse-class boundaries (10%,
// 30%, 50%) under discussion.
//
// Usage: root -l -b -q 'plotFastJet_bkgSubtraction_yieldVsCent_fineCent.C'
// Run from: src/plots/jetPt/pseudoJets/

const char *mixedFile =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/PbPb_MinBias_Part1_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_mixedEventPFClustering_pseudoJetCandPtMin-0.0_2026-8-11_ultraFineCentBins.root";
const char *sameFile =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/PbPb_MinBias_Part1_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_sameEventPFClustering_pseudoJetCandPtMin-0.0_2026-8-6_ultraFineCentBins.root";
const char *histBase = "h_fastJetPt_PF_bkgSub_RC";

const char *outDir  = "../../../../figures/jetPt/pseudoJets/";
const char *outName = "fastJet_bkgSub_yieldVsCent_fineCent.pdf";

const int    NSlice   = 16;
const double centStep = 5.0;
const double ptThreshold = 40.0;   // raw GeV, matches signalJetPtCut

void plotFastJet_bkgSubtraction_yieldVsCent_fineCent()
{
  gStyle->SetOptStat(0);
  gSystem->mkdir(outDir, kTRUE);

  TFile *fM = TFile::Open(mixedFile);
  TFile *fS = TFile::Open(sameFile);
  if(!fM || fM->IsZombie()){ printf("ERROR: cannot open %s\n", mixedFile); return; }
  if(!fS || fS->IsZombie()){ printf("ERROR: cannot open %s\n", sameFile);   return; }

  double centCenter[NSlice], yieldMixed[NSlice], yieldSame[NSlice];
  bool useSlice[NSlice];

  printf("Integrating raw p_{T} > %.0f GeV\n", ptThreshold);
  printf("%-9s %12s %12s\n", "cent [%]", "<N>_mixed", "<N>_same");
  for(int si = 1; si <= NSlice; si++){
    int i = si - 1;
    centCenter[i] = (si - 0.5) * centStep;

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

    double N_mixedEvt = hvzM->Integral(), N_sameEvt = hvzS->Integral();
    int binLo_M = rawMixed->FindBin(ptThreshold + 1e-6);
    int binLo_S = rawSame ->FindBin(ptThreshold + 1e-6);
    yieldMixed[i] = rawMixed->Integral(binLo_M, rawMixed->GetNbinsX()) / N_mixedEvt;
    yieldSame[i]  = rawSame ->Integral(binLo_S, rawSame ->GetNbinsX()) / N_sameEvt;

    printf("%3.0f-%-5.0f %12.5f %12.5f\n",
           (si-1)*centStep, si*centStep, yieldMixed[i], yieldSame[i]);
  }

  TGraph *gMixed = new TGraph();
  TGraph *gSame  = new TGraph();
  double yMin = 1e300, yMax = -1e300;
  int nUsed = 0;
  for(int i = 0; i < NSlice; i++){
    if(!useSlice[i]) continue;
    gMixed->SetPoint(nUsed, centCenter[i], yieldMixed[i]);
    gSame ->SetPoint(nUsed, centCenter[i], yieldSame[i]);
    for(double v : {yieldMixed[i], yieldSame[i]}){
      if(v <= 0.) continue;
      if(v < yMin) yMin = v;
      if(v > yMax) yMax = v;
    }
    nUsed++;
  }

  const double lm = 0.15, rm = 0.05, tm = 0.09, bm = 0.13;
  TCanvas *c = new TCanvas("cFastJetYieldVsCent", "", 780, 620);
  c->SetLeftMargin(lm); c->SetRightMargin(rm);
  c->SetTopMargin(tm);  c->SetBottomMargin(bm);
  c->SetLogy();

  TH1F *frame = c->DrawFrame(0., yMin/2., 80., yMax*3.);
  frame->GetXaxis()->SetTitle("centrality [%]");
  frame->GetYaxis()->SetTitle(Form("mean fake jets / event (raw p_{T} > %.0f GeV)", ptThreshold));
  frame->GetXaxis()->SetTitleSize(0.045); frame->GetXaxis()->SetLabelSize(0.040);
  frame->GetYaxis()->SetTitleSize(0.045); frame->GetYaxis()->SetLabelSize(0.040);
  frame->GetYaxis()->SetTitleOffset(1.55);

  for(double xb : {10., 30., 50.}){
    TLine *l = new TLine(xb, yMin/2., xb, yMax*3.);
    l->SetLineStyle(2); l->SetLineColor(kGray+1);
    l->Draw();
  }

  gSame->SetMarkerStyle(20); gSame->SetMarkerColor(TColor::GetColor("#D55E00"));
  gSame->SetLineColor(TColor::GetColor("#D55E00")); gSame->SetMarkerSize(1.1); gSame->SetLineWidth(2);
  gSame->Draw("lp same");

  gMixed->SetMarkerStyle(24); gMixed->SetMarkerColor(TColor::GetColor("#0072B2"));
  gMixed->SetLineColor(TColor::GetColor("#0072B2")); gMixed->SetMarkerSize(1.1); gMixed->SetLineWidth(2);
  gMixed->Draw("lp same");

  TLegend *leg = new TLegend(0.62, 0.76, 0.93, 0.89);
  leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.032);
  leg->AddEntry(gSame,  "same-event", "lp");
  leg->AddEntry(gMixed, "mixed-event", "lp");
  leg->Draw();

  TLatex lat; lat.SetNDC(); lat.SetTextSize(0.036);
  lat.DrawLatex(lm, 0.955, "FastJet fake-jet yield, raw p_{T}, RC-subtracted -- fine centrality");

  TString out = TString(outDir) + outName;
  c->SaveAs(out);
  printf("\nSaved %s\n", out.Data());
}
