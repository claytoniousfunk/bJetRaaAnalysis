// All-methods twin of plotFastJet_bkgSubtraction_yieldVsCent_fineCent.C:
// same absolute (not ratio) fake-jet-yield-vs-fine-centrality plot, looped
// over every background-subtraction variant available on the raw-pT axis in
// the scan output (RC, RC_geoCorr, RC_geoCorr_etaReflect, dPT, and the two
// dPT sub-variants). One PDF per method, same style throughout so they're
// directly comparable side by side.
//
// dPT is known to over-subtract into negative territory in central events
// (see calculateRAA.C dPT runs, h_C1 integral around -3e9 in no-unfold mode)
// -- log-y silently breaks on a negative/zero point, so each method's yield
// array is checked first; if anything is <= 0, that method's plot falls back
// to a linear y-axis with a y=0 reference line instead, and a NOTE is
// printed rather than crashing or silently mis-rendering.
//
// Usage: root -l -b -q 'plotFastJet_bkgSubtraction_yieldVsCent_fineCent_allMethods.C'
// Run from: src/plots/jetPt/pseudoJets/

const char *mixedFile =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/PbPb_MinBias_Part1_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_mixedEventPFClustering_pseudoJetCandPtMin-0.0_2026-8-11_ultraFineCentBins.root";
const char *sameFile =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/PbPb_MinBias_Part1_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_sameEventPFClustering_pseudoJetCandPtMin-0.0_2026-8-6_ultraFineCentBins.root";

const char *outDir  = "../../../../figures/jetPt/pseudoJets/";

const int    NSlice   = 16;
const double centStep = 5.0;
const double ptThreshold = 40.0;   // raw GeV, matches signalJetPtCut

struct Method { const char *suffix; const char *label; const char *outTag; };
const int NMethod = 5;
Method methods[NMethod] = {
  {"RC",                   "RC",                          ""},
  {"RC_geoCorr",           "RC_{geoCorr}",                "_RCgeoCorr"},
  {"dPT",                  "dPT",                         "_dPT"},
  {"dPT_PFCsPTAbove60",    "dPT (PFCs p_{T}>60)",         "_dPT_PFCsAbove60"},
  {"dPT_dPTAbove0",        "dPT (dPT>0)",                 "_dPT_dPTAbove0"},
};

void plotOneMethod(TFile *fM, TFile *fS, const Method &m)
{
  TString histBase = TString("h_fastJetPt_PF_bkgSub_") + m.suffix;

  double centCenter[NSlice], yieldMixed[NSlice], yieldSame[NSlice];
  bool useSlice[NSlice];

  printf("\n--- %s (raw p_{T} > %.0f GeV) ---\n", m.label, ptThreshold);
  printf("%-9s %12s %12s\n", "cent [%]", "<N>_mixed", "<N>_same");
  for(int si = 1; si <= NSlice; si++){
    int i = si - 1;
    centCenter[i] = (si - 0.5) * centStep;

    TH1D *rawMixed = nullptr, *rawSame = nullptr;
    fM->GetObject(Form("%s_C%d", histBase.Data(), si), rawMixed);
    fS->GetObject(Form("%s_C%d", histBase.Data(), si), rawSame);
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
  bool anyNonPositive = false;
  int nUsed = 0;
  for(int i = 0; i < NSlice; i++){
    if(!useSlice[i]) continue;
    gMixed->SetPoint(nUsed, centCenter[i], yieldMixed[i]);
    gSame ->SetPoint(nUsed, centCenter[i], yieldSame[i]);
    for(double v : {yieldMixed[i], yieldSame[i]}){
      if(v <= 0.) anyNonPositive = true;
      if(v < yMin) yMin = v;
      if(v > yMax) yMax = v;
    }
    nUsed++;
  }

  bool useLogY = !anyNonPositive;
  if(!useLogY)
    printf("NOTE: %s has non-positive yield somewhere -- using linear y-axis\n", m.label);

  double frameYmin, frameYmax;
  if(useLogY){ frameYmin = yMin/2.; frameYmax = yMax*3.; }
  else{
    double pad = (yMax - yMin) * 0.15;
    if(pad <= 0.) pad = std::abs(yMax) > 0. ? std::abs(yMax)*0.15 : 1.;
    frameYmin = yMin - pad; frameYmax = yMax + pad;
  }

  const double lm = 0.15, rm = 0.05, tm = 0.09, bm = 0.13;
  TCanvas *c = new TCanvas(Form("cFastJetYieldVsCent_%s", m.suffix), "", 780, 620);
  c->SetLeftMargin(lm); c->SetRightMargin(rm);
  c->SetTopMargin(tm);  c->SetBottomMargin(bm);
  if(useLogY) c->SetLogy();

  TH1F *frame = c->DrawFrame(0., frameYmin, 80., frameYmax);
  frame->GetXaxis()->SetTitle("centrality [%]");
  frame->GetYaxis()->SetTitle(Form("mean fake jets / event (raw p_{T} > %.0f GeV)", ptThreshold));
  frame->GetXaxis()->SetTitleSize(0.045); frame->GetXaxis()->SetLabelSize(0.040);
  frame->GetYaxis()->SetTitleSize(0.045); frame->GetYaxis()->SetLabelSize(0.040);
  frame->GetYaxis()->SetTitleOffset(1.55);

  for(double xb : {10., 30., 50.}){
    TLine *l = new TLine(xb, frameYmin, xb, frameYmax);
    l->SetLineStyle(2); l->SetLineColor(kGray+1);
    l->Draw();
  }
  if(!useLogY){
    TLine *zero = new TLine(0., 0., 80., 0.);
    zero->SetLineStyle(3); zero->SetLineColor(kGray+2);
    zero->Draw();
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

  TLatex lat; lat.SetNDC(); lat.SetTextSize(0.030);
  lat.DrawLatex(lm, 0.955, Form("FastJet fake-jet yield, raw p_{T}, %s-subtracted -- fine centrality", m.label));

  TString out = TString(outDir) + "fastJet_bkgSub_yieldVsCent_fineCent" + m.outTag + ".pdf";
  c->SaveAs(out);
  printf("Saved %s\n", out.Data());
}

void plotFastJet_bkgSubtraction_yieldVsCent_fineCent_allMethods()
{
  gStyle->SetOptStat(0);
  gSystem->mkdir(outDir, kTRUE);

  TFile *fM = TFile::Open(mixedFile);
  TFile *fS = TFile::Open(sameFile);
  if(!fM || fM->IsZombie()){ printf("ERROR: cannot open %s\n", mixedFile); return; }
  if(!fS || fS->IsZombie()){ printf("ERROR: cannot open %s\n", sameFile);   return; }

  for(int mi = 0; mi < NMethod; mi++)
    plotOneMethod(fM, fS, methods[mi]);
}
