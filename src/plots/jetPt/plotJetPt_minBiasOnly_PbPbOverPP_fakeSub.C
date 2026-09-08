// PbPb/pp jet pT shape ratio in FINE (5%) centrality slices, with the
// combinatorial fake-jet background SUBTRACTED from the PbPb side.
//
// Companion to plotJetPt_minBiasOnly_PbPbOverPP.C, which is the same
// measurement with no subtraction. Run both and compare: the un-subtracted
// version is the reference, this one is the correction, and the difference
// between them is entirely the fake-jet estimate.
//
// WHY A DIFFERENT INPUT FILE THAN THE UN-SUBTRACTED VERSION.
// The subtraction has to be done event-sample-for-event-sample, so the
// measured spectrum and the fake estimate must come from the SAME scan. The
// 2026-8-13 file used by the un-subtracted macro has no FastJet histograms at
// all, so this macro reads both from the 2026-8-14 mixed-event PF clustering
// scan instead. Checked before trusting it: same 5 GeV binning over 0-500 GeV
// on both histograms, same per-slice event counts, and the slice1/slice2 event
// ratio is 1.0009, i.e. the NEW [lo,hi) centrality convention (1.100 would
// mean the OLD one). Running this macro with subtractMethod = mNone
// reproduces the un-subtracted shape from the 8-13 file closely enough
// (7.19/7.67/7.53/9.08 vs 7.47/7.86/8.80/...) that the file swap is not doing
// the work; the residual is scan-to-scan statistics.
//
// PP IS NOT SUBTRACTED, DELIBERATELY. The fake jets being removed here are
// combinatorial pseudo-jets built out of underlying-event PF candidates, which
// is a PbPb-only pathology. pp has no comparable combinatorial background, so
// subtracting anything from it would be inventing a correction.
//
// *** THE METHOD CHOICE DOMINATES THE CENTRAL SLICES -- READ THIS ***
// At 50 GeV the three available fake estimates are wildly different fractions
// of the measured 0-5% yield: RC 51%, RC_geoCorr 85%, dPT 272%. The last two
// drive the most central slices to or below zero, which is unphysical and
// means they cannot be right at this pT. RC is the default here because it is
// the only method that leaves a positive remainder in central events, NOT
// because it has been shown to be correct. The printed sensitivity table at
// the bottom reports the probe point under all three so the size of that
// systematic is visible; do not quote a central-slice number from this plot
// without it.
//
// NEGATIVE BINS are zeroed rather than drawn, and counted per slice in the
// printout. A slice with more than a couple of them is telling you the
// subtraction has overshot, not that the yield is small.
//
// NORMALISATION. Unchanged from the un-subtracted macro: each curve is scaled
// so its own integral above ptNormMin is 1, so the ratio asks purely about
// SHAPE. The fake estimate above 150 GeV is negligible (the raw count above
// 150 in 0-5% moves 238 -> 225 under RC), so the normalisation anchor is
// essentially fake-free and the subtraction acts on the shape, not on the
// scale. The minimum-statistics criterion is evaluated on the UNSUBTRACTED
// count on purpose, so that the set of drawn slices does not change when the
// method changes -- it is a statement about the sample, not about the result.
//
// Usage: root -l -b -q 'plotJetPt_minBiasOnly_PbPbOverPP_fakeSub.C'
// Run from: src/plots/jetPt/

const char *pbpbFile =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/PbPb_MinBias_Part1_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_mixedEventPFClustering_pseudoJetCandPtMin-0.0_2026-8-14_ultraFineCentBins.root";
const char *ppFile =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/pp/pp_MinBias_mu12_pTmu-15to999_tight_deltaR-40_jetTrkMaxFilter_WDecayFilter_2026-7-6.root";

// JEC-corrected pT on both sides, and the JEC variant of the fake estimate to
// match -- mixing raw and corrected between the two would subtract a spectrum
// that is bodily shifted by the ~1.15 JEC factor.
const char *pbpbHistBase = "h_inclRecoJetPt";
const char *ppHistName   = "h_inclRecoJetPt";

enum { mNone = 0, mRC, mRCgeoCorr, mdPT, NMethod };
const char *methHist[NMethod] = {
  "",
  "h_fastJetPt_PF_JEC_bkgSub_RC",
  "h_fastJetPt_PF_JEC_bkgSub_RC_geoCorr",
  "h_fastJetPt_PF_JEC_bkgSub_dPT"
};
const char *methLabel[NMethod] = { "no subtraction", "RC", "RC geoCorr", "dPT" };

const int subtractMethod = mRC;   // see header before changing

const char *outDir  = "../../../figures/jetPt/";
const char *outName = "jetPt_minBiasOnly_PbPbOverPP_fineCent_fakeSub.pdf";

const int    NSlice   = 16;
const double centStep = 5.0;

const double ptNormMin   = 150.;
const int    minNormJets = 25;    // applied to the UNSUBTRACTED count, see header

const double plotPtMin = 30., plotPtMax = 300.;
const double ratioMin  = 0.,  ratioMax  = 11.;   // same axes as the un-subtracted
const double dblMin    = 0.3, dblMax    = 2.4;   // wider: subtraction moves 5-10/0-5 a lot

const double probePt = 52.5;

const int    NEdge = 24;
double       ptEdge[NEdge] = {
  20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,
  110,120,140,160,200,250,300
};

const double PALFRACLO = 0.22;
static void sliceColours(int *col)
{
  TColor::InitializeColors();
  gStyle->SetPalette(kViridis);
  TArrayI pal = TColor::GetPalette();
  int n = pal.GetSize();
  for(int i = 0; i < NSlice; i++){
    double frac = PALFRACLO + (1.0 - PALFRACLO) * i / (NSlice - 1);
    col[i] = pal[(int)(frac * (n - 1))];
  }
}

// Measured minus fake, on the native binning, with negatives zeroed. Returns
// the subtracted spectrum and reports how many bins had to be clipped.
static TH1D* subtracted(TFile *f, int si, int method, const char *name, int &nNeg)
{
  nNeg = 0;
  TH1D *meas = nullptr;
  f->GetObject(Form("%s_C%d", pbpbHistBase, si), meas);
  if(!meas) return nullptr;

  TH1D *w = (TH1D*) meas->Clone(name);
  w->SetDirectory(nullptr);
  if(method == mNone) return w;

  TH1D *fake = nullptr;
  f->GetObject(Form("%s_C%d", methHist[method], si), fake);
  if(!fake){ printf("ERROR: %s_C%d missing\n", methHist[method], si); delete w; return nullptr; }
  if(fake->GetNbinsX() != w->GetNbinsX() ||
     TMath::Abs(fake->GetXaxis()->GetXmax() - w->GetXaxis()->GetXmax()) > 1e-6){
    printf("ERROR: binning mismatch between measured and fake in slice %d\n", si);
    delete w; return nullptr;
  }

  w->Add(fake, -1.);
  for(int b = 1; b <= w->GetNbinsX(); b++){
    if(w->GetBinContent(b) < 0.){ nNeg++; w->SetBinContent(b, 0.); w->SetBinError(b, 0.); }
  }
  return w;
}

static TH1D* buildShape(TH1D *raw, const char *name)
{
  TH1D *h = (TH1D*) raw->Rebin(NEdge-1, name, ptEdge);
  h->SetDirectory(nullptr);
  double norm = h->Integral(h->FindBin(ptNormMin), h->GetNbinsX()+1);
  if(norm <= 0.){ delete h; return nullptr; }
  h->Scale(1./norm);
  h->Scale(1., "width");
  return h;
}

void plotJetPt_minBiasOnly_PbPbOverPP_fakeSub()
{
  gStyle->SetOptStat(0);
  gSystem->mkdir(outDir, kTRUE);

  TFile *fP = TFile::Open(pbpbFile);
  TFile *fp = TFile::Open(ppFile);
  if(!fP || fP->IsZombie()){ printf("ERROR: cannot open %s\n", pbpbFile); return; }
  if(!fp || fp->IsZombie()){ printf("ERROR: cannot open %s\n", ppFile);   return; }

  printf("Fake-jet subtraction: %s\n", methLabel[subtractMethod]);

  // Centrality convention probe, checked not assumed: the convention is not
  // recorded in the file, but slice1/slice2 event counts give it away.
  {
    TH1D *v1 = nullptr, *v2 = nullptr;
    fP->GetObject("h_vz_C1", v1); fP->GetObject("h_vz_C2", v2);
    if(v1 && v2 && v2->Integral() > 0){
      double r = v1->Integral()/v2->Integral();
      printf("Centrality convention probe: slice1/slice2 events = %.4f -> %s\n\n",
             r, r > 1.05 ? "OLD (lo,hi] -- SLICE 1 IS 10%% WIDER, RESULTS WILL BE BIASED"
                         : "NEW [lo,hi), as expected");
    }
  }

  // ---- pp reference, NOT subtracted (see header) ----
  TH1D *hpp0 = nullptr;
  fp->GetObject(ppHistName, hpp0);
  if(!hpp0){ printf("ERROR: %s not found in pp file\n", ppHistName); return; }
  TH1D *hpp = buildShape(hpp0, "hpp");
  if(!hpp){ printf("ERROR: pp has zero integral above %.0f GeV\n", ptNormMin); return; }

  int col[NSlice];
  sliceColours(col);

  TH1D *hRatio[NSlice], *hShape[NSlice], *hDouble[NSlice];
  bool  useSlice[NSlice];
  double probeVal[NSlice], peakVal[NSlice], peakPt[NSlice], fakeFrac[NSlice];
  int    prevOf[NSlice], negBins[NSlice];
  for(int i = 0; i < NSlice; i++){ hShape[i]=nullptr; hDouble[i]=nullptr; prevOf[i]=-1; }

  printf("%-9s %10s %9s %9s %9s %9s %6s\n",
         "cent", "raw>norm", "fake@52", "ratio@52", "peak", "peak pT", "neg");
  for(int si = 1; si <= NSlice; si++){
    int i = si - 1;
    useSlice[i] = false;
    probeVal[i] = peakVal[i] = peakPt[i] = fakeFrac[i] = -1.;
    negBins[i] = 0;

    // statistics criterion on the UNSUBTRACTED spectrum, so the drawn slice
    // set is identical whichever method is selected
    TH1D *measured = nullptr;
    fP->GetObject(Form("%s_C%d", pbpbHistBase, si), measured);
    if(!measured){ printf("WARNING: %s_C%d missing\n", pbpbHistBase, si); continue; }
    double aboveNorm = measured->Integral(measured->FindBin(ptNormMin + 1e-6),
                                          measured->GetNbinsX()+1);

    TH1D *w = subtracted(fP, si, subtractMethod, Form("wPbPb_C%d", si), negBins[i]);
    if(!w) continue;

    // fake fraction of the measured yield at the probe point, the number that
    // says how much of this slice's result is the subtraction's opinion
    if(subtractMethod != mNone){
      TH1D *fake = nullptr;
      fP->GetObject(Form("%s_C%d", methHist[subtractMethod], si), fake);
      int pb = measured->FindBin(probePt);
      if(fake && measured->GetBinContent(pb) > 0.)
        fakeFrac[i] = fake->GetBinContent(pb) / measured->GetBinContent(pb);
    }

    if(aboveNorm < minNormJets){
      printf("%3.0f-%-5.0f %10.0f   SKIPPED (fewer than %d raw jets above %.0f GeV;\n"
             "                        its normalisation would be noise)\n",
             (si-1)*centStep, si*centStep, aboveNorm, minNormJets, ptNormMin);
      delete w; continue;
    }

    TH1D *shaped = buildShape(w, Form("hPbPb_C%d", si));
    delete w;
    if(!shaped){ printf("WARNING: slice %d has zero integral above norm after subtraction\n", si); continue; }

    hShape[i] = shaped;
    hRatio[i] = (TH1D*) shaped->Clone(Form("hRatio_C%d", si));
    hRatio[i]->SetDirectory(nullptr);
    hRatio[i]->Divide(hpp);
    useSlice[i] = true;

    probeVal[i] = hRatio[i]->GetBinContent(hRatio[i]->FindBin(probePt));
    for(int b = 1; b <= hRatio[i]->GetNbinsX(); b++){
      double lo = hRatio[i]->GetXaxis()->GetBinLowEdge(b);
      if(lo < plotPtMin || lo >= 120.) continue;
      double v = hRatio[i]->GetBinContent(b);
      if(v > peakVal[i]){ peakVal[i] = v; peakPt[i] = hRatio[i]->GetXaxis()->GetBinCenter(b); }
    }
    printf("%3.0f-%-5.0f %10.0f %8.0f%% %9.2f %9.2f %9.0f %6d%s\n",
           (si-1)*centStep, si*centStep, aboveNorm, 100.*fakeFrac[i],
           probeVal[i], peakVal[i], peakPt[i], negBins[i],
           negBins[i] > 2 ? "  <-- OVERSUBTRACTED" : "");
  }

  // ---- slice-to-slice double ratio ----
  // The pp reference cancels EXACTLY here, so this evolution is independent of
  // pp and of its statistics. Adjacent slices are disjoint event samples.
  int lastGood = -1;
  printf("\n%-16s %12s %12s\n", "slice pair", "double@52", "meaning");
  for(int i = 0; i < NSlice; i++){
    if(!useSlice[i]) continue;
    if(lastGood >= 0){
      prevOf[i] = lastGood;
      hDouble[i] = (TH1D*) hShape[i]->Clone(Form("hDouble_%d", i));
      hDouble[i]->SetDirectory(nullptr);
      hDouble[i]->Divide(hShape[lastGood]);
      double d = hDouble[i]->GetBinContent(hDouble[i]->FindBin(probePt));
      printf("%3.0f-%-3.0f / %3.0f-%-3.0f %12.3f %12s\n",
             i*centStep, (i+1)*centStep, lastGood*centStep, (lastGood+1)*centStep,
             d, d > 1.0 ? "still rising" : "falling");
    }
    lastGood = i;
  }

  for(int i = 0; i < NSlice; i++){
    if(!useSlice[i]) continue;
    for(int b = 1; b <= hRatio[i]->GetNbinsX(); b++){
      double lo = hRatio[i]->GetXaxis()->GetBinLowEdge(b);
      if(lo < plotPtMin || lo >= plotPtMax) continue;
      double v = hRatio[i]->GetBinContent(b);
      if(v != 0. && (v > ratioMax || v < ratioMin))
        printf("NOTE: slice %d at %.0f-%.0f is %.2f, outside display range [%.1f,%.1f]\n",
               i+1, lo, hRatio[i]->GetXaxis()->GetBinUpEdge(b), v, ratioMin, ratioMax);
    }
  }

  const double lm = 0.13, rm = 0.19, tm = 0.09, bm = 0.28, split = 0.40;
  TCanvas *c = new TCanvas("cPbPbOverPPfakeSub", "", 950, 900);

  // ================= top: PbPb/pp per slice =================
  TPad *pUp = new TPad("pUp", "", 0, split, 1, 1);
  pUp->SetLeftMargin(lm); pUp->SetRightMargin(rm);
  pUp->SetTopMargin(tm);  pUp->SetBottomMargin(0.);
  pUp->Draw(); pUp->cd();

  TH1F *fr = pUp->DrawFrame(plotPtMin, ratioMin, plotPtMax, ratioMax);
  fr->GetYaxis()->SetTitle("PbPb / pp  (shape-normalized)");
  fr->GetYaxis()->SetTitleSize(0.055); fr->GetYaxis()->SetLabelSize(0.048);
  fr->GetYaxis()->SetTitleOffset(1.10);
  fr->GetXaxis()->SetLabelSize(0.);

  TLine *one = new TLine(plotPtMin, 1., plotPtMax, 1.);
  one->SetLineStyle(2); one->SetLineColor(kGray+2); one->Draw();
  TLine *nl = new TLine(ptNormMin, ratioMin, ptNormMin, ratioMax);
  nl->SetLineStyle(3); nl->SetLineColor(kGray+1); nl->Draw();
  TLatex nt; nt.SetTextSize(0.034); nt.SetTextColor(kGray+2); nt.SetTextAngle(90);
  nt.DrawLatex(ptNormMin + 5., ratioMax*0.42, Form("normalized above %.0f", ptNormMin));

  for(int i = 0; i < NSlice; i++){
    if(!useSlice[i]) continue;
    hRatio[i]->SetLineColor(col[i]); hRatio[i]->SetMarkerColor(col[i]);
    hRatio[i]->SetMarkerStyle(20); hRatio[i]->SetMarkerSize(0.8);
    hRatio[i]->SetLineWidth(2);
    hRatio[i]->GetXaxis()->SetRangeUser(plotPtMin, plotPtMax);
    hRatio[i]->Draw("ep same");
  }

  TLatex sl; sl.SetNDC(); sl.SetTextSize(0.040); sl.SetTextColor(kBlack);
  sl.DrawLatex(lm + 0.03, 0.86, Form("fake jets subtracted: %s", methLabel[subtractMethod]));

  // ================= bottom: slice-to-slice double ratio =================
  c->cd();
  TPad *pDn = new TPad("pDn", "", 0, 0, 1, split);
  pDn->SetLeftMargin(lm); pDn->SetRightMargin(rm);
  pDn->SetTopMargin(0.);  pDn->SetBottomMargin(bm);
  pDn->Draw(); pDn->cd();

  const double sc = (1.-split)/split;
  TH1F *fr2 = pDn->DrawFrame(plotPtMin, dblMin, plotPtMax, dblMax);
  fr2->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]  (JEC-corrected)");
  fr2->GetYaxis()->SetTitle("slice / previous slice");
  fr2->GetXaxis()->SetTitleSize(0.055*sc); fr2->GetXaxis()->SetLabelSize(0.048*sc);
  fr2->GetYaxis()->SetTitleSize(0.055*sc); fr2->GetYaxis()->SetLabelSize(0.048*sc);
  fr2->GetYaxis()->SetTitleOffset(1.10/sc);
  fr2->GetYaxis()->SetNdivisions(505);

  TLine *one2 = new TLine(plotPtMin, 1., plotPtMax, 1.);
  one2->SetLineStyle(2); one2->SetLineColor(kGray+2); one2->Draw();

  for(int i = 0; i < NSlice; i++){
    if(!hDouble[i]) continue;
    hDouble[i]->SetLineColor(col[i]); hDouble[i]->SetMarkerColor(col[i]);
    hDouble[i]->SetMarkerStyle(20); hDouble[i]->SetMarkerSize(0.8);
    hDouble[i]->SetLineWidth(2);
    hDouble[i]->GetXaxis()->SetRangeUser(plotPtMin, plotPtMax);
    hDouble[i]->Draw("ep same");
  }

  TLatex dl; dl.SetNDC(); dl.SetTextSize(0.048*sc*0.7); dl.SetTextColor(kGray+3);
  dl.DrawLatex(lm + 0.02, 0.90, "above 1 = ratio still rising with centrality");

  // ================= centrality colour bar, spanning both pads =========
  c->cd();
  double gx0 = 1. - rm + 0.030, gx1 = gx0 + 0.035;
  double gy0 = 0.95, gy1 = 0.12;
  for(int i = 0; i < NSlice; i++){
    double y0 = gy0 - (gy0-gy1) * i     / NSlice;
    double y1 = gy0 - (gy0-gy1) * (i+1) / NSlice;
    TPave *box = new TPave(gx0, y1, gx1, y0, 0, "NDC");
    box->SetFillColor(col[i]); box->SetLineWidth(0); box->Draw();
  }
  TPave *bord = new TPave(gx0, gy1, gx1, gy0, 0, "NDC");
  bord->SetFillStyle(0); bord->SetLineColor(kBlack); bord->SetLineWidth(1); bord->Draw();
  TLatex gl; gl.SetNDC(); gl.SetTextSize(0.022); gl.SetTextAlign(12);
  gl.DrawLatex(gx1 + 0.010, gy0, "0%");
  gl.DrawLatex(gx1 + 0.010, (gy0+gy1)/2., "40%");
  gl.DrawLatex(gx1 + 0.010, gy1, "80%");
  TLatex gt; gt.SetNDC(); gt.SetTextSize(0.024); gt.SetTextAngle(90); gt.SetTextAlign(22);
  gt.DrawLatex(gx1 + 0.072, (gy0+gy1)/2., "centrality");

  TLatex lat; lat.SetNDC(); lat.SetTextSize(0.026);
  lat.DrawLatex(0.04, 0.972,
    "PbPb / pp jet p_{T} shape, MinBias only, 5% centrality slices, fake jets subtracted");

  TString out = TString(outDir) + outName;
  c->SaveAs(out);
  printf("\nSaved %s\n", out.Data());

  // ---- method sensitivity at the probe point ----
  // The point of this table is that the central slices are not a measurement
  // until the spread across these three columns is much smaller than it is.
  printf("\n=== method sensitivity: PbPb/pp at %.0f GeV ===\n", probePt);
  printf("%-9s", "cent");
  for(int m = 0; m < NMethod; m++) printf(" %12s", methLabel[m]);
  printf("\n");
  for(int si = 1; si <= NSlice; si++){
    TH1D *measured = nullptr;
    fP->GetObject(Form("%s_C%d", pbpbHistBase, si), measured);
    if(!measured) continue;
    if(measured->Integral(measured->FindBin(ptNormMin + 1e-6),
                          measured->GetNbinsX()+1) < minNormJets) continue;
    printf("%3.0f-%-5.0f", (si-1)*centStep, si*centStep);
    for(int m = 0; m < NMethod; m++){
      int nn = 0;
      TH1D *w = subtracted(fP, si, m, Form("ms_%d_%d", m, si), nn);
      if(!w){ printf(" %12s", "-"); continue; }
      TH1D *s = buildShape(w, Form("mss_%d_%d", m, si));
      delete w;
      if(!s){ printf(" %12s", "-"); continue; }
      s->Divide(hpp);
      printf(" %12.2f", s->GetBinContent(s->FindBin(probePt)));
      delete s;
    }
    printf("\n");
  }
}
