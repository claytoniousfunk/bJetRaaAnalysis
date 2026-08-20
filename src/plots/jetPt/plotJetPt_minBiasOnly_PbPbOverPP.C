// PbPb/pp jet pT shape ratio, MinBias only, in FINE (5%) centrality slices.
//
// Converted from the coarse 4-class data-vs-MC version (preserved in git at
// commit 373909f9) to answer the smoothness question directly: Olga's
// objection is that the coarse 10-30% class looks identical to 0-10%, and the
// cleanest test is whether the ratio evolves smoothly when the same data is
// sliced 5% at a time instead of merged into four wide classes.
//
// *** DATA ONLY -- THE MC CANNOT FOLLOW ***
// The PYTHIA+HYDJET response file carries only four centrality classes
// (h_matchedRecoJetPt_genJetPt_allJets_C1..C4, checked directly), so there is
// no fine-grained MC counterpart. The data-vs-MC comparison therefore stays
// in the coarse version rather than being approximated here.
//
// NORMALISATION, AND WHY 150 RATHER THAN 200 GeV.
// Each curve is scaled so its own integral above ptNormMin is 1, so the ratio
// asks purely about SHAPE. The coarse version used 200 GeV, which is where
// the ratio has genuinely flattened. Per 5% slice, however, the raw jet
// counts above 200 GeV run from 268 (0-5%) down to 1 (75-80%) -- the
// peripheral slices cannot support a normalisation there at all, and since
// every curve is divided by that integral, a handful of jets would rescale a
// whole curve arbitrarily. Above 150 GeV the same counts run 1042 down to 6,
// which carries most slices.
//
// The cost is real and worth stating: at 100-160 GeV the coarse ratio has not
// fully flattened (0-10% still sits around 1.17 at 150-160), so normalising
// at 150 slightly suppresses the most central curves relative to a 200 GeV
// normalisation. That is a common multiplicative shift per curve, not a
// distortion of its shape, and it does not move the peak POSITION or the
// slice-to-slice ORDERING, which is what this plot is for. Set ptNormMin
// back to 200 to check the effect; slices too thin to normalise are skipped
// with a warning rather than drawn.
//
// Usage: root -l -b -q 'plotJetPt_minBiasOnly_PbPbOverPP.C'
// Run from: src/plots/jetPt/

const char *pbpbFile =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/PbPb_MinBias_Part1_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_2026-8-13_ultraFineCentBins.root";
const char *ppFile =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/pp/pp_MinBias_mu12_pTmu-15to999_tight_deltaR-40_jetTrkMaxFilter_WDecayFilter_2026-7-6.root";

// JEC-corrected pT on both sides, as requested.
const char *pbpbHistBase = "h_inclRecoJetPt";
const char *ppHistName   = "h_inclRecoJetPt";

const char *outDir  = "../../../figures/jetPt/";
const char *outName = "jetPt_minBiasOnly_PbPbOverPP_fineCent.pdf";

const int    NSlice   = 16;
const double centStep = 5.0;

const double ptNormMin   = 150.;  // see header for why not 200
const int    minNormJets = 25;    // skip a slice with fewer raw jets than this
                                  // above ptNormMin -- its normalisation would
                                  // be dominated by a handful of counts

const double plotPtMin = 30., plotPtMax = 300.;
const double ratioMin  = 0.,  ratioMax  = 11.;
const double dblMin    = 0.3, dblMax    = 1.9;   // slice/previous-slice panel

// Probe point for the printed smoothness table. Olga: "at 50, this is
// horrendous."
const double probePt = 52.5;

const int    NEdge = 24;
double       ptEdge[NEdge] = {
  20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,
  110,120,140,160,200,250,300
};

// Centrality is an ordered continuous variable, so it gets a perceptually
// uniform gradient rather than 16 categorical hues.
//
// The bottom PALFRACLO of Viridis is near-black and consecutive slices there
// are indistinguishable, which is exactly where the most central slices sit
// -- the region the plot is about. Sampling starts above it so 0-5% opens at
// a clearly blue tone and successive central slices separate visibly. Costs
// a little of the total range; buys back the discrimination that matters.
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

static TH1D* buildShape(TH1D *raw, const char *name, double &rawAboveNorm)
{
  TH1D *h = (TH1D*) raw->Rebin(NEdge-1, name, ptEdge);
  h->SetDirectory(nullptr);
  int b0 = h->FindBin(ptNormMin);
  rawAboveNorm = raw->Integral(raw->FindBin(ptNormMin + 1e-6), raw->GetNbinsX()+1);
  double norm = h->Integral(b0, h->GetNbinsX()+1);
  if(norm <= 0.) return nullptr;
  h->Scale(1./norm);
  h->Scale(1., "width");   // genuine dN/dpT; cancels in a same-binned ratio
  return h;
}

void plotJetPt_minBiasOnly_PbPbOverPP()
{
  gStyle->SetOptStat(0);
  gSystem->mkdir(outDir, kTRUE);

  TFile *fP = TFile::Open(pbpbFile);
  TFile *fp = TFile::Open(ppFile);
  if(!fP || fP->IsZombie()){ printf("ERROR: cannot open %s\n", pbpbFile); return; }
  if(!fp || fp->IsZombie()){ printf("ERROR: cannot open %s\n", ppFile);   return; }

  // ---- pp reference ----
  TH1D *hpp0 = nullptr;
  fp->GetObject(ppHistName, hpp0);
  if(!hpp0){ printf("ERROR: %s not found in pp file\n", ppHistName); return; }
  double ppAbove = 0.;
  TH1D *hpp = buildShape(hpp0, "hpp", ppAbove);
  if(!hpp){ printf("ERROR: pp has zero integral above %.0f GeV\n", ptNormMin); return; }
  printf("pp reference: %.0f raw jets above %.0f GeV\n\n", ppAbove, ptNormMin);

  int col[NSlice];
  sliceColours(col);

  TH1D *hRatio[NSlice], *hShape[NSlice], *hDouble[NSlice];
  bool  useSlice[NSlice];
  double probeVal[NSlice], peakVal[NSlice], peakPt[NSlice];
  int    prevOf[NSlice];              // previous AVAILABLE slice, for the double ratio
  for(int i = 0; i < NSlice; i++){ hShape[i]=nullptr; hDouble[i]=nullptr; prevOf[i]=-1; }

  printf("%-9s %12s %12s %12s %12s\n",
         "cent", "raw>norm", "ratio@50", "peak", "peak pT");
  for(int si = 1; si <= NSlice; si++){
    int i = si - 1;
    useSlice[i] = false;
    probeVal[i] = peakVal[i] = peakPt[i] = -1.;

    TH1D *raw = nullptr;
    fP->GetObject(Form("%s_C%d", pbpbHistBase, si), raw);
    if(!raw){ printf("WARNING: %s_C%d missing\n", pbpbHistBase, si); continue; }

    double aboveNorm = 0.;
    TH1D *shaped = buildShape(raw, Form("hPbPb_C%d", si), aboveNorm);
    if(!shaped || aboveNorm < minNormJets){
      printf("%3.0f-%-5.0f %12.0f   SKIPPED (fewer than %d raw jets above %.0f GeV;\n"
             "                            its normalisation would be noise)\n",
             (si-1)*centStep, si*centStep, aboveNorm, minNormJets, ptNormMin);
      continue;
    }

    hShape[i] = shaped;
    hRatio[i] = (TH1D*) shaped->Clone(Form("hRatio_C%d", si));
    hRatio[i]->SetDirectory(nullptr);
    hRatio[i]->Divide(hpp);
    useSlice[i] = true;

    probeVal[i] = hRatio[i]->GetBinContent(hRatio[i]->FindBin(probePt));
    for(int b = 1; b <= hRatio[i]->GetNbinsX(); b++){
      double lo = hRatio[i]->GetXaxis()->GetBinLowEdge(b);
      if(lo < plotPtMin || lo >= 120.) continue;   // low-pT bump region only
      double v = hRatio[i]->GetBinContent(b);
      if(v > peakVal[i]){ peakVal[i] = v; peakPt[i] = hRatio[i]->GetXaxis()->GetBinCenter(b); }
    }
    printf("%3.0f-%-5.0f %12.0f %12.2f %12.2f %12.0f\n",
           (si-1)*centStep, si*centStep, aboveNorm, probeVal[i], peakVal[i], peakPt[i]);
  }

  // ---- slice-to-slice double ratio ----
  // R_i / R_(i-1) = (PbPb_i/pp) / (PbPb_(i-1)/pp) = PbPb_i / PbPb_(i-1):
  // the pp reference cancels EXACTLY, so this evolution is independent of pp
  // and of its statistics. Computed from the shaped PbPb spectra directly
  // rather than by dividing the two ratios, so the errors do not double-count
  // the (cancelled) pp uncertainty. Adjacent centrality slices are disjoint
  // event samples, hence statistically independent.
  //
  // Chained on the previous AVAILABLE slice, so a skipped slice bridges rather
  // than breaking the chain; the label says which pair each curve is.
  int lastGood = -1;
  printf("\n%-16s %12s %12s\n", "slice pair", "double@50", "meaning");
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

  // warn rather than silently clip
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
  TCanvas *c = new TCanvas("cPbPbOverPPfine", "", 950, 900);

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
  lat.DrawLatex(0.04, 0.972, "PbPb / pp jet p_{T} shape, MinBias only, 5% centrality slices");

  TString out = TString(outDir) + outName;
  c->SaveAs(out);
  printf("\nSaved %s\n", out.Data());
}
