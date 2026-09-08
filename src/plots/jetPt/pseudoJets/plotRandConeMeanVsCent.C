// Plot the overall mean of each UE-estimator map as a function of centrality,
// comparing same-event and mixed-event running, with a ratio to the mixed-event
// random cone.
//
// The UE maps are TProfile2D in (eta, phi); their z value is the mean estimate
// in each cell. This collapses each map to a single number -- the entry-weighted
// average over all cells, i.e. the mean over the whole acceptance -- and plots
// it against centrality.
//
// Entry weighting (rather than a plain average over cells) is what makes the
// result equal to the mean over all thrown cones / matched jets: cells are not
// populated equally, and an unweighted cell average would over-count sparse
// cells.
//
// The lower panel divides each series by the mixed-event random cone, which is
// the estimator the analysis actually subtracts. Numerator and denominator always
// come from different files -- same- and mixed-event running are separate scans --
// so they are statistically independent and the uncertainties are propagated
// uncorrelated. The denominator's own ratio is identically 1 and is drawn as the
// dashed reference line rather than as a series.
//
// Series that are identically zero are skipped automatically, which matters if
// this is ever pointed at a file written before pseudoJetPt_geoCorr_k was
// accumulated in both branches of PbPb_pfCandAnalyzer.C -- in those files the
// geoCorr map is empty in mixed-event running (but not same-event, which is
// where it is used below).
//
// The mixed-event dPT series is deliberately not plotted. It is filled, but it is
// not a UE estimate there: the PFCs clustering always runs on the current event's
// candidates while the PF clustering runs on the mixed pool, so dPT subtracts a
// real jet from a mixed-pool fake, and it goes negative beyond ~25% centrality.
//
// Usage: root -l -b -q 'plotRandConeMeanVsCent.C'
// Run from: src/plots/jetPt/pseudoJets/

const int NFiles = 2;
const char *inFile[NFiles] = {
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/PbPb_MinBias_Part1_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_mixedEventPFClustering_pseudoJetCandPtMin-0.0_2026-8-11_ultraFineCentBins.root",
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/PbPb_MinBias_Part1_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_sameEventPFClustering_pseudoJetCandPtMin-0.0_2026-8-6_ultraFineCentBins.root"
};
const char *fileLabel[NFiles] = { "mixed", "same" };

const char *outDir =
  "../../../../figures/jetPt/pseudoJets/randConeMeanVsCent/";
const char *outName = "randConeMeanVsCent.pdf";

// ultra-fine slices C1..C16, each 10 hiBin units = 5% centrality
const int    NSlice   = 16;
const double centStep = 5.0;

// One entry per curve: which file, which map.
//
// Colour, marker shape and marker size all differ between series, so no
// distinction relies on colour alone. The two random-cone curves agree to a few
// parts in ten thousand and would otherwise be indistinguishable: the same-event
// one is drawn as an open square slightly larger than the filled mixed-event
// circle, so that circle sits visibly nested inside it wherever they coincide.
//
// Other variants available in the same files:
//   h_dPTEtaPhi_PF_PFCs_dPTAbove0, h_dPTEtaPhi_PF_PFCs_PFCsPTAbove60
const int NSer = 4;
const int   serFile [NSer] = { 0, 1, 1, 1 };
const char *serBase [NSer] = { "h_randConeEtaPhi", "h_randConeEtaPhi",
			       "h_dPTEtaPhi_PF_PFCs", "h_randConeEtaPhi_geoCorr" };
const char *serName [NSer] = { "Random cone", "Random cone",
			       "#Deltap_{T}(PF#minusPFCs)", "RC, geoCorr" };
// Centre-symmetric shapes only -- no triangles, whose visual centre of mass sits
// off the point they mark. geoCorr is same-event like the RC-same and dPT
// series, so it gets an open marker too, matching that convention.
const int   serMark [NSer] = { 20, 25, 27, 28 };
const double serSize[NSer] = { 1.0, 1.5, 1.7, 1.5 };

// Okabe-Ito: blue / vermillion / bluish green / reddish purple. Chosen so all
// four stay mutually distinguishable under deuteranopia, protanopia and
// tritanopia, and so they separate in lightness as well as hue for greyscale
// printing.
const char *serHex [NSer] = { "#0072B2", "#D55E00", "#009E73", "#CC79A7" };

// Ratio denominator: the mixed-event random cone.
const int serDenom = 0;

// Fixed ratio range: wide enough for every series in these files, and stable
// across reruns so plots from different scans stay directly comparable.
const double ratioMin = 0.8, ratioMax = 1.2;

// Entry-weighted mean of a TProfile2D's z values, with the propagated
// uncertainty on that mean. Returns false if the map has no entries.
static bool profileMean(TProfile2D *p, double &mean, double &err)
{
  double sw = 0., swy = 0., var = 0.;
  for(int ix = 1; ix <= p->GetNbinsX(); ix++){
    for(int iy = 1; iy <= p->GetNbinsY(); iy++){
      int b = p->GetBin(ix, iy);
      double n = p->GetBinEntries(b);
      if(n <= 0) continue;
      double y = p->GetBinContent(b);
      double e = p->GetBinError(b);   // error on the cell mean
      sw  += n;
      swy += n * y;
      var += n * n * e * e;
    }
  }
  if(sw <= 0) return false;
  mean = swy / sw;
  err  = TMath::Sqrt(var) / sw;
  return true;
}

void plotRandConeMeanVsCent()
{
  gStyle->SetOptStat(0);
  gSystem->mkdir(outDir, kTRUE);

  TFile *f[NFiles];
  for(int i = 0; i < NFiles; i++){
    f[i] = TFile::Open(inFile[i]);
    if(!f[i] || f[i]->IsZombie()){ printf("ERROR: cannot open %s\n", inFile[i]); return; }
  }

  int serCol[NSer];
  for(int s = 0; s < NSer; s++) serCol[s] = TColor::GetColor(serHex[s]);

  TH1D *h[NSer];
  for(int s = 0; s < NSer; s++){
    h[s] = new TH1D(Form("hMean_%d", s), "", NSlice, 0., NSlice * centStep);
    h[s]->SetDirectory(nullptr);
  }

  printf("\n%-9s", "cent [%]");
  for(int s = 0; s < NSer; s++) printf(" %14s/%-6s", serName[s], fileLabel[serFile[s]]);
  printf("\n");

  for(int si = 1; si <= NSlice; si++){
    printf("%3.0f-%-5.0f", (si-1)*centStep, si*centStep);
    for(int s = 0; s < NSer; s++){
      TProfile2D *p = nullptr;
      f[serFile[s]]->GetObject(Form("%s_C%d", serBase[s], si), p);
      double mean = 0., err = 0.;
      if(!p || !profileMean(p, mean, err)){ printf(" %21s", "--"); continue; }
      h[s]->SetBinContent(si, mean);
      h[s]->SetBinError  (si, err);
      printf(" %21.4f", mean);
    }
    printf("\n");
  }

  bool use[NSer];
  int nUsed = 0;
  for(int s = 0; s < NSer; s++){
    use[s] = false;
    for(int si = 1; si <= NSlice; si++)
      if(h[s]->GetBinContent(si) != 0.){ use[s] = true; break; }
    if(use[s]) nUsed++;
    else printf("NOTE: %s (%s event) is identically zero -- not plotted\n",
		serName[s], fileLabel[serFile[s]]);
  }
  if(nUsed == 0){ printf("ERROR: nothing to plot\n"); return; }

  // --- ratios to the mixed-event random cone ------------------------------
  // Independent files, so the relative errors add in quadrature.
  const bool doRatio = use[serDenom];
  if(!doRatio)
    printf("NOTE: %s (%s event) unusable as denominator -- ratio panel empty\n",
	   serName[serDenom], fileLabel[serFile[serDenom]]);

  TH1D *r[NSer];
  bool useR[NSer];
  for(int s = 0; s < NSer; s++){
    r[s] = new TH1D(Form("hRatio_%d", s), "", NSlice, 0., NSlice * centStep);
    r[s]->SetDirectory(nullptr);
    useR[s] = false;
    if(!doRatio || s == serDenom || !use[s]) continue;
    for(int si = 1; si <= NSlice; si++){
      double a = h[s]       ->GetBinContent(si), ea = h[s]       ->GetBinError(si);
      double b = h[serDenom]->GetBinContent(si), eb = h[serDenom]->GetBinError(si);
      if(a == 0. || b == 0.) continue;
      double ratio = a / b;
      r[s]->SetBinContent(si, ratio);
      r[s]->SetBinError  (si, TMath::Abs(ratio) *
			  TMath::Sqrt((ea/a)*(ea/a) + (eb/b)*(eb/b)));
      useR[s] = true;
    }
  }

  // The axis is anchored at zero so the fall-off with centrality is read against
  // a true baseline rather than a suppressed one. Warn if that would clip a
  // series: every estimator here is positive, but the mixed-event dPT variant
  // (not plotted) goes negative, so a future series swap could hit this.
  double ymin = 0., ymax = -1e30;
  for(int s = 0; s < NSer; s++){
    if(!use[s]) continue;
    for(int si = 1; si <= NSlice; si++){
      double v = h[s]->GetBinContent(si);
      if(v > ymax) ymax = v;
      if(v < 0.)
	printf("NOTE: %s/%s at %.0f-%.0f%% is %.4f, clipped by the zero baseline\n",
	       serName[s], fileLabel[serFile[s]], (si-1)*centStep, si*centStep, v);
    }
  }
  ymax *= 1.12;

  // Warn rather than silently clip if a point falls outside the fixed range.
  for(int s = 0; s < NSer; s++){
    if(!useR[s]) continue;
    for(int si = 1; si <= NSlice; si++){
      double v = r[s]->GetBinContent(si);
      if(v == 0.) continue;
      if(v < ratioMin || v > ratioMax)
	printf("NOTE: %s/%s at %.0f-%.0f%% is %.4f, outside the ratio range [%.2f, %.2f]\n",
	       serName[s], fileLabel[serFile[s]], (si-1)*centStep, si*centStep,
	       v, ratioMin, ratioMax);
    }
  }

  const double lm = 0.15, rm = 0.05, tm = 0.07, bm = 0.28;
  const double split = 0.35;

  TCanvas *c = new TCanvas("cRCMean", "", 700, 800);

  TPad *pUp = new TPad("pUpRCMean", "", 0, split, 1, 1);
  pUp->SetLeftMargin(lm); pUp->SetRightMargin(rm);
  pUp->SetTopMargin(tm);  pUp->SetBottomMargin(0.);
  pUp->Draw(); pUp->cd();

  // Single column, upper right: the RC curve has already fallen well below this
  // region by ~20% centrality, so nothing is hidden and no labels collide.
  TLegend *leg = new TLegend(0.40, 0.60, 0.97, 0.88);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.040);
  leg->SetFillStyle(0);

  bool first = true;
  TH1D *hFirst = nullptr;
  for(int s = 0; s < NSer; s++){
    if(!use[s]) continue;
    h[s]->SetLineColor(serCol[s]);
    h[s]->SetMarkerColor(serCol[s]);
    h[s]->SetMarkerStyle(serMark[s]);
    h[s]->SetMarkerSize(serSize[s]);
    h[s]->SetLineWidth(2);
    h[s]->SetTitle("");
    h[s]->GetYaxis()->SetTitle("mean UE estimate [GeV]");
    h[s]->GetYaxis()->SetTitleSize(0.055);
    h[s]->GetYaxis()->SetLabelSize(0.050);
    h[s]->GetYaxis()->SetTitleOffset(1.20);
    h[s]->SetMinimum(ymin);
    h[s]->SetMaximum(ymax);
    h[s]->Draw(first ? "ep" : "ep same");
    leg->AddEntry(h[s], Form("%s, %s evt", serName[s], fileLabel[serFile[s]]), "lp");
    if(first) hFirst = h[s];
    first = false;
  }

  // The upper pad has no bottom margin, so with the axis anchored at zero the
  // "0" label straddles the pad boundary and is drawn clipped. Suppress it; the
  // ratio panel sits directly below and carries the shared x-axis anyway.
  if(hFirst) hFirst->GetYaxis()->ChangeLabel(1, -1, 0.);

  leg->Draw();

  TLatex lat; lat.SetNDC(); lat.SetTextSize(0.048);
  lat.DrawLatex(0.18, 0.87, "PbPb 5.02 TeV, R = 0.4");

  c->cd();
  TPad *pDn = new TPad("pDnRCMean", "", 0, 0, 1, split);
  pDn->SetLeftMargin(lm); pDn->SetRightMargin(rm);
  pDn->SetTopMargin(0.);  pDn->SetBottomMargin(bm);
  pDn->Draw(); pDn->cd();

  // Font sizes scale with pad height, so the lower pad needs them larger by
  // (1-split)/split to match the upper pad on the page.
  const double sc = (1. - split) / split;

  bool firstR = true;
  for(int s = 0; s < NSer; s++){
    if(!useR[s]) continue;
    r[s]->SetLineColor(serCol[s]);
    r[s]->SetMarkerColor(serCol[s]);
    r[s]->SetMarkerStyle(serMark[s]);
    r[s]->SetMarkerSize(serSize[s]);
    r[s]->SetLineWidth(2);
    r[s]->SetTitle("");
    r[s]->GetXaxis()->SetTitle("centrality [%]");
    // Deliberately shorter and smaller than the size-matched value: spelling the
    // denominator out in full overruns the short pad and runs off the page.
    r[s]->GetYaxis()->SetTitle("ratio to mixed");
    r[s]->GetXaxis()->SetTitleSize(0.050 * sc);
    r[s]->GetXaxis()->SetLabelSize(0.045 * sc);
    r[s]->GetYaxis()->SetTitleSize(0.042 * sc);
    r[s]->GetYaxis()->SetLabelSize(0.045 * sc);
    r[s]->GetYaxis()->SetTitleOffset(1.55 / sc);
    r[s]->GetYaxis()->SetNdivisions(505);
    r[s]->SetMinimum(ratioMin);
    r[s]->SetMaximum(ratioMax);
    r[s]->Draw(firstR ? "ep" : "ep same");
    firstR = false;
  }

  if(!firstR){
    TLine *one = new TLine(0., 1., NSlice * centStep, 1.);
    one->SetLineStyle(2);
    one->SetLineColor(kGray + 1);
    one->Draw();
  }

  // --- ratio table --------------------------------------------------------
  printf("\nratio to %s, %s event\n", serName[serDenom], fileLabel[serFile[serDenom]]);
  printf("%-9s", "cent [%]");
  for(int s = 0; s < NSer; s++){
    if(!useR[s]) continue;
    printf(" %14s/%-6s", serName[s], fileLabel[serFile[s]]);
  }
  printf("\n");
  for(int si = 1; si <= NSlice; si++){
    printf("%3.0f-%-5.0f", (si-1)*centStep, si*centStep);
    for(int s = 0; s < NSer; s++){
      if(!useR[s]) continue;
      double v = r[s]->GetBinContent(si);
      if(v == 0.){ printf(" %21s", "--"); continue; }
      printf(" %15.4f +- %-4.4f", v, r[s]->GetBinError(si));
    }
    printf("\n");
  }

  TString out = TString(outDir) + outName;
  c->SaveAs(out);
  printf("\nSaved %s\n", out.Data());
}
