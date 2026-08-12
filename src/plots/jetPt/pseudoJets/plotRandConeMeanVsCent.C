// Plot the overall mean of each UE-estimator map as a function of centrality,
// comparing same-event and mixed-event running.
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
// Series that are identically zero are skipped automatically, which matters for
// the geoCorr maps: pseudoJetPt_geoCorr_k is only accumulated in the same-event
// branch of PbPb_pfCandAnalyzer.C, so those maps are empty in mixed-event files.
//
// Caveat on the mixed-event dPT curve: it is filled, but it is not a UE estimate
// there. The PFCs clustering always runs on the current event's candidates while
// the PF clustering runs on the mixed pool, so dPT subtracts a real jet from a
// mixed-pool fake, and it goes negative beyond ~25% centrality.
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

// One entry per curve: which file, which map. Colour tracks the estimator,
// filled/open marker tracks same vs mixed event.
// Other variants available in the same files:
//   h_dPTEtaPhi_PF_PFCs_dPTAbove0, h_dPTEtaPhi_PF_PFCs_PFCsPTAbove60
// The geoCorr maps are deliberately left out: they are empty in mixed-event
// files, so the comparison this plot is for cannot be made for them. Add
// "h_randConeEtaPhi_geoCorr" back if looking at same-event files only.
const int NSer = 4;
const int   serFile [NSer] = { 0, 1, 0, 1 };
const char *serBase [NSer] = { "h_randConeEtaPhi",    "h_randConeEtaPhi",
			       "h_dPTEtaPhi_PF_PFCs", "h_dPTEtaPhi_PF_PFCs" };
const char *serName [NSer] = { "Random cone", "Random cone",
			       "#Deltap_{T}(PF#minusPFCs)", "#Deltap_{T}(PF#minusPFCs)" };
const int   serMark [NSer] = { 20, 24, 22, 26 };

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
  serCol[0] = serCol[1] = TColor::GetColor("#0072B2");
  serCol[2] = serCol[3] = TColor::GetColor("#009E73");

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

  double ymin = 1e30, ymax = -1e30;
  for(int s = 0; s < NSer; s++){
    if(!use[s]) continue;
    for(int si = 1; si <= NSlice; si++){
      double v = h[s]->GetBinContent(si);
      if(v < ymin) ymin = v;
      if(v > ymax) ymax = v;
    }
  }
  double pad = 0.12 * (ymax - ymin);
  ymin -= pad; ymax += pad;

  TCanvas *c = new TCanvas("cRCMean", "", 700, 600);
  c->SetLeftMargin(0.15);
  c->SetBottomMargin(0.13);
  c->SetTopMargin(0.07);
  c->SetRightMargin(0.05);

  // Single column, upper right: the RC curve has already fallen well below this
  // region by ~20% centrality, so nothing is hidden and no labels collide.
  TLegend *leg = new TLegend(0.44, 0.66, 0.94, 0.89);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.034);
  leg->SetFillStyle(0);

  bool first = true;
  for(int s = 0; s < NSer; s++){
    if(!use[s]) continue;
    h[s]->SetLineColor(serCol[s]);
    h[s]->SetMarkerColor(serCol[s]);
    h[s]->SetMarkerStyle(serMark[s]);
    h[s]->SetMarkerSize(1.2);
    h[s]->SetLineWidth(2);
    h[s]->SetTitle("");
    h[s]->GetXaxis()->SetTitle("centrality [%]");
    h[s]->GetYaxis()->SetTitle("mean UE estimate [GeV]");
    h[s]->GetXaxis()->SetTitleSize(0.045);
    h[s]->GetYaxis()->SetTitleSize(0.045);
    h[s]->GetXaxis()->SetLabelSize(0.040);
    h[s]->GetYaxis()->SetLabelSize(0.040);
    h[s]->GetYaxis()->SetTitleOffset(1.45);
    h[s]->SetMinimum(ymin);
    h[s]->SetMaximum(ymax);
    h[s]->Draw(first ? "ep" : "ep same");
    leg->AddEntry(h[s], Form("%s, %s evt", serName[s], fileLabel[serFile[s]]), "lp");
    first = false;
  }
  leg->Draw();

  if(ymin < 0.){
    TLine *zero = new TLine(0., 0., NSlice * centStep, 0.);
    zero->SetLineStyle(2);
    zero->SetLineColor(kGray + 1);
    zero->Draw();
  }

  TLatex lat; lat.SetNDC(); lat.SetTextSize(0.038);
  lat.DrawLatex(0.18, 0.89, "PbPb 5.02 TeV, R = 0.4");

  TString out = TString(outDir) + outName;
  c->SaveAs(out);
  printf("\nSaved %s\n", out.Data());
}
