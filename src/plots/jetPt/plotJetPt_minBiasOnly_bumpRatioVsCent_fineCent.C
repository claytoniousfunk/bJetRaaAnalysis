// Why the centrality hierarchy is not visible in the shape-normalized
// PbPb ratio: the observable is FLAT across 0-20%, then falls off a cliff.
//
// Follow-up to plotJetPt_minBiasOnly_rawVsJEC_RCP.C (task #1, 2026-08-18
// meeting), which found 0-10% and 10-30% degenerate in BOTH raw and
// JEC-corrected pT -- ruling out the JEC as the cause.  This plot shows why.
//
// The "bump height" that Olga reads off the slide-2 plot is, in essence,
//   Y(low pT)  /  Y(high pT normalization region)
// i.e. a fake/UE-driven numerator over an Ncoll/hard-scatter-driven
// denominator.  Both fall steeply with centrality, and dividing one steep
// function by another gives something far flatter than either.
//
// Top panel shows both ingredients as raw per-event yields (no rescaling to
// any reference slice), so the absolute rates are directly readable. The
// mechanism is then read off the SLOPES: across 0-20% the two fall at almost
// the same rate, which is what flattens their ratio.
// Bottom panel is their ratio, with the yield-weighted coarse-class averages
// overlaid as horizontal segments -- 0-10% and 10-30% land on essentially
// the same value because 0-10% sits entirely inside the flat region and
// 10-30% is half inside it, with the flat half carrying most of the jets.
//
// Reference window is > 150 GeV rather than the > 200 GeV used for the
// actual shape normalization: same quantity, meaningfully better statistics
// per 5% slice (238 raw jets in the most central slice at >150, versus far
// fewer at >200).  The >200 version is statistically fragile enough per
// slice that it obscures the shape this plot exists to show.  The 75-80%
// slice has zero jets above 150 GeV and is dropped rather than plotted as
// a divide-by-zero.
//
// Usage: root -l -b -q 'plotJetPt_minBiasOnly_bumpRatioVsCent_fineCent.C'
// Run from: src/plots/jetPt/

const char *pbpbFile =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/PbPb_MinBias_Part1_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_sameEventPFClustering_pseudoJetCandPtMin-0.0_2026-8-18_ultraFineCentBins.root";

const char *outDir  = "../../../figures/jetPt/";
const char *outName = "jetPt_minBiasOnly_bumpRatioVsCent_fineCent.pdf";

const int    NSlice   = 16;
const double centStep = 5.0;

// Low-pT probe: the 50-55 GeV bin, where Olga says the anomaly is "in your
// face".  High-pT reference: everything above refPtMin.
const double lowPt    = 52.5;
const double refPtMin = 150.;

// Coarse classes, for the overlaid merged values.
const int   NClass = 4;
const int   sliceLo[NClass] = { 1,  3,  7, 11};
const int   sliceHi[NClass] = { 2,  6, 10, 16};
const char *classLabel[NClass] = { "0-10%", "10-30%", "30-50%", "50-80%" };

const char *lowHex = "#0072B2", *refHex = "#D55E00";

void plotJetPt_minBiasOnly_bumpRatioVsCent_fineCent()
{
  gStyle->SetOptStat(0);
  gSystem->mkdir(outDir, kTRUE);

  TFile *f = TFile::Open(pbpbFile);
  if(!f || f->IsZombie()){ printf("ERROR: cannot open %s\n", pbpbFile); return; }

  double cent[NSlice], y50[NSlice], e50[NSlice], yhi[NSlice], ehi[NSlice], nev[NSlice];
  bool   ok[NSlice];

  printf("%-9s %10s %12s %12s %12s\n", "cent", "N_evt", "Y50/evt", "Yref/evt", "Y50/Yref");
  for(int si = 1; si <= NSlice; si++){
    int i = si - 1;
    cent[i] = (si - 0.5) * centStep;

    TH1D *h = nullptr, *v = nullptr;
    f->GetObject(Form("h_inclRecoJetPt_C%d", si), h);
    f->GetObject(Form("h_vz_C%d", si), v);
    if(!h || !v){ printf("ERROR: missing histograms for slice %d\n", si); return; }

    double N = v->Integral();
    nev[i] = N;

    int b50 = h->FindBin(lowPt);
    y50[i] = h->GetBinContent(b50) / N;
    e50[i] = h->GetBinError(b50)   / N;

    double errHi = 0.;
    int bhi = h->FindBin(refPtMin + 1e-6);
    double intHi = h->IntegralAndError(bhi, h->GetNbinsX()+1, errHi);
    yhi[i] = intHi / N;
    ehi[i] = errHi / N;

    ok[i] = (yhi[i] > 0. && y50[i] > 0.);
    printf("%3.0f-%-5.0f %10.0f %12.5f %12.6f %12.2f\n",
           (si-1)*centStep, si*centStep, N, y50[i], yhi[i], ok[i] ? y50[i]/yhi[i] : -1.);
  }

  // Per-slice ratio, with both ingredients' statistical errors propagated
  // (disjoint pT regions of the same histogram, so treated as independent).
  TGraphErrors *gRatio = new TGraphErrors();
  int np = 0;
  for(int i = 0; i < NSlice; i++){
    if(!ok[i]) continue;
    double r  = y50[i] / yhi[i];
    double er = r * TMath::Sqrt(TMath::Power(e50[i]/y50[i], 2) + TMath::Power(ehi[i]/yhi[i], 2));
    gRatio->SetPoint(np, cent[i], r);
    gRatio->SetPointError(np, centStep/2., er);
    np++;
  }

  // Ingredients as raw per-event yields, not rescaled to any reference slice.
  // The two now sit on genuinely different scales (the 50 GeV yield is ~300x
  // the >150 GeV yield), so the "they fall together" point is read off the
  // slopes rather than off overlapping curves -- and the absolute rates are
  // directly legible, which the normalised version hid.
  TGraphErrors *gLow = new TGraphErrors(), *gRef = new TGraphErrors();
  int nl = 0, nr = 0;
  double yTop = 0., yBot = 1e300;
  for(int i = 0; i < NSlice; i++){
    if(y50[i] > 0.){
      gLow->SetPoint(nl, cent[i], y50[i]);
      gLow->SetPointError(nl, centStep/2., e50[i]);
      nl++;
      if(y50[i] > yTop) yTop = y50[i];
      if(y50[i] < yBot) yBot = y50[i];
    }
    if(yhi[i] > 0.){
      gRef->SetPoint(nr, cent[i], yhi[i]);
      gRef->SetPointError(nr, centStep/2., ehi[i]);
      nr++;
      if(yhi[i] > yTop) yTop = yhi[i];
      if(yhi[i] < yBot) yBot = yhi[i];
    }
  }
  const double frLo = yBot * 0.25, frHi = yTop * 6.;

  // Yield-weighted coarse-class merges -- what sumSlices() actually produces.
  double classRatio[NClass], classEffCent[NClass];
  printf("\n%-9s %12s %12s %12s  %s\n", "class", "Y50/evt", "Yref/evt", "Y50/Yref", "eff.<cent>");
  for(int c = 0; c < NClass; c++){
    double sY50 = 0., sYhi = 0., sN = 0., wsum = 0., wc = 0.;
    for(int si = sliceLo[c]; si <= sliceHi[c]; si++){
      int i = si - 1;
      sY50 += y50[i]*nev[i]; sYhi += yhi[i]*nev[i]; sN += nev[i];
      double w = yhi[i]*nev[i]; wsum += w; wc += w*cent[i];
    }
    classRatio[c]   = sY50 / sYhi;
    classEffCent[c] = wc / wsum;
    printf("%-9s %12.5f %12.6f %12.2f  %9.1f%%\n",
           classLabel[c], sY50/sN, sYhi/sN, classRatio[c], classEffCent[c]);
  }
  printf("\n0-10%% / 10-30%% = %.3f   (a visible hierarchy would give ~1.3-1.5)\n",
         classRatio[0]/classRatio[1]);

  int colLow = TColor::GetColor(lowHex), colRef = TColor::GetColor(refHex);

  const double lm = 0.14, rm = 0.04, tm = 0.08, bm = 0.30, split = 0.42;
  TCanvas *c = new TCanvas("cBumpRatioVsCent", "", 900, 850);

  // ---------------- top: the two ingredients ----------------
  TPad *pUp = new TPad("pUp", "", 0, split, 1, 1);
  pUp->SetLeftMargin(lm); pUp->SetRightMargin(rm);
  pUp->SetTopMargin(tm);  pUp->SetBottomMargin(0.);
  pUp->SetLogy(); pUp->Draw(); pUp->cd();

  TH1F *frUp = pUp->DrawFrame(0., frLo, 80., frHi);
  frUp->GetYaxis()->SetTitle("jets / event");
  frUp->GetYaxis()->SetTitleSize(0.055); frUp->GetYaxis()->SetLabelSize(0.048);
  frUp->GetYaxis()->SetTitleOffset(1.20);
  frUp->GetXaxis()->SetLabelSize(0.);

  for(double xb : {10., 30., 50.}){
    TLine *l = new TLine(xb, frLo, xb, frHi);
    l->SetLineStyle(2); l->SetLineColor(kGray+1); l->Draw();
  }

  gLow->SetMarkerStyle(20); gLow->SetMarkerColor(colLow); gLow->SetLineColor(colLow);
  gLow->SetMarkerSize(1.1); gLow->SetLineWidth(2); gLow->Draw("pz same");

  gRef->SetMarkerStyle(21); gRef->SetMarkerColor(colRef); gRef->SetLineColor(colRef);
  gRef->SetMarkerSize(1.1); gRef->SetLineWidth(2); gRef->Draw("pz same");

  TLegend *leg = new TLegend(0.50, 0.60, 0.95, 0.84);
  leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.048);
  leg->AddEntry(gLow, "low p_{T}: 50-55 GeV", "lp");
  leg->AddEntry(gRef, Form("reference: p_{T} > %.0f GeV", refPtMin), "lp");
  leg->Draw();

  TLatex lat; lat.SetNDC(); lat.SetTextSize(0.052);
  lat.DrawLatex(lm + 0.02, 0.955, "PbPb MinBias, inclusive reco jets, no subtraction");

  // ---------------- bottom: their ratio ----------------
  c->cd();
  TPad *pDn = new TPad("pDn", "", 0, 0, 1, split);
  pDn->SetLeftMargin(lm); pDn->SetRightMargin(rm);
  pDn->SetTopMargin(0.);  pDn->SetBottomMargin(bm);
  pDn->Draw(); pDn->cd();

  const double sc = (1.-split)/split;
  TH1F *frDn = pDn->DrawFrame(0., 0., 80., 545.);
  frDn->GetXaxis()->SetTitle("centrality [%]");
  frDn->GetYaxis()->SetTitle("Y(50 GeV) / Y(>150 GeV)");
  frDn->GetXaxis()->SetTitleSize(0.055*sc); frDn->GetXaxis()->SetLabelSize(0.048*sc);
  frDn->GetYaxis()->SetTitleSize(0.055*sc); frDn->GetYaxis()->SetLabelSize(0.048*sc);
  frDn->GetYaxis()->SetTitleOffset(1.20/sc);
  frDn->GetYaxis()->SetNdivisions(505);

  for(double xb : {10., 30., 50.}){
    TLine *l = new TLine(xb, 0., xb, 545.);
    l->SetLineStyle(2); l->SetLineColor(kGray+1); l->Draw();
  }

  // Coarse-class merged values, drawn across each class's centrality span.
  // Value labels sit on a common row above all the data so they cannot
  // collide with the fine-slice points.
  for(int cl = 0; cl < NClass; cl++){
    double xlo = (sliceLo[cl]-1)*centStep, xhi = sliceHi[cl]*centStep;
    TLine *seg = new TLine(xlo, classRatio[cl], xhi, classRatio[cl]);
    seg->SetLineColor(kBlack); seg->SetLineWidth(3); seg->Draw();
    TLatex cl2; cl2.SetTextSize(0.048*sc); cl2.SetTextAlign(21); cl2.SetTextFont(62);
    cl2.DrawLatex((xlo+xhi)/2., 487., Form("%.0f", classRatio[cl]));
  }

  gRatio->SetMarkerStyle(20); gRatio->SetMarkerColor(kBlack); gRatio->SetLineColor(kBlack);
  gRatio->SetMarkerSize(1.0); gRatio->SetLineWidth(2); gRatio->Draw("pz same");

  TLatex nb; nb.SetNDC(); nb.SetTextSize(0.042*sc);
  nb.DrawLatex(0.46, 0.60, "black = yield-weighted");
  nb.DrawLatex(0.46, 0.52, "coarse-class merge");

  TString out = TString(outDir) + outName;
  c->SaveAs(out);
  printf("\nSaved %s\n", out.Data());
}
