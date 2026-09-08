// The 2026-07-27 TODO this answers: "spectra: check fastJet_bkgSub_mixedEvent
// vs. fastJet_bkgSub_sameEvent per-event distributions."
//
// h_fastJetPt_PF_bkgSub_RC is the full FastJet-reclustered, random-cone-
// subtracted fake-jet spectrum -- not a bare random cone. A cone only tests
// whether the mean/width of the UE density is right; this tests whether the
// actual downstream fake-JET rate is right, which is far more sensitive to the
// fluctuation TAIL (a fake jet only exists because some patch of the event
// fluctuated hard enough for anti-kT to cluster it into something jet-shaped).
// If candidate-level mixed-event pooling flattens that tail -- established
// earlier this session via the RC width comparison (mixed/same RMS ratio
// 0.73-0.85) -- this is where it should show up directly, as a shape
// difference in the fake-jet spectrum itself, not just a width number.
//
// Caveat this plot cannot resolve, same as the TODO's likely reason for never
// getting done: on real MinBias data, same-event FastJet output is a mix of
// genuine combinatorial fakes AND real low-pT jets, with no per-event way to
// separate them. Mixed-event has no such contamination (pool candidates never
// contain a real jet's worth of correlated energy). So a same-event excess
// over mixed-event here is upper-bounded by "mixed under-predicts the fake
// rate" but could partly be "same-event also contains real jets mixed-event
// structurally cannot produce." The HYDJET-only closure test sidesteps this by
// using a sample with no hard scatter at all, where same-event output is
// unambiguously 100% fake; this plot is the real-data cross-check against that
// conclusion, not a replacement for it.
//
// Genuine per-event dN/dpT: both curves come from the same underlying MinBias
// sample (N_events differs by <1.3% between the two scans, confirmed before
// writing this), so each is normalized by its own N_events -- no shape trick,
// no cross-sample stitching, a real absolute comparison.
//
// Usage: root -l -b -q 'plotFastJet_bkgSubtraction_mixedVsSame.C'
// Run from: src/plots/jetPt/pseudoJets/

const char *mixedFile =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/PbPb_MinBias_Part1_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_mixedEventPFClustering_pseudoJetCandPtMin-0.0_2026-8-11_ultraFineCentBins.root";
const char *sameFile =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/PbPb_MinBias_Part1_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_sameEventPFClustering_pseudoJetCandPtMin-0.0_2026-8-6_ultraFineCentBins.root";
const char *histBase = "h_fastJetPt_PF_bkgSub_RC";  // raw axis, random-cone subtracted

const char *outDir  = "../../../../figures/jetPt/pseudoJets/";
const char *outName = "fastJet_bkgSub_mixedVsSame.pdf";

// PbPb coarse class -> [first, last] ultra-fine slice index, same convention
// used throughout this session.
const int   NClass = 4;
const int   sliceLo[NClass] = { 1,  3,  7, 11};
const int   sliceHi[NClass] = { 2,  6, 10, 16};
const char *classLabel[NClass] = { "0-10%", "10-30%", "30-50%", "50-80%" };

const double plotPtMin = 0., plotPtMax = 160.;
// mixed/same is bounded in [0,1] by construction here: mixed peaks near 1 at
// threshold (where the spectrum is still combinatorics-dominated on both
// sides) and falls toward 0 as same-event picks up an increasing admixture of
// real QCD jets mixed-event cannot produce. Reads as "fraction of the
// same-event fake yield mixed-event actually recovers."
const double botMin = 0., botMax = 1.1;

// 5 GeV steps where both curves have real content, coarser beyond ~100 GeV
// where statistics thin out fast (single-digit raw counts per 5 GeV bin).
const int    NEdge = 25;
double       ptEdge[NEdge] = {
  0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,90,100,110,120,130,140,150,160
};

// Okabe-Ito, centre-symmetric filled/open pairs per class: open = mixed event,
// filled = same event -- same convention as plotRandConeMeanVsCent.C.
const char *classHex[NClass]   = { "#D55E00", "#0072B2", "#009E73", "#CC79A7" };
const int   markFilled[NClass] = { 20, 21, 33, 34 };
const int   markOpen[NClass]   = { 24, 25, 27, 28 };

static TH1D* sumSlices(TFile *f, const char *base, int ci, const char *name)
{
  TH1D *sum = nullptr;
  for(int si = sliceLo[ci]; si <= sliceHi[ci]; si++){
    TH1D *h = nullptr;
    f->GetObject(Form("%s_C%d", base, si), h);
    if(!h){ printf("ERROR: %s_C%d not found\n", base, si); return nullptr; }
    if(!sum){ sum = (TH1D*) h->Clone(name); sum->SetDirectory(nullptr); }
    else sum->Add(h);
  }
  return sum;
}

static double sumEvents(TFile *f, int ci)
{
  double N = 0.;
  for(int si = sliceLo[ci]; si <= sliceHi[ci]; si++){
    TH1D *hvz = nullptr;
    f->GetObject(Form("h_vz_C%d", si), hvz);
    N += hvz->Integral();
  }
  return N;
}

void plotFastJet_bkgSubtraction_mixedVsSame()
{
  gStyle->SetOptStat(0);
  gSystem->mkdir(outDir, kTRUE);

  TFile *fM = TFile::Open(mixedFile);
  TFile *fS = TFile::Open(sameFile);
  if(!fM || fM->IsZombie()){ printf("ERROR: cannot open %s\n", mixedFile); return; }
  if(!fS || fS->IsZombie()){ printf("ERROR: cannot open %s\n", sameFile);   return; }

  TH1D *hMixed[NClass], *hSame[NClass], *hRatio[NClass];
  int colClass[NClass];
  double ymax = 0.;

  printf("%-8s %10s %10s\n", "class", "N_mixed", "N_same");
  for(int ci = 0; ci < NClass; ci++){
    colClass[ci] = TColor::GetColor(classHex[ci]);

    TH1D *rawMixed = sumSlices(fM, histBase, ci, Form("rawMixed_%d", ci));
    TH1D *rawSame  = sumSlices(fS, histBase, ci, Form("rawSame_%d", ci));
    if(!rawMixed || !rawSame) return;

    double N_mixed = sumEvents(fM, ci);
    double N_same  = sumEvents(fS, ci);
    printf("%-8s %10.0f %10.0f\n", classLabel[ci], N_mixed, N_same);

    hMixed[ci] = (TH1D*) rawMixed->Rebin(NEdge-1, Form("hMixed_%d", ci), ptEdge);
    hSame[ci]  = (TH1D*) rawSame ->Rebin(NEdge-1, Form("hSame_%d",  ci), ptEdge);
    hMixed[ci]->SetDirectory(nullptr); hSame[ci]->SetDirectory(nullptr);

    hMixed[ci]->Scale(1./N_mixed); hMixed[ci]->Scale(1., "width");
    hSame[ci] ->Scale(1./N_same);  hSame[ci] ->Scale(1., "width");

    hRatio[ci] = (TH1D*) hMixed[ci]->Clone(Form("hRatio_%d", ci));
    hRatio[ci]->Divide(hSame[ci]);

    ymax = TMath::Max(ymax, TMath::Max(hMixed[ci]->GetMaximum(), hSame[ci]->GetMaximum()));
  }

  printf("\n%-12s", "pT [GeV]");
  for(int ci = 0; ci < NClass; ci++) printf(" %20s", classLabel[ci]);
  printf("\n");
  for(int b = 1; b <= hMixed[0]->GetNbinsX(); b++){
    printf("%4.0f-%-6.0f", hMixed[0]->GetXaxis()->GetBinLowEdge(b), hMixed[0]->GetXaxis()->GetBinUpEdge(b));
    for(int ci = 0; ci < NClass; ci++)
      printf(" %9.4g->%-9.4g", hMixed[ci]->GetBinContent(b), hSame[ci]->GetBinContent(b));
    printf("\n");
  }

  for(int ci = 0; ci < NClass; ci++){
    for(int b = 1; b <= hRatio[ci]->GetNbinsX(); b++){
      double lo = hRatio[ci]->GetXaxis()->GetBinLowEdge(b);
      if(lo < plotPtMin || lo >= plotPtMax) continue;
      double r = hRatio[ci]->GetBinContent(b);
      if(r != 0. && (r > botMax || r < botMin))
        printf("NOTE: %s mixed/same at %.0f-%.0f is %.2f, outside bottom-panel range [%.1f,%.1f]\n",
               classLabel[ci], lo, hRatio[ci]->GetXaxis()->GetBinUpEdge(b), r, botMin, botMax);
    }
  }

  const double lm = 0.16, rm = 0.05, tm = 0.07, bm = 0.28, split = 0.40;
  TCanvas *c = new TCanvas("cFastJetMixedVsSame", "", 800, 800);

  TPad *pUp = new TPad("pUp", "", 0, split, 1, 1);
  pUp->SetLeftMargin(lm); pUp->SetRightMargin(rm);
  pUp->SetTopMargin(tm);  pUp->SetBottomMargin(0.);
  pUp->SetLogy(); pUp->Draw(); pUp->cd();

  TLegend *leg = new TLegend(0.42, 0.55, 0.94, 0.90);
  leg->SetBorderSize(0); leg->SetTextSize(0.032); leg->SetFillStyle(0); leg->SetNColumns(2);

  bool first = true;
  for(int ci = 0; ci < NClass; ci++){
    hMixed[ci]->SetLineColor(colClass[ci]); hMixed[ci]->SetMarkerColor(colClass[ci]);
    hMixed[ci]->SetMarkerStyle(markOpen[ci]); hMixed[ci]->SetMarkerSize(1.1);
    hMixed[ci]->SetLineWidth(2); hMixed[ci]->SetLineStyle(2);
    hMixed[ci]->SetTitle("");
    hMixed[ci]->GetYaxis()->SetTitle("1/N_{evt} dN_{fakeJet}/dp_{T}  [GeV^{-1}]");
    hMixed[ci]->GetYaxis()->SetTitleSize(0.046);
    hMixed[ci]->GetYaxis()->SetLabelSize(0.040);
    hMixed[ci]->GetYaxis()->SetTitleOffset(1.55);
    hMixed[ci]->GetXaxis()->SetRangeUser(plotPtMin, plotPtMax);
    hMixed[ci]->SetMinimum(ymax * 1e-5); hMixed[ci]->SetMaximum(ymax * 3.);
    hMixed[ci]->Draw(first ? "ep" : "ep same");
    first = false;

    hSame[ci]->SetLineColor(colClass[ci]); hSame[ci]->SetMarkerColor(colClass[ci]);
    hSame[ci]->SetMarkerStyle(markFilled[ci]); hSame[ci]->SetMarkerSize(1.0);
    hSame[ci]->SetLineWidth(2);
    hSame[ci]->Draw("ep same");

    leg->AddEntry(hMixed[ci], Form("%s, mixed", classLabel[ci]), "lp");
    leg->AddEntry(hSame[ci],  Form("%s, same",  classLabel[ci]), "lp");
  }
  leg->Draw();

  TLatex lat; lat.SetNDC(); lat.SetTextSize(0.036);
  lat.DrawLatex(0.20, 0.955, "FastJet fake-jet spectrum, raw p_{T}, RC-subtracted");

  c->cd();
  TPad *pDn = new TPad("pDn", "", 0, 0, 1, split);
  pDn->SetLeftMargin(lm); pDn->SetRightMargin(rm);
  pDn->SetTopMargin(0.);  pDn->SetBottomMargin(bm);
  pDn->Draw(); pDn->cd();

  const double sc = (1.-split)/split;
  bool firstB = true;
  for(int ci = 0; ci < NClass; ci++){
    hRatio[ci]->SetLineColor(colClass[ci]); hRatio[ci]->SetMarkerColor(colClass[ci]);
    hRatio[ci]->SetMarkerStyle(markFilled[ci]); hRatio[ci]->SetMarkerSize(1.0);
    hRatio[ci]->SetLineWidth(2);
    hRatio[ci]->SetTitle("");
    hRatio[ci]->GetXaxis()->SetTitle("raw p_{T}^{fakeJet} [GeV]");
    hRatio[ci]->GetYaxis()->SetTitle("mixed / same");
    hRatio[ci]->GetXaxis()->SetTitleSize(0.046*sc);
    hRatio[ci]->GetXaxis()->SetLabelSize(0.040*sc);
    hRatio[ci]->GetYaxis()->SetTitleSize(0.046*sc);
    hRatio[ci]->GetYaxis()->SetLabelSize(0.040*sc);
    hRatio[ci]->GetYaxis()->SetTitleOffset(1.55/sc);
    hRatio[ci]->GetYaxis()->SetNdivisions(505);
    hRatio[ci]->GetXaxis()->SetRangeUser(plotPtMin, plotPtMax);
    hRatio[ci]->SetMinimum(botMin); hRatio[ci]->SetMaximum(botMax);
    hRatio[ci]->Draw(firstB ? "ep" : "ep same");
    firstB = false;
  }

  TLine *one = new TLine(plotPtMin, 1., plotPtMax, 1.);
  one->SetLineStyle(2); one->SetLineColor(kGray+1); one->Draw();

  TString out = TString(outDir) + outName;
  c->SaveAs(out);
  printf("\nSaved %s\n", out.Data());
}
