// X-projection of h_muonDR_inclusiveClosestJet_C* (dR vs. closest-jet pT, all
// muons -- not restricted to the muon's own tagged jet) in three jet pT
// windows, merged into the standard 4 coarse centrality classes.
//
// One canvas per coarse class, three overlaid dR distributions (one per jet pT
// window), each divided by the number of events in that coarse class and by bin
// width, giving dN/d(dR) PER EVENT. This is the same normalisation the ptRel
// decomposition uses: sum h_vz over the class's fine slices and divide. Because
// every histogram here is filled at w_resample = w/N_resamples, one real event
// contributes total weight w regardless of the resample count, so the per-event
// rate is directly comparable across scans with different N.
//
// These curves were previously unit-normalised on the grounds that the absolute
// yields spanned two orders of magnitude. They do not: the per-event rates
// across the three pT windows differ by only 1.9x (30-50%) to 4.7x (0-10%), so
// all three sit legibly on one linear axis while now carrying rate information
// that unit normalisation threw away.
//
// Statistics from the 100-resample scan (fastJetResamples-100, 2026-09-08) are
// roughly two orders of magnitude above the earlier single-sample version:
// merged-class entry counts run from 31975 (0-10%, 50-60 GeV) down to 490
// (50-80%, 80-120 GeV), against 364 down to 4 before. The legend still reports
// raw entries per curve.
//
// dR rebinned into 10 uniform 0.05-wide bins. The earlier 0.1-wide binning was
// forced by those old counts and is no longer necessary; the native axis is
// 0.0025 wide, so any multiple of that is a legal edge. Note the peripheral
// panels remain thin at LOW dR specifically: for a flat density the entries go
// as dR, so the first bin holds only ~1% of a curve's total -- about 5 entries
// in the sparsest panel even though the curve as a whole has 490.
//
// Usage: root -l -b -q 'plotMuonDR_inclusiveClosestJet_coarseCent.C'
// Run from: src/plots/muJetDr/

#include "../../../headers/functions/divideByBinwidth.h"

const char *inFile =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/PbPb_SingleMuon_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_mixedEventPFClustering_fastJetResamples-100_2026-9-8_ultraFineCentBins.root";

// Prefer the trigger-gated histogram when the scan has it, so this figure can
// be compared directly against plotMixedMuonPtRel_recoJet_coarseCent.C, whose
// fill is gated. Scans before 2026-09-08 carry only the ungated one; the macro
// falls back to it and says so, and normalises by the matching denominator in
// each case -- h_vz_triggerOn for the gated histogram, h_vz for the ungated.
const char *histBaseGated   = "h_muonDR_inclusiveClosestJet_triggerOn";
const char *histBaseUngated = "h_muonDR_inclusiveClosestJet";

const char *outDir = "../../../figures/muJetDr/";

// coarse class -> [first, last] ultra-fine slice index, standard convention
const int NClass = 4;
const int sliceLo[NClass] = { 1,  3,  7, 11};
const int sliceHi[NClass] = { 2,  6, 10, 16};
const char *classLabel[NClass] = { "0-10%", "10-30%", "30-50%", "50-80%" };

const int    NPt = 3;
const double ptLo[NPt] = {50, 60,  80};
const double ptHi[NPt] = {60, 80, 120};
const char  *ptHex[NPt]   = {"#0072B2", "#D55E00", "#009E73"};   // Okabe-Ito
const int    ptMarker[NPt] = {20, 21, 33};                        // circle, square, diamond -- no triangles

const int    NEdge_dR = 11;
double       edge_dR[NEdge_dR] = {0, 0.05, 0.10, 0.15, 0.20, 0.25,
                                  0.30, 0.35, 0.40, 0.45, 0.50};

const double plotDRmax = 0.5;

// Events in a coarse class, summed over its fine slices. base selects the
// denominator: "h_vz" for all events, "h_vz_triggerOn" for triggered ones.
double classEvents(TFile *f, int ci, const char *base)
{
  double n = 0.;
  for(int si = sliceLo[ci]; si <= sliceHi[ci]; si++){
    TH1D *h = nullptr;
    f->GetObject(Form("%s_C%d", base, si), h);
    if(!h) return -1.;
    n += h->Integral();
  }
  return n;
}

// scanFile overrides inFile, so an alternative scan can be checked without
// editing the constant; outSuffix is appended to each output filename.
void plotMuonDR_inclusiveClosestJet_coarseCent(const char *scanFile = nullptr,
                                               const char *outSuffix = ""){
  gStyle->SetOptStat(0);
  // per-event rates are ~1e-4, so without this the y labels render as
  // 0.0002, 0.0004, ... which are wide enough to push the axis title off the
  // canvas; 3 digits forces a shared "x10^-3" header instead
  TGaxis::SetMaxDigits(3);
  const char *usePath = (scanFile && scanFile[0]) ? scanFile : inFile;
  printf("scan file: %s\n", usePath);
  gSystem->mkdir(outDir, kTRUE);

  TFile *f = TFile::Open(gSystem->ExpandPathName(usePath));
  if(!f || f->IsZombie()){ printf("ERROR: cannot open %s\n", usePath); return; }

  bool gated = (f->GetListOfKeys()->FindObject(Form("%s_C1", histBaseGated)) != nullptr);
  const char *histBase = gated ? histBaseGated : histBaseUngated;
  const char *vzBase   = gated ? "h_vz_triggerOn" : "h_vz";
  printf("histogram : %s\n", histBase);
  printf("normalised: per %s event (%s)\n",
         gated ? "TRIGGERED" : "inclusive", vzBase);
  if(!gated)
    printf("  NOTE: this scan predates %s (added 2026-09-08).\n"
           "        Falling back to the ungated histogram normalised by all events.\n"
           "        Rates from it are NOT comparable to the trigger-gated ptRel plot.\n",
           histBaseGated);

  for(int ci = 0; ci < NClass; ci++){

    // sum the fine slices making up this coarse class
    TH2D *hSum = nullptr;
    for(int si = sliceLo[ci]; si <= sliceHi[ci]; si++){
      TH2D *h = nullptr;
      f->GetObject(Form("%s_C%d", histBase, si), h);
      if(!h){ printf("WARNING: %s_C%d missing\n", histBase, si); continue; }
      if(!hSum){ hSum = (TH2D*) h->Clone(Form("hSum_%d", ci)); hSum->SetDirectory(nullptr); }
      else hSum->Add(h);
    }
    if(!hSum){ printf("WARNING: no histograms found for class %s, skipping\n", classLabel[ci]); continue; }

    double nEvt = classEvents(f, ci, vzBase);
    if(nEvt <= 0.){ printf("WARNING: no %s for class %s, skipping\n", vzBase, classLabel[ci]); continue; }

    TH1D *proj[NPt];
    int   nEntries[NPt];

    for(int p = 0; p < NPt; p++){
      int b1 = hSum->GetYaxis()->FindBin(ptLo[p] + 1e-6);
      int b2 = hSum->GetYaxis()->FindBin(ptHi[p] - 1e-6);
      TH1D *raw = hSum->ProjectionX(Form("proj_%d_%d", ci, p), b1, b2);
      nEntries[p] = (int) raw->GetEntries();

      TH1D *h = (TH1D*) raw->Rebin(NEdge_dR-1, Form("h_dR_%d_%d", ci, p), edge_dR);
      h->SetDirectory(nullptr);
      h->Scale(1./nEvt);      // per-event rate
      divideByBinwidth(h);    // ... per unit dR

      proj[p] = h;
    }

    printf("%-8s N_evt=%.0f", classLabel[ci], nEvt);
    for(int p = 0; p < NPt; p++)
      printf("  |  %.0f-%.0f GeV: %d entries, %.3e/evt",
             ptLo[p], ptHi[p], nEntries[p], proj[p]->Integral("width"));
    printf("\n");

    TCanvas *c = new TCanvas(Form("c_%d", ci), "", 700, 700);
    TPad *pad = new TPad("pad", "", 0, 0, 1, 1);
    pad->SetLeftMargin(0.17);
    pad->SetBottomMargin(0.12);
    pad->Draw(); pad->cd();

    double ymax = 0.;
    for(int p = 0; p < NPt; p++) if(proj[p]->GetMaximum() > ymax) ymax = proj[p]->GetMaximum();

    for(int p = 0; p < NPt; p++){
      int col = TColor::GetColor(ptHex[p]);
      proj[p]->SetLineColor(col);
      proj[p]->SetMarkerColor(col);
      proj[p]->SetMarkerStyle(ptMarker[p]);
      proj[p]->SetMarkerSize(1.0);
      proj[p]->SetLineWidth(2);
      proj[p]->SetTitle("");
      proj[p]->GetXaxis()->SetTitle("#it{#Delta}#it{R}(#it{#mu},jet)");
      proj[p]->GetYaxis()->SetTitle(gated ? "d#it{N}/d(#it{#Delta}#it{R}) per triggered event"
                                          : "d#it{N}/d(#it{#Delta}#it{R}) per event");
      proj[p]->GetYaxis()->SetTitleOffset(1.55);
      proj[p]->GetXaxis()->SetTitleSize(0.045);
      proj[p]->GetYaxis()->SetTitleSize(0.045);
      proj[p]->GetXaxis()->SetLabelSize(0.04);
      proj[p]->GetYaxis()->SetLabelSize(0.04);
      proj[p]->GetXaxis()->SetRangeUser(0., plotDRmax);
      proj[p]->SetMaximum(ymax * 1.6);
      proj[p]->SetMinimum(0.);
      proj[p]->Draw(p == 0 ? "ep" : "ep same");
    }

    TLegend *leg = new TLegend(0.40, 0.60, 0.88, 0.80);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.032);
    for(int p = 0; p < NPt; p++)
      leg->AddEntry(proj[p], Form("%.0f < p_{T}^{jet} < %.0f GeV  (n=%d)", ptLo[p], ptHi[p], nEntries[p]), "lp");
    leg->Draw();

    TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.032);
    // kept clear of the y-axis "x10^-3" exponent, which ROOT draws at the top
    // left of the frame once SetMaxDigits forces scientific notation
    la.DrawLatex(0.21, 0.86,  Form("PbPb SingleMuon (5.02 TeV), %s", classLabel[ci]));
    la.DrawLatex(0.21, 0.815, Form("p_{T}^{#mu} > 15 GeV, |#eta^{#mu}| < 2, |#eta^{jet}| < 1.6%s",
                                   gated ? ", mu12 triggered" : ""));

    TString out = TString(outDir) + Form("muonDR_inclusiveClosestJet_coarseCent_%s%s.pdf",
                    TString(classLabel[ci]).ReplaceAll("%","pct").ReplaceAll("-","to").Data(),
                    outSuffix);
    c->SaveAs(out);
    printf("Saved %s\n\n", out.Data());
  }
}
