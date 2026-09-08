// ptRel of the (mixed-event muon, same-event reco jet) pairing --
// h_mixedMuonPtRel_recoJetPt_C* -- projected in three jet pT windows and merged
// into the standard 4 coarse centrality classes. One canvas per class.
//
// This is the ptRel companion to plotMuonDR_inclusiveClosestJet_coarseCent.C
// and uses the same classes, the same jet pT windows and the same normalisation
// so the two can be read side by side. It is NOT the same event sample though:
// the dR histogram records the closest reco jet at ANY distance, while this one
// is filled only when that jet is within dR < epsilon_mm = 0.4, i.e. only when
// the mixed muon would actually have been tagged. So this covers the dR < 0.4
// part of the dR plot.
//
// NORMALISATION: per event. Sum h_vz over the class's fine slices, divide, then
// divide by bin width, giving dN/d(ptRel) per event -- the same convention as
// plotPtRelTemplateDecomposition.C and the dR plot. Every entry is filled at
// w_resample = w/N_resamples, so one real event contributes total weight w
// regardless of the resample count and the rate is comparable across scans with
// different N.
//
// WHAT THIS IS FOR: this pairing is the (fake mu, real jet) background template
// subtracted from the data ptRel. Two things are worth reading off it.
//
//   1. Where it sits. In 0-10% and 10-30% the mean is ~5.0 GeV, far above the
//      data's ~1.6, so almost none of this template's weight lands in the
//      region the b-purity fit uses. That is the geometry problem documented in
//      plotPtRelTemplateDecomposition.C: the mixed muon is paired against a jet
//      axis computed WITHOUT it, so the opening angle is whatever random
//      geometry gives and, since solid angle grows as dR*d(dR), the pairs pile
//      up near the edge of the cone -- ptRel ~ p_mu * 0.4, about 6 GeV for a
//      15 GeV muon.
//
//   2. That it moves with centrality. The mean drops to ~4.1 in 30-50% and
//      ~1.5 in 50-80%, where it nearly coincides with the data. The template is
//      therefore not a fixed shape being scaled up and down; it changes
//      character across centrality, so a subtraction tuned in one class cannot
//      be assumed to behave in another.
//
// Full 0-10 GeV range is shown deliberately rather than the 0-3 fit window:
// clipping at 3 would cut away most of the distribution and make the template
// look small rather than misplaced.
//
// Usage: root -l -b -q 'plotMixedMuonPtRel_recoJet_coarseCent.C'
// Run from: src/plots/muonPtRel/

#include "../../../headers/functions/divideByBinwidth.h"

const char *inFile =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/"
  "PbPb_SingleMuon_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_"
  "mixedEventPFClustering_fastJetResamples-100_2026-9-8_ultraFineCentBins.root";

const char *histBase = "h_mixedMuonPtRel_recoJetPt";

const char *outDir = "../../../figures/muonPtRel/";

// coarse class -> [first, last] ultra-fine slice index, standard convention
const int NClass = 4;
const int sliceLo[NClass] = { 1,  3,  7, 11};
const int sliceHi[NClass] = { 2,  6, 10, 16};
const char *classLabel[NClass] = { "0-10%", "10-30%", "30-50%", "50-80%" };

// same windows as the dR plot, so the two figures line up
const int    NPt = 3;
const double ptLo[NPt] = {50, 60,  80};
const double ptHi[NPt] = {60, 80, 120};
const char  *ptHex[NPt]    = {"#0072B2", "#D55E00", "#009E73"};   // Okabe-Ito
const int    ptMarker[NPt] = {20, 21, 33};    // circle, square, diamond -- no triangles

// Variable binning: 0.2 GeV through the low-ptRel region where the data lives
// and where the peripheral classes pile up, widening out to 1.0 GeV in the tail
// that only the central classes populate. Native axis is 0.1 wide, so every
// edge below is a legal merge -- check that before editing. Sparsest window
// (30-50%, 80-120 GeV) has 249 entries, about 11 per bin at this granularity.
const double edge_ptRel[] = {0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,
                             2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,
                             7.0,8.0,9.0,10.0};
const int NEdge_ptRel = (int)(sizeof(edge_ptRel)/sizeof(double));

const double plotPtRelMax = 10.0;

// events in a coarse class = sum of h_vz over its fine slices, matching
// classEvents() in plotPtRelTemplateDecomposition.C
double classEvents(TFile *f, int ci)
{
  double n = 0.;
  for(int si = sliceLo[ci]; si <= sliceHi[ci]; si++){
    TH1D *h = nullptr;
    f->GetObject(Form("h_vz_C%d", si), h);
    if(!h) return -1.;
    n += h->Integral();
  }
  return n;
}

// scanFile overrides inFile; outSuffix is appended to each output filename.
void plotMixedMuonPtRel_recoJet_coarseCent(const char *scanFile = nullptr,
                                           const char *outSuffix = "")
{
  gStyle->SetOptStat(0);
  // per-event rates are ~1e-6, so without this the y labels are wide enough to
  // push the axis title off the canvas
  TGaxis::SetMaxDigits(3);

  const char *usePath = (scanFile && scanFile[0]) ? scanFile : inFile;
  printf("scan file: %s\n", usePath);
  gSystem->mkdir(outDir, kTRUE);

  TFile *f = TFile::Open(gSystem->ExpandPathName(usePath));
  if(!f || f->IsZombie()){ printf("ERROR: cannot open %s\n", usePath); return; }

  for(int ci = 0; ci < NClass; ci++){

    TH2D *hSum = nullptr;
    for(int si = sliceLo[ci]; si <= sliceHi[ci]; si++){
      TH2D *h = nullptr;
      f->GetObject(Form("%s_C%d", histBase, si), h);
      if(!h){ printf("WARNING: %s_C%d missing\n", histBase, si); continue; }
      if(!hSum){ hSum = (TH2D*) h->Clone(Form("hSum_%d", ci)); hSum->SetDirectory(nullptr); }
      else hSum->Add(h);
    }
    if(!hSum){ printf("WARNING: no histograms for class %s, skipping\n", classLabel[ci]); continue; }

    double nEvt = classEvents(f, ci);
    if(nEvt <= 0.){ printf("WARNING: no h_vz for class %s, skipping\n", classLabel[ci]); continue; }

    TH1D  *proj[NPt];
    int    nEntries[NPt];
    double meanPtRel[NPt];

    for(int p = 0; p < NPt; p++){
      int b1 = hSum->GetYaxis()->FindBin(ptLo[p] + 1e-6);
      int b2 = hSum->GetYaxis()->FindBin(ptHi[p] - 1e-6);
      TH1D *raw = hSum->ProjectionX(Form("proj_%d_%d", ci, p), b1, b2);
      nEntries[p]  = (int) raw->GetEntries();
      // mean taken BEFORE rebinning, so it is the true mean rather than one
      // computed from bin centres of wide bins
      meanPtRel[p] = raw->GetMean();

      TH1D *h = (TH1D*) raw->Rebin(NEdge_ptRel-1, Form("h_ptRel_%d_%d", ci, p), edge_ptRel);
      h->SetDirectory(nullptr);
      h->Scale(1./nEvt);     // per-event rate
      divideByBinwidth(h);   // ... per unit ptRel
      proj[p] = h;
      delete raw;
    }

    printf("%-8s N_evt=%.0f", classLabel[ci], nEvt);
    for(int p = 0; p < NPt; p++)
      printf("  |  %.0f-%.0f GeV: %d entries, %.3e/evt, <ptRel>=%.2f",
             ptLo[p], ptHi[p], nEntries[p], proj[p]->Integral("width"), meanPtRel[p]);
    printf("\n");

    TCanvas *c = new TCanvas(Form("c_%d", ci), "", 700, 700);
    TPad *pad = new TPad("pad", "", 0, 0, 1, 1);
    pad->SetLeftMargin(0.17);
    pad->SetBottomMargin(0.12);
    pad->Draw(); pad->cd();

    double ymax = 0.;
    for(int p = 0; p < NPt; p++) if(proj[p]->GetMaximum() > ymax) ymax = proj[p]->GetMaximum();
    if(ymax <= 0.){ printf("  class %s: empty, skipping plot\n", classLabel[ci]); delete c; continue; }

    for(int p = 0; p < NPt; p++){
      int col = TColor::GetColor(ptHex[p]);
      proj[p]->SetLineColor(col);
      proj[p]->SetMarkerColor(col);
      proj[p]->SetMarkerStyle(ptMarker[p]);
      proj[p]->SetMarkerSize(1.0);
      proj[p]->SetLineWidth(2);
      proj[p]->SetTitle("");
      proj[p]->GetXaxis()->SetTitle("#it{p}_{T}^{rel}(#it{#mu},jet) [GeV]");
      proj[p]->GetYaxis()->SetTitle("d#it{N}/d#it{p}_{T}^{rel} per event");
      proj[p]->GetXaxis()->SetTitleSize(0.045);
      proj[p]->GetYaxis()->SetTitleSize(0.045);
      proj[p]->GetYaxis()->SetTitleOffset(1.55);
      proj[p]->GetXaxis()->SetLabelSize(0.04);
      proj[p]->GetYaxis()->SetLabelSize(0.04);
      proj[p]->GetXaxis()->SetRangeUser(0., plotPtRelMax);
      proj[p]->SetMaximum(ymax * 1.6);
      proj[p]->SetMinimum(0.);
      proj[p]->Draw(p == 0 ? "ep" : "ep same");
    }

    // wide and low: each entry carries the pT window, the entry count and the
    // mean, which does not fit in the right-hand third of the frame, and the
    // three header lines above run down to 0.77
    TLegend *leg = new TLegend(0.30, 0.50, 0.94, 0.72);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.026);
    for(int p = 0; p < NPt; p++)
      leg->AddEntry(proj[p], Form("%.0f < p_{T}^{jet} < %.0f GeV  (n=%d, #LT#it{p}_{T}^{rel}#GT=%.2f)",
                                  ptLo[p], ptHi[p], nEntries[p], meanPtRel[p]), "lp");
    leg->Draw();

    TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.032);
    la.DrawLatex(0.21, 0.86,  Form("PbPb SingleMuon (5.02 TeV), %s", classLabel[ci]));
    la.DrawLatex(0.21, 0.815, "mixedEvent #mu + sameEvent reco jet, #Delta#it{R} < 0.4");
    la.DrawLatex(0.21, 0.77,  "p_{T}^{#mu} > 15 GeV, |#eta^{#mu}| < 2, |#eta^{jet}| < 1.6");

    TString out = TString(outDir) + Form("mixedMuonPtRel_recoJet_coarseCent_%s%s.pdf",
                    TString(classLabel[ci]).ReplaceAll("%","pct").ReplaceAll("-","to").Data(),
                    outSuffix);
    c->SaveAs(out);
    printf("Saved %s\n\n", out.Data());
  }
}
