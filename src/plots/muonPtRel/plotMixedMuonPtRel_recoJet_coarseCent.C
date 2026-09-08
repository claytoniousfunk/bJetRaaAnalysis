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
// NORMALISATION: per TRIGGERED event, using h_vz_triggerOn summed over the
// class's fine slices. This histogram is filled only when the muon trigger bit
// is set, and the fill of h_mixedMuonPtRel_recoJetPt is likewise gated on
// evtTriggerDecision, so this is the matching denominator.
//
// This deliberately DIFFERS from the dR companion plot, which divides by h_vz
// (all events) because its own fill is ungated. Do not read a rate off one and
// compare it numerically to the other: for 0-10%, 50-60 GeV the two denominators
// differ by a factor ~8 (event trigger fraction 0.12), on top of this plot's
// extra dR < 0.4 requirement.
//
// A caveat on the denominator: h_vz_triggerOn is filled on the raw bit
// HLT_HIL3Mu12_v1 == 1, while evtTriggerDecision additionally requires a
// non-zero prescale via triggerIsOn(). The triggered-event count is therefore a
// slight overestimate of the true denominator, and this rate a slight
// underestimate. Events with the bit set and prescale zero should be rare
// (a disabled trigger should not fire), but the two are not identical by
// construction.
//
// Every entry is filled at w_resample = w/N_resamples, so one real event
// contributes total weight w regardless of the resample count and the rate is
// comparable across scans with different N.
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
// AXIS STOPS AT 5 GeV, matching the decomposition plot, but this distribution
// does not: 55-58% of entries in 0-10% and 10-30% lie ABOVE 5 GeV, 35-45% in
// 30-50%, and only 6-9% in 50-80%. That is the whole point of the template --
// most of its weight is outside the region the b-purity fit uses -- so the
// fraction above the axis maximum is printed and drawn on every panel rather
// than left for the reader to infer from a curve that appears to stop.
//
// The quoted <ptRel> is likewise computed over the FULL 0-10 axis, not the
// plotted range, so it can and does exceed the axis maximum.
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

// Variable binning, identical to plotPtRelTemplateDecomposition.C: 0.2 GeV
// through the peak region, widening to 0.8 GeV by 5 GeV. Native axis is 0.1
// wide, so every edge is a legal merge -- check that before editing.
const double edge_ptRel[] = {0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,
                             2.0,2.4,2.8,3.4,4.2,5.0};
const int NEdge_ptRel = (int)(sizeof(edge_ptRel)/sizeof(double));

const double plotPtRelMax = 5.0;

// Events in a coarse class, summed over its fine slices. base picks the
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

    double nEvt    = classEvents(f, ci, "h_vz_triggerOn");   // matching denominator
    double nEvtAll = classEvents(f, ci, "h_vz");              // for the printout only
    if(nEvt <= 0.){ printf("WARNING: no h_vz_triggerOn for class %s, skipping\n", classLabel[ci]); continue; }

    TH1D  *proj[NPt];
    int    nEntries[NPt];
    double meanPtRel[NPt];
    double fracAbove[NPt];   // fraction of entries beyond plotPtRelMax

    for(int p = 0; p < NPt; p++){
      int b1 = hSum->GetYaxis()->FindBin(ptLo[p] + 1e-6);
      int b2 = hSum->GetYaxis()->FindBin(ptHi[p] - 1e-6);
      TH1D *raw = hSum->ProjectionX(Form("proj_%d_%d", ci, p), b1, b2);
      nEntries[p]  = (int) raw->GetEntries();
      // mean taken BEFORE rebinning, so it is the true mean rather than one
      // computed from bin centres of wide bins
      meanPtRel[p] = raw->GetMean();
      { double tot = raw->Integral(0, raw->GetNbinsX()+1);
        double hi  = raw->Integral(raw->FindBin(plotPtRelMax + 1e-6), raw->GetNbinsX()+1);
        fracAbove[p] = (tot > 0.) ? hi/tot : 0.; }

      TH1D *h = (TH1D*) raw->Rebin(NEdge_ptRel-1, Form("h_ptRel_%d_%d", ci, p), edge_ptRel);
      h->SetDirectory(nullptr);
      h->Scale(1./nEvt);     // per-event rate
      divideByBinwidth(h);   // ... per unit ptRel
      proj[p] = h;
      delete raw;
    }

    printf("%-8s N_trig=%.0f (of %.0f, frac %.3f)", classLabel[ci], nEvt, nEvtAll, nEvt/nEvtAll);
    for(int p = 0; p < NPt; p++)
      printf("  |  %.0f-%.0f GeV: %d entries, %.3e/trig-evt, <ptRel>=%.2f, %.0f%% above %.0f",
             ptLo[p], ptHi[p], nEntries[p], proj[p]->Integral("width"), meanPtRel[p],
             100.*fracAbove[p], plotPtRelMax);
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
      proj[p]->GetYaxis()->SetTitle("d#it{N}/d#it{p}_{T}^{rel} per triggered event");
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
    la.DrawLatex(0.21, 0.815, "mixedEvent #mu + sameEvent reco jet, #Delta#it{R} < 0.4, mu12 triggered");
    la.DrawLatex(0.21, 0.77,  "p_{T}^{#mu} > 15 GeV, |#eta^{#mu}| < 2, |#eta^{jet}| < 1.6");
    { double fLo = 1., fHi = 0.;
      for(int p = 0; p < NPt; p++){ if(fracAbove[p] < fLo) fLo = fracAbove[p];
                                    if(fracAbove[p] > fHi) fHi = fracAbove[p]; }
      TLatex lo; lo.SetNDC(); lo.SetTextFont(42); lo.SetTextSize(0.030);
      lo.SetTextColor(kGray+3);
      lo.DrawLatex(0.21, 0.725, Form("%.0f#minus%.0f%% of entries lie above %.0f GeV",
                                     100.*fLo, 100.*fHi, plotPtRelMax)); }

    TString out = TString(outDir) + Form("mixedMuonPtRel_recoJet_coarseCent_%s%s.pdf",
                    TString(classLabel[ci]).ReplaceAll("%","pct").ReplaceAll("-","to").Data(),
                    outSuffix);
    c->SaveAs(out);
    printf("Saved %s\n\n", out.Data());
  }
}
