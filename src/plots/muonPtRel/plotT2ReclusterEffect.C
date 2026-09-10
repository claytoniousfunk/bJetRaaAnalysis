// What reclustering did to the T2 template.
//
// T2 is the (fake mu, real jet) background: a mixed-event muon paired with a
// jet from a real event. Two ways of building it:
//
//   old  h_mixedMuonPtRel_recoJetPt        the muon is matched by dR < 0.4 to an
//                                          akCs4PF reco jet AFTER clustering, so
//                                          the jet pT does NOT contain the muon
//                                          and the axis was built without it
//   new  h_injMuonPtRel_donorJetPt_inject  the muon is ADDED to the next
//                                          centrality-matched event's PF
//                                          candidates and everything is
//                                          reclustered, so the jet pT contains
//                                          the muon and the axis is pulled onto
//                                          it -- the geometry the data has
//
// Both are drawn SHAPE-ONLY, unit-normalised over the plotted range, because the
// question here is how the distribution moved, not how the yield changed. The
// per-event rates are printed to the terminal for the yield side of it.
//
// TWO OBSERVABLES:
//   ptRel  the fit variable. Plotted over the FULL 0-10 GeV axis, not the 0-5
//          window the decomposition uses -- the old template's weight sits above
//          5 GeV, so cropping there would hide most of the change.
//   dR     muon to jet axis. Once the muon is clustered in it dominates the
//          jet and the axis is pulled onto it, so dR should collapse.
//          Both dR histograms are "nearest jet, no dR cut", so both carry a
//          large tail beyond the 0-0.5 axis (muons with no jet near them);
//          the shapes here are of the in-range part only.
//
// NORMALISATION for the printed rates: the old template is trigger-gated and so
// divided by h_vz_triggerOn; the injection fills are ungated and divided by
// h_vz. The shapes are unaffected by that choice.
//
// Usage: root -l -b -q 'plotT2ReclusterEffect.C'
// Run from: src/plots/muonPtRel/

#include "../../../headers/plotting/coarseCent.h"
#include "../../../headers/plotting/plotStyle.h"
#include "../../../headers/plotting/ratioPanel.h"
#include "../../../headers/functions/divideByBinwidth.h"

const char *scanPath =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/"
  "PbPb_SingleMuon_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_"
  "mixedEventPFClustering_fastJetResamples-100_2026-9-8_ultraFineCentBins.root";

const char *nameOld = "h_mixedMuonPtRel_recoJetPt";
const char *nameNew = "h_injMuonPtRel_donorJetPt_inject";
// dR counterparts of the same two constructions
const char *nameOldDR = "h_muonDR_inclusiveClosestJet";
const char *nameNewDR = "h_injMuonDR_donorJetPt_inject";

const char *outDir = "../../../figures/ptRelDecomposition/";

// full ptRel axis: the old template peaks near 5, so 0-5 would crop the change
const double edgePtRel[] = {0.0,0.4,0.8,1.2,1.6,2.0,2.5,3.0,3.5,4.0,
                            4.5,5.0,5.5,6.0,7.0,8.0,10.0};
const int    nEdgePtRel  = (int)(sizeof(edgePtRel)/sizeof(double)) - 1;

const double edgeDR[] = {0,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50};
const int    nEdgeDR  = (int)(sizeof(edgeDR)/sizeof(double)) - 1;

// one canvas: two shapes overlaid, unit-normalised, with a new/old ratio panel
// mOld/mNew are passed in rather than taken from the drawn histograms: those
// have been rebinned to wide variable bins, so GetMean() on them is computed
// from bin centres and is wrong by ~10%.
void drawShape(TH1D *hOld, TH1D *hNew, double mOld, double mNew,
               const char *xTitle,
               const char *headLine, double xMin, double xMax,
               const char *outName)
{
  TCanvas *c = new TCanvas(Form("c_%s", outName), "", 700, 800);
  TPad *pT = nullptr, *pB = nullptr;
  splitPads(pT, pB);

  pT->cd();
  styleH(hOld, hexBkg1, markOpenCircle);
  styleH(hNew, hexCorrected, markFilledSquare);

  double ymax = TMath::Max(hOld->GetMaximum(), hNew->GetMaximum());
  if(ymax <= 0.){ printf("  %s: empty, skipped\n", outName); delete c; return; }

  hOld->GetXaxis()->SetRangeUser(xMin, xMax);
  hOld->GetXaxis()->SetLabelSize(0);
  hOld->GetYaxis()->SetTitle("normalised to unit area");
  hOld->GetYaxis()->SetTitleSize(0.050); hOld->GetYaxis()->SetTitleOffset(1.55);
  hOld->GetYaxis()->SetLabelSize(0.042);
  hOld->SetTitle(""); hOld->SetMinimum(0.); hOld->SetMaximum(ymax*1.65);
  hOld->Draw("E");
  hNew->Draw("E same");

  TLegend *leg = makeLegend(0.40, 0.62, 0.96, 0.83, 0.036);
  leg->AddEntry(hOld, Form("matched post-cluster  (#mu=%.2f)", mOld), "lp");
  leg->AddEntry(hNew, Form("reclustered  (#mu=%.2f)",          mNew), "lp");
  leg->Draw();

  TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.044);
  la.DrawLatex(0.21, 0.885, headLine);
  // The mean shift used to be spelled out here in bold as
  // "<x>: a -> b (+c%)". Removed 2026-09-10 as too distracting: it competed
  // with the histograms for attention and repeated what the legend's per-curve
  // means already carry. The shift is still printed to the terminal.

  pB->cd();
  TH1D *r = makeRatio(hNew, hOld, Form("r_%s", outName), RatioErr::kBoth);
  styleH(r, hexCorrected, markFilledSquare);
  styleRatioAxes(r, xTitle, "reclust. / old");
  r->GetXaxis()->SetRangeUser(xMin, xMax);
  r->SetMinimum(0.);
  r->SetMaximum(ratioMax(r, hOld, xMin, xMax) * 1.25);
  r->Draw("E");
  TLine *l = unityLine(xMin, xMax); l->Draw();
  r->Draw("E same");

  c->SaveAs(Form("%s%s.pdf", outDir, outName));
  delete c;
}

void plotT2ReclusterEffect()
{
  initPlotStyle();
  gSystem->mkdir(outDir, kTRUE);

  TFile *f = TFile::Open(scanPath);
  if(!f || f->IsZombie()){ printf("ERROR: cannot open %s\n", scanPath); return; }

  printf("T2 template, post-cluster matched vs reclustered\n");
  printf("%-9s %10s %10s %9s | %12s %12s %8s\n",
         "class", "<x> old", "<x> new", "shift",
         "rate old", "rate new", "ratio");

  for(int ci = 0; ci < NCoarse; ci++){

    TH2D *hOld2 = coarseSum2(f, nameOld, ci, "old");
    TH2D *hNew2 = coarseSum2(f, nameNew, ci, "new");
    if(!hOld2 || !hNew2){ printf("  %s: missing, skipped\n", coarseLabel[ci]); continue; }

    double nTrig = coarseEvents(f, ci, "h_vz_triggerOn");
    double nAll  = coarseEvents(f, ci, "h_vz");
    if(nTrig <= 0. || nAll <= 0.){ printf("  %s: no event counts\n", coarseLabel[ci]); continue; }

    // ---- ptRel, integrated over all jet pT (max statistics for a shape) ----
    TH1D *pOldRaw = hOld2->ProjectionX(Form("poR_%d", ci));  pOldRaw->SetDirectory(nullptr);
    TH1D *pNewRaw = hNew2->ProjectionX(Form("pnR_%d", ci));  pNewRaw->SetDirectory(nullptr);

    // means BEFORE rebinning, so they are not computed from wide-bin centres
    double mOld = pOldRaw->GetMean(), mNew = pNewRaw->GetMean();
    // per-event rates, each with its matching denominator
    double rOld = pOldRaw->Integral(0, pOldRaw->GetNbinsX()+1) / nTrig;
    double rNew = pNewRaw->Integral(0, pNewRaw->GetNbinsX()+1) / nAll;

    printf("%-9s %10.3f %10.3f %8.0f%% | %12.3e %12.3e %8.2f   (ptRel)\n",
           coarseLabel[ci], mOld, mNew,
           mOld > 0. ? 100.*(mNew/mOld - 1.) : 0., rOld, rNew,
           rOld > 0. ? rNew/rOld : 0.);

    TH1D *pOld = rebinTo(pOldRaw, nEdgePtRel, edgePtRel, Form("po_%d", ci));
    TH1D *pNew = rebinTo(pNewRaw, nEdgePtRel, edgePtRel, Form("pn_%d", ci));
    delete pOldRaw; delete pNewRaw;
    // unit area, then per unit x so variable widths do not distort the shape
    if(pOld->Integral() > 0.) pOld->Scale(1./pOld->Integral());
    if(pNew->Integral() > 0.) pNew->Scale(1./pNew->Integral());
    divideByBinwidth(pOld); divideByBinwidth(pNew);

    drawShape(pOld, pNew, mOld, mNew, "#it{p}_{T}^{rel} [GeV]",
              Form("PbPb %s, T2 template", coarseLabel[ci]),
              0., 10., Form("t2Recluster_ptRel_%s", coarseTag[ci]));

    // ---- dR, from the matching pair of nearest-jet histograms ----
    TH2D *dOld2 = coarseSum2(f, nameOldDR, ci, "odr");
    TH2D *dNew2 = coarseSum2(f, nameNewDR, ci, "ndr");
    if(dOld2 && dNew2){
      TH1D *dOldRaw = dOld2->ProjectionX(Form("doR_%d", ci)); dOldRaw->SetDirectory(nullptr);
      TH1D *dNewRaw = dNew2->ProjectionX(Form("dnR_%d", ci)); dNewRaw->SetDirectory(nullptr);
      // in-range means: both histograms park "no jet anywhere near" in the
      // overflow, so a mean including it would be meaningless
      double dmOld = dOldRaw->GetMean(), dmNew = dNewRaw->GetMean();

      TH1D *dOld = rebinTo(dOldRaw, nEdgeDR, edgeDR, Form("do_%d", ci));
      TH1D *dNew = rebinTo(dNewRaw, nEdgeDR, edgeDR, Form("dn_%d", ci));
      delete dOldRaw; delete dNewRaw;
      if(dOld->Integral() > 0.) dOld->Scale(1./dOld->Integral());
      if(dNew->Integral() > 0.) dNew->Scale(1./dNew->Integral());
      divideByBinwidth(dOld); divideByBinwidth(dNew);

      printf("%-9s %10.3f %10.3f %8.0f%%   (dR)\n", coarseLabel[ci], dmOld, dmNew,
             dmOld > 0. ? 100.*(dmNew/dmOld - 1.) : 0.);

      drawShape(dOld, dNew, dmOld, dmNew, "#it{#Delta}#it{R}(#mu,jet)",
                Form("PbPb %s, T2 template", coarseLabel[ci]),
                0., 0.5, Form("t2Recluster_dR_%s", coarseTag[ci]));
      delete dOld2; delete dNew2;
    }
    else printf("  %s: dR histograms missing\n", coarseLabel[ci]);

    delete hOld2; delete hNew2;
  }

  printf("\nfigures written to %s\n", outDir);
}
