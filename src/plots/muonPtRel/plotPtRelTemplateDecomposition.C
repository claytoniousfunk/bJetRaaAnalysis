// ptRel template decomposition: isolate the real-muon-real-jet contribution in
// data, then compare the corrected shape to MC.
//
// ---------------------------------------------------------------------------
// THE DECOMPOSITION
// ---------------------------------------------------------------------------
// What data measures is a (tight reco muon, reco jet) pair with dR < 0.4. Two
// independent things can spoil such a pair:
//
//   A = the jet is combinatorial (built out of underlying event, not a hard
//       parton)
//   B = the muon is uncorrelated with the jet (it did not come from that jet's
//       heavy-flavour decay; it just landed inside the cone)
//
// A pair is background if EITHER holds, so by inclusion-exclusion
//
//   Bkg = P(A) + P(B) - P(A and B)
//
// and the three mixed-event templates estimate exactly those three terms:
//
//   T2 = h_realMuonPtRel_mixedFastJetPt   real muon  + mixed-event fastJet  -> P(A)
//   T3 = h_mixedMuonPtRel_recoJetPt       mixed-event PF muon + reco jet    -> P(B)
//   T4 = h_fastJetMuonPtRel_..._RC        mixed muon + mixed fastJet        -> P(A and B)
//
// giving
//
//   S = D - T2 - T3 + T4
//
// The +T4 is NOT optional. T2's real muons sit on top of a fully combinatorial
// jet, and T3's uncorrelated muons sit on top of a reco jet collection that is
// itself part combinatorial, so the doubly-fake configuration is counted once
// in each -- without adding it back the subtraction over-subtracts.
//
// Every term is converted to a PER-EVENT rate (divide by that class's h_vz
// integral) before combining, because that is the only normalisation in which
// the four terms are commensurate. All four also come out of the SAME scan of
// the SAME file, so the event counts cancel exactly rather than approximately;
// do not swap one of them for a histogram from a different scan without
// checking that the event selection is identical.
//
// ---------------------------------------------------------------------------
// CAVEATS -- read before trusting a number out of this
// ---------------------------------------------------------------------------
// 1. T2 and T4 use FastJet anti-kT R=0.4 jets reclustered from mixed-event PF
//    candidates, RC-background-subtracted. Data and T3 use akCs4PF reco jets.
//    These are different jet definitions. The mixed-event fastJet is a PROXY
//    for the fake component of the reco collection, not the same object, and
//    the two need not have the same pT scale or the same ptRel response. This
//    is the single largest systematic in the whole construction and it is not
//    quantified here.
//
// 2. T3's denominator is "all reco jets", i.e. genuine + fake. That is exactly
//    what the inclusion-exclusion above assumes, so it is consistent -- but it
//    also means T3 alone is NOT "fake muons on real jets", despite the name.
//
// 3. S can go negative bin-by-bin where the subtraction overshoots. Negative
//    bins are left in place rather than clipped: clipping would bias the
//    integral upward and hide exactly the failure mode worth seeing.
//
// 4. The corrected data / MC comparison is a SHAPE comparison. Both sides are
//    unit-normalised over the fit range because the MC normalisation carries
//    pThat cross-section weights that have no counterpart in the per-event data
//    rate. Do not read the ratio's overall level as a yield statement.
//
// ---------------------------------------------------------------------------
// Usage: root -l -b -q 'plotPtRelTemplateDecomposition.C'
// Run from: src/plots/muonPtRel/
//
// Runs today with T2/T3 absent (they were added to PbPb_pfCandAnalyzer.C after
// the input scan was produced): missing templates are reported and dropped from
// the subtraction, and every plot is labelled with which terms were actually
// applied so a partial result cannot be mistaken for a complete one.
// ---------------------------------------------------------------------------

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TStyle.h"
#include "TString.h"
#include <cstdio>
#include <vector>

// ---------------------------------------------------------------- inputs ----

// data D, fake-fake T4, and (after the rescan) T2 and T3 all live here
const char *scanPath =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/"
  "PbPb_SingleMuon_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_"
  "mixedEventPFClustering_pseudoJetCandPtMin-0.0_2026-9-1_ultraFineCentBins.root";

// MC reference: PYTHIA+HYDJET, all jets, same muon tag as data
const char *mcPath =
  "/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/scanningOutput/"
  "PYTHIAHYDJET/latest/PYTHIAHYDJET_DiJet_pThat-30_mu12_pTmu-15_tight_"
  "mu12TriggerEfficiencyCorrection_vzReweight_hiBinReweight_hiBinShift-10_"
  "leadingXjetDumpFilter_jetTrkMaxFilter_removeHYDJETjet0p35CutOnGen_WDecayFilter_"
  "weightCut_2026-2-12.root";

const char *nameD  = "h_muptrel_recoJetPt_inclRecoMuonTag_triggerOn";  // real mu + real jet
const char *nameT2 = "h_realMuonPtRel_mixedFastJetPt";                 // real mu + fake jet
const char *nameT3 = "h_mixedMuonPtRel_recoJetPt";                     // fake mu + real jet
// Two candidates for the doubly-fake term. They are the SAME muon on the SAME
// jet; the only difference is the jet pT scale. Prefer the JEC one, because D
// and T3 are both on the JEC scale and a mixed-scale subtraction compares
// different jets in the same nominal pT window. The RC one is the pre-existing
// histogram on the raw rcSub scale, kept only as a fallback for scans produced
// before the JEC twin was added.
const char *nameT4     = "h_mixedMuonPtRel_mixedFastJetPt";            // fake mu + fake jet, JEC scale
const char *nameT4_raw = "h_fastJetMuonPtRel_fastJetPt_PF_bkgSub_RC";  // same, raw rcSub scale
const char *nameMC = "h_muptrel_recoJetPt_inclRecoMuonTag_triggerOn_allJets"; // + "_C%dT0"

const char *outDir = "/home/clayton/Analysis/code/bJetRaaAnalysis/figures/ptRelDecomposition";

// ------------------------------------------------------- binning / classes ----

// coarse class -> [first,last] ultra-fine 5% slice, same convention as
// makeFakeJetFile.C and makeCoarseBins_muonPtRel.C
const int NClass = 4;
const int sliceLo[NClass]      = {1, 3, 7, 11};
const int sliceHi[NClass]      = {2, 6, 10, 16};
const char *classLabel[NClass] = {"0-10%", "10-30%", "30-50%", "50-80%"};

// Jet pT windows. Deliberately reaching down to 40 GeV rather than starting at
// the template fit's 80 GeV: the one background term available before the
// rescan (T4, doubly-fake) is 7.7% of data at 40-60 GeV in 0-10% but only 0.7%
// at 80-100 and 0.2% at 100-120. The whole correction is a low-pT effect, and
// low pT is also where the 10-30/0-10 degeneracy lives, so a decomposition
// evaluated only in 80-150 would be measuring the region where it cannot
// matter. The upper windows are kept because that is where the b-purity fit
// actually runs.
const int NPt = 5;
const double ptLo[NPt] = { 40., 60., 80., 100., 120.};
const double ptHi[NPt] = { 60., 80., 100., 120., 150.};

// ptRel range used for the shape comparison and the unit normalisation --
// matches the 0-3 fit range settled on for the template fit
const double ptRelFitLo = 0.0;
const double ptRelFitHi = 3.0;

// Okabe-Ito
const int colD  = TColor::GetColor("#000000");
const int colT2 = TColor::GetColor("#E69F00");
const int colT3 = TColor::GetColor("#56B4E9");
const int colT4 = TColor::GetColor("#009E73");
const int colS  = TColor::GetColor("#D55E00");
const int colMC = TColor::GetColor("#0072B2");

// ------------------------------------------------------------- helpers ----

// Sum the ultra-fine slices of a TH2 over a coarse class and project the ptRel
// (X) axis in a jet-pT (Y) window. Returns nullptr if ANY slice is missing --
// a partial class sum would be silently mis-normalised against a complete h_vz.
TH1D* projectClass(TFile *f, const char *base, int ci, int pi, const char *tag)
{
  TH2D *sum = nullptr;
  for(int si = sliceLo[ci]; si <= sliceHi[ci]; si++){
    TH2D *h = nullptr;
    f->GetObject(Form("%s_C%d", base, si), h);
    if(!h){ if(sum) delete sum; return nullptr; }
    if(!sum){ sum = (TH2D*) h->Clone(Form("sum_%s_%s_c%d_p%d", base, tag, ci, pi)); sum->SetDirectory(nullptr); }
    else sum->Add(h);
  }
  if(!sum) return nullptr;

  // FindBin on the low edge is correct; on the high edge it returns the bin
  // that STARTS at ptHi, which must be excluded -> -1.
  int by1 = sum->GetYaxis()->FindBin(ptLo[pi] + 1e-6);
  int by2 = sum->GetYaxis()->FindBin(ptHi[pi] - 1e-6);
  TH1D *p = sum->ProjectionX(Form("px_%s_%s_c%d_p%d", base, tag, ci, pi), by1, by2);
  p->SetDirectory(nullptr);
  delete sum;
  return p;
}

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

void styleH(TH1 *h, int col, int mark)
{
  h->SetLineColor(col); h->SetMarkerColor(col);
  h->SetMarkerStyle(mark); h->SetMarkerSize(0.9); h->SetLineWidth(2);
  h->SetStats(0);
}

// integral over the fit range only
double fitRangeIntegral(TH1 *h)
{
  int b1 = h->FindBin(ptRelFitLo + 1e-6);
  int b2 = h->FindBin(ptRelFitHi - 1e-6);
  return h->Integral(b1, b2);
}

// ------------------------------------------------------------- driver ----

void plotPtRelTemplateDecomposition()
{
  gStyle->SetOptStat(0);
  gSystem->mkdir(outDir, kTRUE);

  TFile *fS = TFile::Open(scanPath);
  if(!fS || fS->IsZombie()){ printf("ERROR: cannot open scan file\n  %s\n", scanPath); return; }
  TFile *fM = TFile::Open(mcPath);
  if(!fM || fM->IsZombie()){ printf("ERROR: cannot open MC file\n  %s\n", mcPath); return; }

  // ---- report which templates the input actually carries -------------------
  bool haveT2 = false, haveT3 = false, haveT4 = false, t4IsRaw = false;
  const char *nameT4_jec = nameT4;   // nameT4 is repointed on fallback; keep the JEC name for the message
  { TH2D *h=nullptr;
    fS->GetObject(Form("%s_C1",nameT2),h); haveT2 = (h!=nullptr); h=nullptr;
    fS->GetObject(Form("%s_C1",nameT3),h); haveT3 = (h!=nullptr); h=nullptr;
    fS->GetObject(Form("%s_C1",nameT4),h); haveT4 = (h!=nullptr); h=nullptr;
    if(!haveT4){
      fS->GetObject(Form("%s_C1",nameT4_raw),h);
      if(h){ haveT4 = true; t4IsRaw = true; nameT4 = nameT4_raw; } } }

  printf("\n=== template availability in the input scan ===\n");
  printf("  D  real mu + real jet   %-40s present\n", nameD);
  printf("  T2 real mu + fake jet   %-40s %s\n", nameT2, haveT2 ? "present" : "*** MISSING ***");
  printf("  T3 fake mu + real jet   %-40s %s\n", nameT3, haveT3 ? "present" : "*** MISSING ***");
  printf("  T4 fake mu + fake jet   %-40s %s\n", nameT4, haveT4 ? "present" : "*** MISSING ***");
  if(t4IsRaw)
    printf("\n  WARNING: falling back to the RAW rcSub-scale doubly-fake histogram --\n"
           "           its jet pT is not JEC-corrected while D and T3 are, so the\n"
           "           jet pT windows do not select the same jets on both sides of\n"
           "           the subtraction. Rerun the scan to pick up %s.\n", nameT4_jec);
  if(!haveT2 || !haveT3){
    printf("\n  NOTE: the missing templates were added to PbPb_pfCandAnalyzer.C AFTER\n");
    printf("        this scan was produced. The decomposition below is PARTIAL --\n");
    printf("        rerun the scan to get the full correction.\n");
  }
  // With T2 or T3 absent the inclusion-exclusion is incomplete, and adding T4
  // back would then be correcting a double-count that never happened.
  bool useT4 = haveT4 && haveT2 && haveT3;
  if(haveT4 && !useT4)
    printf("        T4 is present but is being SKIPPED: it only corrects a\n"
           "        double-count between T2 and T3, and both must be applied\n"
           "        for that double-count to exist.\n");

  TString applied = "D";
  if(haveT2) applied += " #minus T2";
  if(haveT3) applied += " #minus T3";
  if(useT4)  applied += " + T4";

  // ---- per-event rate table ------------------------------------------------
  printf("\n=== per-event tagged-muon rate, ptRel in [%.1f,%.1f] ===\n", ptRelFitLo, ptRelFitHi);
  printf("%-8s %-12s %10s %10s %10s %10s %10s %8s\n",
         "class","jet pT","N_evt","D","T2","T3","T4","S/D");

  for(int ci = 0; ci < NClass; ci++){

    double nEvt = classEvents(fS, ci);
    if(nEvt <= 0.){ printf("  class %s: no h_vz, skipping\n", classLabel[ci]); continue; }

    double nEvtMC = 0.;
    { TH1D *h=nullptr; fM->GetObject(Form("h_vz_C%d", ci+1), h); if(h) nEvtMC = h->Integral(); }

    for(int pi = 0; pi < NPt; pi++){

      TH1D *hD  = projectClass(fS, nameD,  ci, pi, "D");
      TH1D *hT2 = haveT2 ? projectClass(fS, nameT2, ci, pi, "T2") : nullptr;
      TH1D *hT3 = haveT3 ? projectClass(fS, nameT3, ci, pi, "T3") : nullptr;
      TH1D *hT4 = useT4  ? projectClass(fS, nameT4, ci, pi, "T4") : nullptr;
      if(!hD){ printf("  class %s pT[%.0f,%.0f]: no data histogram, skipping\n",
                      classLabel[ci], ptLo[pi], ptHi[pi]); continue; }

      // per-event rates
      hD->Scale(1./nEvt);
      if(hT2) hT2->Scale(1./nEvt);
      if(hT3) hT3->Scale(1./nEvt);
      if(hT4) hT4->Scale(1./nEvt);

      // S = D - T2 - T3 + T4
      TH1D *hS = (TH1D*) hD->Clone(Form("hS_c%d_p%d", ci, pi));
      hS->SetDirectory(nullptr);
      if(hT2) hS->Add(hT2, -1.);
      if(hT3) hS->Add(hT3, -1.);
      if(hT4) hS->Add(hT4, +1.);

      double iD  = fitRangeIntegral(hD);
      double iS  = fitRangeIntegral(hS);
      printf("%-8s %4.0f-%-7.0f %10.0f %10.2e %10.2e %10.2e %10.2e %8.3f\n",
             classLabel[ci], ptLo[pi], ptHi[pi], nEvt,
             iD,
             hT2 ? fitRangeIntegral(hT2) : 0.,
             hT3 ? fitRangeIntegral(hT3) : 0.,
             hT4 ? fitRangeIntegral(hT4) : 0.,
             iD > 0. ? iS/iD : 0.);

      // ---- MC reference ----------------------------------------------------
      TH1D *hMC = nullptr;
      { TH2D *h2 = nullptr;
        fM->GetObject(Form("%s_C%dT0", nameMC, ci+1), h2);
        if(h2){
          int by1 = h2->GetYaxis()->FindBin(ptLo[pi] + 1e-6);
          int by2 = h2->GetYaxis()->FindBin(ptHi[pi] - 1e-6);
          hMC = h2->ProjectionX(Form("pxMC_c%d_p%d", ci, pi), by1, by2);
          hMC->SetDirectory(nullptr);
          if(nEvtMC > 0.) hMC->Scale(1./nEvtMC);
        } }

      // ---- shape comparison: unit-normalise over the fit range -------------
      TH1D *hSn  = (TH1D*) hS->Clone(Form("hSn_c%d_p%d", ci, pi));  hSn->SetDirectory(nullptr);
      TH1D *hDn  = (TH1D*) hD->Clone(Form("hDn_c%d_p%d", ci, pi));  hDn->SetDirectory(nullptr);
      TH1D *hMCn = hMC ? (TH1D*) hMC->Clone(Form("hMCn_c%d_p%d", ci, pi)) : nullptr;
      if(hMCn) hMCn->SetDirectory(nullptr);
      { double n;
        n = fitRangeIntegral(hSn);  if(n != 0.) hSn->Scale(1./n);
        n = fitRangeIntegral(hDn);  if(n != 0.) hDn->Scale(1./n);
        if(hMCn){ n = fitRangeIntegral(hMCn); if(n != 0.) hMCn->Scale(1./n); } }

      // ---- draw ------------------------------------------------------------
      TCanvas *c = new TCanvas(Form("c_c%d_p%d", ci, pi), "", 700, 800);
      TPad *pTop = new TPad("pTop","",0,0.34,1,1);
      TPad *pBot = new TPad("pBot","",0,0,1,0.34);
      pTop->SetBottomMargin(0.02); pTop->SetLeftMargin(0.15); pTop->SetTopMargin(0.07);
      pBot->SetTopMargin(0.02);    pBot->SetLeftMargin(0.15); pBot->SetBottomMargin(0.32);
      pTop->Draw(); pBot->Draw();

      // top: the per-event decomposition
      pTop->cd(); pTop->SetLogy();
      styleH(hD, colD, 20); styleH(hS, colS, 21);
      if(hT2) styleH(hT2, colT2, 24);
      if(hT3) styleH(hT3, colT3, 25);
      if(hT4) styleH(hT4, colT4, 27);

      hD->GetXaxis()->SetRangeUser(ptRelFitLo, ptRelFitHi);
      hD->GetXaxis()->SetLabelSize(0);
      hD->GetYaxis()->SetTitle("d#it{N} / d#it{p}_{T}^{rel} per event");
      hD->GetYaxis()->SetTitleSize(0.055); hD->GetYaxis()->SetTitleOffset(1.25);
      hD->GetYaxis()->SetLabelSize(0.045);
      hD->SetTitle("");
      // a log pad silently drops the whole frame if the max is non-positive
      double ymax = hD->GetMaximum();
      if(ymax <= 0.){ printf("  class %s pT[%.0f,%.0f]: empty data, skipping plot\n",
                             classLabel[ci], ptLo[pi], ptHi[pi]); delete c; continue; }
      // Range for the templates when they are there, tight when they are not --
      // a fixed 1e-4 decade floor leaves the pre-rescan plot as one curve
      // floating in four decades of white space.
      hD->SetMinimum(ymax * ((haveT2 || haveT3) ? 1e-4 : 3e-2));
      hD->SetMaximum(ymax*20.);
      hD->Draw("E");
      if(hT2) hT2->Draw("E same");
      if(hT3) hT3->Draw("E same");
      if(hT4) hT4->Draw("E same");
      hS->Draw("E same");

      TLegend *leg = new TLegend(0.55,0.60,0.93,0.90);
      leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.040);
      leg->AddEntry(hD,  "D: real #mu + real jet", "lp");
      if(hT2) leg->AddEntry(hT2, "T2: real #mu + fake jet", "lp");
      if(hT3) leg->AddEntry(hT3, "T3: fake #mu + real jet", "lp");
      if(hT4) leg->AddEntry(hT4, "T4: fake #mu + fake jet", "lp");
      leg->AddEntry(hS,  Form("S = %s", applied.Data()), "lp");
      leg->Draw();

      TLatex tx; tx.SetNDC(); tx.SetTextSize(0.045);
      tx.DrawLatex(0.18, 0.87, Form("PbPb %s", classLabel[ci]));
      tx.DrawLatex(0.18, 0.81, Form("%.0f < jet #it{p}_{T} < %.0f GeV", ptLo[pi], ptHi[pi]));
      if(!haveT2 || !haveT3){
        tx.SetTextColor(kRed+1); tx.SetTextSize(0.038);
        tx.DrawLatex(0.18, 0.74, "PARTIAL: templates missing");
        tx.SetTextColor(kBlack);
      }

      // bottom: shape comparison of corrected data against MC
      pBot->cd();
      TH1D *hRat = (TH1D*) hSn->Clone(Form("hRat_c%d_p%d", ci, pi));
      hRat->SetDirectory(nullptr);
      if(hMCn) hRat->Divide(hMCn);
      styleH(hRat, colS, 21);
      hRat->GetXaxis()->SetRangeUser(ptRelFitLo, ptRelFitHi);
      hRat->SetTitle("");
      hRat->GetXaxis()->SetTitle("#it{p}_{T}^{rel} [GeV]");
      hRat->GetXaxis()->SetTitleSize(0.105); hRat->GetXaxis()->SetTitleOffset(1.25);
      hRat->GetXaxis()->SetLabelSize(0.090);
      hRat->GetYaxis()->SetTitle("corrected / MC");
      hRat->GetYaxis()->SetTitleSize(0.095); hRat->GetYaxis()->SetTitleOffset(0.72);
      hRat->GetYaxis()->SetLabelSize(0.085); hRat->GetYaxis()->SetNdivisions(505);
      hRat->SetMinimum(0.0); hRat->SetMaximum(2.5);
      hRat->Draw("E");

      TLine *l1 = new TLine(ptRelFitLo, 1.0, ptRelFitHi, 1.0);
      l1->SetLineStyle(2); l1->SetLineColor(kGray+2); l1->Draw();
      hRat->Draw("E same");

      if(!hMCn){
        TLatex tw; tw.SetNDC(); tw.SetTextSize(0.09); tw.SetTextColor(kRed+1);
        tw.DrawLatex(0.20, 0.55, "no MC reference for this bin");
      }

      c->SaveAs(Form("%s/ptRelDecomposition_C%d_jetPt%.0f-%.0f.pdf",
                     outDir, ci+1, ptLo[pi], ptHi[pi]));
      delete c;
    }
  }

  printf("\nfigures written to %s\n", outDir);
  printf("S/D is the surviving real-muon-real-jet fraction; 1.000 means nothing\n"
         "was subtracted, which with templates missing is expected, not a result.\n");
}
