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
// MC contains only the (real mu, real jet) pairing, so to compare data against
// it the other pairings are subtracted off the data. Both surviving templates
// are cases where the MUON is uncorrelated with the jet; they differ only in
// whether the jet is real:
//
//   T2 = h_injMuonPtRel_donorJetPt_inject mixed-event PF muon clustered INTO the
//                                         next event and reclustered -- real jets
//   T3 = h_fastJetMuonPtRel_..._RC        mixed muon in a synthetic mixed event --
//                                         combinatorial jets
//
//   S = D - T2 - T3
//
// A (real mu, fake jet) template was built and then dropped: a real tight-ID
// 15 GeV muon is almost always accompanied by a real jet, so that pairing is
// physically rare. Its mixed-event estimator
// (h_realMuonPtRel_mixedFastJetPt, still filled by the scan) measured
// 9.9e-05/event in 0-10% 40-60 GeV, roughly 2% of the data -- but that number
// is inflated by construction, since it pairs a real muon against a mixed event
// full of combinatorial fastJets and so partly counts fake jets we manufactured
// near the muon rather than a physical rate. It is no longer read here.
//
// KNOWN OVERSUBTRACTION: T2's jet is a reco jet, and the reco collection is
// itself part combinatorial, so T2 already contains some (fake mu, fake jet)
// pairs -- the configuration T3 estimates. That piece comes off twice.
//
// NORMALISATION -- two denominators, not one. D and T2 are filled only inside
// if(evtTriggerDecision), so their rate is per TRIGGERED event and they are
// divided by h_vz_triggerOn. The injection fills are ungated, so T2 is divided
// by h_vz. That is not an inconsistency: the injected muon comes from the mixed
// pool and the jets from a donor event, so neither depends on whether THIS
// event fired the trigger, and the per-event rate equals the per-triggered-event
// rate. Dividing all three by h_vz instead would inflate T2 by 1/f_trig -- 8.7x
// in 0-10% rising to 40x in 50-80%, so it would not even cancel across
// centrality classes.
//
// All terms come out of the SAME scan of the SAME file; do not swap one for a
// histogram from a different scan without checking the event selection matches.
//
// ---------------------------------------------------------------------------
// CAVEATS -- read before trusting a number out of this
// ---------------------------------------------------------------------------
// 1. T3 uses FastJet anti-kT R=0.4 jets reclustered from mixed-event PF
//    candidates, RC-background-subtracted. Data and T2 use akCs4PF reco jets.
//    These are different jet definitions. The mixed-event fastJet is a PROXY
//    for the fake component of the reco collection, not the same object, and
//    the two need not have the same ptRel response. All three terms are on the
//    same JEC-corrected, background-subtracted pT scale, so the pT windows do
//    select comparable jets -- but the underlying jet DEFINITION still differs,
//    and that residual is the single largest systematic in the construction.
//    It is not quantified here.
//
// 2. T2's denominator is "all reco jets", i.e. genuine + fake, which is the
//    source of the double subtraction noted above. It also means T2 alone is
//    NOT "fake muons on real jets", despite the name.
//
// 3. GEOMETRY -- LARGELY ADDRESSED for T3, still open for T2. T3 now uses the
//    muon-injection template, where the muon is clustered into the donor event
//    so the axis is pulled onto it exactly as in data. T2 still pairs a muon
//    against an axis built without it. The paragraph below describes the
//    problem the T3 swap fixes, and which T2 still has.
//
//    ORIGINAL NOTE: In D the muon sits inside
//    the jet -- a 15+ GeV muon within dR < 0.4 of a 40-60 GeV jet is clustered
//    into it and the reco jet pT already contains it -- so the axis is pulled
//    toward the muon and the opening angle, hence ptRel, is small. T3 shares
//    that property (its muon IS a constituent). T2 instead pairs a muon
//    against an axis computed WITHOUT it, so the opening angle is whatever
//    random geometry gives; since solid angle grows as dR*d(dR) those pairings
//    pile up near the edge of the cone. Measured in 0-10%, 40-60 GeV, mean
//    ptRel is 1.57 (D), 2.40 (T3), but 5.13 (T2) -- consistent with
//    15 GeV * 0.4 for a muon at the cone edge. T2 therefore puts almost all
//    its weight ABOVE the plotted range and barely subtracts.
//    Fixing this needs the mixed muon inserted into the jet and the axis
//    recomputed before ptRel is formed. Until then these templates do not
//    describe the background they are meant to model.
//
// 4. S can go negative bin-by-bin where the subtraction overshoots. Negative
//    bins are left in place rather than clipped: clipping would bias the
//    integral upward and hide exactly the failure mode worth seeing.
//
// 5. The corrected data / MC comparison is a SHAPE comparison. Both sides are
//    unit-normalised over the fit range because the MC normalisation carries
//    pThat cross-section weights that have no counterpart in the per-event data
//    rate. Do not read the ratio's overall level as a yield statement.
//
// ---------------------------------------------------------------------------
// Usage: root -l -b -q 'plotPtRelTemplateDecomposition.C'
// Run from: src/plots/muonPtRel/
//
// Runs against scans that predate either template: missing terms are reported
// and dropped from the subtraction, and every plot is labelled with which terms
// were actually applied, so a partial result cannot be mistaken for a complete
// one. T2 first appears in scans from 2026-09-08.
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

// D and both templates all come out of this one scan, so the h_vz event counts
// cancel exactly in the per-event normalisation rather than approximately
const char *scanPath =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/"
  "PbPb_SingleMuon_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_"
  "mixedEventPFClustering_fastJetResamples-100_2026-9-8_ultraFineCentBins.root";

// MC reference: PYTHIA+HYDJET, all jets, same muon tag as data
const char *mcPath =
  "/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/scanningOutput/"
  "PYTHIAHYDJET/latest/PYTHIAHYDJET_DiJet_pThat-30_mu12_pTmu-15_tight_"
  "mu12TriggerEfficiencyCorrection_vzReweight_hiBinReweight_hiBinShift-10_"
  "leadingXjetDumpFilter_jetTrkMaxFilter_removeHYDJETjet0p35CutOnGen_WDecayFilter_"
  "weightCut_2026-2-12.root";

const char *nameD  = "h_muptrel_recoJetPt_inclRecoMuonTag_triggerOn";  // measured: all pairings
// T2 is the muon-injection template: a mixed-event muon is added to the NEXT
// centrality-matched event's PF candidates and everything is reclustered. The
// donor is a REAL event, so its jets are real jets -- this is the (fake mu,
// real jet) pairing. Because the muon is clustered in, the jet pT contains it
// and the axis is pulled onto it, which is what the data has and what the
// post-cluster-matched version below does not.
const char *nameT2 = "h_injMuonPtRel_donorJetPt_inject";               // fake mu + real jet, reclustered
// Mixed-event candidates drawn from a 100-event pool and clustered into a
// synthetic event: those jets are combinatorial, so this is the genuine
// (fake mu, fake jet) term. Its muon is already a jet constituent.
const char *nameT3 = "h_fastJetMuonPtRel_fastJetPt_PF_bkgSub_RC";      // fake mu + fake jet
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

// ptRel range shown, and the range over which every normalisation integral is
// taken. Extended past the b-purity fit's 0-3 window because the accidental
// pairing templates peak near 5 GeV (caveat 3): cutting at 3 would hide where
// T2 and T3 actually live and make them look negligible rather than misplaced.
const double ptRelFitLo = 0.0;
const double ptRelFitHi = 5.0;

// Variable binning: 0.2 GeV through the peak region, widening above 2 GeV where
// the data thins out. Uniform 0.2 GeV bins leave the 3-5 GeV range with errors
// larger than the points. Every edge must land on the original 0.1 GeV grid,
// since TH1::Rebin can only merge existing bins -- check that before editing.
//
// Because the widths differ, every drawn histogram is converted to a DENSITY
// (Scale(1,"width")) before plotting, so a wide bin is not read as a tall one.
// That happens after all the integrals and the MC area matching, which are
// counts rather than densities.
const double ptRelBinEdges[] = {0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,
                                2.0,2.4,2.8,3.4,4.2,5.0};
const int nPtRelBins = (int)(sizeof(ptRelBinEdges)/sizeof(double)) - 1;

// Okabe-Ito
const int colD  = TColor::GetColor("#000000");
const int colT2 = TColor::GetColor("#56B4E9");
const int colT3 = TColor::GetColor("#009E73");
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

  // every template and the MC share the original X binning, so rebinning them
  // all onto the same edges here keeps the later Add()/Divide() calls valid.
  // Rebin with an edge array returns a NEW histogram; the original must be
  // released explicitly or it leaks once per panel.
  TH1D *pr = (TH1D*) p->Rebin(nPtRelBins, Form("rb_%s_%s_c%d_p%d", base, tag, ci, pi),
                              ptRelBinEdges);
  pr->SetDirectory(nullptr);
  delete p;
  return pr;
}

// vzBase selects the denominator: "h_vz" for all events, "h_vz_triggerOn" for
// triggered ones. Which is correct depends on the numerator -- see the
// normalisation note in the header.
double classEvents(TFile *f, int ci, const char *vzBase = "h_vz")
{
  double n = 0.;
  for(int si = sliceLo[ci]; si <= sliceHi[ci]; si++){
    TH1D *h = nullptr;
    f->GetObject(Form("%s_C%d", vzBase, si), h);
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

// scanFile: optional override of scanPath, so a single-file or test scan can be
// checked without editing the constant. outSuffix is appended to every output
// filename, to keep a test run from overwriting the nominal figures.
void plotPtRelTemplateDecomposition(const char *scanFile = nullptr,
                                    const char *outSuffix = "")
{
  gStyle->SetOptStat(0);
  gSystem->mkdir(outDir, kTRUE);

  const char *usePath = (scanFile && scanFile[0]) ? scanFile : scanPath;
  TFile *fS = TFile::Open(gSystem->ExpandPathName(usePath));
  if(!fS || fS->IsZombie()){ printf("ERROR: cannot open scan file\n  %s\n", usePath); return; }
  printf("scan file: %s\n", usePath);
  TFile *fM = TFile::Open(mcPath);
  if(!fM || fM->IsZombie()){ printf("ERROR: cannot open MC file\n  %s\n", mcPath); return; }

  // ---- report which templates the input actually carries -------------------
  bool haveT2 = false, haveT3 = false;
  { TH2D *h=nullptr;
    fS->GetObject(Form("%s_C1",nameT2),h); haveT2 = (h!=nullptr); h=nullptr;
    fS->GetObject(Form("%s_C1",nameT3),h); haveT3 = (h!=nullptr); }

  printf("\n=== template availability in the input scan ===\n");
  printf("  D  measured mu + jet     %-40s present\n", nameD);
  printf("  T2 fake mu + real jet   %-40s %s\n", nameT2, haveT2 ? "present" : "*** MISSING ***");
  printf("  T3 fake mu + fake jet   %-40s %s\n", nameT3, haveT3 ? "present" : "*** MISSING ***");
  if(!haveT2)
    printf("\n  NOTE: %s was added to PbPb_pfCandAnalyzer.C on 2026-09-08; scans\n"
           "        older than that do not carry it, and in those the T3 histogram\n"
           "        is also still on the RAW rcSub jet pT scale rather than the JEC\n"
           "        one, so its pT windows do not match D's. The decomposition below\n"
           "        is PARTIAL and every plot is labelled as such.\n", nameT2);

  TString applied = "D";
  if(haveT2) applied += " #minus T2";
  if(haveT3) applied += " #minus T3";

  // ---- per-event rate table ------------------------------------------------
  printf("\n=== per-event tagged-muon rate, ptRel in [%.1f,%.1f] ===\n", ptRelFitLo, ptRelFitHi);
  printf("%-8s %-12s %10s %10s %10s %10s %8s\n",
         "class","jet pT","N_trig","D","T2","T3","S/D");

  for(int ci = 0; ci < NClass; ci++){

    // TWO denominators, because the numerators are not filled over the same
    // event sample. D and T2 are gated on evtTriggerDecision, so their rate is
    // per TRIGGERED event. The injection fills are ungated, so theirs is per
    // event -- and since the injected muon comes from the mixed pool and the
    // jets from a donor event, neither depends on whether THIS event fired,
    // so that rate is also the rate per triggered event. Dividing all three by
    // h_vz would inflate T3 by 1/f_trig, which runs 8.7x in 0-10% to 40x in
    // 50-80% and so would not even cancel across centrality.
    double nEvtTrig = classEvents(fS, ci, "h_vz_triggerOn");
    double nEvtAll  = classEvents(fS, ci, "h_vz");
    if(nEvtTrig <= 0. || nEvtAll <= 0.){
      printf("  class %s: missing h_vz or h_vz_triggerOn, skipping\n", classLabel[ci]); continue; }

    double nEvtMC = 0.;
    { TH1D *h=nullptr; fM->GetObject(Form("h_vz_C%d", ci+1), h); if(h) nEvtMC = h->Integral(); }

    for(int pi = 0; pi < NPt; pi++){

      TH1D *hD  = projectClass(fS, nameD,  ci, pi, "D");
      TH1D *hT2 = haveT2 ? projectClass(fS, nameT2, ci, pi, "T2") : nullptr;
      TH1D *hT3 = haveT3 ? projectClass(fS, nameT3, ci, pi, "T3") : nullptr;
      if(!hD){ printf("  class %s pT[%.0f,%.0f]: no data histogram, skipping\n",
                      classLabel[ci], ptLo[pi], ptHi[pi]); continue; }

      // trigger-gated numerators over triggered events; ungated over all events.
      // Only the injection template (T2) is ungated -- see the header.
      hD->Scale(1./nEvtTrig);
      if(hT2)  hT2->Scale(1./nEvtAll);
      if(hT3)  hT3->Scale(1./nEvtTrig);

      // S = D - T2 - T3 + T4
      TH1D *hS = (TH1D*) hD->Clone(Form("hS_c%d_p%d", ci, pi));
      hS->SetDirectory(nullptr);
      if(hT2) hS->Add(hT2, -1.);
      if(hT3) hS->Add(hT3, -1.);

      double iD  = fitRangeIntegral(hD);
      double iS  = fitRangeIntegral(hS);
      printf("%-8s %4.0f-%-7.0f %10.0f %10.2e %10.2e %10.2e %8.3f\n",
             classLabel[ci], ptLo[pi], ptHi[pi], nEvtTrig,
             iD,
             hT2 ? fitRangeIntegral(hT2) : 0.,
             hT3 ? fitRangeIntegral(hT3) : 0.,
             iD > 0. ? iS/iD : 0.);

      // ---- MC reference ----------------------------------------------------
      TH1D *hMC = nullptr;
      { TH2D *h2 = nullptr;
        fM->GetObject(Form("%s_C%dT0", nameMC, ci+1), h2);
        if(h2){
          int by1 = h2->GetYaxis()->FindBin(ptLo[pi] + 1e-6);
          int by2 = h2->GetYaxis()->FindBin(ptHi[pi] - 1e-6);
          TH1D *hMCraw = h2->ProjectionX(Form("pxMC_c%d_p%d", ci, pi), by1, by2);
          hMCraw->SetDirectory(nullptr);
          hMC = (TH1D*) hMCraw->Rebin(nPtRelBins, Form("rbMC_c%d_p%d", ci, pi), ptRelBinEdges);
          hMC->SetDirectory(nullptr);
          delete hMCraw;
          if(nEvtMC > 0.) hMC->Scale(1./nEvtMC);
        } }

      // MC for the top panel, scaled so its area over the fit range equals the
      // CORRECTED data's. This is a shape overlay: the MC normalisation carries
      // pThat cross-section weights with no counterpart in a per-event data
      // rate, so its absolute level is not meaningful and matching areas is the
      // only way to put the two on one linear axis. Any vertical agreement here
      // is imposed by construction -- only the shape carries information.
      //
      // Scaling to S rather than to D is deliberate: S is the thing the
      // decomposition claims is the real-muon-real-jet contribution, and MC is
      // the prediction for exactly that.
      TH1D *hMCtop = nullptr;
      if(hMC){
        double iMC = fitRangeIntegral(hMC);
        // iS <= 0 means the subtraction removed more than the data contained.
        // Scaling by iS/iMC would then flip or collapse the MC curve and draw a
        // meaningless line, so skip it and say so rather than plot nonsense.
        if(iS > 0. && iMC > 0.){
          hMCtop = (TH1D*) hMC->Clone(Form("hMCtop_c%d_p%d", ci, pi));
          hMCtop->SetDirectory(nullptr);
          hMCtop->Scale(iS/iMC);
        }
        else
          printf("  class %s pT[%.0f,%.0f]: corrected integral %.2e <= 0, "
                 "MC overlay skipped\n", classLabel[ci], ptLo[pi], ptHi[pi], iS);
      }

      // ---- convert to densities for drawing ---------------------------------
      // Must come after every integral above: those are per-event counts, and
      // dividing by bin width first would silently change what they mean. The
      // ratio panel is unaffected either way, since the widths cancel.
      { TH1D *toScale[5] = {hD, hT2, hT3, hS, hMCtop};
        for(int q = 0; q < 5; q++) if(toScale[q]) toScale[q]->Scale(1.0, "width"); }

      // ---- draw ------------------------------------------------------------
      TCanvas *c = new TCanvas(Form("c_c%d_p%d", ci, pi), "", 700, 800);
      TPad *pTop = new TPad("pTop","",0,0.34,1,1);
      TPad *pBot = new TPad("pBot","",0,0,1,0.34);
      pTop->SetBottomMargin(0.02); pTop->SetLeftMargin(0.15); pTop->SetTopMargin(0.07);
      pBot->SetTopMargin(0.02);    pBot->SetLeftMargin(0.15); pBot->SetBottomMargin(0.32);
      pTop->Draw(); pBot->Draw();

      // top: the per-event decomposition
      pTop->cd();
      styleH(hD, colD, 20); styleH(hS, colS, 21);
      if(hT2) styleH(hT2, colT2, 24);
      if(hT3) styleH(hT3, colT3, 25);
      // MC drawn as a line, not points: it is a prediction overlaid on the
      // measurement, and giving it markers makes it read as another dataset.
      if(hMCtop){ styleH(hMCtop, colMC, 1); hMCtop->SetLineWidth(3); }

      hD->GetXaxis()->SetRangeUser(ptRelFitLo, ptRelFitHi);
      hD->GetXaxis()->SetLabelSize(0);
      hD->GetYaxis()->SetTitle("d#it{N} / d#it{p}_{T}^{rel} per event");
      hD->GetYaxis()->SetTitleSize(0.055); hD->GetYaxis()->SetTitleOffset(1.25);
      hD->GetYaxis()->SetLabelSize(0.045);
      hD->SetTitle("");
      double ymax = hD->GetMaximum();
      // area-matched MC can peak above the data if its shape is more sharply
      // peaked, so it has to be allowed to set the ceiling
      if(hMCtop && hMCtop->GetMaximum() > ymax) ymax = hMCtop->GetMaximum();
      if(ymax <= 0.){ printf("  class %s pT[%.0f,%.0f]: empty data, skipping plot\n",
                             classLabel[ci], ptLo[pi], ptHi[pi]); delete c; continue; }

      // Linear axis. Scan every drawn curve for the lowest point in the plotted
      // range so the frame can open downwards when the subtraction drives S
      // negative -- on the old log axis those bins simply vanished, and they are
      // the failure mode most worth seeing. Anchored at 0 otherwise, since a
      // linear axis that does not include zero misreads the relative sizes of
      // the templates, which is the whole point of this panel.
      double ylow = 0.;
      { TH1D *drawn[5] = {hD, hT2, hT3, hS, hMCtop};
        for(int q = 0; q < 5; q++){
          if(!drawn[q]) continue;
          int b1 = drawn[q]->FindBin(ptRelFitLo + 1e-6);
          int b2 = drawn[q]->FindBin(ptRelFitHi - 1e-6);
          for(int b = b1; b <= b2; b++){
            double v = drawn[q]->GetBinContent(b) - drawn[q]->GetBinError(b);
            if(v < ylow) ylow = v;
          } } }
      hD->SetMinimum(ylow < 0. ? ylow*1.15 : 0.);
      hD->SetMaximum(ymax*1.45);   // headroom for the legend
      hD->Draw("E");
      if(hMCtop) hMCtop->Draw("HIST same");
      if(hT2) hT2->Draw("E same");
      if(hT3) hT3->Draw("E same");
      hS->Draw("E same");

      TLegend *leg = new TLegend(0.40,0.575,0.985,0.90);
      leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.031);
      leg->AddEntry(hD,  "D: measured #mu + jet", "lp");
      if(hT2) leg->AddEntry(hT2, "T2: fake #mu + real jet (reclustered)", "lp");
      if(hT3) leg->AddEntry(hT3, "T3: fake #mu + fake jet", "lp");
      leg->AddEntry(hS,  Form("S = %s", applied.Data()), "lp");
      if(hMCtop) leg->AddEntry(hMCtop, "MC (area matched to S)", "l");
      leg->Draw();

      TLatex tx; tx.SetNDC(); tx.SetTextSize(0.045);
      tx.DrawLatex(0.18, 0.87, Form("PbPb %s", classLabel[ci]));
      tx.DrawLatex(0.18, 0.81, Form("%.0f < jet #it{p}_{T} < %.0f GeV", ptLo[pi], ptHi[pi]));
      if(!haveT2 || !haveT3){
        tx.SetTextColor(kRed+1); tx.SetTextSize(0.038);
        tx.DrawLatex(0.18, 0.74, "PARTIAL: templates missing");
        tx.SetTextColor(kBlack);
      }

      // bottom: BOTH data curves against the SAME MC shown on the top panel.
      // Using hMCtop as the common denominator, rather than unit-normalising
      // each ratio separately, is what makes the two curves comparable: S/MC
      // averages to 1 by construction (that is how hMCtop was scaled), so
      // D/MC sits above it by exactly the fraction the subtraction removed.
      // Unit-normalising each would force both to average 1 and hide that.
      pBot->cd();
      TH1D *hRatD = nullptr, *hRatS = nullptr;
      if(hMCtop){
        hRatD = (TH1D*) hD->Clone(Form("hRatD_c%d_p%d", ci, pi)); hRatD->SetDirectory(nullptr);
        hRatS = (TH1D*) hS->Clone(Form("hRatS_c%d_p%d", ci, pi)); hRatS->SetDirectory(nullptr);
        // Divide() would propagate the MC error as if the two were independent
        // measurements of the same thing; here MC is a fixed reference curve,
        // so only the data error belongs on the ratio.
        for(int b = 1; b <= hRatD->GetNbinsX(); b++){
          double m = hMCtop->GetBinContent(b);
          // An empty MC bin makes the ratio undefined, not zero. Parking it at
          // 0 would draw a marker on the axis that reads as a measured value;
          // pushing it below the frame leaves the bin visibly absent instead.
          if(m == 0.){ hRatD->SetBinContent(b,-999); hRatD->SetBinError(b,0);
                       hRatS->SetBinContent(b,-999); hRatS->SetBinError(b,0); continue; }
          hRatD->SetBinContent(b, hRatD->GetBinContent(b)/m);
          hRatD->SetBinError  (b, hRatD->GetBinError(b)/m);
          hRatS->SetBinContent(b, hRatS->GetBinContent(b)/m);
          hRatS->SetBinError  (b, hRatS->GetBinError(b)/m);
        }
        styleH(hRatD, colD, 20);
        styleH(hRatS, colS, 21);
      }

      TH1D *hRatFrame = hRatD ? hRatD : (TH1D*) hD->Clone(Form("hRatF_c%d_p%d", ci, pi));
      if(!hRatD){ hRatFrame->SetDirectory(nullptr); hRatFrame->Reset(); styleH(hRatFrame, colD, 20); }
      hRatFrame->GetXaxis()->SetRangeUser(ptRelFitLo, ptRelFitHi);
      hRatFrame->SetTitle("");
      hRatFrame->GetXaxis()->SetTitle("#it{p}_{T}^{rel} [GeV]");
      hRatFrame->GetXaxis()->SetTitleSize(0.105); hRatFrame->GetXaxis()->SetTitleOffset(1.25);
      hRatFrame->GetXaxis()->SetLabelSize(0.090);
      hRatFrame->GetYaxis()->SetTitle("data / MC");
      hRatFrame->GetYaxis()->SetTitleSize(0.095); hRatFrame->GetYaxis()->SetTitleOffset(0.72);
      hRatFrame->GetYaxis()->SetLabelSize(0.085); hRatFrame->GetYaxis()->SetNdivisions(505);

      // Range set only from bins where the MC denominator is actually populated
      // (>2% of its peak). Out in the tail MC falls to nearly nothing while the
      // data does not, so the ratio there diverges for a reason that says more
      // about MC statistics than about the comparison -- letting those bins set
      // the scale compresses the region the plot is for. They are still drawn,
      // clipped at the frame edge, rather than silently dropped.
      double rmax = 1.5;
      if(hMCtop){
        double mcPeak = hMCtop->GetMaximum();
        TH1D *rr[2] = {hRatD, hRatS};
        for(int q = 0; q < 2; q++){
          if(!rr[q]) continue;
          int b1 = rr[q]->FindBin(ptRelFitLo + 1e-6), b2 = rr[q]->FindBin(ptRelFitHi - 1e-6);
          for(int b = b1; b <= b2; b++){
            if(hMCtop->GetBinContent(b) < 0.02*mcPeak) continue;
            double v = rr[q]->GetBinContent(b) + rr[q]->GetBinError(b);
            if(v > rmax) rmax = v;
          } } }
      hRatFrame->SetMinimum(0.0); hRatFrame->SetMaximum(rmax*1.25);
      hRatFrame->Draw(hRatD ? "E" : "AXIS");

      TLine *l1 = new TLine(ptRelFitLo, 1.0, ptRelFitHi, 1.0);
      l1->SetLineStyle(2); l1->SetLineColor(kGray+2); l1->Draw();
      if(hRatD) hRatD->Draw("E same");
      if(hRatS) hRatS->Draw("E same");

      if(hRatD){
        TLegend *legR = new TLegend(0.18,0.80,0.60,0.97);
        legR->SetBorderSize(0); legR->SetFillStyle(0); legR->SetTextSize(0.080);
        legR->SetNColumns(2);
        legR->AddEntry(hRatD, "D / MC", "lp");
        legR->AddEntry(hRatS, "S / MC", "lp");
        legR->Draw();
      }
      else{
        TLatex tw; tw.SetNDC(); tw.SetTextSize(0.09); tw.SetTextColor(kRed+1);
        tw.DrawLatex(0.20, 0.55, "no MC reference for this bin");
      }

      c->SaveAs(Form("%s/ptRelDecomposition_C%d_jetPt%.0f-%.0f%s.pdf",
                     outDir, ci+1, ptLo[pi], ptHi[pi], outSuffix));
      delete c;
    }
  }

  printf("\nfigures written to %s\n", outDir);
  printf("S/D is the surviving real-muon-real-jet fraction; 1.000 means nothing\n"
         "was subtracted, which with templates missing is expected, not a result.\n");
}
