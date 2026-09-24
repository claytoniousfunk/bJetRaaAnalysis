// Inclusive jets-per-Z, PbPb / pp, with calo jets -- NO UNFOLDING.
//
// This is calculateRAA.C's chain with the unfolding step removed (no calo
// response matrices exist yet) and the calo scan files swapped in. Everything
// else -- stitching, normalization windows, Z window, muon efficiencies -- is
// copied from calculateRAA.C so a calo/PF difference is a jet-collection
// difference and not a method difference. Run it with useCalo=false to get the
// PF answer through this identical chain, for the overlay.
//
// CHAIN (calculateRAA.C, minus unfolding and minus fake-jet subtraction)
//   1. h_inclRecoJetPt per trigger sample; PbPb MinBias scaled by 1/N_evt.
//   2. stitch: MinBias below 150, Jet80 150-200, Jet100 above 200. Jet100 sets
//      the scale; Jet80 is normalized to it over [200,500], Jet60 to Jet80 over
//      [150,200], MinBias to Jet80 (PbPb) or Jet60 (pp) over [130,170].
//   3. rebin, divide by bin width and by the eta range (3.2).
//   4. divide by N_Z (dimuon mass 75-105, corrected for muon efficiency).
//   5. ratio PbPb / pp, per centrality class.
//
// NO FAKE-JET SUBTRACTION. The fake-jet estimate is FastJet run on PF
// candidates; there is no calo equivalent, and subtracting the PF one from a
// calo spectrum would be meaningless. calculateRAA.C's nominal has it on, so
// the PF curve from here is NOT its nominal -- it is the no-subtraction twin,
// which is the like-for-like comparison.
//
// FILE CHOICES (calo). Picked from the file names, all 2026-09 calo scans:
//   pp   jets : pp_HighEGJet_caloJets_Jet{60,80,100}HLT_..._2026-9-15
//        MinBias : pp_MinBias_caloJets_..._2026-9-16
//        N_Z  : pp_SingleMuon_caloJets_..._2026-9-15
//   PbPb jets : PbPb_HardProbes_caloJets_Jet{60,80,100}HLT_..._2026-9-15
//        MinBias : PbPb_MinBias_Part1_caloJets_..._2026-9-14   <- the 4-class
//          file, NOT the 2026-9-10 "_ultraFineCentBins" one. This chain indexes
//          C1..C4 as the coarse classes (centrality_4CentBins.h: hiBin 0-20,
//          20-60, 60-100, 100-160 = 0-10, 10-30, 30-50, 50-80%). In the
//          ultrafine file C1..C4 are 5% slices, so it would silently measure
//          0-5,5-10,10-15,15-20% instead. Event counts in the 9-14 file
//          (4.63e6 : 9.25e6 : 9.26e6 : 1.41e7 = 1:2:2:3) confirm the mapping.
//        N_Z  : no PbPb SingleMuon calo scan exists, and the in-repo PF ones
//          are partial (~1400 Z total vs 19423). The dimuon yield does not
//          depend on the jet collection, so N_Z comes from the same PF
//          SingleMuon file calculateRAA.C uses. Cross-check that this is safe:
//          in pp, where both exist, the calo and PF scans give 97932 and 97760
//          Z, agreeing to 0.2%.
//
// WHY THE AXIS ONCE STARTED AT 100 GeV, AND WHY IT NOW STARTS AT 50. The
// 2026-9-16 pp MinBias calo scan had an 80 GeV jet pT threshold: h_inclRawJetPt
// started at 50 but h_inclRecoJetPt at 80 on the same 1899 entries, i.e. the cut
// was on corrected pT. pp takes everything below 150 GeV from MinBias, so below
// 80 the pp spectrum was empty and 80-100 sat in the threshold turn-on
// (MinBias/Jet60 per-event ratio 0.0006 in [80,90] against a 0.0018 plateau from
// 90 up). The manualJEC rescan (2026-9-17) removed that: it reaches 45 GeV and
// holds 9439 entries in 50-100 against the old 909, so the floor is now 50.
//
// Note the 100 GeV floor was never a PF limitation -- it was a calo file
// limitation applied to a shared axis. pp PF MinBias reaches 20 GeV with 5.4e4
// entries in 50-100 and PbPb PF MinBias reaches 15 GeV, so the PF curve was
// being truncated for no reason of its own. Both sides now start at 50.
//
// The [130,170] MinBias normalization window is untouched by any of this and
// stays well above the extended region.
//
// RAW pT. ppRawPt / pbpbRawPt swap h_inclRecoJetPt for h_inclRawJetPt on that
// side. For calo jets "reco" is the forest's own jtpt (pp_scan.C uses the
// built-in value when useCaloJetsOverride is set), so raw-vs-reco here is the
// forest JEC and nothing else: in the pp Jet60 calo sample the [130,170] yield
// falls from 750k to 460k, a ~15-20% scale shift on a steeply falling spectrum.
//
// Turning it on for ONE side only makes the ratio a comparison between two
// different energy scales, so it is a diagnostic, not a measurement -- the
// macro says so at the top of its output. Two further consequences on the pp
// side: the PF files carry no h_inclRawJetPt at all, so the raw variant is
// calo-only; and the pp MinBias normalization window holds 26 raw entries
// against 203 corrected (20% vs 7%), because that scan's jet cut was applied to
// the corrected pT, so the same jets sit at lower raw values.
//
// Usage, from src/newFractionCalculation/:
//   root -l -b -q 'calculateJetsPerZ_caloJets.C'                    // calo, corrected
//   root -l -b -q 'calculateJetsPerZ_caloJets.C(false)'             // PF, same chain
//   root -l -b -q 'calculateJetsPerZ_caloJets.C(true,true)'         // calo, raw pT on pp
//   root -l -b -q 'calculateJetsPerZ_caloJets.C(true,true,true)'    // calo, raw both sides

#include "../../headers/functions/divideByBinwidth.h"
#include "../../headers/plotting/plotStyle.h"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TSystem.h"
#include <cstdio>

// ------------------------------------------------------------------ inputs ---

const char *repo = "/home/clayton/Analysis/code/bJetRaaAnalysis";
const char *sib  = "/home/clayton/Analysis/code/bJetMuonTaggingAnalysis";  // read-only

struct FileSet {
  const char *tag, *label;
  TString pp_jet60, pp_jet80, pp_jet100, pp_MinBias, pp_SingleMuon;
  TString PbPb_jet60, PbPb_jet80, PbPb_jet100, PbPb_MinBias, PbPb_SingleMuon;
};

FileSet caloSet(){
  FileSet s;
  s.tag = "caloJets"; s.label = "calo jets";
  s.pp_jet60        = Form("%s/rootFiles/scanningOuput/pp/pp_HighEGJet_caloJets_Jet60HLT_mu12_pTmu-15to999_tight_deltaR-40_jetTrkMaxFilter_WDecayFilter_2026-9-15.root", repo);
  s.pp_jet80        = Form("%s/rootFiles/scanningOuput/pp/pp_HighEGJet_caloJets_Jet80HLT_mu12_pTmu-15to999_tight_deltaR-40_jetTrkMaxFilter_WDecayFilter_2026-9-15.root", repo);
  s.pp_jet100       = Form("%s/rootFiles/scanningOuput/pp/pp_HighEGJet_caloJets_Jet100HLT_mu12_pTmu-15to999_tight_deltaR-40_jetTrkMaxFilter_WDecayFilter_2026-9-15.root", repo);
  // the manualJEC rescan: reaches 45 GeV where the 2026-9-16 scan started at 80,
  // with 9439 entries in 50-100 against 909. This is what lets the axis go to 50.
  s.pp_MinBias      = Form("%s/rootFiles/scanningOuput/pp/pp_MinBias_caloJets_manualJEC_mu12_pTmu-15to999_tight_deltaR-40_jetTrkMaxFilter_WDecayFilter_2026-9-17.root", repo);
  s.pp_SingleMuon   = Form("%s/rootFiles/scanningOuput/pp/pp_SingleMuon_caloJets_mu12_pTmu-15to999_tight_deltaR-40_jetTrkMaxFilter_WDecayFilter_2026-9-15.root", repo);
  s.PbPb_jet60      = Form("%s/rootFiles/scanningOuput/PbPb/PbPb_HardProbes_caloJets_Jet60HLT_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_2026-9-15.root", repo);
  s.PbPb_jet80      = Form("%s/rootFiles/scanningOuput/PbPb/PbPb_HardProbes_caloJets_Jet80HLT_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_2026-9-15.root", repo);
  s.PbPb_jet100     = Form("%s/rootFiles/scanningOuput/PbPb/PbPb_HardProbes_caloJets_Jet100HLT_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_2026-9-15.root", repo);
  s.PbPb_MinBias    = Form("%s/rootFiles/scanningOuput/PbPb/PbPb_MinBias_Part1_caloJets_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_2026-9-14.root", repo);
  s.PbPb_SingleMuon = Form("%s/rootFiles/scanningOutput/PbPb/latest/PbPb_SingleMuon_mu12_pTmu-15to999_tight_mu12TriggerEfficiencyCorrection_jetTrkMaxFilter_WDecayFilter_2026-5-4.root", sib);
  return s;
}

// exactly the files calculateRAA.C opens, so this reproduces its spectra
FileSet pfSet(){
  FileSet s;
  s.tag = "PFJets"; s.label = "PF jets";
  s.pp_jet60        = Form("%s/rootFiles/scanningOutput/pp/latest/pp_HighEGJet_Jet60HLT_mu12_pTmu-15to999_tight_deltaR-40_jetTrkMaxFilter_2026-3-10.root", sib);
  s.pp_jet80        = Form("%s/rootFiles/scanningOuput/pp/pp_HighEGJet_Jet80HLT_mu12_pTmu-15to999_tight_deltaR-40_jetTrkMaxFilter_WDecayFilter_2026-7-27.root", repo);
  s.pp_jet100       = Form("%s/rootFiles/scanningOuput/pp/pp_HighEGJet_Jet100HLT_mu12_pTmu-15to999_tight_deltaR-40_jetTrkMaxFilter_WDecayFilter_2026-7-27.root", repo);
  s.pp_MinBias      = Form("%s/rootFiles/scanningOuput/pp/pp_MinBias_mu12_pTmu-15to999_tight_deltaR-40_jetTrkMaxFilter_WDecayFilter_2026-7-6.root", repo);
  s.pp_SingleMuon   = Form("%s/rootFiles/scanningOutput/pp/latest/pp_SingleMuon_mu12_pTmu-15to999_tight_deltaR-40_mu12TriggerEfficiencyCorrection_jetTrkMaxFilter_WDecayFilter_2026-5-4.root", sib);
  s.PbPb_jet60      = Form("%s/rootFiles/scanningOutput/PbPb/latest/PbPb_HardProbes_Jet60HLT_mu12_pTmu-15to999_tight_WDecayFilter_2026-3-11.root", sib);
  s.PbPb_jet80      = Form("%s/rootFiles/scanningOutput/PbPb/latest/PbPb_HardProbes_Jet80HLT_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_2026-5-5.root", sib);
  s.PbPb_jet100     = Form("%s/rootFiles/scanningOutput/PbPb/latest/PbPb_HardProbes_Jet100HLT_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_2026-5-5.root", sib);
  s.PbPb_MinBias    = Form("%s/rootFiles/scanningOutput/PbPb/latest/PbPb_MinBias_mu12_pTmu-15to999_tight_WDecayFilter_2026-4-8.root", sib);
  s.PbPb_SingleMuon = Form("%s/rootFiles/scanningOutput/PbPb/latest/PbPb_SingleMuon_mu12_pTmu-15to999_tight_mu12TriggerEfficiencyCorrection_jetTrkMaxFilter_WDecayFilter_2026-5-4.root", sib);
  return s;
}

// ------------------------------------------------------------------ binning ---
// calculateRAA.C's fine axis (60..500) with 50 and 55 prepended, so this is bin
// for bin identical to the nominal from 60 up and the two extra bins continue
// its 5 GeV spacing. Both 50 and 55 land on boundaries of the input 5 GeV axis.
const int    N_edge = 23;
double newAxis[N_edge] = {50,55,60,65,70,75,80,85,90,95,100,110,120,130,140,150,160,180,200,240,280,350,500};
const int    N_coarse = 9;
double coarseAxis[N_coarse] = {50,60,80,100,120,150,200,300,500};

const double Z_lo = 75., Z_hi = 105.;
// muon reconstruction efficiencies, from calculateRAA.C
const double muEff_pp = 0.9708, muEff_C4 = 0.9778, muEff_C3 = 0.9856, muEff_C2 = 0.9069, muEff_C1 = 0.8627;

// ------------------------------------------------------------------ helpers ---

static bool ok = true;

static TFile* openOrFail(TString path)
{
  TFile *f = TFile::Open(path);
  if(!f || f->IsZombie()){ printf("ERROR: cannot open %s\n", path.Data()); ok = false; return nullptr; }
  return f;
}

static TH1D* getOrFail(TFile *f, const char *name)
{
  if(!f) return nullptr;
  TH1D *h = nullptr;
  f->GetObject(name, h);
  if(!h){ printf("ERROR: %s missing from %s\n", name, f->GetName()); ok = false; return nullptr; }
  h = (TH1D*) h->Clone(Form("%s_%p", name, (void*)f));
  h->SetDirectory(nullptr);
  return h;
}

static double winInt(TH1D *h, double lo, double hi)
{
  const double eps = 0.01;
  return h->Integral(h->FindBin(lo + eps), h->FindBin(hi - eps));
}

// Trigger-only stitch: drop MinBias and start each side at the Jet80 threshold
// (150 GeV here). Set by the main function's useMinBias argument.
//
// Why it matters: below 150 GeV both systems are fed by MinBias, and that is
// where the peripheral calo curve goes ragged -- 50-80% scatters 0.42, 0.52,
// 0.46, 0.66, 0.96 across 100-150 GeV with 11-22% errors, then settles to 2.0%
// errors the moment Jet80 takes over at 150. The scatter is consistent with
// those errors, but it is worthless as a measurement and it dominates the eye.
// Peripheral is worst because every class holds the same number of MinBias
// EVENTS while peripheral collisions make far fewer hard jets.
static bool g_useMinBias = true;

// calculateRAA.C's stitchSamples, with the fake-jet argument dropped
// mbRawEntries is the MinBias count in the normalization window BEFORE any
// per-event scaling, for the printout only: the PbPb MinBias spectrum arrives
// here as a per-event rate, so its own window integral is ~1e-4 and says
// nothing about the statistics behind the normalization.
static TH1D* stitchSamples(TH1D *h_jetMB, TH1D *h_jet60, TH1D *h_jet80, TH1D *h_jet100,
                           bool isPP, const char *name, bool verbose, double mbRawEntries)
{
  TH1D *h_return = (TH1D*) h_jet100->Clone(name);
  h_return->SetDirectory(nullptr);

  double jet60_pTmin = 150, jet80_pTmin = 150, jet100_pTmin = 200;
  if(!isPP) jet60_pTmin = jet80_pTmin;   // PbPb does not use the Jet60 sample

  double mbNormLo = 130., mbNormHi = 170.;

  double N_jet100 = winInt(h_jet100, jet100_pTmin, 500.);
  double N_jet80  = winInt(h_jet80,  jet100_pTmin, 500.);
  double N_jet60  = winInt(h_jet60,  jet80_pTmin, jet100_pTmin);
  double N_jetMB  = winInt(h_jetMB,  mbNormLo, mbNormHi);

  if(N_jet80 <= 0. || N_jet60 <= 0. || (g_useMinBias && N_jetMB <= 0.)){
    printf("ERROR: %s: empty normalization window (jet80 %.0f, jet60 %.0f, MinBias %.0f)\n",
           name, N_jet80, N_jet60, N_jetMB);
    ok = false;
    return h_return;
  }

  TH1D *h_jet80_scaled = (TH1D*) h_jet80->Clone(Form("%s_jet80s", name));
  h_jet80_scaled->Scale(N_jet100 / N_jet80);
  double N_jet80_scaled = winInt(h_jet80_scaled, jet80_pTmin, jet100_pTmin);

  TH1D *h_jet60_scaled = (TH1D*) h_jet60->Clone(Form("%s_jet60s", name));
  h_jet60_scaled->Scale(N_jet80_scaled / N_jet60);

  double N_ref = isPP ? winInt(h_jet60_scaled, mbNormLo, mbNormHi)
                      : winInt(h_jet80_scaled, mbNormLo, mbNormHi);

  TH1D *h_jetMB_scaled = (TH1D*) h_jetMB->Clone(Form("%s_MBs", name));
  if(N_jetMB > 0.) h_jetMB_scaled->Scale(N_ref / N_jetMB);

  if(verbose)
    printf("    stitch %-4s: k(jet80)=%.4g  k(jet60)=%.4g  k(MinBias)=%.4g   MinBias [%.0f,%.0f] entries=%.0f (%.1f%% stat)\n",
           name, N_jet100/N_jet80, N_jet80_scaled/N_jet60, N_ref/N_jetMB,
           mbNormLo, mbNormHi, mbRawEntries, mbRawEntries > 0. ? 100./sqrt(mbRawEntries) : -1.);

  // Iterate over h_return's OWN bins and look each source up by pT.
  //
  // This used to loop over h_jet60's bins, take pT from h_jet60, and then read
  // and write content at that same bin INDEX in the other histograms. The pp
  // Jet60 scan (2026-3-10) has 96 bins over 20-500 while Jet80, Jet100 and
  // MinBias have 100 bins over 0-500 -- same 5 GeV width, offset by four bins.
  // So the pT that chose the source was 20 GeV away from the pT the content was
  // written to: with the handover set at 150 the switch actually happened at
  // 130 in the output frame, which is why pp carried MinBias-derived content at
  // 130-150 that should have come from Jet60/Jet80, roughly a factor two off.
  // Indexing by pT is immune to the binnings differing.
  for(int i = 1; i <= h_return->GetNbinsX(); i++){
    double pT = h_return->GetBinCenter(i);
    TH1D *src = nullptr;
    if(pT < jet60_pTmin)                              src = g_useMinBias ? h_jetMB_scaled : nullptr;
    else if(pT > jet60_pTmin && pT < jet80_pTmin)     src = h_jet60_scaled;
    else if(pT > jet80_pTmin && pT < jet100_pTmin)    src = h_jet80_scaled;
    else if(pT > jet100_pTmin)                        src = h_jet100;
    // No source for this bin: zero it. h_return is CLONED from h_jet100, so
    // skipping would silently leave the raw Jet100 content there -- which is
    // what happened on the first trigger-only attempt, filling the region below
    // 150 GeV with an unscaled Jet100/Jet100 ratio that looked like a
    // measurement (0.28 with 0.0006 errors). This also zeroes bins sitting
    // exactly on a threshold, which the strict inequalities above skip.
    if(!src){ h_return->SetBinContent(i, 0.); h_return->SetBinError(i, 0.); continue; }
    const int j = src->FindBin(pT);
    if(j < 1 || j > src->GetNbinsX()){
      h_return->SetBinContent(i, 0.); h_return->SetBinError(i, 0.); continue; }
    h_return->SetBinContent(i, src->GetBinContent(j));
    h_return->SetBinError  (i, src->GetBinError(j));
  }

  delete h_jet80_scaled; delete h_jet60_scaled; delete h_jetMB_scaled;
  return h_return;
}

// rebin to an axis, then per pT and per eta
static TH1D* toSpectrum(TH1D *h, const char *name, int nEdge, double *axis)
{
  TH1D *r = (TH1D*) h->Rebin(nEdge - 1, name, axis);
  r->SetDirectory(nullptr);
  divideByBinwidth(r);
  r->Scale(1./3.2);   // eta range
  return r;
}

// -------------------------------------------------------------------- main ---

void calculateJetsPerZ_caloJets(bool useCalo = true, bool ppRawPt = false, bool pbpbRawPt = false,
                                bool useMinBias = true)
{
  g_useMinBias = useMinBias;
  initPlotStyle();
  ok = true;

  FileSet s = useCalo ? caloSet() : pfSet();
  const char *ppJetHist   = ppRawPt   ? "h_inclRawJetPt" : "h_inclRecoJetPt";
  const char *pbpbJetHist = pbpbRawPt ? "h_inclRawJetPt" : "h_inclRecoJetPt";

  TString scaleTag = "";
  if(ppRawPt && pbpbRawPt) scaleTag = "_rawPt";
  else if(ppRawPt)         scaleTag = "_ppRawPt";
  else if(pbpbRawPt)       scaleTag = "_PbPbRawPt";

  printf("\n=== inclusive jets per Z, PbPb / pp, %s, no unfolding, no fake-jet subtraction ===\n", s.label);
  printf("  jet pT: pp = %s, PbPb = %s\n", ppJetHist, pbpbJetHist);
  if(ppRawPt != pbpbRawPt)
    printf("  WARNING: the two sides use DIFFERENT energy scales (one raw, one JEC-corrected).\n"
           "           This ratio is a scale diagnostic, not a measurement.\n");

  TFile *f_pp_jet60   = openOrFail(s.pp_jet60);
  TFile *f_pp_jet80   = openOrFail(s.pp_jet80);
  TFile *f_pp_jet100  = openOrFail(s.pp_jet100);
  TFile *f_pp_MB      = openOrFail(s.pp_MinBias);
  TFile *f_pp_mu      = openOrFail(s.pp_SingleMuon);
  TFile *f_PbPb_jet60 = openOrFail(s.PbPb_jet60);
  TFile *f_PbPb_jet80 = openOrFail(s.PbPb_jet80);
  TFile *f_PbPb_jet100= openOrFail(s.PbPb_jet100);
  TFile *f_PbPb_MB    = openOrFail(s.PbPb_MinBias);
  TFile *f_PbPb_mu    = openOrFail(s.PbPb_SingleMuon);
  if(!ok) return;

  // ---- pp ----
  TH1D *h_pp_jet60  = getOrFail(f_pp_jet60,  ppJetHist);
  TH1D *h_pp_jet80  = getOrFail(f_pp_jet80,  ppJetHist);
  TH1D *h_pp_jet100 = getOrFail(f_pp_jet100, ppJetHist);
  TH1D *h_pp_jetMB  = getOrFail(f_pp_MB,     ppJetHist);
  TH1D *h_dimuon_pp = getOrFail(f_pp_mu,     "h_dimuonMass");
  if(!ok) return;

  // ---- PbPb, per class ----
  TH1D *h_jetMB[5], *h_jet60[5], *h_jet80[5], *h_jet100[5], *h_dimuon[5], *h_vzMB[5];
  for(int c = 1; c <= 4; c++){
    h_jetMB[c]  = getOrFail(f_PbPb_MB,     Form("%s_C%d", pbpbJetHist, c));
    h_vzMB[c]   = getOrFail(f_PbPb_MB,     Form("h_vz_C%d", c));
    h_jet60[c]  = getOrFail(f_PbPb_jet60,  Form("%s_C%d", pbpbJetHist, c));
    h_jet80[c]  = getOrFail(f_PbPb_jet80,  Form("%s_C%d", pbpbJetHist, c));
    h_jet100[c] = getOrFail(f_PbPb_jet100, Form("%s_C%d", pbpbJetHist, c));
    h_dimuon[c] = getOrFail(f_PbPb_mu,     Form("h_dimuonMass_C%d", c));
  }
  if(!ok) return;

  // raw MinBias counts in the normalization window, kept for the printout
  // before the per-event scaling turns the spectra into rates
  double mbRaw_pp = winInt(h_pp_jetMB, 130., 170.), mbRaw[5];
  for(int c = 1; c <= 4; c++) mbRaw[c] = winInt(h_jetMB[c], 130., 170.);

  // PbPb MinBias to a per-event rate (pp MinBias is not scaled: the stitch
  // renormalizes it anyway, and calculateRAA.C does the same)
  for(int c = 1; c <= 4; c++) h_jetMB[c]->Scale(1./h_vzMB[c]->Integral());

  // ---- stitch ----
  printf("\n  stitching\n");
  TH1D *h_pp = stitchSamples(h_pp_jetMB, h_pp_jet60, h_pp_jet80, h_pp_jet100, true, "h_pp", true, mbRaw_pp);
  TH1D *h_C[5];
  for(int c = 1; c <= 4; c++)
    h_C[c] = stitchSamples(h_jetMB[c], h_jet60[c], h_jet80[c], h_jet100[c], false, Form("h_C%d", c), true, mbRaw[c]);
  if(!ok) return;

  // ---- spectra ----
  TH1D *sp_pp = toSpectrum(h_pp, "sp_pp", N_edge, newAxis);
  TH1D *sp_C[5];
  for(int c = 1; c <= 4; c++) sp_C[c] = toSpectrum(h_C[c], Form("sp_C%d", c), N_edge, newAxis);

  // ---- Z yields ----
  double NZ_pp = h_dimuon_pp->Integral(h_dimuon_pp->GetXaxis()->FindBin(Z_lo), h_dimuon_pp->GetXaxis()->FindBin(Z_hi));
  double NZ[5];
  for(int c = 1; c <= 4; c++)
    NZ[c] = h_dimuon[c]->Integral(h_dimuon[c]->GetXaxis()->FindBin(Z_lo), h_dimuon[c]->GetXaxis()->FindBin(Z_hi));

  NZ_pp /= (muEff_pp*muEff_pp);
  NZ[1] /= (muEff_C1*muEff_C1);
  NZ[2] /= (muEff_C2*muEff_C2);
  NZ[3] /= (muEff_C3*muEff_C3);
  NZ[4] /= (muEff_C4*muEff_C4);

  printf("\n  N_Z (75-105 GeV, muon-efficiency corrected)\n");
  printf("    pp      %9.1f\n", NZ_pp);
  const char *clsLabel[5] = {"", "0-10%", "10-30%", "30-50%", "50-80%"};
  for(int c = 1; c <= 4; c++) printf("    %-7s %9.1f\n", clsLabel[c], NZ[c]);
  if(NZ_pp <= 0.){ printf("ERROR: no Z in pp\n"); return; }

  sp_pp->Scale(1./NZ_pp);
  for(int c = 1; c <= 4; c++){
    if(NZ[c] <= 0.){ printf("ERROR: no Z in %s\n", clsLabel[c]); return; }
    sp_C[c]->Scale(1./NZ[c]);
  }

  // ---- ratio ----
  TH1D *r[5];
  for(int c = 1; c <= 4; c++){
    r[c] = (TH1D*) sp_C[c]->Clone(Form("r_C%d", c));
    r[c]->SetDirectory(nullptr);
    r[c]->Divide(sp_C[c], sp_pp, 1, 1, "");
  }

  printf("\n  jets per Z, PbPb / pp\n    %-12s", "pT [GeV]");
  for(int c = 1; c <= 4; c++) printf(" %16s", clsLabel[c]);
  printf("\n");
  for(int b = 1; b <= r[1]->GetNbinsX(); b++){
    printf("    %4.0f-%-7.0f", r[1]->GetXaxis()->GetBinLowEdge(b), r[1]->GetXaxis()->GetBinUpEdge(b));
    for(int c = 1; c <= 4; c++) printf(" %8.3f+-%-6.3f", r[c]->GetBinContent(b), r[c]->GetBinError(b));
    printf("\n");
  }

  // ---- coarse-axis ratio, as calculateRAA.C also writes ----
  TH1D *rc[5];
  for(int c = 1; c <= 4; c++){
    TH1D *a = toSpectrum(h_C[c], Form("spc_C%d", c), N_coarse, coarseAxis); a->Scale(1./NZ[c]);
    TH1D *b = toSpectrum(h_pp,   Form("spc_pp_%d", c), N_coarse, coarseAxis); b->Scale(1./NZ_pp);
    rc[c] = (TH1D*) a->Clone(Form("r_C%d_r", c));
    rc[c]->SetDirectory(nullptr);
    rc[c]->Divide(a, b, 1, 1, "");
    delete a; delete b;
  }

  // ---- figure ----
  const char *hexC[5] = {"", okabeHex[5], okabeHex[1], okabeHex[2], okabeHex[0]};
  const int   markC[5] = {0, markFilledCircle, markFilledSquare, markFilledDiamond, markCross};

  // y range is 0-2 for a normal run, but a one-sided raw-pT run inflates the
  // ratio and the 100-150 GeV points carry large errors -- a clipped error bar
  // reads as a smaller uncertainty than it is, so grow the axis to hold them
  double yTop = 2.;
  for(int c = 1; c <= 4; c++)
    for(int b = 1; b <= r[c]->GetNbinsX(); b++){
      if(r[c]->GetXaxis()->GetBinLowEdge(b) < newAxis[0]) continue;
      double v = r[c]->GetBinContent(b) + r[c]->GetBinError(b);
      if(v > yTop) yTop = v;
    }
  if(yTop > 2.) yTop *= 1.08;

  TCanvas *canv = new TCanvas(Form("c_jetsPerZ_%s", s.tag), "", 700, 800);
  canv->SetLeftMargin(0.17); canv->SetBottomMargin(0.12); canv->SetTopMargin(0.14); canv->SetRightMargin(0.05);
  TLegend *leg = makeLegend(0.60, 0.60, 0.92, 0.80, 0.036);
  for(int c = 1; c <= 4; c++){
    styleH(r[c], hexC[c], markC[c]);
    if(c == 1){
      r[c]->SetTitle("");
      r[c]->SetMinimum(0.); r[c]->SetMaximum(yTop);
      r[c]->GetXaxis()->SetRangeUser(newAxis[0], newAxis[N_edge-1]);
      r[c]->GetXaxis()->SetTitle("jet #it{p}_{T} [GeV]");
      r[c]->GetXaxis()->SetTitleSize(0.045); r[c]->GetXaxis()->SetLabelSize(0.040);
      r[c]->GetYaxis()->SetTitle("#frac{1}{#it{N}_{Z}^{PbPb}} #frac{d#it{N}_{jet}^{PbPb}}{d#it{p}_{T}} #scale[2.5]{/} #frac{1}{#it{N}_{Z}^{pp}} #frac{d#it{N}_{jet}^{pp}}{d#it{p}_{T}}");
      r[c]->GetYaxis()->SetTitleSize(0.040); r[c]->GetYaxis()->SetTitleOffset(1.55); r[c]->GetYaxis()->SetLabelSize(0.040);
      r[c]->Draw("E1");
    }
    else r[c]->Draw("E1 same");
    leg->AddEntry(r[c], Form("PbPb %s", clsLabel[c]), "lp");
  }
  TLine *one = new TLine(newAxis[0], 1., newAxis[N_edge-1], 1.);
  one->SetLineStyle(7);
  one->Draw();
  for(int c = 1; c <= 4; c++) r[c]->Draw("E1 same");
  leg->Draw();

  TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.036);
  la.DrawLatex(0.17, 0.94, Form("PbPb / pp 5.02 TeV, inclusive %s, anti-#it{k}_{T} #it{R} = 0.4", s.label));
  // keep this inside the frame width: the longer wording ran off the right edge
  la.DrawLatex(0.17, 0.895, Form("no unfolding, no fake-jet subtraction%s",
                                 scaleTag == "_ppRawPt"   ? ", pp raw #it{p}_{T}" :
                                 scaleTag == "_PbPbRawPt" ? ", PbPb raw #it{p}_{T}" :
                                 scaleTag == "_rawPt"     ? ", raw #it{p}_{T} both sides" : ""));

  TString figDir = Form("%s/figures/JetsPerZ", repo);
  gSystem->mkdir(figDir, kTRUE);
  const char *mbTag = g_useMinBias ? "" : "_trigOnly";
  TString figPath = Form("%s/JetsPerZ_%s%s%s_noUnfold.pdf", figDir.Data(), s.tag, scaleTag.Data(), mbTag);
  canv->SaveAs(figPath);

  TString outDir = "./rootFiles/JetsPerZ";
  gSystem->mkdir(outDir, kTRUE);
  TString outPath = Form("%s/histograms_JetsPerZ_%s%s%s_noUnfold.root", outDir.Data(), s.tag, scaleTag.Data(), mbTag);
  TFile *wf = TFile::Open(outPath, "recreate");
  for(int c = 1; c <= 4; c++){ r[c]->Write(Form("r_C%d_fine", c)); rc[c]->Write(Form("r_C%d_r", c)); }
  sp_pp->Write("spectrum_pp");
  for(int c = 1; c <= 4; c++) sp_C[c]->Write(Form("spectrum_C%d", c));
  wf->Close();

  printf("\n  figure  : %s\n  results : %s\n", figPath.Data(), outPath.Data());
}
