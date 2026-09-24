// Stitched PbPb calo-jet spectra and jets-per-Z in 10%-wide centrality bins.
//
// Built on the manual-JEC PbPb scans of 2026-09-19/20, which are the first calo
// scans whose jet pT comes from the AK4Calo files rather than the forest's jtpt
// (the forest corrects calo jets with the AK4PF payload; see
// headers/config/config_PbPb.h and the 2026-09-18 jet-by-jet check).
//
// CENTRALITY. These scans use centrality_ultraFineCentBins.h: 18 slices of
// hiBin width 10, i.e. 5% each, covering 0-90%. Confirmed against the files
// themselves -- C1..C18 hold 2.31e6-2.47e6 MinBias events each. A 10% bin is
// therefore exactly two adjacent slices, and class k (0-based) is slices
// 2k+1 and 2k+2, giving nine bins: 0-10 ... 80-90%.
//   NOTE coarseCent.h is NOT used here. Its 4-class mapping stops at slice 16
//   (0-80%) and would silently drop 80-90%.
//
// STITCHING. MinBias below the Jet80 threshold, Jet80 up to the Jet100
// threshold, Jet100 above. Jet100 sets the scale, so the result is in Jet100
// trigger units, not per event; that cancels in any class-to-class ratio and in
// PbPb/pp as long as the Jet100 prescale is class-independent.
//
// The two normalization factors are measured in a COMMON pT window and averaged
// over classes, not fitted per class in each class's own overlap window. That
// was established earlier: a per-class window sits where combinatorial jets
// inflate the central-class MinBias yield and drags the factor down ~9%, which
// imports a centrality dependence into a quantity that is a property of the
// datasets. Both factors are printed per class so the spread is visible.
//
// THRESHOLDS are derived here rather than read from
// triggerThresholds_caloJets.txt, because that file was produced on the OLD
// forest-JEC files and the manual JEC moves every threshold up. Jet80 is
// measured ABSOLUTELY against MinBias (which is unbiased and present here), and
// Jet100 by bootstrap off the Jet80 sample. Criterion: on a 10 GeV grid, the
// lowest edge from which the plateau-normalized efficiency stays within
// nSigma of 1 all the way to the ceiling. There is no Jet60 manual-JEC scan, so
// Jet60 is not used at all -- it was already excluded from the spectra.
//
// JETS PER Z. N_Z must come from the SingleMuon dataset, per centrality class.
// The full-statistics PbPb SingleMuon scan (19385 Z) was run with 4 classes
// only, and the ultra-fine SingleMuon scans in this repo are partial (1432 Z,
// a different luminosity), so 10% N_Z does not exist yet. The user chose to
// wait for a SingleMuon rescan with ultraFineCentBins rather than approximate
// it. This macro therefore:
//   - always produces the stitched spectra, which need no Z at all;
//   - fills in jets-per-Z automatically as soon as a SingleMuon file matching
//     singleMuonGlob below appears with >= minZPerClass Z per class.
// Same for the PbPb/pp ratio: the pp trigger samples still carry the forest's
// PF-sized calo correction and are ~40% low, so the ratio is skipped until the
// manual-JEC pp scans exist. Set the paths below when they do.
//
// NOT DONE: no unfolding (no calo response matrices yet) and no fake-jet
// subtraction (the estimate is FastJet on PF candidates, with no calo
// equivalent). Stated on every figure.
//
// Usage: root -l -b -q 'jetsPerZ_caloJets_10pctBins.C'
// Run from: src/newFractionCalculation/

#include "../../headers/functions/triggerStitch.h"
#include "../../headers/functions/divideByBinwidth.h"
#include "../../headers/plotting/plotStyle.h"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TSystem.h"
#include "TMath.h"
#include <cstdio>
#include <vector>

const char *repo = "/home/clayton/Analysis/code/bJetRaaAnalysis";

const char *f_MB_path   = "rootFiles/scanningOuput/PbPb/PbPb_MinBias_Part1_caloJets_manualJEC_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_2026-9-19_ultraFineCentBins.root";
const char *f_j80_path  = "rootFiles/scanningOuput/PbPb/PbPb_HardProbes_caloJets_manualJEC_Jet80HLT_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_2026-9-19_ultraFineCentBins.root";
const char *f_j100_path = "rootFiles/scanningOuput/PbPb/PbPb_HardProbes_caloJets_manualJEC_Jet100HLT_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_2026-9-20_ultraFineCentBins.root";

// The 2026-09-24 SingleMuon rescan in ultra-fine bins (19431 Z over 0-90%).
const char *singleMuonGlob = "rootFiles/scanningOuput/PbPb/PbPb_SingleMuon_*ultraFineCentBins*.root";
const long  minZPerClass   = 400;   // below this a 10% N_Z is too noisy to use

// ---- pp side ---------------------------------------------------------------
// Both manual-JEC: the forest applies the AK4PF payload to calo jets, so a
// forest-JEC pp file mixed with manual-JEC PbPb would shift the ratio by ~20%.
// There is no manual-JEC pp Jet100 scan yet, so pp is MinBias + Jet80 only and
// its reach is whatever Jet80 statistics allow -- the ratio is reported only
// where the pp spectrum still has entries.
const char *f_pp_MB_path  = "rootFiles/scanningOuput/pp/pp_MinBias_caloJets_manualJEC_mu12_pTmu-15to999_tight_deltaR-40_jetTrkMaxFilter_WDecayFilter_2026-9-23.root";
const char *f_pp_j80_path = "rootFiles/scanningOuput/pp/pp_HighEGJet_caloJets_manualJEC_Jet80HLT_mu12_pTmu-15to999_tight_deltaR-40_jetTrkMaxFilter_WDecayFilter_2026-9-2.root";
// N_Z only, so the JEC of this file is irrelevant -- the dimuon mass does not
// depend on the jet calibration.
const char *f_pp_mu_path  = "rootFiles/scanningOuput/pp/pp_SingleMuon_caloJets_mu12_pTmu-15to999_tight_deltaR-40_jetTrkMaxFilter_WDecayFilter_2026-9-15.root";

// Drop MinBias from BOTH stitches and use the triggered samples only, the same
// convention as jetsPerZ_caloJets_stitched_trigOnly_overlay.pdf in the
// 2026-09-17 presentation (stitchSpectra_caloJets.C, useMinBias = false).
//
// Each side then starts at its own Jet80 threshold instead of the reporting
// floor, so the reach is shorter, but the sample whose normalisation is least
// trustworthy is gone. That matters here: the PbPb MinBias calo spectrum has a
// selection turn-on at 65-75 GeV, and MinBias enters the stitch through a
// MinBias->Jet80 factor of order 1e9 on the PbPb side against 6.5e3 on the pp
// side -- the two are per-event rate and raw counts respectively, so the two
// MinBias pieces are not normalised the same way and do not cancel in the
// ratio. Above the Jet80 thresholds both sides are Jet80/Jet80 and Jet100 (PbPb
// only), which is a like-for-like comparison.
//
// MinBias is still USED, but only to derive the Jet80 thresholds absolutely --
// that is a turn-on shape, not a normalisation, so it is unaffected.
const bool useMinBias = false;

// Use the PbPb Jet100 sample as well.
//
// MUST be false while the pp side has no manual-JEC Jet100 scan. The PbPb
// stitch puts everything in Jet100 trigger units by scaling its Jet80 piece by
// k(80->100) = 2.138; pp, having only Jet80, stays in Jet80 units. The ratio of
// two spectra in different trigger units carries that factor and is not a
// measurement. Reproducing the 2026-09-17 trigOnly figure showed the symptom:
// its 50-80% sits at ~0.73 near 100 GeV where this macro gave 0.53, and the
// integrated 0-10% per-Z yields differ by 2.07 (PbPb) against 2.36 (pp) -- the
// same quantity computed two ways should not move by different factors.
// stitchSpectra_caloJets.C avoids this with matchPPToClass, pairing each PbPb
// sample with the SAME pp sample; that needs a pp Jet100 scan, which exists
// only in the forest-JEC generation.
//
// With this false both sides are Jet80/Jet80 throughout, which is like-for-like
// at the cost of reach.
const bool useJet100 = true;

// Put the output on the PF chain's fixed axis instead of deriving per-class
// edges, so a calo/PF ratio can actually be formed bin by bin.
//
// The two chains bin differently for no physics reason: the PF chain
// (calculateJetsPerZ_caloJets.C) uses one hand-chosen axis for every class,
// while this one derives edges per class from stitchDeriveEdges targeting 5%
// statistical error with a 15% resolution floor. Adaptive binning is the better
// choice for a standalone calo figure and the wrong one for a comparison --
// misaligned edges make a per-bin ratio impossible.
//
// Every edge below is a multiple of the underlying 5 GeV grid, so the rebin is
// exact. Edges below the class floor are dropped, which above the calo Jet80
// threshold leaves {120, 150, 200, 300, 500}.
const bool   useFixedAxis = true;
const int    N_fixedAxis  = 9;
const double fixedAxis[N_fixedAxis] = {50, 60, 80, 100, 120, 150, 200, 300, 500};

// pp Jet80 threshold, at file scope so the per-class binning can respect it.
// See the block where it is used for the measured turn-on this comes from.
const double ppThr80 = 120.;

const double Z_lo = 75., Z_hi = 105.;

// Dimuon reconstruction efficiency, squared per Z, following the values used in
// stitchSpectra_caloJets.C and calculateBJetsPerZ.cc.
//
// CAVEAT: those were measured for the FOUR nominal classes. There is no
// per-10%-class determination, so each 10% class takes the value of the nominal
// class containing it: 0-10 -> C1, 10-30 -> C2, 30-50 -> C3, 50-80 -> C4, and
// 80-90 reuses C4 (it lies outside the nominal range entirely). Within 10-30%
// and 30-50% this assumes the efficiency is flat across the two halves, which
// is the main approximation in the ratio's normalisation.
const double muEff_pp = 0.9708;
const double muEff10[9] = {0.8627, 0.9069, 0.9069, 0.9856, 0.9856,
                           0.9778, 0.9778, 0.9778, 0.9778};
// the four nominal classes, as measured (see muonReconstructionEfficiencyCalculator)
const double muEffCoarse[4] = {0.8627, 0.9069, 0.9856, 0.9778};
static const double *g_muEff = muEff10;

const char *outDir  = "/home/clayton/Analysis/code/bJetRaaAnalysis/figures/JetsPerZ";
const char *outFile10     = "./rootFiles/JetsPerZ/jetsPerZ_caloJets_10pctBins.root";
const char *outFileCoarse = "./rootFiles/JetsPerZ/jetsPerZ_caloJets_coarseBins.root";
static const char *outFile = outFile10;

// ------------------------------------------------------------- parameters ---
const int    NCls        = 9;      // 10% bins, 0-10 .. 80-90
// Reporting floor. This was 50 on the assumption that it was the scans' own jet
// pT cut, but the stitched PbPb calo spectrum is EMPTY below 65 GeV and PEAKS
// at 70-75 before falling (0-10%: 0, 0, 0, 2.31e6, 7.25e6, 4.16e6, 2.66e6 in
// the 5 GeV bins from 50), i.e. there is a selection turn-on around 65-75 that
// the stitch cannot repair. Reporting from 50 put a bin spanning that turn-on
// at the bottom of every class and made PbPb/pp read 0.019 there against 0.66
// in the next bin up. 80 GeV clears it with margin and matches the floor the
// b-jet chain quotes from.
const double ptFloor     = 80.;
const double ptCeiling   = 500.;
const double etaRange    = 3.2;    // |eta| < 1.6

// threshold criterion
const double effTarget   = 0.98;
const double nSigmaBelow = 2.0;
const double turnOnRebin = 10.;    // GeV, coarser than the 5 GeV scan grid
// Plateau reference windows. Jet80-vs-MinBias must be read where MinBias still
// has events (it is a per-event sample with a hard 5 GeV grid and empties out
// above ~250 GeV), while the Jet100 bootstrap runs on the jet-triggered sample
// and can use a genuinely high window.
const double platLo_abs  = 130., platHi_abs  = 200.;
const double platLo_boot = 200., platHi_boot = 300.;

// common normalization windows: above the fake-dominated region, inside the
// plateau of both samples in every class
const double kWinMB_lo = 150., kWinMB_hi = 200.;   // MinBias -> Jet80
const double kWin80_lo = 200., kWin80_hi = 280.;   // Jet80   -> Jet100

// binning
const double targetRelErr  = 0.05;
const double minWidthFrac  = 0.15;   // resolution floor; see JER (sigma ~ 0.11-0.22)
const double maxBinWidth   = 80.;

static const char *clsLabel[NCls] = {"0-10%","10-20%","20-30%","30-40%","40-50%",
                                     "50-60%","60-70%","70-80%","80-90%"};
static const char *clsTag[NCls]   = {"0to10","10to20","20to30","30to40","40to50",
                                     "50to60","60to70","70to80","80to90"};

// ---- class map -------------------------------------------------------------
// The scans are 18 slices of 5% covering 0-90%. A class is a contiguous run of
// slices, so the same chain produces either the nine 10% classes or the four
// nominal ones without touching anything downstream.
//
// The nominal map matches coarseCent.h: 0-10 = slices 1-2, 10-30 = 3-6,
// 30-50 = 7-10, 50-80 = 11-16. Slices 17-18 (80-90%) are NOT in any nominal
// class and are dropped there, exactly as coarseCent.h drops them.
static int  g_nCls = NCls;
static int  g_lo[NCls], g_hi[NCls];
static const char *g_label[NCls], *g_tag[NCls];

static const char *coarseLabel[4] = {"0-10%","10-30%","30-50%","50-80%"};
static const char *coarseTag  [4] = {"0to10","10to30","30to50","50to80"};
static const int   coarseLo   [4] = {1, 3, 7, 11};
static const int   coarseHi   [4] = {2, 6, 10, 16};

static void setClassMap(bool coarse)
{
  if(coarse){
    g_nCls = 4;
    for(int k = 0; k < 4; k++){
      g_lo[k] = coarseLo[k]; g_hi[k] = coarseHi[k];
      g_label[k] = coarseLabel[k]; g_tag[k] = coarseTag[k];
    }
  }
  else{
    g_nCls = NCls;
    for(int k = 0; k < g_nCls; k++){
      g_lo[k] = 2*k + 1; g_hi[k] = 2*k + 2;
      g_label[k] = g_label[k]; g_tag[k] = g_tag[k];
    }
  }
}

// sum the 5% slices making up class k; nullptr if any is missing
static TH1D* classSum(TFile *f, const char *fmt, int k, const char *name)
{
  TH1D *s = nullptr;
  for(int i = g_lo[k]; i <= g_hi[k]; i++){
    TH1D *a = nullptr;
    f->GetObject(Form(fmt, i), a);
    if(!a){ printf("ERROR: %s missing\n", Form(fmt, i)); delete s; return nullptr; }
    if(!s){ s = (TH1D*) a->Clone(name); s->SetDirectory(nullptr); }
    else     s->Add(a);
  }
  return s;
}

// Lowest 10-GeV edge from which num/den stays within nSigma of effTarget*plateau
// up to ptCeiling. Returns -1 if no such edge exists.
//
// Two tolerances, both carried over from deriveTriggerThresholds_caloJets.C and
// both necessary here: the MinBias denominator runs out of events well before
// 300 GeV, so (a) only bins measured to maxRelErrForFlat or better can veto an
// edge -- a bin with 30% errors says nothing about flatness -- and (b) a few
// mild downward fluctuations are expected among many bins and are allowed,
// while a single catastrophic one is not. Without these a lone high-pT
// fluctuation pushes the threshold to 400+ GeV, which is what happened first.
static double turnOnThreshold(TH1D *num, TH1D *den, double pLo, double pHi,
                              double &plateau, const char *label, bool verbose)
{
  const double maxRelErrForFlat = 0.05;
  const double nSigmaCatastrophic = 4.0;
  TH1D *n = (TH1D*) num->Clone(Form("n_%s", label)); n->SetDirectory(nullptr);
  TH1D *d = (TH1D*) den->Clone(Form("d_%s", label)); d->SetDirectory(nullptr);
  int group = (int) (turnOnRebin / n->GetBinWidth(1) + 0.5);
  if(group > 1 && n->GetNbinsX() % group == 0){ n->Rebin(group); d->Rebin(group); }

  const double eps = 0.01;
  // The plateau window has to sit where the denominator still has events. In
  // the most peripheral class MinBias is empty above ~130 GeV, so slide the
  // window down until it holds enough to measure a plateau, and say which one
  // was used -- a plateau read off 3 events is not a plateau.
  // A window is only accepted with enough events behind it. Sliding further
  // down would land inside the turn-on itself and return a "plateau" of 0.06,
  // which then yields a meaningless threshold -- better to report the class as
  // unmeasured (-2) and let the caller substitute a measured neighbour.
  const double minPlateauCounts = 200.;
  double nP = 0., dP = 0., usedLo = pLo, usedHi = pHi;
  bool measured = false;
  for(double scale : {1.0, 0.85}){
    usedLo = pLo*scale; usedHi = pHi*scale;
    if(usedLo < ptFloor) break;
    nP = n->Integral(n->FindBin(usedLo + eps), n->FindBin(usedHi - eps));
    dP = d->Integral(d->FindBin(usedLo + eps), d->FindBin(usedHi - eps));
    if(dP >= minPlateauCounts && nP > 0.){ measured = true; break; }
  }
  if(!measured){
    printf("    %-22s NOT MEASURABLE: plateau window holds %.0f events (need %.0f)\n",
           label, dP, minPlateauCounts);
    delete n; delete d; return -2.;
  }
  plateau = nP / dP;

  int bLo = n->FindBin(ptFloor + eps), bHi = n->FindBin(ptCeiling - eps);
  // how many bins above the candidate edge carry enough precision to judge
  int nTestable = 0;
  for(int bb = bLo; bb <= bHi; bb++){
    double dd = d->GetBinContent(bb);
    if(dd <= 0.) continue;
    double e = n->GetBinContent(bb)/dd, se = sqrt(TMath::Max(e*(1.-e)/dd, 1e-12));
    if(e > 0. && se/e <= maxRelErrForFlat) nTestable++;
  }
  const int allowedLow = (int) TMath::Max(1., 0.05*nTestable);

  double best = -1.;
  for(int b = bLo; b <= bHi; b++){
    int nLow = 0; bool catastrophic = false;
    for(int bb = b; bb <= bHi; bb++){
      double nn = n->GetBinContent(bb), dd = d->GetBinContent(bb);
      if(dd <= 0.) continue;                      // no data, not a failure
      double e = nn/dd, se = sqrt(TMath::Max(e*(1.-e)/dd, 1e-12));
      if(e <= 0.) continue;
      if(e < effTarget*plateau - nSigmaCatastrophic*se){ catastrophic = true; break; }
      if(se/e > maxRelErrForFlat) continue;       // too imprecise to veto
      if(e < effTarget*plateau - nSigmaBelow*se) nLow++;
    }
    if(!catastrophic && nLow <= allowedLow){ best = n->GetXaxis()->GetBinLowEdge(b); break; }
  }
  if(verbose) printf("    %-22s plateau %.4g on [%.0f,%.0f]   threshold %.0f GeV\n",
                     label, plateau, usedLo, usedHi, best);
  delete n; delete d;
  return best;
}

// coarse = true -> the four nominal classes (0-10, 10-30, 30-50, 50-80), built
// from the SAME ultra-fine manual-JEC scans, so it is directly comparable to a
// PF result in those classes. coarse = false -> the nine 10% classes.
void jetsPerZ_caloJets_10pctBins(bool coarse = false)
{
  setClassMap(coarse);
  g_muEff = coarse ? muEffCoarse : muEff10;
  outFile = coarse ? outFileCoarse : outFile10;
  initPlotStyle();
  gSystem->mkdir(outDir, kTRUE);
  gSystem->mkdir("./rootFiles/JetsPerZ", kTRUE);

  // ---- input audit ---------------------------------------------------------
  // Every histogram that enters a SPECTRUM must be manual JEC: the forest
  // applies the AK4PF payload to calo jets and its closure is 0.78 in pp but
  // 0.73-0.90 in PbPb, so the error is centrality dependent and does NOT cancel
  // in PbPb/pp. Files used only to count Z are exempt -- the dimuon mass does
  // not depend on the jet calibration.
  {
    struct InRec { const char *role; const char *path; bool spectrum; };
    const InRec in[] = {
      {"PbPb MinBias",    f_MB_path,      true },
      {"PbPb Jet80",      f_j80_path,     true },
      {"PbPb Jet100",     f_j100_path,    true },
      {"pp   MinBias",    f_pp_MB_path,   true },
      {"pp   Jet80",      f_pp_j80_path,  true },
      {"pp   SingleMuon", f_pp_mu_path,   false},
    };
    printf("\n  input audit (spectra must be manualJEC; N_Z files exempt)\n");
    bool bad = false;
    for(const auto &r : in){
      TString b = gSystem->BaseName(r.path);
      const bool mj = b.Contains("manualJEC");
      const bool uf = b.Contains("ultraFineCentBins");
      if(r.spectrum && !mj) bad = true;
      printf("    %-14s %-9s %-13s %s\n", r.role,
             mj ? "manualJEC" : "forestJEC",
             uf ? "ultraFine" : "(no centrality)",
             (r.spectrum && !mj) ? "<-- SPECTRUM ON FOREST JEC" : "");
    }
    // the PbPb Z file is globbed, so it is reported later with N_Z
    if(bad){
      printf("  ERROR: a spectrum is on forest JEC; the ratio would be distorted.\n");
      return;
    }
  }

  TFile *fMB  = TFile::Open(Form("%s/%s", repo, f_MB_path));
  TFile *f80  = TFile::Open(Form("%s/%s", repo, f_j80_path));
  TFile *f100 = TFile::Open(Form("%s/%s", repo, f_j100_path));
  if(!fMB || fMB->IsZombie() || !f80 || f80->IsZombie() || !f100 || f100->IsZombie()){
    printf("ERROR: cannot open one of the PbPb scans\n"); return;
  }
  { TH1D *g = nullptr; fMB->GetObject("h_nEventsNoJetTrigSel", g);
    if(g && g->GetMean() > 0.5){ printf("ERROR: the MinBias scan is jet-triggered; not a valid reference\n"); return; } }

  printf("\n  PbPb calo jets, manual JEC, 10%% centrality bins (two 5%% slices each)\n");

  // ------------------------------------------------------- per-class input ---
  TH1D *spMB[NCls], *sp80[NCls], *sp100[NCls];
  double nEvt[NCls];
  double thr80[NCls], thr100[NCls];

  printf("\n  trigger turn-on, derived on these files (Jet80 absolute vs MinBias, Jet100 off Jet80)\n");
  for(int k = 0; k < g_nCls; k++){
    spMB[k]  = classSum(fMB,  "h_inclRecoJetPt_C%d", k, Form("mb_%d", k));
    sp80[k]  = classSum(f80,  "h_inclRecoJetPt_C%d", k, Form("j80_%d", k));
    sp100[k] = classSum(f100, "h_inclRecoJetPt_C%d", k, Form("j100_%d", k));
    TH1D *vz = classSum(fMB,  "h_vz_C%d", k, Form("vz_%d", k));
    if(!spMB[k] || !sp80[k] || !sp100[k] || !vz) return;
    nEvt[k] = vz->Integral();
    delete vz;

    // MinBias is a rate; the trigger samples stay in counts and are scaled to
    // Jet100 units below
    spMB[k]->Scale(1. / nEvt[k]);

    printf("  %s  (N_evt = %.3g)\n", g_label[k], nEvt[k]);
    TH1D *mbAll = classSum(fMB, "h_leadJetPt_all_C%d",   k, Form("mbAll_%d", k));
    TH1D *mb80  = classSum(fMB, "h_leadJetPt_jet80_C%d", k, Form("mb80_%d", k));
    TH1D *d80   = classSum(f80, "h_leadJetPt_jet80_C%d", k, Form("d80_%d", k));
    TH1D *n100  = classSum(f80, "h_leadJetPt_jet80_and_jet100_C%d", k, Form("n100_%d", k));
    if(!mbAll || !mb80 || !d80 || !n100) return;
    double pl;
    thr80[k]  = turnOnThreshold(mb80, mbAll, platLo_abs,  platHi_abs,  pl, Form("Jet80 abs %s",  g_tag[k]), true);
    thr100[k] = turnOnThreshold(n100, d80,   platLo_boot, platHi_boot, pl, Form("Jet100/Jet80 %s", g_tag[k]), true);
    delete mbAll; delete mb80; delete d80; delete n100;

    // In the most peripheral classes MinBias has too few high-pT jets to show a
    // Jet80 plateau at all. The measured thresholds fall monotonically with
    // peripherality (120 -> 100 GeV), so carrying the last measured value
    // forward is the conservative choice: it sits ABOVE the true turn-on.
    if(thr80[k] == -2.){
      if(k == 0){ printf("ERROR: Jet80 unmeasurable in the most central class\n"); return; }
      thr80[k] = thr80[k-1];
      printf("    Jet80 abs %-12s using %.0f GeV carried over from %s (not measured here)\n",
             g_tag[k], thr80[k], clsLabel[k-1]);
    }
    if(thr100[k] == -2.){
      if(k == 0){ printf("ERROR: Jet100 unmeasurable in the most central class\n"); return; }
      thr100[k] = thr100[k-1];
      printf("    Jet100    %-12s using %.0f GeV carried over from %s (not measured here)\n",
             g_tag[k], thr100[k], clsLabel[k-1]);
    }
    if(thr80[k] < 0. || thr100[k] < 0.){ printf("ERROR: no threshold in %s\n", g_label[k]); return; }
    if(thr100[k] <= thr80[k]){
      printf("ERROR: %s has Jet100 threshold %.0f <= Jet80 %.0f; stitch order broken\n",
             g_label[k], thr100[k], thr80[k]); return;
    }
  }

  // ------------------------------------------ common-window scale factors ---
  // measured on the same summed class histograms used for the spectra
  auto commonK = [&](TH1D **low, TH1D **high, double lo, double hi, const char *label){
    const double eps = 0.01;
    double sw = 0., swk = 0.;
    printf("\n  %s yield ratio in the common window [%.0f,%.0f]\n", label, lo, hi);
    for(int k = 0; k < g_nCls; k++){
      double nl = low[k] ->Integral(low[k] ->FindBin(lo + eps), low[k] ->FindBin(hi - eps));
      double nh = high[k]->Integral(high[k]->FindBin(lo + eps), high[k]->FindBin(hi - eps));
      if(nl <= 0. || nh <= 0.){ printf("    %-7s empty window\n", g_label[k]); continue; }
      // MinBias is already a per-event rate, so its "counts" for the error are
      // recovered from the rate and the event count
      double effLow = (low == spMB) ? nl * nEvt[k] : nl;
      double kk = nh/nl, sk = kk*sqrt(1./TMath::Max(effLow, 1.) + 1./nh);
      printf("    %-7s %12.4g +- %.2g\n", g_label[k], kk, sk);
      sw += 1./(sk*sk); swk += kk/(sk*sk);
    }
    if(sw <= 0.) return -1.;
    double kAvg = swk/sw;
    printf("    weighted average applied to every class: %.4g +- %.2g\n", kAvg, sqrt(1./sw));
    return kAvg;
  };

  double k80to100 = commonK(sp80, sp100, kWin80_lo, kWin80_hi, "Jet80 -> Jet100");
  double kMBto80  = commonK(spMB, sp80,  kWinMB_lo, kWinMB_hi, "MinBias -> Jet80");
  if(k80to100 <= 0. || kMBto80 <= 0.){ printf("ERROR: a normalization window failed\n"); return; }

  // ----------------------------------------------------------- stitching ---
  printf("\n  stitching (Jet100 sets the scale; result is in Jet100 trigger units)\n");
  TH1D *stitched[NCls];
  for(int k = 0; k < g_nCls; k++){
    std::vector<StitchSample> s;
    if(useMinBias){
      s = { {"MinBias", ptFloor,   spMB[k],  kMBto80},
            {"Jet80",   thr80[k],  sp80[k],  k80to100},
            {"Jet100",  thr100[k], sp100[k], 0.} };
      printf("  %s: MinBias %.0f-%.0f, Jet80 %.0f-%.0f, Jet100 %.0f-%.0f\n", g_label[k],
             ptFloor, thr80[k], thr80[k], thr100[k], thr100[k], ptCeiling);
    }
    else if(useJet100){
      s = { {"Jet80",  thr80[k],  sp80[k],  k80to100},
            {"Jet100", thr100[k], sp100[k], 0.} };
      printf("  %s: Jet80 %.0f-%.0f, Jet100 %.0f-%.0f  (MinBias dropped)\n", g_label[k],
             thr80[k], thr100[k], thr100[k], ptCeiling);
    }
    else{
      // Jet80 only, matching what the pp side can offer
      s = { {"Jet80", thr80[k], sp80[k], 0.} };
      printf("  %s: Jet80 %.0f-%.0f  (MinBias and Jet100 dropped)\n",
             g_label[k], thr80[k], ptCeiling);
    }
    stitched[k] = stitchTriggerSamples(s, Form("stitched_%s", g_tag[k]),
                                       useMinBias ? 50. : thr80[k], ptCeiling, true);
    if(!stitched[k]){ printf("ERROR: stitch failed in %s\n", g_label[k]); return; }
  }

  // ------------------------------------------------- per-class rebinning ---
  printf("\n  per-class binning (target %.0f%% stat, resolution floor %.0f%% of pT)\n",
         100.*targetRelErr, 100.*minWidthFrac);
  TH1D *spec[NCls];
  for(int k = 0; k < g_nCls; k++){
    // With MinBias dropped the spectrum starts at the class's own Jet80
    // threshold, so binning from ptFloor would put empty bins at the bottom.
    const double floorK = useMinBias
      ? ptFloor
      : TMath::Max(TMath::Max(ptFloor, thr80[k]), ppThr80);
    std::vector<double> e;
    if(useFixedAxis){
      for(int i = 0; i < N_fixedAxis; i++)
        if(fixedAxis[i] >= floorK - 1e-6 && fixedAxis[i] <= ptCeiling + 1e-6)
          e.push_back(fixedAxis[i]);
    }
    else e = stitchDeriveEdges({stitched[k]}, floorK, ptCeiling,
                               targetRelErr, minWidthFrac, maxBinWidth);
    if(e.size() < 2){ printf("ERROR: no binning in %s\n", g_label[k]); return; }
    spec[k] = rebinTo(stitched[k], (int)e.size() - 1, e.data(), Form("spec_%s", g_tag[k]));
    spec[k]->Scale(1. / etaRange);
    divideByBinwidth(spec[k]);
    printf("  %-7s %2d bins, %.0f-%.0f GeV\n", g_label[k], (int)e.size() - 1, e.front(), e.back());
  }

  // ------------------------------------------------------------ jets per Z ---
  //
  // The Z gate is PER CLASS, not global. It used to demand minZPerClass in every
  // one of the nine classes and abandon the whole calculation otherwise, which
  // threw away six good classes because 80-90% holds 23 Z -- peripheral classes
  // are intrinsically Z-poor and always will be.
  printf("\n  jets per Z\n");
  TString found = gSystem->GetFromPipe(Form("ls -1t %s/%s 2>/dev/null | head -1", repo, singleMuonGlob));
  found = found.Strip(TString::kBoth);

  double nZ[NCls];
  bool   haveZ[NCls];
  int    nGood = 0;
  for(int k = 0; k < g_nCls; k++){ nZ[k] = 0.; haveZ[k] = false; }

  if(found.Length()){
    TFile *fZ = TFile::Open(found);
    if(fZ && !fZ->IsZombie()){
      printf("    SingleMuon file: %s\n", gSystem->BaseName(found.Data()));
      for(int k = 0; k < g_nCls; k++){
        TH1D *m = classSum(fZ, "h_dimuonMass_C%d", k, Form("z_%d", k));
        if(!m) continue;
        nZ[k] = m->Integral(m->FindBin(Z_lo + 0.01), m->FindBin(Z_hi - 0.01));
        // Every class with any Z at all is computed and plotted. minZPerClass
        // no longer gates the calculation -- it only marks a class as
        // statistics-limited, because 60-90% holds 218 / 80 / 22 Z and those
        // points are wanted despite their errors. A class below the threshold
        // is labelled on the figures rather than hidden.
        haveZ[k] = (nZ[k] > 0.);
        if(nZ[k] >= minZPerClass) nGood++;
        delete m;
      }
      fZ->Close();
    }
  }
  printf("    class      N_Z    muEff   status\n");
  for(int k = 0; k < g_nCls; k++)
    printf("    %-7s %7.0f   %.4f   %s\n", g_label[k], nZ[k], g_muEff[k],
           !haveZ[k]              ? "no Z at all" :
           nZ[k] < minZPerClass   ? "statistics-limited" : "ok");

  // ---- pp reference --------------------------------------------------------
  TH1D *ppStitch = nullptr;
  double ppFloorPt = ptFloor;
  double nZ_pp = 0.;
  {
    TFile *fppMB = TFile::Open(Form("%s/%s", repo, f_pp_MB_path));
    TFile *fpp80 = TFile::Open(Form("%s/%s", repo, f_pp_j80_path));
    TFile *fppMu = TFile::Open(Form("%s/%s", repo, f_pp_mu_path));
    if(!fppMB || fppMB->IsZombie() || !fpp80 || fpp80->IsZombie() ||
       !fppMu || fppMu->IsZombie()){
      printf("    pp inputs missing -- PbPb/pp ratio skipped\n");
    }
    else{
      TH1D *hMB = (TH1D*) fppMB->Get("h_inclRecoJetPt");
      TH1D *h80 = (TH1D*) fpp80->Get("h_inclRecoJetPt");
      TH1D *hMu = (TH1D*) fppMu->Get("h_dimuonMass");
      if(!hMB || !h80 || !hMu) printf("    a pp histogram is missing -- ratio skipped\n");
      else{
        hMB = (TH1D*) hMB->Clone("ppMB"); hMB->SetDirectory(nullptr);
        h80 = (TH1D*) h80->Clone("pp80"); h80->SetDirectory(nullptr);
        nZ_pp = hMu->Integral(hMu->FindBin(Z_lo + 0.01), hMu->FindBin(Z_hi - 0.01));

        // pp has no centrality, so one threshold and one normalisation serve
        // every class. commonK cannot be reused: it loops over the nine classes
        // and needs the PbPb per-class event counts, so the pp factor is the
        // plain yield ratio in the same common window.
        // pp Jet80 threshold.
        //
        // turnOnThreshold is NOT used here. Its criterion -- stay within nSigma
        // of the plateau all the way to the ceiling -- cannot be met by pp,
        // because pp MinBias holds only ~318 jets above 200 GeV and the
        // Jet80/MinBias ratio scatters at the top of the range; run on pp it
        // returns 480 GeV, which is the point where only one noisy bin remains
        // rather than a turn-on.
        //
        // Measured directly instead, Jet80/MinBias yield ratio by window:
        //     60-80   452      80-100  3665     100-120  6296
        //    120-150  6316    150-200  6495     200-300  6546
        // so the plateau is ~6500 and the sample reaches 97% of it by 100-120
        // and 99.9% by 150-200. 120 GeV is the first edge at or above 97%.
        printf("    pp Jet80 threshold %.0f GeV (measured turn-on; see comment)\n", ppThr80);

        std::vector<StitchSample> s;
        if(useMinBias){
          const double eps = 0.01;
          double nl = hMB->Integral(hMB->FindBin(kWinMB_lo + eps), hMB->FindBin(kWinMB_hi - eps));
          double nh = h80->Integral(h80->FindBin(kWinMB_lo + eps), h80->FindBin(kWinMB_hi - eps));
          double kpp = (nl > 0.) ? nh/nl : -1.;
          printf("    pp MinBias -> Jet80 in [%.0f,%.0f]: k = %.4g\n", kWinMB_lo, kWinMB_hi, kpp);
          s = { {"MinBias", ptFloor, hMB, kpp}, {"Jet80", ppThr80, h80, 0.} };
          ppStitch = stitchTriggerSamples(s, "ppStitch", 50., ptCeiling, true);
          printf("    pp: N_Z = %.0f, MinBias 50-%.0f, Jet80 %.0f-%.0f\n",
                 nZ_pp, ppThr80, ppThr80, ptCeiling);
        }
        else{
          s = { {"Jet80", ppThr80, h80, 0.} };
          ppStitch = stitchTriggerSamples(s, "ppStitch", ppThr80, ptCeiling, true);
          printf("    pp: N_Z = %.0f, Jet80 %.0f-%.0f  (MinBias dropped)\n",
                 nZ_pp, ppThr80, ptCeiling);
        }
        ppFloorPt = ppThr80;
      }
    }
  }

  // --------------------------------------------------- jets per Z, ratio ---
  //
  // Normalisation follows stitchSpectra_caloJets.C: each stitched spectrum is
  // divided by its own efficiency-corrected N_Z, and the ratio of the two
  // per-Z spectra is the result. Classes without enough Z are left null.
  TH1D *jpzPb[NCls], *jpzPP[NCls], *ratio[NCls];
  for(int k = 0; k < g_nCls; k++){ jpzPb[k] = jpzPP[k] = ratio[k] = nullptr; }

  const double NZ_pp_corr = (nZ_pp > 0.) ? nZ_pp/(muEff_pp*muEff_pp) : 0.;

  printf("\n  PbPb / pp jets per Z\n");
  for(int k = 0; k < g_nCls; k++){
    if(!haveZ[k] || !ppStitch || NZ_pp_corr <= 0.) continue;

    const double NZ_k = nZ[k]/(g_muEff[k]*g_muEff[k]);

    jpzPb[k] = (TH1D*) spec[k]->Clone(Form("jetsPerZ_PbPb_%s", g_tag[k]));
    jpzPb[k]->SetDirectory(nullptr);
    jpzPb[k]->Scale(1./NZ_k);

    // pp rebinned onto THIS class's edges, then treated identically
    std::vector<double> e;
    for(int b = 1; b <= spec[k]->GetNbinsX() + 1; b++) e.push_back(spec[k]->GetBinLowEdge(b));
    jpzPP[k] = rebinTo(ppStitch, (int)e.size() - 1, e.data(), Form("jetsPerZ_pp_%s", g_tag[k]));
    if(!jpzPP[k]) continue;
    jpzPP[k]->SetDirectory(nullptr);
    jpzPP[k]->Scale(1./etaRange);
    divideByBinwidth(jpzPP[k]);
    jpzPP[k]->Scale(1./NZ_pp_corr);

    ratio[k] = (TH1D*) jpzPb[k]->Clone(Form("ratio_%s", g_tag[k]));
    ratio[k]->SetDirectory(nullptr);
    ratio[k]->Divide(jpzPb[k], jpzPP[k], 1., 1., "");

    printf("    %-7s (N_Z = %.0f):", g_label[k], nZ[k]);
    for(int b = 1; b <= ratio[k]->GetNbinsX(); b++){
      // pp runs out before PbPb does: it has no Jet100 scan, so report only
      // bins where the pp reference still has entries
      if(jpzPP[k]->GetBinContent(b) <= 0.){ ratio[k]->SetBinContent(b, 0.); ratio[k]->SetBinError(b, 0.); continue; }
      printf("  %.0f-%.0f: %.3f", ratio[k]->GetBinLowEdge(b),
             ratio[k]->GetBinLowEdge(b+1), ratio[k]->GetBinContent(b));
    }
    printf("\n");
  }

  // --------------------------------------------------------------- output ---
  TFile *out = TFile::Open(outFile, "RECREATE");
  for(int k = 0; k < g_nCls; k++){
    spec[k]->Write(Form("spectrum_%s", g_tag[k]));
    stitched[k]->Write(Form("stitched_fine_%s", g_tag[k]));
    if(jpzPb[k]) jpzPb[k]->Write();
    if(jpzPP[k]) jpzPP[k]->Write();
    if(ratio[k]) ratio[k]->Write();
  }
  TH1D *hThr = new TH1D("thresholds", "Jet80 and Jet100 thresholds;class;GeV", NCls, 0, NCls);
  for(int k = 0; k < g_nCls; k++){ hThr->SetBinContent(k+1, thr80[k]); hThr->GetXaxis()->SetBinLabel(k+1, g_label[k]); }
  hThr->Write();
  out->Close();

  // ---------------------------------------------------------------- figure ---
  {
    TCanvas *c = new TCanvas("c_spec", "", 700, 800);
    c->SetLeftMargin(0.17); c->SetBottomMargin(0.12);
    c->SetTopMargin(0.17);  c->SetRightMargin(0.05);   // three caption lines
    c->SetLogy();   // the spectra span five decades across 50-500 GeV

    TLegend *leg = makeLegend(0.58, 0.40, 0.94, 0.80, 0.030);
    // Nine classes, seven usable Okabe-Ito colors (4 is the unreadable yellow),
    // so color repeats after seven and the marker shape changes with it. No
    // (color, marker) pair repeats, which a plain k%8 did not guarantee -- it
    // gave 0-10% and 80-90% the same black circle.
    const int useHex[7] = {0, 1, 2, 3, 5, 6, 7};
    const int mkCycle[2] = {markFilledCircle, markFilledSquare};
    double yMax = 0., yMin = 1e30;
    for(int k = 0; k < g_nCls; k++){
      for(int b = 1; b <= spec[k]->GetNbinsX(); b++){
        double v = spec[k]->GetBinContent(b);
        if(v > 0.){ yMax = TMath::Max(yMax, v); yMin = TMath::Min(yMin, v); }
      }
    }
    for(int k = 0; k < g_nCls; k++){
      styleH(spec[k], okabeHex[useHex[k % 7]], mkCycle[k / 7], 0.8);
      if(k == 0){
        spec[k]->SetTitle(""); spec[k]->SetStats(0);
        spec[k]->SetMinimum(yMin * 0.4); spec[k]->SetMaximum(yMax * 40.);
        spec[k]->GetXaxis()->SetTitle("jet #it{p}_{T} [GeV]");
        spec[k]->GetXaxis()->SetTitleSize(0.045); spec[k]->GetXaxis()->SetLabelSize(0.040);
        spec[k]->GetYaxis()->SetTitle("d#it{N}_{jet}/d#it{p}_{T}d#eta  [Jet100 units]");
        spec[k]->GetYaxis()->SetTitleSize(0.042); spec[k]->GetYaxis()->SetLabelSize(0.040);
        spec[k]->GetYaxis()->SetTitleOffset(1.65);
        spec[k]->Draw("E1");
      }
      else spec[k]->Draw("E1 same");
      leg->AddEntry(spec[k], g_label[k], "lp");
    }
    leg->Draw();

    TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.036);
    la.DrawLatex(0.17, 0.955, "PbPb 5.02 TeV, calo jets, anti-#it{k}_{T} #it{R} = 0.4");
    la.SetTextSize(0.030);
    la.DrawLatex(0.17, 0.917, "manual AK4Calo JEC, no unfolding, no fake jet subtraction");
    la.DrawLatex(0.17, 0.883, "MinBias + Jet80 + Jet100, 10% centrality bins");

    c->SaveAs(Form("%s/stitchedSpectra_caloJets_10pctBins.pdf", outDir));
    delete c;
  }

  // ------------------------------------------------------- ratio figures ---
  //
  // Nine classes, so the overlay needs nine distinguishable styles. Okabe-Ito
  // gives seven usable colours (index 4, yellow, is unreadable on white), so
  // the two extra classes reuse colours with a different centre-symmetric
  // marker rather than introducing an eighth hue.
  const char *okabe[8] = {"#000000","#E69F00","#56B4E9","#009E73",
                          "#F0E442","#0072B2","#D55E00","#CC79A7"};
  const int colIdx[NCls] = {0, 5, 2, 3, 6, 7, 1, 0, 5};
  const int mrkIdx[NCls] = {20,21,33,24,25,27,20,25,27};

  auto styleRatio = [&](int k){
    const int c0 = TColor::GetColor(okabe[colIdx[k]]);
    ratio[k]->SetLineColor(c0); ratio[k]->SetMarkerColor(c0);
    ratio[k]->SetMarkerStyle(mrkIdx[k]);
    ratio[k]->SetMarkerSize(1.2); ratio[k]->SetLineWidth(2);
    ratio[k]->SetStats(0);
  };

  auto makeFrame = [&](const char *nm, double yHi){
    TH1D *fr = new TH1D(nm, "", 1, ptFloor, ptCeiling);
    fr->SetStats(0);
    fr->GetYaxis()->SetRangeUser(0., yHi);
    fr->GetXaxis()->SetTitle("calo jet p_{T} [GeV]");
    fr->GetYaxis()->SetTitle("(jets per Z)_{PbPb} / (jets per Z)_{pp}");
    fr->GetXaxis()->SetTitleOffset(1.15);
    fr->GetYaxis()->SetTitleOffset(1.35);
    return fr;
  };

  // ---- all classes on one canvas -----------------------------------------
  {
    TCanvas *c = new TCanvas("c_ratio_all", "", 800, 700);
    c->SetLeftMargin(0.15); c->SetRightMargin(0.04);
    c->SetTopMargin(0.09);  c->SetBottomMargin(0.13);
    makeFrame("fr_all", 1.9)->Draw("axis");

    // Upper left: every class sits at 0.4-0.8 across the whole pT range, so the
    // right-hand side and the lower half are occupied. Only the top-left corner
    // is genuinely free.
    TLegend *leg = new TLegend(0.18, 0.55, 0.60, 0.85);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.028);
    leg->SetNColumns(2);

    for(int k = 0; k < g_nCls; k++){
      if(!ratio[k]) continue;
      styleRatio(k);
      ratio[k]->Draw("E1 same");
      leg->AddEntry(ratio[k], nZ[k] < minZPerClass
                    ? Form("%s (%.0f Z)", g_label[k], nZ[k]) : g_label[k], "lp");
    }
    leg->Draw();

    TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.034);
    la.DrawLatex(0.15, 0.955, "calo jets per Z, PbPb / pp, 10% centrality bins");

    c->SaveAs(Form("%s/jetsPerZ_caloJets_PbPbToPP_10pctBins.pdf", outDir));
    delete c;
  }

  // ---- one canvas per class, plus a single multi-page PDF ------------------
  const TString multi = Form("%s/jetsPerZ_caloJets_PbPbToPP_perClass.pdf", outDir);
  int nPages = 0, nLast = -1;
  for(int k = 0; k < g_nCls; k++) if(ratio[k]) { nPages++; nLast = k; }

  int page = 0;
  for(int k = 0; k < g_nCls; k++){
    if(!ratio[k]) continue;

    TCanvas *c = new TCanvas(Form("c_ratio_%s", g_tag[k]), "", 800, 700);
    c->SetLeftMargin(0.15); c->SetRightMargin(0.04);
    c->SetTopMargin(0.09);  c->SetBottomMargin(0.13);
    // the thin classes swing well above 1, so give them their own y range
    makeFrame(Form("fr_%s", g_tag[k]), nZ[k] < minZPerClass ? 2.6 : 1.1)->Draw("axis");

    TLine *one = new TLine(ptFloor, 1., ptCeiling, 1.);
    one->SetLineStyle(2); one->SetLineColor(kGray+2); one->Draw();

    styleRatio(k);
    ratio[k]->Draw("E1 same");

    TLatex la; la.SetNDC(); la.SetTextFont(42);
    la.SetTextSize(0.042);
    la.DrawLatex(0.19, 0.84, Form("PbPb %s", g_label[k]));
    la.SetTextSize(0.034);
    la.DrawLatex(0.15, 0.955, "calo jets per Z, PbPb / pp");
    la.SetTextSize(0.027);
    la.DrawLatex(0.19, 0.795, Form("N_{Z} = %.0f%s", nZ[k],
                 nZ[k] < minZPerClass ? "  (statistics-limited)" : ""));

    c->SaveAs(Form("%s/jetsPerZ_caloJets_PbPbToPP_%s.pdf", outDir, g_tag[k]));

    // same canvases collected into one multi-page file for flipping through
    if(nPages == 1)            c->Print(multi);
    else if(page == 0)         c->Print(multi + "(");
    else if(k == nLast)        c->Print(multi + ")");
    else                       c->Print(multi);
    page++;
    delete c;
  }

  printf("\n  results: %s\n  figures: %s/stitchedSpectra_caloJets_10pctBins.pdf\n"
         "           %s/jetsPerZ_caloJets_PbPbToPP_10pctBins.pdf  (all classes)\n"
         "           %s  (%d pages, one per class)\n"
         "           %s/jetsPerZ_caloJets_PbPbToPP_<class>.pdf  (individual)\n",
         outFile, outDir, outDir, multi.Data(), nPages, outDir);
}
