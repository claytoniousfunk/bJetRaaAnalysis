// Systematic uncertainties of the PF-jet INCLUSIVE jets per Z: PbPb with the
// calo-match jet ID (requireCaloJetMatch, "caloMatchedPF") and pp with plain
// ak4PF ("PFJets"), the pairing of the PF b-jet chain (jetCollectionInputs.h,
// bJetChain_caloJets_PbPb.C("dedicated","caloMatchedPF","PFJets")).
//
// A copy of inclusiveJetsPerZSys_caloJets.C (2026-10-07) with the PF inputs;
// the chain and every source are defined as there -- read that header. What
// differs:
//   INPUTS     the files calculateJetsPerZ.cc opens for each side:
//                PbPb  calculateJetsPerZ(false,-1,-1.,true)  calo-matched PF,
//                      2026-10-06 MinBias Part1 / Jet80 / Jet100; no Jet60
//                      (PbPb never uses it) ; response = the 2026-10-06
//                      calo-matched anaMuTag PYTHIA+HYDJET scan, ultraFine
//                      slices summed into the classes
//                pp    calculateJetsPerZ(false)  plain PF: MinBias 2026-07-06,
//                      Jet60 2026-03-10 (sibling repo), Jet80 / Jet100
//                      2026-07-27; response = plain-PF PYTHIA 2026-08-03.
//                      These pp files predate the 2026-09-08 generation change;
//                      the current-generation plain-PF HighEGJet scans were
//                      configured 2026-10-06 and are not here yet.
//              The reproduction check compares C1..C4 with the caloMatchedPF
//              file of calculateJetsPerZ.cc and pp with its PFJets file.
//   ITERATIONS min MSE of unfoldClosureTest.C(<sys>,15,"weighted","tilt:0.3",   (nIterNom below; the
//              "0-10% 8" in the MODE paragraph further down is the earlier value)
//              80.,500.,...,recoFloor 60), 2026-10-07: pp 3 (8-3 halves),
//              0-10% 9, 10-30% 6, 30-50% 4, 50-80% 3 (calo-matched response
//              merged to coarse classes, train = test: no independent half).
//              0-10% is shallow: N = 6-9 within 4% in MSE.
//   JEU        rootFiles/JEU/unfold_PYTHIA_PFJets_JEUShift.root (pp, Spring18
//              ppRef AK4PF) and unfold_PYTHIAHYDJET_caloMatchedPF_JEUShift.root
//              (PbPb, Autumn18_HI AK4PF), src/unfoldTest.
//   fakeJets   PbPb nominal rootFiles/fakeJets/fakeJets_forestCS.root since
//              2026-11-05 (the exact akCs4PF emulation on the mixed-event PF
//              pool; what calculateJetsPerZ.cc subtracts, via
//              fakeJets_leakSubFit.root), leakage-subtracted and fitted
//              (fakeJetLeakFit.h). No alternative estimate. Until 2026-11-05
//              fakeJets.root (FastJet on PF candidates, random cone) with the
//              variation fakeJets_RCgeoCorr.root. Both are ALL-PF-jet fake
//              estimates, not specific to calo-matched jets.
// MODE "MBcaloMatched" (second argument; 2026-10-07, user: "everything PF
// normal, ONLY PbPb MinBias has the calo match"): calculateJetsPerZ(false,-1,
// -1.,false,true). PbPb Jet80/Jet100 = the plain-PF HardProbes scans of
// calculateJetsPerZ.cc's PF branch (sibling repo 2026-05-05, nominal classes),
// MinBias = the calo-matched 2026-10-06 scan, response = the calo-matched
// 2026-10-06 PYTHIA+HYDJET (fRespAA) and its JEU file (fJEUAA) since 2026-10-08
// (useCaloMatchedRespMB; before: plain-PF 2026-09-21 fRespAAP / fJEUAAP); nominal N
// pp 3, 0-10% 8, 10-30% 6, 30-50% 4, 50-80% 3 (closure on that response,
// tilt:0.3, train = test). fakeJets source: envelope of NO subtraction and the
// RCgeoCorr estimate (the calo match already removes most fakes).
// Outputs inclJetsPerZSys_PF_MBcaloMatched<sfx>.root, figures likewise.
// NOT COVERED: the jet-ID mismatch itself -- PbPb carries the calo-match ID's
// genuine-jet loss and pp does not (no pp ID efficiency correction exists).
//
// Output  rootFiles/systematics/inclJetsPerZSys_caloMatchedPF_ppPFJets<sfx>.root
//         figures/systematics/inclJetsPerZSys_caloMatchedPF_ppPFJets_*<sfx>.pdf
//         (same histogram names as the calo macro)
//
// Usage: root -l -b -q -e 'gSystem->Load("/home/clayton/Programs/RooUnfold/build/libRooUnfold.so");' \
//                         'inclusiveJetsPerZSys_PF.C(true)'   // JEU correlated pp-PbPb
//        JEU not flattened (separate _noFlatJEU file):
//          'inclusiveJetsPerZSys_PF.C(true,"MBcaloMatched","","","","plain",false,true)'
//        + the variation spectra (inclJetsPerZVariations_..._noFlatJEU.root): 9th argument true
//        MinBias only, no stitching (_MBonly): 10th argument true
//        fakeJets from the RCgeoCorr estimate only (_fakeGeoOnly): 11th argument true
//        event-count trigger stitch (_evtCountStitch, header EVENT-COUNT STITCH): 12th argument true;
//          compare with plotInclRZ_compareStitch_PF.C after the deck's dropSysSources passes
// Run from: src/calculateJetsPerZ/

#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#endif
#include "TFile.h"
#include "../../headers/plotting/flattenSystematic.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLine.h"
#include "TBox.h"
#include "TColor.h"
#include "TLegend.h"
#include "TSystem.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TVectorD.h"
#include <vector>
#include "../../headers/plotting/coarseCent.h"
#include "../../headers/plotting/plotStyle.h"
#include "../../headers/functions/caloJetBPurityFitCore.h"   // ptEdges (b-jet axis), ptReportMin
#include "../../headers/functions/jeuShift.h"               // JEU on the data-shaped spectrum (analytic shift)
#include "TGraph.h"
#include "../../headers/functions/dataDrivenPrior.h"        // tiltGenRows, fitPriorTilt (unfoldingBias)
#include "../../headers/functions/fakeJetLeakFit.h"          // fitted, leak-subtracted fake fraction
#include "../../headers/functions/jerSmear.h"                // JER scheme: nominal response x1.10, variation x1.10 (HIN-18-014)
#include "../../headers/publishedData/CMS_ZYield_5p02TeV.h"  // BJET_NZ_SOURCE=published: N_Z from the published yields (pubZ)

namespace ijzpf {

const TString scanDir = "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/";
// the files calculateJetsPerZ.cc opens: PbPb caloMatchedPF branch, pp PF branch
const TString fPPMB    = scanDir + "pp/pp_MinBias_mu12_pTmu-15to999_tight_deltaR-40_jetTrkMaxFilter_WDecayFilter_2026-7-6.root";
const TString fPP60    = "/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/scanningOutput/pp/latest/pp_HighEGJet_Jet60HLT_mu12_pTmu-15to999_tight_deltaR-40_jetTrkMaxFilter_2026-3-10.root";
const TString fPP80    = scanDir + "pp/pp_HighEGJet_Jet80HLT_mu12_pTmu-15to999_tight_deltaR-40_jetTrkMaxFilter_WDecayFilter_2026-7-27.root";
const TString fPP100   = scanDir + "pp/pp_HighEGJet_Jet100HLT_mu12_pTmu-15to999_tight_deltaR-40_jetTrkMaxFilter_WDecayFilter_2026-7-27.root";
const TString fPPZ     = "/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/scanningOutput/pp/latest/pp_SingleMuon_mu12_pTmu-15to999_tight_deltaR-40_mu12TriggerEfficiencyCorrection_jetTrkMaxFilter_WDecayFilter_2026-5-4.root";
const TString fAAMB    = scanDir + "PbPb/PbPb_MinBias_Part1_manualJEC_caloJetMatched_mu12_pTmu-15to999_tight_mu12TriggerEfficiencyCorrection_jetTrkMaxFilter_WDecayFilter_2026-10-6_ultraFineCentBins.root";
const TString fAA80    = scanDir + "PbPb/PbPb_HardProbes_manualJEC_caloJetMatched_Jet80HLT_mu12_pTmu-15to999_tight_mu12TriggerEfficiencyCorrection_jetTrkMaxFilter_WDecayFilter_2026-10-6_ultraFineCentBins.root";
const TString fAA100   = scanDir + "PbPb/PbPb_HardProbes_manualJEC_caloJetMatched_Jet100HLT_mu12_pTmu-15to999_tight_mu12TriggerEfficiencyCorrection_jetTrkMaxFilter_WDecayFilter_2026-10-6_ultraFineCentBins.root";
const TString fAAZ     = scanDir + "PbPb/PbPb_SingleMuon_manualJEC_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_2026-9-24_ultraFineCentBins.root";
const TString fRespPP  = scanDir + "PYTHIA/PYTHIA_DiJet_response_pThat-15_mu12_pTmu-15_tight_jetTrkMaxFilter_doPThatCorrelationFilterTight_2026-8-3.root";
const TString fRespAA  = scanDir + "PYTHIAHYDJET/PYTHIAHYDJET_response_DiJet_manualJEC_caloJetMatched_anaMuTag_pThat-15_mu12_pTmu-15_tight_vzReweight_hiBinReweight_hiBinShift-10_jetTrkMaxFilter_doPThatCorrelationFilterTight_2026-10-6_ultraFineCentBins.root";
// Combinatorial-jet estimate since 2026-11-05 (user: the random cone is deprecated): the EXACT
// akCs4PF emulation run on the mixed-event PF pool (doForestCSJets, 2026-10-29 scan), per event per
// class (plotFakeJetFraction_forestCS_PF.C). It has no alternative estimate (useFakeAlt false); the
// random-cone files (fakeJets.root nominal, fakeJets_RCgeoCorr.root variation) were used until then.
const TString fFakeNom = "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/fakeJets/fakeJets_forestCS.root";
const TString fFakeAlt = "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/fakeJets/fakeJets_RCgeoCorr.root";
bool useFakeAlt = false;
const TString fJEUPP   = "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/JEU/unfold_PYTHIA_PFJets_JEUShift.root";
const TString fJEUAA   = "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/JEU/unfold_PYTHIAHYDJET_caloMatchedPF_JEUShift.root";
const TString fCheckAA = "/home/clayton/Analysis/code/bJetRaaAnalysis/src/calculateJetsPerZ/rootFiles/JetsPerZ/jetsPerZSpectra_PFJets_caloMatchedPF_coarseMuEff_unfold1iter.root";
const TString fCheckPP = "/home/clayton/Analysis/code/bJetRaaAnalysis/src/calculateJetsPerZ/rootFiles/JetsPerZ/jetsPerZSpectra_PFJets_coarseMuEff_unfold1iter.root";
const char *collText = "PF jets (PbPb calo-matched)";
// MBcaloMatched mode
const TString fAA80P   = "/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/scanningOutput/PbPb/latest/PbPb_HardProbes_Jet80HLT_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_2026-5-5.root";
const TString fAA100P  = "/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/scanningOutput/PbPb/latest/PbPb_HardProbes_Jet100HLT_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_2026-5-5.root";
const TString fRespAAP = scanDir + "PYTHIAHYDJET/PYTHIAHYDJET_response_DiJet_manualJEC_pThat-15_mu12_pTmu-15_tight_vzReweight_hiBinReweight_hiBinShift-10_jetTrkMaxFilter_doPThatCorrelationFilterTight_2026-9-21.root";
const TString fJEUAAP  = "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/JEU/unfold_PYTHIAHYDJET_PFJets_JEUShift.root";
const TString fCheckMB = "/home/clayton/Analysis/code/bJetRaaAnalysis/src/calculateJetsPerZ/rootFiles/JetsPerZ/jetsPerZSpectra_PFJets_MBcaloMatched_coarseMuEff_unfold1iter.root";
bool modeMB = false;
// Per-class turn-on stitch edges, the same as calculateJetsPerZ.cc's
// turnOnJ80Lo / turnOnJ100Lo (deriveTriggerThresholds_PFJets.C, 2026-10-08),
// MBcaloMatched mode only. Index: pp, 0-10, 10-30, 30-50, 50-80%. Opts keeps its
// fixed-edge defaults (130/170/200); in this mode they act as OFFSETS from those
// defaults, so the stitchWindow variations (140/180, 190, 220) become +10 GeV on
// the MinBias window, +20 on its upper edge, +20 on the Jet100 edge.
bool turnOnStitch = true;
// PbPb response in MBcaloMatched mode: the calo-matched one (fRespAA, fJEUAA),
// 2026-10-08 -- the plain-PF one biased central R_Z ~6-7% high at 80-100 GeV.
bool useCaloMatchedRespMB = true;
// MinBias only (last argument, 2026-10-20, user: "only use the MinBias spectra, no
// trigger stitching"): the scaled MinBias spectrum at every pT, no Jet60/80/100
// pieces. The MinBias scale factor kMB is kept as it is, since N_Z is counted
// in SingleMuon and only the jet-triggered chain ties MinBias to that
// luminosity; so the stitching then sets ONE number per class (flat in pT) and
// stitchWindow / stitchNormStat become flat normalization sources. No source is
// flattened in this mode (flattenSources). _MBonly in the name.
bool mbOnly = false;
bool fakeGeoOnlyMode = false;          // 11th argument: fakeJets without the "no subtraction" variation
// EVENT-COUNT STITCH (12th argument, 2026-10-21, user; studyStitchByEventCount_PF.C): each
// trigger sample normalized to Jet100-sample units by trigger-bit event counts instead of
// spectrum matching:
//   k_MB = N_evt(HardProbes/HighEGJet Jet100) / N_evt(MinBias with the Jet100 bit)  (h_vz_jet100)
//   k_80 = N_evt(Jet100) / N_evt(Jet80 & Jet100)  (h_leadJetPt_jet100 / _jet80_and_jet100; PbPb
//          from the 2026-10-06 ultraFine HardProbes scans, the only ones that book it -- the
//          fraction is 0.463-0.464 in every class; pp from the Jet100 file, Jet80 unprescaled)
// No Jet60 leg. Switch points (MinBias -> Jet80, Jet80 -> Jet100) where the study finds each
// ratio flat (2 sigma + 1%), at least the trigger turn-on edge, and in 0-10% at least
// 170 / 180 GeV. stitchWindow then moves these switch points (+10 GeV MinBias edge, +20
// Jet100 edge, as the window variations); stitchNormStat = the counts' statistics.
// _evtCountStitch in the name. Only the "" centrality path.
bool evtCountStitch = false;
double cntJ80Lo[5]  = {115., 170., 145., 130., 120.};
double cntJ100Lo[5] = {145., 180., 155., 150., 140.};
double cntKMB[5], cntKMBrel[5], cntK80[5], cntK80rel[5];   // k_MB per MinBias EVENT-normalized spectrum
double tsJ80Lo[5]  = {115., 140., 145., 115., 120.};
double tsJ100Lo[5] = {140., 155., 155., 140., 140.};

TString outStem = "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/systematics/inclJetsPerZSys_caloMatchedPF_ppPFJets";
const TString figDir  = "/home/clayton/Analysis/code/bJetRaaAnalysis/figures/systematics";

const int NS = 5;                                     // pp, C1..C4
const char *sysTag[NS]   = {"pp", "C1", "C2", "C3", "C4"};
const char *sysLabel[NS] = {"pp", "PbPb 0-10%", "PbPb 10-30%", "PbPb 30-50%", "PbPb 50-80%"};

// CENTRALITY CLASSES (fifth argument, 2026-10-09). "" = the code path above:
// nominal classes, Jet80/Jet100 read by class index in MBcaloMatched mode.
// "slices" and "ATLAS" build every PbPb input from the 5% slices instead:
// jets, response, both fake estimates (makeFakeJetFile.C / _RCgeoCorr.C recipe
// on their mixed-event scan), the JEU reco spectra (unfoldPYTHIAHYDJETJetPt_JEUShift.C
// recipe on its scan) and N_Z, each slice with the efficiency (and error) of the
// nominal class it lies in.
//   "slices"  nominal classes; with trig "plain" it must reproduce the "" output
//             exactly (validation of the slice path)
//   "ATLAS"   0-20 / 20-50 / 50-80% (arXiv:2204.13530) in C1..C3, 0-80% in C4
//             (so "R_CP" here is class / 0-80%). Stitch edges = the larger of the
//             overlapping nominal classes' (not re-derived); iterations from a
//             closure test per class (see the class tables in the main function). The b-jet axis is replaced by the pT edges of ATLAS PLB 790
//             (79/89/100/112/125/141/158/177/199/223/251/281/316/398) moved to the
//             5 GeV grid, from 60 to 500 GeV ("atlas"; plotInclRZ_ATLASCent_vsATLAS_PF.C).
//   "10pctA" / "10pctB"  the 10% classes 0-10 .. 30-40% / 40-50 .. 70-80% (four
//             slots each, so two runs), same axis, stitch edges of the nominal class
//             each lies in, per-class closure N_iter; the direct check against
//             ATLAS PLB 790's own classes.
// trig (sixth argument), PbPb Jet80/Jet100 in these modes: "plain" = the
// 2026-05-05 nominal files (nominal classes only; they use the pre-2026-08-12
// hiBin convention), "caloMatched" = the 2026-10-06 ultraFine calo-matched scans,
// "plainUF" = a plain-PF ultraFine rescan (paths below, not scanned yet).
// bjetAxis (seventh argument, 2026-10-15): in the slice modes keep the b-jet axis
// (ptEdges, bJetChain_caloJets_PbPb.C's edges) instead of "atlas", _bjetAxis in the
// name; the b-jet chain's cent "ATLAS" mode reads this file. Every output also
// carries NZ (TVectorD, pp and C1..C4, efficiency corrected) since 2026-10-15.
TString centMode = "";
int clsLo[4] = {1, 3, 7, 11}, clsHi[4] = {2, 6, 10, 16};
// slice modes: the forest-CS mixed-event scan since 2026-11-05 (as fFakeNom); the random-cone
// 2026-8-11 scan (h_fastJetPt_PF_JEC_bkgSub_RC[_geoCorr]) before. [1] is used only with useFakeAlt.
const TString fMixedFake = "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/PbPb_MinBias_Part1_mu12_pTmu-15to999_tight_mu12TriggerEfficiencyCorrection_jetTrkMaxFilter_WDecayFilter_mixedEventPFClustering_fastJetResamples-100_pseudoJetCandPtMin-0.0_caloConstituents_bkgMapOverride_2026-10-29_forestCS_ultraFineCentBins.root";
const char *fakeBase[2] = {"h_csJetPt_JEC", "h_csJetPt_JEC"};   // fFakeNom, fFakeAlt
const TString fJEUScanAA = "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PYTHIAHYDJET/PYTHIAHYDJET_manualJEC_caloJetMatched_DiJet_pThat-15_mu12_pTmu-15to999_tight_vzReweight_hiBinShift-10_jetTrkMaxFilter_doPThatCorrelationFilterTight_removeUnmatchedJets_WDecayFilter_genMuTagFix_2026-10-6_ultraFineCentBins.root";   // the "inputs" of fJEUAA
const double jeuScanFloor = 50.;   // unfoldPYTHIAHYDJETJetPt_JEUShift.C's recoFloor (below the chain's own)
const TString fAA80PlainUF = "", fAA100PlainUF = "";   // plain-PF ultraFine rescan, pending
TString trigAA = "plain";
int        nIterNom[NS]  = {3, 9, 6, 4, 3};      // tilt:0.3 optimum (header); nIterOpt overrides

const double muEff[NS]    = {0.9708, 0.8627, 0.9069, 0.9856, 0.9778};   // calculateJetsPerZ.cc (coarse)
const double muEffErr[NS] = {0.0077, 0.0350, 0.0347, 0.0151, 0.0105};   // delta eps: |SF_ID - 1| (+) sigma_SF (+) T&P bkg (+) sim. N_eff stat., studyMuonEffSF_binned.C (2026-10-21)

// JER (2026-11-05, user): HIN-18-014 scheme, headers/functions/jerSmear.h -- every response is smeared
// to a 10% worse resolution right after loading (loadInputs), the variation is a further 10%;
// BJET_JER_SCHEME=old = the earlier unsmeared nominal with a +20% variation.
const int    nToys    = 500;
const int    toySeed  = 20261005;
bool   flattenSources = !(std::getenv("BJET_SYS_NOFLATTEN") && TString(std::getenv("BJET_SYS_NOFLATTEN")) == "1");   // sources without a pT trend -> their mean (flattenSystematic.h); off with mbOnly,
                                      // where the statistically empty high-pT MinBias bins would set every source's mean
// JEU CROSS-CHECK on the data-shaped spectrum (2026-10-21): true = the measured (stitched,
// fake-subtracted) reco spectrum of each system is shifted analytically by +/- delta(pT), the
// AK4PF JEC uncertainty averaged uniformly over |eta| < 1.6 (pp: Spring18 ppRef, PbPb: Autumn18
// HI), and unfolded with the nominal response and N. It reproduces the nominal (false: the
// simulated reco spectrum shifted jet by jet in the scan) within ~2-5% absolute on the spectra:
// the measured central reco spectrum is as steep as the simulated one (n = 6-7 at 70-100 GeV in
// 0-10%, steepened by the resolution), so the large central low-pT JEU is real. Not nominal: the
// measured spectrum's stitching seams and fluctuations make it noisier, and the analytic shift
// uses an eta-averaged delta. Closure of the analytic shift on the simulated spectrum: 1-2%.
// Output of a true run: rootFiles/systematics/jeuDataShape_crossCheck_2026-10-21/.
bool jeuDataShape = false;
std::vector<TH1D*> jeuVarOut;   // signed JEU shifts per variation, written to the output (2026-10-21)
const char *jeuUncPP = "/home/clayton/Analysis/code/bJetRaaAnalysis/JetEnergyCorrections/Spring18_ppRef5TeV_V6_MC_Uncertainty_AK4PF.txt";
const char *jeuUncAA = "/home/clayton/Analysis/code/bJetRaaAnalysis/JetEnergyCorrections/Autumn18_HI_V8_MC_Uncertainty_AK4PF.txt";
bool noFlattenJEU = false;            // last argument: JEU kept bin by bin (2026-10-20, user: the F-test against a
                                      // line in log pT flattened the 0-10% unfolded JEU, which dips 42 -> 11 -> 23%)

const double atlasPtEdges[17] = {60., 70., 80., 90., 100., 110., 125., 140., 160., 175., 200., 225., 250., 280., 315., 400., 500.};   // centMode "ATLAS": replaces the bjet axis
const int NInclE = 14;
const double inclEdges[NInclE] = {60,70,80,90,100,110,130,150,180,200,240,280,350,500};
struct Axis { const char *tag; int n; const double *edges; };
const int NAx = 2;
Axis axes[NAx] = { {"incl", NInclE - 1, inclEdges}, {"bjet", NPt, ptEdges} };

double nomFloor = 60.;   // set to 70 in MBcaloMatched mode, as calculateJetsPerZ.cc's nominalFloor()
// Real-jet leakage of the mixed-event fake estimate (2026-10-21, user). The estimate carries a flat
// 3% x MinBias of genuine-jet leading fragments drawn into the mixed event (makeFakeJetFile_leakSub.C
// header), and where it is mostly leakage it has slice-to-slice statistical spikes. The subtracted
// rate is the FIT of headers/functions/fakeJetLeakFit.h: (fake - fakeLeak x MinBias)/MinBias fitted
// with exp(a + b pT) over 60-140 GeV per class, times MinBias -- identical to fakeJets_leakSubFit.root
// for the coarse classes (makeFakeJetFile_leakSubFit.C) and done the same way for the slice-mode
// estimates. fakeJets systematic: envelope of the fit's +-1 sigma band, the geometry-corrected
// estimate (fitted the same way), and fakeLeak 0.020 / 0.040 (per-class plateau fits 0.019-0.042).
// fakeLeakNom = 0 restores the pre-2026-10-21 subtraction of the full, unfitted estimate.
// forest-CS (2026-11-05): leakage 0.026 = the flat fake / MinBias plateau of the forest-CS estimate,
// 10-80%, 60-200 GeV in 10 GeV bins (0.0256 +- 0.0003; the random cone gave 0.0275 the same way);
// varied by +-0.010 as before. Random cone: 0.030.
// 2026-11-05 (user): the leakage argument is withdrawn for now (under discussion): nominal is the
// FULL, UNFITTED forest-CS estimate, fakeLeakNom = 0, so no leakage / fit-band variations and,
// with no alternative estimate, no fakeJets systematic in the fakeGeoOnly runs. 0.026 = the
// leak-removed, fitted subtraction of 2026-11-05 (measureFakeLeakage_forestCS_PF.C).
double fakeLeakNom = 0.;
struct Opts {
  int    fakeAlt = 0;   // 0 nominal file, 1 alternative file, 2 no subtraction
  double leak = fakeLeakNom;   // fraction of the MinBias spectrum removed from the fake estimate as leakage
  int    fakeFitVar = 0;       // +1 / -1: the fit's 1 sigma band
  double mbLo = 130., mbHi = 170., j100Lo = 200.;
  double floor = nomFloor;   // 60; 70 in MBcaloMatched mode (2026-10-08)
  double zLo = 60., zHi = 120.;
};

// ---------------------------------------------------------------- inputs ------
struct Raw { TH1D *mb = nullptr, *j60 = nullptr, *j80 = nullptr, *j100 = nullptr, *fake[2] = {nullptr, nullptr}; };
Raw raw[NS];
// MinBias scale factor of the LAST stitch() call per system (read-only diagnostic, 2026-11-09:
// for PbPb = the MinBias-equivalent event count of the Jet100 sample; inclusiveJetRAA_PF.C)
double lastKMB[NS] = {0., 0., 0., 0., 0.};
TH2D *resp[NS];
TH1D *dimu[NS][16];                                   // pp: [0][0]; PbPb: 5% slices of the class
TH1D *dimuSS[NS][16];                                 // same sign, same layout
double dimuEff[NS][16], dimuEffErr[NS][16];           // muon efficiency (error) per entry
int   nDimu[NS];
int   nzSlLo[NS] = {0}, nzSlHi[NS] = {0};              // PbPb: the 5% slices (1-based) of each system, for the published N_Z
TH1D *jeuReco[NS][3];                                 // nominal, up, down

TH1D* own(TH1D *h, const char *name){ if(!h) return nullptr; TH1D *c = (TH1D*) h->Clone(name); c->SetDirectory(nullptr); return c; }

// calculateJetsPerZ.cc's getCentHist: 5% slices summed for ultraFine files, _C<c> otherwise
TH1D* centHist(TFile *f, const char *base, int c, const char *tag){
  TH1D *h = nullptr;
  if(TString(f->GetName()).Contains("ultraFineCentBins")) h = coarseSum1(f, base, c - 1, tag);
  else { f->GetObject(Form("%s_C%d", base, c), h); h = own(h, Form("%s_%s_C%d", base, tag, c)); }
  if(!h) printf("ERROR: %s class %d missing in %s\n", base, c, f->GetName());
  return h;
}

// ---- slice path (centMode "slices" / "ATLAS") ----------------------------------
TH1D* sliceSum1(TFile *f, const char *base, int lo, int hi, const char *name){
  TH1D *s = nullptr;
  for(int si = lo; si <= hi; si++){
    TH1D *h = nullptr; f->GetObject(Form("%s_C%d", base, si), h);
    if(!h){ printf("ERROR: %s_C%d missing in %s\n", base, si, f->GetName()); delete s; return nullptr; }
    if(!s){ s = (TH1D*) h->Clone(name); s->SetDirectory(nullptr); } else s->Add(h);
  }
  return s;
}
TH2D* sliceSum2(TFile *f, const char *base, int lo, int hi, const char *name){
  TH2D *s = nullptr;
  for(int si = lo; si <= hi; si++){
    TH2D *h = nullptr; f->GetObject(Form("%s_C%d", base, si), h);
    if(!h){ printf("ERROR: %s_C%d missing in %s\n", base, si, f->GetName()); delete s; return nullptr; }
    if(!s){ s = (TH2D*) h->Clone(name); s->SetDirectory(nullptr); } else s->Add(h);
  }
  return s;
}
// makeFakeJetFile.C: summed FastJet spectrum / summed events, < 20 GeV and negative bins zeroed
TH1D* sliceFakeRate(TFile *fM, const char *base, int lo, int hi, const char *name){
  TH1D *fj = sliceSum1(fM, base, lo, hi, name), *vz = sliceSum1(fM, "h_vz", lo, hi, Form("%s_vz", name));
  if(!fj || !vz) return nullptr;
  fj->Scale(1./vz->Integral());
  for(int b = 1; b <= fj->GetNbinsX(); b++)
    if(fj->GetXaxis()->GetBinUpEdge(b) <= 20. || fj->GetBinContent(b) < 0.){ fj->SetBinContent(b, 0.); fj->SetBinError(b, 0.); }
  delete vz;
  return fj;
}
int nominalClassOf(int si){ for(int c = 0; c < 4; c++) if(si >= coarseSliceLo[c] && si <= coarseSliceHi[c]) return c + 1; return -1; }

bool loadSliceClasses(){
  const TString p80  = trigAA == "caloMatched" ? fAA80  : trigAA == "plain" ? fAA80P  : fAA80PlainUF;
  const TString p100 = trigAA == "caloMatched" ? fAA100 : trigAA == "plain" ? fAA100P : fAA100PlainUF;
  TFile *aMB = TFile::Open(fAAMB), *aZ = TFile::Open(fAAZ), *rA = TFile::Open(fRespAA), *fM = TFile::Open(fMixedFake), *fJ = TFile::Open(fJEUScanAA);
  TFile *a80 = TFile::Open(p80), *a100 = TFile::Open(p100);
  for(TFile *f : {aMB, aZ, rA, fM, fJ, a80, a100}) if(!f || f->IsZombie()){ printf("ERROR: a slice-path input did not open\n"); return false; }
  const bool ufTrig = p80.Contains("ultraFineCentBins") && p100.Contains("ultraFineCentBins");
  if(!ufTrig && centMode != "slices"){ printf("ERROR: %s classes need ultraFine Jet80/Jet100 files\n", centMode.Data()); return false; }
  for(int c = 1; c <= 4; c++){
    const int lo = clsLo[c-1], hi = clsHi[c-1];
    nzSlLo[c] = lo; nzSlHi[c] = hi;
    Raw &R = raw[c];
    R.mb = sliceSum1(aMB, "h_inclRecoJetPt", lo, hi, Form("sl_mb_%d", c));
    TH1D *vz = sliceSum1(aMB, "h_vz", lo, hi, Form("sl_vz_%d", c));
    if(!R.mb || !vz) return false;
    R.mb->Scale(1./vz->Integral()); delete vz;
    R.j60 = nullptr;
    if(ufTrig){
      R.j80  = sliceSum1(a80,  "h_inclRecoJetPt", lo, hi, Form("sl_j80_%d", c));
      R.j100 = sliceSum1(a100, "h_inclRecoJetPt", lo, hi, Form("sl_j100_%d", c));
    } else {
      R.j80  = centHist(a80,  "h_inclRecoJetPt", c, "sl80");
      R.j100 = centHist(a100, "h_inclRecoJetPt", c, "sl100");
    }
    R.fake[0] = sliceFakeRate(fM, fakeBase[0], lo, hi, Form("sl_fake0_%d", c));
    R.fake[1] = useFakeAlt ? sliceFakeRate(fM, fakeBase[1], lo, hi, Form("sl_fake1_%d", c)) : R.fake[0];
    delete resp[c];
    resp[c] = sliceSum2(rA, "h_matchedRecoJetPt_genJetPt_allJets", lo, hi, Form("sl_resp_%d", c));
    const char *vn[3] = {"nominal", "JEUShiftUp", "JEUShiftDown"};
    for(int k = 0; k < 3; k++){
      TH2D *J = sliceSum2(fJ, Form("h_inclRecoJetPt_flavor_%s", vn[k]), lo, hi, Form("sl_jeu2_%d_%d", k, c));
      if(!J) return false;
      TH1D *h = J->ProjectionX(Form("sl_jeu_%d_%d", k, c), 1, J->GetNbinsY()); h->SetDirectory(nullptr);
      for(int b = 0; b <= h->GetNbinsX() + 1; b++) if(h->GetXaxis()->GetBinCenter(b) < jeuScanFloor){ h->SetBinContent(b, 0.); h->SetBinError(b, 0.); }
      jeuReco[c][k] = h; delete J;
    }
    nDimu[c] = 0;
    TH1D *t = nullptr;
    for(int si = lo; si <= hi; si++){
      const int nc = nominalClassOf(si);
      dimuEff[c][nDimu[c]] = muEff[nc]; dimuEffErr[c][nDimu[c]] = muEffErr[nc];
      aZ->GetObject(Form("h_dimuonMass_sameSign_C%d", si), t); dimuSS[c][nDimu[c]] = own(t, Form("sl_dimuSS_%d_%d", c, si));
      aZ->GetObject(Form("h_dimuonMass_C%d", si), t); dimu[c][nDimu[c]++] = own(t, Form("sl_dimu_%d_%d", c, si));
    }
    if(!R.j80 || !R.j100 || !R.fake[0] || !R.fake[1] || !resp[c]) return false;
    printf("%-12s slices %2d-%2d  stitch MB->Jet80 %g, ->Jet100 %g  N_iter %d\n", sysLabel[c], lo, hi, tsJ80Lo[c], tsJ100Lo[c], nIterNom[c]);
  }
  printf("PbPb Jet80 / Jet100: %s / %s\n", gSystem->BaseName(p80.Data()), gSystem->BaseName(p100.Data()));
  return true;
}

bool loadInputs(){
  TFile *pMB = TFile::Open(fPPMB), *p60 = TFile::Open(fPP60), *p80 = TFile::Open(fPP80), *p100 = TFile::Open(fPP100), *pZ = TFile::Open(fPPZ);
  TFile *aMB = TFile::Open(fAAMB), *a80 = TFile::Open(modeMB ? fAA80P : fAA80), *a100 = TFile::Open(modeMB ? fAA100P : fAA100), *aZ = TFile::Open(fAAZ);
  // MBcaloMatched: calo-matched PbPb response and its JEU file since 2026-10-08
  // (plain-PF fRespAAP / fJEUAAP before; set useCaloMatchedRespMB = false for those)
  const bool plainAA = modeMB && !useCaloMatchedRespMB;
  TFile *rP = TFile::Open(fRespPP), *rA = TFile::Open(plainAA ? fRespAAP : fRespAA), *fk0 = TFile::Open(fFakeNom), *fk1 = TFile::Open(fFakeAlt);
  TFile *jP = TFile::Open(fJEUPP), *jA = TFile::Open(plainAA ? fJEUAAP : fJEUAA);
  for(TFile *f : {pMB, p60, p80, p100, pZ, aMB, a80, a100, aZ, rP, rA, fk0, fk1, jP, jA})
    if(!f || f->IsZombie()){ printf("ERROR: cannot open an input file\n"); return false; }

  TH1D *t = nullptr;
  pMB->GetObject("h_inclRecoJetPt", t);  raw[0].mb   = own(t, "pp_mb");
  p60->GetObject("h_inclRecoJetPt", t);  raw[0].j60  = own(t, "pp_j60");
  p80->GetObject("h_inclRecoJetPt", t);  raw[0].j80  = own(t, "pp_j80");
  p100->GetObject("h_inclRecoJetPt", t); raw[0].j100 = own(t, "pp_j100");
  TH2D *r2 = nullptr;
  rP->GetObject("h_matchedRecoJetPt_genJetPt_allJets", r2);
  resp[0] = r2 ? (TH2D*) r2->Clone("resp_pp") : nullptr;
  pZ->GetObject("h_dimuonMass", t); dimu[0][0] = own(t, "dimu_pp"); nDimu[0] = 1;
  dimuEff[0][0] = muEff[0]; dimuEffErr[0][0] = muEffErr[0];
  pZ->GetObject("h_dimuonMass_sameSign", t); dimuSS[0][0] = own(t, "dimuSS_pp");
  jP->GetObject("reco_nominal", t); jeuReco[0][0] = own(t, "jeu_pp_0");
  jP->GetObject("reco_up", t);      jeuReco[0][1] = own(t, "jeu_pp_1");
  jP->GetObject("reco_down", t);    jeuReco[0][2] = own(t, "jeu_pp_2");

  for(int c = 1; c <= 4; c++){
    Raw &R = raw[c];
    R.mb   = centHist(aMB, "h_inclRecoJetPt", c, "aaMB");
    TH1D *vz = centHist(aMB, "h_vz", c, "aaVz");
    if(!R.mb || !vz) return false;
    R.mb->Scale(1./vz->Integral());                   // per MinBias event, like the fakes
    R.j60  = nullptr;   // no calo-matched Jet60; PbPb is stitched MinBias -> Jet80 -> Jet100
    R.j80  = centHist(a80,  "h_inclRecoJetPt", c, "aa80");
    R.j100 = centHist(a100, "h_inclRecoJetPt", c, "aa100");
    fk0->GetObject(Form("h_fakeJets_C%d", c), t); R.fake[0] = own(t, Form("fake0_C%d", c));
    if(useFakeAlt){ fk1->GetObject(Form("h_fakeJets_C%d", c), t); R.fake[1] = own(t, Form("fake1_C%d", c)); }
    else R.fake[1] = R.fake[0];   // no alternative estimate (forest CS)
    if(TString(rA->GetName()).Contains("ultraFineCentBins")) resp[c] = coarseSum2(rA, "h_matchedRecoJetPt_genJetPt_allJets", c - 1, "respAA");   // ultraFine slices -> class
    else { rA->GetObject(Form("h_matchedRecoJetPt_genJetPt_allJets_C%d", c), r2); resp[c] = r2 ? (TH2D*) r2->Clone(Form("resp_C%d", c)) : nullptr; }
    nDimu[c] = 0;
    nzSlLo[c] = coarseSliceLo[c-1]; nzSlHi[c] = coarseSliceHi[c-1];
    for(int s = coarseSliceLo[c-1]; s <= coarseSliceHi[c-1]; s++){
      aZ->GetObject(Form("h_dimuonMass_sameSign_C%d", s), t); dimuSS[c][nDimu[c]] = own(t, Form("dimuSS_C%d_%d", c, s));
      dimuEff[c][nDimu[c]] = muEff[c]; dimuEffErr[c][nDimu[c]] = muEffErr[c];
      aZ->GetObject(Form("h_dimuonMass_C%d", s), t); dimu[c][nDimu[c]++] = own(t, Form("dimu_C%d_%d", c, s));
    }
    for(int k = 0; k < 3; k++){
      const char *nm[3] = {"reco_nominal", "reco_up", "reco_down"};
      jA->GetObject(Form("%s_C%d", nm[k], c), t); jeuReco[c][k] = own(t, Form("jeu_C%d_%d", c, k));
    }
  }
  if(evtCountStitch){
    TFile *u80 = TFile::Open(fAA80), *u100 = TFile::Open(fAA100);   // 2026-10-06 ultraFine: the Jet80 & Jet100 counts
    if(!u80 || !u100){ printf("ERROR: ultraFine HardProbes files for the Jet80 & Jet100 counts\n"); return false; }
    auto I = [](TFile *f, const char *n){ TH1D *h = nullptr; f->GetObject(n, h); return h ? h->Integral(0, h->GetNbinsX()+1) : -1.; };
    for(int s = 0; s < NS; s++){
      double nMB, nMBJ100, nJ100, nJ100lead, nAnd;
      if(s == 0){
        nMB = I(pMB, "h_vz"); nMBJ100 = I(pMB, "h_vz_jet100"); nJ100 = I(p100, "h_vz_jet100");
        nJ100lead = I(p100, "h_leadJetPt_jet100"); nAnd = I(p100, "h_leadJetPt_jet80");   // Jet100-gated file
      } else {
        TH1D *v = centHist(aMB, "h_vz", s, "cntVz"), *vj = centHist(aMB, "h_vz_jet100", s, "cntVzJ"), *v1 = centHist(a100, "h_vz_jet100", s, "cntV100");
        TH1D *l1 = coarseSum1(u100, "h_leadJetPt_jet100", s - 1, "cntL1"), *la = coarseSum1(u80, "h_leadJetPt_jet80_and_jet100", s - 1, "cntLA");
        if(!v || !vj || !v1 || !l1 || !la){ printf("ERROR: event-count histograms missing for %s\n", sysTag[s]); return false; }
        nMB = v->Integral(0, -1); nMBJ100 = vj->Integral(0, -1); nJ100 = v1->Integral(0, -1);
        nJ100lead = l1->Integral(0, -1); nAnd = la->Integral(0, -1);
        delete v; delete vj; delete v1; delete l1; delete la;
      }
      const double fMB = nMBJ100/nJ100, f80 = nAnd/nJ100lead;   // binomial: the numerator events are a subset
      // raw[s].mb is per MinBias event for PbPb, raw counts for pp (loadInputs above)
      cntKMB[s] = (s == 0 ? 1. : nMB)/fMB;  cntKMBrel[s] = sqrt((1. - fMB)/nMBJ100);
      cntK80[s] = 1./f80;   cntK80rel[s] = sqrt((1. - f80)/nAnd);
      printf("event-count stitch %-12s k_MB %.4g per raw MinBias count +- %.1f%%   k_80 %.4f +- %.2f%%   switches %g / %g GeV\n",
             sysLabel[s], 1./fMB, 100.*cntKMBrel[s], cntK80[s], 100.*cntK80rel[s], cntJ80Lo[s], cntJ100Lo[s]);
    }
  }
  if(centMode != "" && !loadSliceClasses()) return false;
  for(int s = 0; s < NS; s++){
    if(!raw[s].mb || (s == 0 && !raw[s].j60) || !raw[s].j80 || !raw[s].j100 || !resp[s] || !jeuReco[s][0] || !jeuReco[s][1] || !jeuReco[s][2]){
      printf("ERROR: missing input for %s\n", sysTag[s]); return false; }
    if(s > 0 && (!raw[s].fake[0] || !raw[s].fake[1])){ printf("ERROR: missing fakes for %s\n", sysTag[s]); return false; }
    resp[s]->SetDirectory(nullptr);
    { TH2D *sm = jerSmear::nominal(resp[s], Form("%s_jerNom", resp[s]->GetName())); delete resp[s]; resp[s] = sm; }   // JER nominal (jerSmear.h)
    for(int i = 0; i < nDimu[s]; i++) if(!dimu[s][i] || !dimuSS[s][i]){ printf("ERROR: missing dimuon mass for %s\n", sysTag[s]); return false; }
  }
  return true;
}

// ---------------------------------------------------------------- chain -------
double integ(TH1D *h, double lo, double hi, double *err = nullptr){
  double e = 0., v = h->IntegralAndError(h->FindBin(lo + 0.01), h->FindBin(hi - 0.01), e);
  if(err) *err = (v != 0.) ? e/v : 0.;
  return v;
}

// calculateJetsPerZ.cc's stitchSamples with its windows as options. relStat
// (optional) receives, per bin, the relative statistical error of the
// normalisation factors applied to that bin.
TH1D* stitch(int s, const Opts &oIn, const char *name, TH1D *relStat = nullptr){
  const Raw &R = raw[s];
  const bool isPP = (s == 0);
  double j60Lo = 150., j80Lo = 150.;
  const double hiEdge = 500.;
  Opts o = oIn;
  if(modeMB && turnOnStitch){
    j80Lo = j60Lo = tsJ80Lo[s];
    o.mbLo   = j80Lo + (oIn.mbLo - 130.);
    o.mbHi   = j80Lo + 40. + (oIn.mbHi - 170.);
    o.j100Lo = tsJ100Lo[s] + (oIn.j100Lo - 200.);
  }
  TH1D *out = (TH1D*) R.j100->Clone(name); out->SetDirectory(nullptr);

  if(evtCountStitch){
    // switch points moved by the stitchWindow variations as the windows are; no Jet60 leg
    const double sw1 = cntJ80Lo[s] + (oIn.mbLo - 130.), sw2 = TMath::Max(sw1, cntJ100Lo[s] + (oIn.j100Lo - 200.));
    TH1D *c80 = (TH1D*) R.j80->Clone("c80tmp"); c80->Scale(cntK80[s]);
    TH1D *cMB = (TH1D*) R.mb->Clone("cMBtmp");  cMB->Scale(cntKMB[s]);
    lastKMB[s] = cntKMB[s];
    if(!isPP && o.fakeAlt != 2){
      TH1D *fk = o.leak > 0. ? fakeLeakFit(R.fake[o.fakeAlt == 1 ? 1 : 0], R.mb, o.leak, o.fakeFitVar, "fktmp")
                             : (TH1D*) R.fake[o.fakeAlt == 1 ? 1 : 0]->Clone("fktmp");
      fk->SetDirectory(nullptr);
      fk->Scale(cntKMB[s]); cMB->Add(fk, -1.); delete fk;
    }
    for(int i = 0; i < out->GetSize(); i++){
      double pT = out->GetBinCenter(i);
      TH1D *src; double rs;
      if(mbOnly || pT < sw1) { src = cMB;    rs = cntKMBrel[s]; }
      else if(pT < sw2)      { src = c80;    rs = cntK80rel[s]; }
      else                   { src = R.j100; rs = 0.; }
      out->SetBinContent(i, src->GetBinContent(i)); out->SetBinError(i, src->GetBinError(i));
      if(relStat) relStat->SetBinContent(i, rs);
    }
    delete c80; delete cMB;
    return out;
  }

  double e100, e80, e60, e80b, eMB;
  double N100 = integ(R.j100, o.j100Lo, hiEdge, &e100), N80 = integ(R.j80, o.j100Lo, hiEdge, &e80);
  const double k80 = N100/N80, rel80 = sqrt(e100*e100 + e80*e80);
  TH1D *s80 = (TH1D*) R.j80->Clone("s80tmp"); s80->Scale(k80);
  double N80s = integ(s80, j80Lo, o.j100Lo);
  // PbPb has no Jet60 here (and never uses it: j60Lo = j80Lo); s60 is a dummy.
  // Turn-on stitch (MBcaloMatched): no Jet60 bridge in pp either since 2026-11-05
  // (user) -- MinBias is normalized to the scaled Jet80 directly, as in PbPb
  // (the bridge moved the pp normalization by 0.1%, j60bridge check)
  TH1D *j60 = (modeMB && turnOnStitch) ? nullptr : R.j60;
  TH1D *s60 = (TH1D*) (j60 ? j60 : R.j80)->Clone("s60tmp");
  double rel60 = rel80;
  if(j60){
    double N60  = integ(j60, j80Lo, o.j100Lo, &e60);
    integ(R.j80, j80Lo, o.j100Lo, &e80b);
    s60->Scale(N80s/N60);
    rel60 = sqrt(rel80*rel80 + e60*e60 + e80b*e80b);
  }

  double eRef;
  const bool viaJ60 = isPP && j60;
  double Nref = viaJ60 ? integ(s60, o.mbLo, o.mbHi, &eRef) : integ(s80, o.mbLo, o.mbHi, &eRef);
  double NMB  = integ(R.mb, o.mbLo, o.mbHi, &eMB);
  const double kMB = Nref/NMB;
  lastKMB[s] = kMB;
  const double relMB = sqrt((viaJ60 ? rel60*rel60 : rel80*rel80) + eRef*eRef + eMB*eMB);
  TH1D *sMB = (TH1D*) R.mb->Clone("sMBtmp"); sMB->Scale(kMB);
  if(!isPP && o.fakeAlt != 2){
    TH1D *fk = o.leak > 0. ? fakeLeakFit(R.fake[o.fakeAlt == 1 ? 1 : 0], R.mb, o.leak, o.fakeFitVar, "fktmp")   // leak removed, fitted
                           : (TH1D*) R.fake[o.fakeAlt == 1 ? 1 : 0]->Clone("fktmp");
    fk->SetDirectory(nullptr);
    fk->Scale(kMB); sMB->Add(fk, -1.); delete fk;
  }

  for(int i = 0; i < out->GetSize(); i++){
    double pT = out->GetBinCenter(i);
    TH1D *src = nullptr; double rs = 0.;
    if(mbOnly)                           { src = sMB;    rs = relMB; }
    else if(pT < j60Lo)                  { src = sMB;    rs = relMB; }
    else if(pT > j60Lo && pT < j80Lo)    { src = s60;    rs = rel60; }
    else if(pT > j80Lo && pT < o.j100Lo) { src = s80;    rs = rel80; }
    else if(pT > o.j100Lo)               { src = R.j100; rs = 0.; }
    if(!src) continue;
    out->SetBinContent(i, src->GetBinContent(i)); out->SetBinError(i, src->GetBinError(i));
    if(relStat) relStat->SetBinContent(i, rs);
  }
  delete s80; delete s60; delete sMB;
  return out;
}

void zeroBelow(TH1D *h, double floor){
  for(int b = 0; b <= h->GetNbinsX()+1; b++)
    if(h->GetXaxis()->GetBinUpEdge(b) <= floor + 1e-6){ h->SetBinContent(b, 0.); h->SetBinError(b, 0.); }
}
TH2D* floorResp(TH2D *r, double floor, const char *name){
  TH2D *t = (TH2D*) r->Clone(name); t->SetDirectory(nullptr);
  for(int ix = 0; ix <= t->GetNbinsX()+1; ix++){
    if(t->GetXaxis()->GetBinUpEdge(ix) > floor + 1e-6) continue;
    for(int iy = 0; iy <= t->GetNbinsY()+1; iy++){ t->SetBinContent(ix, iy, 0.); t->SetBinError(ix, iy, 0.); }
  }
  return t;
}

// truth from the untruncated matrix, so jets cut by the floor are misses
TH1D* unfold(TH2D *r, TH1D *meas, int nIt, double floor, const char *name){
  TH1D *truth = r->ProjectionY(Form("%s_truth", name), 1, r->GetNbinsX()); truth->SetDirectory(nullptr);
  TH2D *rf = floorResp(r, floor, Form("%s_rf", name));
  TH1D *mMC = rf->ProjectionX(Form("%s_mMC", name), 1, rf->GetNbinsY()); mMC->SetDirectory(nullptr);
  TH1D *m = (TH1D*) meas->Clone(Form("%s_meas", name)); m->SetDirectory(nullptr); zeroBelow(m, floor);
  RooUnfoldResponse ru(mMC, truth, rf, Form("%s_ru", name), "");
  RooUnfoldBayes ub(&ru, m, nIt); ub.SetVerbose(0);
  TH1D *u = (TH1D*) ub.Hunfold(); u->SetName(name); u->SetDirectory(nullptr);
  delete truth; delete rf; delete mMC; delete m;
  return u;
}

// JEU data-shaped shift: headers/functions/jeuShift.h (jeuDelta, densAt, jeuShiftSpectrum)

// our own Z->mumu count (OS, o.zLo-o.zHi, / eps^2), whatever the N_Z source
double NZown(int s, const Opts &o, double *nRaw = nullptr){
  double raw_ = 0., corr = 0.;
  for(int i = 0; i < nDimu[s]; i++){
    TH1D *h = dimu[s][i];
    double n = h->Integral(h->FindBin(o.zLo), h->FindBin(o.zHi - 1e-6));
    raw_ += n; corr += n/(dimuEff[s][i]*dimuEff[s][i]);
  }
  if(nRaw) *nRaw = raw_;
  return corr;
}

// N_Z used by the chain. NOMINAL since 2026-11-09 (user): the published CMS yields
// converted with the luminosity (headers/publishedData/CMS_ZYield_5p02TeV.h); the mass
// window o.zLo / o.zHi then plays no role, and nRaw is that count. BJET_NZ_SOURCE=own =
// our own count (NZown), outputs _ownNZ.
double NZ(int s, const Opts &o, double *nRaw = nullptr){
  if(!pubZ::enabled()) return NZown(s, o, nRaw);
  const double n = s == 0 ? pubZ::nzPP() : pubZ::nzSlices(nzSlLo[s], nzSlHi[s]);
  if(nRaw) *nRaw = n;
  return n;
}

// Fraction of the Z window removed by subtracting the same-sign spectrum bin by
// bin, raw counts (the efficiency cancels). A bin whose OS - SS is negative
// contributes no Z rather than a negative count, so it removes at most its own
// OS content: removed = sum_b min(SS_b, OS_b). PbPb bins are summed over the
// class's 5% slices first, as N_Z is.
double sameSignFrac(int s, const Opts &o){
  TH1D *h0 = dimu[s][0];
  const int bLo = h0->FindBin(o.zLo + 1e-6), bHi = h0->FindBin(o.zHi - 1e-6);
  double os = 0., removed = 0.;
  for(int b = bLo; b <= bHi; b++){
    double osB = 0., ssB = 0.;
    for(int i = 0; i < nDimu[s]; i++){ osB += dimu[s][i]->GetBinContent(b); ssB += dimuSS[s][i]->GetBinContent(b); }
    os += osB;
    removed += TMath::Max(0., TMath::Min(ssB, osB));
  }
  return os > 0. ? removed/os : 0.;
}

// unfolded jets per Z, 5 GeV scan axis, counts
TH1D* jetsPerZ(int s, const Opts &o, int nIt, TH2D *r, const char *name){
  TH1D *m = stitch(s, o, Form("%s_st", name));
  TH1D *u = unfold(r ? r : resp[s], m, nIt, o.floor, name);
  u->Scale(1./NZ(s, o));
  delete m;
  return u;
}

TH1D* rebinAx(TH1D *h, const Axis &a, const char *name){
  TH1D *o = new TH1D(name, h->GetTitle(), a.n, a.edges); o->SetDirectory(nullptr); o->Sumw2();
  for(int i = 1; i <= h->GetNbinsX(); i++){
    int j = o->FindBin(h->GetBinCenter(i)); if(j < 1 || j > a.n) continue;
    o->SetBinContent(j, o->GetBinContent(j) + h->GetBinContent(i));
    o->SetBinError(j, sqrt(pow(o->GetBinError(j), 2) + pow(h->GetBinError(i), 2)));
  }
  return o;
}

// N_opt from the DATA-LIKE closure (2026-11-05, user): per system the data tilt k is fitted (fitPriorTilt,
// as the unfoldingBias block), the response's gen rows are tilted by it, folded into pseudo-data with the
// data's relative errors, and unfolded with the nominal (untilted) response for N = 1..nMax; N_opt is the
// minimum of MSE = sum_i [((u_i - t_i)/t_i)^2 + (sigma_i/t_i)^2] over the reporting bins from 80 GeV, on the
// chain's own reporting axis. k is refitted at N_opt and the scan repeated until N_opt is stable (<= 3 passes).
// The iterations systematic is then N_opt +- 1 and unfoldingBias the residual |u/t - 1| at N_opt.
bool deriveNoptDataLike = true;
// 2026-11-05 (user): the closure truth is gen rows x (g_c/100)^(+closureTilt), g_c clamped to [40, 500] GeV
// (closureTiltRows, dataDrivenPrior.h; 0.3 = harder, nominal); < -90 = fit the data tilt (tiltGenRows, unclamped)
// (the data-like variant of earlier the same day). Used for N_opt and for unfoldingBias alike.
double closureTilt = 0.3;
// STUDY (2026-11-05, user): INCL_CLOSURE_FLOOR=<GeV> sets the reco floor of the closure only (N_opt scan and
// unfoldingBias), the data unfolding keeps nomFloor; unset = the nominal floor. Use with a tag.
inline double closureFloor(double nominal){ static double v = -2.; if(v < -1.){ const char *e = std::getenv("INCL_CLOSURE_FLOOR"); v = e ? atof(e) : -1.;
  if(v >= 0.) printf("STUDY: closure reco floor %.0f GeV (data unfolding unchanged)\n", v); } return v >= 0. ? v : nominal; }
std::vector<TObject*> dlSave;   // the scan: g_dlClosure_{bias2,var,mse}_<sys> vs N, h_dlClosure_uOverT_<sys> at N_opt, h_dlClosure_k
int dataLikeNopt(int s, const Opts &o, int nStart, int nMax = 50){   // 50 since 2026-11-05: the clamped-tilt 0-10% MSE is still falling at 20
  TH1D *m = stitch(s, o, Form("dlM_%s", sysTag[s]));
  int n = nStart;
  for(int pass = 0; pass < 3; pass++){
    auto unfoldFn = [&](TH2D *fine, TH1D *&truthOut){
      truthOut = fine->ProjectionY(Form("dlT_%s", sysTag[s]), 1, fine->GetNbinsX()); truthOut->SetDirectory(nullptr);
      return unfold(fine, m, n, o.floor, Form("dlU_%s", sysTag[s])); };
    const double k = closureTilt > -90. ? closureTilt : fitPriorTilt(resp[s], unfoldFn, 80., 300., 8, false);
    TH2D *P = closureTilt > -90. ? closureTiltRows(resp[s], closureTilt, Form("dlP_%s", sysTag[s])) : tiltGenRows(resp[s], k, Form("dlP_%s", sysTag[s]));
    TH1D *pd = P->ProjectionX(Form("dlPd_%s", sysTag[s]), 1, P->GetNbinsY()); pd->SetDirectory(nullptr);
    for(int b = 1; b <= pd->GetNbinsX(); b++){ const double mc = m->GetBinContent(b); pd->SetBinError(b, mc > 0. ? pd->GetBinContent(b)*m->GetBinError(b)/mc : 0.); }
    TH1D *tt = P->ProjectionY(Form("dlTt_%s", sysTag[s]), 1, P->GetNbinsX()); tt->SetDirectory(nullptr);
    TH1D *Ta = rebinAx(tt, axes[0], "dlTa");
    int best = 1; double bestM = 1e30; std::vector<double> mseV;
    TGraph *gB = new TGraph(), *gV = new TGraph(), *gM = new TGraph();
    printf("  %s data-like k %+.3f (pass %d), MSE vs N:", sysTag[s], k, pass);
    for(int N = 1; N <= nMax; N++){
      TH1D *U = unfold(resp[s], pd, N, closureFloor(o.floor), "dlUn"), *Ua = rebinAx(U, axes[0], "dlUa");
      double b2 = 0., v = 0.;
      for(int i = 1; i <= axes[0].n; i++){
        if(Ta->GetBinLowEdge(i) < 80. - 1e-6) continue;
        const double t = Ta->GetBinContent(i); if(t <= 0.) continue;
        b2 += pow(Ua->GetBinContent(i)/t - 1., 2); v += pow(Ua->GetBinError(i)/t, 2);
      }
      const double mse = b2 + v;
      gB->SetPoint(N-1, N, b2); gV->SetPoint(N-1, N, v); gM->SetPoint(N-1, N, mse);
      printf(" %d:%.4f", N, mse);
      if(mse < bestM){ bestM = mse; }
      mseV.push_back(mse);
      delete U; delete Ua;
    }
    best = plateauN(mseV);   // 5% plateau rule (2026-11-05, user), not the absolute minimum
    { int nMin = 1; for(size_t i = 1; i < mseV.size(); i++) if(mseV[i] < mseV[nMin-1]) nMin = (int) i + 1;
      printf("  -> N_opt %d (5%% plateau; absolute minimum at N = %d)\n", best, nMin); }
    const bool stable = best == n || closureTilt > -90.; n = best;   // fixed tilt: k does not depend on N, one pass
    if(stable || pass == 2){   // keep the final pass
      gB->SetName(Form("g_dlClosure_bias2_%s", sysTag[s])); gV->SetName(Form("g_dlClosure_var_%s", sysTag[s])); gM->SetName(Form("g_dlClosure_mse_%s", sysTag[s]));
      dlSave.push_back(gB); dlSave.push_back(gV); dlSave.push_back(gM);
      TH1D *U = unfold(resp[s], pd, best, closureFloor(o.floor), "dlUb"), *Ua = rebinAx(U, axes[0], Form("h_dlClosure_uOverT_%s", sysTag[s]));
      Ua->Divide(Ta); for(int i = 1; i <= Ua->GetNbinsX(); i++) Ua->SetBinError(i, 0.);
      Ua->SetTitle(Form("data-like closure at N_opt = %d, k = %.3f", best, k)); dlSave.push_back(Ua); delete U;
    } else { delete gB; delete gV; delete gM; }
    delete P; delete pd; delete tt; delete Ta;
    if(stable) break;
  }
  delete m;
  return n;
}

// JER variation of the (nominal-smeared) response: jerSmear.h
TH2D* smearResp(TH2D *r, const char *name){ return jerSmear::variation(r, name); }

// ---------------------------------------------------------------- sources -----
// A source is held as signed relative shifts per variation (corr) or as an
// unsigned relative uncertainty (uncorr), per system, per axis bin.
struct Source {
  TString name, label;
  bool corr = true;          // ratio formed per variation
  bool ppUncorr = false;     // corr, except quadrature between pp and PbPb
  bool ratioOnly = false;    // zero on the spectra
  std::vector<std::vector<std::vector<double>>> d[NAx];   // [ax][var][sys][bin]   (corr)
  std::vector<std::vector<double>> u[NAx];                // [ax][sys][bin]        (uncorr)
};

// quantity q: 0..4 spectra, 5..8 PbPb/pp C1..C4, 9..11 R_CP C1..C3
const int NQ = 12;
void qSys(int q, int &a, int &b){
  if(q < 5){ a = q; b = -1; } else if(q < 9){ a = q - 4; b = 0; } else { a = q - 8; b = 4; }
}
TString qTag(int q){
  if(q < 5) return sysTag[q];
  if(q < 9) return Form("PbPbOverPP_C%d", q - 4);
  return Form("RCP_C%d", q - 8);
}
TString qLabel(int q){
  if(q < 5) return sysLabel[q];
  if(q < 9) return Form("%s / pp", sysLabel[q-4] + 5);
  return Form("%s / %s", sysLabel[q-8] + 5, sysLabel[4] + 5);
}

double relOf(const Source &S, int ax, int q, int i){
  int a, b; qSys(q, a, b);
  if(S.ratioOnly && b < 0) return 0.;
  if(!S.corr){
    double ua = S.u[ax][a][i];
    return b < 0 ? ua : sqrt(ua*ua + pow(S.u[ax][b][i], 2));
  }
  double m = 0., ma = 0., mb = 0.;
  for(auto &v : S.d[ax]){
    double da = v[a][i], db = b < 0 ? 0. : v[b][i];
    m  = TMath::Max(m, fabs((1. + da)/(1. + db) - 1.));
    ma = TMath::Max(ma, fabs(da)); if(b >= 0) mb = TMath::Max(mb, fabs(db));
  }
  if(S.ppUncorr && b == 0) return sqrt(ma*ma + mb*mb);
  return m;
}

// signed shift of each system's binned per-Z spectrum, var vs nom
std::vector<std::vector<double>> shifts(TH1D *const *var, TH1D *const *nom, int ax){
  std::vector<std::vector<double>> d(NS, std::vector<double>(axes[ax].n, 0.));
  for(int s = 0; s < NS; s++){
    if(!var[s]) continue;
    TH1D *V = rebinAx(var[s], axes[ax], "shV"), *N = rebinAx(nom[s], axes[ax], "shN");
    for(int i = 1; i <= axes[ax].n; i++){
      double n = N->GetBinContent(i);
      d[s][i-1] = n > 0. ? V->GetBinContent(i)/n - 1. : 0.;
    }
    delete V; delete N;
  }
  return d;
}

} // namespace ijz

using namespace ijzpf;

void inclusiveJetsPerZSys_PF(bool jeuPPCorr = true,
                                   TString mode = "",         // "" caloMatchedPF PbPb; "MBcaloMatched" (header)
                                   TString nIterOpt = "",     // "pp,C1,C2,C3,C4", e.g. "7,2,3,2,4" (tilt:0.3 optimum)
                                   TString tag = "",          // appended to every output name
                                   TString cent = "",         // "", "slices", "ATLAS", "10pctA", "10pctB" (header, CENTRALITY CLASSES)
                                   TString trig = "plain",    // PbPb Jet80/Jet100 for cent != "": plain, caloMatched, plainUF
                                   bool bjetAxis = false,     // slice modes: keep the b-jet axis (bJetChain's edges) instead of "atlas"; _bjetAxis in the name
                                   bool noFlatJEU = false,    // JEU exempt from flattenSources; _noFlatJEU in the name
                                   bool saveVariations = false,   // also write the nominal and every correlated variation's
                                                                  // jets-per-Z spectra (native binning) to inclJetsPerZVariations_<...>.root
                                   bool mbOnlyOpt = false,        // MinBias spectrum at every pT, no stitching (header); _MBonly
                                   bool fakeGeoOnly = false,      // fakeJets = the RCgeoCorr estimate only, no "no subtraction"
                                                                  // variation (2026-10-20, user); _fakeGeoOnly in the name
                                   bool evtCount = false)         // event-count stitch (header, EVENT-COUNT STITCH); _evtCountStitch
{
  tag = TString(pubZ::tag()) + tag;   // published N_Z nominal since 2026-11-09; BJET_NZ_SOURCE=own: _ownNZ, before every other prefix
  if(!ptBinGuard(tag)) return;   // 7-bin jet-pT axis variant (caloJetBPurityFitCore.h): the b-jet axis of this file
  evtCountStitch = evtCount;
  if(evtCount){
    if(cent != "" || mode != "MBcaloMatched"){ printf("ERROR: the event-count stitch needs mode MBcaloMatched and cent \"\"\n"); return; }
    tag = "_evtCountStitch" + tag;
  }
  fakeGeoOnlyMode = fakeGeoOnly;
  if(fakeGeoOnly) tag = "_fakeGeoOnly" + tag;
  noFlattenJEU = noFlatJEU;
  mbOnly = mbOnlyOpt;
  if(mbOnly){ tag = "_MBonly" + tag; flattenSources = false; }
  if(noFlatJEU) tag = "_noFlatJEU" + tag;
  centMode = cent; trigAA = trig;
  if(cent != "" && cent != "slices" && cent != "ATLAS" && cent != "10pctA" && cent != "10pctB" && cent != "RCP60A" && cent != "RCP60B" && cent != "RCP50"){ printf("ERROR: cent is \"\", slices, ATLAS, 10pctA, 10pctB, RCP60A or RCP60B\n"); return; }
  if(cent != "" && mode != "MBcaloMatched"){ printf("ERROR: cent needs mode MBcaloMatched\n"); return; }
  if(trig != "plain" && trig != "caloMatched" && trig != "plainUF"){ printf("ERROR: trig is plain, caloMatched or plainUF\n"); return; }
  if(trig == "plainUF" && (fAA80PlainUF == "" || fAA100PlainUF == "")){ printf("ERROR: plain-PF ultraFine Jet80/Jet100 paths not set (rescan pending)\n"); return; }
  modeMB = mode == "MBcaloMatched";
  nomFloor = 60.;   // 2026-11-05 (user): 60 in MBcaloMatched mode too (70 from 2026-10-08)
  if(modeMB){
    // 2026-11-05, floor 60: min MSE of unfoldClosureTest.C(<C>,15,"weighted","tilt:0.3",80,500,...,calo-matched
    // response as train = test, recoFloor 60): 7 / 6 / 5 / 5; pp 3 from the even / odd halves (4 train = test).
    // Before (floor 70): 3, 8, 6, 4, 3.
    // 2026-11-05 (user): N_opt from the closure with the truth tilted to each system's DATA-LIKE shape
    // (unfoldClosureTest.C tilt:k with k the fitted data tilt of the unfoldingBias block: pp -0.09 (even /
    // odd halves), PbPb -0.68 / -0.37 / -0.33 / -0.26 on the calo-matched response, floor 60, min MSE
    // over 80-500 GeV): 2 / 9 / 6 / 5 / 4. The iterations systematic is N_opt +- 1 and unfoldingBias is the
    // residual |unfolded / true - 1| of that closure at N_opt. Before (tilt:0.3 closure): 3, 7, 6, 5, 5.
    const int n[NS] = {2, 3, 4, 6, 2}; for(int s = 0; s < NS; s++) nIterNom[s] = n[s];   // starting values; dataLikeNopt() rederives them (2026-11-05: 2 / 3 / 4 / 6 / 2)
    outStem = "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/systematics/inclJetsPerZSys_PF_MBcaloMatched";
    collText = "PF jets (PbPb MinBias calo-matched)";
  } else if(mode != ""){ printf("ERROR: mode is \"\" or MBcaloMatched\n"); return; }
  const TString figTag = modeMB ? "PF_MBcaloMatched" : "caloMatchedPF_ppPFJets";
  if(!nIterOpt.IsNull()){
    TObjArray *tok = nIterOpt.Tokenize(",");
    if(tok->GetEntries() != NS){ printf("nIterOpt needs %d comma-separated values\n", NS); delete tok; return; }
    for(int s = 0; s < NS; s++) nIterNom[s] = ((TObjString*)tok->At(s))->String().Atoi();
    delete tok;
  }
  // class tables of the slice modes: 5% slices, label, stitch edges of the nominal
  // class each lies in (overlapping several: the larger edges), and N_iter = the
  // min-MSE closure optimum of each class (2026-10-09): unfoldClosureTest.C("C<n>",
  // 15,"weighted","tilt:0.3",80.,500.,false,1,<resp>,<resp>,<label>,-1,"","",false,70.)
  // with <resp> = the calo-matched 2026-10-06 response merged into these classes by
  // tools/mergeUltraFineToNominalCent_response.C(<ultraFine file>,"<mode>") (train =
  // test, as for the nominal N). Shallow centrally: 0-20% N = 9-15 and 0-10% 7-15
  // lie within 5% of the minimum MSE. The same test at floor 70 gives 10 / 9 / 6 / 4
  // for the nominal classes, whose nIterNom (8 / 6 / 4 / 3) was set at floor 60
  // (+2-11% in MSE); the "" path is left as it is.
  struct CentSet { const char *mode; int lo[4], hi[4]; const char *lab[4]; double j80[4], j100[4]; int nIt[4]; };
  // RCP60A / RCP60B (2026-10-20, user): 10% classes over a 60-80% baseline in the 4th
  // slot, so the file's RCP_C<n> = class n / 60-80% (n = 1..3); edges and N of each
  // 10% class as in 10pctA / 10pctB, 60-80% as 60-70 / 70-80%
  // RCP50 (2026-11-01, user): 0-20 / 20-50% over a 50-80% baseline in the 4th slot, edges and N as "ATLAS"
  const CentSet centSets[6] = {
    {"ATLAS",   {1, 5, 11, 1}, {4, 10, 16, 16}, {"PbPb 0-20%", "PbPb 20-50%", "PbPb 50-80%", "PbPb 0-80%"},
                {145., 145., 120., 145.}, {155., 155., 140., 155.}, {13, 8, 4, 13}},
    {"10pctA",  {1, 3, 5, 7},  {2, 4, 6, 8},   {"PbPb 0-10%", "PbPb 10-20%", "PbPb 20-30%", "PbPb 30-40%"},
                {140., 145., 145., 115.}, {155., 155., 155., 140.}, {10, 8, 6, 5}},
    {"10pctB",  {9, 11, 13, 15}, {10, 12, 14, 16}, {"PbPb 40-50%", "PbPb 50-60%", "PbPb 60-70%", "PbPb 70-80%"},
                {115., 120., 120., 120.}, {140., 140., 140., 140.}, {4, 3, 3, 3}},
    {"RCP60A",  {1, 3, 5, 13}, {2, 4, 6, 16},  {"PbPb 0-10%", "PbPb 10-20%", "PbPb 20-30%", "PbPb 60-80%"},
                {140., 145., 145., 120.}, {155., 155., 155., 140.}, {10, 8, 6, 3}},
    {"RCP60B",  {7, 9, 11, 13}, {8, 10, 12, 16}, {"PbPb 30-40%", "PbPb 40-50%", "PbPb 50-60%", "PbPb 60-80%"},
                {115., 115., 120., 120.}, {140., 140., 140., 140.}, {5, 4, 3, 3}},
    {"RCP50",   {1, 5, 1, 11}, {4, 10, 16, 16}, {"PbPb 0-20%", "PbPb 20-50%", "PbPb 0-80%", "PbPb 50-80%"},
                {145., 145., 145., 120.}, {155., 155., 155., 140.}, {13, 8, 13, 4}}};
  for(const CentSet &S : centSets){
    if(cent != S.mode) continue;
    for(int c = 0; c < 4; c++){
      clsLo[c] = S.lo[c]; clsHi[c] = S.hi[c]; sysLabel[c+1] = S.lab[c];
      tsJ80Lo[c+1] = S.j80[c]; tsJ100Lo[c+1] = S.j100[c];
      if(nIterOpt.IsNull()) nIterNom[c+1] = S.nIt[c];
    }
    if(!bjetAxis) axes[1] = {"atlas", 16, atlasPtEdges};
  }
  if(cent != ""){ tag = Form("_%sCent_%s%s", cent.Data(), trig.Data(), bjetAxis ? "_bjetAxis" : "") + tag; }
  const TString sfx = TString(jeuPPCorr ? "_JEUppCorr" : "_JEUppUncorr") + tag;
  const TString outRoot = outStem + sfx + ".root";
  initPlotStyle();
  gSystem->mkdir(gSystem->DirName(outRoot), kTRUE);
  gSystem->mkdir(figDir, kTRUE);
  if(!loadInputs()) return;
  const Opts nomO;

  // ---- N_Z: ours and the published-yield equivalent, always printed -----------
  printf("\nN_Z %s: system  ours (OS 60-120 GeV / eps^2)  published x lumi (rel. unc.)  ours / published\n", pubZ::enabled() ? "PUBLISHED (nominal)" : "OWN Z->mumu (BJET_NZ_SOURCE=own)");
  for(int s = 0; s < NS; s++){
    const double own = NZown(s, nomO), pub = s == 0 ? pubZ::nzPP() : pubZ::nzSlices(nzSlLo[s], nzSlHi[s]);
    const double ru = s == 0 ? sqrt(pow(pubZ::relUncPP(), 2) + pow(pubZ::relUncFy(), 2)) : pubZ::relUncSlices(nzSlLo[s], nzSlHi[s]);
    printf("  %-14s %10.1f  %10.1f (%4.1f%%)  %.4f\n", sysLabel[s], own, pub, 100.*ru, own/pub);
  }
  if(pubZ::enabled()) printf("  pp: sigma %.0f pb x f_y %.3f +- %.3f (Pythia 8) x %.2f pb^-1; PbPb: N_MB(0-100%%) %.4e (lumi %.3f nb^-1)\n",
                             pubZ::sigmaPPtot, pubZ::fYpp, pubZ::fYppErr, pubZ::lumiPP, pubZ::nMB100(), pubZ::lumiAA);

  // ---- reproduction check: N = 1 vs calculateJetsPerZ.cc's stored output -----
  // (calculateJetsPerZ.cc divides by our own N_Z: skipped with the published N_Z)
  if(centMode == "" && !pubZ::enabled()){
    TFile *fcAA = TFile::Open(modeMB ? fCheckMB : fCheckAA), *fcPP = TFile::Open(modeMB ? fCheckMB : fCheckPP);
    printf("\ncheck: this chain at N = 1 vs calculateJetsPerZ.cc (pp %s, PbPb %s)\n", gSystem->BaseName(fCheckPP.Data()), gSystem->BaseName(fCheckAA.Data()));
    for(int s = 0; s < NS && fcAA && fcPP; s++){
      TFile *fc = s == 0 ? fcPP : fcAA;
      TH1D *ref = nullptr; fc->GetObject(Form("jetsPerZ_fine_%s", sysTag[s]), ref);
      if(!ref){ printf("  %s: reference missing\n", sysTag[s]); continue; }
      TH1D *u = jetsPerZ(s, nomO, 1, nullptr, Form("chk_%s", sysTag[s]));
      double mx = 0.;
      for(int i = 1; i <= u->GetNbinsX(); i++)
        if(u->GetBinLowEdge(i) >= 60. && ref->GetBinContent(i) != 0.) mx = TMath::Max(mx, fabs(u->GetBinContent(i)/ref->GetBinContent(i) - 1.));
      printf("  %-3s max |rel. diff| above 60 GeV = %.2e\n", sysTag[s], mx);
      delete u;
    }
  }

  // ---- N_opt from the data-like closure (2026-11-05, user) ---------------------
  if(deriveNoptDataLike && modeMB){
    printf("\nN_opt from the data-like closure (MSE over the reporting bins from 80 GeV):\n");
    for(int s = 0; s < NS; s++) nIterNom[s] = dataLikeNopt(s, nomO, nIterNom[s]);
    printf("N_opt:"); for(int s = 0; s < NS; s++) printf(" %s %d", sysTag[s], nIterNom[s]); printf("\n");
  }

  // ---- nominal -----------------------------------------------------------------
  TH1D *nom[NS];
  for(int s = 0; s < NS; s++) nom[s] = jetsPerZ(s, nomO, nIterNom[s], nullptr, Form("nom_%s", sysTag[s]));

  std::vector<Source> srcs;
  std::vector<TH1D*> savedVar;   // saveVariations: clones named h_var_<source>_<k>_<sys>
  std::vector<TH2D*> savedResp;  // saveVariations: h_resp_{nominal,JERsmeared}_<sys>
  std::vector<TH2D*> savedToys;  // saveVariations: h_toys_responseStat_<sys>
  auto corrSource = [&](const char *name, const char *label, const std::vector<std::vector<TH1D*>> &vars){
    Source S; S.name = name; S.label = label; S.corr = true;
    if(saveVariations) for(size_t k = 0; k < vars.size(); k++) for(int s = 0; s < NS; s++) if(vars[k][s]){
      TH1D *c = (TH1D*) vars[k][s]->Clone(Form("h_var_%s_%zu_%s", name, k, sysTag[s])); c->SetDirectory(nullptr); savedVar.push_back(c); }
    for(int ax = 0; ax < NAx; ax++) for(auto &v : vars) S.d[ax].push_back(shifts(v.data(), nom, ax));
    return S;
  };
  auto runAll = [&](const Opts &o, const char *tag){
    std::vector<TH1D*> v(NS);
    for(int s = 0; s < NS; s++) v[s] = jetsPerZ(s, o, nIterNom[s], nullptr, Form("%s_%s", tag, sysTag[s]));
    return v;
  };

  // iterations (uncorr)
  {
    Source S; S.name = "iterations"; S.label = "iterations (#it{N} #pm 1)"; S.corr = false;
    for(int ax = 0; ax < NAx; ax++) S.u[ax].assign(NS, std::vector<double>(axes[ax].n, 0.));
    for(int s = 0; s < NS; s++){
      TH1D *lo = jetsPerZ(s, nomO, nIterNom[s] - 1, nullptr, Form("itLo_%s", sysTag[s]));
      TH1D *hi = jetsPerZ(s, nomO, nIterNom[s] + 1, nullptr, Form("itHi_%s", sysTag[s]));
      if(saveVariations){   // N - 1 (k = 0) and N + 1 (k = 1), for the extraction figures (2026-10-21)
        TH1D *c0 = (TH1D*) lo->Clone(Form("h_var_iterations_0_%s", sysTag[s])), *c1 = (TH1D*) hi->Clone(Form("h_var_iterations_1_%s", sysTag[s]));
        c0->SetDirectory(nullptr); c1->SetDirectory(nullptr); savedVar.push_back(c0); savedVar.push_back(c1);
      }
      TH1D *V[NS] = {}, *W[NS] = {}; V[s] = lo; W[s] = hi;
      for(int ax = 0; ax < NAx; ax++){
        auto a = shifts(V, nom, ax), b = shifts(W, nom, ax);
        for(int i = 0; i < axes[ax].n; i++) S.u[ax][s][i] = TMath::Max(fabs(a[s][i]), fabs(b[s][i]));
      }
      delete lo; delete hi;
    }
    srcs.push_back(S);
  }

  // unfolding bias = residual non-closure at N_opt (uncorr; 2026-10-21, user; N_opt from this same closure
  // since 2026-11-05): the same-sample closure at each system's DATA-LIKE
  // tilt -- k fitted so that the response's gen rows x (g/100)^-k follow the unfolded data (fitPriorTilt,
  // dataDrivenPrior.h, 80-300 GeV; measured only, the prior stays the MC one); pseudo-data = the tilted
  // response folded, with the data's relative errors; unfolded with the nominal response and N;
  // |unfolded / true - 1|. Smooth and without the even/odd halves' statistical noise.
  {
    Source S; S.name = "unfoldingBias"; S.label = "unfolding non-closure (tilted truth)"; S.corr = false;
    for(int ax = 0; ax < NAx; ax++) S.u[ax].assign(NS, std::vector<double>(axes[ax].n, 0.));
    printf("\nunfolding bias: data-like tilt k per system:");
    for(int s = 0; s < NS; s++){
      TH1D *m = stitch(s, nomO, Form("ubM_%s", sysTag[s]));
      auto unfoldFn = [&](TH2D *fine, TH1D *&truthOut){
        truthOut = fine->ProjectionY(Form("ubT_%s", sysTag[s]), 1, fine->GetNbinsX()); truthOut->SetDirectory(nullptr);
        return unfold(fine, m, nIterNom[s], nomO.floor, Form("ubU_%s", sysTag[s])); };
      const double k = closureTilt > -90. ? closureTilt : fitPriorTilt(resp[s], unfoldFn, 80., 300., 8, false);
      printf("  %s %+.3f", sysTag[s], k);
      TH2D *P = closureTilt > -90. ? closureTiltRows(resp[s], closureTilt, Form("ubP_%s", sysTag[s])) : tiltGenRows(resp[s], k, Form("ubP_%s", sysTag[s]));
      TH1D *pd = P->ProjectionX(Form("ubPd_%s", sysTag[s]), 1, P->GetNbinsY()); pd->SetDirectory(nullptr);
      for(int b = 1; b <= pd->GetNbinsX(); b++){ const double mc = m->GetBinContent(b); pd->SetBinError(b, mc > 0. ? pd->GetBinContent(b)*m->GetBinError(b)/mc : 0.); }
      TH1D *tt = P->ProjectionY(Form("ubTt_%s", sysTag[s]), 1, P->GetNbinsX()); tt->SetDirectory(nullptr);
      TH1D *U = unfold(resp[s], pd, nIterNom[s], closureFloor(nomO.floor), Form("ubUn_%s", sysTag[s]));
      for(int ax = 0; ax < NAx; ax++){
        TH1D *Ua = rebinAx(U, axes[ax], "ubUa"), *Ta = rebinAx(tt, axes[ax], "ubTa");
        for(int i = 1; i <= axes[ax].n; i++){ const double t = Ta->GetBinContent(i); S.u[ax][s][i-1] = t > 0. ? fabs(Ua->GetBinContent(i)/t - 1.) : 0.; }
        delete Ua; delete Ta;
      }
      delete m; delete P; delete pd; delete tt; delete U;
    }
    printf("\n");
    srcs.push_back(S);
  }

  // response MC statistics, toys (uncorr)
  {
    Source S; S.name = "responseStat"; S.label = "response MC stat."; S.corr = false;
    for(int ax = 0; ax < NAx; ax++) S.u[ax].assign(NS, std::vector<double>(axes[ax].n, 0.));
    TRandom3 rnd(toySeed);
    for(int s = 0; s < NS; s++){
      TH1D *stt = stitch(s, nomO, "toyStitch");
      const double nz = NZ(s, nomO);
      std::vector<double> s1[NAx], s2[NAx];
      TH1D *nb[NAx];
      for(int ax = 0; ax < NAx; ax++){ s1[ax].assign(axes[ax].n, 0.); s2[ax].assign(axes[ax].n, 0.); nb[ax] = rebinAx(nom[s], axes[ax], Form("toyNom_%d", ax)); }
      // saveVariations (2026-10-21): every toy's toy / nominal on the reporting ("incl") axis,
      // x = reporting bin, y = ratio, for the extraction figures
      TH2D *toyMap = saveVariations ? new TH2D(Form("h_toys_responseStat_%s", sysTag[s]), "response toys;reporting bin;toy / nominal",
                                               axes[0].n, axes[0].edges, 2000, 0.5, 1.5) : nullptr;
      if(toyMap) toyMap->SetDirectory(nullptr);
      for(int t = 0; t < nToys; t++){
        TH2D *r2 = (TH2D*) resp[s]->Clone("toyR"); r2->SetDirectory(nullptr);
        for(int ix = 1; ix <= r2->GetNbinsX(); ix++) for(int iy = 1; iy <= r2->GetNbinsY(); iy++){
          double c = resp[s]->GetBinContent(ix, iy), e = resp[s]->GetBinError(ix, iy);
          if(c <= 0. || e <= 0.) continue;
          double ne = (c/e)*(c/e);
          r2->SetBinContent(ix, iy, rnd.Poisson(ne)*e*e/c);
        }
        TH1D *u = unfold(r2, stt, nIterNom[s], nomO.floor, "toyU"); u->Scale(1./nz);
        for(int ax = 0; ax < NAx; ax++){
          TH1D *ub = rebinAx(u, axes[ax], "toyUb");
          for(int i = 1; i <= axes[ax].n; i++){
            double n = nb[ax]->GetBinContent(i); if(n <= 0.) continue;
            double x = ub->GetBinContent(i)/n; s1[ax][i-1] += x; s2[ax][i-1] += x*x;
            if(ax == 0 && toyMap) toyMap->Fill(axes[0].edges[i-1] + 1e-6, x);
          }
          delete ub;
        }
        delete u; delete r2;
      }
      for(int ax = 0; ax < NAx; ax++){
        for(int i = 0; i < axes[ax].n; i++){
          double m = s1[ax][i]/nToys;
          S.u[ax][s][i] = sqrt(TMath::Max(0., s2[ax][i]/nToys - m*m));
        }
        delete nb[ax];
      }
      delete stt;
      if(toyMap) savedToys.push_back(toyMap);
      printf("  response toys done: %s\n", sysTag[s]);
    }
    srcs.push_back(S);
  }

  // JEU: MC reco spectra shifted, unfolded with this response and N
  {
    std::vector<TH1D*> g[3];
    if(jeuDataShape){
      const std::vector<double> dPP = jeuDelta(jeuUncPP), dAA = jeuDelta(jeuUncAA);
      printf("\nJEU on the data-shaped spectrum: delta(|eta|<1.6 avg) pp %.4f / %.4f, PbPb %.4f / %.4f at 80 / 200 GeV\n", dPP[80], dPP[200], dAA[80], dAA[200]);
      printf("  closure of the analytic shift on the SIMULATED reco spectrum vs the scan's per-jet shift, up, at 80/130/250 GeV:\n");
      for(int s = 0; s < NS; s++){
        const std::vector<double> &dl = s == 0 ? dPP : dAA;
        TH1D *m = stitch(s, nomO, Form("jeuM_%s", sysTag[s]));
        TH1D *mu = jeuShiftSpectrum(m, dl, +1, Form("jeuMu_%s", sysTag[s])), *md = jeuShiftSpectrum(m, dl, -1, Form("jeuMd_%s", sysTag[s]));
        g[0].push_back(unfold(resp[s], m,  nIterNom[s], nomO.floor, Form("jeu0_%s", sysTag[s])));
        g[1].push_back(unfold(resp[s], mu, nIterNom[s], nomO.floor, Form("jeu1_%s", sysTag[s])));
        g[2].push_back(unfold(resp[s], md, nIterNom[s], nomO.floor, Form("jeu2_%s", sysTag[s])));
        // closure: analytic shift of the MC reco spectrum vs the scan's shifted MC spectrum
        TH1D *au = jeuShiftSpectrum(jeuReco[s][0], dl, +1, Form("jeuCl_%s", sysTag[s]));
        printf("    %-8s", sysTag[s]);
        for(double x : {80., 130., 250.}){
          const int i = jeuReco[s][0]->FindBin(x); const double n0 = jeuReco[s][0]->GetBinContent(i);
          printf("  analytic %.3f scan %.3f", n0 > 0. ? au->GetBinContent(i)/n0 : 0., n0 > 0. ? jeuReco[s][1]->GetBinContent(i)/n0 : 0.);
        }
        printf("\n");
        if(saveVariations){
          TH1D *in[3] = {m, mu, md};   // h_jeuData_<k>_<sys>: the measured reco spectrum, nominal / up / down
          for(int k = 0; k < 3; k++){ TH1D *c = (TH1D*) in[k]->Clone(Form("h_jeuData_%d_%s", k, sysTag[s])); c->SetDirectory(nullptr); savedVar.push_back(c); }
          TH1D *c = (TH1D*) au->Clone(Form("h_jeuClosureUp_%s", sysTag[s])); c->SetDirectory(nullptr); savedVar.push_back(c);
        }
        delete m; delete mu; delete md; delete au;
      }
    } else {
      for(int k = 0; k < 3; k++) for(int s = 0; s < NS; s++)
        g[k].push_back(unfold(resp[s], jeuReco[s][k], nIterNom[s], nomO.floor, Form("jeu%d_%s", k, sysTag[s])));
    }
    // saveVariations (2026-10-21): the JEU inputs and outputs for the extraction figures,
    // k = 0 nominal, 1 up, 2 down; h_var_JEU_* unfolded (data-shaped input when jeuDataShape),
    // h_jeuReco_* the simulated reco spectra of the old method (the closure reference)
    if(saveVariations) for(int k = 0; k < 3; k++) for(int s = 0; s < NS; s++){
      TH1D *c = (TH1D*) g[k][s]->Clone(Form("h_var_JEU_%d_%s", k, sysTag[s])); c->SetDirectory(nullptr); savedVar.push_back(c);
      TH1D *r = (TH1D*) jeuReco[s][k]->Clone(Form("h_jeuReco_%d_%s", k, sysTag[s])); r->SetDirectory(nullptr); savedVar.push_back(r);
    }
    Source S; S.name = "JEU"; S.label = jeuPPCorr ? "jet energy scale (pp-PbPb corr.)" : "jet energy scale (pp-PbPb uncorr.)"; S.corr = true; S.ppUncorr = !jeuPPCorr;
    for(int ax = 0; ax < NAx; ax++)
      for(int k = 1; k <= 2; k++) S.d[ax].push_back(shifts(g[k].data(), g[0].data(), ax));
    // the signed shifts per variation, for the correlated b / inclusive JEU of bJetsPerZSys_PF_PbPb.C
    // (2026-10-21): h_varRel_JEU_{up,down}_<sys>_<axis> = var / nom - 1
    for(int ax = 0; ax < NAx; ax++) for(int k = 0; k < 2; k++) for(int s = 0; s < NS; s++){
      TH1D *h = new TH1D(Form("h_varRel_JEU_%s_%s_%s", k == 0 ? "up" : "down", sysTag[s], axes[ax].tag), "JEU, signed shift (var / nom - 1);#it{p}_{T}^{jet} [GeV];var / nom - 1", axes[ax].n, axes[ax].edges);
      h->SetDirectory(nullptr);
      for(int i = 0; i < axes[ax].n; i++) h->SetBinContent(i+1, S.d[ax][k][s][i]);
      jeuVarOut.push_back(h);
    }
    // ... and the same shift evaluated on the MEASURED spectrum (analytic shift, jeuShift.h), the
    // method of the b-jet JEU, so that b / inclusive is formed with one method on both sides
    // (2026-10-21): h_varRel_JEUdata_{up,down}_<sys>_<axis>. Not a source of the inclusive jets.
    {
      const std::vector<double> dPP = jeuDelta(jeuUncPP), dAA = jeuDelta(jeuUncAA);
      std::vector<TH1D*> G[3];
      for(int s = 0; s < NS; s++){
        const std::vector<double> &dl = s == 0 ? dPP : dAA;
        TH1D *m = stitch(s, nomO, Form("jeuDM_%s", sysTag[s]));
        TH1D *mu = jeuShiftMeasured(m, dl, +1, Form("jeuDMu_%s", sysTag[s]), true), *md = jeuShiftMeasured(m, dl, -1, Form("jeuDMd_%s", sysTag[s]), true);   // local smooth fit   // fit-based shape since 2026-11-16 (jeuShift.h; BJET_JEU_SHAPE=measured = before)
        { TH1D *c = (TH1D*) m->Clone(Form("h_jeuDataMeas_%s", sysTag[s])); c->SetDirectory(nullptr); c->SetTitle("measured spectrum the JEUdata shift is evaluated on (counts)"); jeuVarOut.push_back(c); }
        G[0].push_back(unfold(resp[s], m,  nIterNom[s], nomO.floor, Form("jeuD0_%s", sysTag[s])));
        G[1].push_back(unfold(resp[s], mu, nIterNom[s], nomO.floor, Form("jeuD1_%s", sysTag[s])));
        G[2].push_back(unfold(resp[s], md, nIterNom[s], nomO.floor, Form("jeuD2_%s", sysTag[s])));
        delete m; delete mu; delete md;
      }
      for(int ax = 0; ax < NAx; ax++) for(int k = 1; k <= 2; k++){
        std::vector<std::vector<double>> d = shifts(G[k].data(), G[0].data(), ax);
        for(int s = 0; s < NS; s++){
          TH1D *h = new TH1D(Form("h_varRel_JEUdata_%s_%s_%s", k == 1 ? "up" : "down", sysTag[s], axes[ax].tag), "JEU on the measured spectrum, signed shift (var / nom - 1);#it{p}_{T}^{jet} [GeV];var / nom - 1", axes[ax].n, axes[ax].edges);
          h->SetDirectory(nullptr);
          for(int i = 0; i < axes[ax].n; i++) h->SetBinContent(i+1, d[s][i]);
          jeuVarOut.push_back(h);
        }
      }
      for(int k = 0; k < 3; k++) for(auto h : G[k]) delete h;
    }
    srcs.push_back(S);
  }

  // JER: smeared response
  {
    // 2026-11-05 (user): measured on the MC (jerSmear.h): the nominal response's reco spectrum unfolded with
    // the nominal and the varied response, the ratio applied to the data nominal bin by bin (5 GeV axis);
    // BJET_JER_INPUT=data = the data unfolded with the varied response. Both are printed.
    std::vector<TH1D*> v(NS);
    printf("\nJER (%s), var / nom - 1 on the reporting axis, MC | data:\n", jerSmear::inputText());
    for(int s = 0; s < NS; s++){
      TH2D *sm = smearResp(resp[s], Form("jer_%s", sysTag[s]));
      TH1D *vD = jetsPerZ(s, nomO, nIterNom[s], sm, Form("jerUd_%s", sysTag[s]));
      TH1D *mc = resp[s]->ProjectionX(Form("jerMC_%s", sysTag[s]), 1, resp[s]->GetNbinsY()); mc->SetDirectory(nullptr);
      TH1D *u0 = unfold(resp[s], mc, nIterNom[s], nomO.floor, Form("jerMC0_%s", sysTag[s]));
      TH1D *u1 = unfold(sm, mc, nIterNom[s], nomO.floor, Form("jerMC1_%s", sysTag[s]));
      TH1D *vM = (TH1D*) nom[s]->Clone(Form("jerUm_%s", sysTag[s])); vM->SetDirectory(nullptr);
      for(int b = 1; b <= vM->GetNbinsX(); b++){ const double a = u0->GetBinContent(b); const double r = a > 0. ? u1->GetBinContent(b)/a : 1.;
        vM->SetBinContent(b, nom[s]->GetBinContent(b)*r); vM->SetBinError(b, nom[s]->GetBinError(b)*r); }
      { TH1D *VM[NS] = {}, *VD[NS] = {}; VM[s] = vM; VD[s] = vD; auto a = shifts(VM, nom, 0), b = shifts(VD, nom, 0);
        printf("  %-4s", sysTag[s]); for(int i = 0; i < axes[0].n; i++) if(axes[0].edges[i] >= 80. - 1e-6 && axes[0].edges[i+1] <= 200. + 1e-6) printf("  %3.0f: %+.3f | %+.3f", axes[0].edges[i], a[s][i], b[s][i]); printf("\n"); }
      if(saveVariations){ TH1D *c = (TH1D*) (jerSmear::onData() ? vM : vD)->Clone(Form("h_var_JERalt_0_%s", sysTag[s])); c->SetDirectory(nullptr); c->SetTitle(jerSmear::onData() ? "JER measured on the MC" : "JER measured on the data"); savedVar.push_back(c); }
      v[s] = jerSmear::onData() ? vD : vM; delete (jerSmear::onData() ? vM : vD);
      delete mc; delete u0; delete u1;
      if(saveVariations){   // the responses behind the JER source, for the extraction figures (2026-10-21)
        TH2D *r0 = (TH2D*) resp[s]->Clone(Form("h_resp_nominal_%s", sysTag[s])), *r1 = (TH2D*) sm->Clone(Form("h_resp_JERsmeared_%s", sysTag[s]));
        r0->SetDirectory(nullptr); r1->SetDirectory(nullptr); savedResp.push_back(r0); savedResp.push_back(r1);
      }
      delete sm;
    }
    // correlated between pp and PbPb and between the PbPb classes: ratio of the varied spectra per variation
    // (2026-10-22, user; uncorrelated pp-PbPb, quadrature, from 2026-10-21 to 2026-10-22)
    Source J = corrSource("JER", jerSmear::varLabel(), {v});
    // the signed shift per system, for the CORRELATED b / inclusive JER of bJetsPerZSys_PF_PbPb.C
    // (2026-11-12, user): h_varRel_JER_<sys>_<axis> = var / nom - 1
    for(int ax = 0; ax < NAx; ax++) for(int s = 0; s < NS; s++){
      TH1D *h = new TH1D(Form("h_varRel_JER_%s_%s", sysTag[s], axes[ax].tag), "JER, signed shift (var / nom - 1);#it{p}_{T}^{jet} [GeV];var / nom - 1", axes[ax].n, axes[ax].edges);
      h->SetDirectory(nullptr);
      for(int i = 0; i < axes[ax].n; i++) h->SetBinContent(i+1, J.d[ax][0][s][i]);
      jeuVarOut.push_back(h);
    }
    srcs.push_back(J);
  }

  // fake jets: alternative estimate, the leakage fraction varied (and, in MBcaloMatched mode
  // without fakeGeoOnly, no subtraction); PbPb only
  {
    std::vector<std::vector<TH1D*>> vars;
    if(useFakeAlt){   // the alternative estimate (random-cone geometry correction); none for the forest CS
      Opts o; o.fakeAlt = 1;
      std::vector<TH1D*> v = runAll(o, "fake");
      delete v[0]; v[0] = nullptr;
      vars.push_back(v);
    }
    if(fakeLeakNom > 0.){
      for(double lk : {fakeLeakNom - 0.010, fakeLeakNom + 0.010}){
        Opts l; l.leak = lk;
        std::vector<TH1D*> w = runAll(l, Form("fakeLeak%.0f", 1000.*lk));
        delete w[0]; w[0] = nullptr; vars.push_back(w);
      }
      for(int fv : {+1, -1}){
        Opts l; l.fakeFitVar = fv;
        std::vector<TH1D*> w = runAll(l, fv > 0 ? "fakeFitUp" : "fakeFitDn");
        delete w[0]; w[0] = nullptr; vars.push_back(w);
      }
    }
    if(modeMB && !fakeGeoOnlyMode){ Opts n; n.fakeAlt = 2; std::vector<TH1D*> w = runAll(n, "fakeNone"); delete w[0]; w[0] = nullptr; vars.push_back(w); }
    if(!vars.empty()) srcs.push_back(corrSource("fakeJets", "fake-jet subtraction", vars));   // none: forest CS, no leakage (2026-11-05)
  }

  // stitching windows
  {
    Opts a, b, c; a.mbLo = 140.; a.mbHi = 180.; b.mbHi = 190.; c.j100Lo = 220.;
    srcs.push_back(corrSource("stitchWindow", "stitching windows", {runAll(a, "stA"), runAll(b, "stB"), runAll(c, "stC")}));
  }

  // stitching normalisation statistics (uncorr)
  {
    Source S; S.name = "stitchNormStat"; S.label = "stitch normalisation stat."; S.corr = false;
    for(int ax = 0; ax < NAx; ax++) S.u[ax].assign(NS, std::vector<double>(axes[ax].n, 0.));
    for(int s = 0; s < NS; s++){
      TH1D *rs = (TH1D*) raw[s].j100->Clone("relStat"); rs->Reset(); rs->SetDirectory(nullptr);
      TH1D *st = stitch(s, nomO, "stStat", rs);
      // reco-level relative error, weighted by the reco content of each output bin
      for(int ax = 0; ax < NAx; ax++){
        for(int j = 1; j <= axes[ax].n; j++){
          double w = 0., we = 0.;
          for(int i = 1; i <= st->GetNbinsX(); i++){
            double x = st->GetBinCenter(i);
            if(x < axes[ax].edges[j-1] || x > axes[ax].edges[j]) continue;
            double c = TMath::Max(0., st->GetBinContent(i)); w += c; we += c*rs->GetBinContent(i);
          }
          S.u[ax][s][j-1] = w > 0. ? we/w : 0.;
        }
      }
      delete rs; delete st;
    }
    srcs.push_back(S);
  }

  // reco floor
  {
    // MBcaloMatched: 70 -> 80 (2026-10-08), as the b-jet unfolding; else 60 -> 50.
    // Varying down into the fake-dominated 50-60 GeV bin measured the fake
    // subtraction there, not the floor.
    Opts o; o.floor = modeMB ? 70. : 50.;   // MBcaloMatched: 60 -> 70 since 2026-11-05 (70 -> 80 before)
    if(recoFloorSourceOn())   // off since 2026-11-16 (user: the non-closure at the 60 GeV floor, no extra source; purityVariant.h)
    srcs.push_back(corrSource("recoFloor", modeMB ? "reco floor (60 #rightarrow 70 GeV)" : "reco floor (60 #rightarrow 50 GeV)", {runAll(o, "floor")}));
  }

  // N_Z sources. Published N_Z (nominal since 2026-11-09): pubZYield = the
  // published yield (PbPb: Y stat. (+) syst. of the classes' bins; pp: sigma stat. (+)
  // syst.), uncorr; pubZAcc = the pp |y_Z| < 2.1 fraction f_y, pp only (uncorr: it then
  // enters every PbPb / pp ratio in full and cancels in R_CP). Our own count's sources
  // (NZstat, muonEff, zBkg, massWindow) exist only with BJET_NZ_SOURCE=own.
  if(pubZ::enabled()){
    Source A; A.name = "pubZYield"; A.label = "#it{N}_{Z} published yield"; A.corr = false;
    Source F; F.name = "pubZAcc";   F.label = "#it{N}_{Z} pp |#it{y}_{Z}| < 2.1 fraction"; F.corr = false;
    for(int ax = 0; ax < NAx; ax++){
      A.u[ax].assign(NS, std::vector<double>(axes[ax].n, 0.)); F.u[ax].assign(NS, std::vector<double>(axes[ax].n, 0.));
      for(int s = 0; s < NS; s++){
        const double ua = s == 0 ? pubZ::relUncPP() : pubZ::relUncSlices(nzSlLo[s], nzSlHi[s]);
        // spectra are 1/N_Z: a relative N_Z uncertainty u gives the same u to first order
        for(int i = 0; i < axes[ax].n; i++){ A.u[ax][s][i] = ua; F.u[ax][s][i] = s == 0 ? pubZ::relUncFy() : 0.; }
      }
    }
    srcs.push_back(A); srcs.push_back(F);
  } else {
  // N_Z statistics, muon efficiency (uncorr, flat in pT)
  {
    Source A; A.name = "NZstat"; A.label = "#it{N}_{Z} stat."; A.corr = false;
    Source B; B.name = "muonEff"; B.label = "muon efficiency (#it{N}_{Z})"; B.corr = false;
    for(int ax = 0; ax < NAx; ax++){
      A.u[ax].assign(NS, std::vector<double>(axes[ax].n, 0.)); B.u[ax].assign(NS, std::vector<double>(axes[ax].n, 0.));
      for(int s = 0; s < NS; s++){
        double nr; NZ(s, nomO, &nr);
        // N_Z-weighted relative error of the entries' efficiencies (linear: within a
        // nominal class the slices share one efficiency); = muEffErr/muEff for one class
        double w = 0., we = 0.;
        for(int j = 0; j < nDimu[s]; j++){
          TH1D *h = dimu[s][j]; double n = h->Integral(h->FindBin(nomO.zLo), h->FindBin(nomO.zHi - 1e-6))/pow(dimuEff[s][j], 2);
          w += n; we += n*dimuEffErr[s][j]/dimuEff[s][j];
        }
        for(int i = 0; i < axes[ax].n; i++){ A.u[ax][s][i] = 1./sqrt(nr); B.u[ax][s][i] = 2.*we/w; }
      }
    }
    srcs.push_back(A); srcs.push_back(B);
  }

  // Z background from the same-sign count (uncorr, flat in pT)
  {
    Source S; S.name = "zBkg"; S.label = "#it{N}_{Z} background (same sign)"; S.corr = false;
    printf("\nzBkg: removed by OS - SS (bins < 0 -> 0) in %g-%g GeV:", nomO.zLo, nomO.zHi);
    for(int s = 0; s < NS; s++) printf("  %s %.2f%%", sysTag[s], 100.*sameSignFrac(s, nomO));
    printf("\n");
    for(int ax = 0; ax < NAx; ax++){
      S.u[ax].assign(NS, std::vector<double>(axes[ax].n, 0.));
      for(int s = 0; s < NS; s++){ const double f = sameSignFrac(s, nomO); for(int i = 0; i < axes[ax].n; i++) S.u[ax][s][i] = f; }
    }
    srcs.push_back(S);
  }

  // Z mass window, ratios only
  {
    Opts o; o.zLo = 70.; o.zHi = 110.;
    Source S; S.name = "massWindow"; S.label = "#it{N}_{Z} mass window"; S.corr = true; S.ratioOnly = true;
    for(int ax = 0; ax < NAx; ax++){
      std::vector<std::vector<double>> d(NS, std::vector<double>(axes[ax].n, 0.));
      for(int s = 0; s < NS; s++){ double f = NZ(s, nomO)/NZ(s, o) - 1.; for(int i = 0; i < axes[ax].n; i++) d[s][i] = f; }
      S.d[ax].push_back(d);
    }
    srcs.push_back(S);
  }
  }

  // ---- outputs -----------------------------------------------------------------
  TFile *fo = TFile::Open(outRoot, "recreate");
  for(auto ob : dlSave) ob->Write();   // data-like closure scan (2026-11-05)
  for(auto h : jeuVarOut) h->Write();
  const int NSrc = srcs.size();
  for(int ax = 0; ax < NAx; ax++){
    const Axis &A = axes[ax];
    TH1D *nb[NS];
    for(int s = 0; s < NS; s++) nb[s] = rebinAx(nom[s], A, Form("h_nom_%s_%s", sysTag[s], A.tag));
    for(int q = 0; q < NQ; q++){
      int a, b; qSys(q, a, b);
      TH1D *v = (TH1D*) nb[a]->Clone(Form("h_nom_%s_%s", qTag(q).Data(), A.tag)); v->SetDirectory(nullptr);
      if(b >= 0) v->Divide(nb[b]);
      v->SetTitle(Form("%s;#it{p}_{T}^{jet} [GeV];%s", qLabel(q).Data(), b < 0 ? "jets per #it{N}_{Z} per bin" : "ratio of jets per #it{N}_{Z}"));
      TH1D *st = new TH1D(Form("h_statRel_%s_%s", qTag(q).Data(), A.tag), "stat.;#it{p}_{T}^{jet} [GeV];relative uncertainty", A.n, A.edges);
      TH1D *to = new TH1D(Form("h_sysRel_total_%s_%s", qTag(q).Data(), A.tag), "total systematic;#it{p}_{T}^{jet} [GeV];relative uncertainty", A.n, A.edges);
      std::vector<TH1D*> hs;
      for(int k = 0; k < NSrc; k++)
        hs.push_back(new TH1D(Form("h_sysRel_%s_%s_%s", srcs[k].name.Data(), qTag(q).Data(), A.tag), Form("%s;#it{p}_{T}^{jet} [GeV];relative uncertainty", srcs[k].label.Data()), A.n, A.edges));
      for(int i = 1; i <= A.n; i++){
        if(v->GetBinContent(i) <= 0.) continue;
        st->SetBinContent(i, v->GetBinError(i)/v->GetBinContent(i));
        for(int k = 0; k < NSrc; k++) hs[k]->SetBinContent(i, relOf(srcs[k], ax, q, i-1));
      }
      // a source with no pT trend is replaced by its mean over the reported
      // bins (headers/plotting/flattenSystematic.h); bins below ptReportMin
      // are left as they are and do not enter the fit
      if(flattenSources) for(int k = 0; k < NSrc; k++){
        if(noFlattenJEU && srcs[k].name == "JEU") continue;
        double mean, F; int bLo = hs[k]->FindBin(ptReportMin + 1e-6);
        if(flattenIfNoTrend(hs[k], mean, F, bLo)) printf("  flattened %-14s %-24s %-5s mean %5.1f%%  F = %.1f\n", srcs[k].name.Data(), qLabel(q).Data(), A.tag, 100*mean, F);
      }
      for(int i = 1; i <= A.n; i++){
        if(v->GetBinContent(i) <= 0.) continue;
        double t2 = 0.;
        for(int k = 0; k < NSrc; k++){ double r = hs[k]->GetBinContent(i); t2 += r*r; }
        to->SetBinContent(i, sqrt(t2));
      }
      fo->cd(); v->Write(); st->Write(); to->Write(); for(auto h : hs) h->Write();
    }
  }

  // ---- printout (incl axis) ----------------------------------------------------
  for(int q = 0; q < NQ; q++){
    printf("\n=== %s (%%) ===\n%-8s %6s", qLabel(q).Data(), "pT", "stat");
    for(int k = 0; k < NSrc; k++) printf(" %8.8s", srcs[k].name.Data());
    printf(" %7s\n", "TOTAL");
    TH1D *st = (TH1D*) fo->Get(Form("h_statRel_%s_incl", qTag(q).Data())), *to = (TH1D*) fo->Get(Form("h_sysRel_total_%s_incl", qTag(q).Data()));
    for(int i = 1; i <= NInclE - 1; i++){
      if(to->GetBinContent(i) <= 0.) continue;
      printf("%3.0f-%-4.0f %6.1f", inclEdges[i-1], inclEdges[i], 100*st->GetBinContent(i));
      for(int k = 0; k < NSrc; k++) printf(" %8.1f", 100*((TH1D*) fo->Get(Form("h_sysRel_%s_%s_incl", srcs[k].name.Data(), qTag(q).Data())))->GetBinContent(i));
      printf(" %7.1f\n", 100*to->GetBinContent(i));
    }
  }

  // ---- figures -----------------------------------------------------------------
  // breakdown: one panel per quantity, incl axis, reported range only
  const char *srcHex[12]  = {"#E69F00", "#CC79A7", "#0072B2", "#56B4E9", "#D55E00", "#009E73", "#009E73", "#999999", "#E69F00", "#CC79A7", "#0072B2", "#56B4E9"};
  const int   srcStyle[12] = {1, 1, 1, 2, 1, 1, 2, 3, 7, 7, 7, 7};
  // panels on top, one legend strip below shared by all of them
  auto breakdown = [&](int q0, int nq, int nx, int ny, TString fileS){
    const char *file = fileS.Data();   // own copy: Form()'s buffer is reused below
    const double legH = ny > 1 ? 0.13 : 0.24;   // fraction of the canvas
    TCanvas *cv = new TCanvas(Form("cvB_%s", file), "", 520*nx, (int)(460*ny/(1. - legH)));
    TLegend *lg = nullptr;
    for(int p = 0; p < nq; p++){
      int q = q0 + p, col = p % nx, row = p / nx;
      double x1 = (double) col/nx, x2 = (double)(col + 1)/nx;
      double y2 = 1. - (1. - legH)*row/ny, y1 = 1. - (1. - legH)*(row + 1)/ny;
      cv->cd();
      TPad *pd = new TPad(Form("pd_%s_%d", file, p), "", x1, y1, x2, y2); pd->Draw(); pd->cd();
      pd->SetLeftMargin(0.14); pd->SetRightMargin(0.03); pd->SetTopMargin(0.09); pd->SetBottomMargin(0.13);
      TH1D *to = (TH1D*) fo->Get(Form("h_sysRel_total_%s_incl", qTag(q).Data()))->Clone();
      TH1D *st = (TH1D*) fo->Get(Form("h_statRel_%s_incl", qTag(q).Data()))->Clone();
      double ymax = 0.;
      for(int i = 1; i <= to->GetNbinsX(); i++) if(to->GetBinLowEdge(i) >= ptReportMin - 1e-6) ymax = TMath::Max(ymax, TMath::Max(to->GetBinContent(i), st->GetBinContent(i)));
      TH1F *fr = pd->DrawFrame(ptReportMin, 0., inclEdges[NInclE-1], 100.*ymax*1.1);
      fr->GetXaxis()->SetTitle("#it{p}_{T}^{jet} [GeV]"); fr->GetYaxis()->SetTitle("relative uncertainty [%]");
      fr->GetXaxis()->SetTitleSize(0.050); fr->GetYaxis()->SetTitleSize(0.050);
      fr->GetXaxis()->SetLabelSize(0.045); fr->GetYaxis()->SetLabelSize(0.045); fr->GetYaxis()->SetTitleOffset(1.3);
      for(TH1D *h : {to, st}){ h->Scale(100.); for(int i = 0; i <= h->GetNbinsX()+1; i++) h->SetBinError(i, 0.); h->GetXaxis()->SetRangeUser(ptReportMin, inclEdges[NInclE-1]); }
      styleLine(to, "#000000", 5); to->Draw("HIST same");        // first, so the sources stay visible on it
      if(p == 0){ lg = makeLegend(0.03, 0.05, 0.99, 0.95); lg->SetTextFont(43); lg->SetTextSize(19); lg->SetNColumns(4); lg->AddEntry(to, "total syst.", "l"); }
      for(int k = 0; k < NSrc; k++){
        TH1D *h = (TH1D*) fo->Get(Form("h_sysRel_%s_%s_incl", srcs[k].name.Data(), qTag(q).Data()))->Clone();
        h->Scale(100.); for(int i = 0; i <= h->GetNbinsX()+1; i++) h->SetBinError(i, 0.);
        h->GetXaxis()->SetRangeUser(ptReportMin, inclEdges[NInclE-1]);
        styleLine(h, srcHex[k % 12], 2); h->SetLineStyle(srcStyle[k % 12]);
        h->Draw("HIST same");
        if(p == 0) lg->AddEntry(h, srcs[k].label, "l");
      }
      styleLine(st, "#000000", 2); st->SetLineStyle(2); st->Draw("HIST same");
      if(p == 0) lg->AddEntry(st, "stat.", "l");
      TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.050);
      la.DrawLatex(0.17, 0.935, Form("%s, %s: inclusive jets per #it{Z}", qLabel(q).Data(), collText));
    }
    cv->cd();
    TPad *pl = new TPad(Form("pl_%s", file), "", 0., 0., 1., legH); pl->Draw(); pl->cd();
    lg->Draw();
    savePdfTight(cv, Form("%s/%s%s.pdf", figDir.Data(), file, sfx.Data()));
    delete cv;
  };
  breakdown(0, 5, 3, 2, Form("inclJetsPerZSys_%s_spectra", figTag.Data()));
  breakdown(5, 4, 2, 2, Form("inclJetsPerZSys_%s_PbPbOverPP", figTag.Data()));
  breakdown(9, 3, 3, 1, Form("inclJetsPerZSys_%s_RCP", figTag.Data()));

  // ratios with stat. bars and systematic boxes
  auto ratioFig = [&](int q0, int nq, int nx, int ny, const char *ytitle, TString fileS){
    const char *file = fileS.Data();
    TCanvas *cv = new TCanvas(Form("cvR_%s", file), "", 520*nx, 460*ny);
    cv->Divide(nx, ny, 0.001, 0.001);
    for(int p = 0; p < nq; p++){
      int q = q0 + p;
      cv->cd(p + 1);
      gPad->SetLeftMargin(0.14); gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.08); gPad->SetBottomMargin(0.13);
      TH1D *v  = (TH1D*) fo->Get(Form("h_nom_%s_incl", qTag(q).Data()))->Clone();
      TH1D *to = (TH1D*) fo->Get(Form("h_sysRel_total_%s_incl", qTag(q).Data()));
      double ymax = 0.;
      for(int i = 1; i <= v->GetNbinsX(); i++) if(v->GetBinLowEdge(i) >= ptReportMin - 1e-6)
        ymax = TMath::Max(ymax, v->GetBinContent(i)*(1. + to->GetBinContent(i)) + v->GetBinError(i));
      TH1F *fr = gPad->DrawFrame(ptReportMin, 0., inclEdges[NInclE-1], TMath::Max(1.3, 1.35*ymax));
      fr->GetXaxis()->SetTitle("#it{p}_{T}^{jet} [GeV]"); fr->GetYaxis()->SetTitle(ytitle);
      fr->GetXaxis()->SetTitleSize(0.050); fr->GetYaxis()->SetTitleSize(0.050);
      fr->GetXaxis()->SetLabelSize(0.045); fr->GetYaxis()->SetLabelSize(0.045); fr->GetYaxis()->SetTitleOffset(1.3);
      TLine *one = new TLine(ptReportMin, 1., inclEdges[NInclE-1], 1.); one->SetLineStyle(7); one->SetLineColor(kGray+2); one->Draw();
      const int col = TColor::GetColor("#D55E00");
      for(int i = 1; i <= v->GetNbinsX(); i++){
        double y = v->GetBinContent(i); if(y <= 0. || v->GetBinLowEdge(i) < ptReportMin - 1e-6) continue;
        TBox *bx = new TBox(v->GetBinLowEdge(i), y*(1. - to->GetBinContent(i)), v->GetBinLowEdge(i+1), y*(1. + to->GetBinContent(i)));
        bx->SetFillColorAlpha(col, 0.30); bx->SetLineColor(col); bx->Draw();
      }
      for(int i = 1; i <= v->GetNbinsX(); i++) if(v->GetBinLowEdge(i) < ptReportMin - 1e-6){ v->SetBinContent(i, -999.); v->SetBinError(i, 0.); }
      styleH(v, "#000000", markFilledCircle, 1.1); v->Draw("E1 X0 same");
      TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.050);
      la.DrawLatex(0.17, 0.935, Form("%s, %s", qLabel(q).Data(), collText));
      if(p == 0){
        TLegend *lg = makeLegend(0.17, 0.74, 0.70, 0.88, 0.042);
        lg->AddEntry(v, "stat.", "pe");
        TBox *lb = new TBox(); lb->SetFillColorAlpha(col, 0.30); lb->SetLineColor(col);
        lg->AddEntry(lb, jeuPPCorr ? "total syst. (JEU pp-PbPb corr.)" : "total syst. (JEU pp-PbPb uncorr.)", "f"); lg->Draw();
      }
    }
    savePdfTight(cv, Form("%s/%s%s.pdf", figDir.Data(), file, sfx.Data()));
    delete cv;
  };
  // the spectra themselves: one panel per system, stat. bars and systematic
  // boxes, pp repeated as open markers in every PbPb panel; legend in the spare cell
  {
    TCanvas *cv = new TCanvas("cvSpec", "", 520*3, 460*2);
    cv->Divide(3, 2, 0.001, 0.001);
    const int col = TColor::GetColor("#D55E00");
    TH1D *ppD = nullptr;
    for(int q = 0; q < 5; q++){
      cv->cd(q + 1); gPad->SetLogy();
      gPad->SetLeftMargin(0.16); gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.08); gPad->SetBottomMargin(0.13);
      TH1D *v  = (TH1D*) fo->Get(Form("h_nom_%s_incl", qTag(q).Data()))->Clone(Form("spec_%d", q)); v->SetDirectory(nullptr);
      TH1D *to = (TH1D*) fo->Get(Form("h_sysRel_total_%s_incl", qTag(q).Data()));
      v->Scale(1., "width");
      double ymax = 0., ymin = 1e30;
      for(int i = 1; i <= v->GetNbinsX(); i++) if(v->GetBinLowEdge(i) >= ptReportMin - 1e-6 && v->GetBinContent(i) > 0.){
        ymax = TMath::Max(ymax, v->GetBinContent(i)*(1. + to->GetBinContent(i))); ymin = TMath::Min(ymin, v->GetBinContent(i)*(1. - to->GetBinContent(i))); }
      if(q == 0){ ppD = (TH1D*) v->Clone("spec_ppRef"); ppD->SetDirectory(nullptr); }
      else for(int i = 1; i <= ppD->GetNbinsX(); i++) if(ppD->GetBinLowEdge(i) >= ptReportMin - 1e-6 && ppD->GetBinContent(i) > 0.){ ymax = TMath::Max(ymax, ppD->GetBinContent(i)); ymin = TMath::Min(ymin, ppD->GetBinContent(i)); }
      TH1F *fr = gPad->DrawFrame(ptReportMin, 0.3*ymin, inclEdges[NInclE-1], 5.*ymax);
      fr->GetXaxis()->SetTitle("#it{p}_{T}^{jet} [GeV]"); fr->GetYaxis()->SetTitle("jets per #it{N}_{Z} per GeV");
      fr->GetXaxis()->SetTitleSize(0.050); fr->GetYaxis()->SetTitleSize(0.050);
      fr->GetXaxis()->SetLabelSize(0.045); fr->GetYaxis()->SetLabelSize(0.045); fr->GetYaxis()->SetTitleOffset(1.5);
      for(int i = 1; i <= v->GetNbinsX(); i++){
        double y = v->GetBinContent(i); if(y <= 0. || v->GetBinLowEdge(i) < ptReportMin - 1e-6) continue;
        TBox *bx = new TBox(v->GetBinLowEdge(i), y*(1. - to->GetBinContent(i)), v->GetBinLowEdge(i+1), y*(1. + to->GetBinContent(i)));
        bx->SetFillColorAlpha(col, 0.30); bx->SetLineColor(col); bx->Draw();
      }
      if(q > 0){
        TH1D *pr = (TH1D*) ppD->Clone(Form("spec_ppRef_%d", q)); pr->SetDirectory(nullptr);
        for(int i = 1; i <= pr->GetNbinsX(); i++) if(pr->GetBinLowEdge(i) < ptReportMin - 1e-6){ pr->SetBinContent(i, -999.); pr->SetBinError(i, 0.); }
        styleH(pr, "#999999", markOpenSquare, 1.0); pr->Draw("E1 X0 same");
      }
      for(int i = 1; i <= v->GetNbinsX(); i++) if(v->GetBinLowEdge(i) < ptReportMin - 1e-6 || v->GetBinContent(i) <= 0.){ v->SetBinContent(i, -999.); v->SetBinError(i, 0.); }
      styleH(v, "#000000", markFilledCircle, 1.1); v->Draw("E1 X0 same");
      TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.050);
      la.DrawLatex(0.19, 0.935, Form("%s, %s", qLabel(q).Data(), q == 0 || modeMB ? "PF jets" : "calo-matched PF jets"));
    }
    cv->cd(6);
    TLegend *lg = makeLegend(0.08, 0.50, 0.95, 0.85, 0.055);
    TH1D *lv = (TH1D*) ppD->Clone("legSpec"); styleH(lv, "#000000", markFilledCircle, 1.1);
    TH1D *lp = (TH1D*) ppD->Clone("legSpecPP"); styleH(lp, "#999999", markOpenSquare, 1.0);
    TBox *lb = new TBox(); lb->SetFillColorAlpha(col, 0.30); lb->SetLineColor(col);
    lg->AddEntry(lv, "unfolded, stat.", "pe"); lg->AddEntry(lb, "total syst.", "f"); lg->AddEntry(lp, "pp, for reference", "pe"); lg->Draw();
    savePdfTight(cv, Form("%s/inclJetsPerZ_%s_spectra_withSys%s.pdf", figDir.Data(), figTag.Data(), sfx.Data()));
    delete cv;
  }
  ratioFig(5, 4, 2, 2, "jets per #it{N}_{Z}: PbPb / pp", Form("inclJetsPerZ_%s_PbPbOverPP_withSys", figTag.Data()));
  ratioFig(9, 3, 3, 1, Form("jets per #it{N}_{Z}: class / %s", sysLabel[4] + 5), Form("inclJetsPerZ_%s_RCP_withSys", figTag.Data()));

  fo->cd();
  // N_Z per system (pp, C1..C4) as used above (efficiency corrected, nominal window):
  // the b-jet chain reads it in its slice modes (bJetChain_caloJets_PbPb.C cent "ATLAS")
  { TVectorD nzv(NS); for(int s = 0; s < NS; s++) nzv[s] = NZ(s, nomO); nzv.Write("NZ"); }
  TNamed("info", Form("PF inclusive jets per Z systematics, PbPb calo-matched PF / pp plain PF; JEU pp-PbPb %s; nominal N_iter pp/C1-C4 = %d/%d/%d/%d/%d (closure min MSE); "
                      "%d response toys; JER %s; fake nominal %s, variation %s; JEU %s; see the macro header.",
                      jeuPPCorr ? "correlated" : "uncorrelated", nIterNom[0], nIterNom[1], nIterNom[2], nIterNom[3], nIterNom[4], nToys, jerSmear::tagText(),
                      gSystem->BaseName(fFakeNom.Data()), useFakeAlt ? gSystem->BaseName(fFakeAlt.Data()) : "none", !flattenSources ? "NOT flattened, no source (MinBias only)" : noFlattenJEU ? "NOT flattened (bin by bin)" : "flattened if no pT trend")).Write();
  fo->Close();
  if(saveVariations){
    TString vp = outRoot; vp.ReplaceAll("inclJetsPerZSys_", "inclJetsPerZVariations_");
    TFile *fv = TFile::Open(vp, "recreate");
    for(int s = 0; s < NS; s++){ TH1D *c = (TH1D*) nom[s]->Clone(Form("h_nominal_%s", sysTag[s])); c->Write(); }
    for(auto h : savedVar) h->Write();
    for(auto h : savedResp) h->Write();
    for(auto h : savedToys) h->Write();
    TNamed("info", "jets per Z, native binning: h_nominal_<sys> and h_var_<source>_<k>_<sys> for every correlated source of inclusiveJetsPerZSys_PF.C (k = variation index)").Write();
    fv->Close(); printf("written %s\n", vp.Data());
  }
  printf("\nwritten %s\nfigures in %s\n", outRoot.Data(), figDir.Data());
}
