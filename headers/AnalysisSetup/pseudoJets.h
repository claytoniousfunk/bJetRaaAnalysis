#pragma once
// pfCandidateAnalysis variables

#include <cmath>
#include <string>

bool doEventMixing = false;

bool doFastJetClustering = true;      // true = run anti-kT R=0.4 on PF candidates via FastJet (requires -DDO_FASTJET at compile time)

// PF candidate collection fed to the random cones and the FastJet clustering.
//   false -> pfcandAnalyzer   (raw PF), matching PbPb_pfCandAnalyzer.C
//   true  -> pfcandAnalyzerCS (constituent-subtracted)
// Must be false for the HYDJET fake-jet closure test: the data-side estimate is
// built from raw PF, and a CS cone carries roughly a tenth of the energy of a
// raw one, so the two are not comparable.
bool doConstituentSubtraction = false;

bool skipSingleConstituentJets = false;
bool useDeltaPTMapsForBkgSub = false;

bool useGeoCorrForRCMap = false;

int N_mixedEventsInPool = 100;

// The mixed-event candidate pool (N_mixedEventsInPool events' worth of PF
// candidates) is built once per event. This controls how many INDEPENDENT
// random draws of NCandidatesToSample are then taken from that same pool for
// the FastJet clustering step -- each draw is reclustered from scratch and
// filled separately, so this multiplies the statistics of every FastJet
// mixed-event histogram (h_fastJetPt_PF*, h_fastJetMuonPtRel_fastJetPt_PF_
// bkgSub_RC*, h_muonDR_inclusiveClosestJet, the PFCs-matching histograms,
// ...) without requiring more real events to be scanned.
//
// Each draw is filled with weight w/N_fastJetMixedEventResamples, so the
// TOTAL contribution of one real event to any of these histograms is
// unchanged no matter what this is set to -- increasing it only reduces the
// pool-sampling noise (smoother fake-muon/fake-jet shapes at the same
// normalization). No downstream macro needs to change: still just divide by
// h_vz for a per-event rate, exactly as today. This is deliberately DIFFERENT
// from the existing random-cone convention (h_pseudoJetPt etc.), which fills
// every one of its N_mixedEventsInPool cone throws at the FULL weight w and
// requires dividing by N_pool downstream (see makeFakeJetFile.C) -- that
// asymmetry is intentional: the random-cone histograms were designed that
// way already and changing their convention would break every macro that
// reads them, whereas this is a brand-new knob with no consumers yet, so it
// can default to the answer that needs no new bookkeeping.
//
// Has no effect when !doEventMixing: same-event clustering has no randomness
// to resample (it is just the event's own candidates every time), so the
// loop runs exactly once regardless of this value.
int N_fastJetMixedEventResamples = 100;

double pseudoJetCandPt_min = 2.0;

// ---- calo-jet fake-jet study ------------------------------------------------
// caloConstituentsOnly restricts every PF-candidate sum this scan makes -- the
// FastJet inputs (same- and mixed-event), the random cones, the PFCs
// clustering and the injection donors -- to the species a calorimeter jet is
// built from, so the combinatorial-jet estimate is the calo-jet analogue of
// the PF one. It changes ONLY which candidates are clustered: the reco jet
// collection, and so the forest branches required, are exactly those of the
// standard PF scan. No calo-jet branch is needed -- leave useCaloJetsOverride
// off unless you separately want akPu4Calo as the reco jets.
//
// PF ids: 1 charged hadron, 2 electron, 3 muon, 4 photon, 5 neutral hadron,
// 6 HF hadron, 7 HF EM.
//   default                       : all but muons (1,2,4,5,6,7). A muon is a
//                                   MIP in the calorimeter and calo jets exclude
//                                   its momentum; everything else deposits its
//                                   energy in the towers.
//   caloConstituentsIncludeCharged = false
//                                 : photons, neutral hadrons and HF only -- the
//                                   purely calorimetric PF candidates, a lower
//                                   bound on what the calorimeter sees.
//
// CAVEAT. This selects WHICH particles enter; it does not reproduce HOW a
// calorimeter measures them. Charged hadrons enter at their track momentum,
// not the lower calorimeter response; there is no 0.087 tower granularity, no
// 0.3 GeV tower threshold, and no akPu pileup subtraction. With the default
// the selection removes only muons, so the result should sit close to the PF
// one -- the neutral-only mode brackets it from below.
//
// The existing background maps were built from ALL PF candidates, so
// bkgMapFile() returns none in this mode: the first (map-making) run needs
// bkgMapFileOverride, exactly as the 2 GeV maps were bootstrapped.
bool caloConstituentsOnly = false;
bool caloConstituentsIncludeCharged = true;

inline bool isClusteringCand(int pfId)
{
  if(!caloConstituentsOnly) return true;
  if(pfId == 3) return false;                                   // muon
  if(pfId == 1 || pfId == 2) return caloConstituentsIncludeCharged;
  return (pfId == 4 || pfId == 5 || pfId == 6 || pfId == 7);    // gamma, h0, HF
}

// output-name tag for the constituent selection ("" for the standard PF scan)
inline std::string constituentTag()
{
  if(!caloConstituentsOnly) return "";
  return caloConstituentsIncludeCharged ? "_caloConstituents" : "_caloNeutralConstituents";
}

// Background (UE) maps for the RC- and dPT-subtracted FastJet spectra
// (ultraFine centrality). A map is only valid for a scan run with the SAME
// PF-candidate pT cut it was built with: in 0-5% the random-cone map holds
// ~91 GeV per cone without a cut and ~18 GeV with the 2 GeV cut, so a
// mismatched map subtracts the wrong background by tens of GeV.
//
// bkgMapFile() returns the map built with a given cut; PbPb_pfCandAnalyzer.C
// stops if there is none. To use a map not in the table, set
// bkgMapFileOverride to its full path -- the output name then carries
// "_bkgMapOverride" instead of "_matchedBkgMap", and the provenance stamp
// records the file either way.
//
// Each entry is the same-event MinBias scan made with that cut, copied into
// the EOS maps directory under its own file name.
std::string bkgMapFileOverride = "";
std::string bkgMapFileUsed     = "";   // set by the scan; written to provenance

inline std::string bkgMapFile(double candPtMin)
{
  if(!bkgMapFileOverride.empty()) return bkgMapFileOverride;
  // every map below was built from all PF candidates; none is valid for a
  // calo-constituent scan (see caloConstituentsOnly above)
  if(caloConstituentsOnly) return "";
  const std::string dir  = "/eos/cms/store/group/phys_heavyions/cbennett/maps/";
  const std::string stem = "PbPb_MinBias_Part1_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_sameEventPFClustering_";
  if(fabs(candPtMin - 0.0) < 1e-6) return dir + stem + "pseudoJetCandPtMin-0.0_2026-8-17_ultraFineCentBins.root";
  if(fabs(candPtMin - 2.0) < 1e-6) return dir + stem + "pseudoJetCandPtMin-2.0_2026-9-15_ultraFineCentBins.root";
  return "";
}

// Output-name tag saying which map the subtracted spectra used. Needed because
// the 2026-09-15 2 GeV map-making pass subtracted the 0 GeV map, and a rerun
// with the matched map would otherwise get exactly the same name.
inline std::string bkgMapTag()
{
  return bkgMapFileOverride.empty() ? "_matchedBkgMap" : "_bkgMapOverride";
}

double subleadingPFCandPt_min = 15.0;

// Same-event random cone restricted to events that contain a reco jet clearing
// signalJetPtCut -- Olga's multiplicity-matching proposal (2026-08-13): throw
// cones in the literal signal-selected events instead of a separate reference
// sample (MinBias, mixed-event pool, ...), so there is no event-activity
// mismatch between the background estimate and the events it is meant to
// correct, and no reweighting is needed. Filled only when !doEventMixing --
// "does this event pass the signal selection" is meaningless for a candidate
// drawn from a pool of OTHER events.
//
// signalJetPtCutIsRaw: cut on raw jet pT, confirmed by Olga (2026-08-13) --
// not JEC-corrected.
//
// Still open: selection scope is currently just the jet pT cut in isolation,
// not the full downstream analysis chain (eta cut, jetTrkMax filter, muon
// tag, ...). Extend the event-level pre-scan in PbPb_pfCandAnalyzer.C if she
// wants those folded in too.
//
// Runs on PbPb MinBias (~2000 files, fileNames_HIMinimumBias0_Part1_
// withTracksAndPFCandidates.txt) via the existing condor_PbPb_pfCandAnalyzer.py
// setup -- HardProbes/SingleMuon were considered as a higher-statistics source
// of real signal jets, but neither forest carries a PF-candidate (or PFCs)
// branch, so MinBias is the only option despite the rarity of a genuine 40 GeV
// jet there.
bool   doSignalSelectedRC  = true;
double signalJetPtCut      = 40.0;
bool   signalJetPtCutIsRaw = true;
