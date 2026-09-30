#pragma once
// pfCandidateAnalysis variables

#include <cmath>
#include <string>

// true  = mixed-event: FastJet inputs, random cones and the tower-level PU
//         subtraction (doTowerPUSub) are built from a pool of
//         N_mixedEventsInPool same-centrality events -- the combinatorial
//         (fake) jet estimate. Slow.
// false = same-event: each event's own towers; no other event is touched.
// The output name carries _mixedEventPFClustering / _sameEventPFClustering.
//
// 2026-09-29 same-event scan (tower map in hand); 2026-09-30 back to
// mixed-event for the fake-jet rate of the tower-subtracted calo jets.
bool doEventMixing = true;

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

double pseudoJetCandPt_min = 0.0;

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
bool caloConstituentsOnly = true;
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
// MAP-MAKING OVERRIDE. Both the calo-constituent PF scan and, since
// 2026-09-28, PbPb_caloTowerAnalyzer.C are in the same position: no background
// map built from their own inputs exists yet, and the scan refuses to start
// without one. So the first pass subtracts the all-PF map with the same
// candidate cut (0 GeV). What that pass is FOR is the random-cone map it
// writes, which is built from its own clustering and is unaffected by which map
// was subtracted; its own RC- and dPT-subtracted spectra are not usable and
// should not be read. Merge that pass's h_randConeEtaPhi_* into a map file,
// point this override at it, and run again for the subtracted spectra.
//
// The tower map cannot reuse the PF one for real: measured on the first tower
// run, the 0-5% random cone holds 54.9 GeV of tower ET against ~91 GeV of PF
// candidates -- the 0.3 GeV tower threshold and the calorimeter's blindness to
// soft particles are most of the difference.
//
// 2026-09-29: the tower bootstrap pass (doBkgMapMakingPass, 2026-09-28) is
// complete -- PbPb_caloTowerAnalyzer.C, towers, etMin 0.3, etaCluster 3.0,
// mu12/pTmu-15to999/tight/jetTrkMaxFilter/WDecayFilter, MinBias Part1, all 19
// ultraFine slices present with 51-55M cones each (no zero-filled slices, see
// the 2026-09-19 C17/C18 gap above). Its own h_randConeEtaPhi_* IS the map --
// the file needs no separate merge step, it already covers every slice.
// bkgMapFileOverride below now points at that file (copy it to EOS under this
// same name), and doBkgMapMakingPass is turned off so subsequent calo-tower
// scans read it back and produce real RC-subtracted spectra instead of the
// bootstrap's zero-subtracted placeholders.
//
// SCOPE: this override is gated on caloConstituentsOnly, so it applies to the
// calo-constituent PF scan and PbPb_caloTowerAnalyzer.C ONLY. A plain PF scan
// (caloConstituentsOnly = false) still falls through to the dir+stem lookup
// below and is unaffected -- a towers-built map is not a valid background for
// full PF-candidate clustering (see the 54.9 vs 91 GeV gap above), so it must
// not silently become the default for every scan the way the old PF-map
// override effectively did.
std::string bkgMapFileOverride = "/eos/cms/store/group/phys_heavyions/cbennett/maps/PbPb_MinBias_Part1_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_mixedEventPFClustering_fastJetResamples-100_pseudoJetCandPtMin-0.0_towers_etMin-0.30_etaCluster-3.0_bkgMapMaking_2026-9-28_ultraFineCentBins.root";
std::string bkgMapFileUsed     = "";   // set by the scan; written to provenance

inline std::string bkgMapFile(double candPtMin)
{
  if(!bkgMapFileOverride.empty() && caloConstituentsOnly) return bkgMapFileOverride;
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
// BOOTSTRAP PASS. The scan refuses to start unless the background map contains
// every centrality slice, because the FastJet loop dereferences those maps
// without null checks. That is the right default, but it makes the maps
// impossible to bootstrap: no map can exist for a new clustering input (calo
// constituents, towers) until a pass has been run to write one, and that pass
// cannot start without a map.
//
// It also bit the whole ultraFine scheme on 2026-09-19, when commit 53e4effb
// extended NCentralityIndices from 17 to 19 (adding C17 80-85% and C18 85-90%).
// Every map on EOS was built before that, so every one of them is now missing
// C17 and C18 and NO complete map exists for any input.
//
// With this true, a slice absent from the map file is replaced by a zero
// TProfile2D with the scan's own RC binning. Nothing is dereferenced null, and
// every subtraction subtracts zero -- so the RC- and dPT-subtracted spectra
// come out equal to the unsubtracted ones and are NOT results. What the pass is
// for is the h_randConeEtaPhi_* it writes, which is built from its own
// clustering and does not depend on what was subtracted.
//
// The output name carries "_bkgMapMaking" so such a file can never be mistaken
// for a subtracted one.
//
// 2026-09-29: off. The tower bootstrap pass this produced is merged into
// bkgMapFileOverride above, so PbPb_caloTowerAnalyzer.C now has a real map for
// every slice and should run its normal (subtracted) pass.
bool doBkgMapMakingPass = false;

inline std::string bkgMapTag()
{
  if(doBkgMapMakingPass) return "_bkgMapMaking";
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
