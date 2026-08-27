#pragma once
// pfCandidateAnalysis variables

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

double pseudoJetCandPt_min = 0.0;

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
