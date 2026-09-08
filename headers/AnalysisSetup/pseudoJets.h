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
// normalisation). No downstream macro needs to change: still just divide by
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
