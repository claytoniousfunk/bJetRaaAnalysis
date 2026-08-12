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
