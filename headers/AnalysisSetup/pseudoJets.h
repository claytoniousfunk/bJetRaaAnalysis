#pragma once
// pfCandidateAnalysis variables

bool doEventMixing = true;

bool doFastJetClustering = true;      // true = run anti-kT R=0.4 on PF candidates via FastJet (requires -DDO_FASTJET at compile time)

bool skipSingleConstituentJets = false;
bool useDeltaPTMapsForBkgSub = false;

bool useGeoCorrForRCMap = true;

int N_mixedEventsInPool = 100;

double pseudoJetCandPt_min = 0.0;

double subleadingPFCandPt_min = 15.0;
