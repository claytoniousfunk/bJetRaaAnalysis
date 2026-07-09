#pragma once
// pfCandidateAnalysis variables
bool doEventMixing = false;
bool doConstituentSubtraction = false; // true = pfcandAnalyzerCS (CS applied at particle level), false = pfcandAnalyzer + manual rho subtraction
bool doFastJetClustering = true;      // true = run anti-kT R=0.4 on PF candidates via FastJet (requires -DDO_FASTJET at compile time)

int N_mixedEventsInPool = 100;

double pseudoJetCandPt_min = 0.0;

double subleadingPFCandPt_min = 15.0;
