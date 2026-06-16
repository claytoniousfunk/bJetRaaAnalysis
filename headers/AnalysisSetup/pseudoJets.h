// pfCandidateAnalysis variables
bool doEventMixing = false;
bool doConstituentSubtraction = false; // true = pfcandAnalyzerCS (CS applied at particle level), false = pfcandAnalyzer + manual rho subtraction
bool doFastJetClustering = false;      // true = run anti-kT R=0.4 on PF candidates via FastJet

int N_generatedPseudoJets = 100;

double pseudoJetCandPt_min = 0.0;

double subleadingPFCandPt_min = 15.0;

// Path to FastJet shared library (only used when doFastJetClustering = true).
// Set to bare name if libfastjet.so is already on LD_LIBRARY_PATH (e.g. after
// sourcing an LCG view or CMSSW environment on LXPLUS), otherwise supply the
// full absolute path.
const char *fastjetLibPath = "libfastjet";
