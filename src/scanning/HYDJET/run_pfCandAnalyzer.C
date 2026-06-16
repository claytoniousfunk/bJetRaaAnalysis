// Wrapper to compile HYDJET_pfCandAnalyzer.C with ACLiC + FastJet support.
// Usage: root -l -q 'run_pfCandAnalyzer.C(1)'
// Requires the LCG view (or CMSSW) to be sourced so that FastJet headers
// and libfastjet.so are on the compiler/linker search paths.

void run_pfCandAnalyzer(int group = 1){
  // If FASTJET_HOME is set, add its include path explicitly.
  // If not (e.g. LCG view puts FastJet on the default search path),
  // just define DO_FASTJET — the headers are already findable.
  const char *fjHome = gSystem->Getenv("FASTJET_HOME");
  if(fjHome && strlen(fjHome) > 0){
    gSystem->AddIncludePath(Form("-DDO_FASTJET -I%s/include", fjHome));
  } else {
    gSystem->AddIncludePath("-DDO_FASTJET");
  }
  gSystem->AddLinkedLibs("-lfastjet");

  gROOT->LoadMacro("HYDJET_pfCandAnalyzer.C+");
  gROOT->ProcessLine(Form("HYDJET_pfCandAnalyzer(%d)", group));
}
