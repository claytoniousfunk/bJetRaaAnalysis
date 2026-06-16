// Wrapper to compile HYDJET_pfCandAnalyzer.C with ACLiC + FastJet support.
// Usage: root -l -q 'run_pfCandAnalyzer.C(1)'

void run_pfCandAnalyzer(int group = 1){
  const char *fjHome = gSystem->Getenv("FASTJET_HOME");
  if(fjHome && strlen(fjHome) > 0){
    gSystem->AddIncludePath(Form("-DDO_FASTJET -I%s/include", fjHome));
    gSystem->AddLinkerOpts("-lfastjet");
  } else {
    printf("WARNING: FASTJET_HOME not set — compiling without FastJet support\n");
  }
  gROOT->LoadMacro("HYDJET_pfCandAnalyzer.C+");
  gROOT->ProcessLine(Form("HYDJET_pfCandAnalyzer(%d)", group));
}
