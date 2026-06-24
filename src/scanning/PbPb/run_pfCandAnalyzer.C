// Wrapper to compile PbPb_pfCandAnalyzer.C with ACLiC + FastJet support.
// Usage: root -l -q 'run_pfCandAnalyzer.C(1)'
// Requires the LCG view (or CMSSW) to be sourced so that fastjet-config is in PATH.

void run_pfCandAnalyzer(int group = 1){
  TString fjCxxFlags = gSystem->GetFromPipe("fastjet-config --cxxflags 2>/dev/null");
  TString fjLibs     = gSystem->GetFromPipe("fastjet-config --libs 2>/dev/null");

  if(fjCxxFlags.Length() > 0 && fjLibs.Length() > 0){
    gSystem->AddIncludePath(fjCxxFlags + " -DDO_FASTJET");
    gSystem->AddLinkedLibs(fjLibs);
  } else {
    const char *fjHome = gSystem->Getenv("FASTJET_HOME");
    if(fjHome && strlen(fjHome) > 0){
      gSystem->AddIncludePath(Form("-DDO_FASTJET -I%s/include", fjHome));
      gSystem->AddLinkedLibs(Form("-L%s/lib -lfastjet", fjHome));
    } else {
      printf("WARNING: fastjet-config not found and FASTJET_HOME not set — compiling without FastJet support\n");
    }
  }

  gROOT->LoadMacro("PbPb_pfCandAnalyzer.C+");
  gROOT->ProcessLine(Form("PbPb_pfCandAnalyzer(%d)", group));
}
