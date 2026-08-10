// Batch-safe wrapper to compile PbPb_pfCandAnalyzer.C with ACLiC + FastJet support.
//
// Usage: root -l -b -q 'run_pfCandAnalyzer_condor.C(1)'
//
// Differs from run_pfCandAnalyzer.C in two ways that matter when many jobs run
// concurrently out of the same directory:
//
//   1. Private ACLiC build directory. The default build dir is next to the macro,
//      so N concurrent jobs all write PbPb_pfCandAnalyzer_C.so (+ .d + dict pcm)
//      to the same path and corrupt each other. Here each job builds inside its
//      own $_CONDOR_SCRATCH_DIR (falling back to a pid-tagged temp dir), which is
//      job-local and cleaned up automatically.
//
//   2. Hard failure when FastJet is missing. The interactive wrapper warns and
//      then compiles without -DDO_FASTJET, which produces a job that exits 0 and
//      writes empty FastJet histograms -- the worst possible outcome in a batch
//      submission. Here a missing FastJet aborts with a non-zero exit code so
//      condor reports the job as failed.
//
// Requires the LCG view (or CMSSW) to be sourced so that fastjet-config is in PATH.
// The generated condor scripts do this themselves; see condor_PbPb_pfCandAnalyzer.py.
//
// Note: every job recompiles the analyzer (~1 min). That is deliberate -- it keeps
// jobs independent. If the compile time ever dominates, build the .so once and
// gSystem->Load() it here instead of using "+".

void run_pfCandAnalyzer_condor(int group = 1){

  // --- private ACLiC build directory -------------------------------------
  const char *scratch = gSystem->Getenv("_CONDOR_SCRATCH_DIR");
  TString buildDir;
  if(scratch && strlen(scratch) > 0)
    buildDir = Form("%s/aclic", scratch);
  else
    buildDir = Form("%s/aclic_%d_%d", gSystem->TempDirectory(), gSystem->GetPid(), group);

  if(gSystem->mkdir(buildDir, kTRUE) != 0 && gSystem->AccessPathName(buildDir)){
    printf("ERROR: cannot create ACLiC build directory %s\n", buildDir.Data());
    gSystem->Exit(1);
  }
  gSystem->SetBuildDir(buildDir, kTRUE);
  printf("ACLiC build directory: %s\n", buildDir.Data());

  // --- FastJet flags ------------------------------------------------------
  TString fjCxxFlags = gSystem->GetFromPipe("fastjet-config --cxxflags 2>/dev/null");
  TString fjLibs     = gSystem->GetFromPipe("fastjet-config --libs 2>/dev/null");

  if(fjCxxFlags.Length() > 0 && fjLibs.Length() > 0){
    gSystem->AddIncludePath(fjCxxFlags + " -DDO_FASTJET");
    gSystem->AddLinkedLibs(fjLibs);
    printf("FastJet via fastjet-config: %s\n", fjCxxFlags.Data());
  }
  else{
    const char *fjHome = gSystem->Getenv("FASTJET_HOME");
    if(fjHome && strlen(fjHome) > 0){
      gSystem->AddIncludePath(Form("-DDO_FASTJET -I%s/include", fjHome));
      gSystem->AddLinkedLibs(Form("-L%s/lib -lfastjet", fjHome));
      printf("FastJet via FASTJET_HOME: %s\n", fjHome);
    }
    else{
      printf("ERROR: fastjet-config not found and FASTJET_HOME not set.\n");
      printf("       Refusing to compile without FastJet -- the job would run to\n");
      printf("       completion and write empty FastJet histograms.\n");
      printf("       Source the LCG view (or set FASTJET_HOME) before submitting.\n");
      gSystem->Exit(1);
    }
  }

  // --- compile and run ----------------------------------------------------
  if(gROOT->LoadMacro("PbPb_pfCandAnalyzer.C+") != 0){
    printf("ERROR: ACLiC failed to build PbPb_pfCandAnalyzer.C\n");
    gSystem->Exit(1);
  }

  gROOT->ProcessLine(Form("PbPb_pfCandAnalyzer(%d)", group));
}
