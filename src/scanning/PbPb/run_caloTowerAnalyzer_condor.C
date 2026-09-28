// Batch-safe wrapper to compile PbPb_caloTowerAnalyzer.C with ACLiC + FastJet.
//
// Usage: root -l -b -q 'run_caloTowerAnalyzer_condor.C(1)'
//
// This is the FALLBACK path. condor_PbPb_caloTowerAnalyzer.py defaults to the
// standalone binary (caloTowerAnalyzer_PbPb) instead, because on any LCG view
// that ships Delphes the ACLiC route segfaults inside ClusterSequence's
// constructor -- see main_caloTowerAnalyzer.cc for the mechanism. Kept for use
// on a stack without Delphes, and because it is the quickest way to run a
// single job interactively without a make step.
//
// Differs from a plain interactive run in two ways that matter when many jobs
// run concurrently out of the same directory:
//
//   1. Private ACLiC build directory. The default build dir is next to the
//      macro, so N concurrent jobs all write PbPb_caloTowerAnalyzer_C.so (+ .d
//      + dict pcm) to the same path and corrupt each other. Here each job builds
//      inside its own $_CONDOR_SCRATCH_DIR (falling back to a pid-tagged temp
//      dir), which is job-local and cleaned up automatically.
//
//   2. Hard failure when FastJet is missing. Compiling without -DDO_FASTJET
//      produces a job that exits 0 and writes empty FastJet histograms -- the
//      worst possible outcome in a batch submission. Here a missing FastJet
//      aborts with a non-zero exit code so condor reports the job as failed.
//
// It also propagates the analyzer's own completion flag, so a job that bailed
// out early (missing background map, forest with no tower tree, unreadable file
// list) exits non-zero rather than reporting success having written nothing.
//
// Requires the LCG view (or CMSSW) to be sourced so that fastjet-config is in
// PATH. The generated condor scripts do this themselves; see
// condor_PbPb_caloTowerAnalyzer.py.
//
// Note: every job recompiles the analyzer (~1 min). That is deliberate -- it
// keeps jobs independent. If the compile time ever dominates, use the standalone
// binary path, which builds once.

void run_caloTowerAnalyzer_condor(int group = 1,
                                  const char *fileList      = "",
                                  const char *outputBaseDir = ""){

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
  TString fjLibDir   = "";

  if(fjCxxFlags.Length() > 0 && fjLibs.Length() > 0){
    fjLibDir = gSystem->GetFromPipe("fastjet-config --prefix 2>/dev/null");
    fjLibDir += "/lib";
    gSystem->AddIncludePath(fjCxxFlags + " -DDO_FASTJET");
    gSystem->AddLinkedLibs(fjLibs + Form(" -Wl,-rpath,%s", fjLibDir.Data()));
    printf("FastJet via fastjet-config: %s\n", fjCxxFlags.Data());
  }
  else{
    const char *fjHome = gSystem->Getenv("FASTJET_HOME");
    if(fjHome && strlen(fjHome) > 0){
      fjLibDir = Form("%s/lib", fjHome);
      gSystem->AddIncludePath(Form("-DDO_FASTJET -I%s/include", fjHome));
      gSystem->AddLinkedLibs(Form("-L%s/lib -lfastjet -Wl,-rpath,%s/lib", fjHome, fjHome));
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

  // --- preload the real libfastjet ----------------------------------------
  // Delphes bundles its own copy of the FastJet sources and its rootmap claims
  // the fastjet:: namespace. Loading the compiled macro triggers ROOT class
  // autoloading, which pulls in libDelphesDisplay.so; the dynamic linker then
  // resolves fastjet::ClusterSequence etc. to Delphes' embedded build while the
  // macro was compiled against the headers of the standalone one. The layouts
  // differ, so ClusterSequence's ctor segfaults inside SharedPtr<UserInfoBase>.
  // Loading the genuine library first puts its symbols in the global scope
  // ahead of Delphes'. If a site ever loads Delphes even earlier, fall back to
  // LD_PRELOAD in the job script (see condor_PbPb_caloTowerAnalyzer.py).
  TString fjSo = fjLibDir + "/libfastjet.so";
  if(gSystem->Load(fjSo) < 0){
    printf("ERROR: could not load %s\n", fjSo.Data());
    gSystem->Exit(1);
  }
  printf("Preloaded FastJet: %s\n", fjSo.Data());

  // --- compile and run ----------------------------------------------------
  if(gROOT->LoadMacro("PbPb_caloTowerAnalyzer.C+") != 0){
    printf("ERROR: ACLiC failed to build PbPb_caloTowerAnalyzer.C\n");
    gSystem->Exit(1);
  }

  // Report which libraries actually provide the FastJet symbols. If a Delphes
  // library shows up here, the preload lost the race -- use LD_PRELOAD.
  printf("Loaded FastJet/Delphes libraries: %s\n",
         gSystem->GetLibraries("fastjet|Delphes", "D", kFALSE));

  gROOT->ProcessLine(Form("PbPb_caloTowerAnalyzer(%d,\"%s\",\"%s\")",
                          group, fileList, outputBaseDir));

  // The analyzer returns void and bails out with a plain "return" on a missing
  // background map, a forest with no tower tree and a dozen other conditions.
  // Without this the job exits 0 and condor records success with no output.
  Long_t ok = gROOT->ProcessLine("g_caloTowerScanCompletedOK;");
  if(!ok){
    printf("ERROR: PbPb_caloTowerAnalyzer(%d) did not complete -- no output written.\n", group);
    gSystem->Exit(1);
  }
}
