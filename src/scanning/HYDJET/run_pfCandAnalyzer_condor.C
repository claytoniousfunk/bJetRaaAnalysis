// Batch-safe wrapper to compile HYDJET_pfCandAnalyzer.C with ACLiC + FastJet.
//
// Usage: root -l -b -q 'run_pfCandAnalyzer_condor.C(1)'
//
// KNOWN BROKEN on any LCG view that ships Delphes -- it segfaults inside
// fastjet::ClusterSequence's constructor, see main_pfCandAnalyzer.cc. Kept only
// for use on a stack without Delphes; prefer the standalone binary.
//
// Differs from run_pfCandAnalyzer.C in two ways that matter when many jobs run
// concurrently out of the same directory:
//
//   1. Private ACLiC build directory. The default build dir is next to the
//      macro, so N concurrent jobs all write HYDJET_pfCandAnalyzer_C.so to the
//      same path and corrupt each other. Here each job builds inside its own
//      $_CONDOR_SCRATCH_DIR.
//
//   2. Hard failure when FastJet is missing. The interactive wrapper warns and
//      then compiles without -DDO_FASTJET, producing a job that exits 0 and
//      writes empty FastJet histograms.
//
// Mirrors src/scanning/PbPb/run_pfCandAnalyzer_condor.C.

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
      gSystem->Exit(1);
    }
  }

  // --- preload the real libfastjet ----------------------------------------
  // Delphes bundles its own copy of the FastJet sources and its rootmap claims
  // the fastjet:: namespace. Loading the genuine library first puts its symbols
  // in the global scope ahead of Delphes'.
  TString fjSo = fjLibDir + "/libfastjet.so";
  if(gSystem->Load(fjSo) < 0){
    printf("ERROR: could not load %s\n", fjSo.Data());
    gSystem->Exit(1);
  }
  printf("Preloaded FastJet: %s\n", fjSo.Data());

  // --- compile and run ----------------------------------------------------
  if(gROOT->LoadMacro("HYDJET_pfCandAnalyzer.C+") != 0){
    printf("ERROR: ACLiC failed to build HYDJET_pfCandAnalyzer.C\n");
    gSystem->Exit(1);
  }

  printf("Loaded FastJet/Delphes libraries: %s\n",
         gSystem->GetLibraries("fastjet|Delphes", "D", kFALSE));

  gROOT->ProcessLine(Form("HYDJET_pfCandAnalyzer(%d)", group));
}
