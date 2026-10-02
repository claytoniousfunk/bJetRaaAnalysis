// Batch-safe wrapper for caloVsPFJetScale_scan.C. Submitted by
// condor_caloVsPFJetScale.py; one call per condor job.
//
// Usage: root -l -b -q 'run_caloVsPFJetScale_condor.C("list_0.txt", "out_0.root", false)'
//
// Same two differences from running the macro directly as
// src/scanning/PbPb/run_muonTagAndProbe_condor.C:
//
//   1. Private ACLiC build directory per job, so concurrent jobs do not write
//      the same .so next to the macro and corrupt one another.
//
//   2. Non-zero exit on failure. The scan returns early, writing nothing, when
//      its list is unreadable, the chain is empty or a jet tree is missing, and
//      ROOT would still exit 0. This checks the output for h_evtCount and exits
//      1 if it is missing, so condor reports the job as failed.

void run_caloVsPFJetScale_condor(const char *fileList, const char *output, bool isPP = false)
{
  const char *scratch = gSystem->Getenv("_CONDOR_SCRATCH_DIR");
  TString buildDir = (scratch && strlen(scratch) > 0)
                   ? TString::Format("%s/aclic", scratch)
                   : TString::Format("%s/aclic_cvpf_%d", gSystem->TempDirectory(), gSystem->GetPid());
  if(gSystem->mkdir(buildDir, kTRUE) != 0 && gSystem->AccessPathName(buildDir)){
    printf("ERROR: cannot create ACLiC build directory %s\n", buildDir.Data());
    gSystem->Exit(1);
  }
  gSystem->SetBuildDir(buildDir, kTRUE);
  printf("ACLiC build directory: %s\n", buildDir.Data());

  if(gROOT->LoadMacro("caloVsPFJetScale_scan.C+") != 0){
    printf("ERROR: compiling caloVsPFJetScale_scan.C failed\n");
    gSystem->Exit(1);
  }

  gROOT->ProcessLine(Form("caloVsPFJetScale_scan(\"%s\", \"%s\", %s)",
                          fileList, output, isPP ? "true" : "false"));

  TFile *f = TFile::Open(output);
  if(!f || f->IsZombie() || !f->Get("h_evtCount")){
    printf("ERROR: %s missing or has no h_evtCount -- scan did not complete\n", output);
    gSystem->Exit(1);
  }
  f->Close();
  printf("job complete: %s\n", output);
}
