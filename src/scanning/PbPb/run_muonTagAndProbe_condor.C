// Batch-safe wrapper for PbPb_muonTagAndProbe_scan.C. Submitted by
// condor_PbPb_muonTagAndProbe.py; one call per condor job.
//
// Usage: root -l -b -q 'run_muonTagAndProbe_condor.C("list_0.txt", "out_0.root", false)'
//
// Differs from running the macro directly in two ways that matter when many
// jobs run at once out of the same directory:
//
//   1. Private ACLiC build directory. The default puts the .so (+ .d + dict pcm)
//      next to the macro, so concurrent jobs write the same path and corrupt one
//      another. Here each job builds in its own $_CONDOR_SCRATCH_DIR.
//
//   2. Non-zero exit on failure. The scan returns early, writing nothing, when
//      its list is unreadable or the chain is empty, and ROOT would still exit 0.
//      This checks the output for the tnp tree and exits 1 if it is missing, so
//      condor reports the job as failed instead of done.

void run_muonTagAndProbe_condor(const char *fileList, const char *output, bool isMC = false)
{
  const char *scratch = gSystem->Getenv("_CONDOR_SCRATCH_DIR");
  TString buildDir = (scratch && strlen(scratch) > 0)
                   ? TString::Format("%s/aclic", scratch)
                   : TString::Format("%s/aclic_tnp_%d", gSystem->TempDirectory(), gSystem->GetPid());
  if(gSystem->mkdir(buildDir, kTRUE) != 0 && gSystem->AccessPathName(buildDir)){
    printf("ERROR: cannot create ACLiC build directory %s\n", buildDir.Data());
    gSystem->Exit(1);
  }
  gSystem->SetBuildDir(buildDir, kTRUE);
  printf("ACLiC build directory: %s\n", buildDir.Data());

  if(gROOT->LoadMacro("PbPb_muonTagAndProbe_scan.C+") != 0){
    printf("ERROR: compiling PbPb_muonTagAndProbe_scan.C failed\n");
    gSystem->Exit(1);
  }

  gROOT->ProcessLine(Form("PbPb_muonTagAndProbe_scan(\"%s\", \"%s\", %s)",
                          fileList, output, isMC ? "true" : "false"));

  TFile *f = TFile::Open(output);
  if(!f || f->IsZombie() || !f->Get("tnp")){
    printf("ERROR: %s missing or has no tnp tree -- scan did not complete\n", output);
    gSystem->Exit(1);
  }
  f->Close();
  printf("job complete: %s\n", output);
}
