// Load FastJet 3.3.0 + fjcontrib 1.033 (the forest's versions) for ACLiC
// macros that include headers/functions/csJetEmulation.h. FJ330 overrides the
// install prefix (default ~/Programs/fastjet330).
//   root -l -b -q loadFastJet330.C 'validateCSJets.C+(...)'
void loadFastJet330()
{
  TString pre = gSystem->Getenv("FJ330") ? gSystem->Getenv("FJ330")
                                         : Form("%s/Programs/fastjet330", gSystem->HomeDirectory());
  gSystem->AddIncludePath(Form("-I%s/include", pre.Data()));
  for(const char *lib : {"libfastjet.so", "libfastjettools.so", "libfastjetcontribfragile.so"})
    if(gSystem->Load(Form("%s/lib/%s", pre.Data(), lib)) < 0)
      printf("loadFastJet330: could not load %s/lib/%s\n", pre.Data(), lib);
}
