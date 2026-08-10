// Standalone entry point for PbPb_pfCandAnalyzer.C.
//
// Build:  make -f Makefile.pfCandAnalyzer
// Run:    ./pfCandAnalyzer_PbPb <group>
//
// Why this exists rather than running the macro through ROOT:
//
// fastjet::ClusterSequence's constructor is a template defined in the header,
// so it is not exported by libfastjet.so -- every consumer instantiates its own
// weak copy. The LCG views ship Delphes, which embeds its own FastJet build and
// exports that instantiation from libDelphesDisplay.so, and Delphes' rootmap
// claims the fastjet:: namespace. When ROOT parses the ACLiC dictionary for the
// macro, class autoloading dlopens libDelphesDisplay.so, and the dynamic linker
// then binds the weak symbol to Delphes' copy. Its PseudoJet layout differs from
// the standalone build the macro was compiled against, so ClusterSequence's
// constructor segfaults inside SharedPtr<UserInfoBase>.
//
// A standalone binary has no dictionary and no rootmap lookup, so Delphes is
// never loaded and the collision cannot happen. It also removes the per-job
// ACLiC compile, which matters when submitting ~2000 jobs.
//
// The analyzer is a single translation unit with no ClassDef and no cling-only
// constructs, so including the .C directly is safe.

#include <cstdio>
#include <cstdlib>

#include "PbPb_pfCandAnalyzer.C"

int main(int argc, char **argv){

  if(argc < 2){
    fprintf(stderr, "usage: %s <group>\n", argv[0]);
    fprintf(stderr, "  <group> is the 1-based index into the input file list\n");
    return 2;
  }

  char *endp = nullptr;
  long group = strtol(argv[1], &endp, 10);
  if(endp == argv[1] || *endp != '\0' || group < 1){
    fprintf(stderr, "ERROR: group must be a positive integer, got '%s'\n", argv[1]);
    return 2;
  }

  PbPb_pfCandAnalyzer((int)group);

  return 0;
}
