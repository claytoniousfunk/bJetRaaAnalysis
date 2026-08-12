// Standalone entry point for HYDJET_pfCandAnalyzer.C.
//
// Build:  make -f Makefile.pfCandAnalyzer
// Run:    ./pfCandAnalyzer_HYDJET <group>
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
// ACLiC compile.
//
// This mirrors src/scanning/PbPb/main_pfCandAnalyzer.cc; keep the two in step.

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

// rootcling's generated dictionary emits a `using namespace std;` before it
// includes the macro, so ACLiC builds get it for free. The analysis headers
// depend on that: eventMap.h declares bare `vector<Float_t>*`, common.h a bare
// `string`, and several headers use unqualified cout/endl. Replicate it here so
// a plain g++ translation unit sees the same names.
namespace std {}
using namespace std;

#include "HYDJET_pfCandAnalyzer.C"

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

  HYDJET_pfCandAnalyzer((int)group);

  return 0;
}
