// Standalone entry point for PbPb_caloTowerAnalyzer.C.
//
// Build:  make -f Makefile.caloTowerAnalyzer
// Run:    ./caloTowerAnalyzer_PbPb <group> [fileList] [outputBaseDir]
//
// Why this exists rather than running the macro through ROOT -- the same reason
// main_pfCandAnalyzer.cc does:
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
//
// TWO OPTIONAL ARGUMENTS beyond the pfCandAnalyzer version. The tower analyzer
// takes a file-list and an output-directory override so it can be pointed at a
// single local forest for a smoke test. They default to "" here, which is
// exactly the production path, so a job that passes only <group> behaves
// identically to the pfCandAnalyzer one.

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

#include "PbPb_caloTowerAnalyzer.C"

int main(int argc, char **argv){

  if(argc < 2){
    fprintf(stderr, "usage: %s <group> [fileList] [outputBaseDir]\n", argv[0]);
    fprintf(stderr, "  <group>        1-based index into the input file list\n");
    fprintf(stderr, "  [fileList]     optional override of the list chosen from config_PbPb.h\n");
    fprintf(stderr, "  [outputBaseDir] optional override of the EOS output base directory\n");
    return 2;
  }

  char *endp = nullptr;
  long group = strtol(argv[1], &endp, 10);
  if(endp == argv[1] || *endp != '\0' || group < 1){
    fprintf(stderr, "ERROR: group must be a positive integer, got '%s'\n", argv[1]);
    return 2;
  }

  TString fileList  = (argc > 2) ? argv[2] : "";
  TString outputDir = (argc > 3) ? argv[3] : "";

  PbPb_caloTowerAnalyzer((int)group, fileList, outputDir);

  // The analyzer returns void and bails out with a plain "return" on a missing
  // background map, a forest with no tower tree, an unreadable file list and a
  // dozen other conditions. Without this check every one of those exits 0 and
  // condor records the job as successful having produced no output.
  if(!g_caloTowerScanCompletedOK){
    fprintf(stderr, "ERROR: PbPb_caloTowerAnalyzer(%ld) did not complete -- "
                    "no output written. See the messages above.\n", group);
    return 1;
  }

  return 0;
}
