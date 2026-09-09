#pragma once
// Stamp a scan output with the code version and configuration that produced it.
//
// WHY. Several 2026-09-08 commits changed what existing histograms CONTAIN
// without changing their names -- the jet pT axis moved raw -> JEC, and the
// tagging muon moved last-constituent -> leading. Files from either side of
// that look identical by key list, so telling them apart meant sniffing for a
// marker histogram and inferring the rest. That is a workaround for output
// files carrying no record of their own origin. This is the fix.
//
// Read it back with:
//   root -l file.root
//   provenance->GetTitle()
// or from a macro:
//   TNamed *p = nullptr; f->GetObject("provenance", p);
//   if(p) printf("%s\n", p->GetTitle());
//
// A macro that must not mix generations can parse the gitHash line, or simply
// require that two files agree on it.
//
// Requires the config globals listed in the body -- same convention as
// applyJEU_JER and findRecoMuonTag. Include this after the config headers.
//
// The git query is best effort. Condor runs the scan with `cd` into the repo,
// so it normally succeeds; if git is missing or the code was copied out of the
// working tree, the field records "unavailable" rather than failing the job.
// A scan is never worth losing to a bookkeeping call.

#include "TFile.h"
#include "TNamed.h"
#include "TDatime.h"
#include <cstdio>
#include <string>

// one line of output from a shell command, trailing newline stripped
inline std::string shellLine(const char *cmd)
{
  std::string out;
  FILE *p = popen(cmd, "r");
  if(!p) return "unavailable";
  char buf[512];
  if(fgets(buf, sizeof(buf), p)) out = buf;
  pclose(p);
  while(!out.empty() && (out.back() == '\n' || out.back() == '\r')) out.pop_back();
  return out.empty() ? std::string("unavailable") : out;
}

inline void writeProvenance(TFile *wf)
{
  if(!wf) return;

  std::string hash  = shellLine("git rev-parse --short HEAD 2>/dev/null");
  std::string full  = shellLine("git rev-parse HEAD 2>/dev/null");
  std::string brnch = shellLine("git rev-parse --abbrev-ref HEAD 2>/dev/null");

  // A hash alone does not identify the code if the tree had uncommitted edits,
  // and scans are often launched from a dirty tree mid-iteration. Record it,
  // so a file can say "this is NOT reproducible from the hash".
  std::string dirtyCount = shellLine("git status --porcelain --untracked-files=no 2>/dev/null | wc -l");
  bool dirty = (dirtyCount != "unavailable" && dirtyCount != "0");

  TDatime dt;
  std::string s;
  char b[1024];

  snprintf(b, sizeof(b), "scan written    : %04d-%02d-%02d %02d:%02d:%02d\n",
           dt.GetYear(), dt.GetMonth(), dt.GetDay(),
           dt.GetHour(), dt.GetMinute(), dt.GetSecond());               s += b;
  snprintf(b, sizeof(b), "gitHash         : %s\n", hash.c_str());        s += b;
  snprintf(b, sizeof(b), "gitHashFull     : %s\n", full.c_str());        s += b;
  snprintf(b, sizeof(b), "gitBranch       : %s\n", brnch.c_str());       s += b;
  snprintf(b, sizeof(b), "workingTree     : %s\n",
           dirty ? "DIRTY - uncommitted changes, NOT reproducible from the hash"
                 : "clean");                                            s += b;

  s += "--- event mixing / FastJet ---\n";
  snprintf(b, sizeof(b), "doEventMixing               : %d\n", (int)doEventMixing);                 s += b;
  snprintf(b, sizeof(b), "N_mixedEventsInPool         : %d\n", N_mixedEventsInPool);                s += b;
  snprintf(b, sizeof(b), "N_fastJetMixedEventResamples: %d\n", N_fastJetMixedEventResamples);       s += b;
  snprintf(b, sizeof(b), "doFastJetClustering         : %d\n", (int)doFastJetClustering);           s += b;
  snprintf(b, sizeof(b), "doConstituentSubtraction    : %d\n", (int)doConstituentSubtraction);      s += b;
  snprintf(b, sizeof(b), "doSignalSelectedRC          : %d (jetPtCut %.1f, raw %d)\n",
           (int)doSignalSelectedRC, signalJetPtCut, (int)signalJetPtCutIsRaw);                      s += b;

  s += "--- jet collection ---\n";
  snprintf(b, sizeof(b), "useCaloJetsOverride         : %d\n", (int)useCaloJetsOverride);           s += b;
  snprintf(b, sizeof(b), "useFlowJetsOverride         : %d\n", (int)useFlowJetsOverride);           s += b;
  snprintf(b, sizeof(b), "jetTree                     : %s\n",
           useCaloJetsOverride ? "akPu4CaloJetAnalyzer/t"
         : useFlowJetsOverride ? "akFlowPuCs4PFJetAnalyzer/t"
                               : "akCs4PFJetAnalyzer/t");                                           s += b;

  s += "--- selection ---\n";
  snprintf(b, sizeof(b), "muPtCut / muPtMaxCut        : %.1f / %.1f\n", muPtCut, muPtMaxCut);       s += b;
  snprintf(b, sizeof(b), "epsilon_mm                  : %.2f\n", epsilon_mm);                       s += b;
  snprintf(b, sizeof(b), "etaMax / jetPtCut           : %.1f / %.1f\n", etaMax, jetPtCut);          s += b;
  snprintf(b, sizeof(b), "doJetTrkMaxFilter           : %d\n", (int)doJetTrkMaxFilter);             s += b;
  snprintf(b, sizeof(b), "doWDecayFilter              : %d\n", (int)doWDecayFilter);                s += b;
  snprintf(b, sizeof(b), "doEtaPhiMask                : %d\n", (int)doEtaPhiMask);                  s += b;
  snprintf(b, sizeof(b), "applyMu12TrigEffCorrection  : %d\n", (int)applyMu12TriggerEfficiencyCorrection); s += b;

  s += "--- corrections ---\n";
  snprintf(b, sizeof(b), "apply_JEU_shift_up / down   : %d / %d\n",
           (int)apply_JEU_shift_up, (int)apply_JEU_shift_down);                                     s += b;
  snprintf(b, sizeof(b), "apply_JER_smear             : %d\n", (int)apply_JER_smear);               s += b;
  snprintf(b, sizeof(b), "doJERCorrection             : %d\n", (int)doJERCorrection);               s += b;

  s += "--- binning ---\n";
  snprintf(b, sizeof(b), "NCentralityIndices          : %d\n", NCentralityIndices);                 s += b;

  wf->cd();
  TNamed prov("provenance", s.c_str());
  prov.Write();

  printf("\n=== provenance written to output ===\n%s\n", s.c_str());
}
