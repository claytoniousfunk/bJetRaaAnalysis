#pragma once
// Calo jets borrow from the nearest PF jet of the same event what their own jet
// tree does not carry.
//
// The forests' calo jet trees have no hadron-level flavor information: no
// bHadronNumber, so no calo jet can be put in the gluon-splitting b class (17),
// and no jtPartonFlavor / matchedPartonFlavor. The PF jet tree of the same event
// has them. This reads the PF tree a second time and hands each calo jet the
// values of the nearest PF jet in dR.
//
// It cannot go through eventMap: loadJet() attaches the jet tree as a FRIEND of
// evtTree, so a second jet tree would collide on jtpt, jteta and the rest. It is
// read as a plain TTree with its own addresses and stepped with the same entry
// index -- both are per-event trees in the same forest file, so entry i is the
// same event. Call g_pfTree->GetEntry(evi) next to em->getEvent(evi).
//
// How far the match can be trusted, from plotCaloPFMatching_PYTHIA.C (pp) and
// plotCaloPFFlavorMatch_PYTHIAHYDJET.C (PbPb): above 50 GeV, 97-100% of
// calo->PF pairs share one gen jet in every class. At 30-50 GeV in 0-10% only
// ~60% of calo jets find a PF jet and ~40% of those share a gen jet; most of
// those calo jets have no gen jet at all.
//
// Needs eventMap.h (for eventMap::jetMax) included first.

#include <cmath>
#include "TFile.h"
#include "TTree.h"
#include "TVector2.h"

static TTree  *g_pfTree    = nullptr;
static bool    g_pfHasFlav = false;   // jtPartonFlavor attached
static bool    g_pfHasBHad = false;   // bHadronNumber attached
static Int_t   g_pfN       = 0;
static const int g_pfMax = eventMap::jetMax;   // class member, not a global
static Float_t g_pfPt [g_pfMax], g_pfEta[g_pfMax], g_pfPhi[g_pfMax], g_pfFlav[g_pfMax];
static Int_t   g_pfBHad[g_pfMax];
static long    g_nPFMatched   = 0, g_nPFUnmatched   = 0;   // flavor lookups
static long    g_nBHadMatched = 0, g_nBHadUnmatched = 0;   // bHadronNumber lookups

// Attach the PF tree for whichever of the two borrowings is wanted. Anything
// missing from the forest is reported and that borrowing is switched off, so a
// forest without it still scans with the calo jet's own values.
inline void attachCaloPFMatchTree(TFile *f, const char *pfTreeName, bool wantFlavor, bool wantBHad)
{
  g_pfTree = nullptr; g_pfHasFlav = false; g_pfHasBHad = false;
  if(!wantFlavor && !wantBHad) return;

  TTree *t = (TTree*) f->Get(pfTreeName);
  if(!t){
    cout << "\tWARNING: " << pfTreeName << " is absent; calo jets keep their own flavor and bHadronNumber = 0\n";
    return;
  }
  g_pfHasFlav = wantFlavor && t->GetBranch("jtPartonFlavor");
  g_pfHasBHad = wantBHad   && t->GetBranch("bHadronNumber");
  if(wantFlavor && !g_pfHasFlav)
    cout << "\tWARNING: " << pfTreeName << " has no jtPartonFlavor; calo jets keep their own flavor\n";
  if(wantBHad && !g_pfHasBHad)
    cout << "\tWARNING: " << pfTreeName << " has no bHadronNumber; calo jets keep bHadronNumber = 0\n";
  if(!g_pfHasFlav && !g_pfHasBHad) return;

  // read only what is used: GetEntry on the full PF tree every event is not free
  t->SetBranchStatus("*", 0);
  for(const char *b : {"nref", "jtpt", "jteta", "jtphi"}) t->SetBranchStatus(b, 1);
  t->SetBranchAddress("nref",  &g_pfN);
  t->SetBranchAddress("jtpt",   g_pfPt);
  t->SetBranchAddress("jteta",  g_pfEta);
  t->SetBranchAddress("jtphi",  g_pfPhi);
  if(g_pfHasFlav){ t->SetBranchStatus("jtPartonFlavor", 1); t->SetBranchAddress("jtPartonFlavor", g_pfFlav); }
  if(g_pfHasBHad){ t->SetBranchStatus("bHadronNumber",  1); t->SetBranchAddress("bHadronNumber",  g_pfBHad); }
  g_pfTree = t;
}

// index of the nearest PF jet within maxDR, or -1
inline int nearestPFJet(double eta, double phi, double maxDR)
{
  if(!g_pfTree) return -1;
  double best = maxDR; int bestIdx = -1;
  for(int j = 0; j < g_pfN; j++){
    double dEta = eta - g_pfEta[j];
    double dPhi = TVector2::Phi_mpi_pi(phi - g_pfPhi[j]);
    double dr   = sqrt(dEta*dEta + dPhi*dPhi);
    if(dr < best){ best = dr; bestIdx = j; }
  }
  return bestIdx;
}

// matched PF jet's jtPartonFlavor, or fallback (the calo jet's own label)
inline int caloFlavorByPFMatch(double eta, double phi, double maxDR, int fallback)
{
  if(!g_pfHasFlav) return fallback;
  int j = nearestPFJet(eta, phi, maxDR);
  if(j < 0){ g_nPFUnmatched++; return fallback; }
  g_nPFMatched++;
  return (int) g_pfFlav[j];
}

// matched PF jet's bHadronNumber, or 0 (never bGS) when there is no match
inline int caloBHadronNumberByPFMatch(double eta, double phi, double maxDR)
{
  if(!g_pfHasBHad) return 0;
  int j = nearestPFJet(eta, phi, maxDR);
  if(j < 0){ g_nBHadUnmatched++; return 0; }
  g_nBHadMatched++;
  return g_pfBHad[j];
}

inline void printCaloPFMatchSummary()
{
  if(g_pfHasFlav)
    cout << "\tcalo->PF flavor lookups: " << g_nPFMatched << " matched, "
         << g_nPFUnmatched << " kept the calo jet's own flavor\n";
  if(g_pfHasBHad)
    cout << "\tcalo->PF bHadronNumber lookups: " << g_nBHadMatched << " matched, "
         << g_nBHadUnmatched << " held at 0\n";
}
