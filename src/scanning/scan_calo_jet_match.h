#pragma once
// CALO-JET MATCHING as a jet ID for PF jets.
//
// A PF jet is calo-matched when a calo jet (akPu4Calo in PbPb, ak4Calo in pp)
// lies within caloMatchDr of its axis. Combinatorial PF jets are built from
// the UE pedestal, which the calo collection's own pileup subtraction removes
// differently, so a genuine jet should appear in both collections and a fake
// one less often. No pT or eta cut is applied to the calo jet beyond the
// forest's own threshold; raw pT is what the diagnostics keep, because the
// forest calo jtpt carries the PF JEC.
//
// With requireCaloJetMatch (config_PbPb.h, config_PYTHIAHYDJET.h, config_pp.h,
// config_PYTHIA.h) the scans drop every unmatched PF jet at the top of the
// reco-jet loop, so the requirement reaches every reco-jet histogram, the muon
// tag and the response alike; the output name gains _caloJetMatched.
// Gen-level histograms and the event counts (h_vz ...) are untouched. With the
// flag off the tree is still read when present, for the dR / calo-pT
// diagnostics only.
//
// Read as a plain TTree with its own buffers, NOT through eventMap: loadJet()
// attaches the PF jet tree as a friend of evtTree, and the calo tree shares its
// branch names (nref, rawpt, jteta, jtphi). Entry-aligned with evtTree like
// every other forest tree: call g_caloTree->GetEntry(evi) via
// loadCaloJetMatchEntry(evi) next to em->getEvent(evi).
//
// Needs eventMap.h (eventMap::jetMax) and getDr.h included first.

#include <iostream>
#include "TFile.h"
#include "TTree.h"

const double caloMatchDr = 0.2;
const char  *caloJetTreeName   = "akPu4CaloJetAnalyzer/t";   // PbPb forests
const char  *caloJetTreeNamePP = "ak4CaloJetAnalyzer/t";     // pp forests

static TTree  *g_caloTree = nullptr;
static Int_t   g_caloN    = 0;
static Float_t g_caloRawPt[eventMap::jetMax], g_caloEta[eventMap::jetMax], g_caloPhi[eventMap::jetMax];

// Attach the calo tree for a PF scan. Returns false if the tree is missing or
// misaligned; the caller must abort when requireCaloJetMatch is set, since a
// scan with no calo jets would otherwise silently keep every PF jet.
inline bool attachCaloJetMatchTree(TFile *f, Long64_t nEvents, const char *treeName = caloJetTreeName)
{
  g_caloTree = nullptr;
  g_caloN = 0;
  TTree *t = (TTree*) f->Get(treeName);
  if(!t){
    std::cout << "\033[1;31m " << treeName << " not found in " << f->GetName() << "\033[0m" << std::endl;
    return false;
  }
  if(t->GetEntries() != nEvents){
    std::cout << "\033[1;31m " << treeName << " has " << t->GetEntries()
              << " entries against " << nEvents << " events\033[0m" << std::endl;
    return false;
  }
  t->SetBranchStatus("*",0);
  for(const char *b : {"nref", "rawpt", "jteta", "jtphi"}) t->SetBranchStatus(b,1);
  t->SetBranchAddress("nref",  &g_caloN);
  t->SetBranchAddress("rawpt",  g_caloRawPt);
  t->SetBranchAddress("jteta",  g_caloEta);
  t->SetBranchAddress("jtphi",  g_caloPhi);
  g_caloTree = t;
  std::cout << "\tLoaded " << treeName << " for the calo-match jet ID (dR < " << caloMatchDr << ")" << std::endl;
  return true;
}

inline void loadCaloJetMatchEntry(Long64_t evi)
{
  g_caloN = 0;
  if(g_caloTree) g_caloTree->GetEntry(evi);
}

// dR to the closest calo jet (9.9 when the event has none, i.e. overflow in
// the diagnostics) and that calo jet's raw pT (-1 if none).
inline double closestCaloJetDr(double eta, double phi, double &caloRawPt)
{
  double drMin = 9.9;
  caloRawPt = -1.0;
  for(int k = 0; k < g_caloN; k++){
    double dr = getDr(eta, phi, g_caloEta[k], g_caloPhi[k]);
    if(dr < drMin){ drMin = dr; caloRawPt = g_caloRawPt[k]; }
  }
  return drMin;
}
