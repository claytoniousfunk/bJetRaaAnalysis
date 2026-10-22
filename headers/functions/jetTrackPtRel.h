#pragma once
// Tracks for the track-ptRel scans (data-driven light template, ATLAS-style,
// arXiv:2204.13530): src/scanning/PbPb/PbPb_jetTrkMax_scan.C and
// src/scanning/pp/pp_jetTrkPtRel_scan.C. 2026-10-19, moved out of the PbPb scan so
// both systems select and match tracks the same way.
//
// Source, decided per file by openTracks() and reported by trackSource():
//   the forest's track tree (PbPbTracks/trackTree, ppTrack/trackTree) if present,
//   with the cuts below; otherwise the charged hadrons (pfId == 1) of
//   pfcandAnalyzer/pfTree. PF has already cleaned those (track quality, and
//   high-pT tracks without the calorimeter energy to match are not called charged
//   hadrons), but none of the cuts below can be applied to them.
// Cuts (the 2018 PbPb charged-hadron selection): highPurity, sigma(pT)/pT < 0.1,
// |dz/sigma| < 3, |dxy/sigma| < 3, nHit >= 11, chi2/ndof/nLayer < 0.18, and above
// 20 GeV calorimeter matching, (pfEcal + pfHcal)/cosh(eta) > 0.5 pT. A cut whose
// branch is missing is dropped with a warning and left out of cutsApplied().
//
// Requires eventMap.h and getDr.h to be included first (neither has an include guard).
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TMath.h"
#include "TString.h"
#include <vector>
#include <string>
#include <cstring>
#include <cstdio>

namespace jtrk {
const int trkMax = 100000;
int nTrk = 0;
Float_t trkPt[trkMax], trkEta[trkMax], trkPhi[trkMax], trkPtError[trkMax];
Float_t trkDxy1[trkMax], trkDxyError1[trkMax], trkDz1[trkMax], trkDzError1[trkMax];
Float_t trkChi2[trkMax], pfEcal[trkMax], pfHcal[trkMax];
UChar_t trkNHit[trkMax], trkNlayer[trkMax], trkNdof[trkMax];
Bool_t highPurity[trkMax];
bool hasPtErr = false, hasIP = false, hasHits = false, hasChi2 = false, hasCalo = false, hasHP = false;

const char *pfDirName = "pfcandAnalyzer";   // fallback: its pfTree, pfId == 1
bool usePF = false;
TString trkTreeUsed = "";
TTree *trkTree = nullptr;

const double maxPtErrFrac = 0.1, maxIPSig = 3., maxChi2PerLayer = 0.18, caloMatchPtMin = 20., caloMatchFrac = 0.5;
const int minNHit = 11;

// enable and address one branch; false if the forest lacks it
inline bool bind(TTree *t, const char *name, void *addr){
  if(!t->GetBranch(name)){ printf("WARNING: %s has no branch %s -- that cut is dropped\n", trkTreeUsed.Data(), name); return false; }
  t->SetBranchStatus(name, 1); t->SetBranchAddress(name, addr); return true;
}

// returns "" on success, else the reason the tree cannot be read
inline std::string loadTracks(TTree *t){
  t->SetBranchStatus("*", 0);
  for(const char *b : {"nTrk", "trkPt", "trkEta", "trkPhi"}) if(!t->GetBranch(b)) return std::string("no branch ") + b;
  if(strlen(t->GetBranch("trkPt")->GetClassName()) > 0) return std::string("trkPt is a ") + t->GetBranch("trkPt")->GetClassName() + ", expected a C array";
  bind(t, "nTrk", &nTrk); bind(t, "trkPt", trkPt); bind(t, "trkEta", trkEta); bind(t, "trkPhi", trkPhi);
  hasHP = bind(t, "highPurity", highPurity);
  hasPtErr = bind(t, "trkPtError", trkPtError);
  hasIP = bind(t, "trkDxy1", trkDxy1) & bind(t, "trkDxyError1", trkDxyError1) & bind(t, "trkDz1", trkDz1) & bind(t, "trkDzError1", trkDzError1);
  hasHits = bind(t, "trkNHit", trkNHit);
  hasChi2 = bind(t, "trkChi2", trkChi2) & bind(t, "trkNdof", trkNdof) & bind(t, "trkNlayer", trkNlayer);
  hasCalo = bind(t, "pfEcal", pfEcal) & bind(t, "pfHcal", pfHcal);
  return "";
}

// Pick the source for this file: the track tree if it exists, else charged PF
// candidates (loaded into em as a friend of its event tree, read by getEvent).
// Returns "" on success, else the reason no source is usable.
inline std::string openTracks(TFile *f, eventMap *em, const char *trkTreeName, Long64_t nEvents){
  trkTree = nullptr; usePF = false; trkTreeUsed = trkTreeName;
  f->GetObject(trkTreeName, trkTree);
  if(trkTree){
    if(trkTree->GetEntries() != nEvents) return Form("%s has %lld entries, the event tree %lld", trkTreeName, trkTree->GetEntries(), nEvents);
    const std::string e = loadTracks(trkTree);
    return e.empty() ? e : std::string(trkTreeName) + ": " + e;
  }
  if(!f->Get(Form("%s/pfTree", pfDirName))) return Form("neither %s nor %s/pfTree", trkTreeName, pfDirName);
  printf("no %s: tracks = charged PF candidates of %s/pfTree\n", trkTreeName, pfDirName);
  em->loadParticleFlowAnalyzer(pfDirName);
  usePF = true;
  return "";
}

inline bool passTrack(int k){
  const double pt = trkPt[k];
  if(hasHP && !highPurity[k]) return false;
  if(hasPtErr && trkPtError[k] > maxPtErrFrac*pt) return false;
  if(hasIP){
    if(trkDzError1[k] <= 0. || trkDxyError1[k] <= 0.) return false;
    if(TMath::Abs(trkDz1[k]/trkDzError1[k]) > maxIPSig || TMath::Abs(trkDxy1[k]/trkDxyError1[k]) > maxIPSig) return false;
  }
  if(hasHits && trkNHit[k] < minNHit) return false;
  if(hasChi2){
    if(trkNdof[k] == 0 || trkNlayer[k] == 0) return false;
    if(trkChi2[k]/trkNdof[k]/trkNlayer[k] > maxChi2PerLayer) return false;
  }
  if(hasCalo && pt > caloMatchPtMin && (pfEcal[k] + pfHcal[k])/TMath::CosH(trkEta[k]) < caloMatchFrac*pt) return false;
  return true;
}

inline std::string cutsApplied(){
  if(usePF) return std::string("none: charged PF candidates (pfId == 1) of ") + pfDirName + "/pfTree";
  std::string s = "trkPt/Eta/Phi";
  if(hasHP) s += ", highPurity";
  if(hasPtErr) s += ", ptErr/pt";
  if(hasIP) s += ", IP significance";
  if(hasHits) s += ", nHit";
  if(hasChi2) s += ", chi2/ndof/nLayer";
  if(hasCalo) s += ", calo match > 20 GeV";
  return s;
}

inline TString trackSource(){ return usePF ? TString::Format("%s/pfTree, pfId == 1", pfDirName) : trkTreeUsed; }

// This event's selected tracks with pT >= ptFloor, as indices into the arrays.
// Call after em->getEvent(evi). A leading track below the floor would fill
// nothing, so prefiltering does not change which track leads. Returns false if
// the event has more tracks than the arrays hold.
inline bool selectTracks(eventMap *em, Long64_t evi, double ptFloor, std::vector<int> &sel){
  if(usePF){   // charged PF candidates into the track arrays; no cuts apply to them
    nTrk = 0;
    const int nPF = em->pfPt ? (int) em->pfPt->size() : 0;
    for(int l = 0; l < nPF && nTrk < trkMax; l++){
      if(em->pfId->at(l) != 1) continue;
      trkPt[nTrk] = em->pfPt->at(l); trkEta[nTrk] = em->pfEta->at(l); trkPhi[nTrk] = em->pfPhi->at(l);
      nTrk++;
    }
  }
  else trkTree->GetEntry(evi);
  if(nTrk > trkMax){ printf("ERROR: nTrk %d > %d\n", nTrk, trkMax); return false; }
  sel.clear();
  for(int k = 0; k < nTrk; k++) if(trkPt[k] >= ptFloor && passTrack(k)) sel.push_back(k);
  return true;
}

// leading selected track within coneR of (eta, phi); -1 if none
inline int leading(const std::vector<int> &sel, double eta, double phi, double coneR){
  int best = -1;
  for(int k : sel){
    if(getDr(trkEta[k], trkPhi[k], eta, phi) >= coneR) continue;
    if(best < 0 || trkPt[k] > trkPt[best]) best = k;
  }
  return best;
}
}
