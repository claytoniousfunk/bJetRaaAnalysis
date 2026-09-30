#pragma once
// CaloTower clustering knobs for PbPb_caloTowerAnalyzer.C.
//
// WHY THIS EXISTS. pseudoJets.h has a caloConstituentsOnly mode that
// approximates a calorimeter jet by restricting the PF-candidate clustering to
// the species a calo jet is built from, and its own comment names what that
// cannot do:
//
//   "This selects WHICH particles enter; it does not reproduce HOW a
//    calorimeter measures them. Charged hadrons enter at their track momentum,
//    not the lower calorimeter response; there is no 0.087 tower granularity,
//    no 0.3 GeV tower threshold, and no akPu pileup subtraction."
//
// Clustering the towers themselves closes the first three of those. The fourth
// -- pileup subtraction -- is still not done here: these are raw anti-kT jets
// over the towers, the direct analogue of the PF pseudoJets, and the UE is
// handled downstream by the same random-cone and mixed-event machinery the PF
// scan uses. That is deliberate, so the two are comparable.
//
// TOWER KINEMATICS. rechitanalyzerpp/tower stores massless towers: e/cosh(eta)
// reproduces et to 4e-7, and emEt + hadEt equals et exactly (checked on
// HiForestAOD.root, 100 events, 208821 towers). So the FastJet input is
//   (et cos phi, et sin phi, et sinh eta, e)
// with no mass term and no separate energy assumption.

#include <string>

// ---- tower selection --------------------------------------------------------

// Tower ET threshold. CMS HI calo jets are built from towers above a small
// threshold; without one, the thousands of near-zero towers per event dominate
// the input multiplicity and cost time without moving any jet axis. 0.3 GeV is
// the value pseudoJets.h names as the calo-jet convention.
//
// Measured on a peripheral forest event (hiBin 178): 2088 towers per event
// on average, of which 76.4% are above 0.3 GeV -- so this is not a small cut,
// and a scan that changes it cannot be compared with one that did not. The
// value goes in the output file name for that reason.
double towerEtMin = 0.3;

// Clustering acceptance. Towers reach |eta| = 5.05 (HF included). Jets are only
// used to |eta| < 1.6, but the clustering has to see beyond that or jets near
// the edge lose the constituents that belong to them. 3.0 keeps the whole
// barrel+endcap plus a margin and drops HF, whose energy scale is a separate
// problem and which cannot contribute to a |eta| < 1.6 jet with R = 0.4 anyway.
double towerEtaMaxCluster = 3.0;

// Optional: use only the electromagnetic or only the hadronic part of each
// tower. Off by default -- both on is the physical tower. These are here for
// the same reason caloConstituentsIncludeCharged is in pseudoJets.h: to bracket
// the answer when the calo response is in question.
bool towerUseEmOnly  = false;
bool towerUseHadOnly = false;

// ---- forest-style (akPu) tower-level pileup subtraction --------------------
//
// doTowerPUSub runs, on each event's OWN towers, an emulation of the iterative
// pileup subtraction the forest's akPu4CaloJets use (MultipleAlgoIterator):
// per-ieta-ring mean + nSigma * RMS subtracted from every tower, towers that
// go <= 0 dropped, cluster, exclude towers near jets above puPtMin, recompute
// the pedestal from what is left, re-subtract, recluster. Implementation and
// the differences from the forest that remain are documented in
// headers/functions/towerPUSubtraction.h.
//
// Fills the h_fastJetPt_towerPUSub* histograms only. It does not touch the
// existing RC-subtracted spectra, so the two background methods can be compared
// in one scan. Same-event only: in mixed-event running the current event's
// tower arrays are overwritten by the pool loop, and the forest's algorithm is
// defined per real event anyway, so the fill is skipped when doEventMixing.
// The histograms are booked unconditionally and stay empty when this is off;
// read the provenance to tell "off" from "no jets".
bool   doTowerPUSub        = true;
double towerPUSub_nSigma   = 1.0;   // nSigmaPU
double towerPUSub_ptMin    = 8.0;   // puPtMin for R = 0.4
double towerPUSub_radius   = 0.5;   // radiusPU
double towerPUSub_jetPtMin = 1.0;   // jetPtMin of akPu4CaloJets (akCaloJets_cfi.py)

// et of the tower as it should enter the clustering
inline double towerClusterEt(double et, double emEt, double hadEt)
{
  if(towerUseEmOnly)  return emEt;
  if(towerUseHadOnly) return hadEt;
  return et;
}

// true if this tower is clustered at all
inline bool isClusteringTower(double et, double eta)
{
  if(fabs(eta) > towerEtaMaxCluster) return false;
  return et > towerEtMin;
}

// output-name tag, so a scan cannot be confused with one run at a different
// threshold or with a different tower component
inline std::string towerTag()
{
  std::string s = "_towers";
  char buf[64];
  snprintf(buf, sizeof(buf), "_etMin-%.2f", towerEtMin);
  s += buf;
  snprintf(buf, sizeof(buf), "_etaCluster-%.1f", towerEtaMaxCluster);
  s += buf;
  if(towerUseEmOnly)  s += "_emOnly";
  if(towerUseHadOnly) s += "_hadOnly";
  return s;
}
