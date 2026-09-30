#pragma once
// Emulation of the forest's akPu4CaloJets pileup subtraction on the towers of
// rechitanalyzerpp/tower.
//
// SOURCE. CMSSW_10_3_4:
//   RecoJets/JetProducers/plugins/VirtualJetProducer.cc   (produce)
//   RecoJets/JetProducers/src/PileUpSubtractor.cc          (geometry map)
//   RecoHI/HiJetAlgos/src/MultipleAlgoIterator.cc          (the algorithm)
//   HeavyIonsAnalysis/JetAnalysis/python/akCaloJets_cfi.py (akPu4CaloJets)
// with HiCaloJetParameters: nSigmaPU 1, radiusPU 0.5, puPtMin 8 (R = 0.4),
// jetPtMin 1, inputEtMin 0.3, doPVCorrection false, minimumTowersFraction 0
// (fillDescriptions default), dropZeroTowers true (untracked default),
// sumRecHits false.
//
// ALGORITHM, per event, as VirtualJetProducer::produce runs it:
//   1. inputs: every tower with et >= inputEtMin, at ANY eta (HF included --
//      the pedestal is per ieta ring, so HF does not feed the central rings,
//      but a faithful input list keeps them).
//   2. calculatePedestal(inputs): per ring, mean = sum(et)/nt and
//      sigma = nSigma * sqrt(sum(et^2)/nt - mean^2), where nt is the number of
//      GEOMETRIC cells in the ring (72 / 36 / 18), not the number of towers
//      present -- empty cells count as zero.
//   3. subtractPedestal: et -> et - mean - sigma; towers at <= 0 are dropped,
//      the rest rescaled (direction kept).
//   4. cluster anti-kT R, keep jets above jetPtMin.
//   5. calculateOrphanInput: for each of those jets with pT >= puPtMin, every
//      geometric cell within radiusPU of the jet axis is excluded, and its
//      ring's nt is decremented. Orphans are the ORIGINAL (unsubtracted) input
//      towers not in an excluded cell.
//   6. calculatePedestal(orphans), with the reduced nt.
//   7. offsetCorrectJets: sigma *= nSigma AGAIN (rescaleRMS -- a no-op at
//      nSigma = 1, nSigma^2 otherwise; reproduced as written), subtract from the
//      original inputs, drop <= 0, recluster, keep jets above jetPtMin. These
//      are the output jets; their pT is the forest's rawpt.
//   jtpu (the per-jet pileup) is sum over the output jet's constituents of
//   original et - subtracted et.
//
// KNOWN DIFFERENCES FROM THE FOREST that this cannot remove:
//   * The tower tree is vertex-corrected (rechitanalyzerpp: hasVtx true,
//     saveBothVtx false), the forest clustered origin-based towers
//     (doPVCorrection false). et, eta -- and so the 0.3 GeV input cut -- differ
//     slightly per tower, by an amount that grows with |vz|.
//   * Cell positions for the exclusion step are nominal cell centres (below),
//     not CaloGeometry positions. Only the dR < radiusPU decision uses them.
//   * The forest clusters with ghost areas (ClusterSequenceArea); ghosts do
//     not move anti-kT hard jets, so this uses a plain ClusterSequence.
//   * The forest's geomtowers_ is built by walking every HCAL DetId, all
//     depths and subdetectors, resetting the count on each new ieta run. For
//     the standard DetId ordering (iphi fastest) that ends at the iphi count of
//     the ring's last (subdet, depth) -- 72 / 36 / 18, matching the tower
//     grid. That is an inference from the ordering, not checked against a
//     CaloGeometry dump; the validation against the forest's rawpt is the test.
//
// Usage:
//   std::vector<TowerPUInput> in;  (fill from the tower arrays)
//   std::vector<TowerPUJet> jets = clusterTowersWithPUSub(in, 0.4);

#include <vector>
#include <map>
#include <set>
#include <cmath>
#include "TMath.h"
#include "fastjet/ClusterSequence.hh"

struct TowerPUInput { double et, eta, phi; int ieta, iphi; };
struct TowerPUJet   { double pt, eta, phi, pu; int nConst; };

namespace towerPU {

// cells in an ieta ring of the CaloTower grid; checked on the MinBias forest
// (every cell of every ring appears in the tower tree)
inline int nCellsInRing(int ieta)
{
  const int a = std::abs(ieta);
  if(a == 0 || a > 41) return 0;
  if(a <= 20) return 72;
  if(a <= 39) return 36;
  return 18;
}

// iphi values the ring actually uses: all for 72-cell rings, odd for 36-cell
// rings, 3,7,11,... for the 18-cell HF rings (read off the tower tree)
inline bool isValidIphi(int ieta, int iphi)
{
  const int n = nCellsInRing(ieta);
  if(n == 72) return iphi >= 1 && iphi <= 72;
  if(n == 36) return iphi >= 1 && iphi <= 71 && (iphi % 2 == 1);
  if(n == 18) return iphi >= 3 && iphi <= 71 && (iphi % 4 == 3);
  return false;
}

// nominal cell centre in phi, wrapped to (-pi, pi]. Cell edges are multiples
// of 5 deg: a 72-cell iphi covers [iphi-1, iphi]*5 deg, a 36-cell odd iphi
// [iphi-1, iphi+1]*5, an 18-cell iphi [iphi-1, iphi+3]*5 -- centres confirmed
// against the HE/HF tower phi in the tree (barrel tower phi is energy-weighted
// and scatters around the centre).
inline double cellPhi(int ieta, int iphi)
{
  const int n = nCellsInRing(ieta);
  double deg = (n == 72) ? (iphi - 0.5) * 5.
             : (n == 36) ?  iphi * 5.
             :             (iphi + 1) * 5.;
  double p = deg * TMath::Pi() / 180.;
  if(p > TMath::Pi()) p -= 2. * TMath::Pi();
  return p;
}

// nominal ring centre in |eta|: mean tower eta per ring on the MinBias forest
// (vertex-smeared, so a few 1e-3 from the geometric centre -- irrelevant for a
// dR < 0.5 decision)
inline double cellEta(int ieta)
{
  static const double c[42] = { 0.,
    0.044, 0.132, 0.219, 0.307, 0.394, 0.483, 0.570, 0.656, 0.742, 0.832,
    0.919, 1.005, 1.091, 1.181, 1.266, 1.352, 1.439, 1.526, 1.609, 1.698,
    1.787, 1.882, 1.988, 2.111, 2.250, 2.412, 2.579, 2.758, 2.861, 3.044,
    3.219, 3.394, 3.569, 3.745, 3.919, 4.095, 4.270, 4.444, 4.620, 4.796,
    5.033 };
  const int a = std::abs(ieta);
  if(a == 0 || a > 41) return 0.;
  return ieta > 0 ? c[a] : -c[a];
}

struct Cell { int ieta, iphi; double eta, phi; };

inline const std::vector<Cell>& allCells()
{
  static std::vector<Cell> cells;
  if(cells.empty()){
    for(int ie = -41; ie <= 41; ie++){
      if(ie == 0) continue;
      for(int ip = 1; ip <= 72; ip++)
        if(isValidIphi(ie, ip)) cells.push_back({ie, ip, cellEta(ie), cellPhi(ie, ip)});
    }
  }
  return cells;
}

inline double dR(double eta1, double phi1, double eta2, double phi2)
{
  double dphi = std::fabs(phi1 - phi2);
  if(dphi > TMath::Pi()) dphi = 2. * TMath::Pi() - dphi;
  const double deta = eta1 - eta2;
  return std::sqrt(deta * deta + dphi * dphi);
}

// mean and sigma per ring from a set of towers, over (cells - excluded) cells
inline void pedestal(const std::vector<TowerPUInput> &in, const std::vector<int> &use,
                     const std::map<int,int> &nExcluded, double nSigma,
                     std::map<int,double> &mean, std::map<int,double> &sigma)
{
  std::map<int,double> s1, s2;
  for(int k : use){ s1[in[k].ieta] += in[k].et; s2[in[k].ieta] += in[k].et * in[k].et; }
  mean.clear(); sigma.clear();
  for(int ie = -41; ie <= 41; ie++){
    if(ie == 0) continue;
    auto ex = nExcluded.find(ie);
    const int nt = nCellsInRing(ie) - (ex == nExcluded.end() ? 0 : ex->second);
    if(nt > 0){
      const double m = s1[ie] / nt;
      double var = s2[ie] / nt - m * m;
      if(var < 0.) var = 0.;
      mean[ie] = m; sigma[ie] = nSigma * std::sqrt(var);
    }
    else { mean[ie] = 0.; sigma[ie] = 0.; }
  }
}

// subtracted, clipped inputs as PseudoJets; user_index = index into `in`
inline std::vector<fastjet::PseudoJet> subtracted(const std::vector<TowerPUInput> &in,
                                                  const std::map<int,double> &mean,
                                                  const std::map<int,double> &sigma)
{
  std::vector<fastjet::PseudoJet> out;
  out.reserve(in.size());
  for(size_t k = 0; k < in.size(); k++){
    const double etNew = in[k].et - mean.at(in[k].ieta) - sigma.at(in[k].ieta);
    if(etNew <= 0.) continue;                            // dropZeroTowers
    const double eta = in[k].eta, phi = in[k].phi;
    fastjet::PseudoJet p(etNew * std::cos(phi), etNew * std::sin(phi),
                         etNew * std::sinh(eta), etNew * std::cosh(eta));
    p.set_user_index((int)k);
    out.push_back(p);
  }
  return out;
}

} // namespace towerPU

// Towers must already have passed the input et cut (inputEtMin, 0.3 GeV). The
// returned jets are sorted by pT and include every eta -- apply the analysis
// acceptance at the call site.
inline std::vector<TowerPUJet> clusterTowersWithPUSub(const std::vector<TowerPUInput> &in,
                                                      double R,
                                                      double nSigma   = 1.0,
                                                      double puPtMin  = 8.0,
                                                      double radiusPU = 0.5,
                                                      double jetPtMin = 1.0)
{
  using namespace towerPU;
  std::vector<TowerPUJet> result;
  if(in.empty()) return result;

  const fastjet::JetDefinition def(fastjet::antikt_algorithm, R);
  std::vector<int> all(in.size());
  for(size_t k = 0; k < in.size(); k++) all[k] = (int)k;

  // steps 2-4: first pass
  std::map<int,int> nExcluded;
  std::map<int,double> mean, sigma;
  pedestal(in, all, nExcluded, nSigma, mean, sigma);
  std::vector<fastjet::PseudoJet> in1 = subtracted(in, mean, sigma);
  fastjet::ClusterSequence cs1(in1, def);
  std::vector<fastjet::PseudoJet> jets1 = fastjet::sorted_by_pt(cs1.inclusive_jets(jetPtMin));

  // step 5: exclude cells around the hard first-pass jets
  std::set<std::pair<int,int>> excluded;
  for(const auto &j : jets1){
    if(j.perp() < puPtMin) continue;
    const double jeta = j.eta(), jphi = j.phi_std();
    for(const Cell &c : allCells()){
      if(dR(c.eta, c.phi, jeta, jphi) >= radiusPU) continue;
      if(excluded.count({c.ieta, c.iphi})) continue;
      const int already = nExcluded.count(c.ieta) ? nExcluded[c.ieta] : 0;
      if(nCellsInRing(c.ieta) - already > 0){          // minimumTowersFraction = 0
        nExcluded[c.ieta] = already + 1;
        excluded.insert({c.ieta, c.iphi});
      }
    }
  }
  std::vector<int> orphans;
  orphans.reserve(in.size());
  for(size_t k = 0; k < in.size(); k++)
    if(!excluded.count({in[k].ieta, in[k].iphi})) orphans.push_back((int)k);

  // steps 6-7: pedestal from orphans, double sigma scaling, subtract, recluster
  pedestal(in, orphans, nExcluded, nSigma, mean, sigma);
  for(auto &s : sigma) s.second *= nSigma;             // rescaleRMS
  std::vector<fastjet::PseudoJet> in2 = subtracted(in, mean, sigma);
  fastjet::ClusterSequence cs2(in2, def);
  std::vector<fastjet::PseudoJet> jets2 = fastjet::sorted_by_pt(cs2.inclusive_jets(jetPtMin));

  result.reserve(jets2.size());
  for(const auto &j : jets2){
    double pu = 0.;
    const std::vector<fastjet::PseudoJet> cons = j.constituents();
    for(const auto &c : cons) pu += in[c.user_index()].et - c.perp();
    result.push_back({j.perp(), j.eta(), j.phi_std(), pu, (int)cons.size()});
  }
  return result;
}
