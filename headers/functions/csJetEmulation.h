#pragma once
// Emulation of the forest's akCs4PFJets (constituent-subtracted anti-kT R = 0.4
// PF jets) on the PF candidates of pfcandAnalyzer/pfTree -- the PF-jet
// counterpart of towerPUSubtraction.h.
//
// SOURCE. CmsHI/cmssw forest_CMSSW_10_3_1 (the forest is CMSSW_10_3_4-571-...,
// HiForestInfo), which RE-RUNS these producers rather than reading the AOD jets
// (fullJetSequence_pponAA_data_cff -> rerecoRho_cff, rerecoJets_cff):
//   HeavyIonsAnalysis/JetAnalysis/python/akPFJets_cfi.py     (akCs4PFJets)
//   HeavyIonsAnalysis/JetAnalysis/python/ktPFJets_cfi.py     (kt4PFJetsForRho)
//   RecoHI/HiJetAlgos/python/HiRecoPFJets_cff.py              (kt4PFJetsForRho)
//   RecoHI/HiJetAlgos/python/HiPFJetParameters_cff.py         (ghost spec)
//   RecoHI/HiJetAlgos/python/hiFJRhoProducer.py               (rho config)
//   RecoHI/HiJetAlgos/plugins/HiFJRhoProducer.cc              (rho, rhom)
//   RecoJets/JetProducers/plugins/CSJetProducer.cc            (the subtraction)
//   RecoJets/JetProducers/plugins/VirtualJetProducer.cc       (inputs, seed)
// The forest's CSJetProducer.cc and VirtualJetProducer.cc differ from the
// plain CMSSW_10_3_4 release; these are the forest versions.
// Externals: fastjet 3.3.0 (cms-externals fork: thread_local RNG and atomics
// only, ghost placement and RNG identical to vanilla, built -O3 -ffast-math)
// and fjcontrib 1.033 (ConstituentSubtractor byte-identical to hepforge).
// Build against those (~/Programs/fastjet330 here, or the CMSSW_10_3_4
// externals on lxplus): exact reproduction relies on the same ghost positions.
//
// ALGORITHM, per event:
//   0. Inputs: every particleFlow candidate as (px, py, pz, E) WITH its mass
//      (inputEtMin = inputEMin = 0; candidates with pT < 100 eps dropped).
//      Muons included. The tree keeps |eta| < 5 only (pfAbsEtaMax), which can
//      move only kt jets near |eta| = 5 -- not the rho of any bin a
//      |eta| < 1.6 jet touches.
//   1. Seed the FastJet RNG from run/event, exactly as VirtualJetProducer
//      does before every clustering (setForestSeed).
//   2. kt4PFJetsForRho: kt R = 0.4, ClusterSequenceArea with active area,
//      GhostedAreaSpec(6.5, 1, 0.005), jets above 0, sorted by pT.
//   3. HiFJRhoProducer: skip the first 2 kt jets with |eta| < 2 and pT > 20,
//      and the first 1 with 2 < |eta| < 3 and pT > 20; in each bin of
//      {-5, -3, -2.1, -1.3, 1.3, 2.1, 3, 5}, rho = median pT/area and
//      rhom = median sum(mT - pT)/area over the kt jets with eta in
//      [lo + 0.2, hi - 0.2). pT, area and eta pass through float, as in the
//      producer; an empty bin keeps 1e-6.
//   4. Reseed, then anti-kT R = 0.4 with EXPLICIT ghosts, same ghost spec,
//      jets above 0.
//   5. CSJetProducer, per jet with at least one real particle: every ghost
//      gets pT = rho(eta) * area and mass from rhom (no flow modulation,
//      useModulatedRho = false for akCs4PF); fjcontrib ConstituentSubtractor
//      with deltaR distance, max_distance -1 (none), alpha 2, mass
//      subtraction, zero-pT particles removed. The jet is join() of the
//      subtracted particles, kept if pT > 0. Its pT is the forest's rawpt, its
//      direction the forest's jteta / jtphi -- the SUBTRACTED axis.
//
// VALIDATION (src/plots/csJets/validateCSJets.C, MinBias Part 1 forest, 6946
// events, 2026-10-29). Forest akCs4PF jets, rawpt > 10 GeV, |eta| < 1.6:
// 99.4-100% matched, and emulated / forest rawpt within 1e-5 for 95.5 / 97.1
// / 98.0 / 98.1% (0-10 / 10-30 / 30-50 / 50-80%), within 1e-2 for 97-100%;
// emulated jets above 30 GeV raw have a forest partner 1835 of 1836. Matched
// jets are exact to float precision; when one is not, it typically differs by
// a single constituent.
// Two limits, both from the INPUT, not the algorithm:
//   * The forest's pfcandAnalyzer keeps |eta| < 3 only (no HF candidates),
//     while the forest CS ran on HF. Removing more forward input shows how
//     far that reaches: dropping 2.5 < |eta| < 3 as well takes the |eta| < 1
//     jets in 0-30% from 98.5% to 84% exact. The non-exact fraction also grows
//     with |eta| (0-30%: 98% / 89% / 85% within 1e-5 at |eta| < 1.0 / 1.0-1.3
//     / 1.3-1.6, mean emulated/forest - 1 at most 1.2e-3). Most of the
//     remaining mismatch is therefore probably the missing HF (likely, not
//     proven). A forest with pfAbsEtaMax 5 would settle it. CS particles
//     beyond |eta| ~ 1.7 (the [2.1, 3] rho bin) do not reproduce from this
//     tree.
//   * Float storage of the candidates (and the forest's -ffast-math build)
//     can flip near-ties in the CS ghost-particle ordering: one ghost moves to
//     another particle, ~0.05 GeV.
//
// Not emulated: akFlowPuCs4PFJets (hiPuRhoR3Producer on PF towers plus the
// flow fit) -- a different collection.
//
// COST. ~16 300 ghosts per clustering (|y| < 6.5, area 0.005), two clusterings
// per event, plus the subtraction in every jet. ptUnsubMin skips the
// subtraction of jets whose UNSUBTRACTED pT is below it; CS only removes pT,
// so every jet above that threshold after subtraction is still exact.
//
// Usage:
//   std::vector<CSJetInput> in;  (fill from the pfcandAnalyzer vectors)
//   std::vector<CSJet> jets = csJet::forestCSJets(in, run, evt, 10.);

#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/Selector.hh"
#include "fastjet/contrib/ConstituentSubtractor.hh"

struct CSJetInput       { double pt, eta, phi, m; int id; };
struct CSJetConstituent { double pt, eta, phi, m; int index; };   // index into the input vector
struct CSJet {
  double pt, eta, phi, m;     // subtracted jet: the forest's rawpt, jteta, jtphi
  double ptUnsub;             // the same constituents before subtraction
  std::vector<CSJetConstituent> constituents;
};
struct CSRhoMap { std::vector<double> etaEdges, rho, rhom; };

namespace csJet {

const double R            = 0.4;
const double ghostRapMax  = 6.5;     // HiPFJetParameters Ghost_EtaMax
const double ghostArea    = 0.005;
const int    ghostRepeats = 1;
const unsigned int minSeed = 14327;

const double csAlpha       = 2.;
const double csMaxDistance = -1.;

const int    nExcl       = 2;
const double etaMaxExcl  = 2.;
const double ptMinExcl   = 20.;
const int    nExcl2      = 1;
const double etaMaxExcl2 = 3.;
const double ptMinExcl2  = 20.;
const double rhoEdgeRadius = 0.2;    // kt jet must sit this far inside its bin
inline const std::vector<double>& etaRanges()
{
  static const std::vector<double> e = {-5., -3., -2.1, -1.3, 1.3, 2.1, 3., 5.};
  return e;
}

// VirtualJetProducer::produce with useDeterministicSeed. The arithmetic is
// unsigned, then stored in an int vector, as there.
inline void setForestSeed(unsigned int run, unsigned int evt)
{
  fastjet::GhostedAreaSpec gas;
  std::vector<int> seeds(2);
  seeds[0] = std::max(run, minSeed + 3) + 3 * evt;
  seeds[1] = std::max(run, minSeed + 5) + 5 * evt;
  gas.set_random_status(seeds);
}

// any other seed, e.g. per mixed event (keep both positive)
inline void setSeed(int s0, int s1)
{
  fastjet::GhostedAreaSpec gas;
  std::vector<int> seeds = {s0, s1};
  gas.set_random_status(seeds);
}

inline fastjet::AreaDefinition areaDefinition(bool explicitGhosts)
{
  fastjet::GhostedAreaSpec spec(ghostRapMax, ghostRepeats, ghostArea);
  return fastjet::AreaDefinition(explicitGhosts ? fastjet::active_area_explicit_ghosts
                                                : fastjet::active_area, spec);
}

// massive 4-vectors, user_index = position in `in`
inline std::vector<fastjet::PseudoJet> makeInputs(const std::vector<CSJetInput> &in)
{
  std::vector<fastjet::PseudoJet> out;
  out.reserve(in.size());
  for(size_t k = 0; k < in.size(); k++){
    const CSJetInput &c = in[k];
    if(!std::isfinite(c.pt)) continue;
    if(c.pt < 100 * std::numeric_limits<double>::epsilon()) continue;
    const double px = c.pt * std::cos(c.phi), py = c.pt * std::sin(c.phi);
    const double pz = c.pt * std::sinh(c.eta);
    const double E  = std::sqrt(px * px + py * py + pz * pz + c.m * c.m);
    fastjet::PseudoJet p(px, py, pz, E);
    p.set_user_index((int)k);
    out.push_back(p);
  }
  return out;
}

// HiFJRhoProducer::calcMedian
inline double forestMedian(std::vector<double> v)
{
  if(v.empty()) return 0.;
  const size_t n = v.size() / 2;
  std::nth_element(v.begin(), v.begin() + n, v.end());
  double med = v[n];
  if(!(v.size() & 1)){
    const double mx = *std::max_element(v.begin(), v.begin() + n);
    med = (mx + med) / 2.;
  }
  return med;
}

// steps 2-3. Does not seed: call setForestSeed (or seed otherwise) first.
inline CSRhoMap computeRho(const std::vector<fastjet::PseudoJet> &inputs)
{
  fastjet::JetDefinition jetDef(fastjet::kt_algorithm, R);
  fastjet::ClusterSequenceArea cs(inputs, jetDef, areaDefinition(false));
  std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(cs.inclusive_jets(0.));

  const std::vector<double> &edges = etaRanges();
  const int neta = (int)edges.size();
  CSRhoMap map;
  map.etaEdges = edges;
  map.rho.assign(neta - 1, 1e-6);
  map.rhom.assign(neta - 1, 1e-6);

  std::vector<double> rhoVec, rhomVec, etaVec;
  int njetsEx = 0, njetsEx2 = 0;
  for(const fastjet::PseudoJet &jet : jets){
    const double jpt = jet.pt(), jeta = jet.eta();
    if(njetsEx < nExcl && std::fabs(jeta) < etaMaxExcl && jpt > ptMinExcl){ njetsEx++; continue; }
    if(njetsEx2 < nExcl2 && std::fabs(jeta) < etaMaxExcl2 && std::fabs(jeta) > etaMaxExcl
       && jpt > ptMinExcl2){ njetsEx2++; continue; }
    const float pt   = jpt;
    const float area = jet.area();
    const float eta  = jeta;
    if(eta < edges[0] || eta > edges[neta - 1]) continue;
    if(area > 0.){
      double md = 0.;                                    // calcMd
      for(const fastjet::PseudoJet &c : jet.constituents())
        md += std::sqrt(c.m() * c.m() + c.pt() * c.pt()) - c.pt();
      rhoVec.push_back(pt / area);
      rhomVec.push_back(md / area);
      etaVec.push_back(eta);
    }
  }
  for(int ie = 0; ie < neta - 1; ie++){
    std::vector<double> r, rm;
    const double lo = edges[ie] + rhoEdgeRadius, hi = edges[ie + 1] - rhoEdgeRadius;
    for(size_t i = 0; i < etaVec.size(); i++)
      if(etaVec[i] >= lo && etaVec[i] < hi){ r.push_back(rhoVec[i]); rm.push_back(rhomVec[i]); }
    if(!r.empty()){ map.rho[ie] = forestMedian(r); map.rhom[ie] = forestMedian(rm); }
  }
  return map;
}

// CSJetProducer's ghost-to-bin assignment
inline int rhoBin(const CSRhoMap &map, double eta)
{
  const std::vector<double> &e = map.etaEdges;
  if(eta <= e.front()) return 0;
  if(eta >= e.back())  return (int)map.rho.size() - 1;
  for(size_t ie = 0; ie + 1 < e.size(); ie++)
    if(eta >= e[ie] && eta < e[ie + 1]) return (int)ie;
  return -1;
}

// steps 4-5. Does not seed: call setForestSeed (or seed otherwise) first.
inline std::vector<CSJet> clusterCS(const std::vector<fastjet::PseudoJet> &inputs,
                                    const CSRhoMap &map, double ptUnsubMin = 0.)
{
  fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, R);
  fastjet::ClusterSequenceArea cs(inputs, jetDef, areaDefinition(true));
  std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(cs.inclusive_jets(0.));

  fastjet::contrib::ConstituentSubtractor subtractor;
  subtractor.set_distance_type(fastjet::contrib::ConstituentSubtractor::deltaR);
  subtractor.set_max_distance(csMaxDistance);
  subtractor.set_alpha(csAlpha);
  subtractor.set_do_mass_subtraction(true);
  subtractor.set_remove_all_zero_pt_particles(true);

  std::vector<CSJet> out;
  for(const fastjet::PseudoJet &jet : jets){
    if(jet.pt() < ptUnsubMin) break;                     // sorted: the rest are softer
    std::vector<fastjet::PseudoJet> particles, ghosts;
    fastjet::SelectorIsPureGhost().sift(jet.constituents(), ghosts, particles);
    if(particles.empty()) continue;
    for(fastjet::PseudoJet &g : ghosts){
      const int ib = rhoBin(map, g.eta());
      const double pt = map.rho[ib] * g.area();
      const double m2 = std::pow(map.rhom[ib] * g.area() + pt, 2) - std::pow(pt, 2);
      g.reset_momentum_PtYPhiM(pt, g.rap(), g.phi(), m2 > 0 ? std::sqrt(m2) : 0.);
    }
    std::vector<fastjet::PseudoJet> sub = subtractor.do_subtraction(particles, ghosts);
    fastjet::PseudoJet sj = fastjet::join(sub);
    if(!(sj.perp() > 0.)) continue;

    CSJet j;
    j.pt = sj.pt(); j.eta = sj.eta(); j.phi = sj.phi_std(); j.m = sj.m();
    fastjet::PseudoJet unsub = fastjet::join(particles);
    j.ptUnsub = unsub.pt();
    for(const fastjet::PseudoJet &c : sub)
      j.constituents.push_back({c.pt(), c.eta(), c.phi_std(), c.m(), c.user_index()});
    out.push_back(j);
  }
  std::sort(out.begin(), out.end(), [](const CSJet &a, const CSJet &b){ return a.pt > b.pt; });
  return out;
}

// the whole chain for one event, seeded as the forest seeds that event
inline std::vector<CSJet> forestCSJets(const std::vector<CSJetInput> &in,
                                       unsigned int run, unsigned int evt,
                                       double ptUnsubMin = 0., CSRhoMap *rhoOut = nullptr)
{
  std::vector<fastjet::PseudoJet> inputs = makeInputs(in);
  setForestSeed(run, evt);
  CSRhoMap map = computeRho(inputs);
  setForestSeed(run, evt);
  std::vector<CSJet> jets = clusterCS(inputs, map, ptUnsubMin);
  if(rhoOut) *rhoOut = map;
  return jets;
}

} // namespace csJet
