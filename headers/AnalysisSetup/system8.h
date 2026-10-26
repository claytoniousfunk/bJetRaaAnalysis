#pragma once
// System8 and negative-tag inputs for a DATA measurement of the CSVv2 b-tag efficiency of muon-tagged
// jets (2026-10-26, user: "configure for pp, pythia, PbPb, and PYTHIA+HYDJET scans"). Booked and filled by
// pp_scan.C, PYTHIA_scan.C, PbPb_scan.C and PYTHIAHYDJET_scan.C (doSystem8Maps in their configs), next to
// the ptRel x jet pT x CSVv2 map, for the same triggered muon-tagged jets with the same weight and jet pT.
//
// SYSTEM8 (CMS BTV): two samples of muon-tagged jets that differ in b fraction,
//   n   the muon jet has an AWAY jet: the leading other reco jet with |eta| < s8AwayEtaMax,
//       pT > s8AwayPtMin and |Delta phi| > s8AwayDPhiMin (forest jtpt; no jetTrkMax / calo-match / gen-match
//       requirement, so in PYTHIA+HYDJET a HYDJET jet can be the away jet, as in PbPb data)
//   p   ... and the away jet is b-tagged, CSVv2 > s8AwayTagCut
// and two cuts on the muon jet, the CSVv2 working point (z axis) and a ptRel cut (x axis, chosen offline).
// The 8 counts solve for the b and non-b content of n and p and the b / non-b efficiencies of both cuts;
// MC supplies the correlation factors (kappa, alpha, beta) from the same maps per flavor.
//   h_s8n_muptrel_recoJetPt_csvV2_inclRecoMuonTag_triggerOn[_<flavor>][_C<i>]
//   h_s8p_muptrel_recoJetPt_csvV2_inclRecoMuonTag_triggerOn[_<flavor>][_C<i>]
// NEGATIVE TAG: not booked. The forests' ndiscr_csvV2 (and pdiscr_csvV2, ndiscr_ssvHighEff) is 0 for every jet
// in PbPb and PYTHIA+HYDJET and absent in pp (checked 2026-10-26): the negative taggers were not run, so the
// light-mistag method of AN-2020/029 needs a re-forest. eventMap::loadNegCSVv2 is kept for when they are.
// Axes: x muon ptRel 0-5 GeV (0.1 GeV), y jet pT on the CSVv2-map edges (20-500 GeV), z the discriminant in
// WORKING-POINT categories [-11, 0) (not evaluated) [0, 0.54) [0.54, 0.84) [0.84, 0.90) [0.90, 0.96) [0.96, 1]:
// ~10x smaller than the full CSVv2 map (PYTHIA+HYDJET books 19 classes x 9 flavors of each).
#include "TH3D.h"
#include "TMath.h"
#include <vector>
#include <cmath>

double s8AwayEtaMax = 1.6, s8AwayPtMin = 30., s8AwayDPhiMin = 2.5, s8AwayTagCut = 0.84;

inline TH3D* bookS8Map(const char *name, const char *title, const char *tagger)
{
  std::vector<double> x; for(int i = 0; i <= 50; i++) x.push_back(0.1*i);
  const std::vector<double> y = {20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 135, 150, 175, 200, 250, 300, 500};
  const std::vector<double> z = {-11., 0., 0.54, 0.84, 0.90, 0.96, 1.001};
  TH3D *h = new TH3D(name, Form("%s;muon #it{p}_{T}^{rel} [GeV];jet #it{p}_{T} [GeV];%s category", title, tagger),
                     (int) x.size() - 1, x.data(), (int) y.size() - 1, y.data(), (int) z.size() - 1, z.data());
  h->Sumw2();
  return h;
}

// index of the away jet of jet i (leading qualifying other jet), or -1; EM = eventMap (template: no include order)
template <class EM> inline int s8AwayJet(const EM *em, int i)
{
  int best = -1; double bestPt = -1.;
  for(int j = 0; j < em->njet; j++){
    if(j == i) continue;
    if(std::fabs(em->jeteta[j]) > s8AwayEtaMax || em->jetpt[j] < s8AwayPtMin) continue;
    double dphi = std::fabs(em->jetphi[j] - em->jetphi[i]); if(dphi > TMath::Pi()) dphi = 2.*TMath::Pi() - dphi;
    if(dphi < s8AwayDPhiMin) continue;
    if(em->jetpt[j] > bestPt){ bestPt = em->jetpt[j]; best = j; }
  }
  return best;
}
// 0 = no away jet, 1 = away jet not tagged, 2 = away jet tagged (CSVv2 > s8AwayTagCut)
template <class EM> inline int s8AwayState(const EM *em, int i)
{
  const int a = s8AwayJet(em, i);
  if(a < 0) return 0;
  return em->disc_csvV2[a] > s8AwayTagCut ? 2 : 1;
}
