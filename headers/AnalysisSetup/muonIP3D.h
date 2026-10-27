#pragma once
// Muon 3D impact parameter maps for an IP3D template fit of the b purity of muon-tagged jets (2026-10-27,
// user: "i want to try template fitting with muon-jet IP3D for the high pT region"). CMS BTV-12-001 (JINST 8
// P04013) Sec. 5.1: the muon IP separates b from c / light better than ptRel above a jet pT of about 120 GeV;
// they fit ln(|IP| [cm]) with b, c and light templates, separately in the CSVv2-tagged and untagged jets.
// Booked and filled by pp_scan.C, PYTHIA_scan.C, PbPb_scan.C and PYTHIAHYDJET_scan.C (doMuonIP3DMaps in their
// configs), next to the ptRel x jet pT x CSVv2 map, for the same triggered muon-tagged jets, muon, weight and
// jet pT:
//   h_muLnIP3D_recoJetPt_csvV2_inclRecoMuonTag_triggerOn[_<flavor>][_C<i>]   ln|IP3D| x jet pT x CSVv2 category
//   h_muptrel_muLnIP3D_recoJetPt_inclRecoMuonTag_triggerOn[_<flavor>][_C<i>] ptRel x ln|IP3D| x jet pT (joint fit)
// IP3D: the forest's muIP3D (ggHiNtuplizer, cm, unsigned, w.r.t. the primary vertex); eventMap::loadMuonIP3D.
// Axes: ln|IP3D| -11 to 1 in 0.25 (|IP| = 17 um to 2.7 cm; pp data muons > 15 GeV: mean -5.5, RMS 1.5, i.e.
// ~40 um); jet pT on the System8 edges; CSVv2 in the System8 working-point categories (system8.h), so tagged /
// untagged at any of 0.54 / 0.84 / 0.90 / 0.96 and, summed, all jets; ptRel 0-5 GeV in 0.5 GeV (the fit binning).
#include "TH3D.h"
#include "TMath.h"
#include <vector>
#include <cmath>

const double muLnIPMin = -11., muLnIPMax = 1.;
const int muLnIPBins = 48;

inline std::vector<double> muIP3DJetPtEdges(){ return {20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 135, 150, 175, 200, 250, 300, 500}; }

inline TH3D* bookMuIP3DMap(const char *name, const char *title)
{
  std::vector<double> x; for(int i = 0; i <= muLnIPBins; i++) x.push_back(muLnIPMin + i*(muLnIPMax - muLnIPMin)/muLnIPBins);
  const std::vector<double> y = muIP3DJetPtEdges();
  const std::vector<double> z = {-11., 0., 0.54, 0.84, 0.90, 0.96, 1.001};
  TH3D *h = new TH3D(name, Form("%s;muon ln(|IP_{3D}| [cm]);jet #it{p}_{T} [GeV];CSVv2 category", title),
                     (int) x.size() - 1, x.data(), (int) y.size() - 1, y.data(), (int) z.size() - 1, z.data());
  h->Sumw2();
  return h;
}
inline TH3D* bookMuPtRelIP3DMap(const char *name, const char *title)
{
  std::vector<double> x; for(int i = 0; i <= 10; i++) x.push_back(0.5*i);
  std::vector<double> y; for(int i = 0; i <= muLnIPBins; i++) y.push_back(muLnIPMin + i*(muLnIPMax - muLnIPMin)/muLnIPBins);
  const std::vector<double> z = muIP3DJetPtEdges();
  TH3D *h = new TH3D(name, Form("%s;muon #it{p}_{T}^{rel} [GeV];muon ln(|IP_{3D}| [cm]);jet #it{p}_{T} [GeV]", title),
                     (int) x.size() - 1, x.data(), (int) y.size() - 1, y.data(), (int) z.size() - 1, z.data());
  h->Sumw2();
  return h;
}
// ln|IP3D| of muon m, kept inside the axis (under / overflow folded into the edge bins); -99 if unavailable
template <class EM> inline double muLnIP3D(const EM *em, int m)
{
  if(!em->hasMuonIP3D || m < 0 || m >= (int) em->muIP3D->size()) return -99.;
  const double ip = std::fabs(em->muIP3D->at(m));
  double l = ip > 0. ? std::log(ip) : muLnIPMin;
  const double half = 0.5*(muLnIPMax - muLnIPMin)/muLnIPBins;
  if(l < muLnIPMin + half) l = muLnIPMin + half;
  if(l > muLnIPMax - half) l = muLnIPMax - half;
  return l;
}
