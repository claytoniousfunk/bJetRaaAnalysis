#pragma once
// ptRel x muon pT x jet pT histograms (2026-10-09, user: "ptRel vs muon pT per
// flavor"). Booked by pp_scan.C, PbPb_scan.C (per centrality), PYTHIA_scan.C and
// PYTHIAHYDJET_scan.C (per flavor, nominal jet pT T0 only), all from this header
// so data and templates share the binning. Uses:
//   - the pp fit-vs-PYTHIA tension (studyPPFitVsPythia_PF.C): ptRel and muon pT
//     disagree on the b fraction; reweighting the b muon pT (fragmentation,
//     cascade decays) needs ptRel per muon-pT bin
//   - fits with a harder muon cut (any edge below: 15, 20, 25, ...)
//   - a joint ptRel x muon-pT fit, and an axis-smearing toy that draws ptRel and
//     muon pT jointly instead of independently
// Name: h_muptrel_mupt_recoJetPt_inclRecoMuonTag_triggerOn[_<flavor>][_C<i>],
// x = ptRel [GeV], y = muon pT [GeV], z = reco jet pT [GeV]; same selection and
// weight as h_muptrel_recoJetPt_inclRecoMuonTag_triggerOn in each scan.
// muPtLow (2026-10-29): the muon scans pass muPtCut. Below 15 GeV the y axis
// gains 2 GeV edges from muPtLow up (mu5: 7, 9, 11, 13; mu7: 9, 11, 13), which
// the mu5 / mu7 windows would otherwise put entirely in underflow. At the
// default 15 (mu12, and the track scans, whose track-pT axis is divided bin by
// bin by the muon one) the binning is unchanged.
#include "TH3D.h"
#include <vector>

inline TH3D* bookPtRelMuPt3D(const char *name, const char *title, double muPtLow = 15.)
{
  std::vector<double> x; for(int i = 0; i <= 50; i++) x.push_back(0.1*i);   // 0-5 GeV, the fit domain; overflow above
  std::vector<double> y;
  for(double e : {7., 9., 11., 13.}) if(e >= muPtLow - 1e-6) y.push_back(e);
  for(double e : {15., 17.5, 20., 22.5, 25., 27.5, 30., 32.5, 35., 40., 45., 50., 60., 70., 80., 100., 150.}) y.push_back(e);
  const std::vector<double> z = {20, 30, 40, 50, 60, 70, 80, 100, 120, 150, 200, 300, 500};
  TH3D *h = new TH3D(name, Form("%s;muon #it{p}_{T}^{rel} [GeV];muon #it{p}_{T} [GeV];jet #it{p}_{T} [GeV]", title),
                     (int)x.size() - 1, x.data(), (int)y.size() - 1, y.data(), (int)z.size() - 1, z.data());
  h->Sumw2();
  return h;
}
