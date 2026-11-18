#pragma once
// Z+jet pT balance in data, for an in-situ check of the PbPb jet energy scale
// against pp (2026-11-18, user). The Z is not quenched and its pT is measured by
// the muons, so in peripheral PbPb <x_jZ> = pT^jet / pT^Z should equal pp if the
// two jet energy scales agree.
//
// Per event: the highest-pT opposite-sign tight-muon pair with 70 < m < 110 GeV
// (both muons from the scan's Z loop, |eta| < 2.4, its trigger and pair weight),
// then the LEADING jet back to back with it (|dphi| > 7 pi / 8) among the jets
// passing the scan's own jet selection (JEC, calo match if required, trkMax
// filter, eta-phi mask, |eta| < etaMax, pT > jetPtCut), excluding jets within
// dR < 0.4 of either Z muon.
//   h_zJetXjz_zPt   x_jZ (0-2.5, 0.02 bins) x Z pT, leading back-to-back jet
//   h_zPt_zJet      every Z candidate, Z pT: the per-Z normalization
// Jets below jetPtCut are lost, so x_jZ is truncated at jetPtCut / pT^Z: compare
// pp and PbPb at the same Z pT and read the region x_jZ > jetPtCut / pT^Z.

#include "TH1D.h"
#include "TH2D.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TVector2.h"

namespace zjb {

const double massLo = 70., massHi = 110.;
const double dPhiMin = 7.*TMath::Pi()/8.;
const double muJetDrMin = 0.4;
const double muMass = 0.105658;
const int nZPt = 9;
const double zPtEdges[nZPt + 1] = {20., 30., 40., 50., 60., 80., 100., 150., 250., 500.};
const int nX = 125;
const double xMin = 0., xMax = 2.5;

inline TH2D* bookXjz(const char *name, const char *title){
  TH2D *h = new TH2D(name, Form("%s;x_{jZ} = p_{T}^{jet} / p_{T}^{Z};p_{T}^{Z} (GeV)", title), nX, xMin, xMax, nZPt, zPtEdges);
  h->Sumw2(); return h;
}
inline TH1D* bookZPt(const char *name, const char *title){
  TH1D *h = new TH1D(name, Form("%s;p_{T}^{Z} (GeV)", title), nZPt, zPtEdges);
  h->Sumw2(); return h;
}

struct ZCand {
  bool ok = false;
  double pt = -1., phi = 0., w = 0.;
  double eta1 = 0., phi1 = 0., eta2 = 0., phi2 = 0.;
};

// call for every opposite-sign pair; keeps the highest-pT candidate in the window
inline void consider(ZCand &z, double pt1, double eta1, double phi1, double pt2, double eta2, double phi2, double w){
  TLorentzVector a, b; a.SetPtEtaPhiM(pt1, eta1, phi1, muMass); b.SetPtEtaPhiM(pt2, eta2, phi2, muMass);
  TLorentzVector s = a + b;
  if(s.M() < massLo || s.M() > massHi) return;
  if(s.Pt() <= z.pt) return;
  z.ok = true; z.pt = s.Pt(); z.phi = s.Phi(); z.w = w;
  z.eta1 = eta1; z.phi1 = phi1; z.eta2 = eta2; z.phi2 = phi2;
}

inline bool nearZMuon(const ZCand &z, double eta, double phi){
  auto dr = [](double e1, double p1, double e2, double p2){ double dp = TVector2::Phi_mpi_pi(p1 - p2); return sqrt((e1 - e2)*(e1 - e2) + dp*dp); };
  return dr(eta, phi, z.eta1, z.phi1) < muJetDrMin || dr(eta, phi, z.eta2, z.phi2) < muJetDrMin;
}

inline bool backToBack(const ZCand &z, double phi){ return fabs(TVector2::Phi_mpi_pi(phi - z.phi)) > dPhiMin; }

} // namespace zjb
