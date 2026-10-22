#pragma once
// ptRel x jet pT x DeepCSV histograms (2026-10-22, user: "3d maps of muon pTrel, jet pt, and the
// deepcsv b jet tagger"). Booked by PYTHIA_scan.C (per flavor, nominal jet pT T0 only) from this
// header, so any later data or PYTHIA+HYDJET scan that adds them shares the binning.
// Uses: the b / c / light ptRel templates in DeepCSV slices (a tagger-enriched or tagger-depleted
// ptRel fit, the tagger as a second discriminating observable next to ptRel), and the tagger's
// efficiency for muon-tagged jets per flavor vs jet pT.
// Name: h_muptrel_recoJetPt_deepCSV_inclRecoMuonTag_triggerOn[_<flavor>],
// x = muon ptRel [GeV], y = reco jet pT [GeV], z = DeepCSV discriminant of the tagged jet; same
// selection and weight as h_muptrel_recoJetPt_inclRecoMuonTag_triggerOn.
// DeepCSV axis: one bin [-2, 0) for jets the tagger did not evaluate (the forest stores -1),
// then 0.02-wide bins over [0, 1] (the last edge 1.001 so a value of exactly 1 is inside).
// -99 (tree without the branch, eventMap::loadDeepCSV) lands in the underflow.
#include "TH3D.h"
#include <vector>

//
// CSVv2 (2026-10-22, user: "take the csvv2 instead of deepcsv", the tagger of the CMS b-jet shapes
// analyses AN-2017/337 and AN-2020/029, working point 0.9): same axes, but the forest stores -10
// for unevaluated jets, so the first z bin is [-11, 0). PYTHIA_scan.C and pp_scan.C book these
// (h_muptrel_recoJetPt_csvV2_inclRecoMuonTag_triggerOn[_<flavor>]); PYTHIAHYDJET_scan.C still
// books the DeepCSV ones.
inline TH3D* bookPtRelJetPtBTag3D(const char *name, const char *title, const char *taggerName, double unevaluatedLow)
{
  std::vector<double> x; for(int i = 0; i <= 50; i++) x.push_back(0.1*i);   // 0-5 GeV, the fit domain; overflow above
  const std::vector<double> y = {20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 135, 150, 175, 200, 250, 300, 500};
  std::vector<double> z = {unevaluatedLow};
  for(int i = 0; i < 50; i++) z.push_back(0.02*i);
  z.push_back(1.001);
  TH3D *h = new TH3D(name, Form("%s;muon #it{p}_{T}^{rel} [GeV];jet #it{p}_{T} [GeV];%s discriminant", title, taggerName),
                     (int)x.size() - 1, x.data(), (int)y.size() - 1, y.data(), (int)z.size() - 1, z.data());
  h->Sumw2();
  return h;
}

inline TH3D* bookPtRelJetPtDeepCSV3D(const char *name, const char *title){ return bookPtRelJetPtBTag3D(name, title, "DeepCSV", -2.); }
inline TH3D* bookPtRelJetPtCSVv2_3D(const char *name, const char *title){ return bookPtRelJetPtBTag3D(name, title, "CSVv2", -11.); }
