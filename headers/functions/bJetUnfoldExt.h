#pragma once
// b-jet unfolding with the truth axis extended below the analysis range
// (2026-10-08). Used by the PF chain: unfoldBJetSpectrum_caloJets_pp.C,
// bJetChain_caloJets_PbPb.C and bJetsPerZSys_PF_PbPb.C, so the three agree.
//
// Why: the analysis axis starts at 50 GeV, and the older rebinning dropped
// every response entry outside it -- including all b jets with gen pT < 50.
// In 0-10% PbPb 63% of the muon-tagged b jets reconstructed at 50-70 GeV come
// from gen < 50 (pp: 13%). With those rows missing, the unfolding had to put
// that data into truth 50-80 GeV, which in central events feed half of reco
// 80-100 GeV, so it explained away the 80-150 GeV data as migration: central
// b-jet R_Z came out 10-40% low at 80-150 GeV, with a dip at 100-150.
//
// Here both axes are the analysis edges with bExtLow (20, 30, 40) prepended:
// gen < 50 has its own truth bins (never reported), and reco below the floor
// is kept in the truth (as misses) but removed from the measured side. The
// result is returned on the analysis axis.

#include <vector>
#include <cmath>
#include "TH1D.h"
#include "TH2D.h"
#include "purityVariant.h"   // bLowFloorOn
// RooUnfold: unfoldOnExt is a template over the response and Bayes classes, so
// the names are resolved where a macro calls it -- unfoldOnExt<RooUnfoldResponse,
// RooUnfoldBayes>(...). Naming them here makes ROOT's interpreter autoparse
// RooUnfold at the include point, which fails with typedef redefinitions.

const int    bNExtLow = 3;
const double bExtLow[bNExtLow] = {20., 30., 40.};

// UNDER / OVERFLOW AXIS (user, 2026-10-09): results are reported at 80-200 GeV
// only, and the unfolding axis keeps just one bin on either side of that range:
//   [20, 70) [70, 80) | 80 100 120 150 200 | [200, 500]
// The reco side starts at the 70 GeV floor, so reco has ONE underflow bin
// [70, 80); the truth below 80 is split at the floor ([20, 70) is fed only by
// migration up into reco >= 70), so the axis stays square. Above 200 GeV reco
// and truth are one overflow bin. Analysis-axis bins inside a merged bin
// (50-60 ... 70-80, 200-300, 300-500) come back as the merged value spread
// uniformly in pT: finite placeholders, never reported. false: the fine
// extended axis above.
bool bUnfoldUfOf = true;
// BJET_UNFOLD_LOWFLOOR=1 (purityVariant.h, 2026-10-26): [20, 50) [50, 60) | 60 70 80 ... 200 | [200, 500]
// with the reco floor at 50 GeV (bRecoFloorPF), reported from 60 GeV. Nominal: 70 / 80 / floor 70.
// 2026-11-05 (user: 70 is too close to the first reported bin): nominal reco floor 60 GeV, axis
// [20, 60) [60, 80) | 80 100 120 150 200 | [200, 500]; 70 before.
const double bUfSplit = bLowFloorOn() ? 50. : 60., bUfHi = bLowFloorOn() ? 60. : 80., bOfLo = 200.;
inline double bRecoFloorPF(){ return bLowFloorOn() ? 50. : 60.; }

// Bayes iterations of the b-jet unfolding, shared by unfoldBJetSpectrum_caloJets_pp.C,
// bJetChain_caloJets_PbPb.C and bJetsPerZSys_PF_PbPb.C (PFJets). With the
// under/overflow axis: the minimum-MSE N of src/unfoldTest/optimalIterations_bJets_PF.C
// (2026-10-09: pp 3, same-sample data-like closure; PbPb 1 in every class, even/odd
// halves -- with the data-driven prior the bias hardly moves with N while the
// variance grows linearly). Before the under/overflow axis: pp 3, PbPb 4/4/3/3.
// 2026-10-21 (user): Bayes is the nominal b-jet unfolding again, with the MC PRIOR and N chosen as
// for the inclusive jets -- minimum MSE of the tilt-0.3 closure, same sample
// (src/unfoldTest/optimalIterations_bJets_PF_MCprior.C): pp 3, PbPb 2 / 2 / 2 / 1. Before (data-driven
// prior, even/odd data-like closure): pp 3, PbPb 1 / 1 / 1 / 1.
// Closure tilt of the b-jet unfolding non-closure and N_opt (2026-11-05, user): gen rows x (g_c/100)^(+bClosureTilt),
// g_c clamped to [40, 500] (closureTiltRows, dataDrivenPrior.h), 0.3 = harder, as for the inclusive jets;
// < -90 = the data-like fitted tilt (fitPriorTilt + tiltGenRows), used before.
double bClosureTilt = 0.3;
int bNIterPP = bNIterVariant() == 1 ? 6 : bNIterVariant() >= 3 ? 2 : 3;   // BJET_NITER=bias (purityVariant.h, 2026-11-11): 5% plateau on bias^2, pp 6, PbPb 10 / 8 / 6 / 6
int bNIterPbPb[4] = {bNIterVariant() == 1 ? 10 : bNIterVariant() == 4 ? 6 : bNIterVariant() == 3 ? 3 : 2, bNIterVariant() == 1 ? 8 : bNIterVariant() >= 3 ? 3 : 2, bNIterVariant() == 1 ? 6 : bNIterVariant() >= 3 ? 2 : 1, bNIterVariant() == 1 ? 6 : 1};   // BJET_NITER=opt: 3 / 3 / 2 / 1 (purityVariant.h)
//   // 0-10, 10-30, 30-50, 50-80%; 2026-11-05, floor 60, clamped harder tilt (3 / 2 / 2 / 1 unclamped softer; 2 / 2 / 2 / 1 at 70); unchanged on the JER-smeared (x1.10) response, 2026-11-05 (30-50%: N = 2 lowers the MSE by 4.7%, just under the 5% plateau rule)
// Data-driven prior (dataDrivenPrior.h) in the PF b-jet unfolding: false since 2026-10-21 (user, MC
// prior as for the inclusive jets). true restores the 2026-10-07 tilted-response prior and its
// "prior" systematic in bJetsPerZSys_PF_PbPb.C.
bool bUseDataDrivenPrior = false;

inline std::vector<double> bExtEdges(int nA, const double *ana)   // nA = number of edges
{
  std::vector<double> e;
  if(bUnfoldUfOf){
    e.push_back(bExtLow[0]); e.push_back(bUfSplit);
    for(int i = 0; i < nA; i++) if(ana[i] > bUfHi - 1e-6 && ana[i] < bOfLo + 1e-6) e.push_back(ana[i]);
    e.push_back(ana[nA - 1]);
    return e;
  }
  for(int i = 0; i < bNExtLow; i++) if(bExtLow[i] < ana[0] - 1e-6) e.push_back(bExtLow[i]);
  for(int i = 0; i < nA; i++) e.push_back(ana[i]);
  return e;
}

// fine (reco x, gen y) -> e x e, by bin centre; entries outside e are dropped
inline TH2D* rebin2DTo(const TH2D *h, const std::vector<double> &e, const char *name)
{
  const int n = (int) e.size() - 1;
  TH2D *r = new TH2D(name, "", n, &e[0], n, &e[0]); r->SetDirectory(nullptr); r->Sumw2();
  for(int ix = 1; ix <= h->GetNbinsX(); ix++) for(int iy = 1; iy <= h->GetNbinsY(); iy++){
    double c = h->GetBinContent(ix, iy); if(c == 0.) continue;
    int bx = r->GetXaxis()->FindBin(h->GetXaxis()->GetBinCenter(ix)), by = r->GetYaxis()->FindBin(h->GetYaxis()->GetBinCenter(iy));
    if(bx < 1 || bx > n || by < 1 || by > n) continue;
    r->SetBinContent(bx, by, r->GetBinContent(bx, by) + c);
    r->SetBinError(bx, by, std::hypot(r->GetBinError(bx, by), h->GetBinError(ix, iy)));
  }
  return r;
}

// analysis-axis counts -> extended axis (bins below the analysis range empty),
// bins whose upper edge is at or below floor zeroed; analysis bins falling in one
// merged extended bin are summed (errors in quadrature)
inline TH1D* bToExtAxis(const TH1D *hA, const std::vector<double> &e, double floor, const char *name)
{
  TH1D *o = new TH1D(name, "", (int) e.size() - 1, &e[0]); o->SetDirectory(nullptr); o->Sumw2();
  for(int b = 1; b <= hA->GetNbinsX(); b++){
    if(hA->GetXaxis()->GetBinUpEdge(b) <= floor + 1e-6) continue;
    int j = o->FindBin(hA->GetBinCenter(b));
    if(j < 1 || j > o->GetNbinsX()) continue;
    o->SetBinContent(j, o->GetBinContent(j) + hA->GetBinContent(b)); o->SetBinError(j, std::hypot(o->GetBinError(j), hA->GetBinError(b)));
  }
  return o;
}
// extended axis -> analysis axis (nA edges); an analysis bin that is only part of
// a merged extended bin gets that bin's content x (its width / merged width)
inline TH1D* bToAnaAxis(const TH1D *hE, int nA, const double *ana, const char *name)
{
  TH1D *o = new TH1D(name, "", nA - 1, ana); o->SetDirectory(nullptr); o->Sumw2();
  for(int b = 1; b < nA; b++){
    int j = hE->GetXaxis()->FindFixBin(o->GetBinCenter(b));
    if(j < 1 || j > hE->GetNbinsX()) continue;
    const double f = o->GetBinWidth(b)/hE->GetBinWidth(j);
    o->SetBinContent(b, hE->GetBinContent(j)*f); o->SetBinError(b, hE->GetBinError(j)*f);
  }
  return o;
}
// reported analysis bins: 80-200 GeV with the under/overflow axis, else from 80
inline bool bReportedBin(double lo, double hi){ return lo >= bUfHi - 1e-6 && (!bUnfoldUfOf || hi <= bOfLo + 1e-6); }
// upper edge of the data-driven prior fit: inside the reported range
inline double bPriorFitHi(){ return bUnfoldUfOf ? bOfLo : 300.; }

// Unfold measA (counts, analysis axis) with Rext (extended axes). Truth = gen
// projection of the full matrix; reco columns at or below floor removed from the
// response and from the input. Returns counts on the analysis axis; truthA, if
// given, receives the truth on the analysis axis.
template <class RU, class RB>
TH1D* unfoldOnExt(const TH2D *Rext, const TH1D *measA, int nA, const double *ana, double floor, int N,
                  const char *name, TH1D **truthA = nullptr)
{
  const int n = Rext->GetNbinsX();
  std::vector<double> e(n + 1); for(int i = 0; i <= n; i++) e[i] = Rext->GetXaxis()->GetBinLowEdge(i + 1);
  TH1D *tr = Rext->ProjectionY(Form("%s_tr", name), 1, n); tr->SetDirectory(nullptr);
  TH2D *Rt = (TH2D*) Rext->Clone(Form("%s_Rt", name)); Rt->SetDirectory(nullptr);
  for(int ix = 1; ix <= n; ix++) if(e[ix] <= floor + 1e-6) for(int iy = 0; iy <= n + 1; iy++){ Rt->SetBinContent(ix, iy, 0.); Rt->SetBinError(ix, iy, 0.); }
  TH1D *mm = Rt->ProjectionX(Form("%s_mm", name), 1, n); mm->SetDirectory(nullptr);
  TH1D *in = bToExtAxis(measA, e, floor, Form("%s_in", name));
  RU ru(mm, tr, Rt, Form("%s_ru", name), "");
  RB ub(&ru, in, N); ub.SetVerbose(0);
  TH1D *U = (TH1D*) ub.Hunfold(); U->SetDirectory(nullptr);
  TH1D *o = bToAnaAxis(U, nA, ana, name);
  if(truthA) *truthA = bToAnaAxis(tr, nA, ana, Form("%s_truthA", name));
  delete tr; delete Rt; delete mm; delete in; delete U;
  return o;
}
