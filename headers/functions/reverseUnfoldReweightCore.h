#pragma once
// Helpers of the b-jet "reverse unfolding" study (2026-10-08), shared by
// src/calculateBJetsPerZ/reverseUnfoldReweight_PF_PbPb.C and _pp.C: reco-side
// weights (rho = meas / MC reco, fitted), column / row reweighting of a fine
// (reco x, gen y) response, the gen-side forward fit (genFit), and the closure
// pseudo-data shapes. Analysis axis and reco floor as the PF b chain.
#include <functional>
#include <vector>
#include <cmath>
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TString.h"
#include "bJetUnfoldExt.h"

namespace rvu {

#ifdef BJET_PT_7BINS   // 7-bin variant (caloJetBPurityFitCore.h): 80-200 GeV in the scans' native edges
const int NEdge = 13; const double edges[NEdge] = {50, 60, 70, 80, 90, 100, 110, 120, 135, 150, 200, 300, 500};
#else
const int NEdge = 10; const double edges[NEdge] = {50,60,70,80,100,120,150,200,300,500};   // as the chain
#endif
const double ptReportMin = 80., ptDrawMax = 200.;
const double recoFloor = 70.;     // the chain's PFJets floor
const double fitHi = 300.;        // rho fitted on reco [floor, fitHi]
const double leadEWPtLo = 20.;    // scan event weight: pT clamped to [leadEWPtLo, fitHi] (the response reco axis starts at 20)
const int nIterChain[4] = {4, 4, 3, 3};   // nIterPFJets of the chain, for the closure comparison
const char *cls[5] = {"pp", "C1", "C2", "C3", "C4"};
// closure pseudo-data: gen rows x (g/100)^-(k + t log(g/100)); the curved one is
// outside the genFit "pow" family, so it tests that fit's model dependence
// The last one ("data-like") weights the gen rows to the chain's Bayes-unfolded
// data / MC truth, interpolated linearly in (log g, log ratio) between analysis
// bin centres and frozen outside: a shape in neither fit family.
const int NClo = 4;
const double cloK[NClo] = {0.5, 1.0, 1.0, 0.}, cloT[NClo] = {0., 0., 0.5, 0.};
TH2D* dataLikeTilt(const TH2D *R, TH1D *Unom, TH1D *truthMC, const char *name){
  std::vector<double> lx, ly;
  for(int b = 1; b <= Unom->GetNbinsX(); b++){ double u = Unom->GetBinContent(b), t = truthMC->GetBinContent(b);
    if(u > 0. && t > 0.){ lx.push_back(log(Unom->GetBinCenter(b))); ly.push_back(log(u/t)); } }
  TH2D *o = (TH2D*) R->Clone(name); o->SetDirectory(nullptr);
  for(int iy = 1; iy <= o->GetNbinsY(); iy++){
    double x = log(o->GetYaxis()->GetBinCenter(iy)), y;
    if(x <= lx.front()) y = ly.front(); else if(x >= lx.back()) y = ly.back();
    else { int j = 0; while(lx[j+1] < x) j++; y = ly[j] + (ly[j+1] - ly[j])*(x - lx[j])/(lx[j+1] - lx[j]); }
    double w = exp(y);
    for(int ix = 0; ix <= o->GetNbinsX()+1; ix++){ o->SetBinContent(ix, iy, o->GetBinContent(ix, iy)*w); o->SetBinError(ix, iy, o->GetBinError(ix, iy)*w); }
  }
  return o;
}
TH2D* curvedTilt(const TH2D *R, double k, double t, const char *name){
  TH2D *o = (TH2D*) R->Clone(name); o->SetDirectory(nullptr);
  for(int iy = 1; iy <= o->GetNbinsY(); iy++){ double x = log(o->GetYaxis()->GetBinCenter(iy)/100.), w = exp(-(k + t*x)*x);
    for(int ix = 0; ix <= o->GetNbinsX()+1; ix++){ o->SetBinContent(ix, iy, o->GetBinContent(ix, iy)*w); o->SetBinError(ix, iy, o->GetBinError(ix, iy)*w); } }
  return o;
}

// fine MC reco projection summed into the analysis bins (by fine-bin centre)
TH1D* recoOnEdges(const TH2D *R, const char *name){
  TH1D *o = new TH1D(name, "", NEdge-1, edges); o->SetDirectory(nullptr); o->Sumw2();
  for(int ix = 1; ix <= R->GetNbinsX(); ix++){
    int b = o->FindBin(R->GetXaxis()->GetBinCenter(ix)); if(b < 1 || b > NEdge-1) continue;
    double s = 0., e2 = 0.;
    for(int iy = 1; iy <= R->GetNbinsY(); iy++){ s += R->GetBinContent(ix, iy); e2 += pow(R->GetBinError(ix, iy), 2); }
    o->SetBinContent(b, o->GetBinContent(b) + s); o->SetBinError(b, std::hypot(o->GetBinError(b), sqrt(e2)));
  }
  return o;
}

// rho = meas / MC reco; fit on bins with low edge >= floor and high edge <= fitHi
TF1* fitRho(TH1D *rho, const char *form, const char *name){
  TString f = form;
  TF1 *fn = f == "pow2" ? new TF1(name, "[0]*pow(x/100., [1] + [2]*log(x/100.))", 20., 500.)
                        : new TF1(name, "[0]*pow(x/100., [1])", 20., 500.);
  fn->SetParameters(1., -1.); if(f == "pow2") fn->SetParameter(2, 0.);
  TGraphErrors *g = new TGraphErrors();
  for(int b = 1; b < NEdge; b++){
    double lo = edges[b-1], hi = edges[b];
    if(lo < recoFloor - 1e-6 || hi > fitHi + 1e-6 || rho->GetBinContent(b) <= 0.) continue;
    int n = g->GetN(); g->SetPoint(n, 0.5*(lo + hi), rho->GetBinContent(b)); g->SetPointError(n, 0., rho->GetBinError(b));
  }
  if(f == "pow2"){ TF1 *p1 = new TF1("p1tmp", "[0]*pow(x/100., [1])", 20., 500.); p1->SetParameters(1., -1.); g->Fit(p1, "Q0");
    fn->SetParameters(p1->GetParameter(0), p1->GetParameter(1), 0.); delete p1; }
  g->Fit(fn, "Q0");
  delete g;
  return fn;
}

// fine matrix with each reco column ix scaled by wOf(reco pT of ix)
TH2D* scaleColumns(const TH2D *R, std::function<double(double)> wOf, const char *name){
  TH2D *o = (TH2D*) R->Clone(name); o->SetDirectory(nullptr);
  for(int ix = 0; ix <= o->GetNbinsX()+1; ix++){
    double w = (ix < 1 || ix > o->GetNbinsX()) ? 0. : wOf(o->GetXaxis()->GetBinCenter(ix));
    for(int iy = 0; iy <= o->GetNbinsY()+1; iy++){ o->SetBinContent(ix, iy, o->GetBinContent(ix, iy)*w); o->SetBinError(ix, iy, o->GetBinError(ix, iy)*w); }
  }
  return o;
}

// gen projection (all reco columns) summed into the analysis bins
TH1D* genOnEdges(const TH2D *R, const char *name){
  TH1D *o = new TH1D(name, "", NEdge-1, edges); o->SetDirectory(nullptr); o->Sumw2();
  for(int iy = 1; iy <= R->GetNbinsY(); iy++){
    int b = o->FindBin(R->GetYaxis()->GetBinCenter(iy)); if(b < 1 || b > NEdge-1) continue;
    double s = 0.;
    for(int ix = 1; ix <= R->GetNbinsX(); ix++) s += R->GetBinContent(ix, iy);
    o->SetBinContent(b, o->GetBinContent(b) + s);
  }
  return o;
}
// the same, reco columns below the floor only
TH1D* genOnEdgesBelowFloor(const TH2D *R, const char *name){
  TH1D *o = new TH1D(name, "", NEdge-1, edges); o->SetDirectory(nullptr);
  for(int iy = 1; iy <= R->GetNbinsY(); iy++){
    int b = o->FindBin(R->GetYaxis()->GetBinCenter(iy)); if(b < 1 || b > NEdge-1) continue;
    double s = 0.;
    for(int ix = 1; ix <= R->GetNbinsX(); ix++) if(R->GetXaxis()->GetBinUpEdge(ix) <= recoFloor + 1e-6) s += R->GetBinContent(ix, iy);
    o->SetBinContent(b, o->GetBinContent(b) + s);
  }
  return o;
}

// stat. error of the truth from the meas errors through the binned weights:
// var T(g) = sum_b (M(b,g)/M_reco(b))^2 var meas(b), reco bins above the floor
void setTruthErrors(TH1D *T, const TH2D *R, TH1D *Mreco, TH1D *meas){
  std::vector<std::vector<double>> a(NEdge, std::vector<double>(NEdge, 0.));
  for(int ix = 1; ix <= R->GetNbinsX(); ix++){
    int rb = Mreco->FindBin(R->GetXaxis()->GetBinCenter(ix));
    if(rb < 1 || rb > NEdge-1 || edges[rb-1] < recoFloor - 1e-6 || Mreco->GetBinContent(rb) <= 0.) continue;
    for(int iy = 1; iy <= R->GetNbinsY(); iy++){
      int gb = T->FindBin(R->GetYaxis()->GetBinCenter(iy)); if(gb < 1 || gb > NEdge-1) continue;
      a[rb][gb] += R->GetBinContent(ix, iy)/Mreco->GetBinContent(rb);
    }
  }
  for(int gb = 1; gb < NEdge; gb++){
    double e2 = 0.;
    for(int rb = 1; rb < NEdge; rb++) e2 += pow(a[rb][gb]*meas->GetBinError(rb), 2);
    T->SetBinError(gb, sqrt(e2));
  }
}

struct RU { TH1D *rho = nullptr, *Mreco = nullptr, *T = nullptr, *Tbelow = nullptr; TF1 *fn = nullptr; };

// the method: rho from meas, fit, columns scaled, gen projection.
// below: "extrap" or "freeze"
RU reverseUnfold(const TH2D *R0, TH1D *meas, const char *form, const char *below, const char *tag){
  RU u;
  u.Mreco = recoOnEdges(R0, Form("Mreco_%s", tag));
  u.rho = (TH1D*) meas->Clone(Form("rho_%s", tag)); u.rho->SetDirectory(nullptr);
  u.rho->Divide(u.Mreco);
  for(int b = 1; b < NEdge; b++) if(meas->GetBinContent(b) <= 0. || u.Mreco->GetBinContent(b) <= 0.){ u.rho->SetBinContent(b, 0.); u.rho->SetBinError(b, 0.); }
  u.fn = fitRho(u.rho, form, Form("wFit_%s", tag));
  const bool freeze = TString(below) == "freeze";
  TF1 *fn = u.fn;
  auto wOf = [&](double r){ return fn->Eval(freeze && r < recoFloor ? recoFloor : r); };
  TH2D *Rw = scaleColumns(R0, wOf, Form("Rw_%s", tag));
  u.T = genOnEdges(Rw, Form("truthRW_%s", tag));
  u.Tbelow = genOnEdgesBelowFloor(Rw, Form("truthRWbelow_%s", tag));
  setTruthErrors(u.T, R0, u.Mreco, meas);
  delete Rw;
  return u;
}

// binned rho, reco above the floor only, divided by the MC eps(g) = Bayes N=1 with the MC prior.
// On the extended axis (bJetUnfoldExt.h), returned on the analysis axis.
TH1D* binnedOverEps(const TH2D *R0, TH1D *meas, const char *name){
  const std::vector<double> e = bExtEdges(NEdge, edges); const int n = (int) e.size() - 1;
  TH2D *Re = rebin2DTo(R0, e, Form("%s_Re", name));
  TH1D *in = bToExtAxis(meas, e, recoFloor, Form("%s_in", name));
  TH1D *o = new TH1D(Form("%s_ext", name), "", n, &e[0]); o->SetDirectory(nullptr);
  for(int iy = 1; iy <= n; iy++){
    double num = 0., above = 0., all = 0.;
    for(int ix = 1; ix <= n; ix++){
      double m = Re->GetBinContent(ix, iy); all += m;
      if(e[ix-1] < recoFloor - 1e-6) continue;
      above += m;
      double mr = 0.; for(int jy = 1; jy <= n; jy++) mr += Re->GetBinContent(ix, jy);
      if(mr > 0.) num += m*in->GetBinContent(ix)/mr;
    }
    o->SetBinContent(iy, above > 0. ? num*all/above : 0.);
  }
  TH1D *r = bToAnaAxis(o, NEdge, edges, name);
  delete Re; delete in; delete o;
  return r;
}

// The gen-side variant ("genFit"): weight the response's GEN rows by
// w(g) = A (g/100)^-k [x (g/100)^(-t log(g/100)) for pow2], fit A, k (t) so the
// reweighted reco matches meas on [floor, fitHi], truth = gen projection. This is
// a forward-folding fit: the reco still matches the data, but the weight moves
// the gen spectrum itself instead of keeping MC's gen shape per reco column.
// k (t) on a grid with A analytic; stat. error from the Delta chi2 = 1 range of
// k (A profiled) in quadrature with A's error at fixed k.
struct GF { TH1D *T = nullptr; double A = 0., k = 0., t = 0., chi2 = 0., ek = 0.; int ndf = 0; };
GF genFit(const TH2D *R0, TH1D *meas, const char *form, const char *tag){
  const bool two = TString(form) == "pow2";
  const int ny = R0->GetNbinsY();
  // M[rb][iy]: fine gen row iy summed into analysis reco bin rb; tr[iy]: gen projection
  std::vector<std::vector<double>> M(NEdge, std::vector<double>(ny + 1, 0.));
  std::vector<double> tr(ny + 1, 0.), gc(ny + 1, 0.);
  for(int iy = 1; iy <= ny; iy++){ gc[iy] = R0->GetYaxis()->GetBinCenter(iy);
    for(int ix = 1; ix <= R0->GetNbinsX(); ix++){ double c = R0->GetBinContent(ix, iy); tr[iy] += c;
      int rb = TMath::BinarySearch(NEdge, edges, R0->GetXaxis()->GetBinCenter(ix)) + 1;
      if(rb >= 1 && rb < NEdge) M[rb][iy] += c; } }
  std::vector<int> use; for(int b = 1; b < NEdge; b++) if(edges[b-1] >= recoFloor - 1e-6 && edges[b] <= fitHi + 1e-6 && meas->GetBinContent(b) > 0.) use.push_back(b);
  auto wg = [&](double g, double k, double t){ double x = log(g/100.); return exp(-(k + t*x)*x); };
  auto chi2At = [&](double k, double t, double &A, double &eA){
    std::vector<double> f(NEdge, 0.);
    for(int b : use) for(int iy = 1; iy <= ny; iy++) f[b] += M[b][iy]*wg(gc[iy], k, t);
    double sdf = 0., sff = 0.; for(int b : use){ double s2 = pow(meas->GetBinError(b), 2); sdf += meas->GetBinContent(b)*f[b]/s2; sff += f[b]*f[b]/s2; }
    A = sff > 0. ? sdf/sff : 0.; eA = sff > 0. ? 1./sqrt(sff) : 0.;
    double c2 = 0.; for(int b : use) c2 += pow((meas->GetBinContent(b) - A*f[b])/meas->GetBinError(b), 2);
    return c2; };
  GF r; double best = 1e300, A, eA;
  for(double t = two ? -1. : 0.; t <= (two ? 1. + 1e-9 : 0.); t += 0.02)
    for(double k = -2.; k <= 4. + 1e-9; k += 0.005){ double c = chi2At(k, t, A, eA); if(c < best){ best = c; r.k = k; r.t = t; r.A = A; } }
  r.chi2 = best; r.ndf = (int) use.size() - (two ? 3 : 2);
  double kLo = r.k, kHi = r.k;
  while(kLo > -3. && chi2At(kLo, r.t, A, eA) < best + 1.) kLo -= 0.002;
  while(kHi <  5. && chi2At(kHi, r.t, A, eA) < best + 1.) kHi += 0.002;
  r.ek = 0.5*(kHi - kLo);
  chi2At(r.k, r.t, A, eA);
  auto truthAt = [&](double k, double a, const char *nm){
    TH1D *o = new TH1D(nm, "", NEdge-1, edges); o->SetDirectory(nullptr);
    for(int iy = 1; iy <= ny; iy++){ int b = o->FindBin(gc[iy]); if(b >= 1 && b < NEdge) o->SetBinContent(b, o->GetBinContent(b) + a*tr[iy]*wg(gc[iy], k, r.t)); }
    return o; };
  double Alo, Ahi; chi2At(kLo, r.t, Alo, eA); chi2At(kHi, r.t, Ahi, eA); chi2At(r.k, r.t, A, eA);
  r.T = truthAt(r.k, r.A, Form("truthGenFit_%s", tag));
  TH1D *tl = truthAt(kLo, Alo, "gfLo"), *th = truthAt(kHi, Ahi, "gfHi");
  for(int b = 1; b < NEdge; b++){ double v = r.T->GetBinContent(b); if(v <= 0.) continue;
    r.T->SetBinError(b, std::hypot(0.5*fabs(th->GetBinContent(b) - tl->GetBinContent(b)), v*eA/A)); }
  delete tl; delete th;
  return r;
}

TH1D* ratioH(TH1D *n, TH1D *d, const char *name){
  TH1D *r = (TH1D*) n->Clone(name); r->SetDirectory(nullptr);
  for(int b = 1; b <= r->GetNbinsX(); b++){
    double a = n->GetBinContent(b), c = d->GetBinContent(b);
    if(a <= 0. || c <= 0.){ r->SetBinContent(b, 0.); r->SetBinError(b, 0.); continue; }
    r->SetBinContent(b, a/c); r->SetBinError(b, a/c*n->GetBinError(b)/a);   // numerator's relative error only
  }
  return r;
}

void blankUnreported(TH1D *h){ for(int b = 1; b < NEdge; b++) if(edges[b-1] < ptReportMin - 1e-6 || h->GetBinContent(b) <= 0.){ h->SetBinContent(b, -999.); h->SetBinError(b, 0.); } }

} // namespace rvu
using namespace rvu;
