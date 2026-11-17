// Analytic jet-energy-scale shift of a measured spectrum (2026-10-21). Jets at pT move to
// pT (1 +/- delta(pT)), delta the AK4PF JEC uncertainty averaged uniformly over |eta| < 1.6; the
// shifted bin [a, b] holds the old spectrum's [a/(1 +/- d), b/(1 +/- d)], integrated over a
// log-log interpolation of the density between bin centres.
// Users: inclusiveJetsPerZSys_PF.C (jeuDataShape cross-check; closure on the simulated spectrum
// against the scan's per-jet shift: 1-2%), the b-jet JEU of unfoldBJetSpectrum_caloJets_pp.C and
// bJetsPerZSys_PF_PbPb.C (the measured b spectrum shifted and unfolded with the b response).
#pragma once
#include <vector>
#include <cmath>
#include "TH1D.h"
#include "TMath.h"
#include "/home/clayton/Analysis/code/bJetRaaAnalysis/JetEnergyCorrections/JetUncertainty.h"

// delta(pT) on a 1 GeV grid, the JEC uncertainty (upper value) averaged uniformly over |eta| < 1.6
inline std::vector<double> jeuDelta(const char *uncFile){
  JetUncertainty U(uncFile);
  std::vector<double> d(701, 0.);
  for(int pt = 1; pt <= 700; pt++){
    double s = 0.; int n = 0;
    for(double eta = -1.55; eta < 1.6; eta += 0.1){ U.SetJetPT(pt); U.SetJetEta(eta); U.SetJetPhi(0.); s += U.GetUncertainty().second; n++; }
    d[pt] = s/n;
  }
  return d;
}
// Density model of a binned spectrum: each bin's AVERAGE density placed at its Lafferty-Wyatt
// point x_LW (where a local power law with the slope of the neighbouring bins equals the bin
// average; iterated), then log-log interpolation between the points, linear where a neighbour is
// <= 0, power-law extrapolation beyond the ends. On wide bins of a steep spectrum the bin centre
// is far from x_LW, and centre placement overestimated the shift (b-jet bins, 2026-10-21).
struct DensModel { std::vector<double> x, y; };
inline DensModel densModel(TH1D *h){
  DensModel m; const int nb = h->GetNbinsX();
  std::vector<double> a(nb), b(nb), y(nb), x(nb);
  for(int i = 0; i < nb; i++){ a[i] = h->GetXaxis()->GetBinLowEdge(i+1); b[i] = h->GetXaxis()->GetBinUpEdge(i+1); y[i] = h->GetBinContent(i+1)/(b[i] - a[i]); x[i] = 0.5*(a[i] + b[i]); }
  for(int it = 0; it < 5; it++){
    std::vector<double> xn = x;
    for(int i = 0; i < nb; i++){
      const int j0 = i > 0 ? i - 1 : i, j1 = i < nb - 1 ? i + 1 : i;
      if(j0 == j1 || y[j0] <= 0. || y[j1] <= 0. || a[i] <= 0.) continue;
      const double n = -(log(y[j1]) - log(y[j0]))/(log(x[j1]) - log(x[j0]));   // local x^-n
      double avg;   // average of x^-n over [a, b]
      if(fabs(n - 1.) < 1e-6) avg = log(b[i]/a[i])/(b[i] - a[i]);
      else avg = (pow(b[i], 1. - n) - pow(a[i], 1. - n))/((1. - n)*(b[i] - a[i]));
      if(avg > 0. && fabs(n) > 1e-6) xn[i] = pow(avg, -1./n);
    }
    x = xn;
  }
  m.x = x; m.y = y;
  return m;
}
inline double densAt(const DensModel &m, double x){
  const int nb = m.x.size();
  int i = 0; while(i < nb - 2 && x > m.x[i+1]) i++;
  const double x0 = m.x[i], x1 = m.x[i+1], y0 = m.y[i], y1 = m.y[i+1];
  if(y0 > 0. && y1 > 0. && x > 0.) return exp(log(y0) + (log(y1) - log(y0))*(log(x) - log(x0))/(log(x1) - log(x0)));
  if(x < x0 || x > x1) return 0.;
  return y0 + (y1 - y0)*(x - x0)/(x1 - x0);
}
// jets at pT move to pT (1 + sgn delta): the new bin [a, b] holds the old [a/(1+sgn d), b/(1+sgn d)]
inline TH1D* jeuShiftSpectrum(TH1D *m, const std::vector<double> &delta, int sgn, const char *name){
  const DensModel dm = densModel(m);
  TH1D *o = (TH1D*) m->Clone(name); o->SetDirectory(nullptr); o->Reset();
  for(int j = 1; j <= o->GetNbinsX(); j++){
    const double a = o->GetXaxis()->GetBinLowEdge(j), b = o->GetXaxis()->GetBinUpEdge(j);
    const int ip = TMath::Max(1, TMath::Min(700, (int) std::lround(0.5*(a + b))));
    const double f = 1. + sgn*delta[ip], lo = a/f, hi = b/f;
    const int nStep = 200; double sum = 0.;
    for(int k = 0; k < nStep; k++){ const double x = lo + (k + 0.5)*(hi - lo)/nStep; sum += densAt(dm, x); }
    const double cnt = sum*(hi - lo)/nStep;
    const double c0 = m->GetBinContent(j);
    o->SetBinContent(j, cnt);
    o->SetBinError(j, c0 != 0. ? fabs(cnt)*m->GetBinError(j)/fabs(c0) : 0.);
  }
  return o;
}

// ---------------------------------------------------------------------------------------------------------
// FIT-BASED shift (2026-11-16, user: "shouldn't we be fitting ... to isolate the JES effect? else you are folding in
// all the uncertainties of the data spectra"). The analytic shift above follows the LOCAL slope of the measured
// spectrum, so bin-to-bin fluctuations of the measured spectrum (statistics, purity of single bins) enter the shift
// (e.g. +3.5% instead of ~+17% in one 50-80% b-jet bin). Here the measured density is fitted with a smooth
// log-polynomial, ln(dN/dpT) = sum_k p_k L^k, L = ln(pT / 100 GeV) (order 3 by default; fitted to the bin
// AVERAGES, TH1::Fit option "I", over every bin with content > 0), and each bin's shift RATIO is taken from the fit,
//   R_j = int_{a/(1+s d)}^{b/(1+s d)} f  /  int_a^b f,     shifted content = measured content x R_j,
// so the normalisation and the bin-to-bin content stay the data's and only the slope that sets the shift is
// smoothed. The fit keeps the data's slope (quenching included), unlike a shift of the simulated spectrum.
// jeuShiftMeasured() is what the callers use: the fit by default, the analytic shift of the measured spectrum with
// BJET_JEU_SHAPE=measured (the method 2026-10-21 to 2026-11-16).
#include "TF1.h"
#include <cstdlib>
inline bool jeuShapeFit(){ const char *e = std::getenv("BJET_JEU_SHAPE"); return !(e && TString(e) == "measured"); }
struct JeuFitInfo { double chi2 = 0.; int ndf = 0; double lo = 0., hi = 0.; int order = 3; std::vector<double> p; };
inline JeuFitInfo &jeuLastFit(){ static JeuFitInfo i; return i; }
inline double jeuFitDens(const std::vector<double> &p, double x){ const double L = log(x/100.); double s = 0., Lk = 1.; for(double c : p){ s += c*Lk; Lk *= L; } return exp(s); }
inline double jeuFitInt(const std::vector<double> &p, double a, double b){   // Simpson, 200 intervals
  const int n = 200; const double h = (b - a)/n; double s = jeuFitDens(p, a) + jeuFitDens(p, b);
  for(int k = 1; k < n; k++) s += (k % 2 ? 4. : 2.)*jeuFitDens(p, a + k*h);
  return s*h/3.;
}
// the fit of a counts histogram's density; parameters returned, chi2 / ndf in jeuLastFit()
inline std::vector<double> jeuFitSpectrum(TH1D *m, int order = 3){
  static int nCall = 0; nCall++;
  TH1D *d = (TH1D*) m->Clone(Form("jeuFitDens_%d", nCall)); d->SetDirectory(nullptr);
  double lo = 1e9, hi = 0.; int nUse = 0;
  for(int j = 1; j <= d->GetNbinsX(); j++){ const double w = d->GetBinWidth(j), c = d->GetBinContent(j);
    d->SetBinContent(j, c/w); d->SetBinError(j, d->GetBinError(j)/w);
    if(c > 0.){ lo = std::min(lo, d->GetXaxis()->GetBinLowEdge(j)); hi = std::max(hi, d->GetXaxis()->GetBinUpEdge(j)); nUse++; } }
  if(nUse < order + 2) order = std::max(1, nUse - 2);
  TString form = "exp([0]"; for(int k = 1; k <= order; k++) form += Form(" + [%d]*pow(log(x/100.), %d)", k, k); form += ")";
  TF1 f(Form("jeuFitF_%d", nCall), form, lo, hi);
  // start: a straight line in log-log through the two ends
  int j0 = d->FindBin(lo + 1e-3), j1 = d->FindBin(hi - 1e-3);
  while(j1 > j0 && d->GetBinContent(j1) <= 0.) j1--;
  const double x0 = d->GetBinCenter(j0), x1 = d->GetBinCenter(j1), y0 = d->GetBinContent(j0), y1 = d->GetBinContent(j1);
  const double n = (y0 > 0. && y1 > 0. && x1 > x0) ? log(y1/y0)/log(x1/x0) : -5.;
  f.SetParameter(0, y0 > 0. ? log(y0) - n*log(x0/100.) : 0.); f.SetParameter(1, n); for(int k = 2; k <= order; k++) f.SetParameter(k, 0.);
  d->Fit(&f, "IRQN0");
  JeuFitInfo &I = jeuLastFit(); I.chi2 = f.GetChisquare(); I.ndf = f.GetNDF(); I.lo = lo; I.hi = hi; I.order = order;
  I.p.assign(order + 1, 0.); for(int k = 0; k <= order; k++) I.p[k] = f.GetParameter(k);
  delete d;
  return I.p;
}
inline TH1D* jeuShiftSpectrumFit(TH1D *m, const std::vector<double> &delta, int sgn, const char *name, int order = 3){
  const std::vector<double> p = jeuFitSpectrum(m, order);
  TH1D *o = (TH1D*) m->Clone(name); o->SetDirectory(nullptr);
  for(int j = 1; j <= o->GetNbinsX(); j++){
    const double a = o->GetXaxis()->GetBinLowEdge(j), b = o->GetXaxis()->GetBinUpEdge(j);
    const int ip = TMath::Max(1, TMath::Min(700, (int) std::lround(0.5*(a + b))));
    const double f = 1. + sgn*delta[ip];
    const double den = jeuFitInt(p, a, b), R = den > 0. ? jeuFitInt(p, a/f, b/f)/den : 1.;
    o->SetBinContent(j, m->GetBinContent(j)*R); o->SetBinError(j, m->GetBinError(j)*R);
  }
  return o;
}


// LOCAL variant for finely binned spectra over a wide range (the inclusive spectrum, 5 GeV bins, 0-500 GeV, with its
// low-pT turnover and trigger-stitch steps, which no single log-polynomial describes; 2026-11-16): for each bin, a
// quadratic in L = ln(pT / 100 GeV) fitted by least squares to ln(density) of the bins within |ln(pT/pT_j)| < hw
// (relative errors floored at relFloor, so the shape, not a single high-statistics bin, sets it), and the bin's shift
// ratio taken from that local quadratic. Bins with no content keep zero.
inline std::vector<double> jeuLocalQuad(TH1D *m, int j, double hw, double relFloor){
  const double xj = m->GetBinCenter(j);
  double S[5] = {0}, T[3] = {0};   // sums of w L^k (k = 0..4) and w y L^k (k = 0..2)
  int n = 0;
  for(int i = 1; i <= m->GetNbinsX(); i++){
    const double c = m->GetBinContent(i); if(c <= 0.) continue;
    const double x = m->GetBinCenter(i); if(fabs(log(x/xj)) > hw) continue;
    const double y = log(c/m->GetBinWidth(i)), r = m->GetBinError(i)/c, w = 1./(r*r + relFloor*relFloor), L = log(x/100.);
    double Lk = 1.; for(int k = 0; k < 5; k++){ S[k] += w*Lk; if(k < 3) T[k] += w*y*Lk; Lk *= L; } n++;
  }
  std::vector<double> p(3, 0.);
  if(n < 3){ p[0] = log(std::max(1e-30, m->GetBinContent(j)/m->GetBinWidth(j))) ; return p; }   // too few points: flat
  // normal equations, 3 x 3, Cramer
  const double a[3][3] = {{S[0], S[1], S[2]}, {S[1], S[2], S[3]}, {S[2], S[3], S[4]}};
  auto det3 = [](const double M[3][3]){ return M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1]) - M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0]) + M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0]); };
  const double D = det3(a); if(fabs(D) < 1e-300) return p;
  for(int k = 0; k < 3; k++){ double M[3][3]; for(int r = 0; r < 3; r++) for(int c = 0; c < 3; c++) M[r][c] = (c == k) ? T[r] : a[r][c]; p[k] = det3(M)/D; }
  return p;
}
inline TH1D* jeuShiftSpectrumLocal(TH1D *m, const std::vector<double> &delta, int sgn, const char *name, double hw = 0.3, double relFloor = 0.01){
  TH1D *o = (TH1D*) m->Clone(name); o->SetDirectory(nullptr);
  for(int j = 1; j <= o->GetNbinsX(); j++){
    if(m->GetBinContent(j) <= 0.) continue;
    const std::vector<double> p = jeuLocalQuad(m, j, hw, relFloor);
    const double a = o->GetXaxis()->GetBinLowEdge(j), b = o->GetXaxis()->GetBinUpEdge(j);
    const int ip = TMath::Max(1, TMath::Min(700, (int) std::lround(0.5*(a + b))));
    const double f = 1. + sgn*delta[ip];
    const double den = jeuFitInt(p, a, b), R = den > 0. ? jeuFitInt(p, a/f, b/f)/den : 1.;
    o->SetBinContent(j, m->GetBinContent(j)*R); o->SetBinError(j, m->GetBinError(j)*R);
  }
  return o;
}

// what the callers use: b-jet spectra (few wide bins, 50-300 GeV) -> global log-cubic fit; local = true (the finely
// binned inclusive spectrum) -> local quadratic, half-width 0.3 in ln pT; BJET_JEU_SHAPE=measured -> the analytic shift
// of the measured spectrum itself (before 2026-11-16)
inline TH1D* jeuShiftMeasured(TH1D *m, const std::vector<double> &delta, int sgn, const char *name, bool local = false){
  if(!jeuShapeFit()) return jeuShiftSpectrum(m, delta, sgn, name);
  return local ? jeuShiftSpectrumLocal(m, delta, sgn, name) : jeuShiftSpectrumFit(m, delta, sgn, name);
}
