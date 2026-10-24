// Two-dimensional b-purity fit in muon ptRel x CSVv2, pp PF jets (2026-10-24, user: "rerun the purity
// measurements on CSVv2"). Copy of fitPtRelDeepCSV2D_pp_PF.C (2026-10-22) with the tagger swapped for
// CSVv2, the tagger of the CMS b-jet shapes analyses (AN-2017/337, AN-2020/029, working point 0.9).
// Same model, likelihood, templates, axis smearing and closure; only the inputs, the category edges
// and the labels differ.
//
// Inputs: the 2026-10-24 pp SingleMuon scan and the 2026-10-24 PYTHIA DiJet scan, both with
// h_muptrel_recoJetPt_csvV2_inclRecoMuonTag_triggerOn[_<flavor>] (headers/AnalysisSetup/
// ptRelDeepCSV3D.h; first z bin [-11, 0) = not evaluated, the forest stores -10). Same jets and
// weights as the 2026-10-22 DeepCSV scans. Each jet-pT window (80-300 GeV) is fitted on its own.
//
// Binning: ptRel 0.5 GeV over 0-fitMax (the model is defined over 0-5 GeV), CSVv2 in K categories
//   "4cat": [not evaluated, 0.54), [0.54, 0.84), [0.84, 0.96), [0.96, 1]  (~ the CMS loose / medium /
//           tight working points 0.5426 / 0.8484 / 0.9535 on the scan's 0.02 grid)
//   "8cat": edges 0.10, 0.30, 0.54, 0.70, 0.84, 0.90, 0.96
//   "cut<x>": fail / pass at x (fitPtRelCSVv2Cut_pp_PF.C)
//
// Model, per window, bin (k = CSVv2 category, j = ptRel bin):
//   mu_kj = N sum_F f_F w_Fk T_F,kj        F = b (flavor creation + g->bb, natural share), c, light
//   T_F,kj   PYTHIA counts / flavor total, axis-smeared by the nominal per-window pp PF sigma
//            (axisSigmaPPPFFor), applied as the smeared / unsmeared ptRel ratio of the scan's
//            axisSmear grid, the same in every CSVv2 category
//   f_b      b fraction of all tagged jets over 0-5 GeV = the PURITY (free)
//   f_c      (1 - f_b) r: r = the c share of the non-b jets, PYTHIA's or free
//   w_Fk     data / PYTHIA tagger shape of flavor F: exp(a_Fk) / sum_m eps_Fm exp(a_Fm), a_F0 = 0,
//            so the flavor total is unchanged and only its split across CSVv2 moves
//   Barlow-Beeston lite for the templates' MC statistics (Conway, arXiv:1103.0354)
// Models:
//   "ptRel"   1D: CSVv2 summed, b / light+c at PYTHIA's c share (the nominal pp model)
//   "MC"      2D, every tagger shape PYTHIA's, c share PYTHIA's
//   "bFree"   2D, the b tagger shape free (a_bk unconstrained), c and light PYTHIA's
//   "full"    2D, b free, c and light tagger shapes with Gaussian priors (sigma 0.3 on a), r free
// Closure (Asimov): pseudo-data = the PYTHIA templates at the data's size, (i) as they are,
// (ii) with the tagger mis-modelled: b in the top category x 0.85 (graded below), light and c in
// every category above the lowest x 1.5 and x 1.2; fitted with every model against the truth.
//
// Output: rootFiles/bPurityResults/ptRelCSVv2_2D_PF_pp_<cats>_fit<max>.root
//         figures/bPurity/PFJets_pp/csvV2_2D/{purity,tagSF,fit_<lo>to<hi>}_<cats>_fit<max>.pdf
// Usage, from src/plots/bPurity/:  root -l -b -q 'fitPtRelCSVv2_2D_pp_PF.C("4cat", 5.0)'   (0-5 GeV nominal since 2026-10-24, user; 3.0 = the chain's edge)
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TH3D.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLine.h"
#include "TLegend.h"
#include "THStack.h"
#include "TGraphErrors.h"
#include "../../../headers/plotting/plotStyle.h"
#include "../../../headers/plotting/gridPads.h"
#include "../../../headers/functions/caloJetBPurityFit_pp.h"   // axisSigmaPPPFFor, axisSmearSigmaPP, projWin

namespace p2d {
const char *repo = "/home/clayton/Analysis/code/bJetRaaAnalysis";
const char *dataName = "rootFiles/scanningOuput/pp/pp_SingleMuon_manualJEC_mu12_pTmu-15to999_tight_deltaR-40_mu12TriggerEfficiencyCorrection_jetTrkMaxFilter_WDecayFilter_2026-10-24.root";
const char *mcName = "rootFiles/scanningOuput/PYTHIA/PYTHIA_DiJet_manualJEC_pThat-15_mu12_pTmu-15to999_tight_vzReweight_jetTrkMaxFilter_removeHYDJETjet0p35CutOnGen_WDecayFilter_axisSmearScan_genMuTagFix_2026-10-24.root";
const char *h3 = "h_muptrel_recoJetPt_csvV2_inclRecoMuonTag_triggerOn";
const char *h2 = "h_muptrel_recoJetPt_inclRecoMuonTag_triggerOn";
const int NW = 5;
const double wLo[NW] = {80., 100., 120., 150., 200.}, wHi[NW] = {100., 120., 150., 200., 300.};
const int NPR = 10;                    // ptRel 0.5 GeV over 0-5 GeV
const int KMAX = 8;
int K = 4;
double catUp[KMAX];                    // upper CSVv2 edge of each category (category 0 also holds "not evaluated")
const int NF = 3;                      // b, c, light
const char *flLabel[NF] = {"b", "c", "light"};
const std::vector<std::vector<const char*>> flMembers = {{"bJets", "bGSJets"}, {"cJets"}, {"dJets", "uJets", "sJets", "gJets"}};
const double sigmaPrior = 0.3;         // prior on a_ck, a_lk ("full")

struct Tpl { double n[KMAX][NPR] = {}, e2[KMAX][NPR] = {}, tot = 0., eps[KMAX] = {}; };
struct Win { double d[KMAX][NPR] = {}; Tpl t[NF]; double truthB = 0., rc = 0., sigma = 0., nData = 0.; };
Win W[NW];

inline bool setCats(const TString &c){
  std::vector<double> e;
  if(c == "4cat") e = {0.54, 0.84, 0.96, 1.01};
  else if(c == "8cat") e = {0.10, 0.30, 0.54, 0.70, 0.84, 0.90, 0.96, 1.01};
  else if(c.BeginsWith("cut")) e = {TString(c(3, c.Length())).Atof(), 1.01};   // "cut0.90": fail / pass (cut on the 0.02 grid)
  else return false;
  K = (int) e.size(); for(int k = 0; k < K; k++) catUp[k] = e[k];
  return true;
}
inline int catOf(double x){ for(int k = 0; k < K; k++) if(x < catUp[k]) return k; return K - 1; }

// counts of the histograms 'names' in the window, ptRel 0.5 GeV x CSVv2 category; ratio[j]
// (optional) scales ptRel bin j (axis smearing)
bool fillCounts(TFile *f, const std::vector<TString> &names, double lo, double hi, double n[][NPR], double e2[][NPR], const double *ratio = nullptr){
  for(const TString &nm : names){
    TH3D *h = nullptr; f->GetObject(nm, h);
    if(!h){ printf("ERROR: missing %s\n", nm.Data()); return false; }
    const int y1 = h->GetYaxis()->FindBin(lo + 1e-6), y2 = h->GetYaxis()->FindBin(hi - 1e-6);
    for(int ix = 1; ix <= h->GetNbinsX(); ix++){
      const double x = h->GetXaxis()->GetBinCenter(ix); if(x >= 5.) continue;
      const int j = (int)(x/0.5);
      for(int iz = 1; iz <= h->GetNbinsZ(); iz++){
        const int k = catOf(h->GetZaxis()->GetBinCenter(iz));
        double c = 0., e = 0.;
        for(int iy = y1; iy <= y2; iy++){ c += h->GetBinContent(ix, iy, iz); e += pow(h->GetBinError(ix, iy, iz), 2); }
        const double r = ratio ? ratio[j] : 1.;
        n[k][j] += c*r; e2[k][j] += e*r*r;
      }
    }
  }
  return true;
}
// smeared / unsmeared ptRel per 0.5 GeV bin from the scan's axisSmear grid, linear in sigma
bool smearRatio(TFile *f, const std::vector<const char*> &fls, double lo, double hi, double sigma, double *ratio){
  int k = 0; while(k < NAxisSmearPP - 2 && sigma > axisSmearSigmaPP[k+1]) k++;
  const double fr = (sigma - axisSmearSigmaPP[k])/(axisSmearSigmaPP[k+1] - axisSmearSigmaPP[k]);
  double s[NPR] = {}, u[NPR] = {};
  for(const char *fl : fls){
    TH1D *a = projWin(f, Form("%s_%s_axisSmear%d", h2, fl, k), lo, hi), *b = projWin(f, Form("%s_%s_axisSmear%d", h2, fl, k+1), lo, hi);
    TH1D *z = projWin(f, Form("%s_%s_axisSmear0", h2, fl), lo, hi);
    if(!a || !b || !z) return false;
    for(int i = 1; i <= z->GetNbinsX(); i++){
      const double x = z->GetBinCenter(i); if(x >= 5.) continue;
      const int j = (int)(x/0.5);
      s[j] += (1. - fr)*a->GetBinContent(i) + fr*b->GetBinContent(i); u[j] += z->GetBinContent(i);
    }
    delete a; delete b; delete z;
  }
  for(int j = 0; j < NPR; j++) ratio[j] = u[j] > 0. ? s[j]/u[j] : 1.;
  return true;
}

// ---------------------------------------------------------------- model
struct Model { const char *name, *label; bool oneD, freeB, freeCL, freeR; };
const Model models[4] = {
  {"ptRel", "#it{p}_{T}^{rel} only",                         true,  false, false, false},
  {"MC",    "2D, PYTHIA tagger shapes",                       false, false, false, false},
  {"bFree", "2D, b tagger shape free",                        false, true,  false, false},
  {"full",  "2D, all tagger shapes floated",                  false, true,  true,  true}};
const Model *cur = nullptr;
Win *cw = nullptr;                     // window being fitted
const double (*cd)[NPR] = nullptr;     // its data (real or pseudo)
int nFit = 6;                          // ptRel bins in the fit (6 = 0-3 GeV)
int Kf = 4;                            // categories the current model sees (1 for "ptRel")
// parameters: [0] f_b, [1] r, [2 .. K] a_b, [K+1 .. 2K-1] a_c, [2K .. 3K-2] a_l
inline int NPar(){ return 2 + 3*(K - 1); }
inline int aIdx(int F, int k){ return 2 + F*(K - 1) + (k - 1); }

// the model's view of one template: categories merged for "ptRel"
inline double tN(int F, int k, int j){ if(Kf > 1) return cw->t[F].n[k][j]; double s = 0; for(int m = 0; m < K; m++) s += cw->t[F].n[m][j]; return s; }
inline double tE(int F, int k, int j){ if(Kf > 1) return cw->t[F].e2[k][j]; double s = 0; for(int m = 0; m < K; m++) s += cw->t[F].e2[m][j]; return s; }
inline double dN(int k, int j){ if(Kf > 1) return cd[k][j]; double s = 0; for(int m = 0; m < K; m++) s += cd[m][j]; return s; }
// tagger-shape weight w_Fk = exp(a_Fk) / sum_m eps_Fm exp(a_Fm)
inline void weights(const double *p, double w[NF][KMAX]){
  for(int F = 0; F < NF; F++){
    if(Kf == 1){ w[F][0] = 1.; continue; }
    double den = 0., e[KMAX];
    for(int k = 0; k < K; k++){ e[k] = k == 0 ? 1. : exp(p[aIdx(F, k)]); den += cw->t[F].eps[k]*e[k]; }
    for(int k = 0; k < K; k++) w[F][k] = e[k]/den;
  }
}
inline void fracs(const double *p, double f[NF]){ const double r = cur->freeR ? p[1] : cw->rc; f[0] = p[0]; f[1] = (1. - p[0])*r; f[2] = (1. - p[0])*(1. - r); }
// expectation and template variance of every in-range bin, normalised to the in-range data
void expectAll(const double *p, double mu[][NPR], double va[][NPR], double comp[][KMAX][NPR] = nullptr){
  double w[NF][KMAX], f[NF]; weights(p, w); fracs(p, f);
  double sMu = 0., sD = 0.;
  for(int k = 0; k < Kf; k++) for(int j = 0; j < nFit; j++){
    mu[k][j] = 0.; va[k][j] = 0.;
    for(int F = 0; F < NF; F++){
      const double c = cw->t[F].tot > 0. ? f[F]*w[F][k]/cw->t[F].tot : 0.;   // a flavor with no tagged PYTHIA jets (light, tight cut) contributes nothing
      mu[k][j] += c*tN(F, k, j); va[k][j] += c*c*tE(F, k, j);
      if(comp) comp[F][k][j] = c*tN(F, k, j);
    }
    sMu += mu[k][j]; sD += dN(k, j);
  }
  const double sc = sMu > 0. ? sD/sMu : 1.;
  for(int k = 0; k < Kf; k++) for(int j = 0; j < nFit; j++){
    mu[k][j] *= sc; va[k][j] *= sc*sc;
    if(comp) for(int F = 0; F < NF; F++) comp[F][k][j] *= sc;
  }
}
// Poisson -lnL with Barlow-Beeston lite; dev = deviance against the saturated model
double nllCore(const double *p, double *dev = nullptr){
  double mu[KMAX][NPR], va[KMAX][NPR]; expectAll(p, mu, va);
  double L = 0., D = 0.;
  for(int k = 0; k < Kf; k++) for(int j = 0; j < nFit; j++){
    double m = TMath::Max(mu[k][j], 1e-9); const double n = dN(k, j);
    if(va[k][j] > 0.){
      const double s2 = va[k][j]/(m*m), a = 1. - m*s2, beta = 0.5*(a + sqrt(a*a + 4.*n*s2));
      const double mb = TMath::Max(m*beta, 1e-9);
      L += mb - n*log(mb) + (beta - 1.)*(beta - 1.)/(2.*s2);
      D += 2.*(mb - n + (n > 0. ? n*log(n/mb) : 0.)) + (beta - 1.)*(beta - 1.)/s2;   // with the MC-statistics penalty
    } else {
      L += m - n*log(m);
      D += 2.*(m - n + (n > 0. ? n*log(n/m) : 0.));
    }
  }
  if(cur->freeCL) for(int F = 1; F < NF; F++) for(int k = 1; k < K; k++) L += 0.5*pow(p[aIdx(F, k)]/sigmaPrior, 2);
  if(dev) *dev = D;
  return L;
}
double nll(const double *p){ return nllCore(p); }

struct Result { double par[2 + 3*(KMAX-1)] = {}, err[2 + 3*(KMAX-1)] = {}, pur = 0., epur = 0., dev = 0.; int ndf = 0, status = -1, hesse = -1;
                double sf[NF][KMAX] = {}, esf[NF][KMAX] = {}; };
Result fitWin(const Model &M, Win &w, const double (*data)[NPR]){
  cur = &M; cw = &w; cd = data; Kf = M.oneD ? 1 : K;
  ROOT::Math::Minimizer *mn = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
  mn->SetPrintLevel(0); mn->SetStrategy(1); mn->SetMaxFunctionCalls(200000);
  mn->SetErrorDef(0.5); mn->SetTolerance(0.1);
  const int np = NPar(); ROOT::Math::Functor fn(&nll, np); mn->SetFunction(fn);
  mn->SetLimitedVariable(0, "fb", TMath::Min(0.9, TMath::Max(0.1, w.truthB)), 0.02, 0., 1.);
  if(M.freeR) mn->SetLimitedVariable(1, "r", w.rc, 0.02, 0., 1.); else mn->SetFixedVariable(1, "r", w.rc);
  int nFree = 1 + (M.freeR ? 1 : 0);
  for(int F = 0; F < NF; F++) for(int k = 1; k < K; k++){
    const bool fr = !M.oneD && (F == 0 ? M.freeB : M.freeCL);
    if(fr){ mn->SetLimitedVariable(aIdx(F, k), Form("a%d_%d", F, k), 0., 0.05, -3., 3.); nFree++; }
    else mn->SetFixedVariable(aIdx(F, k), Form("a%d_%d", F, k), 0.);
  }
  mn->Minimize();
  Result r; r.status = mn->Status(); mn->Hesse(); r.hesse = mn->Status();
  for(int i = 0; i < np; i++){ r.par[i] = mn->X()[i]; r.err[i] = mn->Errors()[i]; }
  r.pur = r.par[0]; r.epur = r.err[0];
  nllCore(r.par, &r.dev); r.ndf = Kf*nFit - 1 - nFree;
  // data / PYTHIA tagger shapes w_Fk and their errors (linear propagation with the covariance)
  double w0[NF][KMAX]; weights(r.par, w0);
  for(int F = 0; F < NF; F++) for(int k = 0; k < K && !M.oneD; k++){
    r.sf[F][k] = w0[F][k];
    double g[2 + 3*(KMAX-1)] = {};
    for(int i = 0; i < np; i++){
      if(mn->IsFixedVariable(i)) continue;
      double q[2 + 3*(KMAX-1)]; for(int m = 0; m < np; m++) q[m] = r.par[m];
      q[i] += 1e-4; double w1[NF][KMAX]; weights(q, w1); g[i] = (w1[F][k] - w0[F][k])/1e-4;
    }
    double v = 0.; for(int i = 0; i < np; i++) for(int m = 0; m < np; m++) if(g[i] != 0. && g[m] != 0.) v += g[i]*g[m]*mn->CovMatrix(i, m);
    r.esf[F][k] = sqrt(TMath::Max(0., v));
  }
  delete mn;
  return r;
}
// Asimov pseudo-data from the templates at the data's size; mult[F][k] mis-models the tagger.
// Returns the true b fraction over 0-5 GeV.
double asimov(const Win &w, const double mult[NF][KMAX], double out[][NPR]){
  double tot[NF] = {}, all = 0.;
  for(int F = 0; F < NF; F++) for(int k = 0; k < K; k++) for(int j = 0; j < NPR; j++) tot[F] += w.t[F].n[k][j]*mult[F][k];
  for(int F = 0; F < NF; F++) all += tot[F];
  for(int k = 0; k < K; k++) for(int j = 0; j < NPR; j++){
    out[k][j] = 0.; for(int F = 0; F < NF; F++) out[k][j] += w.t[F].n[k][j]*mult[F][k]*w.nData/all;
  }
  return tot[0]/all;
}
// templates and data of every window, in the current categories (setCats)
bool buildWindows(TFile *fD, TFile *fM){
  for(int w = 0; w < NW; w++) W[w] = Win();
  for(int w = 0; w < NW; w++){
    Win &X = W[w];
    X.sigma = axisSigmaPPPFFor(wLo[w], wHi[w]);
    double dE2[KMAX][NPR] = {};   // the likelihood uses the data counts only
    if(!fillCounts(fD, {h3}, wLo[w], wHi[w], X.d, dE2)) return false;
    for(int F = 0; F < NF; F++){
      double ratio[NPR]; if(!smearRatio(fM, flMembers[F], wLo[w], wHi[w], X.sigma, ratio)) return false;
      std::vector<TString> nm; for(const char *fl : flMembers[F]) nm.push_back(Form("%s_%s", h3, fl));
      if(!fillCounts(fM, nm, wLo[w], wHi[w], X.t[F].n, X.t[F].e2, ratio)) return false;
      for(int k = 0; k < K; k++) for(int j = 0; j < NPR; j++){ if(X.t[F].n[k][j] < 0.) X.t[F].n[k][j] = 0.; X.t[F].tot += X.t[F].n[k][j]; }
      for(int k = 0; k < K; k++){ double s = 0; for(int j = 0; j < NPR; j++) s += X.t[F].n[k][j]; X.t[F].eps[k] = s/X.t[F].tot; }
    }
    X.truthB = X.t[0].tot/(X.t[0].tot + X.t[1].tot + X.t[2].tot);
    X.rc = X.t[1].tot/(X.t[1].tot + X.t[2].tot);
  }
  for(int w = 0; w < NW; w++){ double s = 0.; for(int k = 0; k < K; k++) for(int j = 0; j < NPR; j++) s += W[w].d[k][j]; W[w].nData = s; }
  return true;
}
} // namespace p2d

void fitPtRelCSVv2_2D_pp_PF(TString cats = "4cat", double fitMax = 5.0){
  using namespace p2d;
  if(!setCats(cats)){ printf("ERROR: categories are 4cat or 8cat\n"); return; }
  nFit = (int)(fitMax/0.5 + 0.5);
  if(nFit < 2 || nFit > NPR){ printf("ERROR: fitMax in 1-5 GeV\n"); return; }
  setJetCollection("PFJets"); initPlotStyle();
  TFile *fD = TFile::Open(Form("%s/%s", repo, dataName)), *fM = TFile::Open(Form("%s/%s", repo, mcName));
  if(!fD || !fM || fD->IsZombie() || fM->IsZombie()){ printf("ERROR: inputs\n"); return; }
  const TString tag = Form("%s_fit%g", cats.Data(), fitMax), outDir = Form("%s/figures/bPurity/PFJets_pp/csvV2_2D/", repo);
  gSystem->mkdir(outDir, true);

  if(!buildWindows(fD, fM)) return;

  printf("\n2D ptRel x CSVv2 fit, pp PF, %s, fit 0-%.1f GeV (model over 0-5 GeV)\n", cats.Data(), fitMax);
  printf("CSVv2 categories (upper edges):"); for(int k = 0; k < K; k++) printf(" %.2f", catUp[k]); printf("\n");
  for(int w = 0; w < NW; w++){
    printf("%3.0f-%3.0f GeV: data %.0f, axis sigma %.4f, PYTHIA b fraction %.3f, c share of non-b %.3f\n   tagger fractions per category  b:",
           wLo[w], wHi[w], W[w].nData, W[w].sigma, W[w].truthB, W[w].rc);
    for(int k = 0; k < K; k++) printf(" %.3f", W[w].t[0].eps[k]); printf("  c:"); for(int k = 0; k < K; k++) printf(" %.3f", W[w].t[1].eps[k]);
    printf("  light:"); for(int k = 0; k < K; k++) printf(" %.3f", W[w].t[2].eps[k]);
    printf("  data:"); for(int k = 0; k < K; k++){ double s = 0; for(int j = 0; j < NPR; j++) s += W[w].d[k][j]; printf(" %.3f", s/W[w].nData); } printf("\n");
  }

  // ---- Asimov closure
  printf("\nASIMOV CLOSURE (purity fitted - true)\n");
  double one[NF][KMAX], stress[NF][KMAX];
  for(int F = 0; F < NF; F++) for(int k = 0; k < K; k++){
    one[F][k] = 1.;
    const double x = K > 1 ? double(k)/(K - 1) : 0.;
    stress[F][k] = F == 0 ? 1. - 0.15*x*x : (k == 0 ? 1. : (F == 1 ? 1.2 : 1.5));
  }
  TH1D *hClos[2][4];
  for(int s = 0; s < 2; s++){
    printf("  %s\n", s == 0 ? "(i) PYTHIA as is" : "(ii) tagger mis-modelled: b top category x0.85 (graded), c x1.2 and light x1.5 above the lowest category");
    for(int m = 0; m < 4; m++){ hClos[s][m] = new TH1D(Form("h_closure%d_%s", s, models[m].name), "", NW, ptEdges + 3); hClos[s][m]->SetDirectory(nullptr); }
    for(int w = 0; w < NW; w++){
      double pd[KMAX][NPR]; const double tr = asimov(W[w], s == 0 ? one : stress, pd);
      printf("    %3.0f-%3.0f truth %.3f:", wLo[w], wHi[w], tr);
      for(int m = 0; m < 4; m++){
        Result r = fitWin(models[m], W[w], pd);
        printf("  %s %+.3f(+-%.3f)", models[m].name, r.pur - tr, r.epur);
        hClos[s][m]->SetBinContent(w + 1, r.pur - tr); hClos[s][m]->SetBinError(w + 1, r.epur);
      }
      printf("\n");
    }
  }

  // ---- data
  printf("\nDATA\n");
  Result R[4][NW];
  TH1D *hPur[4], *hTruth = new TH1D("h_bPurity_PYTHIA", "PYTHIA b fraction", NW, ptEdges + 3);
  hTruth->SetDirectory(nullptr);
  for(int m = 0; m < 4; m++){ hPur[m] = new TH1D(Form("h_bPurity_%s", models[m].name), models[m].label, NW, ptEdges + 3); hPur[m]->SetDirectory(nullptr); }
  for(int w = 0; w < NW; w++){
    hTruth->SetBinContent(w + 1, W[w].truthB);
    printf("  %3.0f-%3.0f GeV\n", wLo[w], wHi[w]);
    for(int m = 0; m < 4; m++){
      R[m][w] = fitWin(models[m], W[w], W[w].d);
      const Result &r = R[m][w];
      hPur[m]->SetBinContent(w + 1, r.pur); hPur[m]->SetBinError(w + 1, r.epur);
      printf("    %-6s migrad %d hesse %d  purity %.3f +- %.3f  r %.3f  dev/ndf %.1f/%d", models[m].name, r.status, r.hesse, r.pur, r.epur, r.par[1], r.dev, r.ndf);
      if(!models[m].oneD){ printf("  b tag data/PYTHIA:"); for(int k = 0; k < K; k++) printf(" %.3f", r.sf[0][k]); }
      if(models[m].freeCL){ printf("  light:"); for(int k = 0; k < K; k++) printf(" %.2f", r.sf[2][k]); }
      printf("\n");
    }
  }
  TFile *fo = TFile::Open(Form("%s/rootFiles/bPurityResults/ptRelCSVv2_2D_PF_pp_%s.root", repo, tag.Data()), "RECREATE");
  for(int m = 0; m < 4; m++){ hPur[m]->Write(); hClos[0][m]->Write(); hClos[1][m]->Write(); }
  hTruth->Write();
  TNamed info("info", Form("pp PF 2D ptRel x CSVv2 fit, %s, fit 0-%.1f GeV, purity over 0-5 GeV. data: %s. templates: %s", cats.Data(), fitMax, dataName, mcName)); info.Write();
  fo->Close();

  // ---- figure 1: purity vs jet pT, the four models, ratio to ptRel only
  const char *info1 = "anti-#it{k}_{T} #it{R} = 0.4 PF jets, |#it{#eta}^{jet}| < 1.6, #it{p}_{T}^{#mu} > 15 GeV";
  const char *hexM[4] = {hexData, okabeHex[5], okabeHex[1], okabeHex[6]};
  const int markM[4] = {markFilledCircle, markFilledSquare, markFilledDiamond, 34};
  {
    TCanvas *cv = new TCanvas("cvP", "", 1000, 1000);
    GridPads g = makeGridPads(cv, 1, 2, 0.14, 0.04, 0.11, 0.07);
    // top 70%, ratio 30%
    const double split = 0.11 + 0.30*(1. - 0.11 - 0.07);
    g.at(0, 0)->SetPad(0., split, 1., 1.); g.at(0, 0)->SetTopMargin(0.07/(1. - split)); g.at(0, 0)->SetBottomMargin(0.);
    g.at(0, 1)->SetPad(0., 0., 1., split); g.at(0, 1)->SetBottomMargin(0.11/split); g.at(0, 1)->SetTopMargin(0.);
    g.at(0, 0)->cd();
    TH1F *fr = gPad->DrawFrame(80., 0., 300., 1.0); fr->GetXaxis()->SetLabelSize(0.); fr->GetYaxis()->SetLabelSize(0.05); fr->GetYaxis()->SetNdivisions(505);
    styleLine(hTruth, "#8C8C8C", 4); hTruth->SetLineStyle(7); hTruth->Draw("HIST SAME");
    TLegend *lg = makeLegend(0.19, 0.03, 0.95, 0.24, 0.042); lg->SetNColumns(2);
    for(int m = 0; m < 4; m++){ styleH(hPur[m], hexM[m], markM[m], 2.4); hPur[m]->SetLineWidth(3); hPur[m]->Draw("E1 X0 SAME"); lg->AddEntry(hPur[m], models[m].label, "pe"); }
    lg->AddEntry(hTruth, "PYTHIA b fraction", "l"); lg->Draw();
    TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.042); t.DrawLatex(0.19, 0.86, info1);
    t.DrawLatex(0.19, 0.80, Form("fit #it{p}_{T}^{rel} < %.0f GeV, %d CSVv2 categories", fitMax, K));
    g.at(0, 1)->cd();
    TH1F *fr2 = gPad->DrawFrame(80., 0.4, 300., 1.3); fr2->GetXaxis()->SetLabelSize(0.11); fr2->GetYaxis()->SetLabelSize(0.11); fr2->GetYaxis()->SetNdivisions(504);
    fr2->GetXaxis()->SetTickLength(0.07);
    TLine *u = new TLine(80., 1., 300., 1.); u->SetLineStyle(2); u->SetLineWidth(2); u->Draw();
    for(int m = 1; m < 4; m++){
      TH1D *r = (TH1D*) hPur[m]->Clone(Form("rat%d", m)); r->Divide(hPur[0]);
      for(int i = 1; i <= NW; i++) r->SetBinError(i, hPur[m]->GetBinError(i)/hPur[0]->GetBinContent(i));
      r->Draw("E1 X0 SAME");
    }
    cv->cd(); TLatex a; a.SetNDC(); a.SetTextFont(42); a.SetTextSize(0.042);
    a.SetTextAlign(31); a.DrawLatex(0.96, 0.015, "#it{p}_{T}^{jet} (GeV)");
    a.SetTextAngle(90); a.DrawLatex(0.045, 0.98, "b purity"); a.DrawLatex(0.045, split - 0.01, "2D / 1D");
    g.cmsLabel(cv, "Internal", "pp (5.02 TeV)");
    savePdfTight(cv, outDir + "purity_" + tag + ".pdf");
  }
  // ---- figure 2: data / PYTHIA tagger shapes per category from the "full" model, b and light
  {
    TCanvas *cv = new TCanvas("cvS", "", 1600, 800);
    GridPads g = makeGridPads(cv, 2, 1, 0.08, 0.02, 0.13, 0.09);
    const char *hexK[KMAX] = {okabeHex[0], okabeHex[1], okabeHex[2], okabeHex[3], okabeHex[4], okabeHex[5], okabeHex[6], okabeHex[7]};
    const int markK[4] = {markFilledCircle, markFilledSquare, markFilledDiamond, 34};
    TLegend *lg = makeLegend(0.08, 0.58, 0.95, 0.86, 0.050); lg->SetNColumns(2);
    for(int iF = 0; iF < 2; iF++){
      const int F = iF == 0 ? 0 : 2;
      g.at(iF, 0)->cd();
      TH1F *fr = gPad->DrawFrame(80., 0., 300., 6.0); g.styleFrame(fr, iF, 0, 0.05);
      TLine *u = new TLine(80., 1., 300., 1.); u->SetLineStyle(2); u->SetLineWidth(2); u->Draw();
      for(int k = 0; k < K; k++){
        TH1D *h = new TH1D(Form("sf%d_%d", F, k), "", NW, ptEdges + 3); h->SetDirectory(nullptr);
        for(int w = 0; w < NW; w++){ h->SetBinContent(w + 1, R[3][w].sf[F][k]); h->SetBinError(w + 1, R[3][w].esf[F][k]); }
        styleH(h, hexK[k % 8], markK[k % 4], 2.0); h->SetLineWidth(3); h->Draw("E1 X0 SAME");
        if(iF == 0) lg->AddEntry(h, k == 0 ? Form("CSVv2 < %.2f", catUp[0]) : (k == K - 1 ? Form("> %.2f", catUp[k-1]) : Form("%.2f-%.2f", catUp[k-1], catUp[k])), "pe");
      }
      g.text(iF, 0, 0.09, Form("%s jets", flLabel[F]), 0.055, 62);
      if(iF == 1) lg->Draw();
    }
    g.titles(cv, "#it{p}_{T}^{jet} (GeV)", "fitted data / PYTHIA fraction in category", 0.045, 0.025, 0.02);
    g.cmsLabel(cv, "Internal", "pp (5.02 TeV)");
    savePdfTight(cv, outDir + "tagSF_" + tag + ".pdf");
  }
  // ---- figure 3: fit in one window, ptRel per CSVv2 category, "full" model
  for(int w : {1, 3}){
    cur = &models[3]; cw = &W[w]; cd = W[w].d; Kf = K;
    const int nSave = nFit; nFit = NPR;   // draw over 0-5 GeV with the fitted parameters, normalised to the in-range data
    double mu[KMAX][NPR], va[KMAX][NPR], comp[NF][KMAX][NPR];
    // normalisation from the fit range, then shapes over 0-5
    nFit = nSave; expectAll(R[3][w].par, mu, va, comp);
    double scIn = 0., muIn = 0.; for(int k = 0; k < K; k++) for(int j = 0; j < nFit; j++) muIn += mu[k][j];
    nFit = NPR; expectAll(R[3][w].par, mu, va, comp);
    for(int k = 0; k < K; k++) for(int j = 0; j < nSave; j++) scIn += mu[k][j];
    const double rs = muIn/scIn;
    nFit = nSave;
    TCanvas *cv = new TCanvas(Form("cvF%d", w), "", 400*K + 200, 900);
    GridPads g = makeGridPads(cv, K, 2, 0.08, 0.02, 0.12, 0.09);
    const double ratioH = 0.30, inner = 1. - 0.12 - 0.09, yS = 0.12 + ratioH*inner;
    for(int k = 0; k < K; k++){
      TPad *pt = g.at(k, 0), *pb = g.at(k, 1);
      double x1 = pt->GetXlowNDC(), x2 = x1 + pt->GetWNDC();
      pt->SetPad(x1, yS, x2, 1.); pt->SetTopMargin(0.09/(1. - yS)); pt->SetBottomMargin(0.);
      pb->SetPad(x1, 0., x2, yS); pb->SetBottomMargin(0.12/yS); pb->SetTopMargin(0.);
      TH1D *hd = new TH1D(Form("fd%d_%d", w, k), "", NPR, 0., 5.), *hc[NF], *ht = new TH1D(Form("ft%d_%d", w, k), "", NPR, 0., 5.);
      hd->SetDirectory(nullptr); ht->SetDirectory(nullptr);
      for(int F = 0; F < NF; F++){ hc[F] = new TH1D(Form("fc%d_%d_%d", w, k, F), "", NPR, 0., 5.); hc[F]->SetDirectory(nullptr); }
      double ymax = 0.;
      for(int j = 0; j < NPR; j++){
        // drawn in thousands of jets
        hd->SetBinContent(j + 1, 1e-3*W[w].d[k][j]); hd->SetBinError(j + 1, 1e-3*sqrt(W[w].d[k][j]));
        for(int F = 0; F < NF; F++) hc[F]->SetBinContent(j + 1, 1e-3*comp[F][k][j]*rs);
        ht->SetBinContent(j + 1, 1e-3*mu[k][j]*rs);
        ymax = TMath::Max(ymax, 1e-3*TMath::Max(W[w].d[k][j], mu[k][j]*rs));
      }
      pt->cd();
      TH1F *fr = gPad->DrawFrame(0., 0., 5., 1.35*ymax); g.styleFrame(fr, k, 0, 0.05); fr->GetXaxis()->SetLabelSize(0.);
      fr->GetYaxis()->SetLabelSize(k == 0 ? 0.05 : 0.); fr->GetYaxis()->SetMaxDigits(3);
      THStack *st = new THStack(Form("st%d_%d", w, k), "");
      const int colF[NF] = {kRed - 4, kGreen + 2, kBlue - 4};   // b, c, light (user, 2026-10-22)
      for(int F = NF - 1; F >= 0; F--){ hc[F]->SetFillColor(colF[F]); hc[F]->SetLineColor(colF[F]); st->Add(hc[F]); }
      st->Draw("HIST SAME");
      styleH(hd, hexData, markFilledCircle, 1.6); hd->SetLineWidth(3); hd->Draw("E1 X0 SAME");
      TLine *edge = new TLine(fitMax, 0., fitMax, 1.35*ymax); edge->SetLineStyle(7); edge->SetLineWidth(3); edge->SetLineColor(kGray + 2);
      if(fitMax < 5.) edge->Draw();
      gPad->RedrawAxis();
      g.text(k, 0, 0.11, k == 0 ? Form("CSVv2 < %.2f", catUp[0]) : (k == K - 1 ? Form("CSVv2 > %.2f", catUp[k-1]) : Form("CSVv2 %.2f-%.2f", catUp[k-1], catUp[k])), 0.05);
      if(k == 0){
        TLegend *lg = makeLegend(0.66, 0.30, 0.99, 0.66, 0.052);
        lg->AddEntry(hd, "data", "pe"); for(int F = 0; F < NF; F++) lg->AddEntry(hc[F], Form("%s jets", flLabel[F]), "f"); lg->Draw();
      }
      pb->cd();
      TH1F *fr2 = gPad->DrawFrame(0., 0.7, 5., 1.3); g.styleFrame(fr2, k, 1, 0.11); fr2->GetYaxis()->SetNdivisions(503);
      fr2->GetXaxis()->SetLabelSize(0.11); fr2->GetYaxis()->SetLabelSize(k == 0 ? 0.11 : 0.); fr2->GetXaxis()->SetTickLength(0.07);
      TLine *u = new TLine(0., 1., 5., 1.); u->SetLineStyle(2); u->SetLineWidth(2); u->Draw();
      TH1D *r = (TH1D*) hd->Clone(Form("fr%d_%d", w, k)); r->Divide(ht);
      for(int j = 1; j <= NPR; j++) r->SetBinError(j, ht->GetBinContent(j) > 0. ? hd->GetBinError(j)/ht->GetBinContent(j) : 0.);
      r->Draw("E1 X0 SAME");
    }
    cv->cd(); TLatex a; a.SetNDC(); a.SetTextFont(42); a.SetTextSize(0.045);
    a.SetTextAlign(31); a.DrawLatex(0.98, 0.01, "muon #it{p}_{T}^{rel} (GeV)");
    a.SetTextAngle(90); a.DrawLatex(0.022, 0.98, "10^{3} jets per 0.5 GeV"); a.DrawLatex(0.022, yS - 0.01, "data / fit");
    g.cmsLabel(cv, "Internal", Form("pp (5.02 TeV),  %.0f < #it{p}_{T}^{jet} < %.0f GeV", wLo[w], wHi[w]));
    savePdfTight(cv, outDir + Form("fit_%.0fto%.0f_", wLo[w], wHi[w]) + tag + ".pdf");
  }
  printf("\nfigures in %s\n", outDir.Data());
}
