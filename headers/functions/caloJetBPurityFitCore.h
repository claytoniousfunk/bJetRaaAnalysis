// Collision-system-independent machinery for the CALO-jet b-purity ptRel
// template fits: jet-pT windows, ptRel projection, template struct and the
// RooFit fit. Included by caloJetBPurityFit_pp.h and caloJetBPurityFit_PbPb.h,
// which add their own inputs and template construction.
#pragma once
#include <vector>
#include <string>
#include <functional>
#include <cstdlib>
#include "TFile.h"
#include "TString.h"
#include "TSystem.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TBox.h"
#include "TLatex.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooFormulaVar.h"
#include "RooGaussian.h"
#include <memory>

// Jet-pT windows, matched to the unfolding binning. 50-60, 60-70 and 70-80 are
// UNDERFLOW bins kept for the unfolding and not reported; the first reported
// bin starts at ptReportMin. Edges sit on the 5 GeV bins of the input
// histograms.

// 7-BIN VARIANT (2026-11-16, user: "run the b-jet chain in the 7 bins as a variant"; study studyJetPtBinning_2D_PF.C):
// with the preprocessor symbol BJET_PT_7BINS defined BEFORE the macro is loaded,
//   root -l -b -q -e '#define BJET_PT_7BINS 1' 'macro.C(...)'
// the analysis axis has 80-200 GeV in 7 bins (the scans' native edges) instead of 4. Compile-time because NPt sizes
// arrays. bJetChain_caloJets_PbPb.C and bJetsPerZSys_PF_PbPb.C switch their own edges on the same symbol, and every
// chain macro refuses to run unless its output tag contains "7bins" (ptBinGuard, purityVariant.h).
#ifdef BJET_PT_7BINS
const int NPt = 12;
const double ptEdges[NPt+1] = {50, 60, 70, 80, 90, 100, 110, 120, 135, 150, 200, 300, 500};
#else
const int NPt = 9;
const double ptEdges[NPt+1] = {50, 60, 70, 80, 100, 120, 150, 200, 300, 500};
#endif
// output-name tag of the jet-pT axis (2026-11-16): "_7bins" with BJET_PT_7BINS, "" otherwise -- study macros append it to
// files that would otherwise be overwritten by a 7-bin run (closure2D_HF_PF.C, plotFitQuality_1DTagVs2DMC_PF.C,
// fitPtRelCSVv2_2D_ppPbPb_PF.C, plotFit2D_HF_PF.C)
inline const char* ptBinTag(){
#ifdef BJET_PT_7BINS
  return "_7bins";
#else
  return "";
#endif
}
const double ptReportMin = 80.;

// Minimum data entries in 0-5 GeV ptRel for a window to be fitted at all.
// Below it the window is skipped and left empty in every output, rather than
// fitted on nothing.
const double minFitEntries = 50.;

// fitHi is the upper end of the ptRel fit DOMAIN (templates, data, RooFit
// variable, and the jets the purity refers to). Nominal 5 GeV; setPtRelMax(10.)
// (2026-10-08 study) covers the whole scan axis, so the purity refers to every
// tagged jet the chain counts. Call it before building templates.
const double fitLo = 0.0;
double fitHi = 5.0;
inline void setPtRelMax(double m){ fitHi = m; }

int uidP = 0;

// Jet collection and the scan paths that follow it: jetCollectionInputs.h
// (setJetCollection, collTag, collLabel, ppSmearPath ...).
#include "jetCollectionInputs.h"
#include "purityVariant.h"

// PYTHIA templates with the reco jet axis Gaussian-smeared by sigma, the same
// sigma in eta and in phi (PYTHIA_scan.C doAxisSmearScan), histograms
// <name>_<flavor>Jets_axisSmear<k>. Used by pp (axis resolution of the calo
// jets in data, plotAxisSmearTune_caloJets_pp.C) and by PbPb (smeared pp light).
// Calo jets only: ppSmearPath is empty for caloMatchedPF (no such scan).
const int NAxisSmearPP = 9;   // must match config_PYTHIA.h of that scan
const double axisSmearSigmaPP[NAxisSmearPP] = {0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.10};

inline TFile* ppSmearFile()
{
  static TFile *f = nullptr; static TString opened;
  if(!ppSmearPath[0]){ printf("ERROR: no axis-smeared template scan for the %s collection\n", collTag()); return nullptr; }
  if(f && opened != ppSmearPath){ f = nullptr; }
  if(!f){ f = TFile::Open(ppSmearPath); opened = ppSmearPath; if(!f || f->IsZombie()){ printf("ERROR: cannot open %s\n", ppSmearPath); f = nullptr; } }
  return f;
}

// ptRel projection in a jet-pT window. Window edges are on 5 GeV bin
// boundaries; FindBin(hi - eps) keeps the window from picking up the next bin.
TH1D* projWin(TFile *f, const char *name, double lo, double hi)
{
  TH2D *H = nullptr; f->GetObject(name, H);
  if(!H){ printf("ERROR: missing %s\n", name); return nullptr; }
  TH1D *p = H->ProjectionX(Form("pw%d", uidP++), H->GetYaxis()->FindBin(lo + 1e-6),
                           H->GetYaxis()->FindBin(hi - 1e-6));
  p->SetDirectory(nullptr);
  // weighted-MC bins can come out slightly negative; a template cannot
  for(int i = 0; i <= p->GetNbinsX()+1; i++) if(p->GetBinContent(i) < 0.) p->SetBinContent(i, 0.);
  return p;
}

// ptRel fit binning over 0-5 GeV: 0.1 GeV (50 bins) for calo; fitBinPF for
// PFJets (user, 2026-10-07: wider bins for the sparse central and high-pT
// windows; 0.25 GeV, then 0.5 GeV). The scans book 0.1 GeV bins, so the PF
// binning merges them (0.1-wide bin centres never sit on a 0.25 or 0.5 edge).
// The lower-fit-edge scans must use PF bin edges only: 0 and fitBinPF.
double fitBinPF = 0.5;   // NOMINAL 0.5; settable for studies only (setFitBinPF)
// 2026-10-21 (user): studies may refit at the native 0.1 GeV binning, where the
// low-ptRel peak is resolved (0-0.5 GeV is one bin at 0.5). Call before any template
// or data projection; the nominal chain never calls it. Must divide 0.1-GeV scan bins.
inline void setFitBinPF(double w){ fitBinPF = w; printf("PF fit binning: %.2f GeV (%d bins over 0-%.1f GeV)%s\n", w, (int)(fitHi/w + 0.5), fitHi, w == 0.5 ? "" : "  -- NOT the nominal 0.5"); }
inline TString fitBinTag(){ return fitBinPF == 0.5 ? TString("") : TString::Format("_bin%g", fitBinPF).ReplaceAll(".", "p"); }   // "" at the nominal, e.g. _bin0p1
inline int nFit05(){ return jetColl == kPFJets ? (int)(fitHi/fitBinPF + 0.5) : (int)(fitHi/0.1 + 0.5); }
// lower-fit-edge variation for the systematics: one bin up for PF
inline double fitLowStepPF(){ return fitBinPF; }
// display / shape-comparison edges matching the fit binning
const double edgeFit05PF[] = {0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,5.0};
inline TH1D* rebinFit05(TH1D *h, const char *name)
{
  const int nPF = (int)(sizeof(edgeFit05PF)/sizeof(double)) - 1;
  static const double edgeCalo[] = {0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.4,2.8,3.4,4.2,5.0};
  const int nCalo = (int)(sizeof(edgeCalo)/sizeof(double)) - 1;
  TH1D *r = (TH1D*) h->Rebin(jetColl == kPFJets ? nPF : nCalo, name, jetColl == kPFJets ? edgeFit05PF : edgeCalo);
  r->SetDirectory(nullptr);
  return r;
}

// DISPLAY-ONLY override of the range05 bin width (0: the fit binning, nFit05()).
// Set it around a second dataWin / buildTemplatesPbPb call to redraw a fit at
// finer bins (componentsAt below); never around a fit, the fit binning is
// fitBinPF. 0.1-GeV scan bins merge cleanly into 0.2.
double range05BinOverride = 0.;
// copy of the 0-fitHi part (nominal 0-5 GeV), nFit05() bins, as the fit sees it
TH1D* range05(TH1D *h, const char *name)
{
  const int n = range05BinOverride > 0. ? (int)(fitHi/range05BinOverride + 0.5) : nFit05();
  TH1D *o = new TH1D(name, "", n, 0., fitHi); o->SetDirectory(nullptr); o->Sumw2();
  for(int i = 1; i <= h->GetNbinsX(); i++){
    const double x = h->GetBinCenter(i); if(x < 0. || x >= fitHi) continue;
    const int j = o->FindBin(x);
    o->SetBinContent(j, o->GetBinContent(j) + h->GetBinContent(i));
    o->SetBinError(j, std::hypot(o->GetBinError(j), h->GetBinError(i)));
  }
  return o;
}

// Nominal absolute shift of the b template's gluon-splitting (bGS) share above
// PYTHIA's natural share (bGS_multiplier convention: +0.175). PFJets: 0 -- the
// natural share (user, 2026-10-07: the boost pulled the central purity down at
// high pT, where a free fit prefers a LOWER share than PYTHIA's). Calo: +0.175.
// The bGS systematic scans +/-0.025 around this value.
// PURITY VARIANT switch and file suffix: purityVariant.h
inline double bGSShiftNominal(){ if(purityVariantOn()) return 0.175; return jetColl == kPFJets ? 0. : 0.175; }
// relative (x natural) bGS convention: PFJets nominal only
inline bool bGSRelative(){ return jetColl == kPFJets && !purityVariantOn(); }
// histogram-name tag of a shift: 0.150 -> "0p150", -0.025 -> "m0p025"
inline TString bGSTag(double s){ return Form("%s0p%03.0f", s < 0. ? "m" : "", 1000.*fabs(s)); }

// Toy jet-axis smearing of a ptRel template (PFJets; no smeared-template scans
// exist): ptRel ~ pT_mu x DeltaR; draw ptRel from pr and pT_mu from mp (redrawn
// while DeltaR > 0.4), shift the axis by Gaussian (sigma, sigma) in (eta, phi) at
// a random azimuth, drop muons leaving the 0.4 cone, fill pT_mu x DeltaR'. The
// result keeps pr's integral and each bin pr's relative error. Validated against
// the 2025 PF sigma = 0.005 axisSmear scan (caloJetBPurityFit_pp.h).
inline TH1D* axisSmearToyHist(TH1D *pr, TH1D *mp, double sigma, unsigned seed, int N = 200000)
{
  if(!pr || !mp || pr->Integral() <= 0. || mp->Integral() <= 0. || sigma <= 0.) return (TH1D*) pr->Clone(Form("toyc%d", uidP++));
  TH1D *o = (TH1D*) pr->Clone(Form("toy%d", uidP++)); o->SetDirectory(nullptr); o->Reset();
  TRandom3 r(seed);
  for(int i = 0; i < N; i++){
    double prl = pr->GetRandom(&r), pt = mp->GetRandom(&r); int k = 0;
    while((pt <= 0. || prl/pt > 0.4) && k++ < 100) pt = mp->GetRandom(&r);
    if(pt <= 0. || prl/pt > 0.4) continue;
    double dr = prl/pt, a = r.Uniform(0., TMath::TwoPi());
    double x = dr*cos(a) + r.Gaus(0., sigma), y = dr*sin(a) + r.Gaus(0., sigma), d2 = sqrt(x*x + y*y);
    if(d2 < 0.4) o->Fill(pt*d2);
  }
  if(o->Integral() > 0.) o->Scale(pr->Integral()/o->Integral());
  for(int b = 0; b <= o->GetNbinsX()+1; b++){
    double c0 = pr->GetBinContent(b), rel = c0 > 0. ? pr->GetBinError(b)/c0 : 0.;
    o->SetBinError(b, o->GetBinContent(b)*rel);
  }
  return o;
}
// axisToySeedCommon (default false): the same seed for every sigma, so templates smeared
// with different sigma share their random numbers and vary smoothly with sigma (template
// morphing in profileLikelihoodFit_PF_PbPb.C, 2026-10-21)
bool axisToySeedCommon = false;
inline unsigned axisToySeed(const char *tag, double lo, double hi, int idx, double sigma)
{
  unsigned s = 1000003u*(unsigned)(lo + 0.5) + 7919u*(unsigned)(hi + 0.5) + 10007u*(unsigned)idx + (axisToySeedCommon ? 0u : 101u*(unsigned)(1e5*sigma + 0.5));
  for(const char *c = tag; *c; c++) s = s*31u + (unsigned)*c;
  return s;
}

inline bool fittable(TH1D *data05){ return data05 && data05->Integral() >= minFitEntries; }

// Grey band over the underflow bins (ptEdges[0] .. ptReportMin), drawn in the
// current pad's user coordinates between ymin and ymax, with a small label
// starting at labelY (pass labelY below ymin for no label).
inline void drawUnderflowBand(double ymin, double ymax, double labelY, double textSize = 0.030)
{
  TBox *b = new TBox(ptEdges[0], ymin, ptReportMin, ymax);
  b->SetFillColorAlpha(kGray, 0.45); b->SetLineWidth(0); b->Draw();
  if(labelY > ymin){
    // vertical: the band is too narrow for the word written horizontally
    TLatex *t = new TLatex(0.5*(ptEdges[0] + ptReportMin), labelY, "underflow");
    t->SetTextAlign(12); t->SetTextAngle(90); t->SetTextFont(42); t->SetTextSize(textSize);
    t->SetTextColor(kGray+2); t->Draw();
  }
}

// For drawing only: bins whose nominal purity is empty (window not fitted) are
// pushed far below any frame, so they are absent rather than drawn at zero.
inline void blankUnfitted(TH1D *h, TH1D *nominal)
{
  for(int i = 1; i <= h->GetNbinsX(); i++)
    if(nominal->GetBinContent(i) <= 0.){ h->SetBinContent(i, -999.); h->SetBinError(i, 0.); }
}

struct Tpl {
  TH1D *bFC, *bGS, *lc;     // unit-normalised over 0-5 GeV
  TH1D *c = nullptr, *l = nullptr;   // c and light (u+d+s+g) alone, unit-normalised (b/c/light fit)
  double fGSnat;            // PYTHIA's bGS share of all b in this window
  double truthB;            // MC-truth b fraction (b + bGS) over 0-5 GeV
};

// b template with the bGS share set to s (share of all b)
TH1D* bWithGS(const Tpl &t, double s)
{
  TH1D *b = (TH1D*) t.bFC->Clone(Form("bmix%d", uidP++)); b->SetDirectory(nullptr);
  b->Scale(1. - s); b->Add(t.bGS, s);
  return b;
}

struct FitRes { double fb, efb, fGS, efGS; double chi2; int ndf; TH1D *compB, *compBGS, *compL, *total; TH1D *compFake = nullptr; };

// FAKE-ptRel AS A FIXED FIT COMPONENT (2026-10-08, user). When set (PbPb PF:
// dataWin in caloJetBPurityFit_PbPb.h sets it to the normalised T2 + T3 counts of
// the window, nullptr otherwise), doFit fits the RAW data with
//   phi x fake + (1 - phi) x [f x b + (1 - f) x light+c],
// phi FIXED to the fake counts / data counts inside the fit range. This replaces
// the bin-by-bin subtraction, which left empty data bins negative with only the
// template's error (e.g. -0.6 +- 0.2) and broke sparse fits. fb is returned as
// the b fraction of ALL in-range data, (1 - phi) f, so purityFull() is b / all
// tagged jets directly (fakeKeepFrac() is 1 in this mode).
TH1D *activeFakeTpl = nullptr;

// Fit data with N templates (2 or 3) over [lo, hi] (default the nominal range).
// The fractions are fractions of the jets INSIDE the fit range, as in
// templateFitter(). Components come back as counts per 0.1 GeV bin, normalised
// so that inside the fit range they add up to the data.
FitRes doFit(TH1D *data05, std::vector<TH1D*> tpl, double lo = fitLo, double hi = fitHi)
{
  FitRes r{}; const int N = (int) tpl.size();
  RooRealVar x(Form("x%d", uidP), "p_{T}^{rel}", 0., fitHi); x.setBins(nFit05());
  x.setRange("fitRange", lo, hi);
  RooDataHist dh(Form("dh%d", uidP), "", x, RooFit::Import(*data05));
  std::vector<RooDataHist*> th; std::vector<RooHistPdf*> tp;
  for(int i = 0; i < N; i++){
    th.push_back(new RooDataHist(Form("th%d_%d", uidP, i), "", x, RooFit::Import(*tpl[i])));
    tp.push_back(new RooHistPdf(Form("tp%d_%d", uidP, i), "", x, *th.back(), 0));
  }
  const int b1 = data05->FindBin(lo + 1e-6), b2 = data05->FindBin(hi - 1e-6);
  const double nD = data05->Integral(b1, b2);
  // fixed fake component (activeFakeTpl, counts on the data's binning)
  const bool withFake = activeFakeTpl && activeFakeTpl->Integral(b1, b2) > 0.;
  const double phiVal = withFake && nD > 0. ? TMath::Min(0.95, activeFakeTpl->Integral(b1, b2)/nD) : 0.;
  RooDataHist *fh = nullptr; RooHistPdf *fp = nullptr;
  if(withFake){
    TH1D *fpos = (TH1D*) activeFakeTpl->Clone(Form("fpos%d", uidP)); fpos->SetDirectory(nullptr);
    for(int i = 0; i <= fpos->GetNbinsX()+1; i++) if(fpos->GetBinContent(i) < 0.) fpos->SetBinContent(i, 0.);   // a pdf cannot be negative
    fh = new RooDataHist(Form("fh%d", uidP), "", x, RooFit::Import(*fpos)); fp = new RooHistPdf(Form("fp%d", uidP), "", x, *fh, 0);
    delete fpos;
  }
  RooRealVar phi(Form("phi%d", uidP), "", phiVal); phi.setConstant(true);
  RooRealVar f1(Form("f1_%d", uidP), "", 0.5, 0., 1.), f2(Form("f2_%d", uidP), "", 0.1, 0., 1.);
  RooArgList pdfs, coefs;
  if(withFake){ pdfs.add(*fp); coefs.add(phi); }
  for(auto p : tp) pdfs.add(*p);
  coefs.add(f1); if(N == 3) coefs.add(f2);
  // recursive fractions with the fake first: fake phi, then f1 (and f2) of the rest
  RooAddPdf model(Form("model%d", uidP++), "", pdfs, coefs, withFake);
  // Without this, a Range() fit returns the coefficients as fractions over the
  // FULL domain (ROOT 6.40; checked against a fit to templates truncated
  // at hi), and purityFull() and the components below would convert them a
  // second time -- purity biased high by 1-4% at the nominal range, far more
  // for short ranges.
  model.fixCoefRange("fitRange");
  model.fitTo(dh, RooFit::SumW2Error(true), RooFit::Range("fitRange"), RooFit::PrintLevel(-1));

  // fractions of ALL in-range data
  const double keep = 1. - phiVal;
  r.fb = keep*f1.getVal(); r.efb = keep*f1.getError();
  if(N == 3){
    // recursive: f2 is a fraction of what f1 leaves
    if(withFake){ r.fGS = keep*(1. - f1.getVal())*f2.getVal(); r.efGS = keep*(1. - f1.getVal())*f2.getError(); }
    else { r.fGS = f2.getVal(); r.efGS = f2.getError(); }
  }

  // components on the data's normalisation inside the fit range
  double w[3] = {r.fb, N == 3 ? r.fGS : 0., 0.};
  w[N-1] = keep - r.fb - (N == 3 ? r.fGS : 0.);
  r.compB = (TH1D*) tpl[0]->Clone(Form("cb%d", uidP++)); r.compB->SetDirectory(nullptr); r.compB->Scale(w[0]*nD/tpl[0]->Integral(b1, b2));
  if(N == 3){ r.compBGS = (TH1D*) tpl[1]->Clone(Form("cg%d", uidP++)); r.compBGS->SetDirectory(nullptr); r.compBGS->Scale(w[1]*nD/tpl[1]->Integral(b1, b2)); }
  r.compL = (TH1D*) tpl[N-1]->Clone(Form("cl%d", uidP++)); r.compL->SetDirectory(nullptr); r.compL->Scale(w[N-1]*nD/tpl[N-1]->Integral(b1, b2));
  r.total = (TH1D*) r.compB->Clone(Form("ct%d", uidP++)); r.total->SetDirectory(nullptr);
  if(N == 3) r.total->Add(r.compBGS);
  r.total->Add(r.compL);
  if(withFake){ r.compFake = (TH1D*) activeFakeTpl->Clone(Form("cf%d", uidP++)); r.compFake->SetDirectory(nullptr); r.compFake->Scale(phiVal*nD/activeFakeTpl->Integral(b1, b2)); r.total->Add(r.compFake); }

  // chi2 for display, inside the fit range: data errors plus the templates' MC
  // statistical errors. An EMPTY data bin has no error of its own; it takes the
  // model expectation as its variance (Pearson), or a sparse window -- PbPb
  // peripheral, ~150 jets over 50 bins -- shows chi2 in the thousands. Display
  // only: the fit itself is a likelihood and never reads this.
  r.chi2 = 0.; r.ndf = 0;
  for(int i = b1; i <= b2; i++){
    double d = data05->GetBinContent(i), ed = data05->GetBinError(i), m = r.total->GetBinContent(i);
    double e2 = (d > 0. || ed > 0.) ? ed*ed : TMath::Max(m, 0.);
    for(int k = 0; k < N; k++) e2 += pow(w[k]*nD/tpl[k]->Integral(b1, b2)*tpl[k]->GetBinError(i), 2);
    if(withFake) e2 += pow(phiVal*nD/activeFakeTpl->Integral(b1, b2)*activeFakeTpl->GetBinError(i), 2);
    if(e2 > 0.){ r.chi2 += (d - m)*(d - m)/e2; r.ndf++; }
  }
  r.ndf -= (N - 1);

  for(auto p : tp) delete p;
  for(auto h : th) delete h;
  delete fp; delete fh;
  return r;
}

// 3-TEMPLATE FIT WITH A LIMITED g->bb SHARE (2026-10-09, user: the free 3-template
// fit must not be allowed to call every b g->bb just because that fits). The g->bb
// share s of all b (over 0-fitHi, as Tpl::fGSnat) is scanned on [nat x gsShareLimitLo,
// min(1, nat x gsShareLimitHi)] = [1.0, 1.3] x PYTHIA (user, 2026-10-09: never below
// PYTHIA's share, at most 30% above; first version: a factor 2 either way, which
// 12 of 20 PbPb windows ran into); at each s the b template is bWithGS(t, s)
// and doFit() fits b / light+c; the best s minimizes the fit chi2 (parabola through
// the grid minimum; at a limit it stays there and atLimit is set). Returned in the
// 3-template convention of doFit(data05, {bFC, bGS, lc}): fb / fGS = in-range
// fractions of flavor creation / g->bb, so purityFull(r, {bFC, bGS, lc}, 2, ...) and
// the plotting code work unchanged. Errors: the b-fraction fit error at the best s
// (the share's own uncertainty is in info.eShare, not folded in).
// gsShareLimitHi <= 0: plain free fit (doFit with the three templates).
double gsShareLimitLo = 1.0, gsShareLimitHi = 1.3;
struct GSLimitInfo { double share = 0., eShare = 0., lo = 0., hi = 1.; bool atLimit = false; };
inline const char* gsTplTag(){ return gsShareLimitHi > 0. ? "3tpl_limitGS" : "3tpl_freeGS"; }
inline TString gsTplLabel(){ return gsShareLimitHi > 0. ? TString(Form("3 templates, g#rightarrowb#bar{b} share %.3g-%.3g #times PYTHIA", gsShareLimitLo, gsShareLimitHi)) : TString("3 templates, g#rightarrowb#bar{b} share free"); }
// the fit r (fractions of all in-range data) redrawn on another binning: the
// components scaled as in doFit / doFitCL from the same templates and data
// (and activeFakeTpl) rebuilt with range05BinOverride set. The fit's chi2 / ndf
// are carried over; the fit itself is untouched. Display only.
inline FitRes componentsAt(const FitRes &r, TH1D *data, std::vector<TH1D*> tpl, double lo, double hi)
{
  FitRes o = r; o.compB = o.compBGS = o.compL = o.total = o.compFake = nullptr;
  const int N = (int) tpl.size();
  const int b1 = data->FindBin(lo + 1e-6), b2 = data->FindBin(hi - 1e-6);
  const double nD = data->Integral(b1, b2);
  const bool withFake = activeFakeTpl && activeFakeTpl->Integral(b1, b2) > 0.;
  const double phiVal = withFake && nD > 0. ? TMath::Min(0.95, activeFakeTpl->Integral(b1, b2)/nD) : 0.;
  double w[3] = {r.fb, N == 3 ? r.fGS : 0., 0.};
  w[N-1] = 1. - phiVal - r.fb - (N == 3 ? r.fGS : 0.);
  o.compB = (TH1D*) tpl[0]->Clone(Form("fb%d", uidP++)); o.compB->SetDirectory(nullptr); o.compB->Scale(w[0]*nD/tpl[0]->Integral(b1, b2));
  if(N == 3){ o.compBGS = (TH1D*) tpl[1]->Clone(Form("fg%d", uidP++)); o.compBGS->SetDirectory(nullptr); o.compBGS->Scale(w[1]*nD/tpl[1]->Integral(b1, b2)); }
  o.compL = (TH1D*) tpl[N-1]->Clone(Form("fl%d", uidP++)); o.compL->SetDirectory(nullptr); o.compL->Scale(w[N-1]*nD/tpl[N-1]->Integral(b1, b2));
  o.total = (TH1D*) o.compB->Clone(Form("ft%d", uidP++)); o.total->SetDirectory(nullptr);
  if(N == 3) o.total->Add(o.compBGS);
  o.total->Add(o.compL);
  if(withFake){ o.compFake = (TH1D*) activeFakeTpl->Clone(Form("ff%d", uidP++)); o.compFake->SetDirectory(nullptr); o.compFake->Scale(phiVal*nD/activeFakeTpl->Integral(b1, b2)); o.total->Add(o.compFake); }
  return o;
}
inline void deleteComps(FitRes &r){ delete r.compB; delete r.compBGS; delete r.compL; delete r.total; delete r.compFake; r.compB = r.compBGS = r.compL = r.total = r.compFake = nullptr; }
FitRes doFitGSLimited(TH1D *data05, const Tpl &t, double lo = fitLo, double hi = fitHi, GSLimitInfo *info = nullptr)
{
  if(gsShareLimitHi <= 0.) return doFit(data05, {t.bFC, t.bGS, t.lc}, lo, hi);
  const double nat = t.fGSnat, sLo = nat*gsShareLimitLo, sHi = TMath::Min(1., nat*gsShareLimitHi);
  const int nS = 31; std::vector<double> sv(nS), cv(nS); int kb = 0;
  for(int k = 0; k < nS; k++){
    sv[k] = sLo + (sHi - sLo)*k/(nS - 1);
    TH1D *b = bWithGS(t, sv[k]); FitRes r = doFit(data05, {b, t.lc}, lo, hi); cv[k] = r.chi2; deleteComps(r); delete b;
    if(cv[k] < cv[kb]) kb = k;
  }
  double sB = sv[kb], eS = 0.; bool atLim = (kb == 0 || kb == nS - 1);
  if(!atLim){
    const double h = sv[1] - sv[0], A = (cv[kb+1] - 2*cv[kb] + cv[kb-1])/(2*h*h);
    if(A > 0.){ sB = sv[kb] - (cv[kb+1] - cv[kb-1])/(4*A*h); eS = 1./sqrt(A); }
  }
  TH1D *b = bWithGS(t, sB); FitRes r = doFit(data05, {b, t.lc}, lo, hi); delete b;
  // convert to flavor-creation / g->bb in-range fractions
  const int b1 = data05->FindBin(lo + 1e-6), b2 = data05->FindBin(hi - 1e-6);
  const double IF = t.bFC->Integral(b1, b2)/t.bFC->Integral(1, t.bFC->GetNbinsX()), IG = t.bGS->Integral(b1, b2)/t.bGS->Integral(1, t.bGS->GetNbinsX());
  const double Im = (1. - sB)*IF + sB*IG, wF = Im > 0. ? (1. - sB)*IF/Im : 0., wG = Im > 0. ? sB*IG/Im : 0.;
  FitRes o = r; o.fb = r.fb*wF; o.efb = r.efb*wF; o.fGS = r.fb*wG; o.efGS = r.efb*wG; o.ndf = r.ndf - (atLim ? 0 : 1);
  const double nD = data05->Integral(b1, b2);
  o.compB = (TH1D*) t.bFC->Clone(Form("cbL%d", uidP++)); o.compB->SetDirectory(nullptr); o.compB->Scale(o.fb*nD/t.bFC->Integral(b1, b2));
  o.compBGS = (TH1D*) t.bGS->Clone(Form("cgL%d", uidP++)); o.compBGS->SetDirectory(nullptr); o.compBGS->Scale(o.fGS*nD/t.bGS->Integral(b1, b2));
  delete r.compB;   // the mixed-b component, replaced by its two parts; compL, total, compFake kept
  if(info){ info->share = sB; info->eShare = eS; info->lo = sLo; info->hi = sHi; info->atLimit = atLim; }
  return o;
}

// UPPER FIT EDGE, pp and PbPb alike. Nominal since 2026-10-09: a fixed 3.0 GeV
// (kEdgeFixed below). Before: where the light+c template dies out -- the bin
// edge at which its cumulative over 0-5 GeV first reaches q = 0.98 (kEdgeLC98).
// Above that only the b template has support, so any data there -- uncorrelated
// muon-jet pairs, a mismodelled b tail -- could only be absorbed by the b
// fraction. The fit-range systematic moves the upper edge -/+ 0.5 GeV; the lower
// edge stays at 0. Take the edge from the NOMINAL templates and
// keep it fixed for every variation, and report purityFull() so the purity
// always refers to all jets in 0-5 GeV.
const double lcTailQuantileNominal = 0.98;
// Upper-edge rule switch (2026-10-08 study): kEdgeLC98 = the light+c 0.98 quantile
// (nominal, above); kEdgeB99 = where the b template (flavor creation + g->bb at
// the natural share) reaches 0.99 of its integral over the domain.
// kEdgeFixed (NOMINAL since 2026-10-09, user) = the same edge, fitEdgeFixed = 3.0 GeV,
// in every window, pp and PbPb: studyFitEdgeScan_PF.C showed the central 80-100 GeV
// chi2/ndf rising from 1.2 at 3 GeV to 2.9 at the light+c 0.98 edge (3.5) and 6.0 at
// 5 GeV, and above 100 GeV no edge is preferred. The fit-range systematic still
// moves it by -/+ 0.5 GeV.
enum FitEdgeRule { kEdgeLC98 = 0, kEdgeB99 = 1, kEdgeFixed = 2 };
FitEdgeRule fitEdgeRule = kEdgeFixed;
double fitEdgeFixed = purityVariantOn() ? 5.0 : 3.0;   // 5.0 in the purity variant (above)
double bTailQuantile = 0.99;   // kEdgeB99 quantile (0.95 also studied, 2026-10-08)
inline double lightTailEdge(const Tpl &t, double q = lcTailQuantileNominal)
{
  if(fitEdgeRule == kEdgeFixed) return TMath::Min(fitHi, fitEdgeFixed);
  if(fitEdgeRule == kEdgeB99){
    TH1D *b = (TH1D*) t.bFC->Clone("bEdgeTmp"); b->SetDirectory(nullptr); b->Scale(1. - t.fGSnat); b->Add(t.bGS, t.fGSnat);
    const int n = b->GetNbinsX(); const double tot = b->Integral(1, n); double s = 0., e = fitHi;
    for(int i = 1; i <= n; i++){ s += b->GetBinContent(i); if(s >= bTailQuantile*tot){ e = b->GetXaxis()->GetBinUpEdge(i); break; } }
    delete b; return e;
  }
  const int n = t.lc->GetNbinsX();
  const double tot = t.lc->Integral(1, n); double s = 0.;
  for(int i = 1; i <= n; i++){ s += t.lc->GetBinContent(i); if(s >= q*tot) return t.lc->GetXaxis()->GetBinUpEdge(i); }
  return fitHi;
}

// b purity over the full 0-5 GeV ptRel range from a fit over [lo, hi]. The
// fitted fractions are fractions of the jets INSIDE [lo, hi]; each of the first
// nB components (the b templates) is turned into a count over 0-5 GeV by
// dividing by its template's share inside [lo, hi] (templates unit-area over
// 0-5), and the sum is divided by all data in 0-5 GeV. With [lo, hi] = [0, 5]
// this is the fitted b fraction itself. err (optional) scales the fit error the
// same way.
inline double purityFull(const FitRes &r, const std::vector<TH1D*> &tpl, int nB, TH1D *data05,
                         double lo, double hi, double *err = nullptr)
{
  const int b1 = data05->FindBin(lo + 1e-6), b2 = data05->FindBin(hi - 1e-6);
  const double nIn = data05->Integral(b1, b2), nAll = data05->Integral(1, data05->GetNbinsX());
  if(nAll <= 0.){ if(err) *err = 0.; return 0.; }
  const double f[2] = {r.fb, r.fGS}, ef[2] = {r.efb, r.efGS};
  double p = 0., e2 = 0.;
  for(int k = 0; k < nB && k < 2; k++){
    const double share = tpl[k]->Integral(b1, b2)/tpl[k]->Integral(1, tpl[k]->GetNbinsX());
    if(share <= 0.) continue;
    p += f[k]*nIn/share/nAll; e2 += pow(ef[k]*nIn/share/nAll, 2);
  }
  if(err) *err = sqrt(e2);
  return p;
}

// Nominal fit model: 2 templates b / light+c, c at its MC share (cMult varies
// it). fitCLSeparatePF = true makes PFJets fit 3 templates b / c / light with the
// c fraction free (cMult then has no effect); the upper fit edge stays
// lightTailEdge() of light+c either way. OFF: tried 2026-10-07 on pp PF data
// (0.5 GeV bins) and ill-posed -- c and light are nearly degenerate in ptRel, the
// fitted c share of non-b runs to 1.0 or 0 (MC: 0.61-0.71), corr(b, c) reaches
// -1.00, and b at 200-300 GeV came out 0.70 +- 3.0. Not const so a study can turn
// it on (plotBPurity_caloJets_pp.C, clSplit = true).
bool fitCLSeparatePF = false;
inline bool fitCLSeparate(){ return jetColl == kPFJets && fitCLSeparatePF; }
inline std::vector<TH1D*> nominalTpls(TH1D *b, const Tpl &t)
{
  if(fitCLSeparate() && t.c && t.l) return {b, t.c, t.l};
  return {b, t.lc};
}
// b / c / light fit with recursive fractions: fb = b fraction, rc = c share of
// the non-b jets, both in [0, 1] (so the fractions can never sum past 1), rc
// started at its MC value cShareMC. Returned as doFit's 3-template FitRes: fb,
// fGS = c fraction of all = (1 - fb) rc (error propagated with the fit
// covariance), components b, c (compBGS), light. Also reports rc and the
// correlation of (fb, rc).
// cMultSigma > 0 constrains the c weight instead (cMult convention, as the
// c-fraction systematic): a multiplier m on PYTHIA's c weight with a Gaussian
// constraint of mean 1, width cMultSigma, so rc = m r / (m r + 1 - r), r = cShareMC.
struct CLExtra { double rc = 0., erc = 0., corr = 0., m = 1., em = 0.; int status = -1; };
inline FitRes doFitCL(TH1D *data05, TH1D *b, TH1D *c, TH1D *l, double cShareMC, double lo, double hi, CLExtra *ex = nullptr,
                      double cMultSigma = -1.)
{
  FitRes r{}; const int id = uidP++;
  RooRealVar x(Form("xcl%d", id), "p_{T}^{rel}", 0., fitHi); x.setBins(nFit05());
  x.setRange("fitRange", lo, hi);
  RooDataHist dh(Form("dhcl%d", id), "", x, RooFit::Import(*data05));
  RooDataHist hb(Form("hbcl%d", id), "", x, RooFit::Import(*b)), hc(Form("hccl%d", id), "", x, RooFit::Import(*c)), hl(Form("hlcl%d", id), "", x, RooFit::Import(*l));
  RooHistPdf pb(Form("pbcl%d", id), "", x, hb), pc(Form("pccl%d", id), "", x, hc), pl(Form("plcl%d", id), "", x, hl);
  const bool constr = cMultSigma > 0.;
  const double r0 = TMath::Min(TMath::Max(cShareMC, 0.01), 0.99);
  RooRealVar fb(Form("fbcl%d", id), "", 0.5, 0., 1.), rcFree(Form("rccl%d", id), "", r0, 0., 1.);
  RooRealVar mC(Form("mcl%d", id), "", 1., 0., 1. + 6.*TMath::Max(cMultSigma, 0.));
  RooFormulaVar rcM(Form("rcmcl%d", id), Form("@0*%.6f/(@0*%.6f + %.6f)", r0, r0, 1. - r0), RooArgList(mC));
  RooRealVar m0(Form("m0cl%d", id), "", 1.), ms(Form("mscl%d", id), "", TMath::Max(cMultSigma, 1e-3));
  RooGaussian gC(Form("gcl%d", id), "", mC, m0, ms);
  RooAbsReal &rc = constr ? (RooAbsReal&) rcM : (RooAbsReal&) rcFree;
  RooRealVar &pc2 = constr ? mC : rcFree;   // the floating c parameter
  RooAddPdf nonb(Form("nonbcl%d", id), "", RooArgList(pc, pl), RooArgList(rc));
  // fixed fake-ptRel component (activeFakeTpl), as in doFit (2026-10-08)
  const int fb1 = data05->FindBin(lo + 1e-6), fb2 = data05->FindBin(hi - 1e-6);
  const double nDin = data05->Integral(fb1, fb2);
  const bool withFake = activeFakeTpl && activeFakeTpl->Integral(fb1, fb2) > 0.;
  const double phiVal = withFake && nDin > 0. ? TMath::Min(0.95, activeFakeTpl->Integral(fb1, fb2)/nDin) : 0.;
  RooDataHist *fh = nullptr; RooHistPdf *fp = nullptr;
  if(withFake){
    TH1D *fpos = (TH1D*) activeFakeTpl->Clone(Form("fposcl%d", id)); fpos->SetDirectory(nullptr);
    for(int i = 0; i <= fpos->GetNbinsX()+1; i++) if(fpos->GetBinContent(i) < 0.) fpos->SetBinContent(i, 0.);
    fh = new RooDataHist(Form("fhcl%d", id), "", x, RooFit::Import(*fpos)); fp = new RooHistPdf(Form("fpcl%d", id), "", x, *fh, 0); delete fpos;
  }
  RooRealVar phi(Form("phicl%d", id), "", phiVal); phi.setConstant(true);
  RooArgList mPdfs, mCoefs;
  if(withFake){ mPdfs.add(*fp); mCoefs.add(phi); }
  mPdfs.add(pb); mPdfs.add(nonb); mCoefs.add(fb);
  RooAddPdf model(Form("modelcl%d", id), "", mPdfs, mCoefs, withFake);
  nonb.fixCoefRange("fitRange"); model.fixCoefRange("fitRange");   // fractions inside the fit range, as in doFit
  std::unique_ptr<RooFitResult> res(constr
    ? model.fitTo(dh, RooFit::SumW2Error(true), RooFit::Range("fitRange"), RooFit::PrintLevel(-1), RooFit::Strategy(2), RooFit::Save(), RooFit::ExternalConstraints(RooArgSet(gC)))
    : model.fitTo(dh, RooFit::SumW2Error(true), RooFit::Range("fitRange"), RooFit::PrintLevel(-1), RooFit::Strategy(2), RooFit::Save()));
  // rc and its error; for the constrained form propagated from m
  const double rcV = rc.getVal();
  const double dRdM = constr ? r0*(1. - r0)/pow(mC.getVal()*r0 + 1. - r0, 2) : 1.;
  const double erc = constr ? dRdM*mC.getError() : rcFree.getError();
  const double corr = res ? res->correlation(fb, pc2) : 0.;
  const double cov = corr*fb.getError()*erc;
  const double keep = 1. - phiVal;   // fractions of ALL in-range data
  r.fb = keep*fb.getVal(); r.efb = keep*fb.getError();
  r.fGS = keep*(1. - fb.getVal())*rcV;
  r.efGS = keep*sqrt(TMath::Max(0., pow(rcV*fb.getError(), 2) + pow((1. - fb.getVal())*erc, 2) - 2.*rcV*(1. - fb.getVal())*cov));
  if(ex){ ex->rc = rcV; ex->erc = erc; ex->corr = corr; ex->status = res ? res->status() : -1; ex->m = constr ? mC.getVal() : 1.; ex->em = constr ? mC.getError() : 0.; }

  // components and display chi2, as in doFit
  const int b1 = data05->FindBin(lo + 1e-6), b2 = data05->FindBin(hi - 1e-6);
  const double nD = data05->Integral(b1, b2);
  const double w[3] = {r.fb, r.fGS, keep - r.fb - r.fGS};
  TH1D *tp[3] = {b, c, l}, *cp[3];
  for(int k = 0; k < 3; k++){ cp[k] = (TH1D*) tp[k]->Clone(Form("ccl%d_%d", id, k)); cp[k]->SetDirectory(nullptr); cp[k]->Scale(w[k]*nD/tp[k]->Integral(b1, b2)); }
  r.compB = cp[0]; r.compBGS = cp[1]; r.compL = cp[2];
  r.total = (TH1D*) cp[0]->Clone(Form("ctcl%d", id)); r.total->SetDirectory(nullptr); r.total->Add(cp[1]); r.total->Add(cp[2]);
  if(withFake){ r.compFake = (TH1D*) activeFakeTpl->Clone(Form("cfcl%d", id)); r.compFake->SetDirectory(nullptr); r.compFake->Scale(phiVal*nD/activeFakeTpl->Integral(b1, b2)); r.total->Add(r.compFake); }
  r.chi2 = 0.; r.ndf = 0;
  for(int i = b1; i <= b2; i++){
    double d = data05->GetBinContent(i), ed = data05->GetBinError(i), m = r.total->GetBinContent(i);
    double e2 = (d > 0. || ed > 0.) ? ed*ed : TMath::Max(m, 0.);
    for(int k = 0; k < 3; k++) e2 += pow(w[k]*nD/tp[k]->Integral(b1, b2)*tp[k]->GetBinError(i), 2);
    if(withFake) e2 += pow(phiVal*nD/activeFakeTpl->Integral(b1, b2)*activeFakeTpl->GetBinError(i), 2);
    if(e2 > 0.){ r.chi2 += (d - m)*(d - m)/e2; r.ndf++; }
  }
  r.ndf -= 2;
  delete fp; delete fh;
  return r;
}

// ADAPTIVE 2 / 3 TEMPLATES (2026-10-08 study, user): whether c and light are
// fitted as one template (b vs light+c, 1 free parameter) or separately (b / c /
// light, 2 free parameters) is decided per window by how distinct their shapes
// are. clShapeTest compares the unit-normalised c and light templates inside the
// fit range, on the fit binning, with both MC errors: chi2, ndf, p-value. It also
// gives the DATA's power to tell them apart: the expected separation in sigma
// between "all non-b jets are c" and "all are light" for nNonB jets,
// sqrt(nNonB x sum (c_i - l_i)^2 / ((c_i + l_i)/2)).
// clShapeTest also asks whether c is distinct from every MIXTURE of b and light
// (b is free in the fit, so c ~ a b + (1-a) light can be traded for b and light
// at no cost): c is fitted as a b + (1 - a) light, a in [0, 1], over the same
// bins with all three MC errors; chi2Mix / ndfMix / pMix and the best a. A small
// pMix means c is NOT such a mixture, i.e. a 3rd template adds information.
struct CLShape { double chi2 = 0., p = 1., sepData = 0.; int ndf = 0; double chi2Mix = 0., pMix = 1., aMix = 0.; int ndfMix = 0; };
inline CLShape clShapeTest(const Tpl &t, double lo, double hi, double nNonB, const TH1D *bT = nullptr)
{
  CLShape r;
  if(!t.c || !t.l) return r;
  const int b1 = t.c->FindBin(lo + 1e-6), b2 = t.c->FindBin(hi - 1e-6);
  const double ic = t.c->Integral(b1, b2), il = t.l->Integral(b1, b2);
  if(ic <= 0. || il <= 0.) return r;
  double s2 = 0.;
  for(int i = b1; i <= b2; i++){
    const double c = t.c->GetBinContent(i)/ic, l = t.l->GetBinContent(i)/il;
    const double e2 = pow(t.c->GetBinError(i)/ic, 2) + pow(t.l->GetBinError(i)/il, 2);
    if(e2 > 0.){ r.chi2 += (c - l)*(c - l)/e2; r.ndf++; }
    if(c + l > 0.) s2 += (c - l)*(c - l)/(0.5*(c + l));
  }
  r.ndf -= 1;   // both normalised
  r.p = r.ndf > 0 ? TMath::Prob(r.chi2, r.ndf) : 1.;
  r.sepData = sqrt(TMath::Max(0., nNonB*s2));
  if(bT){
    const double ib = bT->Integral(b1, b2);
    if(ib > 0.){
      auto chiAt = [&](double a, int &n){ double c2 = 0.; n = 0;
        for(int i = b1; i <= b2; i++){
          const double c = t.c->GetBinContent(i)/ic, l = t.l->GetBinContent(i)/il, b = bT->GetBinContent(i)/ib;
          const double e2 = pow(t.c->GetBinError(i)/ic, 2) + pow(a*bT->GetBinError(i)/ib, 2) + pow((1. - a)*t.l->GetBinError(i)/il, 2);
          if(e2 > 0.){ c2 += pow(c - a*b - (1. - a)*l, 2)/e2; n++; } }
        return c2; };
      double best = 1e300; int n = 0;
      for(double a = 0.; a <= 1. + 1e-9; a += 0.002){ int nn; double c2 = chiAt(a, nn); if(c2 < best){ best = c2; r.aMix = a; n = nn; } }
      r.chi2Mix = best; r.ndfMix = n - 2;   // normalisation + a
      r.pMix = r.ndfMix > 0 ? TMath::Prob(best, r.ndfMix) : 1.;
    }
  }
  return r;
}
// MC c share of the non-b jets in light+c: least-squares lc = l + r (c - l)
inline double cShareOfLC(const Tpl &t)
{
  double sxx = 0., sxy = 0.;
  for(int k = 1; k <= t.lc->GetNbinsX(); k++){
    double xx = t.c->GetBinContent(k) - t.l->GetBinContent(k), yy = t.lc->GetBinContent(k) - t.l->GetBinContent(k);
    sxx += xx*xx; sxy += xx*yy;
  }
  return sxx > 0. ? sxy/sxx : 0.;
}
inline const char* nomModelLabel(){ return fitCLSeparate() ? "3 templates b / c / light" : "2 templates"; }
// nominal-model fit with b template b; tplOut (optional) gets the templates used
inline FitRes fitNominal(TH1D *data05, TH1D *b, const Tpl &t, double lo, double hi, std::vector<TH1D*> *tplOut = nullptr)
{
  std::vector<TH1D*> v = nominalTpls(b, t);
  if(tplOut) *tplOut = v;
  if(v.size() == 3) return doFitCL(data05, v[0], v[1], v[2], cShareOfLC(t), lo, hi);
  return doFit(data05, v, lo, hi);
}
// nominal-model b purity over 0-5 GeV (only the first template is b)
inline double purityNominal(TH1D *data05, TH1D *b, const Tpl &t, double lo, double hi, double *err = nullptr)
{
  std::vector<TH1D*> v;
  FitRes r = fitNominal(data05, b, t, lo, hi, &v);
  return purityFull(r, v, 1, data05, lo, hi, err);
}

// 2-template (b, light+c) purity over 0-5 GeV from a fit over [lo, hi]
inline double purity2(TH1D *data05, TH1D *b, TH1D *lc, double lo, double hi, double *err = nullptr)
{
  FitRes r = doFit(data05, {b, lc}, lo, hi);
  return purityFull(r, {b, lc}, 1, data05, lo, hi, err);
}

// TEMPLATE MC STATISTICS (2026-10-08, user: a systematic source, pp and PbPb).
// The fit (RooHistPdf) treats the templates as exact, so its error is the
// data's alone. Here each toy fluctuates every template bin -- b (flavor
// creation), g->bb, light+c and, for the 3-template fit, c and light -- by a
// Gaussian of its MC error (templates keep their unsmeared relative errors,
// see axisSmearToy), clips at 0, renormalises, and refits with the nominal
// model, range and g->bb shift. Returns RMS / nominal of the purity. Fixed seed
// per call site, so a rerun reproduces the source.
#include "TRandom3.h"
const int tplStatNToys = 200;
inline TH1D* tplStatFluct(const TH1D *h, TRandom3 &r)
{
  if(!h) return nullptr;
  TH1D *o = (TH1D*) h->Clone(Form("tsf%d", uidP++)); o->SetDirectory(nullptr);
  for(int i = 1; i <= o->GetNbinsX(); i++) o->SetBinContent(i, TMath::Max(0., h->GetBinContent(i) + r.Gaus(0., h->GetBinError(i))));
  if(o->Integral() > 0.) o->Scale(h->Integral()/o->Integral());
  return o;
}
inline unsigned tplStatSeed(int cls, double lo){ return 7919u*(unsigned)(cls + 2) + 104729u*(unsigned)(lo + 0.5); }
inline double templateStatRel(TH1D *data05, const Tpl &t, double gsShift, double lo, double hi, unsigned seed, int nToy = tplStatNToys)
{
  TH1D *b0 = bWithGS(t, TMath::Min(1., t.fGSnat + gsShift));
  const double p0 = purityNominal(data05, b0, t, lo, hi);
  delete b0;
  if(p0 <= 0.) return 0.;
  TRandom3 rng(seed);
  double s1 = 0., s2 = 0.;
  for(int k = 0; k < nToy; k++){
    Tpl tt = t;
    tt.bFC = tplStatFluct(t.bFC, rng); tt.bGS = tplStatFluct(t.bGS, rng); tt.lc = tplStatFluct(t.lc, rng);
    tt.c = tplStatFluct(t.c, rng); tt.l = tplStatFluct(t.l, rng);
    TH1D *b = bWithGS(tt, TMath::Min(1., tt.fGSnat + gsShift));
    const double p = purityNominal(data05, b, tt, lo, hi);
    s1 += p; s2 += p*p;
    delete b; delete tt.bFC; delete tt.bGS; delete tt.lc; delete tt.c; delete tt.l;
  }
  const double mu = s1/nToy;
  return sqrt(TMath::Max(0., s2/nToy - mu*mu))/p0;
}
