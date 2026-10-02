// Muon reco + tight-ID efficiency from Z -> mu mu tag-and-probe, pp or PbPb data.
//
// Reads the tnp TREE written by src/scanning/pp/pp_muonTagAndProbe_scan.C, not
// its histograms: the analysis cuts muons at |eta| < 2.0 (scan_muon_tag.h), and
// 2.0 is not an edge of the scan's |eta| axis. From the tree any binning works
// without a rescan.
//
// STEPS (see the scan header for the definitions)
//   trk  STA probe          -> has an inner track
//   id   inner-track probe  -> analysis tight ID
//   all  any muon object    -> analysis tight ID
// trk x id should reproduce all; the two are printed side by side.
//
// BACKGROUND. Signal = opposite-sign window count minus same-sign window count,
// separately for passing and failing probes. Windows:
//   id, all   75-105
//   trk       71-111
// trk is wider because a failing STA probe's mass uses its standalone momentum,
// whose resolution is poor; a narrow window would lose failing probes faster than
// passing ones and bias the efficiency up.
//
// Same-sign is nominal because it assumes nothing about the background's mass
// shape. The opposite-sign sideband method it replaced assumed a flat background
// across window and sidebands, which fails below ~25 GeV probe pT: there the
// continuum falls roughly exponentially and the sideband average overestimates
// what sits under the peak. Same-sign models the charge-uncorrelated part (fakes,
// random pairings). It misses charge-correlated opposite-sign background (b bbar
// -> mu+ mu-), which is what the sideband variation below is there to bound.
//
// SYSTEMATIC. Largest shift of two variations:
//   sideband   opposite-sign sidebands, scaled to the window width, assuming a
//              flat background:  id/all 61-71 and 111-121,  trk 51-61 and 121-131
//   none       no subtraction at all
// The sideband variation also removes the Z's radiative tail and Drell-Yan as
// "background" although they are real muons, so it overstates the background
// and the systematic is conservative.
//
// STATISTICAL ERRORS. Clopper-Pearson 68% on the raw window counts, placed around
// the subtracted central value. Efficiencies here sit near 1 with few failures,
// where a Gaussian propagation is wrong; with background this small the raw
// counts carry essentially all of the statistical information.
//
// DUPLICATES. Probes with another inner-track muon within dR < 0.3 are dropped,
// as in the scan's histograms (pMinDRTrkMu).
//
// MC. An optional second file -- the same scan run with isMC = true -- adds the
// MC-truth h_gen* ratios for comparison. Truth uses the scan's |eta| edges, so
// in pp its acceptance is |eta| < 2.1, not 2.0 (PbPb has the 2.0 edge); the
// printout says which. The analysis
// value 0.9708 is tight | gen muon, which also counts muons that leave no muon
// object at all -- a step the tag-and-probe cannot see (scan header). It is
// printed for orientation, not as a like-for-like comparison.
//
// PbPb. Output of src/scanning/PbPb/PbPb_muonTagAndProbe_scan.C, recognised by
// its centHiBin branch. centLo/centHi select centHiBin in [centLo, centHi), in
// the scan's 0.5% hiBin units -- any range, since it cuts the tree. When the
// range is one of the four nominal classes the analysis MC value for that class
// is drawn as the reference, and the MC truth uses that class's histograms;
// otherwise neither is shown. Check the same-sign and sideband levels in the
// mass figures per class: the PbPb tag has no isolation cut, so the background
// grows toward central.
//
// Usage (from this directory):
//   root -l -b -q 'muonTagAndProbeEfficiency.C("pp_HighEGJet_tnp.root")'
//   root -l -b -q 'muonTagAndProbeEfficiency.C("pp_HighEGJet_tnp.root", "pp_MC_tnp.root")'
//   root -l -b -q 'muonTagAndProbeEfficiency.C("PbPb_HIHardProbes_tnp.root", "", 0, 20)'   // 0-10%

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TGraphAsymmErrors.h"
#include "TEfficiency.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TBox.h"
#include "TSystem.h"
#include "TMath.h"
#include "TGaxis.h"
#include <vector>
#include <cstdio>
#include <cstring>
#include <cmath>

#include "../../headers/plotting/plotStyle.h"

namespace tpe {

const double probePtMin = 15.;     // analysis muon cut
const double etaMaxAna  = 2.0;     // analysis muon |eta| cut
const double dupDR      = 0.3;
const double effMC_analysis = 0.9708;   // tight | gen, PYTHIA, as applied in the analysis

// PbPb: the analysis MC values (PYTHIA+HYDJET, tight | gen) for the nominal
// classes, and the scan's class index for each, so the matching h_gen*_C<k>
// truth histograms are used
const int    nClsAA = 4;
const int    clsAA_lo[nClsAA]  = {0, 20, 60, 100};
const int    clsAA_hi[nClsAA]  = {20, 60, 100, 160};
const double effMC_AA[nClsAA]  = {0.8627, 0.9069, 0.9856, 0.9778};

const std::vector<double> ptEdges  = {15, 20, 25, 30, 40, 50, 70, 100};
const std::vector<double> etaEdges = {0., 0.4, 0.8, 1.2, 1.6, 2.0};

const char *figDir = "/home/clayton/Analysis/code/bJetRaaAnalysis/figures/muonReconstructionEfficiency";

enum Step { kTrk, kId, kAll, nStep };
const char *stepName [nStep] = {"trk", "id", "all"};
const char *stepLabel[nStep] = {"STA #rightarrow inner track",
                                "inner track #rightarrow tight ID",
                                "any muon #rightarrow tight ID"};

struct MassWin { double lo, hi, sb1lo, sb1hi, sb2lo, sb2hi; };
// per step: counting window, and the sidebands used by the sideband variation
const MassWin win[nStep] = {
  {71, 111, 51, 61, 121, 131}, {75, 105, 61, 71, 111, 121}, {75, 105, 61, 71, 111, 121},
};

// background variants: 0 same-sign (nominal), 1 opposite-sign sideband, 2 none
enum Bkg { kBkgSameSign, kBkgSideband, kBkgNone, nBkg };

struct Pair { float mass, pt, aeta; bool os, sta, trk, tight; float minDR; };

struct Eff {
  double e = 0, lo = 0, hi = 0, syst = 0;   // central, stat down/up, systematic
  double nPassWin = 0, nFailWin = 0;        // raw OS window counts
  bool ok = false;
};

bool isProbe(const Pair &p, int s) { return s == kTrk ? p.sta : (s == kId ? p.trk : true); }
bool isPass (const Pair &p, int s) { return s == kTrk ? p.trk : p.tight; }

// signal = OS window - background estimate, for background variant v (Bkg)
void count(const std::vector<Pair> &P, int s, int v, double ptLo, double ptHi,
           double etaLo, double etaHi, double &sPass, double &sFail,
           double &wPass, double &wFail)
{
  const MassWin &w = win[s];
  const double scale = (w.hi - w.lo) / ((w.sb1hi - w.sb1lo) + (w.sb2hi - w.sb2lo));
  double osWin[2] = {0, 0}, ssWin[2] = {0, 0}, osSB[2] = {0, 0};
  for(const Pair &p : P){
    if(p.minDR < dupDR || !isProbe(p, s)) continue;
    if(p.pt < ptLo || p.pt >= ptHi || p.aeta < etaLo || p.aeta >= etaHi) continue;
    const int k = isPass(p, s) ? 1 : 0;
    const bool inWin = (p.mass >= w.lo && p.mass < w.hi);
    const bool inSB  = (p.mass >= w.sb1lo && p.mass < w.sb1hi) || (p.mass >= w.sb2lo && p.mass < w.sb2hi);
    if(p.os){ if(inWin) osWin[k]++; else if(inSB) osSB[k]++; }
    else if(inWin) ssWin[k]++;
  }
  for(int k = 0; k < 2; k++){
    const double bkg = (v == kBkgSameSign) ? ssWin[k] : (v == kBkgSideband ? scale*osSB[k] : 0.);
    (k ? sPass : sFail) = osWin[k] - bkg;
  }
  wPass = osWin[1]; wFail = osWin[0];
}

Eff measure(const std::vector<Pair> &P, int s, double ptLo, double ptHi, double etaLo, double etaHi)
{
  Eff r;
  double sp, sf, wp, wf;
  count(P, s, kBkgSameSign, ptLo, ptHi, etaLo, etaHi, sp, sf, wp, wf);
  r.nPassWin = wp; r.nFailWin = wf;
  if(wp + wf <= 0) return r;
  sf = std::max(sf, 0.);
  if(sp <= 0) return r;
  r.e = sp / (sp + sf);
  const int N = (int)(wp + wf), K = (int)wp;
  const double eRaw = (double)K / N;
  r.lo = eRaw - TEfficiency::ClopperPearson(N, K, 0.682689, false);
  r.hi = TEfficiency::ClopperPearson(N, K, 0.682689, true) - eRaw;
  for(int v : {kBkgSideband, kBkgNone}){
    double a, b, c, d;
    count(P, s, v, ptLo, ptHi, etaLo, etaHi, a, b, c, d);
    b = std::max(b, 0.);
    if(a > 0) r.syst = std::max(r.syst, std::fabs(a/(a + b) - r.e));
  }
  r.ok = true;
  return r;
}

// MC truth ratio num/den from the scan's h_gen* histograms, pT > ptLo, |eta| < etaHi
bool truth(TFile *f, const char *num, const char *den, double &e, double &err,
           double ptLo = probePtMin, double ptHi = 1e4, double etaHi = 2.1)
{
  TH2D *n = 0, *d = 0;
  f->GetObject(num, n); f->GetObject(den, d);
  if(!n || !d) return false;
  const int x1 = d->GetXaxis()->FindBin(ptLo + 1e-6), x2 = d->GetXaxis()->FindBin(std::min(ptHi, 199.) - 1e-6);
  const int y2 = d->GetYaxis()->FindBin(etaHi - 1e-6);
  double sn = n->Integral(x1, x2, 1, y2), dErr;
  double sd = d->IntegralAndError(x1, x2, 1, y2, dErr);
  if(sd <= 0 || dErr <= 0) return false;
  e = sn / sd;
  err = std::sqrt(std::max(e*(1 - e), 0.) / (sd*sd/(dErr*dErr)));   // binomial on effective entries
  return true;
}

TGraphAsymmErrors* makeGraph(const std::vector<Eff> &E, const std::vector<double> &edges,
                             double dx, bool withSyst)
{
  TGraphAsymmErrors *g = new TGraphAsymmErrors();
  for(size_t b = 0; b < E.size(); b++){
    if(!E[b].ok) continue;
    const double x = 0.5*(edges[b] + edges[b + 1]) + dx;
    const int n = g->GetN();
    g->SetPoint(n, x, E[b].e);
    const double s = withSyst ? E[b].syst : 0.;
    g->SetPointError(n, 0, 0, std::hypot(E[b].lo, s), std::hypot(E[b].hi, s));
  }
  return g;
}

} // namespace tpe

void muonTagAndProbeEfficiency(TString dataFile = "pp_HighEGJet_tnp.root", TString mcFile = "",
                               int centLo = -1, int centHi = -1)
{
  using namespace tpe;
  initPlotStyle();

  gSystem->ExpandPathName(dataFile);
  TFile *fD = TFile::Open(dataFile);
  if(!fD || fD->IsZombie()){ printf("ERROR: cannot open %s\n", dataFile.Data()); return; }
  TTree *t = nullptr; fD->GetObject("tnp", t);
  if(!t){ printf("ERROR: no tnp tree in %s\n", dataFile.Data()); return; }
  TNamed *prov = nullptr; fD->GetObject("provenance", prov);
  if(prov) printf("\n  data provenance:\n%s\n", prov->GetTitle());

  // ---- load the pairs ------------------------------------------------------
  Float_t mass, pPt, pEta, pMinDRTrkMu;
  Int_t pairOS, pIsSTA, pHasTrk, pTightProd;
  t->SetBranchAddress("mass", &mass);       t->SetBranchAddress("pPt", &pPt);
  t->SetBranchAddress("pEta", &pEta);       t->SetBranchAddress("pMinDRTrkMu", &pMinDRTrkMu);
  t->SetBranchAddress("pairOS", &pairOS);   t->SetBranchAddress("pIsSTA", &pIsSTA);
  t->SetBranchAddress("pHasTrk", &pHasTrk); t->SetBranchAddress("pTightProd", &pTightProd);

  // PbPb trees carry centHiBin; pp trees do not
  const bool isPbPb = (t->GetBranch("centHiBin") != nullptr);
  // a scan that required a muon trigger records it in provenance; its probes
  // can be the muon that fired, which biases the efficiency up (PbPb scan header)
  const bool muonTriggered = prov && TString(prov->GetTitle()).Contains("required");
  Int_t centHiBin = -1;
  if(isPbPb) t->SetBranchAddress("centHiBin", &centHiBin);
  const bool centSel = (centLo >= 0 && centHi > centLo);
  if(centSel && !isPbPb){ printf("ERROR: a centrality range was given but %s is a pp tree\n", dataFile.Data()); return; }

  // reference value, label and figure tag for this selection
  double refEff = isPbPb ? -1. : effMC_analysis;
  int    truthCls = isPbPb ? 0 : -1;    // -1: pp names (no suffix); PbPb inclusive is C0 (0-90%)
  // every figure name carries the input file's name, so running on a second
  // dataset (HighEGJet vs SingleMuon) cannot overwrite the first one's figures
  TString dsTag = gSystem->BaseName(dataFile);
  dsTag.ReplaceAll(".root", "");
  TString sample = "pp data", figTag = "_" + dsTag;
  if(isPbPb){
    sample = centSel ? Form("PbPb data, %g-%g%%", centLo/2., centHi/2.) : "PbPb data, 0-90%";
    if(centSel) figTag += Form("_hiBin%dto%d", centLo, centHi);
    if(centSel){
      truthCls = -2;                   // no matching truth class unless found below
      for(int k = 0; k < nClsAA; k++)
        if(centLo == clsAA_lo[k] && centHi == clsAA_hi[k]){ refEff = effMC_AA[k]; truthCls = k + 1; }
    }
  }

  std::vector<Pair> P; P.reserve(t->GetEntries());
  for(Long64_t i = 0; i < t->GetEntries(); i++){
    t->GetEntry(i);
    if(centSel && (centHiBin < centLo || centHiBin >= centHi)) continue;
    P.push_back({mass, pPt, std::fabs(pEta), pairOS != 0, pIsSTA != 0, pHasTrk != 0,
                 pTightProd != 0, pMinDRTrkMu});
  }
  printf("  %zu tag-probe pairs, %s\n", P.size(), sample.Data());

  TFile *fM = nullptr;
  if(mcFile != ""){
    gSystem->ExpandPathName(mcFile);
    fM = TFile::Open(mcFile);
    if(!fM || fM->IsZombie()){ printf("  WARNING: cannot open MC file %s, skipping truth\n", mcFile.Data()); fM = nullptr; }
  }

  // ---- integrated ---------------------------------------------------------
  Eff I[nStep];
  for(int s = 0; s < nStep; s++) I[s] = measure(P, s, probePtMin, 1e4, 0., etaMaxAna);

  printf("\n  integrated: p_T > %.0f GeV, |eta| < %.1f, opposite-sign minus same-sign\n",
         probePtMin, etaMaxAna);
  printf("    step  window pass  fail     efficiency    stat (+/-)       syst\n");
  for(int s = 0; s < nStep; s++)
    printf("    %-4s  %10.0f %5.0f     %.4f      +%.4f -%.4f    %.4f\n", stepName[s],
           I[s].nPassWin, I[s].nFailWin, I[s].e, I[s].hi, I[s].lo, I[s].syst);
  const double prod = I[kTrk].e * I[kId].e;
  const double prodErr = prod * std::hypot(std::max(I[kTrk].lo, I[kTrk].hi)/I[kTrk].e,
                                           std::max(I[kId].lo, I[kId].hi)/I[kId].e);
  printf("    trk x id = %.4f +- %.4f   vs   all = %.4f   (should agree)\n", prod, prodErr, I[kAll].e);
  printf("    by background method:   same-sign   sideband   none\n");
  for(int s = 0; s < nStep; s++){
    printf("      %-4s               ", stepName[s]);
    for(int v : {kBkgSameSign, kBkgSideband, kBkgNone}){
      double a, b, c, d;
      count(P, s, v, probePtMin, 1e4, 0., etaMaxAna, a, b, c, d);
      b = std::max(b, 0.);
      printf("  %.4f  ", a > 0 ? a/(a + b) : 0.);
    }
    printf("\n");
  }
  if(muonTriggered){
    // worst case, failing probes never fire: true odds = measured odds / (2 - t)
    printf("\n    MUON-TRIGGERED SAMPLE: the efficiency above is biased up. Lower bound on the\n"
           "    true value for per-muon trigger efficiency t (failing probes never firing):\n");
    printf("    step      t = 0.95    t = 0.90    t = 0.85\n");
    for(int s = 0; s < nStep; s++){
      if(!I[s].ok || I[s].e >= 1.) continue;
      printf("    %-4s   ", stepName[s]);
      for(double tr : {0.95, 0.90, 0.85}){
        const double odds = I[s].e / (1. - I[s].e) / (2. - tr);
        printf("    %.4f  ", odds / (1. + odds));
      }
      printf("\n");
    }
  }
  if(refEff > 0)
    printf("\n    analysis MC value (tight | gen): %.4f -- not like-for-like, see header\n", refEff);
  if(fM && truthCls == -2){
    printf("    (no nominal class matches hiBin %d-%d: MC truth skipped)\n", centLo, centHi);
    fM = nullptr;
  }
  const TString sfx = truthCls >= 0 ? TString::Format("_C%d", truthCls) : TString("");

  double mcAll = -1, mcAllErr = 0, mcId = -1, mcIdErr = 0, mcAny = -1, mcAnyErr = 0;
  if(fM){
    // the PbPb scan's |eta| axis has an edge at 2.0, the pp one only at 2.1
    const double etaTruth = isPbPb ? 2.0 : 2.1;
    printf("\n    MC truth from %s (|eta| < %.1f, scan edges):\n", mcFile.Data(), etaTruth);
    for(int z = 0; z < 2; z++){
      const char *pop = z ? "FromZ" : "All";
      double e1, r1, e2, r2, e3, r3;
      bool a = truth(fM, Form("h_gen%s_tight%s", pop, sfx.Data()), Form("h_gen%s_anyMu%s", pop, sfx.Data()), e1, r1, probePtMin, 1e4, etaTruth);
      bool b = truth(fM, Form("h_gen%s_tight%s", pop, sfx.Data()), Form("h_gen%s_hasTrk%s", pop, sfx.Data()), e2, r2, probePtMin, 1e4, etaTruth);
      bool c = truth(fM, Form("h_gen%s_anyMu%s", pop, sfx.Data()), Form("h_gen%s_den%s", pop, sfx.Data()), e3, r3, probePtMin, 1e4, etaTruth);
      if(!a){ printf("      %-5s gen muons: empty\n", z ? "Z" : "all"); continue; }
      printf("      %-5s gen muons: tight|anyMu %.4f +- %.4f  tight|hasTrk %.4f +- %.4f  anyMu|gen %.4f +- %.4f\n",
             z ? "Z" : "all", e1, r1, b ? e2 : -1, b ? r2 : 0, c ? e3 : -1, c ? r3 : 0);
      // prefer the Z population, the one the tag-and-probe actually sees
      if(z == 1 || mcAll < 0){ mcAll = e1; mcAllErr = r1; mcId = e2; mcIdErr = r2; mcAny = e3; mcAnyErr = r3; }
    }
    if(mcAny > 0)
      printf("      T&P 'all' x MC anyMu|gen = %.4f  (the closest data analogue of tight | gen)\n",
             I[kAll].e * mcAny);
  }

  // ---- differential -------------------------------------------------------
  std::vector<Eff> Ept[nStep], Eeta[nStep];
  for(int s = 0; s < nStep; s++){
    for(size_t b = 0; b + 1 < ptEdges.size(); b++)
      Ept[s].push_back(measure(P, s, ptEdges[b], ptEdges[b + 1], 0., etaMaxAna));
    for(size_t b = 0; b + 1 < etaEdges.size(); b++)
      Eeta[s].push_back(measure(P, s, probePtMin, 1e4, etaEdges[b], etaEdges[b + 1]));
  }

  printf("\n  vs probe p_T (|eta| < %.1f):\n    bin        ", etaMaxAna);
  for(int s = 0; s < nStep; s++) printf("  %-26s", stepName[s]);
  printf("\n");
  for(size_t b = 0; b + 1 < ptEdges.size(); b++){
    printf("    %3.0f-%-4.0f   ", ptEdges[b], ptEdges[b + 1]);
    for(int s = 0; s < nStep; s++)
      printf("  %.3f +%.3f -%.3f (%4.0f) ", Ept[s][b].e, Ept[s][b].hi, Ept[s][b].lo,
             Ept[s][b].nPassWin + Ept[s][b].nFailWin);
    printf("\n");
  }
  printf("  vs probe |eta| (p_T > %.0f):\n", probePtMin);
  for(size_t b = 0; b + 1 < etaEdges.size(); b++){
    printf("    %.1f-%-4.1f    ", etaEdges[b], etaEdges[b + 1]);
    for(int s = 0; s < nStep; s++)
      printf("  %.3f +%.3f -%.3f (%4.0f) ", Eeta[s][b].e, Eeta[s][b].hi, Eeta[s][b].lo,
             Eeta[s][b].nPassWin + Eeta[s][b].nFailWin);
    printf("\n");
  }

  // ---- figures ------------------------------------------------------------
  gSystem->mkdir(figDir, kTRUE);
  const char *hex[nStep] = {okabeHex[2], okabeHex[6], okabeHex[0]};
  const int   mk [nStep] = {markFilledDiamond, markFilledSquare, markFilledCircle};
  const double ms[nStep] = {1.7, 1.2, 1.2};

  auto drawVs = [&](std::vector<Eff> *E, const std::vector<double> &edges, const char *xTitle,
                    const char *tag, const char *sel, double off){
    TCanvas *c = new TCanvas(Form("c_%s", tag), "", 800, 650);
    c->SetLeftMargin(0.14); c->SetRightMargin(0.04); c->SetTopMargin(0.08); c->SetBottomMargin(0.13);
    TH1D *fr = new TH1D(Form("fr_%s", tag), "", 1, edges.front(), edges.back());
    fr->SetStats(0);
    // 0.88 suits pp; central PbPb can sit lower, so widen to keep every point
    double yLo = 0.88;
    if(refEff > 0) yLo = std::min(yLo, refEff - 0.03);
    for(int s = 0; s < nStep; s++) for(const Eff &x : E[s])
      if(x.ok) yLo = std::min(yLo, x.e - std::hypot(x.lo, x.syst) - 0.01);
    fr->GetYaxis()->SetRangeUser(std::max(yLo, 0.5), 1.02);
    fr->GetXaxis()->SetTitle(xTitle); fr->GetYaxis()->SetTitle("efficiency");
    fr->GetXaxis()->SetTitleOffset(1.15); fr->GetYaxis()->SetTitleOffset(1.25);
    fr->Draw("axis");
    TLine *one = new TLine(edges.front(), 1., edges.back(), 1.);
    one->SetLineStyle(2); one->SetLineColor(kGray + 2); one->Draw();
    TLine *ref = new TLine(edges.front(), refEff, edges.back(), refEff);
    ref->SetLineStyle(7); ref->SetLineWidth(2); ref->SetLineColor(kGray + 1);
    if(refEff > 0) ref->Draw();
    TLegend *leg = makeLegend(0.44, 0.15, 0.95, 0.42, 0.032);
    leg->SetHeader(Form("Z #rightarrow #mu#mu tag & probe, %s", sample.Data()), "L");
    for(int s = 0; s < nStep; s++){
      TGraphAsymmErrors *g = makeGraph(E[s], edges, (s - 1)*off, true);
      g->SetLineColor(TColor::GetColor(hex[s])); g->SetMarkerColor(TColor::GetColor(hex[s]));
      g->SetMarkerStyle(mk[s]); g->SetMarkerSize(ms[s]); g->SetLineWidth(2);
      g->Draw("pz same");
      leg->AddEntry(g, stepLabel[s], "lp");
    }
    if(refEff > 0)
      leg->AddEntry(ref, isPbPb ? "PYTHIA+HYDJET tight | gen (analysis value)"
                                : "PYTHIA tight | gen (analysis value)", "l");
    leg->Draw();
    TLatex tx; tx.SetNDC(); tx.SetTextFont(42); tx.SetTextSize(0.034);
    tx.DrawLatex(0.14, 0.945, Form("muon efficiency, %s", sel));
    tx.SetTextSize(0.028);
    tx.DrawLatex(0.17, 0.86, "errors: stat (Clopper-Pearson) #oplus background syst");
    savePdfTight(c, Form("%s/muonTnP_eff_%s%s.pdf", figDir, tag, figTag.Data()));
  };
  drawVs(Ept,  ptEdges,  "probe muon p_{T} [GeV]", "vsPt",  "|#eta| < 2.0", 1.2);
  drawVs(Eeta, etaEdges, "probe muon |#eta|",      "vsEta", "p_{T} > 15 GeV", 0.04);

  // mass distributions, pass and fail, with the window and sidebands marked.
  // Drawn integrated for 'all' and 'trk', and per probe-pT bin for 'all' as
  // backup: the per-bin efficiencies above are only as good as the background
  // estimate in each bin, and these show it.
  auto drawMass = [&](int s, double ptLo, double ptHi, const TString &name, const TString &ptLabel){
    TCanvas *c = new TCanvas("cm_" + name, "", 1200, 550);
    c->Divide(2, 1);
    double sp, sf, wp, wf, bp, bf;
    count(P, s, kBkgSameSign, ptLo, ptHi, 0., etaMaxAna, sp, sf, wp, wf);
    bp = wp - sp; bf = wf - sf;    // same-sign estimate, pass / fail
    const Eff E = measure(P, s, ptLo, ptHi, 0., etaMaxAna);
    for(int k = 1; k >= 0; k--){
      c->cd(2 - k);
      gPad->SetLeftMargin(0.15); gPad->SetRightMargin(0.04); gPad->SetBottomMargin(0.13);
      TH1D *os = new TH1D("os_" + name + Form("_%d", k), "", 50, 40, 140);
      TH1D *ss = new TH1D("ss_" + name + Form("_%d", k), "", 50, 40, 140);
      os->SetDirectory(nullptr); ss->SetDirectory(nullptr);
      for(const Pair &p : P){
        if(p.minDR < dupDR || !isProbe(p, s) || isPass(p, s) != (k == 1)) continue;
        if(p.pt < ptLo || p.pt >= ptHi || p.aeta >= etaMaxAna) continue;
        (p.os ? os : ss)->Fill(p.mass);
      }
      styleH(os, okabeHex[0], markFilledCircle, 0.9);
      styleH(ss, okabeHex[5], markOpenSquare, 0.9);
      // log: the peak and the continuum under it differ by up to 100x at low pT,
      // and the continuum's shape across the sidebands is the thing to judge
      gPad->SetLogy();
      const double ymin = 0.5, ymax = 20.*std::max(2., os->GetMaximum());
      os->SetMinimum(ymin); os->SetMaximum(ymax);
      os->GetXaxis()->SetTitle("m_{#mu#mu} [GeV]"); os->GetYaxis()->SetTitle("pairs / 2 GeV");
      os->Draw("axis");
      const MassWin &w = win[s];
      TLine *lWin = nullptr; TBox *bSB = nullptr;
      for(double x : {w.lo, w.hi}){ lWin = new TLine(x, ymin, x, ymax); lWin->SetLineColor(kGray + 2); lWin->SetLineStyle(2); lWin->Draw(); }
      for(auto sb : {std::make_pair(w.sb1lo, w.sb1hi), std::make_pair(w.sb2lo, w.sb2hi)}){
        bSB = new TBox(sb.first, ymin, sb.second, ymax);
        bSB->SetFillColorAlpha(TColor::GetColor(okabeHex[1]), 0.25); bSB->SetLineWidth(0); bSB->Draw();
      }
      // as graphs without the empty bins: a marker at zero reads as a measured zero
      auto toG = [](TH1D *h){
        TGraphAsymmErrors *g = new TGraphAsymmErrors();
        for(int b = 1; b <= h->GetNbinsX(); b++){
          const double n = h->GetBinContent(b);
          if(n <= 0) continue;
          const int i = g->GetN();
          g->SetPoint(i, h->GetBinCenter(b), n);
          g->SetPointError(i, 0, 0, n - 0.5*TMath::ChisquareQuantile(0.1586553, 2*n),
                           0.5*TMath::ChisquareQuantile(0.8413447, 2*(n + 1)) - n);   // Garwood
        }
        g->SetLineColor(h->GetLineColor()); g->SetMarkerColor(h->GetMarkerColor());
        g->SetMarkerStyle(h->GetMarkerStyle()); g->SetMarkerSize(h->GetMarkerSize()); g->SetLineWidth(2);
        return g; };
      TGraphAsymmErrors *gOS = toG(os), *gSS = toG(ss);
      gOS->Draw("pz same"); gSS->Draw("pz same");
      TLegend *leg = makeLegend(0.66, 0.70, 0.95, 0.88, 0.030);
      leg->AddEntry(gOS, "opposite sign", "lp"); leg->AddEntry(gSS, "same sign", "lp");
      leg->AddEntry(lWin, Form("window %.0f-%.0f", w.lo, w.hi), "l");
      leg->AddEntry(bSB, "sidebands (syst.)", "f");
      leg->Draw();
      TLatex tx; tx.SetNDC(); tx.SetTextFont(42); tx.SetTextSize(0.040);
      tx.DrawLatex(0.15, 0.93, Form("%s: %s probes, %s", stepLabel[s], k ? "passing" : "failing", ptLabel.Data()));
      tx.SetTextSize(0.030);
      tx.DrawLatex(0.19, 0.85, Form("OS window %.0f, SS window %.0f", k ? wp : wf, k ? bp : bf));
      if(k == 0 && E.ok)
        tx.DrawLatex(0.19, 0.80, Form("#varepsilon = %.4f ^{+%.4f}_{-%.4f} (stat) #pm %.4f (syst)", E.e, E.hi, E.lo, E.syst));
    }
    savePdfTight(c, Form("%s/%s.pdf", figDir, name.Data()));
    delete c;
  };

  // full integers on the count axis: the shared x10^n header lands on the title
  TGaxis::SetMaxDigits(6);
  const TString cls = isPbPb ? TString(", ") + (sample.Data() + strlen("PbPb data, ")) : TString("");
  for(int s : {kAll, kTrk})
    drawMass(s, probePtMin, 1e4, Form("muonTnP_mass_%s%s", stepName[s], figTag.Data()),
             Form("p_{T} > %.0f GeV%s", probePtMin, cls.Data()));

  const TString ptDir = TString(figDir) + "/massPtBins";
  gSystem->mkdir(ptDir, kTRUE);
  for(size_t b = 0; b + 1 < ptEdges.size(); b++)
    drawMass(kAll, ptEdges[b], ptEdges[b + 1],
             Form("massPtBins/muonTnP_mass_all_pt%.0fto%.0f%s", ptEdges[b], ptEdges[b + 1], figTag.Data()),
             Form("%.0f < p_{T} < %.0f GeV%s", ptEdges[b], ptEdges[b + 1], cls.Data()));
  printf("\n  figures in %s (per-pT-bin mass plots in massPtBins/)\n", figDir);
}
