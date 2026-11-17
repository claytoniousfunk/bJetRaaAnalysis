// b-purity systematics of the 2D muon ptRel x CSVv2 fit, PF pp and the four PbPb classes, and the chain's purity files for it
// (2026-11-12, user: "i want to make the 2d fitting the nominal method"; choices: systematics derived for the 2D fit, no
// working-point source, new tag _csvV2_2D). Copy of systematics_CSVv2Cut_PF.C with the purity function swapped:
//   P = f_b of the HistFactory 2D fit (p2dhf::fitHF of fitPtRelCSVv2_2D_ppPbPb_PF.C: one channel per CSVv2 category, 4cat
//       edges 0.54 / 0.84 / 0.96, b and c + light (MC ratio) samples, MC tagger shapes, Barlow-Beeston-lite gammas), the b
//       fraction of ALL muon-tagged jets over ptRel 0-5 GeV, times (fake-subtracted / all muon-tagged jets) in PbPb, so it
//       is b / all muon-tagged jets like the CSVv2-cut purity. No eps_b division.
// Sources, per system and window 80-300 GeV, max |P_var / P_nom - 1| over the variations, as in the cut macro:
//   bGS           g->bb share x 0.94, x 1.41 (nominal natural x 1.175)                       correlated pp-PbPb (h_varRel)
//   cMult         c yield x 0.8, 0.9, 1.1, 1.2 (in every CSVv2 category)                    correlated pp-PbPb (h_varRel)
//   fitRangeHigh  ZERO since 2026-11-16 (user: the 2D fit uses the whole CSVv2 and ptRel range simultaneously, so there
//                 is no fit-range or working-point choice to vary); kept because the pp chain reads it by name
//   tplStat       ZERO, kept because the pp chain reads it: there is no eps_b, and the template MC statistics are in the
//                 2D fit's statistical error through the gammas
//   JEU, JER      template ptRel scaled by T_var / T_0 per flavor (the same in every CSVv2 category, as in the cut macro)
//   fakeMuJet     PbPb: |P without / P with the per-event uncorrelated-pair subtraction - 1|
// NO working-point source: the chain is run without an alternative-WP chain (bJetsPerZSys_PF_PbPb.C 7th argument empty).
// Outputs: rootFiles/systematics/bPuritySys_total_PFJets_pp_csvV2_2D.root, bPuritySys_PFJets_PbPb_csvV2_2D.root (8cat:
//   _csvV2_2D8cat), read by the chain with BJET_PURITY_VARIANT=csvV2__2D (suffix "_csvV2" + "_2D"); buffer bins 50-80 and
//   300-500 GeV = the ptRel-only nominal x the edge-window ratio, as the cut macro; figure
//   figures/systematics/bPurity/CSVv2/bPuritySys_CSVv2_2D.pdf
// Usage, from src/systematics/:  root -l -b -q 'systematics_CSVv2_2D_PF.C("4cat")'
#include "../plots/bPurity/fitPtRelCSVv2_2D_ppPbPb_PF.C"   // includes fitPtRelCSVv2Cut_PF_PbPb.C; p2dhf::fitHF
#include "../../headers/plotting/puritySysFigure.h"

namespace csys2d {
using namespace p2d; using namespace pcut; using namespace pbc;
const int NSY = 5;                                  // pp, C1..C4
const char *syTag[NSY] = {"pp", "C1", "C2", "C3", "C4"}, *syName[NSY] = {"pp", "0-10%", "10-30%", "30-50%", "50-80%"};
const int NSRC = 7;
const char *srcKey[NSRC] = {"bGS", "cMult", "fitRangeHigh", "tplStat", "JEU", "JER", "fakeMuJet"};
const char *srcLab[NSRC] = {"g#rightarrowb#bar{b} share", "c fraction", "fit range (upper edge)", "#varepsilon_{b} MC stat.", "JES", "JER", "uncorrelated #mu-jet pairs"};
TFile *fDpp = nullptr, *fMpp = nullptr, *fDAA = nullptr;
void scaleFl(Win &X, int F, double m){
  X.t[F].tot = 0.;
  for(int k = 0; k < K; k++) for(int j = 0; j < nPR; j++){ X.t[F].n[k][j] *= m; X.t[F].e2[k][j] *= m*m; X.t[F].tot += X.t[F].n[k][j]; }
  X.truthB = X.t[0].tot/(X.t[0].tot + X.t[1].tot + X.t[2].tot); X.rc = X.t[1].tot/(X.t[1].tot + X.t[2].tot);
}
bool build(int sy, int w, Win &X){
  if(sy == 0){ if(!buildWindows(fDpp, fMpp)) return false; X = W[w]; return true; }
  double nEffB = 0.; return buildCls(fDAA, sy - 1, w, X, nEffB);
}
struct Pc { double p = 0., e = 0., eps = 0., epsRel = 0.; };
int hfUid = 0;
// the 2D HistFactory fit of window X: b / all muon-tagged jets (X.nData is the unsubtracted total, X.d the fitted data)
Pc pcorr(Win X){
  double sd = 0.; for(int k = 0; k < K; k++) for(int j = 0; j < nPR; j++) sd += X.d[k][j];
  const double kAll = X.nData > 0. ? sd/X.nData : 1.;
  const p2dhf::HFResult r = p2dhf::fitHF(X, X.d, hfUid++);
  Pc o; o.p = r.pur*kAll; o.e = r.epur*kAll; o.eps = 1.; o.epsRel = 0.;   // no eps_b; template statistics are in the fit
  if(r.status != 0) printf("    NOTE: HistFactory fit status %d\n", r.status);
  return o;
}
// T_var / T_0 of flavor F per 0.5 GeV ptRel bin, in the jet-pT window: pp the PYTHIA map; PbPb the samples and yields
// of buildCls (shapes weighted by effective entries, flavor groups at their DiJet yields)
bool slot(TFile *f, TString nm, double lo, double hi, double *a, double *e2){
  TH2 *h = nullptr; f->GetObject(nm, h);
  if(!h){ printf("ERROR: missing %s\n", nm.Data()); return false; }
  const int y1 = h->GetYaxis()->FindBin(lo + 1e-6), y2 = h->GetYaxis()->FindBin(hi - 1e-6);
  for(int ix = 1; ix <= h->GetNbinsX(); ix++){
    const double x = h->GetXaxis()->GetBinCenter(ix); if(x >= 5.) continue;
    const int j = (int)(x/ptRelBW);
    for(int iy = y1; iy <= y2; iy++){ a[j] += h->GetBinContent(ix, iy); if(e2) e2[j] += pow(h->GetBinError(ix, iy), 2); }
  }
  return true;
}
// explicit window edges (2026-11-16, for studyJetPtBinning_sys_2D_PF.C); the index version below forwards to it
bool tRatio(int sy, double lo, double hi, int F, int T, double *R){
  double v[NPR] = {}, u[NPR] = {};
  if(sy == 0){
    for(const char *fl : flMembers[F]){
      if(!slot(fMpp, Form("%s_%s_T%d", h2, fl, T), lo, hi, v, nullptr)) return false;
      if(!slot(fMpp, Form("%s_%s_T0", h2, fl), lo, hi, u, nullptr)) return false;
    }
  } else {
    const int c = sy - 1;
    const std::vector<int> grp = F == 0 ? std::vector<int>{0, 1} : F == 1 ? std::vector<int>{2} : std::vector<int>{3};
    for(int G : grp){
      double yG = 0.;   // DiJet yield of the group
      for(int s : gSamples[G]){
        double n0[NPR] = {}, nE[NPR] = {};
        for(const char *fl : gMembers[G]) for(int i = slc[c][0]; i <= slc[c][1]; i++) if(!slot(fS[s], Form("%s_%s_C%dT0", h2, fl, i), lo, hi, n0, nE)) return false;
        if(s == 0) for(int j = 0; j < nPR; j++) yG += n0[j];
      }
      double ws = 0.; double mv[NPR] = {}, mu[NPR] = {};
      for(int s : gSamples[G]){
        double n0[NPR] = {}, nE[NPR] = {}, nv[NPR] = {};
        for(const char *fl : gMembers[G]) for(int i = slc[c][0]; i <= slc[c][1]; i++){
          if(!slot(fS[s], Form("%s_%s_C%dT0", h2, fl, i), lo, hi, n0, nE)) return false;
          if(!slot(fS[s], Form("%s_%s_C%dT%d", h2, fl, i, T), lo, hi, nv, nullptr)) return false;
        }
        double tt = 0., ee = 0.; for(int j = 0; j < nPR; j++){ tt += n0[j]; ee += nE[j]; }
        if(tt <= 0. || ee <= 0.) continue;
        const double ne = tt*tt/ee; ws += ne;
        for(int j = 0; j < nPR; j++){ mv[j] += ne*nv[j]/tt; mu[j] += ne*n0[j]/tt; }
      }
      if(ws > 0.) for(int j = 0; j < nPR; j++){ v[j] += yG*mv[j]/ws; u[j] += yG*mu[j]/ws; }
    }
  }
  for(int j = 0; j < nPR; j++) R[j] = u[j] > 0. ? v[j]/u[j] : 1.;
  return true;
}
bool tRatio(int sy, int w, int F, int T, double *R){ return tRatio(sy, wLo[w], wHi[w], F, T, R); }
// P with the templates at T_var: eps_b stays nominal (scaleFl leaves eps alone too)
bool shifted(int sy, int w, const Win &X, int T, double &p){
  Win Y = X;
  for(int F = 0; F < NF; F++){
    double R[NPR]; if(!tRatio(sy, w, F, T, R)) return false;
    Y.t[F].tot = 0.;
    for(int k = 0; k < K; k++) for(int j = 0; j < nPR; j++){ Y.t[F].n[k][j] *= R[j]; Y.t[F].e2[k][j] *= R[j]*R[j]; Y.t[F].tot += Y.t[F].n[k][j]; }
  }
  Y.truthB = Y.t[0].tot/(Y.t[0].tot + Y.t[1].tot + Y.t[2].tot); Y.rc = Y.t[1].tot/(Y.t[1].tot + Y.t[2].tot);
  p = pcorr(Y).p;
  return true;
}
} // namespace csys2d

// outTag (2026-11-12): appended to the output purity files and the figure after fitEngineTag(); empty with fitMax != 5
// gives "_fit<fitMax>", so a shorter fit range never overwrites the 0-5 GeV files (chain: BJET_PURITY_VARIANT=csvV2_<wp>_fit2)
void systematics_CSVv2_2D_PF(TString cats = "4cat", double fitMax = 5.0, TString outTag = ""){
  using namespace csys2d;
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
  setvbuf(stdout, nullptr, _IOLBF, 0);
  if(outTag == "" && fabs(fitMax - 5.) > 1e-6) outTag = Form("_fit%g", fitMax);
  const TString wpTagS = cats == "4cat" ? TString("_2D") : "_2D" + cats; const char *wpN = wpTagS.Data();   // file tag after "csvV2"
  const double wp = 0.;   // unused (labels only)
  if(!setCats(cats)){ printf("ERROR: categories 4cat or 8cat\n"); return; }
  setJetCollection("PFJets"); initPlotStyle();
  fDpp = TFile::Open(Form("%s/%s", repo, dataName)); fMpp = TFile::Open(Form("%s/%s", repo, mcName)); fDAA = TFile::Open(Form("%s/%s", repo, dataPbPb));
  for(int s = 0; s < NS; s++) fS[s] = TFile::Open(Form("%s/%s", repo, Form(phFmt, sName[s])));
  TFile *fNpp = TFile::Open(Form("%s/rootFiles/systematics/bPuritySys_total_PFJets_pp.root", repo)), *fNAA = TFile::Open(Form("%s/rootFiles/systematics/bPuritySys_PFJets_PbPb.root", repo));
  if(!fDpp || !fMpp || !fDAA || !fS[0] || !fS[1] || !fS[2] || !fNpp || !fNAA){ printf("ERROR: inputs\n"); return; }
  const int nNom = (int)(fitMax/ptRelBW + 0.5);

  double P[NSY][NW], E[NSY][NW], S[NSY][NSRC][NW] = {};
  double V[NSY][3][4][NW] = {};   // signed P_var / P_nom - 1 of the pp-PbPb correlated sources (bGS, cMult) per variation
  for(int sy = 0; sy < NSY; sy++){
    printf("\n==== %s, 2D ptRel x CSVv2 (%s, HistFactory), fit 0-%.1f GeV: relative systematic (%%) of P = f_b (b / all muon-tagged jets)\n", syName[sy], cats.Data(), fitMax);
    printf("  window    P_nom         "); for(int s = 0; s < NSRC; s++) printf(" %12s", srcKey[s]); printf("  total\n");
    for(int w = 0; w < NW; w++){
      auto dev = [&](double v){ return fabs(v/P[sy][w] - 1.); };
      nFit = nNom; gsMult = 1.175;
      Win X; if(!build(sy, w, X)) return;
      Pc n = pcorr(X); P[sy][w] = n.p; E[sy][w] = n.e;
      { int k = 0; for(double g : {0.94, 1.41}){ gsMult = g; Win Y; if(!build(sy, w, Y)) return; const double p = pcorr(Y).p; V[sy][0][k++][w] = p/P[sy][w] - 1.; S[sy][0][w] = TMath::Max(S[sy][0][w], dev(p)); } }
      gsMult = 1.175;
      { int k = 0; for(double m : {0.8, 0.9, 1.1, 1.2}){ Win Y = X; scaleFl(Y, 1, m); const double p = pcorr(Y).p; V[sy][1][k++][w] = p/P[sy][w] - 1.; S[sy][1][w] = TMath::Max(S[sy][1][w], dev(p)); } }
      S[sy][2][w] = 0.;   // no fit-range source for the 2D fit (2026-11-16, user)
      S[sy][3][w] = n.epsRel;
      const double c = 0.5*(wLo[w] + wHi[w]);
      double old[2] = {};
      for(int s = 4; s < 6; s++){
        TH1 *h = sy == 0 ? (TH1*) fNpp->Get(Form("h_sysRel_%s", srcKey[s])) : (TH1*) fNAA->Get(Form("h_sysRel_%s_%s", srcKey[s], syTag[sy]));
        if(!h){ printf("ERROR: ptRel-only %s missing\n", srcKey[s]); return; }
        old[s-4] = h->GetBinContent(h->FindBin(c));
      }
      for(int T : {2, 3}){ double p; if(!shifted(sy, w, X, T, p)) return; S[sy][4][w] = TMath::Max(S[sy][4][w], dev(p)); }
      { double p; if(!shifted(sy, w, X, 1, p)) return; S[sy][5][w] = dev(p); }
      if(sy > 0){ fakeSubPerEvt = false; Win Y; if(!build(sy, w, Y)){ fakeSubPerEvt = true; return; } S[sy][6][w] = dev(pcorr(Y).p); fakeSubPerEvt = true;
        build(sy, w, Y); printf("    (%s %3.0f-%3.0f: uncorrelated muon-jet pairs %.2f%% of the muon-tagged jets removed)\n", syName[sy], wLo[w], wHi[w], 100.*lastFakeFrac); }
      double t = 0.; for(int s = 0; s < NSRC; s++) t += S[sy][s][w]*S[sy][s][w];
      printf("  %3.0f-%3.0f   %.3f +- %.3f", wLo[w], wHi[w], n.p, n.e); for(int s = 0; s < NSRC; s++) printf(" %12.1f", 100*S[sy][s][w]); printf("  %5.1f   (old JEU %.1f, JER %.1f)\n", 100*sqrt(t), 100*old[0], 100*old[1]);
    }
    // A fit at the P = 0 limit (e.g. medium, fit 0-2 GeV, 50-80% 200-300 GeV) has no relative uncertainty: every
    // source is inf / nan there and the nan reaches every unfolded bin and the results figure. Such a window takes
    // the sources of the nearest window with a finite purity (2026-11-12).
    for(int w = 0; w < NW; w++){
      bool bad = !(P[sy][w] > 0.); for(int q = 0; q < NSRC; q++) if(!std::isfinite(S[sy][q][w])) bad = true;
      if(!bad) continue;
      int u = -1; for(int d = 1; d < NW && u < 0; d++) for(int k : {w - d, w + d}) if(u < 0 && k >= 0 && k < NW && P[sy][k] > 0.){
        bool ok = true; for(int q = 0; q < NSRC; q++) if(!std::isfinite(S[sy][q][k])) ok = false; if(ok) u = k; }
      if(u < 0){ printf("ERROR: %s: no window with a finite purity\n", syName[sy]); return; }
      printf("  NOTE: %s %3.0f-%3.0f GeV purity %.3f at the fit limit: sources taken from %3.0f-%3.0f GeV\n", syName[sy], wLo[w], wHi[w], P[sy][w], wLo[u], wHi[u]);
      for(int q = 0; q < NSRC; q++) S[sy][q][w] = S[sy][q][u];
      for(int cc = 0; cc < 3; cc++) for(int k = 0; k < 4; k++) V[sy][cc][k][w] = V[sy][cc][k][u];
    }
  }

  // ---- the chain's purity files
  const TString base = Form("%s/rootFiles/systematics/", repo);
  for(int sy = 0; sy < NSY; sy++){
    const bool pp = sy == 0;
    const TString sfx = pp ? "" : Form("_%s", syTag[sy]);
    TFile *fN = pp ? fNpp : fNAA;
    TH1D *pN = (TH1D*) fN->Get("h_bPurity_nominal" + sfx);
    if(!pN){ printf("ERROR: nominal purity %s\n", syTag[sy]); return; }
    TH1D *p = (TH1D*) pN->Clone("h_bPurity_nominal" + sfx); p->SetDirectory(nullptr);
    TH1D *tot = (TH1D*) pN->Clone("h_sysRel_total" + sfx); tot->SetDirectory(nullptr); tot->Reset();
    TH1D *src[NSRC]; for(int s = 0; s < NSRC; s++){ src[s] = (TH1D*) tot->Clone(Form("h_sysRel_%s%s", srcKey[s], sfx.Data())); src[s]->SetDirectory(nullptr); src[s]->SetTitle(srcLab[s]); }
    for(int i = 1; i <= p->GetNbinsX(); i++){
      const double c = p->GetBinCenter(i);
      const int w = c < wLo[0] ? 0 : c > wHi[NW-1] ? NW - 1 : [&]{ int k = 0; while(k < NW - 1 && c > wHi[k]) k++; return k; }();
      const bool inside = c > wLo[0] && c < wHi[NW-1];
      const double r = P[sy][w]/pN->GetBinContent(pN->FindBin(0.5*(wLo[w] + wHi[w])));
      if(inside){ p->SetBinContent(i, P[sy][w]); p->SetBinError(i, E[sy][w]); }
      else { p->SetBinContent(i, pN->GetBinContent(i)*r); p->SetBinError(i, pN->GetBinError(i)*r); }
      double t = 0.; for(int s = 0; s < NSRC; s++){ src[s]->SetBinContent(i, S[sy][s][w]); t += S[sy][s][w]*S[sy][s][w]; }
      tot->SetBinContent(i, sqrt(t));
    }
    // 2026-10-28 (user): bGS, cMult correlated between pp and PbPb -> signed per variation
    // h_varRel_<src>_<k>; the rest (fitRangeHigh, tplStat, JEU, JER, fakeMuJet) in quadrature as h_sysRel_uncorr
    const int nVar[2] = {2, 4};
    TH1D *unc = (TH1D*) tot->Clone("h_sysRel_uncorr" + sfx); unc->SetTitle("uncorrelated purity sources (JEU, JER, fakeMuJet; fitRangeHigh and tplStat are zero)");
    std::vector<TH1D*> vr;
    for(int c = 0; c < 2; c++) for(int k = 0; k < nVar[c]; k++){ TH1D *h = (TH1D*) tot->Clone(Form("h_varRel_%s_%d%s", srcKey[c], k, sfx.Data())); h->SetDirectory(nullptr); h->SetTitle(Form("%s variation %d: P_var / P_nom - 1 (signed)", srcKey[c], k)); vr.push_back(h); }
    for(int i = 1; i <= p->GetNbinsX(); i++){
      const double c = p->GetBinCenter(i);
      const int w = c < wLo[0] ? 0 : c > wHi[NW-1] ? NW - 1 : [&]{ int k = 0; while(k < NW - 1 && c > wHi[k]) k++; return k; }();
      double t = 0.; for(int q = 2; q < NSRC; q++) t += S[sy][q][w]*S[sy][q][w];
      unc->SetBinContent(i, sqrt(t));
      int n = 0; for(int cc = 0; cc < 2; cc++) for(int k = 0; k < nVar[cc]; k++) vr[n++]->SetBinContent(i, V[sy][cc][k][w]);
    }
    TH1D *abs = (TH1D*) tot->Clone("h_sysAbs_total" + sfx), *pa = (TH1D*) p->Clone("h_bPurity_sysAbs" + sfx);
    for(int i = 1; i <= p->GetNbinsX(); i++){ const double a = tot->GetBinContent(i)*p->GetBinContent(i); abs->SetBinContent(i, a); abs->SetBinError(i, 0.); pa->SetBinError(i, a); }
    const TString outName = pp ? base + Form("bPuritySys_total_PFJets_pp_csvV2%s%s.root", wpN, (fitEngineTag() + outTag).Data()) : base + Form("bPuritySys_PFJets_PbPb_csvV2%s%s.root", wpN, (fitEngineTag() + outTag).Data());
    TFile *fo = TFile::Open(outName, pp || sy == 1 ? "RECREATE" : "UPDATE");
    p->Write(); for(int s = 0; s < NSRC; s++) src[s]->Write(); tot->Write(); abs->Write(); pa->Write(); unc->Write(); for(TH1D *h : vr) h->Write();
    if(pp || sy == 1){
      TNamed info("info", Form("%s PF b purity from the 2D muon ptRel x CSVv2 fit (%s categories, HistFactory, MC tagger shapes, c / light at MC; f_b = b / all muon-tagged jets, fit 0-%.1f GeV), 80-300 GeV; buffer bins = ptRel-only nominal x the edge-window ratio. "
                               "Systematics from systematics_CSVv2_2D_PF.C: bGS, cMult, JEU, JER (template ptRel scaled by T_var / T_0) derived with the 2D fit; fitRangeHigh = 0 (whole range fitted, 2026-11-16); tplStat = 0 (template statistics in the fit error); PbPb: uncorrelated muon-jet pairs subtracted per event, fakeMuJet = |without / with - 1|. No working-point source.",
                               pp ? "pp" : "PbPb (C1..C4)", cats.Data(), fitMax));
      info.Write();
    }
    fo->Close();
    printf("wrote %s (%s)\n", outName.Data(), syTag[sy]);
  }

  // ---- figure: relative systematic per source, pp then peripheral to central, 80-200 GeV
  // (drawing shared with plotBPuritySys_CSVv2_PF.C("_2D"), headers/plotting/puritySysFigure.h)
  const TString figDir = Form("%s/figures/systematics/bPurity/CSVv2/", repo); gSystem->mkdir(figDir, true);
  const int order[NSY] = {0, 4, 3, 2, 1};
  TH1D *fTot[5], *fSrc[5][PSF_NMAX] = {{nullptr}};
  for(int ix = 0; ix < NSY; ix++){
    const int sy = order[ix];
    fTot[ix] = new TH1D(Form("ht%d", sy), "", NW, ptEdges + 3); fTot[ix]->SetDirectory(nullptr);
    for(int w = 0; w < NW; w++){ double t = 0; for(int s = 0; s < NSRC; s++) t += S[sy][s][w]*S[sy][s][w]; fTot[ix]->SetBinContent(w + 1, sqrt(t)); }
    for(int s = 0; s < NSRC; s++){
      if(s == 6 && sy == 0) continue;   // fakeMuJet: PbPb only
      fSrc[ix][s] = new TH1D(Form("hs%d_%d", sy, s), "", NW, ptEdges + 3); fSrc[ix][s]->SetDirectory(nullptr);
      for(int w = 0; w < NW; w++) fSrc[ix][s]->SetBinContent(w + 1, S[sy][s][w]);
    }
  }
  drawPuritySysRow(fTot, fSrc, NSRC, srcKey, 80., 200., Form("muon #it{p}_{T}^{rel} #times CSVv2 fit%s", cats == "4cat" ? "" : " (8 categories)"),
                   figDir + Form("bPuritySys_CSVv2%s%s.pdf", wpN, (fitEngineTag() + outTag).Data()));
}
