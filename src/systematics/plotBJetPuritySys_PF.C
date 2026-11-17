// Extraction of the b-PURITY systematic sources, from each variation to every assigned value (2026-11-17, user: figures
// that justify every systematic value, one source at a time; the purity after JES and JER).
//
// The purity sources (systematics_CSVv2_2D_PF.C; 7 bins: studyJetPtBinning_sys_2D_PF.C) and their propagation
// (bJetsPerZSys_PF_PbPb.C):
//   g->bb share  x 0.94, x 1.41 of the nominal (natural x 1.175); signed, CORRELATED pp-PbPb: every variation applied
//                to every system at once, the ratios formed per variation, assigned = max |ratio - 1|
//   c fraction   the charm yield x 0.8, 0.9, 1.1, 1.2; the same treatment
//   other        the purity's JES and JER (taken from the four wide windows in the 7-bin binning) and, in PbPb, the
//                uncorrelated muon-jet pairs, in quadrature; +- that, UNCORRELATED between systems: quadrature in
//                the ratios
// PbPb: the measured b spectrum scaled bin by bin by (1 + purity shift) and unfolded at the nominal N. pp: the
// reco-level shift times the pp chain's unfolded / reco ratio of the total purity systematic (the pp unfolding is
// not redone in bJetsPerZSys_PF_PbPb.C). The inclusive jets have no purity, so the double ratio and the b-jet
// fraction carry exactly the purity uncertainty of R_Z(b) and of the b jets per Z: not drawn twice.
//
// Figure 1 (purity, per system): rows g->bb, c fraction, other; signed P_var / P_nom - 1 per variation (other: JES,
//   JER, muon-jet pairs, magnitudes), band = the purity file's value. Figure 2 (b jets per Z): the same after the
//   unfolding, band = assigned. Figures 3 and 4: R_Z(b) per class and R_CP^Z(b) per class, band = assigned; for
//   the uncorrelated part the class (solid) and the reference (dashed) magnitudes that add in quadrature.
// Every figure has one y range. The macro checks every band (max |variation| for the correlated sources, the
// quadrature sum for the uncorrelated part) and prints the largest difference.
//
// inputs  rootFiles/systematics/bJetsPerZSys_PF_MBcaloMatched_JEUppCorr<bTag>.root (h_bPur_*, saved since 2026-11-17;
//         h_sysRel_bPurity*), the purity files bPuritySys_{total_PFJets_pp,PFJets_PbPb}<purTag>.root
// Usage, from src/systematics/:  root -l -b -q plotBJetPuritySys_PF.C
// Output: figures/systematics/bPurity/bJetPuritySys_{purity,perSystem,RZ,RCP}_PF<outSfx>.pdf
#include "TFile.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TSystem.h"
#include "TVirtualPad.h"
#include <vector>
#include "../../headers/plotting/plotStyle.h"
#include "../../headers/plotting/gridPads.h"
#include "../../headers/plotting/sysStyle.h"

namespace bpur {
const double xLo = 80., xHi = 200.;
double maxDiff = 0.;
double pxText = 34., pxLeg = 30.;
void setTextScale(TCanvas *cv){ const double k = cv->GetWw()/1800.; pxText = 34.*k; pxLeg = 30.*k; }
double padHpx(){ return gPad->GetWh()*gPad->GetAbsHNDC(); }
double padWpx(){ return gPad->GetWw()*gPad->GetAbsWNDC(); }
void ptext(double x, double y, const char *s, double px, int font = 43, int align = 13, Color_t col = kBlack){
  TLatex t; t.SetNDC(); t.SetTextFont(font); t.SetTextSize(px); t.SetTextAlign(align); t.SetTextColor(col); t.DrawLatex(x, y, s); }
void panelLabel(GridPads &g, int ix, int iy, const char *s){ ptext(g.frameLeft(ix) + 20./padWpx(), g.frameTop(iy) - 14./padHpx(), s, pxText, 63, 13); }
double minDim(){ return std::min(padWpx(), padHpx()); }
void sameLabels(TH1F *f, double px){ if(f->GetXaxis()->GetLabelSize() > 0.) f->GetXaxis()->SetLabelSize(px/minDim()); if(f->GetYaxis()->GetLabelSize() > 0.) f->GetYaxis()->SetLabelSize(px/minDim()); }
bool inRange(TH1 *h, int i){ const double lo = h->GetXaxis()->GetBinLowEdge(i); return lo >= xLo - 1e-6 && lo < xHi - 1e-6; }
TH1D* pct(TH1D *h, const char *n, bool absv = false){ TH1D *o = (TH1D*) h->Clone(n); o->SetDirectory(nullptr);
  for(int i = 0; i <= o->GetNbinsX() + 1; i++){ const double v = 100.*h->GetBinContent(i); o->SetBinContent(i, absv ? fabs(v) : v); o->SetBinError(i, 0.); }
  o->GetXaxis()->SetRangeUser(xLo + 1e-3, xHi - 1e-3); return o; }
TH1D* relRatio(TH1D *a, TH1D *b, const char *n){ TH1D *o = (TH1D*) a->Clone(n); o->SetDirectory(nullptr);
  for(int i = 1; i <= o->GetNbinsX(); i++) o->SetBinContent(i, (1. + a->GetBinContent(i))/(1. + b->GetBinContent(i)) - 1.); return o; }
void band(TH1D *assigned){ TH1D *b = pct(assigned, Form("%s_bd", assigned->GetName()));
  for(int i = 1; i <= b->GetNbinsX(); i++){ const double a = fabs(b->GetBinContent(i)); b->SetBinContent(i, 0.); b->SetBinError(i, a); }
  b->SetFillColorAlpha(kGray + 1, 0.35); b->SetLineWidth(0); b->SetMarkerSize(0); b->Draw("E2 SAME"); }
void zeroLine(){ TLine l; l.SetLineColor(kGray + 2); l.SetLineWidth(2); l.DrawLine(xLo, 0., xHi, 0.); }
// check: assigned vs max |lines| (correlated) or the quadrature of two magnitudes (uncorrelated)
void checkMax(TH1D *assigned, const std::vector<TH1D*> &lines){
  for(int i = 1; i <= assigned->GetNbinsX(); i++){ if(!inRange(assigned, i)) continue; double m = 0.; for(TH1D *h : lines) m = std::max(m, fabs(h->GetBinContent(i)));
    maxDiff = std::max(maxDiff, fabs(fabs(assigned->GetBinContent(i)) - m)); } }
void checkQuad(TH1D *assigned, TH1D *a, TH1D *b){
  for(int i = 1; i <= assigned->GetNbinsX(); i++){ if(!inRange(assigned, i)) continue;
    maxDiff = std::max(maxDiff, fabs(fabs(assigned->GetBinContent(i)) - std::hypot(a->GetBinContent(i), b ? b->GetBinContent(i) : 0.))); } }
double yRange(const std::vector<TH1D*> &hs){
  double m = 0.;
  for(TH1D *h : hs) if(h) for(int i = 1; i <= h->GetNbinsX(); i++) if(inRange(h, i)) m = std::max(m, 100.*fabs(h->GetBinContent(i)));
  m *= 1.6; return m > 20. ? 5.*ceil(m/5.) : (m > 8. ? 2.*ceil(m/2.) : std::max(2., ceil(m)));
}
void rowLabel(TCanvas *cv, GridPads &g, int iy, const char *s, double x){
  cv->cd(); TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.027); t.SetTextAngle(90); t.SetTextAlign(21);
  t.DrawLatex(x, 1. - g.tm - (iy + 0.5)*(1. - g.bm - g.tm)/g.ny, s); }
// line styles of the variations
const char *hexGS[2] = {"#0072B2", "#D55E00"};                      // x 0.94, x 1.41
const char *hexC[4]  = {"#0072B2", "#56B4E9", "#E69F00", "#D55E00"}; // x 0.8, 0.9, 1.1, 1.2
const int styC[4] = {1, 2, 2, 1};
const char *labGS[2] = {"#times 0.94", "#times 1.41"}, *labC[4] = {"#times 0.8", "#times 0.9", "#times 1.1", "#times 1.2"};
void legendRow(GridPads &g, int ix, int iy, int n, const char **lab, const char **hex, const int *sty, bool withBand, int nCol = 1){
  g.at(ix, iy)->cd();
  const double top = g.frameTop(iy) - (14. + 1.6*pxText)/padHpx(), rowH = 1.35*pxLeg/padHpx();
  const int nE = n + (withBand ? 1 : 0), nRow = (nE + nCol - 1)/nCol;
  TLegend *lg = makeLegend(g.frameLeft(ix) + 20./padWpx(), top - nRow*rowH, 0.98, top); lg->SetTextFont(43); lg->SetTextSize(pxLeg); lg->SetNColumns(nCol);
  static int nL = 0;
  for(int k = 0; k < n; k++){ TH1D *p = new TH1D(Form("lgP%d", nL++), "", 1, 0, 1); p->SetDirectory(nullptr); styleLine(p, hex[k], 4); p->SetLineStyle(sty ? sty[k] : 1); lg->AddEntry(p, lab[k], "l"); }
  if(withBand){ TH1D *p = new TH1D(Form("lgP%d", nL++), "", 1, 0, 1); p->SetDirectory(nullptr); p->SetFillColorAlpha(kGray + 1, 0.35); p->SetLineWidth(0); lg->AddEntry(p, "assigned", "f"); }
  lg->Draw();
}
}

void plotBJetPuritySys_PF(TString bTag = "_noFlatJEU_fakeGeoOnly_csvV2_2D7bins_nIterOpt_noRecoFloor_normSep_noStitch",
                          TString purTag = "_csvV2_2D7bins", TString outSfx = "_csvV2_2D7bins_nIterOpt"){
  using namespace bpur;
  initPlotStyle();
  const char *repo = "/home/clayton/Analysis/code/bJetRaaAnalysis";
  TFile *fB = TFile::Open(Form("%s/rootFiles/systematics/bJetsPerZSys_PF_MBcaloMatched_JEUppCorr%s.root", repo, bTag.Data()));
  TFile *fPp = TFile::Open(Form("%s/rootFiles/systematics/bPuritySys_total_PFJets_pp%s.root", repo, purTag.Data()));
  TFile *fPa = TFile::Open(Form("%s/rootFiles/systematics/bPuritySys_PFJets_PbPb%s.root", repo, purTag.Data()));
  if(!fB || !fPp || !fPa || fB->IsZombie() || fPp->IsZombie() || fPa->IsZombie()){ printf("ERROR: input missing\n"); return; }
  auto get = [&](TFile *f, TString n) -> TH1D* { TH1D *h = nullptr; f->GetObject(n, h); if(!h) printf("ERROR: %s missing in %s\n", n.Data(), f->GetName()); return h; };
  if(!get(fB, "h_bPur_bGS_0_C1")){ printf("  -> rerun bJetsPerZSys_PF_PbPb.C (2026-11-17 or later saves the purity variations)\n"); return; }
  const char *sys[5] = {"pp", "C4", "C3", "C2", "C1"}, *name[5] = {"pp", "50-80%", "30-50%", "10-30%", "0-10%"};
  const char *src[2] = {"bGS", "cMult"}; const int nv[2] = {2, 4};
  const TString dirOut = Form("%s/figures/systematics/bPurity/", repo); gSystem->mkdir(dirOut, true);
  auto pf = [&](int s, TString n) -> TH1D* { return s == 0 ? get(fPp, n) : get(fPa, n + "_" + sys[s]); };
  const char *rowLab[3] = {"g#rightarrowb#bar{b} fraction", "c fraction", "other (JES, JER, #mu-jet pairs)"};
  const char *otherKey[3] = {"JEU", "JER", "fakeMuJet"};

  // ---- figure 1: the purity
  {
    std::vector<TH1D*> all;
    TH1D *V[5][2][4] = {{{nullptr}}}, *O[5][3] = {{nullptr}}, *Bd[5][3] = {{nullptr}};
    for(int s = 0; s < 5; s++){
      for(int c = 0; c < 2; c++){ for(int k = 0; k < nv[c]; k++){ V[s][c][k] = pf(s, Form("h_varRel_%s_%d", src[c], k)); all.push_back(V[s][c][k]); }
        Bd[s][c] = pf(s, Form("h_sysRel_%s", src[c])); checkMax(Bd[s][c], std::vector<TH1D*>(V[s][c], V[s][c] + nv[c])); }
      for(int o = 0; o < 3; o++){ O[s][o] = (s == 0 && o == 2) ? nullptr : pf(s, Form("h_sysRel_%s", otherKey[o])); if(O[s][o]) all.push_back(O[s][o]); }
      Bd[s][2] = pf(s, "h_sysRel_uncorr");
      { TH1D *q = (TH1D*) Bd[s][2]->Clone("q"); q->Reset(); for(int i = 1; i <= q->GetNbinsX(); i++){ double t = 0.; for(int o = 0; o < 3; o++) if(O[s][o]) t += pow(O[s][o]->GetBinContent(i), 2); q->SetBinContent(i, sqrt(t)); }
        checkQuad(Bd[s][2], q, nullptr); delete q; }
    }
    const double y = yRange(all);
    TCanvas *cv = new TCanvas("cvPur", "", 2200, 1350); setTextScale(cv);
    GridPads g = makeGridPads(cv, 5, 3, 0.095, 0.012, 0.09, 0.055);
    for(int r = 0; r < 3; r++) for(int s = 0; s < 5; s++){
      g.at(s, r)->cd();
      TH1F *f = gPad->DrawFrame(xLo, r < 2 ? -y : 0., xHi, y); g.styleFrame(f, s, r, 0.075);
      if(r != 2) f->GetXaxis()->SetLabelSize(0.);
      sameLabels(f, 30.*cv->GetWw()/1800.);
      if(r < 2){ band(Bd[s][r]); zeroLine();
        for(int k = 0; k < nv[r]; k++){ TH1D *h = pct(V[s][r][k], Form("pv%d_%d_%d", r, s, k)); styleLine(h, r == 0 ? hexGS[k] : hexC[k], 4); h->SetLineStyle(r == 0 ? 1 : styC[k]); h->Draw("HIST SAME"); } }
      else { TH1D *b = pct(Bd[s][2], Form("pb_%d", s)); b->SetFillColorAlpha(kGray + 1, 0.35); b->SetLineWidth(0); b->Draw("HIST SAME");
        for(int o = 0; o < 3; o++) if(O[s][o]){ const SysStyle st = sysStyle(otherKey[o]); TH1D *h = pct(O[s][o], Form("po%d_%d", s, o), true); styleLine(h, st.hex, 4); h->SetLineStyle(st.lineStyle); h->Draw("HIST SAME"); } }
      gPad->RedrawAxis();
      if(r == 0) panelLabel(g, s, r, name[s]);
    }
    for(int r = 0; r < 3; r++) rowLabel(cv, g, r, rowLab[r], 0.042);
    legendRow(g, 0, 0, 2, labGS, hexGS, nullptr, true, 2);
    legendRow(g, 0, 1, 4, labC, hexC, styC, false, 2);
    { const char *lo[3] = {"JES", "JER", "uncorr. #mu-jet pairs"}; const char *ho[3]; int so[3]; for(int o = 0; o < 3; o++){ const SysStyle st = sysStyle(otherKey[o]); ho[o] = st.hex; so[o] = st.lineStyle; }
      legendRow(g, 1, 2, 3, lo, ho, so, true); }
    g.titles(cv, "#it{p}_{T}^{jet} (GeV)", "", 0.036, 0.012, 0.018);
    cv->cd(); TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.027); t.SetTextAngle(90); t.SetTextAlign(21);
    t.DrawLatex(0.012, g.bm + 0.5*(1. - g.bm - g.tm), "#it{P}_{var} / #it{P}_{nom} #minus 1 (%)");
    g.cmsLabel(cv, "Internal", "#it{b} purity (muon #it{p}_{T}^{rel} #times CSVv2 fit),   #sqrt{#it{s}_{NN}} = 5.02 TeV", 0.95);
    savePdfTight(cv, dirOut + "bJetPuritySys_purity_PF" + outSfx + ".pdf");
  }

  // ---- figure 2: b jets per Z
  TH1D *D[5][2][4] = {{{nullptr}}}, *U[5] = {nullptr};
  for(int s = 0; s < 5; s++){ for(int c = 0; c < 2; c++) for(int k = 0; k < nv[c]; k++) D[s][c][k] = get(fB, Form("h_bPur_%s_%d_%s", src[c], k, sys[s])); U[s] = get(fB, Form("h_bPur_uncorr_%s", sys[s])); }
  {
    std::vector<TH1D*> all; TH1D *A[5][3];
    for(int s = 0; s < 5; s++){
      for(int c = 0; c < 2; c++){ for(int k = 0; k < nv[c]; k++) all.push_back(D[s][c][k]); A[s][c] = get(fB, Form("h_sysRel_bPurity_%s_b_%s", src[c], sys[s])); checkMax(A[s][c], std::vector<TH1D*>(D[s][c], D[s][c] + nv[c])); }
      all.push_back(U[s]); A[s][2] = get(fB, Form("h_sysRel_bPurity_b_%s", sys[s])); checkQuad(A[s][2], U[s], nullptr);
    }
    const double y = yRange(all);
    TCanvas *cv = new TCanvas("cvPurB", "", 2200, 1350); setTextScale(cv);
    GridPads g = makeGridPads(cv, 5, 3, 0.095, 0.012, 0.09, 0.055);
    for(int r = 0; r < 3; r++) for(int s = 0; s < 5; s++){
      g.at(s, r)->cd();
      TH1F *f = gPad->DrawFrame(xLo, r < 2 ? -y : 0., xHi, y); g.styleFrame(f, s, r, 0.075);
      if(r != 2) f->GetXaxis()->SetLabelSize(0.);
      sameLabels(f, 30.*cv->GetWw()/1800.);
      if(r < 2){ band(A[s][r]); zeroLine();
        for(int k = 0; k < nv[r]; k++){ TH1D *h = pct(D[s][r][k], Form("dv%d_%d_%d", r, s, k)); styleLine(h, r == 0 ? hexGS[k] : hexC[k], 4); h->SetLineStyle(r == 0 ? 1 : styC[k]); h->Draw("HIST SAME"); } }
      else { TH1D *b = pct(A[s][2], Form("db_%d", s), true); b->SetFillColorAlpha(kGray + 1, 0.35); b->SetLineWidth(0); b->Draw("HIST SAME");
        TH1D *rp = pct(pf(s, "h_sysRel_uncorr"), Form("dr_%d", s), true); styleLine(rp, "#000000", 4); rp->SetLineStyle(2); rp->Draw("HIST SAME");
        TH1D *h = pct(U[s], Form("du_%d", s), true); styleLine(h, "#CC79A7", 5); h->Draw("HIST SAME"); }
      gPad->RedrawAxis();
      if(r == 0) panelLabel(g, s, r, name[s]);
    }
    for(int r = 0; r < 3; r++) rowLabel(cv, g, r, rowLab[r], 0.042);
    legendRow(g, 0, 0, 2, labGS, hexGS, nullptr, true, 2);
    legendRow(g, 0, 1, 4, labC, hexC, styC, false, 2);
    { const char *lo[2] = {"purity", "unfolded"}; const char *ho[2] = {"#000000", "#CC79A7"}; int so[2] = {2, 1}; legendRow(g, 1, 2, 2, lo, ho, so, true); }
    g.titles(cv, "#it{p}_{T}^{jet} (GeV)", "", 0.036, 0.012, 0.018);
    cv->cd(); TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.027); t.SetTextAngle(90); t.SetTextAlign(21);
    t.DrawLatex(0.012, g.bm + 0.5*(1. - g.bm - g.tm), "#it{b} jets per #it{N}_{Z}: varied / nominal #minus 1 (%)");
    g.cmsLabel(cv, "Internal", "#it{b} purity, #it{b} jets per #it{N}_{Z} (unfolded),   #sqrt{#it{s}_{NN}} = 5.02 TeV", 0.95);
    savePdfTight(cv, dirOut + "bJetPuritySys_perSystem_PF" + outSfx + ".pdf");
  }

  // ---- figure 3: ratios, R_Z(b) per class and R_CP^Z(b) per class (columns)
  {
    const int NC = 7; const int sA[NC] = {1, 2, 3, 4, 2, 3, 4}, sR[NC] = {0, 0, 0, 0, 1, 1, 1};
    const char *cls[5] = {"", "C4", "C3", "C2", "C1"};
    std::vector<TH1D*> all; TH1D *L[NC][2][4] = {{{nullptr}}}, *A[NC][3];
    for(int c = 0; c < NC; c++){
      const TString q = sR[c] == 0 ? Form("b_PbPbOverPP_%s", cls[sA[c]]) : Form("b_RCP_%s", cls[sA[c]]);
      for(int r = 0; r < 2; r++){ for(int k = 0; k < nv[r]; k++){ L[c][r][k] = relRatio(D[sA[c]][r][k], D[sR[c]][r][k], Form("rl%d_%d_%d", c, r, k)); all.push_back(L[c][r][k]); }
        A[c][r] = get(fB, Form("h_sysRel_bPurity_%s_%s", src[r], q.Data())); checkMax(A[c][r], std::vector<TH1D*>(L[c][r], L[c][r] + nv[r])); }
      A[c][2] = get(fB, "h_sysRel_bPurity_" + q); checkQuad(A[c][2], U[sA[c]], U[sR[c]]); all.push_back(A[c][2]);
    }
    // two figures, R_Z(b) (4 columns) and R_CP^Z(b) (3 columns), each with its own common y range
    auto drawRatios = [&](int c0, int c1, const char *fig, const char *what){
      std::vector<TH1D*> sub; for(int c = c0; c < c1; c++){ for(int r = 0; r < 2; r++) for(int k = 0; k < nv[r]; k++) sub.push_back(L[c][r][k]); sub.push_back(A[c][2]); }
      const double y = yRange(sub);
      const int nc = c1 - c0;
      TCanvas *cv = new TCanvas(Form("cvPurR_%s", fig), "", 1800, 1350); setTextScale(cv);
      GridPads g = makeGridPads(cv, nc, 3, 0.115, 0.012, 0.09, 0.055);
      for(int r = 0; r < 3; r++) for(int j = 0; j < nc; j++){
        const int c = c0 + j;
        g.at(j, r)->cd();
        TH1F *f = gPad->DrawFrame(xLo, r < 2 ? -y : 0., xHi, y); g.styleFrame(f, j, r, 0.075);
        if(r != 2) f->GetXaxis()->SetLabelSize(0.);
        sameLabels(f, 30.);
        if(r < 2){ band(A[c][r]); zeroLine();
          for(int k = 0; k < nv[r]; k++){ TH1D *h = pct(L[c][r][k], Form("lv%d_%d_%d", r, c, k)); styleLine(h, r == 0 ? hexGS[k] : hexC[k], 4); h->SetLineStyle(r == 0 ? 1 : styC[k]); h->Draw("HIST SAME"); } }
        else { TH1D *b = pct(A[c][2], Form("rb_%d", c), true); b->SetFillColorAlpha(kGray + 1, 0.35); b->SetLineWidth(0); b->Draw("HIST SAME");
          TH1D *aa = pct(U[sA[c]], Form("ra_%d", c), true), *rr = pct(U[sR[c]], Form("rr_%d", c), true);
          styleLine(rr, "#000000", 4); rr->SetLineStyle(2); rr->Draw("HIST SAME"); styleLine(aa, "#CC79A7", 5); aa->Draw("HIST SAME"); }
        gPad->RedrawAxis();
        if(r == 0) panelLabel(g, j, r, sR[c] == 0 ? Form("%s / pp", name[sA[c]]) : Form("%s / 50-80%%", name[sA[c]]));
      }
      for(int r = 0; r < 3; r++) rowLabel(cv, g, r, rowLab[r], 0.058);
      legendRow(g, 0, 0, 2, labGS, hexGS, nullptr, true, 1);
      legendRow(g, 0, 1, 4, labC, hexC, styC, false, 2);
      { const char *lo[2] = {"class", "reference"}; const char *ho[2] = {"#CC79A7", "#000000"}; int so[2] = {1, 2}; legendRow(g, 0, 2, 2, lo, ho, so, true); }
      g.titles(cv, "#it{p}_{T}^{jet} (GeV)", "", 0.036, 0.012, 0.018);
      cv->cd(); TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.027); t.SetTextAngle(90); t.SetTextAlign(21);
      t.DrawLatex(0.028, g.bm + 0.5*(1. - g.bm - g.tm), "varied / nominal #minus 1 (%)");
      g.cmsLabel(cv, "Internal", Form("#it{b} purity on %s,   #sqrt{#it{s}_{NN}} = 5.02 TeV", what), 0.95);
      savePdfTight(cv, dirOut + Form("bJetPuritySys_%s_PF", fig) + outSfx + ".pdf");
    };
    drawRatios(0, 4, "RZ", "#it{R}_{Z}(#it{b})");
    drawRatios(4, 7, "RCP", "#it{R}_{CP}^{Z}(#it{b})");
  }
  printf("largest |assigned - (max |variation| | quadrature)| over every band drawn, 80-200 GeV: %.2e\n", maxDiff);
  // printout: assigned values (%), 80-200 GeV
  auto row = [&](const char *lab, TH1D *h){ printf("  %-36s", lab); for(int i = 1; i <= h->GetNbinsX(); i++) if(inRange(h, i)) printf(" %5.1f", 100.*fabs(h->GetBinContent(i))); printf("\n"); };
  printf("\nassigned purity sources (%%), 80-200 GeV\n");
  const char *grp[3] = {"bPurity_bGS", "bPurity_cMult", "bPurity"};
  for(int gI = 0; gI < 3; gI++){
    for(int s = 0; s < 5; s++) row(Form("%s b per Z %s", grp[gI], name[s]), get(fB, Form("h_sysRel_%s_b_%s", grp[gI], sys[s])));
    for(int c = 4; c >= 1; c--) row(Form("%s R_Z(b) C%d", grp[gI], c), get(fB, Form("h_sysRel_%s_b_PbPbOverPP_C%d", grp[gI], c)));
    for(int c = 3; c >= 1; c--) row(Form("%s R_CP(b) C%d", grp[gI], c), get(fB, Form("h_sysRel_%s_b_RCP_C%d", grp[gI], c)));
  }
  printf("\npurity level (%%), 80-200 GeV: g->bb | c | uncorr\n");
  for(int s = 0; s < 5; s++){ row(Form("bGS %s", name[s]), pf(s, "h_sysRel_bGS")); row(Form("cMult %s", name[s]), pf(s, "h_sysRel_cMult")); row(Form("uncorr %s", name[s]), pf(s, "h_sysRel_uncorr")); }
}
