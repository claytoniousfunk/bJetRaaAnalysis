// Extraction of the b-jet JET ENERGY RESOLUTION systematic, from the variation to every assigned value (2026-11-17,
// user: figures that justify every systematic value, one source at a time; JER after JES).
//
// How the source is made (jerSmear.h; bJetsPerZSys_PF_PbPb.C and unfoldBJetSpectrum_caloJets_pp.C, PFJets):
//   the b response is NOMINALLY smeared to a 10% worse resolution than the simulation (CMS-HIN-18-014), and the
//   variation smears it by a further 10% (1.21 x the simulation): each generator-level column gets an extra Gaussian
//   of width 0.458 x its own sigma(pT^reco / pT^gen). The variation is measured on the SIMULATION: the nominal
//   response's reco spectrum unfolded with the nominal and with the varied response at the nominal N, var / nom - 1,
//   applied to the data (BJET_JER_INPUT=data = the same ratio on the data, drawn for comparison). One variation,
//   correlated between pp and PbPb, between the classes and between b and inclusive jets (the inclusive side is the
//   inclusive chain's own signed JER shift on the b-jet axis): every ratio is formed for the variation, assigned =
//   |ratio - 1|.
//
// Figure 1 (resolution): sigma(pT^reco / pT^gen) of the b response per gen-pT column, simulation, nominal (x 1.10),
//   variation (x 1.21), pp and the four classes.
// Figure 2 (per system, common y): unfolded b jets per Z, var / nom - 1 measured on the simulation (assigned, band) and
//   on the data (comparison); the inclusive shift; b / inclusive (assigned band).
// Figure 3 (ratios, common y): R_Z(b), R_Z(b) / R_Z(incl), R_CP^Z(b), R_CP^Z(b) / R_CP^Z(incl), with the band.
// The macro checks that every band equals |the line it is drawn with| and prints the largest difference.
//
// inputs  rootFiles/systematics/bJetsPerZSys_PF_MBcaloMatched_JEUppCorr<bTag>.root (h_bJER_*, saved since
//         2026-11-17, h_sysRel_JER_*), the inclusive file of the same run (h_varRel_JER_<sys>_bjet), and the pp
//         unfolded file (h_unf_varRel_JERresp[_data], h_bJER_res*)
// Usage, from src/systematics/:  root -l -b -q plotBJetJERSys_PF.C
// Output: figures/systematics/JER/bJetJERSys_{resolution,perSystem,ratios}_PF<outSfx>.pdf
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

namespace bjer {
const double xLo = 80., xHi = 200.;
const char *hexJMC = "#D55E00", *hexJData = "#0072B2";
TFile *fNomFile = nullptr;   // the results file, for the populated bins of each quantity
double maxDiff = 0.; int nFlat = 0, nBin = 0;   // bands flattened to their mean / kept bin by bin
double pxText = 34., pxLeg = 30.;
void setTextScale(TCanvas *cv){ const double k = cv->GetWw()/1800.; pxText = 34.*k; pxLeg = 30.*k; }
double padHpx(){ return gPad->GetWh()*gPad->GetAbsHNDC(); }
double padWpx(){ return gPad->GetWw()*gPad->GetAbsWNDC(); }
void ptext(double x, double y, const char *s, double px, int font = 43, int align = 13, Color_t col = kBlack){
  TLatex t; t.SetNDC(); t.SetTextFont(font); t.SetTextSize(px); t.SetTextAlign(align); t.SetTextColor(col); t.DrawLatex(x, y, s); }
void panelLabel(GridPads &g, int ix, int iy, const char *s){ ptext(g.frameLeft(ix) + 20./padWpx(), g.frameTop(iy) - 14./padHpx(), s, pxText, 63, 13); }
double minDim(){ return std::min(padWpx(), padHpx()); }
void sameLabels(TH1F *f, double px){ if(f->GetXaxis()->GetLabelSize() > 0.) f->GetXaxis()->SetLabelSize(px/minDim()); if(f->GetYaxis()->GetLabelSize() > 0.) f->GetYaxis()->SetLabelSize(px/minDim()); }
TH1D* pct(TH1D *h, const char *n){ TH1D *o = (TH1D*) h->Clone(n); o->SetDirectory(nullptr);
  for(int i = 0; i <= o->GetNbinsX() + 1; i++){ o->SetBinContent(i, 100.*h->GetBinContent(i)); o->SetBinError(i, 0.); } o->GetXaxis()->SetRangeUser(xLo + 1e-3, xHi - 1e-3); return o; }
TH1D* relRatio(TH1D *a, TH1D *b, const char *n){ TH1D *o = (TH1D*) a->Clone(n); o->SetDirectory(nullptr);
  for(int i = 1; i <= o->GetNbinsX(); i++) o->SetBinContent(i, (1. + a->GetBinContent(i))/(1. + b->GetBinContent(i)) - 1.); return o; }
// one panel: zero line, assigned band (+-), the variation (and, optionally, the data-based comparison, dashed)
void panel(TH1D *v, TH1D *assigned, TH1D *cmp = nullptr){
  if(assigned){
    TH1D *b = pct(assigned, Form("%s_band", assigned->GetName()));
    for(int i = 1; i <= b->GetNbinsX(); i++){ const double a = fabs(b->GetBinContent(i)); b->SetBinContent(i, 0.); b->SetBinError(i, a); }
    b->SetFillColorAlpha(kGray + 1, 0.35); b->SetLineWidth(0); b->SetMarkerSize(0); b->Draw("E2 SAME");
    // the source is FLATTENED to its mean over the reported bins where it shows no pT trend (bJetsPerZSys_PF_PbPb.C,
    // flattenSystematic.h): a constant band is checked against the mean |variation|, a varying one bin by bin
    double amin = 1e9, amax = -1e9, msum = 0.; int n = 0;
    for(int i = 1; i <= assigned->GetNbinsX(); i++){ const double lo = assigned->GetXaxis()->GetBinLowEdge(i); if(lo < xLo - 1e-6 || lo >= xHi - 1e-6) continue;
      amin = std::min(amin, fabs(assigned->GetBinContent(i))); amax = std::max(amax, fabs(assigned->GetBinContent(i))); }
    // the flattening (flattenIfNoTrend) averages every bin from 80 GeV up in which the quantity itself is populated: the
    // 200-300 GeV overflow always, 300-500 GeV only where the measured spectrum reaches it (h_nom_<quantity> > 0)
    TString q = assigned->GetName(); q.ReplaceAll("h_sysRel_JER_", "h_nom_"); TH1D *nom = nullptr; if(fNomFile) fNomFile->GetObject(q, nom);
    for(int i = 1; i <= v->GetNbinsX(); i++){ const double lo = v->GetXaxis()->GetBinLowEdge(i); if(lo < xLo - 1e-6) continue;
      if(nom ? nom->GetBinContent(i) <= 0. : lo >= 300. - 1e-6) continue;
      msum += fabs(v->GetBinContent(i)); n++; }
    const bool flat = amax - amin < 1e-9; if(flat) nFlat++; else nBin++;
    for(int i = 1; i <= assigned->GetNbinsX(); i++){ const double lo = assigned->GetXaxis()->GetBinLowEdge(i); if(lo < xLo - 1e-6 || lo >= xHi - 1e-6) continue;
      const double dd = fabs(fabs(assigned->GetBinContent(i)) - (flat ? msum/n : fabs(v->GetBinContent(i))));
      maxDiff = std::max(maxDiff, dd); }
  }
  TLine l; l.SetLineColor(kGray + 2); l.SetLineWidth(2); l.DrawLine(xLo, 0., xHi, 0.);
  TH1D *u = pct(v, Form("%s_pv", v->GetName())); styleLine(u, hexJMC, 5); u->Draw("HIST SAME");
  if(cmp){ TH1D *c = pct(cmp, Form("%s_pc", cmp->GetName())); styleLine(c, hexJData, 4); c->SetLineStyle(2); c->Draw("HIST SAME"); }   // on top: visible where the two agree
}
double yRange(const std::vector<TH1D*> &hs){
  double m = 0.;
  for(TH1D *h : hs) if(h) for(int i = 1; i <= h->GetNbinsX(); i++){ const double lo = h->GetXaxis()->GetBinLowEdge(i); if(lo >= xLo - 1e-6 && lo < xHi - 1e-6) m = std::max(m, 100.*fabs(h->GetBinContent(i))); }
  m *= 1.45; return m > 20. ? 5.*ceil(m/5.) : (m > 8. ? 2.*ceil(m/2.) : std::max(2., ceil(m)));
}
void rowLabel(TCanvas *cv, GridPads &g, int iy, const char *s, double x){
  cv->cd(); TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.027); t.SetTextAngle(90); t.SetTextAlign(21);
  t.DrawLatex(x, 1. - g.tm - (iy + 0.5)*(1. - g.bm - g.tm)/g.ny, s); }
}

void plotBJetJERSys_PF(TString bTag = "_noFlatJEU_fakeGeoOnly_csvV2_2D7bins_nIterOpt_noRecoFloor_normSep_noStitch",
                       TString inclTag = "_noFlatJEU_fakeGeoOnly_pt7bins_noRecoFloor_normSep_noStitch",
                       TString ppVar = "_csvV2_2D7bins_nIterOpt", TString outSfx = "_csvV2_2D7bins_nIterOpt"){
  using namespace bjer;
  initPlotStyle();
  const char *repo = "/home/clayton/Analysis/code/bJetRaaAnalysis";
  TFile *fB = TFile::Open(Form("%s/rootFiles/systematics/bJetsPerZSys_PF_MBcaloMatched_JEUppCorr%s.root", repo, bTag.Data()));
  TFile *fI = TFile::Open(Form("%s/rootFiles/systematics/inclJetsPerZSys_PF_MBcaloMatched_JEUppCorr%s.root", repo, inclTag.Data()));
  TFile *fP = TFile::Open(Form("%s/rootFiles/CorrectedBJetSpectra/Data/unfoldedBJetSpectrum_PFJets_pp%s.root", repo, ppVar.Data()));
  if(!fB || !fI || !fP || fB->IsZombie() || fI->IsZombie() || fP->IsZombie()){ printf("ERROR: input missing\n"); return; }
  fNomFile = fB;
  auto get = [&](TFile *f, TString n) -> TH1D* { TH1D *h = nullptr; f->GetObject(n, h); if(!h) printf("ERROR: %s missing in %s\n", n.Data(), f->GetName()); return h; };
  if(!get(fB, "h_bJER_varMC_C1") || !get(fP, "h_bJER_resNom")){ printf("  -> rerun unfoldBJetSpectrum_caloJets_pp.C and bJetsPerZSys_PF_PbPb.C (2026-11-17 or later)\n"); return; }
  const char *sys[5] = {"pp", "C4", "C3", "C2", "C1"}, *name[5] = {"pp", "50-80%", "30-50%", "10-30%", "0-10%"};
  const TString dirOut = Form("%s/figures/systematics/JER/", repo); gSystem->mkdir(dirOut, true);

  TH1D *bMC[5], *bDat[5], *inc[5], *fr[5], *aB[5], *aF[5], *res[3][5];
  for(int s = 0; s < 5; s++){
    bMC[s]  = s == 0 ? get(fP, "h_unf_varRel_JERresp") : get(fB, Form("h_bJER_varMC_%s", sys[s]));
    bDat[s] = s == 0 ? get(fP, "h_unf_varRel_JERresp_data") : get(fB, Form("h_bJER_varData_%s", sys[s]));
    inc[s]  = get(fI, Form("h_varRel_JER_%s_bjet", sys[s]));
    for(int k = 0; k < 3; k++){ const char *w[3] = {"Sim", "Nom", "Var"}; res[k][s] = s == 0 ? get(fP, Form("h_bJER_res%s", w[k])) : get(fB, Form("h_bJER_res%s_%s", w[k], sys[s])); }
    aB[s] = get(fB, Form("h_sysRel_JER_b_%s", sys[s])); aF[s] = get(fB, Form("h_sysRel_JER_bOverIncl_%s", sys[s]));
    if(!bMC[s] || !bDat[s] || !inc[s] || !aB[s] || !aF[s] || !res[0][s] || !res[1][s] || !res[2][s]) return;
    fr[s] = relRatio(bMC[s], inc[s], Form("fr_%s", sys[s]));
  }

  // ---- figure 1: the resolution
  {
    TCanvas *cv = new TCanvas("cvJERres", "", 2200, 700); setTextScale(cv);
    GridPads g = makeGridPads(cv, 5, 1, 0.075, 0.012, 0.17, 0.10);
    const double gLo = 60., gHi = 300.; double yMax = 0.;
    for(int s = 0; s < 5; s++) for(int k = 0; k < 3; k++) for(int i = 1; i <= res[k][s]->GetNbinsX(); i++){ const double x = res[k][s]->GetBinCenter(i); if(x > gLo && x < gHi) yMax = std::max(yMax, res[k][s]->GetBinContent(i)); }
    yMax = 0.05*ceil(1.35*yMax/0.05);
    const char *hexR[3] = {"#000000", hexJData, hexJMC}; const int styR[3] = {3, 1, 7};
    for(int s = 0; s < 5; s++){
      g.at(s, 0)->cd();
      TH1F *f = gPad->DrawFrame(gLo, 0., gHi, yMax); g.styleFrame(f, s, 0, 0.06); sameLabels(f, 30.*cv->GetWw()/1800.);
      for(int k = 0; k < 3; k++){ TH1D *h = (TH1D*) res[k][s]->Clone(Form("rs%d_%d", k, s)); h->SetDirectory(nullptr); h->GetXaxis()->SetRangeUser(gLo + 1e-3, gHi - 1e-3);
        styleLine(h, hexR[k], 4); h->SetLineStyle(styR[k]); h->Draw("HIST SAME"); }
      gPad->RedrawAxis(); panelLabel(g, s, 0, name[s]);
    }
    g.at(0, 0)->cd();
    { const double top = g.frameTop(0) - (14. + 1.6*pxText)/padHpx(), rowH = 1.4*pxLeg/padHpx();
      TLegend *lg = makeLegend(g.frameLeft(0) + 20./padWpx(), top - 3*rowH, 0.98, top); lg->SetTextFont(43); lg->SetTextSize(pxLeg);
      const char *lab[3] = {"simulation", "nominal (#times 1.10)", "variation (#times 1.21)"};
      for(int k = 0; k < 3; k++){ TH1D *p = new TH1D(Form("lgR%d", k), "", 1, 0, 1); styleLine(p, hexR[k], 4); p->SetLineStyle(styR[k]); lg->AddEntry(p, lab[k], "l"); }
      lg->Draw(); }
    g.titles(cv, "generator-level #it{p}_{T}^{jet} (GeV)", "#sigma(#it{p}_{T}^{reco} / #it{p}_{T}^{gen})", 0.05, 0.022, 0.025);
    g.cmsLabel(cv, "Simulation Internal", "muon-tagged #it{b}-jet response,   anti-#it{k}_{T} #it{R} = 0.4 PF jets,  #sqrt{#it{s}_{NN}} = 5.02 TeV", 1.2);
    savePdfTight(cv, dirOut + "bJetJERSys_resolution_PF" + outSfx + ".pdf");
  }

  // ---- figure 2: per system
  {
    const int NR = 3;
    const char *rowLab[NR] = {"unfolded #it{b} jets", "unfolded inclusive", "#it{b} / inclusive"};
    TCanvas *cv = new TCanvas("cvJERsys", "", 2200, 1350); setTextScale(cv);
    GridPads g = makeGridPads(cv, 5, NR, 0.095, 0.012, 0.09, 0.055);
    std::vector<TH1D*> all; for(int s = 0; s < 5; s++){ all.push_back(bMC[s]); all.push_back(bDat[s]); all.push_back(inc[s]); all.push_back(fr[s]); }
    const double y = yRange(all);
    for(int r = 0; r < NR; r++) for(int s = 0; s < 5; s++){
      g.at(s, r)->cd();
      TH1F *f = gPad->DrawFrame(xLo, -y, xHi, y); g.styleFrame(f, s, r, 0.075);
      if(r != NR - 1) f->GetXaxis()->SetLabelSize(0.);
      sameLabels(f, 30.*cv->GetWw()/1800.);
      if(r == 0) panel(bMC[s], aB[s], bDat[s]);
      if(r == 1) panel(inc[s], nullptr);
      if(r == 2) panel(fr[s], aF[s]);
      gPad->RedrawAxis();
      if(r == 0) panelLabel(g, s, r, name[s]);
    }
    for(int r = 0; r < NR; r++) rowLabel(cv, g, r, rowLab[r], 0.042);
    g.at(0, 1)->cd();   // legend in the pp inclusive panel: no band, values a few % below zero
    { const double top = g.frameTop(1) - 14./padHpx(), rowH = 1.4*pxLeg/padHpx();
      TLegend *lg = makeLegend(g.frameLeft(0) + 20./padWpx(), top - 3*rowH, 0.98, top); lg->SetTextFont(43); lg->SetTextSize(pxLeg);
      TH1D *a = new TH1D("lgA", "", 1, 0, 1), *b = new TH1D("lgB", "", 1, 0, 1), *c = new TH1D("lgC", "", 1, 0, 1);
      styleLine(a, hexJMC, 5); styleLine(b, hexJData, 4); b->SetLineStyle(2); c->SetFillColorAlpha(kGray + 1, 0.35); c->SetLineWidth(0);
      lg->AddEntry(a, "on simulation", "l"); lg->AddEntry(b, "on data", "l"); lg->AddEntry(c, "assigned", "f"); lg->Draw(); }
    g.titles(cv, "#it{p}_{T}^{jet} (GeV)", "", 0.036, 0.012, 0.018);
    cv->cd(); TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.027); t.SetTextAngle(90); t.SetTextAlign(21);
    t.DrawLatex(0.012, g.bm + 0.5*(1. - g.bm - g.tm), "varied / nominal #minus 1 (%)");
    g.cmsLabel(cv, "Internal", "jet energy resolution, #it{b} jets,   anti-#it{k}_{T} #it{R} = 0.4 PF jets,  #sqrt{#it{s}_{NN}} = 5.02 TeV", 0.95);
    savePdfTight(cv, dirOut + "bJetJERSys_perSystem_PF" + outSfx + ".pdf");
  }

  // ---- figure 3: ratios
  {
    const char *cls[4] = {"C4", "C3", "C2", "C1"}, *cName[4] = {"50-80%", "30-50%", "10-30%", "0-10%"};
    const int NR = 4;
    const char *rowLab[NR] = {"#it{R}_{Z}(#it{b})", "#it{R}_{Z}(#it{b}) / #it{R}_{Z}(incl)", "#it{R}_{CP}^{Z}(#it{b})", "#it{R}_{CP}^{Z}(#it{b}) / #it{R}_{CP}^{Z}(incl)"};
    TH1D *L[NR][4] = {{nullptr}}, *A[NR][4] = {{nullptr}};
    std::vector<TH1D*> all;
    for(int r = 0; r < NR; r++) for(int c = 0; c < 4; c++){
      const int s = c + 1; const bool rcp = r >= 2; if(rcp && c == 0) continue;
      const int ref = rcp ? 1 : 0;
      L[r][c] = r % 2 == 0 ? relRatio(bMC[s], bMC[ref], Form("rr%d_%d", r, c)) : relRatio(fr[s], fr[ref], Form("rr%d_%d", r, c));
      A[r][c] = get(fB, Form("h_sysRel_JER_%s_%s_%s", r % 2 == 0 ? "b" : "bOverIncl", rcp ? "RCP" : "PbPbOverPP", cls[c]));
      all.push_back(L[r][c]);
    }
    const double y = yRange(all);
    TCanvas *cv = new TCanvas("cvJERratio", "", 1800, 1700); setTextScale(cv);
    GridPads g = makeGridPads(cv, 4, NR, 0.135, 0.012, 0.075, 0.045);
    for(int r = 0; r < NR; r++) for(int c = 0; c < 4; c++){
      g.at(c, r)->cd();
      TH1F *f = gPad->DrawFrame(xLo, -y, xHi, y); g.styleFrame(f, c, r, 0.075);
      if(r != NR - 1) f->GetXaxis()->SetLabelSize(0.);
      sameLabels(f, 30.);
      if(L[r][c]) panel(L[r][c], A[r][c]);
      gPad->RedrawAxis();
      if(r == 0 || (r == 2 && c > 0)) panelLabel(g, c, r, r < 2 ? Form("%s / pp", cName[c]) : Form("%s / 50-80%%", cName[c]));
      if(r >= 2 && c == 0) ptext(0.5*(g.frameLeft(0) + 1.), 0.5*(g.frameTop(r) + g.frameBottom(r)), "reference", pxText, 43, 22, kGray + 2);
    }
    for(int r = 0; r < NR; r++) rowLabel(cv, g, r, rowLab[r], 0.060);
    g.at(0, 1)->cd();
    { const double top = g.frameTop(1) - 14./padHpx(), rowH = 1.4*pxLeg/padHpx();
      TLegend *lg = makeLegend(g.frameLeft(0) + 20./padWpx(), top - 2*rowH, 0.98, top); lg->SetTextFont(43); lg->SetTextSize(pxLeg);
      TH1D *a = new TH1D("lgA2", "", 1, 0, 1), *c = new TH1D("lgC2", "", 1, 0, 1); styleLine(a, hexJMC, 5); c->SetFillColorAlpha(kGray + 1, 0.35); c->SetLineWidth(0);
      lg->AddEntry(a, "JER variation", "l"); lg->AddEntry(c, "assigned", "f"); lg->Draw(); }
    g.titles(cv, "#it{p}_{T}^{jet} (GeV)", "", 0.034, 0.012, 0.016);
    cv->cd(); TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.027); t.SetTextAngle(90); t.SetTextAlign(21);
    t.DrawLatex(0.030, g.bm + 0.5*(1. - g.bm - g.tm), "varied / nominal #minus 1 (%)");
    g.cmsLabel(cv, "Internal", "jet energy resolution, ratios,   #sqrt{#it{s}_{NN}} = 5.02 TeV", 0.95);
    savePdfTight(cv, dirOut + "bJetJERSys_ratios_PF" + outSfx + ".pdf");
  }
  printf("largest |assigned - |variation|| (flattened bands: - mean |variation|) over every band drawn, 80-200 GeV: %.2e  (%d bands flattened, %d bin by bin)\n", maxDiff, nFlat, nBin);
  printf("\nassigned JER (%%), reported bins 80-200 GeV\n");
  auto row = [&](const char *lab, TH1D *h){ printf("  %-30s", lab); for(int i = 1; i <= h->GetNbinsX(); i++){ const double lo = h->GetXaxis()->GetBinLowEdge(i); if(lo >= xLo - 1e-6 && lo < xHi - 1e-6) printf(" %5.1f", 100.*fabs(h->GetBinContent(i))); } printf("\n"); };
  for(int s = 0; s < 5; s++) row(Form("b jets per Z, %s", name[s]), aB[s]);
  for(int s = 0; s < 5; s++) row(Form("b / inclusive, %s", name[s]), aF[s]);
  for(int s = 0; s < 5; s++) row(Form("data-based b shift, %s", name[s]), bDat[s]);
  for(int c = 4; c >= 1; c--){ row(Form("R_Z(b), C%d", c), get(fB, Form("h_sysRel_JER_b_PbPbOverPP_C%d", c))); row(Form("double ratio, C%d", c), get(fB, Form("h_sysRel_JER_bOverIncl_PbPbOverPP_C%d", c))); }
  for(int c = 3; c >= 1; c--){ row(Form("R_CP(b), C%d", c), get(fB, Form("h_sysRel_JER_b_RCP_C%d", c))); row(Form("R_CP(b)/R_CP(incl), C%d", c), get(fB, Form("h_sysRel_JER_bOverIncl_RCP_C%d", c))); }
  printf("\nresolution sigma(reco/gen) at gen 85 / 105 / 175 GeV: simulation | nominal | variation\n");
  for(int s = 0; s < 5; s++){ printf("  %-7s", name[s]); for(double x : {85., 105., 175.}){ printf("   %.3f %.3f %.3f", res[0][s]->GetBinContent(res[0][s]->FindBin(x)), res[1][s]->GetBinContent(res[1][s]->FindBin(x)), res[2][s]->GetBinContent(res[2][s]->FindBin(x))); } printf("\n"); }
}
