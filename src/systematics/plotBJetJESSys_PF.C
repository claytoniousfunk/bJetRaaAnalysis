// Extraction of the b-jet JET ENERGY SCALE systematic, every step from the variation to the assigned value
// (2026-11-16, user: "figures that justify every single value that goes into the systematic"; b-jet JES first).
//
// How the source is made (bJetsPerZSys_PF_PbPb.C, PFJets mode, and unfoldBJetSpectrum_caloJets_pp.C):
//   the MEASURED (reco-level, purity-corrected) b-jet spectrum is shifted by +/- delta(pT), the AK4PF JEC
//   uncertainty averaged over |eta| < 1.6 (pp: Spring18_ppRef5TeV_V6, PbPb: Autumn18_HI_V8; headers/functions/
//   jeuShift.h), unfolded with the b response at the nominal N, and compared with the nominal unfolded spectrum.
//   The inclusive side gets the same analytic shift on its own measured spectrum (h_varRel_JEUdata_*_bjet of the
//   inclusive file). Correlated between pp and PbPb, between the PbPb classes and between b and inclusive: every
//   ratio is formed per variation, assigned = max over up / down of |ratio - 1|. Not flattened (noFlatJEU).
//
// Figure 1 (per system, pp and the four classes; columns), rows:
//   measured b spectrum shifted / nominal - 1 (reco level)   pp recomputed here with the pp chain's own call
//   unfolded b jets per Z: shifted / nominal - 1, with the assigned band (h_sysRel_JEU_b_<sys>)
//   inclusive jets, the same shift method (b-jet axis)       no band: enters only through b / inclusive
//   b / inclusive (b-jet fraction), with the assigned band (h_sysRel_JEU_bOverIncl_<sys>)
// Figure 2 (ratios between systems): R_Z(b) and R_Z(b) / R_Z(incl) per class; R_CP(b) and its b / inclusive
//   version over 50-80%. Lines = the ratio per variation, band = assigned.
// The macro checks that every band equals max(|up|, |down|) of the lines it is drawn with and prints the largest
// difference (must be 0 up to float precision).
//
// inputs  rootFiles/systematics/bJetsPerZSys_PF_MBcaloMatched_JEUppCorr<bTag>.root   (h_bJEU_*, h_varRel_JEU_*,
//         saved since 2026-11-16, and h_sysRel_JEU_*)
//         rootFiles/MuTagMuTrigBJetSpectra/Data/bJetSpectrum_PFJets_pp<var>.root     (pp measured, reco shift)
// Usage, from src/systematics/:  root -l -b -q plotBJetJESSys_PF.C
// Output: figures/systematics/JEU/bJetJESSys_perSystem_PF.pdf, bJetJESSys_ratios_PF.pdf
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
#include "../../headers/functions/jeuShift.h"

namespace bjes {
const double xLo = 80., xHi = 200.;
const char *hexUp = "#D55E00", *hexDn = "#0072B2";
double maxDiff = 0.;   // assigned vs max(|up|, |down|)

// labels the same printed size in every pad: ROOT scales a label by the pad's SMALLER dimension in pixels
double minDim(TVirtualPad *p){ return std::min(p->GetWw()*p->GetAbsWNDC(), p->GetWh()*p->GetAbsHNDC()); }
void sameLabels(TH1F *f, double px){ const double m = minDim(gPad);
  if(f->GetXaxis()->GetLabelSize() > 0.) f->GetXaxis()->SetLabelSize(px/m);
  if(f->GetYaxis()->GetLabelSize() > 0.) f->GetYaxis()->SetLabelSize(px/m); }
// relative shift (in %) as a histogram on the given histogram's axis
TH1D* pct(TH1D *h, const char *name){ TH1D *o = (TH1D*) h->Clone(name); o->SetDirectory(nullptr);
  for(int i = 0; i <= o->GetNbinsX() + 1; i++){ o->SetBinContent(i, 100.*h->GetBinContent(i)); o->SetBinError(i, 0.); }
  o->GetXaxis()->SetRangeUser(xLo + 1e-3, xHi - 1e-3); return o; }   // drawn over the reported range only (no edge spikes)
// var / nom - 1
TH1D* rel(TH1D *var, TH1D *nom, const char *name){ TH1D *o = (TH1D*) nom->Clone(name); o->SetDirectory(nullptr); o->Reset();
  for(int i = 1; i <= o->GetNbinsX(); i++){ const double n = nom->GetBinContent(i); o->SetBinContent(i, n > 0. ? var->GetBinContent(i)/n - 1. : 0.); } return o; }
// (1 + a) / (1 + b) - 1, per bin
TH1D* relRatio(TH1D *a, TH1D *b, const char *name){ TH1D *o = (TH1D*) a->Clone(name); o->SetDirectory(nullptr);
  for(int i = 1; i <= o->GetNbinsX(); i++) o->SetBinContent(i, (1. + a->GetBinContent(i))/(1. + b->GetBinContent(i)) - 1.); return o; }

// one panel: zero line, assigned band (+-), up and down lines; returns nothing
void panel(TH1D *up, TH1D *dn, TH1D *assigned){
  if(assigned){
    TH1D *b = pct(assigned, Form("%s_band", assigned->GetName())); TH1D *lo = (TH1D*) b->Clone(Form("%s_lo", b->GetName()));
    for(int i = 1; i <= b->GetNbinsX(); i++){ const double v = fabs(b->GetBinContent(i)); b->SetBinContent(i, 0.); b->SetBinError(i, v); }
    b->SetFillColorAlpha(kGray + 1, 0.35); b->SetLineWidth(0); b->SetMarkerSize(0); b->Draw("E2 SAME");
    for(int i = 1; i <= assigned->GetNbinsX(); i++){
      const double lo_ = assigned->GetXaxis()->GetBinLowEdge(i); if(lo_ < xLo - 1e-6 || lo_ >= xHi - 1e-6) continue;
      maxDiff = std::max(maxDiff, fabs(fabs(assigned->GetBinContent(i)) - std::max(fabs(up->GetBinContent(i)), fabs(dn->GetBinContent(i)))));
    }
    delete lo;
  }
  TLine l; l.SetLineStyle(1); l.SetLineColor(kGray + 2); l.SetLineWidth(2); l.DrawLine(xLo, 0., xHi, 0.);
  TH1D *u = pct(up, Form("%s_pct", up->GetName())), *d = pct(dn, Form("%s_pct", dn->GetName()));
  styleLine(u, hexUp, 4); styleLine(d, hexDn, 4); d->SetLineStyle(2);
  u->Draw("HIST SAME"); d->Draw("HIST SAME");
}
// symmetric y range over the reported bins of the given histograms, rounded
double yRange(std::vector<TH1D*> hs){
  double m = 0.;
  for(TH1D *h : hs) if(h) for(int i = 1; i <= h->GetNbinsX(); i++){ const double lo = h->GetXaxis()->GetBinLowEdge(i); if(lo >= xLo - 1e-6 && lo < xHi - 1e-6) m = std::max(m, 100.*fabs(h->GetBinContent(i))); }
  m *= 1.45; return m > 20. ? 5.*ceil(m/5.) : (m > 8. ? 2.*ceil(m/2.) : std::max(2., ceil(m)));
}
// in-panel text in PIXELS (font 43 / 63), so it prints the same size in every pad whatever the pad's shape
// (2026-11-17, user: the "reference" text came out in two sizes with the pad-relative sizes)
double pxText = 34., pxLeg = 30.;   // set per canvas by setTextScale(): the same printed size once the figure is scaled to the page width
void setTextScale(TCanvas *cv){ const double k = cv->GetWw()/1800.; pxText = 34.*k; pxLeg = 30.*k; }
double padHpx(){ return gPad->GetWh()*gPad->GetAbsHNDC(); }
double padWpx(){ return gPad->GetWw()*gPad->GetAbsWNDC(); }
void ptext(double x, double y, const char *s, double px, int font = 43, int align = 13, Color_t col = kBlack){
  TLatex t; t.SetNDC(); t.SetTextFont(font); t.SetTextSize(px); t.SetTextAlign(align); t.SetTextColor(col); t.DrawLatex(x, y, s);
}
// panel label at the frame's top left, 20 px in, 14 px down
void panelLabel(GridPads &g, int ix, int iy, const char *s){ ptext(g.frameLeft(ix) + 20./padWpx(), g.frameTop(iy) - 14./padHpx(), s, pxText, 63, 13); }
void rowLabel(TCanvas *cv, GridPads &g, int iy, const char *s, double size, double x = 0.042){
  cv->cd(); TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(size); t.SetTextAngle(90); t.SetTextAlign(21);
  const double ih = (1. - g.bm - g.tm)/g.ny; t.DrawLatex(x, 1. - g.tm - (iy + 0.5)*ih, s);
}
void legendUpDn(GridPads &g, int ix, int iy, bool band, double y2){
  g.at(ix, iy)->cd();
  static int nL = 0; nL++;
  TH1D *pu = new TH1D(Form("lgU%d", nL), "", 1, 0, 1), *pd = new TH1D(Form("lgD%d", nL), "", 1, 0, 1), *pb = new TH1D(Form("lgB%d", nL), "", 1, 0, 1);
  for(TH1D *h : {pu, pd, pb}) h->SetDirectory(nullptr);
  styleLine(pu, hexUp, 4); styleLine(pd, hexDn, 4); pd->SetLineStyle(2); pb->SetFillColorAlpha(kGray + 1, 0.35); pb->SetLineWidth(0);
  const int n = band ? 3 : 2; const double rowH = 1.4*pxLeg/padHpx();
  TLegend *lg = makeLegend(g.frameLeft(ix) + 20./padWpx(), y2 - n*rowH, 0.98, y2); lg->SetTextFont(43); lg->SetTextSize(pxLeg);
  lg->AddEntry(pu, "JES up", "l"); lg->AddEntry(pd, "JES down", "l"); if(band) lg->AddEntry(pb, "assigned", "f"); lg->Draw();
}
}

void plotBJetJESSys_PF(TString bTag = "_noFlatJEU_fakeGeoOnly_csvV2_2D_noRecoFloor_normSep_noStitch", TString ppVar = "_csvV2_2D", TString outSfx = ""){   // outSfx: appended to the figure names
  using namespace bjes;
  initPlotStyle();
  const char *repo = "/home/clayton/Analysis/code/bJetRaaAnalysis";
  TFile *fB = TFile::Open(Form("%s/rootFiles/systematics/bJetsPerZSys_PF_MBcaloMatched_JEUppCorr%s.root", repo, bTag.Data()));
  TFile *fS = TFile::Open(Form("%s/rootFiles/MuTagMuTrigBJetSpectra/Data/bJetSpectrum_PFJets_pp%s.root", repo, ppVar.Data()));
  if(!fB || fB->IsZombie() || !fS || fS->IsZombie()){ printf("ERROR: input file missing\n"); return; }
  auto get = [&](TFile *f, TString n) -> TH1D* { TH1D *h = nullptr; f->GetObject(n, h); if(!h) printf("ERROR: %s missing in %s\n", n.Data(), f->GetName()); return h; };
  if(!get(fB, "h_varRel_JEU_up_b_C1")) { printf("  -> rerun bJetsPerZSys_PF_PbPb.C (2026-11-16 or later saves the JEU variations)\n"); return; }
  const char *dir[2] = {"up", "down"};
  const char *sys[5] = {"pp", "C4", "C3", "C2", "C1"}, *name[5] = {"pp", "50-80%", "30-50%", "10-30%", "0-10%"};

  // ---- per-system rows
  TH1D *meas[2][5], *bU[2][5], *iU[2][5], *fr[2][5];
  {
    TH1D *mpp = get(fS, "h_bJetPt"); TH1D *mc = (TH1D*) mpp->Clone("ppMeasC"); mc->SetDirectory(nullptr);
    for(int i = 1; i <= mc->GetNbinsX(); i++) mc->SetBinContent(i, mc->GetBinContent(i)*mc->GetBinWidth(i));   // counts, as the pp chain
    const std::vector<double> dPP = jeuDelta(Form("%s/JetEnergyCorrections/Spring18_ppRef5TeV_V6_MC_Uncertainty_AK4PF.txt", repo));
    for(int k = 0; k < 2; k++){ TH1D *m = jeuShiftSpectrum(mc, dPP, k == 0 ? +1 : -1, Form("ppMeas_%s", dir[k])); meas[k][0] = rel(m, mc, Form("measRel_%s_pp", dir[k])); }
  }
  for(int s = 0; s < 5; s++) for(int k = 0; k < 2; k++){
    if(s > 0) meas[k][s] = rel(get(fB, Form("h_bJEU_meas_%s_%s", dir[k], sys[s])), get(fB, Form("h_bJEU_measNom_%s", sys[s])), Form("measRel_%s_%s", dir[k], sys[s]));
    bU[k][s] = get(fB, Form("h_varRel_JEU_%s_b_%s", dir[k], sys[s]));
    iU[k][s] = get(fB, Form("h_varRel_JEU_%s_incl_%s", dir[k], sys[s]));
    fr[k][s] = relRatio(bU[k][s], iU[k][s], Form("frRel_%s_%s", dir[k], sys[s]));
  }
  TH1D *aB[5], *aF[5];
  for(int s = 0; s < 5; s++){ aB[s] = get(fB, Form("h_sysRel_JEU_b_%s", sys[s])); aF[s] = get(fB, Form("h_sysRel_JEU_bOverIncl_%s", sys[s])); }
  {
    const int NR = 4;
    const char *rowLab[NR] = {"measured #it{b} jets", "unfolded #it{b} jets", "unfolded inclusive", "#it{b} / inclusive"};
    TCanvas *cv = new TCanvas("cvJESsys", "", 2200, 1700); setTextScale(cv);
    GridPads g = makeGridPads(cv, 5, NR, 0.095, 0.012, 0.075, 0.045);
    // one y range for every row (2026-11-17, user: a common scale within the figure)
    std::vector<TH1D*> all;
    for(int r = 0; r < NR; r++) for(int s = 0; s < 5; s++) for(int k = 0; k < 2; k++){ TH1D *h = r == 0 ? meas[k][s] : r == 1 ? bU[k][s] : r == 2 ? iU[k][s] : fr[k][s]; all.push_back(h); }
    const double y = yRange(all);
    for(int r = 0; r < NR; r++){
      for(int s = 0; s < 5; s++){
        g.at(s, r)->cd();
        TH1F *f = gPad->DrawFrame(xLo, -y, xHi, y); g.styleFrame(f, s, r, 0.075);
        if(r != NR - 1) f->GetXaxis()->SetLabelSize(0.);
        sameLabels(f, 30.);
        if(r == 0) panel(meas[0][s], meas[1][s], nullptr);
        if(r == 1) panel(bU[0][s], bU[1][s], aB[s]);
        if(r == 2) panel(iU[0][s], iU[1][s], nullptr);
        if(r == 3) panel(fr[0][s], fr[1][s], aF[s]);
        gPad->RedrawAxis();
        if(r == 0) panelLabel(g, s, r, name[s]);
      }
      rowLabel(cv, g, r, rowLab[r], 0.027);
    }
    g.at(0, 1)->cd(); legendUpDn(g, 0, 1, true, g.frameTop(1) - 14./padHpx());
    g.titles(cv, "#it{p}_{T}^{jet} (GeV)", "", 0.034, 0.012, 0.016);
    cv->cd(); TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.027); t.SetTextAngle(90); t.SetTextAlign(21);
    t.DrawLatex(0.012, g.bm + 0.5*(1. - g.bm - g.tm), "shifted / nominal #minus 1 (%)");
    g.cmsLabel(cv, "Internal", "jet energy scale, #it{b} jets,   anti-#it{k}_{T} #it{R} = 0.4 PF jets,  #sqrt{#it{s}_{NN}} = 5.02 TeV", 0.95);
    const TString dirOut = Form("%s/figures/systematics/JEU/", repo); gSystem->mkdir(dirOut, true);
    savePdfTight(cv, dirOut + "bJetJESSys_perSystem_PF" + outSfx + ".pdf");
  }

  // ---- ratios between systems
  {
    const char *cls[4] = {"C4", "C3", "C2", "C1"}, *cName[4] = {"50-80%", "30-50%", "10-30%", "0-10%"};
    const int NR = 4;
    const char *rowLab[NR] = {"#it{R}_{Z}(#it{b})", "#it{R}_{Z}(#it{b}) / #it{R}_{Z}(incl)", "#it{R}_{CP}^{Z}(#it{b})", "#it{R}_{CP}^{Z}(#it{b}) / #it{R}_{CP}^{Z}(incl)"};
    TCanvas *cv = new TCanvas("cvJESratio", "", 1800, 1700); setTextScale(cv);
    GridPads g = makeGridPads(cv, 4, NR, 0.135, 0.012, 0.075, 0.045);
    TH1D *LL[NR][2][4] = {{{nullptr}}}, *AA[NR][4] = {{nullptr}};
    for(int r = 0; r < NR; r++){
      TH1D *(&L)[2][4] = LL[r]; TH1D *(&A)[4] = AA[r];
      for(int c = 0; c < 4; c++){
        const int s = c + 1;   // index into sys[]: C4 .. C1
        const bool rcp = r >= 2; if(rcp && c == 0) continue;   // R_CP: 50-80% is the reference
        for(int k = 0; k < 2; k++){
          const int ref = rcp ? 1 : 0;   // pp or 50-80%
          L[k][c] = r % 2 == 0 ? relRatio(bU[k][s], bU[k][ref], Form("rr%d_%d_%d", r, k, c)) : relRatio(fr[k][s], fr[k][ref], Form("rr%d_%d_%d", r, k, c));
        }
        A[c] = get(fB, Form("h_sysRel_JEU_%s_%s_%s", r % 2 == 0 ? "b" : "bOverIncl", rcp ? "RCP" : "PbPbOverPP", cls[c]));
      }
    }
    // one y range for every row (2026-11-17, user)
    std::vector<TH1D*> all; for(int r = 0; r < NR; r++) for(int c = 0; c < 4; c++) for(int k = 0; k < 2; k++) all.push_back(LL[r][k][c]);
    const double y = yRange(all);
    for(int r = 0; r < NR; r++){
      TH1D *(&L)[2][4] = LL[r]; TH1D *(&A)[4] = AA[r];
      for(int c = 0; c < 4; c++){
        g.at(c, r)->cd();
        TH1F *f = gPad->DrawFrame(xLo, -y, xHi, y); g.styleFrame(f, c, r, 0.075);
        if(r != NR - 1) f->GetXaxis()->SetLabelSize(0.);
        sameLabels(f, 30.);
        if(L[0][c]) panel(L[0][c], L[1][c], A[c]);
        gPad->RedrawAxis();
        if(r == 0 || (r == 2 && c > 0)) panelLabel(g, c, r, r < 2 ? Form("%s / pp", cName[c]) : Form("%s / 50-80%%", cName[c]));
        if(r >= 2 && c == 0) ptext(0.5*(g.frameLeft(0) + 1.), 0.5*(g.frameTop(r) + g.frameBottom(r)), "reference", pxText, 43, 22, kGray + 2);
      }
      rowLabel(cv, g, r, rowLab[r], 0.027, 0.060);
    }
    g.at(0, 1)->cd(); legendUpDn(g, 0, 1, true, g.frameTop(1) - 14./padHpx());
    g.titles(cv, "#it{p}_{T}^{jet} (GeV)", "", 0.034, 0.012, 0.016);
    cv->cd(); TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.027); t.SetTextAngle(90); t.SetTextAlign(21);
    t.DrawLatex(0.030, g.bm + 0.5*(1. - g.bm - g.tm), "shifted / nominal #minus 1 (%)");
    g.cmsLabel(cv, "Internal", "jet energy scale, ratios,   #sqrt{#it{s}_{NN}} = 5.02 TeV", 0.95);
    savePdfTight(cv, Form("%s/figures/systematics/JEU/bJetJESSys_ratios_PF%s.pdf", repo, outSfx.Data()));
  }
  printf("largest |assigned - max(|up|, |down|)| over every band drawn, 80-200 GeV: %.2e\n", maxDiff);
  // printout: the values in the bands, % (pp, 50-80 ... 0-10%)
  printf("\nassigned JES (%%), 80-100 / 100-120 / 120-150 / 150-200 GeV\n");
  auto row = [&](const char *lab, TH1D *h){ printf("  %-34s", lab); for(int i = 1; i <= h->GetNbinsX(); i++){ const double lo = h->GetXaxis()->GetBinLowEdge(i); if(lo >= xLo - 1e-6 && lo < xHi - 1e-6) printf(" %5.1f", 100.*fabs(h->GetBinContent(i))); } printf("\n"); };
  for(int s = 0; s < 5; s++) row(Form("b jets per Z, %s", name[s]), aB[s]);
  for(int s = 0; s < 5; s++) row(Form("b / inclusive, %s", name[s]), aF[s]);
  for(int c = 4; c >= 1; c--){ row(Form("R_Z(b), C%d", c), get(fB, Form("h_sysRel_JEU_b_PbPbOverPP_C%d", c))); row(Form("double ratio, C%d", c), get(fB, Form("h_sysRel_JEU_bOverIncl_PbPbOverPP_C%d", c))); }
  for(int c = 3; c >= 1; c--){ row(Form("R_CP(b), C%d", c), get(fB, Form("h_sysRel_JEU_b_RCP_C%d", c))); row(Form("R_CP(b)/R_CP(incl), C%d", c), get(fB, Form("h_sysRel_JEU_bOverIncl_RCP_C%d", c))); }
}
