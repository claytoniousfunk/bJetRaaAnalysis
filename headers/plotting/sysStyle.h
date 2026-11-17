// One look for every "systematic uncertainty by source" figure (2026-11-16, user: remake the
// systematics plots in the style of the b-purity one).
//
//  * sysStyle(key): a FIXED Okabe-Ito colour + line style + legend label per source key, so a
//    source looks the same in every figure (purity, jets per Z, R_Z, ratios). Colours repeat only
//    with a different line style; the total is a grey band with a thick black edge.
//  * drawSysRow(...): 1 x N shared-border row (gridPads.h), sources over the total band, legend
//    split over the left-most panels INSIDE their frames, below the panel label; the y range is
//    raised just enough that no curve runs into a legend. Header with jet definition and energy.
//
//   SysPanel p[N]; p[i].name = "pp"; p[i].tot = h_total; p[i].src[k] = h_source_k (or nullptr);
//   drawSysRow(N, p, nSrc, keys, labels /*nullptr = sysStyle labels*/, 80., 200.,
//              "inclusive #it{R}_{Z}" /*quantity, in the header*/, "anti-#it{k}_{T} ...", pdfName);
// Histograms are RELATIVE uncertainties (fractions); they are drawn in %, |content|, over
// [xLo, xHi) only. A source that is zero (< 0.05%) or absent in every panel is left out.
#pragma once
#include <vector>
#include <algorithm>
#include "TCanvas.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TLine.h"
#include "TString.h"
#include "TMath.h"
#include "plotStyle.h"
#include "gridPads.h"

const int SYS_NMAX = 24;

struct SysStyle { const char *key, *hex, *label; int lineStyle; };
inline SysStyle sysStyle(const char *key)
{
  // style 11 = medium dashes; ROOT's long-dash 9 shows no gap in a short legend swatch (2026-11-16, user)
  static bool init = false; if(!init){ gStyle->SetLineStyleString(11, "30 18"); init = true; }
  // solid = b purity family, dashed / dotted = unfolding and normalization family.
  // The TABLE ORDER is the reporting order of every figure legend (2026-11-16, user: "a logical order"):
  // jet energy, b purity, unfolding, muon tagging, combinatorial jets / stitching, normalization;
  // printSysSummaryTables_PF.C lists its rows the same way.
  static const SysStyle tab[] = {
    {"JEU",               "#0072B2", "JES",                                   2},   // short dashes: style 7 read as solid in a legend swatch
    {"JER",               "#009E73", "JER",                                   11},  // medium dashes (registered below): style 9 read as solid in a legend swatch
    // b purity (purity figure: own sources; b-jet chain: propagated groups)
    {"bGS",               "#E69F00", "g#rightarrowb#bar{b} fraction",         1},
    {"cMult",             "#56B4E9", "c fraction",                            1},
    {"fitRangeHigh",      "#CC79A7", "fit range",                             1},
    {"tplStat",           "#000000", "MC #varepsilon_{b} stat.",              3},
    {"fakeMuJet",         "#D55E00", "uncorr. #mu-jet pairs",                 1},
    {"bPurity_bGS",       "#E69F00", "purity: g#rightarrowb#bar{b}",          1},
    {"bPurity_cMult",     "#56B4E9", "purity: c fraction",                    1},
    {"bPurity",           "#CC79A7", "purity: other",                         1},
    {"bTagWP",            "#CC79A7", "b-tag working point",                   10},
    // unfolding
    {"iterations",        "#D55E00", "iterations",                            7},
    {"unfoldingBias",     "#CC79A7", "non-closure",                           7},
    {"responseStat",      "#000000", "response stat.",                        3},
    {"prior",             "#009E73", "prior",                                 3},
    {"belowFloorWeight",  "#0072B2", "below-floor weight",                    3},
    {"recoFloor",         "#000000", "reco floor",                            10},
    {"incl_recoFloor",    "#000000", "reco floor",                            10},
    // tagging, combinatorial jets, stitching
    {"muonEff",           "#D55E00", "muon efficiency",                       3},
    {"muonTagFreq",       "#56B4E9", "muon-tag freq.",                        7},
    {"fakeJets",          "#E69F00", "combinatorial jets",                    7},
    {"incl_fakeJets",     "#E69F00", "combinatorial jets",                    7},
    {"stitchWindow",      "#009E73", "trigger stitch window",                 3},
    {"incl_stitchWindow", "#009E73", "trigger stitch window",                 3},
    {"stitchNormStat",    "#0072B2", "trigger stitch norm.",                  10},
    {"incl_stitchNormStat","#0072B2","trigger stitch norm.",                  10},
    // normalization (usually kept apart as a box and not drawn)
    {"NZstat",            "#56B4E9", "#it{N}_{Z} stat.",                      3},
    {"zBkg",              "#E69F00", "#it{Z} background",                     3},
    {"massWindow",        "#CC79A7", "#it{Z} mass window",                    3},
    {"pubZYield",         "#000000", "published #it{Z} yield",                7},
    {"pubZAcc",           "#000000", "pp |#it{y}_{Z}| < 2.1 fraction",        9},
  };
  for(const auto &s : tab) if(TString(s.key) == key) return s;
  return {key, "#000000", key, 1};
}
// position of a source in the reporting order (unknown keys last)
inline int sysRank(const char *key)
{
  static const char *order[] = {"JEU", "JER", "bGS", "cMult", "fitRangeHigh", "tplStat", "fakeMuJet", "bPurity_bGS", "bPurity_cMult", "bPurity", "bTagWP",
    "iterations", "unfoldingBias", "responseStat", "prior", "belowFloorWeight", "recoFloor", "incl_recoFloor", "muonEff", "muonTagFreq",
    "fakeJets", "incl_fakeJets", "stitchWindow", "incl_stitchWindow", "stitchNormStat", "incl_stitchNormStat", "NZstat", "zBkg", "massWindow", "pubZYield", "pubZAcc"};
  for(int i = 0; i < (int)(sizeof(order)/sizeof(order[0])); i++) if(TString(order[i]) == key) return i;
  return 1000;
}

struct SysPanel { TString name; TH1 *tot = nullptr; TH1 *src[SYS_NMAX] = {nullptr}; };

inline void drawSysRow(int nP, SysPanel *P, int nSrc, const char *const *keys, const char *const *labels,
                       double xLo, double xHi, const char *what, const char *header, const TString &pdf,
                       const char *cmsStatus = "Internal", bool binTicks = false)
{
  // binTicks (2026-11-16): thin dotted lines at the internal bin edges of the
  // total, so a fine binning is visible where neighbouring bins carry equal values (e.g. sources inherited
  // from coarser windows, labelled as such by the caller)
  static int nCall = 0; nCall++;
  // in %, |content|, zero outside [xLo, xHi)
  auto clip = [&](TH1 *h, const TString &n) -> TH1D* { if(!h) return nullptr;
    std::vector<double> ed; for(int b = 1; b <= h->GetNbinsX() + 1; b++) ed.push_back(h->GetBinLowEdge(b));
    TH1D *o = new TH1D(n, "", h->GetNbinsX(), ed.data());
    o->SetDirectory(nullptr);
    for(int b = 1; b <= h->GetNbinsX(); b++){ const double lo = h->GetBinLowEdge(b);
      if(lo >= xLo - 1e-6 && lo < xHi - 1e-6) o->SetBinContent(b, 100.*fabs(h->GetBinContent(b))); }
    return o; };
  std::vector<TH1D*> T(nP, nullptr); std::vector<std::vector<TH1D*>> S(nP, std::vector<TH1D*>(nSrc, nullptr));
  std::vector<bool> used(nSrc, false);
  std::vector<double> pMax(nP, 0.);
  for(int i = 0; i < nP; i++){
    T[i] = clip(P[i].tot, Form("sr%d_tot%d", nCall, i)); if(T[i]) pMax[i] = T[i]->GetMaximum();
    for(int s = 0; s < nSrc; s++){ S[i][s] = clip(P[i].src[s], Form("sr%d_s%d_%d", nCall, i, s));
      if(S[i][s]){ const double m = S[i][s]->GetMaximum(); if(m > 0.05) used[s] = true; pMax[i] = std::max(pMax[i], m); } }
  }
  std::vector<int> drawn; for(int s = 0; s < nSrc; s++) if(used[s]) drawn.push_back(s);
  std::stable_sort(drawn.begin(), drawn.end(), [&](int a, int b){ return sysRank(keys[a]) < sysRank(keys[b]); });
  const int nEnt = (int) drawn.size() + 1;   // + total

  const double W = 440.*nP + 160.;
  TCanvas *cv = new TCanvas(Form("cvSysRow%d", nCall), "", (int) W, 760);
  GridPads g = makeGridPads(cv, nP, 1, 0.060*2360./W, 0.020*2360./W, 0.165, 0.095);
  // legend placement: try 1 .. nP legend panels (entries split evenly, at least 3 rows each unless
  // fewer entries), each time in the panels with the SMALLEST values (kept in left-to-right order),
  // and keep the option that needs the lowest y range; ties go to fewer panels. The y range is
  // 1.25 x the largest value and high enough that a legend panel's curves end below its legend.
  const double rowH = 0.080, labelGap = 0.17;
  g.at(0, 0)->cd();
  const double fTop = g.frameTop(0), fBot = g.frameBottom(0);
  double yBase = 0.; for(double m : pMax) yBase = std::max(yBase, 1.25*m);
  std::vector<int> order(nP); for(int i = 0; i < nP; i++) order[i] = i;
  std::stable_sort(order.begin(), order.end(), [&](int a, int b){ return pMax[a] < pMax[b] - 1e-9; });
  double yMax = 1e30; std::vector<int> legPanel, legN;
  for(int nL = 1; nL <= nP; nL++){
    const int per = (nEnt + nL - 1)/nL;
    if(nL > 1 && per < 3 && nEnt >= 3) break;
    std::vector<int> pan(order.begin(), order.begin() + nL); std::sort(pan.begin(), pan.end());
    std::vector<int> ns; int left = nEnt; double y = yBase;
    for(int L = 0; L < nL; L++){
      const int n = std::min(per, left); left -= n; ns.push_back(n);
      const double bottom = fTop - labelGap - rowH*n, frac = (bottom - fBot)/(fTop - fBot) - 0.04;
      y = frac > 0.05 ? std::max(y, pMax[pan[L]]/frac) : 1e29;
    }
    if(y < 0.97*yMax){ yMax = y; legPanel = pan; legN = ns; }
  }
  yMax = yMax > 20. ? 5.*ceil(yMax/5.) : (yMax > 8. ? 2.*ceil(yMax/2.) : ceil(yMax));

  for(int ix = 0; ix < nP; ix++){
    g.at(ix, 0)->cd();
    TH1F *fr = gPad->DrawFrame(xLo, 0., xHi, yMax); g.styleFrame(fr, ix, 0, 0.060);
    fr->GetXaxis()->SetLabelSize(0.060*g.padW(0)/g.padW(ix));   // same printed size in every pad
    fr->GetXaxis()->SetTickLength(0.03*g.padW(0)/g.padW(ix));
    fr->GetYaxis()->SetTickLength(0.015);
    if(binTicks && T[ix]){   // behind everything: thin dotted lines at the internal bin edges
      TLine tk; tk.SetLineColor(kGray + 1); tk.SetLineWidth(1); tk.SetLineStyle(3);
      for(int b = 2; b <= T[ix]->GetNbinsX(); b++){ const double e = T[ix]->GetBinLowEdge(b);
        if(e > xLo + 1e-6 && e < xHi - 1e-6) tk.DrawLine(e, 0., e, std::min(yMax, 1.15*pMax[ix])); }   // up to just above the panel's values: clear of the legends
    }
    if(T[ix]){ T[ix]->SetFillColorAlpha(kGray + 1, 0.35); T[ix]->SetLineColor(kBlack); T[ix]->SetLineWidth(4); T[ix]->Draw("HIST SAME"); }
    for(int s : drawn){
      if(!S[ix][s]) continue;
      const SysStyle st = sysStyle(keys[s]);
      styleLine(S[ix][s], st.hex, 4); S[ix][s]->SetLineStyle(st.lineStyle); S[ix][s]->Draw("HIST SAME");
    }
    gPad->RedrawAxis();
    g.text(ix, 0, 0.10, P[ix].name.Data(), 0.066, 62);
  }
  int next = 0;   // index into [total, drawn...]
  for(size_t L = 0; L < legPanel.size(); L++){
    g.at(legPanel[L], 0)->cd();
    const int n = legN[L];
    const double top = g.frameTop(0) - labelGap;
    TLegend *lg = makeLegend(g.frameLeft(legPanel[L]) + 0.05, top - rowH*n, 0.99, top, 0.050*g.sizeScale());
    for(int e = 0; e < n; e++, next++){
      if(next == 0){
        TH1D *pt = new TH1D(Form("sr%d_legT", nCall), "", 1, 0., 1.); pt->SetDirectory(nullptr);
        pt->SetFillColorAlpha(kGray + 1, 0.35); pt->SetLineColor(kBlack); pt->SetLineWidth(4);
        lg->AddEntry(pt, "total", "f"); continue;
      }
      const int s = drawn[next - 1]; const SysStyle st = sysStyle(keys[s]);
      TH1D *proxy = new TH1D(Form("sr%d_leg%d", nCall, s), "", 1, 0., 1.); proxy->SetDirectory(nullptr);
      styleLine(proxy, st.hex, 4); proxy->SetLineStyle(st.lineStyle);
      lg->AddEntry(proxy, labels && labels[s] ? labels[s] : st.label, "l");
    }
    lg->Draw();
  }
  g.titles(cv, "#it{p}_{T}^{jet} (GeV)", "relative uncertainty (%)", 0.058, 0.022*2360./W, 0.025);
  g.cmsLabel(cv, cmsStatus, Form("%s,   %s", what, header), 1.35);
  savePdfTight(cv, pdf);
}
