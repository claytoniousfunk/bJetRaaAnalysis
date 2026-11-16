// Systematic uncertainties by source vs jet pT, one figure per quantity, for
// talks (2026-10-20, user): the figure form of printSysSummaryTables_PF.C --
// same files, same sources, same 80-200 GeV window, nothing recomputed.
//   bPurity       pp + 0-10 / 10-30 / 30-50 / 50-80%
//   inclJetsPerZ  pp + four classes
//   bJetsPerZ     pp + four classes (nominal reverse unfolding)
//   RZ_incl, RZ_b, doubleRatio   four classes
// Panel order (2026-10-20, user): pp, then 50-80 / 30-50 / 10-30 / 0-10%.
// Every source drawn on its own (not grouped, unlike
// src/calculateJetsPerZ/plotRZSysBreakdown_PF.C); total in black. A source
// that is zero in every panel of a figure is left out, as in the tables.
// Panels share borders (gridPads.h); legend in the empty bottom-right cell
// (pp + classes, 3x2) or to the right of the 2x2 grid (classes only).
//
// inputs  as printSysSummaryTables_PF.C
// Usage: root -l -b -q plotSysBreakdown_PF.C   (from src/systematics/)
//        root -l -b -q 'plotSysBreakdown_PF.C("_noFlatJEU")'   inclusive JEU not flattened
//        root -l -b -q 'plotSysBreakdown_PF.C("_noFlatJEU","_noFlatJEU_noUnfBias")'   b side too
//        (inclusiveJetsPerZSys_PF.C's noFlatJEU file; the b-jet file is unchanged)
// Output: figures/systematics/sysBreakdown_PF_<quantity>.pdf
#include "TFile.h"
#include "../../headers/functions/purityVariant.h"   // purity files of the purity variant (BJET_PURITY_VARIANT; 2026-10-28)
#include "TH1D.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TColor.h"
#include <vector>
#include "../../headers/plotting/plotStyle.h"
#include "../../headers/plotting/gridPads.h"

// inclTag: suffix of the inclusive file (e.g. "_noFlatJEU"), also added to the
// names of the two inclusive figures
// bTag: suffix of the b-jet file after _reverseUnfold (e.g. "_noFlatJEU_noUnfBias"),
// also added to the names of the b-jet, R_Z(b) and double-ratio outputs
// ptMaxOpt (2026-10-21): upper edge of the drawn window (default 200 GeV, the talk window);
// 300 for the AN, where the b-jet results go to 300 GeV; outputs then carry _to<ptMax>
// simple-surface variant (BJET_PURITY_VARIANT=csvV2__surf): its light-yield source is stored under the name cMult
inline bool surfVariantLabels(){ const char *e = std::getenv("BJET_PURITY_VARIANT"); return e && TString(e) == "csvV2__surf"; }
void plotSysBreakdown_PF(TString inclTag = "", TString bTag = "", double ptMaxOpt = 200.){
  const char *repo = "/home/clayton/Analysis/code/bJetRaaAnalysis";
  const double ptMin = 80., ptMax = ptMaxOpt;
  TFile *fPurPP = TFile::Open(Form("%s/rootFiles/systematics/bPuritySys_total_PFJets_pp%s.root", repo, ppPurVar()));
  TFile *fPurAA = TFile::Open(Form("%s/rootFiles/systematics/bPuritySys_PFJets_PbPb%s.root", repo, purVar()));
  TFile *fIncl  = TFile::Open(Form("%s/rootFiles/systematics/inclJetsPerZSys_PF_MBcaloMatched_JEUppCorr%s.root", repo, inclTag.Data()));
  // bTag starting with "_bayes" (or "bayes") (2026-10-21: Bayes, MC prior, is the b-jet nominal again): the Bayes
  // file bJetsPerZSys_PF_MBcaloMatched_JEUppCorr<rest>.root instead of the reverse-unfolding one
  const bool bayesB = bTag.BeginsWith("_bayes") || bTag.BeginsWith("bayes");
  const TString bRest = bayesB ? TString(bTag(bTag.Index("bayes") + 5, bTag.Length())) : bTag;
  TFile *fB     = TFile::Open(Form("%s/rootFiles/systematics/bJetsPerZSys_PF_MBcaloMatched_JEUppCorr%s%s.root", repo, bayesB ? "" : "_reverseUnfold", bRest.Data()));
  if(!fPurPP || !fPurAA || !fIncl || !fB){ printf("ERROR: a systematics file did not open\n"); return; }
  initPlotStyle();

  struct Src { const char *key, *label; };
  struct Col { TString head; TFile *f; TString fmt; };   // fmt: %s = source key (TString: Form() buffers rotate)
  // line styles: six Okabe-Ito colours (no yellow, no black = total) x solid / dashed / dotted
  const char *hex[6] = {okabeHex[1], okabeHex[2], okabeHex[3], okabeHex[5], okabeHex[6], okabeHex[7]};
  const int sty[3] = {1, 7, 3};
  int nFig = 0;

  // relative uncertainty in %, restricted to the window; nullptr when missing or zero there
  auto get = [&](TFile *f, TString name) -> TH1* {
    TH1 *h = nullptr; f->GetObject(name, h);
    if(!h) return nullptr;
    TH1 *o = (TH1*) h->Clone(Form("sb_%d_%s", nFig, name.Data())); o->SetDirectory(nullptr); o->Reset(); o->SetTitle("");
    double mx = 0.;
    for(int i = 1; i <= h->GetNbinsX(); i++){
      if(h->GetBinLowEdge(i) < ptMin - 1e-6 || h->GetBinLowEdge(i+1) > ptMax + 1e-6) continue;
      o->SetBinContent(i, 100.*fabs(h->GetBinContent(i))); mx = std::max(mx, o->GetBinContent(i));
    }
    if(mx < 0.05){ delete o; return nullptr; }
    o->GetXaxis()->SetRangeUser(ptMin, ptMax - 1e-6);
    return o;
  };

  auto figure = [&](const char *tag, const char *what, const std::vector<Src> &srcs, const std::vector<Col> &cols, TString sfx = ""){
    nFig++;
    const int nP = cols.size();
    // keep a source if it is non-zero in any panel; y range from the largest total
    std::vector<int> keep; double yMax = 0.;
    for(size_t k = 0; k < srcs.size(); k++){
      bool any = false;
      for(auto &c : cols){ TH1 *h = get(c.f, Form(c.fmt.Data(), srcs[k].key)); if(h){ any = true; if(TString(srcs[k].key) == "total") yMax = std::max(yMax, h->GetMaximum()); delete h; } }
      if(any && TString(srcs[k].key) != "total") keep.push_back(k);
    }
    yMax = 5.*ceil(1.12*yMax/5.);
    const bool wide = nP == 5;   // pp + classes: 3x2, legend in the sixth pad
    const int nx = wide ? 3 : 2, ny = 2;
    TCanvas *cv = new TCanvas(Form("cv_%s", tag), "", wide ? 1500 : 1400, 950);
    GridPads g = makeGridPads(cv, nx, ny, wide ? 0.08 : 0.085, wide ? 0.02 : 0.33, 0.09, 0.07);
    std::vector<std::pair<TH1*, TString>> ent;   // legend entries, from the first panel
    for(int p = 0; p < nP; p++){
      const int ix = p % nx, iy = p / nx;
      g.at(ix, iy)->cd();
      TH1F *fr = gPad->DrawFrame(ptMin, 0., ptMax, yMax); g.styleFrame(fr, ix, iy, 0.052);
      const bool first = p == 0;
      for(size_t j = 0; j < keep.size(); j++){
        TH1 *h = get(cols[p].f, Form(cols[p].fmt.Data(), srcs[keep[j]].key));
        TH1D *proxy = nullptr;
        if(!h && first){ proxy = new TH1D(Form("lp_%s_%zu", tag, j), "", 1, 0, 1); proxy->SetDirectory(nullptr); }
        TH1 *s = h ? h : (TH1*) proxy;
        if(!s) continue;   // source absent in this panel
        s->SetLineColor(TColor::GetColor(hex[j % 6])); s->SetLineStyle(sty[(j / 6) % 3]); s->SetLineWidth(3);
        if(h) h->Draw("HIST same");
        if(first) ent.push_back({s, srcs[keep[j]].label});
      }
      TH1 *tot = get(cols[p].f, Form(cols[p].fmt.Data(), "total"));
      if(tot){ tot->SetLineColor(kBlack); tot->SetLineWidth(4); tot->Draw("HIST same"); if(first) ent.push_back({tot, "total"}); }
      // a normalization kept apart (h_normRel_*, tools/dropSysSources.C asNorm) is not drawn
      // here (2026-10-20, user): it goes on the results figures as a box at unity
      gPad->RedrawAxis();
      g.text(ix, iy, 0.10, cols[p].head.Data(), 0.060);
    }
    // legend, canvas NDC: the empty bottom-right cell of the 3x2 grid, or the right margin of the 2x2
    cv->cd();
    const double iw = (1. - g.lm - g.rm)/nx, ih = (1. - g.bm - g.tm)/ny, rowH = 0.042;
    double lx1, lx2, ly1, ly2;
    if(wide){ lx1 = g.lm + 2.*iw + 0.02; lx2 = 1. - g.rm; ly2 = g.bm + ih - 0.01; ly1 = std::max(0.005, ly2 - rowH*ent.size()); }
    else    { lx1 = 1. - g.rm + 0.015;   lx2 = 0.995;     ly2 = 1. - g.tm;         ly1 = std::max(g.bm, ly2 - rowH*ent.size()); }
    TLegend *lg = new TLegend(lx1, ly1, lx2, ly2);
    for(auto &e : ent) lg->AddEntry(e.first, e.second, "l");
    lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextFont(42); lg->SetTextSize(0.032); lg->Draw();
    g.titles(cv, "#it{p}_{T}^{jet} (GeV)", Form("relative uncertainty of %s (%%)", what), 0.038, 0.025, 0.022);
    g.cmsLabel(cv, "Internal", "#sqrt{#it{s}_{NN}} = 5.02 TeV");
    savePdfTight(cv, Form("%s/figures/systematics/sysBreakdown_PF_%s%s%s.pdf", repo, tag, sfx.Data(), ptMax != 200. ? Form("_to%.0f", ptMax) : ""));
  };
  const char *cl[4] = {"PbPb 0-10%", "PbPb 10-30%", "PbPb 30-50%", "PbPb 50-80%"};

  // ---- b purity
  {
    std::vector<Col> cols = {{"pp", fPurPP, "h_sysRel_%s"}};
    for(int c = 4; c >= 1; c--) cols.push_back({cl[c-1], fPurAA, Form("h_sysRel_%%s_C%d", c)});
    figure("bPurity", "#it{b} purity", {
      {"JEU", "jet energy scale"}, {"JER", "jet energy resolution"},
      {"bGS", "#it{g}#rightarrow#it{b#bar{b}} share"}, {"cMult", surfVariantLabels() ? "light yield" : "#it{c} fraction"},
      {"fitRangeHigh", "fit range (upper edge)"},
      {"tplStat", "template MC stat."}, {"fakeMuJet", "uncorrelated #mu-jet pairs"}, {"total", "total"}}, cols, TString(purVar()));
  }
  // ---- inclusive jets per Z and R_Z(incl)
  const std::vector<Src> inclSrc = {
    {"JEU", "jet energy scale"}, {"JER", "jet energy resolution"},
    {"iterations", "unfolding iterations"}, {"unfoldingBias", "unfolding non-closure"}, {"responseStat", "response MC stat."},
    {"fakeJets", "combinatorial jets"}, {"stitchWindow", "trigger stitch window"},
    {"stitchNormStat", "trigger stitch norm. stat."}, {"recoFloor", "reco floor"},
    {"muonEff", "muon efficiency"}, {"NZstat", "#it{N}_{Z} stat."},
    {"zBkg", "#it{Z} background"}, {"massWindow", "#it{Z} mass window"}, {"pubZYield", "published #it{Z} yield"}, {"pubZAcc", "pp |#it{y}_{Z}| < 2.1 fraction"}, {"total", "total"}};
  {
    std::vector<Col> cols = {{"pp", fIncl, "h_sysRel_%s_pp_incl"}};
    for(int c = 4; c >= 1; c--) cols.push_back({cl[c-1], fIncl, Form("h_sysRel_%%s_C%d_incl", c)});
    figure("inclJetsPerZ", "inclusive jets per #it{Z}", inclSrc, cols, inclTag);
    std::vector<Col> rz;
    for(int c = 4; c >= 1; c--) rz.push_back({Form("%s / pp", cl[c-1]), fIncl, Form("h_sysRel_%%s_PbPbOverPP_C%d_incl", c)});
    figure("RZ_incl", "inclusive #it{R}_{Z}", inclSrc, rz, inclTag);
  }
  // ---- muon-tagged b jets per Z, R_Z(b), double ratio
  const std::vector<Src> bSrc = {
    {"JEU", "jet energy scale"}, {"JER", "jet energy resolution"}, {"bPurity", "#it{b} purity, uncorr."}, {"bPurity_bGS", "#it{b} purity: #it{g}#rightarrow#it{b#bar{b}} (corr.)"}, {"bPurity_cMult", surfVariantLabels() ? "#it{b} purity: light yield (corr.)" : "#it{b} purity: #it{c} fraction (corr.)"}, {"bTagWP", "#it{b}-tag working point"},
    {"responseStat", "response MC stat."}, {"iterations", "unfolding iterations"}, {"prior", "unfolding prior"}, {"belowFloorWeight", "below-floor weight"},
    {"unfoldingBias", "unfolding non-closure"}, {"muonTagFreq", "muon-tag frequency"},
    {"muonEff", "muon efficiency"}, {"NZstat", "#it{N}_{Z} stat."},
    {"massWindow", "#it{Z} mass window"}, {"incl_fakeJets", "combinatorial jets"},
    {"incl_stitchWindow", "trigger stitch window"}, {"incl_stitchNormStat", "trigger stitch norm. stat."},
    {"incl_recoFloor", "reco floor"}, {"pubZYield", "published #it{Z} yield"}, {"pubZAcc", "pp |#it{y}_{Z}| < 2.1 fraction"}, {"total", "total"}};
  {
    std::vector<Col> cols = {{"pp", fB, "h_sysRel_%s_b_pp"}};
    for(int c = 4; c >= 1; c--) cols.push_back({cl[c-1], fB, Form("h_sysRel_%%s_b_C%d", c)});
    figure("bJetsPerZ", "#mu-tagged #it{b} jets per #it{Z}", bSrc, cols, bTag);
    std::vector<Col> rz, dr;
    for(int c = 4; c >= 1; c--){
      rz.push_back({Form("%s / pp", cl[c-1]), fB, Form("h_sysRel_%%s_b_PbPbOverPP_C%d", c)});
      dr.push_back({Form("%s / pp", cl[c-1]), fB, Form("h_sysRel_%%s_bOverIncl_PbPbOverPP_C%d", c)});
    }
    figure("RZ_b", "#it{b}-jet #it{R}_{Z}", bSrc, rz, bTag);
    figure("doubleRatio", "#it{R}_{Z}^{#it{b}} / #it{R}_{Z}^{incl}", bSrc, dr, bTag);
  }
}
