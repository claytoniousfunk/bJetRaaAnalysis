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
// Since 2026-11-16 one 1 x N row per quantity drawn by drawSysRow() (headers/plotting/sysStyle.h):
// fixed colour + line style per source across all systematics figures, legends inside the panels.
// Also bFraction, RCP_b, RCP_bOverIncl (the appendix figures bJetsPerZSys_PF_PbPb.C draws).
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
#include "../../headers/plotting/sysStyle.h"

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

  struct Src { const char *key, *label; };   // label nullptr = the sysStyle.h label
  struct Col { TString head; TFile *f; TString fmt; };   // fmt: %s = source key (TString: Form() buffers rotate)
  int nFig = 0;
  // one 1 x N row per quantity, drawn by drawSysRow() (headers/plotting/sysStyle.h, 2026-11-16):
  // fixed colour + line style per source in every figure, total as a grey band, legends inside the
  // left-most panels. Panels pp, 50-80 ... 0-10% (ratios: the classes only).
  auto figure = [&](const char *tag, const char *what, const std::vector<Src> &srcs, const std::vector<Col> &cols, TString sfx = ""){
    nFig++;
    const int nP = cols.size();
    std::vector<const char*> keys, labs;
    for(auto &x : srcs) if(TString(x.key) != "total"){ keys.push_back(x.key); labs.push_back(x.label); }
    std::vector<SysPanel> P(nP);
    for(int p = 0; p < nP; p++){
      P[p].name = cols[p].head;
      TH1 *h = nullptr; cols[p].f->GetObject(Form(cols[p].fmt.Data(), "total"), h); P[p].tot = h;
      if(!h) printf("WARNING: %s missing\n", Form(cols[p].fmt.Data(), "total"));
      for(size_t k = 0; k < keys.size(); k++){ TH1 *hs = nullptr; cols[p].f->GetObject(Form(cols[p].fmt.Data(), keys[k]), hs); P[p].src[k] = hs; }
      // a normalization kept apart (h_normRel_*, tools/dropSysSources.C asNorm) is not drawn
      // here (2026-10-20, user): it goes on the results figures as a box at unity
    }
    drawSysRow(nP, P.data(), keys.size(), keys.data(), labs.data(), ptMin, ptMax, what,
               "anti-#it{k}_{T} #it{R} = 0.4 PF jets,  #sqrt{#it{s}_{NN}} = 5.02 TeV",
               Form("%s/figures/systematics/sysBreakdown_PF_%s%s%s.pdf", repo, tag, sfx.Data(), ptMax != 200. ? Form("_to%.0f", ptMax) : ""));
  };
  const char *cl[4] = {"0-10%", "10-30%", "30-50%", "50-80%"};

  // ---- b purity
  {
    std::vector<Col> cols = {{"pp", fPurPP, "h_sysRel_%s"}};
    for(int c = 4; c >= 1; c--) cols.push_back({cl[c-1], fPurAA, Form("h_sysRel_%%s_C%d", c)});
    figure("bPurity", "b purity", {
      {"JEU", nullptr}, {"JER", nullptr}, {"bGS", nullptr}, {"cMult", surfVariantLabels() ? "light-jet fraction" : nullptr},
      {"fitRangeHigh", nullptr}, {"tplStat", nullptr}, {"fakeMuJet", nullptr}, {"total", nullptr}}, cols, TString(purVar()));
  }
  // ---- inclusive jets per Z and R_Z(incl)
  const std::vector<Src> inclSrc = {
    {"JEU", nullptr}, {"JER", nullptr}, {"iterations", nullptr}, {"unfoldingBias", nullptr}, {"responseStat", nullptr},
    {"fakeJets", nullptr}, {"stitchWindow", nullptr}, {"stitchNormStat", nullptr}, {"recoFloor", nullptr},
    {"muonEff", nullptr}, {"NZstat", nullptr}, {"zBkg", nullptr}, {"massWindow", nullptr}, {"pubZYield", nullptr}, {"pubZAcc", nullptr}, {"total", nullptr}};
  {
    std::vector<Col> cols = {{"pp", fIncl, "h_sysRel_%s_pp_incl"}};
    for(int c = 4; c >= 1; c--) cols.push_back({cl[c-1], fIncl, Form("h_sysRel_%%s_C%d_incl", c)});
    figure("inclJetsPerZ", "inclusive jets per #it{N}_{Z}", inclSrc, cols, inclTag);
    std::vector<Col> rz;
    for(int c = 4; c >= 1; c--) rz.push_back({Form("%s / pp", cl[c-1]), fIncl, Form("h_sysRel_%%s_PbPbOverPP_C%d_incl", c)});
    figure("RZ_incl", "inclusive #it{R}_{Z}", inclSrc, rz, inclTag);
  }
  // ---- muon-tagged b jets per Z, R_Z(b), double ratio
  const std::vector<Src> bSrc = {
    {"JEU", nullptr}, {"JER", nullptr}, {"bPurity_bGS", nullptr}, {"bPurity_cMult", surfVariantLabels() ? "purity: light-jet fraction" : nullptr},
    {"bPurity", nullptr}, {"bTagWP", nullptr},
    {"iterations", nullptr}, {"unfoldingBias", nullptr}, {"responseStat", nullptr}, {"prior", nullptr}, {"belowFloorWeight", nullptr},
    {"muonTagFreq", nullptr}, {"muonEff", nullptr}, {"NZstat", nullptr}, {"massWindow", nullptr}, {"incl_fakeJets", nullptr},
    {"incl_stitchWindow", nullptr}, {"incl_stitchNormStat", nullptr}, {"incl_recoFloor", nullptr}, {"pubZYield", nullptr}, {"pubZAcc", nullptr}, {"total", nullptr}};
  {
    std::vector<Col> cols = {{"pp", fB, "h_sysRel_%s_b_pp"}};
    for(int c = 4; c >= 1; c--) cols.push_back({cl[c-1], fB, Form("h_sysRel_%%s_b_C%d", c)});
    figure("bJetsPerZ", "#it{b} jets per #it{N}_{Z}", bSrc, cols, bTag);
    std::vector<Col> rz, dr;
    for(int c = 4; c >= 1; c--){
      rz.push_back({Form("%s / pp", cl[c-1]), fB, Form("h_sysRel_%%s_b_PbPbOverPP_C%d", c)});
      dr.push_back({Form("%s / pp", cl[c-1]), fB, Form("h_sysRel_%%s_bOverIncl_PbPbOverPP_C%d", c)});
    }
    figure("RZ_b", "#it{R}_{Z}(#it{b})", bSrc, rz, bTag);
    figure("doubleRatio", "#it{R}_{Z}(#it{b}) / #it{R}_{Z}(incl)", bSrc, dr, bTag);
    // the appendix figures bJetsPerZSys_PF_PbPb.C draws in its own older style (2026-11-16): b jets over
    // inclusive jets per system, and R_CP over 50-80% of the b jets and of b / inclusive
    std::vector<Col> fr = {{"pp", fB, "h_sysRel_%s_bOverIncl_pp"}}, rcpB, rcpD;
    for(int c = 4; c >= 1; c--) fr.push_back({cl[c-1], fB, Form("h_sysRel_%%s_bOverIncl_C%d", c)});
    for(int c = 3; c >= 1; c--){
      rcpB.push_back({Form("%s / 50-80%%", cl[c-1]), fB, Form("h_sysRel_%%s_b_RCP_C%d", c)});
      rcpD.push_back({Form("%s / 50-80%%", cl[c-1]), fB, Form("h_sysRel_%%s_bOverIncl_RCP_C%d", c)});
    }
    figure("bFraction", "#it{b} / inclusive jets", bSrc, fr, bTag);
    figure("RCP_b", "#it{R}_{CP}^{Z}(#it{b})", bSrc, rcpB, bTag);
    figure("RCP_bOverIncl", "#it{R}_{CP}^{Z}(#it{b}) / #it{R}_{CP}^{Z}(incl)", bSrc, rcpD, bTag);
  }
}
