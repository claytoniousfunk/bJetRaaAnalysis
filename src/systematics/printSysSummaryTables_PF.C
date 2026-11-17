// Systematic uncertainties by source, as LaTeX tables for talks (2026-10-20,
// user): b purity, inclusive jets per Z, muon-tagged b jets per Z, R_Z(incl),
// R_Z(b), and the double ratio R_Z(b) / R_Z(incl). Reads the nominal PF
// systematics files and recomputes nothing. Each cell is the min-max of the
// relative uncertainty (%) over the reported jet-pT bins, 80-200 GeV (as
// plotRZ_bVsInclusive_PF.C / plotBPurity_onePanel_PF.C);
// "--" where the source does not enter. Columns pp, 50-80, 30-50, 10-30, 0-10%
// (peripheral to central, as plotSysBreakdown_PF.C; 2026-10-20, user).
//
// inputs  rootFiles/systematics/bPuritySys_total_PFJets_pp.root   (h_sysRel_<src>)
//         rootFiles/systematics/bPuritySys_PFJets_PbPb.root       (h_sysRel_<src>_C<n>)
//         rootFiles/systematics/inclJetsPerZSys_PF_MBcaloMatched_JEUppCorr.root
//           (h_sysRel_<src>_{pp,C<n>,PbPbOverPP_C<n>}_incl)
//         rootFiles/systematics/bJetsPerZSys_PF_MBcaloMatched_JEUppCorr_reverseUnfold.root
//           (h_sysRel_<src>_b_{pp,C<n>,PbPbOverPP_C<n>}, h_sysRel_<src>_bOverIncl_PbPbOverPP_C<n>) -- the nominal reverse unfolding
// Usage: root -l -b -q printSysSummaryTables_PF.C   (from src/systematics/)
//        root -l -b -q 'printSysSummaryTables_PF.C("_noFlatJEU")'   inclusive JEU not flattened
//        root -l -b -q 'printSysSummaryTables_PF.C("_noFlatJEU","_noFlatJEU_noUnfBias")'   b side too
// Output: figures/systematics/tables/sysSummary_PF_<quantity>.tex, one tabular per file (printed too)
#include "TFile.h"
#include "../../headers/functions/purityVariant.h"   // purity files of the purity variant (BJET_PURITY_VARIANT; 2026-10-28)
#include "TH1.h"
#include "TString.h"
#include "TSystem.h"
#include <vector>
#include <fstream>

// inclTag: suffix of the inclusive file (e.g. "_noFlatJEU"), also added to the
// names of the two inclusive tables
// bTag: suffix of the b-jet file after _reverseUnfold (e.g. "_noFlatJEU_noUnfBias"),
// also added to the names of the b-jet, R_Z(b) and double-ratio outputs
// simple-surface variant (BJET_PURITY_VARIANT=csvV2__surf): its light-yield source is stored under the name cMult
inline bool surfVariantLabels(){ const char *e = std::getenv("BJET_PURITY_VARIANT"); return e && TString(e) == "csvV2__surf"; }
void printSysSummaryTables_PF(TString inclTag = "", TString bTag = ""){
  const char *repo = "/home/clayton/Analysis/code/bJetRaaAnalysis";
  const double ptMin = 80., ptMax = 200.;
  TFile *fPurPP = TFile::Open(Form("%s/rootFiles/systematics/bPuritySys_total_PFJets_pp%s.root", repo, ppPurVar()));
  TFile *fPurAA = TFile::Open(Form("%s/rootFiles/systematics/bPuritySys_PFJets_PbPb%s.root", repo, purVar()));
  TFile *fIncl  = TFile::Open(Form("%s/rootFiles/systematics/inclJetsPerZSys_PF_MBcaloMatched_JEUppCorr%s.root", repo, inclTag.Data()));
  // bTag starting with "_bayes" (or "bayes") (2026-10-21: Bayes, MC prior, is the b-jet nominal again): the Bayes
  // file bJetsPerZSys_PF_MBcaloMatched_JEUppCorr<rest>.root instead of the reverse-unfolding one
  const bool bayesB = bTag.BeginsWith("_bayes") || bTag.BeginsWith("bayes");
  const TString bRest = bayesB ? TString(bTag(bTag.Index("bayes") + 5, bTag.Length())) : bTag;
  TFile *fB     = TFile::Open(Form("%s/rootFiles/systematics/bJetsPerZSys_PF_MBcaloMatched_JEUppCorr%s%s.root", repo, bayesB ? "" : "_reverseUnfold", bRest.Data()));
  if(!fPurPP || !fPurAA || !fIncl || !fB){ printf("ERROR: a systematics file did not open\n"); return; }

  // "min-max" in % over [ptMin, ptMax); "" when the histogram is missing or zero there
  auto cell = [&](TFile *f, TString name) -> TString {
    TH1 *h = nullptr; f->GetObject(name, h);
    if(!h) return "";
    double lo = 1e9, hi = -1.;
    for(int i = 1; i <= h->GetNbinsX(); i++){
      if(h->GetBinLowEdge(i) < ptMin - 1e-6 || h->GetBinLowEdge(i+1) > ptMax + 1e-6) continue;
      const double v = 100.*fabs(h->GetBinContent(i));
      lo = std::min(lo, v); hi = std::max(hi, v);
    }
    if(hi < 0.05) return "";
    return (TString(Form("%.1f", lo)) == Form("%.1f", hi)) ? TString(Form("%.1f", hi)) : TString(Form("%.1f--%.1f", lo, hi));
  };

  struct Src { const char *key, *label; };
  // one column: file + histogram-name pattern, %s = source key
  struct Col { const char *head; TFile *f; TString fmt; };
  auto table = [&](TString tag, const std::vector<Src> &srcs, const std::vector<Col> &cols){
    TString out = Form("\\begin{tabular}{l%s}\n\\hline\n source", TString('c', cols.size()).Data());
    for(auto &c : cols) out += Form(" & %s", c.head);
    out += " \\\\\n\\hline\n";
    for(auto &s : srcs){
      std::vector<TString> v; bool any = false;
      for(auto &c : cols){ v.push_back(cell(c.f, Form(c.fmt.Data(), s.key))); any |= v.back() != ""; }
      if(!any) continue;
      if(TString(s.key) == "total") out += "\\hline\n";
      out += TString(s.key) == "total" ? Form(" \\textbf{%s}", s.label) : Form(" %s", s.label);
      for(auto &x : v) out += Form(" & %s", x == "" ? "--" : (TString(s.key) == "total" ? Form("\\textbf{%s}", x.Data()) : x.Data()));
      out += " \\\\\n";
    }
    // normalization kept apart from the systematics (h_normRel_*, tools/dropSysSources.C asNorm)
    const std::vector<Src> nrm = {{"NZstat", "$N_\\text{Z}$ stat."}, {"zBkg", "$Z$ background"}, {"massWindow", "$Z$ mass window"}, {"muonEff", "muon efficiency"}, {"pubZYield", "published $Z$ yield"}, {"pubZAcc", "pp $|y_\\text{Z}| < 2.1$ fraction"}, {"total", "normalization total"}};
    bool hdr = false;
    for(auto &s : nrm){
      std::vector<TString> v; bool any = false;
      for(auto &c : cols){ TString nf = c.fmt; nf.ReplaceAll("h_sysRel_", "h_normRel_"); v.push_back(cell(c.f, Form(nf.Data(), s.key))); any |= v.back() != ""; }
      if(!any) continue;
      if(!hdr){ out += Form("\\hline\n \\multicolumn{%zu}{l}{\\textit{normalization, separate (flat in $p_\\text{T}$)}} \\\\\n", cols.size() + 1); hdr = true; }
      const bool t = TString(s.key) == "total";
      out += t ? Form(" \\textbf{%s}", s.label) : Form(" %s", s.label);
      for(auto &x : v) out += Form(" & %s", x == "" ? "--" : (t ? Form("\\textbf{%s}", x.Data()) : x.Data()));
      out += " \\\\\n";
    }
    out += "\\hline\n\\end{tabular}\n";
    gSystem->mkdir(Form("%s/figures/systematics/tables", repo), true);
    std::ofstream(Form("%s/figures/systematics/tables/sysSummary_PF_%s.tex", repo, tag.Data())) << out.Data();
    printf("==== %s\n%s\n", tag.Data(), out.Data());
  };
  const char *cl[4] = {"0--10\\%", "10--30\\%", "30--50\\%", "50--80\\%"};

  // ---- b purity
  {
    std::vector<Col> cols = {{"pp", fPurPP, "h_sysRel_%s"}};
    for(int c = 4; c >= 1; c--) cols.push_back({cl[c-1], fPurAA, Form("h_sysRel_%%s_C%d", c)});
    table(TString("bPurity") + purVar(), {
      {"JEU", "jet energy scale"}, {"JER", "jet energy resolution"},
      {"bGS", "$g\\to b\\bar{b}$ share"}, {"cMult", surfVariantLabels() ? "light yield" : "$c$ fraction"},
      {"fitRangeHigh", "fit range (upper edge)"},
      {"tplStat", "template MC stat."}, {"fakeMuJet", "uncorrelated $\\mu$-jet pairs"}, {"total", "total"}}, cols);
  }
  // ---- inclusive jets per Z and R_Z(incl)
  const std::vector<Src> inclSrc = {
    {"JEU", "jet energy scale"}, {"JER", "jet energy resolution"},
    {"iterations", "unfolding iterations"}, {"unfoldingBias", "unfolding non-closure"}, {"responseStat", "response MC stat."},
    {"fakeJets", "combinatorial jets"}, {"stitchWindow", "trigger stitch window"},
    {"stitchNormStat", "trigger stitch norm.\\ stat."}, {"recoFloor", "reco floor"},
    {"muonEff", "muon efficiency"}, {"NZstat", "$N_\\text{Z}$ stat."},
    {"zBkg", "$Z$ background"}, {"massWindow", "$Z$ mass window"}, {"pubZYield", "published $Z$ yield"}, {"pubZAcc", "pp $|y_\\text{Z}| < 2.1$ fraction"}, {"total", "total"}};
  {
    std::vector<Col> cols = {{"pp", fIncl, "h_sysRel_%s_pp_incl"}};
    for(int c = 4; c >= 1; c--) cols.push_back({cl[c-1], fIncl, Form("h_sysRel_%%s_C%d_incl", c)});
    table("inclJetsPerZ" + inclTag, inclSrc, cols);
    std::vector<Col> rz;
    for(int c = 4; c >= 1; c--) rz.push_back({cl[c-1], fIncl, Form("h_sysRel_%%s_PbPbOverPP_C%d_incl", c)});
    table("RZ_incl" + inclTag, inclSrc, rz);
  }
  // ---- muon-tagged b jets per Z and R_Z(b)
  const std::vector<Src> bSrc = {
    // logical order (2026-11-16, user), the same as the figures (sysStyle.h): jet energy, b purity,
    // unfolding, muon tagging, inclusive-side sources, normalization
    {"JEU", "jet energy scale"}, {"JER", "jet energy resolution"},
    {"bPurity_bGS", "$b$ purity: $g\\to b\\bar{b}$ (corr.)"}, {"bPurity_cMult", surfVariantLabels() ? "$b$ purity: light yield (corr.)" : "$b$ purity: $c$ fraction (corr.)"}, {"bPurity", "$b$ purity, other (uncorr.)"}, {"bTagWP", "$b$-tag working point"},
    {"iterations", "unfolding iterations"}, {"unfoldingBias", "unfolding non-closure"}, {"responseStat", "response MC stat."}, {"prior", "unfolding prior"}, {"belowFloorWeight", "below-floor weight"},
    {"muonEff", "muon efficiency"}, {"muonTagFreq", "muon-tag frequency"},
    {"NZstat", "$N_\\text{Z}$ stat."},
    {"massWindow", "$Z$ mass window"}, {"incl_fakeJets", "combinatorial jets"},
    {"incl_stitchWindow", "trigger stitch window"}, {"incl_stitchNormStat", "trigger stitch norm.\\ stat."},
    {"incl_recoFloor", "reco floor"}, {"pubZYield", "published $Z$ yield"}, {"pubZAcc", "pp $|y_\\text{Z}| < 2.1$ fraction"}, {"total", "total"}};
  {
    std::vector<Col> cols = {{"pp", fB, "h_sysRel_%s_b_pp"}};
    for(int c = 4; c >= 1; c--) cols.push_back({cl[c-1], fB, Form("h_sysRel_%%s_b_C%d", c)});
    table("bJetsPerZ" + bTag, bSrc, cols);
    std::vector<Col> rz;
    for(int c = 4; c >= 1; c--) rz.push_back({cl[c-1], fB, Form("h_sysRel_%%s_b_PbPbOverPP_C%d", c)});
    table("RZ_b" + bTag, bSrc, rz);
    std::vector<Col> dr;
    for(int c = 4; c >= 1; c--) dr.push_back({cl[c-1], fB, Form("h_sysRel_%%s_bOverIncl_PbPbOverPP_C%d", c)});
    table("doubleRatio" + bTag, bSrc, dr);
  }
}
