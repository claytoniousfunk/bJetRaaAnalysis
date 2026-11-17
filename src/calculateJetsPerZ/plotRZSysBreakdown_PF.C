// Systematic uncertainties of R_Z by source, inclusive and b jets, two classes,
// CMS style, for talks (2026-10-14, user; the figure form of the backup table).
// Reads the systematics files and recomputes nothing: the chain's sources are
// grouped as in the deck's table, added in quadrature within a group:
//   jet energy scale      JEU
//   jet energy resolution JER
//   b purity              bPurity                       (b only)
//   unfolding             iterations, responseStat, unfoldingBias
//   combinatorial jets    fakeJets                      (incl only)
//   muon tag frequency    muonTagFreq                   (b only)
//   normalization         h_normRel_total (N_Z stat, Z background, mass window, muon efficiency),
//                         a box kept apart from the band (dropSysSources asNorm), dashed grey
// 2x2 panels sharing borders (gridPads.h): top inclusive, bottom b jets; left
// 50-80%, right 0-10% (peripheral to central); total of the band in black. 80-200 GeV.
// 2026-10-21: inputs moved to the deck chain (Bayes b-jet unfolding, JER pp-PbPb uncorrelated);
// stitching and reco floor are not in the deck band.
//
// inputs  rootFiles/systematics/inclJetsPerZSys_PF_MBcaloMatched_JEUppCorr<tag>.root (h_sysRel_<src>_PbPbOverPP_C<n>_incl)
//         rootFiles/systematics/bJetsPerZSys_PF_MBcaloMatched_JEUppCorr<tag>.root (h_sysRel_<src>_b_PbPbOverPP_C<n>)
//         tag = _noFlatJEU_fakeGeoOnly_noRecoFloor_normSep_noStitch by default
// Usage: root -l -b -q plotRZSysBreakdown_PF.C   (from src/calculateJetsPerZ/)
// Output: figures/systematics/RZSysBreakdown_PF.pdf
#include "TFile.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TMath.h"
#include "../../headers/plotting/plotStyle.h"
#include "../../headers/plotting/gridPads.h"
#include "../../headers/plotting/sysStyle.h"   // colours / line styles shared with every systematics figure (2026-11-16)

void plotRZSysBreakdown_PF(TString tag = "_noFlatJEU_fakeGeoOnly_noRecoFloor_normSep_noStitch", TString outSfx = ""){   // outSfx: appended to the figure name (2026-11-16)
  const char *repo = "/home/clayton/Analysis/code/bJetRaaAnalysis";
  TFile *fi = TFile::Open(Form("%s/rootFiles/systematics/inclJetsPerZSys_PF_MBcaloMatched_JEUppCorr%s.root", repo, tag.Data()));
  TFile *fb = TFile::Open(Form("%s/rootFiles/systematics/bJetsPerZSys_PF_MBcaloMatched_JEUppCorr%s.root", repo, tag.Data()));
  if(!fi || !fb){ printf("ERROR: a systematics file did not open\n"); return; }
  initPlotStyle();
  // group colours and line styles = the sysStyle.h style of the group's leading source (2026-11-16)
  struct Grp { const char *label; std::vector<const char*> si, sb; const char *styleKey; };
  const std::vector<Grp> grp = {
    {"JES",                {"JEU"}, {"JEU"}, "JEU"},
    {"JER",                {"JER"}, {"JER"}, "JER"},
    {"#it{b} purity",      {}, {"bPurity", "bPurity_bGS", "bPurity_cMult"}, "bPurity"},   // split sources since 2026-10-28, quadrature
    {"#it{b}-tag working point", {}, {"bTagWP"}, "bTagWP"},   // CSVv2-loose era only; absent in the 2D-fit files
    {"unfolding",          {"iterations","responseStat","unfoldingBias"}, {"iterations","responseStat","unfoldingBias"}, "iterations"},
    {"muon-tag frequency", {}, {"muonTagFreq"}, "muonTagFreq"},
    {"combinatorial jets", {"fakeJets"}, {}, "fakeJets"},
  };
  std::vector<bool> present(grp.size(), false);
  const double ptMin = 80., ptMax = 200., yMax = 40.;
  const int cls[2] = {4, 1}; const char *clsLab[2] = {"50-80%", "0-10%"};
  // a group's relative uncertainty vs pT, in %, quadrature over its sources; nullptr when the lane has no source
  auto group = [&](TFile *f, const std::vector<const char*> &srcs, const char *fmt, int c, const char *name) -> TH1D* {
    if(srcs.empty()) return nullptr;
    TH1D *o = nullptr;
    for(auto s : srcs){
      TH1D *h = nullptr; f->GetObject(Form(fmt, s, c), h);
      if(!h) continue;   // source absent in this file (e.g. bTagWP after 2026-11-16)
      if(!o){ o = (TH1D*) h->Clone(name); o->SetDirectory(nullptr); o->Reset(); }
      for(int b = 1; b <= o->GetNbinsX(); b++) o->SetBinContent(b, hypot(o->GetBinContent(b), h->GetBinContent(b)));
    }
    if(o) for(int b = 1; b <= o->GetNbinsX(); b++){ o->SetBinContent(b, 100.*o->GetBinContent(b)); o->SetBinError(b, 0.); }
    return o;
  };

  TCanvas *cv = new TCanvas("cvSys", "", 1100, 950);
  GridPads g = makeGridPads(cv, 2, 2, 0.10, 0.03, 0.09, 0.07);
  for(int row = 0; row < 2; row++) for(int col = 0; col < 2; col++){
    const bool bLane = row == 1; const int c = cls[col];
    g.at(col, row)->cd();
    TH1F *fr = gPad->DrawFrame(ptMin, 0., ptMax, yMax); g.styleFrame(fr, col, row, 0.052);
    std::vector<TH1D*> drawn;
    for(size_t k = 0; k < grp.size(); k++){
      TH1D *h = bLane ? group(fb, grp[k].sb, "h_sysRel_%s_b_PbPbOverPP_C%d", c, Form("g%zu_%d_%d", k, row, col))
                      : group(fi, grp[k].si, "h_sysRel_%s_PbPbOverPP_C%d_incl", c, Form("g%zu_%d_%d", k, row, col));
      if(!h) continue;
      h->GetXaxis()->SetRangeUser(ptMin, ptMax - 1e-6);
      double mx = 0.; for(int b = 1; b <= h->GetNbinsX(); b++) if(h->GetBinLowEdge(b) >= ptMin - 1e-6 && h->GetBinLowEdge(b) < ptMax - 1e-6) mx = std::max(mx, h->GetBinContent(b));
      if(mx < 0.05) continue;
      present[k] = true;
      const SysStyle st = sysStyle(grp[k].styleKey);
      styleLine(h, st.hex, 4); h->SetLineStyle(st.lineStyle);
      drawn.push_back(h);
    }
    TH1D *tot = bLane ? group(fb, {"total"}, "h_sysRel_%s_b_PbPbOverPP_C%d", c, Form("tot_%d_%d", row, col))
                      : group(fi, {"total"}, "h_sysRel_%s_PbPbOverPP_C%d_incl", c, Form("tot_%d_%d", row, col));
    tot->GetXaxis()->SetRangeUser(ptMin, ptMax - 1e-6); tot->SetFillColorAlpha(kGray + 1, 0.35); tot->SetLineColor(kBlack); tot->SetLineWidth(4); tot->Draw("HIST same");
    for(TH1D *h : drawn) h->Draw("HIST same");
    TH1D *nrm = bLane ? group(fb, {"total"}, "h_normRel_%s_b_PbPbOverPP_C%d", c, Form("nrm_%d_%d", row, col))
                      : group(fi, {"total"}, "h_normRel_%s_PbPbOverPP_C%d_incl", c, Form("nrm_%d_%d", row, col));
    if(nrm){ nrm->GetXaxis()->SetRangeUser(ptMin, ptMax - 1e-6); nrm->SetLineColor(kBlack); nrm->SetLineWidth(4); nrm->SetLineStyle(2); nrm->Draw("HIST same"); }
    for(int b = 1; b <= tot->GetNbinsX(); b++) if(tot->GetBinCenter(b) > ptMin && tot->GetBinCenter(b) < ptMax)
      printf("  %s %-6s %5.0f GeV: total %5.1f%%  norm %4.1f%%\n", bLane ? "b   " : "incl", clsLab[col], tot->GetBinLowEdge(b), tot->GetBinContent(b), nrm ? nrm->GetBinContent(b) : 0.);
    gPad->RedrawAxis();
    g.text(col, row, 0.10, Form("%s, %s / pp", bLane ? "#it{b} jets" : "inclusive jets", clsLab[col]), 0.060, 62);
  }
  // legend after all panels, so it lists only the groups present; inside the inclusive 50-80% panel
  {
    const int row = 0, col = 0; g.at(col, row)->cd();
    TH1D *tot = new TH1D("lpt", "", 1, 0, 1); tot->SetFillColorAlpha(kGray + 1, 0.35); tot->SetLineColor(kBlack); tot->SetLineWidth(4);
    {   // legend in the inclusive 50-80% panel, whose curves stay below ~12%
      TLegend *l2 = makeLegend(g.frameLeft(0) + 0.04, 0.36, 0.98, 0.74, 0.046);
      l2->AddEntry(tot, "total", "f");
      for(size_t k = 0; k < grp.size(); k++){
        if(!present[k]) continue;
        const SysStyle st = sysStyle(grp[k].styleKey);
        TH1D *proxy = new TH1D(Form("lp%zu", k), "", 1, 0, 1); styleLine(proxy, st.hex, 4); proxy->SetLineStyle(st.lineStyle);
        l2->AddEntry(proxy, grp[k].label, "l");
      }
      TH1D *np = new TH1D("lpn", "", 1, 0, 1); np->SetLineColor(kBlack); np->SetLineWidth(4); np->SetLineStyle(2);
      l2->AddEntry(np, "normalization", "l");   // last: kept apart from the band
      l2->SetNColumns(2); l2->Draw();
    }
  }
  g.titles(cv, "#it{p}_{T}^{jet} (GeV)", "relative uncertainty in #it{R}_{Z} (%)", 0.045, 0.024, 0.022);
  g.cmsLabel(cv, "Internal", "anti-#it{k}_{T} #it{R} = 0.4 PF jets,  #sqrt{#it{s}_{NN}} = 5.02 TeV");
  savePdfTight(cv, Form("%s/figures/systematics/RZSysBreakdown_PF%s.pdf", repo, outSfx.Data()));
}
