// b-jet R_Z (top) and the double ratio R_Z(b) / R_Z(incl) (bottom) on the 7-BIN jet-pT axis against the nominal
// 4 bins, four PbPb classes peripheral to central, 80-200 GeV (2026-11-16, user: "run the b-jet chain in the 7 bins
// as a variant"). Copy of the layout of plotBJetRZ_2DvsNominal_PF.C, which needs both files on the same axis; here
// each file is drawn on its own. Stat. bars, total-systematic boxes (normalization not drawn: common to both).
// Also prints both tables and, per nominal bin, the plain mean of the 7-bin points inside it over the nominal value
// (a rough consistency check: the fine points are not averaged with their yields).
// Inputs: rootFiles/systematics/bJetsPerZSys_PF_MBcaloMatched_JEUppCorr<tag>.root for
//   tagA = _noFlatJEU_fakeGeoOnly_csvV2_2D_noRecoFloor_normSep_noStitch       (nominal)
//   tagB = _noFlatJEU_fakeGeoOnly_csvV2_2D7bins_noRecoFloor_normSep_noStitch  (7 bins; chain run with BJET_PT_7BINS)
// Figure: figures/fullChain_PF_MBcaloMatched/14_fileTag/bJetRZ_doubleRatio_7binsVsNominal.pdf
// Usage, from src/calculateBJetsPerZ/:  root -l -b -q plotBJetRZ_7binsVsNominal_PF.C
#include "../../headers/plotting/plotStyle.h"
#include "../../headers/plotting/gridPads.h"

void plotBJetRZ_7binsVsNominal_PF(TString tagA = "_noFlatJEU_fakeGeoOnly_csvV2_2D_noRecoFloor_normSep_noStitch",
                                  TString tagB = "_noFlatJEU_fakeGeoOnly_csvV2_2D7bins_noRecoFloor_normSep_noStitch",
                                  TString labA = "4 bins", TString labB = "7 bins", TString outName = "bJetRZ_doubleRatio_7binsVsNominal"){
  initPlotStyle();
  const char *repo = "/home/clayton/Analysis/code/bJetRaaAnalysis";
  TFile *f[2] = {TFile::Open(Form("%s/rootFiles/systematics/bJetsPerZSys_PF_MBcaloMatched_JEUppCorr%s.root", repo, tagA.Data())),
                 TFile::Open(Form("%s/rootFiles/systematics/bJetsPerZSys_PF_MBcaloMatched_JEUppCorr%s.root", repo, tagB.Data()))};
  if(!f[0] || !f[1]){ printf("ERROR: inputs\n"); return; }
  const int NC = 4; const char *cTag[NC] = {"C4", "C3", "C2", "C1"}, *cName[NC] = {"50-80%", "30-50%", "10-30%", "0-10%"};
  const int NR = 2; const char *q[NR] = {"b_PbPbOverPP", "bOverIncl_PbPbOverPP"}, *qLab[NR] = {"R_Z(b)", "R_Z(b) / R_Z(incl)"};
  const double ptLo = 80., ptHi = 200.;
  TH1D *hv[2][NR][NC], *hs[2][NR][NC], *ht[2][NR][NC];
  for(int v = 0; v < 2; v++) for(int r = 0; r < NR; r++) for(int c = 0; c < NC; c++){
    hv[v][r][c] = (TH1D*) f[v]->Get(Form("h_nom_%s_%s", q[r], cTag[c])); hs[v][r][c] = (TH1D*) f[v]->Get(Form("h_statRel_%s_%s", q[r], cTag[c])); ht[v][r][c] = (TH1D*) f[v]->Get(Form("h_sysRel_total_%s_%s", q[r], cTag[c]));
    if(!hv[v][r][c] || !hs[v][r][c] || !ht[v][r][c]){ printf("ERROR: %s %s missing in file %d\n", q[r], cTag[c], v); return; }
  }
  auto in = [&](TH1D *h, int i){ return h->GetBinLowEdge(i) >= ptLo - 1e-6 && h->GetBinLowEdge(i + 1) <= ptHi + 1e-6; };
  for(int r = 0; r < NR; r++) for(int v = 0; v < 2; v++){
    printf("\n==== %s, %s: value +- stat%% +- syst%%\n  window   ", qLab[r], v == 0 ? labA.Data() : labB.Data()); for(int c = 0; c < NC; c++) printf(" | %-20s", cName[c]); printf("\n");
    TH1D *a = hv[v][r][0];
    for(int i = 1; i <= a->GetNbinsX(); i++){ if(!in(a, i)) continue;
      printf("  %3.0f-%3.0f ", a->GetBinLowEdge(i), a->GetBinLowEdge(i + 1));
      for(int c = 0; c < NC; c++) printf(" | %.2f +- %4.1f +- %4.1f ", hv[v][r][c]->GetBinContent(i), 100.*hs[v][r][c]->GetBinContent(i), 100.*ht[v][r][c]->GetBinContent(i));
      printf("\n"); }
  }
  for(int r = 0; r < NR; r++){
    printf("\n==== %s: plain mean of the %s points inside each %s bin / the %s value\n  window   ", qLab[r], labB.Data(), labA.Data(), labA.Data()); for(int c = 0; c < NC; c++) printf(" | %-8s", cName[c]); printf("\n");
    TH1D *a = hv[0][r][0];
    for(int i = 1; i <= a->GetNbinsX(); i++){ if(!in(a, i)) continue;
      printf("  %3.0f-%3.0f ", a->GetBinLowEdge(i), a->GetBinLowEdge(i + 1));
      for(int c = 0; c < NC; c++){ double s = 0.; int n = 0; TH1D *b = hv[1][r][c];
        for(int j = 1; j <= b->GetNbinsX(); j++) if(b->GetBinCenter(j) > a->GetBinLowEdge(i) && b->GetBinCenter(j) < a->GetBinLowEdge(i + 1)){ s += b->GetBinContent(j); n++; }
        printf(" | %.3f   ", n && hv[0][r][c]->GetBinContent(i) > 0. ? s/n/hv[0][r][c]->GetBinContent(i) : 0.); }
      printf("\n"); }
  }

  // ---- figure
  double yMax[NR] = {};
  for(int r = 0; r < NR; r++) for(int v = 0; v < 2; v++) for(int c = 0; c < NC; c++) for(int i = 1; i <= hv[v][r][c]->GetNbinsX(); i++) if(in(hv[v][r][c], i))
    yMax[r] = TMath::Max(yMax[r], hv[v][r][c]->GetBinContent(i)*(1. + TMath::Max(ht[v][r][c]->GetBinContent(i), hs[v][r][c]->GetBinContent(i))));
  TCanvas *cv = new TCanvas("cv7", "", 2200, 1300);
  GridPads g = makeGridPads(cv, NC, NR, 0.075, 0.01, 0.10, 0.065);
  const char *hex[2] = {hexData, okabeHex[5]}; const int mark[2] = {markOpenCircle, markFilledSquare}; const double msz[2] = {3.0, 2.3};
  TH1D *leg[2] = {};
  for(int r = 0; r < NR; r++) for(int c = 0; c < NC; c++){
    g.at(c, r)->cd();
    TH1F *fr = gPad->DrawFrame(ptLo, 0., ptHi, TMath::Max(r == 0 ? 1.7 : 2.9, 1.25*yMax[r]));   // headroom for the class label and the legend g.styleFrame(fr, c, r, 0.055);
    fr->GetYaxis()->SetLabelFont(43); fr->GetYaxis()->SetLabelSize(c == 0 ? 34 : 0); fr->GetYaxis()->SetNdivisions(504);
    fr->GetXaxis()->SetLabelFont(43); fr->GetXaxis()->SetLabelSize(r == NR - 1 ? 34 : 0); fr->GetXaxis()->SetNdivisions(505);
    TLine *one = new TLine(ptLo, 1., ptHi, 1.); one->SetLineStyle(2); one->SetLineWidth(2); one->SetLineColor(kGray + 2); one->Draw();
    for(int v = 0; v < 2; v++){
      TH1D *h = (TH1D*) hv[v][r][c]->Clone(Form("d%d_%d_%d", v, r, c)); h->SetDirectory(nullptr);
      for(int i = 1; i <= h->GetNbinsX(); i++){
        if(!in(h, i) || h->GetBinContent(i) <= 0.){ h->SetBinContent(i, -999.); h->SetBinError(i, 0.); continue; }
        const double y = h->GetBinContent(i), s = y*ht[v][r][c]->GetBinContent(i);
        TBox *bx = new TBox(h->GetBinLowEdge(i), y - s, h->GetBinLowEdge(i + 1), y + s);
        bx->SetFillColorAlpha(TColor::GetColor(v == 0 ? "#8C8C8C" : hex[v]), v == 0 ? 0.30 : 0.25); bx->SetLineWidth(0); bx->Draw();
        h->SetBinError(i, y*hs[v][r][c]->GetBinContent(i));
      }
      styleH(h, hex[v], mark[v], msz[v]); h->SetLineWidth(3); h->Draw("E1 X0 SAME");
      if(!leg[v]) leg[v] = h;
    }
    gPad->RedrawAxis();
    if(r == 0) g.text(c, r, 0.10, cName[c], 0.065, 62);
    if(r == 0 && c == 0){
      const double fL = gPad->GetLeftMargin(), fW = 1. - fL - gPad->GetRightMargin(), fT = 1. - gPad->GetTopMargin();
      TLegend *lg = new TLegend(fL + 0.52*fW, fT - 0.27, fL + fW, fT - 0.04);
      lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextFont(43); lg->SetTextSize(34);
      lg->AddEntry(leg[0], labA, "pe"); lg->AddEntry(leg[1], labB, "pe"); lg->Draw();
    }
  }
  cv->cd();
  { TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.040); t.SetTextAlign(21);
    t.DrawLatex(g.lm + (1. - g.lm - g.rm)/2., 0.014, "#it{p}_{T}^{jet} (GeV)");
    t.SetTextAngle(90);
    const double mid = g.bm + 0.5*(1. - g.bm - g.tm);
    t.DrawLatex(0.022, mid + 0.25*(1. - g.bm - g.tm), "#it{R}_{Z} (b jets)");
    t.DrawLatex(0.022, mid - 0.25*(1. - g.bm - g.tm), "#it{R}_{Z}(b) / #it{R}_{Z}(incl)"); }
  g.cmsLabel(cv, "Internal", "pp, PbPb (5.02 TeV)   PF jets, |#eta| < 1.6", 1.0);
  const TString dir = Form("%s/figures/fullChain_PF_MBcaloMatched/14_fileTag/", repo); gSystem->mkdir(dir, true);
  savePdfTight(cv, dir + outName + ".pdf");
}
