// Tilted-truth closure of the inclusive-jet unfolding (fixed p_T^-0.3 since 2026-11-05; data-like k earlier that day), as run inside the chain
// (inclusiveJetsPerZSys_PF.C, dataLikeNopt(), 2026-11-05): per system the truth is
// reweighted to the fitted data tilt, folded, and unfolded with the nominal response.
// Top row: unfolded / truth on the reporting axis at N_opt (the residual is the
// "unfolding non-closure" systematic). Bottom row: B^2, V and MSE over the reporting
// bins from 80 GeV vs the number of iterations, N_opt marked. pp, 50-80 ... 0-10%,
// shared-border 2 x 5 grid, CMS style (AN figure, unfolding section).
//
// Input: the tagged inclusive systematics file (g_dlClosure_*, h_dlClosure_uOverT_*).
// Usage, from src/plots/unfolding/:  root -l -b -q plotDataLikeClosure_PF.C
// Writes figures/unfolding/dataLikeClosure_ratio_PF.pdf (top), dataLikeClosure_mse_PF.pdf (bottom, N <= 15),
//        closureRatio_bJets_PF.pdf and closureMSE_bJets_PF.pdf (the same rows for the muon-tagged b jets, means over 80-200 GeV from
//        rootFiles/unfoldTest/optimalIterations_bJets_PF_MCprior.root, scanned to N = 10; 2026-11-06)

#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TSystem.h"
#include "TRegexp.h"
#include "TMath.h"
#include "TColor.h"
#include "TVectorD.h"
#include "../../../headers/plotting/plotStyle.h"
#include "../../../headers/plotting/gridPads.h"

// one shared-border row of B^2, V and MSE vs N (log y), pp then peripheral to central;
// gB/gV/gM(sys) fetch the graphs, nOpt(sys) the chosen N
template<class F, class G>
void drawMSERow(const char *cname, F graphs, G nOptOf, double xMax, const char *yTitle, const char *cmsText, const TString &out)
{
  const char *sys[5] = {"pp", "C4", "C3", "C2", "C1"};
  const char *lab[5] = {"pp", "PbPb 50-80%", "PbPb 30-50%", "PbPb 10-30%", "PbPb 0-10%"};
  TCanvas *c2 = new TCanvas(cname, "", 2000, 650);
  GridPads g2 = makeGridPads(c2, 5, 1, 0.065, 0.012, 0.17, 0.10);
  double yLo = 1e30, yHi = 0.;
  for(int ix = 0; ix < 5; ix++){ TGraph *gs[3]; graphs(sys[ix], gs);
    for(int q = 0; q < 3; q++) for(int i = 0; i < gs[q]->GetN(); i++) if(gs[q]->GetX()[i] <= xMax && gs[q]->GetY()[i] > 0.){
      yLo = TMath::Min(yLo, gs[q]->GetY()[i]); yHi = TMath::Max(yHi, gs[q]->GetY()[i]); } }
  yLo = TMath::Power(10., TMath::Floor(TMath::Log10(yLo))); yHi = TMath::Power(10., TMath::Ceil(TMath::Log10(yHi)) + 1.0);   // room for the labels above the curves
  TGraph *legG[3] = {};
  for(int ix = 0; ix < 5; ix++){
    TGraph *gs[3]; graphs(sys[ix], gs); const int nOpt = nOptOf(sys[ix]);
    g2.at(ix, 0)->cd(); gPad->SetLogy();
    TH1F *f2 = gPad->DrawFrame(0.5, yLo, xMax + 0.5, yHi); g2.styleFrame(f2, ix, 0, 0.075);
    // the SAME x labels in every panel (2026-11-16, user: "I would rather have some clipping than the asymmetry"): the axis'
    // own labels are off (styleFrame hides the last one in every panel but the right-most, and the first used to be dropped
    // in the others) and the labels are drawn by hand at the same values in each panel
    f2->GetXaxis()->SetLabelSize(0.); f2->GetXaxis()->SetNdivisions(505);
    // ... on the CANVAS, not in the pad, so a label next to a panel edge is not clipped by it
    { TPad *pd = (TPad*) gPad; const double step = xMax > 30. ? 10. : xMax > 12. ? 5. : 2., pl = pd->GetLeftMargin(), pw = 1. - pl - pd->GetRightMargin();
      const double x0 = pd->GetXlowNDC(), w = pd->GetWNDC(), y0 = pd->GetYlowNDC(), hgt = pd->GetHNDC();
      TPad *p0 = g2.at(0, 0); const double ts = 0.075*TMath::Min(p0->GetWw()*p0->GetAbsWNDC(), p0->GetWh()*p0->GetAbsHNDC());   // pixels: the y labels' size, the same in every panel
      c2->cd(); TLatex xl; xl.SetNDC(); xl.SetTextFont(43); xl.SetTextSize(ts); xl.SetTextAlign(23);
      for(double x = step; x <= xMax + 1e-6; x += step) xl.DrawLatex(x0 + w*(pl + pw*(x - 0.5)/xMax), y0 + hgt*(pd->GetBottomMargin() - 0.025), Form("%.0f", x));
      pd->cd(); }
    const char *hx[3] = {okabeHex[5], okabeHex[1], hexData}; const int mk[3] = {markOpenSquare, markOpenCircle, markFilledCircle};
    for(int q = 0; q < 3; q++){
      TGraph *gg = new TGraph();   // points beyond xMax dropped, so the lines stop at the frame
      for(int i = 0; i < gs[q]->GetN(); i++) if(gs[q]->GetX()[i] <= xMax) gg->SetPoint(gg->GetN(), gs[q]->GetX()[i], gs[q]->GetY()[i]);
      gg->SetMarkerStyle(mk[q]); gg->SetMarkerSize(1.5);
      gg->SetMarkerColor(TColor::GetColor(hx[q])); gg->SetLineColor(TColor::GetColor(hx[q])); gg->SetLineWidth(3);
      gg->Draw("PL same"); if(ix == 0) legG[q] = gg;
    }
    TLine *lo = new TLine(nOpt, yLo, nOpt, yHi/300.); lo->SetLineStyle(7); lo->SetLineWidth(3); lo->SetLineColor(kGray + 2); lo->Draw();
    g2.text(ix, 0, 0.10, lab[ix], 0.075);
    g2.text(ix, 0, 0.21, Form("#it{N}_{opt} = %d", nOpt), 0.062);   // under the class label, clear of the curves
    gPad->RedrawAxis();
  }
  g2.at(0, 0)->cd();   // pp's upper right is empty in both rows
  TLegend *lg = makeLegend(0.62, 0.60, 0.97, 0.86, 0.060*g2.sizeScale());
  lg->AddEntry(legG[2], "MSE", "pl"); lg->AddEntry(legG[0], "#it{B}^{2}", "pl"); lg->AddEntry(legG[1], "#it{V}", "pl");
  lg->Draw();
  g2.titles(c2, "Iterations #it{N}", yTitle, 0.060, 0.020, 0.008);
  g2.cmsLabel(c2, "Simulation Internal", cmsText, 1.4);
  savePdfTight(c2, out);
}

void plotDataLikeClosure_PF(const char *in = "../../../rootFiles/systematics/inclJetsPerZSys_PF_MBcaloMatched_JEUppCorr_noFlatJEU_fakeGeoOnly.root",
                            const char *bIn = "../../../rootFiles/unfoldTest/optimalIterations_bJets_PF_MCprior.root",
                            TString bTag = "", double bXMax = 10.)   // bXMax: iterations shown in the b-jet MSE row (a longer scan: optimalIterations_bJets_PF_MCprior.C nMax)
   // 2026-11-16: non-empty = only the two b-jet figures, named ..._bJets_PF<bTag>.pdf (the inclusive
                                                 // and the untagged AN figures are left alone), e.g. bIn ..._MCprior_csvV2_2D7bins.root, "_csvV2_2D7bins"
{
  initPlotStyle();
  const bool bOnly = bTag != "";
  TFile *f = TFile::Open(in);
  if(!f || f->IsZombie()){ printf("ERROR: %s did not open\n", in); return; }
  const char *sys[5] = {"pp", "C4", "C3", "C2", "C1"};
  const char *lab[5] = {"pp", "PbPb 50-80%", "PbPb 30-50%", "PbPb 10-30%", "PbPb 0-10%"};

  // unfolded / truth at N_opt
  TCanvas *c1 = new TCanvas("cDL1", "", 2000, 650);
  GridPads g1 = makeGridPads(c1, 5, 1, 0.065, 0.012, 0.17, 0.10);
  const double rLo = 0.79, rHi = 1.21;
  for(int ix = 0; ix < 5; ix++){
    TH1D *r = nullptr; f->GetObject(Form("h_dlClosure_uOverT_%s", sys[ix]), r);
    if(!r){ printf("ERROR: scan missing for %s (rerun inclusiveJetsPerZSys_PF.C)\n", sys[ix]); return; }
    int nOpt = 0; sscanf(strstr(r->GetTitle(), "N_opt = ") + 8, "%d", &nOpt);
    double k = 0.; sscanf(strstr(r->GetTitle(), "k = ") + 4, "%lf", &k);
    printf("%-3s N_opt %d  k %+.3f  u/t - 1 (%%):", sys[ix], nOpt, k);
    for(int b = 1; b <= r->GetNbinsX(); b++) if(r->GetBinLowEdge(b) >= 80. - 1e-6) printf(" %.0f-%.0f:%+.1f", r->GetBinLowEdge(b), r->GetXaxis()->GetBinUpEdge(b), 100.*(r->GetBinContent(b) - 1.));
    printf("\n");
    g1.at(ix, 0)->cd();
    TH1F *fr = gPad->DrawFrame(80., rLo, 500., rHi); g1.styleFrame(fr, ix, 0, 0.075);
    fr->GetXaxis()->SetLabelSize(0.075*g1.padW(0)/g1.padW(ix)); fr->GetXaxis()->SetNdivisions(505);
    TLine *one = new TLine(80., 1., 500., 1.); one->SetLineStyle(2); one->SetLineWidth(2); one->SetLineColor(kGray + 2); one->Draw();
    TH1D *h = (TH1D*) r->Clone(Form("rr_%d", ix));
    for(int b = 1; b <= h->GetNbinsX(); b++) if(h->GetBinLowEdge(b) < 80. - 1e-6) h->SetBinContent(b, -1.);
    styleH(h, hexData, markFilledCircle, 2.0); h->SetLineWidth(3); h->Draw("P same");
    g1.text(ix, 0, 0.10, lab[ix], 0.075);
    g1.text(ix, 0, 0.21, Form("#it{N}_{opt} = %d", nOpt), 0.062);
    gPad->RedrawAxis();
  }
  g1.titles(c1, "Generator-level jet #it{p}_{T} (GeV)", "Unfolded / truth", 0.060, 0.020, 0.030);
  g1.cmsLabel(c1, "Simulation Internal", "Inclusive jets, truth #times #it{p}_{T}^{+0.3} (clamped below 40 GeV), at #it{N}_{opt},  #sqrt{#it{s}_{NN}} = 5.02 TeV", 1.4);
  gSystem->mkdir("../../../figures/unfolding", kTRUE);
  if(!bOnly) savePdfTight(c1, "../../../figures/unfolding/dataLikeClosure_ratio_PF.pdf");

  // B^2, V, MSE vs N up to 15 iterations (user, 2026-11-06): inclusive (sums, from this file) and
  // muon-tagged b jets (means over 80-200 GeV, optimalIterations_bJets_PF_MCprior.C, scanned to N = 10)
  auto inclG = [&](const char *s, TGraph **g){ f->GetObject(Form("g_dlClosure_bias2_%s", s), g[0]); f->GetObject(Form("g_dlClosure_var_%s", s), g[1]); f->GetObject(Form("g_dlClosure_mse_%s", s), g[2]); };
  auto inclN = [&](const char *s){ TH1D *r = nullptr; f->GetObject(Form("h_dlClosure_uOverT_%s", s), r); int n = 0; sscanf(strstr(r->GetTitle(), "N_opt = ") + 8, "%d", &n); return n; };
  if(!bOnly) drawMSERow("cDL2", inclG, inclN, 15., "Sum over bins from 80 GeV",
             "Inclusive jets, truth #times #it{p}_{T}^{+0.3} (clamped below 40 GeV),  #sqrt{#it{s}_{NN}} = 5.02 TeV", "../../../figures/unfolding/dataLikeClosure_mse_PF.pdf");
  TFile *fb = TFile::Open(bIn);
  if(!fb || fb->IsZombie()){ printf("ERROR: %s did not open\n", bIn); return; }
  TVectorD *nb = nullptr; fb->GetObject("optimalN", nb);
  auto bG = [&](const char *s, TGraph **g){ fb->GetObject(Form("g_bias2_%s", s), g[0]); fb->GetObject(Form("g_var_%s", s), g[1]); fb->GetObject(Form("g_mse_%s", s), g[2]); };
  auto bN = [&](const char *s){ const char *t[5] = {"pp", "C1", "C2", "C3", "C4"}; for(int i = 0; i < 5; i++) if(!strcmp(s, t[i])) return (int) (*nb)[i]; return 0; };
  // b-jet unfolded / truth at N_opt on the reported bins (h_uOverT_<sys>, written since 2026-11-06)
  {
    TCanvas *cb = new TCanvas("cDLbR", "", 2000, 650);
    GridPads gb = makeGridPads(cb, 5, 1, 0.065, 0.012, 0.17, 0.10);
    for(int ix = 0; ix < 5; ix++){
      TH1D *r = nullptr; fb->GetObject(Form("h_uOverT_%s", sys[ix]), r);
      if(!r){ printf("ERROR: h_uOverT_%s missing (rerun optimalIterations_bJets_PF_MCprior.C)\n", sys[ix]); return; }
      printf("b %-3s N_opt %d  u/t - 1 (%%):", sys[ix], bN(sys[ix]));
      for(int b = 1; b <= r->GetNbinsX(); b++) printf(" %.0f-%.0f:%+.1f", r->GetBinLowEdge(b), r->GetXaxis()->GetBinUpEdge(b), 100.*(r->GetBinContent(b) - 1.));
      printf("\n");
      gb.at(ix, 0)->cd();
      TH1F *fr = gPad->DrawFrame(80., rLo, 200., rHi); gb.styleFrame(fr, ix, 0, 0.075);
      fr->GetXaxis()->SetLabelSize(0.075*gb.padW(0)/gb.padW(ix)); fr->GetXaxis()->SetNdivisions(505);
      TLine *one = new TLine(80., 1., 200., 1.); one->SetLineStyle(2); one->SetLineWidth(2); one->SetLineColor(kGray + 2); one->Draw();
      TH1D *h = (TH1D*) r->Clone(Form("rb_%d", ix));
      for(int b = 1; b <= h->GetNbinsX(); b++) h->SetBinError(b, 0.);   // the residual only, as the inclusive row
      styleH(h, hexData, markFilledCircle, 2.0); h->SetLineWidth(3); h->Draw("P same");
      gb.text(ix, 0, 0.10, lab[ix], 0.075);
      gb.text(ix, 0, 0.21, Form("#it{N}_{opt} = %d", bN(sys[ix])), 0.062);
      gPad->RedrawAxis();
    }
    gb.titles(cb, "Generator-level jet #it{p}_{T} (GeV)", "Unfolded / truth", 0.060, 0.020, 0.030);
    gb.cmsLabel(cb, "Simulation Internal", "Muon-tagged #it{b} jets, truth #times #it{p}_{T}^{+0.3} (clamped below 40 GeV), at #it{N}_{opt},  #sqrt{#it{s}_{NN}} = 5.02 TeV", 1.4);
    savePdfTight(cb, "../../../figures/unfolding/closureRatio_bJets_PF" + bTag + ".pdf");
  }
  drawMSERow("cDLb", bG, bN, bXMax, "Mean over 80-200 GeV",
             "Muon-tagged #it{b} jets, truth #times #it{p}_{T}^{+0.3} (clamped below 40 GeV),  #sqrt{#it{s}_{NN}} = 5.02 TeV", "../../../figures/unfolding/closureMSE_bJets_PF" + bTag + ".pdf");
}
