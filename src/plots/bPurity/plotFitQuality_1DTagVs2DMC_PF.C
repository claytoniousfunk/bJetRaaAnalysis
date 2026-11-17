// Fit quality of the nominal 1D b-purity fit against the 2D fit, pp and the four PbPb classes in a 1x5
// row, 80-300 GeV (2026-11-12, user: "compare 1d tag to 2d mc in a plot").
//   1D tag  ptRel of the CSVv2 > 0.54 tagged jets, c / light at MC (the nominal cut fit), 25 cells
//   2D MC   ptRel x CSVv2 (4cat: 0.54 / 0.84 / 0.96) of ALL muon-tagged jets, every flavor's tagger
//           shape from MC, c / light at MC (the "MC" model of fitPtRelCSVv2_2D_ppPbPb_PF.C), 100 cells
// Both refitted here with the same code, data, templates and per-event fake-pair subtraction as the
// nominal fits (fitPtRelCSVv2Cut_PF_PbPb.C, buildWindows / buildCls), 0.2 GeV bins, fit 0-5 GeV.
// Top row: deviance / ndf (Poisson likelihood ratio against the saturated model, with the
// Barlow-Beeston penalty for template MC statistics), ndf = cells - 1 (normalization) - 1 (f_b).
// Bottom row: the same as a one-sided Gaussian significance, z = Phi^-1(1 - p), p = P(chi2_ndf > dev),
// because chi2 / ndf alone does not compare across ndf 23 and 98. Asymptotic p-values, no toys.
// Output: rootFiles/bPurityResults/fitQuality_1DTagVs2DMC_PF_fit5.root (h_dev_, h_ndf_ <1Dtag|2DMC>_<sys>)
// Figure: figures/bPurity/PFJets_PbPb/csvV2/fitQuality_1DTagVs2DMC_fit5.pdf
// Usage, from src/plots/bPurity/:  root -l -b -q plotFitQuality_1DTagVs2DMC_PF.C
#include "fitPtRelCSVv2Cut_PF_PbPb.C"
#include "Math/ProbFuncMathCore.h"
#include "Math/QuantFuncMathCore.h"

void plotFitQuality_1DTagVs2DMC_PF(double wp = 0.54){
  using namespace pbc; using namespace pcut;
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
  setvbuf(stdout, nullptr, _IOLBF, 0);
  nFit = nPR;   // fit 0-5 GeV
  setJetCollection("PFJets"); initPlotStyle();
  TFile *fDpp = TFile::Open(Form("%s/%s", repo, dataName)), *fMpp = TFile::Open(Form("%s/%s", repo, mcName));
  TFile *fDpb = TFile::Open(Form("%s/%s", repo, dataPbPb));
  for(int s = 0; s < NS; s++) fS[s] = TFile::Open(Form("%s/%s", repo, Form(phFmt, sName[s])));
  if(!fDpp || !fMpp || !fDpb || !fS[0] || !fS[1] || !fS[2]){ printf("ERROR: inputs\n"); return; }

  const int NCOL = 5; const char *colName[NCOL] = {"pp", "50-80%", "30-50%", "10-30%", "0-10%"}, *colTag[NCOL] = {"pp", "C4", "C3", "C2", "C1"};
  const int colCls[NCOL] = {-1, 3, 2, 1, 0};   // pbc class index (0 = 0-10%)
  const int NV = 2; const char *vTag[NV] = {"1Dtag", "2DMC"};
  TH1D *hDev[NV][NCOL], *hNdf[NV][NCOL], *hRat[NV][NCOL], *hZ[NV][NCOL];
  for(int v = 0; v < NV; v++) for(int c = 0; c < NCOL; c++){
    auto bk = [&](const char *n){ TH1D *h = new TH1D(Form("%s_%s_%s", n, vTag[v], colTag[c]), "", NW, ptEdges + 3); h->SetDirectory(nullptr); return h; };
    hDev[v][c] = bk("h_dev"); hNdf[v][c] = bk("h_ndf"); hRat[v][c] = bk("h_devOverNdf"); hZ[v][c] = bk("h_signif");
  }
  // pass 1: 2D, 4 CSVv2 categories; pass 2: the CSVv2 > wp tagged jets
  for(int pass = 0; pass < 2; pass++){
    const int v = pass == 0 ? 1 : 0;
    if(pass == 0) setCats("4cat"); else setCats(Form("cut%.2f", wp));
    for(int c = 0; c < NCOL; c++){
      if(c == 0 && !buildWindows(fDpp, fMpp)) return;
      for(int w = 0; w < NW; w++){
        Win X; double nEffB = 0.;
        if(c == 0) X = W[w]; else if(!buildCls(fDpb, colCls[c], w, X, nEffB)) return;
        Result r;
        if(pass == 0) r = fitWin(models[1], X, X.d);
        else { Win P = passWin(X); r = fitWin(models[0], P, P.d); }
        hDev[v][c]->SetBinContent(w + 1, r.dev); hNdf[v][c]->SetBinContent(w + 1, r.ndf);
      }
    }
  }
  printf("\ndeviance / ndf and significance, 1D tagged (CSVv2 > %.2f) vs 2D MC tagger shapes, 0.2 GeV bins, fit 0-5 GeV\n", wp);
  for(int c = 0; c < NCOL; c++){
    printf("%-7s", colName[c]);
    for(int w = 1; w <= NW; w++){
      for(int v = 0; v < NV; v++){
        const double d = hDev[v][c]->GetBinContent(w), n = hNdf[v][c]->GetBinContent(w);
        const double p = ROOT::Math::chisquared_cdf_c(d, n), z = ROOT::Math::normal_quantile_c(TMath::Max(p, 1e-300), 1.);
        hRat[v][c]->SetBinContent(w, d/n); hRat[v][c]->SetBinError(w, 0.);
        hZ[v][c]->SetBinContent(w, z); hZ[v][c]->SetBinError(w, 0.);
        printf(" %s %5.1f/%2.0f=%4.2f (%5.1f sigma)", v == 0 ? "|" : "", d, n, d/n, z);
      }
    }
    printf("\n");
  }
  {
    TFile *fo = TFile::Open(Form("%s/rootFiles/bPurityResults/fitQuality_1DTagVs2DMC_PF_fit5%s.root", repo, ptBinTag()), "RECREATE");
    for(int v = 0; v < NV; v++) for(int c = 0; c < NCOL; c++){ hDev[v][c]->Write(); hNdf[v][c]->Write(); }
    TNamed("info", Form("deviance and ndf of the CSVv2 > %.2f tagged-jet ptRel fit (1Dtag) and the 2D ptRel x CSVv2 4cat fit with MC tagger shapes (2DMC), 0.2 GeV bins, fit 0-5 GeV; pp, C4..C1 = 50-80 ... 0-10%%", wp)).Write();
    fo->Close();
  }

  // ---- figure
  const double ptLo = 80., ptHi = 300.;
  const char *vHex[NV] = {okabeHex[5], okabeHex[6]};
  const int vMark[NV] = {markFilledCircle, markFilledSquare};
  const double vSize[NV] = {2.6, 2.4};
  const char *vLeg[NV] = {Form("1D, CSVv2 > %.2f", wp), "2D, MC tagger shapes"};
  double yHi = 0., zLo = 0., zHi = 0.;
  for(int v = 0; v < NV; v++) for(int c = 0; c < NCOL; c++){ yHi = TMath::Max(yHi, hRat[v][c]->GetMaximum()); zHi = TMath::Max(zHi, hZ[v][c]->GetMaximum()); zLo = TMath::Min(zLo, hZ[v][c]->GetMinimum()); }
  yHi = 2.*ceil(0.5*1.2*yHi); zHi = 5.*ceil(1.1*zHi/5.); zLo = -5.*ceil(-1.1*zLo/5. + 0.01);   // round edges, so the axes carry round labels
  TCanvas *cv = new TCanvas("cvFitQ", "", 2300, 1100);
  GridPads g = makeGridPads(cv, NCOL, 2, 0.06, 0.01, 0.11, 0.075);
  const double fs = 0.055;
  for(int c = 0; c < NCOL; c++){
    for(int row = 0; row < 2; row++){
      g.at(c, row)->cd();
      TH1F *fr = gPad->DrawFrame(ptLo, row == 0 ? 0. : zLo, ptHi, row == 0 ? yHi : zHi); g.styleFrame(fr, c, row, fs);
      fr->GetYaxis()->SetLabelFont(43); fr->GetYaxis()->SetLabelSize(c == 0 ? 34 : 0);
      fr->GetXaxis()->SetLabelFont(43); fr->GetXaxis()->SetLabelSize(row == 1 ? 34 : 0);
      if(row == 0) fr->GetYaxis()->SetNdivisions(5, false); else fr->GetYaxis()->SetNdivisions(506);   // top: exactly 0, 2, ... (yHi even)
      if(row == 0) fr->GetYaxis()->ChangeLabel(1, -1, 0.);   // "0" against the lower row's top label
      const double ref = row == 0 ? 1. : 0.;
      TLine *l = new TLine(ptLo, ref, ptHi, ref); l->SetLineStyle(2); l->SetLineColor(kGray + 2); l->SetLineWidth(2); l->Draw();
      for(int v = NV - 1; v >= 0; v--){
        TH1D *h = row == 0 ? hRat[v][c] : hZ[v][c];
        styleH(h, vHex[v], vMark[v], vSize[v]); h->Draw("P X0 SAME");
      }
      gPad->RedrawAxis();
    }
    g.at(c, 0)->cd();
    g.text(c, 0, 0.10, colName[c], 0.065, 62);
    if(c == 1){   // the pp panel's 2D points sit at mid height: the legend goes where the frame is empty
      const double fL = gPad->GetLeftMargin(), fW = 1. - fL - gPad->GetRightMargin(), fT = 1. - gPad->GetTopMargin();
      TLegend *lg = new TLegend(fL + 0.04*fW, fT - 0.45, fL + fW, fT - 0.18);
      lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextFont(43); lg->SetTextSize(29);
      for(int v = 0; v < NV; v++) lg->AddEntry(hRat[v][c], vLeg[v], "p");
      lg->Draw();
    }
  }
  cv->cd();
  { TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.040); t.SetTextAlign(21);
    t.DrawLatex(g.lm + (1. - g.lm - g.rm)/2., 0.015, "#it{p}_{T}^{jet} (GeV)");
    t.SetTextAngle(90);
    const double mid = g.bm + 0.5*(1. - g.bm - g.tm);
    t.DrawLatex(0.018, mid + 0.5*(1. - g.tm - mid), "#chi^{2} / ndf");
    t.DrawLatex(0.018, g.bm + 0.5*(mid - g.bm), "significance (#sigma)"); }
  g.cmsLabel(cv, "Internal", "pp, PbPb (5.02 TeV)   #it{p}_{T}^{rel} fit 0-5 GeV", 1.0);
  const TString out = Form("%s/figures/bPurity/PFJets_PbPb/csvV2/fitQuality_1DTagVs2DMC_fit5%s.pdf", repo, ptBinTag());
  savePdfTight(cv, out);
}
