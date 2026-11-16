// b purity from a fit of the FULL muon ptRel x CSVv2 SURFACE, without the coarse CSVv2 working-point categories
// (2026-11-16, user: "try fitting without using coarse csvv2 classification ... some type of goodness of fit plot,
// specifically the 2d surface fit / data, maybe drawn in colz format"). Copy of the fit loop of
// fitPtRelCSVv2_2D_ppPbPb_PF.C (included: same inputs, templates, g->bb share, PbPb fake-pair subtraction, model
// "MC" = every flavor's CSVv2 shape from simulation, c / light at the MC ratio, f_b the only free parameter), with
// the CSVv2 axis in nCsv UNIFORM bins over [0, 1] instead of the 4 categories: 50 = the scans' own 0.02 grid,
// 25 and 10 for the trend. The first bin also holds the jets without a discriminant value. ptRel 0-5 GeV in the
// nominal 0.2 GeV bins, so 50 x 25 = 1250 cells per window.
// Fit: the hand-written binned Poisson likelihood with Barlow-Beeston lite profiled analytically per cell
// (p2d::fitWinMinuit), NOT HistFactory: one constrained gamma per cell would be 1250 nuisance parameters per fit.
// At 4 categories the two agree within 0.005 except where the model fits worst (see the 4cat macro).
// Reference: the nominal 4-category HistFactory purity, rootFiles/bPurityResults/ptRelCSVv2_2D_PF_ppPbPb_4cat_fit5.root,
// and its hand-written-likelihood twin (h_bPurity2D_MC_<sys>), the like-for-like comparison.
//
// Goodness of fit, per system (one PDF each, columns = the five jet-pT windows 80-300 GeV):
//   top row     fit / data per cell (x = CSVv2, y = ptRel), COLZ, blank where the data cell is empty or <= 0
//   bottom row  pull (fit - data) / sqrt(fit + template variance) (same color sense as fit / data), blank where neither data nor MC
//   The fit is drawn with the Barlow-Beeston factors at 1 (the template prediction itself).
//   Diverging blue - white - vermillion scale (Okabe-Ito blue and vermillion; white = agreement), values beyond the
//   range drawn at the end colors.
// Printed per window: purity, deviance / ndf (with the template-statistics penalty; NOT chi2-distributed where the
// cells are sparse), and sum pull^2 / cells over the cells with fit >= 5 jets, which is.
// Summary figure: purity vs jet pT, 1x5 row, 4 categories against the uniform binnings, sum pull^2 / cells below.
//
// Outputs: rootFiles/bPurityResults/ptRelCSVv2_fineSurface_PF_ppPbPb_fit5.root (h_bPurity_csv<n>_<sys>,
//          h_pull2_csv<n>_<sys>, the maps h_fitOverData_ / h_pull_ / h_data_ / h_fit_csv<n>_<sys>_w<w>)
//          figures/bPurity/PFJets_PbPb/csvV2/fineSurface/fitOverData_csv<n>_<sys>.pdf, purity_fineSurface.pdf
// C / LIGHT FREE (2026-11-16, user: "try the fine-surface fit with the c / light ratio free"): argument true. The c share
// of the non-b jets, r = c / (c + light), floats in [0, 1] next to f_b; the CSVv2 shapes of every flavor stay MC.
// Outputs tagged _cFree (root file, maps, purity_fineSurface_cFree.pdf). Its summary: purity against the c / light-fixed
// fit of the same binning (read from the fixed-mode root file, run that first), and the fitted r against the MC r below.
// Usage, from src/plots/bPurity/:  root -l -b -q 'fitPtRelCSVv2_fineSurface_ppPbPb_PF.C'   then   '...PF.C(true)'
#define P2D_KMAX 50
#include "fitPtRelCSVv2_2D_ppPbPb_PF.C"
#include "../../../headers/plotting/gridPads.h"
#include "TPaletteAxis.h"
#include "TColor.h"
#include "TH2D.h"

namespace fsurf {
using namespace p2d;
// uniform CSVv2 bins; n must divide the scans' 50 bins of 0.02
inline bool setUniform(int n){
  if(n < 2 || n > KMAX || 50 % n != 0) return false;
  K = n; for(int k = 0; k < K; k++) catUp[k] = (k + 1.)/n; catUp[K - 1] = 1.01;
  return true;
}
inline void divergingPalette(){
  double stops[3] = {0., 0.5, 1.}, r[3] = {0./255, 1., 213./255}, g[3] = {114./255, 1., 94./255}, b[3] = {178./255, 1., 0.};
  TColor::CreateGradientColorTable(3, stops, r, g, b, 99); gStyle->SetNumberContours(99);
}
}

void fitPtRelCSVv2_fineSurface_ppPbPb_PF(bool cFree = false){
  using namespace pbc; using namespace fsurf;
  setvbuf(stdout, nullptr, _IOLBF, 0);
  nFit = nPR;   // fit 0-5 GeV
  setJetCollection("PFJets"); initPlotStyle();
  TFile *fDpp = TFile::Open(Form("%s/%s", repo, dataName)), *fMpp = TFile::Open(Form("%s/%s", repo, mcName));
  TFile *fDpb = TFile::Open(Form("%s/%s", repo, dataPbPb));
  for(int s = 0; s < NS; s++) fS[s] = TFile::Open(Form("%s/%s", repo, Form(phFmt, sName[s])));
  TFile *fRef = TFile::Open(Form("%s/rootFiles/bPurityResults/ptRelCSVv2_2D_PF_ppPbPb_4cat_fit5.root", repo));
  if(!fDpp || !fMpp || !fDpb || !fS[0] || !fS[1] || !fS[2] || !fRef){ printf("ERROR: inputs\n"); return; }
  static const Model mCFree = {"cFree", "2D, MC tagger shapes, c / light free", false, false, false, true};
  const Model &M = cFree ? mCFree : models[1];
  const TString mTag = cFree ? "_cFree" : "";
  TFile *fFix = cFree ? TFile::Open(Form("%s/rootFiles/bPurityResults/ptRelCSVv2_fineSurface_PF_ppPbPb_fit5.root", repo)) : nullptr;
  if(cFree && !fFix){ printf("ERROR: run the c / light-fixed mode first\n"); return; }
  const TString outDir = Form("%s/figures/bPurity/PFJets_PbPb/csvV2/fineSurface/", repo); gSystem->mkdir(outDir, true);

  const int NCOL = 5; const char *colName[NCOL] = {"pp", "PbPb 50-80%", "PbPb 30-50%", "PbPb 10-30%", "PbPb 0-10%"}, *colShort[NCOL] = {"pp", "50-80%", "30-50%", "10-30%", "0-10%"}, *colTag[NCOL] = {"pp", "C4", "C3", "C2", "C1"};
  const int colCls[NCOL] = {-1, 3, 2, 1, 0};
  const int NB = 3; const int nCsv[NB] = {10, 25, 50};
  const double rLo = 0., rHi = 2., pMax = 5., yTop = 6.2;   // fit / data and pull color ranges; frame top (labels above 5 GeV)
  TH1D *hPur[NB][NCOL], *hP2[NB][NCOL], *hRefHF[NCOL], *hRefMn[NCOL], *hR_c[NB][NCOL], *hRcMC[NCOL], *hFix[NB][NCOL];
  TFile *fo = TFile::Open(Form("%s/rootFiles/bPurityResults/ptRelCSVv2_fineSurface_PF_ppPbPb_fit5%s.root", repo, mTag.Data()), "RECREATE");
  for(int c = 0; c < NCOL; c++){
    hRefHF[c] = (TH1D*) fRef->Get(Form("h_bPurity2D_MC_HistFactory_%s", colTag[c])); hRefMn[c] = (TH1D*) fRef->Get(Form("h_bPurity2D_MC_%s", colTag[c]));
    if(!hRefHF[c] || !hRefMn[c]){ printf("ERROR: 4cat reference %s\n", colTag[c]); return; }
  }

  for(int ib = 0; ib < NB; ib++){
    if(!setUniform(nCsv[ib])){ printf("ERROR: nCsv %d\n", nCsv[ib]); return; }
    for(int c = 0; c < NCOL; c++){
      hPur[ib][c] = new TH1D(Form("h_bPurity_csv%d_%s", K, colTag[c]), "", NW, ptEdges + 3); hPur[ib][c]->SetDirectory(nullptr);
      hP2[ib][c] = (TH1D*) hPur[ib][c]->Clone(Form("h_pull2_csv%d_%s", K, colTag[c])); hP2[ib][c]->SetDirectory(nullptr);
      hR_c[ib][c] = (TH1D*) hPur[ib][c]->Clone(Form("h_cShare_csv%d_%s", K, colTag[c])); hR_c[ib][c]->SetDirectory(nullptr);
      if(ib == 0){ hRcMC[c] = (TH1D*) hPur[ib][c]->Clone(Form("h_cShareMC_%s", colTag[c])); hRcMC[c]->SetDirectory(nullptr); }
      hFix[ib][c] = cFree ? (TH1D*) fFix->Get(Form("h_bPurity_csv%d_%s", K, colTag[c])) : nullptr;
      if(cFree && !hFix[ib][c]){ printf("ERROR: fixed-mode purity csv%d %s\n", K, colTag[c]); return; }
      if(c == 0 && !buildWindows(fDpp, fMpp)) return;
      printf("\n==== %s, ptRel x CSVv2 surface, %d CSVv2 bins x %d ptRel bins, model %s (hand-written likelihood)\n", colName[c], K, nPR, M.label);
      printf("  window   | purity (stat)   st | 4cat minuit | 4cat HistFactory | deviance / ndf  | sum pull^2 / cells (fit >= 5) | data in cells without MC | cells: data > 0, MC > 0 | c / (c + light): fit, MC | purity c / light fixed, same bins\n");
      TCanvas *cv = new TCanvas(Form("cvS_%d_%d", ib, c), "", 2500, 1100);
      GridPads g = makeGridPads(cv, NW, 2, 0.055, 0.075, 0.11, 0.08);
      divergingPalette();
      TH2D *hR[NW], *hPl[NW];
      for(int w = 0; w < NW; w++){
        static Win X; X = Win(); double fakeFrac = 0.;
        if(c == 0) X = W[w];
        else { double nEffB = 0.; if(!buildCls(fDpb, colCls[c], w, X, nEffB)) return; fakeFrac = lastFakeFrac; }
        Result r = fitWinMinuit(M, X, X.d);
        const double pur = r.pur*(1. - fakeFrac), epur = r.epur*(1. - fakeFrac);
        hPur[ib][c]->SetBinContent(w + 1, pur); hPur[ib][c]->SetBinError(w + 1, epur);
        // the fitted template prediction (Barlow-Beeston factors at 1) and its template variance
        static double mu[KMAX][NPR], va[KMAX][NPR];
        cur = &M; cw = &X; cd = X.d; Kf = K; expectAll(r.par, mu, va);
        const TString id = Form("csv%d_%s_w%d", K, colTag[c], w);
        hR[w] = new TH2D("h_fitOverData_" + id, "", K, 0., 1., nPR, 0., nPR*ptRelBW); hPl[w] = (TH2D*) hR[w]->Clone("h_pull_" + id);
        TH2D *hD = (TH2D*) hR[w]->Clone("h_data_" + id), *hF = (TH2D*) hR[w]->Clone("h_fit_" + id);
        for(TH2D *h : {hR[w], hPl[w], hD, hF}) h->SetDirectory(nullptr);
        double p2 = 0., noMC = 0., sD = 0.; int nP2 = 0, nDat = 0, nMC = 0;
        for(int k = 0; k < K; k++) for(int j = 0; j < nPR; j++){
          const double d = X.d[k][j], m = mu[k][j];
          hD->SetBinContent(k + 1, j + 1, d); hF->SetBinContent(k + 1, j + 1, m);
          sD += d; if(d > 0.) nDat++; if(m > 0.) nMC++;
          if(m <= 0.){ noMC += d; continue; }
          const double pl = (m - d)/sqrt(m + va[k][j]);   // fit - data, so the colors mean the same as in fit / data
          if(m >= 5.){ p2 += pl*pl; nP2++; }
          // clamped just inside the color range so off-scale cells show the end colors instead of white
          hPl[w]->SetBinContent(k + 1, j + 1, TMath::Max(-pMax + 1e-3, TMath::Min(pMax - 1e-3, fabs(pl) < 1e-6 ? 1e-6 : pl)));
          if(d > 0.) hR[w]->SetBinContent(k + 1, j + 1, TMath::Max(rLo + 1e-3, TMath::Min(rHi - 1e-3, m/d)));
        }
        hP2[ib][c]->SetBinContent(w + 1, nP2 > 0 ? p2/nP2 : 0.);
        hR_c[ib][c]->SetBinContent(w + 1, cFree ? r.par[1] : X.rc); hR_c[ib][c]->SetBinError(w + 1, cFree ? r.err[1] : 0.); hRcMC[c]->SetBinContent(w + 1, X.rc);
        const int rb = hRefHF[c]->FindBin(0.5*(wLo[w] + wHi[w]));
        printf("  %3.0f-%3.0f  | %.3f +- %.3f (%d%d) | %.3f       | %.3f            | %7.1f / %-4d  | %7.1f / %-4d = %.2f        | %.2f%%                    | %d, %d of %d | %.3f +- %.3f, %.3f | %s\n", wLo[w], wHi[w], pur, epur, r.status, r.hesse,
               hRefMn[c]->GetBinContent(rb), hRefHF[c]->GetBinContent(rb), r.dev, r.ndf, p2, nP2, nP2 > 0 ? p2/nP2 : 0., sD > 0. ? 100.*noMC/sD : 0., nDat, nMC, K*nPR, hR_c[ib][c]->GetBinContent(w + 1), hR_c[ib][c]->GetBinError(w + 1), X.rc, cFree ? Form("%.3f", hFix[ib][c]->GetBinContent(w + 1)) : "-");
        fo->cd(); hR[w]->Write(); hPl[w]->Write(); hD->Write(); hF->Write();

        for(int row = 0; row < 2; row++){
          g.at(w, row)->cd();
          TH1F *fr = gPad->DrawFrame(0., 0., 1., yTop); g.styleFrame(fr, w, row, 0.055);
          fr->GetXaxis()->SetLabelFont(43); fr->GetXaxis()->SetLabelSize(row == 1 ? 34 : 0); fr->GetYaxis()->SetLabelFont(43); fr->GetYaxis()->SetLabelSize(w == 0 ? 34 : 0);
          fr->GetYaxis()->SetNdivisions(506);
          if(w > 0) fr->GetXaxis()->ChangeLabel(1, -1, 0.);   // the seam's "0" is cut by the pad edge
          TH2D *h = row == 0 ? hR[w] : hPl[w];
          h->SetMinimum(row == 0 ? rLo : -pMax); h->SetMaximum(row == 0 ? rHi : pMax);
          h->Draw("COL SAME");
          TLine *l5 = new TLine(0., 5., 1., 5.); l5->SetLineColor(kGray + 2); l5->SetLineWidth(1); l5->Draw();
          gPad->RedrawAxis();
          TLatex t; t.SetTextFont(43); t.SetTextSize(30); t.SetTextAlign(12);
          if(row == 0){ t.DrawLatex(0.04, 5.83, Form("%.0f-%.0f GeV", wLo[w], wHi[w])); t.DrawLatex(0.04, 5.33, Form("#it{P} = %.3f #pm %.3f", pur, epur)); }
          else t.DrawLatex(0.04, 5.6, Form("#Sigma pull^{2} / cells = %.2f", nP2 > 0 ? p2/nP2 : 0.));
        }
      }
      // one color scale per row, right of the last column
      cv->cd();
      const double ih = (1. - g.bm - g.tm)/2., fx = 5./yTop;   // the maps fill the frames up to 5 GeV
      for(int row = 0; row < 2; row++){
        const double y1 = 1. - g.tm - (row + 1)*ih, y2 = y1 + ih*fx;
        TPaletteAxis *pa = new TPaletteAxis(1. - g.rm + 0.006, y1, 1. - g.rm + 0.020, y2, row == 0 ? hR[NW-1] : hPl[NW-1]);
        pa->SetLabelFont(43); pa->SetLabelSize(30); pa->SetNdivisions(505); pa->Draw();
        TLatex t; t.SetNDC(); t.SetTextFont(43); t.SetTextSize(34); t.SetTextAlign(21); t.SetTextAngle(90);
        t.DrawLatex(0.992, 0.5*(y1 + y2), row == 0 ? "fit / data" : "(fit #minus data) / #sigma");
      }
      g.titles(cv, "CSVv2 discriminant", "muon #it{p}_{T}^{rel} (GeV)", 0.042, 0.016, 0.018);
      g.cmsLabel(cv, "Internal", Form("%s (5.02 TeV)   muon-tagged jets, %d #times %d cells%s", colName[c], K, nPR, cFree ? ", c / light free" : ""), 1.0);
      savePdfTight(cv, outDir + Form("fitOverData%s_csv%d_%s.pdf", mTag.Data(), K, colTag[c]));
      delete cv;
      fo->cd(); hPur[ib][c]->Write(); hP2[ib][c]->Write(); hR_c[ib][c]->Write(); if(ib == 0) hRcMC[c]->Write();
    }
  }

  // ---- summary: purity vs jet pT (top) and sum pull^2 / cells (bottom), 4 categories against the uniform binnings
  {
    const char *bHex[NB] = {okabeHex[1], okabeHex[2], okabeHex[5]}; const int bMark[NB] = {markFilledDiamond, markFilledSquare, markFilledCircle}; const double bSize[NB] = {3.4, 2.4, 1.9};
    TCanvas *cv = new TCanvas("cvSum", "", 2300, 1000);
    GridPads g = makeGridPads(cv, NCOL, 2, 0.055, 0.01, 0.12, 0.075);
    const double split = g.bm + 0.36*(1. - g.bm - g.tm);
    for(int c = 0; c < NCOL; c++){
      TPad *t = g.at(c, 0), *b = g.at(c, 1); double x1, y1, x2, y2;
      t->GetPadPar(x1, y1, x2, y2); t->SetPad(x1, split, x2, y2); t->SetTopMargin(g.tm/(1. - split));
      b->GetPadPar(x1, y1, x2, y2); b->SetPad(x1, 0., x2, split); b->SetBottomMargin(g.bm/split);
    }
    double p2Max = 1.5; for(int ib = 0; ib < NB; ib++) for(int c = 0; c < NCOL; c++) p2Max = TMath::Max(p2Max, hP2[ib][c]->GetMaximum());
    for(int c = 0; c < NCOL; c++){
      g.at(c, 0)->cd();
      TH1F *fr = gPad->DrawFrame(80., 0.2, 300., 0.8); g.styleFrame(fr, c, 0, 0.055);
      fr->GetYaxis()->SetLabelFont(43); fr->GetYaxis()->SetLabelSize(c == 0 ? 34 : 0); fr->GetYaxis()->ChangeLabel(1, -1, 0.); fr->GetYaxis()->SetNdivisions(507);
      // open circles: fixed mode = the 4-category fit; c / light free = the c / light-fixed fit at 50 bins
      TH1D *h4 = (TH1D*) (cFree ? hFix[NB-1][c] : hRefMn[c])->Clone(Form("ref4_%d", c)); h4->SetDirectory(nullptr);
      styleH(h4, hexData, markOpenCircle, 3.6); h4->SetLineWidth(3); h4->Draw("E1 X0 SAME");
      for(int ib = 0; ib < NB; ib++){ styleH(hPur[ib][c], bHex[ib], bMark[ib], bSize[ib]); hPur[ib][c]->SetLineWidth(3); hPur[ib][c]->Draw("E1 X0 SAME"); }
      gPad->RedrawAxis();
      g.text(c, 0, 0.10, colShort[c], 0.065, 62);
      if(c == 0){
        const double fL = gPad->GetLeftMargin(), fW = 1. - fL - gPad->GetRightMargin(), fT = 1. - gPad->GetTopMargin();
        TLegend *lg = new TLegend(fL + 0.04*fW, fT - 0.46, fL + fW, fT - 0.14);
        lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextFont(43); lg->SetTextSize(30);
        lg->AddEntry(h4, cFree ? "50 bins, c / light fixed" : "4 CSVv2 categories", "pe");
        for(int ib = 0; ib < NB; ib++) lg->AddEntry(hPur[ib][c], Form(cFree ? "%d bins, c / light free" : "%d uniform bins", nCsv[ib]), "pe");
        lg->Draw();
      }
      g.at(c, 1)->cd();
      TH1F *fb = gPad->DrawFrame(80., 0., 300., cFree ? 1.1 : 1.15*p2Max); g.styleFrame(fb, c, 1, 0.055);
      fb->GetXaxis()->SetLabelFont(43); fb->GetXaxis()->SetLabelSize(34); fb->GetYaxis()->SetLabelFont(43); fb->GetYaxis()->SetLabelSize(c == 0 ? 34 : 0);
      fb->GetYaxis()->SetNdivisions(505); fb->GetXaxis()->SetTickLength(0.03*0.64/0.36);
      if(!cFree){
        TLine *one = new TLine(80., 1., 300., 1.); one->SetLineStyle(2); one->SetLineColor(kGray + 2); one->SetLineWidth(2); one->Draw();
        for(int ib = 0; ib < NB; ib++){ styleH(hP2[ib][c], bHex[ib], bMark[ib], bSize[ib]); hP2[ib][c]->Draw("P X0 SAME"); }
      } else {   // fitted c share of the non-b jets against the MC value (dashed)
        styleLine(hRcMC[c], "#8C8C8C", 4); hRcMC[c]->SetLineStyle(7); hRcMC[c]->Draw("HIST SAME");
        for(int ib = 0; ib < NB; ib++){ styleH(hR_c[ib][c], bHex[ib], bMark[ib], bSize[ib]); hR_c[ib][c]->SetLineWidth(3); hR_c[ib][c]->Draw("E1 X0 SAME"); }
      }
      gPad->RedrawAxis();
    }
    cv->cd();
    { TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.040); t.SetTextAlign(21);
      t.DrawLatex(g.lm + (1. - g.lm - g.rm)/2., 0.018, "#it{p}_{T}^{jet} (GeV)");
      t.SetTextAngle(90); t.DrawLatex(0.016, split + 0.5*(1. - g.tm - split), "b purity of muon-tagged jets");
      t.SetTextSize(0.030); t.DrawLatex(0.016, g.bm + 0.5*(split - g.bm), cFree ? "c / (c + light)" : "#Sigma pull^{2} / cells"); }
    g.cmsLabel(cv, "Internal", "pp, PbPb (5.02 TeV)   fit of #it{p}_{T}^{rel} #times CSVv2", 1.0);
    savePdfTight(cv, outDir + "purity_fineSurface" + mTag + ".pdf");
  }
  fo->Close();
}
