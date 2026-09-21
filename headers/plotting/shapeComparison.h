#pragma once
// Data-vs-MC shape comparison: two unit-area histograms overlaid, with a
// data/MC ratio panel underneath. Shared by the pp calo-jet and tag-muon
// kinematics macros, which draw the same kind of figure for eight observables.
//
// Returns chi2/ndf of the ratio against 1 over the drawn range, with both
// histograms' errors propagated (they are independent samples). Bins where the
// MC is empty come back from makeRatio() as kAbsent and are skipped.
//
// The x range stops just inside xMax: a range ending exactly on the last bin
// edge makes ROOT pull in the overflow bin and draw it.

#include <cmath>
#include <vector>
#include "TMath.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TString.h"
#include "plotStyle.h"
#include "ratioPanel.h"

// ratioTitle, rMin/rMax and errMode default to the data-vs-MC use. Pass
// RatioErr::kNumerator when the two histograms come from the SAME events (e.g.
// one scan re-run with a different setting): they are then strongly
// correlated, and propagating both errors overstates the ratio's uncertainty.
inline double drawShapeComparison(TH1D *hD, TH1D *hM,
                                  const char *xTitle, const char *yTitle,
                                  bool logY, double xMin, double xMax,
                                  const std::vector<TString> &headLines,
                                  const char *outPath,
                                  const char *labelD = "pp data",
                                  const char *labelM = "PYTHIA",
                                  const char *ratioTitle = "data / MC",
                                  double rMin = 0.5, double rMax = 1.5,
                                  RatioErr::Mode errMode = RatioErr::kBoth)
{
  TCanvas *c = new TCanvas(Form("c_%s", hD->GetName()), "", 700, 800);
  TPad *pT = nullptr, *pB = nullptr;
  splitPads(pT, pB);

  pT->cd();
  if(logY) pT->SetLogy();
  styleH(hD, hexData, markFilledCircle);
  styleH(hM, hexMC, markOpenSquare);

  double ymax = TMath::Max(hD->GetMaximum(), hM->GetMaximum());
  double ymin = 1e30;
  for(int b = 1; b <= hD->GetNbinsX(); b++){
    if(hD->GetBinCenter(b) < xMin || hD->GetBinCenter(b) > xMax) continue;
    if(hD->GetBinContent(b) > 0.) ymin = TMath::Min(ymin, hD->GetBinContent(b));
    if(hM->GetBinContent(b) > 0.) ymin = TMath::Min(ymin, hM->GetBinContent(b));
  }
  hD->SetTitle("");
  hD->GetXaxis()->SetRangeUser(xMin, xMax - 1e-6);
  hD->GetXaxis()->SetLabelSize(0);
  hD->GetYaxis()->SetTitle(yTitle);
  hD->GetYaxis()->SetTitleSize(0.050); hD->GetYaxis()->SetTitleOffset(1.55);
  hD->GetYaxis()->SetLabelSize(0.042);
  if(logY){ hD->SetMinimum(ymin*0.3); hD->SetMaximum(ymax*30.); }
  else    { hD->SetMinimum(0.);       hD->SetMaximum(ymax*1.45); }
  hD->Draw("E");
  hM->Draw("E same");
  hD->Draw("E same");

  TLegend *leg = makeLegend(0.60, 0.70, 0.93, 0.84, 0.040);
  leg->AddEntry(hD, labelD, "lp");
  leg->AddEntry(hM, labelM, "lp");
  leg->Draw();

  TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.038);
  for(size_t i = 0; i < headLines.size(); i++)
    la.DrawLatex(0.21, 0.885 - 0.045*i, headLines[i].Data());

  pB->cd();
  TH1D *r = makeRatio(hD, hM, Form("r_%s", hD->GetName()), errMode);
  styleH(r, hexData, markFilledCircle);
  styleRatioAxes(r, xTitle, ratioTitle);
  r->GetXaxis()->SetRangeUser(xMin, xMax - 1e-6);
  r->SetMinimum(rMin); r->SetMaximum(rMax);
  r->Draw("E");
  unityLine(xMin, xMax)->Draw();
  r->Draw("E same");

  double chi2 = 0.; int ndf = 0;
  for(int b = 1; b <= r->GetNbinsX(); b++){
    if(r->GetBinCenter(b) < xMin || r->GetBinCenter(b) > xMax) continue;
    double v = r->GetBinContent(b), e = r->GetBinError(b);
    if(v == kAbsent || e <= 0.) continue;
    chi2 += (v - 1.)*(v - 1.)/(e*e); ndf++;
  }

  c->SaveAs(outPath);
  delete c;
  return ndf > 0 ? chi2/ndf : -1.;
}
