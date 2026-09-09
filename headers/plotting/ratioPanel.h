#pragma once
// Ratio-panel construction. This exists because the same two bugs were written
// twice, independently, in two macros:
//
//   1. An empty denominator bin was set to 0. ROOT then draws a marker sitting
//      on the axis, which a reader takes as a measured ratio of zero. The
//      ratio there is UNDEFINED, and an undefined point must be absent, not
//      zero. Fixed by pushing it far below the frame.
//
//   2. The y-range was chosen by scanning ratio VALUES with an arbitrary
//      ceiling (6) to reject blow-ups. That works while ratios sit near unity
//      and fails completely otherwise: in the muon-injection study the honest
//      ratio is ~9-25, so the ceiling pushed every real point off-scale and
//      left a panel showing one bin. Blow-ups come from a poorly populated
//      DENOMINATOR, not from a large honest ratio, so gate on the denominator
//      and take the maximum with no ceiling.
//
// ERRORS. Two conventions, and the right one depends on what the denominator
// is:
//   RatioErr::kBoth       independent measurements -- propagate both.
//   RatioErr::kNumerator  the denominator is a fixed reference (an MC curve,
//                         or a baseline built from the SAME events as the
//                         numerator). Propagating it would overstate the
//                         uncertainty, badly when the two are strongly
//                         correlated.
// Most ratios in this analysis are the second kind.

#include "TH1D.h"
#include "TLine.h"
#include "TPad.h"

namespace RatioErr { enum Mode { kBoth, kNumerator }; }

// Bins whose denominator is empty come back at kAbsent, far below any sensible
// frame, so they are simply not drawn.
const double kAbsent = -999.;

inline TH1D* makeRatio(TH1D *num, TH1D *den, const char *name,
                       RatioErr::Mode mode = RatioErr::kNumerator)
{
  TH1D *r = (TH1D*) num->Clone(name);
  r->SetDirectory(nullptr);
  for(int b = 1; b <= r->GetNbinsX(); b++){
    double n = num->GetBinContent(b), d = den->GetBinContent(b);
    if(d == 0.){ r->SetBinContent(b, kAbsent); r->SetBinError(b, 0.); continue; }
    double v = n/d;
    double e = (num->GetBinError(b))/d;
    if(mode == RatioErr::kBoth && n != 0.){
      double rn = num->GetBinError(b)/n, rd = den->GetBinError(b)/d;
      e = fabs(v)*sqrt(rn*rn + rd*rd);
    }
    r->SetBinContent(b, v);
    r->SetBinError(b, e);
  }
  return r;
}

// Upper limit for the ratio frame, from bins where the denominator is actually
// populated (above minDenFrac of its own peak). Out in a tail where the
// denominator vanishes the ratio diverges for reasons that say more about the
// denominator's statistics than about the comparison; letting those bins set
// the scale compresses the region the panel exists for. They are still drawn,
// clipped at the frame edge, rather than dropped.
inline double ratioMax(TH1D *r, TH1D *den, double xMin, double xMax,
                       double floorVal = 1.2, double minDenFrac = 0.02)
{
  double m = floorVal, peak = den->GetMaximum();
  if(peak <= 0.) return floorVal;
  for(int b = 1; b <= r->GetNbinsX(); b++){
    double x = r->GetBinCenter(b);
    if(x < xMin || x > xMax) continue;
    if(den->GetBinContent(b) < minDenFrac*peak) continue;
    double v = r->GetBinContent(b) + r->GetBinError(b);
    if(v > m) m = v;
  }
  return m;
}

// Standard 2-pad split. Top pad gets no x labels; the ratio pad carries them.
inline void splitPads(TPad *&top, TPad *&bot, double split = 0.34,
                      double leftMargin = 0.17)
{
  top = new TPad("padTop", "", 0, split, 1, 1);
  bot = new TPad("padBot", "", 0, 0,     1, split);
  top->SetBottomMargin(0.02); top->SetLeftMargin(leftMargin); top->SetTopMargin(0.07);
  bot->SetTopMargin(0.02);    bot->SetLeftMargin(leftMargin); bot->SetBottomMargin(0.32);
  top->Draw(); bot->Draw();
}

// Axis sizes for the ratio pad, which is ~1/3 the height of the top pad and so
// needs correspondingly larger relative text.
inline void styleRatioAxes(TH1D *r, const char *xTitle, const char *yTitle)
{
  r->SetTitle("");
  r->GetXaxis()->SetTitle(xTitle);
  r->GetXaxis()->SetTitleSize(0.105); r->GetXaxis()->SetTitleOffset(1.25);
  r->GetXaxis()->SetLabelSize(0.090);
  r->GetYaxis()->SetTitle(yTitle);
  r->GetYaxis()->SetTitleSize(0.085); r->GetYaxis()->SetTitleOffset(0.80);
  r->GetYaxis()->SetLabelSize(0.085); r->GetYaxis()->SetNdivisions(505);
}

inline TLine* unityLine(double xMin, double xMax)
{
  TLine *l = new TLine(xMin, 1.0, xMax, 1.0);
  l->SetLineStyle(2); l->SetLineColor(kGray+2);
  return l;
}
