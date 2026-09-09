#pragma once
// Shared plot styling. The Okabe-Ito hex codes were pasted into 52 files
// before this header existed.
//
// PALETTE: Okabe-Ito, colour-blind safe. Use okabe[] in order for series that
// have no natural meaning; use the named constants where the role is fixed, so
// the same thing keeps the same colour across figures.
//
// MARKERS: centre-symmetric shapes only -- circle, square, diamond, cross,
// open variants. No triangles anywhere in this analysis.

#include "TColor.h"
#include "TH1.h"
#include "TGaxis.h"
#include "TStyle.h"
#include "TLegend.h"

// ---------------------------------------------------------------- palette ---
const char *okabeHex[8] = {
  "#000000",  // black
  "#E69F00",  // orange
  "#56B4E9",  // sky blue
  "#009E73",  // bluish green
  "#F0E442",  // yellow   -- poor on white, avoid for lines
  "#0072B2",  // blue
  "#D55E00",  // vermillion
  "#CC79A7"   // reddish purple
};

// Fixed roles, so a reader who has seen one figure can read the next.
const char *hexData     = "#000000";   // the measured thing
const char *hexMC       = "#0072B2";   // prediction
const char *hexCorrected= "#D55E00";   // measurement after a correction
const char *hexBkg1     = "#56B4E9";   // background template
const char *hexBkg2     = "#009E73";   // second background template
const char *hexSame     = "#0072B2";   // same-event
const char *hexMixed    = "#D55E00";   // mixed-event

// centre-symmetric only
const int markFilledCircle = 20, markFilledSquare = 21, markFilledDiamond = 33;
const int markOpenCircle   = 24, markOpenSquare   = 25, markOpenDiamond   = 27;
const int markCross        = 5;

inline void styleH(TH1 *h, const char *hex, int mark, double size = 1.0)
{
  int c = TColor::GetColor(hex);
  h->SetLineColor(c); h->SetMarkerColor(c);
  h->SetMarkerStyle(mark); h->SetMarkerSize(size); h->SetLineWidth(2);
  h->SetStats(0);
}
inline void styleLine(TH1 *h, const char *hex, int width = 3)
{
  int c = TColor::GetColor(hex);
  h->SetLineColor(c); h->SetMarkerColor(c);
  h->SetLineWidth(width); h->SetMarkerSize(0); h->SetStats(0);
}

// Call once at the top of a macro.
//
// SetMaxDigits matters more than it looks: per-event rates here run ~1e-4 to
// ~1e-6, and without it ROOT writes axis labels as 0.0002, 0.0004, ... which
// are wide enough to push the y-axis title clean off the canvas. With it you
// get a shared "x10^-4" header instead -- but ROOT draws that at the top left
// of the frame, so keep any TLatex header lines below about y = 0.87 NDC.
inline void initPlotStyle()
{
  gStyle->SetOptStat(0);
  TGaxis::SetMaxDigits(3);
}

inline TLegend* makeLegend(double x1, double y1, double x2, double y2, double size = 0.038)
{
  TLegend *l = new TLegend(x1, y1, x2, y2);
  l->SetBorderSize(0); l->SetFillStyle(0); l->SetTextSize(size);
  return l;
}

// ------------------------------------------------------------- rebinning ---
// TH1::Rebin with an edge array returns a NEW histogram and leaves the original
// alone -- forgetting to release it leaks one histogram per call, which over a
// loop of panels adds up. Every edge must land on an existing bin boundary of
// the input, since Rebin can only merge existing bins; check that before
// editing an edge list.
//
// Variable widths mean the result MUST be converted to a density before
// drawing, or a wide bin reads as a tall one. Do that AFTER any integrals you
// need, which are counts, not densities -- see divideByBinwidth.h.
inline TH1D* rebinTo(TH1D *h, int nBins, const double *edges, const char *name)
{
  TH1D *r = (TH1D*) h->Rebin(nBins, name, edges);
  r->SetDirectory(nullptr);
  return r;
}

// Binnings used across the muon/ptRel plots. Sharing them is what lets two
// figures be overlaid or divided; changing one silently breaks that.
const double edge_ptRel_0to5[] = {0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,
                                  2.0,2.4,2.8,3.4,4.2,5.0};
const int    nEdge_ptRel_0to5 = (int)(sizeof(edge_ptRel_0to5)/sizeof(double)) - 1;

const double edge_dR_0to0p5[] = {0,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50};
const int    nEdge_dR_0to0p5 = (int)(sizeof(edge_dR_0to0p5)/sizeof(double)) - 1;
