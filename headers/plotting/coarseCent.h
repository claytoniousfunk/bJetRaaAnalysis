#pragma once
// Coarse centrality classes and the fine-slice sums that build them.
//
// This mapping was duplicated in 31 macros before this header existed. Two
// forms were in use and both are kept, because both are correct for their
// callers:
//
//   NCoarse = 4   classes only:              0-10, 10-30, 30-50, 50-80
//   NCoarse5 = 5  inclusive bin prepended:   0-80 (from the input's own C0),
//                                            then the four above
//
// The 5-form's first entry maps to slice 0, i.e. the scan's own inclusive
// histogram, NOT a sum of slices 1-16. Those differ: C0 is filled for every
// event, while slices 1-16 cover only the hiBin range the centrality binning
// spans. Summing 1-16 to make an inclusive bin is a real (and silent) error.
//
// Slice indices refer to the ultraFine (5%) binning, NCentralityIndices = 17.
// A 4CentBins scan has 5 slices and this mapping does NOT apply to it -- check
// how many h_vz_C* keys a file carries before trusting a coarse merge.

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"

const int   NCoarse = 4;
const int   coarseSliceLo[NCoarse] = { 1,  3,  7, 11};
const int   coarseSliceHi[NCoarse] = { 2,  6, 10, 16};
const char *coarseLabel[NCoarse]   = { "0-10%", "10-30%", "30-50%", "50-80%" };
const char *coarseTag[NCoarse]     = { "0to10pct", "10to30pct", "30to50pct", "50to80pct" };
const int   coarseHiBinLo[NCoarse] = {  0,  20,  60, 100};
const int   coarseHiBinHi[NCoarse] = { 20,  60, 100, 160};

const int   NCoarse5 = 5;
const int   coarseSliceLo5[NCoarse5] = { 0,  1,  3,  7, 11};
const int   coarseSliceHi5[NCoarse5] = { 0,  2,  6, 10, 16};
const char *coarseLabel5[NCoarse5]   = { "0-80% (incl.)", "0-10%", "10-30%", "30-50%", "50-80%" };

// Events in a coarse class. vzBase selects the denominator:
//   "h_vz"            all events        -- for ungated numerators
//   "h_vz_triggerOn"  triggered events  -- for numerators gated on
//                                          evtTriggerDecision
//
// Getting this pairing wrong is the single easiest mistake in this codebase:
// h_vz is filled before evtTriggerDecision is even computed, so a gated
// numerator over h_vz understates the rate by the trigger fraction -- which is
// itself strongly centrality dependent (0.120 in 0-10% down to 0.025 in
// 50-80%), so it does not even cancel in a ratio across classes.
//
// Caveat on h_vz_triggerOn: it is filled on the raw bit HLT_HIL3Mu12_v1 == 1,
// while evtTriggerDecision additionally requires a non-zero prescale via
// triggerIsOn(). It is therefore a slight overestimate of the true denominator.
//
// Returns -1 if any slice is missing -- a partial sum would be silently wrong.
inline double coarseEvents(TFile *f, int lo, int hi, const char *vzBase = "h_vz")
{
  double n = 0.;
  for(int si = lo; si <= hi; si++){
    TH1D *h = nullptr;
    f->GetObject(Form("%s_C%d", vzBase, si), h);
    if(!h) return -1.;
    n += h->Integral();
  }
  return n;
}
inline double coarseEvents(TFile *f, int ci, const char *vzBase = "h_vz")
{ return coarseEvents(f, coarseSliceLo[ci], coarseSliceHi[ci], vzBase); }

// Sum a TH1/TH2 over a coarse class. Returns nullptr if ANY slice is missing,
// rather than a partial sum that would be mis-normalised against a complete
// event count. Caller owns the result (SetDirectory(nullptr) is applied).
//
// tag only has to make the returned histogram's name unique within the macro;
// ROOT will warn about duplicate names otherwise.
inline TH2D* coarseSum2(TFile *f, const char *base, int ci, const char *tag)
{
  TH2D *s = nullptr;
  for(int si = coarseSliceLo[ci]; si <= coarseSliceHi[ci]; si++){
    TH2D *h = nullptr;
    f->GetObject(Form("%s_C%d", base, si), h);
    if(!h){ if(s) delete s; return nullptr; }
    if(!s){ s = (TH2D*) h->Clone(Form("cs2_%s_%s_%d", base, tag, ci)); s->SetDirectory(nullptr); }
    else s->Add(h);
  }
  return s;
}
inline TH1D* coarseSum1(TFile *f, const char *base, int ci, const char *tag)
{
  TH1D *s = nullptr;
  for(int si = coarseSliceLo[ci]; si <= coarseSliceHi[ci]; si++){
    TH1D *h = nullptr;
    f->GetObject(Form("%s_C%d", base, si), h);
    if(!h){ if(s) delete s; return nullptr; }
    if(!s){ s = (TH1D*) h->Clone(Form("cs1_%s_%s_%d", base, tag, ci)); s->SetDirectory(nullptr); }
    else s->Add(h);
  }
  return s;
}

// Project the X axis of a summed TH2 in a Y (usually jet pT) window.
//
// FindBin on the low edge is right; on the high edge it returns the bin that
// STARTS at hi, which must be excluded -- hence the epsilon offsets. Getting
// this wrong silently widens every window by one bin.
inline TH1D* projectY(TH2D *h, double lo, double hi, const char *name)
{
  int b1 = h->GetYaxis()->FindBin(lo + 1e-6);
  int b2 = h->GetYaxis()->FindBin(hi - 1e-6);
  TH1D *p = h->ProjectionX(name, b1, b2);
  p->SetDirectory(nullptr);
  return p;
}
