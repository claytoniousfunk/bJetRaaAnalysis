// Relative b-purity systematic uncertainty by source, pp and the four PbPb classes in a 1x5
// shared-border row (pp, 50-80, 30-50, 10-30, 0-10%), drawn over the reported jet-pT range.
// Shared by systematics_CSVv2_2D_PF.C (after deriving the sources) and plotBPuritySys_CSVv2_PF.C
// (from the stored purity files), so the two cannot drift apart (2026-11-16). The drawing itself
// is drawSysRow() of sysStyle.h, the one look of every by-source systematics figure.
//
//   TH1D *tot[5], *src[5][PSF_NMAX];   // relative (fraction, not %), nullptr = not drawn in that panel
//   drawPuritySysRow(tot, src, nSrc, srcKey, 80., 200., "#it{p}_{T}^{rel} #times CSVv2 fit", pdfName);
#pragma once
#include "sysStyle.h"

const int PSF_NMAX = 8;

inline void drawPuritySysRow(TH1D *tot[5], TH1D *src[5][PSF_NMAX], int nSrc, const char *const *srcKey,
                             double xLo, double xHi, const char *fitText, const TString &pdf,
                             const char *const *labels = nullptr, bool binTicks = false)   // labels: nullptr = sysStyle.h; binTicks: drawSysRow
{
  const char *name[5] = {"pp", "50-80%", "30-50%", "10-30%", "0-10%"};
  SysPanel P[5];
  for(int i = 0; i < 5; i++){ P[i].name = name[i]; P[i].tot = tot[i]; for(int s = 0; s < nSrc && s < SYS_NMAX; s++) P[i].src[s] = src[i][s]; }
  drawSysRow(5, P, nSrc, srcKey, labels, xLo, xHi, Form("b purity (%s)", fitText),
             "anti-#it{k}_{T} #it{R} = 0.4 PF jets,  #sqrt{#it{s}_{NN}} = 5.02 TeV", pdf, "Internal", binTicks);
}
