#pragma once
// Trigger-sample stitching and statistics-driven binning, shared by the calo-jet
// measurements in src/newFractionCalculation/.
//
// Pulled out of stitchSpectra_caloJets.C when the R_CP version needed the same
// normalization: the splice below decides the absolute scale of every spectrum,
// and two copies of that drifting apart is exactly the failure CLAUDE.md warns
// about with the coarse-class mapping.
//
// Thresholds come from deriveTriggerThresholds_caloJets.C. Only samples that
// file validates should ever be handed to stitchTriggerSamples(), so a trigger
// whose turn-on could not be measured is never silently trusted.

#include "TH1D.h"
#include "TFile.h"
#include "TString.h"
#include <cstdio>
#include <vector>
#include <map>
#include <algorithm>

struct StitchSample {
  const char *trig;
  double      thr;        // use this sample from here upward
  TH1D       *h;
  // >0 overrides the overlap fit with a supplied RAW ratio to the sample above.
  // Needed when the natural overlap window sits somewhere the yields are biased:
  // normalizing MinBias to Jet80 in each centrality class's own threshold window
  // put the window at 85-150 GeV, where combinatorial jets inflate the MinBias
  // side in central events and drag the factor down by ~9%. Measured instead in
  // a common window above the fake-dominated region the factor is the dataset
  // constant it should be (125.6 +- 5.4, 124.5 +- 5.1, 120.3 +- 8.1 in the three
  // best-measured classes), so the caller can supply that and avoid importing a
  // centrality-dependent scale into the sub-threshold region.
  double      fixedK = 0.;
};

// "system/class/trigger" -> threshold. Comment lines are skipped, so a trigger
// recorded only as a comment (pp Jet60, which has no absolute reference) is
// correctly absent rather than read as validated.
inline std::map<TString, double> stitchReadThresholds(const char *path)
{
  std::map<TString, double> m;
  FILE *fp = fopen(path, "r");
  if(!fp){
    printf("ERROR: cannot read %s -- run deriveTriggerThresholds_caloJets.C first\n", path);
    return m;
  }
  char sys[64], cls[64], trg[64], ref[64], line[512];
  double thr, plateau;
  while(fgets(line, sizeof(line), fp)){
    if(line[0] == '#') continue;
    if(sscanf(line, "%63s %63s %63s %63s %lf %lf", sys, cls, trg, ref, &thr, &plateau) == 6)
      m[Form("%s/%s/%s", sys, cls, trg)] = thr;
  }
  fclose(fp);
  return m;
}

// Normalize each sample to the one above it and splice at the thresholds.
// Samples must be sorted by threshold, ascending.
//
// The highest-threshold sample sets the scale because it is the least
// prescaled, and each lower sample is scaled into its units over the overlap
// window where both are known flat. Two consequences worth keeping in mind:
// the result is in the top sample's trigger units, not per event; and that
// scale is common to every centrality class only when the top trigger's
// prescale is the same in all of them (for HLT_HICsAK4PFJet100Eta1p5 it is
// 1.000 everywhere), which is what lets a class-to-class ratio cancel it.
inline TH1D* stitchTriggerSamples(std::vector<StitchSample> &s, const char *name,
                                  double overlapWidth, double ptCeiling, bool verbose)
{
  if(s.empty()) return nullptr;
  const int top = (int)s.size() - 1;
  const double eps = 0.01;

  TH1D *out = (TH1D*) s[top].h->Clone(name);
  out->SetDirectory(nullptr);
  out->Reset();

  std::vector<double> k(s.size(), 1.);
  for(int i = top - 1; i >= 0; i--){
    if(s[i].fixedK > 0.){
      // supplied as a raw ratio to the sample above, so it still has to be
      // carried up into the top sample's units like any fitted factor
      k[i] = s[i].fixedK * k[i+1];
      if(verbose)
        printf("    %-6s -> %-6s : k = %.4g  (raw %.4g supplied, not fitted here)\n",
               s[i].trig, s[i+1].trig, k[i], s[i].fixedK);
      continue;
    }
    double lo = s[i+1].thr, hi = std::min(s[i+1].thr + overlapWidth, ptCeiling);
    double nAbove = s[i+1].h->Integral(s[i+1].h->FindBin(lo + eps), s[i+1].h->FindBin(hi - eps)) * k[i+1];
    double nThis  = s[i].h  ->Integral(s[i].h  ->FindBin(lo + eps), s[i].h  ->FindBin(hi - eps));
    if(nThis <= 0. || nAbove <= 0.){
      printf("    ERROR: %s: empty overlap [%.0f,%.0f] between %s and %s\n",
             name, lo, hi, s[i].trig, s[i+1].trig);
      return nullptr;
    }
    k[i] = nAbove / nThis;
    if(verbose)
      printf("    %-6s -> %-6s over [%.0f,%.0f]: k = %.4g  (%.0f counts, %.1f%% stat)\n",
             s[i].trig, s[i+1].trig, lo, hi, k[i], nThis, 100./sqrt(nThis));
  }

  for(int b = 1; b <= out->GetNbinsX(); b++){
    double x = out->GetBinCenter(b);
    int use = -1;
    for(size_t i = 0; i < s.size(); i++){
      double hi = (i + 1 < s.size()) ? s[i+1].thr : ptCeiling + 1.;
      if(x >= s[i].thr && x < hi){ use = (int)i; break; }
    }
    if(use < 0) continue;
    out->SetBinContent(b, s[use].h->GetBinContent(b) * k[use]);
    out->SetBinError  (b, s[use].h->GetBinError(b)   * k[use]);
  }
  return out;
}

// Statistics-weighted average over centrality classes 1..nCls of the yield ratio
// fHigh/fLow in [lo,hi], for histograms named Form(fmt, c). A ratio between two
// datasets or trigger paths is a property of them and not of centrality, so it
// is measured per class as a check and the average is what gets applied --
// measuring it in each class's own overlap window let centrality leak in (18%
// for Jet80/MinBias, 1.4% for Jet100/Jet80; see jetsPerZ_RCP_stitched_caloJets.C).
// Returns <= 0 on failure.
inline double stitchCommonScale(TFile *fLow, TFile *fHigh, const char *fmt, int nCls,
                                double lo, double hi, const char *label, bool verbose)
{
  const double eps = 0.01;
  double sw = 0., swk = 0.;
  if(verbose) printf("\n  %s yield ratio in the common window [%.0f,%.0f]\n", label, lo, hi);
  for(int c = 1; c <= nCls; c++){
    TH1D *hl = nullptr, *hh = nullptr;
    fLow ->GetObject(Form(fmt, c), hl);
    fHigh->GetObject(Form(fmt, c), hh);
    if(!hl || !hh){ printf("    ERROR: %s missing\n", Form(fmt, c)); return -1.; }
    double nl = hl->Integral(hl->FindBin(lo + eps), hl->FindBin(hi - eps));
    double nh = hh->Integral(hh->FindBin(lo + eps), hh->FindBin(hi - eps));
    if(nl <= 0. || nh <= 0.){ if(verbose) printf("    class %d: empty window\n", c); continue; }
    double k = nh/nl, sk = k*sqrt(1./nl + 1./nh);
    if(verbose) printf("    class %d  %8.2f +- %.2f   (lower-sample counts %.0f)\n", c, k, sk, nl);
    sw += 1./(sk*sk); swk += k/(sk*sk);
  }
  if(sw <= 0.) return -1.;
  if(verbose) printf("    weighted average applied to every class: %.2f +- %.2f\n", swk/sw, sqrt(1./sw));
  return swk/sw;
}

// Merge bins upward until EVERY histogram passed in carries targetRelErr, and
// never emit a bin narrower than minWidthFrac of its own pT -- the resolution
// floor. Pass one histogram for a single spectrum, two for a ratio, where the
// sparser of the two (usually the peripheral reference) sets the binning.
inline std::vector<double> stitchDeriveEdges(std::vector<TH1D*> hs, double floorPt, double ceilPt,
                                             double targetRelErr, double minWidthFrac, double maxBinWidth)
{
  std::vector<double> edges;
  if(hs.empty() || !hs[0]) return edges;
  const double eps = 0.01;
  TH1D *h0 = hs[0];
  int b0 = h0->FindBin(floorPt + eps), bN = h0->FindBin(ceilPt - eps);
  edges.push_back(h0->GetXaxis()->GetBinLowEdge(b0));

  std::vector<double> sum(hs.size(), 0.), err2(hs.size(), 0.);
  double startEdge = edges.back();
  for(int b = b0; b <= bN; b++){
    bool allPrecise = true;
    for(size_t i = 0; i < hs.size(); i++){
      if(!hs[i]) continue;
      sum[i]  += hs[i]->GetBinContent(b);
      err2[i] += hs[i]->GetBinError(b) * hs[i]->GetBinError(b);
      if(!(sum[i] > 0. && sqrt(err2[i])/sum[i] <= targetRelErr)) allPrecise = false;
    }
    double upper = h0->GetXaxis()->GetBinUpEdge(b);
    bool resolved = (upper - startEdge >= minWidthFrac * 0.5 * (startEdge + upper));
    bool tooWide  = (upper - startEdge >= maxBinWidth);
    if((allPrecise && resolved) || tooWide){
      edges.push_back(upper);
      startEdge = upper;
      std::fill(sum.begin(), sum.end(), 0.);
      std::fill(err2.begin(), err2.end(), 0.);
    }
  }
  if(edges.size() < 2) return std::vector<double>();
  if(edges.back() < ceilPt) edges.back() = ceilPt;   // fold the sparse tail in
  return edges;
}
