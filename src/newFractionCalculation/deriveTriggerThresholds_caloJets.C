// Jet trigger turn-on curves and the pT above which each trigger sample can be
// used, for calo jets. These thresholds set the stitching edges and the binning
// for the jets-per-Z measurement, so this macro runs first and writes them out.
//
// NO MINBIAS IN THE SPECTRUM. MinBias appears here only as a trigger reference
// on the PbPb side: measuring an efficiency needs an unbiased denominator, which
// is a different job from supplying yield to the stitched spectrum.
//
// --------------------------------------------------------------- PbPb ------
// PbPb_scan.C fills the turn-on histograms with the SAME bits it gates on
// (HLT_HICsAK4PFJet{60,80,100}Eta1p5_v1 -- its useCaloJetsOverride only swaps
// the JEC files, there is no calo trigger substitution), so what is recorded is
// what is used:
//   eff(Jet60)         = h_leadJetPt_jet60_C / h_leadJetPt_all_C   in MinBias
//                        (ungated, h_nEventsNoJetTrigSel mean 0) -- ABSOLUTE
//   eff(Jet80 | Jet60) = h_leadJetPt_jet60_and_jet80_C / h_leadJetPt_jet60_C
//                        in the Jet60 file -- bootstrap
//   eff(Jet100| Jet80) = h_leadJetPt_jet80_and_jet100_C / h_leadJetPt_jet80_C
//                        in the Jet80 file -- bootstrap
//
// ----------------------------------------------------------------- pp ------
// pp CANNOT use those histograms: pp_scan.C fills them from HLT_HIAK4PFJet*_v1
// unconditionally (lines 940-981) while the gate substitutes
// HLT_HIAK4CaloJet*_v1 under useCaloJetsOverride (lines 723-763), so the
// recorded bits are not the trigger the sample is gated on -- hence
// h_leadJetPt_all (32.3M) > h_leadJetPt_jet60 (29.6M) in the Jet60 calo file.
//
// Other datasets do not work as denominators either. SingleMuon applies no jet
// trigger but a tight 15 GeV muon selects a different cross-section, and that
// ratio falls monotonically from 2.06 at 68 GeV to 0.70 at 280 GeV -- process
// mismatch, not a turn-on. pp MinBias is unbiased but its scan cut jets at 80
// GeV with 1-268 counts/bin, giving a ratio scattering over 2500-35000.
//
// What works is the SAME dataset: the three pp calo files scanned identical
// events (h_vz = 3.626675e8 in all three, checked at run time) and differ only
// in which gate ran, so
//   eff(JetB)/eff(JetA) = h_leadJetPt_all(JetB file) / h_leadJetPt_all(JetA file)
// is a genuine relative efficiency. Jet100/Jet60 is computed as a closure on
// the other two.
//
// LIMIT OF THE pp RESULT: every pp curve is relative, so Jet60's own absolute
// efficiency is NOT determined here; the pp floor is taken where the curves are
// mutually flat. Fixing that needs a pp MinBias calo rescan at the usual 20 GeV
// jet cut, or calo trigger bits added to the pp_scan.C turn-on fills.
//
// FORWARD JETS. The gated pp files keep large sub-threshold populations (4.1e6
// events with a 20-40 GeV leading jet in the Jet60 file): HLT jets span the
// whole calorimeter while the offline leading jet here is inside |eta| < 1.6,
// so the trigger fires on a forward jet the analysis never sees. Every curve
// therefore starts from a non-zero FLOOR rather than from 0, and rises in two
// stages. For stitching that is fine -- the requirement is that a sample's
// yield be proportional to the truth, i.e. that the curve be FLAT, not that it
// reach unity.
//
// PRESCALES. Each curve plateaus at 1/prescale, not 1 (pp Jet60 sits 6.42x
// below Jet80), so every curve is divided by its own plateau, measured as
// sum(num)/sum(den) over a window above the turn-on -- an efficiency averages
// that way, never as a mean of per-bin ratios. The plateau and its statistical
// error are printed so they can be checked against the recorded prescales.
//
// ------------------------------------------------------------- THRESHOLD ---
// PRIMARY criterion is flatness, which is what stitching actually needs: the
// lowest bin edge above which every point sits within 2 sigma of the target and
// the weighted mean of the points clears it. No shape is assumed.
//
// Only bins with a relative error below 5% take part. Past that the statistics
// are gone and a single low fluctuation in the tail vetoes every edge beneath
// it -- that is what put the 10-30% Jet80/Jet60 threshold at 195 GeV on a curve
// flat from 108, and left 30-50% with none at all. The top of the usable range
// is printed as "eval to".
//
// A single erf is fitted too, but only as an overlaid DIAGNOSTIC, because it
// does not describe these curves: with the forward-jet floor the shape is
// floor + turn-on, and a floorless erf returned chi2/ndf of 70, 61 and 42 on
// the PbPb bootstraps. The fit here carries a floor parameter and its chi2/ndf
// is printed, so a bad model is visible instead of being quoted silently. Where
// the two estimates disagree by more than 20 GeV the curve is flagged.
//
// A bare bin scan ("every bin above threshold exceeds 0.98") was tried first and
// is not usable: one downward fluctuation truncates the scan, which returned
// 230 GeV for a PbPb Jet60 curve that plateaus by 90, and disagreed with the fit
// by up to 150 GeV.
//
// WEAK TEST. On a noisy curve any threshold passes, because 2 sigma spans
// everything -- the peripheral Jet60 anchor has errors up to +-1.9. Those are
// flagged weakTest and no threshold is written for them; they need more MinBias
// statistics, not a smaller number.
//
// Usage, from src/newFractionCalculation/:
//   root -l -b -q 'deriveTriggerThresholds_caloJets.C'

#include "../../headers/plotting/plotStyle.h"
#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TColor.h"
#include "TMath.h"
#include "TSystem.h"
#include <cstdio>
#include <vector>

const char *repo = "/home/clayton/Analysis/code/bJetRaaAnalysis";

// nominal target first: it is the one written out and used downstream
const double targets[2] = {0.98, 0.95};

const double minCountsPerBin = 120.;   // per bin, numerator AND denominator
const double nSigmaBelow     = 2.0;    // a point this far below target counts as low
const double nSigmaCatastrophic = 4.0; // a point this far below target vetoes on its own
const double maxRelErrForPower = 0.25; // median rel. error above which the test has no power

const char *clsLabel[5] = {"", "0-10%", "10-30%", "30-50%", "50-80%"};
const char *clsTag[5]   = {"", "0to10pct", "10to30pct", "30to50pct", "50to80pct"};

struct TrigDef {
  const char *name;
  double nominal;
  double fitLo, fitHi;     // erf diagnostic range
  double platLo, platHi;   // plateau normalization window
};
TrigDef trigs[3] = {
  {"Jet60",   60.,  50., 160.,  110., 200.},
  {"Jet80",   80.,  55., 200.,  150., 260.},
  {"Jet100", 100.,  70., 240.,  190., 300.},
};

static bool ok = true;

static TFile* openOrFail(TString p)
{
  TFile *f = TFile::Open(p);
  if(!f || f->IsZombie()){ printf("ERROR: cannot open %s\n", p.Data()); ok = false; return nullptr; }
  return f;
}

static TH1D* getOrFail(TFile *f, const char *name, const char *newName)
{
  if(!f) return nullptr;
  TH1D *h = nullptr;
  f->GetObject(name, h);
  if(!h){ printf("ERROR: %s missing from %s\n", name, f->GetName()); ok = false; return nullptr; }
  TH1D *c = (TH1D*) h->Clone(newName);
  c->SetDirectory(nullptr);
  return c;
}

// Variable rebinning: fine through the turn-on, wide across the plateau.
//
// A uniform 5 GeV axis left the plateaus scattering by several percent, which is
// what made the threshold scan sensitive to single outliers. The rise still has
// to be resolved -- the threshold is read off it, so its bin width is the
// threshold granularity -- while the plateau only has to be flat, so its bins
// can be several times wider. Both widths stay on the 5 GeV grid of the input
// histograms, so every new edge falls on an existing boundary.
//
// The rise width adapts to the curve: it grows until both numerator and
// denominator carry minCountsPerBin through the turn-on, capped at maxRiseWidth
// so a threshold never becomes coarser than that.
const double maxRiseWidth = 20., maxPlateauWidth = 50.;
const double drawMax = 300.;   // common upper edge, so every curve fills the frame

static std::vector<double> buildEdges(TH1D *num, TH1D *den, const TrigDef &t,
                                      double &riseW, double &platW)
{
  std::vector<double> edges;
  const double grid = den->GetXaxis()->GetBinWidth(1);
  const double lo = 40., hi = t.platHi;

  // counts per grid bin through the rise, whichever of the two is thinner
  double sn = 0., sd = 0.; int nb = 0;
  for(int b = 1; b <= den->GetNbinsX(); b++){
    double x = den->GetBinCenter(b);
    if(x < t.fitLo || x > t.platLo) continue;
    sn += num->GetBinContent(b); sd += den->GetBinContent(b); nb++;
  }
  double perBin = (nb > 0) ? TMath::Min(sn, sd)/nb : 0.;
  int g = 1;
  while(g*grid < maxRiseWidth && perBin*g < minCountsPerBin) g++;
  riseW = g*grid;
  platW = TMath::Min(4.*riseW, maxPlateauWidth);

  for(double x = lo; x < t.platLo - 1e-6; x += riseW) edges.push_back(x);
  // plateau bins, with platHi kept as an edge so the normalization window is
  // exactly [platLo, platHi], then carried on to a common drawMax so the
  // plateaus of all three triggers stay on one frame
  for(double x = t.platLo; x < hi - 1e-6; ){
    edges.push_back(x);
    x = (x + platW > hi - 1e-6) ? hi : x + platW;
  }
  edges.push_back(hi);
  for(double x = hi + platW; x < drawMax + 1e-6; x += platW) edges.push_back(x);
  if(edges.back() < drawMax - 1e-6) edges.push_back(drawMax);
  return edges;
}

struct Curve {
  TH1D  *eff = nullptr;
  TF1   *fit = nullptr;
  double plateau = -1., plateauRelErr = -1.;
  double thr[2]    = {-1., -1.};   // PRIMARY, flatness
  double thrErf[2] = {-1., -1.};   // diagnostic
  double erfChi2ndf = -1.;
  double flatChi2ndf = -1.;
  double evalHi = -1.;            // top of the range where the curve still has power
  double riseW = 0., platW = 0.;  // bin widths through the turn-on and the plateau
  bool   good = false;
  bool   weakTest = false;
  TString note;
};

static Curve makeCurve(TH1D *num, TH1D *den, const TrigDef &t, bool binomial, const char *name)
{
  Curve c;
  if(!num || !den) return c;

  double riseW = 0., platW = 0.;
  std::vector<double> edges = buildEdges(num, den, t, riseW, platW);
  TH1D *n = (TH1D*) num->Rebin((int)edges.size()-1, Form("%s_num", name), &edges[0]);
  TH1D *d = (TH1D*) den->Rebin((int)edges.size()-1, Form("%s_den", name), &edges[0]);
  n->SetDirectory(nullptr); d->SetDirectory(nullptr);
  c.riseW = riseW; c.platW = platW;

  const double eps = 0.01;
  double sn = n->Integral(n->FindBin(t.platLo + eps), n->FindBin(t.platHi - eps));
  double sd = d->Integral(d->FindBin(t.platLo + eps), d->FindBin(t.platHi - eps));
  if(sd <= 0. || sn <= 0.){
    c.note = Form("empty plateau window [%.0f,%.0f]", t.platLo, t.platHi);
    delete n; delete d; return c;
  }
  c.plateau = sn / sd;
  c.plateauRelErr = binomial ? sqrt(TMath::Max(0., c.plateau * (1. - c.plateau) / sd)) / c.plateau
                             : sqrt(1./sn + 1./sd);

  c.eff = (TH1D*) n->Clone(name);
  c.eff->SetDirectory(nullptr);
  c.eff->Divide(n, d, 1, 1, binomial ? "B" : "");
  c.eff->Scale(1. / c.plateau);

  // ---- PRIMARY: lowest edge above which the curve is flat at the target ----
  //
  // Only bins that still carry precision take part. In the tail the statistics
  // die (the 10-30% Jet80/Jet60 curve is flat from ~108 GeV but has a
  // 0.9045 +- 0.0357 point at 192 GeV, 2.7 sigma low), and letting such a point
  // veto every candidate edge beneath it pushed that class to 195 GeV and left
  // 30-50% with no threshold at all. A bin with no power should abstain, not
  // decide.
  const double maxRelErrForFlat = 0.05;
  int bLast = c.eff->FindBin(t.platHi - eps);
  int bEval = 0;
  for(int b = 1; b <= bLast; b++){
    double y = c.eff->GetBinContent(b), e = c.eff->GetBinError(b);
    if(y <= 0. || e <= 0.) continue;
    if(e/y <= maxRelErrForFlat) bEval = b;
  }
  if(bEval > 0) c.evalHi = c.eff->GetXaxis()->GetBinUpEdge(bEval);

  for(int it = 0; it < 2; it++){
    double found = -1.;
    for(int bStart = 1; bStart <= bEval; bStart++){
      double edge = c.eff->GetXaxis()->GetBinLowEdge(bStart);
      if(edge < t.fitLo) continue;

      int    nPts = 0, nLow = 0;
      double sw = 0., swy = 0., chi2 = 0.;
      bool   catastrophic = false;
      for(int b = bStart; b <= bEval; b++){
        double y = c.eff->GetBinContent(b), e = c.eff->GetBinError(b);
        if(y <= 0. || e <= 0.) continue;
        if(e/y > maxRelErrForFlat) continue;      // no power: abstain
        nPts++;
        double w = 1./(e*e);
        sw += w; swy += w*y;
        chi2 += (y - 1.)*(y - 1.)/(e*e);
        double pull = (targets[it] - y)/e;        // sigma BELOW target
        if(pull > nSigmaBelow)     nLow++;
        if(pull > nSigmaCatastrophic) catastrophic = true;
      }
      if(nPts < 3) break;
      double mean = sw > 0. ? swy/sw : -1.;
      // A hard "no point may sit 2 sigma low" veto is wrong: among ~25 points
      // one such fluctuation is expected, and exactly one (0.9045 +- 0.0357 at
      // 192 GeV) dragged the 10-30% Jet80 threshold to 195 GeV on a curve flat
      // from 110. Tolerate the expected rate, but let any single catastrophic
      // point veto -- inside the rise the points are many sigma low, which is
      // what must still be rejected.
      double allowedLow = TMath::Max(1., 0.05 * nPts);
      if(!catastrophic && nLow <= allowedLow && mean >= targets[it]){
        found = edge;
        // Reported, never required: with 0.1% errors a real 2% creep that stays
        // inside the target gives chi2/ndf of tens, and demanding a small chi2
        // would drive the threshold far above where the trigger is usable.
        if(it == 0) c.flatChi2ndf = nPts > 1 ? chi2/(nPts - 1) : -1.;
        break;
      }
    }
    c.thr[it] = found;
  }
  c.good = (c.thr[0] > 0.);

  // A curve whose usable range ends before the turn-on has finished cannot
  // locate a threshold at all -- that is the PbPb Jet60 anchor in every class:
  // the MinBias calo sample gives 0.944 +- 4.2% at 90 GeV in 0-10% and
  // 0.586 +- 15% at 90 GeV in 50-80%, so "consistent with the plateau" is true
  // of almost any edge. Say so instead of printing a bare "none".
  if(c.evalHi > 0. && c.evalHi < t.nominal + 20.){
    c.weakTest = true;
    c.note = Form("usable range ends at %.0f GeV, before the turn-on completes:"
                  " these statistics cannot locate a threshold", c.evalHi);
  }
  else if(c.evalHi <= 0.){
    c.weakTest = true;
    c.note = "no bin reaches the precision needed to test flatness";
  }

  // is the test able to reject anything, or are the errors simply huge?
  {
    std::vector<double> rel;
    for(int b = 1; b <= bLast; b++){
      double y = c.eff->GetBinContent(b), e = c.eff->GetBinError(b);
      if(c.eff->GetXaxis()->GetBinLowEdge(b) < t.nominal) continue;
      if(y <= 0. || e <= 0.) continue;
      rel.push_back(e/y);
    }
    if(!rel.empty()){
      std::sort(rel.begin(), rel.end());
      double med = rel[rel.size()/2];
      if(med > maxRelErrForPower){
        c.weakTest = true;
        c.note = Form("median rel. error %.0f%% above %.0f GeV: the flatness test has no power here",
                      100.*med, t.nominal);
      }
    }
  }

  // ---- DIAGNOSTIC: floor + erf, so the floor is not absorbed into the rise --
  c.fit = new TF1(Form("%s_erf", name),
                  "[2]+(1.-[2])*0.5*(1.+TMath::Erf((x-[0])/(TMath::Sqrt(2.)*[1])))",
                  t.fitLo, t.fitHi);
  c.fit->SetParameters(t.nominal, 10., 0.1);
  c.fit->SetParLimits(1, 0.5, 80.);
  c.fit->SetParLimits(2, 0., 0.9);
  if(c.eff->Fit(c.fit, "QNR") == 0 && c.fit->GetParameter(1) > 0.){
    double mu = c.fit->GetParameter(0), sg = c.fit->GetParameter(1), f0 = c.fit->GetParameter(2);
    if(c.fit->GetNDF() > 0) c.erfChi2ndf = c.fit->GetChisquare()/c.fit->GetNDF();
    for(int i = 0; i < 2; i++){
      double frac = (targets[i] - f0)/(1. - f0);
      if(frac > 0. && frac < 1.)
        c.thrErf[i] = 5.*ceil((mu + TMath::Sqrt(2.)*sg*TMath::ErfInverse(2.*frac - 1.))/5. - 1e-9);
    }
  }

  delete n; delete d;
  return c;
}

static void styleCurve(TH1D *h, int i)
{
  const char *hex[3] = {okabeHex[5], okabeHex[1], okabeHex[2]};
  const int   mk[3]  = {markFilledCircle, markFilledSquare, markFilledDiamond};
  styleH(h, hex[i % 3], mk[i % 3], 0.9);
}
static int curveColor(int i)
{
  const char *hex[3] = {okabeHex[5], okabeHex[1], okabeHex[2]};
  return TColor::GetColor(hex[i % 3]);
}

static void report(const char *tag, const Curve &c)
{
  if(!c.eff){ printf("    %-16s UNUSABLE: %s\n", tag, c.note.Data()); return; }
  printf("    %-16s plateau %.4f +-%4.1f%%  bins %2.0f/%2.0f GeV   %.0f%% at ", tag,
         c.plateau, 100.*c.plateauRelErr, c.riseW, c.platW, 100.*targets[0]);
  if(c.good) printf("%5.0f GeV", c.thr[0]); else printf("  none");
  printf("   %.0f%% at ", 100.*targets[1]);
  if(c.thr[1] > 0.) printf("%5.0f GeV", c.thr[1]); else printf("  none");
  printf("   [eval to %4.0f; erf %5.0f, chi2/ndf %6.1f; flat chi2/ndf %5.2f]\n",
         c.evalHi, c.thrErf[0], c.erfChi2ndf, c.flatChi2ndf);
  if(c.weakTest) printf("      WEAK: %s -- no threshold written\n", c.note.Data());
  else if(c.good && c.thrErf[0] > 0. && fabs(c.thrErf[0] - c.thr[0]) > 20.)
    printf("      NOTE: flatness and erf disagree by %.0f GeV; the erf is a poor model here\n",
           fabs(c.thrErf[0] - c.thr[0]));
}

void deriveTriggerThresholds_caloJets()
{
  initPlotStyle();
  ok = true;

  TString figDir = Form("%s/figures/JetsPerZ", repo);
  gSystem->mkdir(figDir, kTRUE);
  TString outDir = "./rootFiles/JetsPerZ";
  gSystem->mkdir(outDir, kTRUE);

  // ------------------------------------------------------------------ pp ---
  printf("\n=== pp calo jet triggers: RELATIVE efficiency, same dataset, same events ===\n");

  TFile *fpp[3]; TH1D *hpp[3];
  for(int i = 0; i < 3; i++){
    fpp[i] = openOrFail(Form("%s/rootFiles/scanningOuput/pp/pp_HighEGJet_caloJets_%sHLT_mu12_pTmu-15to999_tight_deltaR-40_jetTrkMaxFilter_WDecayFilter_2026-9-15.root", repo, trigs[i].name));
    if(!fpp[i]) return;
    hpp[i] = getOrFail(fpp[i], "h_leadJetPt_all", Form("ppAll_%d", i));
    if(!hpp[i]) return;
  }
  {
    double v0 = -1.;
    for(int i = 0; i < 3; i++){
      TH1D *v = getOrFail(fpp[i], "h_vz", Form("ppVz_%d", i));
      if(!v) return;
      if(i == 0) v0 = v->Integral();
      else if(fabs(v->Integral() - v0) > 1e-6*v0){
        printf("ERROR: the pp calo files are not the same scan (h_vz differs); cross-file ratios invalid\n");
        return;
      }
    }
    printf("  same-events check: h_vz = %.6e in all three files\n", v0);
  }

  Curve ppC[3];
  ppC[1] = makeCurve(hpp[1], hpp[0], trigs[1], false, "effPP_Jet80_relJet60");
  ppC[2] = makeCurve(hpp[2], hpp[1], trigs[2], false, "effPP_Jet100_relJet80");
  Curve ppCross = makeCurve(hpp[2], hpp[0], trigs[2], false, "effPP_Jet100_relJet60");

  report("Jet80 / Jet60",  ppC[1]);
  report("Jet100 / Jet80", ppC[2]);
  report("Jet100 / Jet60", ppCross);
  printf("    Jet60 itself has NO absolute reference in these files; the pp floor below is\n"
         "    where the relative curves are mutually flat, not a measured Jet60 turn-on.\n");

  double ppFloor = -1.;
  if(ppC[1].good && ppC[2].good) ppFloor = TMath::Max(ppC[1].thr[0], ppC[2].thr[0]);

  // ---------------------------------------------------------------- PbPb ---
  printf("\n=== PbPb calo jet triggers, per centrality class ===\n");

  TFile *f_MB  = openOrFail(Form("%s/rootFiles/scanningOuput/PbPb/PbPb_MinBias_Part1_caloJets_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_2026-9-14.root", repo));
  TFile *f_j60 = openOrFail(Form("%s/rootFiles/scanningOuput/PbPb/PbPb_HardProbes_caloJets_Jet60HLT_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_2026-9-15.root", repo));
  TFile *f_j80 = openOrFail(Form("%s/rootFiles/scanningOuput/PbPb/PbPb_HardProbes_caloJets_Jet80HLT_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_2026-9-15.root", repo));
  if(!ok) return;
  { TH1D *g = nullptr; f_MB->GetObject("h_nEventsNoJetTrigSel", g);
    if(g && g->GetMean() > 0.5){ printf("ERROR: the PbPb MinBias scan is jet-triggered, not a valid reference\n"); return; } }

  Curve pbC[5][3];
  for(int c = 1; c <= 4; c++){
    printf("\n  %s\n", clsLabel[c]);
    TH1D *mbAll = getOrFail(f_MB,  Form("h_leadJetPt_all_C%d", c),   Form("mbAll_%d", c));
    TH1D *mb60  = getOrFail(f_MB,  Form("h_leadJetPt_jet60_C%d", c), Form("mb60_%d", c));
    TH1D *d60   = getOrFail(f_j60, Form("h_leadJetPt_jet60_C%d", c), Form("d60_%d", c));
    TH1D *n80   = getOrFail(f_j60, Form("h_leadJetPt_jet60_and_jet80_C%d", c), Form("n80_%d", c));
    TH1D *d80   = getOrFail(f_j80, Form("h_leadJetPt_jet80_C%d", c), Form("d80_%d", c));
    TH1D *n100  = getOrFail(f_j80, Form("h_leadJetPt_jet80_and_jet100_C%d", c), Form("n100_%d", c));
    if(!ok) return;

    pbC[c][0] = makeCurve(mb60, mbAll, trigs[0], true, Form("effPbPb_Jet60_C%d", c));
    pbC[c][1] = makeCurve(n80,  d60,   trigs[1], true, Form("effPbPb_Jet80_C%d", c));
    pbC[c][2] = makeCurve(n100, d80,   trigs[2], true, Form("effPbPb_Jet100_C%d", c));

    report("Jet60 (absolute)", pbC[c][0]);
    report("Jet80 / Jet60",    pbC[c][1]);
    report("Jet100 / Jet80",   pbC[c][2]);

    for(int i = 1; i < 3; i++)
      if(pbC[c][i].good && pbC[c][i-1].good && !pbC[c][i-1].weakTest &&
         pbC[c][i].thr[0] < pbC[c][i-1].thr[0])
        printf("      NOTE: %s clears its target at %.0f GeV, below its own reference %s (%.0f GeV);\n"
               "            a bootstrap only means anything on the reference's plateau\n",
               trigs[i].name, pbC[c][i].thr[0], trigs[i-1].name, pbC[c][i-1].thr[0]);
  }

  // ------------------------------------------------------------- figures ---
  // iStart skips leading curves while keeping the color tied to the trigger,
  // so dropping one does not recolor the others
  auto drawPanel = [&](Curve *cs, const char *const *labels, int n, int iStart,
                       const char *head1, const char *head2, bool small)
  {
    double ts = small ? 0.055 : 0.045;
    TLegend *leg = makeLegend(small ? 0.38 : 0.48, small ? 0.15 : 0.18,
                              small ? 0.97 : 0.94, small ? 0.39 : 0.38, small ? 0.045 : 0.034);
    bool first = true;
    for(int i = iStart; i < n; i++){
      if(!cs[i].eff) continue;
      styleCurve(cs[i].eff, i);
      if(first){
        cs[i].eff->SetTitle("");
        cs[i].eff->SetMinimum(0.); cs[i].eff->SetMaximum(1.35);
        cs[i].eff->GetXaxis()->SetRangeUser(40., 300.);
        cs[i].eff->GetXaxis()->SetTitle("leading jet #it{p}_{T} [GeV]");
        cs[i].eff->GetXaxis()->SetTitleSize(ts); cs[i].eff->GetXaxis()->SetLabelSize(ts*0.9);
        cs[i].eff->GetYaxis()->SetTitle("trigger efficiency (plateau-normalized)");
        cs[i].eff->GetYaxis()->SetTitleSize(ts*0.95); cs[i].eff->GetYaxis()->SetLabelSize(ts*0.9);
        cs[i].eff->GetYaxis()->SetTitleOffset(small ? 1.30 : 1.50);
        cs[i].eff->Draw("E1");
        first = false;
      }
      else cs[i].eff->Draw("E1 same");
    }
    TLine *l1 = new TLine(40., 1., 300., 1.); l1->SetLineStyle(7); l1->Draw();
    TLine *l2 = new TLine(40., targets[0], 300., targets[0]);
    l2->SetLineStyle(3); l2->SetLineColor(kGray+2); l2->Draw();
    for(int i = iStart; i < n; i++){
      if(!cs[i].good || cs[i].weakTest) continue;
      TLine *lt = new TLine(cs[i].thr[0], 0., cs[i].thr[0], targets[0]);
      lt->SetLineStyle(2); lt->SetLineColor(curveColor(i)); lt->Draw();
    }
    for(int i = iStart; i < n; i++){
      if(!cs[i].eff) continue;
      cs[i].eff->Draw("E1 same");
      TString lab = cs[i].weakTest ? Form("%s: too noisy", labels[i])
                  : cs[i].good     ? Form("%s: %.0f GeV", labels[i], cs[i].thr[0])
                                   : Form("%s: none", labels[i]);
      leg->AddEntry(cs[i].eff, lab, "lp");
    }
    leg->Draw();
    TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(small ? 0.055 : 0.036);
    if(head1) la.DrawLatex(small ? 0.21 : 0.15, small ? 0.88 : 0.945, head1);
    if(head2) la.DrawLatex(0.15, 0.900, head2);
  };

  {
    Curve set[3] = {ppC[1], ppC[2], ppCross};
    const char *labs[3] = {"Jet80 / Jet60", "Jet100 / Jet80", "Jet100 / Jet60"};
    TCanvas *cv = new TCanvas("c_turnOn_pp", "", 700, 800);
    cv->SetLeftMargin(0.15); cv->SetBottomMargin(0.12); cv->SetTopMargin(0.13); cv->SetRightMargin(0.05);
    drawPanel(set, labs, 3, 0, "pp 5.02 TeV, calo jets, relative trigger efficiency",
              "same events, the gate is the only difference", false);
    cv->SaveAs(Form("%s/triggerTurnOn_pp_caloJets.pdf", figDir.Data()));
    delete cv;
  }
  {
    const char *labs[3] = {"Jet60 (abs.)", "Jet80 / Jet60", "Jet100 / Jet80"};
    TCanvas *cv = new TCanvas("c_turnOn_PbPb", "", 700, 800);
    cv->Divide(2, 2, 0.001, 0.001);
    for(int c = 1; c <= 4; c++){
      cv->cd(c);
      gPad->SetLeftMargin(0.16); gPad->SetBottomMargin(0.13); gPad->SetTopMargin(0.08); gPad->SetRightMargin(0.04);
      // start at 1: the Jet60 absolute curve is too noisy to locate a threshold
      // and only cluttered the panel (still computed and reported in the printout)
      drawPanel(pbC[c], labs, 3, 1, Form("PbPb %s", clsLabel[c]), nullptr, true);
    }
    cv->SaveAs(Form("%s/triggerTurnOn_PbPb_caloJets.pdf", figDir.Data()));
    delete cv;
  }

  // ------------------------------------------------------------- outputs ---
  TString txtPath = Form("%s/triggerThresholds_caloJets.txt", outDir.Data());
  FILE *fp = fopen(txtPath, "w");
  if(fp){
    fprintf(fp, "# calo jet trigger thresholds from deriveTriggerThresholds_caloJets.C\n");
    fprintf(fp, "# criterion: flat within %.0f sigma of %.2f of the plateau, upward\n", nSigmaBelow, targets[0]);
    fprintf(fp, "# pp values are RELATIVE: no absolute pp reference exists in these files\n");
    fprintf(fp, "# system class trigger reference threshold plateau\n");
    if(ppC[1].good && !ppC[1].weakTest) fprintf(fp, "pp inclusive Jet80 Jet60 %.0f %.6f\n", ppC[1].thr[0], ppC[1].plateau);
    if(ppC[2].good && !ppC[2].weakTest) fprintf(fp, "pp inclusive Jet100 Jet80 %.0f %.6f\n", ppC[2].thr[0], ppC[2].plateau);
    // Jet60 gets NO row: its absolute efficiency was never measured here, and a
    // row would let the stitching macro treat it as validated -- and tie it with
    // Jet100 at the same edge. The mutual-flatness floor is a comment only.
    if(ppFloor > 0.)
      fprintf(fp, "# pp Jet60 has no measured absolute threshold; mutual-flatness floor %.0f GeV\n", ppFloor);
    const char *refs[3] = {"MinBias", "Jet60", "Jet80"};
    for(int c = 1; c <= 4; c++)
      for(int i = 0; i < 3; i++)
        if(pbC[c][i].good && !pbC[c][i].weakTest)
          fprintf(fp, "PbPb %s %s %s %.0f %.6f\n", clsTag[c], trigs[i].name, refs[i],
                  pbC[c][i].thr[0], pbC[c][i].plateau);
    fclose(fp);
  }

  TString rootPath = Form("%s/triggerTurnOn_caloJets.root", outDir.Data());
  TFile *wf = TFile::Open(rootPath, "recreate");
  if(ppC[1].eff)  ppC[1].eff->Write("effPP_Jet80_relJet60");
  if(ppC[2].eff)  ppC[2].eff->Write("effPP_Jet100_relJet80");
  if(ppCross.eff) ppCross.eff->Write("effPP_Jet100_relJet60");
  for(int c = 1; c <= 4; c++)
    for(int i = 0; i < 3; i++)
      if(pbC[c][i].eff) pbC[c][i].eff->Write(Form("effPbPb_%s_C%d", trigs[i].name, c));
  wf->Close();

  printf("\n  pp floor (mutual flatness of the relative curves): %.0f GeV\n", ppFloor);
  printf("  figures   : %s/triggerTurnOn_{pp,PbPb}_caloJets.pdf\n", figDir.Data());
  printf("  thresholds: %s\n  curves    : %s\n", txtPath.Data(), rootPath.Data());
}
