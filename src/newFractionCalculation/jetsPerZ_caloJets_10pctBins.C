// Stitched PbPb calo-jet spectra and jets-per-Z in 10%-wide centrality bins.
//
// Built on the manual-JEC PbPb scans of 2026-09-19/20, which are the first calo
// scans whose jet pT comes from the AK4Calo files rather than the forest's jtpt
// (the forest corrects calo jets with the AK4PF payload; see
// headers/config/config_PbPb.h and the 2026-09-18 jet-by-jet check).
//
// CENTRALITY. These scans use centrality_ultraFineCentBins.h: 18 slices of
// hiBin width 10, i.e. 5% each, covering 0-90%. Confirmed against the files
// themselves -- C1..C18 hold 2.31e6-2.47e6 MinBias events each. A 10% bin is
// therefore exactly two adjacent slices, and class k (0-based) is slices
// 2k+1 and 2k+2, giving nine bins: 0-10 ... 80-90%.
//   NOTE coarseCent.h is NOT used here. Its 4-class mapping stops at slice 16
//   (0-80%) and would silently drop 80-90%.
//
// STITCHING. MinBias below the Jet80 threshold, Jet80 up to the Jet100
// threshold, Jet100 above. Jet100 sets the scale, so the result is in Jet100
// trigger units, not per event; that cancels in any class-to-class ratio and in
// PbPb/pp as long as the Jet100 prescale is class-independent.
//
// The two normalization factors are measured in a COMMON pT window and averaged
// over classes, not fitted per class in each class's own overlap window. That
// was established earlier: a per-class window sits where combinatorial jets
// inflate the central-class MinBias yield and drags the factor down ~9%, which
// imports a centrality dependence into a quantity that is a property of the
// datasets. Both factors are printed per class so the spread is visible.
//
// THRESHOLDS are derived here rather than read from
// triggerThresholds_caloJets.txt, because that file was produced on the OLD
// forest-JEC files and the manual JEC moves every threshold up. Jet80 is
// measured ABSOLUTELY against MinBias (which is unbiased and present here), and
// Jet100 by bootstrap off the Jet80 sample. Criterion: on a 10 GeV grid, the
// lowest edge from which the plateau-normalized efficiency stays within
// nSigma of 1 all the way to the ceiling. There is no Jet60 manual-JEC scan, so
// Jet60 is not used at all -- it was already excluded from the spectra.
//
// JETS PER Z. N_Z must come from the SingleMuon dataset, per centrality class.
// The full-statistics PbPb SingleMuon scan (19385 Z) was run with 4 classes
// only, and the ultra-fine SingleMuon scans in this repo are partial (1432 Z,
// a different luminosity), so 10% N_Z does not exist yet. The user chose to
// wait for a SingleMuon rescan with ultraFineCentBins rather than approximate
// it. This macro therefore:
//   - always produces the stitched spectra, which need no Z at all;
//   - fills in jets-per-Z automatically as soon as a SingleMuon file matching
//     singleMuonGlob below appears with >= minZPerClass Z per class.
// Same for the PbPb/pp ratio: the pp trigger samples still carry the forest's
// PF-sized calo correction and are ~40% low, so the ratio is skipped until the
// manual-JEC pp scans exist. Set the paths below when they do.
//
// NOT DONE: no unfolding (no calo response matrices yet) and no fake-jet
// subtraction (the estimate is FastJet on PF candidates, with no calo
// equivalent). Stated on every figure.
//
// Usage: root -l -b -q 'jetsPerZ_caloJets_10pctBins.C'
// Run from: src/newFractionCalculation/

#include "../../headers/functions/triggerStitch.h"
#include "../../headers/functions/divideByBinwidth.h"
#include "../../headers/plotting/plotStyle.h"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TSystem.h"
#include "TMath.h"
#include <cstdio>
#include <vector>

const char *repo = "/home/clayton/Analysis/code/bJetRaaAnalysis";

const char *f_MB_path   = "rootFiles/scanningOuput/PbPb/PbPb_MinBias_Part1_caloJets_manualJEC_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_2026-9-19_ultraFineCentBins.root";
const char *f_j80_path  = "rootFiles/scanningOuput/PbPb/PbPb_HardProbes_caloJets_manualJEC_Jet80HLT_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_2026-9-19_ultraFineCentBins.root";
const char *f_j100_path = "rootFiles/scanningOuput/PbPb/PbPb_HardProbes_caloJets_manualJEC_Jet100HLT_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_2026-9-20_ultraFineCentBins.root";

// filled in when the rescan lands; the macro checks and reports
const char *singleMuonGlob = "rootFiles/scanningOuput/PbPb/PbPb_SingleMuon_*ultraFineCentBins*.root";
const long  minZPerClass   = 400;   // below this a 10% N_Z is too noisy to use

const char *outDir  = "/home/clayton/Analysis/code/bJetRaaAnalysis/figures/JetsPerZ";
const char *outFile = "./rootFiles/JetsPerZ/jetsPerZ_caloJets_10pctBins.root";

// ------------------------------------------------------------- parameters ---
const int    NCls        = 9;      // 10% bins, 0-10 .. 80-90
const double ptFloor     = 50.;    // the scans' own jet pT cut
const double ptCeiling   = 500.;
const double etaRange    = 3.2;    // |eta| < 1.6

// threshold criterion
const double effTarget   = 0.98;
const double nSigmaBelow = 2.0;
const double turnOnRebin = 10.;    // GeV, coarser than the 5 GeV scan grid
// Plateau reference windows. Jet80-vs-MinBias must be read where MinBias still
// has events (it is a per-event sample with a hard 5 GeV grid and empties out
// above ~250 GeV), while the Jet100 bootstrap runs on the jet-triggered sample
// and can use a genuinely high window.
const double platLo_abs  = 130., platHi_abs  = 200.;
const double platLo_boot = 200., platHi_boot = 300.;

// common normalization windows: above the fake-dominated region, inside the
// plateau of both samples in every class
const double kWinMB_lo = 150., kWinMB_hi = 200.;   // MinBias -> Jet80
const double kWin80_lo = 200., kWin80_hi = 280.;   // Jet80   -> Jet100

// binning
const double targetRelErr  = 0.05;
const double minWidthFrac  = 0.15;   // resolution floor; see JER (sigma ~ 0.11-0.22)
const double maxBinWidth   = 80.;

static const char *clsLabel[NCls] = {"0-10%","10-20%","20-30%","30-40%","40-50%",
                                     "50-60%","60-70%","70-80%","80-90%"};
static const char *clsTag[NCls]   = {"0to10","10to20","20to30","30to40","40to50",
                                     "50to60","60to70","70to80","80to90"};

// sum the two 5% slices making up 10% class k; nullptr if either is missing
static TH1D* classSum(TFile *f, const char *fmt, int k, const char *name)
{
  TH1D *a = nullptr, *b = nullptr;
  f->GetObject(Form(fmt, 2*k + 1), a);
  f->GetObject(Form(fmt, 2*k + 2), b);
  if(!a || !b){ printf("ERROR: %s or %s missing\n", Form(fmt, 2*k+1), Form(fmt, 2*k+2)); return nullptr; }
  TH1D *s = (TH1D*) a->Clone(name);
  s->SetDirectory(nullptr);
  s->Add(b);
  return s;
}

// Lowest 10-GeV edge from which num/den stays within nSigma of effTarget*plateau
// up to ptCeiling. Returns -1 if no such edge exists.
//
// Two tolerances, both carried over from deriveTriggerThresholds_caloJets.C and
// both necessary here: the MinBias denominator runs out of events well before
// 300 GeV, so (a) only bins measured to maxRelErrForFlat or better can veto an
// edge -- a bin with 30% errors says nothing about flatness -- and (b) a few
// mild downward fluctuations are expected among many bins and are allowed,
// while a single catastrophic one is not. Without these a lone high-pT
// fluctuation pushes the threshold to 400+ GeV, which is what happened first.
static double turnOnThreshold(TH1D *num, TH1D *den, double pLo, double pHi,
                              double &plateau, const char *label, bool verbose)
{
  const double maxRelErrForFlat = 0.05;
  const double nSigmaCatastrophic = 4.0;
  TH1D *n = (TH1D*) num->Clone(Form("n_%s", label)); n->SetDirectory(nullptr);
  TH1D *d = (TH1D*) den->Clone(Form("d_%s", label)); d->SetDirectory(nullptr);
  int group = (int) (turnOnRebin / n->GetBinWidth(1) + 0.5);
  if(group > 1 && n->GetNbinsX() % group == 0){ n->Rebin(group); d->Rebin(group); }

  const double eps = 0.01;
  // The plateau window has to sit where the denominator still has events. In
  // the most peripheral class MinBias is empty above ~130 GeV, so slide the
  // window down until it holds enough to measure a plateau, and say which one
  // was used -- a plateau read off 3 events is not a plateau.
  // A window is only accepted with enough events behind it. Sliding further
  // down would land inside the turn-on itself and return a "plateau" of 0.06,
  // which then yields a meaningless threshold -- better to report the class as
  // unmeasured (-2) and let the caller substitute a measured neighbour.
  const double minPlateauCounts = 200.;
  double nP = 0., dP = 0., usedLo = pLo, usedHi = pHi;
  bool measured = false;
  for(double scale : {1.0, 0.85}){
    usedLo = pLo*scale; usedHi = pHi*scale;
    if(usedLo < ptFloor) break;
    nP = n->Integral(n->FindBin(usedLo + eps), n->FindBin(usedHi - eps));
    dP = d->Integral(d->FindBin(usedLo + eps), d->FindBin(usedHi - eps));
    if(dP >= minPlateauCounts && nP > 0.){ measured = true; break; }
  }
  if(!measured){
    printf("    %-22s NOT MEASURABLE: plateau window holds %.0f events (need %.0f)\n",
           label, dP, minPlateauCounts);
    delete n; delete d; return -2.;
  }
  plateau = nP / dP;

  int bLo = n->FindBin(ptFloor + eps), bHi = n->FindBin(ptCeiling - eps);
  // how many bins above the candidate edge carry enough precision to judge
  int nTestable = 0;
  for(int bb = bLo; bb <= bHi; bb++){
    double dd = d->GetBinContent(bb);
    if(dd <= 0.) continue;
    double e = n->GetBinContent(bb)/dd, se = sqrt(TMath::Max(e*(1.-e)/dd, 1e-12));
    if(e > 0. && se/e <= maxRelErrForFlat) nTestable++;
  }
  const int allowedLow = (int) TMath::Max(1., 0.05*nTestable);

  double best = -1.;
  for(int b = bLo; b <= bHi; b++){
    int nLow = 0; bool catastrophic = false;
    for(int bb = b; bb <= bHi; bb++){
      double nn = n->GetBinContent(bb), dd = d->GetBinContent(bb);
      if(dd <= 0.) continue;                      // no data, not a failure
      double e = nn/dd, se = sqrt(TMath::Max(e*(1.-e)/dd, 1e-12));
      if(e <= 0.) continue;
      if(e < effTarget*plateau - nSigmaCatastrophic*se){ catastrophic = true; break; }
      if(se/e > maxRelErrForFlat) continue;       // too imprecise to veto
      if(e < effTarget*plateau - nSigmaBelow*se) nLow++;
    }
    if(!catastrophic && nLow <= allowedLow){ best = n->GetXaxis()->GetBinLowEdge(b); break; }
  }
  if(verbose) printf("    %-22s plateau %.4g on [%.0f,%.0f]   threshold %.0f GeV\n",
                     label, plateau, usedLo, usedHi, best);
  delete n; delete d;
  return best;
}

void jetsPerZ_caloJets_10pctBins()
{
  initPlotStyle();
  gSystem->mkdir(outDir, kTRUE);
  gSystem->mkdir("./rootFiles/JetsPerZ", kTRUE);

  TFile *fMB  = TFile::Open(Form("%s/%s", repo, f_MB_path));
  TFile *f80  = TFile::Open(Form("%s/%s", repo, f_j80_path));
  TFile *f100 = TFile::Open(Form("%s/%s", repo, f_j100_path));
  if(!fMB || fMB->IsZombie() || !f80 || f80->IsZombie() || !f100 || f100->IsZombie()){
    printf("ERROR: cannot open one of the PbPb scans\n"); return;
  }
  { TH1D *g = nullptr; fMB->GetObject("h_nEventsNoJetTrigSel", g);
    if(g && g->GetMean() > 0.5){ printf("ERROR: the MinBias scan is jet-triggered; not a valid reference\n"); return; } }

  printf("\n  PbPb calo jets, manual JEC, 10%% centrality bins (two 5%% slices each)\n");

  // ------------------------------------------------------- per-class input ---
  TH1D *spMB[NCls], *sp80[NCls], *sp100[NCls];
  double nEvt[NCls];
  double thr80[NCls], thr100[NCls];

  printf("\n  trigger turn-on, derived on these files (Jet80 absolute vs MinBias, Jet100 off Jet80)\n");
  for(int k = 0; k < NCls; k++){
    spMB[k]  = classSum(fMB,  "h_inclRecoJetPt_C%d", k, Form("mb_%d", k));
    sp80[k]  = classSum(f80,  "h_inclRecoJetPt_C%d", k, Form("j80_%d", k));
    sp100[k] = classSum(f100, "h_inclRecoJetPt_C%d", k, Form("j100_%d", k));
    TH1D *vz = classSum(fMB,  "h_vz_C%d", k, Form("vz_%d", k));
    if(!spMB[k] || !sp80[k] || !sp100[k] || !vz) return;
    nEvt[k] = vz->Integral();
    delete vz;

    // MinBias is a rate; the trigger samples stay in counts and are scaled to
    // Jet100 units below
    spMB[k]->Scale(1. / nEvt[k]);

    printf("  %s  (N_evt = %.3g)\n", clsLabel[k], nEvt[k]);
    TH1D *mbAll = classSum(fMB, "h_leadJetPt_all_C%d",   k, Form("mbAll_%d", k));
    TH1D *mb80  = classSum(fMB, "h_leadJetPt_jet80_C%d", k, Form("mb80_%d", k));
    TH1D *d80   = classSum(f80, "h_leadJetPt_jet80_C%d", k, Form("d80_%d", k));
    TH1D *n100  = classSum(f80, "h_leadJetPt_jet80_and_jet100_C%d", k, Form("n100_%d", k));
    if(!mbAll || !mb80 || !d80 || !n100) return;
    double pl;
    thr80[k]  = turnOnThreshold(mb80, mbAll, platLo_abs,  platHi_abs,  pl, Form("Jet80 abs %s",  clsTag[k]), true);
    thr100[k] = turnOnThreshold(n100, d80,   platLo_boot, platHi_boot, pl, Form("Jet100/Jet80 %s", clsTag[k]), true);
    delete mbAll; delete mb80; delete d80; delete n100;

    // In the most peripheral classes MinBias has too few high-pT jets to show a
    // Jet80 plateau at all. The measured thresholds fall monotonically with
    // peripherality (120 -> 100 GeV), so carrying the last measured value
    // forward is the conservative choice: it sits ABOVE the true turn-on.
    if(thr80[k] == -2.){
      if(k == 0){ printf("ERROR: Jet80 unmeasurable in the most central class\n"); return; }
      thr80[k] = thr80[k-1];
      printf("    Jet80 abs %-12s using %.0f GeV carried over from %s (not measured here)\n",
             clsTag[k], thr80[k], clsLabel[k-1]);
    }
    if(thr100[k] == -2.){
      if(k == 0){ printf("ERROR: Jet100 unmeasurable in the most central class\n"); return; }
      thr100[k] = thr100[k-1];
      printf("    Jet100    %-12s using %.0f GeV carried over from %s (not measured here)\n",
             clsTag[k], thr100[k], clsLabel[k-1]);
    }
    if(thr80[k] < 0. || thr100[k] < 0.){ printf("ERROR: no threshold in %s\n", clsLabel[k]); return; }
    if(thr100[k] <= thr80[k]){
      printf("ERROR: %s has Jet100 threshold %.0f <= Jet80 %.0f; stitch order broken\n",
             clsLabel[k], thr100[k], thr80[k]); return;
    }
  }

  // ------------------------------------------ common-window scale factors ---
  // measured on the same summed class histograms used for the spectra
  auto commonK = [&](TH1D **low, TH1D **high, double lo, double hi, const char *label){
    const double eps = 0.01;
    double sw = 0., swk = 0.;
    printf("\n  %s yield ratio in the common window [%.0f,%.0f]\n", label, lo, hi);
    for(int k = 0; k < NCls; k++){
      double nl = low[k] ->Integral(low[k] ->FindBin(lo + eps), low[k] ->FindBin(hi - eps));
      double nh = high[k]->Integral(high[k]->FindBin(lo + eps), high[k]->FindBin(hi - eps));
      if(nl <= 0. || nh <= 0.){ printf("    %-7s empty window\n", clsLabel[k]); continue; }
      // MinBias is already a per-event rate, so its "counts" for the error are
      // recovered from the rate and the event count
      double effLow = (low == spMB) ? nl * nEvt[k] : nl;
      double kk = nh/nl, sk = kk*sqrt(1./TMath::Max(effLow, 1.) + 1./nh);
      printf("    %-7s %12.4g +- %.2g\n", clsLabel[k], kk, sk);
      sw += 1./(sk*sk); swk += kk/(sk*sk);
    }
    if(sw <= 0.) return -1.;
    double kAvg = swk/sw;
    printf("    weighted average applied to every class: %.4g +- %.2g\n", kAvg, sqrt(1./sw));
    return kAvg;
  };

  double k80to100 = commonK(sp80, sp100, kWin80_lo, kWin80_hi, "Jet80 -> Jet100");
  double kMBto80  = commonK(spMB, sp80,  kWinMB_lo, kWinMB_hi, "MinBias -> Jet80");
  if(k80to100 <= 0. || kMBto80 <= 0.){ printf("ERROR: a normalization window failed\n"); return; }

  // ----------------------------------------------------------- stitching ---
  printf("\n  stitching (Jet100 sets the scale; result is in Jet100 trigger units)\n");
  TH1D *stitched[NCls];
  for(int k = 0; k < NCls; k++){
    std::vector<StitchSample> s = {
      {"MinBias", ptFloor,   spMB[k],  kMBto80},
      {"Jet80",   thr80[k],  sp80[k],  k80to100},
      {"Jet100",  thr100[k], sp100[k], 0.}
    };
    printf("  %s: MinBias %.0f-%.0f, Jet80 %.0f-%.0f, Jet100 %.0f-%.0f\n", clsLabel[k],
           ptFloor, thr80[k], thr80[k], thr100[k], thr100[k], ptCeiling);
    stitched[k] = stitchTriggerSamples(s, Form("stitched_%s", clsTag[k]), 50., ptCeiling, true);
    if(!stitched[k]){ printf("ERROR: stitch failed in %s\n", clsLabel[k]); return; }
  }

  // ------------------------------------------------- per-class rebinning ---
  printf("\n  per-class binning (target %.0f%% stat, resolution floor %.0f%% of pT)\n",
         100.*targetRelErr, 100.*minWidthFrac);
  TH1D *spec[NCls];
  for(int k = 0; k < NCls; k++){
    std::vector<double> e = stitchDeriveEdges({stitched[k]}, ptFloor, ptCeiling,
                                              targetRelErr, minWidthFrac, maxBinWidth);
    if(e.size() < 2){ printf("ERROR: no binning in %s\n", clsLabel[k]); return; }
    spec[k] = rebinTo(stitched[k], (int)e.size() - 1, e.data(), Form("spec_%s", clsTag[k]));
    spec[k]->Scale(1. / etaRange);
    divideByBinwidth(spec[k]);
    printf("  %-7s %2d bins, %.0f-%.0f GeV\n", clsLabel[k], (int)e.size() - 1, e.front(), e.back());
  }

  // ------------------------------------------------------------ jets per Z ---
  printf("\n  jets per Z\n");
  TString found = gSystem->GetFromPipe(Form("ls -1t %s/%s 2>/dev/null | head -1", repo, singleMuonGlob));
  found = found.Strip(TString::kBoth);
  bool haveZ = false;
  if(found.Length()){
    TFile *fZ = TFile::Open(found);
    if(fZ && !fZ->IsZombie()){
      long zMin = -1;
      for(int k = 0; k < NCls; k++){
        TH1D *m = classSum(fZ, "h_dimuonMass_C%d", k, Form("z_%d", k));
        if(!m){ zMin = -1; break; }
        long nz = (long) m->Integral(m->FindBin(75.01), m->FindBin(104.99));
        zMin = (zMin < 0) ? nz : TMath::Min(zMin, nz);
        delete m;
      }
      printf("    candidate SingleMuon file: %s\n", found.Data());
      printf("    fewest Z in any 10%% class: %ld (need >= %ld)\n", zMin, minZPerClass);
      haveZ = (zMin >= minZPerClass);
      fZ->Close();
    }
  }
  if(!haveZ){
    printf("    SKIPPED: no ultra-fine SingleMuon scan with enough Z per 10%% class.\n");
    printf("    The full-statistics scan (19385 Z) has 4 classes only, and the ultra-fine\n");
    printf("    scans here are partial (1432 Z on a different luminosity). Rerun the PbPb\n");
    printf("    SingleMuon scan with centrality_ultraFineCentBins.h and this macro fills\n");
    printf("    in jets-per-Z automatically -- the spectra below are already final.\n");
  }

  // --------------------------------------------------------------- output ---
  TFile *out = TFile::Open(outFile, "RECREATE");
  for(int k = 0; k < NCls; k++){
    spec[k]->Write(Form("spectrum_%s", clsTag[k]));
    stitched[k]->Write(Form("stitched_fine_%s", clsTag[k]));
  }
  TH1D *hThr = new TH1D("thresholds", "Jet80 and Jet100 thresholds;class;GeV", NCls, 0, NCls);
  for(int k = 0; k < NCls; k++){ hThr->SetBinContent(k+1, thr80[k]); hThr->GetXaxis()->SetBinLabel(k+1, clsLabel[k]); }
  hThr->Write();
  out->Close();

  // ---------------------------------------------------------------- figure ---
  {
    TCanvas *c = new TCanvas("c_spec", "", 700, 800);
    c->SetLeftMargin(0.17); c->SetBottomMargin(0.12);
    c->SetTopMargin(0.17);  c->SetRightMargin(0.05);   // three caption lines
    c->SetLogy();   // the spectra span five decades across 50-500 GeV

    TLegend *leg = makeLegend(0.58, 0.40, 0.94, 0.80, 0.030);
    // Nine classes, seven usable Okabe-Ito colors (4 is the unreadable yellow),
    // so color repeats after seven and the marker shape changes with it. No
    // (color, marker) pair repeats, which a plain k%8 did not guarantee -- it
    // gave 0-10% and 80-90% the same black circle.
    const int useHex[7] = {0, 1, 2, 3, 5, 6, 7};
    const int mkCycle[2] = {markFilledCircle, markFilledSquare};
    double yMax = 0., yMin = 1e30;
    for(int k = 0; k < NCls; k++){
      for(int b = 1; b <= spec[k]->GetNbinsX(); b++){
        double v = spec[k]->GetBinContent(b);
        if(v > 0.){ yMax = TMath::Max(yMax, v); yMin = TMath::Min(yMin, v); }
      }
    }
    for(int k = 0; k < NCls; k++){
      styleH(spec[k], okabeHex[useHex[k % 7]], mkCycle[k / 7], 0.8);
      if(k == 0){
        spec[k]->SetTitle(""); spec[k]->SetStats(0);
        spec[k]->SetMinimum(yMin * 0.4); spec[k]->SetMaximum(yMax * 40.);
        spec[k]->GetXaxis()->SetTitle("jet #it{p}_{T} [GeV]");
        spec[k]->GetXaxis()->SetTitleSize(0.045); spec[k]->GetXaxis()->SetLabelSize(0.040);
        spec[k]->GetYaxis()->SetTitle("d#it{N}_{jet}/d#it{p}_{T}d#eta  [Jet100 units]");
        spec[k]->GetYaxis()->SetTitleSize(0.042); spec[k]->GetYaxis()->SetLabelSize(0.040);
        spec[k]->GetYaxis()->SetTitleOffset(1.65);
        spec[k]->Draw("E1");
      }
      else spec[k]->Draw("E1 same");
      leg->AddEntry(spec[k], clsLabel[k], "lp");
    }
    leg->Draw();

    TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.036);
    la.DrawLatex(0.17, 0.955, "PbPb 5.02 TeV, calo jets, anti-#it{k}_{T} #it{R} = 0.4");
    la.SetTextSize(0.030);
    la.DrawLatex(0.17, 0.917, "manual AK4Calo JEC, no unfolding, no fake jet subtraction");
    la.DrawLatex(0.17, 0.883, "MinBias + Jet80 + Jet100, 10% centrality bins");

    c->SaveAs(Form("%s/stitchedSpectra_caloJets_10pctBins.pdf", outDir));
    delete c;
  }

  printf("\n  results: %s\n  figure : %s/stitchedSpectra_caloJets_10pctBins.pdf\n", outFile, outDir);
}
