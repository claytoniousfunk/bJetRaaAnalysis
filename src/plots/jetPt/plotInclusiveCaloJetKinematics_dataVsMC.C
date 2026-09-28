// Inclusive calo jet kinematics, data vs MC: pT in pp and in each PbPb
// centrality class, eta and phi in pp.
//
// Companion to plotMuTaggedJetKinematics_pp_caloJets_dataVsMC.C, which does the
// same for MUON-TAGGED jets in pp only. This one drops the tag requirement and
// adds PbPb, so it exercises the jet reconstruction itself rather than the
// tagging.
//
// WHAT IS AND IS NOT HERE, and why. The four pieces come from scans that book
// different histogram sets, and the gaps are not symmetric:
//
//   pp    pT, eta, phi   data h_inclRecoJetPt{,_inclRecoJetEta,_inclRecoJetPhi}
//                        MC   the same three names, PYTHIA_DiJet_caloJets
//   PbPb  pT only        data the same h_inclRecoJetPt_C*
//                        MC   h_recoJetPt_all_C*, from the RESPONSE scan
//
//   PbPb eta and phi CANNOT be done from anything on disk. Every PbPb calo MC
//   file is a response scan, and those book no inclusive reco eta or phi at all
//   -- the only angular histogram they carry is
//   h_matchedRecoJetPtOverGenJetPt_genJetEta_*, which is a JES map against GEN
//   eta for matched jets above 100 GeV, not a reco eta distribution. Closing
//   this needs a PYTHIAHYDJET calo scan booking
//   h_inclRecoJetPt_inclRecoJetEta and _inclRecoJetPhi, as the pp PYTHIA scan
//   already does. See the printout at the end.
//
// JET ENERGY SCALE -- READ BEFORE INTERPRETING THE pp pT PANEL.
//   PbPb: both sides carry manual JEC, so the comparison is like for like.
//   pp:   the data scans are manual-JEC, and the only pp calo MC scan booking
//         inclusive histograms (2026-9-15) is NOT -- it uses the forest jtpt,
//         which applies the AK4PF payload to calo jets (closure ~0.78 in pp).
//         There is no manual-JEC pp calo MC scan with these histograms: the
//         manual-JEC PYTHIA calo files are response scans and book none of
//         them. So the pp pT ratio is dominated by that ~20% scale difference
//         and measures the JEC bug, not the modelling. eta and phi are far less
//         sensitive to it -- the pT floor selects a slightly different jet
//         population on each side, but the shapes are set by acceptance.
//
// TRIGGER. Jet100 only, above the threshold deriveTriggerThresholds_caloJets.C
// validated for that sample and class (triggerThresholds_caloJets.txt: pp 90,
// PbPb 140 / 125 / 120 / 115 GeV). Using one trigger above its own measured
// plateau avoids re-implementing the stitch, whose splice normalization would
// otherwise become part of what this figure is testing.
//
// SHAPES ONLY. The MC carries pThat cross-section weights with no counterpart
// in the data counts, so everything is normalized to unit area over the
// compared range.
//
// Usage: root -l -b -q plotInclusiveCaloJetKinematics_dataVsMC.C
// Run from: src/plots/jetPt/

#include <cmath>
#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TSystem.h"
#include "../../../headers/functions/divideByBinwidth.h"
#include "../../../headers/plotting/plotStyle.h"
#include "../../../headers/plotting/ratioPanel.h"
#include "../../../headers/plotting/shapeComparison.h"

namespace {

const char *base = "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/";

const char *ppData = "pp/pp_HighEGJet_caloJets_manualJEC_Jet100HLT_mu12_pTmu-15to999_tight_"
                     "deltaR-40_jetTrkMaxFilter_WDecayFilter_2026-9-25.root";
const char *ppMC   = "PYTHIA/PYTHIA_DiJet_caloJets_pThat-15_mu12_pTmu-15to999_tight_"
                     "vzReweight_jetTrkMaxFilter_removeHYDJETjet0p35CutOnGen_WDecayFilter_2026-9-15.root";

const char *pbData = "PbPb/PbPb_HardProbes_caloJets_manualJEC_Jet100HLT_mu12_pTmu-15to999_tight_"
                     "jetTrkMaxFilter_WDecayFilter_2026-9-20.root";
const char *pbMC   = "PYTHIAHYDJET/PYTHIAHYDJET_response_DiJet_caloJets_manualJEC_pThat-15_"
                     "mu12_pTmu-15_tight_vzReweight_hiBinReweight_hiBinShift-10_jetTrkMaxFilter_"
                     "doPThatCorrelationFilterTight_2026-9-25.root";

const char *outDir = "../../../figures/jetKinematics/inclusiveCaloJets/";

// Jet100 plateau thresholds, triggerThresholds_caloJets.txt
const double ppThr = 90.;
const int    NCls  = 4;
const char  *cls  [NCls] = {"C1", "C2", "C3", "C4"};
const char  *clsLab[NCls] = {"0-10%", "10-30%", "30-50%", "50-80%"};
const double clsThr[NCls] = {140., 125., 120., 115.};

// One PbPb axis for all four classes, starting above the highest threshold
// (C1, 140) so the four panels are directly comparable rather than each
// starting somewhere different. Every edge is on the scans' 5 GeV grid.
const double edgePb[]  = {150, 175, 200, 250, 300, 400, 500};
const int    nEdgePb   = sizeof(edgePb)/sizeof(double) - 1;
const double edgePp[]  = {90, 100, 120, 150, 200, 250, 300, 400, 500};
const int    nEdgePp   = sizeof(edgePp)/sizeof(double) - 1;

const int etaRebin = 4;   // 160 bins of 0.02 -> 0.08
const int phiRebin = 2;   // 100 bins -> 50

double integralRange(TH1D *h, double lo, double hi)
{
  return h->Integral(h->FindBin(lo + 1e-6), h->FindBin(hi - 1e-6));
}

// inclusive pT by projecting the pT-vs-eta 2D over all eta. Used on both sides
// in pp: the MC scan books no 1D h_inclRecoJetPt (only 2D maps against flavour,
// eta, phi and the tag variables), and projecting the same 2D on both sides
// also guarantees the two get an identical selection.
TH1D* ptFromEta2D(TFile *f, const char *name)
{
  TH2D *H = nullptr; f->GetObject("h_inclRecoJetPt_inclRecoJetEta", H);
  if(!H){ printf("  ERROR: h_inclRecoJetPt_inclRecoJetEta missing for %s\n", name); return nullptr; }
  TH1D *p = H->ProjectionX(name, 0, -1);   // all eta, over- and underflow included
  p->SetDirectory(nullptr);
  return p;
}

// unit-area, bin-width-divided copy of h on the given edges
TH1D* shapeOf(TH1D *src, int nEdge, const double *edge, double normLo, const char *name)
{
  // a null here means the caller's histogram was absent; say so rather than
  // dropping the panel with no explanation, which is how the pp pT comparison
  // first came back missing
  if(!src){ printf("  ERROR: no input histogram for %s\n", name); return nullptr; }
  TH1D *r = rebinTo(src, nEdge, edge, name);
  const double n = integralRange(r, normLo, edge[nEdge]);
  if(n <= 0.){ printf("  ERROR: %s has no jets in %.0f-%.0f GeV\n", name, normLo, edge[nEdge]);
               delete r; return nullptr; }
  r->Scale(1./n);
  divideByBinwidth(r);
  return r;
}

// eta or phi for jets above ptMin, unit area
TH1D* angular(TFile *f, const char *hname, int rebin, double ptMin, const char *name)
{
  TH2D *H = nullptr; f->GetObject(hname, H);
  if(!H){ printf("  ERROR: %s missing\n", hname); return nullptr; }
  TH1D *p = H->ProjectionY(name, H->GetXaxis()->FindBin(ptMin + 1e-6), H->GetNbinsX() + 1);
  p->SetDirectory(nullptr);
  p->Rebin(rebin);
  if(p->Integral() <= 0.){ printf("  ERROR: %s empty above %.0f GeV\n", hname, ptMin);
                           delete p; return nullptr; }
  p->Scale(1./p->Integral());
  divideByBinwidth(p);
  return p;
}

} // namespace

void plotInclusiveCaloJetKinematics_dataVsMC()
{
  initPlotStyle();
  gSystem->mkdir(outDir, kTRUE);

  TFile *fPpD = TFile::Open(Form("%s%s", base, ppData));
  TFile *fPpM = TFile::Open(Form("%s%s", base, ppMC));
  TFile *fPbD = TFile::Open(Form("%s%s", base, pbData));
  TFile *fPbM = TFile::Open(Form("%s%s", base, pbMC));
  for(auto *f : {fPpD, fPpM, fPbD, fPbM})
    if(!f || f->IsZombie()){ printf("ERROR: a required file did not open\n"); return; }

  const TString kt = "anti-#it{k}_{T} #it{R} = 0.4 calo jets, |#it{#eta}| < 1.6";

  // ======================= pp: pT, eta, phi =================================
  printf("\n================ pp, inclusive calo jets ================\n");
  printf("  data manual-JEC, MC forest JEC -- the pT ratio below is dominated\n"
         "  by that scale difference, NOT by modelling.\n");

  TH1D *ppPtDraw = ptFromEta2D(fPpD, "ppPt_raw_data");
  TH1D *ppPtMraw = ptFromEta2D(fPpM, "ppPt_raw_mc");
  TH1D *ppPtD = shapeOf(ppPtDraw, nEdgePp, edgePp, ppThr, "ppPt_data");
  TH1D *ppPtM = shapeOf(ppPtMraw, nEdgePp, edgePp, ppThr, "ppPt_mc");

  if(ppPtD && ppPtM){
    const double c = drawShapeComparison(
      ppPtD, ppPtM, "#it{p}_{T}^{jet} [GeV]", "1/N dN/d#it{p}_{T} [GeV^{-1}]",
      true, edgePp[0], edgePp[nEdgePp],
      {Form("pp 5.02 TeV, %s", kt.Data()),
       "HLT_HIAK4CaloJet100, inclusive jets",
       "MC has NO manual JEC: ratio shows the JEC offset"},
      Form("%sppCaloJetPt_dataVsMC.pdf", outDir));
    printf("\n  pT (unit area %.0f-%.0f GeV), chi2/ndf = %.2f\n", ppThr, edgePp[nEdgePp], c);
    printf("  %-11s %10s\n", "jet pT", "data/MC");
    for(int b = 1; b <= ppPtD->GetNbinsX(); b++)
      printf("  %3.0f-%-7.0f %10.3f\n", edgePp[b-1], edgePp[b],
             ppPtM->GetBinContent(b) > 0. ? ppPtD->GetBinContent(b)/ppPtM->GetBinContent(b) : -1.);
  }

  struct Ang { const char *h; int rb; const char *ax; const char *yt; double lo, hi; const char *out; };
  const Ang angs[2] = {
    {"h_inclRecoJetPt_inclRecoJetEta", etaRebin, "#it{#eta}^{jet}", "1/N dN/d#it{#eta}",
     -1.6, 1.6, "ppCaloJetEta_dataVsMC.pdf"},
    {"h_inclRecoJetPt_inclRecoJetPhi", phiRebin, "#it{#phi}^{jet}", "1/N dN/d#it{#phi}",
     -TMath::Pi(), TMath::Pi(), "ppCaloJetPhi_dataVsMC.pdf"}
  };

  for(const auto &a : angs){
    TH1D *dD = angular(fPpD, a.h, a.rb, ppThr, Form("%s_data", a.h));
    TH1D *dM = angular(fPpM, a.h, a.rb, ppThr, Form("%s_mc",   a.h));
    if(!dD || !dM) continue;
    const double c = drawShapeComparison(
      dD, dM, a.ax, a.yt, false, a.lo, a.hi,
      {Form("pp 5.02 TeV, %s", kt.Data()),
       "HLT_HIAK4CaloJet100, inclusive jets",
       Form("#it{p}_{T}^{jet} > %.0f GeV", ppThr)},
      Form("%s%s", outDir, a.out));
    printf("  %-6s mean  data %+.4f  MC %+.4f   RMS  data %.4f  MC %.4f   chi2/ndf %.2f\n",
           a.ax, dD->GetMean(), dM->GetMean(), dD->GetRMS(), dM->GetRMS(), c);
  }

  // ======================= PbPb: pT per class ===============================
  printf("\n================ PbPb, inclusive calo jets ================\n");
  printf("  both sides manual-JEC. Jet100 above its per-class plateau.\n");

  //  Two MC definitions, because the choice turns out to dominate the answer.
  //
  //  allJets  h_recoJetPt_all: every reco jet, matched or not.  Formally the
  //           right counterpart to the data, which also contains combinatorial
  //           jets -- except that the MC's fake component is the broken one.
  //           It escapes the pThat correlation filter (the fill happens before
  //           the filter's continue, see PYTHIAHYDJET_scan_response.C), so it is
  //           built from a handful of high-weight events and is 71% of all MC
  //           jets in C1 at 150-175 GeV.
  //  matched  ProjectionX of the matched response 2D: real jets only, and the
  //           part of the MC the filter actually protects.
  //
  //  The difference is not cosmetic: in C1 the data/MC ratio across 150-500 GeV
  //  goes from a factor 5.6 with allJets to 1.7 with matched.
  const bool useMatched[2] = {false, true};
  const char *mcTag[2]     = {"allJets", "matchedOnly"};
  const char *mcLab[2]     = {"PYTHIA+HYDJET (all reco)", "PYTHIA+HYDJET (matched)"};

  for(int v = 0; v < 2; v++){
    printf("\n  ---- MC = %s ----\n", mcTag[v]);
    for(int i = 0; i < NCls; i++){
      TH1D *rawD = nullptr, *rawM = nullptr;
      fPbD->GetObject(Form("h_inclRecoJetPt_%s", cls[i]), rawD);
      if(useMatched[v]){
        TH2D *H = nullptr;
        fPbM->GetObject(Form("h_matchedRecoJetPt_genJetPt_allJets_%s", cls[i]), H);
        if(H){ rawM = (TH1D*) H->ProjectionX(Form("pbMat_%s", cls[i]), 1, H->GetNbinsY());
               rawM->SetDirectory(nullptr); }
      }
      else fPbM->GetObject(Form("h_recoJetPt_all_%s", cls[i]), rawM);
      if(!rawD || !rawM){ printf("  %s: missing histogram, skipped\n", cls[i]); continue; }

      TH1D *d = shapeOf(rawD, nEdgePb, edgePb, edgePb[0], Form("pbPt_data_%s_%d", cls[i], v));
      TH1D *m = shapeOf(rawM, nEdgePb, edgePb, edgePb[0], Form("pbPt_mc_%s_%d",   cls[i], v));
      if(!d || !m) continue;

      const double c = drawShapeComparison(
        d, m, "#it{p}_{T}^{jet} [GeV]", "1/N dN/d#it{p}_{T} [GeV^{-1}]",
        true, edgePb[0], edgePb[nEdgePb],
        {Form("PbPb 5.02 TeV %s, %s", clsLab[i], kt.Data()),
         Form("HLT_HIAK4CaloJet100 (plateau %.0f GeV), inclusive jets", clsThr[i]),
         useMatched[v] ? "MC: gen-matched jets only"
                       : "MC: all reco jets (fake component is unfiltered)"},
        Form("%sPbPbCaloJetPt_%s_%s_dataVsMC.pdf", outDir, cls[i], mcTag[v]),
        "PbPb data", mcLab[v]);

      printf("\n  %s  chi2/ndf = %.2f\n  %-11s %10s\n", clsLab[i], c, "jet pT", "data/MC");
      for(int b = 1; b <= d->GetNbinsX(); b++)
        printf("  %3.0f-%-7.0f %10.3f\n", edgePb[b-1], edgePb[b],
               m->GetBinContent(b) > 0. ? d->GetBinContent(b)/m->GetBinContent(b) : -1.);
    }
  }

  printf("\n  PbPb eta and phi: NOT PRODUCED. Every PbPb calo MC file on disk is a\n"
         "  response scan, and those book no inclusive reco eta or phi. Needs a\n"
         "  PYTHIAHYDJET calo scan booking h_inclRecoJetPt_inclRecoJetEta and\n"
         "  h_inclRecoJetPt_inclRecoJetPhi, which the pp PYTHIA scan already does.\n");
  printf("\n  figures in %s\n", outDir);
}
