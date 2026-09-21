// R_CP of jets per Z, PbPb only, from MinBias alone -- calo jets.
//
// The companion macro (stitchSpectra_caloJets.C) stitches trigger samples and
// divides by pp. This one drops both: every class comes from the same ungated
// MinBias scan, so there is no trigger turn-on to respect and the spectra reach
// down to the scan's own 50 GeV jet cut instead of the 85-110 GeV trigger
// floors. The price is statistics, which is why this is a ratio of centrality
// classes rather than a measurement against pp.
//
// ------------------------------------------------------- NORMALISATION -----
// For class c,   jets per Z (c)  =  N_jets(c) / N_Z(c)
// with N_jets the RAW MinBias counts and N_Z the raw Z yield of that class,
// muon-efficiency corrected. Both count objects observed in the same centrality
// slice, so
//      N_jets(c)/N_Z(c) = [L_MinBias / L_SingleMuon] x [sigma_jet(c)/sigma_Z(c)]
// and the luminosity factor -- the only thing the two datasets do not share --
// is class independent, so it cancels exactly in
//      R_CP(c) = [jets per Z](c) / [jets per Z](reference class).
//
// Do NOT divide the MinBias counts by N_evt here. That would strip out the
// centrality-slice population which is precisely what pairs with the same
// population inside N_Z(c); the two must be left to cancel against each other.
//
// BINNING is common to all classes, driven by the peripheral reference. That
// class appears in every ratio and has the fewest jets, so it sets the useful
// precision; giving each class its own binning (as the pp measurement does)
// would be pointless when the denominator is shared. Bins merge until the
// reference carries the target relative error and are never narrower than the
// resolution floor.
//
// ------------------------------------------------------------- CAVEATS -----
// NO FAKE-JET SUBTRACTION, and here it matters. Combinatorial jets from the
// underlying event are a central-collision effect concentrated exactly where
// this macro gains its reach: on the PF estimate the fake fraction of the
// MinBias yield is 17.4% at 60-80 GeV in 0-10%, 6.8% at 80-100 and 4.2% at
// 100-130, against ~3% in the peripheral classes. Unsubtracted, that inflates
// the central classes at low pT and biases R_CP UPWARD there. There is no calo
// fake-jet estimate to subtract (the existing one is FastJet on PF candidates),
// so the low-pT end of the central curves should be read as an upper bound
// until one exists.
//
// No unfolding either (no calo response matrices yet).
//
// Usage, from src/newFractionCalculation/:
//   root -l -b -q 'jetsPerZ_RCP_minBias_caloJets.C'

#include "../../headers/functions/divideByBinwidth.h"
#include "../../headers/plotting/plotStyle.h"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TSystem.h"
#include <cstdio>
#include <vector>

const char *repo = "/home/clayton/Analysis/code/bJetRaaAnalysis";
const char *sib  = "/home/clayton/Analysis/code/bJetMuonTaggingAnalysis";  // read-only

const char *clsLabel[5] = {"", "0-10%", "10-30%", "30-50%", "50-80%"};

const int    refClass     = 4;      // 50-80%, the peripheral reference
const double floorPt      = 50.;    // the MinBias scan's own jet cut
// The reference class holds 29 jets in 200-300 GeV and 3 above 300, so a bin
// reaching to 500 would average a steeply falling spectrum across a range where
// the denominator is empty. Stop where the reference still has something to say.
const double ptCeiling    = 300.;
const double targetRelErr = 0.10;   // on the reference class, the bottleneck
const double minWidthFrac = 0.10;   // never narrower than the resolution
const double maxBinWidth  = 200.;

const double Z_lo = 75., Z_hi = 105.;
const double muEff[5] = {1., 0.8627, 0.9069, 0.9856, 0.9778};

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

// merge 5 GeV bins until the reference class carries targetRelErr and the bin
// is at least as wide as the resolution
static std::vector<double> deriveEdges(TH1D *ref)
{
  std::vector<double> edges;
  if(!ref) return edges;
  const double eps = 0.01;
  int b0 = ref->FindBin(floorPt + eps), bN = ref->FindBin(ptCeiling - eps);
  edges.push_back(ref->GetXaxis()->GetBinLowEdge(b0));

  double sum = 0., err2 = 0., startEdge = edges.back();
  for(int b = b0; b <= bN; b++){
    sum  += ref->GetBinContent(b);
    err2 += ref->GetBinError(b)*ref->GetBinError(b);
    double upper = ref->GetXaxis()->GetBinUpEdge(b);
    bool precise  = (sum > 0. && sqrt(err2)/sum <= targetRelErr);
    bool resolved = (upper - startEdge >= minWidthFrac * 0.5 * (startEdge + upper));
    bool tooWide  = (upper - startEdge >= maxBinWidth);
    if((precise && resolved) || tooWide){
      edges.push_back(upper);
      startEdge = upper;
      sum = 0.; err2 = 0.;
    }
  }
  if(edges.size() < 2) return std::vector<double>();
  if(edges.back() < ptCeiling) edges.back() = ptCeiling;   // fold the tail in
  return edges;
}

void jetsPerZ_RCP_minBias_caloJets()
{
  initPlotStyle();
  ok = true;

  TString figDir = Form("%s/figures/JetsPerZ", repo);
  gSystem->mkdir(figDir, kTRUE);
  TString outDir = "./rootFiles/JetsPerZ";
  gSystem->mkdir(outDir, kTRUE);

  TFile *f_MB = openOrFail(Form("%s/rootFiles/scanningOuput/PbPb/PbPb_MinBias_Part1_caloJets_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_2026-9-14.root", repo));
  TFile *f_mu = openOrFail(Form("%s/rootFiles/scanningOutput/PbPb/latest/PbPb_SingleMuon_mu12_pTmu-15to999_tight_mu12TriggerEfficiencyCorrection_jetTrkMaxFilter_WDecayFilter_2026-5-4.root", sib));
  if(!ok) return;

  // the MinBias scan must carry no jet trigger, or the low-pT reach is a lie
  { TH1D *g = nullptr; f_MB->GetObject("h_nEventsNoJetTrigSel", g);
    if(g && g->GetMean() > 0.5){
      printf("ERROR: this MinBias scan has a jet trigger applied; it cannot reach below the turn-on\n");
      return;
    } }

  printf("\n=== R_CP of jets per Z, PbPb MinBias only, calo jets ===\n");
  printf("  reference class: %s   floor: %.0f GeV   no fake-jet subtraction, no unfolding\n",
         clsLabel[refClass], floorPt);

  TH1D *raw[5] = {nullptr,nullptr,nullptr,nullptr,nullptr};
  double NZ[5] = {0.,0.,0.,0.,0.};
  for(int c = 1; c <= 4; c++){
    raw[c] = getOrFail(f_MB, Form("h_inclRecoJetPt_C%d", c), Form("raw_C%d", c));
    TH1D *d = getOrFail(f_mu, Form("h_dimuonMass_C%d", c), Form("dimu_C%d", c));
    if(!ok) return;
    NZ[c] = d->Integral(d->GetXaxis()->FindBin(Z_lo), d->GetXaxis()->FindBin(Z_hi)) / (muEff[c]*muEff[c]);
  }

  printf("\n  %-8s %12s %10s %12s\n", "class", "jets >50 GeV", "N_Z", "jets per Z");
  for(int c = 1; c <= 4; c++){
    double n = raw[c]->Integral(raw[c]->FindBin(floorPt + 0.01), raw[c]->FindBin(ptCeiling - 0.01));
    printf("  %-8s %12.0f %10.1f %12.4f\n", clsLabel[c], n, NZ[c], NZ[c] > 0. ? n/NZ[c] : -1.);
  }

  std::vector<double> edges = deriveEdges(raw[refClass]);
  if(edges.size() < 2){ printf("ERROR: could not derive a binning from the reference class\n"); return; }
  printf("\n  binning (%zu bins, %.0f%% target on %s):", edges.size() - 1, 100.*targetRelErr, clsLabel[refClass]);
  for(size_t e = 0; e < edges.size(); e++) printf(" %.0f", edges[e]);
  printf("\n");

  // jets per Z, per class, on the common binning
  TH1D *spec[5] = {nullptr,nullptr,nullptr,nullptr,nullptr};
  for(int c = 1; c <= 4; c++){
    spec[c] = (TH1D*) raw[c]->Rebin((int)edges.size() - 1, Form("spec_C%d", c), &edges[0]);
    spec[c]->SetDirectory(nullptr);
    divideByBinwidth(spec[c]);
    spec[c]->Scale(1./3.2);        // eta range
    spec[c]->Scale(1./NZ[c]);      // per Z; the luminosity factor cancels in R_CP
  }

  TH1D *rcp[5] = {nullptr,nullptr,nullptr,nullptr,nullptr};
  for(int c = 1; c < refClass; c++){
    rcp[c] = (TH1D*) spec[c]->Clone(Form("RCP_C%d", c));
    rcp[c]->SetDirectory(nullptr);
    rcp[c]->Divide(spec[c], spec[refClass], 1, 1, "");
  }

  printf("\n  R_CP = [jets per Z](class) / [jets per Z](%s)\n", clsLabel[refClass]);
  printf("    %-12s", "pT [GeV]");
  for(int c = 1; c < refClass; c++) printf(" %16s", clsLabel[c]);
  printf("\n");
  for(int b = 1; b <= rcp[1]->GetNbinsX(); b++){
    printf("    %4.0f-%-7.0f", rcp[1]->GetXaxis()->GetBinLowEdge(b), rcp[1]->GetXaxis()->GetBinUpEdge(b));
    for(int c = 1; c < refClass; c++) printf(" %8.3f+-%-6.3f", rcp[c]->GetBinContent(b), rcp[c]->GetBinError(b));
    printf("\n");
  }

  // ------------------------------------------------------------- figures ---
  const char *hexC[5] = {"", okabeHex[5], okabeHex[1], okabeHex[2], okabeHex[0]};
  const int   markC[5] = {0, markFilledCircle, markFilledSquare, markFilledDiamond, markCross};

  // 1. the four jets-per-Z spectra
  {
    TCanvas *cv = new TCanvas("c_spec_RCP", "", 700, 800);
    cv->SetLogy();
    cv->SetLeftMargin(0.17); cv->SetBottomMargin(0.12); cv->SetTopMargin(0.13); cv->SetRightMargin(0.05);
    TLegend *leg = makeLegend(0.58, 0.60, 0.93, 0.80, 0.036);
    for(int c = 1; c <= 4; c++){
      styleH(spec[c], hexC[c], markC[c]);
      if(c == 1){
        spec[c]->SetTitle("");
        spec[c]->GetXaxis()->SetTitle("jet #it{p}_{T} [GeV]");
        spec[c]->GetXaxis()->SetTitleSize(0.045); spec[c]->GetXaxis()->SetLabelSize(0.040);
        spec[c]->GetYaxis()->SetTitle("#frac{1}{#it{N}_{Z}} #frac{d#it{N}_{jet}}{d#it{p}_{T} d#eta}");
        spec[c]->GetYaxis()->SetTitleSize(0.045); spec[c]->GetYaxis()->SetLabelSize(0.040);
        spec[c]->GetYaxis()->SetTitleOffset(1.60);
        spec[c]->Draw("E1");
      }
      else spec[c]->Draw("E1 same");
      leg->AddEntry(spec[c], Form("PbPb %s", clsLabel[c]), "lp");
    }
    leg->Draw();
    TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.036);
    la.DrawLatex(0.17, 0.945, "PbPb 5.02 TeV MinBias, calo jets, anti-#it{k}_{T} #it{R} = 0.4");
    la.DrawLatex(0.17, 0.900, "no trigger stitching, no fake-jet subtraction, no unfolding");
    cv->SaveAs(Form("%s/jetsPerZ_minBias_spectra_caloJets.pdf", figDir.Data()));
    delete cv;
  }

  // 2. R_CP
  {
    TCanvas *cv = new TCanvas("c_RCP", "", 700, 800);
    cv->SetLeftMargin(0.16); cv->SetBottomMargin(0.12); cv->SetTopMargin(0.13); cv->SetRightMargin(0.05);
    TLegend *leg = makeLegend(0.55, 0.62, 0.93, 0.80, 0.036);
    bool first = true;
    for(int c = 1; c < refClass; c++){
      styleH(rcp[c], hexC[c], markC[c]);
      if(first){
        rcp[c]->SetTitle("");
        rcp[c]->SetMinimum(0.); rcp[c]->SetMaximum(1.6);
        rcp[c]->GetXaxis()->SetRangeUser(floorPt, ptCeiling);
        rcp[c]->GetXaxis()->SetTitle("jet #it{p}_{T} [GeV]");
        rcp[c]->GetXaxis()->SetTitleSize(0.045); rcp[c]->GetXaxis()->SetLabelSize(0.040);
        rcp[c]->GetYaxis()->SetTitle(Form("#it{R}_{CP} of jets per #it{Z}   (ref. %s)", clsLabel[refClass]));
        rcp[c]->GetYaxis()->SetTitleSize(0.042); rcp[c]->GetYaxis()->SetLabelSize(0.040);
        rcp[c]->GetYaxis()->SetTitleOffset(1.50);
        rcp[c]->Draw("E1");
        first = false;
      }
      else rcp[c]->Draw("E1 same");
      leg->AddEntry(rcp[c], Form("PbPb %s", clsLabel[c]), "lp");
    }
    TLine *one = new TLine(floorPt, 1., ptCeiling, 1.);
    one->SetLineStyle(7); one->Draw();
    for(int c = 1; c < refClass; c++) rcp[c]->Draw("E1 same");
    leg->Draw();
    TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.036);
    la.DrawLatex(0.16, 0.945, "PbPb 5.02 TeV MinBias, calo jets");
    la.SetTextSize(0.032);
    // kept short deliberately: the longer wordings ran off the right edge
    la.DrawLatex(0.16, 0.902, "no unfolding, no fake jet subtraction");
    cv->SaveAs(Form("%s/jetsPerZ_RCP_minBias_caloJets.pdf", figDir.Data()));
    delete cv;
  }

  TString outPath = Form("%s/jetsPerZ_RCP_minBias_caloJets.root", outDir.Data());
  TFile *wf = TFile::Open(outPath, "recreate");
  for(int c = 1; c <= 4; c++) spec[c]->Write(Form("jetsPerZ_C%d", c));
  for(int c = 1; c < refClass; c++) rcp[c]->Write(Form("RCP_C%d", c));
  wf->Close();

  printf("\n  figures: %s/jetsPerZ_{minBias_spectra,RCP_minBias}_caloJets.pdf\n", figDir.Data());
  printf("  results: %s\n", outPath.Data());
}
