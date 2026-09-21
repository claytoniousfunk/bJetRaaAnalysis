// Jets per Z, PbPb/pp, one canvas per centrality class, three curves each:
//   1. calo jets, this analysis   -- stitchSpectra_caloJets.C (Jet80+Jet100,
//                                    per-class binning and floors)
//   2. PF jets, this analysis     -- calculateJetsPerZ_caloJets.C(false):
//                                    same no-unfolding / no-fake-subtraction
//                                    convention as the calo curve, and bin for
//                                    bin identical to calculateRAA.C's stored
//                                    PF no-unfold output above 130 GeV
//   3. ATLAS R_AA                 -- PLB 790 (2019) 108, HEPData ins1673184
//
// ATLAS is not the same observable. It is R_AA normalized with T_AA and fully
// unfolded; ours normalize by the Z yield and are not unfolded. The comparison
// is of shape and magnitude, not a like-for-like test.
//
// ATLAS centralities do not match ours. The four tables used here (as in
// calculateRAA.C) are 0-10, 10-20, 30-40 and 50-60%, i.e. narrower slices
// inside our 10-30, 30-50 and 50-80% classes. The legend says which.
//
// ATLAS errors (from the HEPData axis titles): e1 = asymmetric systematic,
// e2 = statistical, e3 = T_AA and e4 = luminosity, both global. Drawn here: stat
// as bars, e1 as boxes. The global T_AA and luminosity terms are not drawn.
// Note calculateRAA.C draws e1plus as the error bar, i.e. systematic, not stat.
//
// Usage, from src/newFractionCalculation/ (run stitchSpectra_caloJets.C and
// calculateJetsPerZ_caloJets.C(false) first):
//   root -l -b -q 'plotJetsPerZ_caloVsPFvsATLAS.C'

#include "../../headers/plotting/plotStyle.h"
#include "TFile.h"
#include "TH1D.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TSystem.h"
#include <cstdio>

const char *repo = "/home/clayton/Analysis/code/bJetRaaAnalysis";
const char *hepData = "/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/HEPData";  // read-only

const char *caloPath = "./rootFiles/JetsPerZ/jetsPerZ_caloJets_stitched.root";
const char *pfPath   = "./rootFiles/JetsPerZ/histograms_JetsPerZ_PFJets_noUnfold.root";

const char *clsLabel[5]   = {"", "0-10%", "10-30%", "30-50%", "50-80%"};
const char *clsTag[5]     = {"", "0to10pct", "10to30pct", "30to50pct", "50to80pct"};
const int   atlasTable[5] = {0, 19, 20, 22, 24};
const char *atlasCent[5]  = {"", "0-10%", "10-20%", "30-40%", "50-60%"};

const double xHi = 500.;

struct AtlasCurve {
  TH1 *stat = nullptr;
  TGraphAsymmErrors *sys = nullptr;
};

static AtlasCurve readATLAS(int table, int c)
{
  AtlasCurve a;
  TFile *f = TFile::Open(Form("%s/HEPData-ins1673184-v1-Table_%d.root", hepData, table));
  if(!f || f->IsZombie()){ printf("ERROR: cannot open ATLAS table %d\n", table); return a; }
  TDirectory *d = (TDirectory*) f->Get(Form("Table %d", table));
  if(!d){ printf("ERROR: no 'Table %d' directory\n", table); return a; }
  TH1 *y  = (TH1*) d->Get("Hist1D_y1");
  TH1 *st = (TH1*) d->Get("Hist1D_y1_e2");        // stat
  TH1 *sp = (TH1*) d->Get("Hist1D_y1_e1plus");    // sys +
  TH1 *sm = (TH1*) d->Get("Hist1D_y1_e1minus");   // sys - (stored negative)
  if(!y || !st || !sp || !sm){ printf("ERROR: ATLAS table %d is missing a histogram\n", table); return a; }

  a.stat = (TH1*) y->Clone(Form("atlas_stat_C%d", c));
  a.stat->SetDirectory(nullptr);
  a.sys = new TGraphAsymmErrors();
  a.sys->SetName(Form("atlas_sys_C%d", c));
  for(int b = 1; b <= y->GetNbinsX(); b++){
    a.stat->SetBinError(b, st->GetBinContent(b));
    double x = y->GetBinCenter(b), hw = 0.5 * y->GetBinWidth(b);
    int n = a.sys->GetN();
    a.sys->SetPoint(n, x, y->GetBinContent(b));
    a.sys->SetPointError(n, hw, hw, fabs(sm->GetBinContent(b)), fabs(sp->GetBinContent(b)));
  }
  return a;
}

void plotJetsPerZ_caloVsPFvsATLAS()
{
  initPlotStyle();

  TFile *fCalo = TFile::Open(caloPath);
  TFile *fPF   = TFile::Open(pfPath);
  if(!fCalo || fCalo->IsZombie()){ printf("ERROR: %s missing -- run stitchSpectra_caloJets.C\n", caloPath); return; }
  if(!fPF   || fPF->IsZombie())  { printf("ERROR: %s missing -- run calculateJetsPerZ_caloJets.C(false)\n", pfPath); return; }

  TString figDir = Form("%s/figures/JetsPerZ", repo);
  gSystem->mkdir(figDir, kTRUE);

  for(int c = 1; c <= 4; c++){
    TH1D *calo = nullptr, *pf = nullptr;
    fCalo->GetObject(Form("ratio_C%d", c), calo);
    fPF  ->GetObject(Form("r_C%d_fine", c), pf);
    AtlasCurve atlas = readATLAS(atlasTable[c], c);
    if(!calo || !pf || !atlas.stat){ printf("ERROR: missing input for %s\n", clsLabel[c]); return; }
    calo = (TH1D*) calo->Clone(Form("calo_C%d", c)); calo->SetDirectory(nullptr);
    pf   = (TH1D*) pf  ->Clone(Form("pf_C%d", c));   pf->SetDirectory(nullptr);

    // frame from the lowest pT any of the three reaches, capped at 50
    double xLo = TMath::Min(calo->GetXaxis()->GetBinLowEdge(1), pf->GetXaxis()->GetBinLowEdge(1));
    xLo = TMath::Min(xLo, atlas.stat->GetXaxis()->GetBinLowEdge(1));
    xLo = TMath::Max(xLo, 50.);

    TCanvas *cv = new TCanvas(Form("c_cmp_C%d", c), "", 700, 800);
    cv->SetLeftMargin(0.16); cv->SetBottomMargin(0.12);
    cv->SetTopMargin(0.17);  cv->SetRightMargin(0.05);   // three caption lines

    TH1D *frame = new TH1D(Form("frame_C%d", c), "", 1, xLo, xHi);
    frame->SetStats(0);
    frame->SetMinimum(0.); frame->SetMaximum(1.4);
    frame->GetXaxis()->SetTitle("jet #it{p}_{T} [GeV]");
    frame->GetXaxis()->SetTitleSize(0.045); frame->GetXaxis()->SetLabelSize(0.040);
    frame->GetYaxis()->SetTitle("jets per #it{Z} PbPb / pp   or   #it{R}_{AA}");
    frame->GetYaxis()->SetTitleSize(0.042); frame->GetYaxis()->SetLabelSize(0.040);
    frame->GetYaxis()->SetTitleOffset(1.50);
    frame->Draw("axis");

    // ATLAS systematic boxes underneath everything
    atlas.sys->SetFillColor(TColor::GetColor("#BBBBBB"));
    atlas.sys->SetLineColor(TColor::GetColor("#BBBBBB"));
    atlas.sys->Draw("2 same");

    TLine *one = new TLine(xLo, 1., xHi, 1.);
    one->SetLineStyle(7);
    one->Draw();

    styleH(atlas.stat, hexData,      markFilledDiamond, 1.2);
    styleH(pf,         hexMC,        markOpenCircle);
    styleH(calo,       hexCorrected, markFilledSquare);
    atlas.stat->Draw("E1 X0 same");
    pf->Draw("E1 same");
    calo->Draw("E1 same");

    // upper right, above the unity line: every curve sits below ~0.97.
    // Top edge stays under the frame, which the 0.17 top margin puts at 0.83.
    TLegend *leg = makeLegend(0.42, 0.63, 0.93, 0.81, 0.034);
    leg->AddEntry(calo, "calo jets, this analysis", "lp");
    leg->AddEntry(pf,   "PF jets, this analysis",   "lp");
    leg->AddEntry(atlas.stat, Form("ATLAS #it{R}_{AA}, %s", atlasCent[c]), "p");
    leg->Draw();

    TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.036);
    la.DrawLatex(0.16, 0.955, Form("PbPb %s, 5.02 TeV, anti-#it{k}_{T} #it{R} = 0.4", clsLabel[c]));
    la.SetTextSize(0.030);
    la.DrawLatex(0.16, 0.915, "this analysis: no unfolding, no fake jet subtraction");
    // Without the subtraction the PF ratio below ~100 GeV is fake dominated and
    // leaves the frame in the central classes (0-10% reaches 3.8 at 50-55). Say
    // so rather than let the missing points read as a shorter curve.
    if(c <= 2){
      la.SetTextSize(0.027);
      la.DrawLatex(0.16, 0.879, "PF below #approx100 GeV is fake dominated and runs off scale");
    }

    cv->SaveAs(Form("%s/jetsPerZ_caloVsPFvsATLAS_%s.pdf", figDir.Data(), clsTag[c]));
    delete cv;

    printf("  %-7s calo from %.0f GeV, PF from %.0f, ATLAS %s from %.0f\n", clsLabel[c],
           calo->GetXaxis()->GetBinLowEdge(1), pf->GetXaxis()->GetBinLowEdge(1),
           atlasCent[c], atlas.stat->GetXaxis()->GetBinLowEdge(1));
  }
  printf("  figures: %s/jetsPerZ_caloVsPFvsATLAS_<class>.pdf\n", figDir.Data());
}
