// Gen-level effect of the JEU on the PYTHIA+HYDJET calo-jet spectrum, per
// centrality class: unfold the nominal and JEU-shifted reco spectra with ONE
// response (the nominal one, same class) and take shifted / nominal at gen level.
// The PbPb counterpart of unfoldPYTHIAJetPt_JEUShift.C; read that header for what
// the ratio measures.
//
// Reco spectra  h_inclRecoJetPt_flavor_{nominal,JEUShiftUp,JEUShiftDown}_C<i>
//               from PYTHIAHYDJET_scan.C (booked 2026-10-02; scans before that
//               do not have them), flavors summed. i = 1..4 -> 0-10, 10-30,
//               30-50, 50-80%.
// Response      PYTHIAHYDJET_scan_response, calo + manual JEC, 2026-09-25 even +
//               odd halves: h_matchedRecoJetPt_genJetPt_allJets_C<i>, misses
//               from h_unmatchedGenJetPt_C<i> in the truth. Reco floor 50 GeV,
//               Bayes, 2 iterations; N = 1..4 printed.
// JEU file      Autumn18_HI_V8_MC_Uncertainty_AK4Calo (set in PYTHIAHYDJET_scan.C).
//
// Caveats beyond the pp ones. FAKES ARE NOT INCLUDED in the response, so the
// background-fluctuation jets that populate the low-pT reco spectrum in central
// events are unfolded as if genuine: the first wide bins of 0-10% and 10-30%
// are not reliable. The spectrum scan and the response scan are different runs,
// so only the shifted/nominal ratio is meaningful. Ratio errors are the
// shifted-spectrum statistics only (same jets as nominal, covariance not
// propagated).
//
// Output  figures/unfoldTest/unfold_PYTHIAHYDJET_caloJets_JEUShiftRatio.pdf
//         (one panel per class: gen ratio points, reco ratio dashed)
//
// Usage: root -l -b -q -e 'gSystem->Load("/home/clayton/Programs/RooUnfold/build/libRooUnfold.so");
//                          gInterpreter->AddIncludePath("/home/clayton/Programs/RooUnfold/build");'
//                     'unfoldPYTHIAHYDJETJetPt_JEUShift.C+("<scan.root>")'
// Run from: src/unfoldTest/
#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#endif
#include "../../headers/plotting/plotStyle.h"
#include "../../headers/plotting/ratioPanel.h"
#include "../../headers/functions/divideByBinwidth.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include <cstdio>

namespace {

const char *repo   = "/home/clayton/Analysis/code/bJetRaaAnalysis";
const char *mcDir  = "rootFiles/scanningOuput/PYTHIAHYDJET/";
const char *respFmt = "PYTHIAHYDJET_response_DiJet_caloJets_manualJEC_%sEvents_pThat-15_mu12_pTmu-15_tight_vzReweight_hiBinReweight_hiBinShift-10_jetTrkMaxFilter_doPThatCorrelationFilterTight_2026-9-25.root";
const char *outPdf = "/home/clayton/Analysis/code/bJetRaaAnalysis/figures/unfoldTest/unfold_PYTHIAHYDJET_caloJets_JEUShiftRatio.pdf";

const double recoFloor = 50.;
const int    nIterNominal = 2, nIterMax = 4;
const double xLo = 50., xHi = 400.;

// wide bins for the ratios; every edge on the 5 GeV grid
const double wideEdges[] = {50, 60, 70, 80, 100, 120, 150, 200, 250, 300, 400};
const int    nWide = sizeof(wideEdges)/sizeof(double) - 1;

TH1* getH(TFile *f, const char *name)
{
  TH1 *h = nullptr; f->GetObject(name, h);
  if(!h) printf("ERROR: %s missing in %s\n", name, f->GetName());
  return h;
}

TH1D* wide(TH1D *h, const char *name)
{
  TH1D *r = (TH1D*) h->Rebin(nWide, name, wideEdges);   // new histogram
  r->SetDirectory(nullptr);
  return r;
}

TH1D* flavorSum(TFile *f, const char *name, const char *tag)
{
  TH2D *H = (TH2D*) getH(f, name);
  if(!H) return nullptr;
  TH1D *h = H->ProjectionX(tag, 1, H->GetNbinsY());
  h->SetDirectory(nullptr);
  for(int b = 0; b <= h->GetNbinsX() + 1; b++)
    if(h->GetXaxis()->GetBinCenter(b) < recoFloor){ h->SetBinContent(b, 0.); h->SetBinError(b, 0.); }
  return h;
}

// shifted / nominal in wide bins; error = shifted-spectrum relative error of the
// numerator only (see header)
TH1D* wideRatio(TH1D *num, TH1D *den, const char *name)
{
  TH1D *n = wide(num, Form("%s_n", name)), *d = wide(den, Form("%s_d", name));
  TH1D *r = makeRatio(n, d, name);
  delete n; delete d;
  return r;
}

const char *centLabel[5] = {"", "0-10%", "10-30%", "30-50%", "50-80%"};

} // namespace

void unfoldPYTHIAHYDJETJetPt_JEUShift(const char *scanFile)
{
  initPlotStyle();
  TFile *fe = TFile::Open(Form("%s/%s%s", repo, mcDir, Form(respFmt, "even")));
  TFile *fo = TFile::Open(Form("%s/%s%s", repo, mcDir, Form(respFmt, "odd")));
  TFile *fs = TFile::Open(scanFile);
  if(!fe || fe->IsZombie() || !fo || fo->IsZombie() || !fs || fs->IsZombie()){ printf("ERROR: cannot open inputs\n"); return; }

  const char *var[3] = {"nominal", "JEUShiftUp", "JEUShiftDown"};
  TCanvas *c = new TCanvas("c", "", 1100, 900);
  c->Divide(2, 2, 0.001, 0.001);

  for(int ci = 1; ci <= 4; ci++){
    TH2D *resp2D = (TH2D*) getH(fe, Form("h_matchedRecoJetPt_genJetPt_allJets_C%d", ci));
    TH2D *resp2Do= (TH2D*) getH(fo, Form("h_matchedRecoJetPt_genJetPt_allJets_C%d", ci));
    TH1D *miss   = (TH1D*) getH(fe, Form("h_unmatchedGenJetPt_C%d", ci));
    TH1D *missO  = (TH1D*) getH(fo, Form("h_unmatchedGenJetPt_C%d", ci));
    if(!resp2D || !resp2Do || !miss || !missO) return;
    resp2D = (TH2D*) resp2D->Clone(Form("resp2D_%d", ci)); resp2D->SetDirectory(nullptr); resp2D->Add(resp2Do);
    miss   = (TH1D*) miss->Clone(Form("miss_%d", ci));     miss->SetDirectory(nullptr);   miss->Add(missO);

    TH1D *truth = resp2D->ProjectionY(Form("truth_%d", ci), 1, resp2D->GetNbinsX());
    truth->SetDirectory(nullptr); truth->Add(miss);
    for(int ix = 0; ix <= resp2D->GetNbinsX(); ix++)
      if(resp2D->GetXaxis()->GetBinCenter(ix) < recoFloor)
        for(int iy = 0; iy <= resp2D->GetNbinsY() + 1; iy++){ resp2D->SetBinContent(ix, iy, 0.); resp2D->SetBinError(ix, iy, 0.); }
    TH1D *measMC = resp2D->ProjectionX(Form("measMC_%d", ci), 1, resp2D->GetNbinsY());
    measMC->SetDirectory(nullptr);
    RooUnfoldResponse resp(measMC, truth, resp2D, Form("resp_C%d", ci), "PYTHIAHYDJET calo response");

    TH1D *reco[3], *gen[3][nIterMax + 1];
    for(int v = 0; v < 3; v++){
      reco[v] = flavorSum(fs, Form("h_inclRecoJetPt_flavor_%s_C%d", var[v], ci), Form("reco_%s_C%d", var[v], ci));
      if(!reco[v]){ printf("ERROR: the scan has no JEU histograms\n"); return; }
      RooUnfoldBayes unfold(&resp, reco[v], 1);
      unfold.SetVerbose(0);
      for(int n = 1; n <= nIterMax; n++){
        unfold.SetIterations(n);
        gen[v][n] = (TH1D*) unfold.Hunfold();
        gen[v][n]->SetName(Form("gen_%s_C%d_N%d", var[v], ci, n)); gen[v][n]->SetDirectory(nullptr);
      }
    }

    printf("\n%s: shifted/nominal in wide bins (reco | gen N=1..4)\n pT-bin    up: reco   N1     N2     N3     N4    | down: reco   N1     N2     N3     N4\n", centLabel[ci]);
    TH1D *rReco[2], *rGen[2][nIterMax + 1];
    for(int s = 0; s < 2; s++){
      rReco[s] = wideRatio(reco[s+1], reco[0], Form("rReco_%d_C%d", s, ci));
      for(int n = 1; n <= nIterMax; n++) rGen[s][n] = wideRatio(gen[s+1][n], gen[0][n], Form("rGen_%d_%d_C%d", s, n, ci));
    }
    for(int b = 1; b <= nWide; b++){
      printf("%4.0f-%-4.0f", wideEdges[b-1], wideEdges[b]);
      for(int s = 0; s < 2; s++){
        printf("   %.3f", rReco[s]->GetBinContent(b));
        for(int n = 1; n <= nIterMax; n++) printf(" %.3f", rGen[s][n]->GetBinContent(b));
        printf("  |");
      }
      printf("\n");
    }

    c->cd(ci);
    gPad->SetLeftMargin(0.15); gPad->SetBottomMargin(0.14); gPad->SetTopMargin(0.06); gPad->SetRightMargin(0.04);
    TH1D *fr = rGen[0][nIterNominal];
    fr->SetTitle("");
    fr->GetXaxis()->SetRangeUser(xLo, xHi);
    fr->GetYaxis()->SetRangeUser(0.6, 1.5);
    fr->GetXaxis()->SetTitle("Jet p_{T} (GeV)"); fr->GetYaxis()->SetTitle("Shifted / nominal");
    fr->GetXaxis()->SetTitleSize(0.055); fr->GetYaxis()->SetTitleSize(0.055);
    fr->GetXaxis()->SetLabelSize(0.05);  fr->GetYaxis()->SetLabelSize(0.05);
    fr->GetYaxis()->SetTitleOffset(1.2);
    styleH(fr, okabeHex[6], markFilledSquare, 1.0);
    styleH(rGen[1][nIterNominal], okabeHex[5], markFilledDiamond, 1.2);
    fr->Draw("E1"); rGen[1][nIterNominal]->Draw("E1 SAME");
    for(int s = 0; s < 2; s++){
      rReco[s]->SetLineColor(TColor::GetColor(okabeHex[s ? 5 : 6])); rReco[s]->SetLineStyle(2);
      rReco[s]->SetLineWidth(2); rReco[s]->SetMarkerSize(0);
      rReco[s]->Draw("HIST SAME");
    }
    unityLine(xLo, xHi)->Draw();
    TLegend *l2 = makeLegend(0.18, 0.70, 0.60, 0.88, 0.05);
    l2->AddEntry(fr, "JEU up, gen", "p");
    l2->AddEntry(rGen[1][nIterNominal], "JEU down, gen", "p");
    l2->AddEntry(rReco[0], "reco", "l");
    l2->Draw();
    TLatex t; t.SetNDC(); t.SetTextSize(0.05);
    t.DrawLatex(0.50, 0.86, Form("PYTHIA+HYDJET %s", centLabel[ci]));
  }

  gSystem->mkdir(gSystem->DirName(outPdf), kTRUE);
  c->SaveAs(outPdf);
}
