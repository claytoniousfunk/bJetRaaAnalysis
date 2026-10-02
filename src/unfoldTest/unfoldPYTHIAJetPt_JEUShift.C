// Gen-level effect of the JEU on the PYTHIA calo-jet spectrum: unfold the
// nominal and the JEU-shifted reco spectra with ONE response (the nominal one)
// and take shifted / nominal at gen level.
//
// Reco spectra  h_inclRecoJetPt_flavor_{nominal,JEUShiftUp,JEUShiftDown}
//               from the 2026-10-02 PYTHIA_scan (all flavors summed).
// Response      PYTHIA_scan_response, calo + manual JEC, 2026-10-01 even + odd
//               halves added: h_matchedRecoJetPt_genJetPt_allJets, misses from
//               h_unmatchedGenJetPt in the truth (as unfoldppData_caloJets.C).
//               Reco floor 50 GeV, Bayesian, 2 iterations nominal; N = 1..4 are
//               printed so the iteration dependence of the gen-level ratio can
//               be read off.
//
// What this measures. A JEU shift moves jets in reco pT, so the reco ratio is
// ~ +/-12% on a pT^-5 spectrum. Unfolding with the nominal response asks what
// gen spectrum, folded with the nominal detector, would reproduce each shifted
// reco spectrum; the gen-level ratio is the uncertainty on the corrected
// spectrum. It is ~ the reco ratio when the response is narrow.
//
// Caveats. The two inputs are different scans (the response scan applies the
// pThat-correlation filter; the jet-spectrum scan does not), so the unfolded
// spectra themselves are not a physics result -- only the shifted/nominal
// ratio is meaningful. Errors drawn on the ratios are the measured-spectrum
// statistics of the shifted histogram; the shifted and nominal spectra are
// built from the same jets, so the two are not independent and the unfolded
// covariance is not propagated.
//
// Usage: root -l -b -q -e 'gSystem->Load("/home/clayton/Programs/RooUnfold/build/libRooUnfold.so");
//                          gInterpreter->AddIncludePath("/home/clayton/Programs/RooUnfold/build");'
//                     'unfoldPYTHIAJetPt_JEUShift.C+'
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
const char *mcDir  = "rootFiles/scanningOuput/PYTHIA/";
const char *scanF  = "PYTHIA_DiJet_caloJets_manualJEC_pThat-15_mu12_pTmu-15to999_tight_vzReweight_jetTrkMaxFilter_removeHYDJETjet0p35CutOnGen_WDecayFilter_2026-10-2.root";
const char *respE  = "PYTHIA_DiJet_response_caloJets_PFflavor_manualJEC_evenEvents_pThat-15_mu12_pTmu-15_tight_jetTrkMaxFilter_doPThatCorrelationFilterTight_2026-10-1.root";
const char *respO  = "PYTHIA_DiJet_response_caloJets_PFflavor_manualJEC_oddEvents_pThat-15_mu12_pTmu-15_tight_jetTrkMaxFilter_doPThatCorrelationFilterTight_2026-10-1.root";
const char *outPdf = "/home/clayton/Analysis/code/bJetRaaAnalysis/figures/unfoldTest/unfold_PYTHIA_caloJets_JEUShiftRatio.pdf";

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

} // namespace

void unfoldPYTHIAJetPt_JEUShift()
{
  initPlotStyle();
  TFile *fe = TFile::Open(Form("%s/%s%s", repo, mcDir, respE));
  TFile *fo = TFile::Open(Form("%s/%s%s", repo, mcDir, respO));
  TFile *fs = TFile::Open(Form("%s/%s%s", repo, mcDir, scanF));
  if(!fe || fe->IsZombie() || !fo || fo->IsZombie() || !fs || fs->IsZombie()){ printf("ERROR: cannot open inputs\n"); return; }

  TH2D *resp2D = (TH2D*) getH(fe, "h_matchedRecoJetPt_genJetPt_allJets");
  TH2D *resp2Do= (TH2D*) getH(fo, "h_matchedRecoJetPt_genJetPt_allJets");
  TH1D *miss   = (TH1D*) getH(fe, "h_unmatchedGenJetPt");
  TH1D *missO  = (TH1D*) getH(fo, "h_unmatchedGenJetPt");
  if(!resp2D || !resp2Do || !miss || !missO) return;
  resp2D = (TH2D*) resp2D->Clone("resp2D"); resp2D->SetDirectory(nullptr); resp2D->Add(resp2Do);
  miss   = (TH1D*) miss->Clone("miss");     miss->SetDirectory(nullptr);   miss->Add(missO);

  // truth: every gen jet (matched + missed), before the reco floor
  TH1D *truth = resp2D->ProjectionY("truth", 1, resp2D->GetNbinsX());
  truth->SetDirectory(nullptr); truth->Add(miss);

  for(int ix = 0; ix <= resp2D->GetNbinsX(); ix++)
    if(resp2D->GetXaxis()->GetBinCenter(ix) < recoFloor)
      for(int iy = 0; iy <= resp2D->GetNbinsY() + 1; iy++){ resp2D->SetBinContent(ix, iy, 0.); resp2D->SetBinError(ix, iy, 0.); }
  TH1D *measMC = resp2D->ProjectionX("measMC", 1, resp2D->GetNbinsY());
  measMC->SetDirectory(nullptr);

  RooUnfoldResponse resp(measMC, truth, resp2D, "resp_pp_calo", "PYTHIA calo response");

  const char *var[3]   = {"nominal", "JEUShiftUp", "JEUShiftDown"};
  TH1D *reco[3], *gen[3][nIterMax + 1];
  for(int v = 0; v < 3; v++){
    reco[v] = flavorSum(fs, Form("h_inclRecoJetPt_flavor_%s", var[v]), Form("reco_%s", var[v]));
    if(!reco[v]){ printf("ERROR: the scan has no JEU histograms\n"); return; }
    RooUnfoldBayes unfold(&resp, reco[v], 1);
    unfold.SetVerbose(0);
    for(int n = 1; n <= nIterMax; n++){
      unfold.SetIterations(n);
      gen[v][n] = (TH1D*) unfold.Hunfold();
      gen[v][n]->SetName(Form("gen_%s_N%d", var[v], n)); gen[v][n]->SetDirectory(nullptr);
    }
  }

  // iteration dependence, printed
  printf("\nshifted/nominal in wide bins (reco | gen N=1..4)\n pT-bin    up: reco   N1     N2     N3     N4    | down: reco   N1     N2     N3     N4\n");
  TH1D *rReco[2], *rGen[2][nIterMax + 1];
  for(int s = 0; s < 2; s++){
    rReco[s] = wideRatio(reco[s+1], reco[0], Form("rReco_%d", s));
    for(int n = 1; n <= nIterMax; n++) rGen[s][n] = wideRatio(gen[s+1][n], gen[0][n], Form("rGen_%d_%d", s, n));
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

  // ---- figure ---------------------------------------------------------------
  TH1D *gN = (TH1D*) gen[0][nIterNominal]->Clone("gN"), *gU = (TH1D*) gen[1][nIterNominal]->Clone("gU"),
       *gD = (TH1D*) gen[2][nIterNominal]->Clone("gD");
  divideByBinwidth(gN); divideByBinwidth(gU); divideByBinwidth(gD);
  styleH(gN, okabeHex[0], markFilledCircle, 0.8);
  styleH(gU, okabeHex[6], markFilledSquare, 0.8);
  styleH(gD, okabeHex[5], markFilledDiamond, 1.0);

  TCanvas *c = new TCanvas("c", "", 700, 800);
  TPad *top, *bot; splitPads(top, bot);
  top->cd(); top->SetLogy();
  gN->SetTitle("");
  gN->GetXaxis()->SetRangeUser(xLo, xHi);
  gN->GetYaxis()->SetTitle("dN/dp_{T}^{gen} (arb.)");
  gN->GetYaxis()->SetTitleSize(0.055); gN->GetYaxis()->SetLabelSize(0.050);
  gN->GetXaxis()->SetLabelSize(0);
  gN->Draw("E1"); gU->Draw("E1 SAME"); gD->Draw("E1 SAME");
  TLegend *leg = makeLegend(0.55, 0.66, 0.90, 0.88, 0.048);
  leg->AddEntry(gN, "Unfolded nominal", "lp");
  leg->AddEntry(gU, "Unfolded JEU up", "lp");
  leg->AddEntry(gD, "Unfolded JEU down", "lp");
  leg->Draw();
  TLatex t; t.SetNDC(); t.SetTextSize(0.05);
  t.DrawLatex(0.20, 0.18, "PYTHIA8 calo jets, inclusive");
  t.SetTextSize(0.042);
  t.DrawLatex(0.20, 0.11, Form("Bayes, %d iterations", nIterNominal));

  bot->cd();
  const int colUp = TColor::GetColor(okabeHex[6]), colDn = TColor::GetColor(okabeHex[5]);
  TH1D *fr = rGen[0][nIterNominal];
  styleRatioAxes(fr, "Jet p_{T} (GeV)", "Shifted / nominal");
  fr->GetXaxis()->SetRangeUser(xLo, xHi);
  fr->GetYaxis()->SetRangeUser(0.6, 1.5);
  styleH(fr, okabeHex[6], markFilledSquare, 1.0);
  styleH(rGen[1][nIterNominal], okabeHex[5], markFilledDiamond, 1.2);
  fr->Draw("E1"); rGen[1][nIterNominal]->Draw("E1 SAME");
  for(int s = 0; s < 2; s++){   // reco-level ratio, for comparison
    rReco[s]->SetLineColor(s ? colDn : colUp); rReco[s]->SetLineStyle(2); rReco[s]->SetLineWidth(2);
    rReco[s]->SetMarkerSize(0);
    rReco[s]->Draw("HIST SAME");
  }
  unityLine(xLo, xHi)->Draw();
  TLegend *l2 = makeLegend(0.20, 0.72, 0.60, 0.93, 0.075);
  l2->SetNColumns(2);
  l2->AddEntry(fr, "gen", "lp");
  l2->AddEntry(rReco[0], "reco", "l");
  l2->Draw();

  gSystem->mkdir(gSystem->DirName(outPdf), kTRUE);
  savePdfTight(c, outPdf);
}
