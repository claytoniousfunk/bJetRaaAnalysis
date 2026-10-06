// R_CP of jets per Z, PbPb only, from MinBias alone -- PF jets.
//
// PF counterpart of jetsPerZ_RCP_minBias_caloJets.C; the normalisation argument
// is identical and documented there. In short: jets per Z (c) = N_jets(c)/N_Z(c)
// with raw MinBias jet counts and the efficiency-corrected Z yield of the same
// class; the MinBias/SingleMuon luminosity ratio is class independent and
// cancels in R_CP(c) = [jets per Z](c) / [jets per Z](50-80%).
//
// Two variants from the same 2026-10-5 MinBias Part1 PF scan:
//   all       h_inclRecoJetPt              every PF jet
//   matched   h_inclRecoJetPt_caloMatched  PF jet with an akPu4Calo jet within
//                                          dR < 0.2 -- a jet ID against
//                                          combinatorial PF jets
// Both share one binning, derived from the reference class of the "all"
// variant, so the two can be overlaid bin by bin.
//
// The scan is ultraFine (5%) and is summed to coarse classes; the Z yields come
// from the same SingleMuon file the calo macro uses, so a PF-vs-calo comparison
// differs only in the jets. Z counts carry no jet information, so this does not
// pair calo and PF jets.
//
// UNFOLDING (doUnfolding = true) follows calculateJetsPerZ.cc: D'Agostini
// Bayes, nIter iterations, per-class PYTHIA+HYDJET response
// h_matchedRecoJetPt_genJetPt_allJets_C1-4 (4CentBins file, so C1-4 are the
// coarse classes), reco below recoFloor removed from data and response, truth
// from the UNtruncated response so jets reconstructed below the floor count as
// inefficiency. Unfolded on the scan's 5 GeV bins, then rebinned to the common
// edges. The response contains reco-gen matched jets only -- no fakes -- so
// any combinatorial jets left in the data are unfolded as if they were genuine.
// There is no calo-matched response yet (the PF response scan configured in
// 5c0bde49 has not been run), so both variants use the all-jets response.
//
// CAVEATS: no fake-jet subtraction in "all" (calo matching is the only fake
// suppression here, and its efficiency for genuine jets is not yet measured),
// MinBias Part1 only (a partial file list).
//
// Usage, from src/newFractionCalculation/:
//   root -l -b -q -e 'gSystem->Load("/home/clayton/Programs/RooUnfold/build/libRooUnfold.so");
//                     gInterpreter->AddIncludePath("/home/clayton/Programs/RooUnfold/build");' \
//        'jetsPerZ_RCP_minBias_PFJets.C(true)'
//   argument false: no unfolding (RooUnfold still has to be loadable)

#include "../../headers/functions/divideByBinwidth.h"
#include "../../headers/plotting/plotStyle.h"
#include "../../headers/plotting/coarseCent.h"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TSystem.h"
#include <cstdio>
#include <vector>
#include <algorithm>
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"

const char *repo = "/home/clayton/Analysis/code/bJetRaaAnalysis";
const char *sib  = "/home/clayton/Analysis/code/bJetMuonTaggingAnalysis";  // read-only

const int    refClass     = 3;      // coarse index of 50-80%, the reference
const double floorPt      = 50.;    // kept at the calo macro's floor for comparison
const double ptCeiling    = 300.;
const double targetRelErr = 0.10;
const double minWidthFrac = 0.10;
const double maxBinWidth  = 200.;

const double Z_lo = 75., Z_hi = 105.;
const double muEff[NCoarse] = {0.8627, 0.9069, 0.9856, 0.9778};

const char *responsePath = "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PYTHIAHYDJET/PYTHIAHYDJET_response_DiJet_manualJEC_pThat-15_mu12_pTmu-15_tight_vzReweight_hiBinReweight_hiBinShift-10_jetTrkMaxFilter_doPThatCorrelationFilterTight_2026-9-21.root";
int          nIter     = 1;         // the calculateJetsPerZ.cc nominal; 2nd argument overrides
const double recoFloor = floorPt;

const int NVar = 2;
const char *varBase[NVar]  = {"h_inclRecoJetPt", "h_inclRecoJetPt_caloMatched"};
const char *varTag[NVar]   = {"all", "caloMatched"};
const char *varLabel[NVar] = {"all PF jets", "calo-matched (#DeltaR < 0.2)"};

static std::vector<double> deriveEdges(TH1D *ref)
{
  std::vector<double> edges;
  int b0 = ref->FindBin(floorPt + 1e-6), bN = ref->FindBin(ptCeiling - 1e-6);
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
      startEdge = upper; sum = 0.; err2 = 0.;
    }
  }
  if(edges.size() < 2) return std::vector<double>();
  if(edges.back() < ptCeiling) edges.back() = ptCeiling;
  return edges;
}

static void zeroBelowFloor(TH1D *h)
{
  for(int b = 0; b <= h->GetNbinsX()+1; b++)
    if(h->GetXaxis()->GetBinUpEdge(b) <= recoFloor + 1e-6){ h->SetBinContent(b, 0.); h->SetBinError(b, 0.); }
}

// Unfold one fine (5 GeV) reco spectrum with one class's response.
static TH1D* unfold(TH1D *data, TH2D *resp, int iter, const char *name)
{
  TH1D *truth = resp->ProjectionY(Form("%s_truth", name), 1, resp->GetNbinsX());   // untruncated
  TH2D *r = (TH2D*) resp->Clone(Form("%s_resp", name));
  r->SetDirectory(nullptr);
  for(int ix = 0; ix <= r->GetNbinsX()+1; ix++)
    if(r->GetXaxis()->GetBinUpEdge(ix) <= recoFloor + 1e-6)
      for(int iy = 0; iy <= r->GetNbinsY()+1; iy++){ r->SetBinContent(ix, iy, 0.); r->SetBinError(ix, iy, 0.); }
  TH1D *meas = r->ProjectionX(Form("%s_meas", name), 1, r->GetNbinsY());
  TH1D *d = (TH1D*) data->Clone(Form("%s_in", name));
  d->SetDirectory(nullptr);
  zeroBelowFloor(d);
  RooUnfoldResponse rr(meas, truth, r, Form("%s_rr", name), "");
  RooUnfoldBayes ub(&rr, d, iter);
  TH1D *u = (TH1D*) ub.Hunfold();
  u = (TH1D*) u->Clone(name);
  u->SetDirectory(nullptr);
  return u;
}

void jetsPerZ_RCP_minBias_PFJets(bool doUnfolding = true, int iterations = 1)
{
  nIter = iterations;
  initPlotStyle();

  TString figDir = Form("%s/figures/JetsPerZ", repo);
  gSystem->mkdir(figDir, kTRUE);
  TString outDir = "./rootFiles/JetsPerZ";
  gSystem->mkdir(outDir, kTRUE);

  TFile *f_MB = TFile::Open(Form("%s/rootFiles/scanningOuput/PbPb/PbPb_MinBias_Part1_manualJEC_mu12_pTmu-15to999_tight_mu12TriggerEfficiencyCorrection_jetTrkMaxFilter_WDecayFilter_2026-10-5_ultraFineCentBins.root", repo));
  TFile *f_mu = TFile::Open(Form("%s/rootFiles/scanningOutput/PbPb/latest/PbPb_SingleMuon_mu12_pTmu-15to999_tight_mu12TriggerEfficiencyCorrection_jetTrkMaxFilter_WDecayFilter_2026-5-4.root", sib));
  if(!f_MB || f_MB->IsZombie() || !f_mu || f_mu->IsZombie()){ printf("ERROR: cannot open inputs\n"); return; }

  { TH1D *g = nullptr; f_MB->GetObject("h_nEventsNoJetTrigSel", g);
    if(g && g->GetMean() > 0.5){ printf("ERROR: this MinBias scan has a jet trigger applied\n"); return; } }

  // Z yields; the SingleMuon file is 4CentBins, so its C1-C4 ARE the coarse classes
  double NZ[NCoarse];
  for(int c = 0; c < NCoarse; c++){
    TH1D *d = nullptr; f_mu->GetObject(Form("h_dimuonMass_C%d", c+1), d);
    if(!d){ printf("ERROR: h_dimuonMass_C%d missing\n", c+1); return; }
    NZ[c] = d->Integral(d->FindBin(Z_lo + 1e-6), d->FindBin(Z_hi - 1e-6)) / (muEff[c]*muEff[c]);
  }

  TH1D *raw[NVar][NCoarse];
  for(int v = 0; v < NVar; v++)
    for(int c = 0; c < NCoarse; c++){
      raw[v][c] = coarseSum1(f_MB, varBase[v], c, varTag[v]);
      if(!raw[v][c]){ printf("ERROR: %s missing a slice for %s\n", varBase[v], coarseLabel[c]); return; }
    }

  printf("\n=== R_CP of jets per Z, PbPb MinBias Part1 only, PF jets ===\n");
  printf("  %-8s %10s %14s %14s %10s\n", "class", "N_Z", "jets>50 (all)", "(calo-match)", "match frac");
  for(int c = 0; c < NCoarse; c++){
    double a = raw[0][c]->Integral(raw[0][c]->FindBin(floorPt + 1e-6), raw[0][c]->FindBin(ptCeiling - 1e-6));
    double m = raw[1][c]->Integral(raw[1][c]->FindBin(floorPt + 1e-6), raw[1][c]->FindBin(ptCeiling - 1e-6));
    printf("  %-8s %10.1f %14.0f %14.0f %10.3f\n", coarseLabel[c], NZ[c], a, m, m/a);
  }

  // binning from the reco reference class, so unfolded and folded share edges
  std::vector<double> edges = deriveEdges(raw[0][refClass]);
  if(edges.size() < 2){ printf("ERROR: no binning\n"); return; }
  const int nb = (int)edges.size() - 1;
  printf("\n  binning (%d bins):", nb);
  for(double e : edges) printf(" %.0f", e);
  printf("\n");

  const TString uTag = doUnfolding ? Form("_unfold%diter", nIter) : "_noUnfold";
  if(doUnfolding){
    TFile *f_resp = TFile::Open(responsePath);
    if(!f_resp || f_resp->IsZombie()){ printf("ERROR: cannot open the response\n"); return; }
    printf("\n  unfolding: Bayes, %d iter, reco floor %.0f GeV, response %s\n", nIter, recoFloor, gSystem->BaseName(responsePath));
    printf("  %-8s %-12s %12s %12s %12s\n", "class", "variant", "reco>floor", "unf. 50-300", "3 iter");
    for(int c = 0; c < NCoarse; c++){
      TH2D *resp = nullptr; f_resp->GetObject(Form("h_matchedRecoJetPt_genJetPt_allJets_C%d", c+1), resp);
      if(!resp){ printf("ERROR: response C%d missing\n", c+1); return; }
      for(int v = 0; v < NVar; v++){
        TH1D *u  = unfold(raw[v][c], resp, nIter, Form("unf_%s_C%d", varTag[v], c));
        TH1D *u3 = unfold(raw[v][c], resp, 3,     Form("unf3_%s_C%d", varTag[v], c));
        double r0 = raw[v][c]->Integral(raw[v][c]->FindBin(recoFloor + 1e-6), raw[v][c]->GetNbinsX());
        auto win = [](TH1D *h){ return h->Integral(h->FindBin(floorPt + 1e-6), h->FindBin(ptCeiling - 1e-6)); };
        printf("  %-8s %-12s %12.0f %12.0f %12.0f\n", coarseLabel[c], varTag[v], r0, win(u), win(u3));
        delete u3;
        delete raw[v][c];
        raw[v][c] = u;
      }
    }
  }

  TH1D *spec[NVar][NCoarse], *rcp[NVar][NCoarse];
  for(int v = 0; v < NVar; v++){
    for(int c = 0; c < NCoarse; c++){
      spec[v][c] = rebinTo(raw[v][c], nb, &edges[0], Form("spec_%s_C%d", varTag[v], c));
      divideByBinwidth(spec[v][c]);
      spec[v][c]->Scale(1./(3.2*NZ[c]));     // eta range, per Z
    }
    for(int c = 0; c < refClass; c++){
      rcp[v][c] = (TH1D*) spec[v][c]->Clone(Form("RCP_%s_C%d", varTag[v], c));
      rcp[v][c]->SetDirectory(nullptr);
      rcp[v][c]->Divide(spec[v][refClass]);
    }
  }

  for(int v = 0; v < NVar; v++){
    printf("\n  R_CP%s, %s (ref. %s)\n    %-12s", uTag.Data(), varTag[v], coarseLabel[refClass], "pT [GeV]");
    for(int c = 0; c < refClass; c++) printf(" %16s", coarseLabel[c]);
    printf("\n");
    for(int b = 1; b <= nb; b++){
      printf("    %4.0f-%-7.0f", edges[b-1], edges[b]);
      for(int c = 0; c < refClass; c++) printf(" %8.3f+-%-6.3f", rcp[v][c]->GetBinContent(b), rcp[v][c]->GetBinError(b));
      printf("\n");
    }
  }

  // ------------------------------------------------------------- figures ---
  const char *hexC[NCoarse] = {okabeHex[5], okabeHex[1], okabeHex[2], okabeHex[0]};
  const int markF[NCoarse]  = {markFilledCircle, markFilledSquare, markFilledDiamond, markCross};
  const int markO[NCoarse]  = {markOpenCircle, markOpenSquare, markOpenDiamond, markCross};

  // spectra, one figure per variant
  for(int v = 0; v < NVar; v++){
    TCanvas *cv = new TCanvas(Form("c_spec_%s", varTag[v]), "", 700, 800);
    cv->SetLogy();
    cv->SetLeftMargin(0.17); cv->SetBottomMargin(0.12); cv->SetTopMargin(0.13); cv->SetRightMargin(0.05);
    TLegend *leg = makeLegend(0.58, 0.60, 0.93, 0.80, 0.036);
    for(int c = 0; c < NCoarse; c++){
      styleH(spec[v][c], hexC[c], markF[c]);
      if(c == 0){
        spec[v][c]->SetTitle("");
        spec[v][c]->GetXaxis()->SetTitle("jet #it{p}_{T} [GeV]");
        spec[v][c]->GetYaxis()->SetTitle("#frac{1}{#it{N}_{Z}} #frac{d#it{N}_{jet}}{d#it{p}_{T} d#eta}");
        spec[v][c]->GetYaxis()->SetTitleOffset(1.60);
        spec[v][c]->Draw("E1");
      }
      else spec[v][c]->Draw("E1 same");
      leg->AddEntry(spec[v][c], Form("PbPb %s", coarseLabel[c]), "lp");
    }
    leg->Draw();
    TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.034);
    la.DrawLatex(0.17, 0.945, Form("PbPb 5.02 TeV MinBias, PF jets, %s", varLabel[v]));
    la.DrawLatex(0.17, 0.900, doUnfolding ? Form("no fake-jet subtraction, Bayes unfolded (%d iter)", nIter) : "no fake-jet subtraction, no unfolding");
    cv->SaveAs(Form("%s/jetsPerZ_minBias_spectra_PFJets_%s%s.pdf", figDir.Data(), varTag[v], uTag.Data()));
    delete cv;
  }

  // R_CP: filled = all PF jets, open = calo-matched
  {
    TCanvas *cv = new TCanvas("c_RCP", "", 700, 800);
    cv->SetLeftMargin(0.16); cv->SetBottomMargin(0.12); cv->SetTopMargin(0.13); cv->SetRightMargin(0.05);
    TLegend *leg = makeLegend(0.52, 0.56, 0.93, 0.85, 0.032);
    TH1D *frame = (TH1D*) rcp[0][0]->Clone("frame"); frame->Reset();
    frame->SetTitle(""); double yMax = 0.;
    for(int c = 0; c < refClass; c++) yMax = std::max(yMax, rcp[0][c]->GetMaximum());
    frame->SetMinimum(0.); frame->SetMaximum(1.45*yMax);   // the low-pT fakes reach ~6
    frame->GetXaxis()->SetTitle("jet #it{p}_{T} [GeV]");
    frame->GetYaxis()->SetTitle(Form("#it{R}_{CP} of jets per #it{Z}   (ref. %s)", coarseLabel[refClass]));
    frame->GetYaxis()->SetTitleOffset(1.50);
    frame->Draw("AXIS");
    TLine *one = new TLine(edges.front(), 1., edges.back(), 1.);
    one->SetLineStyle(7); one->Draw();
    for(int c = 0; c < refClass; c++){
      styleH(rcp[0][c], hexC[c], markF[c], c == 2 ? 1.5 : 1.0);   // diamonds read small
      styleH(rcp[1][c], hexC[c], markO[c], c == 2 ? 1.5 : 1.0);
      rcp[0][c]->Draw("E1 same");
      rcp[1][c]->Draw("E1 same");
      leg->AddEntry(rcp[0][c], Form("%s all", coarseLabel[c]), "lp");
      leg->AddEntry(rcp[1][c], Form("%s calo-matched", coarseLabel[c]), "lp");
    }
    leg->Draw();
    TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.034);
    la.DrawLatex(0.16, 0.945, "PbPb 5.02 TeV MinBias Part1, PF jets");
    la.SetTextSize(0.030);
    la.DrawLatex(0.16, 0.902, doUnfolding ? Form("Bayes unfolded (%d iter), no fake jet subtraction", nIter) : "no unfolding, no fake jet subtraction");
    cv->SaveAs(Form("%s/jetsPerZ_RCP_minBias_PFJets%s.pdf", figDir.Data(), uTag.Data()));
    delete cv;
  }

  TString outPath = Form("%s/jetsPerZ_RCP_minBias_PFJets%s.root", outDir.Data(), uTag.Data());
  TFile *wf = TFile::Open(outPath, "recreate");
  for(int v = 0; v < NVar; v++){
    for(int c = 0; c < NCoarse; c++) spec[v][c]->Write(Form("jetsPerZ_%s_C%d", varTag[v], c+1));
    for(int c = 0; c < refClass; c++) rcp[v][c]->Write(Form("RCP_%s_C%d", varTag[v], c+1));
  }
  wf->Close();
  printf("\n  figures: %s/jetsPerZ_{RCP_minBias,minBias_spectra}_PFJets*.pdf\n  results: %s\n", figDir.Data(), outPath.Data());
}
