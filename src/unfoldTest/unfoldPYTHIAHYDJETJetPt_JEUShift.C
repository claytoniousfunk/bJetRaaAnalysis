// Gen-level effect of the JEU on the PYTHIA+HYDJET calo-jet spectrum, per
// centrality class: unfold the nominal and JEU-shifted reco spectra with ONE
// response (the nominal one, same class) and take shifted / nominal at gen level.
// The PbPb counterpart of unfoldPYTHIAJetPt_JEUShift.C; read that header for what
// the ratio measures.
//
// INPUTS are three sets of PYTHIAHYDJET_scan_response scans (even + odd halves
// each): nominal, apply_JEU_shift_up and apply_JEU_shift_down (config_PYTHIAHYDJET.h;
// the output name gains _applyJEUShiftUp / _applyJEUShiftDown). The scan shifts
// the matched reco jet pT before filling h_matchedRecoJetPt_genJetPt_allJets_C<i>,
// so
//   response  = nominal  h_matchedRecoJetPt_genJetPt_allJets_C<i>, with misses
//               from h_unmatchedGenJetPt_C<i> in the truth
//   reco spec = x-projection of each of the three (nominal / up / down)
// Everything is the same sample and selection (the pThat-correlation filter keys
// on the nominal leading jet pT, so it does not move with the shift). i = 1..4
// -> 0-10, 10-30, 30-50, 50-80%. PYTHIAHYDJET_scan.C is NOT used: it has no calo
// support. Uncertainty file: Autumn18_HI_V8_MC_Uncertainty_AK4Calo.
//
// Reco floor 50 GeV, Bayes, 2 iterations nominal; N = 1..4 printed.
//
// Caveats. The response holds no fakes, so the shifted "measured" spectra are
// matched jets only and fakes are not part of this test; the real data carry
// them, and in 0-10% they are largest at low pT. The shift is applied to matched
// reco jets only. Ratio errors are the shifted-spectrum statistics only (same
// jets as nominal; covariance not propagated).
//
// Arguments are substrings that select the files in the PYTHIAHYDJET scan
// directory, e.g. the date of each generation ("2026-9-25", "2026-10-3"); the
// newest match by modification time wins. nominalTag must match a file WITHOUT
// _applyJ..., and shiftTag one WITH each of _applyJEUShiftUp / Down.
//
// Output  figures/unfoldTest/unfold_PYTHIAHYDJET_caloJets_JEUShiftRatio.pdf
//         (one panel per class: gen ratio points, reco ratio dashed)
//
// Usage: root -l -b -q -e 'gSystem->Load("/home/clayton/Programs/RooUnfold/build/libRooUnfold.so");
//                          gInterpreter->AddIncludePath("/home/clayton/Programs/RooUnfold/build");'
//                     'unfoldPYTHIAHYDJETJetPt_JEUShift.C+("2026-9-25","<date of shifted scans>")'
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
#include "TSystemDirectory.h"
#include "TList.h"
#include "TSystemFile.h"
#include <cstdio>
#include <string>

namespace {

const char *repo   = "/home/clayton/Analysis/code/bJetRaaAnalysis";
const char *mcDir  = "rootFiles/scanningOuput/PYTHIAHYDJET/";
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

// shifted / nominal in wide bins; error = shifted-spectrum relative error of the
// numerator only (see header)
TH1D* wideRatio(TH1D *num, TH1D *den, const char *name)
{
  TH1D *n = wide(num, Form("%s_n", name)), *d = wide(den, Form("%s_d", name));
  TH1D *r = makeRatio(n, d, name);
  delete n; delete d;
  return r;
}


// newest file in the scan directory whose name contains every substring of `must`
// (comma-free list) and none of `mustNot`
std::string findScan(const std::vector<std::string> &must, const std::vector<std::string> &mustNot)
{
  std::string dir = Form("%s/%s", repo, mcDir), best; Long64_t bestT = -1;
  TSystemDirectory d("d", dir.c_str());
  TList *l = d.GetListOfFiles();
  if(!l) return best;
  for(TObject *o : *l){
    std::string n = o->GetName();
    if(n.size() < 5 || n.substr(n.size() - 5) != ".root") continue;
    bool ok = true;
    for(auto &m : must)    if(n.find(m) == std::string::npos) ok = false;
    for(auto &m : mustNot) if(n.find(m) != std::string::npos) ok = false;
    if(!ok) continue;
    FileStat_t st;
    if(gSystem->GetPathInfo((dir + n).c_str(), st)) continue;
    if(st.fMtime > bestT){ bestT = st.fMtime; best = dir + n; }
  }
  return best;
}

// summed even + odd response for one variation; shift = "", "_applyJEUShiftUp", "_applyJEUShiftDown"
bool openHalves(const char *genTag, const char *shift, TFile *&fe, TFile *&fo)
{
  fe = fo = nullptr;
  const char *half[2] = {"evenEvents", "oddEvents"};
  TFile **out[2] = {&fe, &fo};
  for(int h = 0; h < 2; h++){
    std::vector<std::string> must = {"PYTHIAHYDJET_response_DiJet_caloJets", "manualJEC", half[h],
                                     "doPThatCorrelationFilterTight", genTag}, mustNot = {"muTaggedJets"};
    if(shift[0]) must.push_back(shift); else mustNot.push_back("_applyJ");
    std::string f = findScan(must, mustNot);
    if(f.empty()){ printf("ERROR: no %s scan matching tag '%s', shift '%s'\n", half[h], genTag, shift); return false; }
    printf("  %s\n", f.c_str());
    *out[h] = TFile::Open(f.c_str());
    if(!*out[h] || (*out[h])->IsZombie()) return false;
  }
  return true;
}

TH2D* sumHalves(TFile *fe, TFile *fo, const char *name, const char *tag)
{
  TH2D *a = (TH2D*) getH(fe, name), *b = (TH2D*) getH(fo, name);
  if(!a || !b) return nullptr;
  TH2D *r = (TH2D*) a->Clone(tag); r->SetDirectory(nullptr); r->Add(b);
  return r;
}

TH1D* recoOf(TH2D *resp, const char *tag)
{
  TH1D *h = resp->ProjectionX(tag, 1, resp->GetNbinsY());
  h->SetDirectory(nullptr);
  for(int b = 0; b <= h->GetNbinsX() + 1; b++)
    if(h->GetXaxis()->GetBinCenter(b) < recoFloor){ h->SetBinContent(b, 0.); h->SetBinError(b, 0.); }
  return h;
}

const char *centLabel[5] = {"", "0-10%", "10-30%", "30-50%", "50-80%"};

} // namespace

void unfoldPYTHIAHYDJETJetPt_JEUShift(const char *nominalTag, const char *shiftTag)
{
  initPlotStyle();
  TFile *fe[3], *fo[3];
  const char *shiftName[3] = {"", "_applyJEUShiftUp", "_applyJEUShiftDown"};
  printf("input scans:\n");
  for(int v = 0; v < 3; v++)
    if(!openHalves(v == 0 ? nominalTag : shiftTag, shiftName[v], fe[v], fo[v])) return;

  TCanvas *c = new TCanvas("c", "", 1100, 900);
  c->Divide(2, 2, 0.001, 0.001);

  for(int ci = 1; ci <= 4; ci++){
    TH2D *resp2D = sumHalves(fe[0], fo[0], Form("h_matchedRecoJetPt_genJetPt_allJets_C%d", ci), Form("resp2D_%d", ci));
    TH1D *miss   = nullptr;
    { TH1D *a = (TH1D*) getH(fe[0], Form("h_unmatchedGenJetPt_C%d", ci)), *b = (TH1D*) getH(fo[0], Form("h_unmatchedGenJetPt_C%d", ci));
      if(!a || !b || !resp2D) return;
      miss = (TH1D*) a->Clone(Form("miss_%d", ci)); miss->SetDirectory(nullptr); miss->Add(b); }

    TH1D *reco[3];
    for(int v = 0; v < 3; v++){
      TH2D *R = (v == 0) ? resp2D : sumHalves(fe[v], fo[v], Form("h_matchedRecoJetPt_genJetPt_allJets_C%d", ci), Form("respShift%d_%d", v, ci));
      if(!R) return;
      reco[v] = recoOf(R, Form("reco_%d_C%d", v, ci));
    }

    TH1D *truth = resp2D->ProjectionY(Form("truth_%d", ci), 1, resp2D->GetNbinsX());
    truth->SetDirectory(nullptr); truth->Add(miss);
    for(int ix = 0; ix <= resp2D->GetNbinsX(); ix++)
      if(resp2D->GetXaxis()->GetBinCenter(ix) < recoFloor)
        for(int iy = 0; iy <= resp2D->GetNbinsY() + 1; iy++){ resp2D->SetBinContent(ix, iy, 0.); resp2D->SetBinError(ix, iy, 0.); }
    TH1D *measMC = resp2D->ProjectionX(Form("measMC_%d", ci), 1, resp2D->GetNbinsY());
    measMC->SetDirectory(nullptr);
    RooUnfoldResponse resp(measMC, truth, resp2D, Form("resp_C%d", ci), "PYTHIAHYDJET calo response");

    TH1D *gen[3][nIterMax + 1];
    for(int v = 0; v < 3; v++){
      RooUnfoldBayes unfold(&resp, reco[v], 1);
      unfold.SetVerbose(0);
      for(int n = 1; n <= nIterMax; n++){
        unfold.SetIterations(n);
        gen[v][n] = (TH1D*) unfold.Hunfold();
        gen[v][n]->SetName(Form("gen_%d_C%d_N%d", v, ci, n)); gen[v][n]->SetDirectory(nullptr);
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
