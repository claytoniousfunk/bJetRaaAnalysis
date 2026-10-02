// Gen-level effect of the JEU on the PYTHIA+HYDJET calo-jet spectrum in the 16
// ultraFine (5%) centrality slices, 0-80%. The fine-binned counterpart of
// unfoldPYTHIAHYDJETJetPt_JEUShift.C (read that header for what is measured and
// the caveats): each slice's nominal / JEU-up / JEU-down reco spectrum is
// unfolded with that SAME slice's nominal response, and shifted / nominal is
// taken at gen level.
//
// Reco spectra  h_inclRecoJetPt_flavor_{nominal,JEUShiftUp,JEUShiftDown}_C<s>
//               of a PYTHIAHYDJET_scan ultraFine output (argument 1), flavors
//               summed. Slices 17, 18 (80-90%, beyond) are not used.
// Response      the full-sample ultraFine PYTHIAHYDJET_scan_response (all jets,
//               not muon-tagged, no even/odd split), found by name tag in
//               rootFiles/scanningOuput/PYTHIAHYDJET/ (argument 2, e.g. the date).
//               h_matchedRecoJetPt_genJetPt_allJets_C<s>, misses from
//               h_unmatchedGenJetPt_C<s> in the truth. Reco floor 50 GeV, Bayes,
//               2 iterations. No fakes in the response; the scan's reco spectrum
//               has them (largest in central slices at low pT), and the two are
//               different runs, so only the ratio is meaningful.
//
// Output  figures/unfoldTest/unfold_PYTHIAHYDJET_caloJets_JEUShift_vsCentrality.pdf
//         shifted / nominal at gen level vs centrality, in three jet-pT windows
//         (up filled, down open), plus a printed table that also gives the
//         iteration dependence (N = 1, 2, 4) and the reco-level ratio.
//         Error bars are the measured-statistics error of the shifted spectrum
//         in the window (shifted and nominal share the same jets; the unfolded
//         covariance is not propagated).
//
// Usage: root -l -b -q -e 'gSystem->Load("/home/clayton/Programs/RooUnfold/build/libRooUnfold.so");
//                          gInterpreter->AddIncludePath("/home/clayton/Programs/RooUnfold/build");'
//                     'unfoldPYTHIAHYDJETJetPt_JEUShift_fineCent.C+("<scan.root>","2026-10-2")'
// Run from: src/unfoldTest/

#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#endif
#include "../../headers/plotting/plotStyle.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TSystemDirectory.h"
#include "TList.h"
#include <cstdio>
#include <string>
#include <vector>
#include <cmath>

namespace {

const char *repo   = "/home/clayton/Analysis/code/bJetRaaAnalysis";
const char *mcDir  = "rootFiles/scanningOuput/PYTHIAHYDJET/";
const char *outPdf = "/home/clayton/Analysis/code/bJetRaaAnalysis/figures/unfoldTest/unfold_PYTHIAHYDJET_caloJets_JEUShift_vsCentrality.pdf";

const double recoFloor = 50.;
const int    nIterNominal = 2;
const int    nIterShow[3] = {1, 2, 4};
const int    nSlice = 16;           // 5% slices 1..16 = 0-80%

const int    nWin = 3;
const double winLo[nWin] = {100, 150, 250};
const double winHi[nWin] = {150, 250, 400};

TH1* getH(TFile *f, const char *name)
{
  TH1 *h = nullptr; f->GetObject(name, h);
  if(!h) printf("ERROR: %s missing in %s\n", name, f->GetName());
  return h;
}

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

TH1D* recoOf(TH2D *H, const char *tag)
{
  TH1D *h = H->ProjectionX(tag, 1, H->GetNbinsY());
  h->SetDirectory(nullptr);
  for(int b = 0; b <= h->GetNbinsX() + 1; b++)
    if(h->GetXaxis()->GetBinCenter(b) < recoFloor){ h->SetBinContent(b, 0.); h->SetBinError(b, 0.); }
  return h;
}

// counts in [lo, hi) on the 5 GeV grid; FindBin(hi) starts at hi, so use hi - eps
double window(TH1D *h, double lo, double hi, double *err = nullptr)
{
  int b1 = h->GetXaxis()->FindBin(lo + 1e-6), b2 = h->GetXaxis()->FindBin(hi - 1e-6);
  double e = 0.;
  double v = h->IntegralAndError(b1, b2, e);
  if(err) *err = e;
  return v;
}

} // namespace

void unfoldPYTHIAHYDJETJetPt_JEUShift_fineCent(const char *scanFile, const char *respTag)
{
  initPlotStyle();
  TFile *fs = TFile::Open(scanFile);
  if(!fs || fs->IsZombie()){ printf("ERROR: cannot open %s\n", scanFile); return; }
  std::string rf = findScan({"PYTHIAHYDJET_response_DiJet_caloJets", "manualJEC", "ultraFineCentBins",
                             "doPThatCorrelationFilterTight", respTag},
                            {"muTaggedJets", "evenEvents", "oddEvents"});
  if(rf.empty()){ printf("ERROR: no ultraFine all-jets response matching '%s'\n", respTag); return; }
  printf("response: %s\nreco:     %s\n", rf.c_str(), scanFile);
  TFile *fr = TFile::Open(rf.c_str());
  if(!fr || fr->IsZombie()) return;

  const char *var[3] = {"nominal", "JEUShiftUp", "JEUShiftDown"};
  // [window][up/down] graphs of gen ratio vs centrality
  TGraphErrors *g[nWin][2];
  for(int w = 0; w < nWin; w++) for(int s = 0; s < 2; s++){ g[w][s] = new TGraphErrors(); }

  printf("\nslice  cent   | window GeV |  up: reco  gen(N=1,2,4)            | down: reco  gen(N=1,2,4)\n");
  for(int si = 1; si <= nSlice; si++){
    TH2D *resp2D = (TH2D*) getH(fr, Form("h_matchedRecoJetPt_genJetPt_allJets_C%d", si));
    TH1D *miss   = (TH1D*) getH(fr, Form("h_unmatchedGenJetPt_C%d", si));
    if(!resp2D || !miss) return;
    resp2D = (TH2D*) resp2D->Clone(Form("resp2D_%d", si)); resp2D->SetDirectory(nullptr);
    TH1D *truth = resp2D->ProjectionY(Form("truth_%d", si), 1, resp2D->GetNbinsX());
    truth->SetDirectory(nullptr); truth->Add(miss);
    for(int ix = 0; ix <= resp2D->GetNbinsX(); ix++)
      if(resp2D->GetXaxis()->GetBinCenter(ix) < recoFloor)
        for(int iy = 0; iy <= resp2D->GetNbinsY() + 1; iy++){ resp2D->SetBinContent(ix, iy, 0.); resp2D->SetBinError(ix, iy, 0.); }
    TH1D *measMC = resp2D->ProjectionX(Form("measMC_%d", si), 1, resp2D->GetNbinsY());
    measMC->SetDirectory(nullptr);
    RooUnfoldResponse resp(measMC, truth, resp2D, Form("resp_S%d", si), "PYTHIAHYDJET calo response");

    TH1D *reco[3], *gen[3][5];
    for(int v = 0; v < 3; v++){
      TH2D *H = (TH2D*) getH(fs, Form("h_inclRecoJetPt_flavor_%s_C%d", var[v], si));
      if(!H){ printf("ERROR: the scan lacks the JEU histograms\n"); return; }
      reco[v] = recoOf(H, Form("reco_%d_S%d", v, si));
      RooUnfoldBayes unfold(&resp, reco[v], 1);
      unfold.SetVerbose(0);
      for(int n = 1; n <= 4; n++){
        unfold.SetIterations(n);
        gen[v][n] = (TH1D*) unfold.Hunfold();
        gen[v][n]->SetName(Form("gen_%d_S%d_N%d", v, si, n)); gen[v][n]->SetDirectory(nullptr);
      }
    }

    const double cent = 5.*(si - 1) + 2.5;
    for(int w = 0; w < nWin; w++){
      double vals[2][3], rr[2];
      for(int s = 0; s < 2; s++){
        double eUp, nR = window(reco[0], winLo[w], winHi[w]), sR = window(reco[s+1], winLo[w], winHi[w], &eUp);
        rr[s] = nR > 0 ? sR / nR : 0.;
        for(int k = 0; k < 3; k++){
          double n = window(gen[0][nIterShow[k]], winLo[w], winHi[w]), sh = window(gen[s+1][nIterShow[k]], winLo[w], winHi[w]);
          vals[s][k] = n > 0 ? sh / n : 0.;
        }
        double relErr = sR > 0 ? eUp / sR : 0.;
        int np = g[w][s]->GetN();
        g[w][s]->SetPoint(np, cent, vals[s][1]);
        g[w][s]->SetPointError(np, 0., vals[s][1] * relErr);
      }
      printf("%4d  %4.1f%%  | %3.0f-%-3.0f    |  %.3f   %.3f %.3f %.3f        |  %.3f   %.3f %.3f %.3f\n",
             si, cent, winLo[w], winHi[w], rr[0], vals[0][0], vals[0][1], vals[0][2], rr[1], vals[1][0], vals[1][1], vals[1][2]);
    }
    delete resp2D;
  }

  // ---- figure ---------------------------------------------------------------
  const char *hex[nWin]  = {okabeHex[1], okabeHex[5], okabeHex[3]};   // orange, blue, green
  const int   mkF[nWin]  = {markFilledCircle, markFilledSquare, markFilledDiamond};
  const int   mkO[nWin]  = {markOpenCircle,   markOpenSquare,   markOpenDiamond};
  TCanvas *c = new TCanvas("c", "", 800, 700);
  c->SetLeftMargin(0.14); c->SetBottomMargin(0.13); c->SetTopMargin(0.06); c->SetRightMargin(0.04);
  TH1F *fr0 = c->DrawFrame(0, 0.6, 80, 1.5);
  fr0->GetXaxis()->SetTitle("Centrality (%)"); fr0->GetYaxis()->SetTitle("Shifted / nominal, unfolded");
  fr0->GetXaxis()->SetTitleSize(0.05); fr0->GetYaxis()->SetTitleSize(0.05);
  fr0->GetXaxis()->SetLabelSize(0.045); fr0->GetYaxis()->SetLabelSize(0.045);
  fr0->GetYaxis()->SetTitleOffset(1.35);
  TLine one; one.SetLineStyle(7); one.SetLineColor(kGray+1); one.DrawLine(0, 1., 80, 1.);
  TLegend *leg = makeLegend(0.50, 0.40, 0.95, 0.62, 0.04);
  for(int w = 0; w < nWin; w++) for(int s = 0; s < 2; s++){
    int col = TColor::GetColor(hex[w]);
    g[w][s]->SetMarkerColor(col); g[w][s]->SetLineColor(col);
    g[w][s]->SetMarkerStyle(s == 0 ? mkF[w] : mkO[w]); g[w][s]->SetMarkerSize(1.2); g[w][s]->SetLineWidth(2);
    g[w][s]->Draw("PE same");
    if(s == 0) leg->AddEntry(g[w][s], Form("%.0f-%.0f GeV", winLo[w], winHi[w]), "p");
  }
  leg->Draw();
  TLatex t; t.SetNDC(); t.SetTextSize(0.045);
  t.DrawLatex(0.18, 0.88, "PYTHIA+HYDJET calo jets");
  t.SetTextSize(0.038);
  t.DrawLatex(0.18, 0.83, Form("Bayes, %d iterations; filled: JEU up, open: JEU down", nIterNominal));

  gSystem->mkdir(gSystem->DirName(outPdf), kTRUE);
  savePdfTight(c, outPdf);
}
