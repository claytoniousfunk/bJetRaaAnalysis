// pp CALO-jet b-jet spectrum: muon-tagged jet spectrum x b purity.
//
// Follows constructBJetSpectra() in calculateBJetsPerZ.cc for the PF chain:
// the muon-tagged, muon-triggered jet spectrum (h_inclRecoJetPt_inclRecoMuonTag_triggerOn)
// is rebinned to the purity's jet-pT bins, divided by bin width, and multiplied
// bin by bin by the b purity. Calo data with calo purity only.
//
// Uncertainties:
//   stat  TH1::Multiply, as in the PF chain: the spectrum's counting error and
//         the purity's fit error added in relative quadrature. Both pieces are
//         also stored separately.
//   syst  two independent sources, added in quadrature:
//     b purity       the spectrum is purity x counts, so every purity source
//                    carries over as the SAME relative uncertainty; their total
//                    is the quadrature sum from systematics_total_caloJets_pp.C.
//     spectrum JEU   the jet-pT spectrum's own dependence on the jet energy
//                    scale, from src/unfoldTest/unfoldPYTHIAJetPt_JEUShift.C:
//                    the PYTHIA calo reco spectrum shifted up/down by the JEU,
//                    unfolded with the nominal response, shifted / nominal at
//                    gen level, nominal number of iterations (N = 2). The
//                    systematic is max(|up - 1|, |down - 1|) per bin.
//   CAVEATS on the spectrum JEU: it is measured on the INCLUSIVE (all-flavor)
//   PYTHIA calo-jet spectrum, not the muon-tagged b one; the shift uses the
//   AK4PF uncertainty file (no pp AK4Calo one exists); and the same JEU also
//   enters the b-purity JEU source, so treating the two as independent ignores
//   their correlation.
// No efficiency / correction factors other than the turn-on below, no unfolding, no per-Z normalisation:
// this is the purity step only.
//
// The purity is fitted on 0 < ptRel < 5 GeV and applied to the full spectrum
// (all ptRel), as in the PF chain.
//
// Underflow bins (50-60, 60-70, 70-80 GeV) are kept for the unfolding and not
// reported. 50-60 and 60-70 are empty: the pp SingleMuon calo data has no jets
// below 70 GeV (the calo forest's jet threshold), so the unfolding floor is 70.
//
// TURN-ON CORRECTION, 70-80 GeV. That threshold is not sharp: 70-80 is a turn-on.
// Per 5 GeV bin the efficiency is the data / PYTHIA ratio of the muon-tagged,
// trigger-on spectrum shape, normalised on the plateau 80-90 GeV:
//   eff_k = [N_data(k) / N_data(80-90)] / [N_MC(k) / N_MC(80-90)]
// (data h_inclRecoJetPt_inclRecoMuonTag_triggerOn; MC the same histogram, all
// flavours, _flavor projected). About 0.14 at 70-75 and 0.88 at 75-80 GeV. The
// 5 GeV counts are divided by it BEFORE rebinning. The b purity is a ratio and
// is not affected. Without this the half-empty bin was unfolded as if fully
// efficient, and since muon-tagged calo b jets reconstruct at ~0.7 of gen pT it
// pulled gen 80-120 GeV low. Systematic (h_bJetPt_sysRel_turnOn, nonzero only in
// 70-80): max |change| of the corrected bin with the plateau window 80-100 or
// 80-120 GeV, in quadrature with the efficiency's statistical error.
//
// Inputs:  pp SingleMuon calo data (dataPath in caloJetBPurityFit_pp.h)
//          rootFiles/systematics/bPuritySys_total_caloJets_pp.root
//          rootFiles/JEU/unfold_PYTHIA_caloJets_JEUShift.root
// Output:  rootFiles/MuTagMuTrigBJetSpectra/Data/bJetSpectrum_caloJets_pp.root
//            h_muTagJetPt_counts            muon-tagged jets per bin (counts)
//            h_muTagJetPt                   muon-tagged dN/dpT
//            h_bPurity                      nominal purity, fit error
//            h_bJetPt                       b-jet dN/dpT, total stat error (counts + purity fit)
//            h_bJetPt_sysAbs                b-jet dN/dpT, error = TOTAL systematic
//            h_bJetPt_statRel_counts        relative stat error from the counts
//            h_bJetPt_statRel_purity        relative stat error from the purity fit
//            h_bJetPt_sysRel                relative TOTAL systematic (purity (+) spectrum JEU)
//            h_bJetPt_sysRel_bPurity        relative systematic from the b purity (all sources)
//            h_bJetPt_sysRel_spectrumJEU    relative systematic from the spectrum JEU
//            h_bJetPt_sysRel_spectrumJEU_up/_down  signed shifted/nominal - 1
//            h_bJetPt_sysRel_<source>       relative systematic per b-purity source
//          figures/bJetSpectra/bJetSpectrum_caloJets_pp.pdf
//          figures/systematics/bJetSpectrum_sysBreakdown_caloJets_pp.pdf
//
// Usage: root -l -b -q constructBJetSpectrum_caloJets_pp.C
// Run from: src/calculateBJetsPerZ/

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include <vector>
#include "TCanvas.h"
#include "TPad.h"
#include "TBox.h"
#include "TLatex.h"
#include "TSystem.h"
#include "TColor.h"
#include "TStyle.h"
#include "../../headers/plotting/plotStyle.h"
#include "../../headers/plotting/ratioPanel.h"
#include "../../headers/functions/divideByBinwidth.h"
#include "../../headers/functions/caloJetBPurityFit_pp.h"   // dataPath, ptEdges, ptReportMin, underflow helpers

const char *sysTotalPath = "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/systematics/bPuritySys_total_caloJets_pp.root";
const char *jeuSpecPath  = "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/JEU/unfold_PYTHIA_caloJets_JEUShift.root";
const int   jeuIter      = 2;   // unfoldPYTHIAJetPt_JEUShift.C's nominal number of iterations
const char *outRoot      = "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/MuTagMuTrigBJetSpectra/Data/bJetSpectrum_caloJets_pp.root";
const char *outFig       = "/home/clayton/Analysis/code/bJetRaaAnalysis/figures/bJetSpectra/bJetSpectrum_caloJets_pp.pdf";
const char *outFigSys    = "/home/clayton/Analysis/code/bJetRaaAnalysis/figures/systematics/bJetSpectrum_sysBreakdown_caloJets_pp.pdf";

const char *specName = "h_inclRecoJetPt_inclRecoMuonTag_triggerOn";
const char *specNameMC = "h_inclRecoJetPt_inclRecoMuonTag_triggerOn_flavor";   // in mcPath, all flavours
const double turnOnLo = 70., turnOnHi = 80.;            // 5 GeV bins corrected
const double plateauHi[3] = {90., 100., 120.};          // nominal, variations; plateau starts at turnOnHi
const int NSrc = 6;
const char *srcKey[NSrc] = {"bGS", "cMult", "fitRangeLow", "fitRangeHigh", "JEU", "JER"};

TH1D* bookRel(const char *name, const char *title)
{
  TH1D *h = new TH1D(name, Form("%s;#it{p}_{T}^{jet} [GeV];relative uncertainty", title), NPt, ptEdges);
  h->SetDirectory(nullptr); return h;
}

// counts of h summed into ptEdges window i (5 GeV input bins)
double windowSum(TH1D *h, int i)
{
  return h->Integral(h->FindBin(ptEdges[i] + 1e-6), h->FindBin(ptEdges[i+1] - 1e-6));
}

void constructBJetSpectrum_caloJets_pp()
{
  initPlotStyle();
  gSystem->mkdir(gSystem->DirName(outRoot), kTRUE);
  gSystem->mkdir(gSystem->DirName(outFig), kTRUE);
  gSystem->mkdir(gSystem->DirName(outFigSys), kTRUE);

  TFile *fD = TFile::Open(dataPath), *fS = TFile::Open(sysTotalPath), *fJ = TFile::Open(jeuSpecPath);
  if(!fD || fD->IsZombie()){ printf("ERROR: cannot open data %s\n", dataPath); return; }
  if(!fS || fS->IsZombie()){ printf("ERROR: cannot open %s -- run systematics_total_caloJets_pp.C first\n", sysTotalPath); return; }
  if(!fJ || fJ->IsZombie()){ printf("ERROR: cannot open %s -- run src/unfoldTest/unfoldPYTHIAJetPt_JEUShift.C first\n", jeuSpecPath); return; }

  TH1D *hSpecIn = nullptr, *hPur = nullptr, *hSysPur = nullptr, *hSysSrc[NSrc];
  fD->GetObject(specName, hSpecIn);
  fS->GetObject("h_bPurity_nominal", hPur);
  fS->GetObject("h_sysRel_total", hSysPur);
  if(!hSpecIn || !hPur || !hSysPur){ printf("ERROR: missing input histogram\n"); return; }
  for(int s = 0; s < NSrc; s++){
    fS->GetObject(Form("h_sysRel_%s", srcKey[s]), hSysSrc[s]);
    if(!hSysSrc[s]){ printf("ERROR: missing h_sysRel_%s\n", srcKey[s]); return; }
  }
  TH1D *gNom = nullptr, *gUp = nullptr, *gDown = nullptr;
  fJ->GetObject(Form("gen_nominal_N%d", jeuIter), gNom);
  fJ->GetObject(Form("gen_up_N%d", jeuIter), gUp);
  fJ->GetObject(Form("gen_down_N%d", jeuIter), gDown);
  if(!gNom || !gUp || !gDown){ printf("ERROR: gen_{nominal,up,down}_N%d missing in %s\n", jeuIter, jeuSpecPath); return; }

  // the purity file must be in exactly the bins used here
  if(hPur->GetNbinsX() != NPt){ printf("ERROR: purity has %d bins, expected %d\n", hPur->GetNbinsX(), NPt); return; }
  for(int i = 0; i <= NPt; i++){
    double edge = (i < NPt) ? hPur->GetXaxis()->GetBinLowEdge(i+1) : hPur->GetXaxis()->GetBinUpEdge(NPt);
    if(fabs(edge - ptEdges[i]) > 1e-6){ printf("ERROR: purity bin edge %d is %.1f, expected %.1f\n", i, edge, ptEdges[i]); return; }
  }

  // ---- turn-on correction, 70-80 GeV (5 GeV bins) -----------------------------
  TFile *fMC = TFile::Open(mcPath);
  TH2D *hMC2 = nullptr; if(fMC && !fMC->IsZombie()) fMC->GetObject(specNameMC, hMC2);
  if(!hMC2){ printf("ERROR: %s missing in %s\n", specNameMC, mcPath); return; }
  TH1D *hMCt = hMC2->ProjectionX("hMCtagged", 0, hMC2->GetNbinsY()+1); hMCt->SetDirectory(nullptr);
  const int kLo = hSpecIn->FindBin(turnOnLo + 1e-6), kHi = hSpecIn->FindBin(turnOnHi - 1e-6);
  const int nK = kHi - kLo + 1;
  TH1D *hEff = new TH1D("h_turnOnEff", "calo-jet turn-on efficiency (data / PYTHIA muon-tagged shape, plateau 80-90 GeV);#it{p}_{T}^{jet} [GeV];efficiency",
                        nK, hSpecIn->GetBinLowEdge(kLo), hSpecIn->GetBinLowEdge(kHi+1));
  hEff->SetDirectory(nullptr);
  std::vector<std::vector<double>> effV(3, std::vector<double>(nK, 1.));
  for(int v = 0; v < 3; v++){
    int r1 = hSpecIn->FindBin(turnOnHi + 1e-6), r2 = hSpecIn->FindBin(plateauHi[v] - 1e-6);
    double nd = hSpecIn->Integral(r1, r2), nm = hMCt->Integral(hMCt->FindBin(turnOnHi + 1e-6), hMCt->FindBin(plateauHi[v] - 1e-6));
    for(int k = 0; k < nK; k++){
      int b = kLo + k, bm = hMCt->FindBin(hSpecIn->GetBinCenter(b));
      double d = hSpecIn->GetBinContent(b), m = hMCt->GetBinContent(bm);
      double e = (d/nd)/(m/nm);
      effV[v][k] = e;
      if(v == 0){
        double re = sqrt(pow(hSpecIn->GetBinError(b)/d, 2) + pow(hMCt->GetBinError(bm)/m, 2));
        hEff->SetBinContent(k+1, e); hEff->SetBinError(k+1, re*e);
      }
    }
  }
  // corrected 5 GeV counts: nominal and the two plateau variations
  TH1D *hSpecCorr = (TH1D*) hSpecIn->Clone("hSpecCorr"); hSpecCorr->SetDirectory(nullptr);
  TH1D *hSpecVar[2];
  for(int v = 0; v < 2; v++){ hSpecVar[v] = (TH1D*) hSpecIn->Clone(Form("hSpecVar%d", v)); hSpecVar[v]->SetDirectory(nullptr); }
  for(int k = 0; k < nK; k++){
    int b = kLo + k;
    hSpecCorr->SetBinContent(b, hSpecIn->GetBinContent(b)/effV[0][k]); hSpecCorr->SetBinError(b, hSpecIn->GetBinError(b)/effV[0][k]);
    for(int v = 0; v < 2; v++) hSpecVar[v]->SetBinContent(b, hSpecIn->GetBinContent(b)/effV[v+1][k]);
    printf("turn-on %3.0f-%-3.0f GeV: eff %.4f +/- %.4f (plateau 80-90), %.4f (80-100), %.4f (80-120)\n",
           hSpecIn->GetBinLowEdge(b), hSpecIn->GetBinLowEdge(b+1), effV[0][k], hEff->GetBinError(k+1), effV[1][k], effV[2][k]);
  }
  TH1D *hCountsRaw = rebinTo(hSpecIn, NPt, ptEdges, "h_muTagJetPt_counts_noTurnOnCorr");
  hCountsRaw->SetTitle("muon-tagged jets, no turn-on correction (counts);#it{p}_{T}^{jet} [GeV];jets per bin");
  TH1D *hCV0 = rebinTo(hSpecVar[0], NPt, ptEdges, "hCV0"), *hCV1 = rebinTo(hSpecVar[1], NPt, ptEdges, "hCV1");
  // the relative stat error of the efficiency, summed into each analysis bin
  TH1D *hEffStatAbs = (TH1D*) hSpecIn->Clone("hEffStatAbs"); hEffStatAbs->SetDirectory(nullptr); hEffStatAbs->Reset();
  for(int k = 0; k < nK; k++) hEffStatAbs->SetBinContent(kLo + k, hSpecCorr->GetBinContent(kLo + k)*hEff->GetBinError(k+1)/effV[0][k]);

  // ---- muon-tagged spectrum in the purity bins ------------------------------
  TH1D *hCounts = rebinTo(hSpecCorr, NPt, ptEdges, "h_muTagJetPt_counts");
  hCounts->SetTitle("muon-tagged jets (counts);#it{p}_{T}^{jet} [GeV];jets per bin");
  TH1D *hSpec = (TH1D*) hCounts->Clone("h_muTagJetPt"); hSpec->SetDirectory(nullptr);
  divideByBinwidth(hSpec);
  hSpec->SetTitle("muon-tagged jets;#it{p}_{T}^{jet} [GeV];d#it{N}/d#it{p}_{T} [GeV^{-1}]");

  TH1D *hP = (TH1D*) hPur->Clone("h_bPurity"); hP->SetDirectory(nullptr);

  // ---- b-jet spectrum ---------------------------------------------------------
  TH1D *hB = (TH1D*) hSpec->Clone("h_bJetPt"); hB->SetDirectory(nullptr);
  hB->Multiply(hP);   // stat: counts and purity-fit errors in relative quadrature
  hB->SetTitle("b jets (muon-tagged x b purity), stat. error;#it{p}_{T}^{jet} [GeV];d#it{N}/d#it{p}_{T} [GeV^{-1}]");

  TH1D *hBsys = (TH1D*) hB->Clone("h_bJetPt_sysAbs"); hBsys->SetDirectory(nullptr);
  hBsys->SetTitle("b jets, error = total systematic;#it{p}_{T}^{jet} [GeV];d#it{N}/d#it{p}_{T} [GeV^{-1}]");

  TH1D *hRelC   = bookRel("h_bJetPt_statRel_counts",          "relative stat. error from the muon-tagged counts");
  TH1D *hRelP   = bookRel("h_bJetPt_statRel_purity",          "relative stat. error from the b-purity fit");
  TH1D *hRelPur = bookRel("h_bJetPt_sysRel_bPurity",          "relative systematic from the b purity (all purity sources)");
  TH1D *hRelJ   = bookRel("h_bJetPt_sysRel_spectrumJEU",      "relative systematic from the spectrum JEU, max(|up-1|,|down-1|)");
  TH1D *hRelJu  = bookRel("h_bJetPt_sysRel_spectrumJEU_up",   "spectrum JEU up: shifted/nominal - 1 (gen level)");
  TH1D *hRelJd  = bookRel("h_bJetPt_sysRel_spectrumJEU_down", "spectrum JEU down: shifted/nominal - 1 (gen level)");
  TH1D *hRelS   = bookRel("h_bJetPt_sysRel",                  "relative total systematic (b purity (+) spectrum JEU (+) turn-on)");
  TH1D *hRelT   = bookRel("h_bJetPt_sysRel_turnOn",           "relative systematic from the 70-80 GeV turn-on correction");
  TH1D *hRelSrc[NSrc];
  for(int s = 0; s < NSrc; s++){
    hRelSrc[s] = (TH1D*) hSysSrc[s]->Clone(Form("h_bJetPt_sysRel_%s", srcKey[s])); hRelSrc[s]->SetDirectory(nullptr);
    hRelSrc[s]->SetTitle(Form("relative systematic on the b-jet spectrum from the b purity, %s;#it{p}_{T}^{jet} [GeV];relative uncertainty", srcKey[s]));
  }

  printf("\n%-9s %10s %8s %12s %8s %8s %8s %8s %8s %8s %8s\n", "jet pT", "muTag N", "purity", "b dN/dpT",
         "stat", "sys pur", "JEU up", "JEU dn", "sys JEU", "turnOn", "sys tot");
  for(int i = 1; i <= NPt; i++){
    double n = hCounts->GetBinContent(i), p = hP->GetBinContent(i), b = hB->GetBinContent(i);
    double rc = n > 0. ? hCounts->GetBinError(i)/n : 0., rp = p > 0. ? hP->GetBinError(i)/p : 0.;
    double rpur = hSysPur->GetBinContent(i);
    double gn = windowSum(gNom, i-1);
    double ju = gn > 0. ? windowSum(gUp, i-1)/gn - 1. : 0., jd = gn > 0. ? windowSum(gDown, i-1)/gn - 1. : 0.;
    double rj = TMath::Max(fabs(ju), fabs(jd));
    double rt = 0.;
    if(n > 0.){
      double ev = TMath::Max(fabs(hCV0->GetBinContent(i)/n - 1.), fabs(hCV1->GetBinContent(i)/n - 1.));
      double es = windowSum(hEffStatAbs, i-1)/n;   // efficiency stat, conservatively added linearly over the 5 GeV bins
      rt = sqrt(ev*ev + es*es);
    }
    hRelT->SetBinContent(i, rt);
    double rs = sqrt(rpur*rpur + rj*rj + rt*rt);
    hRelC->SetBinContent(i, rc); hRelP->SetBinContent(i, rp);
    hRelPur->SetBinContent(i, rpur); hRelJ->SetBinContent(i, rj); hRelJu->SetBinContent(i, ju); hRelJd->SetBinContent(i, jd);
    hRelS->SetBinContent(i, rs);
    hBsys->SetBinError(i, rs*b);
    printf("%3.0f-%-5.0f %10.0f %8.3f %12.4g %7.2f%% %7.2f%% %+7.2f%% %+7.2f%% %7.2f%% %7.2f%% %7.2f%%%s\n", ptEdges[i-1], ptEdges[i], n, p, b,
           b > 0. ? 100*hB->GetBinError(i)/b : 0., 100*rpur, 100*ju, 100*jd, 100*rj, 100*rt, 100*rs,
           ptEdges[i] <= ptReportMin + 1e-6 ? "   (underflow)" : "");
  }

  TFile *fo = TFile::Open(outRoot, "recreate");
  hCounts->Write(); hSpec->Write(); hP->Write(); hB->Write(); hBsys->Write();
  hCountsRaw->Write(); hEff->Write(); hRelT->Write();
  hRelC->Write(); hRelP->Write(); hRelS->Write(); hRelPur->Write(); hRelJ->Write(); hRelJu->Write(); hRelJd->Write();
  for(int s = 0; s < NSrc; s++) hRelSrc[s]->Write();
  TNamed info("info", Form("pp calo jets: b-jet dN/dpT = muon-tagged dN/dpT (%s) x b purity. Stat via TH1::Multiply "
                            "(counts + purity fit). Syst = b purity (quadrature of %d sources, same relative size) (+) "
                            "spectrum JEU (gen-level shifted/nominal, N=%d, from %s; inclusive PYTHIA calo jets, AK4PF "
                            "uncertainty; correlated with the purity JEU source, treated as independent). No efficiency, "
                            "unfolding or per-Z normalisation. 70-80 GeV corrected for the calo-jet turn-on (data/PYTHIA muon-tagged shape, "
                            "plateau 80-90 GeV; syst from plateau 80-100/80-120 and eff. stat). Bins 50-80 GeV are underflow; 50-70 empty in data. data: %s",
                            specName, NSrc, jeuIter, gSystem->BaseName(jeuSpecPath), gSystem->BaseName(dataPath)));
  info.Write();
  fo->Close();
  printf("\nwritten %s\n", outRoot);

  int colSys = TColor::GetColor("#D55E00");

  // ---- figure 1: the spectrum ------------------------------------------------
  // Reported range only: the underflow bins (below ptReportMin) are not shown
  // here; they stay in the output file and in the breakdown figure.
  auto reported = [](int i){ return ptEdges[i-1] >= ptReportMin - 1e-6; };
  {
    TCanvas *c = new TCanvas("c", "", 700, 800);
    TPad *pTop, *pBot; splitPads(pTop, pBot);
    pTop->cd(); pTop->SetLogy();   // log: the spectrum falls by ~3 orders of magnitude and has no negative bins
    double ymax = 0., ymin = 1e30;
    for(int i = 1; i <= NPt; i++){ double v = hSpec->GetBinContent(i); if(v > 0. && reported(i)){ ymax = TMath::Max(ymax, v); ymin = TMath::Min(ymin, hB->GetBinContent(i)); } }
    TH1D *fr = (TH1D*) hB->Clone("frame"); fr->Reset(); fr->SetTitle("");
    fr->SetMinimum(0.3*ymin); fr->SetMaximum(30.*ymax);
    fr->GetXaxis()->SetLabelSize(0);
    fr->GetXaxis()->SetRangeUser(ptReportMin, ptEdges[NPt]);
    fr->GetYaxis()->SetTitle("d#it{N}/d#it{p}_{T} [GeV^{-1}]"); fr->GetYaxis()->SetTitleSize(0.055);
    fr->GetYaxis()->SetTitleOffset(1.35); fr->GetYaxis()->SetLabelSize(0.045);
    fr->Draw("AXIS");
    for(int i = 1; i <= NPt; i++){
      double b = hB->GetBinContent(i); if(b <= 0. || !reported(i)) continue;
      TBox *bx = new TBox(ptEdges[i-1], b - hBsys->GetBinError(i), ptEdges[i], b + hBsys->GetBinError(i));
      bx->SetFillColorAlpha(colSys, 0.30); bx->SetLineColor(colSys); bx->Draw("l"); bx->Draw();
    }
    TH1D *hSd = (TH1D*) hSpec->Clone("hSd"), *hBd = (TH1D*) hB->Clone("hBd");
    styleH(hSd, "#0072B2", markOpenCircle, 1.2);
    styleH(hBd, "#000000", markFilledCircle, 1.2);
    for(TH1D *h : {hSd, hBd}) for(int i = 1; i <= NPt; i++) if(h->GetBinContent(i) <= 0. || !reported(i)){ h->SetBinContent(i, -999.); h->SetBinError(i, 0.); }
    hSd->Draw("E1 X0 same"); hBd->Draw("E1 X0 same");
    TLegend *leg = makeLegend(0.25, 0.05, 0.60, 0.24, 0.036);   // lower left is empty: the spectrum falls
    leg->AddEntry(hSd, "muon-tagged jets", "lp");
    leg->AddEntry(hBd, "b jets (#times b purity), stat.", "lp");
    TBox *lb = new TBox(); lb->SetFillColorAlpha(colSys, 0.30); lb->SetLineColor(colSys);
    leg->AddEntry(lb, "total systematic", "f");
    leg->Draw();
    TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.046);
    la.DrawLatex(0.21, 0.84, "pp 5.02 TeV, calo jets");
    la.SetTextSize(0.036);
    la.DrawLatex(0.21, 0.78, "no efficiency correction, not unfolded");

    pBot->cd();
    TH1D *rS = (TH1D*) hRelS->Clone("rS"), *rT = (TH1D*) hRelS->Clone("rT");
    for(int i = 1; i <= NPt; i++){ double b = hB->GetBinContent(i); rT->SetBinContent(i, b > 0. ? hB->GetBinError(i)/b : 0.); }
    for(TH1D *h : {rS, rT}){ h->Scale(100.); for(int i = 0; i <= NPt+1; i++) h->SetBinError(i, 0.); }
    styleLine(rS, "#D55E00", 3); styleLine(rT, "#0072B2", 3); rT->SetLineStyle(2);
    double m = 0.;
    for(int i = 1; i <= NPt; i++) if(reported(i)) m = TMath::Max(m, TMath::Max(rS->GetBinContent(i), rT->GetBinContent(i)));
    TH1D *rF = (TH1D*) rS->Clone("rF"); rF->Reset();
    rF->SetMinimum(0.); rF->SetMaximum(1.7*m);   // headroom for the legend row
    rF->GetXaxis()->SetRangeUser(ptReportMin, ptEdges[NPt]);
    styleRatioAxes(rF, "#it{p}_{T}^{jet} [GeV]", "rel. unc. [%]");
    rF->Draw("AXIS");
    // HIST lines would start at 50 GeV; restrict them to the reported range too
    rS->GetXaxis()->SetRangeUser(ptReportMin, ptEdges[NPt]); rT->GetXaxis()->SetRangeUser(ptReportMin, ptEdges[NPt]);
    rS->Draw("HIST same"); rT->Draw("HIST same");
    TLegend *legB = makeLegend(0.19, 0.78, 0.95, 0.97, 0.075);
    legB->SetNColumns(2);
    legB->AddEntry(rS, "total syst.", "l"); legB->AddEntry(rT, "total stat.", "l");
    legB->Draw();
    c->SaveAs(outFig);
    delete c;
  }

  // ---- figure 2: systematic sources on the spectrum ---------------------------
  {
    TCanvas *c = new TCanvas("cSys", "", 800, 700);
    c->SetLeftMargin(0.13); c->SetBottomMargin(0.13); c->SetRightMargin(0.04); c->SetTopMargin(0.06);
    TH1D *lT = (TH1D*) hRelS->Clone("lT"), *lP = (TH1D*) hRelPur->Clone("lP"), *lJ = (TH1D*) hRelJ->Clone("lJ");
    TH1D *pU = (TH1D*) hRelJu->Clone("pU"), *pD = (TH1D*) hRelJd->Clone("pD");
    for(TH1D *h : {lT, lP, lJ, pU, pD}){ h->Scale(100.); for(int i = 0; i <= NPt+1; i++) h->SetBinError(i, 0.); }
    // up/down are signed; shown as magnitudes so they sit against the max(|.|) line
    for(TH1D *h : {pU, pD}) for(int i = 1; i <= NPt; i++) h->SetBinContent(i, fabs(h->GetBinContent(i)));
    double ymax = lT->GetMaximum();
    TH1F *fr = c->DrawFrame(ptEdges[0], 0., ptEdges[NPt], 1.8*ymax);   // headroom for header + legend
    fr->GetXaxis()->SetTitle("#it{p}_{T}^{jet} [GeV]");
    fr->GetYaxis()->SetTitle("relative systematic uncertainty [%]");
    fr->GetXaxis()->SetTitleSize(0.048); fr->GetYaxis()->SetTitleSize(0.048);
    fr->GetXaxis()->SetLabelSize(0.040); fr->GetYaxis()->SetLabelSize(0.040);
    drawUnderflowBand(0., 1.8*ymax, 0.05*ymax);
    styleLine(lP, "#009E73", 3); lP->SetLineStyle(2);
    styleLine(lJ, "#0072B2", 3); lJ->SetLineStyle(3);
    styleLine(lT, "#000000", 4);
    styleH(pU, "#D55E00", markOpenDiamond, 1.5);
    styleH(pD, "#CC79A7", markOpenSquare, 1.2);
    lP->Draw("HIST same"); lJ->Draw("HIST same"); lT->Draw("HIST same");
    pU->Draw("P same"); pD->Draw("P same");
    TLegend *leg = makeLegend(0.17, 0.66, 0.95, 0.85, 0.032);
    leg->SetNColumns(2);
    leg->AddEntry(lP, "b purity (all purity sources)", "l");
    leg->AddEntry(pU, "spectrum JEU up, |#Delta|", "p");
    leg->AddEntry(lJ, "spectrum JEU, max(|up|,|down|)", "l");
    leg->AddEntry(pD, "spectrum JEU down, |#Delta|", "p");
    leg->AddEntry(lT, "total (quadrature)", "l");
    leg->Draw();
    TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.038);
    la.DrawLatex(0.17, 0.875, "pp 5.02 TeV, calo jets: b-jet spectrum systematics");
    savePdfTight(c, outFigSys);
    delete c;
  }
  printf("written %s\nwritten %s\n", outFig, outFigSys);
}
