// Validate headers/functions/towerPUSubtraction.h against the forest's own
// akPu4CaloJets, event by event, on a forest that carries both the tower tree
// and akPu4CaloJetAnalyzer/t.
//
// For every forest akPu4Calo jet above ptMin (rawpt = the subtracted jet pT the
// emulation reproduces), find the closest emulated jet within dRMatch and
// compare pT, eta and the per-jet pileup (jtpu). The reverse is also counted:
// every emulated jet above ptMin, and how many have a forest jet within dRMatch.
// Forest->emulated alone cannot see EXTRA emulated jets, and total jet counts
// are what a spectrum is built from. Reported per coarse
// centrality class (hiBin from hiEvtAnalyzer/HiTree): matching efficiency, the
// mean and RMS of emulated/forest pT, and the same for jtpu.
//
// What agreement to expect. The tower tree is vertex-corrected and the forest
// clustered origin-based towers, so per-tower et and eta differ slightly and a
// perfect match is not possible from this tree -- see the header of
// towerPUSubtraction.h. Agreement at the percent level in pT with ~all jets
// matched means the emulation is the forest's algorithm; a systematic offset
// or a centrality trend means it is not.
//
// absVzMin/absVzMax select a |vz| window: the vertex-correction difference
// grows with |vz|, so the ratio RMS should shrink toward vz = 0 if that is what
// limits the agreement.
//
// Figures (figures/towers/): towerPUSub_validation_ptRatio.pdf -- emulated /
// forest pT vs forest pT per class; towerPUSub_validation_pu.pdf -- emulated
// vs forest jtpu.
//
// Needs FastJet. Run with ACLiC:
//   root -l -b -q 'validateTowerPUSub.C+("/path/to/forest.root")'
// (the rootlogon-free compile below adds the FastJet paths itself)
// Run from: src/plots/towers/

#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TSystem.h"
#include "../../../headers/plotting/plotStyle.h"
#include "../../../headers/functions/towerPUSubtraction.h"

void validateTowerPUSub(const char *path = "/home/clayton/Downloads/HiForestAOD_PbPb_MinBias_withCaloAndFlowJets_withPFAndTowers.root",
                        double ptMin = 20., double dRMatch = 0.1, Long64_t maxEvents = -1,
                        double absVzMin = 0., double absVzMax = 15.)
{
  initPlotStyle();
  const char *outDir = "../../../figures/towers/";
  gSystem->mkdir(outDir, kTRUE);

  TFile *f = TFile::Open(path);
  if(!f || f->IsZombie()){ printf("ERROR: %s did not open\n", path); return; }
  TTree *tt = nullptr, *tj = nullptr, *te = nullptr;
  f->GetObject("rechitanalyzerpp/tower", tt);
  f->GetObject("akPu4CaloJetAnalyzer/t", tj);
  f->GetObject("hiEvtAnalyzer/HiTree", te);
  if(!tt || !tj || !te){ printf("ERROR: need rechitanalyzerpp/tower, akPu4CaloJetAnalyzer/t, hiEvtAnalyzer/HiTree\n"); return; }

  static const int tMax = 10000, jMax = 500;
  static Int_t nT, ieta[tMax], iphi[tMax];
  static Float_t et[tMax], eta[tMax], phi[tMax];
  tt->SetBranchAddress("n", &nT);   tt->SetBranchAddress("et", et);
  tt->SetBranchAddress("eta", eta); tt->SetBranchAddress("phi", phi);
  tt->SetBranchAddress("ieta", ieta); tt->SetBranchAddress("iphi", iphi);
  static Int_t nJ; static Float_t rawpt[jMax], jeta[jMax], jphi[jMax], jpu[jMax];
  tj->SetBranchAddress("nref", &nJ); tj->SetBranchAddress("rawpt", rawpt);
  tj->SetBranchAddress("jteta", jeta); tj->SetBranchAddress("jtphi", jphi);
  tj->SetBranchAddress("jtpu", jpu);
  Int_t hiBin = 0; Float_t vz = 0.;
  te->SetBranchAddress("hiBin", &hiBin); te->SetBranchAddress("vz", &vz);

  // coarse classes in hiBin (0.5% units)
  const int nC = 4;
  const int cLo[nC] = {0, 20, 60, 100}, cHi[nC] = {20, 60, 100, 160};
  const char *cLab[nC] = {"0-10%", "10-30%", "30-50%", "50-80%"};
  const char *cHex[nC] = {"#000000", "#0072B2", "#009E73", "#D55E00"};
  const int   cMk[nC]  = {markFilledCircle, markFilledSquare, markFilledDiamond, markOpenCircle};

  TProfile *pRatio[nC];
  TH1D *hRatio[nC];
  TH2D *hPU[nC];
  long nForest[nC] = {0}, nMatched[nC] = {0}, nEmu[nC] = {0}, nEmuMatched[nC] = {0};
  const double edges[] = {20, 30, 40, 50, 60, 80, 100, 130, 170, 250};
  const int nE = sizeof(edges)/sizeof(double) - 1;
  for(int c = 0; c < nC; c++){
    pRatio[c] = new TProfile(Form("pRatio_%d", c), "", nE, edges, "s");
    hRatio[c] = new TH1D(Form("hRatio_%d", c), "", 200, 0.8, 1.2);
    hPU[c]    = new TH2D(Form("hPU_%d", c), "", 60, 0, 150, 60, 0, 150);
    pRatio[c]->SetDirectory(nullptr); hRatio[c]->SetDirectory(nullptr); hPU[c]->SetDirectory(nullptr);
  }

  const Long64_t N = (maxEvents > 0) ? TMath::Min(maxEvents, tt->GetEntries()) : tt->GetEntries();
  printf("  %lld events\n", N);
  for(Long64_t i = 0; i < N; i++){
    tt->GetEntry(i); tj->GetEntry(i); te->GetEntry(i);
    if(std::fabs(vz) > absVzMax || std::fabs(vz) < absVzMin) continue;
    int cls = -1;
    for(int c = 0; c < nC; c++) if(hiBin >= cLo[c] && hiBin < cHi[c]) cls = c;
    if(cls < 0) continue;

    std::vector<TowerPUInput> in;
    in.reserve(nT);
    for(int k = 0; k < nT && k < tMax; k++)
      if(et[k] >= 0.3) in.push_back({et[k], eta[k], phi[k], ieta[k], iphi[k]});
    std::vector<TowerPUJet> emu = clusterTowersWithPUSub(in, 0.4);

    for(int j = 0; j < nJ && j < jMax; j++){
      if(rawpt[j] < ptMin || std::fabs(jeta[j]) > 1.6) continue;
      nForest[cls]++;
      int best = -1; double bestDR = dRMatch;
      for(size_t e = 0; e < emu.size(); e++){
        const double d = towerPU::dR(emu[e].eta, emu[e].phi, jeta[j], jphi[j]);
        if(d < bestDR){ bestDR = d; best = (int)e; }
      }
      if(best < 0) continue;
      nMatched[cls]++;
      const double r = emu[best].pt / rawpt[j];
      pRatio[cls]->Fill(rawpt[j], r);
      hRatio[cls]->Fill(r);
      hPU[cls]->Fill(jpu[j], emu[best].pu);
    }
    for(const auto &e : emu){
      if(e.pt < ptMin || std::fabs(e.eta) > 1.6) continue;
      nEmu[cls]++;
      for(int j = 0; j < nJ && j < jMax; j++)
        if(towerPU::dR(e.eta, e.phi, jeta[j], jphi[j]) < dRMatch){ nEmuMatched[cls]++; break; }
    }
    if(i % 1000 == 0) printf("    event %lld\n", i);
  }

  printf("\n  forest akPu4Calo jets, rawpt > %.0f GeV, |eta| < 1.6, matched within dR < %.2f\n", ptMin, dRMatch);
  printf("  %-8s %8s %8s %8s %10s %10s %10s\n", "class", "forest", "matched", "eff", "<ratio>", "RMS", "median");
  for(int c = 0; c < nC; c++){
    double q = 0.5, med = 0.;
    if(hRatio[c]->GetEntries() > 0) hRatio[c]->GetQuantiles(1, &med, &q);
    printf("  %-8s %8ld %8ld %8.3f %10.4f %10.4f %10.4f\n", cLab[c], nForest[c], nMatched[c],
           nForest[c] ? (double)nMatched[c]/nForest[c] : 0., hRatio[c]->GetMean(), hRatio[c]->GetRMS(), med);
  }
  printf("\n  jet counts above %.0f GeV, |eta| < 1.6 (a spectrum is built from these)\n", ptMin);
  printf("  %-8s %8s %8s %8s %14s\n", "class", "forest", "emulated", "emu/for", "emu w/ match");
  for(int c = 0; c < nC; c++)
    printf("  %-8s %8ld %8ld %8.3f %14ld\n", cLab[c], nForest[c], nEmu[c],
           nForest[c] ? (double)nEmu[c]/nForest[c] : 0., nEmuMatched[c]);
  printf("  (ratio histogram range 0.8-1.2; entries outside are in under/overflow and not in the mean)\n");
  for(int c = 0; c < nC; c++)
    printf("  %-8s ratio underflow %.0f overflow %.0f\n", cLab[c], hRatio[c]->GetBinContent(0), hRatio[c]->GetBinContent(201));

  // ---- pT ratio vs forest pT ----
  TCanvas *c1 = new TCanvas("c1", "", 800, 700);
  c1->SetLeftMargin(0.15); c1->SetBottomMargin(0.13); c1->SetRightMargin(0.05);
  TH1D *fr = new TH1D("fr", "", nE, edges);
  fr->SetStats(0); fr->SetMinimum(0.9); fr->SetMaximum(1.1);
  fr->GetXaxis()->SetTitle("forest akPu4Calo raw #it{p}_{T} [GeV]");
  fr->GetYaxis()->SetTitle("emulated / forest raw #it{p}_{T}");
  fr->GetXaxis()->SetTitleSize(0.05); fr->GetYaxis()->SetTitleSize(0.05);
  fr->GetXaxis()->SetLabelSize(0.042); fr->GetYaxis()->SetLabelSize(0.042);
  fr->GetYaxis()->SetTitleOffset(1.4);
  fr->Draw();
  TLine *one = new TLine(edges[0], 1., edges[nE], 1.); one->SetLineStyle(2); one->SetLineColor(kGray+2); one->Draw();
  TLegend *lg = makeLegend(0.60, 0.16, 0.93, 0.40);
  for(int c = 0; c < nC; c++){
    TH1D *h = pRatio[c]->ProjectionX(Form("pr_%d", c));
    styleH(h, cHex[c], cMk[c], 1.0);
    h->Draw("E SAME");
    lg->AddEntry(h, cLab[c], "pl");
  }
  lg->Draw();
  TLatex l; l.SetNDC(); l.SetTextSize(0.036);
  l.DrawLatex(0.19, 0.87, "Tower-level PU subtraction emulation vs forest");
  l.DrawLatex(0.19, 0.82, Form("PbPb MinBias, |#it{#eta}| < 1.6, #DeltaR < %.2f; error bar = RMS", dRMatch));
  savePdfTight(c1, Form("%stowerPUSub_validation_ptRatio.pdf", outDir));

  // ---- jtpu, emulated vs forest, 0-10% ----
  TCanvas *c2 = new TCanvas("c2", "", 800, 700);
  c2->SetLeftMargin(0.15); c2->SetBottomMargin(0.13); c2->SetRightMargin(0.15);
  gStyle->SetPalette(kViridis);
  hPU[0]->SetStats(0);
  hPU[0]->GetXaxis()->SetTitle("forest jtpu [GeV]");
  hPU[0]->GetYaxis()->SetTitle("emulated pileup [GeV]");
  hPU[0]->GetXaxis()->SetTitleSize(0.05); hPU[0]->GetYaxis()->SetTitleSize(0.05);
  hPU[0]->GetYaxis()->SetTitleOffset(1.4);
  hPU[0]->Draw("COLZ");
  TLine *diag = new TLine(0, 0, 150, 150); diag->SetLineStyle(2); diag->SetLineColor(kGray+2); diag->Draw();
  l.DrawLatex(0.19, 0.87, "0-10%, matched jets");
  savePdfTight(c2, Form("%stowerPUSub_validation_pu.pdf", outDir));
}
