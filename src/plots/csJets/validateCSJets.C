// Validate headers/functions/csJetEmulation.h against the forest's own
// akCs4PFJets, event by event, on a forest that carries pfcandAnalyzer/pfTree,
// pfcandAnalyzerCS/pfTree and akCs4PFJetAnalyzer/t.
//
// Two tests, both meant to agree to float precision, not to a percent:
//   JETS. For every forest akCs4PF jet with rawpt > ptMin and |eta| < 1.6,
//     the closest emulated jet within dRMatch; emulated / forest rawpt, and
//     the eta and phi differences. The reverse (emulated jets with a forest
//     partner) catches EXTRA emulated jets. It has its own threshold,
//     revPtMin: the analyzer keeps only jets with CORRECTED jtpt >= 20 GeV
//     (rawpt down to ~9 GeV), so an emulated jet at raw 10-20 GeV legitimately
//     has no partner. 30 GeV raw is safely above that cut.
//   PARTICLES. pfcandAnalyzerCS holds akCs4PFJets:pfParticlesCs, the subtracted
//     constituents of every CS jet in the event. The emulation runs on every
//     jet (ptUnsubMin = 0) here so the two lists are the same set. Keys do not
//     match the raw PF keys, so particles are matched by phi (the subtraction
//     keeps rapidity and phi) within dPhiMatch, nearest in eta; pT compared.
//     Only |eta| < partEtaMax: the local MinBias forest's pfcandAnalyzer
//     keeps |eta| < 3 only (no HF candidates at all), while the forest CS ran
//     on HF too, so HF particles, and jets straddling |eta| = 3, cannot be
//     reproduced from this tree. The rho of the |eta| > 2.1 bins can be
//     affected for the same reason; the bins a |eta| < 1.6 jet touches
//     (ghosts within 0.4, |eta| < 2.0) are not.
// Reported per coarse class (hiBin). The fraction of jets with
// |emulated/forest - 1| below 1e-5 / 1e-3 / 1e-2 is the headline number: an
// exact emulation puts essentially everything below 1e-5 (float storage, and
// the forest's -ffast-math build, are what is left).
//
// The local forest has no hiFJRhoAnalyzer tree, so rho cannot be compared
// directly. The first printEvents events print the emulated rho map.
//
// Figures (figures/csJets/): csJets_validation_ptRatio.pdf -- log10 |emulated /
// forest rawpt - 1| per class; csJets_validation_particles.pdf -- the same for
// the matched CS particles.
//
// Needs FastJet 3.3.0 + fjcontrib 1.033 (see csJetEmulation.h). Run from
// src/plots/csJets/ with the libraries loaded first:
//   root -l -b -q loadFastJet330.C 'validateCSJets.C+("/path/to/forest.root")'

#include <vector>
#include <algorithm>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLine.h"
#include "TSystem.h"
#include "TMath.h"
#include "../../../headers/plotting/plotStyle.h"
#include "../../../headers/functions/csJetEmulation.h"

namespace {
double dPhi(double a, double b){ double d = std::fabs(a - b); return d > TMath::Pi() ? 2. * TMath::Pi() - d : d; }
double dR(double e1, double p1, double e2, double p2){ const double de = e1 - e2, dp = dPhi(p1, p2); return std::sqrt(de * de + dp * dp); }
double lg(double r){ const double d = std::fabs(r - 1.); return d > 0. ? std::log10(d) : -9.99; }
}

void validateCSJets(const char *path = "/home/clayton/Downloads/HiForestAOD_PbPb_MinBias_withCaloAndFlowJets_withPFAndTowers.root",
                    double ptMin = 10., double dRMatch = 0.05, Long64_t maxEvents = -1,
                    double dPhiMatch = 1e-4, int printEvents = 3,
                    double revPtMin = 30., double partEtaMax = 2.5)
{
  initPlotStyle();
  const char *outDir = "../../../figures/csJets/";
  gSystem->mkdir(outDir, kTRUE);

  TFile *f = TFile::Open(path);
  if(!f || f->IsZombie()){ printf("ERROR: %s did not open\n", path); return; }
  TTree *tp = nullptr, *tc = nullptr, *tj = nullptr, *te = nullptr;
  f->GetObject("pfcandAnalyzer/pfTree", tp);
  f->GetObject("pfcandAnalyzerCS/pfTree", tc);
  f->GetObject("akCs4PFJetAnalyzer/t", tj);
  f->GetObject("hiEvtAnalyzer/HiTree", te);
  if(!tp || !tc || !tj || !te){ printf("ERROR: need pfcandAnalyzer, pfcandAnalyzerCS, akCs4PFJetAnalyzer, hiEvtAnalyzer\n"); return; }

  std::vector<int> *pId = nullptr;
  std::vector<float> *pPt = nullptr, *pEta = nullptr, *pPhi = nullptr, *pM = nullptr;
  tp->SetBranchAddress("pfId", &pId);   tp->SetBranchAddress("pfPt", &pPt);
  tp->SetBranchAddress("pfEta", &pEta); tp->SetBranchAddress("pfPhi", &pPhi);
  tp->SetBranchAddress("pfM", &pM);
  std::vector<float> *cPt = nullptr, *cEta = nullptr, *cPhi = nullptr;
  tc->SetBranchAddress("pfPt", &cPt); tc->SetBranchAddress("pfEta", &cEta); tc->SetBranchAddress("pfPhi", &cPhi);
  static const int jMax = 2000;
  static Int_t nJ; static Float_t rawpt[jMax], jeta[jMax], jphi[jMax];
  tj->SetBranchAddress("nref", &nJ); tj->SetBranchAddress("rawpt", rawpt);
  tj->SetBranchAddress("jteta", jeta); tj->SetBranchAddress("jtphi", jphi);
  UInt_t run = 0; ULong64_t evt = 0; Int_t hiBin = 0;
  te->SetBranchAddress("run", &run); te->SetBranchAddress("evt", &evt); te->SetBranchAddress("hiBin", &hiBin);

  const int nC = 4;
  const int cLo[nC] = {0, 20, 60, 100}, cHi[nC] = {20, 60, 100, 160};
  const char *cLab[nC] = {"0-10%", "10-30%", "30-50%", "50-80%"};
  const char *cHex[nC] = {"#000000", "#0072B2", "#009E73", "#D55E00"};

  TH1D *hJet[nC], *hPart[nC];
  long nForest[nC] = {0}, nMatched[nC] = {0}, nEmu[nC] = {0}, nEmuMatched[nC] = {0};
  long nBelow[nC][3] = {{0}};
  long nCsForest[nC] = {0}, nCsEmu[nC] = {0}, nCsMatched[nC] = {0}, nCsBelow[nC][3] = {{0}};
  double maxDEta[nC] = {0}, maxDPhi[nC] = {0};
  const double thr[3] = {1e-5, 1e-3, 1e-2};
  for(int c = 0; c < nC; c++){
    hJet[c]  = new TH1D(Form("hJet_%d", c),  "", 100, -10, 0);
    hPart[c] = new TH1D(Form("hPart_%d", c), "", 100, -10, 0);
    hJet[c]->SetDirectory(nullptr); hPart[c]->SetDirectory(nullptr);
  }

  const Long64_t N = (maxEvents > 0) ? TMath::Min(maxEvents, tp->GetEntries()) : tp->GetEntries();
  printf("  %lld events\n", N);
  int printed = 0;
  for(Long64_t i = 0; i < N; i++){
    tp->GetEntry(i); tc->GetEntry(i); tj->GetEntry(i); te->GetEntry(i);
    int cls = -1;
    for(int c = 0; c < nC; c++) if(hiBin >= cLo[c] && hiBin < cHi[c]) cls = c;
    if(cls < 0) continue;

    std::vector<CSJetInput> in;
    in.reserve(pPt->size());
    for(size_t k = 0; k < pPt->size(); k++)
      in.push_back({pPt->at(k), pEta->at(k), pPhi->at(k), pM->at(k), pId->at(k)});
    CSRhoMap rho;
    std::vector<CSJet> emu = csJet::forestCSJets(in, run, (unsigned int)evt, 0., &rho);

    if(printed < printEvents){
      printf("  run %u evt %llu hiBin %d: %zu PF, %zu emulated CS jets, rho / rhom per eta bin:\n",
             run, (unsigned long long)evt, hiBin, in.size(), emu.size());
      for(size_t b = 0; b < rho.rho.size(); b++)
        printf("    [%5.1f, %5.1f)  %9.4f  %9.5f\n", rho.etaEdges[b], rho.etaEdges[b + 1], rho.rho[b], rho.rhom[b]);
      printed++;
    }

    // ---- jets ----
    for(int j = 0; j < nJ && j < jMax; j++){
      if(rawpt[j] < ptMin || std::fabs(jeta[j]) > 1.6) continue;
      nForest[cls]++;
      int best = -1; double bestDR = dRMatch;
      for(size_t e = 0; e < emu.size(); e++){
        const double d = dR(emu[e].eta, emu[e].phi, jeta[j], jphi[j]);
        if(d < bestDR){ bestDR = d; best = (int)e; }
      }
      if(best < 0) continue;
      nMatched[cls]++;
      const double r = emu[best].pt / rawpt[j];
      hJet[cls]->Fill(TMath::Max(lg(r), -9.99));
      for(int t = 0; t < 3; t++) if(std::fabs(r - 1.) < thr[t]) nBelow[cls][t]++;
      maxDEta[cls] = TMath::Max(maxDEta[cls], std::fabs(emu[best].eta - jeta[j]));
      maxDPhi[cls] = TMath::Max(maxDPhi[cls], dPhi(emu[best].phi, jphi[j]));
    }
    for(const CSJet &e : emu){
      if(e.pt < revPtMin || std::fabs(e.eta) > 1.6) continue;
      nEmu[cls]++;
      for(int j = 0; j < nJ && j < jMax; j++)
        if(dR(e.eta, e.phi, jeta[j], jphi[j]) < dRMatch){ nEmuMatched[cls]++; break; }
    }

    // ---- CS particles ----
    struct P { double pt, eta, phi; };
    std::vector<P> ep;
    for(const CSJet &e : emu) for(const CSJetConstituent &c : e.constituents)
      if(std::fabs(c.eta) < partEtaMax) ep.push_back({c.pt, c.eta, c.phi});
    std::sort(ep.begin(), ep.end(), [](const P &a, const P &b){ return a.phi < b.phi; });
    nCsEmu[cls] += ep.size();
    for(size_t k = 0; k < cPt->size(); k++){
      const double phi = cPhi->at(k), eta = cEta->at(k);
      if(std::fabs(eta) >= partEtaMax) continue;
      nCsForest[cls]++;
      auto lo = std::lower_bound(ep.begin(), ep.end(), phi - dPhiMatch, [](const P &a, double v){ return a.phi < v; });
      const P *best = nullptr; double bestDEta = 1e9;
      for(auto it = lo; it != ep.end() && it->phi < phi + dPhiMatch; ++it)
        if(std::fabs(it->eta - eta) < bestDEta){ bestDEta = std::fabs(it->eta - eta); best = &*it; }
      if(!best || bestDEta > 1e-2) continue;
      nCsMatched[cls]++;
      const double r = best->pt / cPt->at(k);
      hPart[cls]->Fill(TMath::Max(lg(r), -9.99));
      for(int t = 0; t < 3; t++) if(std::fabs(r - 1.) < thr[t]) nCsBelow[cls][t]++;
    }
    if(i % 500 == 0) printf("    event %lld\n", i);
  }

  printf("\n  forest akCs4PF jets, rawpt > %.0f GeV, |eta| < 1.6, matched within dR < %.2f\n", ptMin, dRMatch);
  printf("  %-8s %8s %8s %8s %9s %9s %9s %10s %10s\n", "class", "forest", "matched", "eff",
         "<1e-5", "<1e-3", "<1e-2", "max|deta|", "max|dphi|");
  for(int c = 0; c < nC; c++){
    const double m = nMatched[c] ? (double)nMatched[c] : 1.;
    printf("  %-8s %8ld %8ld %8.4f %9.4f %9.4f %9.4f %10.2e %10.2e\n", cLab[c], nForest[c], nMatched[c],
           nForest[c] ? nMatched[c] / (double)nForest[c] : 0.,
           nBelow[c][0] / m, nBelow[c][1] / m, nBelow[c][2] / m, maxDEta[c], maxDPhi[c]);
  }
  printf("\n  emulated jets, raw pT > %.0f GeV, |eta| < 1.6, with a forest jet within dR < %.2f\n", revPtMin, dRMatch);
  printf("  %-8s %8s %14s %8s\n", "class", "emulated", "emu w/ match", "frac");
  for(int c = 0; c < nC; c++)
    printf("  %-8s %8ld %14ld %8.4f\n", cLab[c], nEmu[c], nEmuMatched[c],
           nEmu[c] ? nEmuMatched[c] / (double)nEmu[c] : 0.);
  printf("\n  CS particles, |eta| < %.1f (pfcandAnalyzerCS vs emulated, every jet), matched in phi within %.0e\n", partEtaMax, dPhiMatch);
  printf("  %-8s %10s %10s %10s %9s %9s %9s\n", "class", "forest", "emulated", "matched", "<1e-5", "<1e-3", "<1e-2");
  for(int c = 0; c < nC; c++){
    const double m = nCsMatched[c] ? (double)nCsMatched[c] : 1.;
    printf("  %-8s %10ld %10ld %10ld %9.4f %9.4f %9.4f\n", cLab[c], nCsForest[c], nCsEmu[c], nCsMatched[c],
           nCsBelow[c][0] / m, nCsBelow[c][1] / m, nCsBelow[c][2] / m);
  }

  // ---- figures ----
  TH1D **hs[2] = {hJet, hPart};
  const char *xt[2] = {"log_{10} |emulated / forest raw #it{p}_{T} - 1|", "log_{10} |emulated / forest CS particle #it{p}_{T} - 1|"};
  const char *nm[2] = {"csJets_validation_ptRatio.pdf", "csJets_validation_particles.pdf"};
  const char *ti[2] = {Form("akCs4PF jets, raw #it{p}_{T} > %.0f GeV, |#it{#eta}| < 1.6", ptMin),
                       Form("pfcandAnalyzerCS particles, |#it{#eta}| < %.1f", partEtaMax)};
  for(int k = 0; k < 2; k++){
    TCanvas *cv = new TCanvas(Form("cv%d", k), "", 800, 700);
    cv->SetLeftMargin(0.15); cv->SetBottomMargin(0.13); cv->SetRightMargin(0.05);
    double mx = 0.;
    for(int c = 0; c < nC; c++){
      if(hs[k][c]->Integral() > 0) hs[k][c]->Scale(1. / hs[k][c]->Integral());
      mx = TMath::Max(mx, hs[k][c]->GetMaximum());
    }
    TLegend *leg = makeLegend(0.18, 0.55, 0.45, 0.78);
    for(int c = 0; c < nC; c++){
      TH1D *h = hs[k][c];
      h->SetStats(0); h->SetMaximum(1.3 * mx); h->SetMinimum(0.);
      h->SetLineColor(TColor::GetColor(cHex[c])); h->SetLineWidth(3);
      h->GetXaxis()->SetTitle(xt[k]); h->GetYaxis()->SetTitle("fraction");
      h->GetXaxis()->SetTitleSize(0.045); h->GetYaxis()->SetTitleSize(0.05);
      h->Draw(c ? "HIST SAME" : "HIST");
      leg->AddEntry(h, cLab[c], "l");
    }
    leg->Draw();
    TLatex l; l.SetNDC(); l.SetTextSize(0.036);
    l.DrawLatex(0.18, 0.87, "akCs4PF emulation vs forest, PbPb MinBias");
    l.DrawLatex(0.18, 0.82, ti[k]);
    savePdfTight(cv, Form("%s%s", outDir, nm[k]));
  }
}
