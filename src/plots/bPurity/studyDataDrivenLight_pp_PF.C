// pp PF b purity with a DATA-DRIVEN light template, ATLAS-style (arXiv:2204.13530:
// light-jet ptRel from track-jet pairs in inclusive-jet data). 2026-10-09 study, user.
// Light = the leading-track ptRel of pp jet60 data (sibling repo,
// pp_jetTrkMax_trkPt-14_jet60_updatedJEC_2025-11-24.root, h_jetTrkMaxPtRel_J<j>, track
// pT > 14 GeV, our jet-pT windows; NOT reweighted in track pT -- no track-pT axis;
// the mean track pT is close to the light-muon pT). c and b from PYTHIA (nominal
// per-window smearing), c at its MC share of non-b. Fits 0-3 and 0-5 GeV.
// Variant: the UNSMEARED PYTHIA light-muon template x (data / PYTHIA leading-track
// ptRel), PYTHIA = ..._jetTrkMax_pThat-20_..._2025-8-21.root.
// Finding (2026-10-09): data leading tracks have a broader ptRel than PYTHIA's at the
// same track pT (+0.07-0.09 GeV mean, +30-40% at 1.5-3 GeV); replacing the light
// template brings the pp purity to PYTHIA truth within 4-10%; the reweighting variant
// moves it only 2-5%.
// Usage, from src/plots/bPurity/:  root -l -b -q studyDataDrivenLight_pp_PF.C
#include "RooMsgService.h"
#include "../../../headers/functions/caloJetBPurityFit_pp.h"
void studyDataDrivenLight_pp_PF(){
  setJetCollection("PFJets"); RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
  TFile *fD = TFile::Open(dataPath), *fM = TFile::Open(mcPath);
  TFile *fTM = TFile::Open("/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/scanningOutput/PYTHIA/latest/jetTrkMax/PYTHIA_DiJet_jetTrkMax_pThat-20_mu12_pTmu-14_tight_vzReweight_jetPtReweight_jetTrkMaxFilter_removeHYDJETjet0p45_2025-8-21.root");
  TFile *fT = TFile::Open("/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/scanningOutput/pp/latest/jetTrkMax/pp_jetTrkMax_trkPt-14_jet60_updatedJEC_2025-11-24.root");
  const double win[][2] = {{80,100},{100,120},{120,150},{150,200},{200,300}};
  printf("light = data leading-track ptRel (pp jet60, trk pT > 14); c = PYTHIA (smeared) at its MC share; b = PYTHIA (smeared)\n window  truth | nominal (PYTHIA l+c) | data-driven light: 0-3 / 0-5 GeV | <ptRel> 0-3: data-light  PYTHIA-light\n");
  for(int j = 1; j <= 5; j++){
    const double lo = win[j-1][0], hi = win[j-1][1];
    TH1D *d = range05(projWin(fD, hName, lo, hi), "d"); Tpl t = buildTemplates(fM, lo, hi); TH1D *b = bWithGS(t, t.fGSnat);
    TH1D *tr = (TH1D*) fT->Get(Form("h_jetTrkMaxPtRel_J%d", j)); TH1D *L = range05(tr, Form("L%d", j)); L->Scale(1./L->Integral());
    // MC c share of non-b (tagged, the templates' own counts)
    TH1D *pc = projWin(fM, Form("%s_cJets_T0", hName), lo, hi); double nl = 0; for(const char *f : {"u","d","s","g"}){ TH1D *p = projWin(fM, Form("%s_%sJets_T0", hName, f), lo, hi); nl += p->Integral(); delete p; }
    const double cs = pc->Integral()/(pc->Integral() + nl); delete pc;
    std::vector<TH1D*> v; FitRes r0 = fitNominal(d, b, t, fitLo, 3.0, &v); const double p0 = purityFull(r0, v, 1, d, fitLo, 3.0);
    Tpl td = t; td.lc = (TH1D*) t.c->Clone(Form("lcdd%d", j)); td.lc->Scale(cs/td.lc->Integral()); td.lc->Add(L, 1. - cs); td.lc->Scale(1./td.lc->Integral());
    printf(" %3.0f-%-3.0f %.3f |   %.3f (chi2 %4.1f/%d)  |", lo, hi, t.truthB, p0, r0.chi2, r0.ndf);
    for(double e : {3.0, 5.0}){ std::vector<TH1D*> w; FitRes r = fitNominal(d, b, td, fitLo, e, &w); printf("  %.3f (%5.1f/%d)", purityFull(r, w, 1, d, fitLo, e), r.chi2, r.ndf); }
    // variant 2: UNSMEARED PYTHIA light x (data / PYTHIA leading-track ptRel), c and b nominal
    { Tpl tu = buildTemplates(fM, lo, hi, 1.0, 0, -1.); TH1D *Lm = range05((TH1D*) fTM->Get(Form("h_jetTrkMaxPtRel_J%d", j)), Form("Lm%d", j)); Lm->Scale(1./Lm->Integral());
      TH1D *rw = (TH1D*) L->Clone(Form("rw%d", j)); rw->Divide(Lm);
      TH1D *lw = (TH1D*) tu.l->Clone(Form("lw%d", j)); lw->Multiply(rw); lw->Scale(1./lw->Integral());
      Tpl t2 = t; t2.lc = (TH1D*) t.c->Clone(Form("lc2%d", j)); t2.lc->Scale(cs/t2.lc->Integral()); t2.lc->Add(lw, 1. - cs); t2.lc->Scale(1./t2.lc->Integral());
      std::vector<TH1D*> w; FitRes r = fitNominal(d, b, t2, fitLo, 3.0, &w); printf("  | reweighted PYTHIA light: %.3f (%4.1f/%d)", purityFull(r, w, 1, d, fitLo, 3.0), r.chi2, r.ndf); }
    auto mean3 = [](TH1D *h){ double s = 0, n = 0; for(int i = 1; i <= h->FindBin(2.999); i++){ s += h->GetBinCenter(i)*h->GetBinContent(i); n += h->GetBinContent(i); } return s/n; };
    printf("  |   %.3f   %.3f   (MC c share %.2f)\n", mean3(L), mean3(t.l), cs);
    delete d; delete b;
  }
}
