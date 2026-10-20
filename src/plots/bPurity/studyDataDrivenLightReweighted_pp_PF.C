// pp PF b purity with a data-driven light template REWEIGHTED IN TRACK pT, the ATLAS
// recipe (arXiv:2204.13530: track-jet pairs in inclusive-jet data, track pT spectra
// reweighted to the light-jet muon pT of the dijet MC). 2026-10-09, user; follows
// studyDataDrivenLight_pp_PF.C (same idea without the track-pT axis).
// Per jet-pT window:
//   T(ptRel, pT_trk)  leading-track ptRel x track pT of pp jet60 data,
//                     h_jetTrkMaxPtRel_trkPt_recoJetPt of pp_jetTrkMax_scan.C (2026-10-9
//                     version, the legacy HighEGJet skims; binning ptRelMuPt3D.h)
//   M_L(pT_mu)        PYTHIA light (u, d, s, g) tagged-jet muon pT, from the 2026-10-09
//                     PYTHIA scan's h_muptrel_mupt_recoJetPt_inclRecoMuonTag_triggerOn_<fl>Jets
//   light template    sum over track-pT bins of T(ptRel, p) x M_L(p) / T(p) (shapes)
// c and b: the nominal PYTHIA templates (per-window smearing), c at its MC share of
// non-b; fit 0-3 GeV (nominal edge). Printed beside: the nominal fit, the unweighted
// track template, PYTHIA truth, and <pT> of tracks / target to show the reweighting.
// tracksFile = "closure": the PYTHIA light muons themselves stand in for the tracks
// (unsmeared), so the reweighted template must equal PYTHIA's unsmeared light.
// Finding (2026-10-19, full scan, 3.38 M jets): fitted purity vs PYTHIA truth at
// 80-100 / 100-120 / 120-150 / 150-200 / 200-300 GeV: nominal +17/+11/+22/+34/+58%,
// tracks unweighted +10/-4/+2/+6/+3%, reweighted +11/-2/+5/+12/+20%. The data
// track template removes most of the excess; the track-pT reweighting pulls the
// purity back up where the tracks are much harder than the light muons (49 vs 28
// GeV at 200-300). Closure passes (reweighted <ptRel> = PYTHIA light).
// Usage, from src/plots/bPurity/:
//   root -l -b -q 'studyDataDrivenLightReweighted_pp_PF.C("<pp_jetTrkMax output>.root")'
//   root -l -b -q 'studyDataDrivenLightReweighted_pp_PF.C("closure")'
#include "RooMsgService.h"
#include "TH3D.h"
#include "../../../headers/functions/caloJetBPurityFit_pp.h"

namespace ddr {
const char *mc3DFile = "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PYTHIA/PYTHIA_DiJet_manualJEC_pThat-15_mu12_pTmu-15to999_tight_vzReweight_jetTrkMaxFilter_removeHYDJETjet0p35CutOnGen_WDecayFilter_axisSmearScan_genMuTagFix_2026-10-9.root";
const char *defaultTracks = "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/pp/pp_jetTrkMax_trkPtAxis_jet60_ppRefJEC_2026-10-9.root";
int uid = 0;
// ptRel x pT (yx) of a TH3 over a jet-pT window
TH2D* proj2(TH3 *h, double lo, double hi){
  h->GetZaxis()->SetRange(h->GetZaxis()->FindBin(lo + 1e-6), h->GetZaxis()->FindBin(hi - 1e-6));
  TH2D *q = (TH2D*) h->Project3D("yx"); TH2D *p = (TH2D*) q->Clone(Form("ddr2_%d", uid++)); p->SetDirectory(nullptr); delete q; h->GetZaxis()->SetRange(0, 0);
  return p;
}
TH2D* lightMu2(TFile *f, double lo, double hi){
  TH2D *s = nullptr;
  for(const char *fl : {"u", "d", "s", "g"}){ TH3 *h = nullptr; f->GetObject(Form("h_muptrel_mupt_recoJetPt_inclRecoMuonTag_triggerOn_%sJets", fl), h); if(!h) return nullptr;
    TH2D *p = proj2(h, lo, hi); if(!s) s = p; else { s->Add(p); delete p; } }
  return s;
}
double meanY(TH2D *h){ TH1D *p = h->ProjectionY(Form("my%d", uid++)); double m = p->GetMean(); delete p; return m; }
// reweight the pT (y) distribution of T to the target pT spectrum; returns ptRel 1D, unit area
TH1D* reweighted(TH2D *T, TH1D *target){
  TH1D *tp = T->ProjectionY(Form("tp%d", uid++)); tp->Scale(1./tp->Integral(1, tp->GetNbinsX())); TH1D *mt = (TH1D*) target->Clone(Form("mt%d", uid++)); mt->Scale(1./mt->Integral(1, mt->GetNbinsX()));
  TH1D *L = T->ProjectionX(Form("L%d", uid++)); L->Reset();
  for(int j = 1; j <= T->GetNbinsY(); j++){ const double w = tp->GetBinContent(j) > 0. ? mt->GetBinContent(j)/tp->GetBinContent(j) : 0.;
    for(int i = 1; i <= T->GetNbinsX(); i++){ L->SetBinContent(i, L->GetBinContent(i) + w*T->GetBinContent(i, j)); L->SetBinError(i, hypot(L->GetBinError(i), w*T->GetBinError(i, j))); } }
  L->Scale(1./L->Integral()); delete tp; delete mt; return L;
}
double mean3(TH1D *h){ double s = 0, n = 0; for(int i = 1; i <= h->FindBin(2.999); i++){ s += h->GetBinCenter(i)*h->GetBinContent(i); n += h->GetBinContent(i); } return n > 0 ? s/n : 0; }
} // namespace ddr
using namespace ddr;

void studyDataDrivenLightReweighted_pp_PF(const char *tracksFile = ""){
  setJetCollection("PFJets"); RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
  const bool closure = TString(tracksFile) == "closure";
  TFile *fD = TFile::Open(dataPath), *fM = TFile::Open(mcPath), *f3 = TFile::Open(mc3DFile);
  TFile *fT = closure ? nullptr : TFile::Open(tracksFile[0] ? tracksFile : defaultTracks);
  if(!fD || !fM || !f3 || (!closure && !fT)){ printf("ERROR: inputs (tracks file: %s)\n", closure ? "closure" : (tracksFile[0] ? tracksFile : defaultTracks)); return; }
  TH3 *t3 = nullptr; if(!closure){ fT->GetObject("h_jetTrkMaxPtRel_trkPt_recoJetPt", t3); if(!t3){ printf("ERROR: h_jetTrkMaxPtRel_trkPt_recoJetPt missing -- the 2026-10-09 pp_jetTrkMax_scan.C output is needed\n"); return; } }
  const double win[][2] = {{80,100},{100,120},{120,150},{150,200},{200,300}};
  printf("%s\n window  truth | nominal        | tracks unweighted | tracks reweighted | <pT>: tracks  target(light mu) | <ptRel 0-3>: reweighted  PYTHIA light (unsmeared)\n",
         closure ? "CLOSURE: PYTHIA light muons as the 'tracks' (the reweighted light must equal PYTHIA's unsmeared light)" : "light = data leading tracks reweighted in track pT to PYTHIA light-muon pT; c, b = PYTHIA");
  for(auto &w : win){
    const double lo = w[0], hi = w[1];
    TH1D *d = range05(projWin(fD, hName, lo, hi), Form("d%d", uid++)); Tpl t = buildTemplates(fM, lo, hi); TH1D *b = bWithGS(t, t.fGSnat);
    TH2D *Lmu = lightMu2(f3, lo, hi); if(!Lmu){ printf("ERROR: PYTHIA light 3D missing\n"); return; }
    TH2D *T = closure ? (TH2D*) Lmu->Clone(Form("Tc%d", uid++)) : proj2(t3, lo, hi);
    TH1D *target = Lmu->ProjectionY(Form("tg%d", uid++));
    TH1D *Lrw = range05(reweighted(T, target), Form("Lrw%d", uid++)); Lrw->Scale(1./Lrw->Integral());
    TH1D *Tx = T->ProjectionX(Form("Lun%d", uid++)); TH1D *Lun = range05(Tx, Form("Lu%d", uid++)); Lun->Scale(1./Lun->Integral()); delete Tx;
    // MC c share of non-b
    TH1D *pc = projWin(fM, Form("%s_cJets_T0", hName), lo, hi); double nl = 0; for(const char *f : {"u","d","s","g"}){ TH1D *p = projWin(fM, Form("%s_%sJets_T0", hName, f), lo, hi); nl += p->Integral(); delete p; }
    const double cs = pc->Integral()/(pc->Integral() + nl); delete pc;
    auto fitWith = [&](TH1D *L, double &chi2, int &ndf){ Tpl tx = t; tx.lc = (TH1D*) t.c->Clone(Form("lcx%d", uid++)); tx.lc->Scale(cs/tx.lc->Integral()); tx.lc->Add(L, 1. - cs); tx.lc->Scale(1./tx.lc->Integral());
      std::vector<TH1D*> v; FitRes r = fitNominal(d, b, tx, fitLo, lightTailEdge(t), &v); chi2 = r.chi2; ndf = r.ndf; return purityFull(r, v, 1, d, fitLo, lightTailEdge(t)); };
    std::vector<TH1D*> v0; FitRes r0 = fitNominal(d, b, t, fitLo, lightTailEdge(t), &v0); const double p0 = purityFull(r0, v0, 1, d, fitLo, lightTailEdge(t));
    double cu, cr; int nu, nr; const double pu = fitWith(Lun, cu, nu), pr = fitWith(Lrw, cr, nr);
    Tpl tu = buildTemplates(fM, lo, hi, 1.0, 0, -1.);
    printf(" %3.0f-%-3.0f %.3f | %.3f (%4.1f/%d) | %.3f (%4.1f/%d)  | %.3f (%4.1f/%d)  |   %5.1f  %5.1f          |   %.3f   %.3f\n", lo, hi, t.truthB, p0, r0.chi2, r0.ndf, pu, cu, nu, pr, cr, nr, meanY(T), target->GetMean(), mean3(Lrw), mean3(tu.l));
    delete d; delete b; delete Lmu; delete T; delete target;
  }
}
