// PbPb PF b purity with a DATA-DRIVEN light template from ALL tracks, ATLAS-style
// (arXiv:2204.13530: every track-jet pair in inclusive-jet data, track pT reweighted
// to the light-jet muon pT). 2026-10-19, user. The PbPb counterpart of
// studyDataDrivenLightReweighted_pp_PF.C, with all tracks instead of the leading one:
// requiring a track to be the hardest in its jet selects broad jets, so the
// leading-track ptRel exceeds the light-muon ptRel even in PYTHIA (+0.07 to +0.22 GeV
// in pp, rising with jet pT) and no track-pT reweighting removes it. All tracks carry
// no such condition.
// Per coarse class (data summed over its 5% slices) and jet-pT window:
//   T(ptRel, p)    every charged track (pT > 15 GeV, dR < 0.4) of every akCs4PF jet in
//                  PbPb MinBias Part1, h_jetTrkPtRel_trkPt_recoJetPt_C<i> of
//                  PbPb_jetTrkMax_scan.C (charged PF candidates; the forest has no
//                  track tree). The jet axis is the one the muon ptRel uses, so the
//                  underlying event's axis smearing is in the template.
//   T - P          the same minus the leading track of a perpendicular cone per jet
//                  (h_perpTrkMaxPtRel_..._C<i>, same per-jet normalization): the
//                  underlying event's own tracks above 15 GeV, removed
//   M_L(p)         PYTHIA+HYDJET light (u, d, s, g) tagged-jet muon pT of the class,
//                  the 2026-10-19 DiJet pThat-30 scan's
//                  h_muptrel_mupt_recoJetPt_inclRecoMuonTag_triggerOn_<fl>Jets_C<i>;
//                  target0to80 (default) takes it from 0-80% for every class: per
//                  class its effective entries are 9-17 in 0-10% at 80-120 GeV, and
//                  its <pT> differs between classes by less than its noise
//   light template sum over p of T(ptRel, p) M_L(p) / T(p), unit area
// Each replaces the nominal light template (axis-smeared pp light, kLightSmearedPP);
// b, c, the c/light ratio and the fake-ptRel component stay nominal
// (caloJetBPurityFit_PbPb.h), fit 0 to fitEdge (default 5 GeV, as the pp study).
// Printed per window: PYTHIA+HYDJET truth, nominal, all tracks (T, T - P), the
// leading track (for comparison), all reweighted; chi2/ndf; <pT> tracks vs target and
// the target's effective entries (the light-muon statistics are thin per class).
// Figure: figures/bPurity/PFJets_PbPb/dataDrivenLight/bPurity_allTrkLight_fit0to<edge>[_classTarget].pdf
// Usage, from src/plots/bPurity/:  root -l -b -q 'studyDataDrivenLightAllTrk_PF_PbPb.C(5.0, true)'
#include "RooMsgService.h"
#include "TH3D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "../../../headers/plotting/plotStyle.h"
#include "../../../headers/plotting/ratioPanel.h"
#include "../../../headers/functions/caloJetBPurityFit_PbPb.h"

namespace dla {
const char *trkFile = "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/PbPb_MinBias_Part1_caloFlowPFTowers_leadTrkPtRel_trkPtAxis_ultraFineCentBins_2026-10-19.root";
const char *phLightFile = "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PYTHIAHYDJET/PYTHIAHYDJET_manualJEC_DiJet_pThat-30_mu12_pTmu-15to999_tight_vzReweight_hiBinShift-10_jetTrkMaxFilter_doPThatCorrelationFilterTight_removeUnmatchedJets_WDecayFilter_genMuTagFix_2026-10-19_ultraFineCentBins.root";
const char *outDir = "../../../figures/bPurity/PFJets_PbPb/dataDrivenLight";
const double win[][2] = {{80,100},{100,120},{120,150},{150,200},{200,300}};
const int nWin = 5;
int uid = 0;

// slices lo..hi of name_C<i>, summed
TH3* sumSlices(TFile *f, const char *name, int lo, int hi){
  TH3 *s = nullptr;
  for(int i = lo; i <= hi; i++){ TH3 *h = nullptr; f->GetObject(Form("%s_C%d", name, i), h); if(!h) return nullptr;
    if(!s){ s = (TH3*) h->Clone(Form("dla3_%d", uid++)); s->SetDirectory(nullptr); } else s->Add(h); }
  return s;
}
// ptRel x pT (yx) over a jet-pT window
TH2D* proj2(TH3 *h, double lo, double hi){
  h->GetZaxis()->SetRange(h->GetZaxis()->FindBin(lo + 1e-6), h->GetZaxis()->FindBin(hi - 1e-6));
  TH2D *q = (TH2D*) h->Project3D("yx"); TH2D *p = (TH2D*) q->Clone(Form("dla2_%d", uid++)); p->SetDirectory(nullptr); delete q; h->GetZaxis()->SetRange(0, 0);
  return p;
}
// as ddr::reweighted (studyDataDrivenLightReweighted_pp_PF.C): T's pT (y) spectrum
// reweighted to the target's; returns the ptRel 1D, unit area
TH1D* reweighted(TH2D *T, TH1D *target){
  TH1D *tp = T->ProjectionY(Form("tp%d", uid++)); tp->Scale(1./tp->Integral(1, tp->GetNbinsX())); TH1D *mt = (TH1D*) target->Clone(Form("mt%d", uid++)); mt->Scale(1./mt->Integral(1, mt->GetNbinsX()));
  TH1D *L = T->ProjectionX(Form("L%d", uid++)); L->Reset();
  for(int j = 1; j <= T->GetNbinsY(); j++){ const double w = tp->GetBinContent(j) > 0. ? mt->GetBinContent(j)/tp->GetBinContent(j) : 0.;
    for(int i = 1; i <= T->GetNbinsX(); i++){ L->SetBinContent(i, L->GetBinContent(i) + w*T->GetBinContent(i, j)); L->SetBinError(i, hypot(L->GetBinError(i), w*T->GetBinError(i, j))); } }
  L->Scale(1./L->Integral()); delete tp; delete mt; return L;
}
double neff(TH1 *h){ double s = 0, s2 = 0; for(int i = 1; i <= h->GetNbinsX(); i++){ s += h->GetBinContent(i); s2 += pow(h->GetBinError(i), 2); } return s2 > 0 ? s*s/s2 : 0; }
double mean3(TH1D *h, double hi = 3.){ double s = 0, n = 0; for(int i = 1; i <= h->FindBin(hi - 1e-6); i++){ s += h->GetBinCenter(i)*h->GetBinContent(i); n += h->GetBinContent(i); } return n > 0 ? s/n : 0; }
}
using namespace dla;

void studyDataDrivenLightAllTrk_PF_PbPb(double fitEdge = 5.0, bool target0to80 = true){
  setJetCollection("PFJets"); RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL); initPlotStyle();
  fitEdgeFixed = fitEdge;   // kEdgeFixed: every fit below runs 0 to fitEdge
  TFile *fD = TFile::Open(dataPathPbPb), *fT = TFile::Open(trkFile), *fL = TFile::Open(phLightFile); PHFiles ph;
  if(!fD || !fT || !fL || !openPH(ph)){ printf("ERROR: inputs\n"); return; }
  gSystem->mkdir(outDir, kTRUE);
  printf("fit range 0-%.1f GeV; light template: data tracks reweighted to PYTHIA+HYDJET light-muon pT of %s\n", lightTailEdge(Tpl()), target0to80 ? "0-80% (every class)" : "each class");

  double P[NCls][nWin][5], E[NCls][nWin][5], truth[NCls][nWin]; bool ok[NCls][nWin];   // 0 nominal, 1 all, 2 all - perp, 3 leading
  for(int c = 0; c < nCls; c++){
    const int s0 = clsSlice[c][0], s1 = clsSlice[c][1];
    TH3 *tAll = sumSlices(fT, "h_jetTrkPtRel_trkPt_recoJetPt", s0, s1), *tPerp = sumSlices(fT, "h_perpTrkMaxPtRel_trkPt_recoJetPt", s0, s1);
    TH3 *tLead = sumSlices(fT, "h_jetTrkMaxPtRel_trkPt_recoJetPt", s0, s1), *mL = nullptr;
    const int m0 = target0to80 ? 1 : s0, m1 = target0to80 ? 16 : s1;   // light-muon target: 0-80% or the class
    for(const char *fl : {"u", "d", "s", "g"}){ TH3 *h = sumSlices(fL, Form("h_muptrel_mupt_recoJetPt_inclRecoMuonTag_triggerOn_%sJets", fl), m0, m1); if(!h){ mL = nullptr; break; } if(!mL) mL = h; else { mL->Add(h); delete h; } }
    if(!tAll || !tPerp || !tLead || !mL){ printf("ERROR: %s: a track or light-muon histogram is missing\n", clsName[c]); return; }
    printf("\n%s   truth | nominal       | all tracks    | all - perp    | leading       | tracks/mu <pT>  mu neff | <ptRel 0-3> all, nominal light\n", clsName[c]);
    for(int k = 0; k < nWin; k++){
      const double lo = win[k][0], hi = win[k][1]; ok[c][k] = false;
      TH1D *d = dataWin(fD, c, lo, hi); if(!d || !fittable(d)) continue;
      double lLo, lHi; Tpl t = buildTemplatesPbPb(ph, c, lo, hi, 1.0, 0, lightNeffMinNominal, &lLo, &lHi);
      const double e = lightTailEdge(t); TH1D *b = bWithGS(t, TMath::Min(1., t.fGSnat + bGSShiftNominal()));
      // c per light in DiJet over the light template's range, as buildTemplatesPbPb
      TFile *fDJ = phFile(ph, kDiJet, lo);
      TH1D *cW = tplOneSample(fDJ, {"c"}, c, lLo, lHi, 0), *lW = tplOneSample(fDJ, {"u", "d", "s", "g"}, c, lLo, lHi, 0);
      const double rCL = lW->Integral() > 0. ? cW->Integral()/lW->Integral() : 0.; delete cW; delete lW;
      TH2D *Mmu = proj2(mL, lo, hi); TH1D *target = Mmu->ProjectionY(Form("tg%d", uid++));
      auto fitWith = [&](TH1D *L, double &er, double &chi2, int &ndf){
        Tpl tx = t; TH1D *l = range05(L, Form("l%d", uid++)); l->Scale(1./l->Integral()); tx.l = l;
        tx.lc = (TH1D*) l->Clone(Form("lc%d", uid++)); tx.lc->Add(t.c, rCL); tx.lc->Scale(1./tx.lc->Integral());
        std::vector<TH1D*> v; FitRes r = fitNominal(d, b, tx, fitLo, e, &v); chi2 = r.chi2; ndf = r.ndf;
        const double p = purityFull(r, v, 1, d, fitLo, e, &er); er *= fakeKeepFrac(); return p*fakeKeepFrac(); };
      double chi[5]; int nd[5];
      { std::vector<TH1D*> v; FitRes r = fitNominal(d, b, t, fitLo, e, &v); chi[0] = r.chi2; nd[0] = r.ndf; P[c][k][0] = purityFull(r, v, 1, d, fitLo, e, &E[c][k][0])*fakeKeepFrac(); E[c][k][0] *= fakeKeepFrac(); }
      TH2D *Ta = proj2(tAll, lo, hi), *Tp = proj2(tPerp, lo, hi), *Tl = proj2(tLead, lo, hi);
      TH2D *Ts = (TH2D*) Ta->Clone(Form("Ts%d", uid++)); Ts->Add(Tp, -1.);
      for(int i = 0; i <= Ts->GetNbinsX()+1; i++) for(int j = 0; j <= Ts->GetNbinsY()+1; j++) if(Ts->GetBinContent(i, j) < 0.) Ts->SetBinContent(i, j, 0.);
      TH1D *La = reweighted(Ta, target), *Ls = reweighted(Ts, target), *Ll = reweighted(Tl, target);
      P[c][k][1] = fitWith(La, E[c][k][1], chi[1], nd[1]); P[c][k][2] = fitWith(Ls, E[c][k][2], chi[2], nd[2]); P[c][k][3] = fitWith(Ll, E[c][k][3], chi[3], nd[3]);
      truth[c][k] = t.truthB; ok[c][k] = true;
      TH1D *tpy = Ta->ProjectionY(Form("tpy%d", uid++));
      printf(" %3.0f-%-3.0f  %.3f | %.3f (%4.1f/%d) | %.3f (%4.1f/%d) | %.3f (%4.1f/%d) | %.3f (%4.1f/%d) |  %5.1f / %5.1f  %6.0f |  %.3f  %.3f\n", lo, hi, t.truthB,
             P[c][k][0], chi[0], nd[0], P[c][k][1], chi[1], nd[1], P[c][k][2], chi[2], nd[2], P[c][k][3], chi[3], nd[3],
             tpy->GetMean(), target->GetMean(), neff(target), mean3(La), mean3(t.l));
      delete tpy; delete d; delete b; delete Mmu; delete target; delete Ta; delete Tp; delete Tl; delete Ts; delete La; delete Ls; delete Ll;
    }
  }

  // figure: purity vs jet pT per class, fit / truth below
  TCanvas *cv = new TCanvas("cv", "", 1600, 560);
  const double x0 = 0.05, wP = (1. - x0 - 0.01)/nCls, split = 0.34;
  for(int c = 0; c < nCls; c++){
    double x[nWin], ex[nWin], ez[nWin], y[4][nWin], ey[4][nWin], r[4][nWin], er[4][nWin]; int n = 0;
    for(int k = 0; k < nWin; k++){ if(!ok[c][k]) continue; x[n] = 0.5*(win[k][0] + win[k][1]); ex[n] = 0.5*(win[k][1] - win[k][0]); ez[n] = 0.;
      for(int v = 0; v < 4; v++){ const double val = v == 0 ? truth[c][k] : P[c][k][v-1], err = v == 0 ? 0. : E[c][k][v-1];
        y[v][n] = val; ey[v][n] = err; r[v][n] = val/truth[c][k]; er[v][n] = err/truth[c][k]; } n++; }
    const char *hex[4] = {hexData, hexMC, hexCorrected, okabeHex[3]}; const int mk[4] = {1, markOpenSquare, markFilledCircle, markOpenDiamond};
    const double xl = x0 + c*wP - (c == 0 ? x0 : 0.), xh = x0 + (c+1)*wP;
    for(int pane = 0; pane < 2; pane++){
      cv->cd(); TPad *p = new TPad(Form("p%d%d", c, pane), "", xl, pane == 0 ? split : 0., xh, pane == 0 ? 0.89 : split);
      p->SetLeftMargin(c == 0 ? x0/(xh - xl) + 0.12 : 0.02); p->SetRightMargin(0.02);
      if(pane == 0){ p->SetTopMargin(0.09); p->SetBottomMargin(0.015); } else { p->SetTopMargin(0.02); p->SetBottomMargin(0.30); }
      p->Draw(); p->cd();
      TH1D *fr = new TH1D(Form("fr%d%d", c, pane), "", 1, 80, 300); fr->SetStats(0);
      if(pane == 0){ fr->SetMinimum(0.2); fr->SetMaximum(0.9); fr->GetYaxis()->SetTitle(c == 0 ? "b purity" : ""); fr->GetYaxis()->SetTitleSize(0.075); fr->GetYaxis()->SetLabelSize(c == 0 ? 0.065 : 0.); fr->GetXaxis()->SetLabelSize(0.); }
      else { fr->SetMinimum(0.5); fr->SetMaximum(1.8); styleRatioAxes(fr, "p_{T}^{jet} [GeV]", c == 0 ? "fit / truth" : ""); fr->GetYaxis()->SetLabelSize(c == 0 ? 0.085 : 0.); }
      fr->Draw("axis"); if(pane == 1) unityLine(80, 300)->Draw();
      for(int v = (pane == 0 ? 0 : 1); v < 3; v++){   // all - perp equals all to < 0.002 (table): not drawn
        TGraphErrors *g = new TGraphErrors(n, x, pane == 0 ? y[v] : r[v], ex, pane == 0 ? ey[v] : er[v]); int col = TColor::GetColor(hex[v]);
        g->SetLineColor(col); g->SetMarkerColor(col); g->SetMarkerStyle(mk[v]); g->SetMarkerSize(v == 0 ? 0 : 1.2); g->SetLineWidth(v == 0 ? 3 : 2);
        g->Draw(v == 0 ? "z same" : "pz same");
      }
      if(pane == 0){ TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.068); t.DrawLatex(c == 0 ? 0.45 : 0.08, 0.83, clsName[c]); }
    }
  }
  cv->cd(); TLegend *lg = makeLegend(0.30, 0.90, 0.99, 0.99, 0.032); lg->SetNColumns(3);
  { TH1D *dm[4]; const char *hex[4] = {hexData, hexMC, hexCorrected, okabeHex[3]}; const int mk[4] = {1, markOpenSquare, markFilledCircle, markOpenDiamond};
    const char *lab[4] = {"PYTHIA+HYDJET truth", "fit, nominal light", "fit, data all tracks, reweighted", ""};
    for(int v = 0; v < 3; v++){ dm[v] = new TH1D(Form("dm%d", v), "", 1, 0, 1); styleH(dm[v], hex[v], mk[v]); if(v == 0){ dm[v]->SetLineWidth(3); } lg->AddEntry(dm[v], lab[v], v == 0 ? "l" : "p"); } }
  lg->Draw();
  TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.032); t.DrawLatex(0.01, 0.94, Form("PbPb 5.02 TeV, anti-k_{T} R = 0.4 PF jets; fit 0-%g GeV", fitEdge));
  savePdfTight(cv, Form("%s/bPurity_allTrkLight_fit0to%g%s.pdf", outDir, fitEdge, target0to80 ? "" : "_classTarget"));
}
