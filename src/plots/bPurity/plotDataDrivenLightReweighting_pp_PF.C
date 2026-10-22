// Figures for the pp PF data-driven light template and its track-pT reweighting
// (ATLAS recipe, arXiv:2204.13530), the derivation in studyDataDrivenLightReweighted_pp_PF.C
// step by step. 2026-10-19, user. Reuses that macro's helpers (included below), so
// the figures show exactly what the study computes; the purities printed here must
// match its table.
// Per jet-pT window, 80-100 ... 200-300 GeV:
//   1. inputs + weight: the pT spectrum of the data leading tracks T(p) (pp jet60,
//      pp_jetTrkMax_scan.C) against the target, the PYTHIA light (u, d, s, g)
//      tagged-jet muon pT M_L(p); both unit area, per GeV (log y: steeply falling,
//      all positive). Below, the weight w(p) = M_L(p) / T(p).
//   2. templates: light ptRel from the tracks unweighted and reweighted
//      (sum over p of T(ptRel, p) w(p)), against PYTHIA's light template, nominal
//      (per-window axis smearing, what the fit uses) and unsmeared; in the fit's own
//      binning over its range (range05: 0 to fitHi), unit area, per GeV. Below, the
//      ratio to PYTHIA nominal.
//   3. result: fitted b purity vs jet pT with the light template replaced
//      (unweighted, reweighted) and nominal, against PYTHIA truth; below, / truth.
// c and b stay PYTHIA in every fit (c at its MC share of non-b), fit 0 to fitEdge
// (argument, default 5 GeV, 2026-10-19 user; the nominal edge is 3), in the file names.
// Output: figures/bPurity/PFJets_pp/dataDrivenLight/
//   dataDrivenLight_trackPtWeight.pdf, dataDrivenLight_ptRelTemplates.pdf (fit-independent),
//   dataDrivenLight_purity_fit0to<edge>.pdf
// Usage, from src/plots/bPurity/:  root -l -b -q 'plotDataDrivenLightReweighting_pp_PF.C(5.0)'
#include "studyDataDrivenLightReweighted_pp_PF.C"
#include "TCanvas.h"
#include "TPad.h"
#include "TLatex.h"
#include "TLine.h"
#include "TGraphErrors.h"
#include "TGaxis.h"
#include "TColor.h"
#include "TLegend.h"
#include "TStyle.h"
#include "../../../headers/plotting/plotStyle.h"
#include "../../../headers/plotting/ratioPanel.h"

namespace pdl {
const double win[][2] = {{80,100},{100,120},{120,150},{150,200},{200,300}};
const int nWin = 5;
const char *outDir = "../../../figures/bPurity/PFJets_pp/dataDrivenLight";

// unit area, then per GeV (after the integral, which is a count)
void density(TH1D *h){ const double s = h->Integral(1, h->GetNbinsX()); if(s > 0) h->Scale(1./s); h->Scale(1., "width"); }

// one panel of a 1 x nWin row, top (frac of the canvas height above split) or bottom
TPad* cell(TCanvas *c, int k, bool top, double split){
  const double x0 = 0.04, w = (1. - x0 - 0.01)/nWin;
  const double xl = x0 + k*w - (k == 0 ? x0 : 0.), xh = x0 + (k+1)*w;
  TPad *p = new TPad(Form("p%d%d", k, top), "", xl, top ? split : 0., xh, top ? 1. : split);
  p->SetLeftMargin(k == 0 ? x0/(xh - xl) + 0.10 : 0.02); p->SetRightMargin(0.02);
  if(top){ p->SetTopMargin(0.09); p->SetBottomMargin(0.015); } else { p->SetTopMargin(0.02); p->SetBottomMargin(0.30); }
  c->cd(); p->Draw(); p->cd(); return p;
}
void label(double x, double y, const char *s, double size = 0.065, int align = 13){
  TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(size); t.SetTextAlign(align); t.DrawLatex(x, y, s);
}
}
using namespace pdl;

void plotDataDrivenLightReweighting_pp_PF(double fitEdge = 5.0){
  setJetCollection("PFJets"); RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
  fitEdgeFixed = fitEdge;   // kEdgeFixed: every fit below runs 0 to fitEdge
  printf("fit range 0-%.1f GeV\n", lightTailEdge(Tpl()));
  initPlotStyle(); gStyle->SetEndErrorSize(0);
  gSystem->mkdir(outDir, kTRUE);
  TFile *fD = TFile::Open(dataPath), *fM = TFile::Open(mcPath), *f3 = TFile::Open(mc3DFile), *fT = TFile::Open(defaultTracks);
  if(!fD || !fM || !f3 || !fT){ printf("ERROR: inputs\n"); return; }
  TH3 *t3 = nullptr; fT->GetObject("h_jetTrkMaxPtRel_trkPt_recoJetPt", t3);
  if(!t3){ printf("ERROR: h_jetTrkMaxPtRel_trkPt_recoJetPt missing in %s\n", defaultTracks); return; }

  TH1D *tPt[nWin], *mPt[nWin], *wPt[nWin], *lUn[nWin], *lRw[nWin], *lMC[nWin], *lMC0[nWin];
  double tMean[nWin], mMean[nWin];   // <pT>, from the counts (before the density conversion)
  double truth[nWin], pNom[nWin], pUn[nWin], pRw[nWin];
  printf(" window  truth  nominal  unweighted  reweighted\n");
  for(int k = 0; k < nWin; k++){
    const double lo = win[k][0], hi = win[k][1];
    TH1D *d = range05(projWin(fD, hName, lo, hi), Form("d%d", uid++)); Tpl t = buildTemplates(fM, lo, hi); TH1D *b = bWithGS(t, t.fGSnat);
    TH2D *Lmu = lightMu2(f3, lo, hi), *T = proj2(t3, lo, hi);
    if(!Lmu){ printf("ERROR: PYTHIA light 3D missing\n"); return; }
    TH1D *target = Lmu->ProjectionY(Form("tg%d", uid++));
    TH1D *Lrw = range05(reweighted(T, target), Form("Lrw%d", uid++)); Lrw->Scale(1./Lrw->Integral());
    TH1D *Tx = T->ProjectionX(Form("Lun%d", uid++)); TH1D *Lun = range05(Tx, Form("Lu%d", uid++)); Lun->Scale(1./Lun->Integral()); delete Tx;
    // the study's fit, unchanged: c at its MC share of non-b, mixed into the light template
    TH1D *pc = projWin(fM, Form("%s_cJets_T0", hName), lo, hi); double nl = 0; for(const char *f : {"u","d","s","g"}){ TH1D *p = projWin(fM, Form("%s_%sJets_T0", hName, f), lo, hi); nl += p->Integral(); delete p; }
    const double cs = pc->Integral()/(pc->Integral() + nl); delete pc;
    auto fitWith = [&](TH1D *L){ Tpl tx = t; tx.lc = (TH1D*) t.c->Clone(Form("lcx%d", uid++)); tx.lc->Scale(cs/tx.lc->Integral()); tx.lc->Add(L, 1. - cs); tx.lc->Scale(1./tx.lc->Integral());
      std::vector<TH1D*> v; FitRes r = fitNominal(d, b, tx, fitLo, lightTailEdge(t), &v); return purityFull(r, v, 1, d, fitLo, lightTailEdge(t)); };
    std::vector<TH1D*> v0; FitRes r0 = fitNominal(d, b, t, fitLo, lightTailEdge(t), &v0);
    truth[k] = t.truthB; pNom[k] = purityFull(r0, v0, 1, d, fitLo, lightTailEdge(t)); pUn[k] = fitWith(Lun); pRw[k] = fitWith(Lrw);
    printf(" %3.0f-%-3.0f %.3f  %.3f    %.3f       %.3f\n", lo, hi, truth[k], pNom[k], pUn[k], pRw[k]);

    // figure 1: pT spectra and weight
    tPt[k] = T->ProjectionY(Form("tPt%d", k)); tPt[k]->SetDirectory(nullptr);
    mPt[k] = (TH1D*) target->Clone(Form("mPt%d", k)); mPt[k]->SetDirectory(nullptr);
    { const double st = tPt[k]->Integral(1, tPt[k]->GetNbinsX()), sm = mPt[k]->Integral(1, mPt[k]->GetNbinsX());
      TH1D *tn = (TH1D*) tPt[k]->Clone(Form("tn%d", k)), *mn = (TH1D*) mPt[k]->Clone(Form("mn%d", k)); tn->Scale(1./st); mn->Scale(1./sm);
      wPt[k] = makeRatio(mn, tn, Form("wPt%d", k), RatioErr::kBoth); delete tn; delete mn; }   // per-bin weight = reweighted()'s w
    tMean[k] = tPt[k]->GetMean(); mMean[k] = mPt[k]->GetMean();
    density(tPt[k]); density(mPt[k]);
    // figure 2: templates
    TH1D *mc = (TH1D*) t.l->Clone(Form("lmc%d", k)); Tpl tu = buildTemplates(fM, lo, hi, 1.0, 0, -1.); TH1D *mc0 = (TH1D*) tu.l->Clone(Form("lmc0%d", k));
    // in the fit's own binning (range05: 0 to fitHi), what the fit sees
    lUn[k] = (TH1D*) Lun->Clone(Form("lUn%d", k)); lRw[k] = (TH1D*) Lrw->Clone(Form("lRw%d", k));
    lMC[k] = range05(mc, Form("lMC%d", k)); lMC0[k] = range05(mc0, Form("lMC0%d", k));
    for(TH1D *h : {lUn[k], lRw[k]}) h->SetDirectory(nullptr);
    for(TH1D *h : {lUn[k], lRw[k], lMC[k], lMC0[k]}) density(h);
    delete mc; delete mc0; delete d; delete b; delete Lmu; delete T; delete target; delete Lun; delete Lrw;
  }

  const char *hdr = "pp 5.02 TeV, anti-k_{T} R = 0.4 PF jets";
  const double split = 0.36;

  // ---------------- figure 1: track pT vs light-muon pT, and the weight
  {
    TCanvas *c = new TCanvas("cPt", "", 1700, 620);
    TLegend *lg = nullptr;
    for(int k = 0; k < nWin; k++){
      TPad *p = cell(c, k, true, split); p->SetLogy(); p->SetLogx();
      styleH(tPt[k], hexData, markOpenCircle, 0.8); styleH(mPt[k], hexMC, markFilledSquare, 0.8);
      TH1D *fr = (TH1D*) tPt[k]->Clone(Form("fr1_%d", k)); fr->Reset(); fr->SetDirectory(nullptr);
      fr->GetXaxis()->SetRangeUser(15, 150); fr->SetMinimum(2e-5); fr->SetMaximum(0.5);
      fr->GetYaxis()->SetTitle(k == 0 ? "(1/N) dN/dp_{T} [GeV^{-1}]" : ""); fr->GetYaxis()->SetTitleSize(0.075); fr->GetYaxis()->SetTitleOffset(1.25);
      fr->GetYaxis()->SetLabelSize(k == 0 ? 0.065 : 0.); fr->GetXaxis()->SetLabelSize(0.);
      fr->Draw("axis"); mPt[k]->Draw("same e1"); tPt[k]->Draw("same e1");
      label(k == 0 ? 0.42 : 0.06, 0.86, Form("%.0f < p_{T}^{jet} < %.0f GeV", win[k][0], win[k][1]), 0.068);
      label(k == 0 ? 0.42 : 0.06, 0.77, Form("#LTp_{T}#GT: %.1f / %.1f GeV", tMean[k], mMean[k]), 0.060);
      if(k == nWin - 1){ lg = makeLegend(0.06, 0.05, 0.80, 0.25, 0.058); lg->AddEntry(tPt[k], "data leading track", "pe"); lg->AddEntry(mPt[k], "PYTHIA light-jet muon", "pe"); lg->Draw(); }
      TPad *q = cell(c, k, false, split); q->SetLogx();
      styleH(wPt[k], hexCorrected, markFilledCircle, 0.8);
      wPt[k]->GetXaxis()->SetRangeUser(15, 150); wPt[k]->SetMinimum(0.); wPt[k]->SetMaximum(3.5);
      styleRatioAxes(wPt[k], "p_{T} [GeV]", k == 0 ? "w = M_{L} / T" : "");
      wPt[k]->GetXaxis()->SetTitleSize(0.10); wPt[k]->GetXaxis()->SetLabelSize(0.085); wPt[k]->GetXaxis()->SetMoreLogLabels(); wPt[k]->GetXaxis()->SetNoExponent();
      wPt[k]->GetYaxis()->SetTitleSize(0.10); wPt[k]->GetYaxis()->SetTitleOffset(0.75); wPt[k]->GetYaxis()->SetLabelSize(k == 0 ? 0.085 : 0.);
      wPt[k]->Draw("e1"); unityLine(15, 150)->Draw();
    }
    c->cd(); label(0.045, 0.995, hdr, 0.040, 13);
    savePdfTight(c, Form("%s/dataDrivenLight_trackPtWeight.pdf", outDir));
  }

  // ---------------- figure 2: light ptRel templates
  {
    TCanvas *c = new TCanvas("cTpl", "", 1700, 620);
    for(int k = 0; k < nWin; k++){
      TPad *p = cell(c, k, true, split);
      styleH(lUn[k], hexData, markOpenCircle, 0.8); styleH(lRw[k], hexCorrected, markFilledCircle, 0.8);
      styleLine(lMC[k], hexMC, 3); styleLine(lMC0[k], hexBkg2, 3); lMC0[k]->SetLineStyle(2);
      TH1D *fr = (TH1D*) lMC[k]->Clone(Form("fr2_%d", k)); fr->Reset(); fr->SetDirectory(nullptr);
      fr->SetMinimum(0.); fr->SetMaximum(1.75*TMath::Max(TMath::Max(lMC0[k]->GetMaximum(), lMC[k]->GetMaximum()), lRw[k]->GetMaximum()));
      fr->GetYaxis()->SetTitle(k == 0 ? "(1/N) dN/dp_{T}^{rel} [GeV^{-1}]" : ""); fr->GetYaxis()->SetTitleSize(0.075); fr->GetYaxis()->SetTitleOffset(1.25);
      fr->GetYaxis()->SetLabelSize(k == 0 ? 0.065 : 0.); fr->GetXaxis()->SetLabelSize(0.);
      fr->Draw("axis"); lMC0[k]->Draw("same hist"); lMC[k]->Draw("same hist"); lUn[k]->Draw("same e1"); lRw[k]->Draw("same e1");
      label(k == 0 ? 0.42 : 0.06, 0.88, Form("%.0f < p_{T}^{jet} < %.0f GeV", win[k][0], win[k][1]), 0.068);
      if(k == nWin - 1){ TLegend *lg = makeLegend(0.14, 0.52, 0.98, 0.80, 0.050);
        lg->AddEntry(lUn[k], "data, unweighted", "pe"); lg->AddEntry(lRw[k], "data, reweighted", "pe");
        lg->AddEntry(lMC[k], "PYTHIA light (fit)", "l"); lg->AddEntry(lMC0[k], "PYTHIA light, unsmeared", "l"); lg->Draw(); }
      TPad *q = cell(c, k, false, split);
      TH1D *rU = makeRatio(lUn[k], lMC[k], Form("rU%d", k)), *rR = makeRatio(lRw[k], lMC[k], Form("rR%d", k)), *r0 = makeRatio(lMC0[k], lMC[k], Form("r0%d", k));
      for(int i = 1; i <= r0->GetNbinsX(); i++) r0->SetBinError(i, 0.);
      styleH(rU, hexData, markOpenCircle, 0.8); styleH(rR, hexCorrected, markFilledCircle, 0.8); styleLine(r0, hexBkg2, 3); r0->SetLineStyle(2);
      rR->SetMinimum(0.); rR->SetMaximum(3.); styleRatioAxes(rR, "p_{T}^{rel} [GeV]", k == 0 ? "/ PYTHIA (fit)" : "");
      rR->GetXaxis()->SetTitleSize(0.10); rR->GetXaxis()->SetLabelSize(0.085);
      rR->GetYaxis()->SetTitleSize(0.10); rR->GetYaxis()->SetTitleOffset(0.75); rR->GetYaxis()->SetLabelSize(k == 0 ? 0.085 : 0.);
      rR->Draw("e1"); r0->Draw("same hist"); rU->Draw("same e1"); rR->Draw("same e1"); unityLine(0, fitHi)->Draw();
    }
    c->cd(); label(0.045, 0.995, hdr, 0.040, 13);
    savePdfTight(c, Form("%s/dataDrivenLight_ptRelTemplates.pdf", outDir));
  }

  // ---------------- figure 3: b purity vs jet pT
  {
    TCanvas *c = new TCanvas("cPur", "", 700, 700);
    TPad *top, *bot; splitPads(top, bot);
    double x[nWin], ex[nWin], z[nWin];
    for(int k = 0; k < nWin; k++){ x[k] = 0.5*(win[k][0] + win[k][1]); ex[k] = 0.5*(win[k][1] - win[k][0]); z[k] = 0.; }
    auto gr = [&](const double *y, const char *hex, int mark){ TGraphErrors *g = new TGraphErrors(nWin, x, y, ex, z); int col = TColor::GetColor(hex); g->SetLineColor(col); g->SetMarkerColor(col); g->SetMarkerStyle(mark); g->SetMarkerSize(1.3); g->SetLineWidth(2); return g; };
    double rNom[nWin], rUn[nWin], rRw[nWin], one[nWin];
    for(int k = 0; k < nWin; k++){ rNom[k] = pNom[k]/truth[k]; rUn[k] = pUn[k]/truth[k]; rRw[k] = pRw[k]/truth[k]; one[k] = 1.; }
    TGraphErrors *gT = gr(truth, hexData, markFilledSquare), *gN = gr(pNom, hexMC, markOpenSquare), *gU = gr(pUn, okabeHex[1], markOpenCircle), *gR = gr(pRw, hexCorrected, markFilledDiamond);
    gT->SetLineStyle(1); gT->SetMarkerSize(0); gT->SetLineWidth(3);
    top->cd();
    TH1D *fr = new TH1D("frP", "", 1, 80, 300); fr->SetStats(0); fr->SetMinimum(0.25); fr->SetMaximum(0.80);
    fr->GetYaxis()->SetTitle("b purity"); fr->GetYaxis()->SetTitleSize(0.065); fr->GetYaxis()->SetTitleOffset(1.15); fr->GetYaxis()->SetLabelSize(0.055); fr->GetXaxis()->SetLabelSize(0.);
    fr->Draw("axis"); gT->Draw("pz same"); gN->Draw("pz same"); gU->Draw("pz same"); gR->Draw("pz same");
    TLegend *lg = makeLegend(0.20, 0.62, 0.62, 0.88, 0.045);
    lg->AddEntry(gT, "PYTHIA truth", "l"); lg->AddEntry(gN, "fit, PYTHIA light (nominal)", "p");
    lg->AddEntry(gU, "fit, data tracks unweighted", "p"); lg->AddEntry(gR, "fit, data tracks reweighted", "p"); lg->Draw();
    label(0.20, 0.97, hdr, 0.042, 13); label(0.62, 0.86, Form("c, b: PYTHIA; fit 0-%g GeV", fitEdge), 0.040, 13);
    bot->cd();
    TH1D *fb = new TH1D("frR", "", 1, 80, 300); fb->SetStats(0); fb->SetMinimum(0.8); fb->SetMaximum(1.7);
    styleRatioAxes(fb, "p_{T}^{jet} [GeV]", "fit / truth"); fb->Draw("axis");
    unityLine(80, 300)->Draw();
    gr(rNom, hexMC, markOpenSquare)->Draw("pz same"); gr(rUn, okabeHex[1], markOpenCircle)->Draw("pz same"); gr(rRw, hexCorrected, markFilledDiamond)->Draw("pz same");
    savePdfTight(c, Form("%s/dataDrivenLight_purity_fit0to%g.pdf", outDir, fitEdge));
  }
}
