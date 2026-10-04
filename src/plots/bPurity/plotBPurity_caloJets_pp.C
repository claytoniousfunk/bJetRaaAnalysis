// pp b-jet purity from muon ptRel template fits, CALO jets.
//
// Data and templates are both calo jets (never mix with PF -- the two
// collections have different responses):
//   data       pp SingleMuon, ak4Calo, manual JEC                 (2026-9-29)
//   templates  PYTHIA DiJet pThat>15, ak4Calo, manual JEC, flavour and
//              bHadronNumber taken from the matched ak4PF jet      (2026-10-4)
// The PF-matched labels matter: the earlier calo scans used refparton_flavor,
// which keeps only ~35-45% of true b jets as b and has no gluon-splitting
// (bGS) label, so their b template was narrower than the light one at high pT.
//
// Three fit configurations, all RooFit binned likelihood over 0-5 GeV on the
// 0.1 GeV histogram binning, c folded into light at its MC-truth fraction as in
// templateFitter():
//   2tpl_naturalGS   b (+bGS at PYTHIA's own share) vs light+c
//   2tpl_GSplus0p175 as above with the bGS share raised by 0.175 ABSOLUTE,
//                    which is what bGS_multiplier = 1.175 does in
//                    calculateBJetsPerZ.cc
//   3tpl_freeGS      b (flavour creation), bGS and light+c, both b fractions
//                    free; purity = their sum
//
// Outputs (PDF) in figures/bPurity/caloJets_pp/:
//   ptRelFit_<config>_jetPt<lo>-<hi>.pdf   one fit per jet-pT window
//   bPurity_vsJetPt.pdf                    purity vs jet pT, all configs + MC truth
//   bGSShare_vsJetPt.pdf                   fitted vs natural bGS share (3tpl)
//
// Inputs, templates and the fit itself live in
// headers/functions/caloJetBPurityFit_pp.h, shared with the bGS systematic.
//
// The calo data histogram is empty below 70 GeV and 70-80 looks like a turn-on,
// so the windows start at 70 and the first point should be read with care.
//
// Usage: root -l -b -q plotBPurity_caloJets_pp.C
// Run from: src/plots/bPurity/

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TLatex.h"
#include "TSystem.h"
#include "TColor.h"
#include "TGaxis.h"
#include "TStyle.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooAddPdf.h"
#include "RooMsgService.h"
#include "../../../headers/plotting/plotStyle.h"
#include "../../../headers/plotting/ratioPanel.h"
#include "../../../headers/functions/divideByBinwidth.h"
#include "../../../headers/functions/caloJetBPurityFit_pp.h"

const char *outDir = "/home/clayton/Analysis/code/bJetRaaAnalysis/figures/bPurity/caloJets_pp";


// fixed colours per component, so every fit figure reads the same way
const char *hexB    = "#D55E00";   // b (total, or flavour creation in 3tpl)
const char *hexBGS  = "#E69F00";   // gluon-splitting b (3tpl only)
const char *hexL    = "#0072B2";   // light + c
const char *hexFit  = "#009E73";   // total fit

// display copy: shared ptRel binning, converted to a density
TH1D* forDisplay(TH1D *h05, const char *name)
{
  TH1D *r = rebinTo(h05, nEdge_ptRel_0to5, edge_ptRel_0to5, name);
  divideByBinwidth(r);
  return r;
}

void drawFit(TH1D *data05, const FitRes &fr, bool threeTpl, double lo, double hi,
             const char *configLabel, const TString &pdf)
{
  TH1D *hD = forDisplay(data05, Form("dD%d", uidP++));
  TH1D *hB = forDisplay(fr.compB, Form("dB%d", uidP++));
  TH1D *hG = threeTpl ? forDisplay(fr.compBGS, Form("dG%d", uidP++)) : nullptr;
  TH1D *hL = forDisplay(fr.compL, Form("dL%d", uidP++));
  TH1D *hT = forDisplay(fr.total, Form("dT%d", uidP++));
  // the fit total is a fixed reference here: only the data error goes on the ratio
  TH1D *hR = makeRatio(hD, hT, Form("dR%d", uidP++), RatioErr::kNumerator);

  styleH(hD, hexData, markFilledCircle, 1.0);
  styleLine(hB, hexB, 3); styleLine(hL, hexL, 3); styleLine(hT, hexFit, 3);
  hT->SetLineStyle(1); hB->SetLineStyle(2); hL->SetLineStyle(2);
  if(hG){ styleLine(hG, hexBGS, 3); hG->SetLineStyle(3); }
  styleH(hR, hexData, markFilledCircle, 1.0);

  TCanvas *c = new TCanvas(Form("c%d", uidP++), "", 700, 800);
  TPad *pTop, *pBot; splitPads(pTop, pBot);
  pTop->cd();
  double ymax = TMath::Max(hD->GetMaximum(), hT->GetMaximum());
  hD->SetMinimum(0.); hD->SetMaximum(1.45*ymax);
  hD->SetTitle("");
  hD->GetXaxis()->SetLabelSize(0);
  hD->GetYaxis()->SetTitle("d#it{N}/d#it{p}_{T}^{rel} [GeV^{-1}]");
  hD->GetYaxis()->SetTitleSize(0.055); hD->GetYaxis()->SetTitleOffset(1.45);
  hD->GetYaxis()->SetLabelSize(0.045);
  hD->Draw("E");
  hL->Draw("HIST same"); hB->Draw("HIST same"); if(hG) hG->Draw("HIST same");
  hT->Draw("HIST same"); hD->Draw("E same");

  TLegend *leg = makeLegend(0.55, threeTpl ? 0.56 : 0.62, 0.93, 0.89, 0.042);
  leg->AddEntry(hD, "data", "lp");
  leg->AddEntry(hT, "fit", "l");
  if(threeTpl){
    leg->AddEntry(hB,  Form("b, flavor creation (%.2f)", fr.fb), "l");
    leg->AddEntry(hG,  Form("b, gluon splitting (%.2f)", fr.fGS), "l");
    leg->AddEntry(hL,  Form("light + c (%.2f)", 1. - fr.fb - fr.fGS), "l");
  } else {
    leg->AddEntry(hB,  Form("b (%.3f #pm %.3f)", fr.fb, fr.efb), "l");
    leg->AddEntry(hL,  Form("light + c (%.3f)", 1. - fr.fb), "l");
  }
  leg->Draw();

  TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.044);
  la.DrawLatex(0.21, 0.85, "pp 5.02 TeV, calo jets");
  la.DrawLatex(0.21, 0.79, Form("%.0f < #it{p}_{T}^{jet} < %.0f GeV", lo, hi));
  la.SetTextSize(0.038);
  la.DrawLatex(0.21, 0.73, configLabel);
  if(threeTpl) la.DrawLatex(0.21, 0.68, Form("b purity = %.3f", fr.fb + fr.fGS));
  la.DrawLatex(0.21, threeTpl ? 0.63 : 0.68, Form("#chi^{2}/ndf = %.1f/%d", fr.chi2, fr.ndf));

  pBot->cd();
  styleRatioAxes(hR, "#it{p}_{T}^{rel} [GeV]", "data / fit");
  hR->SetMinimum(0.5); hR->SetMaximum(1.5);
  hR->Draw("E");
  unityLine(fitLo, fitHi)->Draw();
  hR->Draw("E same");

  c->SaveAs(pdf);
  delete c;
}

void plotBPurity_caloJets_pp()
{
  initPlotStyle();
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  gSystem->mkdir(outDir, kTRUE);

  TFile *fD = TFile::Open(dataPath), *fM = TFile::Open(mcPath);
  if(!fD || fD->IsZombie() || !fM || fM->IsZombie()){ printf("ERROR: cannot open inputs\n"); return; }

  double x[NPt], ex[NPt];
  double p2n[NPt], e2n[NPt], p2g[NPt], e2g[NPt], p3[NPt], e3[NPt], truth[NPt];
  double gsFit[NPt], egsFit[NPt], gsNat[NPt];
  bool fitted[NPt];

  printf("\n%-9s %-15s %-15s %-15s %-8s | %-12s %-7s\n", "jet pT", "2tpl natGS", "2tpl GS+0.175",
         "3tpl freeGS", "truth", "GS fit", "GS nat");
  for(int i = 0; i < NPt; i++){
    double lo = ptEdges[i], hi = ptEdges[i+1];
    x[i] = 0.5*(lo + hi); ex[i] = 0.5*(hi - lo);

    TH1D *d05 = range05(projWin(fD, hName, lo, hi), Form("d05_%d", i));
    Tpl t = buildTemplates(fM, lo, hi);
    truth[i] = t.truthB; gsNat[i] = t.fGSnat;
    fitted[i] = fittable(d05);
    if(!fitted[i]){
      p2n[i] = e2n[i] = p2g[i] = e2g[i] = p3[i] = e3[i] = gsFit[i] = egsFit[i] = 0.;
      printf("%3.0f-%-5.0f not fitted: %.0f data entries (MC truth %.3f)\n", lo, hi, d05->Integral(), truth[i]);
      continue;
    }

    TH1D *bNat = bWithGS(t, t.fGSnat);
    TH1D *bUp  = bWithGS(t, TMath::Min(1., t.fGSnat + 0.175));

    FitRes rN = doFit(d05, {bNat, t.lc});
    FitRes rU = doFit(d05, {bUp,  t.lc});
    FitRes r3 = doFit(d05, {t.bFC, t.bGS, t.lc});

    p2n[i] = rN.fb; e2n[i] = rN.efb;
    p2g[i] = rU.fb; e2g[i] = rU.efb;
    p3[i]  = r3.fb + r3.fGS;
    // purity error for the 3tpl fit: the two fractions are anticorrelated, so
    // adding in quadrature overstates it; quoted as an upper bound
    e3[i]  = sqrt(r3.efb*r3.efb + r3.efGS*r3.efGS);
    gsFit[i] = p3[i] > 0. ? r3.fGS/p3[i] : 0.;
    egsFit[i] = p3[i] > 0. ? r3.efGS/p3[i] : 0.;

    printf("%3.0f-%-5.0f %.3f+-%.3f(%4.1f) %.3f+-%.3f(%4.1f) %.3f+-%.3f(%4.1f) %.3f   | %.2f+-%.2f  %.2f\n",
           lo, hi, rN.fb, rN.efb, rN.chi2, rU.fb, rU.efb, rU.chi2, p3[i], e3[i], r3.chi2, truth[i],
           gsFit[i], egsFit[i], gsNat[i]);

    drawFit(d05, rN, false, lo, hi, "2 templates, natural g#rightarrowb#bar{b} share",
            Form("%s/ptRelFit_2tpl_naturalGS_jetPt%.0f-%.0f.pdf", outDir, lo, hi));
    drawFit(d05, rU, false, lo, hi, "2 templates, g#rightarrowb#bar{b} share + 0.175",
            Form("%s/ptRelFit_2tpl_GSplus0p175_jetPt%.0f-%.0f.pdf", outDir, lo, hi));
    drawFit(d05, r3, true, lo, hi, "3 templates, g#rightarrowb#bar{b} share free",
            Form("%s/ptRelFit_3tpl_freeGS_jetPt%.0f-%.0f.pdf", outDir, lo, hi));
  }

  // ---- summary: purity vs jet pT ------------------------------------------
  {
    TCanvas *c = new TCanvas("cSum", "", 800, 700);
    c->SetLeftMargin(0.14); c->SetBottomMargin(0.13); c->SetRightMargin(0.04); c->SetTopMargin(0.06);
    TH1F *fr = c->DrawFrame(ptEdges[0], 0., ptEdges[NPt], 1.0);
    fr->GetXaxis()->SetTitle("#it{p}_{T}^{jet} [GeV]");
    fr->GetYaxis()->SetTitle("b-jet purity");
    fr->GetXaxis()->SetTitleSize(0.048); fr->GetYaxis()->SetTitleSize(0.048);
    fr->GetXaxis()->SetLabelSize(0.040); fr->GetYaxis()->SetLabelSize(0.040);

    TGraphErrors *gN = new TGraphErrors(NPt, x, p2n, ex, e2n);
    TGraphErrors *gU = new TGraphErrors(NPt, x, p2g, ex, e2g);
    TGraphErrors *g3 = new TGraphErrors(NPt, x, p3,  ex, e3);
    double zero[NPt] = {0};
    TGraphErrors *gT = new TGraphErrors(NPt, x, truth, ex, zero);
    // unfitted windows: no data point (the MC truth stays, it exists there)
    for(int i = NPt-1; i >= 0; i--) if(!fitted[i]){ gN->RemovePoint(i); gU->RemovePoint(i); g3->RemovePoint(i); }
    drawUnderflowBand(0., 1.0, 0.05);
    auto sty = [](TGraphErrors *g, const char *hex, int m){ int col = TColor::GetColor(hex);
      g->SetLineColor(col); g->SetMarkerColor(col); g->SetMarkerStyle(m); g->SetMarkerSize(1.3); g->SetLineWidth(2); };
    sty(gN, "#D55E00", markFilledCircle);
    sty(gU, "#0072B2", markFilledSquare);
    sty(g3, "#009E73", markFilledDiamond);
    sty(gT, "#000000", markOpenCircle);
    gT->Draw("P same"); gU->Draw("P same"); g3->Draw("P same"); gN->Draw("P same");

    // upper right is empty (purities stay below ~0.7, errors included)
    TLegend *leg = makeLegend(0.45, 0.74, 0.95, 0.93, 0.032);
    leg->AddEntry(gN, "2 templates, natural g#rightarrowb#bar{b} share", "lp");
    leg->AddEntry(gU, "2 templates, g#rightarrowb#bar{b} share + 0.175", "lp");
    leg->AddEntry(g3, "3 templates, g#rightarrowb#bar{b} share free", "lp");
    leg->AddEntry(gT, "PYTHIA truth", "p");
    leg->Draw();
    TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.040);
    la.DrawLatex(0.18, 0.87, "pp 5.02 TeV, calo jets");
    la.SetTextSize(0.032);
    la.DrawLatex(0.18, 0.82, Form("fit range %.0f < #it{p}_{T}^{rel} < %.0f GeV", fitLo, fitHi));
    savePdfTight(c, Form("%s/bPurity_vsJetPt.pdf", outDir));
    delete c;
  }

  // ---- summary: fitted vs natural bGS share --------------------------------
  {
    TCanvas *c = new TCanvas("cGS", "", 800, 700);
    c->SetLeftMargin(0.14); c->SetBottomMargin(0.13); c->SetRightMargin(0.04); c->SetTopMargin(0.06);
    TH1F *fr = c->DrawFrame(ptEdges[0], 0., ptEdges[NPt], 0.8);
    fr->GetXaxis()->SetTitle("#it{p}_{T}^{jet} [GeV]");
    fr->GetYaxis()->SetTitle("g#rightarrowb#bar{b} share of b jets");
    fr->GetXaxis()->SetTitleSize(0.048); fr->GetYaxis()->SetTitleSize(0.048);
    fr->GetXaxis()->SetLabelSize(0.040); fr->GetYaxis()->SetLabelSize(0.040);
    double zero[NPt] = {0}, up[NPt];
    for(int i = 0; i < NPt; i++) up[i] = TMath::Min(1., gsNat[i] + 0.175);
    TGraphErrors *gF = new TGraphErrors(NPt, x, gsFit, ex, egsFit);
    TGraphErrors *gNt = new TGraphErrors(NPt, x, gsNat, ex, zero);
    TGraphErrors *gUp = new TGraphErrors(NPt, x, up, ex, zero);
    for(int i = NPt-1; i >= 0; i--) if(!fitted[i]) gF->RemovePoint(i);
    drawUnderflowBand(0., 0.8, 0.04);
    auto sty = [](TGraphErrors *g, const char *hex, int m){ int col = TColor::GetColor(hex);
      g->SetLineColor(col); g->SetMarkerColor(col); g->SetMarkerStyle(m); g->SetMarkerSize(1.3); g->SetLineWidth(2); };
    sty(gF, "#009E73", markFilledDiamond);
    sty(gNt, "#000000", markOpenCircle);
    sty(gUp, "#0072B2", markOpenSquare);
    gNt->Draw("P same"); gUp->Draw("P same"); gF->Draw("P same");
    TLegend *leg = makeLegend(0.17, 0.70, 0.94, 0.86, 0.034);
    leg->AddEntry(gF, "fitted (3 templates)", "lp");
    leg->AddEntry(gNt, "PYTHIA natural", "p");
    leg->AddEntry(gUp, "natural + 0.175 (bGS_multiplier = 1.175)", "p");
    leg->Draw();
    TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.040);
    la.DrawLatex(0.18, 0.885, "pp 5.02 TeV, calo jets");
    savePdfTight(c, Form("%s/bGSShare_vsJetPt.pdf", outDir));
    delete c;
  }

  printf("\nfigures written to %s\n", outDir);
}
