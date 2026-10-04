// Systematic on the pp CALO-jet b purity from the ptRel fit range.
//
// Nominal fit range 0-5 GeV. As in systematics_lowerBound.C, the LOWER edge is
// varied over 0, 0.1, 0.2, 0.3 GeV with the upper edge fixed, and the
// systematic is the largest |relative deviation| from nominal in each jet-pT
// window.
//
// The UPPER edge is also scanned (4.0, 4.5 GeV) but only stored as a
// sensitivity, not folded into the systematic: the fit range is documented to
// move the baseline purity strongly (CLAUDE.md, open problem 5), so whether and
// how to quote it is a separate decision.
//
// Like templateFitter(), the fitted fraction is the b fraction of the jets
// INSIDE the fit range, so moving an edge also changes which jets the purity
// refers to.
//
// Other settings nominal: bGS share + 0.175, c-multiplier 1.0. Fit, templates
// and inputs: headers/functions/caloJetBPurityFit_pp.h.
//
// Output: rootFiles/systematics/bPuritySys_fitRange_caloJets_pp.root
//   h_bPurity_nominal          purity at 0-5 GeV, fit (stat) error
//   h_bPurity_fitLow0pX        purity at each lower edge (upper 5)
//   h_bPurity_fitHigh4pX       purity at upper edge 4.0 / 4.5 (lower 0)
//   h_sysRel_fitRange          max |relative deviation| over the lower-edge scan  <- the systematic
//   h_sysAbs_fitRange          the same in absolute purity
//   h_devRel_fitLow0p3         signed relative deviation at lower edge 0.3
//   h_sensRel_fitHigh4p0/4p5   signed relative deviation at upper edge 4.0 / 4.5 (sensitivity only)
// and figures/systematics/fitRange_caloJets_pp.pdf
//
// Usage: root -l -b -q systematics_fitRange_caloJets_pp.C
// Run from: src/systematics/

#include "TCanvas.h"
#include "TLatex.h"
#include "TSystem.h"
#include "TColor.h"
#include "TGaxis.h"
#include "TStyle.h"
#include "RooMsgService.h"
#include "../../headers/plotting/plotStyle.h"
#include "../../headers/plotting/ratioPanel.h"
#include "../../headers/functions/caloJetBPurityFit_pp.h"

const int NLow = 4;
const double fitLowVals[NLow] = {0.0, 0.1, 0.2, 0.3};
const int NHigh = 2;
const double fitHighVals[NHigh] = {4.0, 4.5};
const double gsShiftNominal = 0.175;

// Systematics are evaluated in every window, underflow (50-80 GeV) included,
// since the unfolding needs them there too; windows without enough data (see
// fittable() in the header) are skipped and left empty.
const int NPtSys = NPt;

const char *outRoot = "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/systematics/bPuritySys_fitRange_caloJets_pp.root";
const char *outFig  = "/home/clayton/Analysis/code/bJetRaaAnalysis/figures/systematics/fitRange_caloJets_pp.pdf";

TH1D* bookPt(const char *name, const char *title)
{
  TH1D *h = new TH1D(name, title, NPtSys, ptEdges); h->SetDirectory(nullptr); return h;
}
TString tag(double v){ TString s = Form("%.1f", v); s.ReplaceAll(".", "p"); return s; }

void systematics_fitRange_caloJets_pp()
{
  initPlotStyle();
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  gSystem->mkdir(gSystem->DirName(outRoot), kTRUE);
  gSystem->mkdir(gSystem->DirName(outFig), kTRUE);

  TFile *fD = TFile::Open(dataPath), *fM = TFile::Open(mcPath);
  if(!fD || fD->IsZombie() || !fM || fM->IsZombie()){ printf("ERROR: cannot open inputs\n"); return; }

  TH1D *hLow[NLow], *hHigh[NHigh];
  for(int k = 0; k < NLow; k++)
    hLow[k] = bookPt("h_bPurity_fitLow" + tag(fitLowVals[k]), Form("b purity, fit %.1f-%.1f GeV;#it{p}_{T}^{jet} [GeV];b purity", fitLowVals[k], fitHi));
  for(int k = 0; k < NHigh; k++)
    hHigh[k] = bookPt("h_bPurity_fitHigh" + tag(fitHighVals[k]), Form("b purity, fit %.1f-%.1f GeV;#it{p}_{T}^{jet} [GeV];b purity", fitLo, fitHighVals[k]));
  TH1D *hNom = hLow[0];   // lower edge 0 = nominal 0-5 GeV

  TH1D *hRel  = bookPt("h_sysRel_fitRange",  "fit-range systematic (max |relative deviation|, lower edge 0-0.3 GeV);#it{p}_{T}^{jet} [GeV];relative uncertainty");
  TH1D *hAbs  = bookPt("h_sysAbs_fitRange",  "fit-range systematic (absolute);#it{p}_{T}^{jet} [GeV];absolute uncertainty");
  TH1D *hDevL = bookPt("h_devRel_fitLow0p3", "relative deviation at fit 0.3-5 GeV;#it{p}_{T}^{jet} [GeV];(var - nom)/nom");
  TH1D *hSens[NHigh];
  for(int k = 0; k < NHigh; k++)
    hSens[k] = bookPt("h_sensRel_fitHigh" + tag(fitHighVals[k]), Form("relative deviation at fit 0-%.1f GeV (sensitivity, not in systematic);#it{p}_{T}^{jet} [GeV];(var - nom)/nom", fitHighVals[k]));

  printf("\n%-9s", "jet pT");
  for(int k = 0; k < NLow; k++)  printf("  %.1f-5 ", fitLowVals[k]);
  for(int k = 0; k < NHigh; k++) printf("  0-%.1f ", fitHighVals[k]);
  printf("   sysRel   | 0-4.0 / 0-4.5 rel. dev.\n");

  for(int i = 0; i < NPtSys; i++){
    double lo = ptEdges[i], hi = ptEdges[i+1];
    TH1D *d05 = range05(projWin(fD, hName, lo, hi), Form("d05_%d", i));
    if(!fittable(d05)){ printf("%3.0f-%-5.0f  not fitted: %.0f data entries\n", lo, hi, d05->Integral()); continue; }
    Tpl t = buildTemplates(fM, lo, hi);
    TH1D *b = bWithGS(t, TMath::Min(1., t.fGSnat + gsShiftNominal));
    for(int k = 0; k < NLow; k++){
      FitRes r = doFit(d05, {b, t.lc}, fitLowVals[k], fitHi);
      hLow[k]->SetBinContent(i+1, r.fb); hLow[k]->SetBinError(i+1, r.efb);
    }
    for(int k = 0; k < NHigh; k++){
      FitRes r = doFit(d05, {b, t.lc}, fitLo, fitHighVals[k]);
      hHigh[k]->SetBinContent(i+1, r.fb); hHigh[k]->SetBinError(i+1, r.efb);
    }
    double nom = hNom->GetBinContent(i+1), maxDev = 0.;
    for(int k = 0; k < NLow; k++) maxDev = TMath::Max(maxDev, fabs(hLow[k]->GetBinContent(i+1) - nom)/nom);
    hRel->SetBinContent(i+1, maxDev);
    hAbs->SetBinContent(i+1, maxDev*nom);
    hDevL->SetBinContent(i+1, (hLow[NLow-1]->GetBinContent(i+1) - nom)/nom);
    for(int k = 0; k < NHigh; k++) hSens[k]->SetBinContent(i+1, (hHigh[k]->GetBinContent(i+1) - nom)/nom);

    printf("%3.0f-%-5.0f", lo, hi);
    for(int k = 0; k < NLow; k++)  printf("  %.4f", hLow[k]->GetBinContent(i+1));
    for(int k = 0; k < NHigh; k++) printf("  %.4f", hHigh[k]->GetBinContent(i+1));
    printf("   %5.2f%%   | %+5.1f%% / %+5.1f%%\n", 100*maxDev, 100*hSens[0]->GetBinContent(i+1), 100*hSens[1]->GetBinContent(i+1));
  }

  TFile *fo = TFile::Open(outRoot, "recreate");
  hNom->Write("h_bPurity_nominal");
  for(int k = 0; k < NLow; k++)  hLow[k]->Write();
  for(int k = 0; k < NHigh; k++) hHigh[k]->Write();
  hRel->Write(); hAbs->Write(); hDevL->Write();
  for(int k = 0; k < NHigh; k++) hSens[k]->Write();
  TNamed info("info", Form("pp calo jets, 2-template ptRel fit, nominal %.0f-%.0f GeV, bGS share + %.3f, c-multiplier 1.0. "
                            "Systematic: lower edge %.1f-%.1f GeV (upper %.0f), max |relative deviation|. Upper edge 4.0/4.5 stored as "
                            "sensitivity only. data: %s  MC: %s", fitLo, fitHi, gsShiftNominal, fitLowVals[0], fitLowVals[NLow-1], fitHi,
                            gSystem->BaseName(dataPath), gSystem->BaseName(mcPath)));
  info.Write();
  fo->Close();
  printf("\nwritten %s\n", outRoot);

  // ---- figure ---------------------------------------------------------------
  TCanvas *c = new TCanvas("c", "", 700, 800);
  TPad *pTop, *pBot; splitPads(pTop, pBot);
  pTop->cd();
  TH1D *fr = (TH1D*) hNom->Clone("frame"); fr->Reset(); fr->SetTitle("");
  fr->SetMinimum(0.2); fr->SetMaximum(0.75);
  fr->GetXaxis()->SetLabelSize(0);
  fr->GetYaxis()->SetTitle("b-jet purity"); fr->GetYaxis()->SetTitleSize(0.055);
  fr->GetYaxis()->SetTitleOffset(1.35); fr->GetYaxis()->SetLabelSize(0.045);
  fr->Draw("AXIS");
  drawUnderflowBand(fr->GetMinimum(), fr->GetMaximum(), 0.25);
  TH1D *hN = (TH1D*) hNom->Clone("hN"), *hL3 = (TH1D*) hLow[NLow-1]->Clone("hL3"), *hH4 = (TH1D*) hHigh[0]->Clone("hH4");
  styleH(hN,  "#000000", markFilledCircle, 1.2);
  styleH(hL3, "#0072B2", markOpenSquare,   1.2);
  styleH(hH4, "#CC79A7", markCross,        1.4);
  for(int i = 1; i <= NPtSys; i++){ hL3->SetBinError(i, 0); hH4->SetBinError(i, 0); }
  hL3->Draw("P same"); hH4->Draw("P same"); hN->Draw("E1 same");
  TLegend *leg = makeLegend(0.53, 0.68, 0.95, 0.88, 0.036);
  leg->AddEntry(hN,  "nominal: fit 0-5 GeV", "lp");
  leg->AddEntry(hL3, "fit 0.3-5 GeV", "p");
  leg->AddEntry(hH4, "fit 0-4 GeV (not in syst.)", "p");
  leg->Draw();
  TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.046);
  la.DrawLatex(0.21, 0.84, "pp 5.02 TeV, calo jets");
  la.SetTextSize(0.038);
  la.DrawLatex(0.21, 0.78, "2-template #it{p}_{T}^{rel} fit");

  pBot->cd();
  TH1D *rL = (TH1D*) hDevL->Clone("rL"), *rH = (TH1D*) hSens[0]->Clone("rH"), *rS = (TH1D*) hRel->Clone("rS");
  // deviations carry no error of their own; ROOT would otherwise draw sqrt(N)
  for(TH1D *h : {rL, rH, rS}){ h->Scale(100.); for(int i = 0; i <= h->GetNbinsX()+1; i++) h->SetBinError(i, 0.); }
  for(TH1D *h : {rL, rH}) blankUnfitted(h, hNom);
  styleH(rL, "#0072B2", markOpenSquare, 1.2); styleH(rH, "#CC79A7", markCross, 1.4);
  styleLine(rS, "#000000", 2);
  TH1D *rF = (TH1D*) rS->Clone("rF"); rF->Reset();
  double m = TMath::Max(rS->GetMaximum(), 1.0);
  for(int i = 1; i <= NPtSys; i++) if(rH->GetBinContent(i) > -900.) m = TMath::Max(m, fabs(rH->GetBinContent(i)));   // skip blanked bins
  m *= 1.3;
  rF->SetMinimum(-m); rF->SetMaximum(1.8*m);   // headroom for the legend row
  styleRatioAxes(rF, "#it{p}_{T}^{jet} [GeV]", "#Delta / nom. [%]");
  rF->Draw("AXIS");
  drawUnderflowBand(rF->GetMinimum(), rF->GetMaximum(), -1e9);   // no label
  TLine *z = new TLine(ptEdges[0], 0, ptEdges[NPtSys], 0); z->SetLineStyle(2); z->SetLineColor(kGray+2); z->Draw();
  TH1D *rSm = (TH1D*) rS->Clone("rSm"); rSm->Scale(-1.);
  rS->Draw("HIST same"); rSm->Draw("HIST same");
  rL->Draw("P same"); rH->Draw("P same");
  TLegend *legB = makeLegend(0.19, 0.78, 0.95, 0.97, 0.075);
  legB->SetNColumns(3);
  legB->AddEntry(rL, "0.3-5 GeV", "p");
  legB->AddEntry(rH, "0-4 GeV", "p");
  legB->AddEntry(rS, "syst. (max |#Delta|)", "l");
  legB->Draw();

  c->SaveAs(outFig);
  printf("written %s\n", outFig);
}
