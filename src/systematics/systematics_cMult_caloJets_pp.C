// Systematic on the pp CALO-jet b purity from the c-jet fraction.
//
// In the 2-template fit (b vs light+c) the c template is folded into the light
// one at its MC-truth weight times a c-multiplier, as c_multiplier does in
// templateFitter(). Nominal is 1.0; it is scanned over 0.8-1.2 on the grid
// systematics_cMult.C uses, and as there the systematic is the largest
// |relative deviation| from nominal in each jet-pT window. The b template is
// held at the nominal gluon-splitting setting (bGS share + 0.175), so this
// varies the c fraction alone.
//
// Fit, templates and inputs: headers/functions/caloJetBPurityFit_pp.h, shared
// with src/plots/bPurity/plotBPurity_caloJets_pp.C and
// src/systematics/systematics_bGS_caloJets_pp.C.
//
// Output: rootFiles/systematics/bPuritySys_cMult_caloJets_pp.root
//   h_bPurity_nominal        purity at c-multiplier 1.0, fit (stat) error
//   h_bPurity_cMult1pX       purity at each scanned multiplier
//   h_sysRel_cMult           max |relative deviation|   <- the systematic
//   h_sysAbs_cMult           the same in absolute purity
//   h_sysRel_cMult_up/_down  signed relative deviation at 1.2 / 0.8
// and figures/systematics/cMult_caloJets_pp.pdf
//
// Usage: root -l -b -q systematics_cMult_caloJets_pp.C
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

const double cNominal = 1.0;
const int NC = 5;
const double cMultVals[NC] = {0.8, 0.9, 1.0, 1.1, 1.2};
const double gsShiftNominal = 0.175;

// Systematics are evaluated in every window, underflow (50-80 GeV) included,
// since the unfolding needs them there too; windows without enough data (see
// fittable() in the header) are skipped and left empty.
const int NPtSys = NPt;

const char *outRoot = "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/systematics/bPuritySys_cMult_caloJets_pp.root";
const char *outFig  = "/home/clayton/Analysis/code/bJetRaaAnalysis/figures/systematics/cMult_caloJets_pp.pdf";

void systematics_cMult_caloJets_pp()
{
  initPlotStyle();
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  gSystem->mkdir(gSystem->DirName(outRoot), kTRUE);
  gSystem->mkdir(gSystem->DirName(outFig), kTRUE);

  TFile *fD = TFile::Open(dataPath), *fM = TFile::Open(mcPath);
  if(!fD || fD->IsZombie() || !fM || fM->IsZombie()){ printf("ERROR: cannot open inputs\n"); return; }

  TH1D *hP[NC];
  for(int k = 0; k < NC; k++){
    hP[k] = new TH1D(Form("h_bPurity_cMult%.1f", cMultVals[k]),
                     Form("b purity, c-multiplier %.1f;#it{p}_{T}^{jet} [GeV];b purity", cMultVals[k]),
                     NPtSys, ptEdges);
    hP[k]->SetName(TString(hP[k]->GetName()).ReplaceAll(".", "p"));
    hP[k]->SetDirectory(nullptr);
  }
  int kNom = -1; for(int k = 0; k < NC; k++) if(fabs(cMultVals[k] - cNominal) < 1e-9) kNom = k;

  TH1D *hRel  = new TH1D("h_sysRel_cMult",      "c-multiplier systematic (max |relative deviation|);#it{p}_{T}^{jet} [GeV];relative uncertainty", NPtSys, ptEdges);
  TH1D *hAbs  = new TH1D("h_sysAbs_cMult",      "c-multiplier systematic (absolute);#it{p}_{T}^{jet} [GeV];absolute uncertainty", NPtSys, ptEdges);
  TH1D *hUp   = new TH1D("h_sysRel_cMult_up",   "relative deviation at c-multiplier 1.2;#it{p}_{T}^{jet} [GeV];(var - nom)/nom", NPtSys, ptEdges);
  TH1D *hDown = new TH1D("h_sysRel_cMult_down", "relative deviation at c-multiplier 0.8;#it{p}_{T}^{jet} [GeV];(var - nom)/nom", NPtSys, ptEdges);
  for(TH1D *h : {hRel, hAbs, hUp, hDown}) h->SetDirectory(nullptr);

  printf("\n%-9s", "jet pT"); for(int k = 0; k < NC; k++) printf("  c x%.1f ", cMultVals[k]); printf("   sysRel\n");
  for(int i = 0; i < NPtSys; i++){
    double lo = ptEdges[i], hi = ptEdges[i+1];
    TH1D *d05 = range05(projWin(fD, hName, lo, hi), Form("d05_%d", i));
    if(!fittable(d05)){ printf("%3.0f-%-5.0f  not fitted: %.0f data entries\n", lo, hi, d05->Integral()); continue; }
    for(int k = 0; k < NC; k++){
      Tpl t = buildTemplates(fM, lo, hi, cMultVals[k]);
      TH1D *b = bWithGS(t, TMath::Min(1., t.fGSnat + gsShiftNominal));
      FitRes r = doFit(d05, {b, t.lc});
      hP[k]->SetBinContent(i+1, r.fb); hP[k]->SetBinError(i+1, r.efb);
    }
    double nom = hP[kNom]->GetBinContent(i+1), maxDev = 0.;
    for(int k = 0; k < NC; k++) maxDev = TMath::Max(maxDev, fabs(hP[k]->GetBinContent(i+1) - nom)/nom);
    hRel->SetBinContent(i+1, maxDev);
    hAbs->SetBinContent(i+1, maxDev*nom);
    hUp  ->SetBinContent(i+1, (hP[NC-1]->GetBinContent(i+1) - nom)/nom);
    hDown->SetBinContent(i+1, (hP[0]->GetBinContent(i+1) - nom)/nom);

    printf("%3.0f-%-5.0f", lo, hi); for(int k = 0; k < NC; k++) printf("  %.4f ", hP[k]->GetBinContent(i+1));
    printf("   %.2f%%\n", 100*maxDev);
  }

  TFile *fo = TFile::Open(outRoot, "recreate");
  hP[kNom]->Write("h_bPurity_nominal");
  for(int k = 0; k < NC; k++) hP[k]->Write();
  hRel->Write(); hAbs->Write(); hUp->Write(); hDown->Write();
  TNamed info("info", Form("pp calo jets, 2-template ptRel fit %.0f-%.0f GeV. c weight in light+c = MC truth x c-multiplier; "
                            "nominal %.1f, scanned %.1f-%.1f; bGS share + %.3f throughout; systematic = max |relative deviation|. "
                            "data: %s  MC: %s", fitLo, fitHi, cNominal, cMultVals[0], cMultVals[NC-1], gsShiftNominal,
                            gSystem->BaseName(dataPath), gSystem->BaseName(mcPath)));
  info.Write();
  fo->Close();
  printf("\nwritten %s\n", outRoot);

  // ---- figure: purity variations (top), relative deviation (bottom) --------
  TCanvas *c = new TCanvas("c", "", 700, 800);
  TPad *pTop, *pBot; splitPads(pTop, pBot);
  pTop->cd();
  TH1D *fr = (TH1D*) hP[kNom]->Clone("frame"); fr->Reset(); fr->SetTitle("");
  fr->SetMinimum(0.2); fr->SetMaximum(0.75);
  fr->GetXaxis()->SetLabelSize(0);
  fr->GetYaxis()->SetTitle("b-jet purity"); fr->GetYaxis()->SetTitleSize(0.055);
  fr->GetYaxis()->SetTitleOffset(1.35); fr->GetYaxis()->SetLabelSize(0.045);
  fr->Draw("AXIS");
  drawUnderflowBand(fr->GetMinimum(), fr->GetMaximum(), 0.25);
  TH1D *hLo = (TH1D*) hP[0]->Clone("hLo"), *hHi = (TH1D*) hP[NC-1]->Clone("hHi"), *hN = (TH1D*) hP[kNom]->Clone("hN");
  styleH(hN,  "#000000", markFilledCircle, 1.2);
  styleH(hLo, "#0072B2", markOpenSquare,   1.2);
  styleH(hHi, "#D55E00", markOpenDiamond,  1.5);
  for(int i = 1; i <= NPtSys; i++){ hLo->SetBinError(i, 0); hHi->SetBinError(i, 0); }
  hLo->Draw("P same"); hHi->Draw("P same"); hN->Draw("E1 same");
  TLegend *leg = makeLegend(0.53, 0.70, 0.95, 0.88, 0.036);   // purities stay below ~0.52
  leg->AddEntry(hN,  "nominal: c-multiplier 1.0", "lp");
  leg->AddEntry(hLo, "c-multiplier 0.8", "p");
  leg->AddEntry(hHi, "c-multiplier 1.2", "p");
  leg->Draw();
  TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.046);
  la.DrawLatex(0.21, 0.84, "pp 5.02 TeV, calo jets");
  la.SetTextSize(0.038);
  la.DrawLatex(0.21, 0.78, "2-template #it{p}_{T}^{rel} fit");

  pBot->cd();
  TH1D *rU = (TH1D*) hUp->Clone("rU"), *rD = (TH1D*) hDown->Clone("rD"), *rS = (TH1D*) hRel->Clone("rS");
  // deviations carry no error of their own; ROOT would otherwise draw sqrt(N)
  for(TH1D *h : {rU, rD, rS}){ h->Scale(100.); for(int i = 0; i <= h->GetNbinsX()+1; i++) h->SetBinError(i, 0.); }
  for(TH1D *h : {rU, rD}) blankUnfitted(h, hP[kNom]);
  styleH(rD, "#0072B2", markOpenSquare, 1.2); styleH(rU, "#D55E00", markOpenDiamond, 1.5);
  styleLine(rS, "#000000", 2);
  TH1D *rF = (TH1D*) rS->Clone("rF"); rF->Reset();
  double m = 1.3*TMath::Max(rS->GetMaximum(), 1.0);
  rF->SetMinimum(-m); rF->SetMaximum(m);
  styleRatioAxes(rF, "#it{p}_{T}^{jet} [GeV]", "#Delta / nom. [%]");
  rF->Draw("AXIS");
  drawUnderflowBand(rF->GetMinimum(), rF->GetMaximum(), -1e9);   // no label
  TLine *z = new TLine(ptEdges[0], 0, ptEdges[NPtSys], 0); z->SetLineStyle(2); z->SetLineColor(kGray+2); z->Draw();
  TH1D *rSm = (TH1D*) rS->Clone("rSm"); rSm->Scale(-1.);
  rS->Draw("HIST same"); rSm->Draw("HIST same");
  rD->Draw("P same"); rU->Draw("P same");

  c->SaveAs(outFig);
  printf("written %s\n", outFig);
}
