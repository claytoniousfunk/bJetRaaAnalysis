// Systematic on the pp CALO-jet b purity from the UPPER edge of the ptRel fit
// range. Companion to systematics_fitRange_caloJets_pp.C, which varies the
// lower edge (0-0.3 GeV, as systematics_lowerBound.C did).
//
// Nominal fit range 0-5 GeV. The upper edge is varied over 4.0-5.0 GeV in
// 0.2 GeV steps (on the 0.1 GeV histogram binning) with the lower edge fixed
// at 0, and the systematic is the largest |relative deviation| from nominal in
// each jet-pT window. Only downward variations: the fit histograms stop at
// 5 GeV and the data above it are sparse.
//
// Like templateFitter(), the fitted fraction is the b fraction of the jets
// INSIDE the fit range, so moving the edge also changes which jets the purity
// refers to.
//
// Other settings nominal: bGS share + 0.175, c-multiplier 1.0. Fit, templates
// and inputs: headers/functions/caloJetBPurityFit_pp.h.
//
// Output: rootFiles/systematics/bPuritySys_fitRangeHigh_caloJets_pp.root
//   h_bPurity_nominal          purity at 0-5 GeV, fit (stat) error
//   h_bPurity_fitHigh4pX       purity at each upper edge
//   h_sysRel_fitRangeHigh      max |relative deviation|   <- the systematic
//   h_sysAbs_fitRangeHigh      the same in absolute purity
//   h_devRel_fitHigh4p0        signed relative deviation at upper edge 4.0
// and figures/systematics/fitRangeHigh_caloJets_pp.pdf
//
// Usage: root -l -b -q systematics_fitRangeHigh_caloJets_pp.C
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

const int NHigh = 6;
const double fitHighVals[NHigh] = {5.0, 4.8, 4.6, 4.4, 4.2, 4.0};   // [0] = nominal
const double gsShiftNominal = 0.175;

// Systematics are evaluated in every window, underflow (50-80 GeV) included,
// since the unfolding needs them there too; windows without enough data (see
// fittable() in the header) are skipped and left empty.
const int NPtSys = NPt;

const char *outRoot = "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/systematics/bPuritySys_fitRangeHigh_caloJets_pp.root";
const char *outFig  = "/home/clayton/Analysis/code/bJetRaaAnalysis/figures/systematics/fitRangeHigh_caloJets_pp.pdf";

TH1D* bookPt(const char *name, const char *title)
{
  TH1D *h = new TH1D(name, title, NPtSys, ptEdges); h->SetDirectory(nullptr); return h;
}
TString tag(double v){ TString s = Form("%.1f", v); s.ReplaceAll(".", "p"); return s; }

void systematics_fitRangeHigh_caloJets_pp()
{
  initPlotStyle();
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  gSystem->mkdir(gSystem->DirName(outRoot), kTRUE);
  gSystem->mkdir(gSystem->DirName(outFig), kTRUE);

  TFile *fD = TFile::Open(dataPath), *fM = TFile::Open(mcPath);
  if(!fD || fD->IsZombie() || !fM || fM->IsZombie()){ printf("ERROR: cannot open inputs\n"); return; }

  TH1D *hH[NHigh];
  for(int k = 0; k < NHigh; k++)
    hH[k] = bookPt("h_bPurity_fitHigh" + tag(fitHighVals[k]), Form("b purity, fit %.1f-%.1f GeV;#it{p}_{T}^{jet} [GeV];b purity", fitLo, fitHighVals[k]));
  TH1D *hNom = hH[0];

  TH1D *hRel = bookPt("h_sysRel_fitRangeHigh", "fit-range upper-edge systematic (max |relative deviation|, upper edge 4-5 GeV);#it{p}_{T}^{jet} [GeV];relative uncertainty");
  TH1D *hAbs = bookPt("h_sysAbs_fitRangeHigh", "fit-range upper-edge systematic (absolute);#it{p}_{T}^{jet} [GeV];absolute uncertainty");
  TH1D *hDev = bookPt("h_devRel_fitHigh4p0",   "relative deviation at fit 0-4 GeV;#it{p}_{T}^{jet} [GeV];(var - nom)/nom");

  printf("\n%-9s", "jet pT"); for(int k = 0; k < NHigh; k++) printf("  0-%.1f ", fitHighVals[k]); printf("   sysRel\n");
  for(int i = 0; i < NPtSys; i++){
    double lo = ptEdges[i], hi = ptEdges[i+1];
    TH1D *d05 = range05(projWin(fD, hName, lo, hi), Form("d05_%d", i));
    if(!fittable(d05)){ printf("%3.0f-%-5.0f  not fitted: %.0f data entries\n", lo, hi, d05->Integral()); continue; }
    Tpl t = buildTemplates(fM, lo, hi);
    TH1D *b = bWithGS(t, TMath::Min(1., t.fGSnat + gsShiftNominal));
    for(int k = 0; k < NHigh; k++){
      FitRes r = doFit(d05, {b, t.lc}, fitLo, fitHighVals[k]);
      hH[k]->SetBinContent(i+1, r.fb); hH[k]->SetBinError(i+1, r.efb);
    }
    double nom = hNom->GetBinContent(i+1), maxDev = 0.;
    for(int k = 0; k < NHigh; k++) maxDev = TMath::Max(maxDev, fabs(hH[k]->GetBinContent(i+1) - nom)/nom);
    hRel->SetBinContent(i+1, maxDev);
    hAbs->SetBinContent(i+1, maxDev*nom);
    hDev->SetBinContent(i+1, (hH[NHigh-1]->GetBinContent(i+1) - nom)/nom);
    printf("%3.0f-%-5.0f", lo, hi); for(int k = 0; k < NHigh; k++) printf("  %.4f", hH[k]->GetBinContent(i+1));
    printf("   %5.2f%%\n", 100*maxDev);
  }

  TFile *fo = TFile::Open(outRoot, "recreate");
  hNom->Write("h_bPurity_nominal");
  for(int k = 1; k < NHigh; k++) hH[k]->Write();
  hRel->Write(); hAbs->Write(); hDev->Write();
  TNamed info("info", Form("pp calo jets, 2-template ptRel fit, nominal %.0f-%.0f GeV, bGS share + %.3f, c-multiplier 1.0. "
                            "Systematic: upper edge %.1f-%.1f GeV (lower %.0f), max |relative deviation|. Lower edge: see "
                            "bPuritySys_fitRange_caloJets_pp.root. data: %s  MC: %s", fitLo, fitHi, gsShiftNominal,
                            fitHighVals[NHigh-1], fitHighVals[0], fitLo, gSystem->BaseName(dataPath), gSystem->BaseName(mcPath)));
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
  TH1D *hN = (TH1D*) hNom->Clone("hN"), *h45 = (TH1D*) hH[2]->Clone("h46"), *h40 = (TH1D*) hH[NHigh-1]->Clone("h40");
  styleH(hN,  "#000000", markFilledCircle, 1.2);
  styleH(h45, "#0072B2", markOpenSquare,   1.2);
  styleH(h40, "#D55E00", markOpenDiamond,  1.5);
  for(int i = 1; i <= NPtSys; i++){ h45->SetBinError(i, 0); h40->SetBinError(i, 0); }
  h45->Draw("P same"); h40->Draw("P same"); hN->Draw("E1 same");
  TLegend *leg = makeLegend(0.53, 0.70, 0.95, 0.88, 0.036);
  leg->AddEntry(hN,  "nominal: fit 0-5 GeV", "lp");
  leg->AddEntry(h45, Form("fit 0-%.1f GeV", fitHighVals[2]), "p");
  leg->AddEntry(h40, Form("fit 0-%.1f GeV", fitHighVals[NHigh-1]), "p");
  leg->Draw();
  TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.046);
  la.DrawLatex(0.21, 0.84, "pp 5.02 TeV, calo jets");
  la.SetTextSize(0.038);
  la.DrawLatex(0.21, 0.78, "2-template #it{p}_{T}^{rel} fit");

  pBot->cd();
  TH1D *r46 = (TH1D*) hH[2]->Clone("r46"), *r40 = (TH1D*) hDev->Clone("r40"), *rS = (TH1D*) hRel->Clone("rS");
  for(int i = 1; i <= NPtSys; i++){ double n = hNom->GetBinContent(i); r46->SetBinContent(i, (hH[2]->GetBinContent(i) - n)/n); }
  // deviations carry no error of their own; ROOT would otherwise draw sqrt(N)
  for(TH1D *h : {r46, r40, rS}){ h->Scale(100.); for(int i = 0; i <= h->GetNbinsX()+1; i++) h->SetBinError(i, 0.); }
  for(TH1D *h : {r46, r40}) blankUnfitted(h, hNom);
  styleH(r46, "#0072B2", markOpenSquare, 1.2); styleH(r40, "#D55E00", markOpenDiamond, 1.5);
  styleLine(rS, "#000000", 2);
  TH1D *rF = (TH1D*) rS->Clone("rF"); rF->Reset();
  double m = 1.3*TMath::Max(rS->GetMaximum(), 1.0);
  rF->SetMinimum(-m); rF->SetMaximum(1.8*m);   // headroom for the legend row
  styleRatioAxes(rF, "#it{p}_{T}^{jet} [GeV]", "#Delta / nom. [%]");
  rF->Draw("AXIS");
  drawUnderflowBand(rF->GetMinimum(), rF->GetMaximum(), -1e9);   // no label
  TLine *z = new TLine(ptEdges[0], 0, ptEdges[NPtSys], 0); z->SetLineStyle(2); z->SetLineColor(kGray+2); z->Draw();
  TH1D *rSm = (TH1D*) rS->Clone("rSm"); rSm->Scale(-1.);
  rS->Draw("HIST same"); rSm->Draw("HIST same");
  r46->Draw("P same"); r40->Draw("P same");
  TLegend *legB = makeLegend(0.19, 0.78, 0.95, 0.97, 0.075);
  legB->SetNColumns(3);
  legB->AddEntry(r46, Form("0-%.1f GeV", fitHighVals[2]), "p");
  legB->AddEntry(r40, Form("0-%.1f GeV", fitHighVals[NHigh-1]), "p");
  legB->AddEntry(rS, "syst. (max |#Delta|)", "l");
  legB->Draw();

  c->SaveAs(outFig);
  printf("written %s\n", outFig);
}
