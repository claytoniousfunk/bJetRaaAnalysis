// Jet energy scale (JEU) and resolution (JER) systematics on the pp CALO-jet
// b purity.
//
// The data are left alone; the MC templates are re-binned in a varied jet pT,
// which moves jets between jet-pT windows and so changes the template shapes
// and the MC flavour fractions the fit uses. The variants are filled by
// PYTHIA_scan.C into the _T1/_T2/_T3 template histograms (templateIndexNames):
//   T2 / T3   jet pT x (1 +/- JEU uncertainty)
//   T1        jet pT x Gaus(1, 0.663 x JER(pT))
// Systematics, per jet-pT window, relative to nominal (T0):
//   JEU = max(|up - nom|, |down - nom|) / nom
//   JER = |smeared - nom| / nom
//
// *** CAVEAT: PF INPUTS FOR CALO JETS ***
// PYTHIA_scan.C takes the JEU from Spring18_ppRef5TeV_V6_MC_Uncertainty_AK4PF
// for every jet collection (there is no Spring18 pp AK4Calo uncertainty file;
// only Autumn18_HI_V8_MC_Uncertainty_AK4Calo), and the JER from a single PF
// parametrisation. These numbers are therefore the calo-jet fit's response to
// PF-sized shifts, not a calo-specific JEU/JER. Fixing that needs the calo
// uncertainty/resolution in the scan and a rescan.
//
// Other settings nominal: fit 0-5 GeV, bGS share + 0.175, c-multiplier 1.0.
// Fit, templates and inputs: headers/functions/caloJetBPurityFit_pp.h.
//
// Output:
//   rootFiles/systematics/bPuritySys_JEU_caloJets_pp.root
//     h_bPurity_nominal, h_bPurity_JEUShiftUp, h_bPurity_JEUShiftDown,
//     h_sysRel_JEU (the systematic), h_sysAbs_JEU,
//     h_sysRel_JEU_up, h_sysRel_JEU_down (signed)
//   rootFiles/systematics/bPuritySys_JER_caloJets_pp.root
//     h_bPurity_nominal, h_bPurity_JERSmear,
//     h_sysRel_JER (the systematic), h_sysAbs_JER, h_devRel_JER (signed)
//   figures/systematics/JEU_caloJets_pp.pdf, figures/systematics/JER_caloJets_pp.pdf
//
// Usage: root -l -b -q systematics_JEUJER_caloJets_pp.C
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

const double gsShiftNominal = 0.175;

// Systematics are evaluated in every window, underflow (50-80 GeV) included,
// since the unfolding needs them there too; windows without enough data (see
// fittable() in the header) are skipped and left empty.
const int NPtSys = NPt;

const char *outDir   = "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/systematics";
const char *figDir   = "/home/clayton/Analysis/code/bJetRaaAnalysis/figures/systematics";

TH1D* bookPt(const char *name, const char *title)
{
  TH1D *h = new TH1D(name, title, NPtSys, ptEdges); h->SetDirectory(nullptr); return h;
}

// purity in every window with templates from MC variant tplIdx
TH1D* purityVsPt(TFile *fD, TFile *fM, int tplIdx, const char *name, const char *title)
{
  TH1D *h = bookPt(name, title);
  for(int i = 0; i < NPtSys; i++){
    double lo = ptEdges[i], hi = ptEdges[i+1];
    TH1D *d05 = range05(projWin(fD, hName, lo, hi), Form("d05_%d_%d", tplIdx, i));
    if(!fittable(d05)) continue;   // left empty
    Tpl t = buildTemplates(fM, lo, hi, 1.0, tplIdx);
    TH1D *b = bWithGS(t, TMath::Min(1., t.fGSnat + gsShiftNominal));
    FitRes r = doFit(d05, {b, t.lc});
    h->SetBinContent(i+1, r.fb); h->SetBinError(i+1, r.efb);
  }
  return h;
}

TH1D* relDev(TH1D *var, TH1D *nom, const char *name, const char *title)
{
  TH1D *h = bookPt(name, title);
  for(int i = 1; i <= NPtSys; i++) if(nom->GetBinContent(i) > 0.) h->SetBinContent(i, (var->GetBinContent(i) - nom->GetBinContent(i))/nom->GetBinContent(i));
  return h;
}

// two-panel figure: purities (top), signed deviations + symmetric band (bottom)
void drawSys(TH1D *nom, std::vector<TH1D*> vars, std::vector<const char*> varLabels,
             std::vector<TH1D*> devs, TH1D *sysRel, const char *sysLabel, const TString &pdf)
{
  const char *hexV[2] = {"#0072B2", "#D55E00"};
  const int   mkV[2]  = {markOpenSquare, markOpenDiamond};
  TCanvas *c = new TCanvas(Form("c_%s", sysLabel), "", 700, 800);
  TPad *pTop, *pBot; splitPads(pTop, pBot);
  pTop->cd();
  TH1D *fr = (TH1D*) nom->Clone(Form("fr_%s", sysLabel)); fr->Reset(); fr->SetTitle("");
  fr->SetMinimum(0.2); fr->SetMaximum(0.75);
  fr->GetXaxis()->SetLabelSize(0);
  fr->GetYaxis()->SetTitle("b-jet purity"); fr->GetYaxis()->SetTitleSize(0.055);
  fr->GetYaxis()->SetTitleOffset(1.35); fr->GetYaxis()->SetLabelSize(0.045);
  fr->Draw("AXIS");
  drawUnderflowBand(fr->GetMinimum(), fr->GetMaximum(), 0.25);
  TH1D *hN = (TH1D*) nom->Clone(Form("hN_%s", sysLabel)); styleH(hN, "#000000", markFilledCircle, 1.2);
  TLegend *leg = makeLegend(0.53, 0.70, 0.95, 0.88, 0.036);
  leg->AddEntry(hN, "nominal", "lp");
  for(size_t k = 0; k < vars.size(); k++){
    TH1D *v = (TH1D*) vars[k]->Clone(Form("v%zu_%s", k, sysLabel));
    styleH(v, hexV[k], mkV[k], k == 1 ? 1.5 : 1.2);
    for(int i = 1; i <= NPtSys; i++) v->SetBinError(i, 0.);
    v->Draw("P same"); leg->AddEntry(v, varLabels[k], "p");
  }
  hN->Draw("E1 same"); leg->Draw();
  TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.046);
  la.DrawLatex(0.21, 0.84, "pp 5.02 TeV, calo jets");
  la.SetTextSize(0.038);
  la.DrawLatex(0.21, 0.78, "2-template #it{p}_{T}^{rel} fit");
  la.SetTextSize(0.032); la.SetTextColor(kGray+2);
  la.DrawLatex(0.21, 0.07, "MC shifts use the AK4PF uncertainty / resolution");   // purities sit above ~0.42

  pBot->cd();
  TH1D *rS = (TH1D*) sysRel->Clone(Form("rS_%s", sysLabel));
  rS->Scale(100.); for(int i = 0; i <= NPtSys+1; i++) rS->SetBinError(i, 0.);
  styleLine(rS, "#000000", 2);
  std::vector<TH1D*> rD;
  double m = TMath::Max(rS->GetMaximum(), 1.0);
  for(size_t k = 0; k < devs.size(); k++){
    TH1D *r = (TH1D*) devs[k]->Clone(Form("rD%zu_%s", k, sysLabel));
    r->Scale(100.); for(int i = 0; i <= NPtSys+1; i++) r->SetBinError(i, 0.);
    blankUnfitted(r, nom);
    styleH(r, hexV[k], mkV[k], k == 1 ? 1.5 : 1.2); rD.push_back(r);
    for(int i = 1; i <= NPtSys; i++) if(r->GetBinContent(i) > -900.) m = TMath::Max(m, fabs(r->GetBinContent(i)));   // skip blanked bins
  }
  m *= 1.3;
  TH1D *rF = (TH1D*) rS->Clone(Form("rF_%s", sysLabel)); rF->Reset();
  rF->SetMinimum(-m); rF->SetMaximum(1.8*m);   // headroom for the legend row
  styleRatioAxes(rF, "#it{p}_{T}^{jet} [GeV]", "#Delta / nom. [%]");
  rF->Draw("AXIS");
  drawUnderflowBand(rF->GetMinimum(), rF->GetMaximum(), -1e9);   // no label
  TLine *z = new TLine(ptEdges[0], 0, ptEdges[NPtSys], 0); z->SetLineStyle(2); z->SetLineColor(kGray+2); z->Draw();
  TH1D *rSm = (TH1D*) rS->Clone(Form("rSm_%s", sysLabel)); rSm->Scale(-1.);
  rS->Draw("HIST same"); rSm->Draw("HIST same");
  for(auto r : rD) r->Draw("P same");
  TLegend *legB = makeLegend(0.19, 0.78, 0.95, 0.97, 0.075);
  legB->SetNColumns((int) rD.size() + 1);
  for(size_t k = 0; k < rD.size(); k++) legB->AddEntry(rD[k], varLabels[k], "p");
  legB->AddEntry(rS, "syst.", "l");
  legB->Draw();
  c->SaveAs(pdf);
  delete c;
}

void systematics_JEUJER_caloJets_pp()
{
  initPlotStyle();
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  gSystem->mkdir(outDir, kTRUE); gSystem->mkdir(figDir, kTRUE);

  TFile *fD = TFile::Open(dataPath), *fM = TFile::Open(mcPath);
  if(!fD || fD->IsZombie() || !fM || fM->IsZombie()){ printf("ERROR: cannot open inputs\n"); return; }

  const char *ttl = ";#it{p}_{T}^{jet} [GeV];b purity";
  TH1D *hNom  = purityVsPt(fD, fM, 0, "h_bPurity_nominal",      Form("b purity, nominal MC jet pT%s", ttl));
  TH1D *hJER  = purityVsPt(fD, fM, 1, "h_bPurity_JERSmear",     Form("b purity, MC jet pT JER-smeared%s", ttl));
  TH1D *hUp   = purityVsPt(fD, fM, 2, "h_bPurity_JEUShiftUp",   Form("b purity, MC jet pT JEU up%s", ttl));
  TH1D *hDown = purityVsPt(fD, fM, 3, "h_bPurity_JEUShiftDown", Form("b purity, MC jet pT JEU down%s", ttl));

  const char *dt = ";#it{p}_{T}^{jet} [GeV];(var - nom)/nom";
  TH1D *dUp   = relDev(hUp,   hNom, "h_sysRel_JEU_up",   Form("relative deviation, JEU up%s", dt));
  TH1D *dDown = relDev(hDown, hNom, "h_sysRel_JEU_down", Form("relative deviation, JEU down%s", dt));
  TH1D *dJER  = relDev(hJER,  hNom, "h_devRel_JER",      Form("relative deviation, JER smear%s", dt));

  TH1D *sJEU  = bookPt("h_sysRel_JEU", "JEU systematic (max |relative deviation| of up/down);#it{p}_{T}^{jet} [GeV];relative uncertainty");
  TH1D *aJEU  = bookPt("h_sysAbs_JEU", "JEU systematic (absolute);#it{p}_{T}^{jet} [GeV];absolute uncertainty");
  TH1D *sJER  = bookPt("h_sysRel_JER", "JER systematic (|relative deviation|);#it{p}_{T}^{jet} [GeV];relative uncertainty");
  TH1D *aJER  = bookPt("h_sysAbs_JER", "JER systematic (absolute);#it{p}_{T}^{jet} [GeV];absolute uncertainty");

  printf("\n%-9s %8s %8s %8s %8s   %7s %7s   %6s %6s\n", "jet pT", "nominal", "JEU up", "JEU dn", "JER", "dUp", "dDown", "JEU", "JER");
  for(int i = 1; i <= NPtSys; i++){
    double nom = hNom->GetBinContent(i);
    if(nom <= 0.){ printf("%3.0f-%-5.0f  not fitted (too little data)\n", ptEdges[i-1], ptEdges[i]); continue; }
    double sj = TMath::Max(fabs(dUp->GetBinContent(i)), fabs(dDown->GetBinContent(i)));
    double sr = fabs(dJER->GetBinContent(i));
    sJEU->SetBinContent(i, sj); aJEU->SetBinContent(i, sj*nom);
    sJER->SetBinContent(i, sr); aJER->SetBinContent(i, sr*nom);
    printf("%3.0f-%-5.0f %8.4f %8.4f %8.4f %8.4f   %+6.2f%% %+6.2f%%   %5.2f%% %5.2f%%\n", ptEdges[i-1], ptEdges[i],
           nom, hUp->GetBinContent(i), hDown->GetBinContent(i), hJER->GetBinContent(i),
           100*dUp->GetBinContent(i), 100*dDown->GetBinContent(i), 100*sj, 100*sr);
  }

  TString note = Form("pp calo jets, 2-template ptRel fit %.0f-%.0f GeV, bGS share + %.3f, c-multiplier 1.0; data unchanged, "
                      "MC templates in varied jet pT (PYTHIA_scan.C _T1/_T2/_T3). CAVEAT: the MC shifts use the AK4PF JEU "
                      "uncertainty and PF JER parametrisation for calo jets. data: %s  MC: %s",
                      fitLo, fitHi, gsShiftNominal, gSystem->BaseName(dataPath), gSystem->BaseName(mcPath));
  { TFile *fo = TFile::Open(Form("%s/bPuritySys_JEU_caloJets_pp.root", outDir), "recreate");
    hNom->Write(); hUp->Write(); hDown->Write(); sJEU->Write(); aJEU->Write(); dUp->Write(); dDown->Write();
    TNamed info("info", (TString("JEU systematic = max(|up-nom|,|down-nom|)/nom. ") + note).Data()); info.Write(); fo->Close(); }
  { TFile *fo = TFile::Open(Form("%s/bPuritySys_JER_caloJets_pp.root", outDir), "recreate");
    hNom->Write(); hJER->Write(); sJER->Write(); aJER->Write(); dJER->Write();
    TNamed info("info", (TString("JER systematic = |smeared-nom|/nom. ") + note).Data()); info.Write(); fo->Close(); }
  printf("\nwritten %s/bPuritySys_JEU_caloJets_pp.root and bPuritySys_JER_caloJets_pp.root\n", outDir);

  drawSys(hNom, {hDown, hUp}, {"JEU down", "JEU up"}, {dDown, dUp}, sJEU, "JEU", Form("%s/JEU_caloJets_pp.pdf", figDir));
  drawSys(hNom, {hJER}, {"JER smear"}, {dJER}, sJER, "JER", Form("%s/JER_caloJets_pp.pdf", figDir));
  printf("written %s/JEU_caloJets_pp.pdf and JER_caloJets_pp.pdf\n", figDir);
}
