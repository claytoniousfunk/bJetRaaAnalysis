// Generator-level effect of the jet energy scale uncertainty on the unfolded PF inclusive-jet
// spectrum, pp and the four PbPb classes in one 1x5 shared-border row, peripheral to central,
// CMS style (AN figure fig:jeuShift, 2026-11-06; replaces the per-system figures of
// unfoldPYTHIAJetPt_JEUShift.C / unfoldPYTHIAHYDJETJetPt_JEUShift.C, which labelled the PbPb
// PF jets "calo jets"). Draws the saved outputs of those two macros, no unfolding here:
//   rootFiles/JEU/unfold_PYTHIA_PFJets_JEUShift.root        gen_{nominal,up,down}_N<n>
//   rootFiles/JEU/unfold_PYTHIAHYDJET_PFJets_JEUShift.root  gen_{nominal,up,down}_N<n>_C<i>
// (reco pT shifted by +- the AK4PF JEC uncertainty, unfolded with the nominal response;
// 5 GeV bins summed into the analysis axis, ratio JEU-shifted / nominal at the nominal N_opt).
// Ratio errors: the shifted spectrum's statistics only (same jets as the nominal).
//
// Usage, from src/plots/unfolding/:  root -l -b -q plotJEUShiftGen_PF.C
// Output: figures/unfolding/JEUShiftGen_PF.pdf (inclusive), JEUShiftGen_bJets_PF.pdf (muon-tagged b jets,
//         from rootFiles/JEU/jeuShiftMC_bJets_PF.root of src/unfoldTest/jeuShiftMC_bJets_PF.C)
#include "TFile.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TVectorD.h"
#include "TSystem.h"
#include "../../../headers/plotting/plotStyle.h"
#include "../../../headers/plotting/gridPads.h"

void plotJEUShiftGen_PF(TString bTag = "")   // bTag "_7bins" (2026-11-16): only the b-jet figure, from jeuShiftMC_bJets_PF_7bins.root
{
  initPlotStyle();
  const char *dir = "../../../rootFiles/JEU/";
  const char *lab[5] = {"pp", "PbPb 50-80%", "PbPb 30-50%", "PbPb 10-30%", "PbPb 0-10%"};
  if(bTag == ""){   // inclusive figure
  TFile *fpp = TFile::Open(Form("%sunfold_PYTHIA_PFJets_JEUShift.root", dir)), *fpb = TFile::Open(Form("%sunfold_PYTHIAHYDJET_PFJets_JEUShift.root", dir));
  if(!fpp || !fpb){ printf("ERROR: an input did not open\n"); return; }
  const double edges[] = {80, 90, 100, 110, 130, 150, 180, 200, 240, 280, 350, 500};
  const int nE = sizeof(edges)/sizeof(edges[0]) - 1;
  const int cls[5]  = {0, 4, 3, 2, 1};
  const int nOpt[5] = {2, 2, 4, 6, 8};   // inclusive N_opt, Table tab:nIter (pp, 50-80 ... 0-10%)

  auto get = [&](int p, const char *v){
    TH1D *h = nullptr;
    if(cls[p] == 0) fpp->GetObject(Form("gen_%s_N%d", v, nOpt[p]), h);
    else            fpb->GetObject(Form("gen_%s_N%d_C%d", v, nOpt[p], cls[p]), h);
    if(!h){ printf("ERROR: gen_%s N%d for %s missing\n", v, nOpt[p], lab[p]); return (TH1D*) nullptr; }
    TH1D *r = (TH1D*) h->Rebin(nE, Form("rb_%s_%d", v, p), edges); r->SetDirectory(nullptr);
    return r;
  };

  TCanvas *cv = new TCanvas("cJEU", "", 2000, 650);
  GridPads g = makeGridPads(cv, 5, 1, 0.065, 0.012, 0.17, 0.10);
  TH1D *legU = nullptr, *legD = nullptr;
  printf("JEU-shifted / nominal at N_opt, %%:\n");
  for(int p = 0; p < 5; p++){
    TH1D *n = get(p, "nominal"), *u = get(p, "up"), *d = get(p, "down");
    if(!n || !u || !d) return;
    TH1D *ru = (TH1D*) u->Clone(Form("ru_%d", p)), *rd = (TH1D*) d->Clone(Form("rd_%d", p));
    for(int b = 1; b <= nE; b++){
      const double y = n->GetBinContent(b);
      for(TH1D *r : {ru, rd}){ const TH1D *s = r == ru ? u : d;
        r->SetBinContent(b, y > 0. ? s->GetBinContent(b)/y : 0.); r->SetBinError(b, y > 0. ? s->GetBinError(b)/y : 0.); }
    }
    printf("  %-12s up", lab[p]); for(int b = 1; b <= nE; b++) printf(" %+5.1f", 100.*(ru->GetBinContent(b) - 1.));
    printf("\n  %-12s dn", ""); for(int b = 1; b <= nE; b++) printf(" %+5.1f", 100.*(rd->GetBinContent(b) - 1.));
    printf("\n");

    g.at(p, 0)->cd();
    TH1F *fr = gPad->DrawFrame(80., 0.50, 500., 1.60); g.styleFrame(fr, p, 0, 0.075);
    fr->GetXaxis()->SetLabelSize(0.075*g.padW(0)/g.padW(p)); fr->GetXaxis()->SetNdivisions(505);
    TLine *one = new TLine(80., 1., 500., 1.); one->SetLineStyle(2); one->SetLineWidth(2); one->SetLineColor(kGray + 2); one->Draw();
    styleH(ru, okabeHex[6], markFilledSquare, 1.8);  ru->SetLineWidth(3); ru->Draw("E1 X0 same");
    styleH(rd, okabeHex[5], markFilledDiamond, 2.4); rd->SetLineWidth(3); rd->Draw("E1 X0 same");
    if(!legU){ legU = ru; legD = rd; }
    g.text(p, 0, 0.10, lab[p], 0.075);
    { TLatex tn; tn.SetNDC(); tn.SetTextFont(42); tn.SetTextAlign(31); tn.SetTextSize(0.062*g.sizeScale());
      tn.DrawLatex(1. - gPad->GetRightMargin() - 0.05, g.frameBottom(0) + 0.06, Form("#it{N}_{opt} = %d", nOpt[p])); }
    gPad->RedrawAxis();
  }
  g.at(0, 0)->cd();   // pp panel
  TLegend *lg = makeLegend(g.frameLeft(0) + 0.03, 0.625, 0.80, 0.775, 0.055*g.sizeScale());   // between the label and the up points
  lg->AddEntry(legU, "JEU up", "pe"); lg->AddEntry(legD, "JEU down", "pe"); lg->Draw();
  g.titles(cv, "Generator-level jet #it{p}_{T} (GeV)", "Shifted / nominal", 0.060, 0.020, 0.030);
  g.cmsLabel(cv, "Simulation Internal", "Inclusive anti-#it{k}_{T} #it{R} = 0.4 PF jets, |#it{#eta}^{jet}| < 1.6, unfolded,  #sqrt{#it{s}_{NN}} = 5.02 TeV", 1.4);
  gSystem->mkdir("../../../figures/unfolding", kTRUE);
  savePdfTight(cv, "../../../figures/unfolding/JEUShiftGen_PF.pdf");

  }
  // muon-tagged b jets (2026-11-06): jeuShiftMC_bJets_PF.C, the chain's analytic shift applied to the
  // reco projection of the b response, unfolded at the b-jet N_opt, reported bins 80-200 GeV
  TFile *fb = TFile::Open(Form("%sjeuShiftMC_bJets_PF%s.root", dir, bTag.Data()));
  if(!fb){ printf("ERROR: jeuShiftMC_bJets_PF.root missing (run src/unfoldTest/jeuShiftMC_bJets_PF.C)\n"); return; }
  TVectorD *nb = nullptr; fb->GetObject("optimalN", nb);
  const char *tagB[5] = {"pp", "C4", "C3", "C2", "C1"};
  const int idxB[5] = {0, 4, 3, 2, 1};   // optimalN order pp, C1 .. C4
  TCanvas *cb = new TCanvas("cJEUb", "", 2000, 650);
  GridPads gb = makeGridPads(cb, 5, 1, 0.065, 0.012, 0.17, 0.10);
  TH1D *lbU = nullptr, *lbD = nullptr;
  for(int p = 0; p < 5; p++){
    TH1D *ru = nullptr, *rd = nullptr;
    fb->GetObject(Form("h_jeuRatio_up_%s", tagB[p]), ru); fb->GetObject(Form("h_jeuRatio_down_%s", tagB[p]), rd);
    if(!ru || !rd || !nb){ printf("ERROR: b-jet JEU ratios missing for %s\n", tagB[p]); return; }
    const int N = (int) (*nb)[idxB[p]];
    gb.at(p, 0)->cd();
    TH1F *fr = gPad->DrawFrame(80., 0.50, 200., 1.60); gb.styleFrame(fr, p, 0, 0.075);
    fr->GetXaxis()->SetLabelSize(0.075*gb.padW(0)/gb.padW(p)); fr->GetXaxis()->SetNdivisions(505);
    TLine *one = new TLine(80., 1., 200., 1.); one->SetLineStyle(2); one->SetLineWidth(2); one->SetLineColor(kGray + 2); one->Draw();
    styleH(ru, okabeHex[6], markFilledSquare, 1.8);  ru->SetLineWidth(3); ru->Draw("P same");
    styleH(rd, okabeHex[5], markFilledDiamond, 2.4); rd->SetLineWidth(3); rd->Draw("P same");
    if(!lbU){ lbU = ru; lbD = rd; }
    gb.text(p, 0, 0.10, lab[p], 0.075);
    { TLatex tn; tn.SetNDC(); tn.SetTextFont(42); tn.SetTextAlign(31); tn.SetTextSize(0.062*gb.sizeScale());
      tn.DrawLatex(1. - gPad->GetRightMargin() - 0.05, gb.frameBottom(0) + 0.06, Form("#it{N}_{opt} = %d", N)); }
    gPad->RedrawAxis();
  }
  gb.at(0, 0)->cd();
  TLegend *lb = makeLegend(gb.frameLeft(0) + 0.03, 0.625, 0.80, 0.775, 0.055*gb.sizeScale());
  lb->AddEntry(lbU, "JEU up", "p"); lb->AddEntry(lbD, "JEU down", "p"); lb->Draw();
  gb.titles(cb, "Generator-level jet #it{p}_{T} (GeV)", "Shifted / nominal", 0.060, 0.020, 0.030);
  gb.cmsLabel(cb, "Simulation Internal", "Muon-tagged #it{b} jets, anti-#it{k}_{T} #it{R} = 0.4 PF, |#it{#eta}^{jet}| < 1.6, unfolded,  #sqrt{#it{s}_{NN}} = 5.02 TeV", 1.4);
  savePdfTight(cb, "../../../figures/unfolding/JEUShiftGen_bJets_PF" + bTag + ".pdf");
}
