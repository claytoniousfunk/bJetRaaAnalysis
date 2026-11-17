// Iterations systematic of the b-jet unfolding, extraction figure (2026-11-05; b-jet
// analogue of plotIterSys_PF.C). Top row: unfolded b-jet spectrum at N_opt - 1 and
// N_opt + 1 over N_opt, pp and the four PbPb classes (pp, 50-80 ... 0-10%), 80-200 GeV.
// Bottom row: the iterations uncertainty on R_Z(b) per class -- bin by bin (PbPb envelope
// (+) pp, as bJetsPerZSys_PF_PbPb.C combines them, uncorrelated) and as assigned (bin by bin since 2026-11-05; the
// flattening rule, a source with no pT trend replaced by its mean, applied before).
//
// Inputs: the raw b-jet systematics file (h_bIter_{nom,up,dn}_C<i>, saved since 2026-11-05;
// h_sysRel_iterations_b_PbPbOverPP_C<i>) and the pp unfolded file (h_bJetPt_unfolded and the
// N_opt -+ 1 spectra).
// Usage, from src/systematics/:  root -l -b -q plotBJetIterSys_PF.C
// Writes figures/systematics/iterations/bJetIterSys_PF.pdf

#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TSystem.h"
#include "../../headers/plotting/plotStyle.h"
#include "../../headers/plotting/gridPads.h"

void plotBJetIterSys_PF(const char *fSys = "../../rootFiles/systematics/bJetsPerZSys_PF_MBcaloMatched_JEUppCorr_noFlatJEU_fakeGeoOnly_csvV2loose.root",
                        const char *fPP  = "../../rootFiles/CorrectedBJetSpectra/Data/unfoldedBJetSpectrum_PFJets_pp_csvV2loose.root",
                        int nPP = 3, TString outSfx = "")   // outSfx: appended to the figure name (e.g. "_csvV2_2D7bins_nIterOpt", 2026-11-16)
{
  initPlotStyle();
  TFile *fs = TFile::Open(fSys), *fp = TFile::Open(fPP);
  if(!fs || !fp || fs->IsZombie() || fp->IsZombie()){ printf("ERROR: input missing\n"); return; }
  const char *sys[5] = {"pp", "C4", "C3", "C2", "C1"};
  const char *lab[5] = {"pp", "PbPb 50-80%", "PbPb 30-50%", "PbPb 10-30%", "PbPb 0-10%"};
  const double xLo = 80., xHi = 200.;

  // ratios (N -+ 1) / N per system
  TH1D *rUp[5] = {}, *rDn[5] = {}; int nOpt[5] = {nPP, 0, 0, 0, 0};
  for(int k = 0; k < 5; k++){
    TH1D *nom = nullptr, *up = nullptr, *dn = nullptr;
    if(k == 0){ fp->GetObject("h_bJetPt_unfolded", nom); fp->GetObject(Form("h_bJetPt_unfolded_N%d", nPP + 1), up); fp->GetObject(Form("h_bJetPt_unfolded_N%d", nPP - 1), dn); }
    else { fs->GetObject(Form("h_bIter_nom_%s", sys[k]), nom); fs->GetObject(Form("h_bIter_up_%s", sys[k]), up); fs->GetObject(Form("h_bIter_dn_%s", sys[k]), dn);
           if(nom) sscanf(nom->GetTitle(), "N = %d", &nOpt[k]); }
    if(!nom || !up){ printf("ERROR: spectra missing for %s (rerun bJetsPerZSys_PF_PbPb.C)\n", lab[k]); return; }
    rUp[k] = (TH1D*) up->Clone(Form("rUp%d", k)); rUp[k]->SetDirectory(nullptr); rUp[k]->Divide(nom);
    if(dn){ rDn[k] = (TH1D*) dn->Clone(Form("rDn%d", k)); rDn[k]->SetDirectory(nullptr); rDn[k]->Divide(nom); }
    for(TH1D *h : {rUp[k], rDn[k]}) if(h) for(int b = 1; b <= h->GetNbinsX(); b++) h->SetBinError(b, 0.);
  }
  // pp relative envelope per bin
  auto env = [&](int k, int b){ double e = fabs(rUp[k]->GetBinContent(b) - 1.); if(rDn[k]) e = TMath::Max(e, fabs(rDn[k]->GetBinContent(b) - 1.)); return e; };

  printf("iterations, b jets: (N-1)/N, (N+1)/N - 1 (%%) per bin, and on R_Z(b): bin by bin | assigned\n");
  TCanvas *c = new TCanvas("cBIt", "", 2000, 1000);
  GridPads g = makeGridPads(c, 5, 2, 0.07, 0.012, 0.11, 0.07);
  TH1D *lgU = nullptr, *lgD = nullptr, *lgB = nullptr, *lgA = nullptr;
  for(int ix = 0; ix < 5; ix++){
    const int k = ix;
    // top: ratios
    g.at(ix, 0)->cd();
    TH1F *fr = gPad->DrawFrame(xLo, 0.70, xHi, 1.30); g.styleFrame(fr, ix, 0, 0.070); fr->GetXaxis()->SetNdivisions(505);
    TLine *one = new TLine(xLo, 1., xHi, 1.); one->SetLineStyle(2); one->SetLineWidth(2); one->SetLineColor(kGray + 2); one->Draw();
    styleH(rUp[k], okabeHex[5], markFilledSquare, 1.9); rUp[k]->SetLineWidth(3); rUp[k]->Draw("P X0 same");
    if(rDn[k]){ styleH(rDn[k], okabeHex[1], markOpenCircle, 1.9); rDn[k]->SetLineWidth(3); rDn[k]->Draw("P X0 same"); }
    if(ix == 0){ lgU = rUp[k]; lgD = rDn[k]; }
    g.text(ix, 0, 0.10, lab[k], 0.070);
    g.text(ix, 0, 0.20, Form("#it{N}_{opt} = %d", nOpt[k]), 0.060);
    gPad->RedrawAxis();
    // bottom: on R_Z(b)
    g.at(ix, 1)->cd();
    TH1F *f2 = gPad->DrawFrame(xLo, 0., xHi, 0.40); g.styleFrame(f2, ix, 1, 0.070); f2->GetXaxis()->SetNdivisions(505);
    if(k == 0){ g.text(ix, 1, 0.12, "(pp enters each class", 0.055); g.text(ix, 1, 0.21, "in quadrature)", 0.055); gPad->RedrawAxis(); continue; }
    TH1D *bb = (TH1D*) rUp[k]->Clone(Form("bb%d", k)); bb->Reset();
    for(int b = 1; b <= bb->GetNbinsX(); b++){ const double lo = bb->GetBinLowEdge(b); if(lo < xLo - 1e-6 || lo >= xHi - 1e-6) continue;
      const int bp = rUp[0]->FindBin(bb->GetBinCenter(b)); bb->SetBinContent(b, sqrt(pow(env(k, b), 2) + pow(env(0, bp), 2))); }
    TH1D *as = nullptr; fs->GetObject(Form("h_sysRel_iterations_b_PbPbOverPP_%s", sys[k]), as);
    if(!as){ printf("ERROR: assigned iterations source missing for %s\n", lab[k]); return; }
    as = (TH1D*) as->Clone(Form("as%d", k)); as->SetDirectory(nullptr);
    for(int b = 1; b <= as->GetNbinsX(); b++){ const double lo = as->GetBinLowEdge(b); if(lo < xLo - 1e-6 || lo >= xHi - 1e-6){ as->SetBinContent(b, 0.); } as->SetBinError(b, 0.); }
    styleH(bb, hexData, markFilledCircle, 1.9); bb->SetLineWidth(3); bb->Draw("P X0 same");
    as->SetLineColor(TColor::GetColor(okabeHex[5])); as->SetLineWidth(4); as->SetMarkerSize(0); as->Draw("HIST same");
    if(ix == 1){ lgB = bb; lgA = as; }
    printf("  %-12s N %d  (N-1)/N-1:", lab[k], nOpt[k]);
    for(int b = 1; b <= rUp[k]->GetNbinsX(); b++){ const double lo = rUp[k]->GetBinLowEdge(b); if(lo < xLo - 1e-6 || lo >= xHi - 1e-6) continue;
      printf(" %+.1f", rDn[k] ? 100.*(rDn[k]->GetBinContent(b) - 1.) : 0.); }
    printf("  (N+1)/N-1:");
    for(int b = 1; b <= rUp[k]->GetNbinsX(); b++){ const double lo = rUp[k]->GetBinLowEdge(b); if(lo < xLo - 1e-6 || lo >= xHi - 1e-6) continue; printf(" %+.1f", 100.*(rUp[k]->GetBinContent(b) - 1.)); }
    printf("  R_Z bin by bin:");
    for(int b = 1; b <= bb->GetNbinsX(); b++){ const double lo = bb->GetBinLowEdge(b); if(lo < xLo - 1e-6 || lo >= xHi - 1e-6) continue; printf(" %.1f", 100.*bb->GetBinContent(b)); }
    printf("  | assigned %.1f\n", 100.*as->GetBinContent(as->FindBin(90.)));
    gPad->RedrawAxis();
  }
  g.at(0, 0)->cd();
  TLegend *l1 = makeLegend(g.frameLeft(0) + 0.05, 0.08 + g.frameBottom(0), 0.97, 0.30 + g.frameBottom(0), 0.060*g.sizeScale());
  l1->AddEntry(lgU, "#it{N}_{opt} + 1", "p"); if(lgD) l1->AddEntry(lgD, "#it{N}_{opt} #minus 1", "p"); l1->Draw();
  g.at(0, 1)->cd();
  TLegend *l2 = makeLegend(g.frameLeft(0) + 0.05, 0.40, 0.97, 0.62, 0.060*g.sizeScale());
  l2->AddEntry(lgB, "bin by bin", "p"); l2->AddEntry(lgA, "assigned", "l"); l2->Draw();
  c->cd();
  TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.034); t.SetTextAlign(21);
  t.DrawLatex(0.07 + (1. - 0.07 - 0.012)/2., 0.015, "Muon-tagged #it{b}-jet #it{p}_{T} (GeV)");
  t.SetTextAngle(90);
  t.DrawLatex(0.022, 0.11 + 0.82*0.75, "Unfolded / #it{N}_{opt}");
  t.DrawLatex(0.022, 0.11 + 0.82*0.25, "Iterations on #it{R}_{Z}^{b}");
  g.cmsLabel(c, "Internal", "PbPb and pp, #sqrt{#it{s}_{NN}} = 5.02 TeV", 1.2);
  gSystem->mkdir("../../figures/systematics/iterations", kTRUE);
  savePdfTight(c, "../../figures/systematics/iterations/bJetIterSys_PF" + outSfx + ".pdf");
}
