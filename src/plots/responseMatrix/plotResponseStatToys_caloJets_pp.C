// How the response-matrix MC-statistics systematic of the pp CALO-jet b-jet
// unfolding is derived. Reads what unfoldBJetSpectrum_caloJets_pp.C stores:
//   h_response_neff, h_response_relErr   per-cell effective entries (c/sigma)^2
//                                         and sigma/c of the response
//   h_respToys_values                     every toy's unfolded dN/dpT
//   h_bJetPt_unfolded                     the nominal unfolded spectrum
//   h_unf_sysRel_responseStat             the systematic itself (RMS / nominal)
//
// The method (see unfoldBJetSpectrum_caloJets_pp.C): each response cell is
// fluctuated as Poisson(neff) x sigma^2/c, the prior / floor / projections /
// RooUnfoldResponse are rebuilt from the fluctuated matrix, the data are
// unfolded at the nominal iteration count, and the spread over toys is the
// systematic.
//
// Figures, figures/systematics/responseStatToys/:
//   1_responseInputStats.pdf   neff and sigma/c per cell: what sets each cell's
//                              fluctuation, and where the matrix is thin
//   2_toyRatios.pdf            toy / nominal vs jet pT: individual toys, the
//                              +/- RMS band (= the systematic) and the toy mean
//   3_toyDistributions.pdf     per reported bin, the distribution of toy /
//                              nominal over all toys with its RMS and mean
//   4_toyConvergence.pdf       running RMS vs number of toys
//
// Usage: root -l -b -q plotResponseStatToys_caloJets_pp.C
// Run from: src/plots/responseMatrix/

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TBox.h"
#include "TLine.h"
#include "TLatex.h"
#include "TGraph.h"
#include "TSystem.h"
#include "TColor.h"
#include "TStyle.h"
#include "../../../headers/plotting/plotStyle.h"
#include "../../../headers/functions/caloJetBPurityFit_pp.h"   // ptEdges, NPt, ptReportMin

const char *inPath = "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/CorrectedBJetSpectra/Data/unfoldedBJetSpectrum_caloJets_pp.root";
const char *outDir = "/home/clayton/Analysis/code/bJetRaaAnalysis/figures/systematics/responseStatToys";
const int nShowToys = 100;   // individual toy curves drawn in figure 2

bool reported(int i){ return ptEdges[i-1] >= ptReportMin - 1e-6; }

// cell map in bin-index space, one equal-size cell per bin, ranges as labels
TH2D* indexMap(TH2D *h, const char *name)
{
  TH2D *o = new TH2D(name, "", NPt, 0, NPt, NPt, 0, NPt); o->SetDirectory(nullptr);
  for(int ix = 1; ix <= NPt; ix++) for(int iy = 1; iy <= NPt; iy++) o->SetBinContent(ix, iy, h->GetBinContent(ix, iy));
  for(int i = 1; i <= NPt; i++){
    TString lab = Form("%.0f-%.0f", ptEdges[i-1], ptEdges[i]);
    o->GetXaxis()->SetBinLabel(i, lab); o->GetYaxis()->SetBinLabel(i, lab);
  }
  o->GetXaxis()->SetTitle("reco #it{p}_{T}^{jet} [GeV]"); o->GetYaxis()->SetTitle("gen #it{p}_{T}^{jet} [GeV]");
  o->GetXaxis()->SetLabelSize(0.040); o->GetYaxis()->SetLabelSize(0.040);
  o->GetXaxis()->SetTitleSize(0.048); o->GetYaxis()->SetTitleSize(0.048);
  o->GetXaxis()->SetTitleOffset(1.55); o->GetYaxis()->SetTitleOffset(1.85);
  return o;
}

void plotResponseStatToys_caloJets_pp()
{
  initPlotStyle();
  gSystem->mkdir(outDir, kTRUE);

  TFile *f = TFile::Open(inPath);
  if(!f || f->IsZombie()){ printf("ERROR: cannot open %s -- run unfoldBJetSpectrum_caloJets_pp.C\n", inPath); return; }
  TH2D *hNeff = nullptr, *hRelE = nullptr, *hToy = nullptr, *hResp = nullptr;
  TH1D *hNom = nullptr, *hSys = nullptr;
  f->GetObject("h_response_neff", hNeff);
  f->GetObject("h_response_relErr", hRelE);
  f->GetObject("h_respToys_values", hToy);
  f->GetObject("h_bJetPt_unfolded", hNom);
  f->GetObject("h_unf_sysRel_responseStat", hSys);
  f->GetObject("h_response", hResp);
  if(!hNeff || !hRelE || !hToy || !hNom || !hSys || !hResp){ printf("ERROR: missing objects in %s (rerun the unfolding)\n", inPath); return; }
  const int nToys = hToy->GetNbinsX();

  // toy / nominal, per toy and bin
  std::vector<std::vector<double>> r(NPt+1, std::vector<double>(nToys, 0.));
  for(int i = 1; i <= NPt; i++){
    double n = hNom->GetBinContent(i);
    for(int t = 0; t < nToys; t++) r[i][t] = n > 0. ? hToy->GetBinContent(t+1, i)/n : 0.;
  }
  auto meanRms = [&](int i, int nUse, double &m, double &s){
    double s1 = 0., s2 = 0.; for(int t = 0; t < nUse; t++){ s1 += r[i][t]; s2 += r[i][t]*r[i][t]; }
    m = s1/nUse; s = sqrt(TMath::Max(0., s2/nUse - m*m));
  };

  // ---- 1: input statistics -----------------------------------------------------
  {
    TCanvas *c = new TCanvas("c1", "", 1500, 700);
    c->Divide(2, 1);
    gStyle->SetPalette(kViridis);
    TH2D *mN = indexMap(hNeff, "mNeff"), *mE = indexMap(hRelE, "mRelE");
    // only cells inside the reco floor carry information into the unfolding;
    // the sub-floor reco columns are drawn too (they feed the prior) but marked
    c->cd(1); gPad->SetLeftMargin(0.18); gPad->SetRightMargin(0.17); gPad->SetBottomMargin(0.17); gPad->SetTopMargin(0.12);
    gPad->SetLogz();   // log: effective entries span ~4 orders of magnitude, all positive
    mN->GetZaxis()->SetTitle("effective entries (c/#sigma)^{2}"); mN->GetZaxis()->SetTitleOffset(1.25); mN->GetZaxis()->SetTitleSize(0.045);
    mN->SetMinimum(0.5);
    gStyle->SetPaintTextFormat(".3g"); mN->SetMarkerSize(1.0);
    mN->Draw("COLZ TEXT");
    TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.042);
    la.DrawLatex(0.18, 0.94, "response cell statistics: effective entries");
    la.SetTextSize(0.032);
    la.DrawLatex(0.18, 0.895, "toy draw per cell: Poisson(#it{n}_{eff}) #times #sigma^{2}/c");
    c->cd(2); gPad->SetLeftMargin(0.18); gPad->SetRightMargin(0.17); gPad->SetBottomMargin(0.17); gPad->SetTopMargin(0.12);
    mE->Scale(100.);
    mE->GetZaxis()->SetTitle("relative MC error #sigma/c [%]"); mE->GetZaxis()->SetTitleOffset(1.25); mE->GetZaxis()->SetTitleSize(0.045);
    mE->SetMinimum(0.); mE->SetMaximum(TMath::Min(100., mE->GetMaximum()));
    gStyle->SetPaintTextFormat(".1f"); mE->SetMarkerSize(1.0);
    mE->Draw("COLZ TEXT");
    la.SetTextSize(0.042);
    la.DrawLatex(0.18, 0.94, "response cell statistics: relative MC error");
    la.SetTextSize(0.032);
    la.DrawLatex(0.18, 0.895, "= 1/#sqrt{#it{n}_{eff}}; size of each cell's toy fluctuation");
    savePdfTight(c, Form("%s/1_responseInputStats.pdf", outDir));   // landscape canvas
    delete c;
  }

  // ---- 2: toy / nominal vs pT ---------------------------------------------------
  {
    TCanvas *c = new TCanvas("c2", "", 800, 700);
    c->SetLeftMargin(0.14); c->SetBottomMargin(0.13); c->SetRightMargin(0.04); c->SetTopMargin(0.06);
    double lo = 1., hi = 1.;
    for(int i = 1; i <= NPt; i++) if(reported(i)) for(int t = 0; t < nShowToys; t++){ lo = TMath::Min(lo, r[i][t]); hi = TMath::Max(hi, r[i][t]); }
    double pad = 0.25*(hi - lo);
    TH1F *fr = c->DrawFrame(ptReportMin, lo - pad, ptEdges[NPt], hi + 4.0*pad);   // headroom for header + legend
    fr->GetXaxis()->SetTitle("gen #it{p}_{T}^{jet} [GeV]");
    fr->GetYaxis()->SetTitle("toy / nominal unfolded spectrum");
    fr->GetXaxis()->SetTitleSize(0.048); fr->GetYaxis()->SetTitleSize(0.048);
    fr->GetXaxis()->SetLabelSize(0.040); fr->GetYaxis()->SetLabelSize(0.040); fr->GetYaxis()->SetTitleOffset(1.35);
    // +/- RMS band per bin
    int colBand = TColor::GetColor("#D55E00");
    for(int i = 1; i <= NPt; i++){
      if(!reported(i)) continue;
      double m, s; meanRms(i, nToys, m, s);
      TBox *b = new TBox(ptEdges[i-1], 1. - s, ptEdges[i], 1. + s);
      b->SetFillColorAlpha(colBand, 0.30); b->SetLineWidth(0); b->Draw();
    }
    // individual toys as thin steps
    int colToy = TColor::GetColor("#999999");
    for(int t = 0; t < nShowToys; t++){
      TH1D *h = new TH1D(Form("toy%d", t), "", NPt, ptEdges); h->SetDirectory(nullptr);
      for(int i = 1; i <= NPt; i++) h->SetBinContent(i, r[i][t]);
      h->GetXaxis()->SetRangeUser(ptReportMin, ptEdges[NPt]);
      h->SetLineColorAlpha(colToy, 0.35); h->SetLineWidth(1); h->SetStats(0);
      h->Draw("HIST same");
    }
    TH1D *hm = new TH1D("toyMean", "", NPt, ptEdges); hm->SetDirectory(nullptr);
    for(int i = 1; i <= NPt; i++){ double m, s; meanRms(i, nToys, m, s); hm->SetBinContent(i, m); }
    hm->GetXaxis()->SetRangeUser(ptReportMin, ptEdges[NPt]);
    styleLine(hm, "#000000", 3); hm->SetLineStyle(2);
    hm->Draw("HIST same");
    TLine *one = new TLine(ptReportMin, 1., ptEdges[NPt], 1.); one->SetLineColor(kGray+2); one->Draw();
    TLegend *leg = makeLegend(0.17, 0.64, 0.70, 0.80, 0.032);
    TH1D *lt = new TH1D("lt", "", 1, 0, 1); lt->SetLineColor(colToy); lt->SetLineWidth(2);
    leg->AddEntry(lt, Form("%d of %d toys", nShowToys, nToys), "l");
    TBox *lb = new TBox(); lb->SetFillColorAlpha(colBand, 0.30); lb->SetLineWidth(0);
    leg->AddEntry(lb, "#pm RMS over all toys (the systematic)", "f");
    leg->AddEntry(hm, "toy mean (bias check)", "l");
    leg->Draw();
    TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.036);
    la.DrawLatex(0.17, 0.88, "pp calo b jets: response-matrix MC stat.");
    la.SetTextSize(0.030);
    la.DrawLatex(0.17, 0.835, "each toy: fluctuated response, data unfolded at #it{N} = 2");
    savePdfTight(c, Form("%s/2_toyRatios.pdf", outDir));
    delete c;
  }

  // ---- 3: per-bin distributions -------------------------------------------------
  {
    std::vector<int> bins; for(int i = 1; i <= NPt; i++) if(reported(i)) bins.push_back(i);
    const int nc = 3, nr = (bins.size() + nc - 1)/nc;
    TCanvas *c = new TCanvas("c3", "", 1350, 450*nr);
    c->Divide(nc, nr, 0.002, 0.002);
    for(size_t k = 0; k < bins.size(); k++){
      int i = bins[k];
      c->cd(k+1); gPad->SetLeftMargin(0.15); gPad->SetBottomMargin(0.16); gPad->SetTopMargin(0.10); gPad->SetRightMargin(0.04);
      double m, s; meanRms(i, nToys, m, s);
      double w = TMath::Max(5.*s, 0.01);
      TH1D *h = new TH1D(Form("dist%d", i), "", 50, 1. - w, 1. + w); h->SetDirectory(nullptr);
      for(int t = 0; t < nToys; t++) h->Fill(r[i][t]);
      styleH(h, "#000000", markFilledCircle, 0.8);
      h->GetXaxis()->SetTitle("toy / nominal"); h->GetYaxis()->SetTitle("toys");
      h->GetXaxis()->SetTitleSize(0.055); h->GetYaxis()->SetTitleSize(0.055);
      h->GetXaxis()->SetLabelSize(0.045); h->GetYaxis()->SetLabelSize(0.045);
      h->GetXaxis()->SetNdivisions(505);
      h->SetMaximum(1.45*h->GetMaximum()); h->SetMinimum(0.);
      h->Draw("E1");
      TF1 *g = new TF1(Form("g%d", i), "gaus", 1. - w, 1. + w);
      g->SetParameters(nToys*h->GetBinWidth(1)/(sqrt(2*TMath::Pi())*s), m, s);
      g->SetLineColor(TColor::GetColor("#D55E00")); g->SetLineWidth(2); g->Draw("same");
      TLine *l1 = new TLine(1., 0., 1., h->GetMaximum()/1.45); l1->SetLineStyle(2); l1->SetLineColor(kGray+2); l1->Draw();
      TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.055);
      la.DrawLatex(0.19, 0.83, Form("%.0f < #it{p}_{T} < %.0f GeV", ptEdges[i-1], ptEdges[i]));
      la.SetTextSize(0.048);
      la.DrawLatex(0.19, 0.75, Form("RMS = %.2f%%", 100*s));
      la.DrawLatex(0.19, 0.69, Form("mean #minus 1 = %+.2f%%", 100*(m - 1.)));
    }
    c->cd(0);
    savePdfTight(c, Form("%s/3_toyDistributions.pdf", outDir));   // landscape canvas
    delete c;
  }

  // ---- 4: convergence -----------------------------------------------------------
  {
    TCanvas *c = new TCanvas("c4", "", 800, 700);
    c->SetLeftMargin(0.14); c->SetBottomMargin(0.13); c->SetRightMargin(0.04); c->SetTopMargin(0.06);
    std::vector<int> steps; for(int n = 10; n <= nToys; n += 10) steps.push_back(n);
    double ymax = 0.;
    std::vector<TGraph*> gs;
    for(int i = 1; i <= NPt; i++){
      if(!reported(i)) continue;
      TGraph *g = new TGraph();
      for(int n : steps){ double m, s; meanRms(i, n, m, s); g->SetPoint(g->GetN(), n, 100*s); ymax = TMath::Max(ymax, 100*s); }
      gs.push_back(g);
    }
    TGaxis::SetMaxDigits(4);   // 1000 toys, not "1 x10^3"
    TH1F *fr = c->DrawFrame(0., 0., nToys, 1.7*ymax);   // headroom for header + legend
    fr->GetXaxis()->SetTitle("number of toys");
    fr->GetYaxis()->SetTitle("running RMS of toy / nominal [%]");
    fr->GetXaxis()->SetTitleSize(0.048); fr->GetYaxis()->SetTitleSize(0.048);
    fr->GetXaxis()->SetLabelSize(0.040); fr->GetYaxis()->SetLabelSize(0.040); fr->GetYaxis()->SetTitleOffset(1.25);
    TLegend *leg = makeLegend(0.40, 0.62, 0.95, 0.84, 0.032); leg->SetNColumns(2);
    int k = 0;
    for(int i = 1; i <= NPt; i++){
      if(!reported(i)) continue;
      TGraph *g = gs[k];
      const char *hexes[6] = {"#000000", "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00"};   // Okabe-Ito, no yellow
      int col = TColor::GetColor(hexes[k % 6]);
      g->SetLineColor(col); g->SetLineWidth(3); g->SetLineStyle(1 + (k % 3));
      g->Draw("L same");
      leg->AddEntry(g, Form("%.0f-%.0f GeV", ptEdges[i-1], ptEdges[i]), "l");
      k++;
    }
    leg->Draw();
    TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.036);
    la.DrawLatex(0.17, 0.88, "response toys: convergence of the RMS");
    savePdfTight(c, Form("%s/4_toyConvergence.pdf", outDir));
    TGaxis::SetMaxDigits(3);
    delete c;
  }

  // ---- printout -------------------------------------------------------------------
  printf("\n%-9s %10s %10s %12s\n", "jet pT", "RMS", "mean-1", "stored syst");
  for(int i = 1; i <= NPt; i++){
    double m, s; meanRms(i, nToys, m, s);
    printf("%3.0f-%-5.0f %9.2f%% %+9.2f%% %11.2f%%%s\n", ptEdges[i-1], ptEdges[i], 100*s, 100*(m - 1.), 100*hSys->GetBinContent(i),
           reported(i) ? "" : "  (underflow)");
  }
  printf("figures in %s\n", outDir);
}
