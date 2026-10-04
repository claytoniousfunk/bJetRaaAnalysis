// Total systematic uncertainty on the pp CALO-jet b purity.
//
// Adds the relative systematics in quadrature, as systematicsPlots_4CentBins.C
// does. Sources, each from its own macro in src/systematics/:
//   bGS           gluon-splitting share +0.150..+0.200 (nominal +0.175)   systematics_bGS_caloJets_pp.C
//   c fraction    c-multiplier 0.8..1.2                                   systematics_cMult_caloJets_pp.C
//   fit low edge  ptRel lower edge 0..0.3 GeV                             systematics_fitRange_caloJets_pp.C
//   fit high edge ptRel upper edge 4.0..5.0 GeV                           systematics_fitRangeHigh_caloJets_pp.C
//   JEU           MC jet pT up/down by the JEU uncertainty                systematics_JEUJER_caloJets_pp.C
//   JER           MC jet pT smeared                                       systematics_JEUJER_caloJets_pp.C
// Run those first. The lower-edge file also carries old h_sensRel_fitHigh*
// histograms; the upper edge is taken from its own file only, so it is not
// counted twice.
//
// Every input must have the same nominal purity in every bin; the macro checks
// that and stops if not (it would mean the inputs come from different
// settings or scans).
//
// CAVEATS carried in from the inputs: JEU and JER use the AK4PF uncertainty /
// resolution for calo jets, and several sources are partly MC statistical
// fluctuation of the templates (see the individual macros).
//
// Output: rootFiles/systematics/bPuritySys_total_caloJets_pp.root
//   h_bPurity_nominal       nominal purity with fit (stat) error
//   h_sysRel_total          total relative systematic
//   h_sysAbs_total          total absolute systematic
//   h_sysRel_<source>       each input, copied for reference
//   h_bPurity_sysAbs        nominal purity with the total systematic as its error
// and figures/systematics/total_caloJets_pp_breakdown.pdf,
//     figures/systematics/total_caloJets_pp_bPurity.pdf
//
// Usage: root -l -b -q systematics_total_caloJets_pp.C
// Run from: src/systematics/

#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLine.h"
#include "TBox.h"
#include "TSystem.h"
#include "TColor.h"
#include "TGaxis.h"
#include "TStyle.h"
#include "../../headers/plotting/plotStyle.h"
#include "../../headers/functions/caloJetBPurityFit_pp.h"   // ptEdges, ptReportMin, underflow band

const char *sysDir = "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/systematics";
const char *figDir = "/home/clayton/Analysis/code/bJetRaaAnalysis/figures/systematics";

struct Source { const char *file, *hist, *key, *label, *hex; int style; };
const int NSrc = 6;
const Source src[NSrc] = {
  {"bPuritySys_bGS_caloJets_pp.root",          "h_sysRel_bGS",          "bGS",          "g#rightarrowb#bar{b} share", "#E69F00", 1},
  {"bPuritySys_cMult_caloJets_pp.root",        "h_sysRel_cMult",        "cMult",        "c fraction",                  "#56B4E9", 2},
  {"bPuritySys_fitRange_caloJets_pp.root",     "h_sysRel_fitRange",     "fitRangeLow",  "fit range, lower edge",       "#009E73", 1},
  {"bPuritySys_fitRangeHigh_caloJets_pp.root", "h_sysRel_fitRangeHigh", "fitRangeHigh", "fit range, upper edge",       "#0072B2", 2},
  {"bPuritySys_JEU_caloJets_pp.root",          "h_sysRel_JEU",          "JEU",          "jet energy scale",            "#D55E00", 1},
  {"bPuritySys_JER_caloJets_pp.root",          "h_sysRel_JER",          "JER",          "jet energy resolution",       "#CC79A7", 2},
};

void systematics_total_caloJets_pp()
{
  initPlotStyle();
  gSystem->mkdir(figDir, kTRUE);

  TH1D *rel[NSrc], *nomRef = nullptr;
  for(int s = 0; s < NSrc; s++){
    TFile *f = TFile::Open(Form("%s/%s", sysDir, src[s].file));
    if(!f || f->IsZombie()){ printf("ERROR: cannot open %s -- run its macro first\n", src[s].file); return; }
    TH1D *h = nullptr, *n = nullptr;
    f->GetObject(src[s].hist, h); f->GetObject("h_bPurity_nominal", n);
    if(!h || !n){ printf("ERROR: %s lacks %s or h_bPurity_nominal\n", src[s].file, src[s].hist); return; }
    rel[s] = (TH1D*) h->Clone(Form("h_sysRel_%s", src[s].key)); rel[s]->SetDirectory(nullptr);
    if(!nomRef){ nomRef = (TH1D*) n->Clone("h_bPurity_nominal"); nomRef->SetDirectory(nullptr); }
    else {
      if(n->GetNbinsX() != nomRef->GetNbinsX()){ printf("ERROR: %s has different binning\n", src[s].file); return; }
      for(int i = 1; i <= n->GetNbinsX(); i++)
        if(fabs(n->GetBinContent(i) - nomRef->GetBinContent(i)) > 1e-6){
          printf("ERROR: nominal purity differs in %s, bin %d (%.5f vs %.5f) -- inputs from different settings?\n",
                 src[s].file, i, n->GetBinContent(i), nomRef->GetBinContent(i)); return; }
    }
    f->Close();
  }
  const int NB = nomRef->GetNbinsX();

  TH1D *tot  = (TH1D*) rel[0]->Clone("h_sysRel_total"); tot->Reset(); tot->SetTitle("total systematic (quadrature sum);#it{p}_{T}^{jet} [GeV];relative uncertainty");
  TH1D *totA = (TH1D*) tot->Clone("h_sysAbs_total"); totA->SetTitle("total systematic (absolute);#it{p}_{T}^{jet} [GeV];absolute uncertainty");
  TH1D *pSys = (TH1D*) nomRef->Clone("h_bPurity_sysAbs"); pSys->SetTitle("b purity with total systematic as error;#it{p}_{T}^{jet} [GeV];b purity");

  printf("\n%-9s", "jet pT"); for(int s = 0; s < NSrc; s++) printf(" %12s", src[s].key); printf("  %8s %8s %8s\n", "TOTAL", "stat", "purity");
  for(int i = 1; i <= NB; i++){
    double q = 0.; for(int s = 0; s < NSrc; s++) q += pow(rel[s]->GetBinContent(i), 2);
    double t = sqrt(q), p = nomRef->GetBinContent(i);
    tot->SetBinContent(i, t); tot->SetBinError(i, 0.);
    totA->SetBinContent(i, t*p); totA->SetBinError(i, 0.);
    pSys->SetBinError(i, t*p);
    printf("%3.0f-%-5.0f", tot->GetXaxis()->GetBinLowEdge(i), tot->GetXaxis()->GetBinUpEdge(i));
    if(p <= 0.){ printf("  not fitted (too little data)\n"); continue; }
    for(int s = 0; s < NSrc; s++) printf(" %11.2f%%", 100*rel[s]->GetBinContent(i));
    printf("  %7.2f%% %7.2f%%  %.3f\n", 100*t, 100*nomRef->GetBinError(i)/p, p);
  }

  TFile *fo = TFile::Open(Form("%s/bPuritySys_total_caloJets_pp.root", sysDir), "recreate");
  nomRef->Write(); tot->Write(); totA->Write(); pSys->Write();
  for(int s = 0; s < NSrc; s++) rel[s]->Write();
  TString srcList; for(int s = 0; s < NSrc; s++) srcList += Form("%s%s", s ? ", " : "", src[s].key);
  TNamed info("info", Form("pp calo jets b purity. Total = quadrature sum of relative systematics: %s. "
                            "JEU/JER use AK4PF uncertainty/resolution for calo jets.", srcList.Data()));
  info.Write();
  fo->Close();
  printf("\nwritten %s/bPuritySys_total_caloJets_pp.root\n", sysDir);

  const double xLo = tot->GetXaxis()->GetXmin(), xHi = tot->GetXaxis()->GetXmax();

  // ---- breakdown: each source and the total, relative, in % ---------------
  {
    TCanvas *c = new TCanvas("cBreak", "", 800, 700);
    c->SetLeftMargin(0.13); c->SetBottomMargin(0.13); c->SetRightMargin(0.04); c->SetTopMargin(0.06);
    double ymax = 100*tot->GetMaximum();
    TH1F *fr = c->DrawFrame(xLo, 0., xHi, 2.1*ymax);   // headroom for header + legend
    fr->GetXaxis()->SetTitle("#it{p}_{T}^{jet} [GeV]");
    fr->GetYaxis()->SetTitle("relative systematic uncertainty [%]");
    fr->GetXaxis()->SetTitleSize(0.048); fr->GetYaxis()->SetTitleSize(0.048);
    fr->GetXaxis()->SetLabelSize(0.040); fr->GetYaxis()->SetLabelSize(0.040);
    drawUnderflowBand(0., 2.1*ymax, 0.08*ymax);
    TLegend *leg = makeLegend(0.17, 0.62, 0.95, 0.84, 0.032);
    leg->SetNColumns(2);
    for(int s = 0; s < NSrc; s++){
      TH1D *h = (TH1D*) rel[s]->Clone(Form("b%d", s)); h->Scale(100.);
      for(int i = 0; i <= NB+1; i++) h->SetBinError(i, 0.);
      styleLine(h, src[s].hex, 3); h->SetLineStyle(src[s].style);
      h->Draw("HIST same"); leg->AddEntry(h, src[s].label, "l");
    }
    TH1D *ht = (TH1D*) tot->Clone("bt"); ht->Scale(100.); styleLine(ht, "#000000", 4);
    ht->Draw("HIST same"); leg->AddEntry(ht, "total (quadrature)", "l");
    leg->Draw();
    TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.038);
    la.DrawLatex(0.17, 0.875, "pp 5.02 TeV, calo jets: b-jet purity, 2-template #it{p}_{T}^{rel} fit");
    savePdfTight(c, Form("%s/total_caloJets_pp_breakdown.pdf", figDir));
    delete c;
  }

  // ---- purity with stat (bars) and total syst (boxes) ----------------------
  {
    TCanvas *c = new TCanvas("cPur", "", 800, 700);
    c->SetLeftMargin(0.13); c->SetBottomMargin(0.13); c->SetRightMargin(0.04); c->SetTopMargin(0.06);
    TH1F *fr = c->DrawFrame(xLo, 0., xHi, 1.0);
    fr->GetXaxis()->SetTitle("#it{p}_{T}^{jet} [GeV]");
    fr->GetYaxis()->SetTitle("b-jet purity");
    fr->GetXaxis()->SetTitleSize(0.048); fr->GetYaxis()->SetTitleSize(0.048);
    fr->GetXaxis()->SetLabelSize(0.040); fr->GetYaxis()->SetLabelSize(0.040);
    drawUnderflowBand(0., 1.0, 0.05);
    int colSys = TColor::GetColor("#D55E00");
    for(int i = 1; i <= NB; i++){
      if(nomRef->GetBinContent(i) <= 0.) continue;   // not fitted
      double x1 = tot->GetXaxis()->GetBinLowEdge(i), x2 = tot->GetXaxis()->GetBinUpEdge(i);
      double p = nomRef->GetBinContent(i), e = pSys->GetBinError(i);
      TBox *b = new TBox(x1, p - e, x2, p + e);
      b->SetFillColorAlpha(colSys, 0.30); b->SetLineColor(colSys); b->SetLineWidth(1);
      b->Draw("l"); b->Draw();
    }
    TH1D *hp = (TH1D*) nomRef->Clone("hp"); styleH(hp, "#000000", markFilledCircle, 1.3);
    blankUnfitted(hp, nomRef);
    hp->Draw("E1 X0 same");
    TLegend *leg = makeLegend(0.50, 0.78, 0.95, 0.92, 0.036);
    leg->AddEntry(hp, "fit (stat. uncertainty)", "lp");
    TBox *lb = new TBox(); lb->SetFillColorAlpha(colSys, 0.30); lb->SetLineColor(colSys);
    leg->AddEntry(lb, "total systematic", "f");
    leg->Draw();
    TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.040);
    la.DrawLatex(0.17, 0.87, "pp 5.02 TeV, calo jets");
    la.SetTextSize(0.032);
    la.DrawLatex(0.17, 0.82, "2-template #it{p}_{T}^{rel} fit");
    savePdfTight(c, Form("%s/total_caloJets_pp_bPurity.pdf", figDir));
    delete c;
  }
  printf("written %s/total_caloJets_pp_breakdown.pdf and total_caloJets_pp_bPurity.pdf\n", figDir);
}
