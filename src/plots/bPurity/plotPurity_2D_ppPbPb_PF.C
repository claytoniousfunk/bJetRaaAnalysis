// NOMINAL b purity of all muon-tagged jets from the 2D muon ptRel x CSVv2 fit (HistFactory, MC tagger shapes),
// pp and the four PbPb classes in a 1x5 row, 80-200 GeV, with its total systematic as boxes, the MC b fraction,
// and the old nominal CSVv2-loose cut purity (P_tag f_pass / eps_b(MC)) as a cross-check (2026-11-16, user: the 2D
// fit is the nominal in the AN; AN figure fig:bPurity). Draws stored results only, refits nothing:
//   2D:         rootFiles/systematics/bPuritySys_total_PFJets_pp_csvV2_2D.root, bPuritySys_PFJets_PbPb_csvV2_2D.root
//               (systematics_CSVv2_2D_PF.C), h_bPurity_nominal[_C<i>] (bars = fit stat.) and h_sysRel_total
//   CSVv2 cut:  the same names with _csvV2loose (systematics_CSVv2Cut_PF.C(0.54))
//   MC b fraction: h_bFractionMC_<sys> of rootFiles/bPurityResults/ptRelCSVv2_2D_PF_ppPbPb_4cat_fit5.root
//               (fitPtRelCSVv2_2D_ppPbPb_PF.C)
// Figure: figures/bPurity/PFJets_PbPb/csvV2/purity_2D_ppPbPb.pdf
// Usage, from src/plots/bPurity/:  root -l -b -q 'plotPurity_2D_ppPbPb_PF.C'
#include "TFile.h"
#include "TH1D.h"
#include "TBox.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "../../../headers/plotting/plotStyle.h"
#include "../../../headers/plotting/gridPads.h"

void plotPurity_2D_ppPbPb_PF(double ptLo = 80., double ptHi = 200., bool fine = false){   // fine (2026-11-16): the 7-bin nominal, CSVv2 cut stays on its own 4 windows
  initPlotStyle();
  const TString repo = "/home/clayton/Analysis/code/bJetRaaAnalysis";
  const int NCOL = 5; const char *colName[NCOL] = {"pp", "50-80%", "30-50%", "10-30%", "0-10%"}, *colTag[NCOL] = {"pp", "C4", "C3", "C2", "C1"};
  const int NW4 = 4; const double edges4[NW4 + 1] = {80., 100., 120., 150., 200.};
  const int NW7 = 7; const double edges7[NW7 + 1] = {80., 90., 100., 110., 120., 135., 150., 200.};
  const int NW = fine ? NW7 : NW4; const double *edges = fine ? edges7 : edges4;
  TFile *f2pp = TFile::Open(repo + "/rootFiles/systematics/bPuritySys_total_PFJets_pp_csvV2_2D" + (fine ? "7bins" : "") + ".root");
  TFile *f2pb = TFile::Open(repo + "/rootFiles/systematics/bPuritySys_PFJets_PbPb_csvV2_2D" + (fine ? "7bins" : "") + ".root");
  TFile *fCpp = TFile::Open(repo + "/rootFiles/systematics/bPuritySys_total_PFJets_pp_csvV2loose.root");
  TFile *fCpb = TFile::Open(repo + "/rootFiles/systematics/bPuritySys_PFJets_PbPb_csvV2loose.root");
  TFile *fMC  = TFile::Open(repo + "/rootFiles/bPurityResults/ptRelCSVv2_2D_PF_ppPbPb_4cat_fit5" + (fine ? "_7bins" : "") + ".root");
  if(!f2pp || !f2pb || !fCpp || !fCpb || !fMC){ printf("ERROR: inputs\n"); return; }
  auto proj = [&](TFile *f, TString name, TString nm, int NW, const double *edges) -> TH1D* {
    TH1 *h = nullptr; f->GetObject(name, h);
    if(!h){ printf("ERROR: %s missing in %s\n", name.Data(), f->GetName()); return nullptr; }
    TH1D *o = new TH1D(nm, "", NW, edges); o->SetDirectory(nullptr);
    for(int w = 0; w < NW; w++){
      const int b = h->FindBin(0.5*(edges[w] + edges[w+1]));
      if(fabs(h->GetBinLowEdge(b) - edges[w]) > 1e-6 || fabs(h->GetBinLowEdge(b+1) - edges[w+1]) > 1e-6){ printf("ERROR: binning of %s\n", name.Data()); return nullptr; }
      o->SetBinContent(w + 1, h->GetBinContent(b)); o->SetBinError(w + 1, h->GetBinError(b));
    }
    return o;
  };
  TH1D *p2[NCOL], *s2[NCOL], *pc[NCOL], *mc[NCOL];
  for(int c = 0; c < NCOL; c++){
    const TString sfx = c == 0 ? "" : TString("_") + colTag[c];
    p2[c] = proj(c ? f2pb : f2pp, "h_bPurity_nominal" + sfx, Form("p2_%d", c), NW, edges);
    s2[c] = proj(c ? f2pb : f2pp, "h_sysRel_total" + sfx, Form("s2_%d", c), NW, edges);
    pc[c] = proj(c ? fCpb : fCpp, "h_bPurity_nominal" + sfx, Form("pc_%d", c), NW4, edges4);
    mc[c] = proj(fMC, Form("h_bFractionMC_%s", colTag[c]), Form("mc_%d", c), NW, edges);
    if(!p2[c] || !s2[c] || !pc[c] || !mc[c]) return;
  }
  printf("\nb purity of all muon-tagged jets: 2D fit +- stat (syst %%) | CSVv2 > 0.54 cut | MC b fraction\n");
  for(int c = 0; c < NCOL; c++) for(int w = 1; w <= NW; w++)
    printf("  %-7s %3.0f-%3.0f  %.3f +- %.3f (%4.1f%%) | %.3f +- %.3f | %.3f\n", w == 1 ? colName[c] : "", edges[w-1], edges[w],
           p2[c]->GetBinContent(w), p2[c]->GetBinError(w), 100*s2[c]->GetBinContent(w), pc[c]->GetBinContent(pc[c]->FindBin(0.5*(edges[w-1] + edges[w]))), pc[c]->GetBinError(pc[c]->FindBin(0.5*(edges[w-1] + edges[w]))), mc[c]->GetBinContent(w));

  TCanvas *cv = new TCanvas("cv2D", "", 2300, 950);
  GridPads g = makeGridPads(cv, NCOL, 1, 0.055, 0.01, 0.16, 0.10);
  TLegend *lg = nullptr;
  for(int c = 0; c < NCOL; c++){
    g.at(c, 0)->cd();
    TH1F *fr = gPad->DrawFrame(ptLo, 0., ptHi, 0.95); g.styleFrame(fr, c, 0, 0.055);
    fr->GetYaxis()->SetLabelFont(43); fr->GetYaxis()->SetLabelSize(c == 0 ? 34 : 0);
    fr->GetXaxis()->SetLabelFont(43); fr->GetXaxis()->SetLabelSize(34);
    for(int w = 1; w <= NW; w++){
      const double y = p2[c]->GetBinContent(w), d = y*s2[c]->GetBinContent(w);
      TBox *b = new TBox(edges[w-1], y - d, edges[w], y + d); b->SetFillColorAlpha(kGray + 1, 0.45); b->SetLineWidth(0); b->Draw();
    }
    styleLine(mc[c], hexMC, 4); mc[c]->SetLineStyle(2); mc[c]->Draw("HIST SAME");
    styleH(pc[c], okabeHex[6], markOpenCircle, 2.4); pc[c]->SetLineWidth(3); pc[c]->Draw("E1 X0 SAME");
    styleH(p2[c], hexData, markFilledSquare, 2.2); p2[c]->SetLineWidth(3); p2[c]->Draw("E1 X0 SAME");
    gPad->RedrawAxis();
    g.text(c, 0, 0.10, colName[c], 0.065, 62);
    if(c == 0){
      lg = makeLegend(0.14, 0.14, 0.98, 0.36, 0.06);
      TH1D *box = (TH1D*) p2[c]->Clone("legBox"); box->SetFillColorAlpha(kGray + 1, 0.45); box->SetLineWidth(3);
      lg->AddEntry(box, "2D fit (syst.)", "pef"); lg->AddEntry(pc[c], fine ? "CSVv2 > 0.54 cut (4 bins)" : "CSVv2 > 0.54 cut", "pe"); lg->AddEntry(mc[c], "MC b fraction", "l");
      lg->Draw();
    }
  }
  g.titles(cv, "#it{p}_{T}^{jet} (GeV)", "b purity of muon-tagged jets", 0.05, 0.02, 0.02);
  g.cmsLabel(cv, "Internal", "pp, PbPb (5.02 TeV)   2D #it{p}_{T}^{rel} #times CSVv2 fit", 1.0);
  savePdfTight(cv, repo + "/figures/bPurity/PFJets_PbPb/csvV2/purity_2D_ppPbPb" + (fine ? "_7bins" : "") + ".pdf");
}
