// Forest akCs4PF emulation vs the forest's own akCs4PF jets, per coarse class,
// from ONE same-event PbPb_pfCandAnalyzer scan (doForestCSJets, pseudoJets.h).
//
//   emulated : h_csJetPt_JEC_C<i>       (csJetEmulation.h on the event's PF
//                                        candidates, JEC_PF, |eta| < 1.6,
//                                        jetTrkMax filter on the hardest charged
//                                        PF candidate within 0.4)
//   forest   : h_inclRecoJetPt_C<i>     (akCs4PFJetAnalyzer rawpt, JEC_PF,
//                                        |eta| < 1.6, pT > 20, jetTrkMax filter
//                                        on the forest's track-based trackMax)
//
// Both per event (h_vz) and per GeV. The forest stores only jets with jtpt >= 20
// (its own JEC), so the comparison starts at 20 GeV. Same events, and mostly the
// same jets, so the ratio is drawn WITHOUT error bars: either spectrum's
// statistical error would overstate the uncertainty of a ratio of two nearly
// identical, fully correlated jet samples by an order of magnitude.
//
// Log y on the spectra (they fall over four decades), linear ratio.
// Run from: src/plots/csJets/

#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLatex.h"
#include "TLine.h"
#include "TSystem.h"
#include "../../../headers/plotting/plotStyle.h"
#include "../../../headers/plotting/coarseCent.h"
#include "../../../headers/plotting/ratioPanel.h"

const char *inFile = "/home/clayton/Downloads/PbPb_pfCandAnalyzer_output_1.root";
const char *outPdf = "../../../figures/csJets/csJetSpectra_vsForest_sameEvent.pdf";
const double edges[] = {20, 25, 30, 35, 40, 50, 60, 80, 100, 150};
const int    nEdges  = sizeof(edges) / sizeof(double) - 1;

void plotCSJetSpectra_vsForest(const char *path = inFile, const char *pdf = outPdf)
{
  initPlotStyle();
  gSystem->mkdir("../../../figures/csJets/", kTRUE);
  TFile *f = TFile::Open(path);
  if(!f || f->IsZombie()){ printf("ERROR: %s did not open\n", path); return; }

  TCanvas *c = new TCanvas("c", "", 1600, 650);
  const int order[NCoarse] = {3, 2, 1, 0};          // peripheral to central
  for(int p = 0; p < NCoarse; p++){
    const int ci = order[p];
    const double nEvt = coarseEvents(f, ci);
    TH1D *e0 = coarseSum1(f, "h_csJetPt_JEC", ci, "emu");
    TH1D *f0 = coarseSum1(f, "h_inclRecoJetPt", ci, "for");
    if(nEvt <= 0 || !e0 || !f0){ printf("ERROR: class %s incomplete in %s\n", coarseLabel[ci], path); return; }
    TH1D *e = rebinTo(e0, nEdges, edges, Form("emu_%d", ci));
    TH1D *fr = rebinTo(f0, nEdges, edges, Form("for_%d", ci));
    e->Scale(1. / nEvt, "width");
    fr->Scale(1. / nEvt, "width");
    TH1D *r = makeRatio(e, fr, Form("ratio_%d", ci));

    c->cd();
    TPad *pp = new TPad(Form("pp%d", p), "", p / 4., 0, (p + 1) / 4., 1);
    pp->Draw(); pp->cd();
    TPad *top, *bot;
    splitPads(top, bot, 0.34, 0.20);
    top->SetName(Form("top%d", p)); bot->SetName(Form("bot%d", p));
    top->SetRightMargin(0.04); bot->SetRightMargin(0.04);

    top->cd(); top->SetLogy();
    styleH(fr, hexData, markFilledCircle, 1.3);
    styleH(e, hexCorrected, markOpenSquare, 1.5);
    fr->SetTitle("");
    fr->GetYaxis()->SetTitle("1/N_{evt} dN_{jet}/d#it{p}_{T} [GeV^{-1}]");
    fr->GetYaxis()->SetTitleSize(0.060); fr->GetYaxis()->SetTitleOffset(1.45);
    fr->GetYaxis()->SetLabelSize(0.052); fr->GetXaxis()->SetLabelSize(0);
    fr->SetMaximum(5. * TMath::Max(fr->GetMaximum(), e->GetMaximum()));
    double mn = 1e9;
    for(int b = 1; b <= fr->GetNbinsX(); b++) if(fr->GetBinContent(b) > 0) mn = TMath::Min(mn, fr->GetBinContent(b));
    fr->SetMinimum(0.3 * mn);
    fr->Draw("E1");
    e->Draw("E1 SAME");
    TLatex l; l.SetNDC(); l.SetTextSize(0.062);
    l.DrawLatex(p == 0 ? 0.42 : 0.70, 0.84, coarseLabel[ci]);
    if(p == 0){
      TLegend *lg = makeLegend(0.42, 0.50, 0.95, 0.66, 0.055);
      lg->AddEntry(fr, "forest akCs4PF", "pl");
      lg->AddEntry(e, "emulated akCs4PF", "pl");
      lg->Draw();
      TLatex h; h.SetNDC(); h.SetTextSize(0.044);
      h.DrawLatex(0.34, 0.78, "PbPb MinBias, same event");
      h.DrawLatex(0.34, 0.71, "|#it{#eta}| < 1.6, JEC");
    }

    bot->cd();
    styleH(r, hexCorrected, markOpenSquare, 1.5);
    styleRatioAxes(r, "jet #it{p}_{T} [GeV]", "emulated / forest");
    for(int b = 1; b <= r->GetNbinsX(); b++) r->SetBinError(b, 0.);   // see header
    r->SetMinimum(0.90); r->SetMaximum(1.10);
    r->Draw("P");
    unityLine(edges[0], edges[nEdges])->Draw();
    r->Draw("P SAME");
    if(p == 0){
      TLatex nb; nb.SetNDC(); nb.SetTextSize(0.075);
      nb.DrawLatex(0.25, 0.86, "same jets: no error bars");
    }

    printf("  %-7s N_evt %6.0f   emulated/forest per bin:", coarseLabel[ci], nEvt);
    for(int b = 1; b <= r->GetNbinsX(); b++) printf(" %.3f", r->GetBinContent(b));
    printf("\n");
  }
  savePdfTight(c, pdf);
}
