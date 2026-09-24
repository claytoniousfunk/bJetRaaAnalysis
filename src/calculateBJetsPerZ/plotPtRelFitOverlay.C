// The b-purity template fit drawn on top of the data it was fitted to, for
// PbPb 0-10% and pp in the same jet pT window.
//
// WHY. The fit returns a central-PbPb purity of 0.79 +- 0.025 at 80-90 GeV
// against 0.555 +- 0.003 in pp -- 42% higher, where more combinatorial
// background in central events should push it LOWER. The fit is well converged
// (chi2/ndof ~ 1) and an MC closure with the bGS multiplier correctly at 1.0
// shows no significant bias (-0.9 sigma in 0-10%). So the fitter is not at
// fault and the question is whether the TEMPLATES describe the data's ptRel
// shape. That closure could not test this: it fits MC with MC templates.
//
// READING IT. A data/fit ratio that is flat means the two templates span the
// data and the purity is whatever the fit says. Structure in the ratio means
// they do not, and the fitted fraction is absorbing a shape difference.
//
// AND A CAVEAT THAT CUTS BOTH WAYS. The MC has no quenching. Real central-PbPb
// jets lose energy and their fragmentation is modified, so the ptRel of a
// genuine b jet in data need not match the PYTHIA+HYDJET b template even if
// every reconstruction effect were perfect. A template/data mismatch here is
// therefore NOT automatically a bug -- it may be physics the MC does not model.
// What the comparison can say is whether the fitted mixture reproduces the
// data at all, and whether pp (where quenching is absent) behaves better.
//
// The 2-template fit is f = b*h_b + (1-b)*h_l, with b and light-merged
// templates, over ptRel 0-5.
//
// Usage, from src/calculateBJetsPerZ/:
//   root -l -b -q -e 'gSystem->Load("/home/clayton/Programs/RooUnfold/build/libRooUnfold.so");' \
//        -e '.L calculateBJetsPerZ.cc' -e 'openTemplateFitterFiles();' \
//        -e '.x plotPtRelFitOverlay.C'

#include "../../headers/plotting/plotStyle.h"

void plotPtRelFitOverlay()
{
  initPlotStyle();

  const int    do_data = 1, do_mergeB = 1, do_mergeC = 1, do_2t = 1, do_3t = 0;
  const double xlo = 0.0, xhi = 5.0, cm = 1.0, bgs = 1.175;  // data fit: multiplier ON
  const double ptLo = 80., ptHi = 90.;

  struct Sys { const char *name; int pp, c1, c2, c3, c4; };
  Sys sys[2] = { {"PbPb 0-10%", 0,1,0,0,0}, {"pp", 1,0,0,0,0} };

  TCanvas *c = new TCanvas("cFit", "", 1100, 650);

  for(int s = 0; s < 2; s++){

    const double bFrac = templateFitter(do_data, sys[s].pp,sys[s].c1,sys[s].c2,sys[s].c3,sys[s].c4,
                                        ptLo, ptHi, do_mergeC,do_mergeB, do_2t,do_3t,
                                        xlo,xhi, cm,bgs, 1);
    const double bErr  = templateFitter(do_data, sys[s].pp,sys[s].c1,sys[s].c2,sys[s].c3,sys[s].c4,
                                        ptLo, ptHi, do_mergeC,do_mergeB, do_2t,do_3t,
                                        xlo,xhi, cm,bgs, 2);

    // templateFitter leaves the rebinned, unit-normalised data and templates here
    TH1D *dat = (TH1D*) h_draw  ->Clone(Form("dat_%d", s));
    TH1D *tb  = (TH1D*) h_b_draw->Clone(Form("tb_%d",  s));
    TH1D *tl  = (TH1D*) h_l_draw->Clone(Form("tl_%d",  s));
    dat->SetDirectory(nullptr); tb->SetDirectory(nullptr); tl->SetDirectory(nullptr);

    // the fitted mixture, exactly as func_temp_1 builds it
    TH1D *fit = (TH1D*) tb->Clone(Form("fit_%d", s));
    fit->SetDirectory(nullptr);
    fit->Scale(bFrac);
    fit->Add(tl, 1.0 - bFrac);

    TH1D *bComp = (TH1D*) tb->Clone(Form("bc_%d", s)); bComp->Scale(bFrac);
    TH1D *lComp = (TH1D*) tl->Clone(Form("lc_%d", s)); lComp->Scale(1.0 - bFrac);
    bComp->SetDirectory(nullptr); lComp->SetDirectory(nullptr);

    const double x0 = s/2., x1 = (s+1)/2.;
    const double split = 0.32;

    c->cd();
    TPad *pB = new TPad(Form("pB%d", s), "", x0, 0., x1, split);
    pB->SetLeftMargin(0.16); pB->SetRightMargin(0.03);
    pB->SetTopMargin(0.02);  pB->SetBottomMargin(0.32);
    pB->Draw();
    TPad *pT = new TPad(Form("pT%d", s), "", x0, split, x1, 1.);
    pT->SetLeftMargin(0.16); pT->SetRightMargin(0.03);
    pT->SetTopMargin(0.08);  pT->SetBottomMargin(0.02);
    pT->SetLogy();
    pT->Draw(); pT->cd();

    TH1D *fr = new TH1D(Form("frT%d", s), "", 1, xlo, xhi);
    fr->SetStats(0);
    fr->GetYaxis()->SetRangeUser(1e-4, 5.);
    fr->GetYaxis()->SetTitle("normalized / GeV");
    fr->GetXaxis()->SetLabelSize(0.);
    fr->GetYaxis()->SetTitleSize(0.058); fr->GetYaxis()->SetLabelSize(0.052);
    fr->GetYaxis()->SetTitleOffset(1.30);
    fr->Draw("axis");

    styleH(dat,   hexData,     markFilledCircle, 1.1);
    styleLine(fit,   okabeHex[6], 3);   // vermillion: the fitted sum
    styleLine(bComp, okabeHex[5], 2);   // blue: b component
    styleLine(lComp, okabeHex[3], 2);   // green: light component
    lComp->SetLineStyle(2);

    bComp->Draw("hist same");
    lComp->Draw("hist same");
    fit  ->Draw("hist same");
    dat  ->Draw("ep same");

    TLatex t; t.SetNDC(); t.SetTextFont(42);
    t.SetTextSize(0.062); t.DrawLatex(0.22, 0.87, sys[s].name);
    t.SetTextSize(0.044);
    t.DrawLatex(0.22, 0.80, Form("%.0f < jet p_{T} < %.0f GeV", ptLo, ptHi));
    t.DrawLatex(0.22, 0.74, Form("b fraction = %.3f #pm %.3f", bFrac, bErr));

    if(s == 0){
      TLegend *leg = new TLegend(0.55, 0.60, 0.96, 0.89);
      leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.044);
      leg->AddEntry(dat,   "data", "lp");
      leg->AddEntry(fit,   "fit: b + light", "l");
      leg->AddEntry(bComp, "b component", "l");
      leg->AddEntry(lComp, "light component", "l");
      leg->Draw();
    }

    // ---- data / fit -------------------------------------------------------
    pB->cd();
    TH1D *rat = (TH1D*) dat->Clone(Form("rat_%d", s));
    rat->SetDirectory(nullptr);
    rat->Divide(fit);                       // fit treated as fixed
    rat->SetStats(0);
    rat->GetYaxis()->SetRangeUser(0.45, 1.55);
    rat->GetXaxis()->SetTitle("muon p_{T}^{rel} [GeV]");
    rat->GetYaxis()->SetTitle("data / fit");
    rat->GetXaxis()->SetTitleSize(0.125); rat->GetXaxis()->SetLabelSize(0.115);
    rat->GetYaxis()->SetTitleSize(0.110); rat->GetYaxis()->SetLabelSize(0.105);
    rat->GetXaxis()->SetTitleOffset(1.10);
    rat->GetYaxis()->SetTitleOffset(0.60);
    rat->GetYaxis()->SetNdivisions(505);
    styleH(rat, hexData, markFilledCircle, 1.0);
    rat->Draw("ep");

    TLine *one = new TLine(xlo, 1., xhi, 1.);
    one->SetLineStyle(2); one->SetLineColor(kGray+2); one->Draw();

    printf("\n  %s, %.0f-%.0f GeV: b = %.4f +- %.4f\n", sys[s].name, ptLo, ptHi, bFrac, bErr);
    printf("    ptRel      data       fit     data/fit\n");
    for(int i = 1; i <= dat->GetNbinsX(); i++){
      if(dat->GetBinCenter(i) > xhi) break;
      printf("   %4.2f-%4.2f %9.4f %9.4f %9.3f\n",
             dat->GetBinLowEdge(i), dat->GetBinLowEdge(i+1),
             dat->GetBinContent(i), fit->GetBinContent(i), rat->GetBinContent(i));
    }
  }

  c->SaveAs("/home/clayton/Analysis/code/bJetRaaAnalysis/figures/templateFits/"
            "ptRelFitOverlay_0to10_80to90.pdf");
}
