// Jets per Z, calo jets vs PF jets, one pad per centrality class.
//
// Both curves come from calculateJetsPerZ_caloJets.C run on the two file sets,
// so the chain, stitching windows, Z window and muon efficiencies are identical
// and the only difference is the jet collection. Neither is unfolded and
// neither has the fake-jet subtraction (there is no calo fake-jet estimate).
//
// The PF curve is therefore NOT calculateRAA.C's nominal, which unfolds and
// subtracts fake jets. It is the like-for-like twin of the calo curve.
//
// Below 150 GeV both sides come from MinBias, and on the calo side the pp
// MinBias scan has an 80 GeV threshold, so the lowest points carry its 7%
// normalization statistics -- see the header of calculateJetsPerZ_caloJets.C.
//
// Usage, from src/newFractionCalculation/:
//   root -l -b -q 'plotJetsPerZ_caloVsPF_noUnfold.C'

#include "../../headers/plotting/plotStyle.h"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TSystem.h"

// Trigger-only inputs: MinBias dropped from both systems, so each starts at the
// Jet80 threshold (150 GeV). The MinBias-fed region below 150 was pure scatter
// -- 50-80% ran 0.42 / 0.52 / 0.46 / 0.66 / 0.96 with 11-22% errors, against
// 2.0% errors the moment Jet80 takes over. Produced by
//   calculateJetsPerZ_caloJets.C(true ,false,false,false)   calo
//   calculateJetsPerZ_caloJets.C(false,false,false,false)   PF
const char *f_calo_path = "./rootFiles/JetsPerZ/histograms_JetsPerZ_caloJets_trigOnly_noUnfold.root";
const char *f_PF_path   = "./rootFiles/JetsPerZ/histograms_JetsPerZ_PFJets_trigOnly_noUnfold.root";
const char *outPath     = "/home/clayton/Analysis/code/bJetRaaAnalysis/figures/JetsPerZ/JetsPerZ_caloVsPF_trigOnly_noUnfold.pdf";

void plotJetsPerZ_caloVsPF_noUnfold()
{
  initPlotStyle();

  TFile *fC = TFile::Open(f_calo_path);
  TFile *fP = TFile::Open(f_PF_path);
  if(!fC || fC->IsZombie() || !fP || fP->IsZombie()){
    printf("ERROR: run calculateJetsPerZ_caloJets.C for both file sets first\n");
    return;
  }

  const char *clsLabel[5] = {"", "0-10%", "10-30%", "30-50%", "50-80%"};

  TCanvas *c = new TCanvas("c_caloVsPF", "", 700, 800);
  c->Divide(2, 2, 0.001, 0.001);

  printf("\n  jets per Z, calo / PF, per class\n");

  TLine *one = new TLine();
  one->SetLineStyle(7);

  for(int cls = 1; cls <= 4; cls++){
    TH1D *hC = nullptr, *hP = nullptr;
    fC->GetObject(Form("r_C%d_fine", cls), hC);
    fP->GetObject(Form("r_C%d_fine", cls), hP);
    if(!hC || !hP){ printf("ERROR: r_C%d_fine missing\n", cls); return; }

    c->cd(cls);
    gPad->SetLeftMargin(0.17); gPad->SetBottomMargin(0.13); gPad->SetTopMargin(0.09); gPad->SetRightMargin(0.04);

    styleH(hP, hexData,      markOpenCircle);
    styleH(hC, hexCorrected, markFilledSquare);

    hP->SetTitle("");
    hP->SetMinimum(0.); hP->SetMaximum(1.4);
    hP->GetXaxis()->SetTitle("jet #it{p}_{T} [GeV]");
    hP->GetXaxis()->SetTitleSize(0.055); hP->GetXaxis()->SetLabelSize(0.050); hP->GetXaxis()->SetTitleOffset(1.05);
    hP->GetYaxis()->SetTitle("jets per #it{Z}, PbPb / pp");
    hP->GetYaxis()->SetTitleSize(0.050); hP->GetYaxis()->SetLabelSize(0.050); hP->GetYaxis()->SetTitleOffset(1.45);
    hP->Draw("E1");
    one->DrawLine(hP->GetXaxis()->GetXmin(), 1., hP->GetXaxis()->GetXmax(), 1.);
    hP->Draw("E1 same");
    hC->Draw("E1 same");

    TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.055);
    la.DrawLatex(0.22, 0.85, Form("PbPb %s", clsLabel[cls]));

    if(cls == 1){
      TLegend *leg = makeLegend(0.46, 0.68, 0.94, 0.84, 0.050);
      leg->AddEntry(hC, "calo jets", "lp");
      leg->AddEntry(hP, "PF jets", "lp");
      leg->Draw();
    }

    // where both are statistically solid, above the MinBias region
    double a = hC->GetBinContent(hC->FindBin(250.)), b = hP->GetBinContent(hP->FindBin(250.));
    printf("    %-7s at 250 GeV: calo %.3f   PF %.3f   calo/PF %.3f\n", clsLabel[cls], a, b, b > 0. ? a/b : -1.);
  }

  c->SaveAs(outPath);
  printf("\n  figure: %s\n", outPath);
}
