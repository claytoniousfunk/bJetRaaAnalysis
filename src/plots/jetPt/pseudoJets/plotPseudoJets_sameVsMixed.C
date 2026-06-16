// Compare same-event vs. mixed-event pseudo-jet pT spectra for all 16
// ultra-fine centrality bins.  Per-jet normalised (shape comparison).
//
// Fill in f_mixed_path once the hadd is complete.

const char *f_same_path =
  "/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/scanningOutput/HYDJET/canonical/pfCandAnalyzer/"
  "HYDJET_pThat-unweighted_mu12_pTmu-15to999_tight_hiBinShift-0_jetTrkMaxFilter_WDecayFilter_"
  "sameEventPseudoJets_pseudoJetCandPtMin-0.0_subleadingPFCandPtMin-15_2026-6-9_ultraFineCentBins_pfCandAnalyzer.root";

const char *f_mixed_path =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/HYDJET/"
  "HYDJET_pfCandAnalyzer_mixedEventPseudoJets.root";  // <-- fill in after hadd

// centrality labels for C1..C16 (hiBin step = 10 → 5% steps)
const char *centLabel[17] = {
  "",
  "0-5%",   "5-10%",  "10-15%", "15-20%",
  "20-25%", "25-30%", "30-35%", "35-40%",
  "40-45%", "45-50%", "50-55%", "55-60%",
  "60-65%", "65-70%", "70-75%", "75-80%"
};

TString outDir = "/home/clayton/Analysis/code/bJetRaaAnalysis/figures/jetPt/pseudoJets/sameVsMixed/";

void plotPseudoJets_sameVsMixed(){

  gSystem->Exec(Form("mkdir -p %s", outDir.Data()));

  TFile *fSame  = TFile::Open(f_same_path);
  TFile *fMixed = TFile::Open(f_mixed_path);

  if(!fSame  || fSame ->IsZombie()){ printf("ERROR: cannot open same-event file\n");  return; }
  if(!fMixed || fMixed->IsZombie()){ printf("ERROR: cannot open mixed-event file\n"); return; }

  double lm = 0.16, bm = 0.10, tm = 0.12, rm = 0.05;
  double xmin = 0., xmax = 200.;

  TLine *li = new TLine();
  li->SetLineStyle(2); li->SetLineColor(kGray+1); li->SetLineWidth(1);

  for(int iC = 1; iC <= 16; iC++){

    TH1D *hS, *hM;
    fSame ->GetObject(Form("h_pseudoJetPt_C%d", iC), hS);
    fMixed->GetObject(Form("h_pseudoJetPt_C%d", iC), hM);

    if(!hS){ printf("WARNING: h_pseudoJetPt_C%d not found in same-event file\n",  iC); continue; }
    if(!hM){ printf("WARNING: h_pseudoJetPt_C%d not found in mixed-event file\n", iC); continue; }

    // per-jet shape normalisation
    TH1D *hSn = (TH1D*) hS->Clone(Form("hSn_C%d", iC));
    TH1D *hMn = (TH1D*) hM->Clone(Form("hMn_C%d", iC));
    if(hSn->Integral() > 0) hSn->Scale(1./hSn->Integral());
    if(hMn->Integral() > 0) hMn->Scale(1./hMn->Integral());

    // ratio: mixed / same
    TH1D *hR = (TH1D*) hMn->Clone(Form("hR_C%d", iC));
    hR->Divide(hSn);

    // style
    hSn->SetLineColor(kBlack);   hSn->SetLineWidth(2); hSn->SetStats(0);
    hMn->SetLineColor(kRed+1);   hMn->SetLineWidth(2); hMn->SetStats(0);
    hMn->SetLineStyle(2);
    hR ->SetLineColor(kBlack);   hR ->SetLineWidth(2); hR ->SetStats(0);
    hR ->SetMarkerColor(kBlack); hR ->SetMarkerStyle(20); hR->SetMarkerSize(0.7);

    // canvas
    TCanvas *c = new TCanvas(Form("c_C%d",iC),"",700,800);

    TPad *pUp = new TPad(Form("pUp_C%d",iC),"",0,0.35,1,1);
    pUp->SetLeftMargin(lm); pUp->SetRightMargin(rm);
    pUp->SetBottomMargin(0.); pUp->SetTopMargin(tm);
    pUp->SetLogy(); pUp->Draw(); pUp->cd();

    TString yTitle = "#frac{1}{#it{N}^{jet}} #frac{d#it{N}^{jet}}{d#it{p}_{T}} [GeV^{-1}]";

    hSn->SetTitle("");
    hSn->GetXaxis()->SetRangeUser(xmin,xmax);
    hSn->GetYaxis()->SetTitle(yTitle);
    hSn->GetYaxis()->SetTitleSize(0.052);
    hSn->GetYaxis()->SetLabelSize(0.048);
    hSn->GetYaxis()->SetTitleOffset(1.3);
    hSn->Draw("hist");
    hMn->Draw("hist same");

    TLegend *leg = new TLegend(0.42,0.65,0.92,0.82);
    leg->SetBorderSize(0); leg->SetTextSize(0.046);
    leg->AddEntry(hSn, "Same-event pseudo-jets",  "l");
    leg->AddEntry(hMn, "Mixed-event pseudo-jets", "l");
    leg->Draw();

    TLatex *lat = new TLatex();
    lat->SetNDC(); lat->SetTextSize(0.050);
    lat->DrawLatex(lm+0.02, 1.-tm+0.03, Form("PbPb 5.02 TeV  %s", centLabel[iC]));
    lat->DrawLatex(lm+0.02, 1.-tm-0.03, "HYDJET, anti-#it{k}_{T} #it{R} = 0.4 random cone");

    c->cd();
    TPad *pDn = new TPad(Form("pDn_C%d",iC),"",0,0,1,0.35);
    pDn->SetLeftMargin(lm); pDn->SetRightMargin(rm);
    pDn->SetTopMargin(0.); pDn->SetBottomMargin(0.28);
    pDn->Draw(); pDn->cd();

    hR->SetTitle("");
    hR->GetXaxis()->SetRangeUser(xmin,xmax);
    hR->GetXaxis()->SetTitle("Pseudo-jet #it{p}_{T} [GeV]");
    hR->GetYaxis()->SetTitle("Mixed / Same");
    hR->GetYaxis()->SetRangeUser(0.0, 2.0);
    hR->GetXaxis()->SetTitleSize(0.105);
    hR->GetXaxis()->SetLabelSize(0.090);
    hR->GetYaxis()->SetTitleSize(0.095);
    hR->GetYaxis()->SetLabelSize(0.085);
    hR->GetYaxis()->SetTitleOffset(0.65);
    hR->GetYaxis()->SetNdivisions(504);
    hR->Draw("ep");
    li->DrawLine(xmin,1.,xmax,1.);

    TString outFile = outDir + Form("pseudoJetPt_sameVsMixed_C%d.pdf", iC);
    c->SaveAs(outFile);
    printf("Saved C%d (%s)\n", iC, centLabel[iC]);

    delete c;
  }

  printf("\nDone. Output in %s\n", outDir.Data());
}
