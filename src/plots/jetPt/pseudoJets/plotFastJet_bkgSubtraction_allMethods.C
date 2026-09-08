// Compare RC-subtracted vs. unsubtracted FastJet anti-kT R=0.4 pT spectra.
// For each centrality bin, produces a two-pad canvas:
//   upper: h_fastJetPt and h_fastJetPt_bkgSub_RC overlaid (log-y)
//   lower: ratio  h_fastJetPt_bkgSub_RC / h_fastJetPt
//
// Usage (from src/plots/jetPt/pseudoJets/):
//   root -l -b -q plotFastJet_rcSubtraction.C

const char *inFile =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/PbPb_MinBias_Part1_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_sameEventPFClustering_pseudoJetCandPtMin-0.0_2026-8-6_ultraFineCentBins.root";

const char *outDir =
  "../../../../figures/jetPt/pseudoJets/bkgSubtraction/";

// ultraFine scheme: C0 = inclusive 0-80%, C1-C16 = 5% bins
const int NCentBins = 17;
const char *centLabel[NCentBins] = {
  "0-80% (incl.)",
  "0-5%",   "5-10%",  "10-15%", "15-20%",
  "20-25%", "25-30%", "30-35%", "35-40%",
  "40-45%", "45-50%", "50-55%", "55-60%",
  "60-65%", "65-70%", "70-75%", "75-80%"
};
const char *centSuffix[NCentBins] = {
  "C0",
  "C1",  "C2",  "C3",  "C4",
  "C5",  "C6",  "C7",  "C8",
  "C9",  "C10", "C11", "C12",
  "C13", "C14", "C15", "C16"
};

const double pTlo = 0.;
const double pThi = 300.;

const int colRaw  = kBlue  - 4;
const int colSub  = kRed   - 4;
const int colRC = kGray + 1;


static void styleHist(TH1D *h, int col, double lw = 2.)
{
  h->SetLineColor(col);
  h->SetLineWidth(lw);
  h->SetStats(0);
  h->SetTitle("");
  h->GetXaxis()->SetRangeUser(pTlo, pThi);
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetXaxis()->SetLabelSize(0.04);
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetLabelSize(0.04);
  h->GetYaxis()->SetTitleOffset(1.55);
}

void plotFastJet_bkgSubtraction_allMethods()
{
  gSystem->Exec(Form("mkdir -p %s", outDir));

  TFile *f = TFile::Open(inFile);
  if(!f || f->IsZombie()){
    std::cerr << "ERROR: cannot open " << inFile << "\n";
    return;
  }


  TLine *li = new TLine();
  li->SetLineStyle(7);
  li->SetLineColor(kGray + 1);

  for(int ci = 0; ci < NCentBins; ci++){

    
    TH1D *hSub_RC = nullptr;
    TH1D *hSub_RC_geoCorr = nullptr;
    TH1D *hSub_RC_geoCorr_etaReflect = nullptr;
        
    TH1D *hvz = nullptr;
    TH1D *hCS = nullptr;
    
    
    // f->GetObject(Form("h_fastJetPt_PF_JEC_bkgSub_RC_%s", centSuffix[ci]), hSub_RC);
    // f->GetObject(Form("h_fastJetPt_PF_JEC_bkgSub_RC_geoCorr_%s", centSuffix[ci]), hSub_RC_geoCorr);
    // f->GetObject(Form("h_fastJetPt_PF_JEC_bkgSub_RC_geoCorr_etaReflect_%s", centSuffix[ci]), hSub_RC_geoCorr_etaReflect);

    f->GetObject(Form("h_fastJetPt_PF_bkgSub_RC_%s", centSuffix[ci]), hSub_RC);
    f->GetObject(Form("h_fastJetPt_PF_bkgSub_RC_geoCorr_%s", centSuffix[ci]), hSub_RC_geoCorr);
    f->GetObject(Form("h_fastJetPt_PF_bkgSub_RC_geoCorr_etaReflect_%s", centSuffix[ci]), hSub_RC_geoCorr_etaReflect);

    // f->GetObject(Form("h_fastJetPt_PF_bkgSub_dPT_%s", centSuffix[ci]), hSub_RC);
    // f->GetObject(Form("h_fastJetPt_PF_bkgSub_dPT_dPTAbove0_%s", centSuffix[ci]), hSub_RC_geoCorr);
    // f->GetObject(Form("h_fastJetPt_PF_bkgSub_dPT_PFCsPTAbove60_%s", centSuffix[ci]), hSub_RC_geoCorr_etaReflect);

    //f->GetObject(Form("h_inclRecoJetPt_%s", centSuffix[ci]), hCS);
    f->GetObject(Form("h_inclRawJetPt_%s", centSuffix[ci]), hCS);
    f->GetObject(Form("h_vz_%s",centSuffix[ci]),hvz);
    

    
    if(!hSub_RC || !hSub_RC_geoCorr || !hSub_RC_geoCorr_etaReflect){
      std::cerr << "WARNING: histograms not found for " << centSuffix[ci] << " — skipping\n";
      continue;
    }



    if(!hCS){
      std::cerr << "WARNING: h_inclRecoJetPt_" << centSuffix[ci] << " not found — skipping CS comparison\n";
      continue;
    }
    hCS = (TH1D*) hCS->Clone(Form("hCS_%s", centSuffix[ci]));
    hCS->SetDirectory(nullptr);

    
    hSub_RC = (TH1D*) hSub_RC->Clone(Form("hSub_RC_%s", centSuffix[ci]));
    
    hSub_RC->SetDirectory(nullptr);
    
    styleHist(hSub_RC, kRed-4);
    styleHist(hSub_RC_geoCorr, kBlue-4);
    styleHist(hSub_RC_geoCorr_etaReflect, kGreen+2);

  

    // ---- canvas 2: RC-subtracted FastJet vs CS reco jets (shape comparison) --


    // normalize to unity for shape comparison
    TH1D *hSubNorm_RC = (TH1D*) hSub_RC->Clone(Form("hSubNorm_RC_%s", centSuffix[ci]));
    TH1D *hSubNorm_RC_geoCorr = (TH1D*) hSub_RC_geoCorr->Clone(Form("hSubNorm_RC_geoCorr_%s", centSuffix[ci]));
    TH1D *hSubNorm_RC_geoCorr_etaReflect = (TH1D*) hSub_RC_geoCorr_etaReflect->Clone(Form("hSubNorm_RC_geoCorr_etaReflect_%s", centSuffix[ci]));


    
    TH1D *hCSNorm  = (TH1D*) hCS ->Clone(Form("hCSNorm_%s",  centSuffix[ci]));
    hSubNorm_RC->SetDirectory(nullptr);
    hSubNorm_RC_geoCorr->SetDirectory(nullptr);
    hSubNorm_RC_geoCorr_etaReflect->SetDirectory(nullptr);
    hCSNorm ->SetDirectory(nullptr);

    if(hvz->Integral() > 0) {
      hSubNorm_RC->Scale(1./hvz->Integral());
      hSubNorm_RC_geoCorr->Scale(1./hvz->Integral());
      hSubNorm_RC_geoCorr_etaReflect->Scale(1./hvz->Integral());

      hCSNorm->Scale(1./hvz->Integral());
    }
   

    const int colCS = kBlack;
    styleHist(hSubNorm_RC, colSub);
    styleHist(hCSNorm,  colCS);

    TCanvas *c2 = new TCanvas(Form("c2_%s", centSuffix[ci]), "", 700, 700);
    TPad *pUp2 = new TPad(Form("pUp2_%s", centSuffix[ci]), "", 0, 0.35, 1, 1);
    TPad *pDn2 = new TPad(Form("pDn2_%s", centSuffix[ci]), "", 0, 0,    1, 0.35);
    pUp2->SetLeftMargin(0.15);  pUp2->SetRightMargin(0.05);
    pUp2->SetTopMargin(0.10);   pUp2->SetBottomMargin(0.);
    pDn2->SetLeftMargin(0.15);  pDn2->SetRightMargin(0.05);
    pDn2->SetTopMargin(0.);     pDn2->SetBottomMargin(0.25);
    pUp2->SetLogy();
    pUp2->Draw(); pDn2->Draw();

    pUp2->cd();
    double ymax2 = std::max(hSubNorm_RC->GetMaximum(), hCSNorm->GetMaximum());
    double ymin2 = 1e30;
    for(int b = 1; b <= hSubNorm_RC->GetNbinsX(); b++){
      double v = hSubNorm_RC->GetBinContent(b);
      if(v > 0 && v < ymin2) ymin2 = v;
    }
    if(ymin2 > 1e29) ymin2 = 1e-8;

    hSubNorm_RC->GetXaxis()->SetRangeUser(pTlo, pThi);
    hSubNorm_RC->GetYaxis()->SetRangeUser(ymin2 * 0.05, ymax2 * 50.);
    hSubNorm_RC->GetYaxis()->SetTitle("Entries per-event");
    hSubNorm_RC->GetXaxis()->SetLabelSize(0.);
    hSubNorm_RC->Draw("hist");
    hSubNorm_RC_geoCorr->Draw("hist same");
    hSubNorm_RC_geoCorr_etaReflect->Draw("hist same");
    hCSNorm ->Draw("hist same");

    TLatex *la = new TLatex();
    la->SetTextSize(0.055);
    la->DrawLatexNDC(0.68, 0.82, Form("PbPb, %s", centLabel[ci]));
    la->SetTextSize(0.04);
    //la.SetTextColor(colSub);
    //la.DrawLatex(0.18, 0.76, "FastJet (RC subtracted)");
    //la.DrawLatex(0.18, 0.76, "FastJet (RC subtracted + JEC)");
    //la.DrawLatex(0.18, 0.76, "FastJet (#Delta p_{T} subtracted)");
    //la.DrawLatex(0.18, 0.76, "FastJet (#Delta p_{T} subtracted + JEC)");
    //la.SetTextColor(colCS);
    //la.DrawLatex(0.18, 0.68, "CS reco jets (akCs4PF raw)");
    //la.DrawLatex(0.18, 0.68, "CS reco jets (akCs4PF raw + JEC)");

    la->SetTextColor(kBlack);
    la->DrawLatexNDC(0.5,0.7,"ak4PFCs Jets (raw)");
    la->SetTextColor(kRed-4);
    la->DrawLatexNDC(0.5,0.62,"fastJet - RandomCone");
    la->SetTextColor(kBlue-4);
    la->DrawLatexNDC(0.5,0.54,"fastJet - RandomCone+geoCorr");
    la->SetTextColor(kGreen+2);
    la->DrawLatexNDC(0.5,0.46,"fastJet - RandomCone+geoCorr+#it{#eta}-reflect");

    

    pDn2->cd();
    TH1D *hRatio_RC = (TH1D*) hSubNorm_RC->Clone(Form("hRatio_RC_%s", centSuffix[ci]));
    TH1D *hRatio_RC_geoCorr = (TH1D*) hSubNorm_RC_geoCorr->Clone(Form("hRatio_RC_geoCorr_%s", centSuffix[ci]));
    TH1D *hRatio_RC_geoCorr_etaReflect = (TH1D*) hSubNorm_RC_geoCorr_etaReflect->Clone(Form("hRatio_RC_geoCorr_etaReflect_%s", centSuffix[ci]));

    
    hRatio_RC->Divide(hCSNorm);
    hRatio_RC_geoCorr->Divide(hCSNorm);
    hRatio_RC_geoCorr_etaReflect->Divide(hCSNorm);


    
    hRatio_RC->SetDirectory(nullptr);
    hRatio_RC->GetXaxis()->SetRangeUser(pTlo, pThi);
    hRatio_RC->GetXaxis()->SetTitle("Raw Jet #it{p}_{T} [GeV]");
    hRatio_RC->GetYaxis()->SetTitle("Sub. / CS");
    hRatio_RC->GetYaxis()->SetRangeUser(0., 2.);
    hRatio_RC->GetXaxis()->SetTitleSize(0.10);
    hRatio_RC->GetXaxis()->SetLabelSize(0.09);
    hRatio_RC->GetYaxis()->SetTitleSize(0.10);
    hRatio_RC->GetYaxis()->SetLabelSize(0.09);
    hRatio_RC->GetYaxis()->SetTitleOffset(0.65);
    hRatio_RC->GetYaxis()->SetNdivisions(505);
    hRatio_RC->Draw("hist");
    hRatio_RC_geoCorr->Draw("hist same");
    hRatio_RC_geoCorr_etaReflect->Draw("hist same");
    li->DrawLine(pTlo, 1., pThi, 1.);

    c2->SaveAs(Form("%sfastJetPt_rcSub_vs_csJet_%s.pdf", outDir, centSuffix[ci]));
    delete c2;
  }

  f->Close();
  std::cout << "Figures saved to " << outDir << "\n";
}
