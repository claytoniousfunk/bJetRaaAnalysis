// Estimate the per-event fake-jet rate from the mixed-event pseudo-jet
// distribution (CS PF candidates).
//
// Formula:
//   dN_fake/dpT per event = h_raw / (N_events * N_cones * bw) * N_indep
//
// Cone centers are thrown uniformly in eta [-1.6, 1.6], full phi
// (HYDJET_pfCandAnalyzer.C line 1350), so:
//   N_indep = (delta_eta * 2pi) / (pi * R^2) = (3.2 * 2pi) / (pi * 0.16) = 40

const char *fmixed_path =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/HYDJET/"
  "HYDJET_pThat-unweighted_mu12_pTmu-15to999_tight_hiBinShift-0_jetTrkMaxFilter_WDecayFilter_"
  "mixedEventPseudoJets_pfCandCS_pseudoJetCandPtMin-0.0_subleadingPFCandPtMin-15_"
  "2026-6-16_ultraFineCentBins.root";

const int    N_cones   = 100;
const double R         = 0.4;
const double eta_range = 3.2;  // cone centers in [-1.6, 1.6]
const double N_indep   = (eta_range * 2. * TMath::Pi()) / (TMath::Pi() * R * R);

const char *centLabel[17] = {
  "",
  "0-5%",   "5-10%",  "10-15%", "15-20%",
  "20-25%", "25-30%", "30-35%", "35-40%",
  "40-45%", "45-50%", "50-55%", "55-60%",
  "60-65%", "65-70%", "70-75%", "75-80%"
};

TString outDir = "/home/clayton/Analysis/code/bJetRaaAnalysis/figures/jetPt/pseudoJets/fakeJetRate/";

int centColor(int iC){
  float t = (iC - 1) / 15.f;
  return TColor::GetColor(1.f - t, 0.f, t);
}

// Build differential fake-jet rate histogram from a raw histogram
TH1D* makeDiff(TH1D *hRaw, const char *name){
  double N_events = hRaw->GetEntries() / double(N_cones);
  double norm     = N_indep / (N_events * N_cones);
  TH1D *h = (TH1D*) hRaw->Clone(name);
  h->Scale(norm);
  for(int b = 1; b <= h->GetNbinsX(); b++){
    double bw = h->GetBinWidth(b);
    h->SetBinContent(b, h->GetBinContent(b) / bw);
    h->SetBinError  (b, h->GetBinError(b)   / bw);
  }
  h->SetStats(0);
  return h;
}

void styleUp(TH1D *h, double lm, double tm){
  h->SetTitle("");
  h->GetXaxis()->SetRangeUser(0., 100.);
  h->GetYaxis()->SetTitle("d#it{N}_{fake}/d#it{p}_{T} per event [GeV^{-1}]");
  h->GetYaxis()->SetTitleSize(0.052);
  h->GetYaxis()->SetLabelSize(0.048);
  h->GetYaxis()->SetTitleOffset(1.4);
}

void styleDn(TH1D *h){
  h->SetTitle("");
  h->GetXaxis()->SetRangeUser(0., 100.);
  h->GetXaxis()->SetTitle("Pseudo-jet #it{p}_{T} [GeV]");
  h->GetYaxis()->SetTitle("Cent. / peripheral");
  h->GetYaxis()->SetRangeUser(0., 25.);
  h->GetXaxis()->SetTitleSize(0.105);
  h->GetXaxis()->SetLabelSize(0.090);
  h->GetYaxis()->SetTitleSize(0.095);
  h->GetYaxis()->SetLabelSize(0.085);
  h->GetYaxis()->SetTitleOffset(0.65);
  h->GetYaxis()->SetNdivisions(505);
  h->SetStats(0);
}

void plotPseudoJets_fakeJetRate(){

  gSystem->Exec(Form("mkdir -p %s", outDir.Data()));

  TFile *fM = TFile::Open(fmixed_path);
  if(!fM || fM->IsZombie()){ printf("ERROR: cannot open file\n"); return; }

  printf("N_indep = %.0f  (eta_range=%.1f, R=%.1f)\n", N_indep, eta_range, R);

  double lm = 0.16, tm = 0.10, rm = 0.05;

  TLine *unity = new TLine();
  unity->SetLineStyle(2); unity->SetLineColor(kGray+1); unity->SetLineWidth(1);

  // -----------------------------------------------------------------------
  // Plot 1: all 16 fine centrality bins, ratio to C16 (75-80%)
  // -----------------------------------------------------------------------
  {
    // Build all differential histograms first so we have the reference (C16)
    TH1D *hDiff[17] = {};
    for(int iC = 1; iC <= 16; iC++){
      TH1D *hRaw;
      fM->GetObject(Form("h_pseudoJetPt_C%d", iC), hRaw);
      if(!hRaw){ printf("WARNING: missing C%d\n", iC); continue; }
      hDiff[iC] = makeDiff(hRaw, Form("hDiff_C%d", iC));
      hDiff[iC]->SetLineColor(centColor(iC));
      hDiff[iC]->SetLineWidth(2);
    }
    TH1D *hRef = hDiff[16];  // 75-80% reference

    TCanvas *c = new TCanvas("c_diff","",700,900);
    TPad *pUp = new TPad("pUp","",0,0.35,1,1);
    pUp->SetLeftMargin(lm); pUp->SetRightMargin(rm);
    pUp->SetBottomMargin(0.); pUp->SetTopMargin(tm);
    pUp->SetLogy(); pUp->Draw(); pUp->cd();

    styleUp(hDiff[1], lm, tm);
    hDiff[1]->Draw("hist");
    for(int iC = 2; iC <= 16; iC++) if(hDiff[iC]) hDiff[iC]->Draw("hist same");

    TLegend *leg = new TLegend(0.55,0.35,0.92,0.88);
    leg->SetBorderSize(0); leg->SetTextSize(0.040); leg->SetNColumns(2);
    for(int iC = 1; iC <= 16; iC++) if(hDiff[iC]) leg->AddEntry(hDiff[iC], centLabel[iC], "l");
    leg->Draw();

    TLatex *lat = new TLatex(); lat->SetNDC(); lat->SetTextSize(0.050);
    lat->DrawLatex(lm+0.02, 1.-tm+0.02, "PbPb 5.02 TeV  HYDJET");
    lat->DrawLatex(lm+0.02, 1.-tm-0.06,
      Form("Mixed-event, anti-#it{k}_{T} #it{R}=0.4 (#it{N}_{indep}=%.0f)", N_indep));

    c->cd();
    TPad *pDn = new TPad("pDn","",0,0,1,0.35);
    pDn->SetLeftMargin(lm); pDn->SetRightMargin(rm);
    pDn->SetTopMargin(0.); pDn->SetBottomMargin(0.28);
    pDn->Draw(); pDn->cd();

    TH1D *hRatioFrame = nullptr;
    for(int iC = 1; iC <= 15; iC++){
      if(!hDiff[iC] || !hRef) continue;
      TH1D *hr = (TH1D*) hDiff[iC]->Clone(Form("hRatio_C%d", iC));
      hr->Divide(hRef);
      if(!hRatioFrame){
        hRatioFrame = hr;
        styleDn(hRatioFrame);
        hRatioFrame->Draw("hist");
      } else {
        hr->Draw("hist same");
      }
    }
    unity->DrawLine(0., 1., 100., 1.);

    c->SaveAs(outDir + "fakeJetRate_differential.pdf");
    printf("Saved: fakeJetRate_differential.pdf\n");
  }

  // -----------------------------------------------------------------------
  // Plot 2: coarse centrality bins, ratio to 50-80%
  //   0-10%  : C1+C2
  //   10-30% : C3-C6
  //   30-50% : C7-C10
  //   50-80% : C11-C16  (reference)
  // -----------------------------------------------------------------------
  {
    struct CoarseBin { int first, last; const char *label; int color; };
    CoarseBin cb[] = {
      {1,  2,  "0-10%",  kRed+1},
      {3,  6,  "10-30%", kOrange+7},
      {7,  10, "30-50%", kGreen+2},
      {11, 16, "50-80%", kBlue+1}
    };
    const int NCB = 4;

    TH1D *hCoarse[4] = {};
    for(int ic = 0; ic < NCB; ic++){
      TH1D *hSum = nullptr;
      double totalEntries = 0.;
      for(int iC = cb[ic].first; iC <= cb[ic].last; iC++){
        TH1D *hRaw;
        fM->GetObject(Form("h_pseudoJetPt_C%d", iC), hRaw);
        if(!hRaw) continue;
        totalEntries += hRaw->GetEntries();
        if(!hSum) hSum = (TH1D*) hRaw->Clone(Form("hCoarseRaw_%d", ic));
        else      hSum->Add(hRaw);
      }
      if(!hSum) continue;
      // override GetEntries (lost after Add); use totalEntries
      double N_events = totalEntries / double(N_cones);
      double norm     = N_indep / (N_events * N_cones);
      hCoarse[ic] = (TH1D*) hSum->Clone(Form("hCoarse_%d", ic));
      hCoarse[ic]->Scale(norm);
      for(int b = 1; b <= hCoarse[ic]->GetNbinsX(); b++){
        double bw = hCoarse[ic]->GetBinWidth(b);
        hCoarse[ic]->SetBinContent(b, hCoarse[ic]->GetBinContent(b) / bw);
        hCoarse[ic]->SetBinError  (b, hCoarse[ic]->GetBinError(b)   / bw);
      }
      hCoarse[ic]->SetLineColor(cb[ic].color);
      hCoarse[ic]->SetLineWidth(2);
      hCoarse[ic]->SetStats(0);
      delete hSum;
    }
    TH1D *hRef = hCoarse[3];  // 50-80% reference

    TCanvas *c2 = new TCanvas("c_coarse","",700,900);
    TPad *pUp = new TPad("pUp2","",0,0.35,1,1);
    pUp->SetLeftMargin(lm); pUp->SetRightMargin(rm);
    pUp->SetBottomMargin(0.); pUp->SetTopMargin(tm);
    pUp->SetLogy(); pUp->Draw(); pUp->cd();

    TH1D *hFrame = hCoarse[0];
    styleUp(hFrame, lm, tm);
    hFrame->Draw("hist");
    for(int ic = 1; ic < NCB; ic++) if(hCoarse[ic]) hCoarse[ic]->Draw("hist same");

    TLegend *leg2 = new TLegend(0.55,0.68,0.92,0.88);
    leg2->SetBorderSize(0); leg2->SetTextSize(0.048);
    for(int ic = 0; ic < NCB; ic++) if(hCoarse[ic]) leg2->AddEntry(hCoarse[ic], cb[ic].label, "l");
    leg2->Draw();

    TLatex *lat2 = new TLatex(); lat2->SetNDC(); lat2->SetTextSize(0.050);
    lat2->DrawLatex(lm+0.02, 1.-tm+0.02, "PbPb 5.02 TeV  HYDJET");
    lat2->DrawLatex(lm+0.02, 1.-tm-0.06,
      Form("Mixed-event, anti-#it{k}_{T} #it{R}=0.4 (#it{N}_{indep}=%.0f)", N_indep));

    c2->cd();
    TPad *pDn = new TPad("pDn2","",0,0,1,0.35);
    pDn->SetLeftMargin(lm); pDn->SetRightMargin(rm);
    pDn->SetTopMargin(0.); pDn->SetBottomMargin(0.28);
    pDn->Draw(); pDn->cd();

    TH1D *hRatioFrame = nullptr;
    for(int ic = 0; ic < NCB-1; ic++){
      if(!hCoarse[ic] || !hRef) continue;
      TH1D *hr = (TH1D*) hCoarse[ic]->Clone(Form("hCoarseRatio_%d", ic));
      hr->Divide(hRef);
      if(!hRatioFrame){
        hRatioFrame = hr;
        styleDn(hRatioFrame);
        hRatioFrame->Draw("hist");
      } else {
        hr->Draw("hist same");
      }
    }
    unity->DrawLine(0., 1., 100., 1.);

    c2->SaveAs(outDir + "fakeJetRate_differential_coarseCent.pdf");
    printf("Saved: fakeJetRate_differential_coarseCent.pdf\n");
  }

  printf("\nDone. N_indep=%.0f  Output in %s\n", N_indep, outDir.Data());
}
