// Compare mixed-event FastJet anti-kT clusters vs random-cone pseudo-jets
// from the PbPb MinBias pfCandAnalyzer scan.
//
// Both distributions come from the same mixed-event file, so their ratio
// (FastJet / RandomCone) shows how much the anti-kT algorithm finds
// pT-concentrated structure above pure uniform UE sampling.
//
// Normalization: sumW = h_pseudoJetPt->Integral() / N_pool  (same proxy
// used in plotJetPt_fakeJetSubtraction.C and plotPseudoJets_fakeJetRate.C)

const char* inFile =
    "../../../../rootFiles/scanningOuput/PbPb/"
    "PbPb_MinBias_Part1_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_"
    "mixedEventPseudoJets_pfCand_pseudoJetCandPtMin-0.0_2026-7-7.root";

const char* outDir =
    "../../../../figures/jetPt/pseudoJets/fakeJetAnalysis_MinBias/";

const int    N_pool   = 100;   // N_mixedEventsInPool
const double pTlo     = 0.;
const double pThi     = 200.;

// CENT_NOMINAL 4 bins
const int NCent = 5;  // 0=inclusive, 1-4 = centrality bins
const char* centLabel[NCent] = { "0-80%", "0-10%", "10-30%", "30-50%", "50-80%" };

// Okabe-Ito colorblind-safe colors
const Int_t col_fastJet  = TColor::GetColor("#0072B2");  // blue
const Int_t col_randCone = TColor::GetColor("#E69F00");  // orange
const Int_t col_ratio    = TColor::GetColor("#009E73");  // green

void plotOneBin(TFile* f, int ci){

    TH1D* hFJ = nullptr;
    TH1D* hRC = nullptr;
    TH1D* hvz = nullptr;
    f->GetObject(Form("h_fastJetPt_C%d",   ci), hFJ);
    f->GetObject(Form("h_pseudoJetPt_C%d", ci), hRC);
    f->GetObject(Form("h_vz_C%d", ci), hvz);
    if(!hFJ || !hRC){
        printf("ERROR: missing histograms for C%d\n", ci); return;
    }
    hFJ = (TH1D*) hFJ->Clone(Form("hFJ_C%d", ci)); hFJ->SetDirectory(nullptr);
    hRC = (TH1D*) hRC->Clone(Form("hRC_C%d", ci)); hRC->SetDirectory(nullptr);

    // normalize by effective event count and bin width -> dN/dpT per event
    double sumW = hRC->Integral();
    if(sumW <= 0){ printf("WARNING: sumW=0 for C%d, skipping\n", ci); return; }
    hFJ->Scale(1. / hvz->Integral());
    hRC->Scale(1. / (hvz->Integral()*double(N_pool)));
    for(int b = 1; b <= hFJ->GetNbinsX(); b++){
        double bw = hFJ->GetBinWidth(b);
        hFJ->SetBinContent(b, hFJ->GetBinContent(b) / bw);
        hFJ->SetBinError  (b, hFJ->GetBinError(b)   / bw);
        hRC->SetBinContent(b, hRC->GetBinContent(b) / bw);
        hRC->SetBinError  (b, hRC->GetBinError(b)   / bw);
    }

    // ratio FastJet / RandomCone
    TH1D* hRat = (TH1D*) hFJ->Clone(Form("hRat_C%d", ci));
    hRat->Divide(hRC);
    hRat->SetDirectory(nullptr);

    // y range for spectra
    int blo = hFJ->FindBin(pTlo + 0.01);
    int bhi = hFJ->FindBin(pThi - 0.01);
    double ymax = 0.;
    for(int b = blo; b <= bhi; b++){
        ymax = TMath::Max(ymax, hFJ->GetBinContent(b));
        ymax = TMath::Max(ymax, hRC->GetBinContent(b));
    }

    const double lm = 0.20, rm = 0.05;
    const double fUp = 0.65, fDn = 0.35;
    const double sfUp = 1./fUp, sfDn = 1./fDn;

    TCanvas* c = new TCanvas(Form("c_C%d", ci), "", 700, 700);

    TPad* pUp = new TPad("pUp", "", 0., fDn, 1., 1.);
    pUp->SetLeftMargin(lm); pUp->SetRightMargin(rm);
    pUp->SetTopMargin(0.09); pUp->SetBottomMargin(0.02);
    pUp->SetLogy(); pUp->SetTickx(1); pUp->SetTicky(1);
    pUp->Draw(); pUp->cd();

    hFJ->SetLineColor(col_fastJet);  hFJ->SetLineWidth(2); hFJ->SetMarkerColor(col_fastJet);
    hFJ->SetMarkerStyle(20); hFJ->SetMarkerSize(0.6);
    hRC->SetLineColor(col_randCone); hRC->SetLineWidth(2); hRC->SetMarkerColor(col_randCone);
    hRC->SetMarkerStyle(21); hRC->SetMarkerSize(0.6);

    hFJ->SetStats(0); hFJ->SetTitle("");
    hFJ->GetXaxis()->SetRangeUser(pTlo, pThi);
    hFJ->GetXaxis()->SetLabelSize(0); hFJ->GetXaxis()->SetTitleSize(0);
    hFJ->GetYaxis()->SetTitle("dN/dp_{T} per event [GeV^{-1}]");
    hFJ->GetYaxis()->SetTitleSize(0.030 * sfUp);
    hFJ->GetYaxis()->SetTitleOffset(1.80);
    hFJ->GetYaxis()->SetLabelSize(0.028 * sfUp);
    hFJ->SetMaximum(ymax * 5.);
    hFJ->SetMinimum(1e-8);
    hFJ->Draw("ep");
    hRC->Draw("ep same");

    TLegend* lg = new TLegend(0.38, 0.60, 0.93, 0.88);
    lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.030 * sfUp);
    lg->AddEntry(hFJ, "Mixed-event FastJet anti-k_{T} R=0.4", "lp");
    lg->AddEntry(hRC, "Mixed-event random cone R=0.4", "lp");
    lg->Draw();

    TLatex lat; lat.SetNDC(); lat.SetTextSize(0.030 * sfUp);
    lat.DrawLatex(lm + 0.02, 0.94, Form("PbPb 5.02 TeV MinBias  cent. %s", centLabel[ci]));

    c->cd();
    TPad* pDn = new TPad("pDn", "", 0., 0., 1., fDn);
    pDn->SetLeftMargin(lm); pDn->SetRightMargin(rm);
    pDn->SetTopMargin(0.02); pDn->SetBottomMargin(0.28);
    pDn->SetTickx(1); pDn->SetTicky(1);
    pDn->Draw(); pDn->cd();

    hRat->SetLineColor(col_ratio); hRat->SetLineWidth(2);
    hRat->SetMarkerColor(col_ratio); hRat->SetMarkerStyle(20); hRat->SetMarkerSize(0.6);
    hRat->SetStats(0); hRat->SetTitle("");
    hRat->GetXaxis()->SetRangeUser(pTlo, pThi);
    hRat->GetXaxis()->SetTitle("p_{T} [GeV]");
    hRat->GetXaxis()->SetTitleSize(0.038 * sfDn);
    hRat->GetXaxis()->SetTitleOffset(1.1);
    hRat->GetXaxis()->SetLabelSize(0.032 * sfDn);
    hRat->GetYaxis()->SetTitle("FastJet / RandomCone");
    hRat->GetYaxis()->SetTitleSize(0.030 * sfDn);
    hRat->GetYaxis()->SetTitleOffset(1.10);
    hRat->GetYaxis()->SetLabelSize(0.030 * sfDn);
    hRat->GetYaxis()->SetNdivisions(505);
    hRat->SetMaximum(3.0);
    hRat->SetMinimum(0.);
    hRat->Draw("ep");

    TLine* l1 = new TLine(pTlo, 1., pThi, 1.);
    l1->SetLineStyle(2); l1->SetLineColor(kGray+1); l1->SetLineWidth(2); l1->Draw();

    TString outFile = TString(outDir) + Form("fakeJetAnalysis_MinBias_C%d.pdf", ci);
    c->SaveAs(outFile);
    printf("Saved %s\n", outFile.Data());

    delete c; delete hFJ; delete hRC; delete hRat;
}

void plotPseudoJets_fakeJetAnalysis_MinBias(){

    TFile* f = TFile::Open(inFile);
    if(!f || f->IsZombie()){ printf("ERROR: cannot open %s\n", inFile); return; }

    gSystem->Exec(Form("mkdir -p %s", outDir));

    for(int ci = 0; ci < NCent; ci++) plotOneBin(f, ci);

    f->Close();
    printf("\nDone. Output in %s\n", outDir);
}
