// Compare inclusive reco jet pT spectrum from the pp ZeroBias miniAOD test scan
// (single file) against the reference pp MinBias AOD scan.
// Both normalized per inclusive event (1/N_evt * dN/dpT).
// Upper pad: both spectra (log-y).  Lower pad: ratio miniAOD / reference.

#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TStyle.h"

void compareJetPt_ZeroBias_miniAOD_test()
{
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);

    const char* fTestPath =
        "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/pp/"
        "pp_MinBias_mu12_pTmu-15to999_tight_deltaR-40_jetTrkMaxFilter_WDecayFilter_2026-7-6.root";
    const char* fRefPath =
        "/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/"
        "scanningOutput/pp/latest/"
        "pp_MinBias_mu12_pTmu-15to999_tight_deltaR-40_jetTrkMaxFilter_WDecayFilter_2026-4-7.root";
    const char* outPath =
        "/home/clayton/Analysis/code/bJetRaaAnalysis/figures/jetPt/compare_ZeroBias_miniAOD/"
        "compareJetPt_ZeroBias_miniAOD_test.pdf";

    TFile* fTest = TFile::Open(fTestPath);
    TFile* fRef  = TFile::Open(fRefPath);
    if(!fTest || fTest->IsZombie()){ printf("ERROR: cannot open test file\n");      return; }
    if(!fRef  || fRef ->IsZombie()){ printf("ERROR: cannot open reference file\n"); return; }

    TH1D* hJet_test = (TH1D*) fTest->Get("h_inclRecoJetPt"); hJet_test->SetDirectory(nullptr);
    TH1D* hvz_test  = (TH1D*) fTest->Get("h_vz");            hvz_test ->SetDirectory(nullptr);
    TH1D* hJet_ref  = (TH1D*) fRef ->Get("h_inclRecoJetPt"); hJet_ref ->SetDirectory(nullptr);
    TH1D* hvz_ref   = (TH1D*) fRef ->Get("h_vz");            hvz_ref  ->SetDirectory(nullptr);
    fTest->Close(); fRef->Close();

    double N_test = hvz_test->Integral();
    double N_ref  = hvz_ref ->Integral();
    printf("miniAOD test: N_events = %.0f\n", N_test);
    printf("Reference:    N_events = %.0f\n", N_ref);

    hJet_test->Scale(1. / N_test);
    hJet_ref ->Scale(1. / N_ref);

    // Ratio — built by bin-center lookup since the two histograms have
    // different ranges (test: 0–500, ref: 20–500), same 5 GeV bin width.
    TH1D* hRatio = (TH1D*) hJet_ref->Clone("hRatio");
    hRatio->SetDirectory(nullptr);
    hRatio->Reset();
    for(int b = 1; b <= hJet_ref->GetNbinsX(); b++){
        double pt  = hJet_ref->GetXaxis()->GetBinCenter(b);
        double ref = hJet_ref->GetBinContent(b);
	double err_ref = hJet_ref->GetBinError(b);
        double tst = hJet_test->GetBinContent(hJet_test->FindBin(pt));
	double err_tst = hJet_test->GetBinError(hJet_test->FindBin(pt));
        hRatio->SetBinContent(b, ref > 0 ? tst / ref : 0.);
	hRatio->SetBinError(b, ref > 0 ? (tst / ref)*sqrt(pow(err_tst/tst,2) + pow(err_ref/ref,2)) : 0.);
    }

    // Style
    hJet_test->SetLineColor(kAzure+1); hJet_test->SetLineWidth(2); hJet_test->SetLineStyle(1);
    hJet_ref ->SetLineColor(kBlack);   hJet_ref ->SetLineWidth(2); hJet_ref ->SetLineStyle(2);
    hRatio   ->SetLineColor(kBlack);   hRatio   ->SetLineWidth(2);
    hRatio->SetMarkerStyle(20); hRatio->SetMarkerSize(0.7);

    const double xLo = 20., xHi = 200.;
    for(TH1D* h : {hJet_test, hJet_ref, hRatio})
        h->GetXaxis()->SetRangeUser(xLo, xHi);

    // Upper pad axes
    const double fUp = 0.65, fDn = 1. - fUp;
    const double sfUp = 1./fUp, sfDn = 1./fDn;
    const double mL = 0.27, mR = 0.05;

    double ymax = TMath::Max(hJet_test->GetMaximum(), hJet_ref->GetMaximum());
    hJet_test->SetMaximum(ymax * 5.);
    hJet_test->SetMinimum(1e-10);
    hJet_test->SetTitle("");
    hJet_test->GetXaxis()->SetLabelSize(0);
    hJet_test->GetXaxis()->SetTitleSize(0);
    hJet_test->GetYaxis()->SetTitle("#frac{1}{N^{evt}} #frac{dN^{jet}}{dp_{T}} [GeV^{-1}]");
    hJet_test->GetYaxis()->SetTitleSize(0.044 * sfUp);
    hJet_test->GetYaxis()->SetTitleOffset(1.4);
    hJet_test->GetYaxis()->SetLabelSize(0.040 * sfUp);

    // Lower pad axes
    hRatio->SetTitle("");
    hRatio->GetXaxis()->SetTitle("Jet p_{T} [GeV]");
    hRatio->GetXaxis()->SetTitleSize(0.05 * sfDn);
    hRatio->GetXaxis()->SetTitleOffset(1.0);
    hRatio->GetXaxis()->SetLabelSize(0.04 * sfDn);
    hRatio->GetYaxis()->SetTitle("miniAOD / AOD ref.");
    hRatio->GetYaxis()->SetTitleSize(0.035 * sfDn);
    hRatio->GetYaxis()->SetTitleOffset(0.9);
    hRatio->GetYaxis()->SetLabelSize(0.04 * sfDn);
    hRatio->GetYaxis()->SetNdivisions(504);
    hRatio->SetMaximum(1.5);
    hRatio->SetMinimum(0.5);

    TCanvas* c = new TCanvas("c", "", 600, 700);

    TPad* pUp = new TPad("pUp", "", 0., fDn, 1., 1.);
    pUp->SetLeftMargin(mL); pUp->SetRightMargin(mR);
    pUp->SetTopMargin(0.08); pUp->SetBottomMargin(0.02);
    pUp->SetLogy(); pUp->SetTickx(1); pUp->SetTicky(1);
    pUp->Draw(); pUp->cd();

    hJet_test->Draw("hist");
    hJet_ref ->Draw("hist same");

    TLegend* lg = new TLegend(0.38, 0.62, 0.93, 0.84);
    lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.05);
    lg->AddEntry(hJet_test, "pp ZeroBias miniAOD", "l");
    lg->AddEntry(hJet_ref,  "pp ZeroBias AOD",          "l");
    lg->Draw();

    TLatex lat; lat.SetNDC();
    lat.SetTextSize(0.05);
    lat.DrawLatex(0.22, 0.20, "pp 5.02 TeV  anti-k_{T} R=0.4");
    lat.SetTextSize(0.04);
    lat.DrawLatex(0.22, 0.13, Form("miniAOD N_{evt} = %.0f,  AOD N_{evt} = %.0f", N_test, N_ref));

    c->cd();
    TPad* pDn = new TPad("pDn", "", 0., 0., 1., fDn);
    pDn->SetLeftMargin(mL); pDn->SetRightMargin(mR);
    pDn->SetTopMargin(0.02); pDn->SetBottomMargin(0.30);
    pDn->SetTickx(1); pDn->SetTicky(1);
    pDn->Draw(); pDn->cd();

    hRatio->Draw("ep");
    TLine* line = new TLine(xLo, 1., xHi, 1.);
    line->SetLineStyle(2); line->SetLineColor(kGray+1); line->SetLineWidth(1);
    line->Draw();

    c->SaveAs(outPath);
    printf("Saved to %s\n", outPath);
    delete c;
}
