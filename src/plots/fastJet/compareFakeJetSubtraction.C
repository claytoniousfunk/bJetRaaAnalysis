// Compare the effect of using data-driven vs HYDJET fake-jet spectra on the
// PbPb MinBias inclusive jet spectrum.
//
// Per-centrality 2-pad plots:
//   Upper (log-y): raw MinBias, data-subtracted, HYDJET-subtracted (per event)
//   Lower:         ratio (data-subtracted) / (HYDJET-subtracted)
//
// Summary plot: data-sub / HYDJET-sub ratio for all 4 centralities overlaid.

#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TStyle.h"
#include "TROOT.h"

static const double xLo = 20., xHi = 200.;

// Load, clone (with SetDirectory(nullptr)), and return histogram from file.
TH1D* getH(TFile* f, const char* name)
{
    TH1D* h = (TH1D*) f->Get(name);
    if(!h){ printf("WARNING: %s not found in %s\n", name, f->GetName()); return nullptr; }
    TH1D* hc = (TH1D*) h->Clone(Form("%s_clone", name));
    hc->SetDirectory(nullptr);
    return hc;
}

void drawCentComparison(TH1D* hMB, TH1D* hDataSub, TH1D* hHYSub,
                        const char* centLabel, const char* outPath)
{
    const double fUp = 0.65, fDn = 1. - fUp;
    const double sfUp = 1./fUp, sfDn = 1./fDn;
    const double mL = 0.24, mR = 0.05;

    hMB->SetLineColor(kBlack);   hMB->SetLineWidth(2); hMB->SetLineStyle(1);
    hDataSub->SetLineColor(kAzure+1); hDataSub->SetLineWidth(2); hDataSub->SetLineStyle(1);
    hHYSub->SetLineColor(kRed);  hHYSub->SetLineWidth(2); hHYSub->SetLineStyle(2);

    for(TH1D* h : {hMB, hDataSub, hHYSub})
        h->GetXaxis()->SetRangeUser(xLo, xHi);

    double ymax = TMath::Max(hMB->GetMaximum(), TMath::Max(hDataSub->GetMaximum(), hHYSub->GetMaximum()));
    hMB->SetMaximum(ymax * 5.);
    hMB->SetMinimum(1e-9);
    hMB->SetTitle("");
    hMB->GetXaxis()->SetLabelSize(0);
    hMB->GetXaxis()->SetTitleSize(0);
    hMB->GetYaxis()->SetTitle("#frac{1}{N^{evt}} #frac{dN^{jet}}{dp_{T}} [GeV^{-1}]");
    hMB->GetYaxis()->SetTitleSize(0.044 * sfUp);
    hMB->GetYaxis()->SetTitleOffset(1.3);
    hMB->GetYaxis()->SetLabelSize(0.04 * sfUp);

    TH1D* hRatio = (TH1D*) hDataSub->Clone("hRatio_tmp");
    hRatio->SetDirectory(nullptr);
    hRatio->Divide(hHYSub);
    hRatio->SetLineColor(kBlack); hRatio->SetLineWidth(2); hRatio->SetLineStyle(1);
    hRatio->SetMarkerStyle(20);   hRatio->SetMarkerSize(0.7);
    hRatio->GetXaxis()->SetRangeUser(xLo, xHi);
    hRatio->SetTitle("");
    hRatio->GetXaxis()->SetTitle("Jet p_{T} [GeV]");
    hRatio->GetXaxis()->SetTitleSize(0.05 * sfDn);
    hRatio->GetXaxis()->SetTitleOffset(1.0);
    hRatio->GetXaxis()->SetLabelSize(0.04 * sfDn);
    hRatio->GetYaxis()->SetTitle("Data / HYDJET");
    hRatio->GetYaxis()->SetTitleSize(0.045 * sfDn);
    hRatio->GetYaxis()->SetTitleOffset(0.9);
    hRatio->GetYaxis()->SetLabelSize(0.04 * sfDn);
    hRatio->GetYaxis()->SetNdivisions(504);
    hRatio->SetMaximum(2.0);
    hRatio->SetMinimum(0.);

    TCanvas* c = new TCanvas("c_cent", "", 600, 700);

    TPad* pUp = new TPad("pUp", "", 0., fDn, 1., 1.);
    pUp->SetLeftMargin(mL); pUp->SetRightMargin(mR);
    pUp->SetTopMargin(0.08); pUp->SetBottomMargin(0.02);
    pUp->SetLogy(); pUp->SetTickx(1); pUp->SetTicky(1);
    pUp->Draw(); pUp->cd();

    hMB->Draw("hist");
    hDataSub->Draw("hist same");
    hHYSub->Draw("hist same");

    TLegend* lg = new TLegend(0.52, 0.62, 0.93, 0.88);
    lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.038 * sfUp);
    lg->AddEntry(hMB,      "PbPb MinBias (raw)", "l");
    lg->AddEntry(hDataSub, "Data fake-jet sub.", "l");
    lg->AddEntry(hHYSub,   "HYDJET fake-jet sub.", "l");
    lg->Draw();

    TLatex lat;
    lat.SetNDC(); lat.SetTextSize(0.044 * sfUp);
    lat.DrawLatex(0.28, 0.20, centLabel);
    lat.SetTextSize(0.033 * sfUp);
    lat.DrawLatex(0.28, 0.13, "Mixed-event anti-k_{T} R=0.4, CS PF cands");

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
    delete c;
    delete hRatio;
}

void compareFakeJetSubtraction()
{
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);

    const char* fMBPath =
        "/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/scanningOutput/"
        "PbPb/latest/PbPb_MinBias_mu12_pTmu-15to999_tight_WDecayFilter_2026-4-8.root";
    const char* fDataFakePath =
        "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/fakeJets/fakeJets.root";
    const char* fHYFakePath =
        "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/fakeJets/fakeJets_HYDJET.root";
    const char* outDir =
        "/home/clayton/Analysis/code/bJetRaaAnalysis/figures/fakeJetComparison/";

    gSystem->Exec(Form("mkdir -p %s", outDir));

    TFile* fMB     = TFile::Open(fMBPath);
    TFile* fData   = TFile::Open(fDataFakePath);
    TFile* fHY     = TFile::Open(fHYFakePath);
    if(!fMB   || fMB->IsZombie()){ printf("ERROR: cannot open MinBias file\n");    return; }
    if(!fData || fData->IsZombie()){ printf("ERROR: cannot open data fake file\n"); return; }
    if(!fHY   || fHY->IsZombie()){ printf("ERROR: cannot open HYDJET fake file\n"); return; }

    struct Cent { int idx; const char* label; Int_t col; };
    Cent cents[] = {
        {0, "0--80%",   kBlack},
        {1, "0--10%",   TColor::GetColor("#D55E00")},
        {2, "10--30%",  TColor::GetColor("#0072B2")},
        {3, "30--50%",  TColor::GetColor("#009E73")},
        {4, "50--80%",  TColor::GetColor("#E69F00")},
    };

    // Histograms for the summary plot (C1-C4 only)
    TH1D* hRatios[4] = {};

    for(int ci = 0; ci < 5; ci++){
        int idx = cents[ci].idx;

        TH1D* hMB  = getH(fMB, Form("h_inclRecoJetPt_C%d", idx));
        TH1D* hvz  = getH(fMB, Form("h_vz_C%d", idx));
        TH1D* hFD  = getH(fData, Form("h_fakeJets_C%d", idx));
        TH1D* hFHY = getH(fHY,   Form("h_fakeJets_C%d", idx));
        if(!hMB || !hvz || !hFD || !hFHY){
            printf("Skipping C%d due to missing histograms\n", idx);
            continue;
        }

        double N_evt = hvz->Integral();
        hMB->Scale(1. / N_evt);

        printf("C%d: N_MB_events=%.0f  data_fake_integral=%.6f  HY_fake_integral=%.6f\n",
               idx, N_evt, hFD->Integral(), hFHY->Integral());

        TH1D* hDataSub = (TH1D*) hMB->Clone(Form("hDataSub_C%d", idx));
        hDataSub->SetDirectory(nullptr);
        hDataSub->Add(hFD, -1.);

        TH1D* hHYSub = (TH1D*) hMB->Clone(Form("hHYSub_C%d", idx));
        hHYSub->SetDirectory(nullptr);
        hHYSub->Add(hFHY, -1.);

        drawCentComparison(hMB, hDataSub, hHYSub,
                           cents[ci].label,
                           Form("%sfakeJetSub_C%d.pdf", outDir, idx));

        // Keep ratio for summary (C1-C4 only)
        if(idx >= 1){
            int ri = idx - 1;
            hRatios[ri] = (TH1D*) hDataSub->Clone(Form("hRatioSum_C%d", idx));
            hRatios[ri]->SetDirectory(nullptr);
            hRatios[ri]->Divide(hHYSub);
            hRatios[ri]->GetXaxis()->SetRangeUser(xLo, xHi);
            hRatios[ri]->SetLineColor(cents[ci].col);
            hRatios[ri]->SetLineWidth(2);
        }

        delete hMB; delete hvz; delete hFD; delete hFHY;
        delete hDataSub; delete hHYSub;
    }

    // --- Summary: (data-sub)/(HYDJET-sub) ratio for all centralities ---
    TCanvas* cSum = new TCanvas("cSum", "", 900, 650);
    cSum->SetLeftMargin(0.16); cSum->SetRightMargin(0.05);
    cSum->SetBottomMargin(0.14); cSum->SetTopMargin(0.08);
    cSum->SetTickx(1); cSum->SetTicky(1);

    bool first = true;
    for(int ci = 0; ci < 4; ci++){
        if(!hRatios[ci]) continue;
        hRatios[ci]->SetTitle("");
        hRatios[ci]->GetXaxis()->SetTitle("Jet p_{T} [GeV]");
        hRatios[ci]->GetYaxis()->SetTitle("(MinBias #minus data fake) / (MinBias #minus HYDJET fake)");
        hRatios[ci]->GetYaxis()->SetRangeUser(0., 2.5);
        hRatios[ci]->GetXaxis()->SetTitleSize(0.05); hRatios[ci]->GetXaxis()->SetLabelSize(0.045);
        hRatios[ci]->GetYaxis()->SetTitleSize(0.042); hRatios[ci]->GetYaxis()->SetLabelSize(0.042);
        hRatios[ci]->GetYaxis()->SetTitleOffset(1.8);
        if(first){ hRatios[ci]->Draw("hist"); first = false; }
        else        hRatios[ci]->Draw("hist same");
    }

    TLine* lineOne = new TLine(xLo, 1., xHi, 1.);
    lineOne->SetLineStyle(2); lineOne->SetLineColor(kGray+1); lineOne->SetLineWidth(2);
    lineOne->Draw();

    TLegend* lgSum = new TLegend(0.20, 0.68, 0.55, 0.90);
    lgSum->SetBorderSize(0); lgSum->SetTextSize(0.042);
    for(int ci = 0; ci < 4; ci++){
        if(hRatios[ci]) lgSum->AddEntry(hRatios[ci], cents[ci].label, "l");
    }
    lgSum->Draw();

    TLatex latSum;
    latSum.SetNDC(); latSum.SetTextSize(0.038);
    latSum.DrawLatex(0.18, 0.92, "PbPb MinBias (5.02 TeV)  fake jets: mixed-event anti-k_{T} R=0.4 CS PF cands");

    cSum->SaveAs(Form("%sfakeJetSub_ratio_summary.pdf", outDir));
    delete cSum;

    fMB->Close(); fData->Close(); fHY->Close();
    printf("\nPlots saved to %s\n", outDir);
}
