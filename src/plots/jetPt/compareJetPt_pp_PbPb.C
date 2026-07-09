// Compare inclusive reco jet pT spectra: pp ZeroBias vs PbPb MinBias.
// All 16 ultraFine centrality bins shown individually (0-5%, 5-10%, ..., 75-80%).
// Each spectrum pinned to 1 in the 200-300 GeV bin (shape comparison).
// Upper pad: all spectra (log-y).  Lower pad: PbPb/pp ratio (log-y).
// Color gradient: deep blue (most central 0-5%) → dark red (most peripheral 75-80%).
//
// Usage: root -l -q 'compareJetPt_pp_PbPb.C'
// Run from: src/plots/jetPt/

#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TStyle.h"
#include "TColor.h"

void compareJetPt_pp_PbPb()
{
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);

    const char* fppPath =
        "../../../rootFiles/scanningOuput/pp/"
        "pp_MinBias_mu12_pTmu-15to999_tight_deltaR-40_jetTrkMaxFilter_WDecayFilter_2026-7-6.root";
    const char* fPbPbPath =
        "../../../rootFiles/scanningOuput/PbPb/"
        "PbPb_MinBias_Part1_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_2026-7-6_ultraFineCentBins.root";
    const char* outPath =
        "../../../figures/jetPt/compareJetPt_pp_PbPb.pdf";

    TFile* fpp   = TFile::Open(fppPath);
    TFile* fPbPb = TFile::Open(fPbPbPath);
    if(!fpp   || fpp  ->IsZombie()){ printf("ERROR: cannot open pp file\n");   return; }
    if(!fPbPb || fPbPb->IsZombie()){ printf("ERROR: cannot open PbPb file\n"); return; }

    // --- pp ---
    TH1D* hJet_pp = (TH1D*) fpp->Get("h_inclRecoJetPt"); hJet_pp->SetDirectory(nullptr);
    printf("pp ZeroBias:  N_evt = %.0f\n", ((TH1D*)fpp->Get("h_vz"))->Integral());

    // --- PbPb: one histogram per ultraFine bin (C1=0-5% ... C16=75-80%) ---
    const int NCent = 16;
    const char* centLabel[NCent] = {
        "0-5%",  "5-10%",  "10-15%", "15-20%",
        "20-25%","25-30%", "30-35%", "35-40%",
        "40-45%","45-50%", "50-55%", "55-60%",
        "60-65%","65-70%", "70-75%", "75-80%"
    };

    TH1D* hJet[NCent];
    for(int ic = 0; ic < NCent; ic++){
        int ci = ic + 1;  // C1 through C16
        hJet[ic] = (TH1D*) fPbPb->Get(Form("h_inclRecoJetPt_C%d", ci));
        hJet[ic]->SetDirectory(nullptr);
        printf("PbPb %s:  N_evt = %.0f\n", centLabel[ic],
               ((TH1D*)fPbPb->Get(Form("h_vz_C%d",ci)))->Integral());
    }

    fpp->Close(); fPbPb->Close();

    // --- Rebin to variable-width bins; divide by width to get dN/dpT [GeV^-1] ---
    // Finer bins at low pT, coarser at high pT where stats are limited.
    const int NBins = 12;
    double edges[NBins + 1] = {30, 35, 40, 45, 50, 60, 70, 80, 100, 120, 150, 200, 300 };

    auto doRebin = [&](TH1D* h, const char* name) -> TH1D* {
        TH1D* hr = (TH1D*) h->Rebin(NBins, name, edges);
        hr->SetDirectory(nullptr);
        for(int b = 1; b <= hr->GetNbinsX(); b++){
            double w = hr->GetBinWidth(b);
            hr->SetBinContent(b, hr->GetBinContent(b) / w);
            hr->SetBinError(b,   hr->GetBinError(b)   / w);
        }
        return hr;
    };

    TH1D* hJet_pp_r = doRebin(hJet_pp, "hJet_pp_r");
    TH1D* hJet_r[NCent];
    for(int ic = 0; ic < NCent; ic++)
        hJet_r[ic] = doRebin(hJet[ic], Form("hJet_r_%d", ic));

    // --- Pin each spectrum to 1 in the 200-300 GeV bin ---
    auto normalize = [](TH1D* h){
      int S = h->Integral();  // center of the 200-300 bin
      if(S > 0.) h->Scale(1. / S);
    };
    normalize(hJet_pp_r);
    for(int ic = 0; ic < NCent; ic++) normalize(hJet_r[ic]);

    // --- Ratios (PbPb / pp) from pinned rebinned histograms ---
    TH1D* hRatio[NCent];
    for(int ic = 0; ic < NCent; ic++){
        hRatio[ic] = (TH1D*) hJet_r[ic]->Clone(Form("hRatio_%d", ic));
        hRatio[ic]->SetDirectory(nullptr);
        hRatio[ic]->Divide(hJet_pp_r);
    }

    // --- Color gradient: Plasma colormap (colorblind-friendly, perceptually uniform) ---
    // Built via CreateGradientColorTable to avoid named-constant version dependencies.
    // Samples 5%-85% of the full Plasma range: deep blue-purple (central) → orange (peripheral).
    // Avoids the near-white yellow at the bright end of the full Plasma range.
    {
        const int nStop = 5;
        double stops[nStop] = {0.000, 0.250, 0.500, 0.750, 1.000};
        double rr[nStop]    = {0.051, 0.490, 0.796, 0.974, 0.988};
        double gg[nStop]    = {0.031, 0.010, 0.278, 0.541, 0.652};
        double bb[nStop]    = {0.529, 0.530, 0.471, 0.050, 0.020};
        TColor::CreateGradientColorTable(nStop, stops, rr, gg, bb, 256);
    }
    int nPal = TColor::GetNumberOfColors();
    int cols[NCent];
    for(int ic = 0; ic < NCent; ic++){
        float t = 0.05f + (float(ic) / float(NCent - 1)) * 0.80f;
        cols[ic] = TColor::GetColorPalette(int(t * (nPal - 1)));
    }

    // Unique marker per centrality bin (16 distinct ROOT marker styles)
    const int markerStyles[NCent] = {
        20, 21, 22, 23,   // filled: circle, square, tri-up, tri-down
        24, 25, 26, 32,   // open:   circle, square, tri-up, tri-down
        29, 30, 33, 27,   // filled star, open star, filled diamond, open diamond
        34, 28, 43, 45    // filled cross+, open cross+, filled cross×, open cross×
    };

    // --- Style ---
    hJet_pp_r->SetLineColor(kBlack); hJet_pp_r->SetLineWidth(2);
    hJet_pp_r->SetMarkerColor(kBlack); hJet_pp_r->SetMarkerStyle(29); // filled star for pp
    hJet_pp_r->SetMarkerSize(0.8);
    hJet_pp_r->SetLineStyle(2); hJet_pp_r->SetStats(0); hJet_pp_r->SetTitle("");

    for(int ic = 0; ic < NCent; ic++){
        hJet_r[ic]->SetLineColor(cols[ic]); hJet_r[ic]->SetLineWidth(1);
        hJet_r[ic]->SetMarkerColor(cols[ic]); hJet_r[ic]->SetMarkerStyle(markerStyles[ic]);
        hJet_r[ic]->SetMarkerSize(0.6);
        hJet_r[ic]->SetStats(0); hJet_r[ic]->SetTitle("");

        hRatio[ic]->SetLineColor(cols[ic]); hRatio[ic]->SetLineWidth(1);
        hRatio[ic]->SetMarkerColor(cols[ic]); hRatio[ic]->SetMarkerStyle(markerStyles[ic]);
        hRatio[ic]->SetMarkerSize(0.6);
        hRatio[ic]->SetStats(0); hRatio[ic]->SetTitle("");
    }

    const double xLo = 30., xHi = 300.;
    hJet_pp_r->GetXaxis()->SetRangeUser(xLo, xHi);
    for(int ic = 0; ic < NCent; ic++){
        hJet_r[ic]->GetXaxis()->SetRangeUser(xLo, xHi);
        hRatio[ic]->GetXaxis()->SetRangeUser(xLo, xHi);
    }

    // --- Canvas layout ---
    const double fUp = 0.65, fDn = 1. - fUp;
    const double sfUp = 1. / fUp, sfDn = 1. / fDn;
    const double mL = 0.20, mR = 0.05;

    // Upper pad y-range: pin point = 1, so min just below 1 and max captures low-pT peak
    double ymax = 0.;
    for(int ic = 0; ic < NCent; ic++) ymax = TMath::Max(ymax, hJet_r[ic]->GetMaximum());
    // hJet_pp_r->SetMaximum(ymax * 4.);
    // hJet_pp_r->SetMinimum(0.3);
    hJet_pp_r->GetXaxis()->SetLabelSize(0); hJet_pp_r->GetXaxis()->SetTitleSize(0);
    hJet_pp_r->GetYaxis()->SetTitle("1/N dN/dp_{T}");
    hJet_pp_r->GetYaxis()->SetTitleSize(0.030 * sfUp);
    hJet_pp_r->GetYaxis()->SetTitleOffset(1.80);
    hJet_pp_r->GetYaxis()->SetLabelSize(0.034 * sfUp);

    // Lower pad y-range
    hRatio[0]->GetXaxis()->SetTitle("Jet p_{T} [GeV]");
    hRatio[0]->GetXaxis()->SetTitleSize(0.042 * sfDn);
    hRatio[0]->GetXaxis()->SetTitleOffset(1.1);
    hRatio[0]->GetXaxis()->SetLabelSize(0.036 * sfDn);
    hRatio[0]->GetYaxis()->SetTitle("PbPb / pp");
    hRatio[0]->GetYaxis()->SetTitleSize(0.032 * sfDn);
    hRatio[0]->GetYaxis()->SetTitleOffset(1.10);
    hRatio[0]->GetYaxis()->SetLabelSize(0.036 * sfDn);
    hRatio[0]->GetYaxis()->SetNdivisions(505);
    hRatio[0]->SetMaximum(2);
    hRatio[0]->SetMinimum(0.0);

    TCanvas* c = new TCanvas("c", "", 600, 700);

    TPad* pUp = new TPad("pUp", "", 0., fDn, 1., 1.);
    pUp->SetLeftMargin(mL); pUp->SetRightMargin(mR);
    pUp->SetTopMargin(0.08); pUp->SetBottomMargin(0.02);
    pUp->SetLogy(); pUp->SetTickx(1); pUp->SetTicky(1);
    pUp->Draw(); pUp->cd();

    hJet_pp_r->Draw("ep");
    for(int ic = 0; ic < NCent; ic++) hJet_r[ic]->Draw("ep same");

    // Two-column legend pinned to top-right quadrant
    TLegend* lg = new TLegend(0.48, 0.40, 0.93, 0.92);
    lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.028 * sfUp);
    lg->SetNColumns(2);
    lg->AddEntry(hJet_pp_r, "pp ZeroBias", "lp");
    lg->AddEntry((TObject*)nullptr, "", "");  // blank to fill 2nd column
    for(int ic = 0; ic < NCent; ic++)
        lg->AddEntry(hJet_r[ic], Form("PbPb %s", centLabel[ic]), "lp");
    lg->Draw();

    TLatex lat; lat.SetNDC(); lat.SetTextSize(0.035 * sfUp);
    lat.DrawLatex(0.22, 0.93, "5.02 TeV  anti-k_{T} R=0.4");

    c->cd();
    TPad* pDn = new TPad("pDn", "", 0., 0., 1., fDn);
    pDn->SetLeftMargin(mL); pDn->SetRightMargin(mR);
    pDn->SetTopMargin(0.02); pDn->SetBottomMargin(0.30);
    pDn->SetTickx(1); pDn->SetTicky(1);
    pDn->Draw(); pDn->cd();

    hRatio[0]->Draw("ep");
    for(int ic = 1; ic < NCent; ic++) hRatio[ic]->Draw("ep same");
    TLine* line = new TLine(xLo, 1., xHi, 1.);  // ratio = 1 at 200-300 GeV by construction
    line->SetLineStyle(2); line->SetLineColor(kGray+1); line->SetLineWidth(1);
    line->Draw();

    c->SaveAs(outPath);
    printf("Saved to %s\n", outPath);
    delete c;

    // =========================================================================
    // Ratio-to-prior: for each adjacent pair of centrality bins compute
    //   hJet_r[ic] / hJet_r[ic-1]  (more peripheral / less peripheral)
    // Shows where the shape changes most as centrality changes by 5%.
    // =========================================================================
    // Coarser binning for the ratio plot — rebin numerator and denominator
    // independently before dividing so the ratio errors are computed correctly.
    const double rpEdges[] = { 20, 30, 40, 50, 60, 80, 100, 120, 150, 200, 300 };
    const int    nRPEdges  = sizeof(rpEdges) / sizeof(rpEdges[0]);
    const int    nRPBins   = nRPEdges - 1;

    const int NRP = NCent - 1;  // 15 ratios
    TH1D* hRP[NRP];
    for(int ic = 0; ic < NRP; ic++){
        // rebin pinned (but dN/dpT) histograms — bin-width factors cancel in ratio
        TH1D* hNum = (TH1D*) hJet_r[ic + 1]->Rebin(nRPBins, Form("hNum_%d", ic), rpEdges);
        TH1D* hDen = (TH1D*) hJet_r[ic]    ->Rebin(nRPBins, Form("hDen_%d", ic), rpEdges);
        hNum->SetDirectory(nullptr); hDen->SetDirectory(nullptr);
        hRP[ic] = hNum;
        hRP[ic]->Divide(hDen);
        delete hDen;
        // reuse the color and marker of the more-peripheral (upper) bin
        hRP[ic]->SetLineColor(cols[ic + 1]); hRP[ic]->SetLineWidth(1);
        hRP[ic]->SetMarkerColor(cols[ic + 1]); hRP[ic]->SetMarkerStyle(markerStyles[ic + 1]);
        hRP[ic]->SetMarkerSize(0.6);
        hRP[ic]->SetStats(0); hRP[ic]->SetTitle("");
        hRP[ic]->GetXaxis()->SetRangeUser(xLo, xHi);
    }

    // y-range: expect most ratios near 1; find actual extremes in x-range
    double rpMax = 0., rpMin = 1e9;
    for(int ic = 0; ic < NRP; ic++){
        for(int b = hRP[ic]->FindBin(xLo); b <= hRP[ic]->FindBin(xHi - 1); b++){
            double v = hRP[ic]->GetBinContent(b);
            if(v > 0.){ rpMax = TMath::Max(rpMax, v); rpMin = TMath::Min(rpMin, v); }
        }
    }

    hRP[0]->SetMaximum(rpMax * 1.4);
    hRP[0]->SetMinimum(0.);
    hRP[0]->GetXaxis()->SetTitle("Jet p_{T} [GeV]");
    hRP[0]->GetXaxis()->SetTitleSize(0.05);
    hRP[0]->GetXaxis()->SetTitleOffset(1.1);
    hRP[0]->GetXaxis()->SetLabelSize(0.04);
    hRP[0]->GetYaxis()->SetTitle("Cent_{i+1} / Cent_{i}  (adjacent 5% bins)");
    hRP[0]->GetYaxis()->SetTitleSize(0.042);
    hRP[0]->GetYaxis()->SetTitleOffset(1.6);
    hRP[0]->GetYaxis()->SetLabelSize(0.04);
    hRP[0]->GetYaxis()->SetNdivisions(505);

    // Wide canvas: right margin houses the legend outside the frame
    TCanvas* c2 = new TCanvas("c2", "", 950, 520);
    c2->SetLeftMargin(0.13); c2->SetRightMargin(0.38);
    c2->SetBottomMargin(0.15); c2->SetTopMargin(0.08);
    c2->SetTickx(1); c2->SetTicky(1);
    c2->cd();

    hRP[0]->Draw("ep");
    for(int ic = 1; ic < NRP; ic++) hRP[ic]->Draw("ep same");

    TLine* line2 = new TLine(xLo, 1., xHi, 1.);
    line2->SetLineStyle(2); line2->SetLineColor(kGray+1); line2->SetLineWidth(1);
    line2->Draw();

    // Legend placed in the right margin, outside the plot frame
    TLegend* lg2 = new TLegend(0.638, 0.08, 0.995, 0.92);
    lg2->SetBorderSize(0); lg2->SetFillStyle(0); lg2->SetTextSize(0.026);
    lg2->SetNColumns(2);
    for(int ic = 0; ic < NRP; ic++)
        lg2->AddEntry(hRP[ic], Form("%s / %s", centLabel[ic+1], centLabel[ic]), "lp");
    lg2->Draw();

    TLatex lat2; lat2.SetNDC(); lat2.SetTextSize(0.036);
    lat2.DrawLatex(0.14, 0.945, "5.02 TeV PbPb MinBias  anti-k_{T} R=0.4");

    const char* outPath2 = "../../../figures/jetPt/compareJetPt_pp_PbPb_ratioPrior.pdf";
    c2->SaveAs(outPath2);
    printf("Saved to %s\n", outPath2);
    delete c2;

    // =========================================================================
    // One PDF per centrality ratio pair, saved to a dedicated folder.
    // Shared y-range across all plots so they are directly comparable.
    // =========================================================================
    const char* indDir = "../../../figures/jetPt/ratioPrior_individual";
    TCanvas* cInd = new TCanvas("cInd", "", 700, 520);
    cInd->SetLeftMargin(0.14); cInd->SetRightMargin(0.05);
    cInd->SetBottomMargin(0.15); cInd->SetTopMargin(0.08);
    cInd->SetTickx(1); cInd->SetTicky(1);

    TLine* lineInd = new TLine(xLo, 1., xHi, 1.);
    lineInd->SetLineStyle(2); lineInd->SetLineColor(kGray+1); lineInd->SetLineWidth(1);

    for(int ic = 0; ic < NRP; ic++){
        cInd->cd(); cInd->Clear();

        hRP[ic]->SetMaximum(rpMax * 1.4);
        hRP[ic]->SetMinimum(0.);
        hRP[ic]->GetXaxis()->SetTitle("Jet p_{T} [GeV]");
        hRP[ic]->GetXaxis()->SetTitleSize(0.05);
        hRP[ic]->GetXaxis()->SetTitleOffset(1.1);
        hRP[ic]->GetXaxis()->SetLabelSize(0.04);
        hRP[ic]->GetYaxis()->SetTitle("Cent_{i+1} / Cent_{i}");
        hRP[ic]->GetYaxis()->SetTitleSize(0.048);
        hRP[ic]->GetYaxis()->SetTitleOffset(1.25);
        hRP[ic]->GetYaxis()->SetLabelSize(0.04);
        hRP[ic]->GetYaxis()->SetNdivisions(505);
        hRP[ic]->Draw("ep");
        lineInd->Draw();

        TLegend* lgInd = new TLegend(0.18, 0.78, 0.65, 0.90);
        lgInd->SetBorderSize(0); lgInd->SetFillStyle(0); lgInd->SetTextSize(0.038);
        lgInd->AddEntry(hRP[ic], Form("PbPb %s / %s", centLabel[ic+1], centLabel[ic]), "lp");
        lgInd->Draw();

        TLatex latInd; latInd.SetNDC(); latInd.SetTextSize(0.034);
        latInd.DrawLatex(0.15, 0.945, "5.02 TeV PbPb MinBias  anti-k_{T} R=0.4");

        // filename: cent_05to10_over_00to05.pdf  etc.
        TString fname = Form("%s/ratioPrior_%s_over_%s.pdf",
                             indDir, centLabel[ic+1], centLabel[ic]);
        fname.ReplaceAll("%", "");  // remove % signs for clean filenames
        cInd->SaveAs(fname);
    }
    delete cInd;
    delete lineInd;
    printf("Saved %d individual ratio plots to %s/\n", NRP, indDir);
}
