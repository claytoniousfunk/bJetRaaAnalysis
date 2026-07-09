#include "../../headers/functions/divideByBinwidth.h"

// MinBias-only RAA proxy.
// pp  : ZeroBias scan (2026-7-6), no stitching.
// PbPb: ultraFine cent scan (2026-7-6, C1-C16 = 5% bins), merged to nominal.
// No fake-jet subtraction, no unfolding.
// Ratio = (PbPb / N_evt^PbPb) / (T_AA * pp / N_evt^pp).

void calculateRAA_minBiasOnly()
{
    // ------------------------------------------------------------------
    // Input files
    // ------------------------------------------------------------------
    TFile *f_pp = TFile::Open(
        "../../rootFiles/scanningOuput/pp/"
        "pp_MinBias_mu12_pTmu-15to999_tight_deltaR-40_jetTrkMaxFilter_WDecayFilter_2026-7-6.root");
    TFile *f_PbPb = TFile::Open(
        "../../rootFiles/scanningOuput/PbPb/"
        "PbPb_MinBias_Part1_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_"
        "2026-7-6_ultraFineCentBins.root");

    if(!f_pp   || f_pp->IsZombie())  { printf("ERROR: cannot open pp file\n");   return; }
    if(!f_PbPb || f_PbPb->IsZombie()){ printf("ERROR: cannot open PbPb file\n"); return; }

    gSystem->Exec("mkdir -p ../../figures/jetPt/RAA/");

    // ------------------------------------------------------------------
    // pp spectra + event count
    // ------------------------------------------------------------------
    TH1D *h_pp_raw;
    TH1D *h_vz_pp;
    f_pp->GetObject("h_inclRecoJetPt", h_pp_raw);
    f_pp->GetObject("h_vz",            h_vz_pp);
    h_pp_raw->SetDirectory(nullptr);
    h_vz_pp ->SetDirectory(nullptr);

    double N_evt_pp = h_vz_pp->Integral();
    printf("N_evt_pp = %.0f\n", N_evt_pp);

    // ------------------------------------------------------------------
    // ultraFine PbPb spectra -- merge 5% bins into nominal bins
    //
    // ultraFine indices (C1-C16):
    //   C1  = 0-5%    C2  = 5-10%
    //   C3  = 10-15%  C4  = 15-20%  C5  = 20-25%  C6  = 25-30%
    //   C7  = 30-35%  C8  = 35-40%  C9  = 40-45%  C10 = 45-50%
    //   C11 = 50-55%  C12 = 55-60%  C13 = 60-65%  C14 = 65-70%
    //   C15 = 70-75%  C16 = 75-80%
    //
    // Nominal groups (C1-C4):
    //   Nom-C1 (0-10%)  : uF C1,C2
    //   Nom-C2 (10-30%) : uF C3-C6
    //   Nom-C3 (30-50%) : uF C7-C10
    //   Nom-C4 (50-80%) : uF C11-C16
    // ------------------------------------------------------------------

    // ultraFine groups per nominal bin (1-indexed uF bin numbers)
    const int nNom = 4;
    const char *nomLabel[nNom]  = {"0-10%", "10-30%", "30-50%", "50-80%"};
    const int   uF_start[nNom] = {1,  3,  7, 11};
    const int   uF_end  [nNom] = {2,  6, 10, 16};

    // Standard Glauber T_AA values for PbPb 5.02 TeV (mb^{-1})
    const double T_AA[nNom] = {23.05, 11.60, 3.95, 0.7452};

    TH1D *h_PbPb[nNom];
    double N_evt_PbPb[nNom];

    for(int ic = 0; ic < nNom; ic++){
        h_PbPb[ic]     = nullptr;
        N_evt_PbPb[ic] = 0.;

        for(int iuf = uF_start[ic]; iuf <= uF_end[ic]; iuf++){
            TH1D *h_jet_tmp;
            TH1D *h_vz_tmp;
            f_PbPb->GetObject(Form("h_inclRecoJetPt_C%d", iuf), h_jet_tmp);
            f_PbPb->GetObject(Form("h_vz_C%d",            iuf), h_vz_tmp);
            if(!h_jet_tmp || !h_vz_tmp){
                printf("ERROR: missing C%d histograms\n", iuf); return;
            }

            if(!h_PbPb[ic]){
                h_PbPb[ic] = (TH1D*) h_jet_tmp->Clone(Form("h_PbPb_nom%d", ic));
                h_PbPb[ic]->SetDirectory(nullptr);
            } else {
                h_PbPb[ic]->Add(h_jet_tmp);
            }
            N_evt_PbPb[ic] += h_vz_tmp->Integral();
        }
        printf("Nom-C%d (%s): N_evt_PbPb = %.0f\n", ic+1, nomLabel[ic], N_evt_PbPb[ic]);
    }
    f_pp->Close();
    f_PbPb->Close();

    // ------------------------------------------------------------------
    // Per-event normalization
    // ------------------------------------------------------------------
    TH1D *h_pp = (TH1D*) h_pp_raw->Clone("h_pp");
    h_pp->SetDirectory(nullptr);
    h_pp->Scale(1. / N_evt_pp);

    for(int ic = 0; ic < nNom; ic++)
        h_PbPb[ic]->Scale(1. / N_evt_PbPb[ic]);

    // ------------------------------------------------------------------
    // Rebin + divide by binwidth + eta normalization (|eta|<1.6 => width 3.2)
    // ------------------------------------------------------------------
    const int N_edge = 32;
    double newAxis[N_edge] = {
        60,65,70,75,80,85,90,95,
        100,105,110,115,120,125,130,135,140,145,
        150,160,170,180,190,200,
        220,240,260,280,300,350,400,500
    };

    h_pp = (TH1D*) h_pp->Rebin(N_edge-1, "h_pp_rb", newAxis);
    divideByBinwidth(h_pp);
    h_pp->Scale(1./3.2);

    for(int ic = 0; ic < nNom; ic++){
        h_PbPb[ic] = (TH1D*) h_PbPb[ic]->Rebin(N_edge-1,
                                                  Form("h_PbPb_nom%d_rb", ic),
                                                  newAxis);
        divideByBinwidth(h_PbPb[ic]);
        h_PbPb[ic]->Scale(1./3.2);
    }

    // ------------------------------------------------------------------
    // Style
    // ------------------------------------------------------------------
    Int_t col_C1 = TColor::GetColor("#D55E00");  // vermilion  (most central)
    Int_t col_C2 = TColor::GetColor("#0072B2");  // blue
    Int_t col_C3 = TColor::GetColor("#009E73");  // green
    Int_t col_C4 = TColor::GetColor("#E69F00");  // orange     (most peripheral)
    const int cols[nNom]   = {col_C1, col_C2, col_C3, col_C4};
    const int markers[nNom] = {20, 21, 22, 33};
    const double lw = 2.0, ms = 0.9;

    auto style = [&](TH1D* h, int ic){
        h->SetLineColor(cols[ic]); h->SetMarkerColor(cols[ic]);
        h->SetLineWidth(lw); h->SetMarkerStyle(markers[ic]); h->SetMarkerSize(ms);
        h->SetStats(0); h->SetTitle("");
    };
    h_pp->SetLineColor(kBlack); h_pp->SetMarkerColor(kBlack);
    h_pp->SetLineWidth(lw); h_pp->SetMarkerStyle(29); h_pp->SetMarkerSize(1.2);
    h_pp->SetLineStyle(2); h_pp->SetStats(0); h_pp->SetTitle("");
    for(int ic = 0; ic < nNom; ic++) style(h_PbPb[ic], ic);

    const double xLo = newAxis[0], xHi = newAxis[N_edge-1];

    // ------------------------------------------------------------------
    // 1. Spectra comparison: per-event dN/dpT/deta
    // ------------------------------------------------------------------
    {
        const double fUp = 0.60, fDn = 1.-fUp;
        const double mL = 0.18, mR = 0.05;

        // project PbPb onto same axis as pp for ratio
        auto projectOntoAxis = [](TH1D* hSrc, TH1D* hRef, const char* name) -> TH1D* {
            TH1D* hOut = (TH1D*) hRef->Clone(name);
            hOut->SetDirectory(nullptr);
            hOut->Reset();
            for(int b = 1; b <= hOut->GetNbinsX(); b++){
                double pT = hOut->GetBinCenter(b);
                int bS    = hSrc->FindBin(pT);
                hOut->SetBinContent(b, hSrc->GetBinContent(bS));
                hOut->SetBinError  (b, hSrc->GetBinError(bS));
            }
            return hOut;
        };
        TH1D *h_pp_rb = projectOntoAxis(h_pp, h_PbPb[0], "h_pp_proj");

        TH1D *r_nom[nNom];
        for(int ic = 0; ic < nNom; ic++){
            r_nom[ic] = (TH1D*) h_PbPb[ic]->Clone(Form("r_pre_nom%d", ic));
            r_nom[ic]->Divide(h_pp_rb);
            style(r_nom[ic], ic);
        }

        TCanvas *c = new TCanvas("canv_spectra", "", 700, 800);
        TPad *pUp = new TPad("pUp","",0.,fDn,1.,1.);
        pUp->SetLeftMargin(mL); pUp->SetRightMargin(mR);
        pUp->SetTopMargin(0.08); pUp->SetBottomMargin(0.02);
        pUp->SetLogy(); pUp->SetLogx(); pUp->Draw(); pUp->cd();

        h_pp->GetXaxis()->SetRangeUser(xLo, xHi);
        h_pp->GetXaxis()->SetLabelSize(0); h_pp->GetXaxis()->SetTitleSize(0);
        h_pp->GetYaxis()->SetTitle("#frac{1}{N^{evt}} #frac{dN^{jet}}{dp_{T} d#eta} [GeV^{-1}]");
        h_pp->GetYaxis()->SetTitleSize(0.052/fUp);
        h_pp->GetYaxis()->SetTitleOffset(1.2);
        h_pp->GetYaxis()->SetLabelSize(0.040/fUp);
        h_pp->SetMinimum(1e-7);
        h_pp->Draw("hist");
        for(int ic = 0; ic < nNom; ic++) h_PbPb[ic]->Draw("hist same");

        TLegend *lg = new TLegend(0.52, 0.52, 0.92, 0.88);
        lg->SetBorderSize(0); lg->SetTextSize(0.042/fUp);
        lg->AddEntry(h_pp, "pp ZeroBias", "l");
        for(int ic = 0; ic < nNom; ic++)
            lg->AddEntry(h_PbPb[ic], Form("PbPb %s", nomLabel[ic]), "l");
        lg->Draw();

        TLatex lat; lat.SetNDC(); lat.SetTextSize(0.038/fUp);
        lat.DrawLatex(0.20, 0.89, "MinBias only  |#eta_{jet}| < 1.6  anti-#it{k}_{T} R=0.4");

        c->cd();
        TPad *pDn = new TPad("pDn","",0.,0.,1.,fDn);
        pDn->SetLeftMargin(mL); pDn->SetRightMargin(mR);
        pDn->SetTopMargin(0.02); pDn->SetBottomMargin(0.22);
        pDn->SetLogx(); pDn->Draw(); pDn->cd();

        r_nom[0]->GetXaxis()->SetRangeUser(xLo, xHi);
        r_nom[0]->GetXaxis()->SetTitle("Jet #it{p}_{T} [GeV]");
        r_nom[0]->GetXaxis()->SetTitleSize(0.055/fDn);
        r_nom[0]->GetXaxis()->SetTitleOffset(1.0);
        r_nom[0]->GetXaxis()->SetLabelSize(0.040/fDn);
        r_nom[0]->GetYaxis()->SetTitle("PbPb / pp (per evt)");
        r_nom[0]->GetYaxis()->SetTitleSize(0.050/fDn);
        r_nom[0]->GetYaxis()->SetTitleOffset(0.75);
        r_nom[0]->GetYaxis()->SetLabelSize(0.040/fDn);
        r_nom[0]->GetYaxis()->SetNdivisions(505);
        r_nom[0]->SetMinimum(0.); r_nom[0]->SetMaximum(4.);
        r_nom[0]->Draw("ep");
        for(int ic = 1; ic < nNom; ic++) r_nom[ic]->Draw("ep same");
        TLine *l1 = new TLine(xLo,1.,xHi,1.);
        l1->SetLineStyle(2); l1->SetLineColor(kGray+1); l1->Draw();

        c->SaveAs("../../figures/jetPt/RAA/spectra_minBiasOnly.pdf");
        delete c;
    }

    // ------------------------------------------------------------------
    // 2. RAA proxy: (PbPb/N_evt) / (T_AA * pp/N_evt)
    //    T_AA values: 23.05, 11.60, 3.95, 0.7452 mb^{-1} for C1-C4
    // ------------------------------------------------------------------
    {
        TH1D *h_pp_proj = (TH1D*) h_pp->Clone("h_pp_proj2");
        h_pp_proj->SetDirectory(nullptr);

        // Reproject pp onto PbPb axis (they may differ in binning)
        TH1D *h_pp_on_PbPbAxis = (TH1D*) h_PbPb[0]->Clone("h_pp_on_PbPbAxis");
        h_pp_on_PbPbAxis->SetDirectory(nullptr);
        h_pp_on_PbPbAxis->Reset();
        for(int b = 1; b <= h_pp_on_PbPbAxis->GetNbinsX(); b++){
            double pT = h_pp_on_PbPbAxis->GetBinCenter(b);
            int bS    = h_pp_proj->FindBin(pT);
            h_pp_on_PbPbAxis->SetBinContent(b, h_pp_proj->GetBinContent(bS));
            h_pp_on_PbPbAxis->SetBinError  (b, h_pp_proj->GetBinError(bS));
        }

        TH1D *h_RAA[nNom];
        for(int ic = 0; ic < nNom; ic++){
            h_RAA[ic] = (TH1D*) h_PbPb[ic]->Clone(Form("h_RAA_nom%d", ic));
            h_RAA[ic]->SetDirectory(nullptr);
            // divide by T_AA * pp
            TH1D *h_pp_scaled = (TH1D*) h_pp_on_PbPbAxis->Clone(Form("h_pp_taa_%d", ic));
            h_pp_scaled->Scale(T_AA[ic]);
            h_RAA[ic]->Divide(h_pp_scaled);
            style(h_RAA[ic], ic);
            delete h_pp_scaled;
        }
        delete h_pp_on_PbPbAxis;

        TCanvas *c = new TCanvas("canv_RAA", "", 750, 650);
        c->SetLeftMargin(0.18); c->SetRightMargin(0.06);
        c->SetBottomMargin(0.16); c->SetTopMargin(0.08);
        c->SetLogx();

        h_RAA[0]->GetXaxis()->SetRangeUser(xLo, xHi);
        h_RAA[0]->GetXaxis()->SetTitle("Jet #it{p}_{T} [GeV]");
        h_RAA[0]->GetYaxis()->SetTitle("#it{R}_{AA}^{jet} (MinBias proxy)");
        h_RAA[0]->GetXaxis()->SetTitleSize(0.052);
        h_RAA[0]->GetXaxis()->SetLabelSize(0.040);
        h_RAA[0]->GetYaxis()->SetTitleSize(0.052);
        h_RAA[0]->GetYaxis()->SetLabelSize(0.040);
        h_RAA[0]->GetYaxis()->SetTitleOffset(1.6);
        h_RAA[0]->SetMinimum(0.); h_RAA[0]->SetMaximum(2.5);
        h_RAA[0]->Draw("ep");
        for(int ic = 1; ic < nNom; ic++) h_RAA[ic]->Draw("ep same");

        TLine *l1 = new TLine(xLo,1.,xHi,1.);
        l1->SetLineStyle(2); l1->SetLineColor(kGray+1); l1->Draw();

        TLegend *lg = new TLegend(0.52, 0.60, 0.92, 0.88);
        lg->SetBorderSize(0); lg->SetTextSize(0.042);
        for(int ic = 0; ic < nNom; ic++)
            lg->AddEntry(h_RAA[ic], Form("PbPb %s", nomLabel[ic]), "lep");
        lg->Draw();

        TLatex lat; lat.SetNDC(); lat.SetTextSize(0.038);
        lat.DrawLatex(0.20, 0.895, "MinBias only  |#eta_{jet}| < 1.6  anti-#it{k}_{T} R=0.4");
        lat.SetTextSize(0.033);
        lat.DrawLatex(0.20, 0.845, "No unfolding, no fake-jet subtraction");

        c->SaveAs("../../figures/jetPt/RAA/RAA_minBiasOnly.pdf");
        delete c;
    }

    // ------------------------------------------------------------------
    // 3. Per-centrality canvas: PbPb + pp spectra one pad per cent
    // ------------------------------------------------------------------
    {
        TCanvas *c = new TCanvas("canv_spectra_grid", "", 1200, 600);
        c->Divide(4, 1, 0.005, 0.005);

        for(int ic = 0; ic < nNom; ic++){
            c->cd(ic+1);
            gPad->SetLogy(); gPad->SetLogx();
            gPad->SetLeftMargin(0.20); gPad->SetRightMargin(0.05);
            gPad->SetBottomMargin(0.18); gPad->SetTopMargin(0.10);

            h_PbPb[ic]->GetXaxis()->SetRangeUser(xLo, xHi);
            h_PbPb[ic]->GetYaxis()->SetTitle("#frac{1}{N^{evt}} #frac{dN^{jet}}{dp_{T} d#eta}");
            h_PbPb[ic]->GetYaxis()->SetTitleSize(0.064);
            h_PbPb[ic]->GetYaxis()->SetLabelSize(0.054);
            h_PbPb[ic]->GetXaxis()->SetTitleSize(0.060);
            h_PbPb[ic]->GetXaxis()->SetLabelSize(0.050);
            h_PbPb[ic]->GetXaxis()->SetTitle("Jet #it{p}_{T} [GeV]");
            h_PbPb[ic]->GetYaxis()->SetTitleOffset(1.5);
            h_PbPb[ic]->SetMinimum(1e-7);
            h_PbPb[ic]->Draw("ep");
            h_pp->GetXaxis()->SetRangeUser(xLo, xHi);
            h_pp->Draw("hist same");

            TLatex lat; lat.SetNDC(); lat.SetTextSize(0.060);
            lat.DrawLatex(0.24, 0.88, Form("PbPb %s", nomLabel[ic]));

            if(ic == 0){
                TLegend *lg = new TLegend(0.22, 0.22, 0.80, 0.46);
                lg->SetBorderSize(0); lg->SetTextSize(0.055);
                lg->AddEntry(h_pp,       "pp ZeroBias", "l");
                lg->AddEntry(h_PbPb[ic], "PbPb MinBias", "lep");
                lg->Draw();
            }
        }
        c->SaveAs("../../figures/jetPt/RAA/spectra_grid_minBiasOnly.pdf");
        delete c;
    }

    printf("\nDone. Figures saved in figures/jetPt/RAA/\n");
}
