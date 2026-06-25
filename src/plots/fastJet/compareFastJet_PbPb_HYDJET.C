// Compare mixed-event FastJet anti-kT R=0.4 fake-jet spectra between
// PbPb data (MinBias_test) and HYDJET MC, normalised to fake jets per event.
//
// PbPb event count: integral of h_vz_C{i} (inclusive, before any selection).
// HYDJET event count: h_pseudoJetPt_C{i}->Integral() / N_genPJ  (h_vz in HYDJET
//   is gated on inclRecoJet and would undercount; pseudoJetPt gives direct count).
//
// PbPb centrality bins (hiBin / 2 = centrality %):
//   C1: hiBin 0-20   = 0-10%
//   C2: hiBin 20-60  = 10-30%
//   C3: hiBin 60-100 = 30-50%
//   C4: hiBin 100-160= 50-80%
//
// HYDJET ultra-fine bins summed to match:
//   0-10%  : HY C1(0-10)+C2(10-20)
//   10-30% : HY C3..C6(20-60)
//   30-50% : HY C7..C10(60-100)
//   50-80% : HY C11..C16(100-160)

#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TStyle.h"
#include "TROOT.h"

static const int N_genPJ = 100;

// Return PbPb event count from the inclusive vz histogram for centrality index c.
double getEventCount_PbPb(TFile* fPbPb, int c)
{
    TH1D* h = (TH1D*)fPbPb->Get(Form("h_vz_C%d", c));
    if(!h){ printf("WARNING: h_vz_C%d not found in PbPb file\n", c); return 1.; }
    return h->Integral();
}

// Sum HYDJET fine-bin fastJet histograms across cFirst..cLast.
// Returns a cloned, unnormalised sum (caller owns it).
// HYDJET event count comes from pseudoJetPt integral / N_genPJ.
TH1D* sumHYDJET(TFile* fHY, const char* hname_base, int cFirst, int cLast,
                double& nEvents)
{
    TH1D* hSum = nullptr;
    double pjIntegral = 0.;
    for(int c = cFirst; c <= cLast; c++){
        TH1D* hFJ = (TH1D*)fHY->Get(Form("%s_C%d", hname_base, c));
        TH1D* hPJ = (TH1D*)fHY->Get(Form("h_pseudoJetPt_C%d", c));
        if(!hFJ || !hPJ){ printf("WARNING: missing C%d in HYDJET\n", c); continue; }
        pjIntegral += hPJ->Integral();
        if(!hSum) hSum = (TH1D*)hFJ->Clone(Form("hHY_%s_%d_%d", hname_base, cFirst, cLast));
        else       hSum->Add(hFJ);
    }
    nEvents = pjIntegral / N_genPJ;
    if(hSum) hSum->SetDirectory(nullptr);
    return hSum;
}

void drawComparison(TH1D* hPbPb, double nEventsPbPb,
                    TH1D* hHY,   double nEventsHY,
                    const char* centLabel, const char* outPath)
{
    // Pad split fractions
    const double fUp   = 0.68;
    const double fDown = 1. - fUp;
    const double mL = 0.20, mR = 0.05;
    const double mTop = 0.08, mBotUp = 0.02;
    const double mTopDn = 0.02, mBotDn = 0.32;
    const double sfUp   = 1. / fUp;
    const double sfDown = 1. / fDown;

    // Per-event spectra
    TH1D* hP = (TH1D*)hPbPb->Clone("hP_tmp"); hP->SetDirectory(nullptr);
    TH1D* hH = (TH1D*)hHY  ->Clone("hH_tmp"); hH->SetDirectory(nullptr);
    if(nEventsPbPb > 0) hP->Scale(1. / nEventsPbPb);
    if(nEventsHY   > 0) hH->Scale(1. / nEventsHY);

    // Ratio PbPb / HYDJET
    TH1D* hR = (TH1D*)hP->Clone("hR_tmp"); hR->SetDirectory(nullptr);
    hR->Divide(hH);

    // Style
    hP->SetLineColor(kBlack);   hP->SetLineWidth(2);
    hH->SetLineColor(kAzure-4); hH->SetLineWidth(2); hH->SetLineStyle(2);
    hR->SetLineColor(kBlack);   hR->SetLineWidth(2);
    hR->SetMarkerStyle(20);     hR->SetMarkerSize(0.7);

    const double xLo = 20., xHi = 200.;
    hP->GetXaxis()->SetRangeUser(xLo, xHi);
    hH->GetXaxis()->SetRangeUser(xLo, xHi);
    hR->GetXaxis()->SetRangeUser(xLo, xHi);

    // Upper pad axes
    hP->SetTitle("");
    hP->GetXaxis()->SetLabelSize(0);
    hP->GetXaxis()->SetTitleSize(0);
    hP->GetYaxis()->SetTitle("Fake jets / event / bin");
    hP->GetYaxis()->SetTitleSize(0.05 * sfUp);
    hP->GetYaxis()->SetTitleOffset(1.4);
    hP->GetYaxis()->SetLabelSize(0.04 * sfUp);
    double ymax = TMath::Max(hP->GetMaximum(), hH->GetMaximum());
    hP->SetMaximum(ymax * 5.);
    hP->SetMinimum(1e-9);

    // Lower pad axes
    hR->SetTitle("");
    hR->GetXaxis()->SetTitle("Fake jet p_{T} [GeV]");
    hR->GetXaxis()->SetTitleSize(0.05 * sfDown);
    hR->GetXaxis()->SetTitleOffset(1.0);
    hR->GetXaxis()->SetLabelSize(0.04 * sfDown);
    hR->GetYaxis()->SetTitle("PbPb / HYDJET");
    hR->GetYaxis()->SetTitleSize(0.045 * sfDown);
    hR->GetYaxis()->SetTitleOffset(0.75);
    hR->GetYaxis()->SetLabelSize(0.04 * sfDown);
    hR->GetYaxis()->SetNdivisions(504);
    hR->SetMaximum(3.0);
    hR->SetMinimum(0.);

    // Canvas
    TCanvas* c = new TCanvas("c", "", 600, 700);
    c->SetFillColor(0);

    TPad* pUp = new TPad("pUp", "", 0., fDown, 1., 1.);
    pUp->SetLeftMargin(mL);  pUp->SetRightMargin(mR);
    pUp->SetTopMargin(mTop); pUp->SetBottomMargin(mBotUp);
    pUp->SetLogy(); pUp->SetTickx(1); pUp->SetTicky(1);
    pUp->Draw(); pUp->cd();

    hP->Draw("hist");
    hH->Draw("hist same");

    TLegend* lg = new TLegend(0.54, 0.60, 0.93, 0.84);
    lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.04 * sfUp);
    lg->AddEntry(hP, "PbPb data (MinBias)", "l");
    lg->AddEntry(hH, "HYDJET MC",           "l");
    lg->Draw();

    TLatex lat;
    lat.SetNDC(); lat.SetTextSize(0.042 * sfUp);
    lat.DrawLatex(0.18, 0.88, centLabel);
    lat.SetTextSize(0.036 * sfUp);
    lat.DrawLatex(0.18, 0.79, "Mixed-event anti-k_{T} R=0.4, CS PF cands");

    c->cd();
    TPad* pDn = new TPad("pDn", "", 0., 0., 1., fDown);
    pDn->SetLeftMargin(mL);    pDn->SetRightMargin(mR);
    pDn->SetTopMargin(mTopDn); pDn->SetBottomMargin(mBotDn);
    pDn->SetTickx(1); pDn->SetTicky(1);
    pDn->Draw(); pDn->cd();

    hR->Draw("ep");

    TLine* line = new TLine(xLo, 1., xHi, 1.);
    line->SetLineStyle(2); line->SetLineColor(kGray+1); line->SetLineWidth(1);
    line->Draw();

    c->SaveAs(outPath);
    delete c;
    delete hP; delete hH; delete hR;
}

void compareFastJet_PbPb_HYDJET()
{
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);

    const char* fPbPbPath =
        "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/"
        "PbPb_MinBias_test_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_"
        "mixedEventPseudoJets_pfCandCS_pseudoJetCandPtMin-0.0_subleadingPFCandPtMin-15_"
        "2026-6-24_oneFile.root";

    const char* fHYPath =
        "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/HYDJET/"
        "HYDJET_pThat-unweighted_mu12_pTmu-15to999_tight_hiBinShift-0_jetTrkMaxFilter_"
        "WDecayFilter_mixedEventPseudoJets_pfCandCS_pseudoJetCandPtMin-0.0_"
        "subleadingPFCandPtMin-15_2026-6-23_ultraFineCentBins.root";

    const char* outDir =
        "/home/clayton/Documents/nuclear/GroupMeeting/figures/2026-06-30/";

    TFile* fPbPb = TFile::Open(fPbPbPath);
    TFile* fHY   = TFile::Open(fHYPath);
    if(!fPbPb || fPbPb->IsZombie()){ printf("ERROR: cannot open PbPb file\n"); return; }
    if(!fHY   || fHY  ->IsZombie()){ printf("ERROR: cannot open HYDJET file\n"); return; }

    // ---- Inclusive C0 ----
    {
        TH1D* hP = (TH1D*)fPbPb->Get("h_fastJetPt_JEC_C0");
        double nPbPb = getEventCount_PbPb(fPbPb, 0);
        double nHY;  TH1D* hH = sumHYDJET(fHY, "h_fastJetPt_JEC", 1, 16, nHY);
        printf("C0 inclusive: PbPb events = %.0f, HYDJET events = %.0f\n", nPbPb, nHY);
        drawComparison(hP, nPbPb, hH, nHY, "Inclusive  0--80%",
                       Form("%sfastJet_JEC_PbPb_vs_HYDJET_C0.pdf", outDir));
        delete hH;
    }

    struct CentBin { int pbC; int hyFirst; int hyLast; const char* label; };
    CentBin bins[] = {
        {1,  1,  2,  "0--10%"},
        {2,  3,  6,  "10--30%"},
        {3,  7,  10, "30--50%"},
        {4, 11,  16, "50--80%"},
    };

    for(auto& b : bins){
        TH1D* hP  = (TH1D*)fPbPb->Get(Form("h_fastJetPt_JEC_C%d", b.pbC));
        double nPbPb = getEventCount_PbPb(fPbPb, b.pbC);
        double nHY;  TH1D* hH = sumHYDJET(fHY, "h_fastJetPt_JEC", b.hyFirst, b.hyLast, nHY);
        printf("C%d %s: PbPb events = %.0f, HYDJET events = %.0f\n",
               b.pbC, b.label, nPbPb, nHY);
        drawComparison(hP, nPbPb, hH, nHY, b.label,
                       Form("%sfastJet_JEC_PbPb_vs_HYDJET_C%d.pdf", outDir, b.pbC));
        delete hH;
    }

    fPbPb->Close(); fHY->Close();
    printf("\nPlots saved to %s\n", outDir);
}
