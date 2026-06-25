// Compare mixed-event FastJet anti-kT R=0.4 fake-jet spectra between
// PbPb data (MinBias_test) and HYDJET MC, normalised per random cone event.
//
// PbPb centrality bins (hiBin / 2 = centrality %):
//   C1: hiBin 0-20   = 0-10%
//   C2: hiBin 20-60  = 10-30%
//   C3: hiBin 60-100 = 30-50%
//   C4: hiBin 100-160= 50-80%
//
// HYDJET ultra-fine bins are summed to match:
//   0-10%  : HY C1(0-10)+C2(10-20)
//   10-30% : HY C3..C6(20-60)
//   30-50% : HY C7..C10(60-100)
//   50-80% : HY C11..C16(100-160)

#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TROOT.h"

static const int N_genPJ = 100;

// Sum HYDJET fine-bin fastJet histograms and return a clone normalised per event.
// sumPJ accumulates the pseudoJet integral (= N_events * N_genPJ).
TH1D* sumHYDJET(TFile* fHY, const char* hname_base, int cFirst, int cLast,
                double& sumEvents)
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
    sumEvents = pjIntegral / N_genPJ;
    if(hSum) hSum->SetDirectory(nullptr);
    return hSum;
}

void drawComparison(TH1D* hPbPb, TH1D* hHY,
                    const char* title, const char* outPath)
{
    // Normalise per event
    TH1D* hP = (TH1D*)hPbPb->Clone("hP_tmp");
    TH1D* hH = (TH1D*)hHY  ->Clone("hH_tmp");
    hP->SetDirectory(nullptr);
    hH->SetDirectory(nullptr);
    if(hP->Integral() > 0) hP->Scale(1. / hP->Integral());
    if(hH->Integral() > 0) hH->Scale(1. / hH->Integral());

    hP->SetLineColor(kBlack);    hP->SetLineWidth(2);
    hH->SetLineColor(kAzure-4);  hH->SetLineWidth(2); hH->SetLineStyle(2);

    hP->SetTitle("");
    hP->GetXaxis()->SetTitle("Fake jet p_{T} [GeV]");
    hP->GetYaxis()->SetTitle("Normalised entries / bin");
    hP->GetXaxis()->SetTitleSize(0.05); hP->GetYaxis()->SetTitleSize(0.05);
    hP->GetXaxis()->SetLabelSize(0.04); hP->GetYaxis()->SetLabelSize(0.04);
    hP->GetXaxis()->SetRangeUser(20, 200);

    double ymax = TMath::Max(hP->GetMaximum(), hH->GetMaximum());
    hP->SetMaximum(ymax * 5.);
    hP->SetMinimum(1e-7);

    TCanvas* c = new TCanvas("c","",600,600);
    c->SetLogy();
    hP->Draw("hist"); hH->Draw("hist same");

    TLegend* lg = new TLegend(0.55, 0.68, 0.88, 0.88);
    lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.04);
    lg->AddEntry(hP, "PbPb data (MinBias)", "l");
    lg->AddEntry(hH, "HYDJET MC",           "l");
    lg->Draw();

    TLatex lat;
    lat.SetNDC(); lat.SetTextSize(0.04);
    lat.DrawLatex(0.13, 0.92, title);
    lat.SetTextSize(0.035);
    lat.DrawLatex(0.13, 0.87, "Mixed-event anti-k_{T} R=0.4, CS PF cands");

    c->SaveAs(outPath);
    delete c; delete hP; delete hH;
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
        double ev; TH1D* hH = sumHYDJET(fHY, "h_fastJetPt_JEC", 1, 16, ev);
        drawComparison(hP, hH, "Inclusive  0--80%",
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
        TH1D* hP = (TH1D*)fPbPb->Get(Form("h_fastJetPt_JEC_C%d", b.pbC));
        double ev; TH1D* hH = sumHYDJET(fHY, "h_fastJetPt_JEC", b.hyFirst, b.hyLast, ev);
        drawComparison(hP, hH, b.label,
                       Form("%sfastJet_JEC_PbPb_vs_HYDJET_C%d.pdf", outDir, b.pbC));
        delete hH;
    }

    fPbPb->Close(); fHY->Close();
    printf("\nPlots saved to %s\n", outDir);
}
