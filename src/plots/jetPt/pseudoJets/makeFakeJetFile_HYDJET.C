// Build rootFiles/fakeJets/fakeJets_HYDJET.root from the HYDJET mixed-event
// FastJet anti-kT R=0.4 fake-jet spectrum, for comparison with data-driven fake jets.
//
// HYDJET uses ultra-fine centrality bins (C1-C16). Summed to match the coarse bins
// used in calculateRAA.C:
//   C1 (0-10%)  : HY C1+C2    (hiBin 0-20)
//   C2 (10-30%) : HY C3-C6    (hiBin 20-60)
//   C3 (30-50%) : HY C7-C10   (hiBin 60-100)
//   C4 (50-80%) : HY C11-C16  (hiBin 100-160)
//
// Event count: pseudoJetPt integral / N_genPJ (h_vz in HYDJET is gated on
// inclRecoJet and undercounts events).

static const int N_genPJ = 100;
static const double pT_min = 20.;

const char *fHY_path =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/HYDJET/"
  "HYDJET_pThat-unweighted_mu12_pTmu-15to999_tight_hiBinShift-0_jetTrkMaxFilter_"
  "WDecayFilter_mixedEventPseudoJets_pfCandCS_pseudoJetCandPtMin-0.0_"
  "subleadingPFCandPtMin-15_2026-6-23_ultraFineCentBins.root";

const char *outPath_HY =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/fakeJets/fakeJets_HYDJET.root";

struct CoarseBin { int coarseIdx; int hyFirst; int hyLast; };

void makeFakeJetFile_HYDJET()
{
    gSystem->Exec("mkdir -p /home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/fakeJets");

    TFile *fHY = TFile::Open(fHY_path);
    if(!fHY || fHY->IsZombie()){ printf("ERROR: cannot open HYDJET file\n"); return; }

    TFile *fOut = TFile::Open(outPath_HY, "RECREATE");

    CoarseBin bins[] = {
        {0,  1, 16},   // inclusive 0-80%
        {1,  1,  2},
        {2,  3,  6},
        {3,  7,  10},
        {4, 11,  16},
    };

    for(auto& b : bins){
        TH1D *hSum = nullptr;
        double pjIntegral = 0.;
        for(int c = b.hyFirst; c <= b.hyLast; c++){
            TH1D *hFJ = (TH1D*) fHY->Get(Form("h_fastJetPt_JEC_C%d", c));
            TH1D *hPJ = (TH1D*) fHY->Get(Form("h_pseudoJetPt_C%d", c));
            if(!hFJ || !hPJ){ printf("WARNING: missing C%d in HYDJET\n", c); continue; }
            pjIntegral += hPJ->Integral();
            if(!hSum) hSum = (TH1D*) hFJ->Clone(Form("h_fakeJets_C%d", b.coarseIdx));
            else       hSum->Add(hFJ);
        }
        if(!hSum){ printf("ERROR: no histograms for coarse C%d\n", b.coarseIdx); continue; }
        hSum->SetDirectory(nullptr);

        double N_events = pjIntegral / N_genPJ;
        hSum->Scale(1. / N_events);

        for(int bi = 1; bi <= hSum->GetNbinsX(); bi++){
            if(hSum->GetXaxis()->GetBinUpEdge(bi) <= pT_min){
                hSum->SetBinContent(bi, 0.);
                hSum->SetBinError(bi, 0.);
            }
        }

        printf("h_fakeJets_C%d (HYDJET): N_events=%.1f  integral=%.6f (fake jets/event, pT > %.0f GeV)\n",
               b.coarseIdx, N_events, hSum->Integral(), pT_min);

        fOut->cd();
        hSum->Write();
        delete hSum;
    }

    fOut->Close();
    fHY->Close();
    printf("\nWrote %s\n", outPath_HY);
}
