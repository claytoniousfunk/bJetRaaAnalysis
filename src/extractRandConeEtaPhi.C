void extractRandConeEtaPhi()
{
    const char *inPath =
        "../rootFiles/scanningOuput/PbPb/"
        "PbPb_MinBias_Part1_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_"
        "mixedEventPseudoJets_pfCand_pseudoJetCandPtMin-0.0_2026-7-8_ultraFineCentBins.root";

    const char *outPath =
        "../rootFiles/randConeEtaPhi/"
        "PbPb_MinBias_Part1_randConeEtaPhi_2026-7-8_ultraFineCentBins.root";

    TFile *fIn  = TFile::Open(inPath, "READ");
    TFile *fOut = TFile::Open(outPath, "RECREATE");

    // C0 = inclusive 0-80%, C1-C16 = ultraFine 5% bins
    for(int ic = 0; ic <= 16; ic++){
        TProfile2D *h = nullptr;
        fIn->GetObject(Form("h_randConeEtaPhi_C%d", ic), h);
        if(!h){
            std::cerr << "WARNING: h_randConeEtaPhi_C" << ic << " not found\n";
            continue;
        }
        fOut->cd();
        h->Write();
    }

    fOut->Close();
    fIn->Close();
    std::cout << "Saved to " << outPath << "\n";
}
