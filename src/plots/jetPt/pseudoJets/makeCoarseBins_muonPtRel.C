// Merge the ultra-fine (5%) centrality slices of h_fastJetMuonPtRel_fastJetPt_PF_bkgSub_RC_C*
// and h_vz_C* into the standard 4 coarse classes, and write them to a
// separate file. Straight histogram sums -- no per-event normalisation -- so
// the output can still be normalised however the downstream macro needs.
//
// Coarse class -> [first, last] ultra-fine slice index, same convention as
// makeFakeJetFile.C:
//   C0 = 0-80%  (inclusive) : taken directly from the input C0 histogram
//   C1 = 0-10%              : slices  1- 2
//   C2 = 10-30%             : slices  3- 6
//   C3 = 30-50%             : slices  7-10
//   C4 = 50-80%             : slices 11-16
//
// Usage: root -l -b -q 'makeCoarseBins_muonPtRel.C'
// Run from: src/plots/jetPt/pseudoJets/

const char *inPath =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/PbPb_SingleMuon_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_mixedEventPFClustering_pseudoJetCandPtMin-0.0_2026-8-27_ultraFineCentBins_partial.root";

const char *outPath =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/PbPb_SingleMuon_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_mixedEventPFClustering_pseudoJetCandPtMin-0.0_2026-8-27_coarseBins_partial.root";

const char *ptRelHistBase = "h_fastJetMuonPtRel_fastJetPt_PF_bkgSub_RC";
const char *vzHistBase    = "h_vz";

const int NClass = 5;
const int sliceLo[NClass] = {0,  1,  3,  7, 11};
const int sliceHi[NClass] = {0,  2,  6, 10, 16};
const char *classLabel[NClass] = {"0-80% (incl.)", "0-10%", "10-30%", "30-50%", "50-80%"};
const int hiBinLo[NClass] = {0, 0, 20, 60, 100};
const int hiBinHi[NClass] = {160, 20, 60, 100, 160};

void makeCoarseBins_muonPtRel(){
  TFile *fIn = TFile::Open(inPath);
  if(!fIn || fIn->IsZombie()){ printf("ERROR: cannot open %s\n", inPath); return; }

  TFile *fOut = TFile::Open(outPath, "RECREATE");

  printf("%-16s %10s %14s %14s\n", "class", "N_events", "ptRel entries", "ptRel integral");
  for(int ci = 0; ci < NClass; ci++){

    TH2D *hPtRel_sum = nullptr;
    TH1D *hVz_sum    = nullptr;
    bool ok = true;

    for(int si = sliceLo[ci]; si <= sliceHi[ci]; si++){
      TH2D *hPtRel = nullptr;
      TH1D *hVz    = nullptr;
      fIn->GetObject(Form("%s_C%d", ptRelHistBase, si), hPtRel);
      fIn->GetObject(Form("%s_C%d", vzHistBase,    si), hVz);
      if(!hPtRel || !hVz){
        printf("WARNING: missing histograms for slice C%d (coarse class %d), skipping class\n", si, ci);
        ok = false; break;
      }

      if(!hPtRel_sum){
        hPtRel_sum = (TH2D*) hPtRel->Clone(Form("%s_C%d", ptRelHistBase, ci));
        hVz_sum    = (TH1D*) hVz->Clone(Form("%s_C%d", vzHistBase, ci));
        hPtRel_sum->SetDirectory(nullptr);
        hVz_sum->SetDirectory(nullptr);
      }
      else{
        hPtRel_sum->Add(hPtRel);
        hVz_sum->Add(hVz);
      }
    }

    if(!ok || !hPtRel_sum) continue;

    // clone() carries over the FIRST slice's title, which is wrong for a
    // merged range -- reset it to describe the actual coarse hiBin window.
    hPtRel_sum->SetTitle(Form("fastJet muon #it{p}_{T}^{rel} vs fastJet #it{p}_{T}, %d < hiBin < %d",
                               hiBinLo[ci], hiBinHi[ci]));
    hVz_sum->SetTitle(Form("vz, inclusive events, hiBin %d - %d", hiBinLo[ci], hiBinHi[ci]));

    printf("%-16s %10.0f %14.0f %14.2f\n",
           classLabel[ci], hVz_sum->Integral(), hPtRel_sum->GetEntries(), hPtRel_sum->Integral());

    fOut->cd();
    hPtRel_sum->Write();
    hVz_sum->Write();
    delete hPtRel_sum;
    delete hVz_sum;
  }

  fOut->Close();
  fIn->Close();
  printf("\nWrote %s\n", outPath);
}
