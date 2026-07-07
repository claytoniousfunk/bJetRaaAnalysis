// Build rootFiles/fakeJets/fakeJets.root for use in calculateRAA.C.
//
// Fake-jet estimate = FastJet(JEC) - RandomCone, both per event.
//   FastJet  : h_fastJetPt_JEC_C{i} / N_events          (one clustering per event)
//   RandomCone: h_pseudoJetPt_C{i}  / (N_events*N_pool)  (N_pool cones per event)
//
// Subtracting the random-cone baseline removes the uniform UE contribution
// and leaves only the excess from anti-kT clustering, i.e. the jet-like fake rate.
// Bins that go negative after subtraction are zeroed (statistical fluctuations).
//
// Coarse-bin mapping (matches calculateRAA.C):
//   C1 = 0-10%  : hiBin 0-20
//   C2 = 10-30% : hiBin 20-60
//   C3 = 30-50% : hiBin 60-100
//   C4 = 50-80% : hiBin 100-160

const char *fmixed_path =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/"
  "PbPb_MinBias_Part1_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_"
  "mixedEventPseudoJets_pfCand_pseudoJetCandPtMin-0.0_2026-7-7.root";

const char *outPath =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/fakeJets/fakeJets.root";

const int    N_pool  = 100;   // N_mixedEventsInPool
const double pT_min  = 20.;   // zero out bins below this threshold

void makeFakeJetFile(){
  gSystem->Exec("mkdir -p /home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/fakeJets");

  TFile *fM = TFile::Open(fmixed_path);
  if(!fM || fM->IsZombie()){ printf("ERROR: cannot open mixed-event file\n"); return; }

  TFile *fOut = TFile::Open(outPath, "RECREATE");

  for(int ci = 0; ci <= 4; ci++){
    TH1D *hFJ = nullptr, *hRC = nullptr, *hvz = nullptr;
    fM->GetObject(Form("h_fastJetPt_JEC_C%d", ci), hFJ);
    fM->GetObject(Form("h_pseudoJetPt_C%d",   ci), hRC);
    fM->GetObject(Form("h_vz_C%d",             ci), hvz);
    if(!hFJ || !hRC || !hvz){
      printf("WARNING: missing histograms for C%d\n", ci); continue;
    }

    double N_events = hvz->Integral();
    if(N_events <= 0){ printf("WARNING: N_events=0 for C%d, skipping\n", ci); continue; }

    TH1D *hFJ_norm = (TH1D*) hFJ->Clone(Form("hFJ_C%d", ci));
    TH1D *hRC_norm = (TH1D*) hRC->Clone(Form("hRC_C%d", ci));
    hFJ_norm->SetDirectory(nullptr);
    hRC_norm->SetDirectory(nullptr);
    hFJ_norm->Scale(1. / N_events);
    hRC_norm->Scale(1. / (N_events * double(N_pool)));

    // difference: FastJet(JEC) - RandomCone
    TH1D *hOut = (TH1D*) hFJ_norm->Clone(Form("h_fakeJets_C%d", ci));
    hOut->SetDirectory(nullptr);
    hOut->Add(hRC_norm, -1.);

    // zero bins below pT_min or that went negative
    for(int b = 1; b <= hOut->GetNbinsX(); b++){
      if(hOut->GetXaxis()->GetBinUpEdge(b) <= pT_min ||
         hOut->GetBinContent(b) < 0.){
        hOut->SetBinContent(b, 0.);
        hOut->SetBinError(b, 0.);
      }
    }

    printf("h_fakeJets_C%d: N_events=%.0f  FJ/evt=%.4f  RC/evt=%.4f  diff=%.4f\n",
           ci, N_events, hFJ_norm->Integral(), hRC_norm->Integral(), hOut->Integral());

    fOut->cd();
    hOut->Write();
    delete hFJ_norm; delete hRC_norm; delete hOut;
  }

  fOut->Close();
  fM->Close();
  printf("\nWrote %s\n", outPath);
}
