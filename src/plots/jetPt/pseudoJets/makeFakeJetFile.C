// Build rootFiles/fakeJets/fakeJets.root for use in calculateRAA.C.
//
// Fake-jet estimate = FastJet(no JEC, random-cone pT subtracted), per event.
//   FastJet  : h_fastJetPt_PF_bkgSub_RC_C{i} / N_events          (one clustering per event)
//   RandomCone: h_pseudoJetPt_C{i}           / (N_events*N_pool)  (N_pool cones per event)
//
// The random-cone baseline subtraction is currently disabled (see below); the
// per-jet RC pT subtraction is already folded into h_fastJetPt_PF_bkgSub_RC.
// Bins that go negative after subtraction are zeroed (statistical fluctuations).
//
// The input file is binned in ultra-fine centrality slices (10 hiBin units each),
// which do NOT match the coarse classes used by calculateRAA.C. The slices are
// summed here -- raw counts and event counts added separately, then divided --
// so the output is the per-event fake rate of the full coarse class:
//   C1 = 0-10%  : hiBin   0-20  = slices  1- 2
//   C2 = 10-30% : hiBin  20-60  = slices  3- 6
//   C3 = 30-50% : hiBin  60-100 = slices  7-10
//   C4 = 50-80% : hiBin 100-160 = slices 11-16
// C0 (hiBin 0-160, inclusive) is taken directly from the input C0 histogram.
//
// The input sample is unweighted (flat in hiBin). An earlier reweighted-to-Jet80
// scan was used here; within a coarse class that shifted the composition toward
// the central edge and inflated the per-event rate by 1.05 (0-10%) up to 1.71
// (50-80%), which over-subtracted by the same factor. Do not point this at a
// hiBinReweight* file without forming the fake fraction per ultra-fine slice.

// Which FastJet spectrum to use as the fake estimate.
//   h_fastJetPt_PF_bkgSub_RC      : raw pT minus the random-cone UE
//   h_fastJetPt_PF_JEC_bkgSub_RC  : the same, JEC applied
// h_inclRecoJetPt in calculateRAA.C is JEC-corrected, so the JEC variant is
// the like-for-like choice; the raw one places every fake at too low a pT.
const char *fjHistBase = "h_fastJetPt_PF_JEC_bkgSub_RC";

const char *fmixed_path =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/PbPb_MinBias_Part1_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_mixedEventPFClustering_pseudoJetCandPtMin-0.0_2026-8-11_ultraFineCentBins.root";

const char *outPath =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/fakeJets/fakeJets.root";

const int    N_pool  = 100;   // N_mixedEventsInPool
const double pT_min  = 20.;   // zero out bins below this threshold

// coarse class -> [first, last] ultra-fine slice index (inclusive).
// Index 0 is the inclusive hiBin 0-160 histogram, taken as-is.
const int sliceLo[5] = {0,  1,  3,  7, 11};
const int sliceHi[5] = {0,  2,  6, 10, 16};

void makeFakeJetFile(){
  gSystem->Exec("mkdir -p /home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/fakeJets");

  TFile *fM = TFile::Open(fmixed_path);
  if(!fM || fM->IsZombie()){ printf("ERROR: cannot open mixed-event file\n"); return; }

  TFile *fOut = TFile::Open(outPath, "RECREATE");

  for(int ci = 0; ci <= 4; ci++){

    // --- sum the raw (un-normalised) slice histograms and their event counts ---
    TH1D *hFJ_sum = nullptr, *hRC_sum = nullptr;
    double N_events = 0.;
    bool ok = true;

    for(int si = sliceLo[ci]; si <= sliceHi[ci]; si++){
      TH1D *hFJ = nullptr, *hRC = nullptr, *hvz = nullptr;
      fM->GetObject(Form("%s_C%d", fjHistBase, si), hFJ);
      fM->GetObject(Form("h_pseudoJetPt_C%d",             si), hRC);
      fM->GetObject(Form("h_vz_C%d",                      si), hvz);
      if(!hFJ || !hRC || !hvz){
        printf("WARNING: missing histograms for slice C%d (coarse C%d)\n", si, ci); ok = false; break;
      }

      if(!hFJ_sum){
        hFJ_sum = (TH1D*) hFJ->Clone(Form("hFJ_sum_C%d", ci));
        hRC_sum = (TH1D*) hRC->Clone(Form("hRC_sum_C%d", ci));
        hFJ_sum->SetDirectory(nullptr);
        hRC_sum->SetDirectory(nullptr);
      }
      else{
        hFJ_sum->Add(hFJ);
        hRC_sum->Add(hRC);
      }
      N_events += hvz->Integral();
    }

    if(!ok || !hFJ_sum) continue;
    if(N_events <= 0){ printf("WARNING: N_events=0 for C%d, skipping\n", ci); continue; }

    // --- per-event normalisation of the summed class ---
    TH1D *hFJ_norm = hFJ_sum;
    TH1D *hRC_norm = hRC_sum;
    hFJ_norm->Scale(1. / N_events);
    hRC_norm->Scale(1. / (N_events * double(N_pool)));

    // difference: FastJet - RandomCone
    TH1D *hOut = (TH1D*) hFJ_norm->Clone(Form("h_fakeJets_C%d", ci));
    hOut->SetDirectory(nullptr);
    //hOut->Add(hRC_norm, -1.);

    // zero bins below pT_min or that went negative
    for(int b = 1; b <= hOut->GetNbinsX(); b++){
      if(hOut->GetXaxis()->GetBinUpEdge(b) <= pT_min ||
         hOut->GetBinContent(b) < 0.){
        hOut->SetBinContent(b, 0.);
        hOut->SetBinError(b, 0.);
      }
    }

    printf("h_fakeJets_C%d: slices %d-%d  N_events=%.0f  FJ/evt=%.4f  RC/evt=%.4f  diff=%.4f\n",
           ci, sliceLo[ci], sliceHi[ci], N_events,
           hFJ_norm->Integral(), hRC_norm->Integral(), hOut->Integral());

    fOut->cd();
    hOut->Write();
    delete hFJ_norm; delete hRC_norm; delete hOut;
  }

  fOut->Close();
  fM->Close();
  printf("\nWrote %s\n", outPath);
}
