// Build rootFiles/fakeJets/fakeJets_dPT_rawPt.root -- raw-pT axis, dPT
// background subtraction. Combines the two variations makeFakeJetFile.C has
// separate single-axis siblings for:
//   makeFakeJetFile_rawPt.C : raw axis,      RC  background subtraction
//   makeFakeJetFile_dPT.C   : JEC-corrected,  dPT background subtraction
// this file                 : raw axis,      dPT background subtraction
//
// Fake-jet estimate = FastJet, raw (non-JEC) pT, dPT (PF-PFCs) per-jet
// background subtraction:
//   FastJet  : h_fastJetPt_PF_bkgSub_dPT_C{i} / N_events   (one clustering per event)
//   RandomCone: h_pseudoJetPt_C{i}            / (N_events*N_pool)  (N_pool cones per event, for reference only)
//
// Raw axis, so this is the correct estimate to subtract from a genuine raw-pT
// reco-jet spectrum (h_inclRawJetPt) -- see makeFakeJetFile_rawPt.C for why the
// JEC-axis variants are an axis mismatch there.
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

const char *fjHistBase = "h_fastJetPt_PF_bkgSub_dPT";

const char *fmixed_path =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/PbPb_MinBias_Part1_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_mixedEventPFClustering_pseudoJetCandPtMin-0.0_2026-8-11_ultraFineCentBins.root";

const char *outPath =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/fakeJets/fakeJets_dPT_rawPt.root";

const int    N_pool  = 100;   // N_mixedEventsInPool
const double pT_min  = 20.;   // zero out bins below this threshold

const int sliceLo[5] = {0,  1,  3,  7, 11};
const int sliceHi[5] = {0,  2,  6, 10, 16};

void makeFakeJetFile_dPT_rawPt(){
  gSystem->Exec("mkdir -p /home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/fakeJets");

  TFile *fM = TFile::Open(fmixed_path);
  if(!fM || fM->IsZombie()){ printf("ERROR: cannot open mixed-event file\n"); return; }

  TFile *fOut = TFile::Open(outPath, "RECREATE");

  for(int ci = 0; ci <= 4; ci++){

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

    TH1D *hFJ_norm = hFJ_sum;
    TH1D *hRC_norm = hRC_sum;
    hFJ_norm->Scale(1. / N_events);
    hRC_norm->Scale(1. / (N_events * double(N_pool)));

    TH1D *hOut = (TH1D*) hFJ_norm->Clone(Form("h_fakeJets_C%d", ci));
    hOut->SetDirectory(nullptr);
    //hOut->Add(hRC_norm, -1.);

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
