// Build rootFiles/fakeJets/fakeJets.root from the mixed-event FastJet anti-kT
// R=0.4 fake-jet spectrum for use as a fake-jet subtraction in calculateRAA.C.
//
// Normalization: fake jets per bin per event (same units as h_CX_jetMB after
// Scale(1/N_MinBias_events)), so direct subtraction is valid.
//
// Source: data-driven MinBias mixed-event scan (PbPb data, not HYDJET MC).
// The FastJet clustering runs once per real event on a sample of NCandidatesToSample
// PF candidates drawn from a pool of N_mixedEventsInPool=100 mixed events.
// Each event therefore contributes one clustering result to h_fastJetPt_JEC_C{i},
// so N_events = h_vz_C{i}->Integral() is the correct denominator.
//
// Coarse-bin mapping (matches calculateRAA.C):
//   C1 = 0-10%  : hiBin 0-20
//   C2 = 10-30% : hiBin 20-60
//   C3 = 30-50% : hiBin 60-100
//   C4 = 50-80% : hiBin 100-160

const char *fmixed_path =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/"
  "PbPb_MinBias_Part1_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_"
  "mixedEventPseudoJets_pfCandCS_pseudoJetCandPtMin-0.0_subleadingPFCandPtMin-15_"
  "2026-6-30.root";

const char *outPath =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/fakeJets/fakeJets.root";

const double pT_min = 20.;  // zero out fake jets below this threshold

void makeFakeJetFile(){
  gSystem->Exec("mkdir -p /home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/fakeJets");

  TFile *fM = TFile::Open(fmixed_path);
  if(!fM || fM->IsZombie()){ printf("ERROR: cannot open mixed-event file\n"); return; }

  TFile *fOut = TFile::Open(outPath, "RECREATE");

  for(int ci = 0; ci <= 4; ci++){
    TH1D *hFJ = nullptr, *hvz = nullptr;
    fM->GetObject(Form("h_fastJetPt_JEC_C%d", ci), hFJ);
    fM->GetObject(Form("h_vz_C%d", ci), hvz);
    if(!hFJ || !hvz){ printf("WARNING: missing histograms for C%d\n", ci); continue; }

    double N_events = hvz->Integral();
    TH1D *hOut = (TH1D*) hFJ->Clone(Form("h_fakeJets_C%d", ci));
    hOut->SetDirectory(nullptr);
    hOut->Scale(1. / N_events);

    for(int b = 1; b <= hOut->GetNbinsX(); b++){
      if(hOut->GetXaxis()->GetBinUpEdge(b) <= pT_min){
        hOut->SetBinContent(b, 0.);
        hOut->SetBinError(b, 0.);
      }
    }

    printf("h_fakeJets_C%d: N_events=%.0f  integral=%.4f (fake jets/event, pT > %.0f GeV)\n",
           ci, N_events, hOut->Integral(), pT_min);

    fOut->cd();
    hOut->Write();
  }

  fOut->Close();
  fM->Close();
  printf("\nWrote %s\n", outPath);
}
