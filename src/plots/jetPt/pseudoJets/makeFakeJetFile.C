// Build rootFiles/fakeJets/fakeJets.root from the CS mixed-event pseudo-jet
// distribution for use as a fake-jet subtraction in calculateRAA.C.
//
// Normalization: N_fake per bin per event (same units as h_CX_jetMB after
// Scale(1/N_MinBias_events)), so direct subtraction is valid.
//
// Formula per bin:
//   content[b] = hRaw[b] * N_indep / (N_events * N_cones)
// where N_events = total_raw_entries / N_cones  (unit-weight HYDJET sample).
//
// Coarse-bin mapping (matches calculateRAA.C):
//   C1 = 0-10%  : fine C1+C2
//   C2 = 10-30% : fine C3-C6
//   C3 = 30-50% : fine C7-C10
//   C4 = 50-80% : fine C11-C16

const char *fmixed_path =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/HYDJET/"
  "HYDJET_pThat-unweighted_mu12_pTmu-15to999_tight_hiBinShift-0_jetTrkMaxFilter_WDecayFilter_"
  "mixedEventPseudoJets_pfCandCS_pseudoJetCandPtMin-0.0_subleadingPFCandPtMin-15_"
  "2026-6-16_ultraFineCentBins.root";

const char *outPath =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/fakeJets/fakeJets.root";

const double pT_min  = 20.;  // only count fake jets above this threshold

const int    N_cones = 100;
const double R       = 0.4;
const double eta_range = 3.2;
const double N_indep = (eta_range * 2. * TMath::Pi()) / (TMath::Pi() * R * R); // 40

struct CoarseBin { int first, last; int raaIdx; };
const CoarseBin cb[] = {
  {1,  2,  1},   // C1 = 0-10%
  {3,  6,  2},   // C2 = 10-30%
  {7,  10, 3},   // C3 = 30-50%
  {11, 16, 4}    // C4 = 50-80%
};
const int NCB = 4;

// Sum fine-bin raw histograms and normalise to fake-jets per bin per event.
TH1D* makeCoarse(TFile *f, int firstFine, int lastFine, const char *name){
  TH1D *hSum = nullptr;
  double totalEntries = 0.;
  for(int iC = firstFine; iC <= lastFine; iC++){
    TH1D *hRaw; f->GetObject(Form("h_pseudoJetPt_C%d", iC), hRaw);
    if(!hRaw){ printf("WARNING: missing C%d\n", iC); continue; }
    totalEntries += hRaw->GetEntries();
    if(!hSum) hSum = (TH1D*) hRaw->Clone(name);
    else      hSum->Add(hRaw);
  }
  if(!hSum || totalEntries == 0) return nullptr;
  double N_events = totalEntries / double(N_cones);
  double norm = N_indep / (N_events * double(N_cones));
  hSum->Scale(norm);
  // zero out bins below pT_min
  for(int b = 1; b <= hSum->GetNbinsX(); b++){
    if(hSum->GetXaxis()->GetBinUpEdge(b) <= pT_min){
      hSum->SetBinContent(b, 0.);
      hSum->SetBinError(b, 0.);
    }
  }
  hSum->SetDirectory(nullptr);
  return hSum;
}

void makeFakeJetFile(){
  gSystem->Exec("mkdir -p /home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/fakeJets");

  TFile *fM = TFile::Open(fmixed_path);
  if(!fM || fM->IsZombie()){ printf("ERROR: cannot open mixed-event file\n"); return; }

  printf("N_indep = %.0f\n", N_indep);

  TFile *fOut = TFile::Open(outPath, "RECREATE");

  for(int i = 0; i < NCB; i++){
    TH1D *h = makeCoarse(fM, cb[i].first, cb[i].last,
                         Form("h_fakeJets_C%d", cb[i].raaIdx));
    if(!h){ printf("ERROR: failed to build C%d\n", cb[i].raaIdx); continue; }
    double integral = h->Integral();
    printf("h_fakeJets_C%d: integral=%.4f (fake jets/event, pT > %.0f GeV)\n",
           cb[i].raaIdx, integral, pT_min);
    fOut->cd();
    h->Write();
  }

  fOut->Close();
  fM->Close();
  printf("\nWrote %s\n", outPath);
}
