// Build rootFiles/fakeJets/fakeJets.root from the mixed-event FastJet anti-kT
// R=0.4 fake-jet spectrum for use as a fake-jet subtraction in calculateRAA.C.
//
// Normalization: fake jets per bin per event (same units as h_CX_jetMB after
// Scale(1/N_MinBias_events)), so direct subtraction is valid.
//
// The sum of event weights per fine centrality bin is recovered as
//   sumW = h_pseudoJetPt_C{i}->Integral() / N_genPJ
// which equals the number of events (w=1 for pThat-unweighted HYDJET).
//
// Coarse-bin mapping (matches calculateRAA.C):
//   C1 = 0-10%  : fine C1+C2   (hiBin 0-20)
//   C2 = 10-30% : fine C3-C6   (hiBin 20-60)
//   C3 = 30-50% : fine C7-C10  (hiBin 60-100)
//   C4 = 50-80% : fine C11-C16 (hiBin 100-160)

const char *fmixed_path =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/HYDJET/"
  "HYDJET_pThat-unweighted_mu12_pTmu-15to999_tight_hiBinShift-0_jetTrkMaxFilter_WDecayFilter_"
  "mixedEventPseudoJets_pfCandCS_pseudoJetCandPtMin-0.0_subleadingPFCandPtMin-15_"
  "2026-6-23_ultraFineCentBins.root";

const char *outPath =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/fakeJets/fakeJets.root";

const double pT_min = 20.;  // zero out fake jets below this threshold
const int    N_genPJ = 100; // N_generatedPseudoJets (used to recover sumW)

struct CoarseBin { int first, last; int raaIdx; };
const CoarseBin cb[] = {
  {1,  2,  1},   // C1 = 0-10%
  {3,  6,  2},   // C2 = 10-30%
  {7,  10, 3},   // C3 = 30-50%
  {11, 16, 4}    // C4 = 50-80%
};
const int NCB = 4;

// Sum fine-bin FastJet histograms and normalise to fake-jets per bin per event.
TH1D* makeCoarse(TFile *f, int firstFine, int lastFine, const char *name){
  TH1D *hFJsum = nullptr;
  double sumW   = 0.;
  for(int iC = firstFine; iC <= lastFine; iC++){
    TH1D *hFJ = nullptr, *hPJ = nullptr;
    f->GetObject(Form("h_fastJetPt_C%d", iC), hFJ);
    f->GetObject(Form("h_pseudoJetPt_C%d", iC), hPJ);
    if(!hFJ || !hPJ){ printf("WARNING: missing histograms for C%d\n", iC); continue; }
    sumW += hPJ->Integral() / double(N_genPJ);
    if(!hFJsum) hFJsum = (TH1D*) hFJ->Clone(name);
    else        hFJsum->Add(hFJ);
  }
  if(!hFJsum || sumW == 0.) return nullptr;
  hFJsum->Scale(1. / sumW);
  for(int b = 1; b <= hFJsum->GetNbinsX(); b++){
    if(hFJsum->GetXaxis()->GetBinUpEdge(b) <= pT_min){
      hFJsum->SetBinContent(b, 0.);
      hFJsum->SetBinError(b, 0.);
    }
  }
  hFJsum->SetDirectory(nullptr);
  return hFJsum;
}

void makeFakeJetFile(){
  gSystem->Exec("mkdir -p /home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/fakeJets");

  TFile *fM = TFile::Open(fmixed_path);
  if(!fM || fM->IsZombie()){ printf("ERROR: cannot open mixed-event file\n"); return; }

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
