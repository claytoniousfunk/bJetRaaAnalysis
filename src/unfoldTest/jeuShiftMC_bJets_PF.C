// Generator-level effect of the jet energy scale uncertainty on the unfolded MUON-TAGGED b-jet
// spectrum, in SIMULATION (2026-11-06, user: the b-jet analogue of the inclusive figure
// fig:jeuShift). Per system the reco projection of the chain's b response (the MC "measured"
// spectrum) is shifted by +/- the AK4PF JEC uncertainty with the analytic shift the chain applies
// to the data (headers/functions/jeuShift.h; pp Spring18_ppRef5TeV_V6, PbPb Autumn18_HI_V8, |eta|
// < 1.6 average), and the shifted and unshifted spectra are unfolded with the same response at the
// b-jet N_opt (rootFiles/unfoldTest/optimalIterations_bJets_PF_MCprior.root); the ratio is taken on
// the reported bins. Responses and unfolding as optimalIterations_bJets_PF_MCprior.C: the chain's,
// outliers clipped, JER-smeared x1.10 (jerSmear.h), MC prior, reco floor bRecoFloorPF().
// The JEU assigned to the result is the same shift applied to the MEASURED b spectrum
// (bJetsPerZSys_PF_PbPb.C, unfoldBJetSpectrum_caloJets_pp.C); this is its simulation counterpart.
// Ratio errors: none (same jets, analytic shift).
//
// Output: rootFiles/JEU/jeuShiftMC_bJets_PF.root, h_jeuRatio_{up,down}_<sys> (shifted / nominal)
// Usage, from src/unfoldTest/:
//   root -l -b -q -e 'gSystem->Load("/home/clayton/Programs/RooUnfold/build/libRooUnfold.so");
//                     gInterpreter->AddIncludePath("/home/clayton/Programs/RooUnfold/build");' jeuShiftMC_bJets_PF.C
#include "TFile.h"
#include "TSystem.h"
#include "TVectorD.h"
#include "../../headers/functions/jetCollectionInputs.h"
#include "../../headers/plotting/coarseCent.h"
#include "../../headers/plotting/plotStyle.h"
#include "../../headers/functions/dataDrivenPrior.h"
#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#endif
#include "../../headers/functions/bJetUnfoldExt.h"
#include "../../headers/functions/clipResponseOutliers.h"
#include "../../headers/functions/jerSmear.h"
#include "../../headers/functions/reverseUnfoldReweightCore.h"
#include "../../headers/functions/jeuShift.h"
// 7-bin axis (2026-11-16, the 7-bin nominal): with BJET_PT_7BINS the N_opt of the 7-bin closure
// (optimalIterations_bJets_PF_MCprior_csvV2_2D7bins.root) and the output jeuShiftMC_bJets_PF_7bins.root
inline const char* jsBinTag(bool closure){
#ifdef BJET_PT_7BINS
  return closure ? "_csvV2_2D7bins" : "_7bins";
#else
  return "";
#endif
}

void jeuShiftMC_bJets_PF(){
  setJetCollection("PFJets");
  const char *repo = "/home/clayton/Analysis/code/bJetRaaAnalysis";
  TFile *fR = TFile::Open(phRespPath), *fN = TFile::Open(Form("%s/rootFiles/unfoldTest/optimalIterations_bJets_PF_MCprior%s.root", repo, jsBinTag(true)));
  if(!fR || !fN){ printf("ERROR: an input did not open\n"); return; }
  TVectorD *nOpt = nullptr; fN->GetObject("optimalN", nOpt);
  if(!nOpt){ printf("ERROR: optimalN missing\n"); return; }
  const std::vector<double> dPP = jeuDelta(Form("%s/JetEnergyCorrections/Spring18_ppRef5TeV_V6_MC_Uncertainty_AK4PF.txt", repo));
  const std::vector<double> dAA = jeuDelta(Form("%s/JetEnergyCorrections/Autumn18_HI_V8_MC_Uncertainty_AK4PF.txt", repo));
  const double floorB = bRecoFloorPF();
  const std::vector<double> extE = bExtEdges(rvu::NEdge, rvu::edges);
  gSystem->mkdir(Form("%s/rootFiles/JEU", repo), true);
  TFile *out = TFile::Open(Form("%s/rootFiles/JEU/jeuShiftMC_bJets_PF%s.root", repo, jsBinTag(false)), "RECREATE");

  const char *sysTag[5] = {"pp", "C1", "C2", "C3", "C4"};
  for(int s = 0; s < 5; s++){
    TH2D *R = nullptr;
    if(s == 0){
      for(const std::string &fn : ppRespFiles){ TFile *f = TFile::Open(Form("%s%s", ppRespDir, fn.c_str())); TH2D *r = nullptr;
        if(f && !f->IsZombie()) f->GetObject(ppRespName, r);
        if(!r){ printf("ERROR: %s missing\n", ppRespName); return; }
        if(!R){ R = (TH2D*) r->Clone("Rpp"); R->SetDirectory(nullptr); } else R->Add(r); }
    } else {
      TH2D *r = nullptr; fR->GetObject(Form(phRespMuTagFmt, s), r);
      if(!r){ printf("ERROR: response missing for %s\n", sysTag[s]); return; }
      R = (TH2D*) r->Clone(Form("R_%s", sysTag[s])); R->SetDirectory(nullptr);
    }
    clipResponseOutliers(R, 100., false, sysTag[s]);
    { TH2D *sm = jerSmear::nominal(R, Form("%s_jerNom", R->GetName())); delete R; R = sm; }
    const int N = (int) (*nOpt)[s];
    const std::vector<double> &d = s == 0 ? dPP : dAA;
    TH1D *m0 = rvu::recoOnEdges(R, Form("m0_%s", sysTag[s]));
    TH1D *mU = jeuShiftSpectrum(m0, d, +1, Form("mU_%s", sysTag[s])), *mD = jeuShiftSpectrum(m0, d, -1, Form("mD_%s", sysTag[s]));
    TH2D *Re = rebin2DTo(R, extE, Form("Re_%s", sysTag[s]));
    TH1D *u0 = unfoldOnExt<RooUnfoldResponse, RooUnfoldBayes>(Re, m0, rvu::NEdge, rvu::edges, floorB, N, Form("u0_%s", sysTag[s]));
    TH1D *uU = unfoldOnExt<RooUnfoldResponse, RooUnfoldBayes>(Re, mU, rvu::NEdge, rvu::edges, floorB, N, Form("uU_%s", sysTag[s]));
    TH1D *uD = unfoldOnExt<RooUnfoldResponse, RooUnfoldBayes>(Re, mD, rvu::NEdge, rvu::edges, floorB, N, Form("uD_%s", sysTag[s]));
    std::vector<double> e;
    for(int b = 1; b < rvu::NEdge; b++) if(bReportedBin(rvu::edges[b-1], rvu::edges[b])){ if(e.empty()) e.push_back(rvu::edges[b-1]); e.push_back(rvu::edges[b]); }
    TH1D *hU = new TH1D(Form("h_jeuRatio_up_%s", sysTag[s]), Form("JEU up / nominal, N = %d", N), (int) e.size() - 1, e.data());
    TH1D *hD = new TH1D(Form("h_jeuRatio_down_%s", sysTag[s]), Form("JEU down / nominal, N = %d", N), (int) e.size() - 1, e.data());
    printf("%-3s N %d  reco shift up/down at 80-100: %+.1f / %+.1f %%  unfolded up / down (%%):", sysTag[s], N,
           100.*(mU->GetBinContent(mU->FindBin(90.))/m0->GetBinContent(m0->FindBin(90.)) - 1.), 100.*(mD->GetBinContent(mD->FindBin(90.))/m0->GetBinContent(m0->FindBin(90.)) - 1.));
    for(int i = 1; i <= hU->GetNbinsX(); i++){
      const int b = u0->FindBin(hU->GetBinCenter(i)); const double n = u0->GetBinContent(b);
      hU->SetBinContent(i, n > 0. ? uU->GetBinContent(b)/n : 0.); hD->SetBinContent(i, n > 0. ? uD->GetBinContent(b)/n : 0.);
      printf("  %.0f-%.0f %+.1f/%+.1f", hU->GetBinLowEdge(i), hU->GetBinLowEdge(i+1), 100.*(hU->GetBinContent(i) - 1.), 100.*(hD->GetBinContent(i) - 1.));
    }
    printf("\n");
    out->cd(); hU->Write(); hD->Write();
    delete R; delete Re; delete m0; delete mU; delete mD; delete u0; delete uU; delete uD;
  }
  { TVectorD v = *nOpt; out->cd(); v.Write("optimalN"); }
  out->Close();
}
