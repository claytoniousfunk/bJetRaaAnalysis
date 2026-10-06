// general ROOT/C includes
#include <iostream>
#include "TFile.h"
#include "TRandom.h"
#include "TTree.h"
#include "TSystem.h"    // gSystem, for the output-directory check
#include "TRandom2.h"   // TRandom2 is instantiated below; only TRandom.h was included
#include "TVector2.h"   // Phi_mpi_pi, for the calo-PF flavor match
#include "TH1F.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TRandom.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TProfile2D.h"
#include <TF1.h>
#include "assert.h"
#include <fstream>
#include "TMath.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TMath.h"
#include <TNtuple.h>
#include "TChain.h"
#include <TString.h>
#include <TLatex.h>
#include <TCut.h>
#include "TDatime.h"
#include <vector>
#include "TCanvas.h"
#include <dirent.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
// #if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
// #include "RooUnfoldResponse.h"
// #include "RooUnfoldBayes.h"
// #endif

// event map
#include "../../../eventMap/eventMap.h"
// jet corrector
#include "../../../JetEnergyCorrections/JetCorrector.h"
// jet uncertainty
#include "../../../JetEnergyCorrections/JetUncertainty.h"
// general analysis variables
// analysis constants + centrality scheme (classes, getCentBin, CENT_SCHEME_SUFFIX):
// switch CENT_SCHEME in config_centrality.h, as for the data scans. Under
// CENT_ULTRAFINE the classes are 5% slices out to 90%, so a 50-90% response can
// be built by summing slices 11-18.
#include "../../../headers/config/config_centrality.h"

// Nominal class (1 = 0-10, 2 = 10-30, 3 = 30-50, 4 = 50-80%) of a hiBin, for the
// fits that exist only for those four classes (pThat correlation, JER, b-jet
// spectrum reweight). Under CENT_NOMINAL it equals the histogram index; under
// CENT_ULTRAFINE the histogram index is a 5% slice and must not pick the fit.
// 80-90% takes the 50-80% fits, as PbPb_scan.C does for its trigger fits.
inline int nominalCentClass(int hiBin){
  if(hiBin < 20) return 1;
  if(hiBin < 60) return 2;
  if(hiBin < 100) return 3;
  return 4;
}
// vz-fit parameters
//#include "../../../headers/fitParameters/vzFitParams_PH_mu5.h"
//#include "../../../headers/fitParameters/vzFitParams_PH_mu7.h"
//#include "../../../headers/fitParameters/vzFitParams_PH_mu12.h"
#include "../../../headers/fitParameters/vzFitParams_PYTHIAHYDJET.h"
// hiBin-fit parameters
//#include "../../../headers/fitParameters/hiBinFitParams_mu5.h"
//#include "../../../headers/fitParameters/hiBinFitParams_mu7.h"
//#include "../../../headers/fitParameters/hiBinFitParams_mu12.h"
#include "../../../headers/fitParameters/hiBinFitParams_PYTHIAHYDJET.h"
// jetPt-fit parameters
//#include "../../../headers/fitParameters/jetPtFitParams_PYTHIA_mu5.h"
//#include "../../../headers/fitParameters/jetPtFitParams_PYTHIA_mu7.h"
#include "../../../headers/fitParameters/jetPtFitParams_PYTHIA_mu12.h"
// JERCorrection params
#include "../../../headers/fitParameters/JERCorrectionParams_PYTHIA_mu12.h"
TF1 *fitFxn_hiBin, *fitFxn_vz, *fitFxn_jetPt, *fitFxn_PYTHIA_JESb, *fitFxn_PYTHIA_JERCorrection;

// vz-fit function
#include "../../../headers/fitFunctions/fitFxn_vz_PH.h"
// #include "../../../headers/fitFunctions/fitFxn_vz_PH_mu7.h"
// #include "../../../headers/fitFunctions/fitFxn_vz_PH_mu12.h"
// hiBin-fit function
#include "../../../headers/fitFunctions/fitFxn_hiBin.h"
// #include "../../../headers/fitFunctions/fitFxn_hiBin_mu7.h"
// #include "../../../headers/fitFunctions/fitFxn_hiBin_mu12.h"
// jetPt-fit function
#include "../../../headers/fitFunctions/fitFxn_jetPt.h"
// BJetSpectraReweightToData fxn
#include "../../../headers/fitParameters/BJetSpectraReweightToDataFitParams_PYTHIAHYDJET.h"
TF1 *fitFxn_PYTHIAHYDJET_BJetSpectraReweightToData_C4;
TF1 *fitFxn_PYTHIAHYDJET_BJetSpectraReweightToData_C3;
TF1 *fitFxn_PYTHIAHYDJET_BJetSpectraReweightToData_C2;
TF1 *fitFxn_PYTHIAHYDJET_BJetSpectraReweightToData_C1;
#include "../../../headers/fitFunctions/fitFxn_PYTHIAHYDJET_BJetSpectraReweightToData.h"
// JER-correction function
#include "../../../headers/fitFunctions/fitFxn_PYTHIA_JERCorrection.h"
// pThat correlation
#include "../../../headers/fitFunctions/fitFxn_PYTHIAHYDJET_pThatCorrelation.h"
#include "../../../headers/fitFunctions/fitFxn_PYTHIAHYDJET_pThatCorrelation_caloJets.h"
// eta-phi mask function
#include "../../../headers/functions/etaPhiMask.h"
// getDr function
#include "../../../headers/functions/getDr.h"

// dR within which a reco jet counts as matched to a gen jet, for the
// reco-jet-indexed diagnostics (h_recoJetPt_all / _matchedDr / _unmatchedDr).
// Half the jet radius is the usual choice for R = 0.4. This is deliberately
// separate from the refpt == genjetpt test used for the response matrix, so the
// two can be compared rather than conflated.
const double recoGenMatchDr = 0.2;

// getJetPtBin function
#include "../../../headers/functions/getJetPtBin.h"
// getPtRel function
#include "../../../headers/functions/getPtRel.h"
// isQualityMuon_hybridSoft function
#include "../../../headers/functions/isQualityMuon_hybridSoft.h"
// isQualityMuon_tight function
#include "../../../headers/functions/isQualityMuon_tight.h"
// isWDecayMuon function
#include "../../../headers/functions/isWDecayMuon.h"
// triggerIsOn function
#include "../../../headers/functions/triggerIsOn.h"
// pthat filter function
#include "../../../headers/functions/passesLeadingGenJetPthatFilter.h"
// JetTrkMax filter function
#include "../../../headers/functions/jet_filter/passesJetTrkMaxFilter.h"
// print introduction
#include "../../../headers/introductions/printIntroduction_PYTHIAHYDJET_scan_V3p7.h"
// analysis config
#include "../../../headers/config/config_PYTHIAHYDJET.h"
// read config
#include "../../../headers/config/readConfig.h"
// remove HYDJET jets function
#include "../../../headers/functions/jet_filter/remove_HYDJET_jet.h"
// dataset naming functions
#include "../../../headers/functions/getDatasetName/getDatasetName.h"
#include "../../../headers/functions/getInputFileName/getInputFileName.h"
#include "../../../headers/functions/configureOutputDatasetName/configureOutputDatasetName_PYTHIAHYDJET_response.h"


// ---- calo-jet flavor from the matched PF jet --------------------------------
// Each calo jet takes the parton flavor of the nearest PF jet within
// caloPFMatchDR (scan_calo_pf_match.h, shared with the analysis scans), instead
// of its own refparton_flavorForB. ON, so the response's b class is the same
// label the PYTHIAHYDJET_scan.C ptRel templates and tag frequency use (_PFflavor
// there and here). The PbPb akCs4PF trees carry the flavor as
// matchedPartonFlavor, not jtPartonFlavor; the shared header reads either. Calo
// jets with no PF jet that close keep refparton_flavorForB.
//
// The PF collection is akCs4PF (constituent subtracted), the calo one akPu4Calo
// (pileup subtracted), so the same jet can sit further apart than in pp;
// caloPFMatchDR is looser, and the match rate printed at the end of the job
// should be checked.
#include "../scan_calo_pf_match.h"
// calo-jet match jet ID for PF jets (config_PYTHIAHYDJET.h: requireCaloJetMatch)
#include "../scan_calo_jet_match.h"
bool   caloFlavorFromPFMatch = true;
double caloPFMatchDR         = 0.3;   // looser than pp: different subtraction

// ---- muon tag ---------------------------------------------------------------
// true: a jet is muon-tagged with the ANALYSIS reco-muon tag of
// PYTHIAHYDJET_scan.C -- a muon-tree muon passing passesRecoMuonCuts (tight ID,
// muPtCut < pT < muPtMaxCut, |eta| < 2), not a W-decay muon (doWDecayFilter),
// within epsilon_mm of the reco jet axis, each muon used by one jet only.
// false: the jet branch's own muon (mupt/mueta > muPtCut, |eta| < 2), with no
// ID, no W veto and no one-muon-per-jet rule -- what every response scan before
// 2026-10-04 used. Tagged _anaMuTag in the output name.
bool   useAnalysisMuonTag    = true;
#include "../scan_mc_muon_tag.h"

// one place that decides a jet's flavor, so the two fill sites cannot drift
inline int recoJetFlavorFor(bool useCalo, int idx, eventMap *em)
{
  const int fallback = em->refparton_flavorForB[idx];
  if(useCalo && caloFlavorFromPFMatch)
    return caloFlavorByPFMatch(em->jeteta[idx], em->jetphi[idx], caloPFMatchDR, fallback);
  return fallback;
}

// Distorted-truth response matrices for the unfolding closure test
// (unfoldClosureTest.C's "tilt" mode). w(genPt) = (genPt/pTref)^expo, applied
// per event at fill time -- NOT a post-facto reweight of the already-filled
// response histogram -- so P(reco|gen) for each event carries its own weight
// rather than a value averaged over its gen bin. Same formula and pTref as
// unfoldClosureTest.C::distortionWeight, frozen outside [40,500] so the
// extrapolation cannot drive a gen bin to zero or blow it up.
const int    kNDistTilt = 2;
const double kDistTiltExpo[kNDistTilt] = {0.3, 0.5};
const double kDistTiltPtRef = 150.;
const char  *kDistTiltFlavorTag[2] = {"allJets", "bJets"};
const char  *kDistTiltMuTagSuffix[2] = {"", "_muTagged"};
double distTiltWeight(double genPt, double expo){
  double x = genPt;
  if(x <  40.) x =  40.;
  if(x > 500.) x = 500.;
  return TMath::Power(x/kDistTiltPtRef, expo);
}

void PYTHIAHYDJET_scan_response(int group = 1){

  std::cout << "setting pthat cut to 15...\n";
  pthatcut = 15.;
  //std::cout << "turning off hiBin & vz reweights...\n";
  doHiBinReweight = true;
  doVzReweight = true;
  std::cout << "turning off removeHYDJETjet...\n";
  doRemoveHYDJETjet = false;
  std::cout << "turning on pThat correlation filter...\n";
  doPThatCorrelationFilter = true;
  
  // TString inputDataset = "";
  // TString inputFileName = "";

  // //inputDataset = "/eos/user/c/cbennett/skims/output_skim_PH_DiJet_pTjet-5_withJetTriggers/";
  // inputDataset = "/eos/user/c/cbennett/skims/output_skim_PH_DiJet_pTjet-5_withGenNeutrino_withRefPt/";
  // inputFileName = "PYTHIAHYDJET_DiJet_skim_output";
  // TString input = Form("%s%s_%i.root",inputDataset.Data(),inputFileName.Data(),group);

  std::string inputFileList = "";
  inputFileList = "../../../fileNames/fileNames_PH_DiJet_withCaloAndFlowJets_fix2.txt";

  if(group == 0){ cout << "INPUTFILELIST=" << inputFileList << endl; return; } // query mode for condor submit scripts

  std::ifstream instr(inputFileList.c_str(), std::ifstream::in);
  if(!instr.is_open()){
    cout << "filelist not found!! Exiting..." << endl;
    return;
  }
  std::string filename;
  Int_t ifile = 0;

  while(instr>>filename){

    ifile++;

    if(ifile != group) continue;

    std::string input = filename.c_str();

    std::cout << "input dataset = " << input << std::endl;

    TString outputBaseDir = "/eos/cms/store/group/phys_heavyions/cbennett/scanningOutput/";

    TString outputDatasetName = "";

    if(requireCaloJetMatch && useCaloJetsOverride){
      std::cout << "\033[1;31m requireCaloJetMatch is a PF-jet ID; it cannot be combined with useCaloJetsOverride \033[0m" << std::endl;
      return;
    }
    outputDatasetName = configureOutputDatasetName(generator,
						   doDiJetSample,
						   doMuJetSample,
						   doBJetSample,
						   doDiJetSample_batch1,
						   doDiJetSample_batch2,
						   doDiJetSample_batch3,
						   doDiJetSample_batch4,
						   doDiJetSample_batch5,
						   doDiJetSample_batch6,
						   doDiJetSample_batch7,
						   doDiJetSample_batch8,
						   doDiJetSample_batch9,
						   doDiJetSample_batch10,
						   doDiJetSample_batch11,
						   doDiJetSample_batch12,
						   doDiJetSample_batch13,
						   doDiJetSample_batch14,
						   doDiJetSample_batch15,
						   pthatcut,
						   doPThatWeight,
						   doVzReweight,
						   doHiBinReweight,
						   doJetPtReweight,
						   doGenJetPthatFilter,
						   doLeadingXjetDumpFilter,
						   doXdumpReweight,
						   doJetTrkMaxFilter,
						   doRemoveHYDJETjet,
						   doEtaPhiMask,
						   doDRReweight,
						   doWeightCut,
						   doBJetSpectraReweightToData,
						   doHadronPtRelReweight,
						   doBJetEnergyShift,
						   doBJetNeutrinoEnergyShift,						 
						   doJERCorrection,
						   doJESCorrection,						 
						   apply_JER_smear,
						   apply_JEU_shift_up,
						   apply_JEU_shift_down,
						   hiBinShift,
						   applyJet60Trigger,
						   applyJet80Trigger,
						   muPtCut,
						   doPThatCorrelationFilter,
						   useCaloJetsOverride,
						   caloFlavorFromPFMatch,
						   useManualJEC,
						   onlyEvenEvents,
						   onlyOddEvents,
						   onlyMuTaggedJets,
						   useAnalysisMuonTag,
						   requireCaloJetMatch);


    TString suffixEdit = CENT_SCHEME_SUFFIX;   // "_ultraFineCentBins" etc., "" for nominal
    TString output = Form("%s%s%s/PYTHIAHYDJET_scan_output_%i.root",outputBaseDir.Data(),outputDatasetName.Data(),suffixEdit.Data(),group);
    //TString output = Form("%s%s_muTaggedJetsNoTrigger/PYTHIAHYDJET_scan_output_%i.root",outputBaseDir.Data(),outputDatasetName.Data(),group);
    //TString output = Form("%s%s_evenEvents/PYTHIAHYDJET_scan_output_%i.root",outputBaseDir.Data(),outputDatasetName.Data(),group);
    //TString output = Form("%s%s_oddEvents/PYTHIAHYDJET_scan_output_%i.root",outputBaseDir.Data(),outputDatasetName.Data(),group);
    //TString output = Form("%s%s_ultraFineCentBins/PYTHIAHYDJET_scan_output_%i.root",outputBaseDir.Data(),outputDatasetName.Data(),group);

    std::cout << "output dataset = " << output << std::endl;

    if(gSystem->AccessPathName(Form("%s%s%s",outputBaseDir.Data(),outputDatasetName.Data(),suffixEdit.Data()))){
      std::cout << "\033[1;31m Output directory not found: \033[0m " << Form("%s%s%s",outputBaseDir.Data(),outputDatasetName.Data(),suffixEdit.Data()) << std::endl;
      return;
    }

    readConfig();

    // JET ENERGY CORRECTIONS
    vector<string> Files;

    if(useCaloJetsOverride){
      Files.push_back("../../../JetEnergyCorrections/Autumn18_HI_V8_MC_L2Relative_AK4Calo.txt");
    }
    else{
      Files.push_back("../../../JetEnergyCorrections/Autumn18_HI_V8_MC_L2Relative_AK4PF.txt");
    }

    JetCorrector JEC(Files);

    string JEU_path = "";
    
    if(useCaloJetsOverride){
      JEU_path = "../../../JetEnergyCorrections/Autumn18_HI_V8_MC_Uncertainty_AK4Calo.txt";
    }
    else{
      JEU_path = "../../../JetEnergyCorrections/Autumn18_HI_V8_MC_Uncertainty_AK4PF.txt";
    }
    
    JetUncertainty JEU(JEU_path.c_str());

    // WEIGHT FUNCTIONS

    //  initialize histograms
    // JETS

    TH2D *h_matchedRecoJetPt_genJetPt[NCentralityIndices][7];
    //TH2D *h_matchedRecoJetPtWithCut_genJetPt[NCentralityIndices][7];
    TH2D *h_matchedRecoJetPt_genJetPt_var[NCentralityIndices][7];
    TH2D *h_matchedRecoJetPtOverGenJetPt_genJetPt[NCentralityIndices][7];
    TH2D *h_matchedRecoJetPtOverGenJetPt_genJetEta[NCentralityIndices][7];
    TH2D *h_inclGenJetPt_flavor[NCentralityIndices];
    TH2D *h_inclGenJetPt_inclGenMuonTag_flavor[NCentralityIndices];
    TH2D *h_inclGenJetPt_inclRecoMuonTag_flavor[NCentralityIndices];
    TH2D *h_leadingRecoJetPtOverPThat_pThat[NCentralityIndices];
    TH1D *h_unmatchedRecoJetPt[NCentralityIndices][7];

    // Reco-jet-indexed matching diagnostics. Unlike h_matchedRecoJetPt_genJetPt
    // (filled per GEN jet) and h_unmatchedRecoJetPt (filled per reco jet with a
    // different weight), these are all filled in ONE loop over reco jets with the
    // SAME weight and one dR condition, so
    //     all = matchedDr + unmatchedDr
    // holds bin by bin and fractions of reco jets are well defined.
    TH1D *h_recoJetPt_all[NCentralityIndices];
    TH1D *h_recoJetPt_matchedDr[NCentralityIndices];
    TH1D *h_recoJetPt_unmatchedDr[NCentralityIndices];
    // calo-match diagnostics (scan_calo_jet_match.h): dR to the closest calo jet
    // (9.9 = no calo jet in the event, overflow) and that jet's raw pT when it is
    // inside caloMatchDr, vs PF jet pT -- as in PbPb_scan.C. With
    // requireCaloJetMatch the cut is already applied and dR stops at caloMatchDr.
    TH2D *h_caloJetDr_recoJetPt[NCentralityIndices];
    TH2D *h_caloJetRawPt_recoJetPt[NCentralityIndices];
    // For every reco jet, the nearest gen jet regardless of whether it passes the
    // cut: separates pure combinatorial jets (no gen jet anywhere near) from soft
    // gen jets promoted upward by the underlying event (small dR, genPt << recoPt).
    TH2D *h_recoPt_dRnearestGen[NCentralityIndices];
    TH2D *h_recoPt_nearestGenPt[NCentralityIndices];
    TH1D *h_unmatchedGenJetPt[NCentralityIndices];

    // distorted-truth response + misses for unfoldClosureTest.C's "tilt" mode --
    // allJets and bJets, fixed and variable pT binning. muTag=1 additionally
    // requires hasRecoJetMuon, so a muon-tagged-only closure test can be read
    // out of this (onlyMuTaggedJets=false) production instead of needing a
    // second scan run with onlyMuTaggedJets=true.
    TH2D *h_matchedRecoJetPt_genJetPt_distTilt[NCentralityIndices][2][2][kNDistTilt];       // [.][flavor][muTag][tilt]
    TH2D *h_matchedRecoJetPt_genJetPt_var_distTilt[NCentralityIndices][2][2][kNDistTilt];   // same, variable bins
    TH1D *h_unmatchedGenJetPt_distTilt[NCentralityIndices][kNDistTilt];   // misses aren't flavor/muTag-split upstream either

    // plain (undistorted) muon-tagged response -- allJets/bJets, fixed/var bins.
    // The undistorted-response side of unfoldClosureTest.C's closure test (e.g.
    // the odd half) needs this baseline alongside the distTilt histograms above.
    TH2D *h_matchedRecoJetPt_genJetPt_muTagged[NCentralityIndices][2];       // [.][flavor], fixed bins
    TH2D *h_matchedRecoJetPt_genJetPt_var_muTagged[NCentralityIndices][2];   // [.][flavor], variable bins
    // RooUnfoldResponse response_C4(NPtBins,ptMin,ptMax,"response_C4","response_C4");
    // RooUnfoldResponse response_C3(NPtBins,ptMin,ptMax,"response_C3","response_C3");
    // RooUnfoldResponse response_C2(NPtBins,ptMin,ptMax,"response_C2","response_C2");
    // RooUnfoldResponse response_C1(NPtBins,ptMin,ptMax,"response_C1","response_C1");
    // RooUnfoldResponse response_C0(NPtBins,ptMin,ptMax,"response_C0","response_C0");



    // Define histograms
    const int N1 = 10;
    double ptAxis1[N1] = {60,70,80,90,100,120,150,200,300,500};
    const int N2 = 32;
    double ptAxis2[N2] = {0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,500};



    for(int i = 0; i < NCentralityIndices; i++){

      if(i==0) {
	h_unmatchedGenJetPt[i] = new TH1D(Form("h_unmatchedGenJetPt_C%i",i),Form("unmatchedGenJetPt, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),NPtBins,ptMin,ptMax);
	h_unmatchedRecoJetPt[i][0] = new TH1D(Form("h_unmatchedRecoJetPt_allJets_C%i",i),Form("unmatchedRecoJetPt, allJets, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),NPtBins,ptMin,ptMax);
	h_unmatchedRecoJetPt[i][1] = new TH1D(Form("h_unmatchedRecoJetPt_bJets_C%i",i),Form("unmatchedRecoJetPt, bJets, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),NPtBins,ptMin,ptMax);
	h_unmatchedRecoJetPt[i][2] = new TH1D(Form("h_unmatchedRecoJetPt_cJets_C%i",i),Form("unmatchedRecoJetPt, cJets, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),NPtBins,ptMin,ptMax);
	h_unmatchedRecoJetPt[i][3] = new TH1D(Form("h_unmatchedRecoJetPt_udJets_C%i",i),Form("unmatchedRecoJetPt, udJets, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),NPtBins,ptMin,ptMax);
	h_unmatchedRecoJetPt[i][4] = new TH1D(Form("h_unmatchedRecoJetPt_sJets_C%i",i),Form("unmatchedRecoJetPt, sJets, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),NPtBins,ptMin,ptMax);
	h_unmatchedRecoJetPt[i][5] = new TH1D(Form("h_unmatchedRecoJetPt_gJets_C%i",i),Form("unmatchedRecoJetPt, gJets, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),NPtBins,ptMin,ptMax);
	h_unmatchedRecoJetPt[i][6] = new TH1D(Form("h_unmatchedRecoJetPt_xJets_C%i",i),Form("unmatchedRecoJetPt, xJets, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),NPtBins,ptMin,ptMax);
	h_recoJetPt_all[i]         = new TH1D(Form("h_recoJetPt_all_C%i",i),Form("all reco jets, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),NPtBins,ptMin,ptMax);
	h_recoJetPt_matchedDr[i]   = new TH1D(Form("h_recoJetPt_matchedDr_C%i",i),Form("reco jets with a gen jet within dR, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),NPtBins,ptMin,ptMax);
	h_recoJetPt_unmatchedDr[i] = new TH1D(Form("h_recoJetPt_unmatchedDr_C%i",i),Form("reco jets with no gen jet within dR, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),NPtBins,ptMin,ptMax);
	h_recoPt_dRnearestGen[i]   = new TH2D(Form("h_recoPt_dRnearestGen_C%i",i),Form("dR to nearest gen jet vs reco p_{T}, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),NPtBins,ptMin,ptMax,100,0,1.0);
	h_recoPt_nearestGenPt[i]   = new TH2D(Form("h_recoPt_nearestGenPt_C%i",i),Form("nearest gen jet p_{T} vs reco p_{T}, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),NPtBins,ptMin,ptMax,NPtBins,ptMin,ptMax);
	h_matchedRecoJetPt_genJetPt[i][0] = new TH2D(Form("h_matchedRecoJetPt_genJetPt_allJets_C%i",i),Form("genJetPt vs. matchedRecoJetPt, allJets, hiBin %i - %i", centEdges[0]-10,centEdges[NCentralityIndices-1]-10),NPtBins,ptMin,ptMax,NPtBins,ptMin,ptMax);
	h_matchedRecoJetPt_genJetPt[i][1] = new TH2D(Form("h_matchedRecoJetPt_genJetPt_bJets_C%i",i),Form("genJetPt vs. matchedRecoJetPt, bJets, hiBin %i - %i", centEdges[0]-10,centEdges[NCentralityIndices-1]-10),NPtBins,ptMin,ptMax,NPtBins,ptMin,ptMax);
	h_matchedRecoJetPt_genJetPt[i][2] = new TH2D(Form("h_matchedRecoJetPt_genJetPt_cJets_C%i",i),Form("genJetPt vs. matchedRecoJetPt, cJets, hiBin %i - %i", centEdges[0]-10,centEdges[NCentralityIndices-1]-10),NPtBins,ptMin,ptMax,NPtBins,ptMin,ptMax);
	h_matchedRecoJetPt_genJetPt[i][3] = new TH2D(Form("h_matchedRecoJetPt_genJetPt_udJets_C%i",i),Form("genJetPt vs. matchedRecoJetPt, udJets, hiBin %i - %i", centEdges[0]-10,centEdges[NCentralityIndices-1]-10),NPtBins,ptMin,ptMax,NPtBins,ptMin,ptMax);
	h_matchedRecoJetPt_genJetPt[i][4] = new TH2D(Form("h_matchedRecoJetPt_genJetPt_sJets_C%i",i),Form("genJetPt vs. matchedRecoJetPt, sJets, hiBin %i - %i", centEdges[0]-10,centEdges[NCentralityIndices-1]-10),NPtBins,ptMin,ptMax,NPtBins,ptMin,ptMax);
	h_matchedRecoJetPt_genJetPt[i][5] = new TH2D(Form("h_matchedRecoJetPt_genJetPt_gJets_C%i",i),Form("genJetPt vs. matchedRecoJetPt, gJets, hiBin %i - %i", centEdges[0]-10,centEdges[NCentralityIndices-1]-10),NPtBins,ptMin,ptMax,NPtBins,ptMin,ptMax);
	h_matchedRecoJetPt_genJetPt[i][6] = new TH2D(Form("h_matchedRecoJetPt_genJetPt_xJets_C%i",i),Form("genJetPt vs. matchedRecoJetPt, xJets, hiBin %i - %i", centEdges[0]-10,centEdges[NCentralityIndices-1]-10),NPtBins,ptMin,ptMax,NPtBins,ptMin,ptMax);
	h_inclGenJetPt_flavor[i] = new TH2D(Form("h_inclGenJetPt_flavor_C%i",i),Form("JetFlavorID vs incl. gen p_{T}^{jet}, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),NPtBins,ptMin,ptMax,27,-5,22);
	h_inclGenJetPt_inclGenMuonTag_flavor[i] = new TH2D(Form("h_inclGenJetPt_inclGenMuonTag_flavor_C%i",i),Form("JetFlavorID vs incl. gen p_{T}^{jet}, tagged with incl. gen muon, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),NPtBins,ptMin,ptMax,27,-5,22);
	h_inclGenJetPt_inclRecoMuonTag_flavor[i] = new TH2D(Form("h_inclGenJetPt_inclRecoMuonTag_flavor_C%i",i),Form("JetFlavorID vs incl. gen p_{T}^{jet}, tagged with incl. reco muon, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),NPtBins,ptMin,ptMax,27,-5,22);
	h_leadingRecoJetPtOverPThat_pThat[i] = new TH2D(Form("h_leadingRecoJetPtOverPThat_pThat_C%i",i),Form("(leadingRecoJetPt / pThat) vs. pThat, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),500,0,5,500,0,500);
      }
      else {
	h_unmatchedGenJetPt[i] = new TH1D(Form("h_unmatchedGenJetPt_C%i",i),Form("unmatchedGenJetPt, hiBin %i - %i",centEdges[i-1],centEdges[i]),NPtBins,ptMin,ptMax);
	h_unmatchedRecoJetPt[i][0] = new TH1D(Form("h_unmatchedRecoJetPt_allJets_C%i",i),Form("unmatchedRecoJetPt, allJets, hiBin %i - %i",centEdges[i-1],centEdges[i]),NPtBins,ptMin,ptMax);
	h_unmatchedRecoJetPt[i][1] = new TH1D(Form("h_unmatchedRecoJetPt_bJets_C%i",i),Form("unmatchedRecoJetPt, bJets, hiBin %i - %i",centEdges[i-1],centEdges[i]),NPtBins,ptMin,ptMax);
	h_unmatchedRecoJetPt[i][2] = new TH1D(Form("h_unmatchedRecoJetPt_cJets_C%i",i),Form("unmatchedRecoJetPt, cJets, hiBin %i - %i",centEdges[i-1],centEdges[i]),NPtBins,ptMin,ptMax);
	h_unmatchedRecoJetPt[i][3] = new TH1D(Form("h_unmatchedRecoJetPt_udJets_C%i",i),Form("unmatchedRecoJetPt, udJets, hiBin %i - %i",centEdges[i-1],centEdges[i]),NPtBins,ptMin,ptMax);
	h_unmatchedRecoJetPt[i][4] = new TH1D(Form("h_unmatchedRecoJetPt_sJets_C%i",i),Form("unmatchedRecoJetPt, sJets, hiBin %i - %i",centEdges[i-1],centEdges[i]),NPtBins,ptMin,ptMax);
	h_unmatchedRecoJetPt[i][5] = new TH1D(Form("h_unmatchedRecoJetPt_gJets_C%i",i),Form("unmatchedRecoJetPt, gJets, hiBin %i - %i",centEdges[i-1],centEdges[i]),NPtBins,ptMin,ptMax);
	h_unmatchedRecoJetPt[i][6] = new TH1D(Form("h_unmatchedRecoJetPt_xJets_C%i",i),Form("unmatchedRecoJetPt, xJets, hiBin %i - %i",centEdges[i-1],centEdges[i]),NPtBins,ptMin,ptMax);
	h_recoJetPt_all[i]         = new TH1D(Form("h_recoJetPt_all_C%i",i),Form("all reco jets, hiBin %i - %i",centEdges[i-1],centEdges[i]),NPtBins,ptMin,ptMax);
	h_recoJetPt_matchedDr[i]   = new TH1D(Form("h_recoJetPt_matchedDr_C%i",i),Form("reco jets with a gen jet within dR, hiBin %i - %i",centEdges[i-1],centEdges[i]),NPtBins,ptMin,ptMax);
	h_recoJetPt_unmatchedDr[i] = new TH1D(Form("h_recoJetPt_unmatchedDr_C%i",i),Form("reco jets with no gen jet within dR, hiBin %i - %i",centEdges[i-1],centEdges[i]),NPtBins,ptMin,ptMax);
	h_recoPt_dRnearestGen[i]   = new TH2D(Form("h_recoPt_dRnearestGen_C%i",i),Form("dR to nearest gen jet vs reco p_{T}, hiBin %i - %i",centEdges[i-1],centEdges[i]),NPtBins,ptMin,ptMax,100,0,1.0);
	h_recoPt_nearestGenPt[i]   = new TH2D(Form("h_recoPt_nearestGenPt_C%i",i),Form("nearest gen jet p_{T} vs reco p_{T}, hiBin %i - %i",centEdges[i-1],centEdges[i]),NPtBins,ptMin,ptMax,NPtBins,ptMin,ptMax);
	h_matchedRecoJetPt_genJetPt[i][0] = new TH2D(Form("h_matchedRecoJetPt_genJetPt_allJets_C%i",i),Form("genJetPt vs. matchedRecoJetPt, allJets, hiBin %i - %i", centEdges[i-1]-10,centEdges[i]-10),NPtBins,ptMin,ptMax,NPtBins,ptMin,ptMax);
	h_matchedRecoJetPt_genJetPt[i][1] = new TH2D(Form("h_matchedRecoJetPt_genJetPt_bJets_C%i",i),Form("genJetPt vs. matchedRecoJetPt, bJets, hiBin %i - %i", centEdges[i-1]-10,centEdges[i]-10),NPtBins,ptMin,ptMax,NPtBins,ptMin,ptMax);
	h_matchedRecoJetPt_genJetPt[i][2] = new TH2D(Form("h_matchedRecoJetPt_genJetPt_cJets_C%i",i),Form("genJetPt vs. matchedRecoJetPt, cJets, hiBin %i - %i", centEdges[i-1]-10,centEdges[i]-10),NPtBins,ptMin,ptMax,NPtBins,ptMin,ptMax);
	h_matchedRecoJetPt_genJetPt[i][3] = new TH2D(Form("h_matchedRecoJetPt_genJetPt_udJets_C%i",i),Form("genJetPt vs. matchedRecoJetPt, udJets, hiBin %i - %i", centEdges[i-1]-10,centEdges[i]-10),NPtBins,ptMin,ptMax,NPtBins,ptMin,ptMax);
	h_matchedRecoJetPt_genJetPt[i][4] = new TH2D(Form("h_matchedRecoJetPt_genJetPt_sJets_C%i",i),Form("genJetPt vs. matchedRecoJetPt, sJets, hiBin %i - %i", centEdges[i-1]-10,centEdges[i]-10),NPtBins,ptMin,ptMax,NPtBins,ptMin,ptMax);
	h_matchedRecoJetPt_genJetPt[i][5] = new TH2D(Form("h_matchedRecoJetPt_genJetPt_gJets_C%i",i),Form("genJetPt vs. matchedRecoJetPt, gJets, hiBin %i - %i", centEdges[i-1]-10,centEdges[i]-10),NPtBins,ptMin,ptMax,NPtBins,ptMin,ptMax);
	h_matchedRecoJetPt_genJetPt[i][6] = new TH2D(Form("h_matchedRecoJetPt_genJetPt_xJets_C%i",i),Form("genJetPt vs. matchedRecoJetPt, xJets, hiBin %i - %i", centEdges[i-1]-10,centEdges[i]-10),NPtBins,ptMin,ptMax,NPtBins,ptMin,ptMax);
	h_inclGenJetPt_flavor[i] = new TH2D(Form("h_inclGenJetPt_flavor_C%i",i),Form("JetFlavorID vs incl. gen p_{T}^{jet}, hiBin %i - %i",centEdges[i-1],centEdges[i]),NPtBins,ptMin,ptMax,27,-5,22);
	h_inclGenJetPt_inclGenMuonTag_flavor[i] = new TH2D(Form("h_inclGenJetPt_inclGenMuonTag_flavor_C%i",i),Form("JetFlavorID vs incl. gen p_{T}^{jet}, tagged with incl. gen muon, hiBin %i - %i",centEdges[i-1],centEdges[i]),NPtBins,ptMin,ptMax,27,-5,22);
	h_inclGenJetPt_inclRecoMuonTag_flavor[i] = new TH2D(Form("h_inclGenJetPt_inclRecoMuonTag_flavor_C%i",i),Form("JetFlavorID vs incl. gen p_{T}^{jet}, tagged with incl. reco muon, hiBin %i - %i",centEdges[i-1],centEdges[i]),NPtBins,ptMin,ptMax,27,-5,22);
	h_leadingRecoJetPtOverPThat_pThat[i] = new TH2D(Form("h_leadingRecoJetPtOverPThat_pThat_C%i",i),Form("(leadingRecoJetPt / pThat) vs. pThat, hiBin %i - %i",centEdges[i-1],centEdges[i]),500,0,5,100,0,500);
      }

      h_unmatchedGenJetPt[i]->Sumw2();
      h_unmatchedRecoJetPt[i][0]->Sumw2();
      h_unmatchedRecoJetPt[i][1]->Sumw2();
      h_unmatchedRecoJetPt[i][2]->Sumw2();
      h_unmatchedRecoJetPt[i][3]->Sumw2();
      h_unmatchedRecoJetPt[i][4]->Sumw2();
      h_unmatchedRecoJetPt[i][5]->Sumw2();
      h_unmatchedRecoJetPt[i][6]->Sumw2();
      {
	TString centTitle = (i==0) ? Form("hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1])
	                           : Form("hiBin %i - %i",centEdges[i-1],centEdges[i]);
	h_caloJetDr_recoJetPt[i] = new TH2D(Form("h_caloJetDr_recoJetPt_C%i",i),Form("#DeltaR(PF jet, closest calo jet) vs p_{T}^{jet}, %s; p_{T}^{jet} [GeV]; #DeltaR",centTitle.Data()),NPtBins,ptMin,ptMax,100,0.0,1.0);
	h_caloJetRawPt_recoJetPt[i] = new TH2D(Form("h_caloJetRawPt_recoJetPt_C%i",i),Form("closest calo raw p_{T} vs PF p_{T}^{jet}, #DeltaR < %.2f, %s; p_{T}^{jet} [GeV]; calo raw p_{T} [GeV]",caloMatchDr,centTitle.Data()),NPtBins,ptMin,ptMax,NPtBins,ptMin,ptMax);
      }
      h_recoJetPt_all[i]->Sumw2();
      h_recoJetPt_matchedDr[i]->Sumw2();
      h_recoJetPt_unmatchedDr[i]->Sumw2();
      h_caloJetDr_recoJetPt[i]->Sumw2();
      h_caloJetRawPt_recoJetPt[i]->Sumw2();
      h_recoPt_dRnearestGen[i]->Sumw2();
      h_recoPt_nearestGenPt[i]->Sumw2();
      h_matchedRecoJetPt_genJetPt[i][0]->Sumw2();
      h_matchedRecoJetPt_genJetPt[i][1]->Sumw2();
      h_matchedRecoJetPt_genJetPt[i][2]->Sumw2();
      h_matchedRecoJetPt_genJetPt[i][3]->Sumw2();
      h_matchedRecoJetPt_genJetPt[i][4]->Sumw2();
      h_matchedRecoJetPt_genJetPt[i][5]->Sumw2();
      h_matchedRecoJetPt_genJetPt[i][6]->Sumw2();
      h_inclGenJetPt_flavor[i]->Sumw2();
      h_inclGenJetPt_inclGenMuonTag_flavor[i]->Sumw2();
      h_inclGenJetPt_inclRecoMuonTag_flavor[i]->Sumw2();
      h_leadingRecoJetPtOverPThat_pThat[i]->Sumw2();

      for(int f = 0; f < 2; f++){
	for(int m = 0; m < 2; m++){
	  for(int t = 0; t < kNDistTilt; t++){
	    h_matchedRecoJetPt_genJetPt_distTilt[i][f][m][t] = new TH2D(Form("h_matchedRecoJetPt_genJetPt_%s%s_C%i_distTilt%.1f",kDistTiltFlavorTag[f],kDistTiltMuTagSuffix[m],i,kDistTiltExpo[t]),Form("genJetPt vs. matchedRecoJetPt, %s%s, C%i, gen reweighted w#propto p_{T,gen}^{%.1f}",kDistTiltFlavorTag[f],kDistTiltMuTagSuffix[m],i,kDistTiltExpo[t]),NPtBins,ptMin,ptMax,NPtBins,ptMin,ptMax);
	    h_matchedRecoJetPt_genJetPt_distTilt[i][f][m][t]->Sumw2();
	    h_matchedRecoJetPt_genJetPt_var_distTilt[i][f][m][t] = new TH2D(Form("h_matchedRecoJetPt_genJetPt_var_%s%s_C%i_distTilt%.1f",kDistTiltFlavorTag[f],kDistTiltMuTagSuffix[m],i,kDistTiltExpo[t]),Form("genJetPt vs. matchedRecoJetPt, var bins, %s%s, C%i, gen reweighted w#propto p_{T,gen}^{%.1f}",kDistTiltFlavorTag[f],kDistTiltMuTagSuffix[m],i,kDistTiltExpo[t]),N1-1,ptAxis1,N1-1,ptAxis1);
	    h_matchedRecoJetPt_genJetPt_var_distTilt[i][f][m][t]->Sumw2();
	  }
	}
	h_matchedRecoJetPt_genJetPt_muTagged[i][f] = new TH2D(Form("h_matchedRecoJetPt_genJetPt_%s_C%i_muTagged",kDistTiltFlavorTag[f],i),Form("genJetPt vs. matchedRecoJetPt, %s, C%i, muon-tagged",kDistTiltFlavorTag[f],i),NPtBins,ptMin,ptMax,NPtBins,ptMin,ptMax);
	h_matchedRecoJetPt_genJetPt_muTagged[i][f]->Sumw2();
	h_matchedRecoJetPt_genJetPt_var_muTagged[i][f] = new TH2D(Form("h_matchedRecoJetPt_genJetPt_var_%s_C%i_muTagged",kDistTiltFlavorTag[f],i),Form("genJetPt vs. matchedRecoJetPt, var bins, %s, C%i, muon-tagged",kDistTiltFlavorTag[f],i),N1-1,ptAxis1,N1-1,ptAxis1);
	h_matchedRecoJetPt_genJetPt_var_muTagged[i][f]->Sumw2();
      }
      for(int t = 0; t < kNDistTilt; t++){
	h_unmatchedGenJetPt_distTilt[i][t] = new TH1D(Form("h_unmatchedGenJetPt_C%i_distTilt%.1f",i,kDistTiltExpo[t]),Form("unmatchedGenJetPt, C%i, gen reweighted w#propto p_{T,gen}^{%.1f}",i,kDistTiltExpo[t]),NPtBins,ptMin,ptMax);
	h_unmatchedGenJetPt_distTilt[i][t]->Sumw2();
      }


      if(i==0) {

	h_matchedRecoJetPt_genJetPt_var[i][0] = new TH2D(Form("h_matchedRecoJetPt_genJetPt_var_allJets_C%i",i),Form("genJetPt vs. matchedRecoJetPt, hiBin %i - %i, allJets, var bins", centEdges[0]-10,centEdges[NCentralityIndices-1]-10),N1-1,ptAxis1,N1-1,ptAxis1);
	h_matchedRecoJetPt_genJetPt_var[i][1] = new TH2D(Form("h_matchedRecoJetPt_genJetPt_var_bJets_C%i",i),Form("genJetPt vs. matchedRecoJetPt, hiBin %i - %i, bJets, var bins", centEdges[0]-10,centEdges[NCentralityIndices-1]-10),N1-1,ptAxis1,N1-1,ptAxis1);
	h_matchedRecoJetPt_genJetPt_var[i][2] = new TH2D(Form("h_matchedRecoJetPt_genJetPt_var_cJets_C%i",i),Form("genJetPt vs. matchedRecoJetPt, hiBin %i - %i, cJets, var bins", centEdges[0]-10,centEdges[NCentralityIndices-1]-10),N1-1,ptAxis1,N1-1,ptAxis1);
	h_matchedRecoJetPt_genJetPt_var[i][3] = new TH2D(Form("h_matchedRecoJetPt_genJetPt_var_udJets_C%i",i),Form("genJetPt vs. matchedRecoJetPt, hiBin %i - %i, udJets, var bins", centEdges[0]-10,centEdges[NCentralityIndices-1]-10),N1-1,ptAxis1,N1-1,ptAxis1);
	h_matchedRecoJetPt_genJetPt_var[i][4] = new TH2D(Form("h_matchedRecoJetPt_genJetPt_var_sJets_C%i",i),Form("genJetPt vs. matchedRecoJetPt, hiBin %i - %i, sJets, var bins", centEdges[0]-10,centEdges[NCentralityIndices-1]-10),N1-1,ptAxis1,N1-1,ptAxis1);
	h_matchedRecoJetPt_genJetPt_var[i][5] = new TH2D(Form("h_matchedRecoJetPt_genJetPt_var_gJets_C%i",i),Form("genJetPt vs. matchedRecoJetPt, hiBin %i - %i, gJets, var bins", centEdges[0]-10,centEdges[NCentralityIndices-1]-10),N1-1,ptAxis1,N1-1,ptAxis1);
	h_matchedRecoJetPt_genJetPt_var[i][6] = new TH2D(Form("h_matchedRecoJetPt_genJetPt_var_xJets_C%i",i),Form("genJetPt vs. matchedRecoJetPt, hiBin %i - %i, xJets, var bins", centEdges[0]-10,centEdges[NCentralityIndices-1]-10),N1-1,ptAxis1,N1-1,ptAxis1);
      
	h_matchedRecoJetPtOverGenJetPt_genJetPt[i][0] = new TH2D(Form("h_matchedRecoJetPtOverGenJetPt_genJetPt_allJets_C%i",i),Form("matchedRecoJetPt/genJetPt vs genJetPt, all flavors, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),500,0,5,NPtBins,ptMin,ptMax);
	h_matchedRecoJetPtOverGenJetPt_genJetPt[i][1] = new TH2D(Form("h_matchedRecoJetPtOverGenJetPt_genJetPt_bJets_C%i",i),Form("matchedRecoJetPt/genJetPt vs genJetPt, bJets, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),500,0,5,NPtBins,ptMin,ptMax);
	h_matchedRecoJetPtOverGenJetPt_genJetPt[i][2] = new TH2D(Form("h_matchedRecoJetPtOverGenJetPt_genJetPt_cJets_C%i",i),Form("matchedRecoJetPt/genJetPt vs genJetPt, cJets, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),500,0,5,NPtBins,ptMin,ptMax);
	h_matchedRecoJetPtOverGenJetPt_genJetPt[i][3] = new TH2D(Form("h_matchedRecoJetPtOverGenJetPt_genJetPt_udJets_C%i",i),Form("matchedRecoJetPt/genJetPt vs genJetPt, udJets, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),500,0,5,NPtBins,ptMin,ptMax);
	h_matchedRecoJetPtOverGenJetPt_genJetPt[i][4] = new TH2D(Form("h_matchedRecoJetPtOverGenJetPt_genJetPt_sJets_C%i",i),Form("matchedRecoJetPt/genJetPt vs genJetPt, sJets, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),500,0,5,NPtBins,ptMin,ptMax);
	h_matchedRecoJetPtOverGenJetPt_genJetPt[i][5] = new TH2D(Form("h_matchedRecoJetPtOverGenJetPt_genJetPt_gJets_C%i",i),Form("matchedRecoJetPt/genJetPt vs genJetPt, gJets, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),500,0,5,NPtBins,ptMin,ptMax);
	h_matchedRecoJetPtOverGenJetPt_genJetPt[i][6] = new TH2D(Form("h_matchedRecoJetPtOverGenJetPt_genJetPt_xJets_C%i",i),Form("matchedRecoJetPt/genJetPt vs genJetPt, xJets, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),500,0,5,NPtBins,ptMin,ptMax);

	h_matchedRecoJetPtOverGenJetPt_genJetEta[i][0] = new TH2D(Form("h_matchedRecoJetPtOverGenJetPt_genJetEta_allJets_C%i",i),Form("matchedRecoJetPt/genJetPt vs genJetEta, all flavors, hiBin %i - %i, p_{T}^{jet} > 50 GeV",centEdges[0],centEdges[NCentralityIndices-1]),500,0,5,NEtaBins,etaMin,etaMax);
	h_matchedRecoJetPtOverGenJetPt_genJetEta[i][1] = new TH2D(Form("h_matchedRecoJetPtOverGenJetPt_genJetEta_bJets_C%i",i),Form("matchedRecoJetPt/genJetPt vs genJetEta, bJets, hiBin %i - %i, p_{T}^{jet} > 50 GeV",centEdges[0],centEdges[NCentralityIndices-1]),500,0,5,NEtaBins,etaMin,etaMax);
	h_matchedRecoJetPtOverGenJetPt_genJetEta[i][2] = new TH2D(Form("h_matchedRecoJetPtOverGenJetPt_genJetEta_cJets_C%i",i),Form("matchedRecoJetPt/genJetPt vs genJetEta, cJets, hiBin %i - %i, p_{T}^{jet} > 50 GeV",centEdges[0],centEdges[NCentralityIndices-1]),500,0,5,NEtaBins,etaMin,etaMax);
	h_matchedRecoJetPtOverGenJetPt_genJetEta[i][3] = new TH2D(Form("h_matchedRecoJetPtOverGenJetPt_genJetEta_udJets_C%i",i),Form("matchedRecoJetPt/genJetPt vs genJetEta, udJets, hiBin %i - %i, p_{T}^{jet} > 50 GeV",centEdges[0],centEdges[NCentralityIndices-1]),500,0,5,NEtaBins,etaMin,etaMax);
	h_matchedRecoJetPtOverGenJetPt_genJetEta[i][4] = new TH2D(Form("h_matchedRecoJetPtOverGenJetPt_genJetEta_sJets_C%i",i),Form("matchedRecoJetPt/genJetPt vs genJetEta, sJets, hiBin %i - %i, p_{T}^{jet} > 50 GeV",centEdges[0],centEdges[NCentralityIndices-1]),500,0,5,NEtaBins,etaMin,etaMax);
	h_matchedRecoJetPtOverGenJetPt_genJetEta[i][5] = new TH2D(Form("h_matchedRecoJetPtOverGenJetPt_genJetEta_gJets_C%i",i),Form("matchedRecoJetPt/genJetPt vs genJetEta, gJets, hiBin %i - %i, p_{T}^{jet} > 50 GeV",centEdges[0],centEdges[NCentralityIndices-1]),500,0,5,NEtaBins,etaMin,etaMax);
	h_matchedRecoJetPtOverGenJetPt_genJetEta[i][6] = new TH2D(Form("h_matchedRecoJetPtOverGenJetPt_genJetEta_xJets_C%i",i),Form("matchedRecoJetPt/genJetPt vs genJetEta, xJets, hiBin %i - %i, p_{T}^{jet} > 50 GeV",centEdges[0],centEdges[NCentralityIndices-1]),500,0,5,NEtaBins,etaMin,etaMax);      
      }
      else{

	h_matchedRecoJetPt_genJetPt_var[i][0] = new TH2D(Form("h_matchedRecoJetPt_genJetPt_var_allJets_C%i",i),Form("genJetPt vs. matchedRecoJetPt, hiBin %i - %i, allJets, var bins", centEdges[i-1]-10,centEdges[i]-10),N1-1,ptAxis1,N1-1,ptAxis1) ;
	h_matchedRecoJetPt_genJetPt_var[i][1] = new TH2D(Form("h_matchedRecoJetPt_genJetPt_var_bJets_C%i",i),Form("genJetPt vs. matchedRecoJetPt, hiBin %i - %i, bJets, var bins", centEdges[i-1]-10,centEdges[i]-10),N1-1,ptAxis1,N1-1,ptAxis1) ;
	h_matchedRecoJetPt_genJetPt_var[i][2] = new TH2D(Form("h_matchedRecoJetPt_genJetPt_var_cJets_C%i",i),Form("genJetPt vs. matchedRecoJetPt, hiBin %i - %i, cJets, var bins", centEdges[i-1]-10,centEdges[i]-10),N1-1,ptAxis1,N1-1,ptAxis1) ;
	h_matchedRecoJetPt_genJetPt_var[i][3] = new TH2D(Form("h_matchedRecoJetPt_genJetPt_var_udJets_C%i",i),Form("genJetPt vs. matchedRecoJetPt, hiBin %i - %i, udJets, var bins", centEdges[i-1]-10,centEdges[i]-10),N1-1,ptAxis1,N1-1,ptAxis1) ;
	h_matchedRecoJetPt_genJetPt_var[i][4] = new TH2D(Form("h_matchedRecoJetPt_genJetPt_var_sJets_C%i",i),Form("genJetPt vs. matchedRecoJetPt, hiBin %i - %i, sJets, var bins", centEdges[i-1]-10,centEdges[i]-10),N1-1,ptAxis1,N1-1,ptAxis1) ;
	h_matchedRecoJetPt_genJetPt_var[i][5] = new TH2D(Form("h_matchedRecoJetPt_genJetPt_var_gJets_C%i",i),Form("genJetPt vs. matchedRecoJetPt, hiBin %i - %i, gJets, var bins", centEdges[i-1]-10,centEdges[i]-10),N1-1,ptAxis1,N1-1,ptAxis1) ;
	h_matchedRecoJetPt_genJetPt_var[i][6] = new TH2D(Form("h_matchedRecoJetPt_genJetPt_var_xJets_C%i",i),Form("genJetPt vs. matchedRecoJetPt, hiBin %i - %i, xJets, var bins", centEdges[i-1]-10,centEdges[i]-10),N1-1,ptAxis1,N1-1,ptAxis1) ;

      
	h_matchedRecoJetPtOverGenJetPt_genJetPt[i][0] = new TH2D(Form("h_matchedRecoJetPtOverGenJetPt_genJetPt_allJets_C%i",i),Form("matchedRecoJetPt/genJetPt vs genJetPt, all flavors, hiBin %i - %i",centEdges[i-1],centEdges[i]),500,0,5,NPtBins,ptMin,ptMax);
	h_matchedRecoJetPtOverGenJetPt_genJetPt[i][1] = new TH2D(Form("h_matchedRecoJetPtOverGenJetPt_genJetPt_bJets_C%i",i),Form("matchedRecoJetPt/genJetPt vs genJetPt, bJets, hiBin %i - %i",centEdges[i-1],centEdges[i]),500,0,5,NPtBins,ptMin,ptMax);
	h_matchedRecoJetPtOverGenJetPt_genJetPt[i][2] = new TH2D(Form("h_matchedRecoJetPtOverGenJetPt_genJetPt_cJets_C%i",i),Form("matchedRecoJetPt/genJetPt vs genJetPt, cJets, hiBin %i - %i",centEdges[i-1],centEdges[i]),500,0,5,NPtBins,ptMin,ptMax);
	h_matchedRecoJetPtOverGenJetPt_genJetPt[i][3] = new TH2D(Form("h_matchedRecoJetPtOverGenJetPt_genJetPt_udJets_C%i",i),Form("matchedRecoJetPt/genJetPt vs genJetPt, udJets, hiBin %i - %i",centEdges[i-1],centEdges[i]),500,0,5,NPtBins,ptMin,ptMax);
	h_matchedRecoJetPtOverGenJetPt_genJetPt[i][4] = new TH2D(Form("h_matchedRecoJetPtOverGenJetPt_genJetPt_sJets_C%i",i),Form("matchedRecoJetPt/genJetPt vs genJetPt, sJets, hiBin %i - %i",centEdges[i-1],centEdges[i]),500,0,5,NPtBins,ptMin,ptMax);
	h_matchedRecoJetPtOverGenJetPt_genJetPt[i][5] = new TH2D(Form("h_matchedRecoJetPtOverGenJetPt_genJetPt_gJets_C%i",i),Form("matchedRecoJetPt/genJetPt vs genJetPt, gJets, hiBin %i - %i",centEdges[i-1],centEdges[i]),500,0,5,NPtBins,ptMin,ptMax);
	h_matchedRecoJetPtOverGenJetPt_genJetPt[i][6] = new TH2D(Form("h_matchedRecoJetPtOverGenJetPt_genJetPt_xJets_C%i",i),Form("matchedRecoJetPt/genJetPt vs genJetPt, xJets, hiBin %i - %i",centEdges[i-1],centEdges[i]),500,0,5,NPtBins,ptMin,ptMax); 

	h_matchedRecoJetPtOverGenJetPt_genJetEta[i][0] = new TH2D(Form("h_matchedRecoJetPtOverGenJetPt_genJetEta_allJets_C%i",i),Form("matchedRecoJetPt/genJetPt vs genJetEta, all flavors, hiBin %i - %i, p_{T}^{jet} > 50 GeV",centEdges[i-1],centEdges[i]),500,0,5,NEtaBins,etaMin,etaMax);
	h_matchedRecoJetPtOverGenJetPt_genJetEta[i][1] = new TH2D(Form("h_matchedRecoJetPtOverGenJetPt_genJetEta_bJets_C%i",i),Form("matchedRecoJetPt/genJetPt vs genJetEta, bJets, hiBin %i - %i, p_{T}^{jet} > 50 GeV",centEdges[i-1],centEdges[i]),500,0,5,NEtaBins,etaMin,etaMax);
	h_matchedRecoJetPtOverGenJetPt_genJetEta[i][2] = new TH2D(Form("h_matchedRecoJetPtOverGenJetPt_genJetEta_cJets_C%i",i),Form("matchedRecoJetPt/genJetPt vs genJetEta, cJets, hiBin %i - %i, p_{T}^{jet} > 50 GeV",centEdges[i-1],centEdges[i]),500,0,5,NEtaBins,etaMin,etaMax);
	h_matchedRecoJetPtOverGenJetPt_genJetEta[i][3] = new TH2D(Form("h_matchedRecoJetPtOverGenJetPt_genJetEta_udJets_C%i",i),Form("matchedRecoJetPt/genJetPt vs genJetEta, udJets, hiBin %i - %i, p_{T}^{jet} > 50 GeV",centEdges[i-1],centEdges[i]),500,0,5,NEtaBins,etaMin,etaMax);
	h_matchedRecoJetPtOverGenJetPt_genJetEta[i][4] = new TH2D(Form("h_matchedRecoJetPtOverGenJetPt_genJetEta_sJets_C%i",i),Form("matchedRecoJetPt/genJetPt vs genJetEta, sJets, hiBin %i - %i, p_{T}^{jet} > 50 GeV",centEdges[i-1],centEdges[i]),500,0,5,NEtaBins,etaMin,etaMax);
	h_matchedRecoJetPtOverGenJetPt_genJetEta[i][5] = new TH2D(Form("h_matchedRecoJetPtOverGenJetPt_genJetEta_gJets_C%i",i),Form("matchedRecoJetPt/genJetPt vs genJetEta, gJets, hiBin %i - %i, p_{T}^{jet} > 50 GeV",centEdges[i-1],centEdges[i]),500,0,5,NEtaBins,etaMin,etaMax);
	h_matchedRecoJetPtOverGenJetPt_genJetEta[i][6] = new TH2D(Form("h_matchedRecoJetPtOverGenJetPt_genJetEta_xJets_C%i",i),Form("matchedRecoJetPt/genJetPt vs genJetEta, xJets, hiBin %i - %i, p_{T}^{jet} > 50 GeV",centEdges[i-1],centEdges[i]),500,0,5,NEtaBins,etaMin,etaMax); 
      }

      h_matchedRecoJetPt_genJetPt_var[i][0]->Sumw2();
      h_matchedRecoJetPt_genJetPt_var[i][1]->Sumw2();
      h_matchedRecoJetPt_genJetPt_var[i][2]->Sumw2();
      h_matchedRecoJetPt_genJetPt_var[i][3]->Sumw2();
      h_matchedRecoJetPt_genJetPt_var[i][4]->Sumw2();
      h_matchedRecoJetPt_genJetPt_var[i][5]->Sumw2();
      h_matchedRecoJetPt_genJetPt_var[i][6]->Sumw2();

      h_matchedRecoJetPtOverGenJetPt_genJetPt[i][0]->Sumw2();
      h_matchedRecoJetPtOverGenJetPt_genJetPt[i][1]->Sumw2();
      h_matchedRecoJetPtOverGenJetPt_genJetPt[i][2]->Sumw2();
      h_matchedRecoJetPtOverGenJetPt_genJetPt[i][3]->Sumw2();
      h_matchedRecoJetPtOverGenJetPt_genJetPt[i][4]->Sumw2();
      h_matchedRecoJetPtOverGenJetPt_genJetPt[i][5]->Sumw2();
      h_matchedRecoJetPtOverGenJetPt_genJetPt[i][6]->Sumw2();

      h_matchedRecoJetPtOverGenJetPt_genJetEta[i][0]->Sumw2();
      h_matchedRecoJetPtOverGenJetPt_genJetEta[i][1]->Sumw2();
      h_matchedRecoJetPtOverGenJetPt_genJetEta[i][2]->Sumw2();
      h_matchedRecoJetPtOverGenJetPt_genJetEta[i][3]->Sumw2();
      h_matchedRecoJetPtOverGenJetPt_genJetEta[i][4]->Sumw2();
      h_matchedRecoJetPtOverGenJetPt_genJetEta[i][5]->Sumw2();
      h_matchedRecoJetPtOverGenJetPt_genJetEta[i][6]->Sumw2();

    }



    TFile *f = TFile::Open(input.c_str());
    cout << "	File opened!" << endl;
    auto em = new eventMap(f);
    em->isMC = isMC_status;
    em->AASetup = AASetup_status;
    cout << "	Initializing variables ... " << endl;
    em->init();
    cout << "	Loading jet..." << endl;
    if(useCaloJetsOverride){
      em->loadJet("akPu4CaloJetAnalyzer/t");
    }
    else{
      em->loadJet("akCs4PFJetAnalyzer/t");
    }

    // second, un-friended read of the PF jet tree, for calo flavor by dR match.
    // It cannot go through eventMap: loadJet() attaches the jet tree as a FRIEND
    // of evtTree, so a second one would collide on jtpt, jteta and the rest.
    if(useCaloJetsOverride && caloFlavorFromPFMatch){
      attachCaloPFMatchTree(f, "akCs4PFJetAnalyzer/t", true, false);
      if(g_pfHasFlav) cout << "	calo flavor from the PF jet within dR < " << caloPFMatchDR << "\n";
    }
    else g_pfTree = nullptr;

    // calo jets for the PF-jet calo-match ID (scan_calo_jet_match.h)
    g_caloTree = nullptr;
    if(!useCaloJetsOverride){
      if(!attachCaloJetMatchTree(f, em->evtTree->GetEntries()) && requireCaloJetMatch){
        std::cout << "\033[1;31m requireCaloJetMatch is set but the calo jets cannot be read; aborting \033[0m" << std::endl;
        return;
      }
    }
    if(useAnalysisMuonTag){
      cout << "	Loading muon..." << endl;
      em->loadMuon("ggHiNtuplizerGED/EventTree");
    }
    cout << "Loading muon triggers..." << endl;
    em->loadHLT("hltanalysis/HltTree");
    cout << "	Loading gen particles..." << endl;
    em->loadGenParticle("HiGenParticleAna/hi");
    cout << "	Variables initilized!" << endl << endl ;
    int NEvents = em->evtTree->GetEntries();
    cout << "	Number of events = " << NEvents << endl;


    // define event filters
    em->regEventFilter(NeventFilters, eventFilters);

    TRandom *randomGenerator = new TRandom2();

    // jet-energy resolution fit function
    TF1 *JER_fxn[NCentralityIndices];
    // pp
    JER_fxn[0] = new TF1("JER_fxn_pp","sqrt([0]*[0] + [1]*[1]/x + [2]*[2]/(x*x))",30,500);
    JER_fxn[0]->SetParameter(0,0.0640995);
    JER_fxn[0]->SetParameter(1,0.917851);
    JER_fxn[0]->SetParameter(2,-0.00211695);
    // C4
    JER_fxn[4] = new TF1("JER_fxn_C4","sqrt([0]*[0] + [1]*[1]/x + [2]*[2]/(x*x))",30,500);
    JER_fxn[4]->SetParameter(0,0.0606347);
    JER_fxn[4]->SetParameter(1,1.08087);
    JER_fxn[4]->SetParameter(2,-0.374138);
    // C3
    JER_fxn[3] = new TF1("JER_fxn_C3","sqrt([0]*[0] + [1]*[1]/x + [2]*[2]/(x*x))",30,500);
    JER_fxn[3]->SetParameter(0,0.0576787);
    JER_fxn[3]->SetParameter(1,1.1762);
    JER_fxn[3]->SetParameter(2,-5.67268);
    // C2
    JER_fxn[2] = new TF1("JER_fxn_C2","sqrt([0]*[0] + [1]*[1]/x + [2]*[2]/(x*x))",30,500);
    JER_fxn[2]->SetParameter(0,0.0547384);
    JER_fxn[2]->SetParameter(1,1.306);
    JER_fxn[2]->SetParameter(2,-11.1249);
    // C1
    JER_fxn[1] = new TF1("JER_fxn_C1","sqrt([0]*[0] + [1]*[1]/x + [2]*[2]/(x*x))",30,500);
    JER_fxn[1]->SetParameter(0,0.0598659);
    JER_fxn[1]->SetParameter(1,1.30631);
    JER_fxn[1]->SetParameter(2,-16.893);

    // define vz & hiBin reweighting functions
    if(fillMu5){
      loadFitFxn_vz_mu5();
      loadFitFxn_hiBin_mu5();
    }
    else if(fillMu7){
      loadFitFxn_vz_mu7();
      loadFitFxn_hiBin_mu7();
    }
    else if(fillMu12){
      loadFitFxn_vz_mu12();
      loadFitFxn_hiBin_mu12();
    }
    else{};
    // load JER correction fit fxn
    loadFitFxn_PYTHIA_JERCorrection();
    loadFitFxn_PYTHIAHYDJET_pThatCorrelation();
    loadFitFxn_PYTHIAHYDJET_pThatCorrelation_caloJets();
    loadFitFxn_PYTHIAHYDJET_BJetSpectraReweightToData();

    TFile *f_neutrino_energy_fraction_map = TFile::Open("/eos/cms/store/group/phys_heavyions/cbennett/maps/neutrino_energy_fraction_map.root");
    TH2D *neutrino_energy_fraction_map;
    TH1D *neutrino_energy_fraction_map_proj;
    f_neutrino_energy_fraction_map->GetObject("neutrino_energy_fraction_map",neutrino_energy_fraction_map);

    TFile *f_neutrino_energy_map = TFile::Open("/eos/cms/store/group/phys_heavyions/cbennett/maps/neutrino_energy_map.root");
    TH2D *neutrino_energy_map;
    TH1D *neutrino_energy_map_proj;
    f_neutrino_energy_map->GetObject("neutrino_energy_map",neutrino_energy_map);

    TFile *f_neutrino_tag_fraction = TFile::Open("/eos/cms/store/group/phys_heavyions/cbennett/maps/neutrino_tag_fraction.root");
    TH1D *neutrino_tag_fraction;
    f_neutrino_tag_fraction->GetObject("neutrino_tag_fraction",neutrino_tag_fraction);


  
    // event loop
    int evi_frac = 0;
    for(int evi = 0; evi < NEvents ; evi++){


      if(evi==0) cout << "Processing events..." << endl;

      // // only take even events
      if(onlyEvenEvents){
	if(evi % 2 == 1) continue;
      }

      // only take odd events
      if(onlyOddEvents){
	if(evi % 2 == 0) continue;
      }
    
      em->getEvent(evi);

      // step the un-friended PF tree to the same event
      if(g_pfTree) g_pfTree->GetEntry(evi);
      // and the calo tree for the calo-match ID
      loadCaloJetMatchEntry(evi);

      if((100*evi / NEvents) % 5 == 0 && 100*evi / NEvents > evi_frac) cout << "evt frac: " << evi_frac << "%" << endl;
      evi_frac = 100 * evi/NEvents;

      // global event cuts
      //cout << "Applying global event cuts..." << endl;
      if(em->pthat <= pthatcut) continue;
      if(fabs(em->vz) > 15.0) continue;
      if(em->hiBin > 190) continue;
      if(em->checkEventFilter()) continue;
      //cout << "Event #" << evi << " passed the global cuts!" << endl;

      if(applyJet60Trigger){
	if(em->HLT_HICsAK4PFJet60Eta1p5_v1 == 0) continue;
      }
      if(applyJet80Trigger){
	//std::cout << "jet80 (event "<< evi << ") = " << em->HLT_HICsAK4PFJet80Eta1p5_v1 << std::endl;
	if(em->HLT_HICsAK4PFJet80Eta1p5_v1 == 0) continue;
      }

      // apply HLT
      int triggerDecision = em->HLT_HIL3Mu12_v1;
      int triggerDecision_Prescl = em->HLT_HIL3Mu12_v1_Prescl;
      //if(!triggerIsOn(triggerDecision,triggerDecision_Prescl)) continue;

      // RECO VARIABLES
	
      int matchFlag[10] = {0,0,0,0,0,0,0,0,0,0};

      int CentralityIndex = getCentBin(em->hiBin-hiBinShift);
    
      if(CentralityIndex < 0) continue;

      const int nominalClass = nominalCentClass(em->hiBin-hiBinShift);   // for the 4-class fits


    
      double w_reweight_hiBin = fitFxn_hiBin->Eval(em->hiBin-hiBinShift);

	
      //double w_reweight_vz = fitFxn_vz->Eval(em->vz);
      double w_reweight_vz = 1.0;
	
      double w_pthat = 1.0;
      if(doPThatWeight){
	w_pthat = em->weight;
      }
    
      double w = w_pthat * w_reweight_vz * w_reweight_hiBin;

double leadingMatchedRecoJetPt = -999.0;

      // RESPONSE FROM THE RECO SIDE, with the nominal scan's reco selection
      // (PYTHIAHYDJET_scan.C reco-jet loop), so the measured axis of the matrix
      // is the same population the unfolded spectrum is.  Per reco jet:
      //   pThat-correlation cut  reco pT / pThat against the nominal-class fit,
      //                          PER JET (was: per event, on the leading reco jet)
      //   jetTrkMax filter       doJetTrkMaxFilter (was: in the file name only)
      //   eta-phi mask           doEtaPhiMask
      //   |eta| < etaMax, pT > jetPtCut
      // A reco jet failing any of these is dropped; if it had a gen match, that
      // gen jet becomes a miss.  A reco jet passing them with refpt > 0 enters the
      // matrix at (reco pT, refpt) WHATEVER the gen jet's eta: matched = "has a
      // forest reference", the nominal scan's definition (its flavor 18 is
      // refpt < 0).  refpt <= 0 is a fake -> h_unmatchedRecoJetPt.  Gen jets with
      // |eta| < etaMax and no accepted reco match are the misses.
      //
      // The old event-level pThat filter ran after the reco-jet loop, so the fakes
      // and the reco diagnostics were filled from events it then rejected (the
      // KNOWN BUG of 2026-09-25).  With the cut per jet there is no event-level
      // filter left and every histogram sees the same jets.
      std::vector<int> genMatchedReco(em->ngj, -1);   // accepted reco jet matched to gen jet j
      // analysis muon tag: each muon tags one jet (PYTHIAHYDJET_scan.C's matchFlagR)
      std::vector<int> matchFlagR(useAnalysisMuonTag && em->nMu > 0 ? em->nMu : 1, 0);

      // RECO JET LOOP -- matrix, fakes, reco diagnostics
      for(int i = 0; i < em->njet; i++){
	JEC.SetJetPT(em->rawpt[i]);
	JEC.SetJetEta(em->jeteta[i]);
	JEC.SetJetPhi(em->jetphi[i]);
	// manual JEC on rawpt (AK4Calo or AK4PF, matching the collection), or the
	// forest jtpt (config_PYTHIAHYDJET.h: useManualJEC)
	double recoJetPt_i = useManualJEC ? JEC.GetCorrectedPT() : em->jetpt[i];
	double recoJetEta_i = em->jeteta[i];
	double recoJetPhi_i = em->jetphi[i];
	double refJetPt_i = em->refpt[i];
	int recoJetFlavor_i = recoJetFlavorFor(useCaloJetsOverride, i, em);

	// calo-match jet ID, ahead of every reco-jet histogram, the muon tag and the
	// response: a rejected jet leaves its gen jet unmatched, i.e. a miss
	double caloRawPtMatched = -1.0;
	double caloDrMin = closestCaloJetDr(recoJetEta_i, recoJetPhi_i, caloRawPtMatched);
	bool caloMatched = (caloDrMin < caloMatchDr);
	if(requireCaloJetMatch && !caloMatched) continue;

	// ---- nominal reco selection, same order as PYTHIAHYDJET_scan.C
	if(doPThatCorrelationFilter){
	  TF1 *fPThatCorr = nullptr;
	  switch(nominalClass){
	  case 1: fPThatCorr = useCaloJetsOverride ? fitFxn_PYTHIAHYDJET_pThatCorrelation_caloJets_C1 : fitFxn_PYTHIAHYDJET_pThatCorrelation_C1; break;
	  case 2: fPThatCorr = useCaloJetsOverride ? fitFxn_PYTHIAHYDJET_pThatCorrelation_caloJets_C2 : fitFxn_PYTHIAHYDJET_pThatCorrelation_C2; break;
	  case 3: fPThatCorr = useCaloJetsOverride ? fitFxn_PYTHIAHYDJET_pThatCorrelation_caloJets_C3 : fitFxn_PYTHIAHYDJET_pThatCorrelation_C3; break;
	  default: fPThatCorr = useCaloJetsOverride ? fitFxn_PYTHIAHYDJET_pThatCorrelation_caloJets_C4 : fitFxn_PYTHIAHYDJET_pThatCorrelation_C4; break;
	  }
	  if((recoJetPt_i / em->pthat) > fPThatCorr->Eval(em->pthat)) continue;
	}
	if(doJetTrkMaxFilter){
	  if(!passesJetTrkMaxFilter(em->jetTrkMax[i],recoJetPt_i)) continue;
	}
	if(doEtaPhiMask){
	  if(etaPhiMask(recoJetEta_i,recoJetPhi_i)) continue;
	}
	if(TMath::Abs(recoJetEta_i) > etaMax || recoJetPt_i < jetPtCut) continue;
	if(doRemoveHYDJETjet && refJetPt_i > 0){
	  if(remove_HYDJET_jet(em->pthat, refJetPt_i)) continue;
	}

	// analysis reco-muon tag, on every accepted reco jet in forest order as in
	// PYTHIAHYDJET_scan.C, so muons are consumed by the same jets there and here
	bool anaMuTag_i = false;
	if(useAnalysisMuonTag){
	  for(int m = 0; m < em->nMu; m++){
	    if(!passesRecoMuonCuts(em, m, matchFlagR.data())) continue;
	    if(doWDecayFilter && isWDecayMuon(em->muPt->at(m), recoJetPt_i)) continue;
	    if(getDr(em->muEta->at(m), em->muPhi->at(m), recoJetEta_i, recoJetPhi_i) < epsilon_mm){
	      matchFlagR[m] = 1;
	      anaMuTag_i = true;
	    }
	  }
	}

	// nearest gen jet to this reco jet, by dR, regardless of any cut; and the
	// gen jet the forest reference points at (refpt == genjetpt)
	double dRnearest_i   = 999.0;
	double nearestGenPt_i = -1.0;
	int    refGen_i = -1;
	for(int j = 0; j < em->ngj ; j++){
	  if(refJetPt_i > 0 && refJetPt_i == em->genjetpt[j]) refGen_i = j;
	  double dR_ij = getDr(recoJetEta_i, recoJetPhi_i, em->genjeteta[j], em->genjetphi[j]);
	  if(dR_ij < dRnearest_i){
	    dRnearest_i   = dR_ij;
	    nearestGenPt_i = em->genjetpt[j];
	  }
	}

	// Reco-jet-indexed diagnostics: one loop, one weight, one condition, so
	// all = matchedDr + unmatchedDr holds bin by bin.
	bool hasGenJetMatchDr_i = (dRnearest_i < recoGenMatchDr);

	h_recoJetPt_all[0]->Fill(recoJetPt_i,w);
	h_recoJetPt_all[CentralityIndex]->Fill(recoJetPt_i,w);

	if(hasGenJetMatchDr_i){
	  h_recoJetPt_matchedDr[0]->Fill(recoJetPt_i,w);
	  h_recoJetPt_matchedDr[CentralityIndex]->Fill(recoJetPt_i,w);
	}
	else{
	  h_recoJetPt_unmatchedDr[0]->Fill(recoJetPt_i,w);
	  h_recoJetPt_unmatchedDr[CentralityIndex]->Fill(recoJetPt_i,w);
	}

	// calo-match diagnostics
	if(g_caloTree){
	  h_caloJetDr_recoJetPt[0]->Fill(recoJetPt_i,caloDrMin,w);
	  h_caloJetDr_recoJetPt[CentralityIndex]->Fill(recoJetPt_i,caloDrMin,w);
	  if(caloMatched){
	    h_caloJetRawPt_recoJetPt[0]->Fill(recoJetPt_i,caloRawPtMatched,w);
	    h_caloJetRawPt_recoJetPt[CentralityIndex]->Fill(recoJetPt_i,caloRawPtMatched,w);
	  }
	}

	// dR is capped at the histogram range so overflow does not hide the
	// "no gen jet anywhere near" population; genPt = -1 lands in underflow
	// when the event has no gen jets at all.
	h_recoPt_dRnearestGen[0]->Fill(recoJetPt_i, TMath::Min(dRnearest_i,0.9999), w);
	h_recoPt_dRnearestGen[CentralityIndex]->Fill(recoJetPt_i, TMath::Min(dRnearest_i,0.9999), w);
	h_recoPt_nearestGenPt[0]->Fill(recoJetPt_i, nearestGenPt_i, w);
	h_recoPt_nearestGenPt[CentralityIndex]->Fill(recoJetPt_i, nearestGenPt_i, w);

	// ---- fakes
	if(refJetPt_i <= 0){
	  h_unmatchedRecoJetPt[0][0]->Fill(recoJetPt_i,w);
	  h_unmatchedRecoJetPt[CentralityIndex][0]->Fill(recoJetPt_i,w);
	  if(fabs(recoJetFlavor_i) == 5){
	    h_unmatchedRecoJetPt[0][1]->Fill(recoJetPt_i,w);
	    h_unmatchedRecoJetPt[CentralityIndex][1]->Fill(recoJetPt_i,w);
	  }
	  else if(fabs(recoJetFlavor_i) == 4){
	    h_unmatchedRecoJetPt[0][2]->Fill(recoJetPt_i,w);
	    h_unmatchedRecoJetPt[CentralityIndex][2]->Fill(recoJetPt_i,w);
	  }
	  else if(fabs(recoJetFlavor_i) == 1 || fabs(recoJetFlavor_i) == 2){
	    h_unmatchedRecoJetPt[0][3]->Fill(recoJetPt_i,w);
	    h_unmatchedRecoJetPt[CentralityIndex][3]->Fill(recoJetPt_i,w);
	  }
	  else if(fabs(recoJetFlavor_i) == 3){
	    h_unmatchedRecoJetPt[0][4]->Fill(recoJetPt_i,w);
	    h_unmatchedRecoJetPt[CentralityIndex][4]->Fill(recoJetPt_i,w);
	  }
	  else if(recoJetFlavor_i == 21){
	    h_unmatchedRecoJetPt[0][5]->Fill(recoJetPt_i,w);
	    h_unmatchedRecoJetPt[CentralityIndex][5]->Fill(recoJetPt_i,w);
	  }
	  else if(recoJetFlavor_i == 0){
	    h_unmatchedRecoJetPt[0][6]->Fill(recoJetPt_i,w);
	    h_unmatchedRecoJetPt[CentralityIndex][6]->Fill(recoJetPt_i,w);
	  }
	  else{};
	  continue;
	}

	// ---- matched: response entry at (reco pT, refpt)
	double w_jet = w;
	double x = refJetPt_i;
	const bool hasGenEta = (refGen_i >= 0);   // gen eta known only via refpt == genjetpt
	double y = hasGenEta ? em->genjeteta[refGen_i] : -999.;
	if(refGen_i >= 0) genMatchedReco[refGen_i] = i;

	const bool hasRecoJetMuon = useAnalysisMuonTag ? anaMuTag_i
	                                               : (em->mupt[i] > muPtCut && fabs(em->mueta[i]) < 2.);
	const bool hasRecoMuon = hasRecoJetMuon;
	int jetFlavorInt = recoJetFlavor_i;

	double matchedRecoJetPt = recoJetPt_i;
	if(matchedRecoJetPt > leadingMatchedRecoJetPt) leadingMatchedRecoJetPt = matchedRecoJetPt;

	JEU.SetJetPT(matchedRecoJetPt);
	JEU.SetJetEta(em->jeteta[i]);
	JEU.SetJetPhi(em->jetphi[i]);

	if(apply_JEU_shift_up){
	  matchedRecoJetPt = matchedRecoJetPt * (1 + JEU.GetUncertainty().second);
	}
	else if(apply_JEU_shift_down){
	  matchedRecoJetPt = matchedRecoJetPt * (1 - JEU.GetUncertainty().first);
	}

	if(apply_JER_smear){
	  double sigma = 0.663*JER_fxn[nominalClass]->Eval(matchedRecoJetPt); // apply a 20% smear
	  matchedRecoJetPt = matchedRecoJetPt * randomGenerator->Gaus(1.0,sigma);
	}

	if(doJERCorrection){
	  double k_JERCorrection = TMath::Sqrt(fitFxn_PYTHIA_JERCorrection->Eval(x)*fitFxn_PYTHIA_JERCorrection->Eval(x) - 1.);
	  double sigma_JERCorrection = k_JERCorrection*JER_fxn[nominalClass]->Eval(matchedRecoJetPt);
	  matchedRecoJetPt = matchedRecoJetPt * randomGenerator->Gaus(1.0,sigma_JERCorrection);
	}

	// neutrino energy shift for a fraction of jets; a failed dice roll leaves the
	// jet unshifted but still in the response (as before the restructure)
	if(doBJetNeutrinoEnergyShift){
	  double diceRoll = randomGenerator->Rndm();
	  if(diceRoll <= neutrino_tag_fraction->GetBinContent(neutrino_tag_fraction->FindBin(matchedRecoJetPt))){
	    neutrino_energy_map_proj = (TH1D*) neutrino_energy_map->ProjectionX("neutrino_energy_map_proj", neutrino_energy_map->GetYaxis()->FindBin(matchedRecoJetPt),neutrino_energy_map->GetYaxis()->FindBin(matchedRecoJetPt)+1);
	    matchedRecoJetPt += neutrino_energy_map_proj->GetRandom();
	  }
	}

	if(!onlyMuTaggedJets || hasRecoJetMuon) {
	  if(doBJetSpectraReweightToData){
	    if(nominalClass == 4) w_jet = w_jet * fitFxn_PYTHIAHYDJET_BJetSpectraReweightToData_C4->Eval(matchedRecoJetPt);
	    else if(nominalClass == 3) w_jet = w_jet * fitFxn_PYTHIAHYDJET_BJetSpectraReweightToData_C3->Eval(matchedRecoJetPt);
	    else if(nominalClass == 2) w_jet = w_jet * fitFxn_PYTHIAHYDJET_BJetSpectraReweightToData_C2->Eval(matchedRecoJetPt);
	    else if(nominalClass == 1) w_jet = w_jet * fitFxn_PYTHIAHYDJET_BJetSpectraReweightToData_C1->Eval(matchedRecoJetPt);
	    else{};
	  }

	  // response_C0.Fill(matchedRecoJetPt,x,w);
	  // if(CentralityIndex == 4) response_C4.Fill(matchedRecoJetPt,x,w);
	  // else if(CentralityIndex == 3) response_C3.Fill(matchedRecoJetPt,x,w);
	  // else if(CentralityIndex == 2) response_C2.Fill(matchedRecoJetPt,x,w);
	  // else if(CentralityIndex == 1) response_C1.Fill(matchedRecoJetPt,x,w);
	  // else{};
	
	  h_matchedRecoJetPt_genJetPt[0][0]->Fill(matchedRecoJetPt,x,w_jet);
	  h_matchedRecoJetPt_genJetPt[CentralityIndex][0]->Fill(matchedRecoJetPt,x,w_jet);
	  for(int t = 0; t < kNDistTilt; t++){
	    double wd = w_jet*distTiltWeight(x,kDistTiltExpo[t]);
	    h_matchedRecoJetPt_genJetPt_distTilt[0][0][0][t]->Fill(matchedRecoJetPt,x,wd);
	    h_matchedRecoJetPt_genJetPt_distTilt[CentralityIndex][0][0][t]->Fill(matchedRecoJetPt,x,wd);
	    h_matchedRecoJetPt_genJetPt_var_distTilt[0][0][0][t]->Fill(matchedRecoJetPt,x,wd);
	    h_matchedRecoJetPt_genJetPt_var_distTilt[CentralityIndex][0][0][t]->Fill(matchedRecoJetPt,x,wd);
	    if(hasRecoJetMuon){
	      h_matchedRecoJetPt_genJetPt_distTilt[0][0][1][t]->Fill(matchedRecoJetPt,x,wd);
	      h_matchedRecoJetPt_genJetPt_distTilt[CentralityIndex][0][1][t]->Fill(matchedRecoJetPt,x,wd);
	      h_matchedRecoJetPt_genJetPt_var_distTilt[0][0][1][t]->Fill(matchedRecoJetPt,x,wd);
	      h_matchedRecoJetPt_genJetPt_var_distTilt[CentralityIndex][0][1][t]->Fill(matchedRecoJetPt,x,wd);
	    }
	  }
	  if(hasRecoJetMuon){
	    h_matchedRecoJetPt_genJetPt_muTagged[0][0]->Fill(matchedRecoJetPt,x,w_jet);
	    h_matchedRecoJetPt_genJetPt_muTagged[CentralityIndex][0]->Fill(matchedRecoJetPt,x,w_jet);
	    h_matchedRecoJetPt_genJetPt_var_muTagged[0][0]->Fill(matchedRecoJetPt,x,w_jet);
	    h_matchedRecoJetPt_genJetPt_var_muTagged[CentralityIndex][0]->Fill(matchedRecoJetPt,x,w_jet);
	  }

	  h_matchedRecoJetPt_genJetPt_var[0][0]->Fill(matchedRecoJetPt,x,w_jet);
	  h_matchedRecoJetPt_genJetPt_var[CentralityIndex][0]->Fill(matchedRecoJetPt,x,w_jet);

	  h_matchedRecoJetPtOverGenJetPt_genJetPt[0][0]->Fill(matchedRecoJetPt/x,x,w_jet);
	  h_matchedRecoJetPtOverGenJetPt_genJetPt[CentralityIndex][0]->Fill(matchedRecoJetPt/x,x,w_jet);
	
	  if(x>100 && hasGenEta){
	    h_matchedRecoJetPtOverGenJetPt_genJetEta[0][0]->Fill(matchedRecoJetPt/x,y,w_jet);
	    h_matchedRecoJetPtOverGenJetPt_genJetEta[CentralityIndex][0]->Fill(matchedRecoJetPt/x,y,w_jet);
	  }


	  if(fabs(jetFlavorInt) == 5){
	    h_matchedRecoJetPt_genJetPt_var[0][1]->Fill(matchedRecoJetPt,x,w_jet);
	    h_matchedRecoJetPt_genJetPt_var[CentralityIndex][1]->Fill(matchedRecoJetPt,x,w_jet);
	    h_matchedRecoJetPt_genJetPt[0][1]->Fill(matchedRecoJetPt,x,w_jet);
	    h_matchedRecoJetPt_genJetPt[CentralityIndex][1]->Fill(matchedRecoJetPt,x,w_jet);
	    for(int t = 0; t < kNDistTilt; t++){
	      double wd = w_jet*distTiltWeight(x,kDistTiltExpo[t]);
	      h_matchedRecoJetPt_genJetPt_distTilt[0][1][0][t]->Fill(matchedRecoJetPt,x,wd);
	      h_matchedRecoJetPt_genJetPt_distTilt[CentralityIndex][1][0][t]->Fill(matchedRecoJetPt,x,wd);
	      h_matchedRecoJetPt_genJetPt_var_distTilt[0][1][0][t]->Fill(matchedRecoJetPt,x,wd);
	      h_matchedRecoJetPt_genJetPt_var_distTilt[CentralityIndex][1][0][t]->Fill(matchedRecoJetPt,x,wd);
	      if(hasRecoJetMuon){
		h_matchedRecoJetPt_genJetPt_distTilt[0][1][1][t]->Fill(matchedRecoJetPt,x,wd);
		h_matchedRecoJetPt_genJetPt_distTilt[CentralityIndex][1][1][t]->Fill(matchedRecoJetPt,x,wd);
		h_matchedRecoJetPt_genJetPt_var_distTilt[0][1][1][t]->Fill(matchedRecoJetPt,x,wd);
		h_matchedRecoJetPt_genJetPt_var_distTilt[CentralityIndex][1][1][t]->Fill(matchedRecoJetPt,x,wd);
	      }
	    }
	    if(hasRecoJetMuon){
	      h_matchedRecoJetPt_genJetPt_muTagged[0][1]->Fill(matchedRecoJetPt,x,w_jet);
	      h_matchedRecoJetPt_genJetPt_muTagged[CentralityIndex][1]->Fill(matchedRecoJetPt,x,w_jet);
	      h_matchedRecoJetPt_genJetPt_var_muTagged[0][1]->Fill(matchedRecoJetPt,x,w_jet);
	      h_matchedRecoJetPt_genJetPt_var_muTagged[CentralityIndex][1]->Fill(matchedRecoJetPt,x,w_jet);
	    }
	    h_matchedRecoJetPtOverGenJetPt_genJetPt[0][1]->Fill(matchedRecoJetPt/x,x,w_jet);
	    h_matchedRecoJetPtOverGenJetPt_genJetPt[CentralityIndex][1]->Fill(matchedRecoJetPt/x,x,w_jet);
				
	    if(x>100 && hasGenEta){
	      h_matchedRecoJetPtOverGenJetPt_genJetEta[0][1]->Fill(matchedRecoJetPt/x,y,w_jet);
	      h_matchedRecoJetPtOverGenJetPt_genJetEta[CentralityIndex][1]->Fill(matchedRecoJetPt/x,y,w_jet);
	    }
	  } 
	  if(fabs(jetFlavorInt) == 4){
	    h_matchedRecoJetPt_genJetPt_var[0][2]->Fill(matchedRecoJetPt,x,w_jet);
	    h_matchedRecoJetPt_genJetPt_var[CentralityIndex][2]->Fill(matchedRecoJetPt,x,w_jet);
	    h_matchedRecoJetPt_genJetPt[0][2]->Fill(matchedRecoJetPt,x,w_jet);
	    h_matchedRecoJetPt_genJetPt[CentralityIndex][2]->Fill(matchedRecoJetPt,x,w_jet);
	    h_matchedRecoJetPtOverGenJetPt_genJetPt[0][2]->Fill(matchedRecoJetPt/x,x,w_jet);
	    h_matchedRecoJetPtOverGenJetPt_genJetPt[CentralityIndex][2]->Fill(matchedRecoJetPt/x,x,w_jet);

	    if(x>100 && hasGenEta){
	      h_matchedRecoJetPtOverGenJetPt_genJetEta[0][2]->Fill(matchedRecoJetPt/x,y,w_jet);
	      h_matchedRecoJetPtOverGenJetPt_genJetEta[CentralityIndex][2]->Fill(matchedRecoJetPt/x,y,w_jet);
	    }

	  } 
	  if(fabs(jetFlavorInt) == 1 || fabs(jetFlavorInt) == 2){
	    h_matchedRecoJetPt_genJetPt_var[0][3]->Fill(matchedRecoJetPt,x,w_jet);
	    h_matchedRecoJetPt_genJetPt_var[CentralityIndex][3]->Fill(matchedRecoJetPt,x,w_jet);
	    h_matchedRecoJetPt_genJetPt[0][3]->Fill(matchedRecoJetPt,x,w_jet);
	    h_matchedRecoJetPt_genJetPt[CentralityIndex][3]->Fill(matchedRecoJetPt,x,w_jet);
	    h_matchedRecoJetPtOverGenJetPt_genJetPt[0][3]->Fill(matchedRecoJetPt/x,x,w_jet);
	    h_matchedRecoJetPtOverGenJetPt_genJetPt[CentralityIndex][3]->Fill(matchedRecoJetPt/x,x,w_jet);

	    if(x>100 && hasGenEta){
	      h_matchedRecoJetPtOverGenJetPt_genJetEta[0][3]->Fill(matchedRecoJetPt/x,y,w_jet);
	      h_matchedRecoJetPtOverGenJetPt_genJetEta[CentralityIndex][3]->Fill(matchedRecoJetPt/x,y,w_jet);
	    }
	  } 
	  if(fabs(jetFlavorInt) == 3){
	    h_matchedRecoJetPt_genJetPt_var[0][4]->Fill(matchedRecoJetPt,x,w_jet);
	    h_matchedRecoJetPt_genJetPt_var[CentralityIndex][4]->Fill(matchedRecoJetPt,x,w_jet);
	    h_matchedRecoJetPt_genJetPt[0][4]->Fill(matchedRecoJetPt,x,w_jet);
	    h_matchedRecoJetPt_genJetPt[CentralityIndex][4]->Fill(matchedRecoJetPt,x,w_jet);
	    h_matchedRecoJetPtOverGenJetPt_genJetPt[0][4]->Fill(matchedRecoJetPt/x,x,w_jet);
	    h_matchedRecoJetPtOverGenJetPt_genJetPt[CentralityIndex][4]->Fill(matchedRecoJetPt/x,x,w_jet);

	    if(x>100 && hasGenEta){
	      h_matchedRecoJetPtOverGenJetPt_genJetEta[0][4]->Fill(matchedRecoJetPt/x,y,w_jet);
	      h_matchedRecoJetPtOverGenJetPt_genJetEta[CentralityIndex][4]->Fill(matchedRecoJetPt/x,y,w_jet);
	    }

	  }  
	  if(jetFlavorInt == 21){
	    h_matchedRecoJetPt_genJetPt_var[0][5]->Fill(matchedRecoJetPt,x,w_jet);
	    h_matchedRecoJetPt_genJetPt_var[CentralityIndex][5]->Fill(matchedRecoJetPt,x,w_jet);
	    h_matchedRecoJetPt_genJetPt[0][5]->Fill(matchedRecoJetPt,x,w_jet);
	    h_matchedRecoJetPt_genJetPt[CentralityIndex][5]->Fill(matchedRecoJetPt,x,w_jet);
	    h_matchedRecoJetPtOverGenJetPt_genJetPt[0][5]->Fill(matchedRecoJetPt/x,x,w_jet);
	    h_matchedRecoJetPtOverGenJetPt_genJetPt[CentralityIndex][5]->Fill(matchedRecoJetPt/x,x,w_jet);

	    if(x>100 && hasGenEta){
	      h_matchedRecoJetPtOverGenJetPt_genJetEta[0][5]->Fill(matchedRecoJetPt/x,y,w_jet);
	      h_matchedRecoJetPtOverGenJetPt_genJetEta[CentralityIndex][5]->Fill(matchedRecoJetPt/x,y,w_jet);
	    }

	  }  
	  if(jetFlavorInt == 0){
	    h_matchedRecoJetPt_genJetPt_var[0][6]->Fill(matchedRecoJetPt,x,w_jet);
	    h_matchedRecoJetPt_genJetPt_var[CentralityIndex][6]->Fill(matchedRecoJetPt,x,w_jet);
	    h_matchedRecoJetPt_genJetPt[0][6]->Fill(matchedRecoJetPt,x,w_jet);
	    h_matchedRecoJetPt_genJetPt[CentralityIndex][6]->Fill(matchedRecoJetPt,x,w_jet);
	    h_matchedRecoJetPtOverGenJetPt_genJetPt[0][6]->Fill(matchedRecoJetPt/x,x,w_jet);
	    h_matchedRecoJetPtOverGenJetPt_genJetPt[CentralityIndex][6]->Fill(matchedRecoJetPt/x,x,w_jet);

	    if(x>100 && hasGenEta){
	      h_matchedRecoJetPtOverGenJetPt_genJetEta[0][6]->Fill(matchedRecoJetPt/x,y,w_jet);
	      h_matchedRecoJetPtOverGenJetPt_genJetEta[CentralityIndex][6]->Fill(matchedRecoJetPt/x,y,w_jet);
	    }

	  }

	  if(hasRecoMuon){
	    h_inclGenJetPt_inclRecoMuonTag_flavor[0]->Fill(x,jetFlavorInt,w_jet);
	    h_inclGenJetPt_inclRecoMuonTag_flavor[CentralityIndex]->Fill(x,jetFlavorInt,w_jet);
	  }
	}
      }
      // END RECO JET LOOP

      // GEN JET LOOP -- misses and gen-level spectra, gen acceptance
      for(int i = 0; i < em->ngj ; i++){

	double w_jet = w;
	double x = em->genjetpt[i];
	double y = em->genjeteta[i];
	double z = em->genjetphi[i];

	if(TMath::Abs(y) > etaMax) continue;

	// flavor from the accepted reco jet matched to it; 19 = no accepted match
	const int k = genMatchedReco[i];
	int jetFlavorInt = (k >= 0) ? recoJetFlavorFor(useCaloJetsOverride, k, em) : 19;

	if(k < 0){
	  h_unmatchedGenJetPt[0]->Fill(x,w_jet);
	  h_unmatchedGenJetPt[CentralityIndex]->Fill(x,w_jet);
	  for(int t = 0; t < kNDistTilt; t++){
	    h_unmatchedGenJetPt_distTilt[0][t]->Fill(x,w_jet*distTiltWeight(x,kDistTiltExpo[t]));
	    h_unmatchedGenJetPt_distTilt[CentralityIndex][t]->Fill(x,w_jet*distTiltWeight(x,kDistTiltExpo[t]));
	  }
	}

	h_inclGenJetPt_flavor[0]->Fill(x,jetFlavorInt,w_jet);
	h_inclGenJetPt_flavor[CentralityIndex]->Fill(x,jetFlavorInt,w_jet);
	// begin gen-muon loop

	bool hasGenMuon = false;
      
	for(int j = 0; j < em->gpptp->size(); j++){

	  if(hasGenMuon) continue;

	  if(TMath::Abs(em->gppdgIDp->at(j)) != 13) continue;

	  if(isWDecayMuon(em->gpptp->at(j),x)) continue; // skip if "WDecay" muon (has majority of jet pt)

	  double genMuonPt_j = em->gpptp->at(j);
	  double genMuonEta_j = em->gpetap->at(j);
	  double genMuonPhi_j = em->gpphip->at(j);

	  if(genMuonPt_j < muPtCut || fabs(genMuonEta_j) > 2.0) continue;

	  if(getDr(genMuonEta_j,genMuonPhi_j,y,z) < deltaRCut){
	    hasGenMuon = true;
	    h_inclGenJetPt_inclGenMuonTag_flavor[0]->Fill(x,jetFlavorInt,w_jet);
	    h_inclGenJetPt_inclGenMuonTag_flavor[CentralityIndex]->Fill(x,jetFlavorInt,w_jet);
	  }

	} // end gen-muon loop

      }
      // END GEN JET LOOP

      if(leadingMatchedRecoJetPt > 0){
	h_leadingRecoJetPtOverPThat_pThat[0]->Fill(leadingMatchedRecoJetPt / em->pthat, em->pthat,w);
	h_leadingRecoJetPtOverPThat_pThat[CentralityIndex]->Fill(leadingMatchedRecoJetPt / em->pthat, em->pthat,w);
      }
	

    } // END EVENT LOOP
    delete f;
    // WRITE
    auto wf = TFile::Open(output,"recreate");

    for(int j = 0; j < NCentralityIndices; j++){

      h_inclGenJetPt_flavor[j]->Write();
      h_inclGenJetPt_inclGenMuonTag_flavor[j]->Write();
      h_inclGenJetPt_inclRecoMuonTag_flavor[j]->Write();

      h_unmatchedGenJetPt[j]->Write();
      for(int f = 0; f < 2; f++){
	for(int m = 0; m < 2; m++){
	  for(int t = 0; t < kNDistTilt; t++){
	    h_matchedRecoJetPt_genJetPt_distTilt[j][f][m][t]->Write();
	    h_matchedRecoJetPt_genJetPt_var_distTilt[j][f][m][t]->Write();
	  }
	}
      }
      for(int t = 0; t < kNDistTilt; t++) h_unmatchedGenJetPt_distTilt[j][t]->Write();
      for(int f = 0; f < 2; f++){
	h_matchedRecoJetPt_genJetPt_muTagged[j][f]->Write();
	h_matchedRecoJetPt_genJetPt_var_muTagged[j][f]->Write();
      }
      h_unmatchedRecoJetPt[j][0]->Write();
      h_recoJetPt_all[j]->Write();
      h_recoJetPt_matchedDr[j]->Write();
      h_recoJetPt_unmatchedDr[j]->Write();
      h_caloJetDr_recoJetPt[j]->Write();
      h_caloJetRawPt_recoJetPt[j]->Write();
      h_recoPt_dRnearestGen[j]->Write();
      h_recoPt_nearestGenPt[j]->Write();
      // h_unmatchedRecoJetPt[j][1]->Write();
      // h_unmatchedRecoJetPt[j][2]->Write();
      // h_unmatchedRecoJetPt[j][3]->Write();
      // h_unmatchedRecoJetPt[j][4]->Write();
      // h_unmatchedRecoJetPt[j][5]->Write();
      // h_unmatchedRecoJetPt[j][6]->Write();
    
      h_matchedRecoJetPt_genJetPt[j][0]->Write();
      h_matchedRecoJetPt_genJetPt[j][1]->Write();
      h_matchedRecoJetPt_genJetPt[j][2]->Write();
      h_matchedRecoJetPt_genJetPt[j][3]->Write();
      h_matchedRecoJetPt_genJetPt[j][4]->Write();
      h_matchedRecoJetPt_genJetPt[j][5]->Write();
      h_matchedRecoJetPt_genJetPt[j][6]->Write();



      h_matchedRecoJetPt_genJetPt_var[j][0]->Write();
      h_matchedRecoJetPt_genJetPt_var[j][1]->Write();
      h_matchedRecoJetPt_genJetPt_var[j][2]->Write();
      h_matchedRecoJetPt_genJetPt_var[j][3]->Write();
      h_matchedRecoJetPt_genJetPt_var[j][4]->Write();
      h_matchedRecoJetPt_genJetPt_var[j][5]->Write();
      h_matchedRecoJetPt_genJetPt_var[j][6]->Write();
    
      h_matchedRecoJetPtOverGenJetPt_genJetPt[j][0]->Write();
      h_matchedRecoJetPtOverGenJetPt_genJetPt[j][1]->Write();
      h_matchedRecoJetPtOverGenJetPt_genJetPt[j][2]->Write();
      h_matchedRecoJetPtOverGenJetPt_genJetPt[j][3]->Write();
      h_matchedRecoJetPtOverGenJetPt_genJetPt[j][4]->Write();
      h_matchedRecoJetPtOverGenJetPt_genJetPt[j][5]->Write();
      h_matchedRecoJetPtOverGenJetPt_genJetPt[j][6]->Write();

      h_matchedRecoJetPtOverGenJetPt_genJetEta[j][0]->Write();
      // h_matchedRecoJetPtOverGenJetPt_genJetEta[j][1]->Write();
      // h_matchedRecoJetPtOverGenJetPt_genJetEta[j][2]->Write();
      // h_matchedRecoJetPtOverGenJetPt_genJetEta[j][3]->Write();
      // h_matchedRecoJetPtOverGenJetPt_genJetEta[j][4]->Write();
      // h_matchedRecoJetPtOverGenJetPt_genJetEta[j][5]->Write();
      // h_matchedRecoJetPtOverGenJetPt_genJetEta[j][6]->Write();

      h_leadingRecoJetPtOverPThat_pThat[j]->Write();


    }

    // response_C4.Write();
    // response_C3.Write();
    // response_C2.Write();
    // response_C1.Write();
    // response_C0.Write();




    printCaloPFMatchSummary();
    wf->Close();
    return;
    // END WRITE

  }

}







