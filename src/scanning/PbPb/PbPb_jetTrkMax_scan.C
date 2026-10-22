// Leading-track ptRel in PbPb MinBias data, per ultraFine centrality slice, for a
// data-driven light template (ATLAS-style, arXiv:2204.13530: track-jet pairs in
// inclusive-jet data, track pT reweighted to the light-jet muon pT). The PbPb
// counterpart of src/scanning/pp/pp_jetTrkMax_scan.C. 2026-10-19 (user).
//
// Input: the HIMinimumBias0 Part1 forest with calo and flow jets, PF candidates and
// towers (the PbPb_caloTowerAnalyzer.C production, CRAB 260926), one file per index of
//   fileNames/fileNames_HIMinimumBias0_Part1_withCaloAndFlowJets_withPFAndTowers_fresh.txt
// These do not store the jet trackMax direction, so the leading track is found here:
// the highest-pT SELECTED track within dR < 0.4 of the akCs4PF jet axis. (The forests
// that did store trackMaxEta/Phi, the 2025-01-31 HIHardProbes trackMax production and
// its skims, are gone from EOS, 2026-10-19.)
//
// Track source, decided per file and recorded in the provenance:
//   PbPbTracks/trackTree if the forest has it, with the cuts below;
//   otherwise the charged hadrons (pfId == 1) of pfcandAnalyzer/pfTree. PF has
//   already cleaned these (track quality, and high-pT tracks without the calorimeter
//   energy to match are not called charged hadrons), but none of the cuts below can
//   be applied to them, and PF candidates are the akCs4PF jets' own constituents.
//
// Track selection (the 2018 PbPb charged-hadron cuts): highPurity, sigma(pT)/pT < 0.1,
// |dz/sigma| < 3, |dxy/sigma| < 3, nHit >= 11, chi2/ndof/nLayer < 0.18, and above
// 20 GeV calorimeter matching, (pfEcal + pfHcal)/cosh(eta) > 0.5 pT, against fake
// high-pT tracks. A cut whose branch is missing from the forest is dropped with a
// warning and listed in the provenance.
//
// Histograms, per class C<i> (C0 = every slice), x = ptRel, y = track pT, z = jet pT
// on the muon-template binning (headers/AnalysisSetup/ptRelMuPt3D.h, track pT from
// 15 GeV), the same as the data's h_muptrel_mupt_recoJetPt_inclRecoMuonTag_triggerOn_C<i>:
//   h_jetTrkMaxPtRel_trkPt_recoJetPt_C<i>    leading selected track in the jet
//   h_jetTrkPtRel_trkPt_recoJetPt_C<i>       every selected track in the jet (ATLAS)
//   h_perpTrkMaxPtRel_trkPt_recoJetPt_C<i>   leading selected track in a dR < 0.4 cone
//                                            at the jet's eta, phi + pi/2, ptRel against
//                                            that cone axis; z = the jet's pT. The rate
//                                            at which the underlying event alone supplies
//                                            a track this hard: compare its integral with
//                                            the jet histogram's, per slice
//   h_inclRecoJetPt_C<i>                     every jet, |eta| < 1.6 (per-jet normalization)
//   h_jetTrkMax{Pt,PtOverJetPt,Eta,Phi,DR,PtRel}_C<i>J<j>, h_jetPt_C<i>
//                                            the 2025 1D set (jet pT > 80, track pT > 14)
// Inclusive jets are ~97% light, so the jet histograms are a light template up to
// that b + c share. The jet axis is the akCs4PF axis the muon ptRel is measured
// against, so the template carries the underlying event's smearing of that axis --
// the effect the PYTHIA+HYDJET templates have to model.
//
// JEC: Autumn18_HI_V8 AK4PF L2Relative + L2L3Residual, as PbPb_scan.C. Events: the
// five filters PbPb_scan.C and PbPb_pfCandAnalyzer.C apply, |vz| < 15, hiHF < 6000,
// no trigger bit (MinimumBias dataset, as the pfCand analyzer with
// applyMinBiasTrigger = false). hiBin: half-open ultraFine classes (getCentBinImpl.h).
//
// Run per forest file, from src/scanning/PbPb/:
//   root -l -b -q 'PbPb_jetTrkMax_scan.C(<index>)'   1-based line of the list
//   root -l -b -q 'PbPb_jetTrkMax_scan.C(0)'         query mode: list and output dir
// Condor: condor_PbPb_jetTrkMax_scan.py. Merge with hadd.
#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1D.h"
#include "TH3D.h"
#include "TMath.h"
#include "TNamed.h"
#include "TString.h"
#include "TSystem.h"
#include "TDatime.h"
#include "TVector2.h"
#include <vector>
#include <string>
#include <cstdio>

// event map
#include "../../../eventMap/eventMap.h"
// jet corrector
#include "../../../JetEnergyCorrections/JetCorrector.h"
// centrality scheme and analysis constants (common.h)
#include "../../../headers/config/config_centrality.h"
// ptRel x (track) pT x jet pT binning, shared with the muon templates
#include "../../../headers/AnalysisSetup/ptRelMuPt3D.h"
// getDr function
#include "../../../headers/functions/getDr.h"
// getJetPtBin function
#include "../../../headers/functions/getJetPtBin.h"
// getPtRel function
#include "../../../headers/functions/getPtRel.h"
// track source, selection and matching, shared with pp_jetTrkPtRel_scan.C
#include "../../../headers/functions/jetTrackPtRel.h"


namespace jtm {
const char *inputFileList = "../../../fileNames/fileNames_HIMinimumBias0_Part1_withCaloAndFlowJets_withPFAndTowers_fresh.txt";
const char *outputBaseDir = "/eos/cms/store/group/phys_heavyions/cbennett/scanningOutput/";
const char *outputName = "output_PbPb_MinBias_Part1_caloFlowPFTowers_leadTrkPtRel_trkPtAxis";
const char *jetTreeName = "akCs4PFJetAnalyzer/t";
const char *trkTreeName = "PbPbTracks/trackTree";
// as PbPb_scan.C and PbPb_pfCandAnalyzer.C
const int NeventFilters = 5;
std::string eventFilters[NeventFilters] = {"pprimaryVertexFilter", "HBHENoiseFilterResultRun2Loose", "collisionEventSelectionAODv2", "phfCoincFilter2Th4", "pclusterCompatibilityFilter"};
const double coneR = 0.4;          // track-jet match
const double trkPtFloor = 14.;     // lowest track pT any histogram uses (the 1D set)
const double trkPtCut3D = 15.;     // the 3D axis starts here
const double trkMaxCut1D = 14.;    // the 2025 1D set
const double jetPtCut1D = 80.;

TString outputDir(){
  // SCAN_DATE (set by myProcesses/condor/scan_condor.py) pins the date to the submission
  TString date;
  const char *scanDate = gSystem->Getenv("SCAN_DATE");
  if(scanDate && scanDate[0]) date = scanDate;
  else { TDatime dt; date = Form("%i-%i-%i", dt.GetYear(), dt.GetMonth(), dt.GetDay()); }
  return Form("%s%s%s_%s", outputBaseDir, outputName, CENT_SCHEME_SUFFIX, date.Data());
}

std::string shellLine(const char *cmd){
  std::string out; FILE *p = popen(cmd, "r"); if(!p) return "unavailable";
  char buf[512]; if(fgets(buf, sizeof(buf), p)) out = buf; pclose(p);
  while(!out.empty() && (out.back() == '\n' || out.back() == '\r')) out.pop_back();
  return out.empty() ? "unavailable" : out;
}
}


TH3D *h_jetTrkMaxPtRel_trkPt_recoJetPt[NCentralityIndices];    // leading track in the jet
TH3D *h_jetTrkPtRel_trkPt_recoJetPt[NCentralityIndices];       // every track in the jet
TH3D *h_perpTrkMaxPtRel_trkPt_recoJetPt[NCentralityIndices];   // leading track in the perpendicular cone
TH1D *h_inclRecoJetPt[NCentralityIndices];
TH1D *h_jetPt[NCentralityIndices];
TH1D *h_jetTrkMaxPt[NCentralityIndices][NJetPtIndices];
TH1D *h_jetTrkMaxPtOverJetPt[NCentralityIndices][NJetPtIndices];
TH1D *h_jetTrkMaxEta[NCentralityIndices][NJetPtIndices];
TH1D *h_jetTrkMaxPhi[NCentralityIndices][NJetPtIndices];
TH1D *h_jetTrkMaxDR[NCentralityIndices][NJetPtIndices];
TH1D *h_jetTrkMaxPtRel[NCentralityIndices][NJetPtIndices];


void PbPb_jetTrkMax_scan(int group = 1){

  const TString outputDir = jtm::outputDir();
  if(group == 0){ // query mode for condor submit scripts
    cout << "INPUTFILELIST=" << jtm::inputFileList << endl;
    cout << "OUTPUTDIR=" << outputDir << endl;
    return;
  }

  std::ifstream instr(jtm::inputFileList, std::ifstream::in);
  if(!instr.is_open()){ cout << "ERROR: file list not found: " << jtm::inputFileList << endl; return; }
  std::string filename, input;
  int ifile = 0;
  while(instr >> filename){ ifile++; if(ifile == group){ input = filename; break; } }
  if(input.empty()){ cout << "ERROR: no line " << group << " in " << jtm::inputFileList << endl; return; }
  if(gSystem->AccessPathName(outputDir)){ cout << "ERROR: output directory not found: " << outputDir << endl; return; }
  const TString output = Form("%s/PbPb_scan_output_%i.root", outputDir.Data(), group);
  cout << "input  = " << input << endl << "output = " << output << endl;

  // JET ENERGY CORRECTIONS (as PbPb_scan.C, PF)
  vector<string> Files;
  Files.push_back("../../../JetEnergyCorrections/Autumn18_HI_V8_DATA_L2Relative_AK4PF.txt");
  Files.push_back("../../../JetEnergyCorrections/Autumn18_HI_V8_DATA_L2L3Residual_AK4PF.txt");
  JetCorrector JEC(Files);

  for(int i = 0; i < NCentralityIndices; i++){

    const int cLo = (i == 0) ? centEdges[0] : centEdges[i-1];
    const int cHi = (i == 0) ? centEdges[NCentralityIndices-1] : centEdges[i];

    h_jetTrkMaxPtRel_trkPt_recoJetPt[i] = bookPtRelMuPt3D(Form("h_jetTrkMaxPtRel_trkPt_recoJetPt_C%i",i), Form("leading track in jet: p_{T}^{rel} x track p_{T} x jet p_{T}, hiBin %i - %i",cLo,cHi));
    h_jetTrkPtRel_trkPt_recoJetPt[i] = bookPtRelMuPt3D(Form("h_jetTrkPtRel_trkPt_recoJetPt_C%i",i), Form("every track in jet: p_{T}^{rel} x track p_{T} x jet p_{T}, hiBin %i - %i",cLo,cHi));
    h_perpTrkMaxPtRel_trkPt_recoJetPt[i] = bookPtRelMuPt3D(Form("h_perpTrkMaxPtRel_trkPt_recoJetPt_C%i",i), Form("leading track in perpendicular cone: p_{T}^{rel} x track p_{T} x jet p_{T}, hiBin %i - %i",cLo,cHi));

    h_inclRecoJetPt[i] = new TH1D(Form("h_inclRecoJetPt_C%i",i),Form("every jet, |#eta| < %.1f, hiBin %i - %i",etaMax,cLo,cHi),NPtBins,ptMin,ptMax);
    h_inclRecoJetPt[i]->Sumw2();
    h_jetPt[i] = new TH1D(Form("h_jetPt_C%i",i),Form("jetPt, hiBin %i - %i",cLo,cHi),NPtBins,ptMin,ptMax);
    h_jetPt[i]->Sumw2();

    for(int j = 0; j < NJetPtIndices; j++){

      const double jLo = (j == 0) ? jetPtEdges[0] : jetPtEdges[j-1];
      const double jHi = (j == 0) ? jetPtEdges[NJetPtIndices-1] : jetPtEdges[j];
      const TString tagS = TString::Format("hiBin %i - %i, p_{T}^{jet} %3.0f - %3.0f",cLo,cHi,jLo,jHi);
      const char *tag = tagS.Data();

      h_jetTrkMaxPt[i][j] = new TH1D(Form("h_jetTrkMaxPt_C%iJ%i",i,j),Form("jetTrkMaxPt, %s",tag),500,0,500);
      h_jetTrkMaxPtOverJetPt[i][j] = new TH1D(Form("h_jetTrkMaxPtOverJetPt_C%iJ%i",i,j),Form("jetTrkMaxPt / jetPt, %s",tag),100,0,1);
      h_jetTrkMaxEta[i][j] = new TH1D(Form("h_jetTrkMaxEta_C%iJ%i",i,j),Form("jetTrkMaxEta, %s",tag),NTrkEtaBins,trkEtaMin,trkEtaMax);
      h_jetTrkMaxPhi[i][j] = new TH1D(Form("h_jetTrkMaxPhi_C%iJ%i",i,j),Form("jetTrkMaxPhi, %s",tag),NPhiBins,phiMin,phiMax);
      h_jetTrkMaxDR[i][j] = new TH1D(Form("h_jetTrkMaxDR_C%iJ%i",i,j),Form("jetTrkMaxDR, %s",tag),NdRBins,dRBinMin,dRBinMax);
      h_jetTrkMaxPtRel[i][j] = new TH1D(Form("h_jetTrkMaxPtRel_C%iJ%i",i,j),Form("jetTrkMaxPtRel, %s",tag),NMuRelPtBins,muRelPtMin,muRelPtMax);

      h_jetTrkMaxPt[i][j]->Sumw2();
      h_jetTrkMaxPtOverJetPt[i][j]->Sumw2();
      h_jetTrkMaxEta[i][j]->Sumw2();
      h_jetTrkMaxPhi[i][j]->Sumw2();
      h_jetTrkMaxDR[i][j]->Sumw2();
      h_jetTrkMaxPtRel[i][j]->Sumw2();
    }
  }

  TFile *f = TFile::Open(input.c_str());
  if(!f || f->IsZombie()){ cout << "ERROR: cannot open " << input << endl; return; }
  auto em = new eventMap(f);
  em->isMC = 0;
  em->AASetup = 1;
  em->init();                       // hiEvtAnalyzer/HiTree
  em->loadJet(jtm::jetTreeName);    // (its trackMaxEta/Phi branches are absent here and unused)
  em->loadHLT("hltanalysis/HltTree");
  em->regEventFilter(jtm::NeventFilters, jtm::eventFilters);
  const Long64_t NEvents = em->evtTree->GetEntries();
  cout << "Number of events = " << NEvents << endl;

  const std::string trkErr = jtrk::openTracks(f, em, jtm::trkTreeName, NEvents);
  if(!trkErr.empty()){ cout << "ERROR: " << trkErr << " in " << input << endl; return; }
  cout << "track source: " << jtrk::trackSource() << "; cuts: " << jtrk::cutsApplied() << endl;

  long nEvt = 0, nJet = 0, nLead = 0, nPerp = 0;
  std::vector<int> sel;
  for(Long64_t evi = 0; evi < NEvents; evi++){

    em->getEvent(evi);

    // global event cuts, as PbPb_scan.C
    if(fabs(em->vz) > 15.0) continue;
    if(em->checkEventFilter()) continue;
    if(em->hiHF > 6000) continue;

    const int CentralityIndex = getCentBin(em->hiBin);
    if(CentralityIndex < 0) continue;
    nEvt++;

    const double w = 1.0;   // data

    // selected tracks hard enough for any histogram (jetTrackPtRel.h)
    if(!jtrk::selectTracks(em, evi, jtm::trkPtFloor, sel)) return;

    for(int i = 0; i < em->njet; i++){

      JEC.SetJetPT(em->rawpt[i]);
      JEC.SetJetEta(em->jeteta[i]);
      JEC.SetJetPhi(em->jetphi[i]);
      const double recoJetPt_i = JEC.GetCorrectedPT();
      const double recoJetEta_i = em->jeteta[i];
      const double recoJetPhi_i = em->jetphi[i];
      if(TMath::Abs(recoJetEta_i) > etaMax) continue;
      nJet++;

      h_inclRecoJetPt[0]->Fill(recoJetPt_i,w);
      h_inclRecoJetPt[CentralityIndex]->Fill(recoJetPt_i,w);
      if(sel.empty()) continue;

      // every selected track in the jet (ATLAS)
      for(int k : sel){
        if(jtrk::trkPt[k] < jtm::trkPtCut3D) continue;
        if(getDr(jtrk::trkEta[k], jtrk::trkPhi[k], recoJetEta_i, recoJetPhi_i) >= jtm::coneR) continue;
        const double pr = getPtRel(jtrk::trkPt[k], jtrk::trkEta[k], jtrk::trkPhi[k], recoJetPt_i, recoJetEta_i, recoJetPhi_i);
        h_jetTrkPtRel_trkPt_recoJetPt[0]->Fill(pr, jtrk::trkPt[k], recoJetPt_i, w);
        h_jetTrkPtRel_trkPt_recoJetPt[CentralityIndex]->Fill(pr, jtrk::trkPt[k], recoJetPt_i, w);
      }

      // leading track in the perpendicular cone (eta, phi + pi/2): underlying event only
      const double perpPhi = TVector2::Phi_mpi_pi(recoJetPhi_i + 0.5*TMath::Pi());
      const int kp = jtrk::leading(sel, recoJetEta_i, perpPhi, jtm::coneR);
      if(kp >= 0 && jtrk::trkPt[kp] >= jtm::trkPtCut3D){
        const double pr = getPtRel(jtrk::trkPt[kp], jtrk::trkEta[kp], jtrk::trkPhi[kp], recoJetPt_i, recoJetEta_i, perpPhi);
        h_perpTrkMaxPtRel_trkPt_recoJetPt[0]->Fill(pr, jtrk::trkPt[kp], recoJetPt_i, w);
        h_perpTrkMaxPtRel_trkPt_recoJetPt[CentralityIndex]->Fill(pr, jtrk::trkPt[kp], recoJetPt_i, w);
        nPerp++;
      }

      // leading track in the jet
      const int k = jtrk::leading(sel, recoJetEta_i, recoJetPhi_i, jtm::coneR);
      if(k < 0) continue;
      const double jetTrkMax_i = jtrk::trkPt[k];
      const double jetTrkMaxEta_i = jtrk::trkEta[k];
      const double jetTrkMaxPhi_i = jtrk::trkPhi[k];
      const double jetTrkMaxDR_i = getDr(jetTrkMaxEta_i, jetTrkMaxPhi_i, recoJetEta_i, recoJetPhi_i);
      const double jetTrkMaxPtRel_i = getPtRel(jetTrkMax_i, jetTrkMaxEta_i, jetTrkMaxPhi_i, recoJetPt_i, recoJetEta_i, recoJetPhi_i);

      if(jetTrkMax_i >= jtm::trkPtCut3D){
        h_jetTrkMaxPtRel_trkPt_recoJetPt[0]->Fill(jetTrkMaxPtRel_i, jetTrkMax_i, recoJetPt_i, w);
        h_jetTrkMaxPtRel_trkPt_recoJetPt[CentralityIndex]->Fill(jetTrkMaxPtRel_i, jetTrkMax_i, recoJetPt_i, w);
        nLead++;
      }

      // the 2025 1D set: jet pT > 80, track pT > 14
      if(recoJetPt_i < jtm::jetPtCut1D || jetTrkMax_i < jtm::trkMaxCut1D) continue;
      const int jetPtIndex = getJetPtBin(recoJetPt_i);
      if(jetPtIndex < 0) continue;

      h_jetPt[0]->Fill(recoJetPt_i,w);
      h_jetPt[CentralityIndex]->Fill(recoJetPt_i,w);
      for(int c : {0, CentralityIndex}){
        for(int j : {0, jetPtIndex}){
          h_jetTrkMaxPt[c][j]->Fill(jetTrkMax_i,w);
          h_jetTrkMaxPtOverJetPt[c][j]->Fill(jetTrkMax_i/recoJetPt_i,w);
          h_jetTrkMaxEta[c][j]->Fill(jetTrkMaxEta_i,w);
          h_jetTrkMaxPhi[c][j]->Fill(jetTrkMaxPhi_i,w);
          h_jetTrkMaxDR[c][j]->Fill(jetTrkMaxDR_i,w);
          h_jetTrkMaxPtRel[c][j]->Fill(jetTrkMaxPtRel_i,w);
        }
      }
    }
  }
  cout << "events " << nEvt << ", jets " << nJet << ", jets with a leading track > " << jtm::trkPtCut3D << " GeV " << nLead
       << ", perpendicular cones with one " << nPerp << endl;

  delete f;
  auto wf = TFile::Open(output,"recreate");
  if(!wf || wf->IsZombie()){ cout << "ERROR: cannot write " << output << endl; return; }

  const TString trkSource = jtrk::trackSource();
  TNamed prov("provenance", Form("macro: PbPb_jetTrkMax_scan.C\ngitHash: %s\ndirty: %s\ninput list: %s\ninput file: %s\njets: %s\ntracks: %s\ntrack cuts: %s\ncentrality: %s\nJEC: Autumn18_HI_V8 AK4PF L2Relative + L2L3Residual\nwritten: %s",
                                 jtm::shellLine("git rev-parse --short HEAD 2>/dev/null").c_str(),
                                 jtm::shellLine("git status --porcelain --untracked-files=no 2>/dev/null | head -1 | grep -q . && echo yes || echo no").c_str(),
                                 jtm::inputFileList, input.c_str(), jtm::jetTreeName, trkSource.Data(), jtrk::cutsApplied().c_str(),
                                 CENT_SCHEME_SUFFIX, TDatime().AsSQLString()));
  prov.Write();

  for(int i = 0; i < NCentralityIndices; i++){
    h_jetTrkMaxPtRel_trkPt_recoJetPt[i]->Write();
    h_jetTrkPtRel_trkPt_recoJetPt[i]->Write();
    h_perpTrkMaxPtRel_trkPt_recoJetPt[i]->Write();
    h_inclRecoJetPt[i]->Write();
    h_jetPt[i]->Write();
    for(int j = 0; j < NJetPtIndices; j++){
      h_jetTrkMaxPt[i][j]->Write();
      h_jetTrkMaxPtOverJetPt[i][j]->Write();
      h_jetTrkMaxEta[i][j]->Write();
      h_jetTrkMaxPhi[i][j]->Write();
      h_jetTrkMaxDR[i][j]->Write();
      h_jetTrkMaxPtRel[i][j]->Write();
    }
  }
  wf->Close();
}
