// Leading-track ptRel in PbPb jet-triggered data, per ultraFine centrality slice, for
// a data-driven light template (ATLAS-style, arXiv:2204.13530: track-jet pairs in
// inclusive-jet data, track pT reweighted to the light-jet muon pT). The PbPb
// counterpart of src/scanning/pp/pp_jetTrkMax_scan.C. 2026-10-19 (user).
//
// Input: the 2025-01-31 HIHardProbes forest that stores the jet trackMax DIRECTION
// (trackMaxEta / trackMaxPhi; the production forests store only its pT), read
// directly, one forest file per index of
//   fileNames/fileNames_PbPb_DiJet_onlyJets_withTrackMaxInfo_allFiles_partial.txt
// The 2025-08-04 version read skims of the shorter _partial list, and those skims
// name their trees evtTree / jetTree / hltTree, which eventMap no longer opens.
//
// Histograms, per class C<i> (C0 = every slice):
//   h_jetTrkMaxPtRel_trkPt_recoJetPt_C<i>         leading track ptRel x track pT x jet pT,
//                                                  Jet60, every jet with |eta| < 1.6
//   h_jetTrkMaxPtRel_trkPt_recoJetPt_jet80_C<i>   the same, Jet80: trigger-bias check
//   h_jetTrkMax{Pt,PtOverJetPt,Eta,Phi,DR,PtRel}_C<i>J<j>   the 2025 1D set (Jet60,
//                                                  jet pT > 80, track pT > 14 GeV)
// The 3D uses the muon-template binning (headers/AnalysisSetup/ptRelMuPt3D.h, track
// pT from 15 GeV, the muon threshold), the same as the data's
// h_muptrel_mupt_recoJetPt_inclRecoMuonTag_triggerOn_C<i> in PbPb_scan.C. Inclusive
// jets in a Jet60 sample are ~97% light, so it is a light template up to that b + c
// share.
//
// The jet axis is the same akCs4PF axis the muon ptRel is measured against, so the
// template carries the underlying event's smearing of that axis -- the effect the
// PYTHIA+HYDJET templates have to model (the central b-purity degeneracy). The
// forest's trackMax is the highest-pT track in the jet, so in central events it can
// be an underlying-event track; above 15 GeV that is rare, and comparing slices shows
// how rare.
//
// JEC: Autumn18_HI_V8 AK4PF L2Relative + L2L3Residual, as PbPb_scan.C. Event filters:
// the five PbPb_scan.C applies to every sample. hiBin: half-open ultraFine classes
// (getCentBinImpl.h), the current convention.
//
// Run per forest file, from src/scanning/PbPb/:
//   root -l -b -q 'PbPb_jetTrkMax_scan.C(<index>)'   1-based line of the list
//   root -l -b -q 'PbPb_jetTrkMax_scan.C(0)'         query mode: list and output dir
// Condor: condor_PbPb_jetTrkMax_scan.py. Merge with hadd.
#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH3D.h"
#include "TMath.h"
#include "TNamed.h"
#include "TString.h"
#include "TSystem.h"
#include "TDatime.h"
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
// getJetPtBin function
#include "../../../headers/functions/getJetPtBin.h"
// getPtRel function
#include "../../../headers/functions/getPtRel.h"


namespace jtm {
const char *inputFileList = "../../../fileNames/fileNames_PbPb_DiJet_onlyJets_withTrackMaxInfo_allFiles_partial.txt";
const char *outputBaseDir = "/eos/cms/store/group/phys_heavyions/cbennett/scanningOutput/";
const char *outputName = "output_PbPb_HardProbes_jetTrkMax_trkPtAxis_jet60";
// as PbPb_scan.C, which applies the SingleMuon list to every sample
const int NeventFilters = 5;
std::string eventFilters[NeventFilters] = {"pprimaryVertexFilter", "HBHENoiseFilterResultRun2Loose", "collisionEventSelectionAODv2", "phfCoincFilter2Th4", "pclusterCompatibilityFilter"};
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


TH3D *h_jetTrkMaxPtRel_trkPt_recoJetPt[NCentralityIndices];         // Jet60
TH3D *h_jetTrkMaxPtRel_trkPt_recoJetPt_jet80[NCentralityIndices];   // Jet80
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

    h_jetTrkMaxPtRel_trkPt_recoJetPt[i] = bookPtRelMuPt3D(Form("h_jetTrkMaxPtRel_trkPt_recoJetPt_C%i",i), Form("leading track p_{T}^{rel} x track p_{T} x jet p_{T}, Jet60, hiBin %i - %i",cLo,cHi));
    h_jetTrkMaxPtRel_trkPt_recoJetPt_jet80[i] = bookPtRelMuPt3D(Form("h_jetTrkMaxPtRel_trkPt_recoJetPt_jet80_C%i",i), Form("leading track p_{T}^{rel} x track p_{T} x jet p_{T}, Jet80, hiBin %i - %i",cLo,cHi));

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
  em->init();                                   // hiEvtAnalyzer/HiTree
  em->loadJet("akCs4PFJetAnalyzer/t");          // trackMax, trackMaxEta, trackMaxPhi, trackMaxDR
  em->loadHLT("hltanalysis/HltTree");
  em->regEventFilter(jtm::NeventFilters, jtm::eventFilters);
  const int NEvents = em->evtTree->GetEntries();
  cout << "Number of events = " << NEvents << endl;

  long nJet60 = 0, nJet80 = 0, nFilled3D = 0;
  for(int evi = 0; evi < NEvents; evi++){

    em->getEvent(evi);

    // global event cuts, as PbPb_scan.C
    if(fabs(em->vz) > 15.0) continue;
    if(em->checkEventFilter()) continue;
    if(em->hiHF > 6000) continue;

    const bool jet60 = (em->HLT_HICsAK4PFJet60Eta1p5_v1 != 0);
    const bool jet80 = (em->HLT_HICsAK4PFJet80Eta1p5_v1 != 0);
    if(!jet60 && !jet80) continue;
    if(jet60) nJet60++;
    if(jet80) nJet80++;

    const int CentralityIndex = getCentBin(em->hiBin);
    if(CentralityIndex < 0) continue;

    const double w = 1.0;   // data

    for(int i = 0; i < em->njet; i++){

      JEC.SetJetPT(em->rawpt[i]);
      JEC.SetJetEta(em->jeteta[i]);
      JEC.SetJetPhi(em->jetphi[i]);
      const double recoJetPt_i = JEC.GetCorrectedPT();
      const double recoJetEta_i = em->jeteta[i];
      const double recoJetPhi_i = em->jetphi[i];
      const double jetTrkMax_i = em->jetTrkMax[i];
      const double jetTrkMaxEta_i = em->jetTrkMaxEta[i];
      const double jetTrkMaxPhi_i = em->jetTrkMaxPhi[i];
      const double jetTrkMaxDR_i = em->jetTrkMaxDR[i];

      if(TMath::Abs(recoJetEta_i) > etaMax) continue;
      if(jetTrkMax_i <= 0.) continue;   // no track in the jet

      const double jetTrkMaxPtRel_i = getPtRel(jetTrkMax_i,jetTrkMaxEta_i,jetTrkMaxPhi_i,recoJetPt_i,recoJetEta_i,recoJetPhi_i);

      // 3D: every jet; track pT and jet pT windows are chosen downstream from the axes
      if(jet60){
        h_jetTrkMaxPtRel_trkPt_recoJetPt[0]->Fill(jetTrkMaxPtRel_i, jetTrkMax_i, recoJetPt_i, w);
        h_jetTrkMaxPtRel_trkPt_recoJetPt[CentralityIndex]->Fill(jetTrkMaxPtRel_i, jetTrkMax_i, recoJetPt_i, w);
        nFilled3D++;
      }
      if(jet80){
        h_jetTrkMaxPtRel_trkPt_recoJetPt_jet80[0]->Fill(jetTrkMaxPtRel_i, jetTrkMax_i, recoJetPt_i, w);
        h_jetTrkMaxPtRel_trkPt_recoJetPt_jet80[CentralityIndex]->Fill(jetTrkMaxPtRel_i, jetTrkMax_i, recoJetPt_i, w);
      }

      // the 2025 1D set: Jet60, jet pT > 80, track pT > 14
      if(!jet60) continue;
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
  cout << "events: Jet60 " << nJet60 << ", Jet80 " << nJet80 << "; Jet60 jets in the 3D: " << nFilled3D << endl;

  delete f;
  auto wf = TFile::Open(output,"recreate");
  if(!wf || wf->IsZombie()){ cout << "ERROR: cannot write " << output << endl; return; }

  TNamed prov("provenance", Form("macro: PbPb_jetTrkMax_scan.C\ngitHash: %s\ndirty: %s\ninput list: %s\ninput file: %s\ncentrality: %s\nJEC: Autumn18_HI_V8 AK4PF L2Relative + L2L3Residual\nwritten: %s",
                                 jtm::shellLine("git rev-parse --short HEAD 2>/dev/null").c_str(),
                                 jtm::shellLine("git status --porcelain --untracked-files=no 2>/dev/null | head -1 | grep -q . && echo yes || echo no").c_str(),
                                 jtm::inputFileList, input.c_str(), CENT_SCHEME_SUFFIX, TDatime().AsSQLString()));
  prov.Write();

  for(int i = 0; i < NCentralityIndices; i++){
    h_jetTrkMaxPtRel_trkPt_recoJetPt[i]->Write();
    h_jetTrkMaxPtRel_trkPt_recoJetPt_jet80[i]->Write();
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
