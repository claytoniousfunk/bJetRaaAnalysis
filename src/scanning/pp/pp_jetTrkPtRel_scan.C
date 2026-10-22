// Track ptRel in pp jet60 data for a data-driven light template, ATLAS-style
// (arXiv:2204.13530: ALL track-jet pairs in inclusive-jet data, track pT reweighted
// to the light-jet muon pT). 2026-10-19 (user). The leading-track version
// (pp_jetTrkMax_scan.C, the 2025 skims) biases the jets: requiring a 20 GeV track to
// be the hardest one selects broad jets with no harder charged particle, so its
// ptRel exceeds the light-muon ptRel even in PYTHIA (+0.07 to +0.22 GeV, rising with
// jet pT), and no track-pT reweighting removes that. Every track above threshold
// carries no such condition. The leading track is filled too, for the comparison.
//
// Input: the HighEGJet forests pp_scan.C reads, one file per index of
//   fileNames/fileNames_pp_HighEGJet.txt
// Tracks (headers/functions/jetTrackPtRel.h, shared with PbPb_jetTrkMax_scan.C):
// ppTrack/trackTree with the selection there if the forest has it, else the
// charged PF candidates of pfcandAnalyzer/pfTree; the provenance says which.
//
// Histograms, x = ptRel, y = track pT, z = jet pT on the muon-template binning
// (headers/AnalysisSetup/ptRelMuPt3D.h, track pT from 15 GeV):
//   h_jetTrkPtRel_trkPt_recoJetPt       every selected track within dR < 0.4 (ATLAS)
//   h_jetTrkMaxPtRel_trkPt_recoJetPt    leading selected track within dR < 0.4
//   h_perpTrkMaxPtRel_trkPt_recoJetPt   leading selected track in a dR < 0.4 cone at
//                                       the jet's eta, phi + pi/2 (underlying event
//                                       and pileup), ptRel against that axis
//   h_inclRecoJetPt                     every jet, |eta| < 1.6
// ak4PF jets, JEC Spring18 ppRef V6 L2Relative + L2L3Residual (as
// pp_scan.C), |eta| < 1.6. Events: HLT_HIAK4PFJet60, |vz| < 15, the filters
// pp_scan.C applies to HighEGJet.
//
// Run per forest file, from src/scanning/pp/:
//   root -l -b -q 'pp_jetTrkPtRel_scan.C(<index>)'   1-based line of the list
//   root -l -b -q 'pp_jetTrkPtRel_scan.C(0)'         query mode: list and output dir
// Condor: condor_pp_jetTrkPtRel_scan.py. Merge with hadd.
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
#include "TVector2.h"
#include <vector>
#include <string>
#include <cstdio>

// event map
#include "../../../eventMap/eventMap.h"
// jet corrector
#include "../../../JetEnergyCorrections/JetCorrector.h"
// analysis constants (etaMax, ptMin, ptMax, NPtBins)
#include "../../../headers/AnalysisSetup/common.h"
// ptRel x (track) pT x jet pT binning, shared with the muon templates
#include "../../../headers/AnalysisSetup/ptRelMuPt3D.h"
// getDr function
#include "../../../headers/functions/getDr.h"
// getPtRel function
#include "../../../headers/functions/getPtRel.h"
// track source, selection and matching, shared with PbPb_jetTrkMax_scan.C
#include "../../../headers/functions/jetTrackPtRel.h"


namespace jpp {
const char *inputFileList = "../../../fileNames/fileNames_pp_HighEGJet.txt";
const char *outputBaseDir = "/eos/cms/store/group/phys_heavyions/cbennett/scanningOutput/";
const char *outputName = "output_pp_HighEGJet_jet60_trkPtRel_trkPtAxis";
const char *jetTreeName = "ak4PFJetAnalyzer/t";
const char *trkTreeName = "ppTrack/trackTree";
// as pp_scan.C for HighEGJet (config_pp.h eventFilters_SingleMuon)
const int NeventFilters = 3;
std::string eventFilters[NeventFilters] = {"pPAprimaryVertexFilter", "HBHENoiseFilterResultRun2Loose", "pBeamScrapingFilter"};
const double coneR = 0.4;          // track-jet match
const double trkPtCut = 15.;       // the 3D axis starts here

TString outputDir(){
  // SCAN_DATE (set by myProcesses/condor/scan_condor.py) pins the date to the submission
  TString date;
  const char *scanDate = gSystem->Getenv("SCAN_DATE");
  if(scanDate && scanDate[0]) date = scanDate;
  else { TDatime dt; date = Form("%i-%i-%i", dt.GetYear(), dt.GetMonth(), dt.GetDay()); }
  return Form("%s%s_%s", outputBaseDir, outputName, date.Data());
}

std::string shellLine(const char *cmd){
  std::string out; FILE *p = popen(cmd, "r"); if(!p) return "unavailable";
  char buf[512]; if(fgets(buf, sizeof(buf), p)) out = buf; pclose(p);
  while(!out.empty() && (out.back() == '\n' || out.back() == '\r')) out.pop_back();
  return out.empty() ? "unavailable" : out;
}
}


void pp_jetTrkPtRel_scan(int group = 1){

  const TString outputDir = jpp::outputDir();
  if(group == 0){ // query mode for condor submit scripts
    cout << "INPUTFILELIST=" << jpp::inputFileList << endl;
    cout << "OUTPUTDIR=" << outputDir << endl;
    return;
  }

  std::ifstream instr(jpp::inputFileList, std::ifstream::in);
  if(!instr.is_open()){ cout << "ERROR: file list not found: " << jpp::inputFileList << endl; return; }
  std::string filename, input;
  int ifile = 0;
  while(instr >> filename){ ifile++; if(ifile == group){ input = filename; break; } }
  if(input.empty()){ cout << "ERROR: no line " << group << " in " << jpp::inputFileList << endl; return; }
  if(gSystem->AccessPathName(outputDir)){ cout << "ERROR: output directory not found: " << outputDir << endl; return; }
  const TString output = Form("%s/pp_scan_output_%i.root", outputDir.Data(), group);
  cout << "input  = " << input << endl << "output = " << output << endl;

  // JET ENERGY CORRECTIONS (as pp_scan.C, PF)
  vector<string> Files;
  Files.push_back("../../../JetEnergyCorrections/Spring18_ppRef5TeV_V6_DATA_L2Relative_AK4PF.txt");
  Files.push_back("../../../JetEnergyCorrections/Spring18_ppRef5TeV_V6_DATA_L2L3Residual_AK4PF.txt");
  JetCorrector JEC(Files);

  TH3D *hAll  = bookPtRelMuPt3D("h_jetTrkPtRel_trkPt_recoJetPt", "every track in jet: p_{T}^{rel} x track p_{T} x jet p_{T}");
  TH3D *hLead = bookPtRelMuPt3D("h_jetTrkMaxPtRel_trkPt_recoJetPt", "leading track in jet: p_{T}^{rel} x track p_{T} x jet p_{T}");
  TH3D *hPerp = bookPtRelMuPt3D("h_perpTrkMaxPtRel_trkPt_recoJetPt", "leading track in perpendicular cone: p_{T}^{rel} x track p_{T} x jet p_{T}");
  TH1D *hJet  = new TH1D("h_inclRecoJetPt", Form("every jet, |#eta| < %.1f", etaMax), NPtBins, ptMin, ptMax); hJet->Sumw2();

  TFile *f = TFile::Open(input.c_str());
  if(!f || f->IsZombie()){ cout << "ERROR: cannot open " << input << endl; return; }
  auto em = new eventMap(f);
  em->isMC = 0;
  em->AASetup = 0;
  em->init();                       // hiEvtAnalyzer/HiTree
  em->loadJet(jpp::jetTreeName);    // (its trackMaxEta/Phi branches may be absent; unused)
  em->loadHLT("hltanalysis/HltTree");
  em->regEventFilter(jpp::NeventFilters, jpp::eventFilters);
  const Long64_t NEvents = em->evtTree->GetEntries();
  cout << "Number of events = " << NEvents << endl;

  const std::string trkErr = jtrk::openTracks(f, em, jpp::trkTreeName, NEvents);
  if(!trkErr.empty()){ cout << "ERROR: " << trkErr << " in " << input << endl; return; }
  cout << "track source: " << jtrk::trackSource() << "; cuts: " << jtrk::cutsApplied() << endl;

  long nEvt = 0, nJet = 0, nAll = 0, nLead = 0, nPerp = 0;
  std::vector<int> sel;
  for(Long64_t evi = 0; evi < NEvents; evi++){

    em->getEvent(evi);
    if(fabs(em->vz) > 15.0) continue;
    if(em->checkEventFilter()) continue;
    if(em->HLT_HIAK4PFJet60_v1 == 0) continue;
    nEvt++;

    const double w = 1.0;   // data
    if(!jtrk::selectTracks(em, evi, jpp::trkPtCut, sel)) return;

    for(int i = 0; i < em->njet; i++){

      JEC.SetJetPT(em->rawpt[i]);
      JEC.SetJetEta(em->jeteta[i]);
      JEC.SetJetPhi(em->jetphi[i]);
      const double jetPt = JEC.GetCorrectedPT(), jetEta = em->jeteta[i], jetPhi = em->jetphi[i];
      if(TMath::Abs(jetEta) > etaMax) continue;
      nJet++;
      hJet->Fill(jetPt, w);
      if(sel.empty()) continue;

      // every selected track in the jet (ATLAS)
      for(int k : sel){
        if(getDr(jtrk::trkEta[k], jtrk::trkPhi[k], jetEta, jetPhi) >= jpp::coneR) continue;
        hAll->Fill(getPtRel(jtrk::trkPt[k], jtrk::trkEta[k], jtrk::trkPhi[k], jetPt, jetEta, jetPhi), jtrk::trkPt[k], jetPt, w);
        nAll++;
      }
      // leading selected track in the jet
      const int kl = jtrk::leading(sel, jetEta, jetPhi, jpp::coneR);
      if(kl >= 0){ hLead->Fill(getPtRel(jtrk::trkPt[kl], jtrk::trkEta[kl], jtrk::trkPhi[kl], jetPt, jetEta, jetPhi), jtrk::trkPt[kl], jetPt, w); nLead++; }
      // leading selected track in the perpendicular cone
      const double perpPhi = TVector2::Phi_mpi_pi(jetPhi + 0.5*TMath::Pi());
      const int kp = jtrk::leading(sel, jetEta, perpPhi, jpp::coneR);
      if(kp >= 0){ hPerp->Fill(getPtRel(jtrk::trkPt[kp], jtrk::trkEta[kp], jtrk::trkPhi[kp], jetPt, jetEta, perpPhi), jtrk::trkPt[kp], jetPt, w); nPerp++; }
    }
  }
  cout << "events " << nEvt << ", jets " << nJet << ", tracks > " << jpp::trkPtCut << " GeV in jets " << nAll
       << ", jets with a leading one " << nLead << ", perpendicular cones with one " << nPerp << endl;

  delete f;
  auto wf = TFile::Open(output,"recreate");
  if(!wf || wf->IsZombie()){ cout << "ERROR: cannot write " << output << endl; return; }
  const TString trkSource = jtrk::trackSource();
  TNamed prov("provenance", Form("macro: pp_jetTrkPtRel_scan.C\ngitHash: %s\ndirty: %s\ninput list: %s\ninput file: %s\njets: %s\ntracks: %s\ntrack cuts: %s\ntrigger: HLT_HIAK4PFJet60_v1\nJEC: Spring18_ppRef5TeV_V6 AK4PF L2Relative + L2L3Residual\nwritten: %s",
                                 jpp::shellLine("git rev-parse --short HEAD 2>/dev/null").c_str(),
                                 jpp::shellLine("git status --porcelain --untracked-files=no 2>/dev/null | head -1 | grep -q . && echo yes || echo no").c_str(),
                                 jpp::inputFileList, input.c_str(), jpp::jetTreeName, trkSource.Data(), jtrk::cutsApplied().c_str(), TDatime().AsSQLString()));
  prov.Write();
  hAll->Write(); hLead->Write(); hPerp->Write(); hJet->Write();
  wf->Close();
}
