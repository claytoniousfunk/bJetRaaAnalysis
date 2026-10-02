// Z -> mu mu tag-and-probe scan for the muon reco + tight-ID efficiency in PbPb.
//
// The PbPb counterpart of src/scanning/pp/pp_muonTagAndProbe_scan.C. Same steps,
// same tree layout, same mass windows, so the two outputs can go through one
// efficiency macro. Read that file's header first; only the differences are
// described here.
//
// WHY PbPb NEEDS ITS OWN MEASUREMENT. The analysis applies MC efficiencies of
// 0.8627 / 0.9069 / 0.9856 / 0.9778 for 0-10 / 10-30 / 30-50 / 50-80%, from
// PYTHIA+HYDJET. The 0-10% value is 11% below pp, and it is the largest
// centrality-dependent correction in the R_AA after the jets themselves. Whether
// the data show that drop is exactly what this measures.
//
// DATASET. HISingleMuon by default -- the analysis's own PbPb sample, and far
// more Z than HIHardProbes. That costs the clean trigger argument the pp scan
// has: the forests carry no hltobject tree, so the tag cannot be matched to the
// muon that fired HLT_HIL3Mu12, and some events are in the sample only because
// the PROBE fired it. Those probes are reconstructed L3 muons, so they pass more
// often, and the measured efficiency is biased UP.
//
// The bias is bounded. requireTrigger keeps only events with the trigger bit set
// (every SingleMuon event should have one of its paths, not necessarily this
// one), and the tag is tight with pT > 20, on the trigger plateau. With t the
// per-muon trigger efficiency and R = N_pass / N_fail, assuming failing probes
// never fire the trigger (the worst case):
//     eps_true / (1 - eps_true)  >=  R / (2 - t)
// i.e. the measured pass/fail odds are inflated by at most (2 - t), and the
// inefficiency is low by a fraction ~(1 - t). At eps = 0.87 that shifts the
// efficiency up by 0.005 for t = 0.95, and by 0.015 if t falls to 0.85 in central
// events. muonTagAndProbeEfficiency.C prints this bound. The clean cross-check
// is the same scan on HIHardProbes (jet/photon triggers, no muon paths): set
// requireTrigger = false and point condor_PbPb_muonTagAndProbe.py at that list.
//
// CENTRALITY. hiBin is stored per pair. Data use it as is; MC is shifted by
// hiBinShift = 10 (config_PYTHIAHYDJET.h), the analysis convention, and the tree
// stores both. Histograms are booked per class:
//   C0 0-90%   C1 0-10%   C2 10-30%   C3 30-50%   C4 50-80%   C5 80-90%
// C0 is its own fill, not the sum of the others (events above 90% are dropped).
//
// TAG. Tight ID and pT > tagPtMin as in pp, but NO isolation cut. PF isolation in
// PbPb is dominated by the underlying event and would remove tags preferentially
// in central events -- harmless for the probe's efficiency but a needless loss of
// statistics exactly where they are scarcest. tagRelIso is still stored. The
// background this lets in is handled by the same-sign and sideband estimate in
// the efficiency macro; check it per class, since it grows toward central.
//
// EVENT FILTERS. The analysis's PbPb SingleMuon set (config_PbPb.h).
//
// MC. The PYTHIA+HYDJET dijet forests have essentially no Z; on them only the
// h_gen* truth histograms are useful, now per centrality class. No centrality
// reweighting is applied to MC -- within a class the truth ratio barely depends
// on it, but do not integrate across classes without it.
//
// OUTPUT
//   tnp                                 TTree, as in pp plus hiBin, centHiBin, trigFired
//   h_<step>_<pass|fail>_<OS|SS>_C<k>   TH3D mass x probe pT x probe |eta|
//   h_gen<All|FromZ>_<level>_C<k>       (MC) TH2D gen pT x |eta|
//   h_evtCount, h_hiBin, provenance
//
// Usage (from this directory):
//   root -l -b -q 'PbPb_muonTagAndProbe_scan.C+("files.txt", "out.root", false)'
//   root -l -b -q 'PbPb_muonTagAndProbe_scan.C+("HiForestAOD.root", "out_MC.root", true)'

#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TNamed.h"
#include "TString.h"
#include "TSystem.h"
#include "TVector2.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include <vector>
#include <string>
#include <fstream>
#include <cstdio>
#include <cmath>

#include "../../../headers/functions/isQualityMuon_tight.h"

namespace tnpAA {

const double vzMax        = 15.0;
const double muEtaMax     = 2.4;
const double tagPtMin     = 20.0;
const double probePtMin   = 10.0;
const double pairDRMin    = 0.3;
const double dupDR        = 0.3;
const double genMatchDR   = 0.1;
const double massMuon     = 0.1056583755;
const int    hiBinShiftMC = 10;   // config_PYTHIAHYDJET.h

// the analysis's PbPb muon trigger (config_PbPb.h, fillMu12). See the header
// for why it is required and what it costs. Set false for HIHardProbes.
const char  *trigName       = "HLT_HIL3Mu12_v1";
const bool   requireTrigger = true;

const std::vector<std::string> filters = {"pprimaryVertexFilter",
                                          "HBHENoiseFilterResultRun2Loose",
                                          "collisionEventSelectionAODv2",
                                          "phfCoincFilter2Th4",
                                          "pclusterCompatibilityFilter"};

// centrality classes in hiBin (0.5% each); C0 is the inclusive 0-90%
const int  nCls = 6;
const int  clsLo[nCls] = {0,   0, 20, 60, 100, 160};
const int  clsHi[nCls] = {180, 20, 60, 100, 160, 180};
const char *clsName[nCls] = {"0-90%", "0-10%", "10-30%", "30-50%", "50-80%", "80-90%"};

const std::vector<double> ptEdges  = {10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 100, 200};
const std::vector<double> etaEdges = {0., 0.3, 0.6, 0.9, 1.2, 1.6, 2.0, 2.1, 2.4};
const int    nMass = 100;
const double massLo = 40., massHi = 140.;

const char *stepName[3] = {"trk", "id", "all"};

double dR(double e1, double p1, double e2, double p2)
{
  double dp = TVector2::Phi_mpi_pi(p1 - p2);
  return std::sqrt((e1 - e2)*(e1 - e2) + dp*dp);
}

bool tightFloat(float chi2, float d0, float dz, int muHits, int pixHits,
                int isGlobal, int isPF, int stations, int trkLayers)
{
  return chi2 >= 0 && chi2 < 10. && std::fabs(d0) < 0.2 && std::fabs(dz) < 0.5 &&
         muHits > 0 && pixHits > 0 && isGlobal && isPF && stations > 1 && trkLayers > 5;
}

} // namespace tnpAA

void PbPb_muonTagAndProbe_scan(TString input, TString output, bool isMC = false,
                               long maxEvents = -1)
{
  using namespace tnpAA;

  TChain evt("hiEvtAnalyzer/HiTree"), mu("ggHiNtuplizerGED/EventTree"),
         skim("skimanalysis/HltTree"), hlt("hltanalysis/HltTree");
  gSystem->ExpandPathName(input);
  if(input.EndsWith(".txt")){
    std::ifstream in(input.Data());
    if(!in){ printf("ERROR: cannot open file list %s\n", input.Data()); return; }
    std::string f;
    while(in >> f){ evt.Add(f.c_str()); mu.Add(f.c_str()); skim.Add(f.c_str()); hlt.Add(f.c_str()); }
  }
  else{ evt.Add(input); mu.Add(input); skim.Add(input); hlt.Add(input); }
  if(evt.GetNtrees() == 0 || evt.GetEntries() == 0){
    printf("ERROR: no events in %s (%d files)\n", input.Data(), evt.GetNtrees()); return; }

  // A forest without the GED muon block (several PbPb productions lack it)
  // would otherwise run with null vectors and crash on the first event.
  for(const char *b : {"muPt", "muIsSTA", "muInnerPt", "muPFChIso", "muIDTight"})
    if(!mu.GetBranch(b)){
      printf("ERROR: ggHiNtuplizerGED/EventTree has no %s in %s -- wrong forest?\n", b, input.Data());
      return;
    }
  if(!evt.GetBranch("hiBin")){ printf("ERROR: no hiBin in hiEvtAnalyzer/HiTree\n"); return; }
  // without the bit, a SingleMuon scan cannot bound its trigger bias; in MC the
  // trigger only matters if the MC is meant to mimic it, so a warning suffices
  const bool haveTrig = (hlt.GetBranch(trigName) != nullptr);
  if(requireTrigger && !haveTrig){
    if(!isMC){ printf("ERROR: %s not in hltanalysis/HltTree, cannot require it\n", trigName); return; }
    printf("  WARNING: %s absent in this MC, trigger NOT required\n", trigName);
  }
  evt.AddFriend(&mu); evt.AddFriend(&skim);
  if(haveTrig) evt.AddFriend(&hlt);

  // ---- branches -------------------------------------------------------------
  Float_t vz = 0, weight = 1; Int_t hiBin = -1;
  UInt_t run = 0, lumi = 0; ULong64_t event = 0;
  std::vector<float> *muPt = 0, *muEta = 0, *muPhi = 0, *muInnerPt = 0, *muChi2 = 0,
                     *muD0 = 0, *muDz = 0, *chIso = 0, *phoIso = 0, *neuIso = 0, *puIso = 0;
  std::vector<int>   *muCharge = 0, *isGlobal = 0, *isTracker = 0, *isSTA = 0, *isPF = 0,
                     *muHits = 0, *pixHits = 0, *stations = 0, *trkLayers = 0, *idTight = 0;
  std::vector<int>   *mcPID = 0, *mcStatus = 0, *mcMomPID = 0;
  std::vector<float> *mcPt = 0, *mcEta = 0, *mcPhi = 0;

  evt.SetBranchStatus("*", 0);
  auto on = [&](const char *b, void *addr){
    evt.SetBranchStatus(b, 1); evt.SetBranchAddress(b, addr); };
  on("vz", &vz); on("run", &run); on("lumi", &lumi); on("evt", &event); on("hiBin", &hiBin);
  if(isMC) on("weight", &weight);
  Int_t trigBit = 1;
  if(haveTrig) on(trigName, &trigBit);
  const bool applyTrig = requireTrigger && haveTrig;
  on("muPt", &muPt);           on("muEta", &muEta);         on("muPhi", &muPhi);
  on("muInnerPt", &muInnerPt); on("muCharge", &muCharge);
  on("muChi2NDF", &muChi2);    on("muInnerD0", &muD0);      on("muInnerDz", &muDz);
  on("muMuonHits", &muHits);   on("muPixelHits", &pixHits);
  on("muIsGlobal", &isGlobal); on("muIsTracker", &isTracker); on("muIsSTA", &isSTA);
  on("muIsPF", &isPF);         on("muStations", &stations); on("muTrkLayers", &trkLayers);
  on("muIDTight", &idTight);
  on("muPFChIso", &chIso);     on("muPFPhoIso", &phoIso);
  on("muPFNeuIso", &neuIso);   on("muPFPUIso", &puIso);
  if(isMC){
    on("mcPID", &mcPID); on("mcStatus", &mcStatus); on("mcMomPID", &mcMomPID);
    on("mcPt", &mcPt);   on("mcEta", &mcEta);       on("mcPhi", &mcPhi);
  }

  std::vector<Int_t> fval(filters.size(), 1);
  for(size_t i = 0; i < filters.size(); i++){
    if(!skim.GetBranch(filters[i].c_str())){
      printf("  filter %s absent, NOT applied\n", filters[i].c_str()); continue; }
    on(filters[i].c_str(), &fval[i]);
  }

  // ---- output ---------------------------------------------------------------
  TFile *wf = TFile::Open(output, "RECREATE");
  if(!wf || wf->IsZombie()){ printf("cannot open %s\n", output.Data()); return; }

  TTree *tree = new TTree("tnp", "Z->mumu tag-and-probe pairs, PbPb");
  UInt_t t_run, t_lumi; ULong64_t t_evt; Float_t t_w;
  Int_t   t_hiBin, centHiBin, t_trig;
  Float_t tagPt, tagEta, tagPhi, tagRelIso, mass;
  Int_t   tagCharge, nTags, pairOS;
  Float_t pPt, pEta, pPhi, pInnerPt, pRelIso, pDRtag, pMinDRTrkMu;
  Int_t   pCharge, pIsGlobal, pIsTracker, pIsSTA, pIsPF, pHasTrk;
  Int_t   pTightProd, pTightFloat, pIDTightCMSSW;
  Int_t   pGenMatch = -1, pGenMomPID = 0; Float_t pGenPt = -1;
  tree->Branch("run", &t_run);           tree->Branch("lumi", &t_lumi);
  tree->Branch("evt", &t_evt);           tree->Branch("w", &t_w);
  tree->Branch("hiBin", &t_hiBin);       // as stored in the forest
  tree->Branch("centHiBin", &centHiBin); // analysis convention: hiBin - 10 in MC
  tree->Branch("trigFired", &t_trig);    // trigName; always 1 when requireTrigger
  tree->Branch("nTags", &nTags);
  tree->Branch("tagPt", &tagPt);         tree->Branch("tagEta", &tagEta);
  tree->Branch("tagPhi", &tagPhi);       tree->Branch("tagCharge", &tagCharge);
  tree->Branch("tagRelIso", &tagRelIso);
  tree->Branch("mass", &mass);           tree->Branch("pairOS", &pairOS);
  tree->Branch("pPt", &pPt);             tree->Branch("pEta", &pEta);
  tree->Branch("pPhi", &pPhi);           tree->Branch("pInnerPt", &pInnerPt);
  tree->Branch("pCharge", &pCharge);     tree->Branch("pRelIso", &pRelIso);
  tree->Branch("pDRtag", &pDRtag);       tree->Branch("pMinDRTrkMu", &pMinDRTrkMu);
  tree->Branch("pIsGlobal", &pIsGlobal); tree->Branch("pIsTracker", &pIsTracker);
  tree->Branch("pIsSTA", &pIsSTA);       tree->Branch("pIsPF", &pIsPF);
  tree->Branch("pHasTrk", &pHasTrk);
  tree->Branch("pTightProd", &pTightProd);
  tree->Branch("pTightFloat", &pTightFloat);
  tree->Branch("pIDTightCMSSW", &pIDTightCMSSW);
  if(isMC){
    tree->Branch("pGenMatch", &pGenMatch);
    tree->Branch("pGenPt", &pGenPt);
    tree->Branch("pGenMomPID", &pGenMomPID);
  }

  std::vector<double> mEdges(nMass + 1);
  for(int i = 0; i <= nMass; i++) mEdges[i] = massLo + i*(massHi - massLo)/nMass;

  TH3D *h[3][2][2][nCls];   // step, fail/pass, OS/SS, class
  for(int s = 0; s < 3; s++) for(int p = 0; p < 2; p++) for(int q = 0; q < 2; q++)
    for(int c = 0; c < nCls; c++){
      h[s][p][q][c] = new TH3D(Form("h_%s_%s_%s_C%d", stepName[s], p ? "pass" : "fail", q ? "SS" : "OS", c),
                               Form("%s step, %s probes, %s, %s; m_{#mu#mu} [GeV]; probe p_{T} [GeV]; probe |#eta|",
                                    stepName[s], p ? "passing" : "failing", q ? "same-sign" : "opposite-sign", clsName[c]),
                               nMass, mEdges.data(), ptEdges.size() - 1, ptEdges.data(),
                               etaEdges.size() - 1, etaEdges.data());
      h[s][p][q][c]->Sumw2();
    }

  const char *genLevel[5] = {"den", "anyMu", "STA", "hasTrk", "tight"};
  TH2D *hGen[2][5][nCls] = {};
  if(isMC)
    for(int g = 0; g < 2; g++) for(int l = 0; l < 5; l++) for(int c = 0; c < nCls; c++){
      hGen[g][l][c] = new TH2D(Form("h_gen%s_%s_C%d", g ? "FromZ" : "All", genLevel[l], c),
                               Form("gen muons%s, %s, %s; gen p_{T} [GeV]; gen |#eta|",
                                    g ? " from Z" : "", genLevel[l], clsName[c]),
                               ptEdges.size() - 1, ptEdges.data(), etaEdges.size() - 1, etaEdges.data());
      hGen[g][l][c]->Sumw2();
    }

  const int nCut = 8;
  const char *cutName[nCut] = {"read", "|vz|<15", "filters", "0-90%", "trigger", "#geq2 #mu", "#geq1 tag", "#geq1 pair"};
  TH1D *hCount = new TH1D("h_evtCount", "event cut flow", nCut, 0, nCut);
  for(int i = 0; i < nCut; i++) hCount->GetXaxis()->SetBinLabel(i + 1, cutName[i]);
  TH1D *hHiBin = new TH1D("h_hiBin", "centHiBin, events passing filters; centHiBin", 200, 0, 200);

  // ---- event loop -----------------------------------------------------------
  long nEvt = evt.GetEntries();
  if(maxEvents > 0 && maxEvents < nEvt) nEvt = maxEvents;
  printf("  %ld events, %s\n", nEvt, isMC ? "MC" : "data");
  long nPairs = 0;

  for(long i = 0; i < nEvt; i++){
    if(i % 200000 == 0) printf("  event %ld / %ld\n", i, nEvt);
    evt.GetEntry(i);
    const double w = isMC ? weight : 1.;
    hCount->Fill(0., w);
    if(std::fabs(vz) > vzMax) continue;
    hCount->Fill(1., w);
    bool bad = false; for(int v : fval) if(!v) bad = true;
    if(bad) continue;
    hCount->Fill(2., w);

    const int cb = isMC ? hiBin - hiBinShiftMC : hiBin;
    hHiBin->Fill(cb, w);
    if(cb < clsLo[0] || cb >= clsHi[0]) continue;
    hCount->Fill(3., w);
    int cls = -1;
    for(int c = 1; c < nCls; c++) if(cb >= clsLo[c] && cb < clsHi[c]) cls = c;
    const int fillCls[2] = {0, cls};

    const int n = muPt->size();
    auto hasTrk = [&](int j){ return isGlobal->at(j) || isTracker->at(j); };
    auto relIso = [&](int j){
      double neu = neuIso->at(j) + phoIso->at(j) - 0.5*puIso->at(j);
      return (chIso->at(j) + (neu > 0 ? neu : 0.)) / muPt->at(j); };
    auto tightProd = [&](int j){
      return isQualityMuon_tight(muChi2->at(j), muD0->at(j), muDz->at(j), muHits->at(j),
                                 pixHits->at(j), isGlobal->at(j), isPF->at(j),
                                 stations->at(j), trkLayers->at(j)); };

    std::vector<int> genOf(n, -1);
    if(isMC){
      for(int j = 0; j < n; j++){
        double best = genMatchDR;
        for(size_t g = 0; g < mcPID->size(); g++){
          if(std::abs(mcPID->at(g)) != 13 || mcStatus->at(g) != 1) continue;
          double d = dR(muEta->at(j), muPhi->at(j), mcEta->at(g), mcPhi->at(g));
          if(d < best){ best = d; genOf[j] = g; }
        }
      }
      for(size_t g = 0; g < mcPID->size(); g++){
        if(std::abs(mcPID->at(g)) != 13 || mcStatus->at(g) != 1) continue;
        if(std::fabs(mcEta->at(g)) > muEtaMax || mcPt->at(g) < ptEdges.front()) continue;
        bool lv[5] = {true, false, false, false, false};
        for(int j = 0; j < n; j++){
          if(genOf[j] != (int) g) continue;
          lv[1] = true;
          if(isSTA->at(j)) lv[2] = true;
          if(hasTrk(j))    lv[3] = true;
          if(tightProd(j)) lv[4] = true;
        }
        const bool fromZ = (mcMomPID->at(g) == 23);
        for(int l = 0; l < 5; l++){
          if(!lv[l]) continue;
          for(int c : fillCls){
            if(c < 0) continue;
            hGen[0][l][c]->Fill(mcPt->at(g), std::fabs(mcEta->at(g)), w);
            if(fromZ) hGen[1][l][c]->Fill(mcPt->at(g), std::fabs(mcEta->at(g)), w);
          }
        }
      }
    }

    // the trigger gates the tag-and-probe only, not the MC truth above: truth is
    // the efficiency of every gen muon, and gating it on a muon trigger would
    // select events whose muons were already reconstructed
    if(applyTrig && !trigBit) continue;
    hCount->Fill(4., w);
    if(n < 2) continue;
    hCount->Fill(5., w);

    std::vector<int> tags;
    for(int j = 0; j < n; j++)
      if(muPt->at(j) > tagPtMin && std::fabs(muEta->at(j)) < muEtaMax && tightProd(j))
        tags.push_back(j);
    if(tags.empty()) continue;
    hCount->Fill(6., w);

    bool anyPair = false;
    for(int t : tags){
      TLorentzVector vt; vt.SetPtEtaPhiM(muPt->at(t), muEta->at(t), muPhi->at(t), massMuon);
      for(int p = 0; p < n; p++){
        if(p == t) continue;
        if(muPt->at(p) < probePtMin || std::fabs(muEta->at(p)) > muEtaMax) continue;
        const double dtp = dR(muEta->at(t), muPhi->at(t), muEta->at(p), muPhi->at(p));
        if(dtp < pairDRMin) continue;
        TLorentzVector vp; vp.SetPtEtaPhiM(muPt->at(p), muEta->at(p), muPhi->at(p), massMuon);
        const double m = (vt + vp).M();
        if(m < massLo || m > massHi) continue;

        double minDR = 99.;
        for(int k = 0; k < n; k++){
          if(k == t || k == p || !hasTrk(k)) continue;
          double d = dR(muEta->at(p), muPhi->at(p), muEta->at(k), muPhi->at(k));
          if(d < minDR) minDR = d;
        }

        t_run = run; t_lumi = lumi; t_evt = event; t_w = w;
        t_hiBin = hiBin; centHiBin = cb; t_trig = haveTrig ? trigBit : -1;
        nTags = tags.size();
        tagPt = muPt->at(t); tagEta = muEta->at(t); tagPhi = muPhi->at(t);
        tagCharge = muCharge->at(t); tagRelIso = relIso(t);
        mass = m; pairOS = (muCharge->at(t)*muCharge->at(p) < 0);
        pPt = muPt->at(p); pEta = muEta->at(p); pPhi = muPhi->at(p);
        pInnerPt = muInnerPt->at(p); pCharge = muCharge->at(p); pRelIso = relIso(p);
        pDRtag = dtp; pMinDRTrkMu = minDR;
        pIsGlobal = isGlobal->at(p); pIsTracker = isTracker->at(p);
        pIsSTA = isSTA->at(p); pIsPF = isPF->at(p); pHasTrk = hasTrk(p);
        pTightProd = tightProd(p);
        pTightFloat = tightFloat(muChi2->at(p), muD0->at(p), muDz->at(p), muHits->at(p),
                                 pixHits->at(p), isGlobal->at(p), isPF->at(p),
                                 stations->at(p), trkLayers->at(p));
        pIDTightCMSSW = idTight->at(p);
        if(isMC){
          pGenMatch = genOf[p];
          pGenPt = pGenMatch >= 0 ? mcPt->at(pGenMatch) : -1.;
          pGenMomPID = pGenMatch >= 0 ? mcMomPID->at(pGenMatch) : 0;
        }
        tree->Fill();
        nPairs++; anyPair = true;

        if(minDR < dupDR) continue;
        const int q = pairOS ? 0 : 1;
        const double ae = std::fabs(pEta);
        for(int c : fillCls){
          if(c < 0) continue;
          if(pIsSTA)  h[0][pHasTrk   ][q][c]->Fill(m, pPt, ae, w);
          if(pHasTrk) h[1][pTightProd][q][c]->Fill(m, pPt, ae, w);
                      h[2][pTightProd][q][c]->Fill(m, pPt, ae, w);
        }
      }
    }
    if(anyPair) hCount->Fill(7., w);
  }

  // ---- summary --------------------------------------------------------------
  printf("\n  cut flow:\n");
  for(int i = 1; i <= nCut; i++)
    printf("    %-12s %12.0f\n", cutName[i - 1], hCount->GetBinContent(i));
  printf("  pairs written: %ld\n", nPairs);
  printf("  opposite-sign, 81 < m < 101 GeV, probe pT > 15 (raw counts, no fit):\n");
  for(int c = 0; c < nCls; c++)
    for(int s = 0; s < 3; s++){
      auto win = [&](TH3D *x){
        return x->Integral(x->GetXaxis()->FindBin(81. + 1e-6), x->GetXaxis()->FindBin(101. - 1e-6),
                           x->GetYaxis()->FindBin(15. + 1e-6), x->GetNbinsY(), 1, x->GetNbinsZ()); };
      double P = win(h[s][1][0][c]), F = win(h[s][0][0][c]);
      printf("    %-7s %-4s pass %8.0f  fail %6.0f  eff %.4f\n", s ? "" : clsName[c], stepName[s],
             P, F, P + F > 0 ? P/(P + F) : 0.);
    }

  // ---- provenance -----------------------------------------------------------
  auto shell = [](const char *cmd){
    std::string out; FILE *p = popen(cmd, "r"); if(!p) return std::string("unavailable");
    char buf[512]; if(fgets(buf, sizeof(buf), p)) out = buf; pclose(p);
    while(!out.empty() && (out.back() == '\n' || out.back() == '\r')) out.pop_back();
    return out.empty() ? std::string("unavailable") : out; };
  // shell() maps empty output to "unavailable", so a clean tree reads that way
  const std::string hash  = shell("git rev-parse --short HEAD 2>/dev/null");
  const std::string dirty = shell("git status --porcelain 2>/dev/null | head -1");
  TString prov = Form("scan: PbPb_muonTagAndProbe_scan.C\ngitHash: %s\ndirty: %s\n"
                      "input: %s\nisMC: %d\nhiBinShiftMC: %d\ntrigger: %s %s\ntagPtMin: %g\ntagIsolation: none\n"
                      "probePtMin: %g\nmuEtaMax: %g\npairDRMin: %g\ndupDR: %g\nvzMax: %g\nevents: %ld\n",
                      hash.c_str(),
                      hash == "unavailable" ? "unknown" : (dirty == "unavailable" ? "no" : "yes"),
                      input.Data(), isMC, isMC ? hiBinShiftMC : 0, trigName,
                      applyTrig ? "required" : (haveTrig ? "stored, not required" : "absent"),
                      tagPtMin, probePtMin, muEtaMax,
                      pairDRMin, dupDR, vzMax, nEvt);

  wf->cd();
  TNamed("provenance", prov.Data()).Write();
  tree->Write();
  hCount->Write();
  hHiBin->Write();
  for(int s = 0; s < 3; s++) for(int p = 0; p < 2; p++) for(int q = 0; q < 2; q++)
    for(int c = 0; c < nCls; c++) h[s][p][q][c]->Write();
  if(isMC) for(int g = 0; g < 2; g++) for(int l = 0; l < 5; l++) for(int c = 0; c < nCls; c++)
    hGen[g][l][c]->Write();
  wf->Close();
  printf("  wrote %s\n", output.Data());
}
