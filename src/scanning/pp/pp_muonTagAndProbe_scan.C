// Z -> mu mu tag-and-probe scan for the muon reco + tight-ID efficiency in pp.
//
// WHY. The efficiency the analysis applies (0.9708 for pp) is pure MC: gen muons
// matched to tight reco muons in PYTHIA. Nothing checks that against data. This
// scan writes what a data-driven measurement needs, and runs unchanged on MC so
// the two can be compared as a data/MC scale factor.
//
// INPUT. The HiForest directly, not the skims -- the skims drop muIsSTA and the
// isolation branches, and both are needed here. Trees read:
//   hiEvtAnalyzer/HiTree        vz, run/lumi/evt, weight (MC)
//   ggHiNtuplizerGED/EventTree  the muon block (ggHiNtuplizer holds no muons in
//                               the pp forests) and, in MC, the mc* block
//   skimanalysis/HltTree        event filters
//
// NO TRIGGER MATCHING, AND WHY THAT IS FINE HERE. The pp forests carry no
// hltobject tree, so the tag cannot be matched to the object that fired the
// trigger. In HighEGJet the triggers are photon/jet paths, which do not use
// muons at all, so the probe is not biased by the trigger and the tag needs no
// matching. This does NOT hold for SingleMuon: there the probe may be the muon
// that fired the event, which biases the passing fraction upward. Do not run
// this on SingleMuon without solving that first.
//
// TAG. Analysis tight ID (isQualityMuon_tight, the production definition with
// its int chi2), pT > tagPtMin, |eta| < 2.4, Delta-beta relative PF isolation
// < tagRelIsoMax. Every muon passing is used as a tag, so an event with two
// tags contributes two pairs -- the standard, unbiased choice.
//
// PROBES AND STEPS. There is no general-track collection in these forests, so
// the usual tracker-track probe -- the only probe that is independent of the
// muon system -- is not available. What is available is the reco::Muon
// collection, which holds global, tracker and standalone (STA) muons merged
// into one object per muon. Three steps are written:
//
//   trk  probe = any STA muon            pass = has an inner track
//                                               (isGlobal || isTracker)
//   id   probe = has an inner track      pass = analysis tight ID
//   all  probe = any muon object         pass = analysis tight ID
//
// trk x id is the factorized efficiency; `all` is the same thing in one step and
// serves as a consistency check on the factorization.
//
// MC NEEDS A Z SAMPLE. The dijet PYTHIA forests contain essentially no Z, so the
// MC side of the scale factor needs a Drell-Yan / Z->mumu pp sample. On the dijet
// sample only the h_gen* truth histograms are useful.
//
// WHAT NONE OF THESE MEASURE. Every probe already has a signature in the muon
// system (an STA track or matched segments). A muon leaving no muon-system
// signature at all is never a probe, so eps(any muon object | true muon) is not
// measured. In MC the gen-level histograms below give it directly; use them to
// size that factor.
//
// ALSO NOT MEASURED: muons inside jets. Z muons are isolated. The analysis
// tags muons in b jets, so this checks the ID and reconstruction on clean
// muons only; the in-jet environment is a separate question.
//
// DUPLICATES. An STA fragment can sit next to a muon that does have a track and
// appear as a second, failing probe. pMinDRTrkMu stores the distance to the
// nearest other inner-track muon (tag and probe excluded). The histograms drop
// probes with pMinDRTrkMu < dupDR; the tree keeps everything so the cut can be
// varied without a rescan.
//
// MASS. Computed from muPt/muEta/muPhi. For an STA-only probe that is the
// standalone momentum, whose resolution is poor -- hence the wide 40-140 window.
// Fit the passing and failing mass distributions separately (signal + background);
// the same-sign histograms give a background shape.
//
// OUTPUT
//   tnp                        TTree, one row per (tag, probe) pair
//   h_<step>_<pass|fail>_<OS|SS>  TH3D mass x probe pT x probe |eta|
//   h_gen_*  (MC only)         TH2D gen pT x |eta| for gen muons, prompt Z muons
//                              separately, with what each one was reconstructed as
//   h_evtCount                 cut flow (weighted in MC)
//   provenance                 TNamed
//
// Usage (from this directory):
//   root -l -b -q 'pp_muonTagAndProbe_scan.C+("files.txt", "out.root", false)'
//   root -l -b -q 'pp_muonTagAndProbe_scan.C+("HiForestAOD.root", "out_MC.root", true)'
// A path ending in .txt is read as a list of forest files, one per line.

#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TNamed.h"
#include "TString.h"
#include "TVector2.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include <vector>
#include <string>
#include <fstream>
#include <cstdio>
#include <cmath>

#include "../../../headers/functions/isQualityMuon_tight.h"

namespace tnp {

const double vzMax        = 15.0;
const double muEtaMax     = 2.4;
const double tagPtMin     = 20.0;
const double tagRelIsoMax = 0.15;
const double probePtMin   = 10.0;
const double pairDRMin    = 0.3;   // tag and probe must be distinct muons
const double dupDR        = 0.3;   // histogram-level duplicate veto, see header
const double genMatchDR   = 0.1;
const double massMuon     = 0.1056583755;

const std::vector<std::string> filters = {"pPAprimaryVertexFilter",
                                          "HBHENoiseFilterResultRun2Loose",
                                          "pBeamScrapingFilter"};

// mass is 1 GeV over 40-140; pT and |eta| edges are deliberately finer than any
// one plot will use, so they can be merged later. Every edge here is a boundary
// a later Rebin may need.
const std::vector<double> ptEdges  = {10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 100, 200};
const std::vector<double> etaEdges = {0., 0.3, 0.6, 0.9, 1.2, 1.6, 2.1, 2.4};
const int    nMass = 100;
const double massLo = 40., massHi = 140.;

const char *stepName[3] = {"trk", "id", "all"};

double dR(double e1, double p1, double e2, double p2)
{
  double dp = TVector2::Phi_mpi_pi(p1 - p2);
  return std::sqrt((e1 - e2)*(e1 - e2) + dp*dp);
}

// the CMS tight muon with a FLOAT chi2/ndf. isQualityMuon_tight takes chi2 as
// an int, so 10 <= chi2 < 11 passes there; both are stored.
bool tightFloat(float chi2, float d0, float dz, int muHits, int pixHits,
                int isGlobal, int isPF, int stations, int trkLayers)
{
  return chi2 >= 0 && chi2 < 10. && std::fabs(d0) < 0.2 && std::fabs(dz) < 0.5 &&
         muHits > 0 && pixHits > 0 && isGlobal && isPF && stations > 1 && trkLayers > 5;
}

} // namespace tnp

void pp_muonTagAndProbe_scan(TString input, TString output, bool isMC = false,
                             long maxEvents = -1)
{
  using namespace tnp;

  TChain evt("hiEvtAnalyzer/HiTree"), mu("ggHiNtuplizerGED/EventTree"),
         skim("skimanalysis/HltTree");
  if(input.EndsWith(".txt")){
    std::ifstream in(input.Data()); std::string f;
    while(in >> f){ evt.Add(f.c_str()); mu.Add(f.c_str()); skim.Add(f.c_str()); }
  }
  else{ evt.Add(input); mu.Add(input); skim.Add(input); }
  evt.AddFriend(&mu); evt.AddFriend(&skim);

  // ---- branches -------------------------------------------------------------
  Float_t vz = 0, weight = 1;
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
  on("vz", &vz); on("run", &run); on("lumi", &lumi); on("evt", &event);
  if(isMC) on("weight", &weight);
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

  TTree *tree = new TTree("tnp", "Z->mumu tag-and-probe pairs");
  UInt_t t_run, t_lumi; ULong64_t t_evt; Float_t t_w;
  Float_t tagPt, tagEta, tagPhi, tagRelIso, mass;
  Int_t   tagCharge, nTags, pairOS;
  Float_t pPt, pEta, pPhi, pInnerPt, pRelIso, pDRtag, pMinDRTrkMu;
  Int_t   pCharge, pIsGlobal, pIsTracker, pIsSTA, pIsPF, pHasTrk;
  Int_t   pTightProd, pTightFloat, pIDTightCMSSW;
  Int_t   pGenMatch = -1, pGenMomPID = 0; Float_t pGenPt = -1;
  tree->Branch("run", &t_run);           tree->Branch("lumi", &t_lumi);
  tree->Branch("evt", &t_evt);           tree->Branch("w", &t_w);
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
  tree->Branch("pTightProd", &pTightProd);       // analysis definition, int chi2
  tree->Branch("pTightFloat", &pTightFloat);     // same, float chi2
  tree->Branch("pIDTightCMSSW", &pIDTightCMSSW); // muIDTight from the forest
  if(isMC){
    tree->Branch("pGenMatch", &pGenMatch);       // index of matched gen muon, -1 if none
    tree->Branch("pGenPt", &pGenPt);
    tree->Branch("pGenMomPID", &pGenMomPID);
  }

  TH3D *h[3][2][2];   // step, fail/pass, OS/SS
  for(int s = 0; s < 3; s++)
    for(int p = 0; p < 2; p++)
      for(int q = 0; q < 2; q++){
        std::vector<double> mEdges(nMass + 1);
        for(int i = 0; i <= nMass; i++) mEdges[i] = massLo + i*(massHi - massLo)/nMass;
        h[s][p][q] = new TH3D(Form("h_%s_%s_%s", stepName[s], p ? "pass" : "fail", q ? "SS" : "OS"),
                              Form("%s step, %s probes, %s; m_{#mu#mu} [GeV]; probe p_{T} [GeV]; probe |#eta|",
                                   stepName[s], p ? "passing" : "failing", q ? "same-sign" : "opposite-sign"),
                              nMass, mEdges.data(), ptEdges.size() - 1, ptEdges.data(),
                              etaEdges.size() - 1, etaEdges.data());
        h[s][p][q]->Sumw2();
      }

  // MC truth: what each gen muon in acceptance was reconstructed as. [0] all
  // status-1 muons, [1] those with mcMomPID == 23 (the population the T&P sees)
  const char *genLevel[5] = {"den", "anyMu", "STA", "hasTrk", "tight"};
  TH2D *hGen[2][5] = {};
  if(isMC)
    for(int g = 0; g < 2; g++)
      for(int l = 0; l < 5; l++){
        hGen[g][l] = new TH2D(Form("h_gen%s_%s", g ? "FromZ" : "All", genLevel[l]),
                              Form("gen muons%s, %s; gen p_{T} [GeV]; gen |#eta|",
                                   g ? " from Z" : "", genLevel[l]),
                              ptEdges.size() - 1, ptEdges.data(), etaEdges.size() - 1, etaEdges.data());
        hGen[g][l]->Sumw2();
      }

  const char *cutName[6] = {"read", "|vz|<15", "filters", "#geq2 #mu", "#geq1 tag", "#geq1 pair"};
  TH1D *hCount = new TH1D("h_evtCount", "event cut flow", 6, 0, 6);
  for(int i = 0; i < 6; i++) hCount->GetXaxis()->SetBinLabel(i + 1, cutName[i]);

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

    const int n = muPt->size();

    auto hasTrk = [&](int j){ return isGlobal->at(j) || isTracker->at(j); };
    auto relIso = [&](int j){
      double neu = neuIso->at(j) + phoIso->at(j) - 0.5*puIso->at(j);
      return (chIso->at(j) + (neu > 0 ? neu : 0.)) / muPt->at(j); };
    auto tightProd = [&](int j){
      return isQualityMuon_tight(muChi2->at(j), muD0->at(j), muDz->at(j), muHits->at(j),
                                 pixHits->at(j), isGlobal->at(j), isPF->at(j),
                                 stations->at(j), trkLayers->at(j)); };

    // gen matching, closest gen muon within genMatchDR
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
      // truth: for each gen muon, the best reco level among the muons matched to it
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
          hGen[0][l]->Fill(mcPt->at(g), std::fabs(mcEta->at(g)), w);
          if(fromZ) hGen[1][l]->Fill(mcPt->at(g), std::fabs(mcEta->at(g)), w);
        }
      }
    }

    if(n < 2) continue;
    hCount->Fill(3., w);

    std::vector<int> tags;
    for(int j = 0; j < n; j++)
      if(muPt->at(j) > tagPtMin && std::fabs(muEta->at(j)) < muEtaMax &&
         tightProd(j) && relIso(j) < tagRelIsoMax)
        tags.push_back(j);
    if(tags.empty()) continue;
    hCount->Fill(4., w);

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
        if(pIsSTA)  h[0][pHasTrk   ][q]->Fill(m, pPt, ae, w);
        if(pHasTrk) h[1][pTightProd][q]->Fill(m, pPt, ae, w);
                    h[2][pTightProd][q]->Fill(m, pPt, ae, w);
      }
    }
    if(anyPair) hCount->Fill(5., w);
  }

  // ---- summary --------------------------------------------------------------
  printf("\n  cut flow:\n");
  for(int i = 1; i <= 6; i++)
    printf("    %-12s %12.0f\n", cutName[i - 1], hCount->GetBinContent(i));
  printf("  pairs written: %ld\n", nPairs);
  printf("  opposite-sign, 81 < m < 101 GeV, probe pT > 15 (raw counts, no fit):\n");
  for(int s = 0; s < 3; s++){
    auto win = [&](TH3D *x){
      return x->Integral(x->GetXaxis()->FindBin(81. + 1e-6), x->GetXaxis()->FindBin(101. - 1e-6),
                         x->GetYaxis()->FindBin(15. + 1e-6), x->GetNbinsY(), 1, x->GetNbinsZ()); };
    double P = win(h[s][1][0]), F = win(h[s][0][0]);
    printf("    %-4s pass %8.0f  fail %6.0f  eff %.4f\n", stepName[s], P, F, P + F > 0 ? P/(P + F) : 0.);
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
  TString prov = Form("scan: pp_muonTagAndProbe_scan.C\ngitHash: %s\ndirty: %s\n"
                      "input: %s\nisMC: %d\ntagPtMin: %g\ntagRelIsoMax: %g\nprobePtMin: %g\n"
                      "muEtaMax: %g\npairDRMin: %g\ndupDR: %g\nvzMax: %g\nevents: %ld\n",
                      hash.c_str(),
                      hash == "unavailable" ? "unknown" : (dirty == "unavailable" ? "no" : "yes"),
                      input.Data(), isMC, tagPtMin, tagRelIsoMax, probePtMin, muEtaMax,
                      pairDRMin, dupDR, vzMax, nEvt);

  wf->cd();
  TNamed("provenance", prov.Data()).Write();
  tree->Write();
  hCount->Write();
  for(int s = 0; s < 3; s++) for(int p = 0; p < 2; p++) for(int q = 0; q < 2; q++) h[s][p][q]->Write();
  if(isMC) for(int g = 0; g < 2; g++) for(int l = 0; l < 5; l++) hGen[g][l]->Write();
  wf->Close();
  printf("  wrote %s\n", output.Data());
}
