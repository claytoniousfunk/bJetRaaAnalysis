// Calo vs PF jet energy scale, event by event, in PbPb and pp data.
//
// WHY. With T_AA normalization the calo-jet R_AA from calculateJetsPerZ.cc sits
// 25-35% below the published inclusive-jet R_AA in 50-80% even above 200 GeV,
// where the PbPb spectrum is pure Jet100 (no stitching), fakes are negligible
// and the N_evt / T_AA / L_pp inputs all check out against the data. With a
// spectrum falling like pT^-5, a ~6-7% lower calo scale in PbPb than in pp would
// produce that. The PYTHIA vs PYTHIA+HYDJET responses differ by only 2-3% and
// the unfolding removes that, so the suspect is a scale difference that exists
// in DATA only: the akPu4Calo pileup subtraction and the manual calo JEC
// behaving differently on real PbPb than on HYDJET.
//
// WHAT. In each event, every calo jet is matched to the PF jet of the same
// event, both with the manual JEC the analysis applies (scan_jet_corrections /
// PbPb_scan.C, pp_scan.C). PF is the reference: it is the collection the
// published CMS R_AA and the analysis's PF scans use. The observable is
//     double ratio = <calo/PF>_PbPb / <calo/PF>_pp
// in bins of jet pT, per centrality class. A double ratio below 1 is a calo
// scale that is lower in PbPb than in pp; R_AA moves by roughly (ratio)^(n-1)
// for a local spectral index n. The PF scale itself cancels in the double
// ratio as long as PF behaves the same in PbPb and pp -- which is itself not
// guaranteed (CS subtraction), so this tests the calo-vs-PF consistency, not
// either scale absolutely.
//
// COLLECTIONS
//   PbPb  akPu4CaloJetAnalyzer/t  vs  akCs4PFJetAnalyzer/t   (Autumn18_HI_V8 AK4Calo / AK4PF)
//   pp    ak4CaloJetAnalyzer/t    vs  ak4PFJetAnalyzer/t     (Spring18_ppRef5TeV_V6 AK4Calo / AK4PF)
// Same JEC files as PbPb_scan.C / pp_scan.C with useManualJEC.
//
// MATCHING. Jets with |eta| < etaMax (common.h, 1.6) and pT > ptMinMatch.
// Greedy one-to-one by smallest dR, dR < matchDR. Both conventions are stored:
//   calo/PF  vs PF pT    and    PF/calo  vs calo pT
// Binning in one collection's pT biases the ratio of a steep spectrum toward
// that collection (the "resolution bias"); a genuine scale shift moves the two
// conventions in opposite directions, a resolution difference moves them the
// same way. Compare the two before reading a number off either.
// The ratio is also stored for raw pT (no JEC: isolates the background
// subtraction) and for forest jtpt (the forest's own JEC, which for calo jets
// is the PF L2Relative -- see the memory note on forest calo jets).
//
// EVENTS. The analysis filters (config_PbPb.h HardProbes scan uses the
// SingleMuon set; config_pp.h HighEGJet likewise), |vz| < 15, and the OR of
// the PF and calo Jet100 triggers, so that neither collection is preferred at
// selection:
//   PbPb  HLT_HICsAK4PFJet100Eta1p5_v1 || HLT_HIPuAK4CaloJet100Eta5p1_v1
//   pp    HLT_HIAK4PFJet100_v1         || HLT_HIAK4CaloJet100_v1
// Both triggers are fully efficient well below 200 GeV; read results above
// ~150 GeV. Each bit is also stored per event in the cut-flow histogram.
//
// CENTRALITY (PbPb; hiBin in 0.5% steps). C0 is its own fill, 0-90%.
//   C0 0-90%  C1 0-10%  C2 10-30%  C3 30-50%  C4 50-80%  C5 80-90%
// pp fills C0 only.
//
// OUTPUT, per class k:
//   h_caloOverPF_vsPFPt_C<k>          TH2D  PF pT x calo/PF          (manual JEC)
//   h_PFOverCalo_vsCaloPt_C<k>        TH2D  calo pT x PF/calo        (manual JEC)
//   h_rawCaloOverRawPF_vsPFPt_C<k>    TH2D  PF pT x raw calo/raw PF
//   h_forestCaloOverPF_vsPFPt_C<k>    TH2D  PF pT x forest-jtpt calo/PF
//   h_pfPt_C<k>, h_pfPt_matched_C<k>  TH1D  PF jets, all / with a calo match
//   h_caloPt_C<k>, h_caloPt_matched_C<k>
//   h_matchDR_C<k>                    TH1D
//   h_evtCount, h_hiBin, provenance
//
// Usage (from this directory):
//   root -l -b -q 'caloVsPFJetScale_scan.C+("files.txt", "out_PbPb.root", false)'
//   root -l -b -q 'caloVsPFJetScale_scan.C+("HiForestAOD.root", "out_pp.root", true)'
//   root -l -b -q 'caloVsPFJetScale_scan.C+("mb.root", "out.root", false, -1, false)'  # no trigger
// Batch: condor_caloVsPFJetScale.py.

#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TNamed.h"
#include "TString.h"
#include "TSystem.h"
#include "TVector2.h"
#include "TMath.h"
#include <vector>
#include <string>
#include <fstream>
#include <cstdio>
#include <cmath>

#include "../../../JetEnergyCorrections/JetCorrector.h"

namespace cvpf {

const double vzMax       = 15.0;
const double etaMax      = 1.6;    // common.h
const double ptMinMatch  = 30.0;   // both collections, after JEC
const double matchDR     = 0.2;
const int    jetMax      = 500;

// centrality classes in hiBin (0.5% each); C0 is the inclusive 0-90%
const int  nCls = 6;
const int  clsLo[nCls] = {0,   0, 20, 60, 100, 160};
const int  clsHi[nCls] = {180, 20, 60, 100, 160, 180};
const char *clsName[nCls] = {"0-90%", "0-10%", "10-30%", "30-50%", "50-80%", "80-90%"};

// pT edges: fine enough to follow the trigger region, coarse above 300
const std::vector<double> ptEdges = {30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140,
                                     150, 160, 180, 200, 220, 250, 280, 320, 400, 500, 700, 1000};
const int    nRatio = 300;
const double ratioLo = 0., ratioHi = 3.;

const std::vector<std::string> filtersPbPb = {"pprimaryVertexFilter",
                                              "HBHENoiseFilterResultRun2Loose",
                                              "collisionEventSelectionAODv2",
                                              "phfCoincFilter2Th4",
                                              "pclusterCompatibilityFilter"};
const std::vector<std::string> filtersPP   = {"pPAprimaryVertexFilter",
                                              "HBHENoiseFilterResultRun2Loose",
                                              "pBeamScrapingFilter"};

const char *trigPFPbPb   = "HLT_HICsAK4PFJet100Eta1p5_v1";
const char *trigCaloPbPb = "HLT_HIPuAK4CaloJet100Eta5p1_v1";
const char *trigPFPP     = "HLT_HIAK4PFJet100_v1";
const char *trigCaloPP   = "HLT_HIAK4CaloJet100_v1";

struct Jets {
  Int_t   n = 0;
  Float_t raw[jetMax], pt[jetMax], eta[jetMax], phi[jetMax];
};

double dR(double e1, double p1, double e2, double p2)
{
  double dp = TVector2::Phi_mpi_pi(p1 - p2);
  return std::sqrt((e1 - e2)*(e1 - e2) + dp*dp);
}

} // namespace cvpf

// requireTrigger = false keeps every event (MinBias forests, or to see the
// scale below the trigger turn-on); the trigger bits are still counted.
void caloVsPFJetScale_scan(TString input, TString output, bool isPP = false,
                           long maxEvents = -1, bool requireTrigger = true)
{
  using namespace cvpf;

  const char *caloTree = isPP ? "ak4CaloJetAnalyzer/t" : "akPu4CaloJetAnalyzer/t";
  const char *pfTree   = isPP ? "ak4PFJetAnalyzer/t"   : "akCs4PFJetAnalyzer/t";
  const char *trigPF   = isPP ? trigPFPP   : trigPFPbPb;
  const char *trigCalo = isPP ? trigCaloPP : trigCaloPbPb;
  const std::vector<std::string> &filters = isPP ? filtersPP : filtersPbPb;

  TChain evt("hiEvtAnalyzer/HiTree"), skim("skimanalysis/HltTree"), hlt("hltanalysis/HltTree"),
         calo(caloTree), pf(pfTree);
  gSystem->ExpandPathName(input);
  auto addAll = [&](const char *f){ evt.Add(f); skim.Add(f); hlt.Add(f); calo.Add(f); pf.Add(f); };
  if(input.EndsWith(".txt")){
    std::ifstream in(input.Data());
    if(!in){ printf("ERROR: cannot open file list %s\n", input.Data()); return; }
    std::string f;
    while(in >> f) addAll(f.c_str());
  }
  else addAll(input.Data());
  if(evt.GetNtrees() == 0 || evt.GetEntries() == 0){
    printf("ERROR: no events in %s (%d files)\n", input.Data(), evt.GetNtrees()); return; }

  // both collections must be in the forest; a forest with only one would
  // otherwise run with an empty chain and silently match nothing
  for(TChain *c : {&calo, &pf})
    if(!c->GetBranch("rawpt") || !c->GetBranch("jtpt")){
      printf("ERROR: %s has no rawpt/jtpt in %s -- wrong forest?\n", c->GetName(), input.Data());
      return;
    }
  if(calo.GetEntries() != evt.GetEntries() || pf.GetEntries() != evt.GetEntries()){
    printf("ERROR: entry mismatch: evt %lld, calo %lld, pf %lld\n",
           evt.GetEntries(), calo.GetEntries(), pf.GetEntries());
    return;
  }
  const bool havePF   = (hlt.GetBranch(trigPF)   != nullptr);
  const bool haveCalo = (hlt.GetBranch(trigCalo) != nullptr);
  if(requireTrigger && !havePF && !haveCalo){ printf("ERROR: neither %s nor %s in hltanalysis/HltTree\n", trigPF, trigCalo); return; }
  if(!havePF)   printf("  WARNING: %s absent, selecting on %s alone\n", trigPF, trigCalo);
  if(!haveCalo) printf("  WARNING: %s absent, selecting on %s alone\n", trigCalo, trigPF);

  // ---- branches -------------------------------------------------------------
  Float_t vz = 0; Int_t hiBin = -1;
  evt.SetBranchStatus("*", 0);
  evt.SetBranchStatus("vz", 1);    evt.SetBranchAddress("vz", &vz);
  if(!isPP){
    if(!evt.GetBranch("hiBin")){ printf("ERROR: no hiBin in hiEvtAnalyzer/HiTree\n"); return; }
    evt.SetBranchStatus("hiBin", 1); evt.SetBranchAddress("hiBin", &hiBin);
  }

  Int_t bitPF = 0, bitCalo = 0;
  hlt.SetBranchStatus("*", 0);
  if(havePF)  { hlt.SetBranchStatus(trigPF, 1);   hlt.SetBranchAddress(trigPF, &bitPF); }
  if(haveCalo){ hlt.SetBranchStatus(trigCalo, 1); hlt.SetBranchAddress(trigCalo, &bitCalo); }

  std::vector<Int_t> fval(filters.size(), 1);
  skim.SetBranchStatus("*", 0);
  for(size_t i = 0; i < filters.size(); i++){
    if(!skim.GetBranch(filters[i].c_str())){
      printf("  filter %s absent, NOT applied\n", filters[i].c_str()); continue; }
    skim.SetBranchStatus(filters[i].c_str(), 1);
    skim.SetBranchAddress(filters[i].c_str(), &fval[i]);
  }

  Jets jc, jp;
  for(auto pr : {std::make_pair(&calo, &jc), std::make_pair(&pf, &jp)}){
    TChain *c = pr.first; Jets *j = pr.second;
    c->SetBranchStatus("*", 0);
    for(const char *b : {"nref", "rawpt", "jtpt", "jteta", "jtphi"}) c->SetBranchStatus(b, 1);
    c->SetBranchAddress("nref", &j->n);
    c->SetBranchAddress("rawpt", j->raw);
    c->SetBranchAddress("jtpt", j->pt);
    c->SetBranchAddress("jteta", j->eta);
    c->SetBranchAddress("jtphi", j->phi);
  }

  // ---- JEC, as in PbPb_scan.C / pp_scan.C (paths relative to this directory)
  const std::string jecDir = "../../../JetEnergyCorrections/";
  std::vector<std::string> fCalo, fPF;
  if(isPP){
    fCalo = {jecDir + "Spring18_ppRef5TeV_V6_DATA_L2Relative_AK4Calo.txt",
             jecDir + "Spring18_ppRef5TeV_V6_DATA_L2L3Residual_AK4Calo.txt"};
    fPF   = {jecDir + "Spring18_ppRef5TeV_V6_DATA_L2Relative_AK4PF.txt",
             jecDir + "Spring18_ppRef5TeV_V6_DATA_L2L3Residual_AK4PF.txt"};
  }
  else{
    fCalo = {jecDir + "Autumn18_HI_V8_DATA_L2Relative_AK4Calo.txt",
             jecDir + "Autumn18_HI_V8_DATA_L2L3Residual_AK4Calo.txt"};
    fPF   = {jecDir + "Autumn18_HI_V8_DATA_L2Relative_AK4PF.txt",
             jecDir + "Autumn18_HI_V8_DATA_L2L3Residual_AK4PF.txt"};
  }
  for(auto &f : fCalo) if(gSystem->AccessPathName(f.c_str())){ printf("ERROR: JEC file %s missing\n", f.c_str()); return; }
  for(auto &f : fPF)   if(gSystem->AccessPathName(f.c_str())){ printf("ERROR: JEC file %s missing\n", f.c_str()); return; }
  JetCorrector jecCalo(fCalo), jecPF(fPF);
  auto corr = [](JetCorrector &J, double raw, double eta, double phi){
    J.SetJetPT(raw); J.SetJetEta(eta); J.SetJetPhi(phi); return J.GetCorrectedPT(); };

  // ---- output ---------------------------------------------------------------
  TFile *wf = TFile::Open(output, "RECREATE");
  if(!wf || wf->IsZombie()){ printf("cannot open %s\n", output.Data()); return; }

  const int nPt = ptEdges.size() - 1;
  const int nUse = isPP ? 1 : nCls;
  TH2D *hCP[nCls], *hPC[nCls], *hRaw[nCls], *hFor[nCls];
  TH1D *hPF[nCls], *hPFm[nCls], *hCa[nCls], *hCam[nCls], *hDR[nCls];
  for(int c = 0; c < nUse; c++){
    const char *cn = isPP ? "pp" : clsName[c];
    hCP[c]  = new TH2D(Form("h_caloOverPF_vsPFPt_C%d", c),
                       Form("%s; PF jet p_{T} [GeV]; p_{T}^{calo} / p_{T}^{PF}", cn),
                       nPt, ptEdges.data(), nRatio, ratioLo, ratioHi);
    hPC[c]  = new TH2D(Form("h_PFOverCalo_vsCaloPt_C%d", c),
                       Form("%s; calo jet p_{T} [GeV]; p_{T}^{PF} / p_{T}^{calo}", cn),
                       nPt, ptEdges.data(), nRatio, ratioLo, ratioHi);
    hRaw[c] = new TH2D(Form("h_rawCaloOverRawPF_vsPFPt_C%d", c),
                       Form("%s, no JEC; PF jet p_{T} [GeV]; raw p_{T}^{calo} / raw p_{T}^{PF}", cn),
                       nPt, ptEdges.data(), nRatio, ratioLo, ratioHi);
    hFor[c] = new TH2D(Form("h_forestCaloOverPF_vsPFPt_C%d", c),
                       Form("%s, forest jtpt; forest PF jet p_{T} [GeV]; jtpt^{calo} / jtpt^{PF}", cn),
                       nPt, ptEdges.data(), nRatio, ratioLo, ratioHi);
    hPF[c]  = new TH1D(Form("h_pfPt_C%d", c),          Form("%s, PF jets; p_{T} [GeV]", cn), nPt, ptEdges.data());
    hPFm[c] = new TH1D(Form("h_pfPt_matched_C%d", c),  Form("%s, PF jets with a calo match; p_{T} [GeV]", cn), nPt, ptEdges.data());
    hCa[c]  = new TH1D(Form("h_caloPt_C%d", c),        Form("%s, calo jets; p_{T} [GeV]", cn), nPt, ptEdges.data());
    hCam[c] = new TH1D(Form("h_caloPt_matched_C%d", c),Form("%s, calo jets with a PF match; p_{T} [GeV]", cn), nPt, ptEdges.data());
    hDR[c]  = new TH1D(Form("h_matchDR_C%d", c),       Form("%s, matched pairs; #DeltaR", cn), 40, 0., matchDR);
  }

  const char *cutName[] = {"all", "filters", "|vz|<15", "centrality", "trigger OR / all", "PF bit", "calo bit"};
  const int nCut = 7;
  TH1D *hCount = new TH1D("h_evtCount", "event cut flow", nCut, 0, nCut);
  for(int i = 0; i < nCut; i++) hCount->GetXaxis()->SetBinLabel(i + 1, cutName[i]);
  TH1D *hHiBin = new TH1D("h_hiBin", "hiBin, selected events", 200, 0, 200);

  // ---- loop -----------------------------------------------------------------
  long nEvt = evt.GetEntries();
  if(maxEvents > 0 && maxEvents < nEvt) nEvt = maxEvents;
  long nPairs = 0;

  struct JetRec { double pt, raw, forest, eta, phi; };
  std::vector<JetRec> vc, vp;

  for(long ev = 0; ev < nEvt; ev++){
    if(ev % 100000 == 0) printf("  event %ld / %ld\n", ev, nEvt);
    evt.GetEntry(ev); skim.GetEntry(ev); hlt.GetEntry(ev);
    hCount->Fill(0.);

    bool pass = true;
    for(size_t i = 0; i < filters.size(); i++) if(fval[i] != 1) pass = false;
    if(!pass) continue;
    hCount->Fill(1.);
    if(std::fabs(vz) > vzMax) continue;
    hCount->Fill(2.);

    std::vector<int> fillCls;
    if(isPP) fillCls.push_back(0);
    else{
      if(hiBin < 0 || hiBin >= 180) continue;
      for(int c = 0; c < nCls; c++) if(hiBin >= clsLo[c] && hiBin < clsHi[c]) fillCls.push_back(c);
    }
    hCount->Fill(3.);

    const bool fPF = havePF && bitPF == 1, fCa = haveCalo && bitCalo == 1;
    if(requireTrigger && !fPF && !fCa) continue;
    hCount->Fill(4.);
    if(fPF) hCount->Fill(5.);
    if(fCa) hCount->Fill(6.);
    if(!isPP) hHiBin->Fill(hiBin);

    calo.GetEntry(ev); pf.GetEntry(ev);
    if(jc.n > jetMax || jp.n > jetMax){ printf("ERROR: nref %d / %d > jetMax %d\n", jc.n, jp.n, jetMax); break; }

    auto build = [&](Jets &src, JetCorrector &jec, std::vector<JetRec> &out){
      out.clear();
      for(int i = 0; i < src.n; i++){
        if(std::fabs(src.eta[i]) > etaMax) continue;
        double x = corr(jec, src.raw[i], src.eta[i], src.phi[i]);
        if(x < ptMinMatch) continue;
        out.push_back({x, src.raw[i], src.pt[i], src.eta[i], src.phi[i]});
      }
    };
    build(jc, jecCalo, vc);
    build(jp, jecPF, vp);

    // greedy one-to-one: repeatedly take the closest remaining pair
    std::vector<int> mc(vc.size(), -1), mp(vp.size(), -1);
    while(true){
      double best = matchDR; int bi = -1, bj = -1;
      for(size_t i = 0; i < vc.size(); i++){
        if(mc[i] >= 0) continue;
        for(size_t j = 0; j < vp.size(); j++){
          if(mp[j] >= 0) continue;
          double d = dR(vc[i].eta, vc[i].phi, vp[j].eta, vp[j].phi);
          if(d < best){ best = d; bi = i; bj = j; }
        }
      }
      if(bi < 0) break;
      mc[bi] = bj; mp[bj] = bi;
    }

    for(int c : fillCls){
      for(size_t j = 0; j < vp.size(); j++){
        hPF[c]->Fill(vp[j].pt);
        if(mp[j] >= 0) hPFm[c]->Fill(vp[j].pt);
      }
      for(size_t i = 0; i < vc.size(); i++){
        hCa[c]->Fill(vc[i].pt);
        if(mc[i] < 0) continue;
        hCam[c]->Fill(vc[i].pt);
        const JetRec &a = vc[i], &b = vp[mc[i]];
        hCP[c]->Fill(b.pt, a.pt/b.pt);
        hPC[c]->Fill(a.pt, b.pt/a.pt);
        if(b.raw > 0)    hRaw[c]->Fill(b.pt, a.raw/b.raw);
        if(b.forest > 0) hFor[c]->Fill(b.forest, a.forest/b.forest);
        hDR[c]->Fill(dR(a.eta, a.phi, b.eta, b.phi));
      }
    }
    for(size_t i = 0; i < vc.size(); i++) if(mc[i] >= 0) nPairs++;
  }

  // ---- summary --------------------------------------------------------------
  printf("\n  cut flow:\n");
  for(int i = 1; i <= nCut; i++)
    printf("    %-12s %12.0f\n", cutName[i - 1], hCount->GetBinContent(i));
  printf("  matched pairs: %ld\n", nPairs);
  printf("  mean calo/PF (manual JEC), PF pT 200-320 GeV:\n");
  for(int c = 0; c < nUse; c++){
    TH1D *p = hCP[c]->ProjectionY(Form("tmp_%d", c), hCP[c]->GetXaxis()->FindBin(200. + 1e-6),
                                  hCP[c]->GetXaxis()->FindBin(320. - 1e-6));
    printf("    %-7s %8.0f pairs  mean %.4f +- %.4f\n", isPP ? "pp" : clsName[c],
           p->GetEntries(), p->GetMean(), p->GetMeanError());
    delete p;
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
  TString prov = Form("scan: caloVsPFJetScale_scan.C\ngitHash: %s\ndirty: %s\n"
                      "input: %s\nsystem: %s\ncalo: %s (%s)\npf: %s (%s)\n"
                      "trigger: %s %s || %s %s, %s\netaMax: %g\nptMinMatch: %g\nmatchDR: %g\nvzMax: %g\nevents: %ld\n",
                      hash.c_str(),
                      hash == "unavailable" ? "unknown" : (dirty == "unavailable" ? "no" : "yes"),
                      input.Data(), isPP ? "pp" : "PbPb",
                      caloTree, fCalo[0].c_str(), pfTree, fPF[0].c_str(),
                      trigPF, havePF ? "" : "(absent)", trigCalo, haveCalo ? "" : "(absent)", requireTrigger ? "required" : "NOT required",
                      etaMax, ptMinMatch, matchDR, vzMax, nEvt);

  wf->cd();
  TNamed("provenance", prov.Data()).Write();
  hCount->Write();
  if(!isPP) hHiBin->Write();
  for(int c = 0; c < nUse; c++){
    hCP[c]->Write(); hPC[c]->Write(); hRaw[c]->Write(); hFor[c]->Write();
    hPF[c]->Write(); hPFm[c]->Write(); hCa[c]->Write(); hCam[c]->Write(); hDR[c]->Write();
  }
  wf->Close();
  printf("  wrote %s\n", output.Data());
}
