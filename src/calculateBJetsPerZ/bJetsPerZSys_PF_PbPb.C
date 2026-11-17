// Systematic uncertainties of the PF b jets per Z and of b / inclusive: PbPb
// with the calo-match jet ID (caloMatchedPF), pp with plain ak4PF (PFJets).
// A copy of bJetsPerZSys_caloJets_PbPb.C (2026-10-07) for the chain
//   bJetChain_caloJets_PbPb.C("dedicated","caloMatchedPF","PFJets")
// (run it first); sources and their correlations are as in the calo macro --
// read its header -- except:
//   inputs      PbPb response = the 2026-10-06 calo-matched anaMuTag PYTHIA+HYDJET
//               scan, h_matchedRecoJetPt_genJetPt_bJets_C<slice>_muTagged; purity
//               systematic bPuritySys_caloMatchedPF_PbPb.root; pp chain
//               bJetSpectrum_inclusive_PFJets_pp.root; JEU unfold_PYTHIA_PFJets /
//               unfold_PYTHIAHYDJET_caloMatchedPF; inclusive side
//               inclJetsPerZSys_caloMatchedPF_ppPFJets_JEUppCorr.root.
//   unfolding   the chain's: ID-FREE response (reco rows of the calo-matched
//               matrix x 1/eps_MC, measureCaloMatchEff.C), floor 70, N = 2 in
//               pp and every PbPb class; inclusive N 3/9/6/4/3.
//   caloMatchID      NEW. The ID-free response's row scaling at +/-50%
//               (f = 0.5, 1.5) on the PbPb side; corr within PbPb, in quadrature
//               with pp, whose own +/-50% comes from the pp chain's turnOn slot
//               (unfoldBJetSpectrum_caloJets_pp.C carries it there for PFJets;
//               the calo macro's separate turnOn source is dropped).
//   caloMatchIDdata  NEW, PbPb only. The PbPb data divided by eps_data, the
//               calo-matched / plain-PF ratio of the muon-tagged SingleMuon
//               spectra (normalised over 150-300 GeV), unfolded, full
//               difference: the most genuine-jet loss the ID can cause in data
//               (it also holds fake rejection). corr within PbPb.
//   muonTagFreq the chain divides every system by the SAME f_pp (x a muon-
//               efficiency ratio), so f_pp's uncertainty (pp chain file) is
//               applied to every system as one correlated shift: it stays on
//               b per Z and on the b-jet fraction and cancels in PbPb/pp, R_CP
//               and b/incl PbPb/pp.
// MODE "PFJets" (argument; 2026-10-07, user: everything plain PF, only the PbPb
// MinBias leg of the inclusive calo-matched): chain bJetChain_caloJets_PbPb.C(
// "dedicated","PFJets",""); PbPb response = sibling-repo 2026-02-18
// muTaggedJetsNoTrigger (nominal classes, h_matchedRecoJetPt_genJetPt_bJets_C<n>),
// plain (no ID scaling, no caloMatchID sources); floor 70 with the truth axis
// extended to 20 GeV and the response's weight outliers clipped (2026-10-08,
// bJetUnfoldExt.h / clipResponseOutliers.h; floor 50, truth from 50 before);
// N pp 3, PbPb 4/4/3/3;
// inclusive inclJetsPerZSys_PF_MBcaloMatched (N 3/9/6/4/3); JEU
// unfold_PYTHIAHYDJET_PFJets_JEUShift.root; purity bPuritySys_PFJets_PbPb.root.
// The PbPb data carry no mu12 trigger-efficiency correction (user's choice; not a
// source here). Outputs bJetsPerZSys_PF_MBcaloMatched<sfx>.root.
// Data-driven prior (2026-10-07, dataDrivenPrior.h): the PbPb response gen rows
// are weighted by (g/100)^-k with the chain's k (priorTilt_<cls>), so the nominal
// unfolding is the chain's. Sources added in this mode:
//   prior       MC prior (k = 0) instead, full difference; uncorr (the k differ
//               per class); pp: h_unf_sysRel_prior of unfoldBJetSpectrum_caloJets_pp.C
//   recoFloor   reco floor 80 instead of 70 GeV (was 100 vs 50); corr within
//               PbPb, quadrature with pp (h_unf_sysRel_recoFloor)
//   iterations  PbPb: envelope of N - 1 and N + 1 (was N + 1 only)
// NOT INCLUDED, beyond the calo list: the PbPb purity's light-template MC
// statistics (lightSource = PYTHIA+HYDJET for caloMatchedPF), and the mixed scan
// generations of the pp side (jetCollectionInputs.h, PFJets).
//
// Output  rootFiles/systematics/bJetsPerZSys_caloMatchedPF_ppPFJets<sfx>.root, same names
//         figures/systematics/bJetsPerZSys_caloMatchedPF_ppPFJets_*<sfx>.pdf
//         figures/bJetsPerZ_caloMatchedPF_dedicatedPurity_ppPFJets/{bJetsPerZ,bOverIncl,bJetFraction}_PF_*<sfx>.pdf
//
// Usage, from src/calculateBJetsPerZ/:
//   root -l -b -q -e 'gSystem->Load("/home/clayton/Programs/RooUnfold/build/libRooUnfold.so"); gInterpreter->AddIncludePath("/home/clayton/Programs/RooUnfold/build");' \
//        'bJetsPerZSys_PF_PbPb.C(true)'      // JEU pp-PbPb correlated
//        'bJetsPerZSys_PF_PbPb.C(true,"PFJets",true)'   // JEU not flattened, _noFlatJEU file
//        'bJetsPerZSys_PF_PbPb.C(true,"PFJets",true,"_noFlatJEU_fakeGeoOnly","_fakeGeoOnly")'   // + that inclusive file
//        event-count trigger stitch (bJetChain with outTag "_evtCountStitch" first):
//        'bJetsPerZSys_PF_PbPb.C(true,"PFJets",true,"_noFlatJEU_fakeGeoOnly_evtCountStitch","_fakeGeoOnly_evtCountStitch","_evtCountStitch")'

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TBox.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TSystem.h"
#include "TColor.h"
#include <vector>
#include <cmath>
#include "../../headers/plotting/coarseCent.h"
#include "../../headers/plotting/plotStyle.h"
#include "../../headers/plotting/flattenSystematic.h"
#include "../../headers/functions/dataDrivenPrior.h"
#include "../../headers/functions/jeuShift.h"
#include "TVectorD.h"
#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#endif
#include "../../headers/functions/bJetUnfoldExt.h"
#include "../../headers/functions/purityVariant.h"   // purVar(): purity-variant file suffix
#include "../../headers/functions/clipResponseOutliers.h"
#include "../../headers/functions/jerSmear.h"   // JER scheme (PFJets): nominal x1.10, variation x1.10
#include "../../headers/publishedData/CMS_ZYield_5p02TeV.h"   // BJET_NZ_SOURCE=published: N_Z from the published yields

namespace bjspf {

TString chainPath = "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/BJetPerZ_caloMatchedPF/bJetChain_caloMatchedPF_PbPb_dedicatedPurity_ppPFJets.root";
const char *respPath  = "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PYTHIAHYDJET/PYTHIAHYDJET_response_DiJet_manualJEC_caloJetMatched_anaMuTag_pThat-15_mu12_pTmu-15_tight_vzReweight_hiBinReweight_hiBinShift-10_jetTrkMaxFilter_doPThatCorrelationFilterTight_2026-10-6_ultraFineCentBins.root";
const char *respFmt   = "h_matchedRecoJetPt_genJetPt_bJets_C%d_muTagged";
bool modePF = false;   // "PFJets" mode (header)
const char *purSysPath= "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/systematics/bPuritySys_caloMatchedPF_PbPb.root";
const char *ppBPath   = "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/CorrectedBJetSpectra/Data/bJetSpectrum_inclusive_PFJets_pp.root";
const char *jeuPbPbPath = "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/JEU/unfold_PYTHIAHYDJET_caloMatchedPF_JEUShift.root";
const char *jeuPPPath   = "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/JEU/unfold_PYTHIA_PFJets_JEUShift.root";
const char *inclSysPath = "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/systematics/inclJetsPerZSys_caloMatchedPF_ppPFJets_JEUppCorr.root";
const char *idEffPath   = "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/correctionFactors/caloMatchEff.root";
const double idEpsMin   = 0.1;   // as bJetChain_caloJets_PbPb.C
// data: calo-matched and plain-PF PbPb SingleMuon (measureCaloMatchEff.C's eps_data)
const char *dataIDM = "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/PbPb_SingleMuon_manualJEC_caloJetMatched_mu12_pTmu-15to999_tight_mu12TriggerEfficiencyCorrection_jetTrkMaxFilter_WDecayFilter_2026-10-6_ultraFineCentBins.root";
// plain PF: the chain's PbPb data (jetCollectionInputs.h dataPathPbPbPFJets, set by setJetCollection): the unweighted
// 2026-9-24 scan since 2026-11-10 (user), the weighted 2026-10-26 one with BJET_PBPB_MU12WEIGHT=1 (2026-11-04 to 11-10)
const char *dataIDP = [](){ const char *e = std::getenv("BJET_PBPB_MU12WEIGHT"); return (e && TString(e) == "1")
  ? "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/PbPb_SingleMuon_manualJEC_mu12_pTmu-15to999_tight_mu12TriggerEfficiencyCorrection_jetTrkMaxFilter_WDecayFilter_2026-10-26_ultraFineCentBins.root"
  : "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/PbPb_SingleMuon_manualJEC_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_2026-9-24_ultraFineCentBins.root"; }();
const char *collText = "PF jets";   // PbPb calo-matched; short for the panel titles
const char *outDir    = "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/systematics/";
const char *figSysDir = "/home/clayton/Analysis/code/bJetRaaAnalysis/figures/systematics/";
TString figDir    = "/home/clayton/Analysis/code/bJetRaaAnalysis/figures/bJetsPerZ_caloMatchedPF_dedicatedPurity_ppPFJets/";

#ifdef BJET_PT_7BINS   // 7-bin variant (caloJetBPurityFitCore.h): 80-200 GeV in the scans' native edges
const int NEdge = 13; const double edges[NEdge] = {50, 60, 70, 80, 90, 100, 110, 120, 135, 150, 200, 300, 500};
#else
const int NEdge = 10; const double edges[NEdge] = {50,60,70,80,100,120,150,200,300,500};
#endif
const int NB = NEdge - 1;
const double ptReportMin = 80., ptDrawMax = 200.;
double recoFloor = 70.;   // the chain's, with the ID-free response; PFJets mode 50
const int NS = 5;
const char *sysTag[NS]   = {"pp", "C1", "C2", "C3", "C4"};
const char *sysLabel[NS] = {"pp", "PbPb 0-10%", "PbPb 10-30%", "PbPb 30-50%", "PbPb 50-80%"};
int          nIterB[NS]    = {2, 2, 2, 2, 2};           // b chains: pp (unfoldBJetSpectrum_caloJets_pp.C), PbPb (bJetChain, caloMatchedPF)
int          nIterIncl[NS] = {3, 9, 6, 4, 3};           // inclusiveJetsPerZSys_PF.C (tilt:0.3 closure optimum)
const double muEff[NS]     = {0.9708, 0.8627, 0.9069, 0.9856, 0.9778};
const double muEffErr[NS]  = {0.0077, 0.0350, 0.0347, 0.0151, 0.0105};   // delta eps: |SF_ID - 1| (+) sigma_SF (+) T&P bkg (+) sim. N_eff stat., studyMuonEffSF_binned.C (2026-10-21)
// JER: jerSmear.h (2026-11-05: nominal x1.10, variation x1.10 on it; before: unsmeared, +20%)
const int    nToys = 200, toySeed = 20261005;
// 2026-11-05 (user): the purity sources (bPurity, bPurity_bGS, bPurity_cMult) are exempt from the flattening,
// as the JEU is with noFlatJEU: in 0-10% the trend test missed a step (5% to 150 GeV, 22% at 150-200, 98%
// above) and the mean, 20.5%, exceeded every reported bin but one.
const bool   noFlattenPurity = true;
// 2026-11-05 (user): the iterations source too -- in 0-10% it is 1.4-2.9% up to 150 GeV and 16% at 150-200 GeV
// on R_Z(b), and the mean (over every bin from 80 GeV, 200-500 included) was 19.3%.
const bool   noFlattenIterations = true;
// 2026-11-17 (user): the JER too -- kept bin by bin on every quantity; the flattened mean took in the unreported
// 200-300 GeV overflow (and 300-500 GeV where populated). BJET_SYS_FLATJER=1 = the flattening as before.
const bool   noFlattenJER = !(std::getenv("BJET_SYS_FLATJER") && TString(std::getenv("BJET_SYS_FLATJER")) == "1");
// 2026-11-05 (user): the b-tag working-point source is assigned flat in pT, the mean of its bin-by-bin values over
// the reported bins (not the trend-tested mean over every bin from 80 GeV).
const double bWPMeanLo = 80., bWPMeanHi = 200.;
std::vector<TH1D*> iterSave;   // b-jet N_opt, N_opt +- 1 unfolded spectra per PbPb class, and the JEU inputs / signed shifts (extraction figures)
const bool   flattenSources = !(std::getenv("BJET_SYS_NOFLATTEN") && TString(std::getenv("BJET_SYS_NOFLATTEN")) == "1");   // sources without a pT trend -> their mean (flattenSystematic.h); BJET_SYS_NOFLATTEN=1 (2026-10-28, study) keeps every source bin by bin
bool noFlattenJEU = false;            // last argument: JEU kept bin by bin (2026-10-20, user; as inclusiveJetsPerZSys_PF.C)

bool reported(int b){ return edges[b-1] >= ptReportMin - 1e-6; }   // b = 1..NB

// ---- pieces of the chain, copied so the nominal is reproduced exactly ---------
TH2D* sliceSum2(TFile *f, const char *base, int ci, const char *tag){
  TH2D *s = nullptr;
  for(int si = coarseSliceLo[ci]; si <= coarseSliceHi[ci]; si++){
    TH2D *h = nullptr; f->GetObject(Form("%s_C%d", base, si), h);
    if(!h) continue;
    if(!s){ s = (TH2D*) h->Clone(Form("%s_%s_%d", base, tag, ci)); s->SetDirectory(nullptr); } else s->Add(h);
  }
  return s;
}
TH2D* sliceSumFmt2(TFile *f, const char *fmt, int ci, const char *tag){
  TH2D *s = nullptr;
  if(!TString(f->GetName()).Contains("ultraFineCentBins")){   // nominal classes C1..C4
    f->GetObject(Form(fmt, ci + 1), s);
    if(s){ s = (TH2D*) s->Clone(Form("fmt_%s_%d", tag, ci)); s->SetDirectory(nullptr); }
    return s;
  }
  for(int si = coarseSliceLo[ci]; si <= coarseSliceHi[ci]; si++){
    TH2D *h = nullptr; f->GetObject(Form(fmt, si), h);
    if(!h) continue;
    if(!s){ s = (TH2D*) h->Clone(Form("fmt_%s_%d", tag, ci)); s->SetDirectory(nullptr); } else s->Add(h);
  }
  return s;
}
// bJetChain_caloJets_PbPb.C's idScaleRows: reco rows x (1 + frac (1/eps - 1))
TH2D* idScaleRows(TH2D *src, TH1D *eff, double frac, const char *name){
  TH2D *o = (TH2D*) src->Clone(name); o->SetDirectory(nullptr);
  for(int ix = 1; ix <= o->GetNbinsX(); ix++){
    double e = eff->GetBinContent(eff->FindBin(o->GetXaxis()->GetBinCenter(ix)));
    if(e < idEpsMin) continue;
    double w = 1. + frac*(1./e - 1.);
    for(int iy = 0; iy <= o->GetNbinsY()+1; iy++){ o->SetBinContent(ix, iy, o->GetBinContent(ix, iy)*w); o->SetBinError(ix, iy, o->GetBinError(ix, iy)*w); }
  }
  return o;
}
// PFJets mode (2026-10-08, as bJetChain_caloJets_PbPb.C): both axes extended to
// 20 GeV (bJetUnfoldExt.h), so every variation below unfolds like the chain
TH2D* rebin2D(TH2D *h, const char *name){
  if(modePF) return rebin2DTo(h, bExtEdges(NB + 1, edges), name);
  TH2D *r = new TH2D(name, "", NB, edges, NB, edges); r->SetDirectory(nullptr); r->Sumw2();
  for(int ix = 1; ix <= h->GetNbinsX(); ix++) for(int iy = 1; iy <= h->GetNbinsY(); iy++){
    double c = h->GetBinContent(ix, iy); if(c == 0.) continue;
    int bx = r->GetXaxis()->FindBin(h->GetXaxis()->GetBinCenter(ix)), by = r->GetYaxis()->FindBin(h->GetYaxis()->GetBinCenter(iy));
    if(bx < 1 || bx > NB || by < 1 || by > NB) continue;
    r->SetBinContent(bx, by, r->GetBinContent(bx, by) + c);
    r->SetBinError(bx, by, std::hypot(r->GetBinError(bx, by), h->GetBinError(ix, iy)));
  }
  return r;
}
// the chain's unfolding step: floor on the reco side of R and of the input
TH1D* unfoldB(TH2D *R, TH1D *meas, int N, const char *name){
  if(modePF) return unfoldOnExt<RooUnfoldResponse, RooUnfoldBayes>(R, meas, NB + 1, edges, recoFloor, N, name);
  TH1D *truth = R->ProjectionY(Form("%s_t", name), 1, NB); truth->SetDirectory(nullptr);
  TH2D *Rt = (TH2D*) R->Clone(Form("%s_Rt", name)); Rt->SetDirectory(nullptr);
  for(int ix = 1; ix <= NB; ix++) if(edges[ix] <= recoFloor + 1e-6) for(int iy = 1; iy <= NB; iy++){ Rt->SetBinContent(ix, iy, 0.); Rt->SetBinError(ix, iy, 0.); }
  TH1D *measMC = Rt->ProjectionX(Form("%s_m", name), 1, NB); measMC->SetDirectory(nullptr);
  TH1D *in = (TH1D*) meas->Clone(Form("%s_in", name)); in->SetDirectory(nullptr);
  for(int b = 1; b <= NB; b++) if(edges[b] <= recoFloor + 1e-6){ in->SetBinContent(b, 0.); in->SetBinError(b, 0.); }
  RooUnfoldResponse resp(measMC, truth, Rt, Form("%s_resp", name), "");
  RooUnfoldBayes unf(&resp, in, N);
  unf.SetVerbose(0);
  TH1D *U = (TH1D*) unf.Hunfold()->Clone(name); U->SetDirectory(nullptr);
  delete truth; delete Rt; delete measMC; delete in;
  return U;
}
// extra Gaussian smearing of the reco axis of the FINE matrix (inclusive macro)
TH2D* smearResp(TH2D *r, const char *name){ return modePF ? jerSmear::variation(r, name) : jerSmear::smear(r, name, 0.663); }   // calo: the old +20%
// sum a fine histogram into the b-jet edges
TH1D* rebinEdges(TH1D *h, const char *name){
  TH1D *o = new TH1D(name, "", NB, edges); o->SetDirectory(nullptr); o->Sumw2();
  for(int i = 1; i <= h->GetNbinsX(); i++){
    int j = o->FindBin(h->GetBinCenter(i)); if(j < 1 || j > NB) continue;
    o->SetBinContent(j, o->GetBinContent(j) + h->GetBinContent(i));
    o->SetBinError(j, sqrt(pow(o->GetBinError(j), 2) + pow(h->GetBinError(i), 2)));
  }
  return o;
}
// signed relative shift var/nom - 1 per bin (0 where the nominal is empty)
std::vector<double> shift(TH1D *var, TH1D *nom){
  std::vector<double> d(NB, 0.);
  for(int i = 1; i <= NB; i++){ double n = nom->GetBinContent(i); d[i-1] = n > 0. ? var->GetBinContent(i)/n - 1. : 0.; }
  return d;
}
std::vector<double> relFromFile(TFile *f, const char *name, bool required = true){
  std::vector<double> u(NB, 0.);
  TH1D *h = nullptr; f->GetObject(name, h);
  if(!h){ if(required) printf("WARNING: %s missing in %s\n", name, f->GetName()); return u; }
  if(h->GetNbinsX() != NB){ printf("ERROR: %s has %d bins, expected %d\n", name, h->GetNbinsX(), NB); return u; }
  for(int i = 1; i <= NB; i++) u[i-1] = fabs(h->GetBinContent(i));
  return u;
}

// ---- sources ---------------------------------------------------------------------
typedef std::vector<std::vector<double>> VV;          // [sys][bin]
typedef std::vector<VV> VVV;                          // [var][sys][bin]
VV zeroVV(){ return VV(NS, std::vector<double>(NB, 0.)); }
struct Source {
  TString name, label;
  bool corr = true, ppUncorr = false, ratioOnly = false, nzOnly = false;
  VVV d;                 // b side, signed per variation (corr)
  VV  u = zeroVV();      // b side, unsigned (uncorr)
  VVV di;                // inclusive side, signed per variation (JEU)
  VV  ui = zeroVV();     // inclusive side, unsigned
  VV  ratioVal;          // [q 5..11][bin], ratioOnly (massWindow)
};

// quantities: 0..4 b per Z spectra; 5..8 b PbPb/pp C1..C4; 9..11 b RCP C1..C3;
// 12..15 b/incl PbPb/pp C1..C4; 16..18 b/incl RCP C1..C3;
// 19..23 the b-jet FRACTION b/incl per system (pp, C1..C4)
const int NQ = 24;
void qSys(int q, int &a, int &b, bool &bi){
  if(q >= 19){ bi = true; a = q - 19; b = -1; return; }
  bi = q >= 12;
  int r = bi ? q - 12 + 5 : q;
  if(r < 5){ a = r; b = -1; } else if(r < 9){ a = r - 4; b = 0; } else { a = r - 8; b = 4; }
}
TString qTag(int q){
  int a, b; bool bi; qSys(q, a, b, bi);
  TString p = bi ? "bOverIncl" : "b";
  if(b < 0) return Form("%s_%s", p.Data(), sysTag[a]);
  return Form("%s_%s_C%d", p.Data(), b == 0 ? "PbPbOverPP" : "RCP", a);
}
TString qLabel(int q){
  int a, b; bool bi; qSys(q, a, b, bi);
  if(b < 0) return bi ? Form("#it{b}-jet fraction, %s", sysLabel[a]) : TString(sysLabel[a]);
  return Form("%s%s / %s", bi ? "#it{b} / inclusive, " : "", sysLabel[a] + 5, b == 0 ? "pp" : "50-80%");
}
double relOf(const Source &S, int q, int i){
  int a, b; bool bi; qSys(q, a, b, bi);
  if(S.ratioOnly){ if(bi || b < 0) return 0.; return S.ratioVal[q-5][i]; }
  if(bi && S.nzOnly) return 0.;
  if(!S.corr){
    double ua = S.u[a][i], ub = b < 0 ? 0. : S.u[b][i];
    double ia = bi ? S.ui[a][i] : 0., ib = (bi && b >= 0) ? S.ui[b][i] : 0.;
    return sqrt(ua*ua + ub*ub + ia*ia + ib*ib);
  }
  double m = 0., ma = 0., mb = 0.;
  for(size_t v = 0; v < S.d.size(); v++){
    double ra = 1. + S.d[v][a][i], rb = b < 0 ? 1. : 1. + S.d[v][b][i];
    if(bi && v < S.di.size()){ ra /= 1. + S.di[v][a][i]; if(b >= 0) rb /= 1. + S.di[v][b][i]; }
    m  = TMath::Max(m, fabs(ra/rb - 1.));
    ma = TMath::Max(ma, fabs(ra - 1.)); if(b >= 0) mb = TMath::Max(mb, fabs(rb - 1.));
  }
  double r = (S.ppUncorr && b == 0) ? sqrt(ma*ma + mb*mb) : m;
  if(bi){ double ia = S.ui[a][i], ib = b >= 0 ? S.ui[b][i] : 0.; r = sqrt(r*r + ia*ia + ib*ib); }
  return r;
}

} // namespace bjs
using namespace bjspf;

// chainTag: the chain variant to read, e.g. "_dedicatedPurity" (bJetChain_caloJets_PbPb.C
// purityMode "dedicated"); every output name carries it in front of the JEU suffix
void bJetsPerZSys_PF_PbPb(bool jeuPPCorr = true, TString mode = "",
                          bool noFlatJEU = false,   // JEU exempt from flattenSources; _noFlatJEU in every output name
                          TString inclTag = "",     // PFJets mode: suffix of the inclusive file read (e.g. "_noFlatJEU_fakeGeoOnly")
                          TString outTag = "",      // appended to every output name after _noFlatJEU (e.g. "_fakeGeoOnly")
                          TString chainTag = "",    // PFJets mode: suffix of the chain file read (bJetChain's outTag, e.g. "_evtCountStitch")
                          TString altWPChainTag = "")   // PFJets mode: chain file of the alternative b-tag working point ("_csvV2medium"); adds source bTagWP
{
  noFlattenJEU = noFlatJEU;
  initPlotStyle();
  modePF = mode == "PFJets";
  // N_Z from the PUBLISHED yields, nominal since 2026-11-09 (user), PFJets mode: the N_Z
  // sources are pubZYield and pubZAcc (block "N_Z sources"). BJET_NZ_SOURCE=own = our own
  // count: reads and writes _ownNZ files (chain and inclusive run with the same environment).
  // The legacy calo-matched mode always uses its own count.
  if(!ptBinGuard(outTag)) return;   // 7-bin jet-pT axis variant (caloJetBPurityFitCore.h)
  const bool pubNZ = modePF && pubZ::enabled();
  if(modePF){
    chainTag += pubZ::tag(); inclTag += pubZ::tag(); outTag += pubZ::tag();
    if(altWPChainTag != "") altWPChainTag += pubZ::tag();
    printf("N_Z: %s; chain %s, inclusive %s, output %s\n", pubNZ ? "published yields" : "own Z->mumu count", chainTag.Data(), inclTag.Data(), outTag.Data());
  }
  TString outName = "bJetsPerZSys_caloMatchedPF_ppPFJets";
  if(modePF){
    static TString chainPathS; chainPathS = "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/BJetPerZ_PFJets/bJetChain_PFJets_PbPb_dedicatedPurity" + chainTag + ".root";
    chainPath   = chainPathS.Data();
    respPath    = "/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/scanningOutput/PYTHIAHYDJET/latest/response/PYTHIAHYDJET_response_DiJet_pThat-15_mu12_pTmu-15_tight_vzReweight_hiBinReweight_hiBinShift-10_leadingXjetDumpFilter_jetTrkMaxFilter_removeHYDJETjet0p35CutOnGen_weightCut_2026-2-18_muTaggedJetsNoTrigger.root";
    respFmt     = "h_matchedRecoJetPt_genJetPt_bJets_C%d";
    purSysPath  = strdup(Form("/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/systematics/bPuritySys_PFJets_PbPb%s.root", purVar()));
    ppBPath     = strdup(Form("/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/CorrectedBJetSpectra/Data/bJetSpectrum_inclusive_PFJets_pp%s.root", ppPurVar()));
    jeuPbPbPath = "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/JEU/unfold_PYTHIAHYDJET_PFJets_JEUShift.root";
    static TString inclSysPathS; inclSysPathS = Form("/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/systematics/inclJetsPerZSys_PF_MBcaloMatched_JEUppCorr%s.root", inclTag.Data());
    inclSysPath = inclSysPathS.Data();   // owned: Form() buffers rotate
    figDir      = "/home/clayton/Analysis/code/bJetRaaAnalysis/figures/bJetsPerZ_PFJets_dedicatedPurity/";
    recoFloor   = bRecoFloorPF();   // 60 since 2026-11-05 (70 from 2026-10-08, with the extended truth axis; 50 before)
    const int nb[NS] = {bNIterPP, bNIterPbPb[0], bNIterPbPb[1], bNIterPbPb[2], bNIterPbPb[3]}, ni[NS] = {2, 7, 5, 3, 2};   // b: bJetUnfoldExt.h; inclusive N as inclusiveJetsPerZSys_PF.C nIterNom (0-10%: 9, was 8 until 2026-10-21)
    for(int s = 0; s < NS; s++){ nIterB[s] = nb[s]; nIterIncl[s] = ni[s]; }
    outName = "bJetsPerZSys_PF_MBcaloMatched";
  } else if(mode != ""){ printf("ERROR: mode is \"\" or PFJets\n"); return; }
  gSystem->mkdir(outDir, true); gSystem->mkdir(figSysDir, true); gSystem->mkdir(figDir, true);
  // _noFlatJEU: JEU is computed here from the JEU files, not read from the
  // inclusive file; the inclusive file read is chosen by inclTag (default nominal)
  const TString sfx = TString(jeuPPCorr ? "_JEUppCorr" : "_JEUppUncorr") + (noFlatJEU ? "_noFlatJEU" : "") + outTag;
  TFile *fC = TFile::Open(chainPath), *fR = TFile::Open(respPath), *fPS = TFile::Open(purSysPath), *fP = TFile::Open(ppBPath),
        *fJA = TFile::Open(jeuPbPbPath), *fJP = TFile::Open(jeuPPPath), *fI = TFile::Open(inclSysPath), *fE = TFile::Open(idEffPath),
        *fDM = TFile::Open(dataIDM), *fDP = TFile::Open(dataIDP);
  for(TFile *f : {fC, fR, fPS, fP, fJA, fJP, fI, fE, fDM, fDP}) if(!f || f->IsZombie()){ printf("ERROR: an input did not open\n"); return; }

  // ---- nominal pieces per PbPb class -----------------------------------------------
  TH1D *meas[NS] = {nullptr}, *unfNom[NS] = {nullptr}, *purity[NS] = {nullptr}, *tagFreq[NS] = {nullptr}, *tagged[NS] = {nullptr};
  TH2D *Rfine[NS] = {nullptr}, *R[NS] = {nullptr}, *RfineRaw[NS] = {nullptr};
  TH1D *epsMC[NS] = {nullptr};
  double maxDev = 0.;
  for(int s = 1; s < NS; s++){
    fC->GetObject(Form("bJetPt_meas_%s", sysTag[s]), meas[s]);
    fC->GetObject(Form("bPurity_%s", sysTag[s]), purity[s]);
    fC->GetObject(Form("muTagFreq_%s", sysTag[s]), tagFreq[s]);
    fC->GetObject(Form("muTagJetPt_%s", sysTag[s]), tagged[s]);
    TH1D *uChain = nullptr; fC->GetObject(Form("bJetPt_unf_%s", sysTag[s]), uChain);
    if(!meas[s] || !purity[s] || !tagFreq[s] || !tagged[s] || !uChain){ printf("ERROR: chain output missing for %s\n", sysTag[s]); return; }
    RfineRaw[s] = sliceSumFmt2(fR, respFmt, s-1, "resp");
    if(modePF && RfineRaw[s]) clipResponseOutliers(RfineRaw[s], 100., false);   // as the chain
    if(modePF && RfineRaw[s]){
      iterSave.push_back(jerSmear::resolution(RfineRaw[s], Form("h_bJER_resSim_%s", sysTag[s])));   // simulation width (JER figure, 2026-11-17)
      TH2D *sm = jerSmear::nominal(RfineRaw[s], Form("RfineJERnom_%s", sysTag[s])); delete RfineRaw[s]; RfineRaw[s] = sm; }   // JER nominal, as the chain
    fE->GetObject(Form("h_eff_muTagB_%s", sysTag[s]), epsMC[s]);
    if(!RfineRaw[s] || !epsMC[s]){ printf("ERROR: response or eps_MC missing for %s\n", sysTag[s]); return; }
    double kTilt = 0.;
    if(modePF){ TVectorD *kv = (TVectorD*) fC->Get(Form("priorTilt_%s", sysTag[s])); if(!kv){ printf("ERROR: priorTilt_%s missing in the chain file\n", sysTag[s]); return; } kTilt = (*kv)[0]; }
    Rfine[s] = modePF ? tiltGenRows(RfineRaw[s], kTilt, Form("Rfine_%s", sysTag[s]))
                      : idScaleRows(RfineRaw[s], epsMC[s], 1., Form("RfineID_%s", sysTag[s]));   // ID-free response, as the chain
    Rfine[s]->SetDirectory(nullptr);
    R[s] = rebin2D(Rfine[s], Form("R_%s", sysTag[s]));
    unfNom[s] = unfoldB(R[s], meas[s], nIterB[s], Form("unfNom_%s", sysTag[s]));
    for(int b = 1; b <= NB; b++) if(reported(b) && uChain->GetBinContent(b) > 0.) maxDev = TMath::Max(maxDev, fabs(unfNom[s]->GetBinContent(b)/uChain->GetBinContent(b) - 1.));
  }
  printf("nominal unfolding rebuilt: largest deviation from the chain's bJetPt_unf %.2e\n", maxDev);
  if(maxDev > 1e-6){ printf("ERROR: the chain's unfolding is not reproduced; inputs changed?\n"); return; }

  std::vector<Source> srcs;
  auto addPPFile = [&](Source &S, const char *hist){ S.u[0] = relFromFile(fP, hist); };

  // b purity. SPLIT (2026-10-28, user), when the purity file carries h_sysRel_uncorr_C<n> (systematics_CSVv2Cut_PF.C):
  // bGS and cMult are the same template assumptions in pp and PbPb -> one CORRELATED source each, every
  // signed variation h_varRel_<src>_<k> applied to all systems at once (PbPb: measured spectrum scaled, unfolded;
  // pp: the reco-level signed shift x the pp chain's own unfolded / reco ratio of the total purity systematic,
  // h_incl_sysRel_bPurity / h_sysRel_total, since the pp unfolding is not redone here); the rest of the purity
  // systematic (fitRangeHigh, tplStat, JEU, JER) stays "bPurity", uncorrelated, scaled the same way in pp.
  TH1D *relUncTest = nullptr; fPS->GetObject("h_sysRel_uncorr_C1", relUncTest);
  const bool purSplit = relUncTest != nullptr;
  std::vector<double> ppGenOverReco(NB, 1.);
  TFile *fPSpp = nullptr;
  if(purSplit){
    fPSpp = TFile::Open(Form("/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/systematics/bPuritySys_total_PFJets_pp%s.root", ppPurVar()));
    if(!fPSpp || fPSpp->IsZombie()){ printf("ERROR: pp purity file for the split purity sources\n"); return; }
    std::vector<double> g = relFromFile(fP, "h_incl_sysRel_bPurity"), r = relFromFile(fPSpp, "h_sysRel_total");
    for(int i = 0; i < NB; i++) ppGenOverReco[i] = r[i] > 0. ? g[i]/r[i] : 1.;
    printf("purity split: pp gen/reco ratio of the purity systematic:"); for(int i = 0; i < NB; i++) printf(" %.2f", ppGenOverReco[i]); printf("\n");
    const char *cs[2] = {"bGS", "cMult"}; const char *cl[2] = {"#it{b} purity: g#rightarrowb#bar{b} share", "#it{b} purity: c fraction"};
    const int nv[2] = {2, 4};
    for(int c = 0; c < 2; c++){
      Source S; S.name = Form("bPurity_%s", cs[c]); S.label = cl[c]; S.corr = true;
      for(int k = 0; k < nv[c]; k++){
        VV d = zeroVV();
        std::vector<double> vp(NB, 0.);
        { TH1D *h = nullptr; fPSpp->GetObject(Form("h_varRel_%s_%d", cs[c], k), h); if(!h){ printf("ERROR: pp h_varRel_%s_%d\n", cs[c], k); return; }
          for(int i = 0; i < NB; i++) d[0][i] = h->GetBinContent(i + 1)*ppGenOverReco[i]; }
        for(int s = 1; s < NS; s++){
          TH1D *h = nullptr; fPS->GetObject(Form("h_varRel_%s_%d_%s", cs[c], k, sysTag[s]), h);
          if(!h){ printf("ERROR: h_varRel_%s_%d_%s\n", cs[c], k, sysTag[s]); return; }
          TH1D *mv = (TH1D*) meas[s]->Clone("mv"); mv->SetDirectory(nullptr);
          for(int b = 1; b <= NB; b++){ const double f = 1. + h->GetBinContent(b); mv->SetBinContent(b, meas[s]->GetBinContent(b)*f); mv->SetBinError(b, meas[s]->GetBinError(b)*f); }
          TH1D *u = unfoldB(R[s], mv, nIterB[s], "uv"); d[s] = shift(u, unfNom[s]);
          delete u; delete mv;
        }
        S.d.push_back(d);
        // saved for the purity extraction figure (plotBJetPuritySys_PF.C, 2026-11-17): the signed unfolded shift per system
        for(int s = 0; s < NS; s++){ TH1D *h = new TH1D(Form("h_bPur_%s_%d_%s", cs[c], k, sysTag[s]), "", NB, edges); h->SetDirectory(nullptr);
          for(int b = 1; b <= NB; b++) h->SetBinContent(b, d[s][b-1]); iterSave.push_back(h); }
      }
      srcs.push_back(S);
    }
  }
  // b purity (uncorr): the total, or with purSplit the uncorrelated remainder
  {
    Source S; S.name = "bPurity"; S.label = purSplit ? "#it{b} purity (uncorr. part)" : "#it{b} purity"; S.corr = false;
    addPPFile(S, "h_incl_sysRel_bPurity");
    if(purSplit){ std::vector<double> r = relFromFile(fPSpp, "h_sysRel_uncorr"); for(int i = 0; i < NB; i++) S.u[0][i] = r[i]*ppGenOverReco[i]; }
    for(int s = 1; s < NS; s++){
      TH1D *rel = nullptr; fPS->GetObject(Form(purSplit ? "h_sysRel_uncorr_%s" : "h_sysRel_total_%s", sysTag[s]), rel);
      if(!rel){ printf("ERROR: purity systematic missing for %s\n", sysTag[s]); return; }
      double m[NB] = {0.};
      for(int k = -1; k <= 1; k += 2){
        TH1D *mv = (TH1D*) meas[s]->Clone("mv"); mv->SetDirectory(nullptr);
        for(int b = 1; b <= NB; b++){
          // bins without a fit keep the chain's buffer treatment: the sys file has 0 there
          double f = 1. + k*rel->GetBinContent(b);
          mv->SetBinContent(b, meas[s]->GetBinContent(b)*f); mv->SetBinError(b, meas[s]->GetBinError(b)*f);
        }
        TH1D *u = unfoldB(R[s], mv, nIterB[s], "uv"); std::vector<double> d = shift(u, unfNom[s]);
        for(int i = 0; i < NB; i++) m[i] = TMath::Max(m[i], fabs(d[i]));
        delete u; delete mv;
      }
      for(int i = 0; i < NB; i++) S.u[s][i] = m[i];
    }
    // saved for the purity extraction figure (2026-11-17): the unfolded magnitude of the uncorrelated part per system
    for(int s = 0; s < NS; s++){ TH1D *h = new TH1D(Form("h_bPur_uncorr_%s", sysTag[s]), "", NB, edges); h->SetDirectory(nullptr);
      for(int b = 1; b <= NB; b++) h->SetBinContent(b, S.u[s][b-1]); iterSave.push_back(h); }
    srcs.push_back(S);
  }

  // b-tag working point (2026-10-28, user: CSVv2 loose nominal, medium "an alternative for systematic
  // variation"): the whole chain rerun with the medium-cut purity; signed shift of every system's b per Z,
  // correlated across systems (the same choice everywhere), so it partly cancels in PbPb / pp
  if(modePF && altWPChainTag != ""){
    TFile *fW = TFile::Open(Form("/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/BJetPerZ_PFJets/bJetChain_PFJets_PbPb_dedicatedPurity%s.root", altWPChainTag.Data()));
    if(!fW || fW->IsZombie()){ printf("ERROR: alternative working-point chain file %s\n", altWPChainTag.Data()); return; }
    Source S; S.name = "bTagWP"; S.label = "#it{b}-tag working point"; S.corr = true;
    VV d = zeroVV();
    for(int s = 0; s < NS; s++){
      TH1D *hn = nullptr, *ha = nullptr; fC->GetObject(Form("bJetsPerZ_%s", sysTag[s]), hn); fW->GetObject(Form("bJetsPerZ_%s", sysTag[s]), ha);
      if(!hn || !ha){ printf("ERROR: bJetsPerZ_%s missing for bTagWP\n", sysTag[s]); return; }
      for(int i = 0; i < NB; i++) if(reported(i + 1) && hn->GetBinContent(i + 1) > 0.) d[s][i] = ha->GetBinContent(i + 1)/hn->GetBinContent(i + 1) - 1.;
    }
    S.d.push_back(d);
    srcs.push_back(S);
  }

  // JEU: gen-level shifted/nominal from the JEU files, at the b chain's N (b side)
  // and the inclusive's N (inclusive side); corr, pp-PbPb per jeuPPCorr
  {
    Source S; S.name = "JEU"; S.label = jeuPPCorr ? "jet energy scale (pp-PbPb corr.)" : "jet energy scale (pp-PbPb uncorr.)"; S.corr = true; S.ppUncorr = !jeuPPCorr;
    const char *dir[2] = {"up", "down"};
    // PFJets (2026-10-21, user): the b side is the MEASURED b spectrum shifted by +/- the JEC
    // uncertainty (headers/functions/jeuShift.h) and unfolded with the b response at its N -- PbPb here,
    // pp from the pp chain (h_unf_varRel_JEUshift_*); the inclusive side is the inclusive chain's own
    // signed shift on its MEASURED spectrum, with the same analytic method
    // (h_varRel_JEUdata_*_bjet of the inclusive file), so a common shift cancels in b / incl. Before: both sides from the all-flavour
    // simulated gen-level ratios at N_b and N_incl, so b / incl measured the iteration difference.
    const bool jeuOnB = modePF;
    std::vector<double> dAA;
    TFile *fPUj = nullptr;
    if(jeuOnB){
      dAA = jeuDelta("/home/clayton/Analysis/code/bJetRaaAnalysis/JetEnergyCorrections/Autumn18_HI_V8_MC_Uncertainty_AK4PF.txt");
      fPUj = TFile::Open(Form("/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/CorrectedBJetSpectra/Data/unfoldedBJetSpectrum_PFJets_pp%s.root", ppPurVar()));
      if(!fPUj || fPUj->IsZombie()){ printf("ERROR: pp unfolded file missing\n"); return; }
    }
    auto signedFromFile = [&](TFile *f, const char *name){
      std::vector<double> v(NB, 0.); TH1D *h = nullptr; f->GetObject(name, h);
      if(!h || h->GetNbinsX() != NB){ printf("ERROR: %s missing or binned wrongly in %s\n", name, f->GetName()); return v; }
      for(int i = 1; i <= NB; i++) v[i-1] = h->GetBinContent(i);
      return v; };
    for(int k = 0; k < 2; k++){
      VV db = zeroVV(), di = zeroVV();
      if(jeuOnB){
        for(int s = 0; s < NS; s++){
          if(s == 0) db[s] = signedFromFile(fPUj, Form("h_unf_varRel_JEUshift_%s", dir[k]));
          else {
            TH1D *m = jeuShiftMeasured(meas[s], dAA, k == 0 ? +1 : -1, "jeuMeas");   // fit-based shape since 2026-11-16 (jeuShift.h; BJET_JEU_SHAPE=measured = before)
            TH1D *u = unfoldB(R[s], m, nIterB[s], "jeuU"); db[s] = shift(u, unfNom[s]);
            // saved for the extraction figure (plotBJetJESSys_PF.C, 2026-11-16): the measured b spectrum,
            // its shifted version, and the unfolded shifted spectrum (the nominal is h_bIter_nom_<sys>)
            if(k == 0){ TH1D *c = (TH1D*) meas[s]->Clone(Form("h_bJEU_measNom_%s", sysTag[s])); c->SetDirectory(nullptr); iterSave.push_back(c); }
            { TH1D *c = (TH1D*) m->Clone(Form("h_bJEU_meas_%s_%s", dir[k], sysTag[s])); c->SetDirectory(nullptr); iterSave.push_back(c);
              TH1D *c2 = (TH1D*) u->Clone(Form("h_bJEU_unf_%s_%s", dir[k], sysTag[s])); c2->SetDirectory(nullptr); iterSave.push_back(c2); }
            delete m; delete u;
          }
          di[s] = signedFromFile(fI, Form("h_varRel_JEUdata_%s_%s_bjet", dir[k], sysTag[s]));   // same method as the b side
          // signed b and inclusive shifts on the reported axis, every system (pp from its chain), for the figure
          for(int side = 0; side < 2; side++){
            TH1D *h = new TH1D(Form("h_varRel_JEU_%s_%s_%s", dir[k], side == 0 ? "b" : "incl", sysTag[s]), "", NB, edges); h->SetDirectory(nullptr);
            for(int i = 1; i <= NB; i++) h->SetBinContent(i, side == 0 ? db[s][i-1] : di[s][i-1]);
            iterSave.push_back(h);
          }
        }
        S.d.push_back(db); S.di.push_back(di);
        continue;
      }
      for(int s = 0; s < NS; s++){
        TFile *fj = s == 0 ? fJP : fJA;
        auto get = [&](const char *what, int N){ TH1D *h = nullptr; fj->GetObject(s == 0 ? Form("gen_%s_N%d", what, N) : Form("gen_%s_N%d_%s", what, N, sysTag[s]), h); return h; };
        TH1D *n0 = get("nominal", nIterB[s]), *v0 = get(dir[k], nIterB[s]), *n1 = get("nominal", nIterIncl[s]), *v1 = get(dir[k], nIterIncl[s]);
        if(!n0 || !v0 || !n1 || !v1){ printf("ERROR: JEU inputs missing for %s\n", sysTag[s]); return; }
        TH1D *N0 = rebinEdges(n0, "n0"), *V0 = rebinEdges(v0, "v0"), *N1 = rebinEdges(n1, "n1"), *V1 = rebinEdges(v1, "v1");
        db[s] = shift(V0, N0); di[s] = shift(V1, N1);
        delete N0; delete V0; delete N1; delete V1;
      }
      S.d.push_back(db); S.di.push_back(di);
    }
    srcs.push_back(S);
  }

  // JER: PbPb response smeared (corr within PbPb); pp from the pp chain's smeared b response
  // (h_unf_varRel_JERresp, 2026-10-21; before, the pp b spectrum had none); pp-PbPb CORRELATED, ratio per
  // variation (2026-10-22, user; uncorrelated 2026-10-21 to 10-22); inclusive side unsigned from the inclusive file
  {
    Source S; S.name = "JER"; S.label = modePF ? jerSmear::varLabel() : "jet energy resolution (+20%)"; S.corr = true;
    VV d = zeroVV();
    if(modePF){
      TFile *fPUj = TFile::Open(Form("/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/CorrectedBJetSpectra/Data/unfoldedBJetSpectrum_PFJets_pp%s.root", ppPurVar()));
      TH1D *hj = fPUj ? (TH1D*) fPUj->Get("h_unf_varRel_JERresp") : nullptr;
      if(!hj || hj->GetNbinsX() != NB){ printf("ERROR: h_unf_varRel_JERresp missing in the pp unfolded file\n"); return; }
      for(int i = 1; i <= NB; i++) d[0][i-1] = hj->GetBinContent(i);
    }
    // PFJets, 2026-11-05 (user): measured on the MC (jerSmear.h) -- the nominal response's reco spectrum
    // unfolded with the nominal and the varied response; BJET_JER_INPUT=data = the data unfolded with the
    // varied response. Both printed.
    if(modePF) printf("JER (PbPb b, %s), var / nom - 1, reported bins, MC | data:\n", jerSmear::inputText());
    for(int s = 1; s < NS; s++){
      TH2D *sm = smearResp(Rfine[s], "sm"); TH2D *Rs = rebin2D(sm, "Rs");
      TH1D *u = unfoldB(Rs, meas[s], nIterB[s], "ujer"); std::vector<double> dD = shift(u, unfNom[s]);
      d[s] = dD;
      if(modePF){
        TH1D *mx = R[s]->ProjectionX("jerMCx", 1, R[s]->GetNbinsY()); mx->SetDirectory(nullptr);
        TH1D *mA = bToAnaAxis(mx, NB + 1, edges, "jerMCa");
        TH1D *u0 = unfoldB(R[s], mA, nIterB[s], "ujerMC0"), *u1 = unfoldB(Rs, mA, nIterB[s], "ujerMC1");
        std::vector<double> dM = shift(u1, u0);
        // saved for the JER extraction figure (plotBJetJERSys_PF.C, 2026-11-17): widths and both signed shifts
        iterSave.push_back(jerSmear::resolution(Rfine[s], Form("h_bJER_resNom_%s", sysTag[s])));
        iterSave.push_back(jerSmear::resolution(sm, Form("h_bJER_resVar_%s", sysTag[s])));
        for(int side = 0; side < 2; side++){ TH1D *h = new TH1D(Form("h_bJER_%s_%s", side == 0 ? "varMC" : "varData", sysTag[s]), "", NB, edges); h->SetDirectory(nullptr);
          for(int b = 1; b <= NB; b++) h->SetBinContent(b, side == 0 ? dM[b-1] : dD[b-1]); iterSave.push_back(h); }
        printf("  %-4s", sysTag[s]); for(int b = 1; b <= NB; b++) if(reported(b) && edges[b] <= 200. + 1e-6) printf("  %3.0f: %+.3f | %+.3f", edges[b-1], dM[b-1], dD[b-1]); printf("\n");
        if(!jerSmear::onData()) d[s] = dM;
        delete mx; delete mA; delete u0; delete u1;
      }
      delete u; delete Rs; delete sm;
    }
    S.d.push_back(d);
    // inclusive side. PFJets since 2026-11-12 (user): the inclusive chain's SIGNED shift per system
    // (h_varRel_JER_<sys>_bjet), so in b / incl the ratio is formed per variation: correlated between b and
    // inclusive and between pp and PbPb. BJET_JER_DR=uncorr (and the calo mode) = the earlier unsigned
    // per-system magnitudes added in quadrature (b vs incl and incl pp vs PbPb uncorrelated).
    const char *jerDR = gSystem->Getenv("BJET_JER_DR");
    if(modePF && !(jerDR && TString(jerDR) == "uncorr")){
      VV di = zeroVV();
      for(int s = 0; s < NS; s++){
        TH1D *h = nullptr; fI->GetObject(Form("h_varRel_JER_%s_bjet", sysTag[s]), h);
        if(!h || h->GetNbinsX() != NB){
          printf("ERROR: h_varRel_JER_%s_bjet missing or binned wrongly in %s -- rerun inclusiveJetsPerZSys_PF.C"
                 " (2026-11-12 or later), or set BJET_JER_DR=uncorr for the old quadrature\n", sysTag[s], fI->GetName());
          return;
        }
        for(int i = 1; i <= NB; i++) di[s][i-1] = h->GetBinContent(i);
      }
      S.di.push_back(di);
      printf("JER on b / incl: correlated (signed inclusive shift)\n");
    } else {
      for(int s = 0; s < NS; s++) S.ui[s] = relFromFile(fI, Form("h_sysRel_JER_%s_bjet", sysTag[s]));
      printf("JER on b / incl: inclusive side in quadrature (BJET_JER_DR=uncorr)\n");
    }
    srcs.push_back(S);
  }

  // iterations (uncorr)
  {
    Source S; S.name = "iterations"; S.label = "iterations"; S.corr = false;
    addPPFile(S, "h_incl_sysRel_iterations");
    for(int s = 1; s < NS; s++){
      TH1D *u = unfoldB(R[s], meas[s], nIterB[s] + 1, "uit"); std::vector<double> d = shift(u, unfNom[s]);
      for(int i = 0; i < NB; i++) S.u[s][i] = fabs(d[i]);
      // saved for the extraction figure (src/systematics/plotBJetIterSys_PF.C, 2026-11-05)
      { TH1D *c = (TH1D*) u->Clone(Form("h_bIter_up_%s", sysTag[s])); c->SetDirectory(nullptr); c->SetTitle(Form("N = %d", nIterB[s] + 1)); iterSave.push_back(c);
        TH1D *c0 = (TH1D*) unfNom[s]->Clone(Form("h_bIter_nom_%s", sysTag[s])); c0->SetDirectory(nullptr); c0->SetTitle(Form("N = %d", nIterB[s])); iterSave.push_back(c0); }
      delete u;
      if(modePF && nIterB[s] > 1){   // full N +/- 1 envelope
        TH1D *v = unfoldB(R[s], meas[s], nIterB[s] - 1, "uitm"); std::vector<double> e = shift(v, unfNom[s]);
        for(int i = 0; i < NB; i++) S.u[s][i] = TMath::Max(S.u[s][i], fabs(e[i]));
        TH1D *c = (TH1D*) v->Clone(Form("h_bIter_dn_%s", sysTag[s])); c->SetDirectory(nullptr); c->SetTitle(Form("N = %d", nIterB[s] - 1)); iterSave.push_back(c);
        delete v;
      }
    }
    for(int s = 0; s < NS; s++) S.ui[s] = relFromFile(fI, Form("h_sysRel_iterations_%s_bjet", sysTag[s]));
    srcs.push_back(S);
  }

  // unfolding bias (uncorr; 2026-10-21, user): the same-sample closure at the class's DATA-LIKE tilt
  // (k fitted to the unfolded data, fitPriorTilt; measured only, the prior stays the MC one), as for the
  // inclusive jets: pseudo-data = the response with gen rows x (g/100)^-k, the data's relative errors,
  // unfolded with the nominal response and N, |unfolded / true - 1|. pp from the pp chain
  // (h_unf_sysRel_unfoldingBias); b / incl adds the inclusive part (h_sysRel_unfoldingBias_<sys>_bjet).
  if(modePF){
    Source S; S.name = "unfoldingBias"; S.label = "unfolding non-closure (tilted truth)"; S.corr = false;
    TFile *fPUb = TFile::Open(Form("/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/CorrectedBJetSpectra/Data/unfoldedBJetSpectrum_PFJets_pp%s.root", ppPurVar()));
    if(!fPUb || fPUb->IsZombie()){ printf("ERROR: pp unfolded file missing\n"); return; }
    S.u[0] = relFromFile(fPUb, "h_unf_sysRel_unfoldingBias");
    printf("\nunfolding bias: closure tilt k per class:");
    for(int s = 1; s < NS; s++){
      auto unfoldFine = [&](TH2D *fine, TH1D *&truthOut){
        TH2D *r2 = rebin2D(fine, "ubR");
        TH1D *U = unfoldOnExt<RooUnfoldResponse, RooUnfoldBayes>(r2, meas[s], NB + 1, edges, recoFloor, nIterB[s], "ubU", &truthOut);
        delete r2; return U; };
      const double k = bClosureTilt > -90. ? bClosureTilt : fitPriorTilt(Rfine[s], unfoldFine, 80., bPriorFitHi(), 8, false);   // 0.3 since 2026-11-05
      printf("  %s %+.3f", sysTag[s], k);
      TH2D *P = bClosureTilt > -90. ? closureTiltRows(Rfine[s], bClosureTilt, "ubP") : tiltGenRows(Rfine[s], k, "ubP");
      // pseudo measurement and truth on the analysis edges (all gen / all reco), data's relative errors
      TH1D *pd = (TH1D*) meas[s]->Clone("ubPd"); pd->Reset(); pd->SetDirectory(nullptr);
      TH1D *tt = (TH1D*) meas[s]->Clone("ubTt"); tt->Reset(); tt->SetDirectory(nullptr);
      for(int ix = 1; ix <= P->GetNbinsX(); ix++) for(int iy = 1; iy <= P->GetNbinsY(); iy++){
        const double c = P->GetBinContent(ix, iy); if(c == 0.) continue;
        pd->Fill(P->GetXaxis()->GetBinCenter(ix), c); tt->Fill(P->GetYaxis()->GetBinCenter(iy), c);
      }
      for(int b = 1; b <= NB; b++){ const double mc = meas[s]->GetBinContent(b); pd->SetBinError(b, mc > 0. ? pd->GetBinContent(b)*meas[s]->GetBinError(b)/mc : 0.); }
      TH1D *U = unfoldB(R[s], pd, nIterB[s], "ubUn");
      // reported bins only: outside 80-200 GeV the under/overflow axis returns placeholders
      for(int i = 1; i <= NB; i++){ const double t = tt->GetBinContent(i);
        S.u[s][i-1] = (t > 0. && bReportedBin(edges[i-1], edges[i])) ? fabs(U->GetBinContent(i)/t - 1.) : 0.; }
      delete P; delete pd; delete tt; delete U;
    }
    printf("\n");
    for(int s = 0; s < NS; s++) S.ui[s] = relFromFile(fI, Form("h_sysRel_unfoldingBias_%s_bjet", sysTag[s]));
    srcs.push_back(S);
  }

  // data-driven prior and reco floor (PFJets mode)
  if(modePF){
    TFile *fPU = TFile::Open(Form("/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/CorrectedBJetSpectra/Data/unfoldedBJetSpectrum_PFJets_pp%s.root", ppPurVar()));
    if(!fPU || fPU->IsZombie()){ printf("ERROR: pp unfolded file missing\n"); return; }
    Source P; P.name = "prior"; P.label = "unfolding prior (data-driven vs MC)"; P.corr = false;
    P.u[0] = relFromFile(fPU, "h_unf_sysRel_prior");
    const double floorAlt = 70.;   // 60 -> 70 since 2026-11-05 (70 -> 80 from 2026-10-08; 50 -> 100 before)
    Source F; F.name = "recoFloor"; F.label = "reco floor (60 #rightarrow 70 GeV)"; F.corr = true; F.ppUncorr = true;
    VV dF = zeroVV(); dF[0] = relFromFile(fPU, "h_unf_sysRel_recoFloor");
    for(int s = 1; s < NS; s++){
      TH2D *Rm = rebin2D(RfineRaw[s], "Rmc"); TH1D *u = unfoldB(Rm, meas[s], nIterB[s], "uPrior");
      std::vector<double> d = shift(u, unfNom[s]); for(int i = 0; i < NB; i++) P.u[s][i] = fabs(d[i]);
      delete u; delete Rm;
      if(!recoFloorSourceOn()) continue;   // no reco-floor source since 2026-11-16 (purityVariant.h)
      const double fl0 = recoFloor; recoFloor = floorAlt;
      TH1D *w = unfoldB(R[s], meas[s], nIterB[s], "uFloor"); recoFloor = fl0;
      dF[s] = shift(w, unfNom[s]); for(int i = 0; i < NB; i++) if(edges[i] < floorAlt - 1e-6) dF[s][i] = 0.;
      delete w;
    }
    F.d.push_back(dF);
    if(bUseDataDrivenPrior) srcs.push_back(P);   // with the MC prior (bJetUnfoldExt.h) there is no prior variation: nominal = MC
    if(recoFloorSourceOn()) srcs.push_back(F);   // off since 2026-11-16 (user): the non-closure at the 60 GeV floor stands for it
  }

  // response MC statistics, toys on the fine matrix (uncorr)
  {
    Source S; S.name = "responseStat"; S.label = "response MC stat."; S.corr = false;
    addPPFile(S, "h_incl_sysRel_responseStat");
    TRandom3 rnd(toySeed);
    for(int s = 1; s < NS; s++){
      std::vector<double> s1(NB, 0.), s2(NB, 0.);
      for(int t = 0; t < nToys; t++){
        TH2D *r2 = (TH2D*) Rfine[s]->Clone("toyR"); r2->SetDirectory(nullptr);
        for(int ix = 1; ix <= r2->GetNbinsX(); ix++) for(int iy = 1; iy <= r2->GetNbinsY(); iy++){
          double c = Rfine[s]->GetBinContent(ix, iy), e = Rfine[s]->GetBinError(ix, iy);
          if(c <= 0. || e <= 0.) continue;
          double ne = (c/e)*(c/e);
          r2->SetBinContent(ix, iy, rnd.Poisson(ne)*e*e/c);
        }
        TH2D *Rt = rebin2D(r2, "toyRb"); TH1D *u = unfoldB(Rt, meas[s], nIterB[s], "toyU");
        std::vector<double> d = shift(u, unfNom[s]);
        for(int i = 0; i < NB; i++){ s1[i] += 1. + d[i]; s2[i] += (1. + d[i])*(1. + d[i]); }
        delete u; delete Rt; delete r2;
      }
      for(int i = 0; i < NB; i++){ double m = s1[i]/nToys; S.u[s][i] = sqrt(TMath::Max(0., s2[i]/nToys - m*m)); }
      printf("  response toys done: %s\n", sysTag[s]);
    }
    for(int s = 0; s < NS; s++) S.ui[s] = relFromFile(fI, Form("h_sysRel_responseStat_%s_bjet", sysTag[s]));
    srcs.push_back(S);
  }

  // muon-tag frequency: the same f_pp on every system (chain), so its pp-chain
  // uncertainty is one correlated shift; cancels in PbPb/pp and R_CP
  {
    Source S; S.name = "muonTagFreq"; S.label = "muon-tag frequency (common f_{pp})"; S.corr = true;
    std::vector<double> u = relFromFile(fP, "h_incl_sysRel_muonTagFrequency");
    for(int k = -1; k <= 1; k += 2){ VV d = zeroVV(); for(int s = 0; s < NS; s++) for(int i = 0; i < NB; i++) d[s][i] = k*u[i]; S.d.push_back(d); }
    srcs.push_back(S);
  }

  // calo-match ID, MC: the ID-free response's row scaling at f = 0.5 / 1.5
  // (corr within PbPb, quadrature with pp; pp from its chain's turnOn slot)
  if(!modePF){
    Source S; S.name = "caloMatchID"; S.label = "calo-match ID: response (#pm50%)"; S.corr = true; S.ppUncorr = true;
    std::vector<double> ppU = relFromFile(fP, "h_incl_sysRel_turnOn");
    for(double fr : {0.5, 1.5}){
      VV d = zeroVV(); d[0] = ppU;
      for(int s = 1; s < NS; s++){
        TH2D *rf = idScaleRows(RfineRaw[s], epsMC[s], fr, "rfID"), *rb = rebin2D(rf, "rbID");
        TH1D *u = unfoldB(rb, meas[s], nIterB[s], "uID"); d[s] = shift(u, unfNom[s]);
        delete u; delete rb; delete rf;
      }
      S.d.push_back(d);
    }
    srcs.push_back(S);
  }
  // calo-match ID, data: PbPb data / eps_data (corr within PbPb)
  if(!modePF){
    Source S; S.name = "caloMatchIDdata"; S.label = "calo-match ID: data efficiency"; S.corr = true;
    VV d = zeroVV();
    const char *dn = "h_inclRecoJetPt_inclRecoMuonTag_triggerOn";
    printf("eps_data (calo-matched / plain muon-tagged data, normalised 150-300 GeV) on the b-jet bins:\n");
    for(int s = 1; s < NS; s++){
      TH1D *m = coarseSum1(fDM, dn, s-1, "dIDm"), *p = coarseSum1(fDP, dn, s-1, "dIDp");
      if(!m || !p){ printf("ERROR: data ID spectra missing for %s\n", sysTag[s]); return; }
      auto I = [](TH1D *h, double lo, double hi){ return h->Integral(h->FindBin(lo + 1e-6), h->FindBin(hi - 1e-6)); };
      const double nm = I(m, 150., 300.), np = I(p, 150., 300.);
      TH1D *mv = (TH1D*) meas[s]->Clone("mvID"); mv->SetDirectory(nullptr);
      printf("  %s:", sysTag[s]);
      for(int b = 1; b <= NB; b++){
        double ip = I(p, edges[b-1], edges[b]), e = ip > 0. ? (I(m, edges[b-1], edges[b])/nm)/(ip/np) : 1.;
        e = TMath::Min(1., TMath::Max(e, 0.05));
        mv->SetBinContent(b, meas[s]->GetBinContent(b)/e); mv->SetBinError(b, meas[s]->GetBinError(b)/e);
        if(edges[b-1] >= recoFloor - 1e-6) printf("  %3.0f-%3.0f %.3f", edges[b-1], edges[b], e);
      }
      printf("\n");
      TH1D *u = unfoldB(R[s], mv, nIterB[s], "uIDd"); d[s] = shift(u, unfNom[s]);
      delete u; delete mv;
    }
    S.d.push_back(d);
    srcs.push_back(S);
  }

  // single-muon efficiency (uncorr)
  {
    Source S; S.name = "muonEff"; S.label = "muon efficiency"; S.corr = false;
    for(int s = 0; s < NS; s++) for(int i = 0; i < NB; i++) S.u[s][i] = muEffErr[s]/muEff[s];
    srcs.push_back(S);
  }

  // N_Z sources. Published N_Z (nominal): pubZYield (published yield) and pubZAcc (pp
  // |y_Z| < 2.1 fraction, pp only), from the inclusive file, N_Z-only like NZstat
  if(pubNZ){
    for(const char *nm : {"pubZYield", "pubZAcc"}){
      Source A; A.name = nm; A.label = TString(nm) == "pubZYield" ? "#it{N}_{Z} published yield" : "#it{N}_{Z} pp |#it{y}_{Z}| < 2.1 fraction";
      A.corr = false; A.nzOnly = true;
      for(int s = 0; s < NS; s++) A.u[s] = relFromFile(fI, Form("h_sysRel_%s_%s_bjet", nm, sysTag[s]));
      srcs.push_back(A);
    }
  } else
  // N_Z statistics (uncorr, b per Z only) and mass window (ratios of b per Z only)
  {
    Source A; A.name = "NZstat"; A.label = "#it{N}_{Z} stat."; A.corr = false; A.nzOnly = true;
    for(int s = 0; s < NS; s++) A.u[s] = relFromFile(fI, Form("h_sysRel_NZstat_%s_bjet", sysTag[s]));
    srcs.push_back(A);
    // same-sign Z background (2026-10-21): N_Z-only like NZstat, cancels in b / inclusive;
    // it was in the inclusive chain but never propagated here
    Source Z; Z.name = "zBkg"; Z.label = "#it{N}_{Z} background (same sign)"; Z.corr = false; Z.nzOnly = true;
    for(int s = 0; s < NS; s++) Z.u[s] = relFromFile(fI, Form("h_sysRel_zBkg_%s_bjet", sysTag[s]));
    srcs.push_back(Z);
    Source M; M.name = "massWindow"; M.label = "#it{N}_{Z} mass window"; M.ratioOnly = true;
    for(int q = 5; q < 12; q++){
      int a, b; bool bi; qSys(q, a, b, bi);
      M.ratioVal.push_back(relFromFile(fI, Form("h_sysRel_massWindow_%s_C%d_bjet", b == 0 ? "PbPbOverPP" : "RCP", a)));
    }
    srcs.push_back(M);
  }

  // inclusive-only sources, b / inclusive only (uncorr across the two systems)
  {
    const char *nm[4] = {"fakeJets", "stitchWindow", "stitchNormStat", "recoFloor"};
    const char *lb[4] = {"incl.: fake-jet subtraction", "incl.: stitching windows", "incl.: stitch normalization stat.", "incl.: reco floor"};
    for(int k = 0; k < 4; k++){
      if(k == 3 && !recoFloorSourceOn()) continue;   // no reco-floor source since 2026-11-16 (purityVariant.h)
      TH1D *probe = nullptr; fI->GetObject(Form("h_sysRel_%s_%s_bjet", nm[k], sysTag[1]), probe);
      if(!probe && k == 0){ printf("  no inclusive fakeJets source (forest CS, no leakage, 2026-11-05): skipped\n"); continue; }
      Source S; S.name = Form("incl_%s", nm[k]); S.label = lb[k]; S.corr = false;
      for(int s = 0; s < NS; s++) S.ui[s] = relFromFile(fI, Form("h_sysRel_%s_%s_bjet", nm[k], sysTag[s]));
      srcs.push_back(S);
    }
  }

  // ---- outputs ---------------------------------------------------------------------
  TFile *fo = TFile::Open(Form("%s%s%s.root", outDir, outName.Data(), sfx.Data()), "recreate");
  for(auto h : iterSave) h->Write();
  const int NSrc = srcs.size();
  for(int q = 0; q < NQ; q++){
    int a, b; bool bi; qSys(q, a, b, bi);
    TString nomName = b < 0 ? Form("bJetsPerZ_%s", sysTag[a]) : Form("r_%s_%s_C%d", bi ? "bOverIncl" : "b", b == 0 ? "PbPbOverPP" : "CP", a);
    TH1D *n0 = nullptr; fC->GetObject(nomName, n0);
    if(!n0){ printf("ERROR: %s missing in the chain file\n", nomName.Data()); return; }
    if(b < 0 && bi){   // the fraction: b per Z over inclusive per Z, same system, both per GeV
      TH1D *inc = nullptr; fC->GetObject(Form("jetsPerZ_%s", sysTag[a]), inc);
      if(!inc){ printf("ERROR: jetsPerZ_%s missing in the chain file\n", sysTag[a]); return; }
      TH1D *fr = (TH1D*) n0->Clone(Form("frac_%s", sysTag[a])); fr->SetDirectory(nullptr); fr->Divide(inc);
      for(int i = 1; i <= NB; i++) if(n0->GetBinContent(i) <= 0. || inc->GetBinContent(i) <= 0.){ fr->SetBinContent(i, 0.); fr->SetBinError(i, 0.); }
      n0 = fr;
    }
    TH1D *v = (TH1D*) n0->Clone(Form("h_nom_%s", qTag(q).Data())); v->SetDirectory(nullptr);
    v->SetTitle(Form("%s;#it{p}_{T}^{jet} [GeV];%s", qLabel(q).Data(), b < 0 ? (bi ? "#it{b}-jet fraction" : "#it{b} jets per #it{N}_{Z} per GeV") : "ratio"));
    TH1D *st = new TH1D(Form("h_statRel_%s", qTag(q).Data()), "stat.;#it{p}_{T}^{jet} [GeV];relative uncertainty", NB, edges);
    TH1D *to = new TH1D(Form("h_sysRel_total_%s", qTag(q).Data()), "total systematic;#it{p}_{T}^{jet} [GeV];relative uncertainty", NB, edges);
    std::vector<TH1D*> hs;
    for(int k = 0; k < NSrc; k++)
      hs.push_back(new TH1D(Form("h_sysRel_%s_%s", srcs[k].name.Data(), qTag(q).Data()), Form("%s;#it{p}_{T}^{jet} [GeV];relative uncertainty", srcs[k].label.Data()), NB, edges));
    for(int i = 1; i <= NB; i++){
      if(v->GetBinContent(i) <= 0. || !reported(i)) continue;
      st->SetBinContent(i, v->GetBinError(i)/v->GetBinContent(i));
      for(int k = 0; k < NSrc; k++) hs[k]->SetBinContent(i, relOf(srcs[k], q, i-1));
    }
    // a source with no pT trend is replaced by its mean (headers/plotting/flattenSystematic.h)
    if(flattenSources) for(int k = 0; k < NSrc; k++){
      if(noFlattenJEU && srcs[k].name == "JEU") continue;
      if(noFlattenPurity && srcs[k].name.BeginsWith("bPurity")) continue;   // purity sources kept bin by bin (2026-11-05, user)
      if(noFlattenIterations && srcs[k].name == "iterations") continue;   // iterations kept bin by bin (2026-11-05, user)
      if(noFlattenJER && srcs[k].name == "JER") continue;   // JER kept bin by bin (2026-11-17, user)
      if(srcs[k].name == "bTagWP"){   // 2026-11-05 (user): always flat, the mean of the bin-by-bin values over the reported 80-200 GeV
        double sum = 0.; int n = 0;
        for(int b = 1; b <= hs[k]->GetNbinsX(); b++){ const double lo = hs[k]->GetBinLowEdge(b); if(lo >= bWPMeanLo - 1e-6 && lo < bWPMeanHi - 1e-6 && hs[k]->GetBinContent(b) > 0.){ sum += hs[k]->GetBinContent(b); n++; } }
        if(n > 0){ const double m = sum/n; for(int b = 1; b <= hs[k]->GetNbinsX(); b++) if(hs[k]->GetBinContent(b) > 0.) hs[k]->SetBinContent(b, m);
          printf("  bTagWP flat     %-36s mean over %.0f-%.0f GeV %5.1f%%\n", qLabel(q).Data(), bWPMeanLo, bWPMeanHi, 100*m); }
        continue;
      }
      double mean, F; if(flattenIfNoTrend(hs[k], mean, F)) printf("  flattened %-14s %-36s mean %5.1f%%  F = %.1f\n", srcs[k].name.Data(), qLabel(q).Data(), 100*mean, F);
    }
    for(int i = 1; i <= NB; i++){
      if(v->GetBinContent(i) <= 0. || !reported(i)) continue;
      double t2 = 0.;
      for(int k = 0; k < NSrc; k++){ double r = hs[k]->GetBinContent(i); t2 += r*r; }
      to->SetBinContent(i, sqrt(t2));
    }
    fo->cd(); v->Write(); st->Write(); to->Write(); for(auto h : hs) h->Write();
  }

  // ---- printout --------------------------------------------------------------------
  for(int q = 0; q < NQ; q++){
    printf("\n=== %s (%%) ===\n%-8s %6s", qLabel(q).Data(), "pT", "stat");
    for(int k = 0; k < NSrc; k++) printf(" %8.8s", srcs[k].name.Data());
    printf(" %7s\n", "TOTAL");
    TH1D *st = (TH1D*) fo->Get(Form("h_statRel_%s", qTag(q).Data())), *to = (TH1D*) fo->Get(Form("h_sysRel_total_%s", qTag(q).Data()));
    for(int i = 1; i <= NB; i++){
      if(to->GetBinContent(i) <= 0.) continue;
      printf("%3.0f-%-4.0f %6.1f", edges[i-1], edges[i], 100*st->GetBinContent(i));
      for(int k = 0; k < NSrc; k++) printf(" %8.1f", 100*((TH1D*) fo->Get(Form("h_sysRel_%s_%s", srcs[k].name.Data(), qTag(q).Data())))->GetBinContent(i));
      printf(" %7.1f\n", 100*to->GetBinContent(i));
    }
  }

  // ---- figures ---------------------------------------------------------------------
  const char *srcHex[16]  = {"#E69F00", "#D55E00", "#CC79A7", "#0072B2", "#56B4E9", "#009E73", "#999999", "#000000", "#E69F00", "#CC79A7", "#0072B2", "#56B4E9", "#009E73", "#D55E00", "#999999", "#E69F00"};
  const int   srcStyle[16] = {1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 7, 7, 7, 7, 7};
  auto breakdown = [&](int q0, int nq, int nx, int ny, const char *file){
    const double legH = ny > 1 ? 0.16 : (nx == 1 ? 0.34 : 0.28);
    TCanvas *cv = new TCanvas(Form("cvB_%s", file), "", 520*nx, (int)(460*ny/(1. - legH)));
    TLegend *lg = nullptr;
    for(int p = 0; p < nq; p++){
      int q = q0 + p, col = p % nx, row = p / nx;
      double x1 = (double) col/nx, x2 = (double)(col + 1)/nx;
      double y2 = 1. - (1. - legH)*row/ny, y1 = 1. - (1. - legH)*(row + 1)/ny;
      cv->cd();
      TPad *pd = new TPad(Form("pd_%s_%d", file, p), "", x1, y1, x2, y2); pd->Draw(); pd->cd();
      pd->SetLeftMargin(0.14); pd->SetRightMargin(0.03); pd->SetTopMargin(0.09); pd->SetBottomMargin(0.13);
      TH1D *to = (TH1D*) fo->Get(Form("h_sysRel_total_%s", qTag(q).Data()))->Clone();
      TH1D *st = (TH1D*) fo->Get(Form("h_statRel_%s", qTag(q).Data()))->Clone();
      double ymax = 0.;
      for(int i = 1; i <= NB; i++) if(reported(i) && edges[i-1] < ptDrawMax - 1e-6) ymax = TMath::Max(ymax, TMath::Max(to->GetBinContent(i), st->GetBinContent(i)));
      TH1F *fr = pd->DrawFrame(ptReportMin, 0., ptDrawMax, 100.*ymax*1.1);
      fr->GetXaxis()->SetTitle("#it{p}_{T}^{jet} [GeV]"); fr->GetYaxis()->SetTitle("relative uncertainty [%]");
      fr->GetXaxis()->SetTitleSize(0.050); fr->GetYaxis()->SetTitleSize(0.050);
      fr->GetXaxis()->SetLabelSize(0.045); fr->GetYaxis()->SetLabelSize(0.045); fr->GetYaxis()->SetTitleOffset(1.3);
      for(TH1D *h : {to, st}){ h->Scale(100.); for(int i = 0; i <= NB+1; i++) h->SetBinError(i, 0.); h->GetXaxis()->SetRangeUser(ptReportMin, ptDrawMax); }
      styleLine(to, "#000000", 5); to->Draw("HIST same");
      if(p == 0){ lg = makeLegend(0.03, 0.05, 0.99, 0.95); lg->SetTextFont(43); lg->SetTextSize(nx == 1 ? 15 : 18); lg->SetNColumns(nx == 1 ? 2 : 4); lg->AddEntry(to, "total syst.", "l"); }
      for(int k = 0; k < NSrc; k++){
        TH1D *h = (TH1D*) fo->Get(Form("h_sysRel_%s_%s", srcs[k].name.Data(), qTag(q).Data()))->Clone();
        bool any = false; for(int i = 1; i <= NB; i++) if(h->GetBinContent(i) > 0.) any = true;
        h->Scale(100.); for(int i = 0; i <= NB+1; i++) h->SetBinError(i, 0.);
        h->GetXaxis()->SetRangeUser(ptReportMin, ptDrawMax);
        styleLine(h, srcHex[k % 16], 2); h->SetLineStyle(srcStyle[k % 16]);
        if(any) h->Draw("HIST same");
        if(p == 0 && any) lg->AddEntry(h, srcs[k].label, "l");
      }
      styleLine(st, "#000000", 2); st->SetLineStyle(2); st->Draw("HIST same");
      if(p == 0) lg->AddEntry(st, "stat.", "l");
      TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.050);
      la.DrawLatex(0.17, 0.935, Form("%s, %s", qLabel(q).Data(), collText));
    }
    cv->cd();
    TPad *pl = new TPad(Form("pl_%s", file), "", 0., 0., 1., legH); pl->Draw(); pl->cd();
    lg->Draw();
    savePdfTight(cv, Form("%s%s%s.pdf", figSysDir, file, sfx.Data()));
    delete cv;
  };
  breakdown(0, 5, 3, 2, (outName + "_b_spectra").Data());
  breakdown(1, 1, 1, 1, (outName + "_b_C1").Data());        // single panel for talks
  breakdown(5, 4, 2, 2, (outName + "_b_PbPbOverPP").Data());
  breakdown(9, 3, 3, 1, (outName + "_b_RCP").Data());
  breakdown(12, 4, 2, 2, (outName + "_bOverIncl_PbPbOverPP").Data());
  breakdown(16, 3, 3, 1, (outName + "_bOverIncl_RCP").Data());
  breakdown(19, 5, 3, 2, (outName + "_bFraction").Data());

  // ratios with stat. bars and systematic boxes, the inclusive panels' layout
  auto ratioFig = [&](int q0, int nq, int nx, int ny, const char *ytitle, const char *file){
    TCanvas *cv = new TCanvas(Form("cvR_%s", file), "", 520*nx, 460*ny);
    cv->Divide(nx, ny, 0.001, 0.001);
    for(int p = 0; p < nq; p++){
      int q = q0 + p;
      cv->cd(p + 1);
      gPad->SetLeftMargin(0.14); gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.08); gPad->SetBottomMargin(0.13);
      TH1D *v  = (TH1D*) fo->Get(Form("h_nom_%s", qTag(q).Data()))->Clone();
      TH1D *to = (TH1D*) fo->Get(Form("h_sysRel_total_%s", qTag(q).Data()));
      double ymax = 0.;
      for(int i = 1; i <= NB; i++) if(reported(i) && edges[i-1] < ptDrawMax - 1e-6 && v->GetBinContent(i) > 0.)
        ymax = TMath::Max(ymax, v->GetBinContent(i)*(1. + to->GetBinContent(i)) + v->GetBinError(i));
      TH1F *fr = gPad->DrawFrame(ptReportMin, 0., ptDrawMax, TMath::Max(1.3, 1.35*ymax));
      fr->GetXaxis()->SetTitle("#it{p}_{T}^{jet} [GeV]"); fr->GetYaxis()->SetTitle(ytitle);
      fr->GetXaxis()->SetTitleSize(0.050); fr->GetYaxis()->SetTitleSize(0.050);
      fr->GetXaxis()->SetLabelSize(0.045); fr->GetYaxis()->SetLabelSize(0.045); fr->GetYaxis()->SetTitleOffset(1.3);
      TLine *one = new TLine(ptReportMin, 1., ptDrawMax, 1.); one->SetLineStyle(7); one->SetLineColor(kGray+2); one->Draw();
      const int col = TColor::GetColor("#D55E00");
      for(int i = 1; i <= NB; i++){
        double y = v->GetBinContent(i); if(y <= 0. || !reported(i) || edges[i-1] >= ptDrawMax - 1e-6) continue;
        TBox *bx = new TBox(edges[i-1], y*(1. - to->GetBinContent(i)), edges[i], y*(1. + to->GetBinContent(i)));
        bx->SetFillColorAlpha(col, 0.30); bx->SetLineColor(col); bx->Draw();
      }
      for(int i = 1; i <= NB; i++) if(!reported(i) || v->GetBinContent(i) <= 0.){ v->SetBinContent(i, -999.); v->SetBinError(i, 0.); }
      styleH(v, "#000000", markFilledCircle, 1.1); v->Draw("E1 X0 same");
      TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.050);
      la.DrawLatex(0.17, 0.935, Form("%s, %s", qLabel(q).Data(), collText));
      if(p == 0){
        TLegend *lg = makeLegend(0.17, 0.74, 0.70, 0.88, 0.042);
        lg->AddEntry(v, "stat.", "pe");
        TBox *lb = new TBox(); lb->SetFillColorAlpha(col, 0.30); lb->SetLineColor(col);
        lg->AddEntry(lb, jeuPPCorr ? "total syst. (JEU pp-PbPb corr.)" : "total syst. (JEU pp-PbPb uncorr.)", "f"); lg->Draw();
      }
    }
    savePdfTight(cv, Form("%s%s%s.pdf", figDir.Data(), file, sfx.Data()));
    delete cv;
  };
  // the b spectra themselves: one panel per system (already per Z per GeV),
  // stat. bars and systematic boxes, pp repeated as open markers in the PbPb panels
  {
    TCanvas *cv = new TCanvas("cvSpecB", "", 520*3, 460*2);
    cv->Divide(3, 2, 0.001, 0.001);
    const int col = TColor::GetColor("#D55E00");
    TH1D *ppD = nullptr;
    for(int q = 0; q < 5; q++){
      cv->cd(q + 1); gPad->SetLogy();
      gPad->SetLeftMargin(0.16); gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.08); gPad->SetBottomMargin(0.13);
      TH1D *v  = (TH1D*) fo->Get(Form("h_nom_%s", qTag(q).Data()))->Clone(Form("specB_%d", q)); v->SetDirectory(nullptr);
      TH1D *to = (TH1D*) fo->Get(Form("h_sysRel_total_%s", qTag(q).Data()));
      auto shown = [&](TH1D *h, int i){ return reported(i) && edges[i-1] < ptDrawMax - 1e-6 && h->GetBinContent(i) > 0.; };
      double ymax = 0., ymin = 1e30;
      for(int i = 1; i <= NB; i++) if(shown(v, i)){ ymax = TMath::Max(ymax, v->GetBinContent(i)*(1. + to->GetBinContent(i))); ymin = TMath::Min(ymin, v->GetBinContent(i)*(1. - to->GetBinContent(i))); }
      if(q == 0){ ppD = (TH1D*) v->Clone("specB_ppRef"); ppD->SetDirectory(nullptr); }
      else for(int i = 1; i <= NB; i++) if(shown(ppD, i)){ ymax = TMath::Max(ymax, ppD->GetBinContent(i)); ymin = TMath::Min(ymin, ppD->GetBinContent(i)); }
      TH1F *fr = gPad->DrawFrame(ptReportMin, 0.3*ymin, ptDrawMax, 5.*ymax);
      fr->GetXaxis()->SetTitle("#it{p}_{T}^{jet} [GeV]"); fr->GetYaxis()->SetTitle("#it{b} jets per #it{N}_{Z} per GeV");
      fr->GetXaxis()->SetTitleSize(0.050); fr->GetYaxis()->SetTitleSize(0.050);
      fr->GetXaxis()->SetLabelSize(0.045); fr->GetYaxis()->SetLabelSize(0.045); fr->GetYaxis()->SetTitleOffset(1.5);
      for(int i = 1; i <= NB; i++){
        if(!shown(v, i)) continue; double y = v->GetBinContent(i);
        TBox *bx = new TBox(edges[i-1], y*(1. - to->GetBinContent(i)), edges[i], y*(1. + to->GetBinContent(i)));
        bx->SetFillColorAlpha(col, 0.30); bx->SetLineColor(col); bx->Draw();
      }
      if(q > 0){
        TH1D *pr = (TH1D*) ppD->Clone(Form("specB_ppRef_%d", q)); pr->SetDirectory(nullptr);
        for(int i = 1; i <= NB; i++) if(!shown(ppD, i)){ pr->SetBinContent(i, -999.); pr->SetBinError(i, 0.); }
        styleH(pr, "#999999", markOpenSquare, 1.0); pr->Draw("E1 X0 same");
      }
      for(int i = 1; i <= NB; i++) if(!shown(v, i)){ v->SetBinContent(i, -999.); v->SetBinError(i, 0.); }
      styleH(v, "#000000", markFilledCircle, 1.1); v->Draw("E1 X0 same");
      TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.050);
      la.DrawLatex(0.19, 0.935, Form("%s, #it{b} jets, %s", qLabel(q).Data(), q == 0 || modePF ? "PF" : "calo-matched PF"));
    }
    cv->cd(6);
    TLegend *lg = makeLegend(0.08, 0.50, 0.95, 0.85, 0.055);
    TH1D *lv = (TH1D*) ppD->Clone("legSpecB"); styleH(lv, "#000000", markFilledCircle, 1.1);
    TH1D *lp = (TH1D*) ppD->Clone("legSpecBPP"); styleH(lp, "#999999", markOpenSquare, 1.0);
    TBox *lb = new TBox(); lb->SetFillColorAlpha(col, 0.30); lb->SetLineColor(col);
    lg->AddEntry(lv, "stat.", "pe"); lg->AddEntry(lb, "total syst.", "f"); lg->AddEntry(lp, "pp, for reference", "pe"); lg->Draw();
    savePdfTight(cv, Form("%sbJetsPerZ_PF_spectra_withSys%s.pdf", figDir.Data(), sfx.Data()));
    delete cv;
  }
  ratioFig(5,  4, 2, 2, "#it{b} jets per #it{N}_{Z}: PbPb / pp",     "bJetsPerZ_PF_PbPbOverPP_withSys");
  ratioFig(9,  3, 3, 1, "#it{b} jets per #it{N}_{Z}: class / 50-80%", "bJetsPerZ_PF_RCP_withSys");
  ratioFig(12, 4, 2, 2, "#it{b} / inclusive: PbPb / pp",             "bOverIncl_PF_PbPbOverPP_withSys");
  ratioFig(16, 3, 3, 1, "#it{b} / inclusive: class / 50-80%",        "bOverIncl_PF_RCP_withSys");

  // the b-jet fraction per system: a grid with bands, and one panel with every
  // system overlaid (stat. only, one color and marker per system)
  {
    TCanvas *cv = new TCanvas("cvFrac", "", 520*3, 460*2);
    cv->Divide(3, 2, 0.001, 0.001);
    const int col = TColor::GetColor("#D55E00");
    double yHi = 0.;
    for(int q = 19; q < 24; q++){
      TH1D *v = (TH1D*) fo->Get(Form("h_nom_%s", qTag(q).Data())), *to = (TH1D*) fo->Get(Form("h_sysRel_total_%s", qTag(q).Data()));
      for(int i = 1; i <= NB; i++) if(reported(i) && edges[i-1] < ptDrawMax - 1e-6 && v->GetBinContent(i) > 0.) yHi = TMath::Max(yHi, v->GetBinContent(i)*(1. + to->GetBinContent(i)) + v->GetBinError(i));
    }
    yHi *= 135.;   // in percent, with headroom
    for(int q = 19; q < 24; q++){
      cv->cd(q - 18);
      gPad->SetLeftMargin(0.14); gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.08); gPad->SetBottomMargin(0.13);
      TH1D *v  = (TH1D*) fo->Get(Form("h_nom_%s", qTag(q).Data()))->Clone(Form("frac_%d", q)); v->SetDirectory(nullptr); v->Scale(100.);
      TH1D *to = (TH1D*) fo->Get(Form("h_sysRel_total_%s", qTag(q).Data()));
      TH1F *fr = gPad->DrawFrame(ptReportMin, 0., ptDrawMax, yHi);
      fr->GetXaxis()->SetTitle("#it{p}_{T}^{jet} [GeV]"); fr->GetYaxis()->SetTitle("#it{b}-jet fraction [%]");
      fr->GetXaxis()->SetTitleSize(0.050); fr->GetYaxis()->SetTitleSize(0.050);
      fr->GetXaxis()->SetLabelSize(0.045); fr->GetYaxis()->SetLabelSize(0.045); fr->GetYaxis()->SetTitleOffset(1.3);
      for(int i = 1; i <= NB; i++){
        double y = v->GetBinContent(i); if(y <= 0. || !reported(i) || edges[i-1] >= ptDrawMax - 1e-6) continue;
        TBox *bx = new TBox(edges[i-1], y*(1. - to->GetBinContent(i)), edges[i], y*(1. + to->GetBinContent(i)));
        bx->SetFillColorAlpha(col, 0.30); bx->SetLineColor(col); bx->Draw();
      }
      for(int i = 1; i <= NB; i++) if(!reported(i) || v->GetBinContent(i) <= 0.){ v->SetBinContent(i, -999.); v->SetBinError(i, 0.); }
      styleH(v, "#000000", markFilledCircle, 1.1); v->Draw("E1 X0 same");
      TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.050);
      la.DrawLatex(0.17, 0.935, Form("%s, %s", sysLabel[q - 19], q == 19 || modePF ? "PF jets" : "calo-matched PF jets"));
      if(q == 19){
        TLegend *lg = makeLegend(0.17, 0.74, 0.70, 0.88, 0.042);
        lg->AddEntry(v, "stat.", "pe");
        TBox *lb = new TBox(); lb->SetFillColorAlpha(col, 0.30); lb->SetLineColor(col);
        lg->AddEntry(lb, "total syst.", "f"); lg->Draw();
      }
    }
    savePdfTight(cv, Form("%sbJetFraction_PF_withSys%s.pdf", figDir.Data(), sfx.Data()));
    delete cv;

    TCanvas *co = new TCanvas("cvFracOv", "", 700, 650);
    co->SetLeftMargin(0.15); co->SetRightMargin(0.04); co->SetTopMargin(0.05); co->SetBottomMargin(0.14);
    TH1F *fr = co->DrawFrame(ptReportMin, 0., ptDrawMax, yHi);
    fr->GetXaxis()->SetTitle("#it{p}_{T}^{jet} [GeV]"); fr->GetYaxis()->SetTitle("#it{b}-jet fraction [%]");
    fr->GetXaxis()->SetTitleSize(0.055); fr->GetYaxis()->SetTitleSize(0.055);
    fr->GetXaxis()->SetLabelSize(0.048); fr->GetYaxis()->SetLabelSize(0.048);
    fr->GetXaxis()->SetTitleOffset(1.1); fr->GetYaxis()->SetTitleOffset(1.25); fr->GetYaxis()->SetNdivisions(505);
    const int ovColor[5] = {7, 6, 5, 3, 0}, ovMark[5] = {markOpenSquare, markFilledCircle, markFilledSquare, markFilledDiamond, markCross};
    TLegend *lg = makeLegend(0.58, 0.55, 0.93, 0.90, 0.048);
    for(int q = 19; q < 24; q++){
      TH1D *v = (TH1D*) fo->Get(Form("h_nom_%s", qTag(q).Data()))->Clone(Form("fracOv_%d", q)); v->SetDirectory(nullptr); v->Scale(100.);
      for(int i = 1; i <= NB; i++) if(!reported(i) || v->GetBinContent(i) <= 0.){ v->SetBinContent(i, -999.); v->SetBinError(i, 0.); }
      styleH(v, okabeHex[ovColor[q-19]], ovMark[q-19], ovMark[q-19] == markFilledDiamond ? 2.3 : 1.6); v->SetLineWidth(2); v->Draw("E1 X0 same");
      lg->AddEntry(v, sysLabel[q-19], "pe");
    }
    lg->Draw();
    TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.050);
    la.DrawLatex(0.19, 0.88, modePF ? "PF jets, stat. only" : "PF jets (PbPb calo-matched), stat. only");
    savePdfTight(co, Form("%sbJetFraction_PF_overlay%s.pdf", figDir.Data(), sfx.Data()));
    delete co;
  }

  fo->Close();
  printf("\nwrote %s%s%s.root and figures in %s, %s\n", outDir, outName.Data(), sfx.Data(), figSysDir, figDir.Data());
}
