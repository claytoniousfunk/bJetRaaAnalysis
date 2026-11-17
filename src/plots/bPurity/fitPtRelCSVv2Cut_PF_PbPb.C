// PbPb PF b purity from a ptRel fit AFTER a CSVv2 cut, in the nominal centrality classes (2026-10-26,
// user: "run the CSVv2+pTrel template fit at the medium working point and derive the b-purity in the
// nominal centrality classes"). PbPb counterpart of fitPtRelCSVv2Cut_pp_PF.C (included: same
// likelihood, Barlow-Beeston lite, tagged-jet window), at one working point (default medium, 0.84).
//
// Per class (0-10 / 10-30 / 30-50 / 50-80% = ultraFine slices 1-2 / 3-6 / 7-10 / 11-16) and window:
//   P_tag    b fraction of the TAGGED jets: ptRel fit (0-fitMax GeV) of the jets above the cut, b vs
//            c + light at PYTHIA+HYDJET's c share among the tagged jets
//   eps_b    PYTHIA+HYDJET b-tag efficiency of muon-tagged b jets (tagged / all, ptRel 0-5 GeV)
//   f_pass   tagged / all muon-tagged jets in data
//   P_all^corr = P_tag f_pass / eps_b   the b purity of ALL muon-tagged jets the cut implies
//            (data / PYTHIA+HYDJET b-tag efficiency taken as 1), i.e. what the b-jet chain would use
//   compared with the nominal chain purity (rootFiles/systematics/bPuritySys_PFJets_PbPb.root) and with
//   a ptRel-only fit of all jets with the SAME templates as here (reference, no fake component).
//
// DATA: the 2026-10-26 PbPb SingleMuon scan (plain PF, mu12 with the trigger-efficiency correction),
//   h_muptrel_recoJetPt_csvV2_inclRecoMuonTag_triggerOn_C<i>, summed over the class's slices.
// TEMPLATES: the 2026-10-26 PYTHIA+HYDJET pThat-30 DiJet, BJet and MuJet scans (CSVv2 maps per flavor
//   and slice), combined as the nominal PbPb fit does (caloJetBPurityFit_PbPb.h):
//   SHAPES    (ptRel x pass/fail together) b flavor creation and g->bb from DiJet + BJet + MuJet, c from
//             DiJet + MuJet, light u+d+s+g from DiJet + BJet + MuJet: each sample normalised to unit
//             area, averaged with weights = its effective entries in the window
//   FRACTIONS DiJet only: the natural g->bb share, the b fraction, the c share of the non-b jets
//   Each combined shape is scaled to DiJet's yield of that flavor, so tag efficiencies come from the
//   combined shapes and flavor fractions from DiJet.
// DIFFERENCES FROM THE NOMINAL PbPb FIT, all deliberate: light from PYTHIA+HYDJET (not smeared pp; after
//   the cut light is a few per mille of the tagged jets), no fake-ptRel fit component; instead, since
//   2026-11-09 (user), the uncorrelated muon-jet pairs (forest-CS mixed events) are SUBTRACTED from the data
//   per triggered event in buildCls (applyFakeSub below; fakeSubPerEvt = false switches it off), fit 0-fitMax GeV
//   (default 5, the CSVv2 convention; the nominal edge is the light+c 0.98 quantile).
// CLOSURE: the DiJet all-flavor map itself as pseudo-data (fluctuating, not Asimov; DiJet is also in
//   the templates), against DiJet's true tagged and all-jet b fractions.
//
// Output: rootFiles/bPurityResults/ptRelCSVv2Cut_PF_PbPb_wp<cut>_fit<max>.root
//         figures/bPurity/PFJets_PbPb/csvV2/{purity,fits_<cls>}_wp<cut>_fit<max>.pdf
// Usage, from src/plots/bPurity/:  root -l -b -q 'fitPtRelCSVv2Cut_PF_PbPb.C(0.84, 5.0)'
// include guard (2026-11-16): lets a study include both fitPtRelCSVv2Cut_lowPt_PF.C and fitPtRelCSVv2_2D_ppPbPb_PF.C
#ifndef FITPTRELCSVV2CUT_PF_PBPB_C
#define FITPTRELCSVV2CUT_PF_PBPB_C
#include "fitPtRelCSVv2Cut_pp_PF.C"

namespace pbc {
using namespace p2d;
const char *dataPbPb = "rootFiles/scanningOuput/PbPb/PbPb_SingleMuon_manualJEC_mu12_pTmu-15to999_tight_mu12TriggerEfficiencyCorrection_jetTrkMaxFilter_WDecayFilter_2026-10-26_ultraFineCentBins.root";
const char *phFmt = "rootFiles/scanningOuput/PYTHIAHYDJET/PYTHIAHYDJET_manualJEC_%s_pThat-30_mu12_pTmu-15to999_tight_vzReweight_hiBinShift-10_jetTrkMaxFilter_doPThatCorrelationFilterTight_removeUnmatchedJets_WDecayFilter_genMuTagFix_2026-10-26_ultraFineCentBins.root";
const char *nomPurPath = "rootFiles/systematics/bPuritySys_PFJets_PbPb.root";
const int NC = 4;
const int slc[NC][2] = {{1, 2}, {3, 6}, {7, 10}, {11, 16}};
const char *cName[NC] = {"0-10%", "10-30%", "30-50%", "50-80%"}, *cTag[NC] = {"C1", "C2", "C3", "C4"};
const int NS = 3; const char *sName[NS] = {"DiJet", "BJet", "MuJet"};
// flavor groups: b flavor creation, b from g->bb, c, light
const int NG = 4;
const std::vector<std::vector<const char*>> gMembers = {{"bJets"}, {"bGSJets"}, {"cJets"}, {"uJets", "dJets", "sJets", "gJets"}};
const std::vector<std::vector<int>> gSamples = {{0, 1, 2}, {0, 1, 2}, {0, 2}, {0, 1, 2}};
TFile *fS[NS] = {};

struct Arr { double n[KMAX][NPR] = {}, e2[KMAX][NPR] = {}; double tot() const { double s = 0; for(int k = 0; k < K; k++) for(int j = 0; j < nPR; j++) s += n[k][j]; return s; }
             double neff() const { double s = 0, v = 0; for(int k = 0; k < K; k++) for(int j = 0; j < nPR; j++){ s += n[k][j]; v += e2[k][j]; } return v > 0. ? s*s/v : 0.; } };
std::vector<TString> names(const char *base, const std::vector<const char*> &fls, int c){
  std::vector<TString> v; for(const char *fl : fls) for(int i = slc[c][0]; i <= slc[c][1]; i++) v.push_back(fl[0] ? Form("%s_%s_C%d", base, fl, i) : Form("%s_C%d", base, i)); return v;
}
bool fillArr(TFile *f, const std::vector<TString> &nm, double lo, double hi, Arr &a){ return fillCounts(f, nm, lo, hi, a.n, a.e2); }
// combined unit-area shape of flavor group G, samples weighted by their effective entries
bool combined(int G, int c, double lo, double hi, Arr &out){
  out = Arr(); double wsum = 0.; Arr a[NS]; double ne[NS] = {}, tt[NS] = {};
  for(int s : gSamples[G]){ if(!fillArr(fS[s], names(h3, gMembers[G], c), lo, hi, a[s])) return false; tt[s] = a[s].tot(); ne[s] = tt[s] > 0. ? a[s].neff() : 0.; wsum += ne[s]; }
  if(wsum <= 0.) return true;
  for(int s : gSamples[G]){ if(tt[s] <= 0.) continue; const double w = ne[s]/wsum;
    for(int k = 0; k < K; k++) for(int j = 0; j < nPR; j++){ out.n[k][j] += w*a[s].n[k][j]/tt[s]; out.e2[k][j] += w*w*a[s].e2[k][j]/(tt[s]*tt[s]); } }
  return true;
}
// UNCORRELATED MUON-JET PAIRS, subtracted PER EVENT from the data (NOMINAL since 2026-11-09, user:
// "Change the nominal flow to use this correction"; study: studyFakeMuJetCS_CSVv2_PF_PbPb.C).
// fake ptRel = h_csJetMuonPtRel_csJetPt_C<s> of the forest-CS mixed-event SingleMuon pfCand scan (a
// mixed-event PF muon > 15 GeV inside an akCs4PF-emulated jet of the mixed event), per triggered event of
// that scan; data per triggered event of the data scan; both trigger-gated, so
//   data_k(ptRel) -= N_trig(data) / N_trig(mixed) x fake(ptRel) x (share of the window's data in CSVv2
//   category k)   [the emulated jets carry no CSVv2: the fakes pass the cut as the data do, eps_fake = f_pass]
// bin by bin, floored at 0. X.nData stays the UNSUBTRACTED total, so every consumer's
// P = P_tag (P.nData / X.nData) / eps_b is b / ALL muon-tagged jets, which is what the chain multiplies.
// fakeSubPerEvt = false: the fit before (the bPurity source "fakeMuJet" = |with / without - 1|).
bool fakeSubPerEvt = true;
const char *fakeMixPath = "rootFiles/scanningOuput/PbPb/PbPb_SingleMuon_mu12_pTmu-15to999_tight_mu12TriggerEfficiencyCorrection_jetTrkMaxFilter_WDecayFilter_mixedEventPFClustering_fastJetResamples-100_pseudoJetCandPtMin-0.0_caloConstituents_bkgMapOverride_2026-11-6_forestCS_ultraFineCentBins.root";
double lastFakeFrac = 0.;   // fake / all muon-tagged jets (0-5 GeV) removed by the last applyFakeSub
// subtract the fake ptRel from X.d (class c, window [lo, hi)); returns the unsubtracted total, or -1 on error
double applyFakeSub(TFile *fD, int c, double lo, double hi, Win &X){
  double raw = 0., cat[KMAX] = {};
  for(int k = 0; k < K; k++) for(int j = 0; j < nPR; j++){ raw += X.d[k][j]; cat[k] += X.d[k][j]; }
  lastFakeFrac = 0.;
  if(!fakeSubPerEvt || raw <= 0.) return raw;
  static TFile *fF = nullptr;
  if(!fF){ fF = TFile::Open(Form("%s/%s", repo, fakeMixPath)); if(!fF || fF->IsZombie()){ printf("ERROR: cannot open %s\n", fakeMixPath); fF = nullptr; return -1.; } }
  double nTD = 0., nTM = 0., f[NPR] = {};
  for(int s = slc[c][0]; s <= slc[c][1]; s++){
    TH1 *a = (TH1*) fD->Get(Form("h_vz_triggerOn_C%d", s)), *b = (TH1*) fF->Get(Form("h_vz_triggerOn_C%d", s));
    TH2 *h = (TH2*) fF->Get(Form("h_csJetMuonPtRel_csJetPt_C%d", s));
    if(!a || !b || !h){ printf("ERROR: fake subtraction inputs missing for slice %d\n", s); return -1.; }
    nTD += a->Integral(0, a->GetNbinsX() + 1); nTM += b->Integral(0, b->GetNbinsX() + 1);
    const int y1 = h->GetYaxis()->FindBin(lo + 1e-6), y2 = h->GetYaxis()->FindBin(hi - 1e-6);   // the template's own jet-pT axis
    for(int ix = 1; ix <= h->GetNbinsX(); ix++){
      const int j = (int)(h->GetXaxis()->GetBinCenter(ix)/ptRelBW); if(j >= nPR) continue;
      for(int iy = y1; iy <= y2; iy++) f[j] += h->GetBinContent(ix, iy);
    }
  }
  if(nTM <= 0.) return -1.;
  double removed = 0.;
  for(int k = 0; k < K; k++) for(int j = 0; j < nPR; j++){
    const double sub = f[j]*nTD/nTM*cat[k]/raw, v = TMath::Max(0., X.d[k][j] - sub);
    removed += X.d[k][j] - v; X.d[k][j] = v;
  }
  lastFakeFrac = removed/raw;
  return raw;
}
// one window of one class: data, combined templates scaled to DiJet yields; nEffB = combined b statistics
bool buildCls(TFile *fD, int c, int w, Win &X, double &nEffB){
  X = Win(); X.sigma = 0.;
  double dE2[KMAX][NPR] = {};
  if(!fillCounts(fD, names(h3, {""}, c), wLo[w], wHi[w], X.d, dE2)) return false;
  const double nRaw = applyFakeSub(fD, c, wLo[w], wHi[w], X); if(nRaw < 0.) return false;
  Arr dj[NG], sh[NG];
  for(int G = 0; G < NG; G++){ if(!fillArr(fS[0], names(h3, gMembers[G], c), wLo[w], wHi[w], dj[G])) return false; if(!combined(G, c, wLo[w], wHi[w], sh[G])) return false; }
  const double yFC = dj[0].tot(), yGS = dj[1].tot(), yB = yFC + yGS, gs = yB > 0. ? gsShare(yGS/yB) : 0.;   // natural share x gsMult (gsShare)
  const double yF[NF] = {yB, dj[2].tot(), dj[3].tot()};
  for(int k = 0; k < K; k++) for(int j = 0; j < nPR; j++){
    // b = flavor creation + g->bb at gsShare(DiJet's natural share), then every flavor scaled to DiJet's yield
    const double bn = (1. - gs)*sh[0].n[k][j] + gs*sh[1].n[k][j], be = (1. - gs)*(1. - gs)*sh[0].e2[k][j] + gs*gs*sh[1].e2[k][j];
    X.t[0].n[k][j] = yF[0]*bn; X.t[0].e2[k][j] = yF[0]*yF[0]*be;
    X.t[1].n[k][j] = yF[1]*sh[2].n[k][j]; X.t[1].e2[k][j] = yF[1]*yF[1]*sh[2].e2[k][j];
    X.t[2].n[k][j] = yF[2]*sh[3].n[k][j]; X.t[2].e2[k][j] = yF[2]*yF[2]*sh[3].e2[k][j];
  }
  for(int F = 0; F < NF; F++){ for(int k = 0; k < K; k++) for(int j = 0; j < nPR; j++){ if(X.t[F].n[k][j] < 0.) X.t[F].n[k][j] = 0.; X.t[F].tot += X.t[F].n[k][j]; }
    for(int k = 0; k < K; k++){ double s = 0; for(int j = 0; j < nPR; j++) s += X.t[F].n[k][j]; X.t[F].eps[k] = X.t[F].tot > 0. ? s/X.t[F].tot : 0.; } }
  X.truthB = X.t[0].tot/(X.t[0].tot + X.t[1].tot + X.t[2].tot);
  X.rc = X.t[1].tot/(X.t[1].tot + X.t[2].tot);
  X.nData = nRaw;   // the UNSUBTRACTED total (see applyFakeSub)
  // statistics of the combined b template (FC and GS samples together)
  nEffB = 0.; for(int s = 0; s < NS; s++){ Arr a; fillArr(fS[s], names(h3, {"bJets", "bGSJets"}, c), wLo[w], wHi[w], a); nEffB += a.neff(); }
  return true;
}
} // namespace pbc

void fitPtRelCSVv2Cut_PF_PbPb(double wp = 0.84, double fitMax = 5.0){
  using namespace pbc; using namespace pcut;
  nFit = (int)(fitMax/ptRelBW + 0.5);
  if(nFit < 2 || nFit > nPR){ printf("ERROR: fitMax in 1-5 GeV\n"); return; }
  setCats(Form("cut%.2f", wp));
  setJetCollection("PFJets"); initPlotStyle();
  TFile *fD = TFile::Open(Form("%s/%s", repo, dataPbPb)), *fN = TFile::Open(Form("%s/%s", repo, nomPurPath));
  for(int s = 0; s < NS; s++) fS[s] = TFile::Open(Form("%s/%s", repo, Form(phFmt, sName[s])));
  if(!fD || !fN || !fS[0] || !fS[1] || !fS[2]){ printf("ERROR: inputs\n"); return; }
  const TString wpTag = TString::Format("wp%s", TString::Format("%.2f", wp).ReplaceAll(".", "p").Data()), tag = wpTag + Form("_fit%g", fitMax) + inputTag + fitEngineTag();   // inputTag (p2d, default "") names a driver's variant, e.g. "_gsNat"
  const TString outDir = Form("%s/figures/bPurity/PFJets_PbPb/csvV2/", repo); gSystem->mkdir(outDir, true);
  const Model &M1 = models[0];

  TH1D *hTag[NC], *hTagMC[NC], *hCorr[NC], *hAll[NC], *hAllMC[NC], *hNom[NC], *hNomSys[NC], *hEps[NC], *hPass[NC];
  Win Pw[NC][NW]; Result Rt[NC][NW];
  for(int c = 0; c < NC; c++){
    auto bk = [&](const char *n){ TH1D *h = new TH1D(Form("%s_%s", n, cTag[c]), "", NW, ptEdges + 3); h->SetDirectory(nullptr); return h; };
    hTag[c] = bk("h_purTag"); hTagMC[c] = bk("h_purTagMC"); hCorr[c] = bk("h_purCorr"); hAll[c] = bk("h_purAllSameTpl"); hAllMC[c] = bk("h_purAllMC");
    hNom[c] = bk("h_purNominal"); hNomSys[c] = bk("h_purNominalSys"); hEps[c] = bk("h_epsB"); hPass[c] = bk("h_fPass");
    TH1 *nP = (TH1*) fN->Get(Form("h_bPurity_nominal_%s", cTag[c])), *nS = (TH1*) fN->Get(Form("h_sysAbs_total_%s", cTag[c]));
    printf("\n==== %s, CSVv2 > %.2f, fit 0-%.1f GeV\n", cName[c], wp, fitMax);
    printf("  window   | f_pass eps_b  eps_c  eps_l | b Neff | P_tag fit      P_tag P+H | P_all^corr     P+H truth | nominal (sys)    | P_all same tpl | implied SF_b vs nominal | dev/ndf tag\n");
    for(int w = 0; w < NW; w++){
      Win X; double nEffB = 0.;
      if(!buildCls(fD, c, w, X, nEffB)) return;
      Win P = passWin(X);
      P.truthB = X.t[0].tot > 0. ? X.t[0].eps[1]*X.t[0].tot/(X.t[0].eps[1]*X.t[0].tot + X.t[1].eps[1]*X.t[1].tot + X.t[2].eps[1]*X.t[2].tot) : 0.;
      Result rT = fitWin(M1, P, P.d), rA = fitWin(M1, X, X.d);
      const double fPass = P.nData/X.nData, eB = X.t[0].eps[1];
      const double corr = rT.pur*fPass/eB, ecorr = rT.epur*fPass/eB;
      const int nb = nP ? nP->FindBin(0.5*(wLo[w] + wHi[w])) : 0;
      const double pN = nP ? nP->GetBinContent(nb) : 0., sN = nS ? nS->GetBinContent(nb) : 0.;   // h_sysAbs_total: the absolute systematic is the CONTENT
      printf("  %3.0f-%3.0f  | %.3f  %.3f  %.3f  %.3f | %6.0f | %.3f +- %.3f  %.3f     | %.3f +- %.3f  %.3f     | %.3f +- %.3f    | %.3f +- %.3f  | %.3f                   | %.1f/%d\n",
             wLo[w], wHi[w], fPass, eB, X.t[1].eps[1], X.t[2].eps[1], nEffB, rT.pur, rT.epur, P.truthB, corr, ecorr, X.truthB, pN, sN, rA.pur, rA.epur, pN > 0. ? corr/pN : 0., rT.dev, rT.ndf);
      hTag[c]->SetBinContent(w + 1, rT.pur); hTag[c]->SetBinError(w + 1, rT.epur); hTagMC[c]->SetBinContent(w + 1, P.truthB);
      hCorr[c]->SetBinContent(w + 1, corr); hCorr[c]->SetBinError(w + 1, ecorr);
      hAll[c]->SetBinContent(w + 1, rA.pur); hAll[c]->SetBinError(w + 1, rA.epur); hAllMC[c]->SetBinContent(w + 1, X.truthB);
      hNom[c]->SetBinContent(w + 1, pN); hNomSys[c]->SetBinContent(w + 1, pN); hNomSys[c]->SetBinError(w + 1, sN);
      hEps[c]->SetBinContent(w + 1, eB); hPass[c]->SetBinContent(w + 1, fPass);
      Pw[c][w] = P; Rt[c][w] = rT;
    }
    // closure: DiJet all flavors as pseudo-data
    printf("  closure, DiJet as pseudo-data (fit - truth): tagged purity | corrected all-jet purity | DiJet effective jets, tagged / all\n");
    for(int w = 0; w < NW; w++){
      Win X; double nEffB = 0.; if(!buildCls(fD, c, w, X, nEffB)) return;
      double ps[KMAX][NPR] = {}, pe[KMAX][NPR] = {};
      fillCounts(fS[0], names(h3, {"allJets"}, c), wLo[w], wHi[w], ps, pe);
      double nAll = 0., nPass = 0., vAll = 0., vPass = 0.; for(int j = 0; j < nPR; j++){ nAll += ps[0][j] + ps[1][j]; nPass += ps[1][j]; vAll += pe[0][j] + pe[1][j]; vPass += pe[1][j]; }
      for(int k = 0; k < K; k++) for(int j = 0; j < nPR; j++) X.d[k][j] = ps[k][j]*X.nData/nAll;   // the data's size, DiJet's shape
      Win P = passWin(X);
      const double trTag = X.t[0].eps[1]*X.t[0].tot/(X.t[0].eps[1]*X.t[0].tot + X.t[1].eps[1]*X.t[1].tot + X.t[2].eps[1]*X.t[2].tot);
      Result rT = fitWin(M1, P, P.d);
      const double corr = rT.pur*(nPass/nAll)/X.t[0].eps[1];
      // pull: the fit error rescaled from the data's size to DiJet's effective size (DiJet's own statistics)
      const double nEffPass = vPass > 0. ? nPass*nPass/vPass : 0., eDJ = nEffPass > 0. ? rT.epur*sqrt(P.nData/nEffPass) : 0.;
      printf("    %3.0f-%3.0f  %+.3f (truth %.3f, pull %+.1f) | %+.3f (truth %.3f) | %.0f / %.0f\n", wLo[w], wHi[w], rT.pur - trTag, trTag, eDJ > 0. ? (rT.pur - trTag)/eDJ : 0., corr - X.truthB, X.truthB, nEffPass, vAll > 0. ? nAll*nAll/vAll : 0.);
    }
  }
  {
    TFile *fo = TFile::Open(Form("%s/rootFiles/bPurityResults/ptRelCSVv2Cut_PF_PbPb_%s.root", repo, tag.Data()), "RECREATE");
    for(int c = 0; c < NC; c++) for(TH1D *h : {hTag[c], hTagMC[c], hCorr[c], hAll[c], hAllMC[c], hNom[c], hNomSys[c], hEps[c], hPass[c]}) h->Write();
    TNamed info("info", Form("PbPb PF b purity from a ptRel fit after CSVv2 > %.2f, fit 0-%.1f GeV, per class (C1..C4 = 0-10, 10-30, 30-50, 50-80%%). h_purCorr = P_tag f_pass / eps_b(P+H). data: %s. templates: %s (DiJet, BJet, MuJet)", wp, fitMax, dataPbPb, phFmt));
    info.Write(); fo->Close();
  }

  // ---- figure 1: all-jet purity per class, peripheral (left) to central (right)
  {
    TCanvas *cv = new TCanvas("cvPur", "", 2000, 700);
    GridPads g = makeGridPads(cv, NC, 1, 0.06, 0.01, 0.25, 0.10);
    TLegend *lgT = nullptr;
    for(int i = 0; i < NC; i++){
      const int c = NC - 1 - i;
      g.at(i, 0)->cd();
      TH1F *fr = gPad->DrawFrame(80., 0., 300., 1.15); g.styleFrame(fr, i, 0, 0.055);
      hNomSys[c]->SetFillColorAlpha(kGray + 1, 0.45); hNomSys[c]->SetLineWidth(0); hNomSys[c]->SetMarkerSize(0); hNomSys[c]->Draw("E2 SAME");
      styleLine(hAllMC[c], "#8C8C8C", 4); hAllMC[c]->SetLineStyle(7); hAllMC[c]->Draw("HIST SAME");
      styleH(hNom[c], hexData, markFilledCircle, 2.4); hNom[c]->SetLineWidth(3); hNom[c]->Draw("P X0 SAME");
      styleH(hCorr[c], okabeHex[1], markFilledDiamond, 2.8); hCorr[c]->SetLineWidth(3); hCorr[c]->Draw("E1 X0 SAME");
      g.text(i, 0, 0.09, cName[c], 0.065, 62);
      if(i == 0){   // one legend row in the top margin; the tagged-jet purity is on the fits_<cls> figures
        cv->cd(); lgT = makeLegend(0.08, 0.005, 0.98, 0.075, 0.045); lgT->SetNColumns(3); TLegend *lg = lgT;
        lg->AddEntry(hNom[c], "nominal #it{p}_{T}^{rel} fit (syst.)", "p");
        lg->AddEntry(hCorr[c], Form("CSVv2 > %.2f: #it{P}_{tag} #it{f}_{pass} / #varepsilon_{b}^{P+H}", wp), "pe");
        lg->AddEntry(hAllMC[c], "PYTHIA+HYDJET b fraction", "l");
      }
    }
    cv->cd(); lgT->Draw();
    g.titles(cv, "#it{p}_{T}^{jet} (GeV)", "b purity", 0.050, 0.020, 0.10);
    g.cmsLabel(cv, "Internal", Form("PbPb (5.02 TeV),  CSVv2 > %.2f,  fit #it{p}_{T}^{rel} < %.0f GeV", wp, fitMax), 1.2);
    savePdfTight(cv, outDir + "purity_" + tag + ".pdf");
  }
  // ---- figure 2: the tagged-jet ptRel fits, one figure per class, every window; fit / data below
  const int colF[NF] = {kRed - 4, kGreen + 2, kBlue - 4};   // b, c, light (user, 2026-10-22)
  for(int c = 0; c < NC; c++){
    TCanvas *cv = new TCanvas(Form("cvF%d", c), "", 2200, 950);
    const double lm = 0.06, rm = 0.02, bm = 0.13, tm = 0.09;
    GridPads g = makeGridPads(cv, NW, 2, lm, rm, bm, tm);
    const double inner = 1. - bm - tm, yS = bm + 0.30*inner;
    for(int w = 0; w < NW; w++){
      cur = &M1; cw = &Pw[c][w]; cd = Pw[c][w].d; Kf = 1;
      double mu[KMAX][NPR], va[KMAX][NPR], comp[NF][KMAX][NPR];
      const int nSave = nFit; expectAll(Rt[c][w].par, mu, va, comp);
      double muIn = 0.; for(int j = 0; j < nSave; j++) muIn += mu[0][j];
      nFit = nPR; expectAll(Rt[c][w].par, mu, va, comp);
      double muIn2 = 0.; for(int j = 0; j < nSave; j++) muIn2 += mu[0][j];
      nFit = nSave; const double rs = muIn2 > 0. ? muIn/muIn2 : 1.;
      TH1D *hd = new TH1D(Form("pd%d_%d", c, w), "", nPR, 0., 5.), *ht = new TH1D(Form("pt%d_%d", c, w), "", nPR, 0., 5.), *hc[NF];
      hd->SetDirectory(nullptr); ht->SetDirectory(nullptr);
      for(int F = 0; F < NF; F++){ hc[F] = new TH1D(Form("pc%d_%d_%d", c, w, F), "", nPR, 0., 5.); hc[F]->SetDirectory(nullptr); }
      double ymax = 0.;
      for(int j = 0; j < nPR; j++){
        const double d = Pw[c][w].d[0][j];
        hd->SetBinContent(j + 1, d); hd->SetBinError(j + 1, sqrt(d));
        for(int F = 0; F < NF; F++) hc[F]->SetBinContent(j + 1, comp[F][0][j]*rs);
        ht->SetBinContent(j + 1, mu[0][j]*rs); ymax = TMath::Max(ymax, TMath::Max(d, mu[0][j]*rs));
      }
      TPad *pt = g.at(w, 0), *pb = g.at(w, 1);
      const double x1 = pt->GetXlowNDC(), x2 = x1 + pt->GetWNDC();
      pt->SetPad(x1, yS, x2, 1.); pt->SetTopMargin(tm/(1. - yS)); pt->SetBottomMargin(0.);
      pb->SetPad(x1, 0., x2, yS); pb->SetBottomMargin(bm/yS); pb->SetTopMargin(0.);
      pt->cd();
      TH1F *fr = gPad->DrawFrame(0., 0., 5., 1.35*ymax); g.styleFrame(fr, w, 0, 0.05);
      fr->GetXaxis()->SetLabelSize(0.); fr->GetYaxis()->SetLabelSize(0.050); fr->GetYaxis()->SetMaxDigits(3);
      THStack *st = new THStack(Form("pst%d_%d", c, w), "");
      for(int F = NF - 1; F >= 0; F--){ hc[F]->SetFillColor(colF[F]); hc[F]->SetLineColor(colF[F]); st->Add(hc[F]); }
      st->Draw("HIST SAME");
      styleH(hd, hexData, markFilledCircle, 1.8); hd->SetLineWidth(3); hd->Draw("E1 X0 SAME");
      if(fitMax < 5.){ TLine *e = new TLine(fitMax, 0., fitMax, 1.35*ymax); e->SetLineStyle(7); e->SetLineWidth(3); e->SetLineColor(kGray + 2); e->Draw(); }
      gPad->RedrawAxis();
      g.text(w, 0, 0.10, Form("%.0f-%.0f GeV", wLo[w], wHi[w]), 0.060);
      g.text(w, 0, 0.19, Form("tagged purity %.2f", Rt[c][w].pur), 0.055);
      if(w == NW - 1){ TLegend *lg = makeLegend(0.50, 0.36, 0.97, 0.72, 0.060); lg->AddEntry(hd, "data", "pe"); for(int F = 0; F < NF; F++) lg->AddEntry(hc[F], Form("%s jets", flLabel[F]), "f"); lg->Draw(); }
      pb->cd();
      TH1F *fr2 = gPad->DrawFrame(0., 0.4, 5., 1.6); g.styleFrame(fr2, w, 1, 0.12);
      fr2->GetYaxis()->SetNdivisions(203); fr2->GetXaxis()->SetTickLength(0.07);
      fr2->GetYaxis()->SetLabelSize(w == 0 ? 0.12 : 0.); fr2->GetXaxis()->SetLabelSize(0.12);
      TLine *u = new TLine(0., 1., 5., 1.); u->SetLineStyle(2); u->SetLineWidth(2); u->Draw();
      TH1D *r = (TH1D*) ht->Clone(Form("pr%d_%d", c, w)); r->SetDirectory(nullptr);
      for(int j = 1; j <= nPR; j++){ const double d = hd->GetBinContent(j);
        if(d > 0.){ r->SetBinContent(j, ht->GetBinContent(j)/d); r->SetBinError(j, ht->GetBinContent(j)/d*hd->GetBinError(j)/d); } else { r->SetBinContent(j, -999.); r->SetBinError(j, 0.); } }
      styleH(r, hexData, markFilledCircle, 1.6); r->SetLineWidth(3); r->Draw("E1 X0 SAME");
    }
    cv->cd(); TLatex a; a.SetNDC(); a.SetTextFont(42); a.SetTextSize(0.045);
    a.SetTextAlign(31); a.DrawLatex(1. - rm, 0.012, "muon #it{p}_{T}^{rel} (GeV)");
    a.SetTextAngle(90); a.DrawLatex(0.016, 1. - tm, Form("jets per %g GeV", ptRelBW)); a.DrawLatex(0.016, yS - 0.01, "fit / data");
    g.cmsLabel(cv, "Internal", Form("PbPb (5.02 TeV), %s,  CSVv2 > %.2f", cName[c], wp), 1.3);
    savePdfTight(cv, outDir + Form("fits_%s_", cTag[c]) + tag + ".pdf");
  }
  printf("\nfigures in %s\n", outDir.Data());
}
#endif // FITPTRELCSVV2CUT_PF_PBPB_C
