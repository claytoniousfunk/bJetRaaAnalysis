// PbPb CALO-jet b-jet chain per coarse centrality class, end to end, on the
// 2026-10-04 calo SingleMuon scan -- the PbPb counterpart of the pp calo chain
// (constructBJetSpectrum / unfoldBJetSpectrum / correctMuonRecoEff /
// correctMuonTagFrequency _caloJets_pp.C), kept in one macro so the four
// classes are treated identically.
//
// Per class:
//   1. b purity: 2-template RooFit (b vs light+c, c at its MC share) of the
//      tagging-muon ptRel on 0-5 GeV in 0.1 GeV bins, per jet-pT window, with
//      PYTHIA+HYDJET calo templates (PF-matched flavor) summed over the 5%
//      slices of the class. The calo PbPb scan has no separate g->bb label
//      (bGSJets is empty), so the b template is the merged b and there is no
//      bGS enhancement, unlike the pp nominal (+0.175). TEMPLATE STATISTICS:
//      the calo PYTHIA+HYDJET scan is pThat > 15 with no weight cut, and below
//      120 GeV its light+c templates are dominated by a few high-weight
//      low-pThat entries (effective entries 16 / 3 in 0-10% at 80-100 /
//      100-120 GeV). So (i) the two central classes share templates summed
//      over 0-30% and the two peripheral ones over 30-80%, and (ii) for the
//      80-100 and 100-120 GeV windows the light+c template SHAPE is taken
//      from 120-150 GeV of the same group. Both are stop-gaps until a
//      weight-cut or higher-pThat calo template scan exists.
//      No fake-ptRel (T2+T3)
//      subtraction: the tower-level mixed-event muon templates for calo jets
//      are not scanned yet. Both are differences from the pp chain, to be
//      closed.
//   2. b spectrum = muon-tagged, muon-triggered spectrum x purity (counts) on
//      ptEdges; 50-80 GeV are unfolding buffer bins (purity from the purity file's buffer values since 2026-11-10;
//      BJET_PBPB_BUFFER_PURITY=first = the 80-100 GeV purity, as before).
//   3. unfolding with the class's muon-tagged b response (odd-event half of
//      the 2026-10-02 calo response scan), reco floor 70 GeV, truth from the
//      untruncated matrix, Bayes N per class (nIter below).
//   4. / muon reco+ID efficiency (simulation per class, as the inclusive N_Z)
//      / muon-tag frequency f(gen pT) = gen-muon-tagged / all gen b jets from
//      the response scan (|flavor| 5 and 17), per class.
//   5. per Z with the inclusive chain's N_Z (60-120 GeV, same efficiencies),
//      so the b / inclusive ratio has one N_Z on both sides.
// pp: the pp calo chain's final spectrum (h_bJetPt_inclusive) over N_Z(pp).
// Inclusive: h_nom_<sys>_bjet of inclusiveJetsPerZSys_caloJets.C run at the
// tilt:0.3 closure optimum (pp 7, 0-10% 2, 10-30% 3, 30-50% 2, 50-80% 4), already
// on ptEdges, counts per Z. N_Z itself still comes from calculateJetsPerZ.cc's
// file; that macro reproduces it exactly.
//
// ITERATIONS  minimum MSE of
//   unfoldClosureTest.C("C<n>",15,"weighted","tiltScan:0.3",80.,500.,true,1,
//     <respPath without _ultraFineCentBins>, <same>, "_PbPb_bMuTagChain", -1,
//     "h_matchedRecoJetPt_genJetPt_bJets_C<n>_muTagged",
//     "h_matchedRecoJetPt_genJetPt_bJets_muTagged_C<n>", true, 70.,
//     "50,60,70,80,100,120,150,200,300,500")
// run 2026-10-05: N = 1 in all four classes. No even-half muon-tagged C1-C4
// response exists, so the tilted pseudo-data and the response come from the
// same (odd) half: the prior bias is tested properly, the variance term is not
// independent. The minimum is shallow: N = 2 is within 1-2% in MSE in every
// class and has half the bias (max|r-1| 2.4-4.1% vs 7.2-8.1% at N = 1). At
// tiltScan:0.5 the minimum is N = 2 in all four classes. pp (unfoldBJetSpectrum_
// caloJets_pp.C) is N = 2 at both tilts.
//
// Uncertainties: statistical only (counts (+) purity fit, through the
// unfolding covariance). Systematics follow the pp machinery, not yet run.
//
// Outputs
//   rootFiles/BJetPerZ_caloJets/bJetChain_caloJets_PbPb.root
//     bPurity_C*, bJetPt_meas_C*, bJetPt_unf_C*, bJetsPerZ_{pp,C*},
//     jetsPerZ_{pp,C*}, r_b_PbPbOverPP_C*, r_b_CP_C*, r_bOverIncl_PbPbOverPP_C*,
//     r_bOverIncl_CP_C*, muTagFreq_C*
//   figures/bJetsPerZ_caloJets/
//     ptRelFit_<class>_<lo>-<hi>.pdf, bPurity_vsJetPt_caloJets.pdf,
//     bJetsPerZ_PbPbOverPP_caloJets.pdf, bJetsPerZ_CP_caloJets.pdf,
//     bOverInclusive_PbPbOverPP_caloJets.pdf, bOverInclusive_CP_caloJets.pdf,
//     bJetUnfold_<class>.pdf
//
// Usage, from src/calculateBJetsPerZ/:
//   root -l -b -q -e 'gSystem->Load("/home/clayton/Programs/RooUnfold/build/libRooUnfold.so");
//                     gInterpreter->AddIncludePath("/home/clayton/Programs/RooUnfold/build");' bJetChain_caloJets_PbPb.C
#include "../../headers/functions/jetCollectionInputs.h"
#include "../../headers/functions/purityVariant.h"   // purVar(): purity-variant file suffix   // setJetCollection, dataPathPbPb, phRespPath, phRespMuTagFmt
#include "../../headers/plotting/coarseCent.h"
#include "../../headers/plotting/plotStyle.h"
#include "../../headers/plotting/ptRelFitPanel.h"
#include "../../headers/functions/divideByBinwidth.h"
#include "../../headers/functions/dataDrivenPrior.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooMsgService.h"
#include "TVectorD.h"
#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#endif
#include "../../headers/functions/bJetUnfoldExt.h"         // PFJets: truth axis from 20 GeV
#include "../../headers/functions/clipResponseOutliers.h"   // PFJets: single-event weight outliers
#include "../../headers/functions/jerSmear.h"   // JER scheme (PFJets): nominal response x1.10
#include "../../headers/publishedData/CMS_ZYield_5p02TeV.h"   // BJET_NZ_SOURCE=published: N_Z from the published yields

namespace bcal {

// PbPb data and response follow the jet collection: dataPathPbPb, phRespPath,
// phRespMuTagFmt (jetCollectionInputs.h). tplPath is the calo scan of the
// "chain" purity mode only. The paths below are rebuilt per collection in the
// macro (collTag() in the names); the literals are the caloJets defaults.
const char *dataPathCaloLegacy = "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/PbPb_SingleMuon_caloJets_manualJEC_mu12_pTmu-15to999_tight_mu12TriggerEfficiencyCorrection_jetTrkMaxFilter_WDecayFilter_2026-10-4_ultraFineCentBins.root";
const char *tplPath  = "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PYTHIAHYDJET/PYTHIAHYDJET_manualJEC_DiJet_caloJets_pThat-15_mu12_pTmu-15to999_tight_vzReweight_hiBinShift-10_jetTrkMaxFilter_doPThatCorrelationFilterTight_WDecayFilter_2026-10-2_ultraFineCentBins.root";
const char *respPathCaloLegacy = "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PYTHIAHYDJET/PYTHIAHYDJET_response_DiJet_caloJets_manualJEC_muTaggedJets_oddEvents_pThat-15_mu12_pTmu-15_tight_vzReweight_hiBinReweight_hiBinShift-10_jetTrkMaxFilter_doPThatCorrelationFilterTight_2026-10-2_ultraFineCentBins.root";
const char *inclPath = "/home/clayton/Analysis/code/bJetRaaAnalysis/src/calculateJetsPerZ/rootFiles/JetsPerZ/jetsPerZSpectra_caloJets_coarseMuEff_unfold1iter.root";   // N_Z only
const char *inclSpecPath = "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/systematics/inclJetsPerZSys_caloJets_JEUppUncorr_tilt0p3.root";
const char *ppBPath  = "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/CorrectedBJetSpectra/Data/bJetSpectrum_inclusive_caloJets_pp.root";
const char *purSysPath = "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/systematics/bPuritySys_caloJets_PbPb.root";   // purityMode "dedicated"
TString outRoot = "../../rootFiles/BJetPerZ_caloJets/bJetChain_caloJets_PbPb.root";   // + tag, set in the macro
TString figDir  = "../../figures/bJetsPerZ_caloJets/";

const int NE = 10; const double ptEdges[NE] = {50,60,70,80,90,100,120,150,200,300};
// NOTE: 9 bins 50..300 would drop 300-500; keep the full pp axis instead:
#ifdef BJET_PT_7BINS   // 7-bin variant (caloJetBPurityFitCore.h): 80-200 GeV in the scans' native edges
const int NEdge = 13; const double edges[NEdge] = {50, 60, 70, 80, 90, 100, 110, 120, 135, 150, 200, 300, 500};
#else
const int NEdge = 10; const double edges[NEdge] = {50,60,70,80,100,120,150,200,300,500};
#endif
const double ptReportMin = 80.;
double recoFloor = 70.;   // caloJets 70 (as pp); caloMatchedPF 50: the PbPb data reach 20 GeV and every window from 50 has a purity
const double ptDrawMax   = 200.;   // ratio figures: the bins above have no purity fit in every class
int nIter[4] = {1, 1, 1, 1};   // 0-10, 10-30, 30-50, 50-80%; see ITERATIONS in the header
// caloMatchedPF: N = 2 in every class, the minimum MSE of the same closure test on
// the 2026-10-06 calo-matched anaMuTag response merged to coarse classes
// (tools/mergeUltraFineToNominalCent_response.C), tiltScan:0.3, floor 50, run 2026-10-07
const int nIterCaloMatchedPF[4] = {2, 2, 2, 2};
// PFJets (plain PF both sides, sibling-repo 2026-02-18 muTaggedJetsNoTrigger response,
// nominal classes): min MSE of unfoldClosureTest.C("C<n>",15,"weighted","tilt:0.3",80.,500.,
// false,1,<resp>,<resp>,...,"h_matchedRecoJetPt_genJetPt_bJets_C<n>","",true,50.,<b edges>),
// 2026-10-07 (train = test). calculateBJetsPerZ.cc uses 5/4/3/2 with the same response.
// Since 2026-10-08 the PFJets unfolding runs on a truth axis extended to 20 GeV,
// floor 70, with single-event weight outliers of the fine response clipped
// (bJetUnfoldExt.h, clipResponseOutliers.h, K = 100). These N were tuned on the
// old axis (truth from 50, floor 50) and still need a closure re-check.
const int nIterPFJets[4] = {4, 4, 3, 3};   // superseded 2026-10-09: bNIterPbPb (bJetUnfoldExt.h)
// CALO-MATCH ID (caloMatchedPF with a plain-PF pp reference, ppColl = PFJets).
// In PYTHIA+HYDJET the ID keeps only ~50% of muon-tagged b jets at reco 80-90
// GeV, so the calo-matched response's reco axis is depleted and its truth lacks
// those jets -- but the PbPb DATA lose far fewer (calo-matched / plain-PF
// SingleMuon: 0.77-0.94 at 80-90 GeV, >= 0.95 above 110, and in central partly
// fake rejection): the MC ID does not describe the data (measureCaloMatchEff.C).
// So the PbPb data are taken as ID-free and unfolded, like the plain-PF pp data,
// with an ID-FREE response: each reco row of the calo-matched matrix scaled by
// w = 1 + f (1/eps_MC - 1), eps_MC = h_eff_muTagB_C<n> at the row's reco pT,
// f = idCorrFrac = 1 (rows with eps_MC < idEpsMin unscaled); matrix and truth
// both from the scaled matrix. Reco floor 70 (eps_MC >= 0.27 there; 50 would
// put weights up to 40 on a few surviving jets). bJetsPerZSys_PF_PbPb.C varies
// f to 0.5 / 1.5 and divides the data by eps_data (h_effData_muTag_C<n>).
const char  *idEffPath = "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/correctionFactors/caloMatchEff.root";
const double idEpsMin  = 0.1;
// MUON-TAG FREQUENCY FROM pp (user decision 2026-10-05; applied here for the
// caloMatchedPF / PFJets pairing). PYTHIA+HYDJET's gen-muon-tag f sits 1.4-1.6x
// above PYTHIA's at 80 GeV in EVERY class, 50-80% included -- a sample artifact --
// while its class dependence follows the muon efficiency. So
//   f(c) = f_pp x muEff[c] / 0.9708,
// f_pp = the pp chain's h_muTagFrequency (analysis reco-muon tag, which already
// holds PYTHIA's 0.9708), and the spectrum is divided by f(c) alone (muEff[c]
// is inside it). f_pp's uncertainty is common to pp and PbPb and cancels in
// PbPb/pp. The PYTHIA+HYDJET f is kept as muTagFreqPH_<cls>, a cross-check.
const double effMuPP = 0.9708;
double idCorrFrac = 1.;
// fine matrix with reco rows scaled by the ID weight (see above)
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
double muEffSim[4] = {0.8627, 0.9069, 0.9856, 0.9778};   // tight | gen, as calculateJetsPerZ.cc; cent "ATLAS": Z-weighted mix
const double fitLo = 0., fitHi = 5.;

// sum the T0 template over the class GROUP's 5% slices: 0-30% for the two
// central classes, 30-80% for the two peripheral ones (statistics, see header)
TH2D* tplSum(TFile *f, const char *fl, int ci){
  TH2D *s = nullptr;
  const int gLo = ci < 2 ? coarseSliceLo[0] : coarseSliceLo[2], gHi = ci < 2 ? coarseSliceHi[1] : coarseSliceHi[3];
  for(int si = gLo; si <= gHi; si++){
    TH2D *h = nullptr; f->GetObject(Form("h_muptrel_recoJetPt_inclRecoMuonTag_triggerOn_%sJets_C%dT0", fl, si), h);
    if(!h) continue;
    if(!s){ s = (TH2D*) h->Clone(Form("tpl_%s_%d", fl, ci)); s->SetDirectory(nullptr); } else s->Add(h);
  }
  return s;
}
// nameFmt takes the slice index, e.g. "h_matchedRecoJetPt_genJetPt_bJets_C%d_muTagged"
TH2D* sliceSumFmt2(TFile *f, const char *nameFmt, int ci, const char *tag){
  TH2D *s = nullptr;
  if(!TString(f->GetName()).Contains("ultraFineCentBins")){   // nominal classes C1..C4
    if(coarseSchemeATLAS) return nominalMix<TH2D>(f, nameFmt, ci, Form("%s_fmt", tag));   // ATLAS classes: mixed (coarseCent.h)
    f->GetObject(Form(nameFmt, ci + 1), s);
    if(s){ s = (TH2D*) s->Clone(Form("%s_fmt_%d", tag, ci)); s->SetDirectory(nullptr); }
    return s;
  }
  for(int si = coarseSliceLo[ci]; si <= coarseSliceHi[ci]; si++){
    TH2D *h = nullptr; f->GetObject(Form(nameFmt, si), h);
    if(!h) continue;
    if(!s){ s = (TH2D*) h->Clone(Form("%s_fmt_%d", tag, ci)); s->SetDirectory(nullptr); } else s->Add(h);
  }
  return s;
}
TH2D* sliceSum2(TFile *f, const char *base, int ci, const char *tag){
  TH2D *s = nullptr;
  if(!TString(f->GetName()).Contains("ultraFineCentBins")){   // nominal classes C1..C4
    if(coarseSchemeATLAS) return nominalMix<TH2D>(f, Form("%s_C%%d", base), ci, Form("%s_%s", base, tag));
    f->GetObject(Form("%s_C%d", base, ci + 1), s);
    if(s){ s = (TH2D*) s->Clone(Form("%s_%s_%d", base, tag, ci)); s->SetDirectory(nullptr); }
    return s;
  }
  for(int si = coarseSliceLo[ci]; si <= coarseSliceHi[ci]; si++){
    TH2D *h = nullptr; f->GetObject(Form("%s_C%d", base, si), h);
    if(!h) continue;
    if(!s){ s = (TH2D*) h->Clone(Form("%s_%s_%d", base, tag, ci)); s->SetDirectory(nullptr); } else s->Add(h);
  }
  return s;
}
// 0-5 GeV part of the scan's ptRel axis (100 bins over 0-10) in nFitBins bins:
// 0.2 GeV (two scan bins each; an uneven merge makes a sawtooth), coarser than
// the pp chain's 0.1 GeV because the PbPb calo templates are too thin for that
const int nFitBins = 25;
TH1D* range05(TH1D *h, const char *name){
  TH1D *r = new TH1D(name, "", nFitBins, 0., 5.); r->SetDirectory(nullptr); r->Sumw2();
  for(int b = 1; b <= 50; b++){
    int rb = r->FindBin(h->GetBinCenter(b));
    r->SetBinContent(rb, r->GetBinContent(rb) + h->GetBinContent(b));
    r->SetBinError(rb, std::hypot(r->GetBinError(rb), h->GetBinError(b)));
  }
  return r;
}
double effEntries(TH1D *h){ double w = 0., w2 = 0.; for(int b = 1; b <= h->GetNbinsX(); b++){ w += h->GetBinContent(b); w2 += h->GetBinError(b)*h->GetBinError(b); } return w2 > 0 ? w*w/w2 : 0.; }

struct Fit { double fb = -1., efb = 0., chi2 = 0.; int ndf = 0; bool ok = false; TH1D *d = nullptr, *tb = nullptr, *tl = nullptr; };

Fit fitWindow(TH2D *D, TH2D *Tb, TH2D *Tc, TH2D *Tu, TH2D *Td, TH2D *Ts, TH2D *Tg, double lo, double hi, const char *tag){
  Fit r;
  TH1D *d  = range05(projectY(D,  lo, hi, Form("d_%s", tag)),  Form("d05_%s", tag));
  TH1D *tb = range05(projectY(Tb, lo, hi, Form("b_%s", tag)),  Form("b05_%s", tag));
  // light+c shape: own window above 120 GeV, else the 120-150 GeV window (header)
  const double sLo = hi <= 120. + 1e-6 ? 120. : lo, sHi = hi <= 120. + 1e-6 ? 150. : hi;
  TH1D *tl = range05(projectY(Tu, sLo, sHi, Form("u_%s", tag)),  Form("l05_%s", tag));
  for(TH2D *T : {Td, Ts, Tg, Tc}) if(T) tl->Add(range05(projectY(T, sLo, sHi, Form("x_%s", tag)), Form("x05_%s", tag)));
  r.d = d; r.tb = tb; r.tl = tl;
  if(d->Integral() < 50 || effEntries(tb) < 30 || effEntries(tl) < 30){
    printf("    %5.0f-%5.0f: not fitted (data %.0f, b eff %.0f, light eff %.0f)\n", lo, hi, d->Integral(), effEntries(tb), effEntries(tl));
    return r;
  }
  // the likelihood is -inf if the data populate a bin where both templates are empty
  for(int b = 1; b <= nFitBins; b++) if(d->GetBinContent(b) > 0 && tb->GetBinContent(b) <= 0 && tl->GetBinContent(b) <= 0){
    printf("    %5.0f-%5.0f: not fitted, data in an empty-template bin at %.1f GeV\n", lo, hi, d->GetBinCenter(b)); return r; }
  RooRealVar x("x", "p_{T}^{rel}", fitLo, fitHi);
  RooDataHist dh("dh", "", x, d);
  RooDataHist hb("hb", "", x, tb), hl("hl", "", x, tl);
  RooHistPdf pb("pb", "", x, hb), pl("pl", "", x, hl);
  RooRealVar fb("fb", "b fraction", 0.5, 0., 1.);
  RooAddPdf model("model", "", RooArgList(pb, pl), RooArgList(fb));
  RooFitResult *res = model.fitTo(dh, RooFit::Save(true), RooFit::PrintLevel(-1), RooFit::Verbose(false), RooFit::SumW2Error(false));
  if(!res || res->status() != 0){ printf("    %5.0f-%5.0f: fit did not converge\n", lo, hi); return r; }
  r.fb = fb.getVal(); r.efb = fb.getError();
  // chi2/ndf by hand on the fit model
  double nd = d->Integral(), nb = tb->Integral(), nl = tl->Integral(), c2 = 0.; int ndf = 0;
  for(int b = 1; b <= nFitBins; b++){
    double m = nd*(r.fb*tb->GetBinContent(b)/nb + (1 - r.fb)*tl->GetBinContent(b)/nl);
    double e2 = d->GetBinError(b)*d->GetBinError(b);
    if(e2 <= 0) continue;
    c2 += (d->GetBinContent(b) - m)*(d->GetBinContent(b) - m)/e2; ndf++;
  }
  r.chi2 = c2; r.ndf = ndf - 1; r.ok = true;
  return r;
}

TH2D* rebin2D(TH2D *h, const char *name){
  TH2D *r = new TH2D(name, "", NEdge-1, edges, NEdge-1, edges); r->SetDirectory(nullptr); r->Sumw2();
  for(int ix = 1; ix <= h->GetNbinsX(); ix++) for(int iy = 1; iy <= h->GetNbinsY(); iy++){
    double c = h->GetBinContent(ix, iy); if(c == 0.) continue;
    int bx = r->GetXaxis()->FindBin(h->GetXaxis()->GetBinCenter(ix)), by = r->GetYaxis()->FindBin(h->GetYaxis()->GetBinCenter(iy));
    if(bx < 1 || bx > NEdge-1 || by < 1 || by > NEdge-1) continue;
    r->SetBinContent(bx, by, r->GetBinContent(bx, by) + c);
    r->SetBinError(bx, by, std::hypot(r->GetBinError(bx, by), h->GetBinError(ix, iy)));
  }
  return r;
}
TH1D* flavorPt(TH2D *H, const char *name){   // |flavor| 5 and 17, on H's X (gen pT) axis, rebinned to edges
  TH1D *p = new TH1D(Form("%s_fine", name), "", H->GetNbinsX(), H->GetXaxis()->GetXmin(), H->GetXaxis()->GetXmax()); p->SetDirectory(nullptr); p->Sumw2();
  for(int fv : {-5, 5, 17}){
    int by = H->GetYaxis()->FindBin(fv + 0.1);
    for(int ix = 1; ix <= H->GetNbinsX(); ix++){ p->SetBinContent(ix, p->GetBinContent(ix) + H->GetBinContent(ix, by)); p->SetBinError(ix, std::hypot(p->GetBinError(ix), H->GetBinError(ix, by))); }
  }
  TH1D *r = (TH1D*) p->Rebin(NEdge-1, name, edges); r->SetDirectory(nullptr); delete p; return r;
}
bool reported(int b){ return edges[b-1] >= ptReportMin - 1e-6; }

} // namespace bcal
using namespace bcal;

// purityMode  "chain"      the stop-gap fit above (shared class groups, borrowed light shape)
//             "dedicated"  the purity of systematics_caloJets_PbPb.C (h_bPurity_nominal_C*):
//                          smeared-pp light template, b and c shapes from the
//                          class's OWN PYTHIA+HYDJET templates, no sharing; fit to
//                          the light+c 0.98 quantile, purity over 0-5 GeV. Outputs
//                          carry the tag _dedicatedPurity.
// coll    PbPb jet collection (jetCollectionInputs.h; "" = BJET_JETCOLL env or caloJets).
//          caloMatchedPF forces purityMode "dedicated" (its stop-gap fit needs the calo
//          template scan) and takes the purity from plotBPurity_caloJets_PbPb.C's
//          bPurity_caloMatchedPF_PbPb.root (h_bPurity_2tpl_GSplus0p175_C*) when no
//          bPuritySys_caloMatchedPF_PbPb.root exists; reco floor 50.
// ppColl   collection of the pp reference ("" = coll). When its pp chain output
//          (bJetSpectrum_inclusive_<ppColl>_pp.root) does not exist, the pp side,
//          the PbPb/pp ratios and the b / inclusive ratios are skipped and only the
//          per-class b jets per Z and R_CP are produced. The inclusive jets per Z on
//          the b-jet binning come from inclSpecPath (inclusiveJetsPerZSys_caloJets.C,
//          closure-optimum N) when it exists, else from calculateJetsPerZ.cc's own
//          file (jetsPerZ_fine_<cls>: unfolded at ITS N, 1, counts per 5 GeV bin per
//          Z) summed into the b-jet bins; the source is printed and stored.
// Outputs carry _pp<ppColl> when ppColl differs from coll.
// cent     "" nominal classes; "ATLAS" (2026-10-15, PFJets only): C1..C4 = 0-20 / 20-50 /
//          50-80 / 0-80% for the comparison with arXiv:2204.13530 and LIDO. Data
//          exact (slices); the nominal-class February response MIXED from the
//          overlapping nominal classes (coarseCent.h setCoarseSchemeATLAS /
//          nominalMix); purity bPuritySys_PFJets_PbPb_ATLASCent.root
//          (systematics_caloJets_PbPb.C("PFJets","ATLAS"), templates mixed the same
//          way); N_Z and inclusive from inclusiveJetsPerZSys_PF.C(true,"MBcaloMatched",
//          "","","ATLAS","caloMatched",true) (calo-matched ultraFine Jet80/Jet100);
//          muEff = Z-weighted mix of the nominal classes'. Outputs carry _ATLASCent.
//          "R_CP" is class / 0-80% here.
// inclTag (PFJets, 2026-10-21): read inclJetsPerZSys_PF_MBcaloMatched_JEUppCorr<inclTag>.root for the
// inclusive side (b / inclusive) instead of the plain file, and append outTag to the output names, e.g.
// ("dedicated","PFJets","","","_noFlatJEU_fakeGeoOnly_evtCountStitch","_evtCountStitch") for the
// event-count trigger stitch. The b-jet spectra do not depend on it.
void bJetChain_caloJets_PbPb(const char *purityMode = "chain", const char *coll = "", const char *ppColl = "", TString cent = "",
                             TString inclTag = "", TString outTag = ""){
  setJetCollection(coll);
  if(cent != "" && (cent != "ATLAS" || jetColl != kPFJets)){ printf("ERROR: cent is \"\" or ATLAS (PFJets only)\n"); return; }
  const bool atlas = cent == "ATLAS";
  // N_Z from the published CMS yields x luminosity, NOMINAL since 2026-11-09 (user;
  // headers/publishedData/CMS_ZYield_5p02TeV.h). BJET_NZ_SOURCE=own = our own count from
  // the inclusive file, _ownNZ appended to outTag. cent ATLAS: setCoarseSchemeATLAS sets
  // coarseSliceLo/Hi to the ATLAS slices before N_Z is read, so the override follows them.
  // Plain-PF chain only; the legacy calo / calo-matched chains keep their own count.
  const bool pubNZ = jetColl == kPFJets && pubZ::enabled();
  if(jetColl == kPFJets) outTag += pubZ::tag();
  if(!ptBinGuard(outTag)) return;   // 7-bin jet-pT axis variant (caloJetBPurityFitCore.h)
  const bool pf = jetColl == kCaloMatchedPF;
  TString pm = purityMode;
  if(jetColl != kCaloJets && pm == "chain"){ printf("NOTE: %s has no stop-gap template scan: purityMode -> dedicated\n", collTag()); pm = "dedicated"; }
  const bool dedicated = pm == "dedicated";
  if(!dedicated && pm != "chain"){ printf("ERROR: purityMode is chain or dedicated\n"); return; }
  const TString ppTag = TString(ppColl) == "" ? collTag() : ppColl;
  const TString tag = TString(dedicated ? "_dedicatedPurity" : "") + (ppTag != collTag() ? Form("_pp%s", ppTag.Data()) : "") + (atlas ? "_ATLASCent" : "") + outTag;
  const char *base = "/home/clayton/Analysis/code/bJetRaaAnalysis";
  inclPath     = strdup(Form("%s/src/calculateJetsPerZ/rootFiles/JetsPerZ/%s", base, ppInclJetsPerZFile));   // the PbPb collection's file: N_Z and PbPb inclusive
  // the pp collection's jets-per-Z file for the pp inclusive side (same N_Z: one SingleMuon scan)
  const char *ppInclPath = strdup(Form("%s/src/calculateJetsPerZ/rootFiles/JetsPerZ/%s", base,
    ppTag == "PFJets" ? "jetsPerZSpectra_PFJets_coarseMuEff_unfold1iter.root" : ppTag == "caloMatchedPF" ? "jetsPerZSpectra_PFJets_caloMatchedPF_coarseMuEff_unfold1iter.root" : "jetsPerZSpectra_caloJets_coarseMuEff_unfold1iter.root"));
  inclSpecPath = strdup(Form("%s/rootFiles/systematics/inclJetsPerZSys_%s_JEUppUncorr_tilt0p3.root", base, collTag()));
  // calo-matched PbPb with plain-PF pp: inclusiveJetsPerZSys_PF.C (closure-optimum N per system)
  if(pf && ppTag == "PFJets") inclSpecPath = strdup(Form("%s/rootFiles/systematics/inclJetsPerZSys_caloMatchedPF_ppPFJets_JEUppCorr.root", base));
  if(pf) for(int c = 0; c < 4; c++) nIter[c] = nIterCaloMatchedPF[c];
  const bool plainPF = jetColl == kPFJets;
  if(plainPF){
    for(int c = 0; c < 4; c++) nIter[c] = bNIterPbPb[c];   // bJetUnfoldExt.h
    // floor 70 and truth from 20 GeV (bJetUnfoldExt.h), 2026-10-08: with truth
    // from 50 and floor 50, reco 50-70 (63% gen < 50 in 0-10%) had no truth bin
    // to go to and pulled central R_Z 10-40% low at 80-150 GeV
    recoFloor = bRecoFloorPF();   // 70; 50 with BJET_UNFOLD_LOWFLOOR=1 (bJetUnfoldExt.h)
    inclSpecPath = strdup(Form("%s/rootFiles/systematics/inclJetsPerZSys_PF_MBcaloMatched_JEUppCorr%s.root", base, inclTag.Data()));
  }
  const bool idCorr = pf && ppTag == "PFJets";
  if(idCorr) recoFloor = 70.;
  const bool tagFreqPP = idCorr || plainPF;
  TFile *fIDE = idCorr ? TFile::Open(idEffPath) : nullptr;
  if(idCorr && (!fIDE || fIDE->IsZombie())){ printf("ERROR: %s missing -- run measureCaloMatchEff.C\n", idEffPath); return; }
  ppBPath      = strdup(Form("%s/rootFiles/CorrectedBJetSpectra/Data/bJetSpectrum_inclusive_%s_pp%s.root", base, ppTag.Data(), ppPurVar()));
  TH1D *fPP = nullptr;
  if(tagFreqPP){
    TFile *fPf = TFile::Open(ppBPath);
    if(fPf && !fPf->IsZombie()) fPf->GetObject("h_muTagFrequency", fPP);
    if(!fPP || fPP->GetNbinsX() != NEdge-1){ printf("ERROR: h_muTagFrequency missing in %s or on another axis\n", ppBPath); return; }
    fPP->SetDirectory(nullptr);
    printf("muon-tag frequency: f(c) = f_pp x muEff[c] / %.4f (f_pp from %s)\n", effMuPP, gSystem->BaseName(ppBPath));
  }
  purSysPath   = strdup(Form("%s/rootFiles/systematics/bPuritySys_%s_PbPb%s%s.root", base, collTag(), atlas ? "_ATLASCent" : "", purVar()));
  if(atlas){   // N_Z and inclusive in the ATLAS classes, b-jet axis (header)
    inclSpecPath = strdup(Form("%s/rootFiles/systematics/inclJetsPerZSys_PF_MBcaloMatched_JEUppCorr_ATLASCent_caloMatched_bjetAxis.root", base));
    inclPath = inclSpecPath;
    TFile *fZ = TFile::Open(dataPathPbPb);
    if(!fZ || fZ->IsZombie() || !setCoarseSchemeATLAS(fZ)) return;
    const double mE[4] = {0.8627, 0.9069, 0.9856, 0.9778};
    for(int c = 0; c < 4; c++){ muEffSim[c] = 0.; for(int k = 0; k < 4; k++) muEffSim[c] += coarseMixW[c][k]*mE[k]; }
    printf("muon efficiency per ATLAS class: %.4f %.4f %.4f %.4f\n", muEffSim[0], muEffSim[1], muEffSim[2], muEffSim[3]);
  }
  const char *purFilePath = strdup(Form("%s/rootFiles/bPurityResults/bPurity_%s_PbPb.root", base, collTag()));
  outRoot = Form("%s/rootFiles/BJetPerZ_%s/bJetChain_%s_PbPb%s.root", base, collTag(), collTag(), tag.Data());
  figDir  = Form("%s/figures/bJetsPerZ_%s%s/", base, collTag(), tag.Data());
  if(pf && !(ppTag == "PFJets")) recoFloor = 50.;
  initPlotStyle();
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  gSystem->mkdir(figDir.Data(), true); gSystem->mkdir(gSystem->DirName(outRoot.Data()), true);
  const bool havePP = !gSystem->AccessPathName(ppBPath), haveInclSys = !gSystem->AccessPathName(inclSpecPath);
  const bool haveIncl = true;   // inclSpecPath, else calculateJetsPerZ.cc's fine spectra (below)
  if(!havePP) printf("NOTE: no pp chain output for the %s collection (%s): pp side and the PbPb/pp ratios skipped\n", ppTag.Data(), ppBPath);
  if(!haveInclSys) printf("NOTE: no %s: inclusive jets per Z from calculateJetsPerZ.cc's jetsPerZ_fine_<cls> (its own N = 1 unfolding) summed into the b-jet bins\n", inclSpecPath);
  printf("data %s\nresponse %s (%s), reco floor %.0f, purity mode %s\n", gSystem->BaseName(dataPathPbPb), gSystem->BaseName(phRespPath), phRespMuTagFmt, recoFloor, pm.Data());
  TFile *fD = TFile::Open(dataPathPbPb), *fT = dedicated ? nullptr : TFile::Open(tplPath), *fR = TFile::Open(phRespPath), *fI = TFile::Open(inclPath);
  TFile *fP = havePP ? TFile::Open(ppBPath) : nullptr, *fIS = haveInclSys ? TFile::Open(inclSpecPath) : nullptr;
  TFile *fPI = TString(ppInclPath) == inclPath ? nullptr : TFile::Open(ppInclPath);
  if(TString(ppInclPath) != inclPath){ if(!fPI || fPI->IsZombie()){ printf("ERROR: %s did not open\n", ppInclPath); return; } printf("pp inclusive jets per Z from %s\n", gSystem->BaseName(ppInclPath)); }
  for(TFile *f : {fD, fR, fI}) if(!f || f->IsZombie()){ printf("ERROR: an input did not open\n"); return; }
  if(!dedicated && (!fT || fT->IsZombie())){ printf("ERROR: the template scan did not open\n"); return; }
  if(havePP && (!fP || fP->IsZombie())){ printf("ERROR: %s did not open\n", ppBPath); return; }
  if(haveInclSys && (!fIS || fIS->IsZombie())){ printf("ERROR: %s did not open\n", inclSpecPath); return; }
  TVectorD *NZ = (TVectorD*) fI->Get("NZ");
  if(!NZ){ printf("ERROR: NZ missing from the inclusive file\n"); return; }
  printf("N_Z (inclusive chain): pp %.0f  %s %.0f  %s %.0f  %s %.0f  %s %.0f\n", (*NZ)[0], coarseLabel[0], (*NZ)[1], coarseLabel[1], (*NZ)[2], coarseLabel[2], (*NZ)[3], coarseLabel[3], (*NZ)[4]);
  const TVectorD nzOwn(*NZ);   // the N_Z of the inclusive files read (our own count unless they are _pubNZ files)
  if(pubNZ){
    (*NZ)[0] = pubZ::nzPP();
    for(int c = 0; c < 4; c++) (*NZ)[c+1] = pubZ::nzSlices(coarseSliceLo[c], coarseSliceHi[c]);
    printf("N_Z PUBLISHED x lumi (nominal): pp %.0f  %s %.0f  %s %.0f  %s %.0f  %s %.0f  (pp f_y %.3f)\n",
           (*NZ)[0], coarseLabel[0], (*NZ)[1], coarseLabel[1], (*NZ)[2], coarseLabel[2], (*NZ)[3], coarseLabel[3], (*NZ)[4], pubZ::fYpp);
  }

  // dedicated purity: the systematics file (h_bPurity_nominal_C*) when it
  // exists, else plotBPurity_caloJets_PbPb.C's results (h_bPurity_2tpl_GSplus0p175_C*)
  const bool purFromSys = dedicated && !gSystem->AccessPathName(purSysPath);
  const char *purSrc = purFromSys ? purSysPath : purFilePath, *purHist = purFromSys ? "h_bPurity_nominal_%s" : "h_bPurity_2tpl_GSplus0p175_%s";
  TFile *fPS = dedicated ? TFile::Open(purSrc) : nullptr;
  if(dedicated && (!fPS || fPS->IsZombie())){ printf("ERROR: %s did not open -- run systematics_caloJets_PbPb.C or plotBPurity_caloJets_PbPb.C first\n", purSrc); return; }
  TFile *out = TFile::Open(outRoot, "RECREATE");
  const char *cls[5] = {"pp", "C1", "C2", "C3", "C4"};
  TH1D *bPerZ[5] = {nullptr}, *iPerZ[5] = {nullptr}, *purity[4] = {nullptr};

  // ---- pp b spectrum (pp chain, dN/dpT) per Z; inclusive per Z on the b axis
  if(havePP){ TH1D *h = nullptr; fP->GetObject("h_bJetPt_inclusive", h);
    if(!h){ printf("ERROR: h_bJetPt_inclusive missing in %s\n", ppBPath); return; }
    bPerZ[0] = (TH1D*) h->Clone("bJetsPerZ_pp"); bPerZ[0]->SetDirectory(nullptr); bPerZ[0]->Scale(1./(*NZ)[0]); }
  for(int j = 0; j < 5; j++){
    TH1D *r = nullptr;
    if(haveInclSys){
      TH1D *h = nullptr; fIS->GetObject(Form("h_nom_%s_bjet", cls[j]), h);
      if(!h){ printf("ERROR: h_nom_%s_bjet missing\n", cls[j]); return; }
      for(int b = 0; b < NEdge; b++) if(fabs(h->GetXaxis()->GetBinUpEdge(b) - edges[b]) > 1e-6){ printf("ERROR: h_nom_%s_bjet axis differs from edges\n", cls[j]); return; }
      r = (TH1D*) h->Clone(Form("jetsPerZ_%s", cls[j])); r->SetDirectory(nullptr);
    } else {
      // counts per 5 GeV bin per Z, summed into edges by bin centre (as inclusiveJetsPerZSys_caloJets.C's rebinAx)
      TH1D *h = nullptr; (j == 0 && fPI ? fPI : fI)->GetObject(Form("jetsPerZ_fine_%s", cls[j]), h);
      if(!h){ printf("ERROR: jetsPerZ_fine_%s missing in %s\n", cls[j], j == 0 && fPI ? ppInclPath : inclPath); return; }
      r = new TH1D(Form("jetsPerZ_%s", cls[j]), "", NEdge-1, edges); r->SetDirectory(nullptr); r->Sumw2();
      for(int i = 1; i <= h->GetNbinsX(); i++){
        int k = r->FindBin(h->GetBinCenter(i)); if(k < 1 || k > NEdge-1) continue;
        r->SetBinContent(k, r->GetBinContent(k) + h->GetBinContent(i));
        r->SetBinError(k, std::hypot(r->GetBinError(k), h->GetBinError(i)));
      }
    }
    divideByBinwidth(r); iPerZ[j] = r;
    if(pubNZ){   // inclusive per Z: divided by its file's N_Z, moved to the published N_Z
      double src = nzOwn[j];
      TVectorD *v = haveInclSys ? (TVectorD*) fIS->Get("NZ") : nullptr;
      if(v) src = (*v)[j];
      r->Scale(src/(*NZ)[j]);
    }
  }

  // ---- PbPb classes
  for(int ci = 0; ci < NCoarse; ci++){
    printf("\n===== PbPb %s =====\n", coarseLabel[ci]);
    TH2D *D = nullptr, *Tb = nullptr, *Tg2 = nullptr, *Tc = nullptr, *Tu = nullptr, *Td = nullptr, *Ts = nullptr, *Tg = nullptr;
    if(!dedicated){
    D  = coarseSum2(fD, "h_muptrel_recoJetPt_inclRecoMuonTag_triggerOn", ci, "data");
    Tb = tplSum(fT, "b", ci); Tg2 = tplSum(fT, "bGS", ci); Tc = tplSum(fT, "c", ci);
    Tu = tplSum(fT, "u", ci); Td = tplSum(fT, "d", ci); Ts = tplSum(fT, "s", ci); Tg = tplSum(fT, "g", ci);
    if(!D || !Tb || !Tu){ printf("ERROR: inputs missing for %s\n", coarseLabel[ci]); return; }
    if(Tg2 && Tg2->GetEntries() > 0) Tb->Add(Tg2);
    }

    // purity per window
    purity[ci] = new TH1D(Form("bPurity_%s", cls[ci+1]), "", NEdge-1, edges); purity[ci]->SetDirectory(nullptr);
    if(dedicated){
      TH1D *hp = nullptr; fPS->GetObject(Form(purHist, cls[ci+1]), hp);
      if(!hp || hp->GetNbinsX() != NEdge-1){ printf("ERROR: %s missing or on another axis\n", Form(purHist, cls[ci+1])); return; }
      printf("  purity from %s (own-class b and c; light: %s):\n", gSystem->BaseName(purSrc), pf ? "PYTHIA+HYDJET" : "smeared pp");
      for(int b = 1; b < NEdge; b++){
        if(!reported(b) || hp->GetBinContent(b) <= 0.) continue;
        purity[ci]->SetBinContent(b, hp->GetBinContent(b)); purity[ci]->SetBinError(b, hp->GetBinError(b));
        printf("    %5.0f-%5.0f: f_b = %.3f +- %.3f\n", edges[b-1], edges[b], hp->GetBinContent(b), hp->GetBinError(b));
      }
    }
    printf("  purity (b vs light+c, 0-5 GeV):\n");
    for(int b = 1; b < NEdge && !dedicated; b++){
      if(!reported(b)) continue;
      Fit f = fitWindow(D, Tb, Tc, Tu, Td, Ts, Tg, edges[b-1], edges[b], Form("%s_%d", cls[ci+1], b));
      if(!f.ok) continue;
      purity[ci]->SetBinContent(b, f.fb); purity[ci]->SetBinError(b, f.efb);
      printf("    %5.0f-%5.0f: f_b = %.3f +- %.3f   chi2/ndf = %.1f/%d   (data %.0f)\n", edges[b-1], edges[b], f.fb, f.efb, f.chi2, f.ndf, f.d->Integral());
      // fit figure: shared panel (headers/plotting/ptRelFitPanel.h), densities
      { double nd = f.d->Integral();
        TH1D *sb = (TH1D*) f.tb->Clone(Form("sb_%s_%d", cls[ci+1], b)); sb->SetDirectory(nullptr); sb->Scale(nd*f.fb/sb->Integral());
        TH1D *sl = (TH1D*) f.tl->Clone(Form("sl_%s_%d", cls[ci+1], b)); sl->SetDirectory(nullptr); sl->Scale(nd*(1-f.fb)/sl->Integral());
        TH1D *sum = (TH1D*) sb->Clone(Form("sum_%s_%d", cls[ci+1], b)); sum->SetDirectory(nullptr); sum->Add(sl);
        TH1D *dd = (TH1D*) f.d->Clone(Form("dd_%s_%d", cls[ci+1], b)); dd->SetDirectory(nullptr);
        for(TH1D *h : {dd, sb, sl, sum}) h->Scale(1., "width");
        drawPtRelFitPanel(dd, sum, {{sb, Form("b (%.3f #pm %.3f)", f.fb, f.efb), hexPtRelB, 2}, {sl, "light + c", hexPtRelL, 2}},
                          Form("PbPb %s, %s", coarseLabel[ci], collLabel()), edges[b-1], edges[b], f.chi2, f.ndf, fitHi,
                          Form("%sptRelFit_%s_%.0f-%.0f.pdf", figDir.Data(), coarseTag[ci], edges[b-1], edges[b])); }
    }
    // buffer bins (reco below the first reported window, which feed the unfolding): NOMINAL since 2026-11-10 (user) the
    // purity file's own buffer values, as the pp chain does (they track the low-pT CSVv2 fits). Before, they took the first
    // reported purity, which overstated the 50-80% purity at reco 70-80 GeV by 18% and raised its 80-100 GeV double ratio
    // by 5%; BJET_PBPB_BUFFER_PURITY=first restores that (and is used when the file has no buffer value)
    int firstRep = 1; while(firstRep < NEdge && !reported(firstRep)) firstRep++;
    static const bool bufFromFile = [](){ const char *e = std::getenv("BJET_PBPB_BUFFER_PURITY"); const bool first = e && TString(e) == "first";
                                          if(first) printf("PbPb BUFFER PURITY = the first reported window (BJET_PBPB_BUFFER_PURITY=first, pre-2026-11-10)\n"); return !first; }();
    TH1D *hpBuf = nullptr; if(dedicated && bufFromFile) fPS->GetObject(Form(purHist, cls[ci+1]), hpBuf);
    for(int b = 1; b < firstRep; b++){
      if(hpBuf && hpBuf->GetBinContent(b) > 0.){ purity[ci]->SetBinContent(b, hpBuf->GetBinContent(b)); purity[ci]->SetBinError(b, hpBuf->GetBinError(b)); continue; }
      purity[ci]->SetBinContent(b, purity[ci]->GetBinContent(firstRep)); purity[ci]->SetBinError(b, purity[ci]->GetBinError(firstRep)); }
    out->cd(); purity[ci]->Write();

    // measured b spectrum, counts
    TH1D *tag0 = coarseSum1(fD, "h_inclRecoJetPt_inclRecoMuonTag_triggerOn", ci, "tag");
    TH1D *tag = (TH1D*) tag0->Rebin(NEdge-1, Form("muTagJetPt_%s", cls[ci+1]), edges); tag->SetDirectory(nullptr);
    TH1D *meas = (TH1D*) tag->Clone(Form("bJetPt_meas_%s", cls[ci+1])); meas->SetDirectory(nullptr); meas->Multiply(purity[ci]);
    // windows without a purity fit are empty here and are not reported below
    out->cd(); tag->Write(); meas->Write();

    // response, floor, unfold
    TH2D *R0 = sliceSumFmt2(fR, phRespMuTagFmt, ci, "resp");
    if(!R0){ printf("ERROR: response missing for %s\n", coarseLabel[ci]); return; }
    if(plainPF) clipResponseOutliers(R0, 100., true, cls[ci+1]);
    if(plainPF){   // JER nominal (2026-11-05, user): the response smeared to a 10% worse resolution (jerSmear.h);
                   // bJetsPerZSys_PF_PbPb.C does the same before rebuilding this unfolding
      TH2D *sm = jerSmear::nominal(R0, Form("respJERnom_%s", cls[ci+1])); delete R0; R0 = sm;
      if(ci == 0) jerSmear::announce("PbPb b response");
    }
    const std::vector<double> extE = bExtEdges(NEdge, edges);
    // PFJets: data-driven prior (dataDrivenPrior.h) -- the fine response's gen rows
    // x (g/100)^-k, k fitted so the prior follows the unfolded data over 80-300 GeV;
    // stored as priorTilt_<cls> for bJetsPerZSys_PF_PbPb.C, which takes the MC
    // prior (k = 0) as the "prior" systematic
    double priorTilt = 0.;
    if(plainPF && !bUseDataDrivenPrior){ out->cd(); TVectorD kv(1); kv[0] = 0.; kv.Write(Form("priorTilt_%s", cls[ci+1])); }   // MC prior (bJetUnfoldExt.h)
    if(plainPF && bUseDataDrivenPrior){
      auto unfoldFine = [&](TH2D *fine, TH1D *&truthOut){
        TH2D *r2 = rebin2DTo(fine, extE, "pfR");
        TH1D *U = unfoldOnExt<RooUnfoldResponse, RooUnfoldBayes>(r2, meas, NEdge, edges, recoFloor, nIter[ci], "pfU", &truthOut);
        delete r2; return U; };
      priorTilt = fitPriorTilt(R0, unfoldFine, 80., bPriorFitHi(), 8, true);
      TH2D *t = tiltGenRows(R0, priorTilt, Form("respTilted_%s", cls[ci+1])); delete R0; R0 = t;
      printf("  data-driven prior: response gen rows x (g/100)^%+.3f\n", -priorTilt);
      out->cd(); TVectorD kv(1); kv[0] = priorTilt; kv.Write(Form("priorTilt_%s", cls[ci+1]));
    }
    TH2D *R = plainPF ? rebin2DTo(R0, extE, Form("resp_%s", cls[ci+1])) : rebin2D(R0, Form("resp_%s", cls[ci+1]));
    TH1D *truth = plainPF ? nullptr : R->ProjectionY(Form("truth_%s", cls[ci+1]), 1, NEdge-1);
    if(truth) truth->SetDirectory(nullptr);
    TH1D *U = nullptr;
    if(plainPF){   // extended truth axis (both axes from 20 GeV), result on the analysis axis
      U = unfoldOnExt<RooUnfoldResponse, RooUnfoldBayes>(R, meas, NEdge, edges, recoFloor, nIter[ci], Form("bJetPt_unf_%s", cls[ci+1]));
    } else {
    if(idCorr){   // ID-free response (CALO-MATCH ID)
      TH1D *eps = nullptr; fIDE->GetObject(Form("h_eff_muTagB_%s", cls[ci+1]), eps);
      if(!eps){ printf("ERROR: h_eff_muTagB_%s missing\n", cls[ci+1]); return; }
      TH2D *R0w = idScaleRows(R0, eps, idCorrFrac, Form("respIDw0_%s", cls[ci+1])), *Rw = rebin2D(R0w, Form("respIDw_%s", cls[ci+1]));
      TH1D *tw = Rw->ProjectionY(Form("truthIDw_%s", cls[ci+1]), 1, NEdge-1); tw->SetDirectory(nullptr);
      TH1D *eg = (TH1D*) truth->Clone(Form("caloMatchIDEffGen_MC_%s", cls[ci+1])); eg->SetDirectory(nullptr); eg->Divide(tw);
      printf("  ID-free response; the MC calo-match ID efficiency it undoes, at gen level:");
      for(int b = 1; b < NEdge; b++) if(reported(b)) printf("  %3.0f-%3.0f %.3f", edges[b-1], edges[b], eg->GetBinContent(b));
      printf("\n");
      out->cd(); eg->Write();
      delete truth; truth = tw; delete R; R = Rw; R->SetName(Form("resp_%s", cls[ci+1])); delete R0w;
    }
    TH2D *Rt = (TH2D*) R->Clone(Form("respTrunc_%s", cls[ci+1])); Rt->SetDirectory(nullptr);
    for(int ix = 1; ix < NEdge; ix++) if(edges[ix] <= recoFloor + 1e-6) for(int iy = 1; iy < NEdge; iy++){ Rt->SetBinContent(ix, iy, 0.); Rt->SetBinError(ix, iy, 0.); }
    TH1D *measMC = Rt->ProjectionX(Form("measMC_%s", cls[ci+1]), 1, NEdge-1); measMC->SetDirectory(nullptr);
    TH1D *in = (TH1D*) meas->Clone(Form("in_%s", cls[ci+1])); in->SetDirectory(nullptr);
    for(int b = 1; b < NEdge; b++) if(edges[b] <= recoFloor + 1e-6){ in->SetBinContent(b, 0.); in->SetBinError(b, 0.); }
    RooUnfoldResponse resp(measMC, truth, Rt, Form("resp_%s", cls[ci+1]), "");
    RooUnfoldBayes unf(&resp, in, nIter[ci]);
    U = (TH1D*) unf.Hunfold(); U->SetName(Form("bJetPt_unf_%s", cls[ci+1])); U->SetDirectory(nullptr);
    }
    printf("  unfolded / measured (N = %d, floor %.0f):", nIter[ci], recoFloor);
    for(int b = 1; b < NEdge; b++) if(reported(b)) printf("  %3.0f-%3.0f %.2f", edges[b-1], edges[b], meas->GetBinContent(b) > 0 ? U->GetBinContent(b)/meas->GetBinContent(b) : -1.);
    printf("\n");
    out->cd(); U->Write(); R->Write();
    { TCanvas *c = new TCanvas(Form("cunf_%s", cls[ci+1]), "", 700, 650); c->SetLeftMargin(0.15); c->SetBottomMargin(0.13); c->SetLogy();
      TH1D *m = (TH1D*) meas->Clone(); divideByBinwidth(m); TH1D *u = (TH1D*) U->Clone(); divideByBinwidth(u);
      styleH(m, okabeHex[5], markOpenCircle, 1.0); styleH(u, hexData, markFilledCircle, 1.0); m->SetTitle(""); m->GetXaxis()->SetTitle("jet #it{p}_{T} [GeV]"); m->GetYaxis()->SetTitle("b jets / GeV");
      m->GetXaxis()->SetRangeUser(50., 500.); m->GetYaxis()->SetRangeUser(0.3*std::min(m->GetMinimum(1e-9), u->GetMinimum(1e-9)), 30.*std::max(m->GetMaximum(), u->GetMaximum())); m->Draw("e1"); u->Draw("e1 same");
      TLegend *lg = makeLegend(0.50, 0.68, 0.92, 0.84, 0.038); lg->AddEntry(m, "measured (reco), #times purity", "pe"); lg->AddEntry(u, Form("unfolded, N = %d", nIter[ci]), "pe"); lg->Draw();
      TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.042); t.DrawLatex(0.19, 0.86, Form("PbPb %s, b jets, %s", coarseLabel[ci], collLabel()));
      savePdfTight(c, Form("%sbJetUnfold_%s.pdf", figDir.Data(), coarseTag[ci])); delete c; }

    // corrections: muon efficiency, tag frequency, per Z
    TH2D *Gall = sliceSum2(fR, "h_inclGenJetPt_flavor", ci, "gall"), *Gtag = sliceSum2(fR, "h_inclGenJetPt_inclGenMuonTag_flavor", ci, "gtag");
    TH1D *fAll = flavorPt(Gall, Form("fAll_%s", cls[ci+1])), *fTag = flavorPt(Gtag, Form("fTag_%s", cls[ci+1]));
    TH1D *fr = (TH1D*) fTag->Clone(Form(tagFreqPP ? "muTagFreqPH_%s" : "muTagFreq_%s", cls[ci+1])); fr->SetDirectory(nullptr); fr->Divide(fTag, fAll, 1., 1., "B");
    printf("  muon-tag frequency f%s:", tagFreqPP ? " (PYTHIA+HYDJET, cross-check only)" : ""); for(int b = 1; b < NEdge; b++) printf("  %3.0f-%3.0f %.3f", edges[b-1], edges[b], fr->GetBinContent(b)); printf("\n");
    out->cd(); fr->Write();
    TH1D *corr = (TH1D*) U->Clone(Form("bJetsPerZ_%s", cls[ci+1])); corr->SetDirectory(nullptr);
    if(tagFreqPP){
      TH1D *fc = (TH1D*) fPP->Clone(Form("muTagFreq_%s", cls[ci+1])); fc->SetDirectory(nullptr);
      fc->Scale(muEffSim[ci]/effMuPP);
      printf("  muon-tag frequency f = f_pp x %.4f/%.4f:", muEffSim[ci], effMuPP); for(int b = 1; b < NEdge; b++) printf("  %3.0f-%3.0f %.3f", edges[b-1], edges[b], fc->GetBinContent(b)); printf("\n");
      out->cd(); fc->Write();
      for(int b = 1; b < NEdge; b++){ double f = fc->GetBinContent(b); if(f > 0.){ corr->SetBinContent(b, corr->GetBinContent(b)/f); corr->SetBinError(b, corr->GetBinError(b)/f); } }   // f_pp stat: a common systematic, not here
    } else {
      corr->Scale(1./muEffSim[ci]);
      corr->Divide(fr);
    }
    corr->Scale(1./(*NZ)[ci+1]);
    divideByBinwidth(corr);
    for(int b = 1; b < NEdge; b++) if(purity[ci]->GetBinContent(b) <= 0.){ corr->SetBinContent(b, 0.); corr->SetBinError(b, 0.); }   // no fit: not reported
    bPerZ[ci+1] = corr; out->cd(); corr->Write();
  }
  for(int j = 0; j < 5 && haveIncl; j++){ out->cd(); iPerZ[j]->Write(); } out->cd(); if(havePP) bPerZ[0]->Write();
  { TNamed collN("jetCollection", Form("PbPb %s, pp %s%s", collTag(), ppTag.Data(), havePP ? "" : " (no pp chain output: pp side skipped)")); collN.Write();
    TNamed purN("puritySource", dedicated ? purSrc : "stop-gap fit in bJetChain_caloJets_PbPb.C"); purN.Write();
    TNamed incN("inclusiveSource", haveInclSys ? inclSpecPath : Form("%s jetsPerZ_fine_<cls> (N = 1), summed into the b-jet bins", inclPath)); incN.Write(); }

  // ---- ratios
  auto ratio = [&](TH1D *n, TH1D *d, const char *name){ TH1D *r = (TH1D*) n->Clone(name); r->SetDirectory(nullptr); r->Divide(d);
    for(int b = 1; b < NEdge; b++) if(n->GetBinContent(b) <= 0. || d->GetBinContent(b) <= 0.){ r->SetBinContent(b, 0.); r->SetBinError(b, 0.); } return r; };
  TH1D *rB_pp[4] = {nullptr}, *rB_cp[3] = {nullptr}, *rBI_pp[4] = {nullptr}, *rBI_cp[3] = {nullptr}, *bOverIncl[4] = {nullptr};
  printf("\n  pT bin     class     b PbPb/pp        b C/P         (b/incl) PbPb/pp   (b/incl) C/P\n");
  for(int ci = 0; ci < 4; ci++){
    if(havePP) rB_pp[ci] = ratio(bPerZ[ci+1], bPerZ[0], Form("r_b_PbPbOverPP_%s", cls[ci+1]));
    if(ci < 3) rB_cp[ci] = ratio(bPerZ[ci+1], bPerZ[4], Form("r_b_CP_%s", cls[ci+1]));
    if(haveIncl){
      TH1D *bi = ratio(bPerZ[ci+1], iPerZ[ci+1], Form("bOverIncl_%s", cls[ci+1]));
      bOverIncl[ci] = bi; out->cd(); bi->Write();
      if(havePP){ TH1D *bi0 = ratio(bPerZ[0], iPerZ[0], "bOverIncl_pp"); rBI_pp[ci] = ratio(bi, bi0, Form("r_bOverIncl_PbPbOverPP_%s", cls[ci+1])); if(ci == 0){ out->cd(); bi0->Write(); } }
      if(ci < 3){ TH1D *bi4 = ratio(bPerZ[4], iPerZ[4], "bOverIncl_C4"); rBI_cp[ci] = ratio(bi, bi4, Form("r_bOverIncl_CP_%s", cls[ci+1])); }
    }
    for(int b = 1; b < NEdge; b++){ if(!reported(b)) continue;
      printf("  %3.0f-%3.0f   %-7s", edges[b-1], edges[b], coarseLabel[ci]);
      if(havePP) printf("   %.3f +- %.3f", rB_pp[ci]->GetBinContent(b), rB_pp[ci]->GetBinError(b)); else printf("        --       ");
      if(ci < 3) printf("   %.3f +- %.3f", rB_cp[ci]->GetBinContent(b), rB_cp[ci]->GetBinError(b)); else printf("        --       ");
      if(haveIncl){
        if(havePP) printf("   %.3f +- %.3f", rBI_pp[ci]->GetBinContent(b), rBI_pp[ci]->GetBinError(b)); else printf("        --       ");
        if(ci < 3) printf("   %.3f +- %.3f", rBI_cp[ci]->GetBinContent(b), rBI_cp[ci]->GetBinError(b));
      }
      printf("\n"); }
    out->cd(); if(havePP) rB_pp[ci]->Write(); if(ci < 3) rB_cp[ci]->Write();
    if(haveIncl){ if(havePP) rBI_pp[ci]->Write(); if(ci < 3) rBI_cp[ci]->Write(); }
  }
  if(haveIncl){
    printf("\n  b / inclusive jets per Z (inclusive: %s)\n  pT bin  ", haveInclSys ? gSystem->BaseName(inclSpecPath) : Form("%s jetsPerZ_fine, N = 1", gSystem->BaseName(inclPath)));
    for(int ci = 0; ci < 4; ci++) printf("   %-16s", coarseLabel[ci]); if(havePP) printf("   %-16s", "pp"); printf("\n");
    for(int b = 1; b < NEdge; b++){ if(!reported(b)) continue;
      printf("  %3.0f-%-4.0f", edges[b-1], edges[b]);
      for(int ci = 0; ci < 4; ci++) printf("   %.4f +- %.4f", bOverIncl[ci]->GetBinContent(b), bOverIncl[ci]->GetBinError(b));
      if(havePP){ TH1D *bi0 = (TH1D*) out->Get("bOverIncl_pp"); if(bi0) printf("   %.4f +- %.4f", bi0->GetBinContent(b), bi0->GetBinError(b)); }
      printf("\n"); }
  }

  // ---- figures
  const int mk[4] = {markFilledCircle, markFilledSquare, markFilledDiamond, markCross};
  // ratio figures in the layout of inclusiveJetsPerZSys_caloJets.C's ratioFig:
  // one panel per class, 2x2 for PbPb/pp and 3x1 for C/P. No systematic band
  // yet -- the b-jet systematics have not been run for PbPb.
  auto drawSet = [&](TH1D **h, int n, int nx, int ny, const char *ytitle, const char *name, const char *note){
    TCanvas *cv = new TCanvas(name, "", 520*nx, 460*ny);
    cv->Divide(nx, ny, 0.001, 0.001);
    for(int p = 0; p < n; p++){
      cv->cd(p + 1);
      gPad->SetLeftMargin(0.14); gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.08); gPad->SetBottomMargin(0.13);
      TH1D *v = (TH1D*) h[p]->Clone(); v->SetDirectory(nullptr);
      double ymax = 0.;
      for(int b = 1; b < NEdge; b++){
        if(!reported(b) || v->GetBinContent(b) <= 0.){ v->SetBinContent(b, -999.); v->SetBinError(b, 0.); continue; }
        if(edges[b-1] < ptDrawMax - 1e-6) ymax = TMath::Max(ymax, v->GetBinContent(b) + v->GetBinError(b));
      }
      TH1F *fr = gPad->DrawFrame(ptReportMin, 0., ptDrawMax, TMath::Max(1.3, 1.35*ymax));
      fr->GetXaxis()->SetTitle("#it{p}_{T}^{jet} [GeV]"); fr->GetYaxis()->SetTitle(ytitle);
      fr->GetXaxis()->SetTitleSize(0.050); fr->GetYaxis()->SetTitleSize(0.050);
      fr->GetXaxis()->SetLabelSize(0.045); fr->GetYaxis()->SetLabelSize(0.045); fr->GetYaxis()->SetTitleOffset(1.3);
      TLine *one = new TLine(ptReportMin, 1., ptDrawMax, 1.); one->SetLineStyle(7); one->SetLineColor(kGray+2); one->Draw();
      styleH(v, "#000000", markFilledCircle, 1.1); v->Draw("E1 X0 same");
      TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.050);
      la.DrawLatex(0.17, 0.935, n == 4 ? Form("%s / pp, %s", coarseLabel[p], collLabel()) : Form("%s / 50-80%%, %s", coarseLabel[p], collLabel()));
      if(p == 0){
        TLegend *lg = makeLegend(0.17, 0.78, 0.70, 0.90, 0.042);
        lg->AddEntry(v, "stat.", "pe"); lg->AddEntry((TObject*)nullptr, note, ""); lg->Draw();
      }
    }
    savePdfTight(cv, Form("%s%s.pdf", figDir.Data(), name)); delete cv;
  };
  if(havePP) drawSet(rB_pp, 4, 2, 2, "#it{b} jets per #it{N}_{Z}: PbPb / pp",      Form("bJetsPerZ_PbPbOverPP_%s", collTag()),      "no syst. yet");
  drawSet(rB_cp,  3, 3, 1, "#it{b} jets per #it{N}_{Z}: class / 50-80%",      Form("bJetsPerZ_CP_%s", collTag()),              "no syst. yet");
  if(haveIncl && havePP) drawSet(rBI_pp, 4, 2, 2, "#it{b} / inclusive: PbPb / pp", Form("bOverInclusive_PbPbOverPP_%s", collTag()), "#it{N}_{Z} cancels; no syst. yet");
  if(haveIncl) drawSet(rBI_cp, 3, 3, 1, "#it{b} / inclusive: class / 50-80%",   Form("bOverInclusive_CP_%s", collTag()),         "#it{N}_{Z} cancels; no syst. yet");

  // the same ratios as single panels with the classes overlaid, for talks:
  // one color and marker per class (as bPurity_vsJetPt), heavier lines and
  // larger text, drawn over ptReportMin-ptDrawMax; <name>_overlay.pdf
  auto drawOverlay = [&](TH1D **h, int n, const char *ytitle, const char *name, double yHi, const char *head){
    TCanvas *cv = new TCanvas(Form("%s_ov", name), "", 700, 650);
    cv->SetLeftMargin(0.15); cv->SetRightMargin(0.04); cv->SetTopMargin(0.05); cv->SetBottomMargin(0.14);
    TH1F *fr = cv->DrawFrame(ptReportMin, 0., ptDrawMax, yHi);
    fr->GetXaxis()->SetTitle("#it{p}_{T}^{jet} [GeV]"); fr->GetYaxis()->SetTitle(ytitle);
    fr->GetXaxis()->SetTitleSize(0.055); fr->GetYaxis()->SetTitleSize(0.055);
    fr->GetXaxis()->SetLabelSize(0.048); fr->GetYaxis()->SetLabelSize(0.048);
    fr->GetXaxis()->SetTitleOffset(1.1); fr->GetYaxis()->SetTitleOffset(1.25); fr->GetYaxis()->SetNdivisions(505);
    TLine *one = new TLine(ptReportMin, 1., ptDrawMax, 1.); one->SetLineStyle(7); one->SetLineWidth(2); one->SetLineColor(kGray+2); one->Draw();
    TLegend *lg = makeLegend(0.58, 0.90 - 0.075*n, 0.93, 0.90, 0.050);
    const int ovColor[4] = {6, 5, 3, 0};   // reddish purple, vermilion, blue, black: central to peripheral
    for(int p = 0; p < n; p++){
      TH1D *v = (TH1D*) h[p]->Clone(Form("%s_ov_%d", name, p)); v->SetDirectory(nullptr);
      for(int b = 1; b < NEdge; b++) if(!reported(b) || v->GetBinContent(b) <= 0.){ v->SetBinContent(b, -999.); v->SetBinError(b, 0.); }
      styleH(v, okabeHex[ovColor[p]], mk[p], mk[p] == markFilledDiamond ? 2.3 : 1.6); v->SetLineWidth(2); v->Draw("E1 X0 same");
      lg->AddEntry(v, coarseLabel[p], "pe");   // the y title says what the denominator is
    }
    lg->Draw();
    TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.050);
    la.DrawLatex(0.19, 0.88, head);
    la.SetTextSize(0.040); la.DrawLatex(0.19, 0.82, Form("%s, stat. only", collLabel()));
    savePdfTight(cv, Form("%s%s_overlay.pdf", figDir.Data(), name)); delete cv;
  };
  if(havePP) drawOverlay(rB_pp, 4, "#it{b} jets per #it{N}_{Z}: PbPb / pp", Form("bJetsPerZ_PbPbOverPP_%s", collTag()),      1.6, "#it{b} jets per #it{N}_{Z}");
  drawOverlay(rB_cp,  3, "#it{b} jets per #it{N}_{Z}: class / 50-80%", Form("bJetsPerZ_CP_%s", collTag()),              2.2, "#it{b} jets per #it{N}_{Z}");
  if(haveIncl && havePP) drawOverlay(rBI_pp, 4, "#it{b} / inclusive: PbPb / pp", Form("bOverInclusive_PbPbOverPP_%s", collTag()), 2.0, "#it{b} / inclusive");
  if(haveIncl) drawOverlay(rBI_cp, 3, "#it{b} / inclusive: class / 50-80%", Form("bOverInclusive_CP_%s", collTag()),         2.2, "#it{b} / inclusive");
  if(haveIncl){   // the fraction itself, classes overlaid
    TCanvas *cv = new TCanvas("cBI", "", 700, 650);
    cv->SetLeftMargin(0.15); cv->SetRightMargin(0.04); cv->SetTopMargin(0.05); cv->SetBottomMargin(0.14);
    double yHi = 0.; for(int ci = 0; ci < 4; ci++) for(int b = 1; b < NEdge; b++) if(reported(b) && edges[b-1] < ptDrawMax - 1e-6) yHi = TMath::Max(yHi, bOverIncl[ci]->GetBinContent(b) + bOverIncl[ci]->GetBinError(b));
    TH1F *fr = cv->DrawFrame(ptReportMin, 0., ptDrawMax, 1.4*yHi);
    fr->GetXaxis()->SetTitle("#it{p}_{T}^{jet} [GeV]"); fr->GetYaxis()->SetTitle("#it{b} jets / inclusive jets (per #it{N}_{Z})");
    fr->GetXaxis()->SetTitleSize(0.055); fr->GetYaxis()->SetTitleSize(0.055); fr->GetXaxis()->SetLabelSize(0.048); fr->GetYaxis()->SetLabelSize(0.048);
    fr->GetXaxis()->SetTitleOffset(1.1); fr->GetYaxis()->SetTitleOffset(1.25); fr->GetYaxis()->SetNdivisions(505);
    TLegend *lg = makeLegend(0.58, 0.60, 0.93, 0.90, 0.045);
    const int ovColor[4] = {6, 5, 3, 0};
    for(int ci = 0; ci < 4; ci++){
      TH1D *v = (TH1D*) bOverIncl[ci]->Clone(Form("bi_ov_%d", ci)); v->SetDirectory(nullptr);
      for(int b = 1; b < NEdge; b++) if(!reported(b) || v->GetBinContent(b) <= 0.){ v->SetBinContent(b, -999.); v->SetBinError(b, 0.); }
      styleH(v, okabeHex[ovColor[ci]], mk[ci], mk[ci] == markFilledDiamond ? 2.3 : 1.6); v->SetLineWidth(2); v->Draw("E1 X0 same");
      lg->AddEntry(v, coarseLabel[ci], "pe");
    }
    lg->Draw();
    TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.050); la.DrawLatex(0.19, 0.88, "#it{b} / inclusive");
    la.SetTextSize(0.040); la.DrawLatex(0.19, 0.82, Form("%s, stat. only", collLabel()));
    savePdfTight(cv, Form("%sbOverInclusive_%s.pdf", figDir.Data(), collTag())); delete cv;
  }
  { TCanvas *c = new TCanvas("cpur", "", 700, 650); c->SetLeftMargin(0.15); c->SetBottomMargin(0.13);
    TH1F *fr = new TH1F("fr_pur", "", 1, ptReportMin, 500.); fr->GetYaxis()->SetRangeUser(0., 1.); fr->GetXaxis()->SetTitle("jet #it{p}_{T} [GeV]"); fr->GetYaxis()->SetTitle("#it{b}-jet purity");
    fr->GetXaxis()->SetTitleSize(0.048); fr->GetYaxis()->SetTitleSize(0.048); fr->GetXaxis()->SetLabelSize(0.042); fr->GetYaxis()->SetLabelSize(0.042); fr->GetYaxis()->SetTitleOffset(1.4); fr->Draw("axis");
    TLegend *lg = makeLegend(0.55, 0.66, 0.92, 0.86, 0.038);
    for(int ci = 0; ci < 4; ci++){ TH1D *d = (TH1D*) purity[ci]->Clone(); d->SetDirectory(nullptr); for(int b = 1; b < NEdge; b++) if(!reported(b) || d->GetBinContent(b) <= 0){ d->SetBinContent(b, -9.); d->SetBinError(b, 0.); }
      styleH(d, okabeHex[ci == 0 ? 6 : (ci == 1 ? 5 : (ci == 2 ? 3 : 0))], mk[ci], 1.1); d->Draw("e1 same"); lg->AddEntry(d, Form("PbPb %s", coarseLabel[ci]), "pe"); }
    lg->Draw(); TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.040); t.DrawLatex(0.19, 0.86, Form("5.02 TeV, %s", collLabel())); t.SetTextSize(0.032); t.DrawLatex(0.19, 0.81, "2-template fit, b vs light+c, 0-5 GeV");
    savePdfTight(c, Form("%sbPurity_vsJetPt_%s.pdf", figDir.Data(), collTag())); delete c; }
  out->Close();
  printf("\nwrote %s and figures in %s\n", outRoot.Data(), figDir.Data());
}
