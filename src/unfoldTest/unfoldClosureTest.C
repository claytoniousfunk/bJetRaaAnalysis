//////////////////////////////////////////////////////////////////////////////////////
//
//  unfoldClosureTest.C
//
//  Split-sample closure test of the iterative Bayesian (D'Agostini) unfolding,
//  run for ONE centrality bin (or pp) at a time.
//
//    - the response matrix and the prior are built from one half of the MC
//      (even events by default)
//    - the other half (odd events) is unfolded and compared to its own
//      gen-level spectrum, so the test is statistically independent
//    - bias^2, variance, MSE and chi^2/ndf vs. N_iterations all end up on a
//      single summary plot, starting from N = 0 (the measured spectrum with no
//      unfolding applied), which is the reference the iterations improve on
//
//  Three figures per run, all in figures/unfoldTest/:
//    closure_<tag>.pdf          every iteration, for watching them walk away
//                               from the truth
//    closure_reduced_<tag>.pdf  gen truth, N = 0 and the min-MSE iteration only,
//                               with that iteration's uncertainty -- the version
//                               to show when the question is what the unfolding
//                               delivers rather than how it got there
//    closure_summary_<tag>.pdf  bias^2 / variance / MSE / chi^2 vs. N
//
//  The third argument selects which response is used to unfold the (always
//  pThat-weighted) spectrum; the fourth distorts the truth away from the prior so
//  that the optimal number of iterations is actually measurable.  See the comment
//  on the function itself.
//
//  usage:  root -l 'unfoldClosureTest.C("C1")'
//          root -l 'unfoldClosureTest.C("pp",12)'
//          root -l 'unfoldClosureTest.C("C1",10,"unweighted")'
//          root -l 'unfoldClosureTest.C("C1",10,"unweightedMatrix")'
//          root -l 'unfoldClosureTest.C("C1",15,"weighted","dataMC")'
//          root -l 'unfoldClosureTest.C("C1",15,"weighted","tilt:0.5")'
//          root -l 'unfoldClosureTest.C("C3",10,"weighted","none",100.,300.)'
//          root -l 'unfoldClosureTest.C("C1",10,"weighted","tilt:0.5",80.,300.,true)'   // calo jets
//          root -l 'unfoldClosureTest.C("pp",15,"weighted","none",80.,300.,true,4)'     // 20 GeV bins
//
//  REBINNING (last argument).  Merges the native 5 GeV bins by rebinFactor on
//  BOTH axes of both response matrices and in all four unmatched spectra, right
//  after loading, so every spectrum the test builds (all projections of those)
//  shares the coarser binning.  The factor must divide the native bin count.
//  Wider bins average down the statistical fluctuations that single high-weight
//  events put into individual bins -- the point of the option.
//  Since 2026-09-26 bias^2, variance and MSE are RELATIVE -- each bin divided by
//  its own truth before squaring -- so they no longer carry a counts^2 scale and
//  are broadly comparable between rebin factors, as chi2/ndf and max
//  |unfolded/truth - 1| always were.  MSE_abs in the table is the old counts^2
//  sum, kept only for comparison with runs made before that date.
//
//  FILE OVERRIDES (trainOverride, testOverride, outLabel).  Replace the training
//  file (response + prior) and/or the test file (unfolded spectrum + truth) with
//  explicit paths, for pairings the sample configs do not cover -- e.g. the full
//  pThat-weighted response unfolding a pThat-unweighted half.  The histogram
//  names are unchanged.  outLabel is appended to the output tag so the figures
//  do not overwrite the standard ones; with overrides it is required.
//
//  CALO JETS (last argument).  Uses the manual-JEC calo response scans' even and
//  odd halves from this repo -- the newest files matching
//    PYTHIA/PYTHIA_DiJet_response_caloJets_manualJEC*{even,odd}Events*.root
//    PYTHIAHYDJET/PYTHIAHYDJET_response_DiJet_caloJets_manualJEC*{even,odd}Events*.root
//  produced with onlyEvenEvents / onlyOddEvents in config_PYTHIA.h and
//  config_PYTHIAHYDJET.h.  Differences from PF, all forced by the calo samples:
//    - the response is built from allJets, NOT the sum of named flavours.  Calo
//      flavour comes from refparton_flavor, which leaves ~17-21% of real jets
//      unassigned (x); dropping them would discard a fifth of the sample, and
//      pp calo has no xJets histogram to subtract anyway.
//      The last argument, removeXJetsOpt, overrides that (-1 default, 0 keep,
//      1 drop) and tags the output _noX / _withX.  Tried 2026-09-25: it changes
//      nothing.  The ~20% is a whole-spectrum average sitting at gen pT 0-20;
//      above 80 GeV x is 0.02% of the C1 calo response, and dropping it moves
//      chi2/ndf at iteration 1 by ~1% (C1 11.74 -> 11.87, C4 7.99 -> 7.97) and
//      max|unfolded/truth - 1| by less than 0.1%.  x jets are not what breaks
//      the central-PbPb closure; the unfiltered fakes are.
//    - no pThat-unweighted calo production exists ("unweighted" modes refuse).
//    - "dataMC" is refused: its data files are the PF-era scans, and the
//      manual-JEC PbPb calo data use 5% classes that do not map onto C1..C4.
//      "none" and "tilt" work unchanged.
//    - outputs carry "_caloJets" so they never overwrite the PF results.
//
//////////////////////////////////////////////////////////////////////////////////////

#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include <iostream>

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#endif


//====================================================================================
//  sample configuration -- everything that depends on pp vs. centrality bin
//====================================================================================

struct SampleConfig {
  bool    ok    = false;
  bool    isPP  = false;
  bool    isCalo = false;
  TString label;        // short tag used in output file names
  TString title;        // pretty title drawn on the plots
  TString trainFile;    // even half, pThat-weighted   -> response + prior
  TString trainFileUnw; // even half, pThat-unweighted -> response + prior
                        // (empty when no such sample exists)
  TString testFile;     // odd half, pThat-weighted -> measured + truth, always
  TString responseName; // name of the 2D reco-vs-gen histogram

  // data, used only to derive the data/MC distortion.  Same files as the
  // original unfoldTest.C.
  TString dataMB, dataJet60, dataJet80, dataJet100;
  TString dataName;     // name of the 1D reco jet pT histogram
  double  dataMinPt;    // below this the stitched spectrum is not usable
  TString unmRecoName;  // reco jets with no gen match  -> fakes
  TString unmGenName;   // gen  jets with no reco match -> inefficiency
};

// newest file matching a shell glob, or "" if none.  `exclude`, when given, drops
// any match containing that substring -- needed because the PbPb calo data scans
// exist in both nominal and _ultraFineCentBins flavours under otherwise identical
// names, and the ultraFine ones carry 5% slices instead of C0..C4.
TString newestMatching(const TString &pattern, const TString &exclude = ""){
  TString cmd = Form("ls -1t %s 2>/dev/null", pattern.Data());
  if(!exclude.IsNull()) cmd += Form(" | grep -v '%s'", exclude.Data());
  cmd += " | head -1";
  TString s = gSystem->GetFromPipe(cmd);
  return s.Strip(TString::kBoth);
}

// Point an already-built config at the calo half-sample response scans. The
// histogram names are identical between the calo and PF response scans, so only
// the files, the title and the flags change.
void applyCaloConfig(SampleConfig &c){
  const TString dir  = "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/";
  const TString stem = c.isPP ? "PYTHIA/PYTHIA_DiJet_response_caloJets_manualJEC"
                              : "PYTHIAHYDJET/PYTHIAHYDJET_response_DiJet_caloJets_manualJEC";
  c.isCalo       = true;
  // pThat-WEIGHTED halves only ("_pThat-<number>"): the pThat-unweighted halves
  // share the stem and would otherwise be picked whenever they are newer
  c.trainFile    = newestMatching(dir + stem + "*evenEvents_pThat-[0-9]*.root");
  c.testFile     = newestMatching(dir + stem + "*oddEvents_pThat-[0-9]*.root");
  c.trainFileUnw = "";
  c.title       += ", calo jets";

  //  Calo data, for the dataMC distortion.  These replace the PF-era scans the
  //  config set up; without them "dataMC" would fit a PF data spectrum against a
  //  calo MC one and call the difference a prior error.
  //
  //  The PbPb globs end in [0-9].root via the `exclude`: every one of these
  //  scans also exists as *_ultraFineCentBins.root with 5% slices in place of
  //  C0..C4, and `ls -1t` would happily hand back whichever was written last.
  //
  //  There is no manual-JEC PbPb Jet60 scan, and none is needed: PbPb stitching
  //  sets jet60_pTmin = jet80_pTmin, so the Jet60 window closes to zero width and
  //  that histogram is never read as a source.  It is still dereferenced, so the
  //  Jet80 file stands in -- deliberately NOT the 2026-9-15 calo Jet60 scan,
  //  which has no manual JEC and would silently mix jet energy scales if the
  //  stitch logic ever changed.
  const TString dDir = "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/";
  if(c.isPP){
    c.dataMB     = newestMatching(dDir + "pp/pp_MinBias_caloJets_manualJEC_*.root");
    c.dataJet60  = newestMatching(dDir + "pp/pp_HighEGJet_caloJets_manualJEC_Jet60HLT_*.root");
    c.dataJet80  = newestMatching(dDir + "pp/pp_HighEGJet_caloJets_manualJEC_Jet80HLT_*.root");
    c.dataJet100 = newestMatching(dDir + "pp/pp_HighEGJet_caloJets_manualJEC_Jet100HLT_*.root");
  }
  else{
    const TString ex = "ultraFineCentBins";
    c.dataMB     = newestMatching(dDir + "PbPb/PbPb_MinBias_Part1_caloJets_manualJEC_*.root", ex);
    c.dataJet80  = newestMatching(dDir + "PbPb/PbPb_HardProbes_caloJets_manualJEC_Jet80HLT_*.root",  ex);
    c.dataJet100 = newestMatching(dDir + "PbPb/PbPb_HardProbes_caloJets_manualJEC_Jet100HLT_*.root", ex);
    c.dataJet60  = c.dataJet80;   // never read for PbPb; see above
  }

  if(c.trainFile.IsNull() || c.testFile.IsNull()){
    std::cout << "unfoldClosureTest: calo half-sample response scans not found:\n"
              << "   even: " << (c.trainFile.IsNull() ? TString("MISSING") : c.trainFile) << "\n"
              << "   odd : " << (c.testFile.IsNull()  ? TString("MISSING") : c.testFile)  << "\n"
              << "   looked for " << dir << stem << "*{even,odd}Events_pThat-<number>*.root\n"
              << "   Produce them with useCaloJetsOverride = true and onlyEvenEvents / onlyOddEvents\n"
              << "   in " << (c.isPP ? "config_PYTHIA.h (PYTHIA_scan_response.C)"
                                     : "config_PYTHIAHYDJET.h (PYTHIAHYDJET_scan_response.C)") << ".\n";
    c.ok = false;
    return;
  }
  std::cout << "calo response halves:\n   even: " << c.trainFile << "\n   odd : " << c.testFile << "\n";
}

SampleConfig getSampleConfig(TString sample){

  // pp: no pThat-unweighted production exists
  const TString ppDir = "../../rootFiles/scanningOuput/PYTHIA/";
  const TString ppFile = "PYTHIA_DiJet_response_pThat-15_mu12_pTmu-15_tight_jetTrkMaxFilter_doPThatCorrelationFilterTight_2026-8-3";

  // PbPb: same selection either way, only the pThat cross-section weight differs
  const TString PHDir = "/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/scanningOutput/PYTHIAHYDJET/latest/response/";
  const TString PHFile = "PYTHIAHYDJET_response_DiJet_pThat-15_mu12_pTmu-15_tight_vzReweight_hiBinReweight_hiBinShift-10_jetTrkMaxFilter_doPThatCorrelationFilterTight_2026-4-6";

  const TString PHDirUnw = "../../rootFiles/scanningOuput/PYTHIAHYDJET/";
  const TString PHFileUnw = "PYTHIAHYDJET_response_DiJet_pThat-unweighted_pThat-15_mu12_pTmu-15_tight_vzReweight_hiBinReweight_hiBinShift-10_jetTrkMaxFilter_doPThatCorrelationFilterTight_2026-7-28";

  // data, for the data/MC distortion -- the files used by the original unfoldTest.C
  const TString datDir = "/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/scanningOutput/";

  SampleConfig c;
  sample.ToLower();

  if(sample == "pp"){
    c.ok           = true;
    c.isPP         = true;
    c.label        = "pp";
    c.title        = "PYTHIA, pp";
    c.trainFile    = ppDir + ppFile + "_evenEvents.root";
    c.trainFileUnw = "";
    c.testFile     = ppDir + ppFile + "_oddEvents.root";
    c.responseName = "h_matchedRecoJetPt_genJetPt_allJets";
    c.dataMB       = datDir + "pp/latest/pp_MinBias_mu12_pTmu-14_tight_jetTrkMaxFilter_2025-10-15.root";
    c.dataJet60    = datDir + "pp/latest/pp_HighEGJet_Jet60HLT_mu12_pTmu-15to999_tight_deltaR-40_jetTrkMaxFilter_2026-3-10.root";
    c.dataJet80    = datDir + "pp/latest/pp_HighEGJet_Jet80HLT_mu12_pTmu-15to999_tight_deltaR-40_jetTrkMaxFilter_2026-3-10.root";
    c.dataJet100   = datDir + "pp/latest/pp_HighEGJet_Jet100HLT_mu12_pTmu-15to999_tight_deltaR-40_jetTrkMaxFilter_WDecayFilter_2026-3-11.root";
    c.dataName     = "h_inclRecoJetPt";
    c.dataMinPt    = 100.;   // MinBias carries the spectrum below this, poor stats
    c.unmRecoName  = "h_unmatchedRecoJetPt_allJets";
    c.unmGenName   = "h_unmatchedGenJetPt";
    return c;
  }

  TString centTitle = "";
  if     (sample == "c1") centTitle = "0-10%";
  else if(sample == "c2") centTitle = "10-30%";
  else if(sample == "c3") centTitle = "30-50%";
  else if(sample == "c4") centTitle = "50-80%";
  else{
    std::cout << "unfoldClosureTest: unknown sample \"" << sample << "\"\n"
              << "                   choose one of: pp, C1, C2, C3, C4\n";
    return c;
  }

  c.ok           = true;
  c.isPP         = false;
  c.label        = sample; c.label.ToUpper();
  c.title        = Form("PYTHIA+HYDJET %s", centTitle.Data());
  c.trainFile    = PHDir    + PHFile    + "_evenEvents.root";
  c.trainFileUnw = PHDirUnw + PHFileUnw + "_evenEvents.root";
  c.testFile     = PHDir    + PHFile    + "_oddEvents.root";
  c.responseName = Form("h_matchedRecoJetPt_genJetPt_allJets_%s", c.label.Data());
  c.dataMB       = datDir + "PbPb/latest/PbPb_MinBias_mu12_pTmu-15to999_tight_WDecayFilter_2026-4-8.root";
  c.dataJet60    = datDir + "PbPb/latest/PbPb_HardProbes_Jet60HLT_mu12_pTmu-15to999_tight_WDecayFilter_2026-3-11.root";
  c.dataJet80    = datDir + "PbPb/latest/PbPb_HardProbes_Jet80HLT_mu12_pTmu-15to999_tight_WDecayFilter_2026-4-8.root";
  c.dataJet100   = datDir + "PbPb/latest/PbPb_HardProbes_Jet100HLT_mu12_pTmu-15to999_tight_WDecayFilter_2026-4-8.root";
  c.dataName     = Form("h_inclRecoJetPt_%s", c.label.Data());
  // PbPb stitching sets jet60_pTmin = jet80_pTmin, so the MinBias normalisation
  // integral is an inverted range and the spectrum below 130 GeV is not usable
  c.dataMinPt    = 130.;
  c.unmRecoName  = Form("h_unmatchedRecoJetPt_allJets_%s", c.label.Data());
  c.unmGenName   = Form("h_unmatchedGenJetPt_%s",          c.label.Data());

  return c;

}


//====================================================================================
//  response loading, with the option to drop unassigned-flavour ("x") jets
//
//  x jets carry no parton flavour assignment; in the embedded PbPb sample they
//  are the background-origin category.  They are 5% of matched jets at 58 GeV in
//  C1, falling to <0.5% above 100 GeV, and are ordered by centrality
//  (C1 > C2 > C3 > C4), i.e. they behave like background rather than signal.
//
//  Rather than subtracting the xJets histogram, the response is rebuilt as the
//  sum of the five named parton flavours.  That works uniformly: in PbPb
//  allJets is exactly the sum of six categories including xJets, while pp has no
//  xJets histogram at all yet its five named flavours still fall 1.5% short of
//  allJets -- so summing removes the unassigned component in both cases.
//====================================================================================

TH2D* loadResponse(TFile *f, const TString &allName, bool removeX, const char *cloneName){

  TH2D *all = 0;
  f->GetObject(allName,all);
  if(!all) return 0;

  if(!removeX){
    TH2D *out = (TH2D*) all->Clone(cloneName);
    out->SetDirectory(0);
    return out;
  }

  const char *fl[5] = {"bJets","cJets","udJets","sJets","gJets"};
  TH2D *sum = 0;
  for(int i = 0; i < 5; i++){
    TString nm = allName; nm.ReplaceAll("allJets",fl[i]);
    TH2D *h = 0;
    f->GetObject(nm,h);
    if(!h){
      std::cout << "unfoldClosureTest: missing \"" << nm
                << "\", falling back to allJets (x jets NOT removed)\n";
      TH2D *out = (TH2D*) all->Clone(cloneName);
      out->SetDirectory(0);
      return out;
    }
    if(!sum){ sum = (TH2D*) h->Clone(cloneName); sum->SetDirectory(0); }
    else      sum->Add(h);
  }

  return sum;

}


//====================================================================================
//  projections
//
//  TH2::ProjectionX/Y default to firstbin=0, lastbin=-1, which INCLUDES the
//  under- and overflow of the other axis.  RooUnfold ignores over/underflow
//  (Overflow() is 0 by default), so the projections have to as well -- otherwise
//  the truth we compare against contains gen jets whose reco pT left the axis,
//  which the unfolded result can never reproduce.
//====================================================================================

TH1D* projX(TH2D *h, const char *name){
  return (TH1D*) h->ProjectionX(name,1,h->GetNbinsY());
}

TH1D* projY(TH2D *h, const char *name){
  return (TH1D*) h->ProjectionY(name,1,h->GetNbinsX());
}


//====================================================================================
//  manual variable-bin rebin, straight from the native 5 GeV response
//
//  Same edges as the scan's own "_var_" booking (PYTHIA_scan_response.C /
//  PYTHIAHYDJET_scan_response.C's N1/ptAxis1), but built here from the native
//  fixed-bin histogram instead of trusting the scan's separately-booked "_var_"
//  histogram -- a cross-check for whether that booking itself is behaving, and
//  the only option where a plain "_var_" histogram for a given
//  flavor/muTag combination was never booked at all.
//
//  Content below the lowest edge (60 GeV) has nowhere to go and is dropped --
//  FindBin returns 0 (underflow) for it, same as it would fall outside the
//  scan's own var-binned axis.
//====================================================================================

const int    kVarBinN = 10;
const double kVarBinEdges[kVarBinN] = {60,70,80,90,100,120,150,200,300,500};

TH2D* rebinToVarBins(TH2D *h, const char *name){
  TH2D *out = new TH2D(name,h->GetTitle(),kVarBinN-1,kVarBinEdges,kVarBinN-1,kVarBinEdges);
  out->SetDirectory(0);
  out->Sumw2();
  for(int ix = 1; ix <= h->GetNbinsX(); ix++){
    int jx = out->GetXaxis()->FindBin(h->GetXaxis()->GetBinCenter(ix));
    if(jx < 1 || jx > out->GetNbinsX()) continue;
    for(int iy = 1; iy <= h->GetNbinsY(); iy++){
      int jy = out->GetYaxis()->FindBin(h->GetYaxis()->GetBinCenter(iy));
      if(jy < 1 || jy > out->GetNbinsY()) continue;
      double c  = h->GetBinContent(ix,iy);
      double e  = h->GetBinError(ix,iy);
      double c0 = out->GetBinContent(jx,jy);
      double e0 = out->GetBinError(jx,jy);
      out->SetBinContent(jx,jy, c0 + c);
      out->SetBinError  (jx,jy, std::sqrt(e0*e0 + e*e));
    }
  }
  return out;
}


//  Okabe-Ito qualitative palette: eight hues chosen to stay distinguishable
//  under deuteranopia, protanopia and tritanopia.  Used for the summary curves,
//  which are categorical rather than ordered.
int cbBlue()   { return TColor::GetColor("#0072B2"); }
int cbVermil() { return TColor::GetColor("#D55E00"); }
int cbGreen()  { return TColor::GetColor("#009E73"); }

//====================================================================================
//  data spectrum -> distortion function
//
//  For the distorted-truth closure we need a realistic guess at how wrong the MC
//  prior might be.  We take it from the data: stitch the trigger samples the same
//  way the original unfoldTest.C does, divide by the MC reco spectrum, and fit the
//  shape with a power law.  Only the exponent is kept -- the overall data/MC
//  normalisation is meaningless here, and a constant weight cancels in the prior.
//====================================================================================

TH1D* stitchTriggers(TH1D *mb, TH1D *j60, TH1D *j80, TH1D *j100,
                     bool isPP, const char *name){

  TH1D *out = (TH1D*) j100->Clone(name);

  const double jet80_pTmin  = 130.;
  const double jet100_pTmin = 200.;
  const double e = 0.01;
  // for PbPb the jet60 sample is not used: its window closes to zero width
  const double jet60_pTmin = isPP ? 100. : jet80_pTmin;

  double N100 = j100->Integral(j100->FindBin(jet100_pTmin+e), j100->FindBin(500.-e));
  double N80  = j80 ->Integral(j80 ->FindBin(jet100_pTmin+e), j80 ->FindBin(500.-e));

  TH1D *s80 = (TH1D*) j80->Clone(Form("%s_s80",name));
  if(N80 > 0.) s80->Scale(N100/N80);

  double N80s = s80->Integral(s80->FindBin(jet80_pTmin+e), s80->FindBin(jet100_pTmin-e));
  double N60  = j60->Integral(j60 ->FindBin(jet80_pTmin+e), j60 ->FindBin(jet100_pTmin-e));

  TH1D *s60 = (TH1D*) j60->Clone(Form("%s_s60",name));
  if(N60 > 0.) s60->Scale(N80s/N60);

  double N60s = s60->Integral(s60->FindBin(jet60_pTmin+e), s60->FindBin(jet80_pTmin-e));
  double NMB  = mb ->Integral(mb ->FindBin(jet60_pTmin+e), mb ->FindBin(jet80_pTmin-e));

  TH1D *sMB = (TH1D*) mb->Clone(Form("%s_sMB",name));
  if(NMB > 0. && N60s > 0.) sMB->Scale(N60s/NMB);

  for(int i = 0; i < out->GetSize(); i++){
    double pt = out->GetBinCenter(i);
    TH1D *src = 0;
    if     (pt < jet60_pTmin)   src = sMB;
    else if(pt < jet80_pTmin)   src = s60;
    else if(pt < jet100_pTmin)  src = s80;
    else                        src = j100;
    out->SetBinContent(i,src->GetBinContent(i));
    out->SetBinError  (i,src->GetBinError(i));
  }

  return out;

}

TH1D* loadDataSpectrum(const SampleConfig &cfg){

  TFile *fMB  = TFile::Open(cfg.dataMB);
  TFile *f60  = TFile::Open(cfg.dataJet60);
  TFile *f80  = TFile::Open(cfg.dataJet80);
  TFile *f100 = TFile::Open(cfg.dataJet100);
  if(!fMB || !f60 || !f80 || !f100 ||
     fMB->IsZombie() || f60->IsZombie() || f80->IsZombie() || f100->IsZombie()){
    std::cout << "unfoldClosureTest: could not open the data files\n";
    return 0;
  }

  TH1D *hMB, *h60, *h80, *h100;
  fMB ->GetObject(cfg.dataName,hMB);
  f60 ->GetObject(cfg.dataName,h60);
  f80 ->GetObject(cfg.dataName,h80);
  f100->GetObject(cfg.dataName,h100);
  if(!hMB || !h60 || !h80 || !h100){
    std::cout << "unfoldClosureTest: could not find \"" << cfg.dataName << "\" in the data files\n";
    return 0;
  }

  TH1D *h = stitchTriggers(hMB,h60,h80,h100,cfg.isPP,"h_data_stitched");
  h->SetDirectory(0);
  return h;

}

//  Power-law fit to (data/MC), plus the data's relative statistical error mapped
//  onto the MC binning.  The data and response histograms do not share a binning
//  (pp data is 96 bins from 20 GeV, the response 100 from 0), so everything is
//  matched by bin centre and compared as densities rather than by bin index.
//  Show where the distortion comes from: the normalized data and MC reco
//  spectra on top, their ratio and the fitted power law underneath.  The fit is
//  drawn solid inside the fit range and dashed where it is extrapolated, since
//  the weight is applied over the whole gen axis but only constrained here.
void drawDistortionDerivation(TH1D *h_data, TH1D *h_mcReco, TGraphErrors *g_ratio, TF1 *f,
                              double fitLow, double fitHigh, double pTref,
                              double expo, double expoErr, double chi2ndf,
                              const char *titleStr, const char *outPath){

  const double xLo = 40., xHi = 400.;

  TCanvas *canv = new TCanvas("canv_dist","canv_dist",750,750);
  canv->cd();
  TPad *up = new TPad("pad_dist_up","",0,0.38,1,1);
  TPad *dn = new TPad("pad_dist_dn","",0,0,1,0.38);
  up->SetLeftMargin(0.15); dn->SetLeftMargin(0.15);
  up->SetRightMargin(0.05); dn->SetRightMargin(0.05);
  up->SetBottomMargin(0.);  dn->SetBottomMargin(0.30);
  up->SetTopMargin(0.14);   dn->SetTopMargin(0.);
  up->SetLogy(); up->SetTickx(1); up->SetTicky(1);
  dn->SetTickx(1); dn->SetTicky(1);
  up->Draw(); dn->Draw();

  //---- spectra ----------------------------------------------------------------
  up->cd();
  h_data->SetStats(0); h_data->SetTitle("");
  h_data->SetLineColor(kBlack); h_data->SetMarkerColor(kBlack);
  h_data->SetMarkerStyle(20); h_data->SetMarkerSize(0.8); h_data->SetLineWidth(2);
  h_data->GetXaxis()->SetRangeUser(xLo,xHi);
  h_data->GetYaxis()->SetTitle("normalized  (1/N) d#it{N}/d#it{p}_{T}");
  h_data->GetYaxis()->SetTitleSize(0.055); h_data->GetYaxis()->SetTitleOffset(1.25);
  h_data->GetYaxis()->SetLabelSize(0.045);
  h_data->Draw("p");
  h_mcReco->SetLineColor(cbBlue()); h_mcReco->SetLineWidth(3);
  h_mcReco->Draw("hist same");
  h_data->Draw("p same");

  TLegend *leg = new TLegend(0.58,0.62,0.93,0.84);
  leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.045);
  leg->AddEntry(h_data,"data (stitched)","lp");
  leg->AddEntry(h_mcReco,"MC reco (prior)","l");
  leg->Draw();

  TLatex lat; lat.SetNDC(); lat.SetTextSize(0.055);
  lat.DrawLatex(0.15,0.90,titleStr);

  //---- ratio and fit ----------------------------------------------------------
  dn->cd();
  TH1F *fr = dn->DrawFrame(xLo,0.,xHi,2.4);
  fr->GetXaxis()->SetTitle("Jet #it{p}_{T} [GeV]");
  fr->GetYaxis()->SetTitle("data / MC");
  fr->GetXaxis()->SetTitleSize(0.105); fr->GetXaxis()->SetTitleOffset(1.15);
  fr->GetYaxis()->SetTitleSize(0.095); fr->GetYaxis()->SetTitleOffset(0.70);
  fr->GetXaxis()->SetLabelSize(0.090); fr->GetYaxis()->SetLabelSize(0.090);
  fr->GetYaxis()->SetNdivisions(505);

  g_ratio->SetMarkerStyle(20); g_ratio->SetMarkerSize(0.8);
  g_ratio->SetMarkerColor(kBlack); g_ratio->SetLineColor(kBlack);
  g_ratio->Draw("P same");

  // extrapolated part, dashed
  TF1 *fx = (TF1*) f->Clone("f_dataMC_extrap");
  fx->SetRange(xLo,xHi);
  fx->SetLineColor(cbVermil()); fx->SetLineStyle(2); fx->SetLineWidth(2);
  fx->Draw("l same");
  // constrained part, solid
  TF1 *fi = (TF1*) f->Clone("f_dataMC_infit");
  fi->SetRange(fitLow,fitHigh);
  fi->SetLineColor(cbVermil()); fi->SetLineStyle(1); fi->SetLineWidth(3);
  fi->Draw("l same");

  TLine *li = new TLine(); li->SetLineStyle(7); li->SetLineColor(kGray+2);
  li->DrawLine(xLo,1.,xHi,1.);
  li->DrawLine(fitLow,0.,fitLow,2.4);
  li->DrawLine(fitHigh,0.,fitHigh,2.4);

  TLatex l2; l2.SetNDC(); l2.SetTextSize(0.075);
  l2.DrawLatex(0.19,0.88,Form("#it{w} #propto (#it{p}_{T}/%.0f)^{%.3f #pm %.3f}",pTref,expo,expoErr));
  l2.SetTextSize(0.062); l2.SetTextColor(kGray+2);
  l2.DrawLatex(0.19,0.78,Form("fit %.0f-%.0f GeV,  #chi^{2}/ndf = %.2f",fitLow,fitHigh,chi2ndf));

  canv->SaveAs(outPath);

}

bool buildDataMCDistortion(TH1D *h_data, TH1D *h_mcReco,
                           double fitLow, double fitHigh, double pTref, double dataMinPt,
                           double &expo, std::vector<double> &relErr,
                           const char *titleStr = 0, const char *outPath = 0){

  int b1 = h_data ->FindBin(fitLow), b2 = h_data ->FindBin(fitHigh);
  int m1 = h_mcReco->FindBin(fitLow), m2 = h_mcReco->FindBin(fitHigh);
  double sD = h_data ->Integral(b1,b2);
  double sM = h_mcReco->Integral(m1,m2);
  if(sD <= 0. || sM <= 0.){
    std::cout << "unfoldClosureTest: empty data or MC spectrum over the fit range\n";
    return false;
  }

  TH1D *r = (TH1D*) h_mcReco->Clone("h_dataMC_ratio");   // MC binning
  r->SetDirectory(0);
  r->Reset();

  //  normalized copies, kept only for the derivation plot
  TH1D *h_dataN = (TH1D*) h_data  ->Clone("h_data_norm");   h_dataN->SetDirectory(0);
  TH1D *h_mcN   = (TH1D*) h_mcReco->Clone("h_mcReco_norm"); h_mcN  ->SetDirectory(0);
  for(int b = 1; b <= h_dataN->GetNbinsX(); b++){
    double w = h_dataN->GetBinWidth(b);
    h_dataN->SetBinContent(b,h_dataN->GetBinContent(b)/sD/w);
    h_dataN->SetBinError  (b,h_dataN->GetBinError(b)  /sD/w);
  }
  for(int b = 1; b <= h_mcN->GetNbinsX(); b++){
    double w = h_mcN->GetBinWidth(b);
    h_mcN->SetBinContent(b,h_mcN->GetBinContent(b)/sM/w);
    h_mcN->SetBinError  (b,h_mcN->GetBinError(b)  /sM/w);
  }
  TGraphErrors *g_ratio = new TGraphErrors();

  relErr.assign(h_mcReco->GetNbinsX()+2, 0.);
  double relRef = 0.;

  for(int b = 1; b <= h_mcReco->GetNbinsX(); b++){

    double x  = h_mcReco->GetBinCenter(b);
    int    db = h_data->FindBin(x);
    if(db < 1 || db > h_data->GetNbinsX()) continue;

    double cD = h_data ->GetBinContent(db), eD = h_data ->GetBinError(db);
    double cM = h_mcReco->GetBinContent(b), eM = h_mcReco->GetBinError(b);
    if(cD <= 0. || cM <= 0.) continue;

    double rel = eD/cD;
    if(x >= dataMinPt){ relErr[b] = rel; if(relRef == 0.) relRef = rel; }

    // densities, so unequal bin widths do not bias the shape
    double dD = (cD/sD)/h_data ->GetBinWidth(db);
    double dM = (cM/sM)/h_mcReco->GetBinWidth(b);
    double ratio = dD/dM;

    double ratioErr = ratio*std::sqrt(rel*rel + (eM/cM)*(eM/cM));
    r->SetBinContent(b,ratio);
    r->SetBinError  (b,ratioErr);

    int ip = g_ratio->GetN();
    g_ratio->SetPoint(ip,x,ratio);
    g_ratio->SetPointError(ip,0.,ratioErr);

  }

  // bins with no usable data (below dataMinPt, or off the data axis) inherit the
  // precision of the lowest usable bin
  for(int b = 1; b <= h_mcReco->GetNbinsX(); b++)
    if(relErr[b] == 0.) relErr[b] = relRef;

  TF1 *f = new TF1("f_dataMC","[0]*pow(x/[2],[1])",fitLow,fitHigh);
  f->FixParameter(2,pTref);
  f->SetParameters(1.,0.3);
  TFitResultPtr res = r->Fit(f,"QNRS");
  if(!res.Get() || !res->IsValid()){
    std::cout << "unfoldClosureTest: data/MC power-law fit failed\n";
    return false;
  }

  expo = f->GetParameter(1);
  double chi2ndf = (f->GetNDF() > 0) ? f->GetChisquare()/f->GetNDF() : 0.;
  printf("  data/MC power-law fit over %.0f-%.0f GeV: exponent = %.3f +/- %.3f  (chi2/ndf = %.2f)\n",
         fitLow,fitHigh,expo,f->GetParError(1),chi2ndf);

  if(titleStr && outPath)
    drawDistortionDerivation(h_dataN,h_mcN,g_ratio,f,fitLow,fitHigh,pTref,
                             expo,f->GetParError(1),chi2ndf,titleStr,outPath);

  return true;

}

// w(pT) = (pT/pTref)^expo, frozen outside [40,500] so the extrapolation cannot
// drive a gen row to zero or blow it up
double distortionWeight(double pt, double expo, double pTref){
  double x = pt;
  if(x <  40.) x =  40.;
  if(x > 500.) x = 500.;
  return std::pow(x/pTref, expo);
}

//  Spectrum panels only: with variable bins, counts per bin make a wide bin read
//  as a tall one, so divide by width for display.  Metrics and ratio panels keep
//  using the unscaled histograms (counts), where the width cancels anyway.
bool hasVariableBins(const TH1 *h){ return h->GetXaxis()->GetXbins()->GetSize() > 0; }

TH1D* displayClone(const TH1D *h, const char *name){
  TH1D *c = (TH1D*) h->Clone(name);
  c->SetDirectory(0);
  if(hasVariableBins(h)) c->Scale(1.,"width");
  return c;
}

const char* spectrumYTitle(const TH1 *h){
  return hasVariableBins(h) ? "counts / GeV" : "counts";
}

//  Original vs. distorted spectrum, with a ratio panel underneath. Used for
//  both the reco-level input to the unfolding (h_meas*) and the gen-level
//  truth the distortion is actually defined on (h_truth*) -- see
//  drawTiltDerivation below, which calls this once for each.
void drawSpectrumPair(TH1D *hOrig, TH1D *hDist, const char *xTitle,
                      const char *titleStr, const char *outPath){

  const double xLo = 40., xHi = 400.;

  TCanvas *canv = new TCanvas("canv_tilt","canv_tilt",750,750);
  canv->cd();
  TPad *up = new TPad("pad_tilt_up","",0,0.38,1,1);
  TPad *dn = new TPad("pad_tilt_dn","",0,0,1,0.38);
  up->SetLeftMargin(0.15); dn->SetLeftMargin(0.15);
  up->SetRightMargin(0.05); dn->SetRightMargin(0.05);
  up->SetBottomMargin(0.);  dn->SetBottomMargin(0.30);
  up->SetTopMargin(0.10);   dn->SetTopMargin(0.);
  up->SetLogy(); up->SetTickx(1); up->SetTicky(1);
  dn->SetTickx(1); dn->SetTicky(1);
  up->Draw(); dn->Draw();

  //---- spectrum, original and distorted ------------------------------------------
  up->cd();
  TH1D *mO = displayClone(hOrig,"tilt_mO");
  TH1D *mD = displayClone(hDist,"tilt_mD");
  mO->SetStats(0); mO->SetTitle("");
  mO->GetXaxis()->SetRangeUser(xLo,xHi);
  mO->GetYaxis()->SetTitle(spectrumYTitle(hOrig));
  mO->GetYaxis()->SetTitleSize(0.055); mO->GetYaxis()->SetTitleOffset(1.25);
  mO->GetYaxis()->SetLabelSize(0.045);
  mO->SetLineColor(kBlack);    mO->SetLineWidth(3);
  mD->SetLineColor(cbVermil()); mD->SetLineWidth(3);
  mO->Draw("hist");
  mD->Draw("hist same");

  TLegend *leg = new TLegend(0.55,0.70,0.93,0.86);
  leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.045);
  leg->AddEntry(mO,"original","l");
  leg->AddEntry(mD,"distorted","l");
  leg->Draw();

  TLatex lat; lat.SetNDC(); lat.SetTextSize(0.050);
  lat.DrawLatex(0.15,0.925,titleStr);

  //---- ratio --------------------------------------------------------------------
  dn->cd();
  TH1D *rM = (TH1D*) mD->Clone("tilt_rM"); rM->SetDirectory(0); rM->Divide(mO);
  double yLo = 1e9, yHi = -1e9;
  for(int b = rM->FindBin(xLo); b <= rM->FindBin(xHi - 1e-6); b++){
    if(rM->GetBinContent(b) <= 0.) continue;
    yLo = std::min(yLo, rM->GetBinContent(b)); yHi = std::max(yHi, rM->GetBinContent(b));
  }
  yLo = std::min(yLo, 1.) - 0.1; yHi = std::max(yHi, 1.) + 0.1;
  TH1F *fr = dn->DrawFrame(xLo,yLo,xHi,yHi);
  fr->GetXaxis()->SetTitle(xTitle);
  fr->GetYaxis()->SetTitle("distorted / original");
  fr->GetXaxis()->SetTitleSize(0.105); fr->GetXaxis()->SetTitleOffset(1.15);
  fr->GetYaxis()->SetTitleSize(0.085); fr->GetYaxis()->SetTitleOffset(0.78);
  fr->GetXaxis()->SetLabelSize(0.090); fr->GetYaxis()->SetLabelSize(0.090);
  fr->GetYaxis()->SetNdivisions(505);

  rM->SetLineColor(cbVermil()); rM->SetLineWidth(3);
  rM->Draw("hist same");

  TLine *li = new TLine(); li->SetLineStyle(7); li->SetLineColor(kGray+2);
  li->DrawLine(xLo,1.,xHi,1.);

  canv->SaveAs(outPath);

}

//  How much the tilt distorts (a) the reco spectrum that actually gets
//  unfolded and (b) the gen-level truth the distortion is defined on -- two
//  separate figures, since the reco-level shift is smeared by the response
//  and is not the same curve as w(genPt) itself.
void drawTiltDerivation(TH1D *h_truthOrig, TH1D *h_truthDist,
                        TH1D *h_measOrig,  TH1D *h_measDist,
                        double expo, double pTref, double ptLow, double ptHigh,
                        const char *titleStr, const char *outPath){

  drawSpectrumPair(h_measOrig,h_measDist,"Reco jet #it{p}_{T} [GeV]",titleStr,outPath);

  TString truthPath(outPath);
  truthPath.ReplaceAll(".pdf","_genTruth.pdf");
  drawSpectrumPair(h_truthOrig,h_truthDist,"Gen jet #it{p}_{T} [GeV]",titleStr,truthPath);

}


//====================================================================================
//  closure metrics
//====================================================================================

//  RELATIVE METRICS, adopted 2026-09-26.  bias2/variance/mse are sums of
//  squared COUNTS, so they are dominated by whichever bin holds the most yield
//  -- on a steeply falling spectrum, always the lowest one in the window.  In
//  C3 the 80-100 GeV bin alone supplied 85-99% of bias2 at nearly every
//  iteration (mean share 0.90; C4 0.83), and the "minimum MSE at N = 9" it
//  produced was simply where that one bin's deviation crossed zero on its way
//  from positive to negative -- the other ten bins stayed inside 1.5% the whole
//  time and max|unfolded/truth - 1| had no minimum there at all.  On a log axis
//  the crossing looked like a spike in the bias curve.
//
//  Dividing each bin by its own truth removes the yield weighting, so every bin
//  in the window counts equally.  The minima become broad and consistent
//  (N ~ 8-11 in all four classes, where the absolute metric gave 3 to 10), and
//  the spurious structure disappears.  The absolute MSE is kept in the table for
//  continuity with earlier runs.
struct ClosureMetrics {
  double chi2     = 0.;   // sum (unfolded - truth)^2 / sigma^2
  int    ndf      = 0;    // number of bins entering the sum
  double bias2    = 0.;   // sum ((unfolded - truth)/truth)^2
  double variance = 0.;   // sum (sigma/truth)^2
  double mse      = 0.;   // bias^2 + variance
  double bias2Abs = 0.;   // sum (unfolded - truth)^2, counts^2
  double varAbs   = 0.;   // sum sigma^2, counts^2
  double mseAbs   = 0.;   // bias2Abs + varAbs
  double maxDev   = 0.;   // max |unfolded/truth - 1|: bin-width independent
};

ClosureMetrics computeMetrics(TH1D *h_unfold, TH1D *h_truth, double ptLow, double ptHigh){

  ClosureMetrics m;

  for(int i = 1; i <= h_unfold->GetNbinsX(); i++){

    double pt = h_unfold->GetBinCenter(i);
    if(pt < ptLow || pt > ptHigh) continue;

    double u   = h_unfold->GetBinContent(i);
    double t   = h_truth ->GetBinContent(i);
    double sig = h_unfold->GetBinError(i);

    if(t == 0.) continue;

    const double d = (u - t) / t;      // relative deviation
    const double s = sig / t;          // relative uncertainty

    m.bias2    += d * d;
    m.variance += s * s;
    m.bias2Abs += (u - t) * (u - t);
    m.varAbs   += sig * sig;
    m.maxDev    = std::max(m.maxDev, std::fabs(d));

    if(sig > 0.){
      m.chi2 += (u - t) * (u - t) / (sig * sig);   // scale-free already
      m.ndf++;
    }

  }

  m.mse    = m.bias2    + m.variance;
  m.mseAbs = m.bias2Abs + m.varAbs;

  return m;

}


//====================================================================================
//  plotting helpers
//====================================================================================

//  The iterations are an ordered series, so they get a sequential palette.
//  kCividis is perceptually uniform and built for deuteranopia; the top of it is
//  a pale yellow that washes out on white, so only the lower 85% is used.
int iterColor(int iter, int N_iter_max){
  // N = 0 (the measured spectrum, no unfolding applied) is kept grey and dashed
  if(iter <= 0) return kGray+2;
  if(N_iter_max <= 1) return TColor::GetColorPalette(0);
  const double frac = 0.85;
  return TColor::GetColorPalette(int((iter-1) * frac * (TColor::GetNumberOfColors() - 1)
                                     / double(N_iter_max - 1)));
}

const char* iterLabel(int iter){
  return (iter == 0) ? "#it{N} = 0 (meas.)" : Form("#it{N} = %i",iter);
}

// spectra (top) and unfolded/truth ratio (bottom) for N = 0 ... N_iter_max
void drawClosure(int N_iter_max, TH1D **h_unfold, TH1D *h_truth,
                 double drawPtLow, double drawPtHigh,
                 double ptLow, double ptHigh,
                 double ratioLow, double ratioHigh,
                 const char *titleStr, const char *subTitleStr,
                 const char *subTitleStr2, const char *outPath){

  TLine *li = new TLine();
  li->SetLineStyle(7);

  TCanvas *canv = new TCanvas("canv_closure","canv_closure",750,750);
  canv->cd();
  TPad *pad_upper = new TPad("pad_closure_upper","",0,0.35,1,1);
  TPad *pad_lower = new TPad("pad_closure_lower","",0,0,1,0.35);
  pad_upper->SetLeftMargin(0.15);  pad_lower->SetLeftMargin(0.15);
  pad_upper->SetRightMargin(0.05); pad_lower->SetRightMargin(0.05);
  pad_upper->SetBottomMargin(0.);  pad_lower->SetBottomMargin(0.30);
  pad_upper->SetTopMargin(0.16);   pad_lower->SetTopMargin(0.);
  pad_upper->SetLogy();
  pad_upper->SetTickx(1); pad_upper->SetTicky(1);
  pad_lower->SetTickx(1); pad_lower->SetTicky(1);
  pad_upper->Draw(); pad_lower->Draw();

  //  entries = gen truth + N = 0 ... N_iter_max.  A long scan overruns a
  //  two-column legend, so widen and shrink it rather than let the rows overlap.
  const int nLegEntries = N_iter_max + 2;
  TLegend *leg = new TLegend(0.60,0.32,0.93,0.84);
  leg->SetBorderSize(0); leg->SetFillStyle(0);
  leg->SetTextSize(nLegEntries > 14 ? 0.024 : 0.036);
  leg->SetNColumns(nLegEntries > 14 ? 3 : 2);

  //---- spectra ----------------------------------------------------------------
  pad_upper->cd();

  TH1D *td = displayClone(h_truth,"disp_closure_truth");
  td->SetStats(0);
  td->SetTitle("");
  td->SetLineColor(kBlack); td->SetLineWidth(3);
  td->GetXaxis()->SetRangeUser(drawPtLow,drawPtHigh);
  td->GetYaxis()->SetTitle(spectrumYTitle(h_truth));
  td->GetYaxis()->SetTitleSize(0.055); td->GetYaxis()->SetTitleOffset(1.3);
  td->GetYaxis()->SetLabelSize(0.045);
  td->Draw("hist");

  std::vector<TH1D*> ud(N_iter_max+1);
  for(int n = 0; n <= N_iter_max; n++){
    ud[n] = displayClone(h_unfold[n],Form("disp_closure_iter%i",n));
    ud[n]->Draw("hist same");
  }
  td->Draw("hist same");

  leg->AddEntry(td,"gen truth","l");
  for(int n = 0; n <= N_iter_max; n++) leg->AddEntry(ud[n],iterLabel(n),"l");
  leg->Draw();

  TLatex lat; lat.SetNDC(); lat.SetTextSize(0.050);
  lat.DrawLatex(0.15,0.945,titleStr);
  lat.SetTextSize(0.030); lat.SetTextColor(kGray+2);
  lat.DrawLatex(0.15,0.885,
                (subTitleStr2 && subTitleStr2[0])
                ? Form("%s   |   %s",subTitleStr,subTitleStr2) : subTitleStr);

  //---- ratio to truth ---------------------------------------------------------
  pad_lower->cd();

  for(int n = 0; n <= N_iter_max; n++){

    TH1D *r = (TH1D*) h_unfold[n]->Clone(Form("r_closure_iter%i",n));
    //  binomial for the unfolded iterations, which share their events with the
    //  truth (see the note in drawClosureReduced); plain for N = 0, whose ratio
    //  to truth runs to 8-14 before the fakes are unfolded away.  Neither is
    //  drawn on this figure -- ten curves leave no room for error bars -- but
    //  the two figures should not disagree about what the ratio means.
    r->Divide(h_unfold[n],h_truth,1,1, n > 0 ? "B" : "");

    if(n > 0){ r->Draw("hist same"); continue; }

    r->SetStats(0);
    r->SetTitle("");
    r->GetXaxis()->SetRangeUser(drawPtLow,drawPtHigh);
    r->GetYaxis()->SetRangeUser(ratioLow,ratioHigh);
    r->GetXaxis()->SetTitle("Jet #it{p}_{T} [GeV]");
    r->GetYaxis()->SetTitle("unfolded / truth");
    r->GetXaxis()->SetTitleSize(0.100); r->GetXaxis()->SetTitleOffset(1.15);
    r->GetYaxis()->SetTitleSize(0.090); r->GetYaxis()->SetTitleOffset(0.75);
    r->GetXaxis()->SetLabelSize(0.085);
    r->GetYaxis()->SetLabelSize(0.085);
    r->GetYaxis()->SetNdivisions(505);
    r->Draw("hist");

  }

  li->DrawLine(drawPtLow,1.00,drawPtHigh,1.00);
  li->DrawLine(drawPtLow,0.95,drawPtHigh,0.95);
  li->DrawLine(drawPtLow,1.05,drawPtHigh,1.05);

  // window used for the bias/variance/chi2 metrics
  li->SetLineColor(kGray+2);
  li->DrawLine(ptLow, ratioLow,ptLow, ratioHigh);
  li->DrawLine(ptHigh,ratioLow,ptHigh,ratioHigh);
  li->SetLineColor(kBlack);

  canv->SaveAs(outPath);

}

// Same figure reduced to three curves: gen truth, N = 0 (the measured spectrum,
// no unfolding) and the one iteration the MSE scan picks out.  The full version
// stacks eleven curves, which is the right plot for watching the iterations walk
// away from the truth but the wrong one for showing what the unfolding actually
// delivers -- there the two references and the answer are all that matter.
//
// The optimal N is only meaningful for a distorted-truth run.  In a plain
// split-sample test the prior IS the truth up to the even/odd split, so MSE is
// minimised at N = 1 by construction; nOptMeaningful says which case this is and
// the legend is labelled accordingly.
void drawClosureReduced(TH1D *h_meas, TH1D *h_opt, int nOpt, bool nOptMeaningful,
                        TH1D *h_truth,
                        double drawPtLow, double drawPtHigh,
                        double ptLow, double ptHigh,
                        double ratioLow, double ratioHigh,
                        const char *titleStr,
                        const std::vector<TString> &infoLines,
                        const char *outPath){

  TLine *li = new TLine();
  li->SetLineStyle(7);

  TCanvas *canv = new TCanvas("canv_closure_red","canv_closure_red",750,750);
  canv->cd();
  TPad *pad_upper = new TPad("pad_red_upper","",0,0.35,1,1);
  TPad *pad_lower = new TPad("pad_red_lower","",0,0,1,0.35);
  pad_upper->SetLeftMargin(0.15);  pad_lower->SetLeftMargin(0.15);
  pad_upper->SetRightMargin(0.05); pad_lower->SetRightMargin(0.05);
  pad_upper->SetBottomMargin(0.);  pad_lower->SetBottomMargin(0.30);
  pad_upper->SetTopMargin(0.16);   pad_lower->SetTopMargin(0.);
  pad_upper->SetLogy();
  pad_upper->SetTickx(1); pad_upper->SetTicky(1);
  pad_lower->SetTickx(1); pad_lower->SetTicky(1);
  pad_upper->Draw(); pad_lower->Draw();

  // clones, so restyling here cannot disturb the full eleven-curve figure
  TH1D *t = (TH1D*) h_truth->Clone("red_truth");
  TH1D *m = (TH1D*) h_meas ->Clone("red_meas");
  TH1D *u = (TH1D*) h_opt  ->Clone("red_opt");
  t->SetDirectory(0); m->SetDirectory(0); u->SetDirectory(0);

  t->SetLineColor(kBlack);    t->SetLineWidth(3); t->SetLineStyle(1);
  m->SetLineColor(kGray+2);   m->SetLineWidth(2); m->SetLineStyle(7);
  u->SetLineColor(cbVermil()); u->SetLineWidth(3); u->SetLineStyle(1);
  u->SetMarkerColor(cbVermil());

  //---- spectra ----------------------------------------------------------------
  pad_upper->cd();

  TH1D *td = displayClone(t,"red_disp_truth");
  TH1D *md = displayClone(m,"red_disp_meas");
  TH1D *ud = displayClone(u,"red_disp_opt");
  td->SetStats(0);
  td->SetTitle("");
  td->GetXaxis()->SetRangeUser(drawPtLow,drawPtHigh);
  td->GetYaxis()->SetTitle(spectrumYTitle(t));
  td->GetYaxis()->SetTitleSize(0.055); td->GetYaxis()->SetTitleOffset(1.3);
  td->GetYaxis()->SetLabelSize(0.045);
  td->Draw("hist");
  md->Draw("hist same");
  ud->Draw("hist same");
  td->Draw("hist same");   // back on top: it is the reference being judged

  TLegend *leg = new TLegend(0.55,0.60,0.93,0.84);
  leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.040);
  leg->AddEntry(td,"gen truth","l");
  leg->AddEntry(md,"#it{N} = 0 (meas.)","l");
  leg->AddEntry(ud, nOptMeaningful ? Form("#it{N} = %i (min MSE)",nOpt)
                                   : Form("#it{N} = %i",nOpt), "l");
  leg->Draw();

  TLatex lat; lat.SetNDC(); lat.SetTextSize(0.050);
  lat.DrawLatex(0.15,0.945,titleStr);

  //  Provenance goes bottom left, not under the title: the spectra fall by five
  //  decades across the pad, so that corner is the only empty one, and a reader
  //  checking which half produced which curve is looking at the plot, not at
  //  the header.
  lat.SetTextSize(0.038); lat.SetTextColor(kBlack); lat.SetTextFont(42);
  const double lineStep = 0.052;
  for(size_t i = 0; i < infoLines.size(); i++)
    lat.DrawLatex(0.19, 0.075 + lineStep*(infoLines.size()-1-i), infoLines[i].Data());

  //---- ratio to truth ---------------------------------------------------------
  pad_lower->cd();

  TH1D *rm = (TH1D*) m->Clone("red_r_meas");
  TH1D *ru = (TH1D*) u->Clone("red_r_opt");
  rm->SetDirectory(0); ru->SetDirectory(0);
  //  "B" -- binomial errors on unfolded/truth.  The unfolded spectrum and the
  //  truth it is compared against are built from the SAME events (the odd half,
  //  or the training half for a distorted run), so treating them as independent
  //  and adding the two relative errors in quadrature overstates the spread.
  //
  //  Caveat, since it matters if these bars are ever quoted: the binomial form
  //  assumes the numerator is a subset of the denominator, and an unfolded
  //  spectrum is not -- it is a different estimator on the same events, and the
  //  ratio sits above 1 in plenty of bins.  ROOT's weighted-binomial expression
  //  |((1-2w)e1^2 + w^2 e2^2)/(c2 b2)^2| goes outside its domain there and is
  //  saved only by the absolute value.  In practice it makes almost no
  //  difference: across 80-300 GeV in C1 the unfolded error dominates and the
  //  binomial and uncorrelated bars agree to better than 15% (worst bin
  //  160-180, 0.057 vs 0.066), so this is a presentational choice, not a
  //  correction.
  //
  //  rm stays uncorrelated -- measured/truth runs to 8-14 before unfolding,
  //  where a binomial error would be meaningless.  Its bars are not drawn.
  rm->Divide(m,t,1,1,"");
  ru->Divide(u,t,1,1,"B");
  // Bars, not a filled band.  Outside the metric window the relative error runs
  // to ~100% (the 40-60 bin, and everything above ~220 GeV where the C1 yield is
  // a fraction of a count), and a shaded band there covers the whole pad and
  // reads as "nothing is known".  Bars carry the same numbers per bin and stay
  // out of the way.
  ru->SetMarkerStyle(kFullCircle);
  ru->SetMarkerSize(0.7);

  rm->SetStats(0);
  rm->SetTitle("");
  rm->GetXaxis()->SetRangeUser(drawPtLow,drawPtHigh);
  rm->GetYaxis()->SetRangeUser(ratioLow,ratioHigh);
  rm->GetXaxis()->SetTitle("Jet #it{p}_{T} [GeV]");
  rm->GetYaxis()->SetTitle("unfolded / truth");
  rm->GetXaxis()->SetTitleSize(0.100); rm->GetXaxis()->SetTitleOffset(1.15);
  rm->GetYaxis()->SetTitleSize(0.090); rm->GetYaxis()->SetTitleOffset(0.75);
  rm->GetXaxis()->SetLabelSize(0.085);
  rm->GetYaxis()->SetLabelSize(0.085);
  rm->GetYaxis()->SetNdivisions(505);
  rm->Draw("hist");
  // with one curve instead of ten there is room for its uncertainty, which is
  // the whole point of picking an iteration: it trades bias for variance
  ru->Draw("hist same");
  ru->Draw("e1 same");

  // unity only -- no +-5% guides.  This figure carries the iteration's own
  // uncertainty, and fixed guide lines invite reading the closure against 5%
  // rather than against the error bars, which are what actually bound it.
  li->DrawLine(drawPtLow,1.00,drawPtHigh,1.00);

  li->SetLineColor(kGray+2);
  li->DrawLine(ptLow, ratioLow,ptLow, ratioHigh);
  li->DrawLine(ptHigh,ratioLow,ptHigh,ratioHigh);
  li->SetLineColor(kBlack);

  canv->SaveAs(outPath);

}

// bias^2 / variance / MSE (left, log) and chi^2/ndf (right, log) on one canvas.
// The arrays hold N_iter_max+1 points, starting at N = 0 (no unfolding).
void drawSummary(int N_iter_max, double *x, double *bias2, double *var, double *mse,
                 double *chi2ndf, const char *titleStr, const char *subTitleStr,
                 const char *subTitleStr2, bool nOptMeaningful, const char *outPath){

  const int    N     = N_iter_max + 1;
  const double xLow  = -0.5;
  const double xHigh = N_iter_max + 0.5;
  const double dataFrac = 0.75;   // curves occupy the lower 75% of the frame,
                                  // the rest is headroom for the legend

  // ---- left axis range: only positive values survive the log scale ----------
  double yLow = 1e30, yHigh = -1e30;
  for(int i = 0; i < N; i++){
    for(double v : {bias2[i], var[i], mse[i]}){
      if(v <= 0.) continue;
      if(v < yLow)  yLow  = v;
      if(v > yHigh) yHigh = v;
    }
  }
  if(yLow > yHigh){ yLow = 1e-3; yHigh = 1.; }
  yLow  *= 0.5;
  yHigh *= 2.;
  yHigh *= std::pow(10., std::log10(yHigh/yLow) * (1./dataFrac - 1.));

  TCanvas *canv = new TCanvas("canv_summary","canv_summary",750,650);
  canv->cd();
  TPad *pad = new TPad("pad_summary","",0,0,1,1);
  pad->SetLeftMargin(0.15); pad->SetRightMargin(0.15);
  pad->SetBottomMargin(0.14); pad->SetTopMargin(0.13);
  pad->SetLogy(); pad->SetTickx(1);
  pad->Draw(); pad->cd();

  TGraph *gB = new TGraph(N,x,bias2);
  TGraph *gV = new TGraph(N,x,var);
  TGraph *gM = new TGraph(N,x,mse);
  gB->SetLineColor(cbBlue());   gB->SetMarkerColor(cbBlue());
  gV->SetLineColor(cbVermil()); gV->SetMarkerColor(cbVermil());
  gM->SetLineColor(kBlack);     gM->SetMarkerColor(kBlack);
  for(auto *g : {gB,gV,gM}){ g->SetLineWidth(2); g->SetMarkerSize(1.0); }
  // bias^2 and MSE often sit on top of each other; separate them by dash and
  // marker as well as hue
  gB->SetLineStyle(2); gB->SetMarkerStyle(24);   // open circle, dashed
  gV->SetMarkerStyle(21);                        // filled square
  gM->SetLineWidth(3); gM->SetMarkerStyle(20);   // filled circle, thick

  TH1F *frame = pad->DrawFrame(xLow,yLow,xHigh,yHigh);
  frame->GetXaxis()->SetTitle("#it{N}_{iterations}");
  frame->GetYaxis()->SetTitle("relative bias^{2}, variance, MSE");
  frame->GetXaxis()->SetTitleSize(0.048); frame->GetXaxis()->SetTitleOffset(1.20);
  frame->GetYaxis()->SetTitleSize(0.048); frame->GetYaxis()->SetTitleOffset(1.45);
  frame->GetXaxis()->SetLabelSize(0.042); frame->GetYaxis()->SetLabelSize(0.042);
  frame->GetXaxis()->SetNdivisions(N < 12 ? 100 + N : 510);

  // MSE is the thick line, so draw it first and let the dashed bias^2 sit on top
  // of it -- the two coincide whenever bias dominates
  gM->Draw("LP same");
  gV->Draw("LP same");
  gB->Draw("LP same");

  // ---- mark the MSE minimum (N = 0 is the no-unfolding reference, not a
  //      candidate regularisation strength, so start the search at N = 1) ----
  int iBest = 1;
  for(int i = 2; i < N; i++) if(mse[i] < mse[iBest]) iBest = i;

  //  Stop the marker at the MSE point rather than running it to the top of the
  //  frame: it is pointing at one value on one curve, and a full-height line
  //  crosses bias^2, variance and the chi2/ndf overlay on the way, implying it
  //  marks something about all four.
  TLine *li = new TLine();
  li->SetLineStyle(7); li->SetLineColor(kGray+2);
  li->DrawLine(x[iBest],yLow,x[iBest],mse[iBest]);

  // ---- chi2/ndf on a transparent overlay pad with its own right-hand axis.
  //      Goes log if the N = 0 point drags the range over a couple of decades,
  //      which would otherwise flatten every unfolded point onto one line. ----
  double cLow = 1e30, cHigh = -1e30;
  for(int i = 0; i < N; i++){
    if(chi2ndf[i] <= 0.) continue;
    if(chi2ndf[i] < cLow)  cLow  = chi2ndf[i];
    if(chi2ndf[i] > cHigh) cHigh = chi2ndf[i];
  }
  if(cLow > cHigh){ cLow = 0.1; cHigh = 1.; }

  bool logChi2 = (cHigh / cLow) > 20.;
  if(logChi2){
    cLow  *= 0.5;
    cHigh *= 2.;
    cHigh *= std::pow(10., std::log10(cHigh/cLow) * (1./dataFrac - 1.));
  }
  else{
    cLow  -= 0.10 * (cHigh - cLow);
    cHigh  = cLow + (cHigh - cLow) / dataFrac;
  }

  canv->cd();
  TPad *overlay = new TPad("pad_summary_overlay","",0,0,1,1);
  overlay->SetFillStyle(4000); overlay->SetFrameFillStyle(0);
  overlay->SetLeftMargin(0.15); overlay->SetRightMargin(0.15);
  overlay->SetBottomMargin(0.14); overlay->SetTopMargin(0.13);
  overlay->SetLogy(logChi2);
  overlay->Draw(); overlay->cd();

  TH1F *frameOv = overlay->DrawFrame(xLow,cLow,xHigh,cHigh);
  frameOv->GetXaxis()->SetLabelOffset(999); frameOv->GetXaxis()->SetTickLength(0);
  frameOv->GetYaxis()->SetLabelOffset(999); frameOv->GetYaxis()->SetTickLength(0);

  TGraph *gC = new TGraph(N,x,chi2ndf);
  gC->SetLineColor(cbGreen()); gC->SetMarkerColor(cbGreen());
  // diamond, not a triangle: markers here are centre-symmetric by convention
  gC->SetLineWidth(2); gC->SetLineStyle(3); gC->SetMarkerStyle(33);
  gC->SetMarkerSize(1.1);
  gC->Draw("LP same");

  TGaxis *axisC = new TGaxis(xHigh,cLow,xHigh,cHigh,cLow,cHigh,510,
                             logChi2 ? "+LG" : "+L");
  axisC->SetTitle("#it{#chi}^{2} / ndf");
  axisC->SetLineColor(cbGreen()); axisC->SetLabelColor(cbGreen()); axisC->SetTitleColor(cbGreen());
  axisC->SetTitleSize(0.048); axisC->SetTitleOffset(1.35);
  axisC->SetLabelSize(0.042);
  axisC->SetLabelFont(42);  axisC->SetTitleFont(42);
  axisC->Draw();

  TLegend *leg = new TLegend(0.19,0.64,0.86,0.80);
  leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.033);
  leg->SetNColumns(2);
  leg->SetHeader(nOptMeaningful
                 ? Form("min MSE at #it{N} = %.0f",x[iBest])
                 : "no distortion: prior #approx truth, #it{N}_{opt} not meaningful");
  leg->AddEntry(gB,"bias^{2}","lp");
  leg->AddEntry(gV,"variance","lp");
  leg->AddEntry(gM,"MSE = bias^{2}+var","lp");
  leg->AddEntry(gC,"#it{#chi}^{2}/ndf (right)","lp");
  leg->Draw();

  TLatex lat; lat.SetNDC(); lat.SetTextSize(0.042);
  lat.DrawLatex(0.15,0.955,titleStr);
  lat.SetTextSize(0.031); lat.SetTextColor(kGray+2);
  lat.DrawLatex(0.15,0.905,
                (subTitleStr2 && subTitleStr2[0])
                ? Form("%s   |   %s",subTitleStr,subTitleStr2) : subTitleStr);

  canv->SaveAs(outPath);

}


//====================================================================================
//  main
//====================================================================================

void unfoldClosureTest(TString sample = "C1",              // pp, C1, C2, C3, C4
                       int    N_iter_max = 30,             // iterations to scan
                       TString responseMode = "weighted",  // see below
                       TString distortion   = "none",      // see below
                       double ptLow  =  80.,               // metric window, gen pT
                       double ptHigh = 300.,
                       bool   caloJets = false,            // calo response halves
                       int    rebinFactor = 1,             // merge native bins by this
                       TString trainOverride = "",         // explicit response/prior file
                       TString testOverride  = "",         // explicit unfolded/truth file
                       TString outLabel      = "",         // appended to the output tag
                       int    removeXJetsOpt = -1,         // -1 default, 0 keep x, 1 drop x
                       TString responseNameOverride = "",  // e.g. a bJets / muTagged / var-bin response
                       TString distNameOverride = "",      // tiltScan: base name for "<name>_distTilt<expo>",
                                                            // defaults to responseNameOverride/cfg.responseName.
                                                            // Lets the (even) distorted pseudo-data come from a
                                                            // more specific selection (e.g. muTagged) than the
                                                            // (odd) response actually being unfolded with.
                       bool   varRebin = false){           // manually rebin the native (fixed-bin) response(s)
                                                            // to kVarBinEdges here, instead of reading the scan's
                                                            // separately-booked "_var_" histogram

  //  The spectrum that gets unfolded, and the truth it is compared to, are
  //  ALWAYS the pThat-weighted ones.  Only the response changes:
  //
  //   "weighted"         response matrix and prior from the pThat-weighted MC.
  //                      The standard closure test.
  //   "unweighted"       response matrix and prior from the pThat-unweighted MC.
  //                      Unfolds the weighted spectrum with an unweighted
  //                      response, exactly as it comes out of the files.
  //   "unweightedMatrix" migration probabilities P(reco|gen) from the
  //                      pThat-unweighted MC, but each gen row rescaled to the
  //                      weighted truth so the Bayesian prior is the physical
  //                      one.  Separates the effect of the matrix from the
  //                      effect of the prior, which "unweighted" conflates.
  //
  //  pThat-unweighted samples exist for PbPb only.
  //
  //  "distortion" controls what is unfolded:
  //
  //   "none"     the odd MC half is unfolded and compared to its own truth.
  //              A technical closure test.  Note that the prior is then the
  //              truth up to the even/odd split, so there is essentially no
  //              bias for the iterations to remove -- MSE is minimised at
  //              N = 1 by construction and says nothing about the optimal
  //              number of iterations.
  //   "dataMC"   the training truth is reweighted by a power law fitted to the
  //              data/MC ratio, folded through the nominal response to make a
  //              noiseless pseudo-measurement, and unfolded with the NOMINAL
  //              (unreweighted) prior.  The prior is now wrong by a realistic
  //              amount, so bias^2 genuinely falls with iterations and the MSE
  //              minimum is meaningful.
  //   "tilt"     same, but with a hand-set exponent: "tilt" (0.3) or "tilt:0.5".

  responseMode.ToLower();
  const bool useUnwMatrix     = (responseMode == "unweighted" ||
                                 responseMode == "unweightedmatrix");
  const bool reweightToTruthW = (responseMode == "unweightedmatrix");

  if(responseMode != "weighted" && !useUnwMatrix){
    std::cout << "unfoldClosureTest: unknown responseMode \"" << responseMode << "\"\n"
              << "                   choose one of: weighted, unweighted, unweightedMatrix\n";
    return;
  }

  SampleConfig cfg = getSampleConfig(sample);
  if(cfg.ok && caloJets) applyCaloConfig(cfg);
  if(!cfg.ok) return;

  if(!trainOverride.IsNull() || !testOverride.IsNull()){
    if(outLabel.IsNull()){
      std::cout << "unfoldClosureTest: file overrides need an outLabel, so the figures do not\n"
                << "                   overwrite the standard ones.\n";
      return;
    }
    if(!trainOverride.IsNull()) cfg.trainFile = trainOverride;
    if(!testOverride.IsNull())  cfg.testFile  = testOverride;
    std::cout << "file overrides:\n   train: " << cfg.trainFile << "\n   test : " << cfg.testFile << "\n";
  }

  if(!responseNameOverride.IsNull()){
    cfg.responseName = responseNameOverride;
    std::cout << "response name override: " << cfg.responseName << "\n";
  }

  if(useUnwMatrix && cfg.trainFileUnw.IsNull()){
    std::cout << "unfoldClosureTest: no pThat-unweighted response exists for \""
              << cfg.label << "\" -- that production is PbPb only.\n";
    return;
  }

  const TString trainFile = useUnwMatrix ? cfg.trainFileUnw : cfg.trainFile;

  //  Drop unassigned-flavour ("x") jets from the response -- see loadResponse.
  //  Default is off for calo, where refparton_flavor leaves ~17-21% of jets
  //  unassigned.  That number is honest but it is a whole-spectrum average
  //  dominated by the softest bins: measured on the 2026-9-25 C1 calo response
  //  the x fraction is 0.1997 at gen pT 0-20, 0.0337 at 20-40, 0.0032 at 40-60
  //  and 0.0002 everywhere above 80.  So x is not a contaminant of the
  //  measurement region at all -- but those low bins are where the buffer
  //  migration comes from, so it is not obviously irrelevant either, which is
  //  what removeXJetsOpt is for.
  //
  //  (With includeUnmatched off, below, x removal now applies to everything the
  //  test uses -- the matched response is the whole input.)
  const bool removeXJets = (removeXJetsOpt < 0) ? !cfg.isCalo : (removeXJetsOpt == 1);

  //  MATCHED JETS ONLY.  With this off, measured and truth are just the X and Y
  //  projections of the matched 2D, so the response has efficiency 1 and no
  //  fakes: the test asks whether the unfolding recovers a real jet spectrum,
  //  and nothing else.  Fake-jet subtraction is a separate problem and is not
  //  handled here.
  //
  //  It was on until 2026-09-25, and it was the dominant term in the closure
  //  deviation rather than a realistic addition.  The fake histogram escapes the
  //  pThat correlation filter (see PYTHIAHYDJET_scan_response.C), leaving it
  //  built from a handful of events -- even and odd disagree by 31% / 79% /
  //  133% / 293% in C2 / C1 / C3 / C4 over reco pT 80-200, against 0.4-2.0% for
  //  the matched response over the same range.  Turning it off drops
  //  max|unfolded/truth - 1| at iteration 1 from 12.4% to 2.6% in C1, 26.2% to
  //  4.1% in C2, 37.9% to 2.4% in C3 and 15.2% to 1.8% in C4.
  const bool includeUnmatched = false;

  TString respTitle = "response: pThat-weighted";
  TString distTitle = "";
  TString outTag    = cfg.label;
  if(includeUnmatched) outTag += "_unm";
  if(responseMode == "unweighted"){
    respTitle = "response: pThat-unweighted";
    outTag    = Form("%s_unwResp",cfg.label.Data());
  }
  else if(reweightToTruthW){
    respTitle = "response: unwgt. matrix, wgt. prior";
    outTag    = Form("%s_unwMatrix",cfg.label.Data());
  }
  if(cfg.isCalo) outTag += "_caloJets";   // never overwrite the PF outputs
  //  only tag when the x handling differs from this collection's default, so
  //  the existing filenames keep meaning what they meant
  if(removeXJets != !cfg.isCalo) outTag += removeXJets ? "_noX" : "_withX";
  if(!outLabel.IsNull()) outTag += outLabel;
  // with a train override the mode name no longer says what the response is
  if(!trainOverride.IsNull())
    respTitle = trainOverride.Contains("unweighted") ? "response: pThat-unweighted (override)"
                                                     : "response: pThat-weighted (override)";
  if(rebinFactor > 1){
    outTag    += Form("_rebin%d", rebinFactor);
    respTitle += Form(", %d GeV bins", 5*rebinFactor);
  }
  if(varRebin){
    outTag    += "_varRebin";
    respTitle += ", manual variable bins";
  }

  //---- distortion selection ---------------------------------------------------
  distortion.ToLower();
  const bool doDataMC   = (distortion == "datamc");
  const bool doTiltScan = distortion.BeginsWith("tiltscan");
  const bool doTilt     = distortion.BeginsWith("tilt") && !doTiltScan;
  const bool doDistort  = doDataMC || doTilt || doTiltScan;

  if(distortion != "none" && !doDistort){
    std::cout << "unfoldClosureTest: unknown distortion \"" << distortion << "\"\n"
              << "                   choose one of: none, dataMC, tilt, tilt:<exponent>,\n"
              << "                   tiltScan, tiltScan:<exponent>\n";
    return;
  }

  //  "dataMC" works for calo since 2026-09-26; applyCaloConfig points the four
  //  data slots at the manual-JEC calo scans.  Guard against a missing one
  //  rather than letting TFile::Open fail on an empty path deep in the fit.
  if(doDataMC &&
     (cfg.dataMB.IsNull() || cfg.dataJet80.IsNull() || cfg.dataJet100.IsNull())){
    std::cout << "unfoldClosureTest: \"dataMC\" needs the data scans and at least one is missing:\n"
              << "   MinBias: " << (cfg.dataMB    .IsNull() ? TString("MISSING") : cfg.dataMB)     << "\n"
              << "   Jet80  : " << (cfg.dataJet80 .IsNull() ? TString("MISSING") : cfg.dataJet80)  << "\n"
              << "   Jet100 : " << (cfg.dataJet100.IsNull() ? TString("MISSING") : cfg.dataJet100) << "\n";
    return;
  }

  double tiltExpo = 0.3;
  if((doTilt || doTiltScan) && distortion.Contains(":"))
    tiltExpo = TString(distortion(distortion.Index(":")+1,distortion.Length())).Atof();

  const double drawPtLow  =  40.;   // x-range of the closure plot
  const double drawPtHigh = 400.;
  const double ratioLow   = doDistort ? 0.5 : 0.8;  // closure ratio panel range
  const double ratioHigh  = doDistort ? 1.5 : 1.2;
  const double pTref      = 150.;   // w(pTref) = 1
  const TString outDir = "../../figures/unfoldTest/";

  gStyle->SetPalette(kCividis);
  gSystem->mkdir(outDir,kTRUE);

  //---- input ------------------------------------------------------------------
  TFile *file_train = TFile::Open(trainFile);
  TFile *file_test  = TFile::Open(cfg.testFile);
  if(!file_train || file_train->IsZombie() || !file_test || file_test->IsZombie()){
    std::cout << "unfoldClosureTest: could not open\n  " << trainFile
              << "\n  " << cfg.testFile << "\n";
    return;
  }

  TH2D *h_response_train, *h_response_test;
  h_response_train = loadResponse(file_train,cfg.responseName,removeXJets,"h_resp_train");
  h_response_test  = loadResponse(file_test ,cfg.responseName,removeXJets,"h_resp_test");

  TH1D *h_unmReco_train = 0, *h_unmGen_train = 0;
  TH1D *h_unmReco_test  = 0, *h_unmGen_test  = 0;
  if(includeUnmatched){
    file_train->GetObject(cfg.unmRecoName,h_unmReco_train);
    file_train->GetObject(cfg.unmGenName ,h_unmGen_train);
    file_test ->GetObject(cfg.unmRecoName,h_unmReco_test);
    file_test ->GetObject(cfg.unmGenName ,h_unmGen_test);
    if(!h_unmReco_train || !h_unmGen_train || !h_unmReco_test || !h_unmGen_test){
      std::cout << "unfoldClosureTest: could not find the unmatched histograms ("
                << cfg.unmRecoName << ", " << cfg.unmGenName << ")\n";
      return;
    }
  }
  if(!h_response_train || !h_response_test){
    std::cout << "unfoldClosureTest: could not find \"" << cfg.responseName << "\"\n";
    return;
  }

  //  Coarser binning, applied to every input before anything is derived from it.
  auto divides = [&](TH1 *h){ return h && h->GetNbinsX() % rebinFactor == 0; };
  if(rebinFactor < 1 || !divides(h_response_train) ||
     h_response_train->GetNbinsY() % rebinFactor != 0 ||
     (includeUnmatched && (!divides(h_unmReco_train) || !divides(h_unmGen_train)))){
    std::cout << "unfoldClosureTest: rebinFactor " << rebinFactor << " does not divide the native "
              << h_response_train->GetNbinsX() << "-bin axis\n";
    return;
  }
  if(rebinFactor > 1){
    h_response_train->Rebin2D(rebinFactor,rebinFactor);
    h_response_test ->Rebin2D(rebinFactor,rebinFactor);
    if(includeUnmatched){
      h_unmReco_train->Rebin(rebinFactor); h_unmGen_train->Rebin(rebinFactor);
      h_unmReco_test ->Rebin(rebinFactor); h_unmGen_test ->Rebin(rebinFactor);
    }
    std::cout << "rebinned by " << rebinFactor << ": "
              << h_response_train->GetXaxis()->GetBinWidth(1) << " GeV bins\n";
  }

  //  Manual variable-bin rebin, straight from the native response -- see
  //  rebinToVarBins above for why this exists alongside the scan's own "_var_"
  //  booking.
  if(varRebin){
    h_response_train = rebinToVarBins(h_response_train,"h_resp_train_varRebin");
    h_response_test  = rebinToVarBins(h_response_test ,"h_resp_test_varRebin");
    if(includeUnmatched){
      std::cout << "unfoldClosureTest: varRebin does not support includeUnmatched\n";
      return;
    }
    std::cout << "manually rebinned to " << kVarBinN-1 << " variable-width bins ("
              << kVarBinEdges[0] << "-" << kVarBinEdges[kVarBinN-1] << " GeV)\n";
  }

  //  Rescale every gen row of the unweighted matrix by (weighted truth) /
  //  (unweighted truth).  P(reco|gen) is untouched -- only the gen marginal,
  //  i.e. the Bayesian prior, is swapped for the physical one.
  //  (RooUnfoldBayes::SetPriors would be the obvious route, but in this build
  //  it leaves _N0C at zero, which blows up the error propagation.)
  if(reweightToTruthW){

    TFile *file_trainW = TFile::Open(cfg.trainFile);
    TH2D  *h_response_trainW = 0;
    if(file_trainW && !file_trainW->IsZombie())
      h_response_trainW = loadResponse(file_trainW,cfg.responseName,removeXJets,"h_resp_trainW");
    if(!h_response_trainW){
      std::cout << "unfoldClosureTest: could not read the pThat-weighted response from\n  "
                << cfg.trainFile << "\n";
      return;
    }
    if(rebinFactor > 1) h_response_trainW->Rebin2D(rebinFactor,rebinFactor);

    TH1D *h_truth_unw = projY(h_response_train ,"h_truth_unw");
    TH1D *h_truth_w   = projY(h_response_trainW,"h_truth_w");

    TH2D *h_hybrid = (TH2D*) h_response_train->Clone("h_response_hybrid");
    h_hybrid->SetDirectory(0);

    for(int iy = 1; iy <= h_hybrid->GetNbinsY(); iy++){
      double u = h_truth_unw->GetBinContent(iy);
      double s = (u > 0.) ? h_truth_w->GetBinContent(iy) / u : 0.;
      for(int ix = 1; ix <= h_hybrid->GetNbinsX(); ix++){
        h_hybrid->SetBinContent(ix,iy, s * h_hybrid->GetBinContent(ix,iy));
        h_hybrid->SetBinError  (ix,iy, s * h_hybrid->GetBinError  (ix,iy));
      }
    }

    h_response_train = h_hybrid;

  }

  // training half -> response matrix and prior
  TH1D *h_meas_train  = projX(h_response_train,"h_meas_train");
  TH1D *h_truth_train = projY(h_response_train,"h_truth_train");

  // testing half -> what gets unfolded, and what we compare against
  TH1D *h_meas_test  = projX(h_response_test,"h_meas_test");
  TH1D *h_truth_test = projY(h_response_test,"h_truth_test");

  //  unmatched reco -> fakes, unmatched gen -> inefficiency.  Added to the
  //  measured and truth spectra; RooUnfold derives both from the difference
  //  against the projections of the matched 2D.
  if(includeUnmatched){
    h_meas_train ->Add(h_unmReco_train);
    h_truth_train->Add(h_unmGen_train);
    h_meas_test  ->Add(h_unmReco_test);
    h_truth_test ->Add(h_unmGen_test);
  }

  //  WHICH HALF UNFOLDS.  Both modes keep the pseudo-data and the response in
  //  different halves, so the two are statistically independent:
  //
  //    split-sample   measured = odd, truth = odd, response = EVEN
  //    distorted      even gen reweighted, folded through the EVEN response to
  //                   make the pseudo-measurement; unfolded with the ODD
  //                   response, whose prior is the undistorted odd truth
  //
  //  The distorted mode used the even response on both sides until 2026-09-25.
  //  That folded the pseudo-data through the very matrix that then unfolded it,
  //  so the fluctuations were common to both sides and cancelled -- hence its
  //  chi2/ndf of 0.04, which measured correlation rather than closure.  Taking
  //  the response from the odd half removes that, at the cost of the prior
  //  being wrong by the even/odd difference as well as by the tilt, which is
  //  what it is supposed to be wrong about anyway.
  const bool unfoldWithOddHalf = doDistort;
  RooUnfoldResponse response(unfoldWithOddHalf ? h_meas_test     : h_meas_train,
                             unfoldWithOddHalf ? h_truth_test    : h_truth_train,
                             unfoldWithOddHalf ? h_response_test : h_response_train,
                             "response",cfg.title);

  //---- distorted truth --------------------------------------------------------
  //  Scale gen row t of the training response by w(t).  That leaves P(reco|gen)
  //  untouched -- w cancels in the row normalisation -- so the only thing the
  //  unfolding has to recover from is the wrong prior.  Folding through the same
  //  matrix that will be used to unfold makes the pseudo-measurement noiseless,
  //  which keeps bias^2 free of test-sample statistical scatter.
  double usedExpo = 0.;   // kept for the figure label
  if(doDistort){

    double expo = tiltExpo;

    //  relative statistical precision to hand the pseudo-measurement.  Empty
    //  means "keep the MC errors"; for dataMC it is filled from the data, since
    //  otherwise the variance term reflects MC statistics and the bias/variance
    //  trade-off -- and therefore the MSE minimum -- is meaningless.
    std::vector<double> relErr;

    if(doDataMC){
      TH1D *h_data = loadDataSpectrum(cfg);
      if(!h_data) return;
      double fitLow = std::max(ptLow, cfg.dataMinPt);
      if(fitLow >= ptHigh){
        std::cout << "unfoldClosureTest: no usable data range for the fit\n";
        return;
      }
      if(!buildDataMCDistortion(h_data,h_meas_train,fitLow,ptHigh,pTref,
                                cfg.dataMinPt,expo,relErr,cfg.title,
                                //  outTag, not cfg.label: the label is just "C1",
                                //  so the calo run overwrote the PF derivation
                                //  figure for the same class.  Harmless while
                                //  dataMC was PF-only; not once calo can run it.
                                Form("%sdistortion_%s.pdf",outDir.Data(),outTag.Data())))
        return;
    }

    //  "tilt"/"dataMC" scale each gen row of the already-filled training
    //  response by w evaluated at the row's bin CENTER -- a discretised
    //  approximation of the per-event weight.  "tiltScan" instead reads the
    //  response PYTHIA_scan_response.C / PYTHIAHYDJET_scan_response.C already
    //  built with w(genPt) applied to each event's own continuous gen pT at
    //  fill time, so there is no bin-center approximation. Either way
    //  P(reco|gen) is untouched, so the only thing the unfolding has to
    //  recover from is the wrong prior.
    TH2D *h_dist = 0;

    if(doTiltScan){
      TString distBase = distNameOverride.IsNull() ? cfg.responseName : distNameOverride;
      TString distName = distBase + Form("_distTilt%.1f",expo);
      h_dist = loadResponse(file_train,distName,false,"h_response_distorted");
      if(!h_dist){
        std::cout << "unfoldClosureTest: could not find \"" << distName << "\" in\n  "
                  << trainFile << "\n";
        return;
      }
      if(rebinFactor > 1) h_dist->Rebin2D(rebinFactor,rebinFactor);
      if(varRebin) h_dist = rebinToVarBins(h_dist,"h_response_distorted_varRebin");
    }
    else{
      h_dist = (TH2D*) h_response_train->Clone("h_response_distorted");
      h_dist->SetDirectory(0);
      for(int iy = 1; iy <= h_dist->GetNbinsY(); iy++){
        double w = distortionWeight(h_dist->GetYaxis()->GetBinCenter(iy),expo,pTref);
        for(int ix = 1; ix <= h_dist->GetNbinsX(); ix++){
          h_dist->SetBinContent(ix,iy, w * h_dist->GetBinContent(ix,iy));
          h_dist->SetBinError  (ix,iy, w * h_dist->GetBinError  (ix,iy));
        }
      }
    }

    h_meas_test  = projX(h_dist,"h_meas_distorted");
    h_truth_test = projY(h_dist,"h_truth_distorted");

    //  The pseudo-measurement is folded from the training half, so its
    //  unmatched jets come from there too.  Misses are a gen-level quantity and
    //  follow the reweighting; fakes are underlying-event background and do not
    //  depend on the signal truth shape, so they are left unweighted.
    //  (tiltScan has no per-flavor/muTag miss histogram to read -- dead code
    //  either way while includeUnmatched is hardcoded off, above.)
    if(includeUnmatched && !doTiltScan){
      TH1D *h_unmGen_w = (TH1D*) h_unmGen_train->Clone("h_unmGen_distorted");
      h_unmGen_w->SetDirectory(0);
      for(int b = 1; b <= h_unmGen_w->GetNbinsX(); b++){
        double w = distortionWeight(h_unmGen_w->GetBinCenter(b),expo,pTref);
        h_unmGen_w->SetBinContent(b, w * h_unmGen_w->GetBinContent(b));
        h_unmGen_w->SetBinError  (b, w * h_unmGen_w->GetBinError(b));
      }
      h_meas_test ->Add(h_unmReco_train);
      h_truth_test->Add(h_unmGen_w);
    }

    if(!relErr.empty())
      for(int b = 1; b <= h_meas_test->GetNbinsX() && b < (int)relErr.size(); b++)
        h_meas_test->SetBinError(b, relErr[b] * h_meas_test->GetBinContent(b));

    distTitle = doTiltScan
      ? Form("distorted truth (scan-level, per-event #it{w}): #it{w} #propto #it{p}_{T,gen}^{%.2f}",expo)
      : Form("distorted truth: #it{w} #propto #it{p}_{T}^{%.2f}",expo);
    usedExpo  = expo;
    outTag   += doDataMC ? "_distDataMC" : Form("_distTilt%s%.2f", doTiltScan ? "Scan" : "", expo);

    //  the tilt's own derivation figure (dataMC draws its own, above)
    if(doTilt || doTiltScan)
      drawTiltDerivation(h_truth_train, h_truth_test, h_meas_train, h_meas_test,
                         expo, pTref, ptLow, ptHigh, cfg.title,
                         Form("%sdistortion_%s.pdf",outDir.Data(),outTag.Data()));

  }

  //---- scan N = 0 ... N_iter_max ----------------------------------------------
  //  N = 0 is the measured spectrum itself: the Bayesian iteration has not been
  //  applied at all, so it is the reference the unfolding has to improve on.
  const int N_points = N_iter_max + 1;

  std::vector<TH1D*>  h_unfold(N_points);
  std::vector<double> x(N_points), bias2(N_points), var(N_points),
                      mse(N_points), chi2(N_points), chi2ndf(N_points),
                      maxdev(N_points), mseAbs(N_points);
  std::vector<int>    ndf(N_points);

  RooUnfoldBayes unfold(&response,h_meas_test,1);
  unfold.SetVerbose(0);
  if(includeUnmatched) unfold.HandleFakes(true);   // off by default in RooUnfold

  for(int n = 0; n <= N_iter_max; n++){

    TH1D *h = 0;

    if(n == 0){
      h = (TH1D*) h_meas_test->Clone(Form("h_unfold_%s_iter0",cfg.label.Data()));
      h->SetLineStyle(2);
    }
    else{
      unfold.SetIterations(n);
      h = (TH1D*) unfold.Hunfold();
      h->SetName(Form("h_unfold_%s_iter%i",cfg.label.Data(),n));

      // errors from the diagonal of the full unfolding covariance matrix
      TMatrixD cov = unfold.Eunfold();
      for(int b = 1; b <= h->GetNbinsX(); b++){
        double v = cov(b-1,b-1);
        h->SetBinError(b, v > 0. ? std::sqrt(v) : 0.);
      }
    }

    h->SetDirectory(0);
    h->SetLineColor(iterColor(n,N_iter_max));
    h->SetLineWidth(2);
    h->SetStats(0);
    h_unfold[n] = h;

    ClosureMetrics m = computeMetrics(h,h_truth_test,ptLow,ptHigh);
    x[n]       = n;
    bias2[n]   = m.bias2;
    var[n]     = m.variance;
    mse[n]     = m.mse;
    mseAbs[n]  = m.mseAbs;
    chi2[n]    = m.chi2;
    ndf[n]     = m.ndf;
    chi2ndf[n] = (m.ndf > 0) ? m.chi2 / m.ndf : 0.;
    maxdev[n]  = m.maxDev;

  }

  //---- plots ------------------------------------------------------------------
  //  the iteration the MSE scan picks, needed by the reduced figure below as
  //  well as by the table
  int iBest = 1;
  for(int n = 2; n < N_points; n++) if(mse[n] < mse[iBest]) iBest = n;

  //  A minimum sitting on the last point is not a minimum -- the scan stopped
  //  before MSE turned over, and the number quoted is wherever it was cut off.
  //  This is silent otherwise, and it bit a 6-iteration scan where C3 and C4
  //  both "optimised" at 6 while still falling (they turn over at 9).
  if(doDistort && iBest == N_iter_max)
    printf("\n  WARNING: MSE is still falling at N = %d, the last iteration scanned.\n"
           "           The minimum is beyond the scan -- rerun with a larger N_iter_max\n"
           "           before quoting an optimal iteration count.\n", N_iter_max);

  drawClosure(N_iter_max,h_unfold.data(),h_truth_test,
              drawPtLow,drawPtHigh,ptLow,ptHigh,ratioLow,ratioHigh,
              cfg.title,respTitle,distTitle,
              Form("%sclosure_%s.pdf",outDir.Data(),outTag.Data()));

  //  What the reduced figure says about itself, bottom left of the spectra pad.
  //  The jet collection moves out of the title and in here; which half is which
  //  is not cosmetic, because the two modes do not use the split the same way
  //  and the figures are otherwise indistinguishable.
  //
  //  Both axes are |eta| < 1.6: the scans cut reco jets at 1.6 and gen jets at
  //  etaMax, which common.h sets to the same value.
  std::vector<TString> infoLines;
  infoLines.push_back(Form("#it{R} = 0.4 %s, |#eta| < 1.6, matched jets only",
                           cfg.isCalo ? "calo jets" : "PF jets"));
  if(doDistort){
    //  The weight is a function of GEN pT alone -- distortionWeight is evaluated
    //  on the Y (gen) axis and applied to a whole gen row at once -- so the reco
    //  side moves only by being folded through an unchanged P(reco|gen).  The
    //  fakes are added unweighted, being background that does not track the
    //  signal truth shape.  And h_dist clones the TRAINING response, so the
    //  pseudo-data, its truth and the response all come from the even half:
    //  the odd file is read but nothing from it enters this mode.
    infoLines.push_back(Form("even gen reweighted #it{w} #propto #it{p}_{T,gen}^{%.2f}, folded to reco", usedExpo));
    infoLines.push_back("folded with even response, unfolded with odd");
  }
  else{
    infoLines.push_back("response: even half, unfolded: odd half");
  }

  //  the collection now lives in the info block, so strip it from the title
  TString redTitle = cfg.title;
  redTitle.ReplaceAll(", calo jets","");

  drawClosureReduced(h_unfold[0],h_unfold[iBest],iBest,doDistort,h_truth_test,
                     drawPtLow,drawPtHigh,ptLow,ptHigh,ratioLow,ratioHigh,
                     redTitle,infoLines,
                     Form("%sclosure_reduced_%s.pdf",outDir.Data(),outTag.Data()));

  drawSummary(N_iter_max,x.data(),bias2.data(),var.data(),mse.data(),chi2ndf.data(),
              cfg.title,respTitle,distTitle,doDistort,
              Form("%sclosure_summary_%s.pdf",outDir.Data(),outTag.Data()));

  //---- table ------------------------------------------------------------------
  printf("\n=== %s [%s] : %s closure, %.0f < gen pT < %.0f GeV ===\n",
         cfg.title.Data(), outTag.Data(),
         doDistort ? "distorted-truth" : "split-sample", ptLow, ptHigh);
  printf("  train: %s\n  test:  %s\n\n",trainFile.Data(),cfg.testFile.Data());
  printf("  bias^2, variance and MSE are RELATIVE: each bin divided by its own\n"
         "  truth before squaring, so no single high-yield bin dominates the sum.\n"
         "  MSE_abs is the old counts^2 version, kept for comparison only.\n\n");
  printf("%-6s  %12s  %12s  %12s  %12s  %12s  %5s  %10s  %11s\n",
         "iter","bias^2","variance","MSE","MSE_abs","chi2","ndf","chi2/ndf","max|r-1|");
  for(int n = 0; n < N_points; n++){
    printf("%-6.0f  %12.4g  %12.4g  %12.4g  %12.4g  %12.4g  %5d  %10.3f  %10.2f%%%s\n",
           x[n],bias2[n],var[n],mse[n],mseAbs[n],chi2[n],ndf[n],chi2ndf[n],100.*maxdev[n],
           n == 0     ? "   (measured, no unfolding)" :
           (n == iBest && doDistort) ? "   <- min MSE" : "");
  }
  printf("\n");

}
