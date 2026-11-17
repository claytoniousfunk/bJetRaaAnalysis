// Muon-tag frequency correction: bring the pp CALO-jet b-jet spectrum from
// muon-tagged b jets back to all (inclusive) b jets.
//
//   inclusive b spectrum = muon-tagged b spectrum / f(gen jet pT)
//
// f is the fraction of gen b jets (|flavor| 5 plus gluon splitting 17) carrying
// a muon tag, on the GEN jet pT axis (the spectrum is unfolded to gen level).
// Sources, chosen by fSource:
//
//   "scanRecoTag"  (current default) gen b jets tagged with the ANALYSIS reco
//                  muon (tight ID, pT > 15 GeV, |eta| < 2, W-decay veto,
//                  deltaR < 0.4) / all gen b jets, from the calo PYTHIA scan:
//                  h_inclGenJetPt_inclRecoMuonTag_flavor / h_inclGenJetPt_flavor.
//                  One factor for branching fraction x acceptance x MC muon
//                  efficiency, including what a flat efficiency misses: muons
//                  migrating across the 15 GeV threshold and the tight ID inside
//                  jets. Relative to the gen tag it is 0.97 at 80-100 GeV and
//                  0.89 at 300-500 GeV, against a flat 0.97 (MC) / 0.966 (data
//                  tag-and-probe). Pairs with the data/MC muon scale factor of
//                  correctMuonRecoEff_caloJets_pp.C, not the absolute efficiency.
//                  The deltaR is to the gen-jet axis (data: reco-jet axis).
//
//   "scanGenTag"   (comparison; stored as h_muTagFrequency_genTag) gen-muon-tagged / all gen b jets from the calo
//                  PYTHIA scan: h_inclGenJetPt_inclGenMuonTag_flavor /
//                  h_inclGenJetPt_flavor. Gen-muon tag (PYTHIA_scan.C): gen muon
//                  pT > 15 GeV in the tracker eta acceptance, deltaR < 0.4, not a
//                  W-decay muon, each muon used once.
//                  Needs a scan tagged _genMuTagFix: earlier scans reused the
//                  reco-jet loop's per-event matchFlag / matchFlagR in the gen-jet
//                  loop without resetting them, so a muon that had tagged a reco
//                  jet could not tag its gen jet (b-jet f 0.05-1% instead of
//                  ~5-15%). Fixed 2026-10-04; mcPath is the fixed calo scan. It
//                  agrees with responseRecoTag to 1-4%. Gen-jet axis, gen pT, to
//                  match the unfolded spectrum: the reco-axis fraction
//                  (calculateBJetsPerZ.cc corrFactor_1) is ~3x smaller for calo
//                  jets because muon-tagged calo b jets reconstruct at ~0.7 of gen
//                  pT, a shift the unfolding already corrects.
//
//   "responseRecoTag"  (cross-check; was the stopgap before the rescan) from the calo response
//                  scans (PF flavor, even + odd, 2026-10-01): the gen projection
//                  of h_matchedRecoJetPt_genJetPt_bJets_muTagged over that of
//                  h_matchedRecoJetPt_genJetPt_bJets, i.e. the fraction of
//                  matched gen b jets with a RECO muon tag (tight, pT > 15 GeV,
//                  |eta| < 2, no trigger), divided by the MC muon efficiency
//                  effMuMC = 0.9708 (tight | gen) to approximate the gen-muon
//                  fraction -- the spectrum is already corrected for the muon
//                  reconstruction efficiency (tag-and-probe). Matched gen jets
//                  only, as in the unfolding. The 0.9708 is a single number and
//                  not the same quantity as the tag-and-probe efficiency; replace
//                  this source with scanGenTag once the rescan is in.
//
// Uncertainty: the MC statistical error on f (binomial, TH1::Divide "B"),
// added as a new systematic in quadrature with the existing ones.
//
// GLUON SPLITTING. f uses PYTHIA's natural g->bbbar share. g->bbbar jets hold
// two b hadrons and are tagged more often. The b-purity fit assumed the tagged
// bGS share raised by +0.175; the f that implies is stored for comparison
// (h_muTagFrequency_GSshift) but NOT used.
//
// Input   rootFiles/CorrectedBJetSpectra/Data/bJetSpectrum_muRecoCorr_caloJets_pp.root
// Output  rootFiles/CorrectedBJetSpectra/Data/bJetSpectrum_inclusive_caloJets_pp.root
//           h_muTagFrequency           f vs gen jet pT (MC stat error)
//           h_muTagFrequency_FC/_GS    f for flavor-creation / gluon-splitting b jets
//           h_muTagFrequency_GSshift   f with the tagged bGS share + 0.175 (not used)
//           h_bJetPt_inclusive         corrected / f, stat error
//           h_bJetPt_inclusive_sysAbs  same, error = total systematic
//           h_incl_statRel, h_incl_sysRel, h_incl_sysRel_<source>
//         figures/bJetSpectra/muonTagFrequency_caloJets_pp.pdf
//         figures/bJetSpectra/bJetSpectrum_inclusive_caloJets_pp.pdf
//         ... and *_simple.pdf versions of both figures for talks
//         figures/systematics/bJetSpectrum_inclusive_sysBreakdown_caloJets_pp.pdf
//
// Usage: root -l -b -q correctMuonTagFrequency_caloJets_pp.C
// Run from: src/calculateBJetsPerZ/

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TBox.h"
#include "TLatex.h"
#include "TSystem.h"
#include "TColor.h"
#include "TStyle.h"
#include "../../headers/plotting/plotStyle.h"
#include "../../headers/plotting/ratioPanel.h"
#include "../../headers/functions/caloJetBPurityFit_pp.h"   // ptEdges, NPt, ptReportMin, mcPath

// paths follow the jet collection (collTag() in the names, set in the macro)
const char *inPath = nullptr, *outRoot = nullptr, *outFigF = nullptr, *outFig = nullptr, *outFigS = nullptr;

const double gsShift = 0.175;   // the purity fit's tagged-bGS shift, for the comparison only

// see the header: "responseRecoTag" (provisional) until the scans are rerun.
// PFJets: "responseTagOverGen" (set in the macro), see RESPONSE-CONSISTENT f below.
TString fSource = "scanRecoTag";
// RESPONSE-CONSISTENT f (PFJets, 2026-10-07). f = [gen projection of the muon-
// tagged b response, h_matchedRecoJetPt_genJetPt_bJets of the plain-PF
// muTaggedJetsNoTrigger scan] / [all gen b jets, h_inclGenJetPt_flavor |flavor|
// 5,17, SAME file]: the unfolding returns exactly the matched muon-tagged gen b
// jets, so dividing by this f closes by construction (same events, weights and
// tag definition). The scanRecoTag f of the calo-matched PYTHIA scan does NOT:
// its gen-jet tag histograms count only gen jets whose reco jet passes the calo
// match, which loses most low-pT muon-tagged b jets (measureCaloMatchEff.C), so
// that f was 0.009 / 0.041 / 0.064 at 60-70 / 80-100 / 100-120 GeV against
// 0.048 / 0.067 / 0.081 here, and the pp b fraction came out 6.6% at 80-100 GeV
// against PYTHIA's 3.5% (MC closure of unfold / f: 1.70 at 80-100). Reco-muon
// tag (MC efficiency inside), as scanRecoTag: pairs with the data/MC scale
// factor of correctMuonRecoEff_caloJets_pp.C. FC and GS are not separable in
// the response (b merged), so f_FC / f_GS stay empty.
const double  effMuMC = 0.9708;   // MC muon efficiency, tight | gen (responseRecoTag only)
// response scans: ppRespDir / ppRespFiles (jetCollectionInputs.h), summed

const int NIn = 6;
const char *inKey[NIn]   = {"bPurity", "spectrumJEU", "iterations", "responseStat", "turnOn", "muonRecoEff"};
const char *inLabel[NIn] = {"b purity", "jet energy scale", "iterations",
                            "response-matrix MC stat.", "70-80 GeV turn-on correction", "muon efficiency"};
const char *inHex[NIn]   = {"#009E73", "#0072B2", "#E69F00", "#999999", "#56B4E9", "#CC79A7"};
const int   inStyle[NIn] = {2, 3, 7, 5, 6, 9};

bool reported(int i){ return ptEdges[i-1] >= ptReportMin - 1e-6; }

// gen jet pT projection of a (pT, flavor) map for the listed flavors, in ptEdges
TH1D* flavorPt(TH2D *H, std::vector<int> flavors, const char *name)
{
  TH1D *o = nullptr;
  for(int fl : flavors){
    int b = H->GetYaxis()->FindBin(fl + 0.1);   // flavor bins are [v, v+1)
    TH1D *p = H->ProjectionX(Form("%s_f%d", name, fl), b, b); p->SetDirectory(nullptr);
    if(!o){ o = p; o->SetName(Form("%s_raw", name)); } else { o->Add(p); delete p; }
  }
  TH1D *r = rebinTo(o, NPt, ptEdges, name);
  delete o;
  return r;
}

// gen projection of a (reco, gen) response, summed over ppRespFiles, in ptEdges
TH1D* respGen(const char *hname, const char *name)
{
  TH2D *h = nullptr;
  for(const std::string &fn : ppRespFiles){
    TFile *f = TFile::Open(Form("%s%s", ppRespDir, fn.c_str())); TH2D *r = nullptr;
    if(f && !f->IsZombie()) f->GetObject(hname, r);
    if(!r){ printf("ERROR: %s missing in %s\n", hname, fn.c_str()); return nullptr; }
    if(!h){ h = (TH2D*) r->Clone(Form("%s_2d", name)); h->SetDirectory(nullptr); } else h->Add(r);
  }
  TH1D *p = h->ProjectionY(Form("%s_raw", name), 0, h->GetNbinsX()+1); p->SetDirectory(nullptr);
  TH1D *r = rebinTo(p, NPt, ptEdges, name);
  delete p; delete h;
  return r;
}

TH1D* ratioB(TH1D *num, TH1D *den, const char *name)
{
  TH1D *r = (TH1D*) num->Clone(name); r->SetDirectory(nullptr);
  r->Divide(num, den, 1., 1., "B");
  return r;
}

// coll: jet collection (jetCollectionInputs.h; "" = BJET_JETCOLL env or caloJets)
void correctMuonTagFrequency_caloJets_pp(const char *coll = "")
{
  initPlotStyle();
  setJetCollection(coll);
  if(!ptBinGuard(purVar())) return;   // 7-bin jet-pT axis variant (caloJetBPurityFitCore.h)
  const char *base = "/home/clayton/Analysis/code/bJetRaaAnalysis";
  inPath  = strdup(Form("%s/rootFiles/CorrectedBJetSpectra/Data/bJetSpectrum_muRecoCorr_%s_pp%s.root", base, collTag(), purVar()));
  outRoot = strdup(Form("%s/rootFiles/CorrectedBJetSpectra/Data/bJetSpectrum_inclusive_%s_pp%s.root", base, collTag(), purVar()));
  outFigF = strdup(Form("%s/figures/bJetSpectra/muonTagFrequency_%s_pp%s.pdf", base, collTag(), purVar()));
  outFig  = strdup(Form("%s/figures/bJetSpectra/bJetSpectrum_inclusive_%s_pp%s.pdf", base, collTag(), purVar()));
  outFigS = strdup(Form("%s/figures/systematics/bJetSpectrum_inclusive_sysBreakdown_%s_pp%s.pdf", base, collTag(), purVar()));
  for(const char *p : {outRoot, outFigF, outFigS}) gSystem->mkdir(gSystem->DirName(p), kTRUE);

  if(jetColl == kPFJets) fSource = "responseTagOverGen";
  // ---- f ---------------------------------------------------------------------
  TH1D *allFC = nullptr, *allGS = nullptr, *tagFC = nullptr, *tagGS = nullptr;
  TH1D *fGenTag = nullptr;   // gen-muon-tag fraction, for comparison
  if(fSource == "scanGenTag" || fSource == "scanRecoTag"){
    TFile *fM = TFile::Open(ppTagFreqPath);   // the template scan, except PFJets (jetCollectionInputs.h)
    if(!fM || fM->IsZombie()){ printf("ERROR: cannot open %s\n", ppTagFreqPath); return; }
    printf("tag-frequency scan: %s\n", gSystem->BaseName(ppTagFreqPath));
    TH2D *Hall = nullptr, *Htag = nullptr, *Hgen = nullptr;
    fM->GetObject("h_inclGenJetPt_flavor", Hall);
    fM->GetObject(fSource == "scanRecoTag" ? "h_inclGenJetPt_inclRecoMuonTag_flavor" : "h_inclGenJetPt_inclGenMuonTag_flavor", Htag);
    fM->GetObject("h_inclGenJetPt_inclGenMuonTag_flavor", Hgen);
    if(!Hall || !Htag || !Hgen){ printf("ERROR: gen-jet flavor histograms missing in %s\n", ppTagFreqPath); return; }
    allFC = flavorPt(Hall, {-5, 5}, "allFC"); allGS = flavorPt(Hall, {17}, "allGS");
    tagFC = flavorPt(Htag, {-5, 5}, "tagFC"); tagGS = flavorPt(Htag, {17}, "tagGS");
    TH1D *gB = flavorPt(Hgen, {-5, 5, 17}, "genTagB"), *aB = flavorPt(Hall, {-5, 5, 17}, "allBforGen");
    fGenTag = ratioB(gB, aB, "h_muTagFrequency_genTag");
    fGenTag->SetTitle("fraction of gen b jets with a GEN muon tag (comparison);gen #it{p}_{T}^{jet} [GeV];#it{f}_{#mu-tag}");
  }
  else if(fSource == "responseRecoTag"){
    // the response scans book b (flavor creation + splitting together) only;
    // FC and GS are not separable there, so f_FC / f_GS / f_shift are left empty
    TH1D *tagB0 = respGen("h_matchedRecoJetPt_genJetPt_bJets_muTagged", "tagB0");
    TH1D *allB0 = respGen("h_matchedRecoJetPt_genJetPt_bJets", "allB0");
    if(!tagB0 || !allB0) return;
    tagB0->Scale(1./effMuMC);   // reco-muon tag -> approximate gen-muon tag
    allFC = allB0; tagFC = tagB0;
    allGS = (TH1D*) allB0->Clone("allGS"); allGS->Reset();
    tagGS = (TH1D*) tagB0->Clone("tagGS"); tagGS->Reset();
  }
  else if(fSource == "responseTagOverGen"){
    TH1D *tagB0 = respGen(ppRespName, "tagB0");
    // DENOMINATOR. February muTaggedJetsNoTrigger scan (one file, muon-tagged jets only): all gen b
    // jets, h_inclGenJetPt_flavor |flavor| 5, 17. 2026-11-12 anaMuTag rescan (ppRespName ends in
    // _muTagged; even + odd summed): ALL MATCHED gen b jets, the gen projection of the untagged
    // h_matchedRecoJetPt_genJetPt_bJets (2026-11-11). Its h_inclGenJetPt_flavor is filled for every
    // gen jet and the pThat-correlation filter acts on the leading RECO jet only, so low-pThat events
    // with a mis-reconstructed high-pT gen jet (no reco partner) leave single entries carrying 6-32%
    // of a 20 GeV bin (effective entries 6-138 at 140-250 GeV, against 20-40k elsewhere; the February
    // scan's removeHYDJETjet0p35CutOnGen removed them): f fell 23% at 150-200 GeV. In the clean bins
    // matched b / all gen b = 0.992-1.000 (1.000 above 120 GeV), so the matched denominator differs
    // from all gen b by < 1% and carries none of the outliers; the unfolding returns matched gen b.
    const TString nm = ppRespName;
    const bool matchedDen = nm.EndsWith("_muTagged");
    TH1D *allB0 = nullptr;
    if(matchedDen){
      const TString denName = TString(nm(0, nm.Length() - 9));   // strip "_muTagged" (copy: TSubString::Data() is not terminated)
      allB0 = respGen(denName.Data(), "allB0");
      printf("muon-tag frequency denominator: all matched gen b jets (%s)\n", denName.Data());
    } else {
      TH2D *G = nullptr;
      for(const std::string &fn : ppRespFiles){
        TFile *fR = TFile::Open(Form("%s%s", ppRespDir, fn.c_str())); TH2D *g = nullptr;
        if(fR && !fR->IsZombie()) fR->GetObject("h_inclGenJetPt_flavor", g);
        if(!g){ printf("ERROR: h_inclGenJetPt_flavor missing in %s\n", fn.c_str()); return; }
        if(!G){ G = (TH2D*) g->Clone("G_inclGenFlavor"); G->SetDirectory(nullptr); } else G->Add(g);
      }
      if(G) allB0 = flavorPt(G, {-5, 5, 17}, "allB0");
      printf("muon-tag frequency denominator: all gen b jets (h_inclGenJetPt_flavor)\n");
    }
    if(!tagB0 || !allB0){ printf("ERROR: responseTagOverGen inputs missing (%s and its denominator)\n", ppRespName); return; }
    allFC = allB0; tagFC = tagB0;
    allGS = (TH1D*) allB0->Clone("allGS"); allGS->Reset();
    tagGS = (TH1D*) tagB0->Clone("tagGS"); tagGS->Reset();
    // the old scan f, for comparison only (h_muTagFrequency_genTag slot)
    TFile *fM = TFile::Open(ppTagFreqPath); TH2D *Hall = nullptr, *Htag = nullptr;
    if(fM && !fM->IsZombie()){ fM->GetObject("h_inclGenJetPt_flavor", Hall); fM->GetObject("h_inclGenJetPt_inclRecoMuonTag_flavor", Htag); }
    if(Hall && Htag){ fGenTag = ratioB(flavorPt(Htag, {-5, 5, 17}, "oTag"), flavorPt(Hall, {-5, 5, 17}, "oAll"), "h_muTagFrequency_genTag");
      fGenTag->SetTitle("scanRecoTag f of the calo-matched scan (comparison; biased by the calo match);gen #it{p}_{T}^{jet} [GeV];#it{f}_{#mu-tag}"); }
  }
  else { printf("ERROR: unknown fSource %s\n", fSource.Data()); return; }
  printf("\nmuon-tag frequency source: %s%s\n", fSource.Data(),
         fSource == "responseRecoTag" ? Form(" (PROVISIONAL: reco-muon tag / %.4f; gen-tag scans need a rescan)", effMuMC) : "");

  TH1D *allB = (TH1D*) allFC->Clone("allB"); allB->Add(allGS);
  TH1D *tagB = (TH1D*) tagFC->Clone("tagB"); tagB->Add(tagGS);

  TH1D *fB  = ratioB(tagB,  allB,  "h_muTagFrequency");
  TH1D *fFC = ratioB(tagFC, allFC, "h_muTagFrequency_FC");
  TH1D *fGS = ratioB(tagGS, allGS, "h_muTagFrequency_GS");
  fB->SetTitle("fraction of gen b jets with a muon tag;gen #it{p}_{T}^{jet} [GeV];#it{f}_{#mu-tag}");
  fFC->SetTitle("flavor-creation b jets;gen #it{p}_{T}^{jet} [GeV];#it{f}_{#mu-tag}");
  fGS->SetTitle("gluon-splitting b jets;gen #it{p}_{T}^{jet} [GeV];#it{f}_{#mu-tag}");
  const bool haveFCGS = fSource.BeginsWith("scan");

  // f implied by the purity fit's tagged-bGS shift (comparison only):
  //   1/f = (1 - s)/f_FC + s/f_GS,  s = tagged GS share (natural + shift)
  TH1D *fShift = (TH1D*) fB->Clone("h_muTagFrequency_GSshift"); fShift->SetDirectory(nullptr);
  fShift->SetTitle(Form("f with the tagged bGS share + %.3f (comparison only);gen #it{p}_{T}^{jet} [GeV];#it{f}_{#mu-tag}", gsShift));
  for(int i = 1; i <= NPt; i++){
    double tb = tagB->GetBinContent(i), a = fFC->GetBinContent(i), g = fGS->GetBinContent(i);
    if(tb <= 0. || a <= 0. || g <= 0.){ fShift->SetBinContent(i, 0.); fShift->SetBinError(i, 0.); continue; }
    double s = TMath::Min(1., tagGS->GetBinContent(i)/tb + gsShift);
    fShift->SetBinContent(i, 1./((1. - s)/a + s/g)); fShift->SetBinError(i, 0.);
  }

  // ---- correct the spectrum ---------------------------------------------------
  TFile *fI = TFile::Open(inPath);
  if(!fI || fI->IsZombie()){ printf("ERROR: cannot open %s -- run correctMuonRecoEff_caloJets_pp.C\n", inPath); return; }
  TH1D *hC = nullptr, *rStat = nullptr, *rIn[NIn];
  fI->GetObject("h_bJetPt_corrected", hC);
  fI->GetObject("h_corr_statRel", rStat);
  for(int k = 0; k < NIn; k++) fI->GetObject(Form("h_corr_sysRel_%s", inKey[k]), rIn[k]);
  if(!hC || !rStat){ printf("ERROR: missing histograms in %s\n", inPath); return; }
  for(int k = 0; k < NIn; k++) if(!rIn[k]){ printf("ERROR: missing h_corr_sysRel_%s\n", inKey[k]); return; }

  TH1D *hI = (TH1D*) hC->Clone("h_bJetPt_inclusive"); hI->SetDirectory(nullptr);
  hI->SetTitle("b jets (inclusive), unfolded, muon-reco and tag-frequency corrected, stat.;#it{p}_{T}^{jet} [GeV];d#it{N}/d#it{p}_{T} [GeV^{-1}]");
  TH1D *hIsys = (TH1D*) hI->Clone("h_bJetPt_inclusive_sysAbs"); hIsys->SetDirectory(nullptr);
  TH1D *iStat = (TH1D*) rStat->Clone("h_incl_statRel"); iStat->SetDirectory(nullptr);
  TH1D *iSrc[NIn]; for(int k = 0; k < NIn; k++){ iSrc[k] = (TH1D*) rIn[k]->Clone(Form("h_incl_sysRel_%s", inKey[k])); iSrc[k]->SetDirectory(nullptr); }
  TH1D *iF   = (TH1D*) rStat->Clone("h_incl_sysRel_muonTagFrequency"); iF->Reset(); iF->SetDirectory(nullptr);
  iF->SetTitle("muon-tag frequency (MC stat.);#it{p}_{T}^{jet} [GeV];relative uncertainty");
  TH1D *iTot = (TH1D*) iF->Clone("h_incl_sysRel"); iTot->SetTitle("total systematic;#it{p}_{T}^{jet} [GeV];relative uncertainty");

  printf("\n%-9s %8s %8s %8s %9s %9s | %11s %11s | %8s %8s %8s\n", "gen pT", "f", "f_FC", "f_GS", "GS share", "f_shift",
         "corrected", "inclusive", "f syst", "TOTAL", "stat");
  for(int i = 1; i <= NPt; i++){
    double f = fB->GetBinContent(i), c = hC->GetBinContent(i);
    double rf = f > 0. ? fB->GetBinError(i)/f : 0.;
    double v = (f > 0.) ? c/f : 0.;
    hI->SetBinContent(i, v); hI->SetBinError(i, f > 0. ? hC->GetBinError(i)/f : 0.);
    double q = rf*rf; for(int k = 0; k < NIn; k++) q += pow(rIn[k]->GetBinContent(i), 2);
    double t = (v > 0.) ? sqrt(q) : 0.;
    iF->SetBinContent(i, v > 0. ? rf : 0.); iTot->SetBinContent(i, t);
    hIsys->SetBinError(i, t*v);
    double tb = tagB->GetBinContent(i);
    printf("%3.0f-%-5.0f %8.4f %8.4f %8.4f %9.3f %9.4f | %11.4g %11.4g | %7.2f%% %7.2f%% %7.2f%%%s\n", ptEdges[i-1], ptEdges[i],
           f, fFC->GetBinContent(i), fGS->GetBinContent(i), tb > 0. ? tagGS->GetBinContent(i)/tb : 0., fShift->GetBinContent(i),
           c, v, 100*rf, 100*t, 100*iStat->GetBinContent(i), reported(i) ? "" : "  (underflow)");
  }

  TFile *fo = TFile::Open(outRoot, "recreate");
  fB->Write(); fFC->Write(); fGS->Write(); fShift->Write();
  if(fGenTag) fGenTag->Write();
  hI->Write(); hIsys->Write(); iStat->Write(); iTot->Write(); iF->Write();
  for(int k = 0; k < NIn; k++) iSrc[k]->Write();
  TNamed info("info", Form("pp calo b jets, inclusive: muon-SF-corrected unfolded spectrum / f, f (%s) = %s / all gen b jets "
                            "(|flavor| 5 + bGS 17) vs gen jet pT, PYTHIA natural bGS share, from %s. Syst adds f's MC stat (binomial). "
                            "Gen-muon-tag fraction stored as h_muTagFrequency_genTag; f with the purity fit's tagged-bGS shift (+%.3f) "
                            "stored as h_muTagFrequency_GSshift; neither used.",
                            fSource.Data(), fSource == "scanRecoTag" ? "analysis-reco-muon-tagged" : "muon-tagged",
                            gSystem->BaseName(mcPath), gsShift));
  info.Write();
  fo->Close();
  printf("\nwritten %s\n", outRoot);

  int colSys = TColor::GetColor("#D55E00");

  // ---- figure: f vs gen pT ----------------------------------------------------
  {
    TCanvas *c = new TCanvas("cF", "", 800, 700);
    c->SetLeftMargin(0.14); c->SetBottomMargin(0.13); c->SetRightMargin(0.04); c->SetTopMargin(0.06);
    double ymax = 0.;
    for(int i = 1; i <= NPt; i++) for(TH1D *h : {fB, fFC, fGS}) ymax = TMath::Max(ymax, h->GetBinContent(i) + h->GetBinError(i));
    if(fGenTag && fSource != "scanGenTag") for(int i = 1; i <= NPt; i++) ymax = TMath::Max(ymax, fGenTag->GetBinContent(i));
    TH1F *fr = c->DrawFrame(ptEdges[0], 0., ptEdges[NPt], 1.6*ymax);
    fr->GetXaxis()->SetTitle("gen #it{p}_{T}^{jet} [GeV]");
    fr->GetYaxis()->SetTitle(fSource == "scanGenTag" ? "fraction of b jets with a gen muon tag" : "fraction of b jets with a muon tag");
    fr->GetXaxis()->SetTitleSize(0.048); fr->GetYaxis()->SetTitleSize(0.048);
    fr->GetXaxis()->SetLabelSize(0.040); fr->GetYaxis()->SetLabelSize(0.040); fr->GetYaxis()->SetTitleOffset(1.35);
    drawUnderflowBand(0., 1.6*ymax, 0.05*ymax);
    TH1D *dB = (TH1D*) fB->Clone("dB"), *dF = (TH1D*) fFC->Clone("dF"), *dG = (TH1D*) fGS->Clone("dG"), *dS = (TH1D*) fShift->Clone("dS");
    styleH(dB, "#000000", markFilledCircle, 1.3);
    styleH(dF, "#0072B2", markOpenSquare, 1.2);
    styleH(dG, "#E69F00", markOpenDiamond, 1.5);
    styleH(dS, "#CC79A7", markCross, 1.4);
    for(int i = 1; i <= NPt; i++) if(dS->GetBinContent(i) <= 0.) dS->SetBinContent(i, -999.);
    if(haveFCGS){ dF->Draw("E1 X0 same"); dG->Draw("E1 X0 same"); dS->Draw("P same"); }
    TH1D *dGT = nullptr;
    if(fGenTag && fSource != "scanGenTag"){
      dGT = (TH1D*) fGenTag->Clone("dGT"); dGT->SetDirectory(nullptr);
      for(int i = 0; i <= NPt+1; i++) dGT->SetBinError(i, 0.);
      styleLine(dGT, "#999999", 2); dGT->SetLineStyle(2); dGT->Draw("HIST same");
    }
    dB->Draw("E1 X0 same");
    TLegend *leg = haveFCGS ? makeLegend(0.42, 0.18, 0.95, dGT ? 0.445 : 0.40, 0.034) : makeLegend(0.55, 0.18, 0.95, 0.25, 0.034);   // lower right is empty: f rises with pT
    leg->AddEntry(dB, "all b jets (used)", "lp");
    if(haveFCGS){
      leg->AddEntry(dF, "flavor creation", "lp");
      leg->AddEntry(dG, "gluon splitting", "lp");
      leg->AddEntry(dS, Form("tagged g#rightarrowb#bar{b} share + %.3f (not used)", gsShift), "p");
    }
    if(dGT) leg->AddEntry(dGT, fSource == "responseTagOverGen" ? "calo-matched scan (comp.)" : "all b, gen-muon tag (comparison)", "l");
    leg->Draw();
    TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.038);
    la.DrawLatex(0.18, 0.88, "PYTHIA, pp 5.02 TeV");
    la.SetTextSize(0.030);
    if(fSource == "scanRecoTag") la.DrawLatex(0.18, 0.83, "reco tight muon #it{p}_{T} > 15 GeV, |#eta| < 2, #DeltaR < 0.4 (gen jets)");
    else if(haveFCGS) la.DrawLatex(0.18, 0.83, "gen muon #it{p}_{T} > 15 GeV, #DeltaR < 0.4");
    else if(fSource == "responseTagOverGen") la.DrawLatex(0.18, 0.83, "muon-tagged / all b jets of the unfolding response");
    else {
      la.DrawLatex(0.18, 0.83, Form("PROVISIONAL: reco-#mu-tagged / all matched b jets, #div %.4f", effMuMC));
      la.DrawLatex(0.18, 0.79, "(gen-muon tag needs the fixed-scan rerun)");
    }
    savePdfTight(c, outFigF);
    delete c;
  }
  // the same for talks: all b jets only, reported range, no band
  {
    TCanvas *c = new TCanvas("cFs", "", 700, 800);   // the shape of the two spectrum figures it sits between
    c->SetLeftMargin(0.17); c->SetBottomMargin(0.12); c->SetRightMargin(0.04); c->SetTopMargin(0.05);
    double ymax = 0.;
    for(int i = 1; i <= NPt; i++) if(reported(i)) ymax = TMath::Max(ymax, fB->GetBinContent(i) + fB->GetBinError(i));
    TH1F *fr = c->DrawFrame(ptReportMin, 0., ptEdges[NPt], 1.5*ymax);
    fr->GetXaxis()->SetTitle("gen #it{p}_{T}^{jet} [GeV]"); fr->GetYaxis()->SetTitle("fraction of #it{b} jets with a muon tag");
    fr->GetXaxis()->SetTitleSize(0.050); fr->GetYaxis()->SetTitleSize(0.050);
    fr->GetXaxis()->SetLabelSize(0.042); fr->GetYaxis()->SetLabelSize(0.042); fr->GetYaxis()->SetTitleOffset(1.6);
    TH1D *dB = (TH1D*) fB->Clone("dBs");
    for(int i = 1; i <= NPt; i++) if(!reported(i)){ dB->SetBinContent(i, -999.); dB->SetBinError(i, 0.); }
    styleH(dB, "#000000", markFilledCircle, 1.4); dB->Draw("E1 X0 same");
    TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.046);
    la.DrawLatex(0.18, 0.88, "PYTHIA, pp");
    savePdfTight(c, TString(outFigF).ReplaceAll(".pdf", "_simple.pdf"));
    delete c;
  }

  // ---- figure: inclusive b-jet spectrum ---------------------------------------
  {
    TCanvas *c = new TCanvas("cI", "", 700, 800);
    TPad *pTop, *pBot; splitPads(pTop, pBot);
    pTop->cd(); pTop->SetLogy();   // log: falls by ~3 orders of magnitude, no negative bins
    TH1D *fr = (TH1D*) hI->Clone("frI"); fr->Reset(); fr->SetTitle("");
    double ymax = 0., ymin = 1e30;
    for(int i = 1; i <= NPt; i++) if(reported(i)){ ymax = TMath::Max(ymax, hI->GetBinContent(i)); ymin = TMath::Min(ymin, hI->GetBinContent(i)); }
    fr->SetMinimum(0.3*ymin); fr->SetMaximum(30.*ymax);
    fr->GetXaxis()->SetRangeUser(ptReportMin, ptEdges[NPt]); fr->GetXaxis()->SetLabelSize(0);
    fr->GetYaxis()->SetTitle("d#it{N}/d#it{p}_{T} [GeV^{-1}]"); fr->GetYaxis()->SetTitleSize(0.055);
    fr->GetYaxis()->SetTitleOffset(1.35); fr->GetYaxis()->SetLabelSize(0.045);
    fr->Draw("AXIS");
    for(int i = 1; i <= NPt; i++){
      double v = hI->GetBinContent(i); if(v <= 0. || !reported(i)) continue;
      TBox *bx = new TBox(ptEdges[i-1], v - hIsys->GetBinError(i), ptEdges[i], v + hIsys->GetBinError(i));
      bx->SetFillColorAlpha(colSys, 0.30); bx->SetLineColor(colSys); bx->Draw("l"); bx->Draw();
    }
    TH1D *d = (TH1D*) hI->Clone("dI");
    for(int i = 1; i <= NPt; i++) if(!reported(i)){ d->SetBinContent(i, -999.); d->SetBinError(i, 0.); }
    styleH(d, "#000000", markFilledCircle, 1.2);
    d->Draw("E1 X0 same");
    TLegend *leg = makeLegend(0.25, 0.05, 0.60, 0.20, 0.036);
    leg->AddEntry(d, "inclusive b jets, stat.", "lp");
    TBox *lb = new TBox(); lb->SetFillColorAlpha(colSys, 0.30); lb->SetLineColor(colSys);
    leg->AddEntry(lb, "total systematic", "f");
    leg->Draw();
    TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.046);
    la.DrawLatex(0.21, 0.84, Form("pp 5.02 TeV, %s: inclusive b jets", collLabel()));
    la.SetTextSize(0.034);
    la.DrawLatex(0.21, 0.785, "unfolded; muon reco and muon-tag frequency corrected");
    if(fSource == "responseRecoTag"){ la.SetTextColor(kRed+1); la.DrawLatex(0.21, 0.74, "tag frequency PROVISIONAL (reco-#mu tag / MC #varepsilon_{#mu})"); la.SetTextColor(kBlack); }

    pBot->cd();
    TH1D *rS = (TH1D*) iTot->Clone("rSI"), *rT = (TH1D*) iStat->Clone("rTI");
    for(TH1D *h : {rS, rT}){ h->Scale(100.); for(int i = 0; i <= NPt+1; i++) h->SetBinError(i, 0.); h->GetXaxis()->SetRangeUser(ptReportMin, ptEdges[NPt]); }
    styleLine(rS, "#D55E00", 3); styleLine(rT, "#0072B2", 3); rT->SetLineStyle(2);
    double m = 0.; for(int i = 1; i <= NPt; i++) if(reported(i)) m = TMath::Max(m, TMath::Max(rS->GetBinContent(i), rT->GetBinContent(i)));
    TH1D *rF = (TH1D*) rS->Clone("rFI"); rF->Reset();
    rF->SetMinimum(0.); rF->SetMaximum(1.7*m);   // headroom for the legend row
    rF->GetXaxis()->SetRangeUser(ptReportMin, ptEdges[NPt]);
    styleRatioAxes(rF, "#it{p}_{T}^{jet} [GeV]", "rel. unc. [%]");
    rF->Draw("AXIS");
    rS->Draw("HIST same"); rT->Draw("HIST same");
    TLegend *legB = makeLegend(0.19, 0.78, 0.95, 0.97, 0.075);
    legB->SetNColumns(2);
    legB->AddEntry(rS, "total syst.", "l"); legB->AddEntry(rT, "stat.", "l");
    legB->Draw();
    savePdfTight(c, outFig);
    delete c;
  }
  // the same for talks: spectrum with its band, shorter text
  {
    TCanvas *c = new TCanvas("cIs", "", 700, 800);
    TPad *pTop, *pBot; splitPads(pTop, pBot);
    pTop->cd(); pTop->SetLogy();
    TH1D *fr = (TH1D*) hI->Clone("frIs"); fr->Reset(); fr->SetTitle("");
    double ymax = 0., ymin = 1e30;
    for(int i = 1; i <= NPt; i++) if(reported(i)){ ymax = TMath::Max(ymax, hI->GetBinContent(i)); ymin = TMath::Min(ymin, hI->GetBinContent(i)); }
    fr->SetMinimum(0.3*ymin); fr->SetMaximum(30.*ymax);
    fr->GetXaxis()->SetRangeUser(ptReportMin, ptEdges[NPt]); fr->GetXaxis()->SetLabelSize(0);
    fr->GetYaxis()->SetTitle("d#it{N}/d#it{p}_{T} [GeV^{-1}]"); fr->GetYaxis()->SetTitleSize(0.055);
    fr->GetYaxis()->SetTitleOffset(1.35); fr->GetYaxis()->SetLabelSize(0.045);
    fr->Draw("AXIS");
    for(int i = 1; i <= NPt; i++){
      double v = hI->GetBinContent(i); if(v <= 0. || !reported(i)) continue;
      TBox *bx = new TBox(ptEdges[i-1], v - hIsys->GetBinError(i), ptEdges[i], v + hIsys->GetBinError(i));
      bx->SetFillColorAlpha(colSys, 0.30); bx->SetLineColor(colSys); bx->Draw();
    }
    TH1D *d = (TH1D*) hI->Clone("dIs");
    for(int i = 1; i <= NPt; i++) if(!reported(i)){ d->SetBinContent(i, -999.); d->SetBinError(i, 0.); }
    styleH(d, "#000000", markFilledCircle, 1.3); d->Draw("E1 X0 same");
    TLegend *leg = makeLegend(0.25, 0.08, 0.60, 0.24, 0.045);
    leg->AddEntry(d, "stat.", "lp");
    TBox *lb = new TBox(); lb->SetFillColorAlpha(colSys, 0.30); lb->SetLineColor(colSys);
    leg->AddEntry(lb, "total syst.", "f"); leg->Draw();
    TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.050);
    la.DrawLatex(0.21, 0.84, Form("pp, #it{b} jets, %s", collLabel()));
    pBot->cd();
    TH1D *rS = (TH1D*) iTot->Clone("rSIs"), *rT = (TH1D*) iStat->Clone("rTIs");
    for(TH1D *h : {rS, rT}){ h->Scale(100.); for(int i = 0; i <= NPt+1; i++) h->SetBinError(i, 0.); h->GetXaxis()->SetRangeUser(ptReportMin, ptEdges[NPt]); }
    styleLine(rS, "#D55E00", 3); styleLine(rT, "#0072B2", 3); rT->SetLineStyle(2);
    double m = 0.; for(int i = 1; i <= NPt; i++) if(reported(i)) m = TMath::Max(m, TMath::Max(rS->GetBinContent(i), rT->GetBinContent(i)));
    TH1D *rF = (TH1D*) rS->Clone("rFIs"); rF->Reset();
    rF->SetMinimum(0.); rF->SetMaximum(1.7*m);
    rF->GetXaxis()->SetRangeUser(ptReportMin, ptEdges[NPt]);
    styleRatioAxes(rF, "#it{p}_{T}^{jet} [GeV]", "rel. unc. [%]");
    rF->Draw("AXIS"); rS->Draw("HIST same"); rT->Draw("HIST same");
    TLegend *legB = makeLegend(0.19, 0.78, 0.95, 0.97, 0.075); legB->SetNColumns(2);
    legB->AddEntry(rS, "total syst.", "l"); legB->AddEntry(rT, "stat.", "l"); legB->Draw();
    savePdfTight(c, TString(outFig).ReplaceAll(".pdf", "_simple.pdf"));
    delete c;
  }

  // ---- figure: systematic breakdown -------------------------------------------
  {
    TCanvas *c = new TCanvas("cS", "", 800, 700);
    c->SetLeftMargin(0.13); c->SetBottomMargin(0.13); c->SetRightMargin(0.04); c->SetTopMargin(0.06);
    TH1D *lT = (TH1D*) iTot->Clone("lT"), *lF = (TH1D*) iF->Clone("lF"), *lS[NIn];
    for(int k = 0; k < NIn; k++) lS[k] = (TH1D*) iSrc[k]->Clone(Form("lS%d", k));
    std::vector<TH1D*> all = {lT, lF}; for(int k = 0; k < NIn; k++) all.push_back(lS[k]);
    for(TH1D *h : all){ h->Scale(100.); for(int i = 0; i <= NPt+1; i++) h->SetBinError(i, 0.); h->GetXaxis()->SetRangeUser(ptReportMin, ptEdges[NPt]); }
    double ymax = 0.; for(int i = 1; i <= NPt; i++) if(reported(i)) ymax = TMath::Max(ymax, lT->GetBinContent(i));
    TH1F *fr = c->DrawFrame(ptReportMin, 0., ptEdges[NPt], 2.1*ymax);   // headroom for header + legend
    fr->GetXaxis()->SetTitle("#it{p}_{T}^{jet} [GeV]");
    fr->GetYaxis()->SetTitle("relative systematic uncertainty [%]");
    fr->GetXaxis()->SetTitleSize(0.048); fr->GetYaxis()->SetTitleSize(0.048);
    fr->GetXaxis()->SetLabelSize(0.040); fr->GetYaxis()->SetLabelSize(0.040);
    TLegend *leg = makeLegend(0.40, 0.55, 0.95, 0.86, 0.030);
    // the 70-80 GeV turn-on correction stays in the total but is too small to earn a line
    for(int k = 0; k < NIn; k++){ if(TString(inKey[k]) == "turnOn") continue; styleLine(lS[k], inHex[k], 3); lS[k]->SetLineStyle(inStyle[k]); lS[k]->Draw("HIST same"); leg->AddEntry(lS[k], inLabel[k], "l"); }
    styleLine(lF, "#56B4E9", 3); lF->SetLineStyle(10); lF->Draw("HIST same"); leg->AddEntry(lF, "muon-tag frequency", "l");
    styleLine(lT, "#000000", 4); lT->Draw("HIST same"); leg->AddEntry(lT, "total", "l");
    leg->Draw();
    TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.038);
    la.DrawLatex(0.17, 0.885, Form("pp 5.02 TeV, %s: inclusive b-jet spectrum systematics", collLabel()));
    savePdfTight(c, outFigS);
    delete c;
  }
  printf("written %s\nwritten %s\nwritten %s\n", outFigF, outFig, outFigS);
}
