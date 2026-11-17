// Systematic uncertainties by source of the pp PF b-jet spectrum after each step of the pp b-jet chain, in one
// 1x4 row (measured, unfolded, muon-efficiency corrected, tag-frequency corrected), 80-200 GeV (2026-11-16, user:
// remake the systematics plots in the style of the b-purity one). Drawn from the four steps' output files, nothing
// recomputed; replaces the four single-panel step figures the chain macros write (bJetSpectrum_sysBreakdown_*,
// unfoldedBJetSpectrum_sysBreakdown_*, bJetSpectrum_muRecoCorr_sysBreakdown_*, bJetSpectrum_inclusive_sysBreakdown_*),
// which stay as they are. Style: drawSysRow() of headers/plotting/sysStyle.h.
//   b purity = all purity sources together (the step-1 file also stores its parts, not drawn here)
//
// inputs  rootFiles/MuTagMuTrigBJetSpectra/Data/bJetSpectrum_PFJets_pp<var>.root            (h_bJetPt_sysRel_*)
//         rootFiles/CorrectedBJetSpectra/Data/unfoldedBJetSpectrum_PFJets_pp<var>.root       (h_unf_sysRel_*)
//         rootFiles/CorrectedBJetSpectra/Data/bJetSpectrum_muRecoCorr_PFJets_pp<var>.root    (h_corr_sysRel_*)
//         rootFiles/CorrectedBJetSpectra/Data/bJetSpectrum_inclusive_PFJets_pp<var>.root     (h_incl_sysRel_*)
//         <var> = purity-variant suffix, e.g. "_csvV2_2D" (the nominal)
// Usage, from src/calculateBJetsPerZ/:  root -l -b -q 'plotBJetChainSys_pp_PF.C("_csvV2_2D")'
// Output: figures/systematics/bJetChainSys_PFJets_pp<var>.pdf
#include "TFile.h"
#include "TH1.h"
#include "TSystem.h"
#include "../../headers/plotting/sysStyle.h"

void plotBJetChainSys_pp_PF(TString var = "_csvV2_2D")
{
  initPlotStyle();
  const TString base = "../../rootFiles/";
  struct Step { const char *name, *file, *pre; };
  const Step steps[4] = {
    {"measured",         "MuTagMuTrigBJetSpectra/Data/bJetSpectrum_PFJets_pp",            "h_bJetPt_sysRel"},
    {"unfolded",         "CorrectedBJetSpectra/Data/unfoldedBJetSpectrum_PFJets_pp",      "h_unf_sysRel"},
    {"#mu eff. corr.",   "CorrectedBJetSpectra/Data/bJetSpectrum_muRecoCorr_PFJets_pp",   "h_corr_sysRel"},
    {"tag freq. corr.",  "CorrectedBJetSpectra/Data/bJetSpectrum_inclusive_PFJets_pp",    "h_incl_sysRel"},
  };
  // stored name -> style key (sysStyle.h) and label
  const int NS = 10;
  const char *stored[NS] = {"spectrumJEU", "bPurity", "iterations", "unfoldingBias", "responseStat", "prior", "recoFloor", "turnOn", "muonRecoEff", "muonTagFrequency"};
  const char *keys[NS]   = {"JEU",         "bPurity", "iterations", "unfoldingBias", "responseStat", "prior", "recoFloor", "turnOn", "muonEff",     "muonTagFreq"};
  const char *labs[NS]   = {nullptr,       "b purity", nullptr,     nullptr,         nullptr,        nullptr, nullptr,     "trigger turn-on", nullptr, nullptr};
  SysPanel P[4];
  for(int i = 0; i < 4; i++){
    TFile *f = TFile::Open(base + steps[i].file + var + ".root");
    if(!f || f->IsZombie()){ printf("ERROR: %s%s.root missing\n", steps[i].file, var.Data()); return; }
    P[i].name = steps[i].name;
    TH1 *t = nullptr; f->GetObject(steps[i].pre, t); P[i].tot = t;
    if(!t){ printf("ERROR: %s missing in step %d\n", steps[i].pre, i); return; }
    for(int s = 0; s < NS; s++){ TH1 *h = nullptr; f->GetObject(Form("%s_%s", steps[i].pre, stored[s]), h); P[i].src[s] = h; }
  }
  const TString dir = "../../figures/systematics/"; gSystem->mkdir(dir, true);
  drawSysRow(4, P, NS, keys, labs, 80., 200., "pp #it{b} jets, chain steps",
             "anti-#it{k}_{T} #it{R} = 0.4 PF jets,  #sqrt{#it{s}} = 5.02 TeV", dir + "bJetChainSys_PFJets_pp" + var + ".pdf");
}
