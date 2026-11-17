// b-purity systematics of the CSVv2-cut (or, with "_2D", the nominal 2D-fit) purity by source, pp and the four PbPb classes
// (pp, 50-80 ... 0-10%), drawn from the stored purity files without refitting (2026-11-05).
// Same layout, sources and styles as the figure of systematics_CSVv2Cut_PF.C, shown over the
// reported 80-200 GeV only (user: the 200-300 GeV window, 98% in 0-10%, set the y scale).
//
// Inputs: rootFiles/systematics/bPuritySys_total_PFJets_pp_csvV2<wp>.root and
//         bPuritySys_PFJets_PbPb_csvV2<wp>.root (h_sysRel_<src>[_C<i>]), from systematics_CSVv2Cut_PF.C
// Usage, from src/systematics/:  root -l -b -q 'plotBPuritySys_CSVv2_PF.C("loose")'   ("medium", "_2D")
// Drawing (style of 2026-11-16) in headers/plotting/puritySysFigure.h, shared with systematics_CSVv2_2D_PF.C.
// Writes figures/systematics/bPurity/CSVv2/bPuritySys_CSVv2_<wp>.pdf ("_2D": bPuritySys_CSVv2_2D.pdf)

#include "TFile.h"
#include "TH1D.h"
#include "TSystem.h"
#include "../../headers/plotting/puritySysFigure.h"

void plotBPuritySys_CSVv2_PF(const char *wpN = "loose")
{
  initPlotStyle();
  const TString base = "../../rootFiles/systematics/", w = wpN;
  const bool is2D = w.BeginsWith("_");                      // "_2D": the nominal 2D fit (systematics_CSVv2_2D_PF.C)
  TFile *fp = TFile::Open(base + Form("bPuritySys_total_PFJets_pp_csvV2%s.root", wpN)), *fa = TFile::Open(base + Form("bPuritySys_PFJets_PbPb_csvV2%s.root", wpN));
  if(!fp || !fa || fp->IsZombie() || fa->IsZombie()){ printf("ERROR: purity systematics files for %s missing\n", wpN); return; }
  const double wp = w == "medium" ? 0.84 : 0.54;
  const int NSRC = 7;
  const char *key[NSRC] = {"bGS", "cMult", "fitRangeHigh", "tplStat", "JEU", "JER", "fakeMuJet"};
  const char *sfx[5] = {"", "_C4", "_C3", "_C2", "_C1"};   // pp, 50-80 ... 0-10%
  TH1D *tot[5], *src[5][PSF_NMAX] = {{nullptr}};
  for(int ix = 0; ix < 5; ix++){
    TFile *f = ix ? fa : fp;
    tot[ix] = nullptr; f->GetObject(Form("h_sysRel_total%s", sfx[ix]), tot[ix]);
    if(!tot[ix]){ printf("ERROR: h_sysRel_total%s missing\n", sfx[ix]); return; }
    for(int s = 0; s < NSRC; s++){
      if(!is2D && s == 4 && ix != 0) continue;   // as systematics_CSVv2Cut_PF.C: JES drawn for pp only
      f->GetObject(Form("h_sysRel_%s%s", key[s], sfx[ix]), src[ix][s]);
    }
  }
  // 7-bin variant (2026-11-16, user): its JES and JER are taken from the nominal 80-100-120-150-200 GeV windows
  // (studyJetPtBinning_sys_2D_PF.C, _mergedJESJER), so neighbouring fine bins repeat them: labelled as such
  const bool fine = w.Contains("7bins");
  const char *labFine[NSRC] = {nullptr, nullptr, nullptr, nullptr, "JES (nominal bins)", "JER (nominal bins)", nullptr};
  const TString dir = "../../figures/systematics/bPurity/CSVv2/"; gSystem->mkdir(dir, true);
  drawPuritySysRow(tot, src, NSRC, key, 80., 200.,
                   is2D ? "muon #it{p}_{T}^{rel} #times CSVv2 fit" : Form("CSVv2 > %.2f, muon #it{p}_{T}^{rel} fit", wp),
                   dir + (is2D ? Form("bPuritySys_CSVv2%s.pdf", wpN) : Form("bPuritySys_CSVv2_%s.pdf", wpN)),
                   nullptr, false);   // plain legend labels (2026-11-16, user: no "(nominal bins)"; the caption says it)   // no bin-edge lines (2026-11-16, user: removed from the AN figure)
}
