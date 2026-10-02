// PYTHIA reco jet pT spectrum with the JEU shifted up and down, as shifted / nominal.
//
// Reads h_inclRecoJetPt[_inclRecoMuonTag]_flavor_{nominal,JEUShiftUp,JEUShiftDown},
// booked in PYTHIA_scan.C (the JEU variants were added 2026-10-02; scans before
// that do not have them). Four figures: {all jets, b jets} x {inclusive,
// muon-tagged}. Flavor is the signed JetFlavorID; b jets are |flavor| = 5.
//
// The shift is applied to the jet pT only: the jetPtCut, jetTrkMax filter and
// weights all follow the NOMINAL pT, so below ~jetPtCut/(1-JEU) the down ratio
// is not a clean shift of a complete spectrum. Read the lowest bins with that in
// mind. The uncertainty file is the Spring18 ppRef AK4PF one, also used for calo
// jets; there is no calo-specific one in the repo.
//
// The ratio errors are the numerator's alone (RatioErr::kNumerator): the shifted
// and nominal spectra are built from the SAME jets.
//
// Usage: root -l -b -q 'plotPYTHIAJetPt_JEUShiftRatio.C("<scan.root>")'
// Run from: src/plots/jetPt/jetEnergyUncertainty/

#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TSystem.h"
#include "../../../../headers/functions/divideByBinwidth.h"
#include "../../../../headers/plotting/plotStyle.h"
#include "../../../../headers/plotting/ratioPanel.h"

namespace {

const char *outDir = "../../../../figures/jetKinematics/jetEnergyUncertainty/";

const double ptLo = 30, ptHi = 400;
const int    ptRebin = 4;   // 5 GeV -> 20 GeV bins

// flavor axis of the *_flavor histograms: 27 bins on [-5,22), signed, so a b jet
// is flavor +5 OR -5 (both populated: the bbar bin holds as many as the b bin)
TH1D* spectrum(TFile *f, const char *base, const char *var, bool bOnly, const char *name)
{
  auto *H = (TH2D*) f->Get(Form("%s_%s", base, var));
  if(!H){ printf("missing %s_%s -- scan predates the JEU histograms?\n", base, var); return nullptr; }
  TH1D *h;
  if(bOnly){
    int bm = H->GetYaxis()->FindBin(-4.5), bp = H->GetYaxis()->FindBin(5.5);
    h = H->ProjectionX(name, bm, bm);
    TH1D *hp = H->ProjectionX(Form("%s_bbar", name), bp, bp);
    h->Add(hp); delete hp;
  }
  else h = H->ProjectionX(name, 1, H->GetNbinsY());
  h->SetDirectory(nullptr);
  h->Rebin(ptRebin);
  return h;
}

void oneFigure(TFile *f, const char *base, bool bOnly, const char *label, const char *tag)
{
  TH1D *nom = spectrum(f, base, "nominal",     bOnly, Form("nom_%s", tag));
  TH1D *up  = spectrum(f, base, "JEUShiftUp",  bOnly, Form("up_%s",  tag));
  TH1D *dn  = spectrum(f, base, "JEUShiftDown",bOnly, Form("dn_%s",  tag));
  if(!nom || !up || !dn) return;

  // ratios from counts, then densities for the spectrum panel
  TH1D *rUp = makeRatio(up, nom, Form("rUp_%s", tag));
  TH1D *rDn = makeRatio(dn, nom, Form("rDn_%s", tag));
  divideByBinwidth(nom); divideByBinwidth(up); divideByBinwidth(dn);

  styleH(nom, okabeHex[0], markFilledCircle, 1.0);
  styleH(up,  okabeHex[6], markFilledSquare, 1.0);
  styleH(dn,  okabeHex[5], markFilledDiamond, 1.2);
  styleH(rUp, okabeHex[6], markFilledSquare, 1.0);
  styleH(rDn, okabeHex[5], markFilledDiamond, 1.2);

  TCanvas *c = new TCanvas(Form("c_%s", tag), "", 700, 800);
  TPad *top, *bot; splitPads(top, bot);

  top->cd(); top->SetLogy();   // steeply falling spectrum; no negative bins to hide here
  nom->GetXaxis()->SetRangeUser(ptLo, ptHi);
  nom->SetTitle("");
  nom->GetYaxis()->SetTitle("dN/dp_{T} (arb.)");
  nom->GetYaxis()->SetTitleSize(0.055); nom->GetYaxis()->SetLabelSize(0.050);
  nom->GetXaxis()->SetLabelSize(0);
  nom->Draw("E1");
  up->Draw("E1 SAME"); dn->Draw("E1 SAME");
  TLegend *leg = makeLegend(0.55, 0.66, 0.90, 0.88, 0.048);
  leg->AddEntry(nom, "Nominal", "lp");
  leg->AddEntry(up,  "JEU up", "lp");
  leg->AddEntry(dn,  "JEU down", "lp");
  leg->Draw();
  TLatex t; t.SetNDC(); t.SetTextSize(0.05);
  t.DrawLatex(0.20, 0.30, "PYTHIA8 calo jets");
  t.DrawLatex(0.20, 0.23, label);

  bot->cd();
  rUp->GetXaxis()->SetRangeUser(ptLo, ptHi);
  styleRatioAxes(rUp, "Reco jet p_{T} (GeV)", "Shifted / nominal");
  rUp->GetYaxis()->SetRangeUser(0.6, 1.5);
  rUp->Draw("E1"); rDn->Draw("E1 SAME");
  unityLine(ptLo, ptHi)->Draw();

  gSystem->mkdir(outDir, kTRUE);
  savePdfTight(c, Form("%sPYTHIAJetPt_JEUShiftRatio_%s.pdf", outDir, tag));
  delete c;
}

} // namespace

void plotPYTHIAJetPt_JEUShiftRatio(const char *scanFile)
{
  initPlotStyle();
  TFile *f = TFile::Open(scanFile);
  if(!f || f->IsZombie()){ printf("cannot open %s\n", scanFile); return; }

  oneFigure(f, "h_inclRecoJetPt_flavor",                false, "Inclusive jets",         "inclusive_allJets");
  oneFigure(f, "h_inclRecoJetPt_flavor",                true,  "Inclusive b jets",       "inclusive_bJets");
  oneFigure(f, "h_inclRecoJetPt_inclRecoMuonTag_flavor",false, "#mu-tagged jets",        "muTagged_allJets");
  oneFigure(f, "h_inclRecoJetPt_inclRecoMuonTag_flavor",true,  "#mu-tagged b jets",      "muTagged_bJets");
}
