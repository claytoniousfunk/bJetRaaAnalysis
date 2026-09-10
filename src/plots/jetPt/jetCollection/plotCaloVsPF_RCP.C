// RCP for akPu4Calo vs akCs4PF, on the high-statistics MinBias scans.
//
// Supersedes the calo/flow/PF comparison in plotJetsByCollection.cc for this
// purpose: that one runs on ~185k-event noRhoModification scans in 4CentBins,
// these are 30M (calo) and 37M (PF) events in ultraFine bins. Flow is dropped
// here; add it back by pointing a third file at the same machinery.
//
// PT FLOOR = 70 GeV, and this is not arbitrary. The calo scan's RAW jets start
// at 35-40 GeV, but the calo JEC is large at low pT (~1.7x at 25 GeV against
// ~1.14x for PF), so on the corrected scale the first populated calo bin is
// 65-70. Below that calo is structurally empty while PF is not, and a ratio
// there compares a real spectrum against nothing. PF alone reaches down to 20.
//
// RCP = (class per-event yield) / (50-80% per-event yield), each collection
// against its OWN peripheral bin. No N_coll scaling, so the absolute value is
// not a suppression measurement -- what is meaningful is calo vs PF at fixed
// centrality, and the pT dependence.
//
// Two variants per class:
//   caloVsPF_RCP_<class>.pdf         as measured
//   caloVsPF_RCPpinned_<class>.pdf   each curve divided by its own value in the
//                                    highest pT bin, so both pass through 1
//                                    there and only the pT dependence remains
//
// Usage: root -l -b -q 'plotCaloVsPF_RCP.C'
// Run from: src/plots/jetPt/jetCollection/

#include "../../../../headers/plotting/coarseCent.h"
#include "../../../../headers/plotting/plotStyle.h"
#include "../../../../headers/plotting/ratioPanel.h"
#include "../../../../headers/functions/divideByBinwidth.h"

const char *caloFile =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/"
  "PbPb_MinBias_Part1_caloJets_mu12_pTmu-15to999_tight_jetTrkMaxFilter_"
  "WDecayFilter_2026-9-10_ultraFineCentBins.root";
const char *pfFile =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/"
  "PbPb_MinBias_Part1_mu12_pTmu-15to999_tight_jetTrkMaxFilter_"
  "WDecayFilter_2026-8-13_ultraFineCentBins.root";

const char *histBase = "h_inclRecoJetPt";   // JEC-corrected inclusive reco jets
const char *outDir   = "../../../../figures/jetCollection/";

// 70 GeV floor, widening with pT to keep the 50-80% denominator populated --
// it is the thinnest sample here (2247 calo jets in 60-80 against 43052 in
// 0-10%), and RCP is only as good as its denominator.
const double edgePt[]  = {70, 85, 100, 130, 180, 300};
const int    nEdgePt   = (int)(sizeof(edgePt)/sizeof(double)) - 1;
const double pinAtPt   = 180.;   // pin in the last bin, 180-300

const int PERIPH = 3;   // index of the 50-80% class, the RCP denominator

// per-event, bin-width-normalised spectrum for one class
TH1D* spectrum(TFile *f, int ci, const char *tag)
{
  TH1D *s = coarseSum1(f, histBase, ci, tag);
  if(!s) return nullptr;
  double n = coarseEvents(f, ci, "h_vz");   // these fills are not trigger-gated
  if(n <= 0.){ delete s; return nullptr; }
  TH1D *r = rebinTo(s, nEdgePt, edgePt, Form("sp_%s_%d", tag, ci));
  delete s;
  r->Scale(1./n);
  divideByBinwidth(r);
  return r;
}

void plotCaloVsPF_RCP()
{
  initPlotStyle();
  gSystem->mkdir(outDir, kTRUE);

  TFile *fc = TFile::Open(caloFile);
  TFile *fp = TFile::Open(pfFile);
  if(!fc || fc->IsZombie()){ printf("ERROR: cannot open calo file\n"); return; }
  if(!fp || fp->IsZombie()){ printf("ERROR: cannot open PF file\n");   return; }

  // peripheral reference, one per collection
  TH1D *pC = spectrum(fc, PERIPH, "cRef");
  TH1D *pP = spectrum(fp, PERIPH, "pRef");
  if(!pC || !pP){ printf("ERROR: missing 50-80%% reference\n"); return; }

  printf("RCP = class / %s, per-event, no N_coll scaling\n", coarseLabel[PERIPH]);
  printf("%-9s %-6s", "class", "coll");
  for(int b = 0; b < nEdgePt; b++) printf(" %11.0f-%-4.0f", edgePt[b], edgePt[b+1]);
  printf("\n");

  for(int ci = 0; ci < PERIPH; ci++){

    TH1D *hC = spectrum(fc, ci, Form("c%d", ci));
    TH1D *hP = spectrum(fp, ci, Form("p%d", ci));
    if(!hC || !hP){ printf("  %s: missing, skipped\n", coarseLabel[ci]); continue; }

    // errors on both: numerator and denominator are different centrality
    // classes, so they are independent samples -- unlike the injection study,
    // where the two collections shared events.
    TH1D *rC = makeRatio(hC, pC, Form("rcpC_%d", ci), RatioErr::kBoth);
    TH1D *rP = makeRatio(hP, pP, Form("rcpP_%d", ci), RatioErr::kBoth);
    styleH(rC, hexBkg2, markFilledSquare);    // calo
    styleH(rP, hexData, markFilledCircle);    // PF

    for(int k = 0; k < 2; k++){
      TH1D *r = k ? rP : rC;
      printf("%-9s %-6s", coarseLabel[ci], k ? "PF" : "calo");
      for(int b = 1; b <= r->GetNbinsX(); b++) printf(" %16.3f", r->GetBinContent(b));
      printf("\n");
    }

    // ---- absolute ----------------------------------------------------------
    for(int pass = 0; pass < 2; pass++){

      TH1D *dC = (TH1D*) rC->Clone(Form("dC_%d_%d", ci, pass)); dC->SetDirectory(nullptr);
      TH1D *dP = (TH1D*) rP->Clone(Form("dP_%d_%d", ci, pass)); dP->SetDirectory(nullptr);

      if(pass == 1){
        // pin: divide each by its OWN value in the top bin
        int bPin = dC->FindBin(pinAtPt + 1e-6);
        double vC = dC->GetBinContent(bPin), vP = dP->GetBinContent(bPin);
        if(vC <= 0. || vP <= 0.){
          printf("  %s: no yield in the pinning bin, pinned variant skipped\n", coarseLabel[ci]);
          continue;
        }
        dC->Scale(1./vC); dP->Scale(1./vP);
      }

      TCanvas *c = new TCanvas(Form("c_%d_%d", ci, pass), "", 700, 700);
      TPad *pad = new TPad("pad","",0,0,1,1);
      pad->SetLeftMargin(0.19); pad->SetBottomMargin(0.13);
      pad->Draw(); pad->cd();

      double ymax = 0.;
      for(int b = 1; b <= dC->GetNbinsX(); b++){
        ymax = TMath::Max(ymax, dC->GetBinContent(b) + dC->GetBinError(b));
        ymax = TMath::Max(ymax, dP->GetBinContent(b) + dP->GetBinError(b));
      }
      dC->SetTitle("");
      dC->GetXaxis()->SetTitle("Jet #it{p}_{T} [GeV]");
      dC->GetXaxis()->SetTitleSize(0.045); dC->GetXaxis()->SetLabelSize(0.038);
      dC->GetYaxis()->SetTitle(pass ? Form("(%s / %s), pinned at %.0f GeV",
                                           coarseLabel[ci], coarseLabel[PERIPH], pinAtPt)
                                    : Form("%s / %s", coarseLabel[ci], coarseLabel[PERIPH]));
      dC->GetYaxis()->SetTitleSize(0.045); dC->GetYaxis()->SetTitleOffset(1.85);
      dC->GetYaxis()->SetLabelSize(0.038);
      dC->SetMinimum(0.); dC->SetMaximum(ymax*1.45);
      dC->Draw("E");
      dP->Draw("E same");

      if(pass == 1){
        TLine *l = unityLine(edgePt[0], edgePt[nEdgePt]);
        l->Draw();
        dC->Draw("E same"); dP->Draw("E same");
      }

      TLegend *leg = makeLegend(0.55, 0.74, 0.93, 0.88, 0.036);
      leg->AddEntry(dP, "akCs4PF", "lp");
      leg->AddEntry(dC, "akPu4Calo", "lp");
      leg->Draw();

      TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.036);
      la.DrawLatex(0.22, 0.92, Form("PbPb MinBias, %s / %s%s",
                                    coarseLabel[ci], coarseLabel[PERIPH],
                                    pass ? "" : "   (no N_{coll} scaling)"));

      c->SaveAs(Form("%scaloVsPF_%s_%s.pdf", outDir,
                     pass ? "RCPpinned" : "RCP", coarseTag[ci]));
      delete c;
    }
  }

  printf("\nfigures written to %s\n", outDir);
}
