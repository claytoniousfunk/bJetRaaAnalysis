// Cross-check of our PbPb MinBias per-event jet spectrum against an
// independent scan produced by a colleague with separate code.
//
// This is the "ask Rhagu for per-event spectra with independent code" item
// from the 2026-08-18 meeting, run against whatever external file is handed
// over. It deliberately compares the simplest possible observable -- jets per
// event vs jet pT, no Z, no centrality split, no subtraction -- because that
// is the one thing two unrelated frameworks can be expected to agree on.
//
// *** READ THE CAVEATS BEFORE BELIEVING ANY DISAGREEMENT ***
// The external file inspected here (PbPb_Outfile_pTHat30p0_JetpT50p0_...)
// carries three problems that block a like-for-like comparison:
//
//   1. hCent is entirely zero. All 17181 entries sit in hiBin bin [0,1) with
//      RMS exactly 0, i.e. the centrality variable is never filled. hZvtx is
//      healthy (mean 0.24, RMS 4.9 cm, spanning +/-15), so event reading
//      works and it is centrality specifically that is broken. Consequence:
//      the external sample's centrality range is unknown, so it cannot be
//      matched to our 0-80% (hiBin < 160) selection.
//
//   2. The filename says JetpT50p0 but hJetpT is populated from 20 GeV, so
//      the filename parameters do not all correspond to applied cuts. It is
//      therefore not safe to assume the leading/subleading dijet cuts in the
//      name (LdpT50p0, SldpT30p0) were applied either -- and indeed they
//      cannot have been, since requiring a leading jet above 50 and a
//      subleading above 30 would force >= 2 jets/event, whereas the file has
//      1.43.
//
//   3. pTHat30p0 in the name suggests a generator-level cut, i.e. MC rather
//      than MinBias data. Not confirmable from the file contents.
//
// pT DEFINITION. The external scan uses RAW jet pT, so this comparison uses
// h_inclRawJetPt on our side, NOT h_inclRecoJetPt. That matters a lot here:
// our measured JEC factor <corrected>/<raw> is 1.15-1.19, so comparing our
// corrected spectrum against their raw one shifts us bodily to higher pT and,
// on a steeply falling spectrum, inflates our apparent per-event yield at any
// fixed pT. An earlier version of this macro made exactly that mistake and
// showed a spurious 1.4-1.7x excess.
//
// Known selection differences even setting those aside:
//   ours   |eta| < 1.6, raw jet pT > 20, jetTrkMax filter applied, hiBin < 160
//   theirs |eta| < 2.0 (per filename), raw jet pT > 20 (per contents), no
//          known trkMax filter, centrality range unknown
// The wider eta acceptance and absent trkMax filter would both push the
// external yield UP relative to ours. An unrestricted centrality range
// including 80-100% -- which we reject and which contributes almost no jets
// -- would push it DOWN.
//
// Usage: root -l -b -q 'plotJetPt_externalCrossCheck.C'
// Run from: src/plots/jetPt/

const char *ourFile =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/PbPb_MinBias_Part1_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_sameEventPFClustering_pseudoJetCandPtMin-0.0_2026-8-18_ultraFineCentBins.root";
const char *extFile =
  "/home/clayton/Downloads/PbPb_Outfile_pTHat30p0_JetpT50p0_LdpT50p0_SldpT30p0_JetEta2p0_LdEtaN1p0_1p0_0.root";

const char *extJetHist = "Event_Hist/hJetpT";
const char *extEvtHist = "Event_Hist/hCent";   // entry count = events passing

const char *outDir  = "../../../figures/jetPt/";
const char *outName = "jetPt_externalCrossCheck.pdf";

const double plotPtMin = 20., plotPtMax = 120.;

// Coarser than the native 5 GeV so the external file's thin high-pT tail
// stays visible rather than dissolving into single-count bins.
const int    NEdge = 14;
double       ptEdge[NEdge] = { 20,25,30,35,40,45,50,55,60,70,80,95,110,130 };

const char *ourHex = "#0072B2", *extHex = "#D55E00";

void plotJetPt_externalCrossCheck()
{
  gStyle->SetOptStat(0);
  gSystem->mkdir(outDir, kTRUE);

  TFile *fo = TFile::Open(ourFile);
  TFile *fe = TFile::Open(extFile);
  if(!fo || fo->IsZombie()){ printf("ERROR: cannot open %s\n", ourFile); return; }
  if(!fe || fe->IsZombie()){ printf("ERROR: cannot open %s\n", extFile); return; }

  // ---- ours: sum all 16 fine slices (0-80%), normalise by its own events ----
  TH1D *ourRaw = nullptr; double Nour = 0.;
  for(int si = 1; si <= 16; si++){
    TH1D *h = nullptr, *v = nullptr;
    fo->GetObject(Form("h_inclRawJetPt_C%d", si), h);   // RAW, to match the external scan
    fo->GetObject(Form("h_vz_C%d", si), v);
    if(!h || !v){ printf("ERROR: our slice %d missing\n", si); return; }
    if(!ourRaw){ ourRaw = (TH1D*) h->Clone("ourRaw"); ourRaw->SetDirectory(nullptr); }
    else ourRaw->Add(h);
    Nour += v->Integral();
  }

  // ---- theirs ----
  TH1D *extRaw = nullptr, *extEvt = nullptr;
  fe->GetObject(extJetHist, extRaw);
  fe->GetObject(extEvtHist, extEvt);
  if(!extRaw || !extEvt){ printf("ERROR: external histograms missing\n"); return; }
  double Next = extEvt->Integral();

  // centrality sanity check on the external file, reported not assumed
  printf("External centrality check: entries=%.0f mean=%.4f RMS=%.4f\n",
         extEvt->Integral(), extEvt->GetMean(), extEvt->GetRMS());
  if(extEvt->GetRMS() < 1e-6)
    printf("  *** hCent is degenerate (RMS = 0): centrality never filled.\n"
           "      The external centrality range is therefore UNKNOWN and cannot be\n"
           "      matched to our hiBin < 160 (0-80%%) selection. Any normalisation\n"
           "      difference below may simply be a different centrality mix. ***\n");

  TH1D *ourH = (TH1D*) ourRaw->Rebin(NEdge-1, "ourH", ptEdge);
  TH1D *extH = (TH1D*) extRaw->Rebin(NEdge-1, "extH", ptEdge);
  ourH->SetDirectory(nullptr); extH->SetDirectory(nullptr);
  ourH->Scale(1./Nour); ourH->Scale(1., "width");
  extH->Scale(1./Next); extH->Scale(1., "width");

  TH1D *ratio = (TH1D*) ourH->Clone("ratio");
  ratio->Divide(extH);

  printf("\nours   : N_evt=%10.0f  jets=%12.0f  jets/evt=%.4f\n", Nour, ourRaw->Integral(), ourRaw->Integral()/Nour);
  printf("theirs : N_evt=%10.0f  jets=%12.0f  jets/evt=%.4f\n", Next, extRaw->Integral(), extRaw->Integral()/Next);
  printf("  ratio of inclusive jets/event = %.3f\n", (ourRaw->Integral()/Nour)/(extRaw->Integral()/Next));
  printf("  ratio above 50 GeV            = %.3f\n",
         (ourRaw->Integral(ourRaw->FindBin(50.001),ourRaw->GetNbinsX()+1)/Nour) /
         (extRaw->Integral(extRaw->FindBin(50.001),extRaw->GetNbinsX()+1)/Next));

  printf("\n%-12s %14s %14s %9s\n", "pT [GeV]", "ours/evt/GeV", "theirs/evt/GeV", "ours/theirs");
  for(int b = 1; b <= ourH->GetNbinsX(); b++){
    double lo = ourH->GetXaxis()->GetBinLowEdge(b);
    if(lo < plotPtMin || lo >= plotPtMax) continue;
    printf("%4.0f-%-7.0f %14.6f %14.6f %9.3f\n", lo, ourH->GetXaxis()->GetBinUpEdge(b),
           ourH->GetBinContent(b), extH->GetBinContent(b), ratio->GetBinContent(b));
  }

  int colOur = TColor::GetColor(ourHex), colExt = TColor::GetColor(extHex);

  const double lm = 0.15, rm = 0.05, tm = 0.09, bm = 0.30, split = 0.36;
  TCanvas *c = new TCanvas("cExtXCheck", "", 850, 800);

  TPad *pUp = new TPad("pUp", "", 0, split, 1, 1);
  pUp->SetLeftMargin(lm); pUp->SetRightMargin(rm);
  pUp->SetTopMargin(tm);  pUp->SetBottomMargin(0.);
  pUp->SetLogy(); pUp->Draw(); pUp->cd();

  ourH->SetLineColor(colOur); ourH->SetMarkerColor(colOur);
  ourH->SetMarkerStyle(20); ourH->SetMarkerSize(1.0); ourH->SetLineWidth(2);
  ourH->SetTitle("");
  ourH->GetYaxis()->SetTitle("1/N_{evt} dN_{jet}/dp_{T}  [GeV^{-1}]");
  ourH->GetYaxis()->SetTitleSize(0.052); ourH->GetYaxis()->SetLabelSize(0.046);
  ourH->GetYaxis()->SetTitleOffset(1.35);
  ourH->GetXaxis()->SetRangeUser(plotPtMin, plotPtMax);
  ourH->SetMinimum(2e-5); ourH->SetMaximum(1.0);
  ourH->Draw("ep");

  extH->SetLineColor(colExt); extH->SetMarkerColor(colExt);
  extH->SetMarkerStyle(21); extH->SetMarkerSize(1.0); extH->SetLineWidth(2);
  extH->Draw("ep same");

  TLegend *leg = new TLegend(0.45, 0.62, 0.94, 0.86);
  leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.042);
  leg->AddEntry(ourH, "ours (PbPb MinBias, 0-80%, raw p_{T})", "lp");
  leg->AddEntry(extH, "external scan (raw p_{T})", "lp");
  leg->Draw();

  TLatex lat; lat.SetNDC(); lat.SetTextSize(0.045);
  lat.DrawLatex(lm + 0.02, 0.945, "PbPb per-event inclusive jet spectrum, raw p_{T}");

  c->cd();
  TPad *pDn = new TPad("pDn", "", 0, 0, 1, split);
  pDn->SetLeftMargin(lm); pDn->SetRightMargin(rm);
  pDn->SetTopMargin(0.);  pDn->SetBottomMargin(bm);
  pDn->Draw(); pDn->cd();

  const double sc = (1.-split)/split;
  ratio->SetLineColor(kBlack); ratio->SetMarkerColor(kBlack);
  ratio->SetMarkerStyle(20); ratio->SetMarkerSize(0.9); ratio->SetLineWidth(2);
  ratio->SetTitle("");
  ratio->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
  ratio->GetYaxis()->SetTitle("ours / external");
  ratio->GetXaxis()->SetTitleSize(0.052*sc); ratio->GetXaxis()->SetLabelSize(0.046*sc);
  ratio->GetYaxis()->SetTitleSize(0.052*sc); ratio->GetYaxis()->SetLabelSize(0.046*sc);
  ratio->GetYaxis()->SetTitleOffset(1.35/sc);
  ratio->GetYaxis()->SetNdivisions(505);
  ratio->GetXaxis()->SetRangeUser(plotPtMin, plotPtMax);
  ratio->SetMinimum(0.4); ratio->SetMaximum(2.1);
  ratio->Draw("ep");

  for(double yl : {1.0, 1.1, 0.9}){
    TLine *l = new TLine(plotPtMin, yl, plotPtMax, yl);
    l->SetLineStyle(yl == 1.0 ? 2 : 3); l->SetLineColor(yl == 1.0 ? kGray+2 : kGray+1);
    l->Draw();
  }

  TString out = TString(outDir) + outName;
  c->SaveAs(out);
  printf("\nSaved %s\n", out.Data());
}
