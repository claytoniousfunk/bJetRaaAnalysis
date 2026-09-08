// Task #3 from the 2026-08-18 meeting: fakeJet/event vs realJet/event at
// 50 GeV, both as a function of centrality in 5% slices.
//
// Olga: "at 50, where we are at, it is dominated by fakes.  The fakes should
// scale per event, so binary scattering will not be most important... You can
// plot your one jet pT bin, projected fakeJet/event, at exactly 50, as a
// function of centrality, along with the realJet/event at 50 as a fxn of
// centrality.  And these two should closely resemble each other for this
// small jet pT window."
//
// This is the discriminator for the flat-then-cliff shape found in
// plotJetPt_minBiasOnly_bumpRatioVsCent_fineCent.C: if the measured
// ("real") jet rate at 50 GeV tracks the fake-jet estimate across
// centrality, the fake description is consistent and the flatness is
// structural.  If they diverge -- especially in the most central slices --
// that isolates a genuine problem in either the fake estimate or the
// measured spectrum.
//
// Both curves come from the SAME mixed-event pfCandAnalyzer file, so N_evt
// (from h_vz) is common to numerator and denominator and cancels exactly in
// the ratio.
//   "real" = h_inclRecoJetPt / h_inclRawJetPt -- every reconstructed jet in
//            the event, which at 50 GeV in MinBias is largely fakes already
//   "fake" = h_fastJetPt_PF_[JEC_]bkgSub_RC -- the mixed-event, random-cone-
//            subtracted FastJet estimate, i.e. the quantity actually
//            subtracted downstream
//
// Shown for both pT definitions side by side, since task #1 established that
// the centrality degeneracy is present in raw pT as well as JEC-corrected pT.
// Note the fake spectrum's JEC and raw variants are genuinely different
// histograms (not a relabelling), so "50 GeV" means the same thing on each
// side of a given panel but not across panels.
//
// CENTRALITY CONVENTION: both files must come from the same side of the
// 2026-08-12 convention change (commit 7a2b47ad, "standardise the centrality
// binning convention").  Before it, getCentBin used (lo,hi] with the first
// bin inclusive at both ends, making slice 1 eleven hiBin units wide instead
// of ten; after it, [lo,hi) throughout.  The convention is not recorded in
// the output files, so mixing eras silently compares different hiBin ranges
// -- verified by integrating h_hiBin against h_vz per slice, which matches
// each file's own convention exactly (0/16 mismatches).  The mixed-event
// file here is 2026-8-14 rather than 2026-8-11 precisely so that both sides
// are post-change.
//
// Statistical errors on the ratio are combined in quadrature.  The two are
// not strictly independent -- the mixed-event pool is drawn from the same
// run -- but the pool aggregates candidates across many other events, so
// treating them as independent is a good approximation here.
//
// Usage: root -l -b -q 'plotFakeVsRealJet_at50GeV_vsCent.C'
// Run from: src/plots/jetPt/pseudoJets/

const char *mixedFile =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/PbPb_MinBias_Part1_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_mixedEventPFClustering_pseudoJetCandPtMin-0.0_2026-8-14_ultraFineCentBins.root";
// Same-event fake estimate, as the upper bracket: known to be contaminated by
// real jets (see plotFastJet_bkgSubtraction_mixedVsSame.C), so mixed and same
// bound the true fake rate from below and above respectively.
const char *sameFile =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/PbPb_MinBias_Part1_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_sameEventPFClustering_pseudoJetCandPtMin-0.0_2026-8-18_ultraFineCentBins.root";

const char *outDir  = "../../../../figures/jetPt/pseudoJets/";
const char *outName = "fakeVsRealJet_at50GeV_vsCent.pdf";

const int    NSlice   = 16;
const double centStep = 5.0;
const double probePt  = 52.5;   // centre of the 50-55 GeV bin

// Two pT definitions: {real-jet hist, fake-jet hist, label}
const int   NDef = 2;
const char *realHist[NDef] = { "h_inclRecoJetPt",              "h_inclRawJetPt" };
const char *fakeHist[NDef] = { "h_fastJetPt_PF_JEC_bkgSub_RC", "h_fastJetPt_PF_bkgSub_RC" };
const char *defLabel[NDef] = { "JEC-corrected p_{T}",          "raw p_{T}" };

const int   NClass = 4;
const int   sliceLo[NClass] = { 1,  3,  7, 11};
const int   sliceHi[NClass] = { 2,  6, 10, 16};
const char *classLabel[NClass] = { "0-10%", "10-30%", "30-50%", "50-80%" };

const char *realHex = "#0072B2", *fakeHex = "#D55E00", *fakeSameHex = "#009E73";

void plotFakeVsRealJet_at50GeV_vsCent()
{
  gStyle->SetOptStat(0);
  gSystem->mkdir(outDir, kTRUE);

  TFile *f = TFile::Open(mixedFile);
  if(!f || f->IsZombie()){ printf("ERROR: cannot open %s\n", mixedFile); return; }
  TFile *fs = TFile::Open(sameFile);
  if(!fs || fs->IsZombie()){ printf("ERROR: cannot open %s\n", sameFile); return; }

  int colReal = TColor::GetColor(realHex), colFake = TColor::GetColor(fakeHex);
  int colFakeSame = TColor::GetColor(fakeSameHex);

  TGraphErrors *gReal[NDef], *gFake[NDef], *gRatio[NDef];
  TGraphErrors *gFakeS[NDef], *gRatioS[NDef];
  double ymin_all = 1e300, ymax_all = 0., rmax_all = 0.;

  for(int di = 0; di < NDef; di++){
    gReal[di]   = new TGraphErrors();
    gFake[di]   = new TGraphErrors();
    gRatio[di]  = new TGraphErrors();
    gFakeS[di]  = new TGraphErrors();
    gRatioS[di] = new TGraphErrors();

    printf("\n=== %s ===\n", defLabel[di]);
    printf("%-9s %10s %13s %13s %10s %13s %10s\n", "cent", "N_evt", "real/evt", "fakeMix/evt", "mix/real", "fakeSame/evt", "same/real");

    double sR[NClass] = {0}, sF[NClass] = {0}, sN[NClass] = {0}, sFS[NClass] = {0};

    int np = 0;
    for(int si = 1; si <= NSlice; si++){
      TH1D *hr = nullptr, *hf = nullptr, *hv = nullptr, *hfs = nullptr, *hvs = nullptr;
      f->GetObject(Form("%s_C%d", realHist[di], si), hr);
      f->GetObject(Form("%s_C%d", fakeHist[di], si), hf);
      f->GetObject(Form("h_vz_C%d", si), hv);
      fs->GetObject(Form("%s_C%d", fakeHist[di], si), hfs);
      fs->GetObject(Form("h_vz_C%d", si), hvs);
      if(!hr || !hf || !hv || !hfs || !hvs){ printf("ERROR: missing hists for slice %d\n", si); return; }

      double N  = hv->Integral(), Ns = hvs->Integral();
      int    br = hr->FindBin(probePt), bf = hf->FindBin(probePt);
      double R  = hr->GetBinContent(br)/N, eR = hr->GetBinError(br)/N;
      double F  = hf->GetBinContent(bf)/N, eF = hf->GetBinError(bf)/N;
      double FS = hfs->GetBinContent(bf)/Ns, eFS = hfs->GetBinError(bf)/Ns;

      double cent = (si - 0.5)*centStep;
      printf("%3.0f-%-5.0f %10.0f %13.6f %13.6f %10.3f %13.6f %10.3f\n",
             (si-1)*centStep, si*centStep, N, R, F, R>0 ? F/R : -1., FS, R>0 ? FS/R : -1.);

      for(int c = 0; c < NClass; c++)
        if(si >= sliceLo[c] && si <= sliceHi[c]){ sR[c] += R*N; sF[c] += F*N; sFS[c] += FS*N; sN[c] += N; }

      if(R > 0.){
        gReal[di]->SetPoint(np, cent, R);
        gReal[di]->SetPointError(np, centStep/2., eR);
        if(R < ymin_all) ymin_all = R;
        if(R > ymax_all) ymax_all = R;
      }
      if(F > 0.){
        gFake[di]->SetPoint(np, cent, F);
        gFake[di]->SetPointError(np, centStep/2., eF);
        if(F < ymin_all) ymin_all = F;
        if(F > ymax_all) ymax_all = F;
      }
      if(FS > 0.){
        gFakeS[di]->SetPoint(np, cent, FS);
        gFakeS[di]->SetPointError(np, centStep/2., eFS);
        if(FS < ymin_all) ymin_all = FS;
        if(FS > ymax_all) ymax_all = FS;
      }
      if(R > 0. && F > 0.){
        double rr = F/R;
        double er = rr*TMath::Sqrt(TMath::Power(eF/F,2) + TMath::Power(eR/R,2));
        gRatio[di]->SetPoint(np, cent, rr);
        gRatio[di]->SetPointError(np, centStep/2., er);
        if(rr > rmax_all) rmax_all = rr;
      }
      if(R > 0. && FS > 0.){
        double rr = FS/R;
        double er = rr*TMath::Sqrt(TMath::Power(eFS/FS,2) + TMath::Power(eR/R,2));
        gRatioS[di]->SetPoint(np, cent, rr);
        gRatioS[di]->SetPointError(np, centStep/2., er);
        if(rr > rmax_all) rmax_all = rr;
      }
      np++;
    }

    printf("%-9s %10s %13s %13s %10s %13s %10s\n", "class", "", "real/evt", "fakeMix/evt", "mix/real", "fakeSame/evt", "same/real");
    for(int c = 0; c < NClass; c++)
      printf("%-9s %10s %13.6f %13.6f %10.3f %13.6f %10.3f\n",
             classLabel[c], "", sR[c]/sN[c], sF[c]/sN[c], sF[c]/sR[c], sFS[c]/sN[c], sFS[c]/sR[c]);
  }

  const double rMax = rmax_all*1.25;

  TCanvas *c = new TCanvas("cFakeVsReal50", "", 1300, 800);

  const double split = 0.40;
  const double sc = (1. - split)/split;

  for(int di = 0; di < NDef; di++){
    double x0 = di*0.5, x1 = (di+1)*0.5;

    c->cd();
    TPad *outer = new TPad(Form("outer_%d", di), "", x0, 0., x1, 1.);
    outer->SetMargin(0., 0., 0., 0.); outer->Draw(); outer->cd();

    TPad *pUp = new TPad(Form("pUp_%d", di), "", 0, split, 1, 1);
    pUp->SetLeftMargin(0.19); pUp->SetRightMargin(0.05);
    pUp->SetTopMargin(0.12);  pUp->SetBottomMargin(0.);
    pUp->SetLogy(); pUp->Draw();

    TPad *pDn = new TPad(Form("pDn_%d", di), "", 0, 0, 1, split);
    pDn->SetLeftMargin(0.19); pDn->SetRightMargin(0.05);
    pDn->SetTopMargin(0.);    pDn->SetBottomMargin(0.28);
    pDn->Draw();

    pUp->cd();
    TH1F *frUp = pUp->DrawFrame(0., ymin_all*0.4, 80., ymax_all*4.);
    frUp->GetYaxis()->SetTitle("jets / event at 50-55 GeV");
    frUp->GetYaxis()->SetTitleSize(0.058); frUp->GetYaxis()->SetLabelSize(0.050);
    frUp->GetYaxis()->SetTitleOffset(1.35);
    frUp->GetXaxis()->SetLabelSize(0.);
    for(double xb : {10., 30., 50.}){
      TLine *l = new TLine(xb, ymin_all*0.4, xb, ymax_all*4.);
      l->SetLineStyle(2); l->SetLineColor(kGray+1); l->Draw();
    }

    gReal[di]->SetMarkerStyle(20); gReal[di]->SetMarkerColor(colReal);
    gReal[di]->SetLineColor(colReal); gReal[di]->SetMarkerSize(1.0);
    gReal[di]->SetLineWidth(2); gReal[di]->Draw("pz same");

    gFakeS[di]->SetMarkerStyle(33); gFakeS[di]->SetMarkerColor(colFakeSame);
    gFakeS[di]->SetLineColor(colFakeSame); gFakeS[di]->SetMarkerSize(1.3);
    //gFakeS[di]->SetLineWidth(2); gFakeS[di]->Draw("pz same");

    gFake[di]->SetMarkerStyle(21); gFake[di]->SetMarkerColor(colFake);
    gFake[di]->SetLineColor(colFake); gFake[di]->SetMarkerSize(1.0);
    gFake[di]->SetLineWidth(2); gFake[di]->Draw("pz same");

    TLatex dl; dl.SetNDC(); dl.SetTextSize(0.062); dl.SetTextFont(62); dl.SetTextAlign(13);
    dl.DrawLatex(0.22, 0.96, defLabel[di]);

    if(di == 0){
      TLegend *leg = new TLegend(0.40, 0.11, 0.95, 0.36);
      leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.049);
      leg->AddEntry(gReal[di],  "real (all reco jets)", "lp");
      //leg->AddEntry(gFakeS[di], "fake, same-event", "lp");
      leg->AddEntry(gFake[di],  "fake, mixed-event", "lp");
      leg->Draw();
    }

    pDn->cd();
    TH1F *frDn = pDn->DrawFrame(0., 0., 80., rMax);
    frDn->GetXaxis()->SetTitle("centrality [%]");
    frDn->GetYaxis()->SetTitle("fake / real");
    frDn->GetXaxis()->SetTitleSize(0.058*sc); frDn->GetXaxis()->SetLabelSize(0.050*sc);
    frDn->GetYaxis()->SetTitleSize(0.058*sc); frDn->GetYaxis()->SetLabelSize(0.050*sc);
    frDn->GetYaxis()->SetTitleOffset(1.35/sc);
    frDn->GetYaxis()->SetNdivisions(505);
    for(double xb : {10., 30., 50.}){
      TLine *l = new TLine(xb, 0., xb, rMax);
      l->SetLineStyle(2); l->SetLineColor(kGray+1); l->Draw();
    }
    TLine *one = new TLine(0., 1., 80., 1.);
    one->SetLineStyle(2); one->SetLineColor(kGray+2); one->Draw();

    gRatioS[di]->SetMarkerStyle(33); gRatioS[di]->SetMarkerColor(colFakeSame);
    gRatioS[di]->SetLineColor(colFakeSame); gRatioS[di]->SetMarkerSize(1.1);
    //gRatioS[di]->SetLineWidth(2); gRatioS[di]->Draw("pz same");

    gRatio[di]->SetMarkerStyle(21); gRatio[di]->SetMarkerColor(colFake);
    gRatio[di]->SetLineColor(colFake); gRatio[di]->SetMarkerSize(0.9);
    gRatio[di]->SetLineWidth(2); gRatio[di]->Draw("pz same");
  }

  c->cd();
  TLatex lat; lat.SetNDC(); lat.SetTextSize(0.023);
  lat.DrawLatex(0.05, 0.982, "PbPb MinBias, per-event, at 50-55 GeV -- measured jets vs mixed-event fake-jet estimate");

  TString out = TString(outDir) + outName;
  c->SaveAs(out);
  printf("\nSaved %s\n", out.Data());
}
