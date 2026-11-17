// Is the b-purity JER (and JES) a resolution effect or the statistical structure of the template ratio?
// (2026-11-17, user: "revisit purity JER", after the 0-10% 80-100 GeV purity JER of 21% turned out to be the largest
// purity entry.)
//
// The purity JER of systematics_CSVv2_2D_PF.C scales each flavor's template, bin by bin in ptRel, by R = T1 / T0, the
// template of the JER-smeared simulation over the nominal one, and refits. Smearing the jet pT changes the ptRel shape
// only through the jets it moves between pT windows (ptRel does not depend on the jet-pT scale), and the ptRel shape
// changes slowly with jet pT, so a resolution effect on R is SMOOTH in ptRel. Statistical structure of R (thin,
// pThat-weighted PYTHIA+HYDJET migrations) is not. Per system and window the fit is redone with
//   raw      R as used (the assigned value)
//   norm     R replaced by its template-weighted mean per flavor (only the flavor yields change)
//   linear   R replaced by a template-weighted straight line in ptRel per flavor (smooth shape change kept)
//   A        that line per flavor, rescaled so every flavor's yield changes by the b + c weighted mean of R
//   B        one R for all flavors summed, straight line in ptRel
// and the bin-to-bin scatter of R around that line (template-weighted RMS) is printed. Also for the JES (T2, T3).
// Windows: the four wide ones of the nominal binning (the 7-bin purity inherits them) and 200-300 GeV.
// Usage, from src/plots/bPurity/:  root -l -b -q studyPurityJER_PF.C
// Output: printout; rootFiles/bPurityResults/studyPurityJER_PF.root (R per flavor, raw and linear, per system and window)
#include "../../systematics/systematics_CSVv2_2D_PF.C"

void studyPurityJER_PF(){
  using namespace csys2d;
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
  if(!setCats("4cat")){ printf("ERROR: categories\n"); return; }
  setJetCollection("PFJets");
  fDpp = TFile::Open(Form("%s/%s", repo, dataName)); fMpp = TFile::Open(Form("%s/%s", repo, mcName)); fDAA = TFile::Open(Form("%s/%s", repo, dataPbPb));
  for(int s = 0; s < NS; s++) fS[s] = TFile::Open(Form("%s/%s", repo, Form(phFmt, sName[s])));
  nFit = (int)(5./ptRelBW + 0.5); gsMult = 1.175;
  const char *fl[3] = {"b", "c", "light"};
  TFile *fo = TFile::Open(Form("%s/rootFiles/bPurityResults/studyPurityJER_PF.root", repo), "RECREATE");
  printf("purity shift P_var / P_nom - 1 (%%): raw (assigned) | norm only | linear in ptRel;  scatter of R around the line (%%, b / c / light)\n");
  for(int sy = 0; sy < NSY; sy++){
    printf("\n==== %s\n", syName[sy]);
    for(int w = 0; w < NW; w++){
      Win X; if(!build(sy, w, X)) return;
      const double P0 = pcorr(X).p;
      for(int T : {1, 2, 3}){
        double R[NF][NPR], Rn[NF][NPR], Rl[NF][NPR], scat[NF];
        for(int F = 0; F < NF; F++){
          if(!tRatio(sy, w, F, T, R[F])) return;
          // template weights: the flavor's content per ptRel bin, summed over the CSVv2 categories
          double wt[NPR] = {}, sw = 0., sx = 0., sy2 = 0., sxx = 0., sxy = 0.;
          for(int j = 0; j < nFit; j++){ for(int k = 0; k < K; k++) wt[j] += X.t[F].n[k][j]; }
          for(int j = 0; j < nFit; j++){ const double x = (j + 0.5)*ptRelBW, v = wt[j]; sw += v; sx += v*x; sy2 += v*R[F][j]; sxx += v*x*x; sxy += v*x*R[F][j]; }
          const double mean = sw > 0. ? sy2/sw : 1., den = sw*sxx - sx*sx;
          const double b1 = fabs(den) > 0. ? (sw*sxy - sx*sy2)/den : 0., b0 = sw > 0. ? (sy2 - b1*sx)/sw : 1.;
          double s2 = 0.;
          for(int j = 0; j < NPR; j++){ const double x = (j + 0.5)*ptRelBW; Rn[F][j] = mean; Rl[F][j] = b0 + b1*x; }
          for(int j = 0; j < nFit; j++) s2 += wt[j]*pow(R[F][j] - Rl[F][j], 2);
          scat[F] = sw > 0. ? sqrt(s2/sw) : 0.;
          if(T == 1){ fo->cd();
            TH1D *hr = new TH1D(Form("h_R_%s_%s_w%d", fl[F], syTag[sy], w), Form("%s, %.0f-%.0f GeV: T1 / T0;p_{T}^{rel} (GeV);R", syName[sy], wLo[w], wHi[w]), nFit, 0., nFit*ptRelBW);
            TH1D *hl = (TH1D*) hr->Clone(Form("h_Rlin_%s_%s_w%d", fl[F], syTag[sy], w)), *hw = (TH1D*) hr->Clone(Form("h_wt_%s_%s_w%d", fl[F], syTag[sy], w));
            for(int j = 0; j < nFit; j++){ hr->SetBinContent(j + 1, R[F][j]); hl->SetBinContent(j + 1, Rl[F][j]); hw->SetBinContent(j + 1, wt[j]); }
            hr->Write(); hl->Write(); hw->Write(); }
        }
        auto refit = [&](double Rv[NF][NPR]){ Win Y = X;
          for(int F = 0; F < NF; F++){ Y.t[F].tot = 0.; for(int k = 0; k < K; k++) for(int j = 0; j < nPR; j++){ Y.t[F].n[k][j] *= Rv[F][j]; Y.t[F].e2[k][j] *= Rv[F][j]*Rv[F][j]; Y.t[F].tot += Y.t[F].n[k][j]; } }
          Y.truthB = Y.t[0].tot/(Y.t[0].tot + Y.t[1].tot + Y.t[2].tot); Y.rc = Y.t[1].tot/(Y.t[1].tot + Y.t[2].tot);
          return pcorr(Y).p/P0 - 1.; };
        const double dRaw = refit(R), dN = refit(Rn), dL = refit(Rl);
        // A: per-flavor linear shape, common yield factor = the b + c template-weighted mean of R
        // B: one common R for every flavor (all flavors summed), linear in ptRel
        double RA[NF][NPR], RB[NF][NPR];
        { double num = 0., den = 0., numB[NPR] = {}, denB[NPR] = {};
          for(int F = 0; F < NF; F++) for(int j = 0; j < nFit; j++){ double v = 0.; for(int k = 0; k < K; k++) v += X.t[F].n[k][j]; if(F < 2){ num += v*R[F][j]; den += v; } numB[j] += v*R[F][j]; denB[j] += v; }
          const double gY = den > 0. ? num/den : 1.;
          for(int F = 0; F < NF; F++){ double mF = 0., wF = 0.; for(int j = 0; j < nFit; j++){ double v = 0.; for(int k = 0; k < K; k++) v += X.t[F].n[k][j]; mF += v*Rl[F][j]; wF += v; }
            mF = wF > 0. ? mF/wF : 1.; for(int j = 0; j < NPR; j++) RA[F][j] = Rl[F][j]/mF*gY; }
          double sw = 0., sx = 0., sy2 = 0., sxx = 0., sxy = 0.;
          for(int j = 0; j < nFit; j++){ const double x = (j + 0.5)*ptRelBW, v = denB[j], r = denB[j] > 0. ? numB[j]/denB[j] : 1.; sw += v; sx += v*x; sy2 += v*r; sxx += v*x*x; sxy += v*x*r; }
          const double dd = sw*sxx - sx*sx, b1 = fabs(dd) > 0. ? (sw*sxy - sx*sy2)/dd : 0., b0 = sw > 0. ? (sy2 - b1*sx)/sw : 1.;
          for(int F = 0; F < NF; F++) for(int j = 0; j < NPR; j++) RB[F][j] = b0 + b1*(j + 0.5)*ptRelBW; }
        const double dA = refit(RA), dB = refit(RB);
        printf("  %3.0f-%3.0f  %s   raw %+6.1f   norm %+6.1f   linear %+6.1f   A %+6.1f   B %+6.1f   scatter %5.1f %5.1f %5.1f\n", wLo[w], wHi[w], T == 1 ? "JER " : (T == 2 ? "JES+" : "JES-"),
               100*dRaw, 100*dN, 100*dL, 100*dA, 100*dB, 100*scat[0], 100*scat[1], 100*scat[2]);
      }
    }
  }
  fo->Close();
}
