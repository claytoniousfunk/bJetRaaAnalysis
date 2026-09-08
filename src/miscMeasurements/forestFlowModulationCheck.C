// Did the constituent-subtraction flow modulation actually fire?
//
// Compares two jet collections in the SAME forest jet-by-jet and decomposes
// the per-jet pT difference into a phi-independent offset and a cos(2 dphi)
// modulation. Use it to confirm that useModulatedRho is doing something before
// committing to a production scan.
//
// ---------------------------------------------------------------------------
// WHY THIS EXISTS
//
// In CMSSW_10_3_4, CSJetProducer.cc:171-172 reads
//
//     desc.add<double> ("minFlowChi2Prob", false);
//     desc.add<double> ("maxFlowChi2Prob", false);
//
// -- double parameters given a bool default, so both default to 0.0. The gate
// at line 98-99 is then
//
//     bool minProb = val > minFlowChi2Prob_;   // val > 0.0  -> true
//     bool maxProb = val < maxFlowChi2Prob_;   // val < 0.0  -> ALWAYS FALSE
//
// and since val is a p-value in (0,1), the modulation never applies. Neither
// akCs4PFJets nor akFlowPuCs4PFJets in HeavyIonsAnalysis/JetAnalysis/python/
// akPFJets_cfi.py sets these, so setting useModulatedRho = True by itself does
// NOTHING. Only RecoJets/JetProducers/python/akCs4PFJets_cfi.py sets them
// (0.05/0.95), and the forest does not use that file.
//
// Verified empirically on a 1000-event MinBias forest: with the stock config
// the covariance below was 0.09 +- 0.23; after setting minFlowChi2Prob and
// maxFlowChi2Prob it was -3.79, with the correct sign and pi-periodicity.
//
// ---------------------------------------------------------------------------
// *** THE ESTIMATOR MATTERS -- DO NOT USE <D cos2dphi> ***
//
// The naive <D cos(2 dphi)> is BIASED and will show a confident 4-5 sigma
// "signal" even when the modulation is entirely off. Two things multiply:
//
//   - a large phi-independent offset <D> (e.g. -1.9 GeV if the two collections
//     also differ in their rho source), and
//   - the jet sample's own azimuthal anisotropy <cos 2 dphi>, which is large
//     (0.35-0.75) mostly from autocorrelation with the event plane.
//
// Their product reproduced the entire apparent signal in the first test. The
// unbiased estimator is the COVARIANCE
//
//     Cov = <D cos2dphi> - <D><cos2dphi>
//
// which is what this macro reports. The tell for the biased case is the dphi
// profile: a genuine flow correction subtracts MORE in-plane, so D must be
// most negative near dphi = 0 and pi. If the profile is flat, or most negative
// out-of-plane, there is no modulation and you are looking at the bias.
//
// EXPECTED Cov if the modulation is live: about -v2 * rho * A, i.e. roughly
// -1 to -4 GeV. Zero within errors means it is not firing.
//
// The event plane here is rebuilt exactly as HiFJRhoFlowModulationProducer
// builds it (charged PF candidates, |eta|<1, 0.3<pT<3, particleId==1) so it
// matches the plane the producer actually used. That still leaves some
// dilution because the producer uses its FITTED plane, not the raw Q-vector.
// To remove that entirely, enable hiPuRhoR3Analyzer in the forest config -- it
// is the only analyzer configured with useModulatedRho = True and so is the
// only one that writes rhoFlowFitParams to the tree.
//
// Usage:
//   root -l -b -q 'forestFlowModulationCheck.C("forest.root")'
//   root -l -b -q 'forestFlowModulationCheck.C("forest.root","akCs4PFJetAnalyzer","akFlowPuCs4PFJetAnalyzer")'
// Run from: src/miscMeasurements/

// Jets are taken from the FIRST collection and matched into the second, so the
// first should be the nominal one.
const char *defNominal  = "akCs4PFJetAnalyzer";
const char *defModulated = "akFlowPuCs4PFJetAnalyzer";

const double jetEtaMax = 1.3;      // stay inside the flat central rho bin
const double jetPtLo = 20., jetPtHi = 40.;   // narrow, so composition cannot bias Cov
const double matchDR = 0.2;
const int    epMinCand = 100;      // the producer's own nFill requirement

const int    NC = 4;
const int    centEdge[NC+1] = {0, 20, 60, 100, 140};
const char  *centLabel[NC]  = {"0-10%", "10-30%", "30-50%", "50-70%"};

const int NB = 6;                  // dphi profile bins over [0, pi)

void forestFlowModulationCheck(const char *forestFile,
                               const char *nomTree = defNominal,
                               const char *modTree = defModulated)
{
  TFile *f = TFile::Open(forestFile);
  if(!f || f->IsZombie()){ printf("ERROR: cannot open %s\n", forestFile); return; }

  TTree *ev = (TTree*) f->Get("hiEvtAnalyzer/HiTree");
  TTree *a  = (TTree*) f->Get(Form("%s/t", nomTree));
  TTree *b  = (TTree*) f->Get(Form("%s/t", modTree));
  TTree *pf = (TTree*) f->Get("pfcandAnalyzer/pfTree");
  TTree *rh = (TTree*) f->Get("hiFJRhoAnalyzer/t");
  if(!ev || !a || !b || !pf || !rh){
    printf("ERROR: missing tree (ev=%p nom=%p mod=%p pf=%p rho=%p)\n",
           (void*)ev,(void*)a,(void*)b,(void*)pf,(void*)rh);
    return;
  }
  if(a->GetEntries() != ev->GetEntries() || b->GetEntries() != ev->GetEntries()){
    printf("ERROR: entry counts differ; trees are read by index and must correspond.\n");
    return;
  }
  printf("Forest   : %s\n", forestFile);
  printf("nominal  : %s\nmodulated: %s\n\n", nomTree, modTree);

  int hiBin; ev->SetBranchAddress("hiBin", &hiBin);
  int na, nb;
  float pa[10000], ea[10000], fa[10000], pb[10000], eb[10000], fb[10000];
  a->SetBranchAddress("nref",&na); a->SetBranchAddress("jtpt",pa);
  a->SetBranchAddress("jteta",ea); a->SetBranchAddress("jtphi",fa);
  b->SetBranchAddress("nref",&nb); b->SetBranchAddress("jtpt",pb);
  b->SetBranchAddress("jteta",eb); b->SetBranchAddress("jtphi",fb);
  std::vector<float> *pt=nullptr,*et=nullptr,*ph=nullptr; std::vector<int> *id=nullptr;
  pf->SetBranchAddress("pfPt",&pt); pf->SetBranchAddress("pfEta",&et);
  pf->SetBranchAddress("pfPhi",&ph); pf->SetBranchAddress("pfId",&id);
  std::vector<double> *R=nullptr; rh->SetBranchAddress("rho",&R);

  const double A = TMath::Pi()*0.4*0.4;
  double sD[NC]={0}, sC[NC]={0}, sDC[NC]={0}, sD2[NC]={0}, sC2[NC]={0}, srho[NC]={0};
  long   n[NC]={0}, nEv[NC]={0}, nFired[NC]={0}, nFillOK[NC]={0}, nTry[NC]={0};
  double pD[NC][NB]={{0}}; long pN[NC][NB]={{0}};

  for(Long64_t i = 0; i < ev->GetEntries(); i++){
    ev->GetEntry(i); a->GetEntry(i); b->GetEntry(i); pf->GetEntry(i); rh->GetEntry(i);
    int c = -1;
    for(int k = 0; k < NC; k++) if(hiBin >= centEdge[k] && hiBin < centEdge[k+1]) c = k;
    if(c < 0) continue;
    nEv[c]++;

    bool same = (na == nb);
    if(same) for(int j = 0; j < na; j++) if(fabs(pa[j]-pb[j]) > 1e-6){ same = false; break; }
    if(!same) nFired[c]++;

    // the producer's own event-plane definition
    double C=0, S=0; int nFill=0;
    for(size_t p = 0; p < pt->size(); p++){
      if(id->at(p) != 1) continue;
      double q = pt->at(p);
      if(q < 0.3 || q > 3.0) continue;
      if(fabs(et->at(p)) > 1.0) continue;
      C += cos(2*ph->at(p)); S += sin(2*ph->at(p)); nFill++;
    }
    if(nFill >= epMinCand) nFillOK[c]++;
    if(nFill < epMinCand) continue;
    double Psi = atan2(S, C)/2.;

    for(int j = 0; j < na; j++){
      if(fabs(ea[j]) > jetEtaMax || pa[j] < jetPtLo || pa[j] >= jetPtHi) continue;
      nTry[c]++;
      int best = -1; double bdr = matchDR;
      for(int k = 0; k < nb; k++){
        double de = ea[j]-eb[k], dp = TVector2::Phi_mpi_pi(fa[j]-fb[k]);
        double dr = sqrt(de*de + dp*dp);
        if(dr < bdr){ bdr = dr; best = k; }
      }
      if(best < 0) continue;
      double D = pb[best] - pa[j], cc = cos(2*(fa[j]-Psi));
      sD[c]+=D; sC[c]+=cc; sDC[c]+=D*cc; sD2[c]+=D*D; sC2[c]+=cc*cc; n[c]++;
      if(R->size() > 3) srho[c] += R->at(3);
      double d = TVector2::Phi_mpi_pi(fa[j]-Psi); if(d < 0) d += TMath::Pi();
      int ib = (int)(d/(TMath::Pi()/NB)); if(ib >= NB) ib = NB-1;
      pD[c][ib] += D; pN[c][ib]++;
    }
  }

  printf("Jets %.0f-%.0f GeV, |eta|<%.1f, matched dR<%.2f\n\n", jetPtLo, jetPtHi, jetEtaMax, matchDR);
  printf("%-9s %6s %8s %8s %7s %7s %8s %9s %10s %8s %11s\n",
         "cent","Nevt","EPok%","fired%","njet","match%","<D>","RMS(D)","Cov","signif","v2 applied");
  for(int c = 0; c < NC; c++){
    if(!nEv[c]) continue;
    if(n[c] < 15){
      printf("%-9s %6ld %7.0f%% %7.0f%% %7ld  (too few matched jets)\n",
             centLabel[c], nEv[c], 100.*nFillOK[c]/nEv[c], 100.*nFired[c]/nEv[c], n[c]);
      continue;
    }
    double N=n[c], mD=sD[c]/N, mC=sC[c]/N, mDC=sDC[c]/N;
    double vD=sD2[c]/N-mD*mD, vC=sC2[c]/N-mC*mC;
    double cov = mDC - mD*mC;                 // unbiased; see header
    double ecov = sqrt(vD*vC/N);
    double rhoA = (srho[c]/N)*A;
    printf("%-9s %6ld %7.0f%% %7.0f%% %7ld %6.0f%% %8.3f %9.3f %10.3f %8.1f %11.4f\n",
           centLabel[c], nEv[c], 100.*nFillOK[c]/nEv[c], 100.*nFired[c]/nEv[c],
           n[c], 100.*n[c]/std::max(1L,nTry[c]), mD, sqrt(vD), cov,
           ecov>0?fabs(cov)/ecov:0, rhoA>0?-cov/rhoA:0);
  }

  printf("\n--- <D> vs (phi_jet - Psi2). In-plane is first and last row. ---\n");
  printf("A live modulation is MOST NEGATIVE in-plane. Flat, or most negative\n");
  printf("in the middle rows, means it is not firing.\n\n");
  printf("%-9s", "dphi");
  for(int c = 0; c < NC; c++) printf(" %14s", centLabel[c]);
  printf("\n");
  for(int ib = 0; ib < NB; ib++){
    printf("%3.0f-%-5.0f", 180.*ib/NB, 180.*(ib+1)/NB);
    for(int c = 0; c < NC; c++){
      if(pN[c][ib] < 5){ printf(" %14s", "-"); continue; }
      printf(" %9.3f(%3ld)", pD[c][ib]/pN[c][ib], pN[c][ib]);
    }
    printf("\n");
  }

  printf("\nIf 'fired%%' is well below 100 with EPok%% at 100, the chi2 gate is\n"
         "rejecting the fit -- set minFlowChi2Prob=0.0 and maxFlowChi2Prob=1.0.\n");
}
