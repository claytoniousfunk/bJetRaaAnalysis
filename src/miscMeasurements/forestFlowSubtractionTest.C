// Test of the hypothesis that the PbPb low-pT jet excess is driven by an
// AZIMUTHALLY UNMODULATED background subtraction.
//
// ---------------------------------------------------------------------------
// WHAT THE CMSSW RELEASE SAYS (CMSSW_10_3_4, read 2026-08-22)
//
// The MinBias forest was produced by runForestAOD_pponAA_DATA_103X.py, whose
// jet sequence is fullJetSequence_pponAA_data_cff.py. In that sequence:
//
//     akCs4PFJets +          <-- this one runs
//     #aakFlowPuCs4PFJets    <-- commented out
//
// and akCs4PFJets is configured in HeavyIonsAnalysis/JetAnalysis/python/
// akPFJets_cfi.py with
//
//     useModulatedRho = cms.bool(False)
//
// The flow-corrected clone akFlowPuCs4PFJets (useModulatedRho = True) exists
// in the same file and is never scheduled. hiFJRhoFlowModulationProducer IS
// run in rerecoRho_cff.py, so the event-by-event v2/v3 fit is computed and
// then thrown away: CSJetProducer.cc only reads rhoFlowFitParams when
// useModulatedRho_ is true.
//
// Consequence: every ghost in the constituent subtraction is assigned
// pt = rho(eta) * area with NO cos(2(phi - Psi2)) term. The true underlying
// event has one. So each jet carries a residual
//
//     delta(phi) = 2 * v2 * rho * A * cos(2(phi - Psi2))
//
// that was never subtracted -- positive in-plane (under-subtracted), negative
// out-of-plane (over-subtracted).
//
// ---------------------------------------------------------------------------
// WHY THAT BREAKS THE CENTRALITY HIERARCHY
//
// On a steeply falling spectrum N ~ exp(-pT/T), a residual delta shifts the
// yield at fixed measured pT by exp(delta/T). Averaged over phi that is a
// Bessel factor I0(2 v2 rho A / T) > 1: a net INFLATION of the low-pT yield.
//
// The size of that inflation is controlled by the product v2 * rho. rho falls
// monotonically toward peripheral while v2 RISES, so the product is
// non-monotonic and peaks in mid-central. The inflation is therefore LARGER
// in 10-30% than in 0-10% -- which lifts the mid-central classes up toward
// the central ones and destroys the expected ordering.
//
// ---------------------------------------------------------------------------
// WHAT THIS MACRO MEASURES, AND WHAT IT CANNOT
//
// It builds an event plane from forward charged PF candidates (2.2<|eta|<3.0)
// in two sub-events, and measures <cos2(phi - Psi_A)> for
//   (a) midrapidity charged hadrons  -- a CONTROL with a known answer
//   (b) jets, in pT bins             -- the quantity of interest
// both against Psi_A alone, corrected by R_A = sqrt(<cos2(Psi_A - Psi_B)>).
//
// *** THE ABSOLUTE SCALE IS NOT TRUSTWORTHY. *** The hadron control comes out
// around 0.12-0.20 where the literature says 0.05-0.12, i.e. this simple
// resolution correction over-corrects by roughly a factor 2 (no recentering,
// no flattening, and eta-decorrelation between sub-events separated by 4.4
// units of eta biases R_A low). Do NOT quote the absolute jet v2 from here.
//
// What IS robust is the RATIO jet/hadron, because both are divided by the
// same R_A and it cancels exactly. That ratio is the result: low-pT jets are
// modulated several times more strongly than the bulk, by an amount that
// tracks v2 * rho * A / T.
//
// An earlier version of this test put the event plane at |eta|<1.0, the same
// region the jets live in. That gave <cos2dphi> ~ 0.8, which is pure
// autocorrelation -- a 30-45 GeV "jet" in central PbPb is largely built from
// the very 0.3-3 GeV candidates that define the plane. The forward gap here
// (>= 0.8 units between the jet cone edge and the nearest EP candidate) is
// what makes the number mean anything.
//
// Usage: root -l -b -q 'forestFlowSubtractionTest.C("/path/to/HiForestAOD.root")'
// Run from: src/miscMeasurements/

const char *defaultForest = "/home/clayton/Downloads/HiForestAOD_100.root";

// Event plane sub-events: forward, symmetric, well clear of the jet cones.
const double epEtaMin = 2.2, epEtaMax = 3.0;
const int    epMinCand = 30;      // per sub-event, else the plane is noise

// Jets restricted well inside the flat central rho bin (|eta|<1.3) AND far
// enough from the EP region that the cone cannot touch it.
const double jetEtaMax = 1.0;

const int    NC = 4;
const int    centEdge[NC+1] = {0, 20, 60, 100, 160};
const char  *centLabel[NC]  = {"0-10%", "10-30%", "30-50%", "50-80%"};

const int    NP = 3;
const double ptLo[NP] = {30, 45, 70}, ptHi[NP] = {45, 70, 400};

void forestFlowSubtractionTest(const char *forestFile = defaultForest)
{
  TFile *f = TFile::Open(forestFile);
  if(!f || f->IsZombie()){ printf("ERROR: cannot open %s\n", forestFile); return; }
  printf("Forest: %s\n\n", forestFile);

  TTree *ev = (TTree*) f->Get("hiEvtAnalyzer/HiTree");
  TTree *jt = (TTree*) f->Get("akCs4PFJetAnalyzer/t");
  TTree *pf = (TTree*) f->Get("pfcandAnalyzer/pfTree");
  if(!ev || !jt || !pf){ printf("ERROR: a required tree is missing\n"); return; }
  if(ev->GetEntries() != jt->GetEntries() || ev->GetEntries() != pf->GetEntries()){
    printf("ERROR: tree entry counts differ (%lld/%lld/%lld); they are read by\n"
           "       index here and must correspond event for event.\n",
           ev->GetEntries(), jt->GetEntries(), pf->GetEntries());
    return;
  }

  int hiBin; ev->SetBranchAddress("hiBin", &hiBin);
  int nref; float jtpt[10000], jteta[10000], jtphi[10000];
  jt->SetBranchAddress("nref", &nref);
  jt->SetBranchAddress("jtpt", jtpt);
  jt->SetBranchAddress("jteta", jteta);
  jt->SetBranchAddress("jtphi", jtphi);
  std::vector<float> *pfPt = nullptr, *pfEta = nullptr, *pfPhi = nullptr;
  std::vector<int> *pfId = nullptr;
  pf->SetBranchAddress("pfPt",  &pfPt);
  pf->SetBranchAddress("pfEta", &pfEta);
  pf->SetBranchAddress("pfPhi", &pfPhi);
  pf->SetBranchAddress("pfId",  &pfId);

  double sumJet[NC][NP] = {{0}}, sumHad[NC] = {0}, sumRes[NC] = {0};
  long   nJet[NC][NP]   = {{0}}, nHad[NC]   = {0}, nRes[NC]   = {0};

  for(Long64_t i = 0; i < ev->GetEntries(); i++){
    ev->GetEntry(i); jt->GetEntry(i); pf->GetEntry(i);
    int c = -1;
    for(int k = 0; k < NC; k++) if(hiBin >= centEdge[k] && hiBin < centEdge[k+1]) c = k;
    if(c < 0) continue;

    // two forward sub-event planes, charged hadrons only, as the CMSSW flow
    // producer selects them (particleId == 1, 0.3 < pT < 3)
    double CA = 0, SA = 0, CB = 0, SB = 0; int nA = 0, nB = 0;
    for(size_t p = 0; p < pfPt->size(); p++){
      if(pfId->at(p) != 1) continue;
      double pt = pfPt->at(p);
      if(pt < 0.3 || pt > 3.0) continue;
      double e = pfEta->at(p);
      if(e > epEtaMin && e < epEtaMax){ CA += cos(2*pfPhi->at(p)); SA += sin(2*pfPhi->at(p)); nA++; }
      else if(e < -epEtaMin && e > -epEtaMax){ CB += cos(2*pfPhi->at(p)); SB += sin(2*pfPhi->at(p)); nB++; }
    }
    if(nA < epMinCand || nB < epMinCand) continue;
    double PsiA = atan2(SA, CA)/2., PsiB = atan2(SB, CB)/2.;
    sumRes[c] += cos(2*(PsiA - PsiB)); nRes[c]++;

    // control: midrapidity charged hadrons, known answer
    for(size_t p = 0; p < pfPt->size(); p++){
      if(pfId->at(p) != 1) continue;
      double pt = pfPt->at(p);
      if(pt < 1.0 || pt > 3.0) continue;
      if(fabs(pfEta->at(p)) > 1.0) continue;
      sumHad[c] += cos(2*(pfPhi->at(p) - PsiA)); nHad[c]++;
    }

    for(int j = 0; j < nref; j++){
      if(fabs(jteta[j]) > jetEtaMax) continue;
      double d = jtphi[j] - PsiA;
      for(int k = 0; k < NP; k++)
        if(jtpt[j] >= ptLo[k] && jtpt[j] < ptHi[k]){ sumJet[c][k] += cos(2*d); nJet[c][k]++; }
    }
  }

  printf("Modulation w.r.t. a forward event plane (%.1f<|eta|<%.1f), jets |eta|<%.1f\n",
         epEtaMin, epEtaMax, jetEtaMax);
  printf("Absolute scale is NOT calibrated -- read the RATIO column.\n\n");
  printf("%-9s %7s %10s", "cent", "R_A", "hadron");
  for(int k = 0; k < NP; k++) printf(" %16s", Form("%.0f-%.0f GeV", ptLo[k], ptHi[k]));
  printf(" %10s\n", "ratio");
  for(int c = 0; c < NC; c++){
    if(!nRes[c] || !nHad[c]) continue;
    double R = sqrt(std::max(sumRes[c]/nRes[c], 1e-9));
    double had = (sumHad[c]/nHad[c])/R;
    printf("%-9s %7.3f %10.4f", centLabel[c], R, had);
    double first = -1.;
    for(int k = 0; k < NP; k++){
      if(nJet[c][k] < 25){ printf(" %16s", Form("(n=%ld)", nJet[c][k])); continue; }
      double m = (sumJet[c][k]/nJet[c][k])/R, e = sqrt(0.5/nJet[c][k])/R;
      if(first < 0) first = m;
      printf(" %7.3f+-%-7.3f", m, e);
    }
    if(first > 0 && had > 0) printf(" %9.1fx\n", first/had); else printf(" %10s\n", "-");
  }

  // Expected inflation of the phi-integrated low-pT yield, for reference.
  // rho values are the per-class averages measured from hiFJRhoAnalyzer in
  // this same forest; v2 is bulk charged-hadron v2 from the literature; T is
  // the local inverse slope fit to our own PbPb jet spectrum over 30-60 GeV.
  printf("\nExpected phi-integrated yield inflation I0(2 v2 rho A / T):\n");
  printf("%-9s %8s %7s %6s %12s %12s\n", "cent", "rho", "v2", "T", "inflation", "rel 0-10%");
  double rho[NC] = {166.2, 91.5, 35.2, 8.0};
  double v2b[NC] = {0.050, 0.090, 0.115, 0.090};
  double Tsl[NC] = {10.73, 8.82, 6.27, 7.59};
  double A = TMath::Pi()*0.4*0.4;
  double ref = 0.;
  for(int c = 0; c < NC; c++){
    double x = 2*v2b[c]*rho[c]*A/Tsl[c];
    double I0 = TMath::BesselI0(x);
    if(c == 0) ref = I0 - 1.;
    printf("%-9s %8.1f %7.3f %6.2f %12.3f %12.2f\n",
           centLabel[c], rho[c], v2b[c], Tsl[c], I0, ref > 0 ? (I0-1.)/ref : 0.);
  }
  printf("\nNon-monotonic: rho falls toward peripheral while v2 rises, so the\n"
         "inflation peaks mid-central and 10-30%% is lifted relative to 0-10%%.\n");
}
