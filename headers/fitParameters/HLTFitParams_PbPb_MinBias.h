// HLT_HIL3Mu12 single-muon efficiency fit, f(pT) = p0 - (p0 - p1) exp(-p2 pT),
// 0 <= p1 <= p0 enforced, fit for 15 < pT < 100 GeV, tight ID, |eta| < 2.0.
// PbPb HIMinimumBias0 Parts 1-4 (2026-09-19..21 scans), h_muTrigEff histograms
// (prescale-1 events, other-muon veto), 5% slices merged into the nominal
// classes by src/plots/triggerEfficiency/mergeMuTrigEff_coarseCent.C and fit by
// src/plots/triggerEfficiency/plotMuonTriggerEfficiency.C.
// C1 = 0-10%, C2 = 10-30%, C3 = 30-50%, C4 = 50-80%.

double HLTFitParam_C4_0 = 0.974438;
double HLTFitParam_C4_1 = 9.71528e-08;
double HLTFitParam_C4_2 = 0.287125;

double e_HLTFitParam_C4_0 = 0.0118152;
double e_HLTFitParam_C4_1 = 0.488776;
double e_HLTFitParam_C4_2 = 0.096725;

////////////////////////////////////////

double HLTFitParam_C3_0 = 0.982905;
double HLTFitParam_C3_1 = 0.907164;
double HLTFitParam_C3_2 = 0.047198;

double e_HLTFitParam_C3_0 = 0.0409154;
double e_HLTFitParam_C3_1 = 0.0717132;
double e_HLTFitParam_C3_2 = 0.0880277;

////////////////////////////////////////

double HLTFitParam_C2_0 = 0.934924;
double HLTFitParam_C2_1 = 1.83983e-05;
double HLTFitParam_C2_2 = 0.191132;

double e_HLTFitParam_C2_0 = 0.00668975;
double e_HLTFitParam_C2_1 = 0.746734;
double e_HLTFitParam_C2_2 = 0.0142427;

////////////////////////////////////////

double HLTFitParam_C1_0 = 0.910733;
double HLTFitParam_C1_1 = 0.499749;
double HLTFitParam_C1_2 = 0.0880751;

double e_HLTFitParam_C1_0 = 0.0237961;
double e_HLTFitParam_C1_1 = 0.247626;
double e_HLTFitParam_C1_2 = 0.0459619;

////////////////////////////////////////

