
TF1 *fitFxn_PYTHIA_pThatCorrelation_caloJets;

void loadFitFxn_PYTHIA_pThatCorrelation_caloJets(){

  cout << "~~ loading pThatCorrelation fit function ..." << endl;

  fitFxn_PYTHIA_pThatCorrelation_caloJets = new TF1("fitFxn_PYTHIA_pThatCorrelation_caloJets","[0] + [1]*exp(-[2]*x) + [3]*exp(-[4]*x) + [5]*exp(-[6]*x)");


  fitFxn_PYTHIA_pThatCorrelation->SetParameter(0,1.77678);
  fitFxn_PYTHIA_pThatCorrelation->SetParameter(1,-11.894);
  fitFxn_PYTHIA_pThatCorrelation->SetParameter(2,0.342266);
  fitFxn_PYTHIA_pThatCorrelation->SetParameter(3,8.18003);
  fitFxn_PYTHIA_pThatCorrelation->SetParameter(4,0.0650943);
  fitFxn_PYTHIA_pThatCorrelation->SetParameter(5,0.660459);
  fitFxn_PYTHIA_pThatCorrelation->SetParameter(6,0.00462879);


}
