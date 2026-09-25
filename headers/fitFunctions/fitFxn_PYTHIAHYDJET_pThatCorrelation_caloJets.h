
TF1 *fitFxn_PYTHIAHYDJET_pThatCorrelation_caloJets_C4;
TF1 *fitFxn_PYTHIAHYDJET_pThatCorrelation_caloJets_C3;
TF1 *fitFxn_PYTHIAHYDJET_pThatCorrelation_caloJets_C2;
TF1 *fitFxn_PYTHIAHYDJET_pThatCorrelation_caloJets_C1;

void loadFitFxn_PYTHIAHYDJET_pThatCorrelation_caloJets(){

  cout << "~~ loading pThatCorrelation fit function ..." << endl;

  fitFxn_PYTHIAHYDJET_pThatCorrelation_caloJets_C4 = new TF1("fitFxn_PYTHIAHYDJET_pThatCorrelation_caloJets_C4","[0] + [1]*exp(-[2]*x) + [3]*exp(-[4]*x) + [5]*exp(-[6]*x)");
  fitFxn_PYTHIAHYDJET_pThatCorrelation_caloJets_C3 = new TF1("fitFxn_PYTHIAHYDJET_pThatCorrelation_caloJets_C3","[0] + [1]*exp(-[2]*x) + [3]*exp(-[4]*x) + [5]*exp(-[6]*x)");
  fitFxn_PYTHIAHYDJET_pThatCorrelation_caloJets_C2 = new TF1("fitFxn_PYTHIAHYDJET_pThatCorrelation_caloJets_C2","[0] + [1]*exp(-[2]*x) + [3]*exp(-[4]*x) + [5]*exp(-[6]*x)");
  fitFxn_PYTHIAHYDJET_pThatCorrelation_caloJets_C1 = new TF1("fitFxn_PYTHIAHYDJET_pThatCorrelation_caloJets_C1","[0] + [1]*exp(-[2]*x) + [3]*exp(-[4]*x) + [5]*exp(-[6]*x)");

  fitFxn_PYTHIAHYDJET_pThatCorrelation_caloJets_C4->SetParameter(0,-0.740806);
  fitFxn_PYTHIAHYDJET_pThatCorrelation_caloJets_C4->SetParameter(1,4.52383);
  fitFxn_PYTHIAHYDJET_pThatCorrelation_caloJets_C4->SetParameter(2,0.0694607);
  fitFxn_PYTHIAHYDJET_pThatCorrelation_caloJets_C4->SetParameter(3,16.7354);
  fitFxn_PYTHIAHYDJET_pThatCorrelation_caloJets_C4->SetParameter(4,0.275494);
  fitFxn_PYTHIAHYDJET_pThatCorrelation_caloJets_C4->SetParameter(5,3.02867);
  fitFxn_PYTHIAHYDJET_pThatCorrelation_caloJets_C4->SetParameter(6,0.00040674);

  fitFxn_PYTHIAHYDJET_pThatCorrelation_caloJets_C3->SetParameter(0,3.36004);
  fitFxn_PYTHIAHYDJET_pThatCorrelation_caloJets_C3->SetParameter(1,3.21552);
  fitFxn_PYTHIAHYDJET_pThatCorrelation_caloJets_C3->SetParameter(2,0.0719387);
  fitFxn_PYTHIAHYDJET_pThatCorrelation_caloJets_C3->SetParameter(3,2.29909);
  fitFxn_PYTHIAHYDJET_pThatCorrelation_caloJets_C3->SetParameter(4,0.0450385);
  fitFxn_PYTHIAHYDJET_pThatCorrelation_caloJets_C3->SetParameter(5,-1.16322);
  fitFxn_PYTHIAHYDJET_pThatCorrelation_caloJets_C3->SetParameter(6,-0.000618867);

  fitFxn_PYTHIAHYDJET_pThatCorrelation_caloJets_C2->SetParameter(0,1.59362);
  fitFxn_PYTHIAHYDJET_pThatCorrelation_caloJets_C2->SetParameter(1,0.794807);
  fitFxn_PYTHIAHYDJET_pThatCorrelation_caloJets_C2->SetParameter(2,0.00279327);
  fitFxn_PYTHIAHYDJET_pThatCorrelation_caloJets_C2->SetParameter(3,0.543028);
  fitFxn_PYTHIAHYDJET_pThatCorrelation_caloJets_C2->SetParameter(4,0.0586818);
  fitFxn_PYTHIAHYDJET_pThatCorrelation_caloJets_C2->SetParameter(5,5.58825);
  fitFxn_PYTHIAHYDJET_pThatCorrelation_caloJets_C2->SetParameter(6,0.0545184);

  fitFxn_PYTHIAHYDJET_pThatCorrelation_caloJets_C1->SetParameter(0,1.85442);
  fitFxn_PYTHIAHYDJET_pThatCorrelation_caloJets_C1->SetParameter(1,-50.9639);
  fitFxn_PYTHIAHYDJET_pThatCorrelation_caloJets_C1->SetParameter(2,0.0220517);
  fitFxn_PYTHIAHYDJET_pThatCorrelation_caloJets_C1->SetParameter(3,50.9561);
  fitFxn_PYTHIAHYDJET_pThatCorrelation_caloJets_C1->SetParameter(4,0.0246736);
  fitFxn_PYTHIAHYDJET_pThatCorrelation_caloJets_C1->SetParameter(5,5.84145);
  fitFxn_PYTHIAHYDJET_pThatCorrelation_caloJets_C1->SetParameter(6,0.0126487);



}
