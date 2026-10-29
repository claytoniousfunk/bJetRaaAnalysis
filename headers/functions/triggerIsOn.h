// A path is on when it fired and its prescale is non-zero. Both must be > 0:
// hltanalysis writes -1 to the decision and the prescale of a path that is not
// in the menu (e.g. HLT_HIL3Mu5_v1 in the 2018 PbPb and PYTHIA+HYDJET forests,
// -1 in every event), which a != 0 test would count as fired.
bool triggerIsOn(int triggerValue, int prescaleValue){

	bool result = true;
	if(triggerValue <= 0 || prescaleValue <= 0) result = false;
	return result;

}
