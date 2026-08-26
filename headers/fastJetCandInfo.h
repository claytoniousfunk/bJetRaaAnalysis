#pragma once
#include "fastjet/PseudoJet.hh"


class CandInfo : public fastjet::PseudoJet::UserInfoBase{

 public:

 CandInfo(int idx, int id) : _idx(idx), _id(id) {}

  int getIndex() const {return _idx;}
  
  int getId() const {return _id;}

  bool isCharged()   const {return _id == 1 || _id == 2 || _id == 3; }

 private:
  
  int _idx, _id;

};
