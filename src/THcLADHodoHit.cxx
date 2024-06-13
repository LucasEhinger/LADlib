#include "THcLADHodoHit.h"

#include <iostream>

//___________________________________________________________________
THcLADHodoHit::THcLADHodoHit(Int_t postdc, Int_t negtdc, Double_t posadc, Double_t negadc, Int_t ipad,
                             THcLADHodoPlane *hp)
    : fPosTDC(postdc), fNegTDC(negtdc), fPosADC_Ped(posadc), fNegADC_Ped(negadc), fPaddleNumber(ipad), fPlane(hp) {}

//___________________________________________________________________
THcLADHodoHit::~THcLADHodoHit() {}

void THcLADHodoHit::SetCorrectedTimes(Double_t pos, Double_t neg, Double_t postof, Double_t negtof, Double_t timeave) {
  fPosCorrectedTime = pos;
  fNegCorrectedTime = neg;
  fPosTOFCorrectedTime = postof;
  fNegTOFCorrectedTime = negtof;
  fScinCorrectedTime = timeave;
  fHasCorrectedTimes = kTRUE;
}

void THcLADHodoHit::SetCorrectedTimes(Double_t pos, Double_t neg) {
  fPosCorrectedTime = pos;
  fNegCorrectedTime = neg;
  fHasCorrectedTimes = kFALSE;
}

ClassImp(THcLADHodoHit)
