#ifndef THcLADHodoHit_h
#define THcLADHodoHit_h

#include "THcLADHodoPlane.h"
#include "THcRawHodoHit.h"
#include "TObject.h"

class THcLADHodoHit : public TObject {
public:
  THcLADHodoHit(Int_t postdc, Int_t negtdc, Double_t posadc, Double_t negadc, Int_t ipad, THcLADHodoPlane *hp);
  virtual ~THcLADHodoHit();

  Double_t GetPosADC() const { return fPosADC_Ped; }
  Double_t GetNegADC() const { return fNegADC_Ped; }
  Double_t GetPosADCpeak() const { return fPosADC_Peak; }
  Double_t GetNegADCpeak() const { return fNegADC_Peak; }
  Double_t GetPosADCtime() const { return fPosADC_Time; }
  Double_t GetNegADCtime() const { return fNegADC_Time; }
  Int_t GetPosTDC() const { return fPosTDC; }
  Int_t GetNegTDC() const { return fNegTDC; }
  Int_t GetPaddleNumber() const { return fPaddleNumber; }
  Int_t GetPaddleCenter() const { return fPaddleCenter; }

  void SetPaddleCenter(Double_t padcenter) { fPaddleCenter = padcenter; }
  void SetPosADCpeak(Double_t adc) { fPosADC_Peak = adc; }
  void SetNegADCpeak(Double_t adc) { fNegADC_Peak = adc; }
  void SetPosADCtime(Double_t ptime) { fPosADC_Time = ptime; }
  void SetNegADCtime(Double_t ptime) { fNegADC_Time = ptime; }
  void SetCorrectedTimes(Double_t pos, Double_t neg);
  void SetCorrectedTimes(Double_t pos, Double_t neg, Double_t postof, Double_t negtof, Double_t timeave);

protected:
  Int_t fPosTDC;
  Int_t fNegTDC;
  Double_t fPosADC_Ped; // Pedestal subtracted ADC
  Double_t fNegADC_Ped;
  Double_t fPosADC_Peak; // ADC peak amplitude
  Double_t fNegADC_Peak; // ADC peak amplitude
  Double_t fPosADC_Time;
  Double_t fNegADC_Time;

  Double_t fPosCorrectedTime;	// Pulse height corrected time
  Double_t fNegCorrectedTime;	// Pulse height corrected time
  Double_t fScinCorrectedTime;  // Time average corrected for position
                                // based on ADCs.
  Double_t fPosTOFCorrectedTime; // Times corrected for z position
  Double_t fNegTOFCorrectedTime; // using nominal beta

  Int_t fPaddleNumber;
  Double_t fPaddleCenter;

  Bool_t fHasCorrectedTimes;

  THcLADHodoPlane *fPlane;

private:
  THcLADHodoHit(const THcLADHodoHit &);
  THcLADHodoHit &operator=(const THcLADHodoHit &);

  ClassDef(THcLADHodoHit, 0)
};

#endif
