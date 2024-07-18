#ifndef THcLADHodoHit_h
#define THcLADHodoHit_h

#include "THcLADHodoPlane.h"
#include "THcRawHodoHit.h"
#include "TObject.h"

class THcLADHodoHit : public TObject {
public:
  THcLADHodoHit(Int_t toptdc, Int_t btmtdc, Double_t topadc, Double_t btmadc, Int_t ipad, THcLADHodoPlane *hp);
  virtual ~THcLADHodoHit();

  Double_t GetTopADC() const { return fTopADC_Ped; }
  Double_t GetBtmADC() const { return fBtmADC_Ped; }
  Double_t GetTopADCpeak() const { return fTopADC_Peak; }
  Double_t GetBtmADCpeak() const { return fBtmADC_Peak; }
  Double_t GetTopADCtime() const { return fTopADC_Time; }
  Double_t GetBtmADCtime() const { return fBtmADC_Time; }
  Int_t GetTopTDC() const { return fTopTDC; }
  Int_t GetBtmTDC() const { return fBtmTDC; }
  Int_t GetPaddleNumber() const { return fPaddleNumber; }
  Int_t GetPaddleCenter() const { return fPaddleCenter; }

  void SetPaddleCenter(Double_t padcenter) { fPaddleCenter = padcenter; }
  void SetTopADCpeak(Double_t adc) { fTopADC_Peak = adc; }
  void SetBtmADCpeak(Double_t adc) { fBtmADC_Peak = adc; }
  void SetTopADCtime(Double_t ptime) { fTopADC_Time = ptime; }
  void SetBtmADCtime(Double_t ptime) { fBtmADC_Time = ptime; }
  void SetCorrectedTimes(Double_t top, Double_t btm);
  void SetCorrectedTimes(Double_t top, Double_t btm, Double_t toptof, Double_t btmtof, Double_t timeave);

protected:
  Int_t fTopTDC;
  Int_t fBtmTDC;
  Double_t fTopADC_Ped; // Pedestal subtracted ADC
  Double_t fBtmADC_Ped;
  Double_t fTopADC_Peak; // ADC peak amplitude
  Double_t fBtmADC_Peak; // ADC peak amplitude
  Double_t fTopADC_Time;
  Double_t fBtmADC_Time;

  Double_t fTopCorrectedTime;	// Pulse height corrected time
  Double_t fBtmCorrectedTime;	// Pulse height corrected time
  Double_t fScinCorrectedTime;  // Time average corrected for position
                                // based on ADCs.
  Double_t fTopTOFCorrectedTime; // Times corrected for z position
  Double_t fBtmTOFCorrectedTime; // using nominal beta

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
