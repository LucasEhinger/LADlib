#ifndef THcLADGEMCluster_h
#define THcLADGEMCluster_h

#include "TObject.h"

class THcLADGEMCluster : public TObject {
 public:
  THcLADGEMCluster();
  virtual ~THcLADGEMCluster() {}
  virtual void Clear( Option_t* opt );
  virtual void Print( Option_t* opt="" ) const;
  /*
  void AddHit(THcLADGEMHit* hit) {
    fHits.push_back(hit);
    fNHits++;
  }

  THcLADGEMHits* GetHit(Int_t ihit) { return fHits[ihit]; }
  */

  Int_t    GetNHits() { return fNHits; }
  Int_t    GetNStrips() { return fNStrips; }
  Int_t    GetStripLow() { return fStripLow;}
  Int_t    GetStripHigh() { return fStripHigh; }
  Int_t    GetStripMax() { return fStripMax; }

  Double_t GetADCsum() { return fADCsum; }
  Double_t GetPos()     { return fPos; }
  Double_t GetPosMax()     { return fPosMax; }
  Double_t GetMoments()     { return fMoments; }
  Double_t GetPosDiff()     { return fPosDiff; }
  Double_t GetE()     { return fE; }
  Double_t GetTime()     { return fTime; }
  Double_t GetTimeFit()     { return fTimeFit; }

  Double_t GetMPD() { return fMPD; }
  Double_t GetAPV() { return fAPV; }
  Double_t GetLayer() { return fLayer; }
  Double_t GetAxis() { return fAxis; }

  void SetMode(int this_value) { fClusteringFlag = this_value; }
  void SetNStrips(int this_value) { fNStrips = this_value; }
  void SetStripLow(int this_value) { fStripLow = this_value; }
  void SetStripHigh(int this_value) { fStripHigh = this_value; }
  void SetStripMax(int this_value) { fStripMax = this_value; }  
  void SetStrips(int nstrip, int striplo, int striphi, int stripmax)
  { fNStrips = nstrip, fStripLow = striplo, fStripHigh = striphi, fStripMax = stripmax; }

  void SetPosition(double this_value) { fPos = this_value; }
  void SetPosMaxStrip(double this_value) { fPosMax = this_value; }
  void SetPosSigma(double this_value) { fPosSigma = this_value; }
  void SetMoments(double this_value) { fMoments = this_value; }

  void SetTime(double this_value) { fTime = this_value; }
  void SetTimeFit(double this_value) { fTimeFit = this_value; }
  void SetTime(double value1, double value2) { fTime = value1; fTSigma = value2; }
  void SetADCsum(double this_value) { fADCsum = this_value; }
  void SetSampMax(int this_value) { fSampMax = this_value; }

  void SetLayer(int this_value ) { fLayer = this_value; }
  void SetMPD(int this_value ) { fMPD = this_value; }
  void SetAPV(int this_value ) { fAPV = this_value; }
  void SetAxis(int this_value) { fAxis = this_value; }

 protected:

  Int_t    fClusteringFlag;

  Int_t    fNStrips;
  Int_t    fStripLow;
  Int_t    fStripHigh;
  Int_t    fStripMax;

  Int_t    fLayer;
  Int_t    fMPD;
  Int_t    fAPV; // APV adc id
  Int_t    fAxis; // U/V, X/Y 

  Double_t fPos; // ADC weighted mean coordinate along the direction measured by the strip
  Double_t fPosMax;  // Max strip position
  Double_t fPosDiff;
  Double_t fPosSigma; // RMS coordinate deviation from the mean 
  Double_t fMoments;      // clustermoments
  Int_t    fSampMax; // time sample in which peak of cluster summed ADC values occurs
  Double_t fADCsum;

  Int_t    fNHits;
  Double_t fE;
  Double_t fTime;   // mean time
  Double_t fTSigma; 
  Double_t fTimeFit;

  //  std::vector<THcLADGEMHit*> fHits;

  ClassDef(THcLADGEMCluster,0)
};

#endif
