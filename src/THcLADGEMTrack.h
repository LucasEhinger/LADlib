#ifndef THcLADGEMTrack_h
#define THcLADGEMTrack_h

// LAD GEM Track Object

#include "TVector3.h"
#include "TObject.h"
#include "THcLADGEM.h"

class GEM2DHits {
  public:
  GEM2DHits() : 
    layer(-1), posX(kBig), posY(kBig), posZ(kBig),
    TimeMean(kBig), TimeDiff(kBig), TimeCorr(kBig),
    IsGoodHit(kFALSE), ADCMean(kBig), ADCasym(kBig),
    trackID(-1), clusID{-1,-1} {}
  virtual ~GEM2DHits() = default;
  void Set(Int_t _layer, Double_t _x, Double_t _y, Double_t _z,
	   Double_t _t, Double_t _dt, Double_t _ct, Bool_t _flag,
	   Double_t _adc, Double_t _asym) {
    layer = _layer; posX = _x; posY = _y; posZ = _z;
    TimeMean = _t; TimeDiff = _dt; TimeCorr = _ct;
    IsGoodHit = _flag; ADCMean = _adc; ADCasym = _asym;
  }
  void SetClusterIDs(Int_t id1, Int_t id2) {
    clusID[0] = id1;
    clusID[1] = id2;
  }
  void SetTrackID(Int_t id) { trackID = id; }

  Int_t    layer;
  Double_t posX;
  Double_t posY;
  Double_t posZ;
  Double_t TimeMean; // average time
  Double_t TimeDiff;
  Double_t TimeCorr;
  Bool_t   IsGoodHit;
  Double_t ADCMean; // average adc sum
  Double_t ADCasym;
  Int_t    trackID;   // associated Track ID
  Int_t    clusID[2]; // associated cluster IDs
};

class THcLADGEMTrack : public TObject {
 public:
  THcLADGEMTrack(Int_t nlayers);
  virtual ~THcLADGEMTrack();
  virtual void Clear( Option_t* opt="" );

  // THcLADGEMHit* GetHit(Int_t ihit) {};  

  Int_t    GetTrackID()      const { return fTrackID; }
  // Space points
  Int_t    GetNSpacePoints() const { return fNSp; }
  Int_t*   GetClusID_Sp1()   const { return fSp[0].clusID; }
  Int_t*   GetClusID_Sp2()   const { return fSp[1].clusID; }
  Double_t GetADCMean_Sp1()  const { return fSp[0].ADCMean; }
  Double_t GetADCMean_Sp2()  const { return fSp[1].ADCMean; }
  Double_t GetADCasym_Sp1()  const { return fSp[0].ADCasym; }
  Double_t GetADCasym_Sp2()  const { return fSp[1].ADCasym; }
  Double_t GetX1()           const { return fSp[0].posX; }
  Double_t GetY1()           const { return fSp[0].posY; }
  Double_t GetZ1()           const { return fSp[0].posZ; }
  Double_t GetX2()           const { return fSp[1].posX; }
  Double_t GetY2()           const { return fSp[1].posY; }
  Double_t GetZ2()           const { return fSp[1].posZ; }
  // Track quantities
  Double_t GetProjVz()       const { return fProjVz; } // projected z-vertex
  Double_t GetD0()           const { return fD0; }
  Double_t GetT()            const { return fT; }
  Double_t GetdT()           const { return fTdiff; }

  void     SetTrackID(int itrk) { fTrackID = itrk; }
  void     SetD0(Double_t d0) { fD0 = d0; }
  void     SetZVertex(Double_t vz) { fProjVz = vz; }
  void     SetTime(Double_t t, Double_t dt) { fT = t; fTdiff = dt; }
    
  GEM2DHits GetSpacePoint(int isp) { return fSp[isp]; }
  virtual void AddSpacePoint(GEM2DHits &sp);

 protected:

  Int_t    fNSp;
  Double_t fProjVz; 
  Double_t fD0;
  Int_t    fTrackID;
  Double_t fT;
  Double_t fTdiff;

  GEM2DHits* fSp;

 private:
  THcLADGEMTrack( const THcLADGEMTrack& );
  THcLADGEMTrack& operator=( const THcLADGEMTrack& );

  ClassDef(THcLADGEMTrack,0)
};

#endif
