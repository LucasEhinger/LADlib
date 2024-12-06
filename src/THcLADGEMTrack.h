#ifndef THcLADGEMTrack_h
#define THcLADGEMTrack_h

#include "TVector3.h"
#include "TObject.h"
#include "THcLADGEM.h"

// LAD GEM Track Object

class GEM2DHits {
  public:
   Int_t    layer;
   Double_t posX;
   Double_t posY;
   Double_t posZ;
   Double_t TimeMean; // average time
   Double_t TimeDiff;
   Double_t TimeCorr;
   Bool_t   IsGoodHit;
   //    Bool_t   Filtered;
   Double_t ADCMean; // average adc sum
   Double_t ADCasym;
};

class THcLADGEMTrack : public TObject {
 public:
  THcLADGEMTrack(Int_t nlayers);
  virtual ~THcLADGEMTrack();
  virtual void Clear( Option_t* opt="" );

  // THcLADGEMHit* GetHit(Int_t ihit) {};  

  Int_t    GetTrackID()      const { return fTrackID; }
  Int_t    GetNSpacePoints() const { return fNSp; }
  Double_t GetX1()   const { return fSp[0].posX; }
  Double_t GetY1()   const { return fSp[0].posY; }
  Double_t GetZ1()   const { return fSp[0].posZ; }
  Double_t GetX2()   const { return fSp[1].posX; }
  Double_t GetY2()   const { return fSp[1].posY; }
  Double_t GetZ2()   const { return fSp[1].posZ; }
  Double_t GetProjVz()       const { return fProjVz; } // projected z-vertex
  Double_t GetD0()           const { return fD0; }
  Double_t GetT()         const { return fT; }
  Double_t GetdT()         const { return fTdiff; }

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
