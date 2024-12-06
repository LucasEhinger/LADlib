#ifndef THcLADSpectrometer_h
#define THcLADSpectrometer_h

/////////////////////////////////////////////////
//
// Skeleton class for HallC LAD spectrometer
//
/////////////////////////////////////////////////

#include "THaAnalysisObject.h"
#include "THaApparatus.h"
#include "TClonesArray.h"
#include "TMath.h"

class THcLADSpectrometer : public THaApparatus {

 public:

  THcLADSpectrometer( const char* name, const char* description );
  virtual ~THcLADSpectrometer();
  
  virtual Int_t CoarseReconstruct();
  virtual Int_t Reconstruct();

  // Mass of nominal detected particle type
  Double_t GetParticleMass() const { return fPartMass; }
  Double_t GetBetaAtPcentral() const { return fPcentral / TMath::Sqrt(fPcentral * fPcentral + fPartMass * fPartMass); }

  virtual Int_t DefineVariables( EMode mode = kDefine );
  virtual Int_t ReadDatabase( const TDatime& date );

  enum EStagesDone {
    kCoarseTrack = BIT(0),
    kCoarseRecon = BIT(1),
    kTracking    = BIT(2),
    kReconstruct = BIT(3)
  };

 protected:
  Double_t fPartMass;
  Double_t fPcentral;
  Int_t fNtracks;
  //  THcLADHodoscope *fHodo; // we don't quite need these
  //  THcLADGEM *fGEM;

  TClonesArray* fTracks;

  UInt_t fStagesDone;

  TList* fNonTrackingDetectors;
  //  TList* fTrackingDetectors;

 private:
  Bool_t fListInit;
  void   ListInit();

  THcLADSpectrometer();
  ClassDef(THcLADSpectrometer,0)
    
};

#endif /* THcLADSpectrometer_h */
