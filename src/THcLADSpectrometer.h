#ifndef THcLADSpectrometer_h
#define THcLADSpectrometer_h

/////////////////////////////////////////////////
//
// Skeleton class for HallC LAD spectrometer
//
/////////////////////////////////////////////////

#include "THaSpectrometer.h"
#include "THcLADGEM.h"
#include "THcLADHodoscope.h"

class THcLADSpectrometer : public THaSpectrometer {

public:
  THcLADSpectrometer(const char *name, const char *description);
  virtual ~THcLADSpectrometer();

  virtual void Clear(Option_t *opt = "");

  // In case we need to redefine these functions
  /*
  virtual Int_t CoarseReconstruct();
  virtual Int_t CoarseTrack();
  virtual Int_t Reconstruct();
  virtual Int_t Track();
  */

  // Mass of nominal detected particle type
  Double_t GetParticleMass() const { return fPartMass; }
  Double_t GetBetaAtPcentral() const { return fPcentral / TMath::Sqrt(fPcentral * fPcentral + fPartMass * fPartMass); }

  virtual Int_t FindVertices(TClonesArray &tracks);
  virtual Int_t TrackCalc();
  virtual Int_t Decode(const THaEvData &);

  virtual Int_t DefineVariables(EMode mode = kDefine);
  virtual Int_t ReadDatabase(const TDatime &date);
  virtual Int_t ReadRunDatabase(const TDatime &date);

protected:
  Double_t fPartMass;
  Int_t fNtracks;
  THcLADHodoscope *fHodo;
  THcLADGEM *fGEM;

  ClassDef(THcLADSpectrometer, 0)
};

#endif /* THcLADSpectrometer_h */
