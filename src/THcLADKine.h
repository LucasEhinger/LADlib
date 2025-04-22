#ifndef ROOT_THcLADKine
#define ROOT_THcLADKine
//////////////////////////////////////////////////////////////////////////

#include "THaSpectrometer.h" // Include for THaSpectrometer
#include "THcLADGEM.h"       // Include for THcLADGEM
#include "THcLADHodoscope.h" // Include for THcLADHodoscope
#include "THcPrimaryKine.h"
#include "THcReactionPoint.h" // Include for THcReactionPoint

class THcReactionPoint;
class THcLADHodoscope;
class THCLADGEM;

class THcLADKine : public THcPrimaryKine {
public:
  THcLADKine(const char *name, const char *description = "", const char *spectro = "", const char *primary_kine = "",
             const char *vertex_module = "");

  virtual ~THcLADKine();
  virtual EStatus Init(const TDatime &run_time);
  virtual Int_t ReadDatabase(const TDatime &date);

  // not sure if we ever want to override this..
  virtual Int_t Process(const THaEvData &);
  virtual void Clear(Option_t *opt = "");
  void SetApparatus(const char *name);

protected:
  TString fSpecName;
  THaSpectrometer *fSpec;
  THcLADGEM *fGEM;
  THcLADHodoscope *fHodoscope;
  TString fVertexModuleName;
  THcReactionPoint *fVertexModule;

  Int_t MAXGOODHITS = 500;
  Int_t goodhit_n;
  TClonesArray *fGoodLADHits;
  TClonesArray *fGEMTracks;
  Double_t fD0Cut_wVertex;
  Double_t fD0Cut_noVertex;
  Double_t fmax_dTrans_match;
  Double_t fmax_dLong_match;
  Int_t fNfixed_z;
  Double_t *fFixed_z;
  Double_t fglobal_time_offset;
  virtual Int_t DefineVariables(EMode mode = kDefine);

  ClassDef(THcLADKine, 0)
};

#endif /* ROOT_THcLADKine */