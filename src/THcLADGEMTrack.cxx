#include "THcLADGEMTrack.h"

//_________________________________________________________
THcLADGEMTrack::THcLADGEMTrack(Int_t nlayers)
{
  fProjVz = -999.;
  fProjVy = -999.;
  fProjVx = -999.;
  fNSp = 0;
  fD0 = -999;
  fhasGoodD0 = kFALSE;
  fSp = new GEM2DHits[nlayers];
  chisq = -999.;
  ftheta = -999.;
  fphi = -999.;

}

//_________________________________________________________
THcLADGEMTrack::~THcLADGEMTrack()
{
  delete [] fSp; fSp = nullptr;
}

//_________________________________________________________
void THcLADGEMTrack::AddSpacePoint(GEM2DHits sp)
{
  fSp[fNSp] = sp;
  fNSp++;

}

//_________________________________________________________

void THcLADGEMTrack::Clear( Option_t* opt )
{
  fNSp = 0;
}

ClassImp(THcLADGEMTrack)
