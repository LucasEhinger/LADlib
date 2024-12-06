#include "THcLADSpectrometer.h"
#include "THaNonTrackingDetector.h"
#include "TDatime.h"
#include "TList.h"
#include "THaTrack.h"
#include "THcGlobals.h"
#include "THcParmList.h"


#ifdef WITH_DEBUG
#include <iostream>
#endif

using namespace std;

THcLADSpectrometer::THcLADSpectrometer( const char* name, const char* description ) :
  THaApparatus( name, description )
{
  // default constructor 
  fTracks = 0;
  fNtracks = 0;
  fStagesDone = 0;
  fNonTrackingDetectors = new TList;
  fListInit = kFALSE;
}

//____________________________________________________________________
THcLADSpectrometer::~THcLADSpectrometer()
{
  // Destructor

  DefineVariables( kDelete );
  
}

//____________________________________________________________________
void THcLADSpectrometer::ListInit()
{

  fNonTrackingDetectors->Clear();
  TIter next(fDetectors);
  while( THaDetector* theDetector =
	 static_cast<THaDetector*>( next() )) {

    // We don't really have a tracking detector 
    // fTrackingDetectors->Add( theDetector );

    fNonTrackingDetectors->Add( theDetector );
  }

  fListInit = kTRUE;  
}

//____________________________________________________________________
Int_t THcLADSpectrometer::DefineVariables( EMode mode )
{
  if (mode == kDefine && fIsSetup) return kOK;
  fIsSetup = (mode == kDefine);

  return kOK;
}

//____________________________________________________________________
Int_t THcLADSpectrometer::ReadDatabase( const TDatime& date )
{

  DBRequest list[]={
    {"lpartmass",    &fPartMass, kDouble },
    {"lad_pcentral", &fPcentral, kDouble },
    {0}    
  };
  
  gHcParms->LoadParmValues((DBRequest*)&list);

  return kOK;

}

//____________________________________________________________________
/*
Int_t THcLADSpectrometer::CoarseTrack()
{

  //  THaSpectrometer::CoarseTrack();
  TIter next( fTrackingDetectors );
  while( THaTrackingDetector* theTrackDetector =
	 static_cast<THaTrackingDetector*>( next() )) {
#ifdef WITH_DEBUG
    if( fDebug >1 ) cout << "Call CoarseProcess() for "
			 << theTrackDetector->GetName() << "... ";
#endif
    theTrackDetector->CoarseTrack( *fTracks );
#ifdef WITH_DEBUG
    if ( fDebug>1 ) cout << "done.\n";
#endif
  }

  fStatesDone |= kCoarseTrack;
  return 0;

}
*/
//____________________________________________________________________
Int_t THcLADSpectrometer::Reconstruct()
{

  // Fine Process

  TIter next( fNonTrackingDetectors );
  while( THaNonTrackingDetector* theNonTrackDetector =
	 static_cast<THaNonTrackingDetector*>( next() )) {
#ifdef WITH_DEBUG
    if( fDebug > 1 ) cout << "Call FineProcess() for"
			  << theNonTrackDetector->GetName() << "... ";
#endif
    theNonTrackDetector->FineProcess( *fTracks );
#ifdef WITH_DEBUG
    if( fDebug > 1 ) cout << "done.\n";
#endif
  }

  fStagesDone |= kReconstruct;
  return 0;
}

//____________________________________________________________________
Int_t THcLADSpectrometer::CoarseReconstruct()
{

  if( !fListInit )
    ListInit();

  TIter next( fNonTrackingDetectors );
  while( THaNonTrackingDetector* theNonTrackDetector =
	 static_cast<THaNonTrackingDetector*>( next() )) {
#ifdef WITH_DEBUG
    if( fDebug >1 ) cout << "Call CoarseProcess() for "
			 << theNonTrackDetector->GetName() << "... ";
#endif
    theNonTrackDetector->CoarseProcess( *fTracks );
#ifdef WITH_DEBUG
    if ( fDebug>1 ) cout << "done.\n";
#endif
  }
  
  fStagesDone |= kCoarseRecon;
  return 0;

}

//____________________________________________________________________
THcLADSpectrometer::THcLADSpectrometer() {}

ClassImp(THcLADSpectrometer)

