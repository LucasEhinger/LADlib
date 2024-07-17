#ifndef THcLADHodoPlane_h
#define THcLADHodoPlane_h

#include "TClonesArray.h"
#include "THaSubDetector.h"
#include <vector>

using namespace std;

class THaEvData;
class THaSignalHit;

// rough class outline for LAD Hodoscope plane
// contains xx number of paddles per plane
// For each hodoscope detector, we will define it as a separate detector
// since they are not really placed together

class THcLADHodoPlane : public THaSubDetector {

public:
  THcLADHodoPlane(const char *name, const char *description, const Int_t planenum, THaDetectorBase *parent = nullptr);
  virtual ~THcLADHodoPlane();

  virtual void Clear(Option_t *opt = "");
  virtual Int_t Decode(const THaEvData &);
  virtual EStatus Init(const TDatime &date);

  virtual Int_t CoarseProcess(TClonesArray &tracks);
  virtual Int_t FineProcess(TClonesArray &tracks);
  virtual Int_t ProcessHits(TClonesArray *rawhits, Int_t nexthit);

  // these are called in THcLADHodoscope::Decode, to analyze pedestal events
  virtual Int_t AccumulatePedestals(TClonesArray *rawhits, Int_t nexthit);
  virtual void CalculatePedestals();

  // Some getter/setter functions here

  Int_t GetNelem() { return fNelem; };         // return number of paddles in this plane
  Int_t GetNScinHits() { return fNScinHits; }; // Get # hits in plane (that pass min/max TDC cuts)
  Int_t GetNGoodHits() { return fNGoodHits; }; // Get # hits in plane (used in determining focal plane time)

  TClonesArray *GetHits() { return fHodoHits; }

protected:
  TClonesArray *fHodoHits;

  TClonesArray *frPosTdcTimeRaw;
  TClonesArray *frPosTdcTime;
  TClonesArray *frNegTdcTimeRaw;
  TClonesArray *frNegTdcTime;

  Double_t fPosTdcRefTime;
  Double_t fNegTdcRefTime;
  Double_t fPosTdcRefDiffTime;
  Double_t fNegTdcRefDiffTime;

  Double_t fPosAdcRefTime;
  Double_t fNegAdcRefTime;
  Double_t fPosAdcRefDiffTime;
  Double_t fNegAdcRefDiffTime;

  // Counters
  Int_t fTotNumTdcHits;
  Int_t fTotNumPosTdcHits;
  Int_t fTotNumNegTdcHits;

  Int_t fTotNumAdcHits;
  Int_t fTotNumPosAdcHits;
  Int_t fTotNumNegAdcHits;

  Int_t fTotNumGoodPosAdcHits;
  Int_t fTotNumGoodNegAdcHits;
  Int_t fTotNumGoodAdcHits;

  Int_t fTotNumGoodPosTdcHits;
  Int_t fTotNumGoodNegTdcHits;
  Int_t fTotNumGoodTdcHits;

  vector<Int_t> fNumGoodPosAdcHits;
  vector<Int_t> fNumGoodNegAdcHits;
  vector<Int_t> fNumGoodPosTdcHits;
  vector<Int_t> fNumGoodNegTdcHits;

  vector<Double_t> fGoodPosAdcPed;
  vector<Double_t> fGoodPosAdcMult;
  vector<Double_t> fGoodPosAdcHitUsed;
  vector<Double_t> fGoodPosAdcPulseInt;
  vector<Double_t> fGoodPosAdcPulseAmp;
  vector<Double_t> fGoodPosAdcPulseTime;
  vector<Double_t> fGoodPosAdcTdcDiffTime;

  vector<Double_t> fGoodNegAdcPed;
  vector<Double_t> fGoodNegAdcMult;
  vector<Double_t> fGoodNegAdcHitUsed;
  vector<Double_t> fGoodNegAdcPulseInt;
  vector<Double_t> fGoodNegAdcPulseAmp;
  vector<Double_t> fGoodNegAdcPulseTime;
  vector<Double_t> fGoodNegAdcTdcDiffTime;

  TClonesArray *frPosAdcPedRaw;
  TClonesArray *frPosAdcPed;
  TClonesArray *frPosAdcPulseIntRaw;
  TClonesArray *frPosAdcPulseInt;
  TClonesArray *frPosAdcPulseAmpRaw;
  TClonesArray *frPosAdcPulseAmp;
  TClonesArray *frPosAdcPulseTimeRaw;
  TClonesArray *frPosAdcPulseTime;
  TClonesArray *frPosAdcErrorFlag;

  TClonesArray *frPosAdcSampPedRaw;
  TClonesArray *frPosAdcSampPed;
  TClonesArray *frPosAdcSampPulseIntRaw;
  TClonesArray *frPosAdcSampPulseInt;
  TClonesArray *frPosAdcSampPulseAmpRaw;
  TClonesArray *frPosAdcSampPulseAmp;
  TClonesArray *frPosAdcSampPulseTimeRaw;
  TClonesArray *frPosAdcSampPulseTime;

  TClonesArray *frNegAdcPedRaw;
  TClonesArray *frNegAdcPed;
  TClonesArray *frNegAdcPulseIntRaw;
  TClonesArray *frNegAdcPulseInt;
  TClonesArray *frNegAdcPulseAmpRaw;
  TClonesArray *frNegAdcPulseAmp;
  TClonesArray *frNegAdcPulseTimeRaw;
  TClonesArray *frNegAdcPulseTime;
  TClonesArray *frNegAdcErrorFlag;

  TClonesArray *frNegAdcSampPedRaw;
  TClonesArray *frNegAdcSampPed;
  TClonesArray *frNegAdcSampPulseIntRaw;
  TClonesArray *frNegAdcSampPulseInt;
  TClonesArray *frNegAdcSampPulseAmpRaw;
  TClonesArray *frNegAdcSampPulseAmp;
  TClonesArray *frNegAdcSampPulseTimeRaw;
  TClonesArray *frNegAdcSampPulseTime;

  vector<Double_t> fPosAdcSampWaveform;
  vector<Double_t> fNegAdcSampWaveform;

  // Hodoscopoe "GOOD" TDC Variables
  vector<Double_t> fGoodPosTdcTimeUnCorr;
  vector<Double_t> fGoodPosTdcTimeCorr;
  vector<Double_t> fGoodPosTdcTimeTOFCorr;

  vector<Double_t> fGoodNegTdcTimeUnCorr;
  vector<Double_t> fGoodNegTdcTimeCorr;
  vector<Double_t> fGoodNegTdcTimeTOFCorr;

  // Time Walk Corrected
  vector<Double_t> fGoodPosTdcTimeWalkCorr;
  vector<Double_t> fGoodNegTdcTimeWalkCorr;
  vector<Double_t> fGoodDiffDistTrack;

  TClonesArray *frPosTDCHits;
  TClonesArray *frPosADCHits;
  TClonesArray *frPosADCSums;
  TClonesArray *frPosADCPeds;
  TClonesArray *frNegTDCHits;
  TClonesArray *frNegADCHits;
  TClonesArray *frNegADCSums;
  TClonesArray *frNegADCPeds;

  Int_t fPlaneNum;
  Int_t fNelem;

  Int_t fNScinHits;     /* number of hits in plane (that pass min/max TDC cuts) */
  Int_t fNGoodHits;     /* number of hits in plane (used in determining focal plane time) */
  Int_t fNScinGoodHits; /* number of hits for which both ends of the paddle fired in time -- diff from fNGoodHits? */
  Double_t fpTime;

  Double_t fHitDistance;
  Double_t fScinXPos;
  Double_t fScinYPos;
  Double_t fTrackXPosition;
  Double_t fTrackYPosition;
  Double_t *fPosCenter; /* array with centers for all scintillators in the plane */

  Int_t fDebugAdc;
  Int_t fADCMode;
  Int_t fADCDiagCut;
  Int_t fCosmicFlag;
  Double_t fSampThreshold;
  Int_t fSampNSA;
  Int_t fSampNSAT;
  Int_t fSampNSB;
  Int_t fOutputSampWaveform;
  Int_t fUseSampWaveform;

  Double_t fDzpos;
  Double_t fHodoSlop;
  Double_t fPosLeft;  /* NOTE: "left" = "bottom" for a Y scintillator */
  Double_t fPosRight; /* NOTE: "right" = "top" for a Y scintillator */
  Int_t fTdcOffset;
  Double_t fAdcTdcOffset;
  Double_t fScinTdcMin;
  Double_t fScinTdcMax;
  Double_t fScinTdcToTime;
  Double_t fBetaNominal;

  Double_t *fHodoPosAdcTimeWindowMin;
  Double_t *fHodoPosAdcTimeWindowMax;
  Double_t *fHodoNegAdcTimeWindowMin;
  Double_t *fHodoNegAdcTimeWindowMax;

  // Hodoscope Calib Parameters
  Int_t fTofUsingInvAdc;
  Double_t *fHodoVelLight;
  Double_t *fHodoPosInvAdcOffset;
  Double_t *fHodoNegInvAdcOffset;
  Double_t *fHodoPosInvAdcAdc;
  Double_t *fHodoNegInvAdcAdc;
  Double_t *fHodoPosInvAdcLinear;
  Double_t *fHodoNegInvAdcLinear;

  // Time-Walk Parameters
  Double_t *fHodoVelFit;
  Double_t *fHodoCableFit;
  Double_t *fHodo_LCoeff;
  Double_t *fHodoPos_c1;
  Double_t *fHodoNeg_c1;
  Double_t *fHodoPos_c2;
  Double_t *fHodoNeg_c2;
  Double_t fTdc_Thrs;

  Double_t tw_corr_pos;
  Double_t tw_corr_neg;

  // Pedestal calculations
  Int_t fNPedestalEvents; /* Number of pedestal events */
  Int_t fMinPeds;         /* Only analyze/update if num events > */
  Int_t *fPosPedSum;      /* Accumulators for pedestals */
  Int_t *fPosPedSum2;
  Int_t *fPosPedLimit;
  Int_t *fPosPedCount;
  Int_t *fNegPedSum;
  Int_t *fNegPedSum2;
  Int_t *fNegPedLimit;
  Int_t *fNegPedCount;

  Double_t *fPosPed;
  Double_t *fPosSig;
  Double_t *fPosThresh;
  Double_t *fNegPed;
  Double_t *fNegSig;
  Double_t *fNegThresh;

  enum { kADCStandard = 0, kADCDynamicPedestal, kADCSampleIntegral, kADCSampIntDynPed };

  // Geometry parameters
  Double_t fSpacing; /* paddle spacing */
  Double_t fSize;    /* paddle size */
  Double_t fZpos;    /* z position */

  virtual Int_t ReadDatabase(const TDatime &date);
  virtual Int_t DefineVariables(EMode mode = kDefine);
  virtual void InitializePedestals();

  class HodoCluster {
  public:
    Int_t fClusterNumber; // cluster index
    Double_t fClusterPos;
    Double_t fClusterSize;
  };

  vector<HodoCluster> fCluster;

public:
  vector<HodoCluster> GetClusters() { return fCluster; }

  ClassDef(THcLADHodoPlane, 0);
};

#endif /* THcLADHodPlane_h */
