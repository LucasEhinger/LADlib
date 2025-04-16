#ifndef THcLADHodoscope_h
#define THcLADHodoscope_h

#include "TClonesArray.h"
#include "TH1F.h"
#include "THaNonTrackingDetector.h"
#include "THaSpectrometer.h"
#include "THcHitList.h"
#include "THcLADGEM.h"
#include "THcLADHodoHit.h"
#include "THcLADHodoPlane.h"
#include "THcRawHodoHit.h"

class THcGoodLADHit : public TObject {
public:
  Int_t plane[2];
  Int_t paddle[2];
  Int_t track_id[2];
  Double_t beta[2];
  Double_t delta_pos_trans[2];
  Double_t delta_pos_long[2];
  Double_t hit_time[2];
  Double_t hit_theta[2];
  Double_t hit_phi[2];
  Double_t hit_edep[2];

  static constexpr Double_t kbig = 1.0e38;

  THcGoodLADHit() {
    for (int i = 0; i < 2; ++i) {
      plane[i] = paddle[i] = track_id[i] =-1;
      beta[i] = delta_pos_trans[i] = delta_pos_long[i] = hit_time[i] = kbig;
      hit_theta[i] = hit_phi[i] = hit_edep[i] = kbig;
    }
  }
  virtual ~THcGoodLADHit() = default;

  void SetPlane(Int_t hit, Int_t value) { plane[hit] = value; }
  void SetPaddle(Int_t hit, Int_t value) { paddle[hit] = value; }
  void SetTrackID(Int_t hit, Int_t value) { track_id[hit] = value; }
  void SetBeta(Int_t hit, Double_t value) { beta[hit] = value; }
  void SetDeltaPosTrans(Int_t hit, Double_t value) { delta_pos_trans[hit] = value; }
  void SetDeltaPosLong(Int_t hit, Double_t value) { delta_pos_long[hit] = value; }
  void SetHitTime(Int_t hit, Double_t value) { hit_time[hit] = value; }
  void SetHitTheta(Int_t hit, Double_t value) { hit_theta[hit] = value; }
  void SetHitPhi(Int_t hit, Double_t value) { hit_phi[hit] = value; }
  void SetHitEdep(Int_t hit, Double_t value) { hit_edep[hit] = value; }

  // Int_t GetPlane(Int_t hit) const { return plane[hit]; }
  // Int_t GetPaddle(Int_t hit) const { return paddle[hit]; }
  // Int_t GetTrackID(Int_t hit) const { return track_id[hit]; }
  // Double_t GetBeta(Int_t hit) const { return beta[hit]; }
  // Double_t GetDeltaPosTrans(Int_t hit) const { return delta_pos_trans[hit]; }
  // Double_t GetDeltaPosLong(Int_t hit) const { return delta_pos_long[hit]; }
  // Double_t GetHitTime(Int_t hit) const { return hit_time[hit]; }
  // Double_t GetHitTheta(Int_t hit) const { return hit_theta[hit]; }
  // Double_t GetHitPhi(Int_t hit) const { return hit_phi[hit]; }
  // Double_t GetHitEdep(Int_t hit) const { return hit_edep[hit]; }

  Int_t GetPlaneHit0() const { return plane[0]; }
  Int_t GetPlaneHit1() const { return plane[1]; }

  Int_t GetPaddleHit0() const { return paddle[0]; }
  Int_t GetPaddleHit1() const { return paddle[1]; }

  Int_t GetTrackIDHit0() const { return track_id[0]; }
  Int_t GetTrackIDHit1() const { return track_id[1]; }

  Double_t GetBetaHit0() const { return beta[0]; }
  Double_t GetBetaHit1() const { return beta[1]; }

  Double_t GetDeltaPosTransHit0() const { return delta_pos_trans[0]; }
  Double_t GetDeltaPosTransHit1() const { return delta_pos_trans[1]; }

  Double_t GetDeltaPosLongHit0() const { return delta_pos_long[0]; }
  Double_t GetDeltaPosLongHit1() const { return delta_pos_long[1]; }

  Double_t GetHitTimeHit0() const { return hit_time[0]; }
  Double_t GetHitTimeHit1() const { return hit_time[1]; }

  Double_t GetHitThetaHit0() const { return hit_theta[0]; }
  Double_t GetHitThetaHit1() const { return hit_theta[1]; }

  Double_t GetHitPhiHit0() const { return hit_phi[0]; }
  Double_t GetHitPhiHit1() const { return hit_phi[1]; }

  Double_t GetHitEdepHit0() const { return hit_edep[0]; }
  Double_t GetHitEdepHit1() const { return hit_edep[1]; }
  ClassDef(THcGoodLADHit, 1)
};

class THcLADHodoscope : public THaNonTrackingDetector, public THcHitList {
public:
  THcLADHodoscope(const char *name, const char *description = "", THaApparatus *apparatus = nullptr);
  virtual ~THcLADHodoscope();

  virtual Int_t Decode(const THaEvData &);
  virtual EStatus Init(const TDatime &date);
  virtual void Clear(Option_t *opt = "");
  virtual Int_t End(THaRunBase *run = 0);

  virtual Int_t CoarseProcess(TClonesArray &tracks);
  virtual Int_t FineProcess(TClonesArray &tracks);

  Double_t DetermineTimePeak(Int_t FillFlag) { return 0.0; }; // TODO: Implement this function

  Int_t GetScinIndex(Int_t nPlane, Int_t nPaddle) { return fNPlanes * nPaddle + nPlane; }
  Int_t GetScinIndex(Int_t nSide, Int_t nPlane, Int_t nPaddle) {
    return nSide * fMaxHodoScin + fNPlanes * nPaddle + nPlane - 1;
  }

  THcLADHodoPlane *GetPlane(Int_t iplane) { return fPlanes[iplane]; }

  Int_t GetNPlanes() { return fNPlanes; }
  Int_t GetTdcOffset(Int_t iplane) const { return fTdcOffset[iplane]; }
  Double_t GetHodoSlop(Int_t iplane) const { return fHodoSlop[iplane]; }
  Double_t GetAdcTdcOffset(Int_t iplane) const { return fAdcTdcOffset[iplane]; }
  Double_t GetTdcMin() const { return fScinTdcMin; }
  Double_t GetTdcMax() const { return fScinTdcMax; }
  Double_t GetTdcToTime() const { return fScinTdcToTime; }

  Double_t GetHodoTopAdcTimeWindowMax(Int_t iii) const { return fHodoTopAdcTimeWindowMax[iii]; }
  Double_t GetHodoTopAdcTimeWindowMin(Int_t iii) const { return fHodoTopAdcTimeWindowMin[iii]; }
  Double_t GetHodoBtmAdcTimeWindowMax(Int_t iii) const { return fHodoBtmAdcTimeWindowMax[iii]; }
  Double_t GetHodoBtmAdcTimeWindowMin(Int_t iii) const { return fHodoBtmAdcTimeWindowMin[iii]; }

  Double_t GetHodoVelLight(Int_t iii) const { return fHodoVelLight[iii]; }

  // Time walk
  Double_t GetHodoVelFit(Int_t iii) const { return fHodoVelFit[iii]; }
  Double_t GetHodoCableFit(Int_t iii) const { return fHodoCableFit[iii]; }
  Double_t GetHodoLCoeff(Int_t iii) const { return fHodo_LCoeff[iii]; }
  Double_t GetHodoTop_c1(Int_t iii) const { return fHodoTop_c1[iii]; }
  Double_t GetHodoBtm_c1(Int_t iii) const { return fHodoBtm_c1[iii]; }
  Double_t GetHodoTop_c2(Int_t iii) const { return fHodoTop_c2[iii]; }
  Double_t GetHodoBtm_c2(Int_t iii) const { return fHodoBtm_c2[iii]; }
  Double_t GetTDCThrs() const { return fTdc_Thrs; }

  TClonesArray *GetLADGoodHits() { return fGoodLADHits; }

protected:
  TH1F *hTime;

  Int_t fNPlanes;
  Int_t fNHits;
  Int_t *fNPaddle;
  Double_t fTofTolerance;
  Int_t fNumPlanesBetaCalc; // Number of planes to use in beta calculation
  Int_t fMaxHodoScin;
  Int_t fCosmicFlag;
  Double_t fScinTdcMin, fScinTdcMax;
  Double_t fScinTdcToTime;
  Int_t *fTdcOffset;
  Double_t *fAdcTdcOffset;
  Double_t *fHodoSlop;
  Int_t fIsMC;
  Double_t fTrackToleranceLong;
  Double_t fTrackToleranceTrans;

  // Output variables
  static const Int_t MAXGOODHITS = 500;
  Int_t goodhit_n;
  TClonesArray *fGoodLADHits;
  Int_t num_unique_hits;
  Int_t num_unique_good_hits;
  Double_t *fHodoVelLight;

  Double_t fStartTime;
  Double_t fFPTimeAll;
  Double_t *fFPTime; // [fNPlanes] Array

  THaApparatus *fSpectro;
  THcLADGEM *fGEM;

  struct TOFPInfo {
    Bool_t onTrack;
    Bool_t keep_pos;
    Bool_t keep_neg;
    Double_t time_pos;      // Times also corrected for particle
    Double_t time_neg;      // flight time
    Double_t scin_pos_time; // Times corrected for position on
    Double_t scin_neg_time; // the bar
    Double_t pathp;
    Double_t pathn;
    Double_t zcor;
    Double_t scinTrnsCoord;
    Double_t scinLongCoord;
    Int_t planeIndex;
    Int_t hitNumInPlane;
    THcLADHodoHit *hit;
    TOFPInfo()
        : onTrack(kFALSE), keep_pos(kFALSE), keep_neg(kFALSE), time_pos(-99.0), time_neg(-99.0), scin_pos_time(0.0),
          scin_neg_time(0.0) {}
  };
  std::vector<TOFPInfo> fTOFPInfo;

  struct TOFCalc {
    Int_t hit_paddle;
    Int_t pindex; // Plane index
    Int_t good_raw_pad;
    Bool_t good_scin_time;
    Bool_t good_tdc_pos;
    Bool_t good_tdc_neg;
    Double_t scin_time;
    Double_t scin_time_fp;
    Double_t scin_sigma;
    Double_t dedx;
    TOFCalc() : good_scin_time(kFALSE), good_tdc_pos(kFALSE), good_tdc_neg(kFALSE) {}
  };
  std::vector<TOFCalc> fTOFCalc;

  std::vector<std::vector<Double_t>> fdEdX; // Vector over track #
  std::vector<Int_t> fNScinHit;             // # scins hit for the track

  Double_t *fSumPlaneTime; // [fNPlanes]
  Int_t *fNScinHits;       // [fNPlanes]
  Int_t *fNPlaneTime;      // [fNPlanes]
  Bool_t *fGoodPlaneTime;  // [fNPlanes]

  struct GoodFlags {
    Bool_t onTrack;
    Bool_t goodScinTime;
    Bool_t goodTdcNeg;
    Bool_t goodTdcPos;
    GoodFlags() : onTrack(false), goodScinTime(false), goodTdcNeg(false), goodTdcPos(false) {}
  };
  std::vector<std::vector<std::vector<GoodFlags>>> fGoodFlags;
  //

  // Time walk
  Double_t *fHodoVelFit;
  Double_t *fHodoCableFit;
  Double_t *fHodo_LCoeff;
  Double_t *fHodoTop_c1;
  Double_t *fHodoBtm_c1;
  Double_t *fHodoTop_c2;
  Double_t *fHodoBtm_c2;
  Double_t fTdc_Thrs;

  Int_t fAnalyzePedestals;

  Double_t *fHodoBtmAdcTimeWindowMin; // per element? per plane?
  Double_t *fHodoBtmAdcTimeWindowMax;
  Double_t *fHodoTopAdcTimeWindowMin;
  Double_t *fHodoTopAdcTimeWindowMax;

  // FIXME. Neither of these are used or initialized.
  Double_t fPartMass;    // Nominal particle mass
  Double_t fBetaNominal; // Beta for central ray of nominal particle type

  THcLADHodoPlane **fPlanes;
  char **fPlaneNames;

  Bool_t *fPresentP;

  virtual Int_t DefineVariables(EMode mode = kDefine);
  virtual Int_t ReadDatabase(const TDatime &date);
  void Setup(const char *name, const char *description);

  ClassDef(THcLADHodoscope, 0)
};

#endif /* THcLADHodoscope_h */
