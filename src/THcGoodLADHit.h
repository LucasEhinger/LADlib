#ifndef THcGoodLADHit_H
#define THcGoodLADHit_H
#include "TClonesArray.h"
#include "TObject.h"

class THcGoodLADHit : public TObject {
public:
  THcGoodLADHit() {
    for (int i = 0; i < 2; ++i) {
      plane[i] = paddle[i] = track_id[i] = -1;
      hit_tof[i] = delta_pos_trans[i] = delta_pos_long[i] = hit_time[i] = hit_theta[i] = hit_phi[i] = hit_edep[i] =
          hit_yPos[i] = hit_yPos[i] = hit_alpha[i] = hit_tof[i] = 1e30;
    }
  };
  virtual ~THcGoodLADHit() = default;

  void CheckHitIndex(Int_t hit) const {
    if (hit < 0 || hit > 1) {
      throw std::out_of_range("Invalid hit index. Must be 0 or 1.");
    }
  }

  void SetPlane(Int_t hit, Int_t value) {
    CheckHitIndex(hit);
    plane[hit] = value;
  }
  void SetPaddle(Int_t hit, Int_t value) {
    CheckHitIndex(hit);
    paddle[hit] = value;
  }
  void SetTrackID(Int_t hit, Int_t value) {
    CheckHitIndex(hit);
    track_id[hit] = value;
  }
  void SetDeltaPosTrans(Int_t hit, Double_t value) {
    CheckHitIndex(hit);
    delta_pos_trans[hit] = value;
  }
  void SetDeltaPosLong(Int_t hit, Double_t value) {
    CheckHitIndex(hit);
    delta_pos_long[hit] = value;
  }
  void SetHitTime(Int_t hit, Double_t value) {
    CheckHitIndex(hit);
    hit_time[hit] = value;
  }
  void SetHitTheta(Int_t hit, Double_t value) {
    CheckHitIndex(hit);
    hit_theta[hit] = value;
  }
  void SetHitPhi(Int_t hit, Double_t value) {
    CheckHitIndex(hit);
    hit_phi[hit] = value;
  }
  void SetHitEdep(Int_t hit, Double_t value) {
    CheckHitIndex(hit);
    hit_edep[hit] = value;
  }
  void SetHitAlpha(Int_t hit, Double_t value) {
    CheckHitIndex(hit);
    hit_alpha[hit] = value;
  }
  void SetHitYPos(Int_t hit, Double_t value) {
    CheckHitIndex(hit);
    hit_yPos[hit] = value;
  }
  void SetHitTOF(Int_t hit, Double_t value) {
    CheckHitIndex(hit);
    hit_tof[hit] = value;
  }

  void SetHitBeta(Int_t hit, Double_t value) {
    CheckHitIndex(hit);
    hit_beta[hit] = value;
  }

  void CopyHit(Int_t this_plane, Int_t copy_plane, THcGoodLADHit *copyhit) {
    if (copy_plane == 0) {
      SetPlane(this_plane, copyhit->GetPlaneHit0());
      SetPaddle(this_plane, copyhit->GetPaddleHit0());
      SetTrackID(this_plane, copyhit->GetTrackIDHit0());
      SetDeltaPosTrans(this_plane, copyhit->GetDeltaPosTransHit0());
      SetDeltaPosLong(this_plane, copyhit->GetDeltaPosLongHit0());
      SetHitTime(this_plane, copyhit->GetHitTimeHit0());
      SetHitTheta(this_plane, copyhit->GetHitThetaHit0());
      SetHitPhi(this_plane, copyhit->GetHitPhiHit0());
      SetHitEdep(this_plane, copyhit->GetHitEdepHit0());
      SetHitAlpha(this_plane, copyhit->GetHitAlphaHit0());
      SetHitYPos(this_plane, copyhit->GetHitYPosHit0());
      SetHitTOF(this_plane, copyhit->GetHitTOFHit0());
    } else if (copy_plane == 1) {
      SetPlane(this_plane, copyhit->GetPlaneHit1());
      SetPaddle(this_plane, copyhit->GetPaddleHit1());
      SetTrackID(this_plane, copyhit->GetTrackIDHit1());
      SetDeltaPosTrans(this_plane, copyhit->GetDeltaPosTransHit1());
      SetDeltaPosLong(this_plane, copyhit->GetDeltaPosLongHit1());
      SetHitTime(this_plane, copyhit->GetHitTimeHit1());
      SetHitTheta(this_plane, copyhit->GetHitThetaHit1());
      SetHitPhi(this_plane, copyhit->GetHitPhiHit1());
      SetHitEdep(this_plane, copyhit->GetHitEdepHit1());
      SetHitAlpha(this_plane, copyhit->GetHitAlphaHit1());
      SetHitYPos(this_plane, copyhit->GetHitYPosHit1());
      SetHitTOF(this_plane, copyhit->GetHitTOFHit1());
    } else {
      throw std::out_of_range("Invalid copy_plane index. Must be 0 or 1.");
    }
  }

  // Declaring two different methods for each hit is dumb, but RDefVars can't take TClonesArray
  //  objects calling a method with a parameter (or I haven't been able to figure it out)
  //  so we have to do this
  Int_t GetPlaneHit0() const { return plane[0]; }
  Int_t GetPlaneHit1() const { return plane[1]; }

  Int_t GetPaddleHit0() const { return paddle[0]; }
  Int_t GetPaddleHit1() const { return paddle[1]; }

  Int_t GetTrackIDHit0() const { return track_id[0]; }
  Int_t GetTrackIDHit1() const { return track_id[1]; }

  Double_t GetBetaHit0() const { return hit_beta[0]; }
  Double_t GetBetaHit1() const { return hit_beta[1]; }

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

  Double_t GetHitAlphaHit0() const { return hit_alpha[0]; }
  Double_t GetHitAlphaHit1() const { return hit_alpha[1]; }

  Double_t GetHitYPosHit0() const { return hit_yPos[0]; }
  Double_t GetHitYPosHit1() const { return hit_yPos[1]; }

  Double_t GetHitTOFHit0() const { return hit_tof[0]; }
  Double_t GetHitTOFHit1() const { return hit_tof[1]; }

protected:
  Int_t plane[2];
  Int_t paddle[2];
  Int_t track_id[2];
  Double_t delta_pos_trans[2];
  Double_t delta_pos_long[2];
  Double_t hit_time[2];
  Double_t hit_beta[2];
  Double_t hit_theta[2];
  Double_t hit_phi[2];
  Double_t hit_edep[2];
  Double_t hit_tof[2];
  Double_t hit_yPos[2];
  Double_t hit_alpha[2];

private:
  THcGoodLADHit(const THcGoodLADHit &);            // Prevent copy constructor
  THcGoodLADHit &operator=(const THcGoodLADHit &); // Prevent assignment operator
  ClassDef(THcGoodLADHit, 0)                       // Class for good hodo hits
};

#endif
//_____________________________________________________________________________