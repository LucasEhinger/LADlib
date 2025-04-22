
#include "THcLADKine.h"
#include "THcGlobals.h"
#include "THcHallCSpectrometer.h"
#include "THcLADGEM.h"
#include "THcLADHodoscope.h"
#include "THcParmList.h"
#include "THcPrimaryKine.h"
#include "THcReactionPoint.h"
#include "TVector3.h"
#include "VarDef.h"
#include "VarType.h"

ClassImp(THcLADKine)
    //_____________________________________________________________________________
    THcLADKine::THcLADKine(const char *name, const char *description, const char *spectro, const char *primary_kine,
                           const char *vertex_module)
    : THcPrimaryKine(name, description, spectro, primary_kine, 0.938), fSpecName(spectro), fGEM(nullptr),
      fHodoscope(nullptr), fVertexModuleName(vertex_module), fVertexModule(nullptr) {
  // Constructor
  fGoodLADHits = new TClonesArray("THcLADHodoHit", MAXGOODHITS);
  goodhit_n    = 0;
  fFixed_z     = nullptr;
}
//_____________________________________________________________________________
THcLADKine::~THcLADKine() {
  // Destructor
  DefineVariables(kDelete);
  delete fGoodLADHits;
  delete[] fFixed_z;
}
//_____________________________________________________________________________
void THcLADKine::Clear(Option_t *opt) {
  // Clear the object
  THcPrimaryKine::Clear(opt);
  fGoodLADHits->Clear();
  goodhit_n = 0;
}
//_____________________________________________________________________________
THaAnalysisObject::EStatus THcLADKine::Init(const TDatime &run_time) {

  // Initialize the object
  if (fSpectro == nullptr) {
    fSpectro = dynamic_cast<THcHallCSpectrometer *>(FindModule(fSpecName.Data(), "THcHallCSpectrometer"));
    if (fSpectro == nullptr) {
      Error("Init", "Cannot find spectrometer %s", fSpecName.Data());
      return kInitError;
    }
  }

  if (fHodoscope == nullptr) {
    if (fSpectro) {
      fHodoscope = dynamic_cast<THcLADHodoscope *>(fSpectro->GetDetector("ladhod"));
    }
    if (fHodoscope == nullptr) {
      Error("Init", "No Hodoscope module found in spectrometer");
      return kInitError;
    }
  }
  if (fGEM == nullptr) {
    if (fSpectro) {
      fGEM = dynamic_cast<THcLADGEM *>(fSpectro->GetDetector("gem"));
    }
    if (fGEM == nullptr) {
      Error("Init", "No GEM module found in spectrometer");
      return kInitError;
    }
  }

  if (fVertexModuleName != "") {
    fVertexModule = static_cast<THcReactionPoint *>(FindModule(fVertexModuleName, "THcReactionPoint"));
    if (fVertexModule == nullptr) {
      Error("Init", "No Vertex module found");
      return kInitError;
    }
  }

  return THaPhysicsModule::Init(run_time);
}
//_____________________________________________________________________________
Int_t THcLADKine::ReadDatabase(const TDatime &date) {
  // Read the database
  Int_t err = THcPrimaryKine::ReadDatabase(date);
  if (err)
    return err;

  char prefix[2];
  prefix[0] = 'l';
  prefix[1] = '\0';

  fD0Cut_wVertex      = 0.0;
  fD0Cut_noVertex     = 0.0;
  fmax_dTrans_match   = 0.0;
  fmax_dLong_match    = 0.0;
  fNfixed_z           = 0;
  fglobal_time_offset = 0.0;

  DBRequest list[] = {{"d0_cut_wVertex", &fD0Cut_wVertex, kDouble, 0, 1},
                      {"d0_cut_noVertex", &fD0Cut_noVertex, kDouble, 0, 1},
                      {"max_dTrans_hitMatch", &fmax_dTrans_match, kDouble, 0, 1},
                      {"max_dLong_hitMatch", &fmax_dLong_match, kDouble, 0, 1},
                      {"nfixed_z", &fNfixed_z, kInt, 0, 1},
                      {"global_time_offset", &fglobal_time_offset, kDouble, 0, 1},
                      {0}};
  gHcParms->LoadParmValues((DBRequest *)&list, prefix);

  delete[] fFixed_z;
  if (fNfixed_z > 0) {
    fFixed_z          = new Double_t[fNfixed_z];
    DBRequest list2[] = {{"fixed_z_pos", fFixed_z, kDouble, fNfixed_z}, {0}};
    gHcParms->LoadParmValues((DBRequest *)&list2, prefix);
  } else {
    fFixed_z = nullptr;
  }

  return err;
}
//_____________________________________________________________________________
Int_t THcLADKine::Process(const THaEvData &evdata) {

  // Check track projection to vertex
  fGEMTracks    = fGEM->GetTracks();
  Int_t ntracks = fGEMTracks->GetLast() + 1;
  std::vector<bool> isGoodTrack(ntracks, false);
  for (Int_t i = 0; i < ntracks; i++) {
    THcLADGEMTrack *track = static_cast<THcLADGEMTrack *>(fGEMTracks->At(i));
    if (track == nullptr)
      continue;

    // Get the hit positions
    TVector3 v_hit1(track->GetX1(), track->GetY1(), track->GetZ1());
    TVector3 v_hit2(track->GetX2(), track->GetY2(), track->GetZ2());

    TVector3 vertex;
    if (fVertexModule->HasVertex()) {
      vertex = fVertexModule->GetVertex();
      track->SetGoodD0(kTRUE);
    } else {
      vertex.SetXYZ(0, 0, 0);
    }

    // Fix track vertex (for improved resolution on multifoils)
    if (fNfixed_z > 0 && track->GetGoodD0()) {
      std::vector<double> distances;
      for (Int_t j = 0; j < fNfixed_z; j++) {
        distances.push_back(abs(vertex.Z() - fFixed_z[j]));
      }
      auto min_it     = std::min_element(distances.begin(), distances.end());
      Int_t min_index = std::distance(distances.begin(), min_it);
      vertex.SetZ(fFixed_z[min_index]);
    }

    // Calculate d0
    double numer = ((vertex - v_hit1).Cross((vertex - v_hit2))).Mag();
    double denom = (v_hit2 - v_hit1).Mag();
    // here we can put a range/fiducial cut on d0 taking into account the target size
    double d0 = numer / denom;
    track->SetD0(d0);

    if (track->GetGoodD0()) {
      if (d0 > fD0Cut_wVertex) {
        track->SetGoodD0(kFALSE);
      }
    } else {
      if (d0 > fD0Cut_noVertex) {
        track->SetGoodD0(kFALSE);
      }
    }
  }

  //////////////////////////////////////////////////////////////////////////////
  // Remove hits with bad tracks and duplicate hits
  TClonesArray *LADHits_unfiltered = fHodoscope->GetLADGoodHits();
  Int_t nHits                      = LADHits_unfiltered->GetLast() + 1;
  for (Int_t i = 0; i < nHits; i++) {
    THcGoodLADHit *hit = static_cast<THcGoodLADHit *>(LADHits_unfiltered->At(i));
    if (hit == nullptr)
      continue;

    // Check if the hit has a track that points back to the origin.
    Int_t track_id = hit->GetTrackIDHit0();
    if (track_id < 0 || track_id >= ntracks || !isGoodTrack[track_id])
      continue;

    // Get the track parameters
    Int_t plane  = hit->GetPlaneHit0();
    Int_t paddle = hit->GetPaddleHit0();

    Int_t plane_loc; // Front vs back plane
    if (plane == 0 || plane == 2 || plane == 4)
      plane_loc = 0;
    else
      plane_loc = 1;

    Int_t matching_hit_index = -1;
    Int_t partner_hit_index  = -1;
    for (Int_t j = 0; j < goodhit_n; j++) {
      THcGoodLADHit *goodhit = static_cast<THcGoodLADHit *>(fGoodLADHits->At(j));
      if (goodhit == nullptr)
        continue;

      // Check if the hit is already in the list
      Int_t good_hit_plane  = plane_loc ? goodhit->GetPlaneHit0() : goodhit->GetPlaneHit1();
      Int_t good_hit_paddle = plane_loc ? goodhit->GetPaddleHit0() : goodhit->GetPaddleHit1();
      if (good_hit_plane == plane && good_hit_paddle == paddle) {
        matching_hit_index = j;
        break;
      }
      // Check if the hit has a partner hit
      Int_t partner_hit_track_id = plane_loc ? goodhit->GetTrackIDHit1() : goodhit->GetTrackIDHit0();
      if (partner_hit_track_id == track_id) {
        partner_hit_index = j;
      }
    }
    if (matching_hit_index == -1 && partner_hit_index == -1) {
      THcGoodLADHit *goodhit = new ((*fGoodLADHits)[goodhit_n]) THcGoodLADHit();
      goodhit_n++;
      goodhit->SetPlane(plane_loc, plane);
      goodhit->SetPaddle(plane_loc, paddle);
      goodhit->SetTrackID(plane_loc, track_id);
      goodhit->SetBeta(plane_loc, hit->GetBetaHit0());
      goodhit->SetDeltaPosTrans(plane_loc, hit->GetDeltaPosTransHit0());
      goodhit->SetDeltaPosLong(plane_loc, hit->GetDeltaPosLongHit0());
      goodhit->SetHitTime(plane_loc, hit->GetHitTimeHit0());
      goodhit->SetHitTheta(plane_loc, hit->GetHitThetaHit0());
      goodhit->SetHitPhi(plane_loc, hit->GetHitPhiHit0());
      goodhit->SetHitEdep(plane_loc, hit->GetHitEdepHit0());
      // Could probably overload an operator to do this, but that's too much work...
    }
    if (matching_hit_index != -1) {
      THcGoodLADHit *goodhit = static_cast<THcGoodLADHit *>(fGoodLADHits->At(matching_hit_index));
      if (abs(hit->GetDeltaPosTransHit0()) <
          (plane_loc ? abs(goodhit->GetDeltaPosTransHit0()) : abs(goodhit->GetDeltaPosTransHit1()))) {
        goodhit->SetTrackID(plane_loc, track_id);
        goodhit->SetBeta(plane_loc, hit->GetBetaHit0());
        goodhit->SetDeltaPosTrans(plane_loc, hit->GetDeltaPosTransHit0());
        goodhit->SetDeltaPosLong(plane_loc, hit->GetDeltaPosLongHit0());
        goodhit->SetHitTime(plane_loc, hit->GetHitTimeHit0());
        goodhit->SetHitTheta(plane_loc, hit->GetHitThetaHit0());
        goodhit->SetHitPhi(plane_loc, hit->GetHitPhiHit0());
        goodhit->SetHitEdep(plane_loc, hit->GetHitEdepHit0());
      }
    }
    // Todo: Calculate beta, alpha, etc. for the hit
    // FIXME. One track can still have multiple (distinct) hodo hits that are counted as good
  }

  //////////////////////////////////////////////////////////////////////////////
  // Find matching front + back plane hits
  // Create a vector of pairs to store matching hits for each track
  std::vector<std::pair<Int_t, Int_t>> matchingHits(ntracks, {-1, -1});
  std::vector<double> deltaPosTransValues;

  //TODO: This can probably be incorporated into the loop above
  for (Int_t i = 0; i < goodhit_n; i++) {
    THcGoodLADHit *goodhit = static_cast<THcGoodLADHit *>(fGoodLADHits->At(i));
    if (goodhit == nullptr)
      continue;

    // Add delta_pos_trans values for both planes (if valid) to the vector
    if (goodhit->GetDeltaPosTransHit0() != 0) {
      deltaPosTransValues.push_back(goodhit->GetDeltaPosTransHit0());
    }
    if (goodhit->GetDeltaPosTransHit1() != 0) {
      deltaPosTransValues.push_back(goodhit->GetDeltaPosTransHit1());
    }

    Int_t trackID0 = goodhit->GetTrackIDHit0();
    Int_t trackID1 = goodhit->GetTrackIDHit1();

    if (trackID0 >= 0 && trackID0 < ntracks) {
      // Check if the track ID is already in the vector
      if (matchingHits[trackID0].first == -1) {
        // If it is, set the second hit index to the current index
        matchingHits[trackID0].first = i;
      } else if (deltaPosTransValues[matchingHits[trackID0].first] > goodhit->GetDeltaPosTransHit0()) {
        matchingHits[trackID0].first = i;
      }
    }
    if (trackID1 >= 0 && trackID1 < ntracks) {
      // Check if the track ID is already in the vector
      if (matchingHits[trackID1].second == -1) {
        // If it is, set the second hit index to the current index
        matchingHits[trackID1].second = i;
      } else if (deltaPosTransValues[matchingHits[trackID1].second] > goodhit->GetDeltaPosTransHit1()) {
        matchingHits[trackID1].second = i;
      }
    }
  }

  // Loop over the matching hits and set the partner hit index
  for (auto &match : matchingHits) {
    if (match.first != -1 && match.second != -1) {
      THcGoodLADHit *firstHit  = static_cast<THcGoodLADHit *>(fGoodLADHits->At(match.first));
      THcGoodLADHit *secondHit = static_cast<THcGoodLADHit *>(fGoodLADHits->At(match.second));
      if (firstHit && secondHit) {
        // Fill the second hit into the first for plane 1
        firstHit->SetPlane(1, secondHit->GetPlaneHit1());
        firstHit->SetPaddle(1, secondHit->GetPaddleHit1());
        firstHit->SetTrackID(1, secondHit->GetTrackIDHit1());
        firstHit->SetBeta(1, secondHit->GetBetaHit1());
        firstHit->SetDeltaPosTrans(1, secondHit->GetDeltaPosTransHit1());
        firstHit->SetDeltaPosLong(1, secondHit->GetDeltaPosLongHit1());
        firstHit->SetHitTime(1, secondHit->GetHitTimeHit1());
        firstHit->SetHitTheta(1, secondHit->GetHitThetaHit1());
        firstHit->SetHitPhi(1, secondHit->GetHitPhiHit1());
        firstHit->SetHitEdep(1, secondHit->GetHitEdepHit1());

        // Remove the second hit from the good hits array
        fGoodLADHits->RemoveAt(match.second);
        goodhit_n--;
      }
    }
  }
  fGoodLADHits->Compress(); // Compress the array to remove null entries

  // LADHits_unfiltered->Clear(); // Clear the unfiltered hits
  // fGEMTracks->Clear();      // Clear the tracks
  return kOK;
}
//_____________________________________________________________________________
Int_t THcLADKine::DefineVariables(EMode mode) {
  if (mode == kDefine && fIsSetup)
    return kOK;
  fIsSetup = (mode == kDefine);

  // Define variables for the analysis tree
  if (mode == kDefine) {
    RVarDef vars[] = {
        {"n_goodhits", "Number of good hits", "goodhit_n"},
        {"goodhit_plane_0", "Good hit plane", "fGoodLADHits.THcGoodLADHit.GetPlaneHit0()"},
        {"goodhit_paddle_0", "Good hit paddle", "fGoodLADHits.THcGoodLADHit.GetPaddleHit0()"},
        {"goodhit_trackid_0", "Good hit track ID", "fGoodLADHits.THcGoodLADHit.GetTrackIDHit0()"},
        {"goodhit_beta_0", "Good hit beta", "fGoodLADHits.THcGoodLADHit.GetBetaHit0()"},
        {"goodhit_deltapostrans_0", "Good hit delta pos trans", "fGoodLADHits.THcGoodLADHit.GetDeltaPosTransHit0()"},
        {"goodhit_deltaposlong_0", "Good hit delta pos long", "fGoodLADHits.THcGoodLADHit.GetDeltaPosLongHit0()"},
        {"goodhit_hittime_0", "Good hit time", "fGoodLADHits.THcGoodLADHit.GetHitTimeHit0()"},
        {"goodhit_hittheta_0", "Good hit theta", "fGoodLADHits.THcGoodLADHit.GetHitThetaHit0()"},
        {"goodhit_hitphi_0", "Good hit phi", "fGoodLADHits.THcGoodLADHit.GetHitPhiHit0()"},
        {"goodhit_hitedep_0", "Good hit energy deposition", "fGoodLADHits.THcGoodLADHit.GetHitEdepHit0()"},
        {"goodhit_plane_1", "Good hit plane (second plane)", "fGoodLADHits.THcGoodLADHit.GetPlaneHit1()"},
        {"goodhit_paddle_1", "Good hit paddle (second plane)", "fGoodLADHits.THcGoodLADHit.GetPaddleHit1()"},
        {"goodhit_trackid_1", "Good hit track ID (second plane)", "fGoodLADHits.THcGoodLADHit.GetTrackIDHit1()"},
        {"goodhit_beta_1", "Good hit beta (second plane)", "fGoodLADHits.THcGoodLADHit.GetBetaHit1()"},
        {"goodhit_deltapostrans_1", "Good hit delta pos trans (second plane)",
         "fGoodLADHits.THcGoodLADHit.GetDeltaPosTransHit1()"},
        {"goodhit_deltaposlong_1", "Good hit delta pos long (second plane)",
         "fGoodLADHits.THcGoodLADHit.GetDeltaPosLongHit1()"},
        {"goodhit_hittime_1", "Good hit time (second plane)", "fGoodLADHits.THcGoodLADHit.GetHitTimeHit1()"},
        {"goodhit_hittheta_1", "Good hit theta (second plane)", "fGoodLADHits.THcGoodLADHit.GetHitThetaHit1()"},
        {"goodhit_hitphi_1", "Good hit phi (second plane)", "fGoodLADHits.THcGoodLADHit.GetHitPhiHit1()"},
        {"goodhit_hitedep_1", "Good hit energy deposition (second plane)",
         "fGoodLADHits.THcGoodLADHit.GetHitEdepHit1()"},
        {0}};
    return DefineVarsFromList(vars, mode);
  }
  return kOK;
}
//_____________________________________________________________________________
