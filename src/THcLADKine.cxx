
#include "THcLADKine.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"
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
#include <algorithm>
#include <cmath>
#include <limits>

ClassImp(THcLADKine)
    //_____________________________________________________________________________
    THcLADKine::THcLADKine(const char *name, const char *description, const char *spectro, const char *primary_kine,
                           const char *vertex_module)
    : THcPrimaryKine(name, description, spectro, primary_kine, 0.938), fSpecName(spectro), fGEM(nullptr),
      fHodoscope(nullptr), fVertexModuleName(vertex_module), fVertexModule(nullptr), fTrack(nullptr),
      fMinimizer(nullptr) {
  // Constructor
  // fGoodLADHits    = nullptr;
  goodhit_n       = 0;
  fTVertex        = kBig;
  fRFTime         = kBig;
  fTVertex_RFcorr = kBig;
  fFixed_z        = nullptr;
  rf_offset       = nullptr;
  n_rf_offsets    = 0;
  rf_period       = 4.00801; // Default RF period in ns

  fFrontPlaneEdepCut = 100.0; // MeV
  fBackPlaneDtMin    = 2.8;   // ns
  fBackPlaneDtMax    = 10.0;  // ns
  fBackPlaneYDiag    = 35.0;  // MeV
  fBackPlaneXDiag    = 5.0;   // ns

  fZCellMin     = -15.0;                     // cm
  fZCellMax     = +15.0;                     // cm
  fThetaMin     = 60.0 * TMath::DegToRad();  // rad
  fThetaMax     = 170.0 * TMath::DegToRad(); // rad
  fPhiMin       = -50.0 * TMath::DegToRad(); // rad
  fPhiMax       = +50.0 * TMath::DegToRad(); // rad
  fchisq_cut[0] = 20.0;                      // default chi2 cut for tracks with 2 hodo hits
  fchisq_cut[1] = 15.0;                      // default chi2 cut for tracks with 1 hodo hit
  fSigma_GEM    = 0.1;                       // default GEM resolution in cm
  fSigma_Hodo   = 10;                        // default Hodoscope resolution in cm

  fDoNoVertexTracking    = 1; // run the no-vertex tracking pass by default
  fDoXZTracking          = 0; // x-z (no-y) tracking off by default
  fDo1DClusterTracking   = 0; // 1D-cluster tracking off by default
  fMax1DClusterPerSlot   = 3;    // max clusters per (layer,axis) slot
  fMax1DCandidates       = 2000; // max projective fits attempted per event
  fSigma_Hodo_x          = 6.0;  // cm, x-z / paddle-direction resolution
  fSigma_Hodo_y          = 3.0;  // cm, along-paddle (y) resolution
  fHodoXZZeroHW          = 0.0;  // cm, x-z hodo residual zeroing half-width (off)
}
//_____________________________________________________________________________
THcLADKine::~THcLADKine() {
  // Destructor
  DefineVariables(kDelete);
  // fGoodLADHits->Delete();
  // delete fGoodLADHits;
  delete[] fFixed_z;
  delete[] rf_offset;
  delete fMinimizer;
}
//_____________________________________________________________________________
void THcLADKine::Clear(Option_t *opt) {
  // Clear the object
  THcPrimaryKine::Clear(opt);
  // fGoodLADHits->Delete();
  goodhit_n = 0;
  // NB: rf_offset is the DB-loaded offset table (3*n_rf_offsets entries) and must
  // NOT be zeroed per-event here -- doing so corrupted the run lookup.
  fTVertex        = kBig;
  fRFTime         = kBig;
  fTVertex_RFcorr = kBig;
  vertex.SetXYZ(0, 0, 0);

  // 1D-cluster projective tracking writes its results directly onto the
  // THcGoodLADHit objects (reset in their own Clear), so nothing to clear here.
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

  TString fTrigDetName;
  // Get the apparatus prefix
  char aparatus_prefix[2];
  aparatus_prefix[0] = tolower(fSpecName[0]);
  aparatus_prefix[1] = '\0';
  if (aparatus_prefix[0] == 'h') {
    fTrigDetName = "T.hms";
  } else if (aparatus_prefix[0] == 'p') {
    fTrigDetName = "T.shms";
  } else {
    Error("Init", "Invalid apparatus prefix: %s", aparatus_prefix);
    return kInitError;
  }

  fTrigDet = dynamic_cast<THcTrigDet *>(FindModule(fTrigDetName.Data(), "THcTrigDet"));
  if (!fTrigDet) {
    cout << "THcCoinTime module  Cannnot find TrigDet = " << fTrigDetName.Data() << endl;
    fStatus = kInitError;
    return fStatus;
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

  fD0Cut_wVertex       = 0.0;
  fD0Cut_noVertex      = 0.0;
  fNfixed_z            = 0;
  fglobal_time_offset  = 0.0;
  fTrk_dtCut           = 10.0;
  fDoNoVertexTracking  = 1;    // default: no-vertex tracking on
  fDoXZTracking        = 0;    // default: x-z (no-y) tracking off
  fDo1DClusterTracking = 0;    // default: 1D-cluster tracking off
  fMax1DClusterPerSlot = 3;    // default: keep 3 clusters per (layer,axis) slot
  fMax1DCandidates     = 2000; // default: cap projective fits per event
  fSigma_Hodo_x        = 6.0;  // default: x-z / paddle-direction resolution (cm)
  fSigma_Hodo_y        = 3.0;  // default: along-paddle (y) resolution (cm)
  fHodoXZZeroHW        = 0.0;  // default: no x-z hodo residual zeroing

  cout << "Reading LAD Kinematics parameters from database..." << endl;

  DBRequest list[] = {{"d0_cut_wVertex", &fD0Cut_wVertex, kDouble, 0, 1},
                      {"d0_cut_noVertex", &fD0Cut_noVertex, kDouble, 0, 1},
                      {"nfixed_z", &fNfixed_z, kInt, 0, 1},
                      {"global_time_offset", &fglobal_time_offset, kDouble, 0, 1},
                      {"trk_dt_cut", &fTrk_dtCut, kDouble, 0, 1},
                      {"do_noVertex_tracking", &fDoNoVertexTracking, kInt, 0, 1},
                      {"do_xz_tracking", &fDoXZTracking, kInt, 0, 1},
                      {"do_1Dcluster_tracking", &fDo1DClusterTracking, kInt, 0, 1},
                      {"max_1Dcluster_per_slot", &fMax1DClusterPerSlot, kInt, 0, 1},
                      {"max_1Dcluster_candidates", &fMax1DCandidates, kInt, 0, 1},
                      {"sigma_hodo_x", &fSigma_Hodo_x, kDouble, 0, 1},
                      {"sigma_hodo_y", &fSigma_Hodo_y, kDouble, 0, 1},
                      {"hodo_xz_zero_hw", &fHodoXZZeroHW, kDouble, 0, 1},
                      {"_rf_period", &rf_period, kDouble, 0, 1},
                      {"proton_front_edep_cut", &fFrontPlaneEdepCut, kDouble, 0, 1},
                      {"proton_back_dt_min", &fBackPlaneDtMin, kDouble, 0, 1},
                      {"proton_back_dt_max", &fBackPlaneDtMax, kDouble, 0, 1},
                      {"proton_back_y_diag", &fBackPlaneYDiag, kDouble, 0, 1},
                      {"proton_back_x_diag", &fBackPlaneXDiag, kDouble, 0, 1},
                      {0}};
  gHcParms->LoadParmValues((DBRequest *)&list, prefix);

  DBRequest list_rf[] = {{"_nr_rf_offsets", &n_rf_offsets, kInt, 0, 1}, {0}};
  gHcParms->LoadParmValues((DBRequest *)&list_rf, prefix);
  rf_offset = new Double_t[n_rf_offsets * 3];
  for (Int_t i = 0; i < n_rf_offsets * 3; i++) {
    rf_offset[i] = 0.0;
  }
  DBRequest list_rf_offsets[] = {{"_rf_offset", &rf_offset[0], kDouble, (UInt_t)(n_rf_offsets * 3)}, {0}};
  gHcParms->LoadParmValues((DBRequest *)&list_rf_offsets, prefix);

  delete[] fFixed_z;
  if (fNfixed_z > 0) {
    fFixed_z          = new Double_t[fNfixed_z];
    DBRequest list2[] = {{"fixed_z_pos", fFixed_z, kDouble, static_cast<UInt_t>(fNfixed_z)}, {0}};
    gHcParms->LoadParmValues((DBRequest *)&list2, prefix);
  } else {
    fFixed_z = nullptr;
  }

  // NB: LoadParmValues prepends the prefix ("l"), so these keys are read as
  // "lz_cell_min", "ltheta_min", ... in the database (matching the other keys
  // above). The previous "ladkin." names resolved to "lladkin.*" and never matched.
  DBRequest track_constraints[] = {
      {"z_cell_min", &fZCellMin, kDouble, 0, 1},   {"z_cell_max", &fZCellMax, kDouble, 0, 1},
      {"theta_min", &fThetaMin, kDouble, 0, 1},    {"theta_max", &fThetaMax, kDouble, 0, 1},
      {"phi_min", &fPhiMin, kDouble, 0, 1},        {"phi_max", &fPhiMax, kDouble, 0, 1},
      {"chisq_cut", fchisq_cut, kDouble, 2, 1},    {"sigma_gem", &fSigma_GEM, kDouble, 0, 1},
      {"sigma_hodo", &fSigma_Hodo, kDouble, 0, 1}, {0}};
  gHcParms->LoadParmValues((DBRequest *)&track_constraints, prefix);

  if (fZCellMin > fZCellMax) {
    Double_t tmp = fZCellMin;
    fZCellMin    = fZCellMax;
    fZCellMax    = tmp;
  }
  if (fThetaMin > fThetaMax) {
    Double_t tmp = fThetaMin;
    fThetaMin    = fThetaMax;
    fThetaMax    = tmp;
  }
  if (fPhiMin > fPhiMax) {
    Double_t tmp = fPhiMin;
    fPhiMin      = fPhiMax;
    fPhiMax      = tmp;
  }

  fThetaMin = TMath::Max(0.0, TMath::Min(TMath::Pi(), fThetaMin));
  fThetaMax = TMath::Max(0.0, TMath::Min(TMath::Pi(), fThetaMax));
  fPhiMin   = TMath::Max(-TMath::Pi(), TMath::Min(TMath::Pi(), fPhiMin));
  fPhiMax   = TMath::Max(-TMath::Pi(), TMath::Min(TMath::Pi(), fPhiMax));

  return err;
}
//_____________________________________________________________________________
Int_t THcLADKine::Process(const THaEvData &evdata) {
  // Get Golden Track (for tof calculation later)
  fTrack = fSpectro->GetGoldenTrack();
  CalculateTVertex();

  // Get LADGoodHits for track matching
  TClonesArray *LADHits_unfiltered = fHodoscope->GetLADGoodHits();
  Int_t nHits                      = LADHits_unfiltered->GetLast() + 1;
  goodhit_n                        = 0;
  //////////////////////////////////////////////////////////////////////////////
  // Check track projection to vertex
  fGEMTracks    = fGEM->GetTracks();
  Int_t ntracks = fGEMTracks->GetLast() + 1;
  std::vector<bool> isGoodTrack(ntracks, false);
  for (Int_t iTrack = 0; iTrack < ntracks; iTrack++) {
    THcLADGEMTrack *track = static_cast<THcLADGEMTrack *>(fGEMTracks->At(iTrack));
    if (track == nullptr)
      continue;

    // Get the hit positions
    TVector3 v_hit1(track->GetX1(), track->GetY1(), track->GetZ1());
    TVector3 v_hit2(track->GetX2(), track->GetY2(), track->GetZ2());

    if (fVertexModule && fVertexModule->HasVertex()) {
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

    Double_t gemdir[3];
    gemdir[0] = TMath::ACos((v_hit2.Z() - v_hit1.Z()) / (v_hit2 - v_hit1).Mag());
    gemdir[1] = TMath::ATan2((v_hit2.Y() - v_hit1.Y()), (v_hit2.X() - v_hit1.X()));
    gemdir[2] = vertex.Z(); // Hopefully the GEM track will be the same as elecctron vertex Z

    Double_t bestchisq[3][3] = {{kBig, kBig, kBig},
                                {kBig, kBig, kBig},
                                {kBig, kBig, kBig}}; // no hodo hits, 1 hodo hit, 2 hodo hits+ (F-only, B-only, F&B)
    bool usedHodoHit[3][3]   = {{kFALSE, kFALSE, kFALSE}, {kFALSE, kFALSE, kFALSE}, {kFALSE, kFALSE, kFALSE}};
    Double_t best_gemdir[3][3][3] = {}; // zero-initialized so an unselected entry can never be read uninitialized
    int bestHodoHitIndex[3][3]    = {{-1, -1, -1}, {-1, -1, -1}, {-1, -1, -1}};
    // fit GEM only first, if chisq is negative, skip the track and mark as bad
    Double_t dir[3];
    dir[0] = gemdir[0];
    dir[1] = gemdir[1];
    dir[2] = gemdir[2]; // if the track doesn't have a good hodoscope match, just fit to the GEM hits and vertex
    bestchisq[0][0]   = FitTrack(vertex, {v_hit1, v_hit2}, {fSigma_GEM, fSigma_GEM}, dir);
    usedHodoHit[0][0] = kTRUE;
    if (bestchisq[0][0] < 0) {
      track->SetIsGoodTrack(kFALSE);
      track->SetChisq(bestchisq[0][0]);
      track->SetAngles(-kBig, -kBig);
      track->SetProjVertex(-kBig, -kBig, -kBig);
      track->SetD0(-kBig);
      track->SetHasHodoHit(kFALSE);
      track->SetBestHodoHit(nullptr);
      continue;
    }
    best_gemdir[0][0][0]         = dir[0];
    best_gemdir[0][0][1]         = dir[1];
    best_gemdir[0][0][2]         = dir[2];
    int hit_index_for_best_chisq = -1; // (0-2)*3 for (no hodo hits, 1 hodo hit, 2 hodo hits) + (F-only, B-only))
    // loop over hodoscope hits to find the best match to the track projection
    for (int iHodoHit = 0; iHodoHit < nHits; iHodoHit++) {
      THcGoodLADHit *goodHit = static_cast<THcGoodLADHit *>(LADHits_unfiltered->At(iHodoHit));
      if (goodHit == nullptr)
        continue;
      // Check the number hits in goohit
      int num_hodo_hits = 0;
      if (goodHit->GetPlaneHit0() >= 0 && goodHit->GetPlaneHit0() < 999) {
        num_hodo_hits++;
      }
      if (goodHit->GetPlaneHit1() >= 0 && goodHit->GetPlaneHit1() < 999) {
        num_hodo_hits++;
      }

      std::vector<TVector3> v_hits;
      std::vector<Double_t> v_resolutions;
      v_hits.push_back(v_hit1);
      v_hits.push_back(v_hit2);
      v_resolutions.push_back(fSigma_GEM);
      v_resolutions.push_back(fSigma_GEM);

      if (num_hodo_hits == 1) {
        if (goodHit->GetPlaneHit0() >= 0 && goodHit->GetPlaneHit0() < 999) {
          TVector3 hodo_hit_pos = fHodoscope->GetHitPositionLab(goodHit->GetPlaneHit0(), goodHit->GetPaddleHit0(),
                                                                goodHit->GetHitYPosHit0());
          v_hits.push_back(hodo_hit_pos);
          v_resolutions.push_back(fSigma_Hodo);
        } else if (goodHit->GetPlaneHit1() >= 0 && goodHit->GetPlaneHit1() < 999) {
          TVector3 hodo_hit_pos = fHodoscope->GetHitPositionLab(goodHit->GetPlaneHit1(), goodHit->GetPaddleHit1(),
                                                                goodHit->GetHitYPosHit1());
          v_hits.push_back(hodo_hit_pos);
          v_resolutions.push_back(fSigma_Hodo);
        }
      } else if (num_hodo_hits == 2) {
        TVector3 hodo_hit_pos1 =
            fHodoscope->GetHitPositionLab(goodHit->GetPlaneHit0(), goodHit->GetPaddleHit0(), goodHit->GetHitYPosHit0());
        TVector3 hodo_hit_pos2 =
            fHodoscope->GetHitPositionLab(goodHit->GetPlaneHit1(), goodHit->GetPaddleHit1(), goodHit->GetHitYPosHit1());
        v_hits.push_back(hodo_hit_pos1);
        v_hits.push_back(hodo_hit_pos2);
        v_resolutions.push_back(fSigma_Hodo);
        v_resolutions.push_back(fSigma_Hodo);
      }
      Double_t dir[3];
      dir[0]         = gemdir[0];
      dir[1]         = gemdir[1];
      dir[2]         = gemdir[2];
      Double_t chisq = FitTrack(vertex, v_hits, v_resolutions, dir);
      if (num_hodo_hits == 1) {
        if (chisq >= 0 && chisq < bestchisq[1][0]) {
          bestchisq[1][0]        = chisq;
          best_gemdir[1][0][0]   = dir[0];
          best_gemdir[1][0][1]   = dir[1];
          best_gemdir[1][0][2]   = dir[2];
          bestHodoHitIndex[1][0] = iHodoHit;
          usedHodoHit[1][0]      = kTRUE;
        }
      } else if (num_hodo_hits == 2) {
        if (chisq >= 0 && chisq < bestchisq[2][0]) {
          bestchisq[2][0]        = chisq;
          best_gemdir[2][0][0]   = dir[0];
          best_gemdir[2][0][1]   = dir[1];
          best_gemdir[2][0][2]   = dir[2];
          bestHodoHitIndex[2][0] = iHodoHit;
          usedHodoHit[2][0]      = kTRUE;
        }
        v_hits.pop_back();
        v_hits.pop_back();
        v_resolutions.pop_back();
        dir[0] = gemdir[0];
        dir[1] = gemdir[1];
        dir[2] = gemdir[2];
        TVector3 hodo_hit_pos1 =
            fHodoscope->GetHitPositionLab(goodHit->GetPlaneHit0(), goodHit->GetPaddleHit0(), goodHit->GetHitYPosHit0());
        TVector3 hodo_hit_pos2 =
            fHodoscope->GetHitPositionLab(goodHit->GetPlaneHit1(), goodHit->GetPaddleHit1(), goodHit->GetHitYPosHit1());
        v_hits.push_back(hodo_hit_pos1);
        Double_t chisq1 = FitTrack(vertex, v_hits, v_resolutions, dir);
        if (chisq1 >= 0 && chisq1 < bestchisq[2][1]) {
          bestchisq[2][1]        = chisq1;
          best_gemdir[2][1][0]   = dir[0];
          best_gemdir[2][1][1]   = dir[1];
          best_gemdir[2][1][2]   = dir[2];
          bestHodoHitIndex[2][1] = iHodoHit;
          usedHodoHit[2][1]      = kTRUE;
        }
        v_hits.pop_back();
        v_hits.push_back(hodo_hit_pos2);
        dir[0]          = gemdir[0];
        dir[1]          = gemdir[1];
        dir[2]          = gemdir[2];
        Double_t chisq2 = FitTrack(vertex, v_hits, v_resolutions, dir);
        if (chisq2 >= 0 && chisq2 < bestchisq[2][2]) {
          bestchisq[2][2]        = chisq2;
          best_gemdir[2][2][0]   = dir[0];
          best_gemdir[2][2][1]   = dir[1];
          best_gemdir[2][2][2]   = dir[2];
          bestHodoHitIndex[2][2] = iHodoHit;
          usedHodoHit[2][2]      = kTRUE;
        }
      }
      v_hits.clear();
      v_resolutions.clear();
    }
    // bestchisq should always be positive?
    // check between f-only, b-only, and single, which has the best chisq
    hit_index_for_best_chisq = 2 * 3 + 0; // default to 2 hodo hit
    if (usedHodoHit[2][1] && bestchisq[2][1] < bestchisq[2][2] && bestchisq[2][1] < bestchisq[1][0]) {
      hit_index_for_best_chisq = 2 * 3 + 1;
    } else if (usedHodoHit[2][2] && bestchisq[2][2] < bestchisq[2][1] && bestchisq[2][2] < bestchisq[1][0]) {
      hit_index_for_best_chisq = 2 * 3 + 2;
    } else if (usedHodoHit[1][0]) {
      hit_index_for_best_chisq = 1 * 3 + 0;
    }
    if (usedHodoHit[2][0] &&
        (bestchisq[2][0] - bestchisq[hit_index_for_best_chisq / 3][hit_index_for_best_chisq % 3] <= fchisq_cut[0])) {
      // check if the 2 hodo hit fit is better than the best single hodo hit fit by more than the chisq cut, if so, use
      // the 2 hodo hit fit
      hit_index_for_best_chisq = 2 * 3 + 0;
    } else if (usedHodoHit[0][0] && (bestchisq[hit_index_for_best_chisq / 3][hit_index_for_best_chisq % 3] -
                                     bestchisq[0][0]) > fchisq_cut[1]) {
      // check if the GEM only fit is better than the best hodo hit fit by more than the chisq cut, if so, use the GEM
      // only fit
      hit_index_for_best_chisq = 0 * 3 + 0;
    }
    track->SetChisq(bestchisq[hit_index_for_best_chisq / 3][hit_index_for_best_chisq % 3]);
    track->SetAngles(best_gemdir[hit_index_for_best_chisq / 3][hit_index_for_best_chisq % 3][0],
                     best_gemdir[hit_index_for_best_chisq / 3][hit_index_for_best_chisq % 3][1]);
    track->SetProjVertex(vertex.X(), vertex.Y(),
                         best_gemdir[hit_index_for_best_chisq / 3][hit_index_for_best_chisq % 3][2]);
    Double_t d0 = fabs(vertex.Z() - best_gemdir[hit_index_for_best_chisq / 3][hit_index_for_best_chisq % 3][2]);
    track->SetD0(d0);

    track->SetHasHodoHit(hit_index_for_best_chisq / 3);
    if (hit_index_for_best_chisq >= 1) {
      THcGoodLADHit *besthit = static_cast<THcGoodLADHit *>(
          LADHits_unfiltered->At(bestHodoHitIndex[hit_index_for_best_chisq / 3][hit_index_for_best_chisq % 3]));
      // THcGoodLADHit *newhit  = static_cast<THcGoodLADHit *>(fGoodLADHits->ConstructedAt(goodhit_n));
      goodhit_n++;
      if (besthit == nullptr) {
        track->SetBestHodoHit(nullptr);
        continue;
      }
      if (goodhit_n <= MAXGOODHITS) {
        if (hit_index_for_best_chisq / 3 == 2) {
          int ntmp = 0;
          if (hit_index_for_best_chisq % 3 == 1 || hit_index_for_best_chisq % 3 == 0) {
            ntmp++;
          }
          if (hit_index_for_best_chisq % 3 == 2 || hit_index_for_best_chisq % 3 == 0) {
            ntmp++;
          }
          track->SetHasHodoHit(ntmp);
          besthit->SetTrackID(track->GetTrackID());
          besthit->SetTrkChiSqr(track->GetChisq());
          // dof = 3*(#GEM + #hodo) - 3 params (full 3D residual, vertex-constrained)
          besthit->SetTrkNdof(3 * (2 + ntmp) - 3);
        } else if (hit_index_for_best_chisq / 3 == 1) {
          besthit->SetTrackID(track->GetTrackID());
          besthit->SetTrkChiSqr(track->GetChisq());
          besthit->SetTrkNdof(3 * (2 + 1) - 3);
          track->SetHasHodoHit(1);
        }
        if (hit_index_for_best_chisq / 3 != 0) {
          track->SetBestHodoHit(besthit);
          isGoodTrack[iTrack] = true;
        }
      }else{
        cout<<"THcLADKine::Process: Exceeded maximum number of good hits ("<<goodhit_n<<"/"<<MAXGOODHITS<<"), skipping hit."<<endl;
      }
    } else {
      track->SetBestHodoHit(nullptr);
      track->SetHasHodoHit(0);
      isGoodTrack[iTrack] = true;
    }

    if (track->GetGoodD0()) {
      if (!(track->GetD0() > 0 && track->GetD0() < fD0Cut_wVertex)) {
        isGoodTrack[iTrack] = false; // reject tracks whose D0 is outside the allowed range
      }
    } else {
      if (!(track->GetD0() > 0 && track->GetD0() < fD0Cut_noVertex)) {
        isGoodTrack[iTrack] = false;
      }
    }
    if (track->GetChisq() < 0) {
      isGoodTrack[iTrack] = false;
    }
    if (track->GetdT() > fTrk_dtCut) {
      isGoodTrack[iTrack] = false;
    }
    track->SetIsGoodTrack(isGoodTrack[iTrack]);
  }

  //////////////////////////////////////////////////////////////////////////////
  // No-vertex tracking
  // Identical in spirit to the loop above, but the fit does NOT use the target
  // vertex. The track is a free 3D line through GEM1, GEM2 and the hodoscope
  // hit(s). Because two points (GEM1, GEM2) already define a line exactly, a
  // no-vertex fit is only meaningful when there is at least one hodoscope hit,
  // so tracks with no hodo hit get no valid no-vertex result.
  for (Int_t iTrack = 0; iTrack < ntracks; iTrack++) {
    // no-vertex tracking pass; enabled via the ldo_noVertex_tracking DB flag (default on).
    // The flag is loop-invariant, so break skips the whole pass (tracks are freshly
    // constructed each event, so their _noTrackVertex members keep the invalid defaults).
    if (!fDoNoVertexTracking)
      break;
    THcLADGEMTrack *track = static_cast<THcLADGEMTrack *>(fGEMTracks->At(iTrack));
    if (track == nullptr)
      continue;

    // Default (invalid) values: overwritten below only if a valid fit is found
    track->SetChisq_noTrackVertex(-kBig);
    track->SetAngles_noTrackVertex(-kBig, -kBig);
    track->SetProjVertex_noTrackVertex(-kBig, -kBig, -kBig);
    track->SetD0_noTrackVertex(-kBig);
    track->SetHasHodoHit_noTrackVertex(0);
    track->SetBestHodoHit_noTrackVertex(nullptr);
    track->SetIsGoodTrack_noTrackVertex(kFALSE);

    TVector3 v_hit1(track->GetX1(), track->GetY1(), track->GetZ1());
    TVector3 v_hit2(track->GetX2(), track->GetY2(), track->GetZ2());

    // Direction seed from the GEM space points
    Double_t gem_theta = TMath::ACos((v_hit2.Z() - v_hit1.Z()) / (v_hit2 - v_hit1).Mag());
    Double_t gem_phi   = TMath::ATan2((v_hit2.Y() - v_hit1.Y()), (v_hit2.X() - v_hit1.X()));

    // Best chi-square per hypothesis, indexed [category][subtype]:
    //   category 1 -> a single hodoscope plane; category 2 -> two hodoscope planes
    //   for category 2: subtype 0 = Front&Back, 1 = Front-only, 2 = Back-only
    Double_t bestchisq[3][3] = {{kBig, kBig, kBig}, {kBig, kBig, kBig}, {kBig, kBig, kBig}};
    bool usedHodoHit[3][3]   = {{kFALSE, kFALSE, kFALSE}, {kFALSE, kFALSE, kFALSE}, {kFALSE, kFALSE, kFALSE}};
    Double_t best_dir[3][3][2];    // theta, phi
    Double_t best_anchor[3][3][2]; // y0, z0 (line y,z at x = GEM1 x)
    int bestHodoHitIndex[3][3] = {{-1, -1, -1}, {-1, -1, -1}, {-1, -1, -1}};

    // Loop over hodoscope hits to find the best match to the GEM track (no vertex)
    for (int iHodoHit = 0; iHodoHit < nHits; iHodoHit++) {
      THcGoodLADHit *goodHit = static_cast<THcGoodLADHit *>(LADHits_unfiltered->At(iHodoHit));
      if (goodHit == nullptr)
        continue;
      int num_hodo_hits = 0;
      if (goodHit->GetPlaneHit0() >= 0 && goodHit->GetPlaneHit0() < 999)
        num_hodo_hits++;
      if (goodHit->GetPlaneHit1() >= 0 && goodHit->GetPlaneHit1() < 999)
        num_hodo_hits++;
      if (num_hodo_hits == 0)
        continue; // no-vertex tracking requires a hodoscope hit

      if (num_hodo_hits == 1) {
        TVector3 hodo_hit_pos;
        if (goodHit->GetPlaneHit0() >= 0 && goodHit->GetPlaneHit0() < 999) {
          hodo_hit_pos = fHodoscope->GetHitPositionLab(goodHit->GetPlaneHit0(), goodHit->GetPaddleHit0(),
                                                       goodHit->GetHitYPosHit0());
        } else {
          hodo_hit_pos = fHodoscope->GetHitPositionLab(goodHit->GetPlaneHit1(), goodHit->GetPaddleHit1(),
                                                       goodHit->GetHitYPosHit1());
        }
        std::vector<TVector3> v_hits = {v_hit1, v_hit2, hodo_hit_pos};
        std::vector<Double_t> v_res  = {fSigma_GEM, fSigma_GEM, fSigma_Hodo};
        Double_t dir[3]              = {gem_theta, gem_phi, 0.0};
        Double_t anchor[2]           = {v_hit1.Y(), v_hit1.Z()};
        Double_t chisq               = FitTrack_noTrackVertex(v_hits, v_res, dir, anchor);
        if (chisq >= 0 && chisq < bestchisq[1][0]) {
          bestchisq[1][0]        = chisq;
          best_dir[1][0][0]      = dir[0];
          best_dir[1][0][1]      = dir[1];
          best_anchor[1][0][0]   = anchor[0];
          best_anchor[1][0][1]   = anchor[1];
          bestHodoHitIndex[1][0] = iHodoHit;
          usedHodoHit[1][0]      = kTRUE;
        }
      } else { // num_hodo_hits == 2
        TVector3 hodo_hit_posF =
            fHodoscope->GetHitPositionLab(goodHit->GetPlaneHit0(), goodHit->GetPaddleHit0(), goodHit->GetHitYPosHit0());
        TVector3 hodo_hit_posB =
            fHodoscope->GetHitPositionLab(goodHit->GetPlaneHit1(), goodHit->GetPaddleHit1(), goodHit->GetHitYPosHit1());

        // Front & Back -> [2][0]
        {
          std::vector<TVector3> v_hits = {v_hit1, v_hit2, hodo_hit_posF, hodo_hit_posB};
          std::vector<Double_t> v_res  = {fSigma_GEM, fSigma_GEM, fSigma_Hodo, fSigma_Hodo};
          Double_t dir[3]              = {gem_theta, gem_phi, 0.0};
          Double_t anchor[2]           = {v_hit1.Y(), v_hit1.Z()};
          Double_t chisq               = FitTrack_noTrackVertex(v_hits, v_res, dir, anchor);
          if (chisq >= 0 && chisq < bestchisq[2][0]) {
            bestchisq[2][0]        = chisq;
            best_dir[2][0][0]      = dir[0];
            best_dir[2][0][1]      = dir[1];
            best_anchor[2][0][0]   = anchor[0];
            best_anchor[2][0][1]   = anchor[1];
            bestHodoHitIndex[2][0] = iHodoHit;
            usedHodoHit[2][0]      = kTRUE;
          }
        }
        // Front only -> [2][1]
        {
          std::vector<TVector3> v_hits = {v_hit1, v_hit2, hodo_hit_posF};
          std::vector<Double_t> v_res  = {fSigma_GEM, fSigma_GEM, fSigma_Hodo};
          Double_t dir[3]              = {gem_theta, gem_phi, 0.0};
          Double_t anchor[2]           = {v_hit1.Y(), v_hit1.Z()};
          Double_t chisq               = FitTrack_noTrackVertex(v_hits, v_res, dir, anchor);
          if (chisq >= 0 && chisq < bestchisq[2][1]) {
            bestchisq[2][1]        = chisq;
            best_dir[2][1][0]      = dir[0];
            best_dir[2][1][1]      = dir[1];
            best_anchor[2][1][0]   = anchor[0];
            best_anchor[2][1][1]   = anchor[1];
            bestHodoHitIndex[2][1] = iHodoHit;
            usedHodoHit[2][1]      = kTRUE;
          }
        }
        // Back only -> [2][2]
        {
          std::vector<TVector3> v_hits = {v_hit1, v_hit2, hodo_hit_posB};
          std::vector<Double_t> v_res  = {fSigma_GEM, fSigma_GEM, fSigma_Hodo};
          Double_t dir[3]              = {gem_theta, gem_phi, 0.0};
          Double_t anchor[2]           = {v_hit1.Y(), v_hit1.Z()};
          Double_t chisq               = FitTrack_noTrackVertex(v_hits, v_res, dir, anchor);
          if (chisq >= 0 && chisq < bestchisq[2][2]) {
            bestchisq[2][2]        = chisq;
            best_dir[2][2][0]      = dir[0];
            best_dir[2][2][1]      = dir[1];
            best_anchor[2][2][0]   = anchor[0];
            best_anchor[2][2][1]   = anchor[1];
            bestHodoHitIndex[2][2] = iHodoHit;
            usedHodoHit[2][2]      = kTRUE;
          }
        }
      }
    }

    // Pick the best hypothesis. First find the best single-plane fit, then
    // prefer the two-plane (F&B) fit if it is no worse than the best single
    // fit by more than fchisq_cut[0].
    int hit_index      = -1;
    Double_t bc_single = kBig;
    if (usedHodoHit[1][0] && bestchisq[1][0] < bc_single) {
      bc_single = bestchisq[1][0];
      hit_index = 1 * 3 + 0;
    }
    if (usedHodoHit[2][1] && bestchisq[2][1] < bc_single) {
      bc_single = bestchisq[2][1];
      hit_index = 2 * 3 + 1;
    }
    if (usedHodoHit[2][2] && bestchisq[2][2] < bc_single) {
      bc_single = bestchisq[2][2];
      hit_index = 2 * 3 + 2;
    }
    if (usedHodoHit[2][0] && (hit_index < 0 || (bestchisq[2][0] - bc_single) <= fchisq_cut[0])) {
      hit_index = 2 * 3 + 0;
    }
    if (hit_index < 0)
      continue; // no valid no-vertex hypothesis for this track -> leave invalid defaults

    int cat        = hit_index / 3;
    int sub        = hit_index % 3;
    Double_t theta = best_dir[cat][sub][0];
    Double_t phi   = best_dir[cat][sub][1];
    Double_t y0    = best_anchor[cat][sub][0];
    Double_t z0    = best_anchor[cat][sub][1];
    Double_t x_ref = v_hit1.X(); // anchor plane used by FitTrack_noTrackVertex

    Double_t sx = TMath::Sin(theta) * TMath::Cos(phi);
    Double_t sy = TMath::Sin(theta) * TMath::Sin(phi);
    Double_t sz = TMath::Cos(theta);

    track->SetChisq_noTrackVertex(bestchisq[cat][sub]);
    track->SetAngles_noTrackVertex(theta, phi);

    // Reproject the fitted line back toward the beamline (target region).
    // Projected vertex = point where the track crosses the x = 0 plane.
    // D0 = distance of closest approach of the track to the beam (z) axis.
    if (TMath::Abs(sx) > 1e-9) {
      Double_t t0    = -x_ref / sx;
      Double_t projY = y0 + t0 * sy;
      Double_t projZ = z0 + t0 * sz;
      Double_t d0    = TMath::Abs(x_ref * sy - y0 * sx) / TMath::Sqrt(sx * sx + sy * sy);
      track->SetProjVertex_noTrackVertex(0.0, projY, projZ);
      track->SetD0_noTrackVertex(d0);
    }

    int nplanes = (cat == 2 && sub == 0) ? 2 : 1;
    track->SetHasHodoHit_noTrackVertex(nplanes);
    THcGoodLADHit *besthit = static_cast<THcGoodLADHit *>(LADHits_unfiltered->At(bestHodoHitIndex[cat][sub]));
    track->SetBestHodoHit_noTrackVertex(besthit);
    if (besthit) {
      // per-hit association from the no-vertex fit (parallel to the vertex-based track_id / trk_chiSqr)
      besthit->SetTrackID_noTrackVertex(track->GetTrackID());
      besthit->SetTrkChiSqr_noTrackVertex(track->GetChisq_noTrackVertex());
      // dof = 3*(#GEM + #hodo) - 4 params (full 3D residual, free line, no vertex)
      besthit->SetTrkNdof_noTrackVertex(3 * (2 + nplanes) - 4);
    }

    bool good     = kTRUE;
    Double_t d0nv = track->GetD0_noTrackVertex();
    if (!(d0nv > 0 && d0nv < fD0Cut_noVertex))
      good = kFALSE;
    if (bestchisq[cat][sub] < 0)
      good = kFALSE;
    if (track->GetdT() > fTrk_dtCut)
      good = kFALSE;
    track->SetIsGoodTrack_noTrackVertex(good);
  }

  //////////////////////////////////////////////////////////////////////////////
  // x-z tracking (vertex-constrained)
  // Identical to the vertex-constrained loop above, but the chi-square ignores
  // the y position (FitTrack called with use_y = false). Results are stored in
  // the parallel _xz members so they can be compared against the full 3D fit.
  std::vector<bool> isGoodTrack_xz(ntracks, false);
  for (Int_t iTrack = 0; iTrack < ntracks; iTrack++) {
    // vertex-constrained x-z tracking pass; enabled via the ldo_xz_tracking DB flag (default off).
    if (!fDoXZTracking)
      break;
    THcLADGEMTrack *track = static_cast<THcLADGEMTrack *>(fGEMTracks->At(iTrack));
    if (track == nullptr)
      continue;

    // Get the hit positions
    TVector3 v_hit1(track->GetX1(), track->GetY1(), track->GetZ1());
    TVector3 v_hit2(track->GetX2(), track->GetY2(), track->GetZ2());

    if (fVertexModule && fVertexModule->HasVertex()) {
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

    Double_t gemdir[3];
    gemdir[0] = TMath::ACos((v_hit2.Z() - v_hit1.Z()) / (v_hit2 - v_hit1).Mag());
    gemdir[1] = TMath::ATan2((v_hit2.Y() - v_hit1.Y()), (v_hit2.X() - v_hit1.X()));
    gemdir[2] = vertex.Z();

    Double_t bestchisq[3][3]      = {{kBig, kBig, kBig}, {kBig, kBig, kBig}, {kBig, kBig, kBig}};
    bool usedHodoHit[3][3]        = {{kFALSE, kFALSE, kFALSE}, {kFALSE, kFALSE, kFALSE}, {kFALSE, kFALSE, kFALSE}};
    Double_t best_gemdir[3][3][3] = {}; // zero-initialized so an unselected entry can never be read uninitialized
    int bestHodoHitIndex[3][3]    = {{-1, -1, -1}, {-1, -1, -1}, {-1, -1, -1}};
    Double_t dir[3];
    dir[0]            = gemdir[0];
    dir[1]            = gemdir[1];
    dir[2]            = gemdir[2];
    bestchisq[0][0]   = FitTrack(vertex, {v_hit1, v_hit2}, {fSigma_GEM, fSigma_GEM}, dir, false);
    usedHodoHit[0][0] = kTRUE;
    if (bestchisq[0][0] < 0) {
      track->SetIsGoodTrack_xz(kFALSE);
      track->SetChisq_xz(bestchisq[0][0]);
      track->SetAngles_xz(-kBig, -kBig);
      track->SetProjVertex_xz(-kBig, -kBig, -kBig);
      track->SetD0_xz(-kBig);
      track->SetHasHodoHit_xz(kFALSE);
      track->SetBestHodoHit_xz(nullptr);
      continue;
    }
    best_gemdir[0][0][0]         = dir[0];
    best_gemdir[0][0][1]         = dir[1];
    best_gemdir[0][0][2]         = dir[2];
    int hit_index_for_best_chisq = -1;
    // loop over hodoscope hits to find the best match to the track projection
    for (int iHodoHit = 0; iHodoHit < nHits; iHodoHit++) {
      THcGoodLADHit *goodHit = static_cast<THcGoodLADHit *>(LADHits_unfiltered->At(iHodoHit));
      if (goodHit == nullptr)
        continue;
      int num_hodo_hits = 0;
      if (goodHit->GetPlaneHit0() >= 0 && goodHit->GetPlaneHit0() < 999) {
        num_hodo_hits++;
      }
      if (goodHit->GetPlaneHit1() >= 0 && goodHit->GetPlaneHit1() < 999) {
        num_hodo_hits++;
      }

      std::vector<TVector3> v_hits;
      std::vector<Double_t> v_resolutions;
      v_hits.push_back(v_hit1);
      v_hits.push_back(v_hit2);
      v_resolutions.push_back(fSigma_GEM);
      v_resolutions.push_back(fSigma_GEM);

      if (num_hodo_hits == 1) {
        if (goodHit->GetPlaneHit0() >= 0 && goodHit->GetPlaneHit0() < 999) {
          TVector3 hodo_hit_pos = fHodoscope->GetHitPositionLab(goodHit->GetPlaneHit0(), goodHit->GetPaddleHit0(),
                                                                goodHit->GetHitYPosHit0());
          v_hits.push_back(hodo_hit_pos);
          v_resolutions.push_back(fSigma_Hodo);
        } else if (goodHit->GetPlaneHit1() >= 0 && goodHit->GetPlaneHit1() < 999) {
          TVector3 hodo_hit_pos = fHodoscope->GetHitPositionLab(goodHit->GetPlaneHit1(), goodHit->GetPaddleHit1(),
                                                                goodHit->GetHitYPosHit1());
          v_hits.push_back(hodo_hit_pos);
          v_resolutions.push_back(fSigma_Hodo);
        }
      } else if (num_hodo_hits == 2) {
        TVector3 hodo_hit_pos1 =
            fHodoscope->GetHitPositionLab(goodHit->GetPlaneHit0(), goodHit->GetPaddleHit0(), goodHit->GetHitYPosHit0());
        TVector3 hodo_hit_pos2 =
            fHodoscope->GetHitPositionLab(goodHit->GetPlaneHit1(), goodHit->GetPaddleHit1(), goodHit->GetHitYPosHit1());
        v_hits.push_back(hodo_hit_pos1);
        v_hits.push_back(hodo_hit_pos2);
        v_resolutions.push_back(fSigma_Hodo);
        v_resolutions.push_back(fSigma_Hodo);
      }
      Double_t dir[3];
      dir[0]         = gemdir[0];
      dir[1]         = gemdir[1];
      dir[2]         = gemdir[2];
      Double_t chisq = FitTrack(vertex, v_hits, v_resolutions, dir, false);
      if (num_hodo_hits == 1) {
        if (chisq >= 0 && chisq < bestchisq[1][0]) {
          bestchisq[1][0]        = chisq;
          best_gemdir[1][0][0]   = dir[0];
          best_gemdir[1][0][1]   = dir[1];
          best_gemdir[1][0][2]   = dir[2];
          bestHodoHitIndex[1][0] = iHodoHit;
          usedHodoHit[1][0]      = kTRUE;
        }
      } else if (num_hodo_hits == 2) {
        if (chisq >= 0 && chisq < bestchisq[2][0]) {
          bestchisq[2][0]        = chisq;
          best_gemdir[2][0][0]   = dir[0];
          best_gemdir[2][0][1]   = dir[1];
          best_gemdir[2][0][2]   = dir[2];
          bestHodoHitIndex[2][0] = iHodoHit;
          usedHodoHit[2][0]      = kTRUE;
        }
        v_hits.pop_back();
        v_hits.pop_back();
        v_resolutions.pop_back();
        dir[0] = gemdir[0];
        dir[1] = gemdir[1];
        dir[2] = gemdir[2];
        TVector3 hodo_hit_pos1 =
            fHodoscope->GetHitPositionLab(goodHit->GetPlaneHit0(), goodHit->GetPaddleHit0(), goodHit->GetHitYPosHit0());
        TVector3 hodo_hit_pos2 =
            fHodoscope->GetHitPositionLab(goodHit->GetPlaneHit1(), goodHit->GetPaddleHit1(), goodHit->GetHitYPosHit1());
        v_hits.push_back(hodo_hit_pos1);
        Double_t chisq1 = FitTrack(vertex, v_hits, v_resolutions, dir, false);
        if (chisq1 >= 0 && chisq1 < bestchisq[2][1]) {
          bestchisq[2][1]        = chisq1;
          best_gemdir[2][1][0]   = dir[0];
          best_gemdir[2][1][1]   = dir[1];
          best_gemdir[2][1][2]   = dir[2];
          bestHodoHitIndex[2][1] = iHodoHit;
          usedHodoHit[2][1]      = kTRUE;
        }
        v_hits.pop_back();
        v_hits.push_back(hodo_hit_pos2);
        dir[0]          = gemdir[0];
        dir[1]          = gemdir[1];
        dir[2]          = gemdir[2];
        Double_t chisq2 = FitTrack(vertex, v_hits, v_resolutions, dir, false);
        if (chisq2 >= 0 && chisq2 < bestchisq[2][2]) {
          bestchisq[2][2]        = chisq2;
          best_gemdir[2][2][0]   = dir[0];
          best_gemdir[2][2][1]   = dir[1];
          best_gemdir[2][2][2]   = dir[2];
          bestHodoHitIndex[2][2] = iHodoHit;
          usedHodoHit[2][2]      = kTRUE;
        }
      }
      v_hits.clear();
      v_resolutions.clear();
    }
    // check between f-only, b-only, and single, which has the best chisq
    hit_index_for_best_chisq = 2 * 3 + 0; // default to 2 hodo hit
    if (usedHodoHit[2][1] && bestchisq[2][1] < bestchisq[2][2] && bestchisq[2][1] < bestchisq[1][0]) {
      hit_index_for_best_chisq = 2 * 3 + 1;
    } else if (usedHodoHit[2][2] && bestchisq[2][2] < bestchisq[2][1] && bestchisq[2][2] < bestchisq[1][0]) {
      hit_index_for_best_chisq = 2 * 3 + 2;
    } else if (usedHodoHit[1][0]) {
      hit_index_for_best_chisq = 1 * 3 + 0;
    }
    if (usedHodoHit[2][0] &&
        (bestchisq[2][0] - bestchisq[hit_index_for_best_chisq / 3][hit_index_for_best_chisq % 3] <= fchisq_cut[0])) {
      hit_index_for_best_chisq = 2 * 3 + 0;
    } else if (usedHodoHit[0][0] && (bestchisq[hit_index_for_best_chisq / 3][hit_index_for_best_chisq % 3] -
                                     bestchisq[0][0]) > fchisq_cut[1]) {
      hit_index_for_best_chisq = 0 * 3 + 0;
    }
    track->SetChisq_xz(bestchisq[hit_index_for_best_chisq / 3][hit_index_for_best_chisq % 3]);
    track->SetAngles_xz(best_gemdir[hit_index_for_best_chisq / 3][hit_index_for_best_chisq % 3][0],
                        best_gemdir[hit_index_for_best_chisq / 3][hit_index_for_best_chisq % 3][1]);
    track->SetProjVertex_xz(vertex.X(), vertex.Y(),
                            best_gemdir[hit_index_for_best_chisq / 3][hit_index_for_best_chisq % 3][2]);
    Double_t d0_xz = fabs(vertex.Z() - best_gemdir[hit_index_for_best_chisq / 3][hit_index_for_best_chisq % 3][2]);
    track->SetD0_xz(d0_xz);

    track->SetHasHodoHit_xz(hit_index_for_best_chisq / 3);
    if (hit_index_for_best_chisq >= 1) {
      THcGoodLADHit *besthit = static_cast<THcGoodLADHit *>(
          LADHits_unfiltered->At(bestHodoHitIndex[hit_index_for_best_chisq / 3][hit_index_for_best_chisq % 3]));
      if (besthit == nullptr) {
        track->SetBestHodoHit_xz(nullptr);
        continue;
      }
      if (hit_index_for_best_chisq / 3 == 2) {
        int ntmp = 0;
        if (hit_index_for_best_chisq % 3 == 1 || hit_index_for_best_chisq % 3 == 0) {
          ntmp++;
        }
        if (hit_index_for_best_chisq % 3 == 2 || hit_index_for_best_chisq % 3 == 0) {
          ntmp++;
        }
        track->SetHasHodoHit_xz(ntmp);
        besthit->SetTrackID_xz(track->GetTrackID());
        besthit->SetTrkChiSqr_xz(track->GetChisq_xz());
        // dof = 2*(#GEM + #hodo) - 3 params (x-z residual only, vertex-constrained)
        besthit->SetTrkNdof_xz(2 * (2 + ntmp) - 3);
      } else if (hit_index_for_best_chisq / 3 == 1) {
        besthit->SetTrackID_xz(track->GetTrackID());
        besthit->SetTrkChiSqr_xz(track->GetChisq_xz());
        besthit->SetTrkNdof_xz(2 * (2 + 1) - 3);
        track->SetHasHodoHit_xz(1);
      }
      if (hit_index_for_best_chisq / 3 != 0) {
        track->SetBestHodoHit_xz(besthit);
        isGoodTrack_xz[iTrack] = true;
      }
    } else {
      track->SetBestHodoHit_xz(nullptr);
      track->SetHasHodoHit_xz(0);
      isGoodTrack_xz[iTrack] = true;
    }

    if (track->GetGoodD0()) {
      if (!(track->GetD0_xz() > 0 && track->GetD0_xz() < fD0Cut_wVertex)) {
        isGoodTrack_xz[iTrack] = false;
      }
    } else {
      if (!(track->GetD0_xz() > 0 && track->GetD0_xz() < fD0Cut_noVertex)) {
        isGoodTrack_xz[iTrack] = false;
      }
    }
    if (track->GetChisq_xz() < 0) {
      isGoodTrack_xz[iTrack] = false;
    }
    if (track->GetdT() > fTrk_dtCut) {
      isGoodTrack_xz[iTrack] = false;
    }
    track->SetIsGoodTrack_xz(isGoodTrack_xz[iTrack]);
  }

  //////////////////////////////////////////////////////////////////////////////
  // x-z tracking, no vertex
  // Identical to the no-vertex loop above, but the chi-square ignores the y
  // position (FitTrack_noTrackVertex called with use_y = false). Results are
  // stored in the parallel _noTrackVertex_xz members.
  for (Int_t iTrack = 0; iTrack < ntracks; iTrack++) {
    // no-vertex x-z tracking pass; needs BOTH ldo_noVertex_tracking and ldo_xz_tracking.
    if (!(fDoNoVertexTracking && fDoXZTracking))
      break;
    THcLADGEMTrack *track = static_cast<THcLADGEMTrack *>(fGEMTracks->At(iTrack));
    if (track == nullptr)
      continue;

    // Default (invalid) values: overwritten below only if a valid fit is found
    track->SetChisq_noTrackVertex_xz(-kBig);
    track->SetAngles_noTrackVertex_xz(-kBig, -kBig);
    track->SetProjVertex_noTrackVertex_xz(-kBig, -kBig, -kBig);
    track->SetD0_noTrackVertex_xz(-kBig);
    track->SetHasHodoHit_noTrackVertex_xz(0);
    track->SetBestHodoHit_noTrackVertex_xz(nullptr);
    track->SetIsGoodTrack_noTrackVertex_xz(kFALSE);

    TVector3 v_hit1(track->GetX1(), track->GetY1(), track->GetZ1());
    TVector3 v_hit2(track->GetX2(), track->GetY2(), track->GetZ2());

    // Direction seed from the GEM space points
    Double_t gem_theta = TMath::ACos((v_hit2.Z() - v_hit1.Z()) / (v_hit2 - v_hit1).Mag());
    Double_t gem_phi   = TMath::ATan2((v_hit2.Y() - v_hit1.Y()), (v_hit2.X() - v_hit1.X()));

    Double_t bestchisq[3][3] = {{kBig, kBig, kBig}, {kBig, kBig, kBig}, {kBig, kBig, kBig}};
    bool usedHodoHit[3][3]   = {{kFALSE, kFALSE, kFALSE}, {kFALSE, kFALSE, kFALSE}, {kFALSE, kFALSE, kFALSE}};
    Double_t best_dir[3][3][2];    // theta, phi
    Double_t best_anchor[3][3][2]; // y0, z0 (line y,z at x = GEM1 x)
    int bestHodoHitIndex[3][3] = {{-1, -1, -1}, {-1, -1, -1}, {-1, -1, -1}};

    // Loop over hodoscope hits to find the best match to the GEM track (no vertex)
    for (int iHodoHit = 0; iHodoHit < nHits; iHodoHit++) {
      THcGoodLADHit *goodHit = static_cast<THcGoodLADHit *>(LADHits_unfiltered->At(iHodoHit));
      if (goodHit == nullptr)
        continue;
      int num_hodo_hits = 0;
      if (goodHit->GetPlaneHit0() >= 0 && goodHit->GetPlaneHit0() < 999)
        num_hodo_hits++;
      if (goodHit->GetPlaneHit1() >= 0 && goodHit->GetPlaneHit1() < 999)
        num_hodo_hits++;
      if (num_hodo_hits == 0)
        continue; // no-vertex tracking requires a hodoscope hit

      if (num_hodo_hits == 1) {
        TVector3 hodo_hit_pos;
        if (goodHit->GetPlaneHit0() >= 0 && goodHit->GetPlaneHit0() < 999) {
          hodo_hit_pos = fHodoscope->GetHitPositionLab(goodHit->GetPlaneHit0(), goodHit->GetPaddleHit0(),
                                                       goodHit->GetHitYPosHit0());
        } else {
          hodo_hit_pos = fHodoscope->GetHitPositionLab(goodHit->GetPlaneHit1(), goodHit->GetPaddleHit1(),
                                                       goodHit->GetHitYPosHit1());
        }
        std::vector<TVector3> v_hits = {v_hit1, v_hit2, hodo_hit_pos};
        std::vector<Double_t> v_res  = {fSigma_GEM, fSigma_GEM, fSigma_Hodo};
        Double_t dir[3]              = {gem_theta, gem_phi, 0.0};
        Double_t anchor[2]           = {v_hit1.Y(), v_hit1.Z()};
        Double_t chisq               = FitTrack_noTrackVertex(v_hits, v_res, dir, anchor, false);
        if (chisq >= 0 && chisq < bestchisq[1][0]) {
          bestchisq[1][0]        = chisq;
          best_dir[1][0][0]      = dir[0];
          best_dir[1][0][1]      = dir[1];
          best_anchor[1][0][0]   = anchor[0];
          best_anchor[1][0][1]   = anchor[1];
          bestHodoHitIndex[1][0] = iHodoHit;
          usedHodoHit[1][0]      = kTRUE;
        }
      } else { // num_hodo_hits == 2
        TVector3 hodo_hit_posF =
            fHodoscope->GetHitPositionLab(goodHit->GetPlaneHit0(), goodHit->GetPaddleHit0(), goodHit->GetHitYPosHit0());
        TVector3 hodo_hit_posB =
            fHodoscope->GetHitPositionLab(goodHit->GetPlaneHit1(), goodHit->GetPaddleHit1(), goodHit->GetHitYPosHit1());

        // Front & Back -> [2][0]
        {
          std::vector<TVector3> v_hits = {v_hit1, v_hit2, hodo_hit_posF, hodo_hit_posB};
          std::vector<Double_t> v_res  = {fSigma_GEM, fSigma_GEM, fSigma_Hodo, fSigma_Hodo};
          Double_t dir[3]              = {gem_theta, gem_phi, 0.0};
          Double_t anchor[2]           = {v_hit1.Y(), v_hit1.Z()};
          Double_t chisq               = FitTrack_noTrackVertex(v_hits, v_res, dir, anchor, false);
          if (chisq >= 0 && chisq < bestchisq[2][0]) {
            bestchisq[2][0]        = chisq;
            best_dir[2][0][0]      = dir[0];
            best_dir[2][0][1]      = dir[1];
            best_anchor[2][0][0]   = anchor[0];
            best_anchor[2][0][1]   = anchor[1];
            bestHodoHitIndex[2][0] = iHodoHit;
            usedHodoHit[2][0]      = kTRUE;
          }
        }
        // Front only -> [2][1]
        {
          std::vector<TVector3> v_hits = {v_hit1, v_hit2, hodo_hit_posF};
          std::vector<Double_t> v_res  = {fSigma_GEM, fSigma_GEM, fSigma_Hodo};
          Double_t dir[3]              = {gem_theta, gem_phi, 0.0};
          Double_t anchor[2]           = {v_hit1.Y(), v_hit1.Z()};
          Double_t chisq               = FitTrack_noTrackVertex(v_hits, v_res, dir, anchor, false);
          if (chisq >= 0 && chisq < bestchisq[2][1]) {
            bestchisq[2][1]        = chisq;
            best_dir[2][1][0]      = dir[0];
            best_dir[2][1][1]      = dir[1];
            best_anchor[2][1][0]   = anchor[0];
            best_anchor[2][1][1]   = anchor[1];
            bestHodoHitIndex[2][1] = iHodoHit;
            usedHodoHit[2][1]      = kTRUE;
          }
        }
        // Back only -> [2][2]
        {
          std::vector<TVector3> v_hits = {v_hit1, v_hit2, hodo_hit_posB};
          std::vector<Double_t> v_res  = {fSigma_GEM, fSigma_GEM, fSigma_Hodo};
          Double_t dir[3]              = {gem_theta, gem_phi, 0.0};
          Double_t anchor[2]           = {v_hit1.Y(), v_hit1.Z()};
          Double_t chisq               = FitTrack_noTrackVertex(v_hits, v_res, dir, anchor, false);
          if (chisq >= 0 && chisq < bestchisq[2][2]) {
            bestchisq[2][2]        = chisq;
            best_dir[2][2][0]      = dir[0];
            best_dir[2][2][1]      = dir[1];
            best_anchor[2][2][0]   = anchor[0];
            best_anchor[2][2][1]   = anchor[1];
            bestHodoHitIndex[2][2] = iHodoHit;
            usedHodoHit[2][2]      = kTRUE;
          }
        }
      }
    }

    // Pick the best hypothesis (same logic as the 3D no-vertex loop).
    int hit_index      = -1;
    Double_t bc_single = kBig;
    if (usedHodoHit[1][0] && bestchisq[1][0] < bc_single) {
      bc_single = bestchisq[1][0];
      hit_index = 1 * 3 + 0;
    }
    if (usedHodoHit[2][1] && bestchisq[2][1] < bc_single) {
      bc_single = bestchisq[2][1];
      hit_index = 2 * 3 + 1;
    }
    if (usedHodoHit[2][2] && bestchisq[2][2] < bc_single) {
      bc_single = bestchisq[2][2];
      hit_index = 2 * 3 + 2;
    }
    if (usedHodoHit[2][0] && (hit_index < 0 || (bestchisq[2][0] - bc_single) <= fchisq_cut[0])) {
      hit_index = 2 * 3 + 0;
    }
    if (hit_index < 0)
      continue; // no valid no-vertex hypothesis for this track -> leave invalid defaults

    int cat        = hit_index / 3;
    int sub        = hit_index % 3;
    Double_t theta = best_dir[cat][sub][0];
    Double_t phi   = best_dir[cat][sub][1];
    Double_t y0    = best_anchor[cat][sub][0];
    Double_t z0    = best_anchor[cat][sub][1];
    Double_t x_ref = v_hit1.X(); // anchor plane used by FitTrack_noTrackVertex

    Double_t sx = TMath::Sin(theta) * TMath::Cos(phi);
    Double_t sy = TMath::Sin(theta) * TMath::Sin(phi);
    Double_t sz = TMath::Cos(theta);

    track->SetChisq_noTrackVertex_xz(bestchisq[cat][sub]);
    track->SetAngles_noTrackVertex_xz(theta, phi);

    // Reproject the fitted line back toward the beamline (same formula as the 3D
    // no-vertex fit). NB: with y dropped from the chi-square the y0/phi (hence sy)
    // are only weakly constrained, so projY and this D0 carry the GEM-seeded y.
    if (TMath::Abs(sx) > 1e-9) {
      Double_t t0    = -x_ref / sx;
      Double_t projY = y0 + t0 * sy;
      Double_t projZ = z0 + t0 * sz;
      Double_t d0    = TMath::Abs(x_ref * sy - y0 * sx) / TMath::Sqrt(sx * sx + sy * sy);
      track->SetProjVertex_noTrackVertex_xz(0.0, projY, projZ);
      track->SetD0_noTrackVertex_xz(d0);
    }

    int nplanes = (cat == 2 && sub == 0) ? 2 : 1;
    track->SetHasHodoHit_noTrackVertex_xz(nplanes);
    THcGoodLADHit *besthit = static_cast<THcGoodLADHit *>(LADHits_unfiltered->At(bestHodoHitIndex[cat][sub]));
    track->SetBestHodoHit_noTrackVertex_xz(besthit);
    if (besthit) {
      // per-hit association from the no-vertex x-z fit
      besthit->SetTrackID_noTrackVertex_xz(track->GetTrackID());
      besthit->SetTrkChiSqr_noTrackVertex_xz(track->GetChisq_noTrackVertex_xz());
      // dof = 2*(#GEM + #hodo) - 4 params (x-z residual only, free line, no vertex)
      besthit->SetTrkNdof_noTrackVertex_xz(2 * (2 + nplanes) - 4);
    }

    bool good     = kTRUE;
    Double_t d0nv = track->GetD0_noTrackVertex_xz();
    if (!(d0nv > 0 && d0nv < fD0Cut_noVertex))
      good = kFALSE;
    if (bestchisq[cat][sub] < 0)
      good = kFALSE;
    if (track->GetdT() > fTrk_dtCut)
      good = kFALSE;
    track->SetIsGoodTrack_noTrackVertex_xz(good);
  }

  // Refresh the event vertex for the hit ToF RF correction. The per-track loops
  // above leave the vertex member in a track-dependent state; use the current
  // event's reaction-point vertex here instead.
  if (fVertexModule && fVertexModule->HasVertex())
    vertex = fVertexModule->GetVertex();
  else
    vertex.SetXYZ(0, 0, 0);

  //////////////////////////////////////////////////////////////////////////////
  // Calculate ToF and RF-corrected ToF for good LAD hits

  for (Int_t i = 0; i < nHits; i++) {
    THcGoodLADHit *goodhit = static_cast<THcGoodLADHit *>(LADHits_unfiltered->At(i));
    if (goodhit == nullptr)
      continue;

    Double_t tof, tof_rfcorr, beta, alpha;
    if (goodhit->GetPlaneHit0() >= 0 && goodhit->GetPlaneHit0() < 999) {
      tof = CalculateToF(goodhit->GetHitTimeHit0());
      goodhit->SetHitTOF(0, tof);
      tof_rfcorr = CalculateTOFRFcorr(goodhit->GetHitTimeHit0());
      goodhit->SetHitTOFRFcorr(0, tof_rfcorr);
      // Calculate beta, alpha, etc. for the first plane
    }
    if (goodhit->GetPlaneHit1() >= 0 && goodhit->GetPlaneHit1() < 999) {
      tof = CalculateToF(goodhit->GetHitTimeHit1());
      goodhit->SetHitTOF(1, tof);
      tof_rfcorr = CalculateTOFRFcorr(goodhit->GetHitTimeHit1());
      goodhit->SetHitTOFRFcorr(1, tof_rfcorr);
      // Calculate beta, alpha, etc. for the second plane
    }
  }

  //////////////////////////////////////////////////////////////////////////////
  // Calculate tof for raw planes (not good lad hits)
  // Not sure this is necessary. TODO: revisit and possibly remove

  Int_t nPlanes = fHodoscope->GetNPlanes();

  // Loop through the planes and get LAD hits
  for (Int_t plane = 0; plane < nPlanes; plane++) {
    TClonesArray *planeHits = fHodoscope->GetLADHits(plane);
    Int_t nPlaneHits        = planeHits->GetLast() + 1;

    for (Int_t i = 0; i < nPlaneHits; i++) {
      THcLADHodoHit *hit = static_cast<THcLADHodoHit *>(planeHits->At(i));
      if (hit == nullptr)
        continue;

      // Calculate time of flight for each hit
      Double_t tof_avg = CalculateToF(hit->GetScinCorrectedTime());
      Double_t tof_top = CalculateToF(hit->GetTopCorrectedTime());
      Double_t tof_btm = CalculateToF(hit->GetBtmCorrectedTime());
      hit->SetTOFCorrectedTimes(tof_top, tof_btm, tof_avg);

      tof_avg = CalculateTOFRFcorr(hit->GetScinCorrectedTime());
      tof_top = CalculateTOFRFcorr(hit->GetTopCorrectedTime());
      tof_btm = CalculateTOFRFcorr(hit->GetBtmCorrectedTime());
      hit->SetTOF_RF_CorrectedTimes(tof_top, tof_btm, tof_avg);
    }
  }

  MakeProtonCut(LADHits_unfiltered);

  // Optional alternative tracking that seeds from individual GEM strip clusters
  // (unpaired U/V) instead of 2D space points. Enabled with ldo_1Dcluster_tracking.
  if (fDo1DClusterTracking)
    Do1DClusterTracking();

  return kOK;
}
//_____________________________________________________________________________
void THcLADKine::MakeProtonCut(TClonesArray *hits) {
  // Front plane (index 0): is_proton if edep_MeV > 100
  // Back plane  (index 1): is_proton if point (dt, edep_MeV) is above the line
  //   through (0,300) and (5,0), i.e. edep_MeV > 300 - 60*dt
  //   where dt = hit_time[1] - hit_time[0] in ns

  Int_t nHits = hits->GetLast() + 1;
  for (Int_t i = 0; i < nHits; i++) {
    THcGoodLADHit *hit = static_cast<THcGoodLADHit *>(hits->At(i));
    if (hit == nullptr)
      continue;

    // Front plane
    if (hit->GetPlaneHit0() >= 0 && hit->GetPlaneHit0() < 999) {
      hit->SetIsProton(0, hit->GetHitEdepMeVHit0() > fFrontPlaneEdepCut);
    }

    // Back plane: piecewise boundary
    //   dt < dt_min or dt > dt_max  -> not proton
    //   dt_min <= dt <= x_diag      -> proton if edep > diagonal (y_diag scaled linearly to 0)
    //   x_diag < dt < dt_max        -> proton if edep > 0
    if (hit->GetPlaneHit1() >= 0 && hit->GetPlaneHit1() < 999) {
      Double_t edep    = hit->GetHitEdepMeVHit1();
      Double_t dt      = hit->GetHitTimeHit1() - hit->GetHitTimeHit0();
      Bool_t is_proton = kFALSE;
      if (dt > fBackPlaneDtMin && dt < fBackPlaneDtMax) {
        Double_t cut_y =
            TMath::Max(0.0, fBackPlaneYDiag * (fBackPlaneXDiag - dt) / (fBackPlaneXDiag - fBackPlaneDtMin));
        is_proton = edep > cut_y;
      }
      hit->SetIsProton(1, is_proton);
    }
  }
}
//_____________________________________________________________________________
void THcLADKine::CalculateTVertex() {
  // Calculate the time of flight
  if (fTrack == nullptr) {
    // Error("CalculateToF", "No track found");
    return;
  }

  Double_t HMScentralPathLen  = 22.0 * 100.;
  Double_t SHMScentralPathLen = 18.1 * 100.;
  Double_t elecMass           = 0.510998 / 1000.0; // electron mass in GeV/c^2
  Double_t lightSpeed         = 29.9792;           // in cm/ns

  Double_t fPtime = fTrack->GetFPTime();
  Double_t xptar  = fTrack->GetTTheta();
  Double_t dP     = fTrack->GetDp();
  Double_t elec_P = fTrack->GetP();
  if (fPtime == -2000 || fPtime == -1000) {
    return;
  }

  // Leave these here for now, but don't think I need them
  // Double_t pSHMS_TdcTime_ROC1 = fTrigDet->Get_CT_Trigtime(0); // SHMS
  // Double_t pHMS_TdcTime_ROC1  = fTrigDet->Get_CT_Trigtime(1); // HMS
  // Double_t pSHMS_TdcTime_ROC2 = fTrigDet->Get_CT_Trigtime(2); // SHMS pTrig1
  // Double_t pHMS_TdcTime_ROC2  = fTrigDet->Get_CT_Trigtime(3); // HMS pTrig3

  char aparatus_prefix[2];
  aparatus_prefix[0] = tolower(fSpecName[0]);
  aparatus_prefix[1] = '\0';

  Double_t DeltaPathLength;
  if (aparatus_prefix[0] == 'h') {
    DeltaPathLength = (.12 * xptar * 1000 + 0.17 * dP / 100.);
    DeltaPathLength += HMScentralPathLen;
  } else if (aparatus_prefix[0] == 'p') {
    DeltaPathLength = (.11 * xptar * 1000 + 0.057 * dP / 100.);
    DeltaPathLength += SHMScentralPathLen;
  } else {
    Error("CalculateToF", "Invalid apparatus prefix: %s", aparatus_prefix);
    return;
  }

  Double_t beta_calc      = elec_P / sqrt(elec_P * elec_P + elecMass * elecMass);
  Double_t PathLengthCorr = DeltaPathLength / (lightSpeed * beta_calc);

  fTVertex = fPtime - PathLengthCorr;

  // Look up the run number; default to a safe value (0 -> first offset row) if unavailable
  Int_t run_number = 0;
  THaVar *varptr   = gHcParms->Find("gen_run_number");
  if (varptr) {
    Int_t *runnum = (Int_t *)varptr->GetValuePointer(); // Assume correct type
    if (runnum)
      run_number = *runnum;
  }

  int fRFSpecIndex = (aparatus_prefix[0] == 'p') ? 0 : 1;
  fRFTime          = fTrigDet->Get_RF_TrigTime(fRFSpecIndex);

  // The RF offset table is stored as n_rf_offsets triplets: {run_start, SHMS, HMS}
  if (n_rf_offsets > 0) {
    int fRFTimeIndex = -1;
    int run_idx      = 0;
    while (fRFTimeIndex < 0) {
      if (run_idx >= 3 * n_rf_offsets || rf_offset[run_idx] > run_number) {
        fRFTimeIndex = run_idx - 3; // Go back 1 step (3 indices) to the applicable row
        break;
      }
      run_idx += 3;
    }
    if (fRFTimeIndex < 0)
      fRFTimeIndex = 0; // run precedes the first offset row -> use the first row
    fRFTimeIndex += (1 + fRFSpecIndex);
    if (fRFTimeIndex >= 3 * n_rf_offsets)
      fRFTimeIndex = 3 * n_rf_offsets - 1; // safety clamp within bounds

    // std::remainder maps the value into [-rf_period/2, +rf_period/2] directly.
    // (fmod would instead return a full-period range with the sign of the dividend.)
    double tmp      = std::remainder(fTVertex - fRFTime + rf_offset[fRFTimeIndex], rf_period);
    fTVertex_RFcorr = fTVertex - tmp;
  } else {
    fTVertex_RFcorr = fTVertex; // no RF offsets available
  }

  // cout << "Offset" << rf_offset[fRFTimeIndex] << " Offset Index " << fRFTimeIndex << " RF Time " << fRFTime
  //      << " TVertex before RF corr " << fTVertex << " TVertex after RF corr " << fTVertex_RFcorr << endl;
  return;
}
//_____________________________________________________________________________
Double_t THcLADKine::CalculateToF(Double_t t_raw) {

  Double_t tof = t_raw - fTVertex + fglobal_time_offset;

  return tof;
}
Double_t THcLADKine::CalculateTOFRFcorr(Double_t t_raw) {

  Double_t tof = t_raw - fTVertex_RFcorr + fglobal_time_offset;

  tof -= vertex.Z() / lightSpeed;

  return tof;
}
//_____________________________________________________________________________
Int_t THcLADKine::DefineVariables(EMode mode) {
  if (mode == kDefine && fIsSetup)
    return kOK;
  fIsSetup = (mode == kDefine);

  // Define variables for the analysis tree
  if (mode == kDefine) {
    RVarDef vars[] = {
        {"t_vertex", "Calculated vertex time", "fTVertex"},
        {"t_vertex_RFcorr", "Calculated vertex time (RF corrected)", "fTVertex_RFcorr"},
        // All below variables are defined and identical to those in ladhodo. Just use H/P.ladhod.var

        // {"n_goodhits", "Number of good hits", "goodhit_n"},
        // {"goodhit_plane_0", "Good hit plane", "fGoodLADHits.THcGoodLADHit.GetPlaneHit0()"},
        // {"goodhit_paddle_0", "Good hit paddle", "fGoodLADHits.THcGoodLADHit.GetPaddleHit0()"},
        // {"goodhit_trackid_0", "Good hit track ID", "fGoodLADHits.THcGoodLADHit.GetTrackIDHit0()"},
        // {"goodhit_beta_0", "Good hit beta", "fGoodLADHits.THcGoodLADHit.GetBetaHit0()"},
        // {"goodhit_dTrkHoriz_0", "Good hit horizontal trk proj - hit position",
        //  "fGoodLADHits.THcGoodLADHit.GetdTrkHorizHit0()"},
        // {"goodhit_dTrkVert_0", "Good hit vertical trk proj - hit position",
        //  "fGoodLADHits.THcGoodLADHit.GetdTrkVertHit0()"},
        // {"goodhit_hittime_0", "Good hit time", "fGoodLADHits.THcGoodLADHit.GetHitTimeHit0()"},
        // {"goodhit_hittheta_0", "Good hit theta", "fGoodLADHits.THcGoodLADHit.GetHitThetaHit0()"},
        // {"goodhit_hitphi_0", "Good hit phi", "fGoodLADHits.THcGoodLADHit.GetHitPhiHit0()"},
        // {"goodhit_hitedep_0", "Good hit energy deposition", "fGoodLADHits.THcGoodLADHit.GetHitEdepHit0()"},
        // {"goodhit_hitedep_amp_0", "Good hit energy deposition (amplitude)",
        //  "fGoodLADHits.THcGoodLADHit.GetHitEdepAmpHit0()"},
        // {"goodhit_plane_1", "Good hit plane (second plane)", "fGoodLADHits.THcGoodLADHit.GetPlaneHit1()"},
        // {"goodhit_paddle_1", "Good hit paddle (second plane)", "fGoodLADHits.THcGoodLADHit.GetPaddleHit1()"},
        // {"goodhit_trackid_1", "Good hit track ID (second plane)", "fGoodLADHits.THcGoodLADHit.GetTrackIDHit1()"},
        // {"goodhit_beta_1", "Good hit beta (second plane)", "fGoodLADHits.THcGoodLADHit.GetBetaHit1()"},
        // {"goodhit_dTrkHoriz_1", "Good hit horizontal trk proj - hit position (second plane)",
        //  "fGoodLADHits.THcGoodLADHit.GetdTrkHorizHit1()"},
        // {"goodhit_dTrkVert_1", "Good hit vertical trk proj - hit position (second plane)",
        //  "fGoodLADHits.THcGoodLADHit.GetdTrkVertHit1()"},
        // {"goodhit_hittime_1", "Good hit time (second plane)", "fGoodLADHits.THcGoodLADHit.GetHitTimeHit1()"},
        // {"goodhit_hittheta_1", "Good hit theta (second plane)", "fGoodLADHits.THcGoodLADHit.GetHitThetaHit1()"},
        // {"goodhit_hitphi_1", "Good hit phi (second plane)", "fGoodLADHits.THcGoodLADHit.GetHitPhiHit1()"},
        // {"goodhit_hitedep_1", "Good hit energy deposition (second plane)",
        //  "fGoodLADHits.THcGoodLADHit.GetHitEdepHit1()"},
        // {"goodhit_hitedep_amp_1", "Good hit energy deposition (amplitude) (second plane)",
        //  "fGoodLADHits.THcGoodLADHit.GetHitEdepAmpHit1()"},
        // {"goodhit_ypos_0", "Good hit y position", "fGoodLADHits.THcGoodLADHit.GetHitYPosHit0()"},
        // {"goodhit_ypos_1", "Good hit y position (second plane)", "fGoodLADHits.THcGoodLADHit.GetHitYPosHit1()"},
        // {"goodhit_tof_0", "Good hit time of flight", "fGoodLADHits.THcGoodLADHit.GetHitTOFHit0()"},
        // {"goodhit_tof_1", "Good hit time of flight (second plane)", "fGoodLADHits.THcGoodLADHit.GetHitTOFHit1()"},
        // {"goodhit_alpha_0", "Good hit alpha", "fGoodLADHits.THcGoodLADHit.GetHitAlphaHit0()"},
        // {"goodhit_alpha_1", "Good hit alpha (second plane)", "fGoodLADHits.THcGoodLADHit.GetHitAlphaHit1()"},

        // 1D-cluster projective tracking (ldo_1Dcluster_tracking) writes its
        // per-hodo-hit chi-squares onto the THcGoodLADHit objects; see the
        // goodhit_chiSquare_1D_* leaves defined by the hodoscope.
        {0}};
    return DefineVarsFromList(vars, mode);
  }
  return kOK;
}
//_____________________________________________________________________________

Double_t THcLADKine::FitTrack(TVector3 vertex, std::vector<TVector3> sp_positions, std::vector<double> sp_resolutions,
                              double dir[3], bool use_y) {
  if (dir == nullptr || sp_positions.empty() || sp_resolutions.empty()) {
    return -1; // Invalid input
  }
  int nPoints = sp_positions.size();
  if (nPoints < 2 || nPoints > 4) {
    return -1; // Not enough points to fit a track
  }
  if (nPoints != sp_resolutions.size()) {
    return -2; // Mismatch in number of points and resolutions
  }
  for (int i = 0; i < nPoints; i++) {
    if (sp_resolutions[i] <= 0) {
      return -3; // Invalid resolution value
    }
  }
  if (dir[0] < 0 || dir[0] > TMath::Pi() || dir[1] < -TMath::Pi() || dir[1] > TMath::Pi()) {
    return -4; // Invalid initial direction values
  }
  // check if the gem track is pointing to the hodoscope hits if there are any, return negative chisq if not
  if (nPoints > 2) {
    TVector3 dir_vec = sp_positions[1] - sp_positions[0];
    dir_vec          = dir_vec.Unit();
    for (int i = 2; i < nPoints; i++) {
      // closest approach of the track to the point
      double t         = ((sp_positions[i].X() - sp_positions[0].X()) * dir_vec.X() +
                          (sp_positions[i].Y() - sp_positions[0].Y()) * dir_vec.Y() +
                          (sp_positions[i].Z() - sp_positions[0].Z()) * dir_vec.Z());
      double x_closest = sp_positions[0].X() + t * dir_vec.X();
      double y_closest = sp_positions[0].Y() + t * dir_vec.Y();
      double z_closest = sp_positions[0].Z() + t * dir_vec.Z();
      double dx        = sp_positions[i].X() - x_closest;
      double dy        = sp_positions[i].Y() - y_closest;
      double dz        = sp_positions[i].Z() - z_closest;
      double dist2     = dx * dx + dz * dz;
      if (dist2 > (22 * 22)) { // if the track is more than 22cm away from the hodoscope hit, return negative chisq
        return -5;             // Track does not point to the hodoscope hit
      }
    }
  }

  // requireing the track to originate from vertex (x,y), and only fitting for the track direction (theta, phi) and z
  // vertex position.
  double chi2 = -kBig;
  if (!fMinimizer)
    fMinimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
  fMinimizer->Clear(); // reset variables/state from the previous fit (reused across calls)
  fMinimizer->SetMaxFunctionCalls(1000);
  fMinimizer->SetTolerance(1e-6);
  ROOT::Math::Minimizer *minimizer = fMinimizer;
  ROOT::Math::Functor f(
      [=](const double *params) {
        double theta      = params[0]; //  in radians
        double phi        = params[1]; //  in radians
        double z          = params[2];
        double chi2_local = 0;
        double sx         = TMath::Sin(theta) * TMath::Cos(phi);
        double sy         = TMath::Sin(theta) * TMath::Sin(phi);
        double sz         = TMath::Cos(theta);
        for (int i = 0; i < nPoints; i++) {
          // closest approach of the track to the point
          double t         = ((sp_positions[i].X() - vertex.X()) * sx + (sp_positions[i].Y() - vertex.Y()) * sy +
                              (sp_positions[i].Z() - z) * sz);
          double x_closest = vertex.X() + t * sx;
          double y_closest = vertex.Y() + t * sy;
          double z_closest = z + t * sz;
          double dx        = sp_positions[i].X() - x_closest;
          double dy        = sp_positions[i].Y() - y_closest;
          double dz        = sp_positions[i].Z() - z_closest;
          double dist2     = dx * dx + dz * dz;
          if (i >= 2) {
            // Assume these are hodoscope hit, so no chisq penalty if the hit is within the width of the paddle
            double paddle_width = 22; // in cm, TODO: get actual paddle width from database
            if (dist2 < (paddle_width / 2.0) * (paddle_width / 2.0)) {
              dist2 = 0;
            }
          }
          // use_y = false drops the (dy)^2 term so only the x-z residual is fit
          chi2_local += (dist2 + (use_y ? dy * dy : 0.0)) /
                        (sp_resolutions[i] * sp_resolutions[i]); // do we want to weight the vertical and horizontal
                                                                 // residuals differently based on detector performance?
        }
        return chi2_local;
      },
      3);
  minimizer->SetFunction(f);

  double initial_params[3] = {dir[0], dir[1], vertex.Z()};
  minimizer->SetLimitedVariable(0, "theta", initial_params[0], 0.1, fThetaMin,
                                fThetaMax); // Limit theta  to avoid unphysical solutions (tracks going backwards)
  minimizer->SetLimitedVariable(
      1, "phi", initial_params[1], 0.1, fPhiMin,
      fPhiMax); // Limit phi as the GEMs are only on the left side of the target in beam direction
  minimizer->SetLimitedVariable(2, "z_vertex", initial_params[2], 0.1, fZCellMin,
                                fZCellMax); // the target length is 20cm
  minimizer->Minimize();
  if (minimizer->Status() != 0) {
    return -6; // Fit did not converge (minimizer is reused, do not delete)
  }
  const double *best_params = minimizer->X();
  dir[0]                    = best_params[0]; // theta in radians
  dir[1]                    = best_params[1]; // phi in radians
  dir[2]                    = best_params[2]; // z vertex position
  chi2                      = minimizer->MinValue();
  // minimizer is reused across calls; do not delete here

  // return the chi2 of the fit
  return chi2;
}
//_____________________________________________________________________________

Double_t THcLADKine::FitTrack_noTrackVertex(std::vector<TVector3> sp_positions, std::vector<double> sp_resolutions,
                                            double dir[3], double anchor[2], bool use_y) {
  // Free 3D-line fit that does NOT use the target vertex.
  //
  // A 3D line has 4 degrees of freedom. Here it is parametrized by its
  // direction (theta, phi) and by the (y, z) coordinates where it crosses a
  // fixed reference plane x = x_ref = sp_positions[0].X() (the GEM1 x). Over the
  // LAD acceptance the along-track x-component sx = sin(theta)cos(phi) is always
  // positive (theta in [60,170] deg, phi in [-50,50] deg), so every candidate
  // line crosses x = x_ref exactly once and the parametrization is non-singular.
  //
  // All space points (GEM1, GEM2 and the hodoscope hit(s)) contribute weighted
  // residuals; unlike the vertex-constrained FitTrack no point is forced exact.
  if (dir == nullptr || anchor == nullptr || sp_positions.empty() || sp_resolutions.empty()) {
    return -1; // Invalid input
  }
  int nPoints = sp_positions.size();
  // Require GEM1, GEM2 + at least one hodoscope hit (3 points); up to GEM1,
  // GEM2, front, back (4 points). With only the 2 GEM points the fit is
  // degenerate (the line is exactly determined), so we reject it here.
  if (nPoints < 3 || nPoints > 4) {
    return -1;
  }
  if (nPoints != (int)sp_resolutions.size()) {
    return -2; // Mismatch in number of points and resolutions
  }
  for (int i = 0; i < nPoints; i++) {
    if (sp_resolutions[i] <= 0) {
      return -3; // Invalid resolution value
    }
  }
  if (dir[0] < 0 || dir[0] > TMath::Pi() || dir[1] < -TMath::Pi() || dir[1] > TMath::Pi()) {
    return -4; // Invalid initial direction values
  }

  const double x_ref = sp_positions[0].X(); // anchor plane at GEM1 x

  // Check that the GEM line (GEM1 -> GEM2) points near the hodoscope hits;
  // if not, reject the track. When use_y is false the check uses only the x-z
  // separation (consistent with the x-z chi-square), otherwise the full 3D one.
  {
    TVector3 dir_vec = (sp_positions[1] - sp_positions[0]).Unit();
    for (int i = 2; i < nPoints; i++) {
      double t         = (sp_positions[i] - sp_positions[0]).Dot(dir_vec);
      TVector3 closest = sp_positions[0] + t * dir_vec;
      double dx        = sp_positions[i].X() - closest.X();
      double dy        = sp_positions[i].Y() - closest.Y();
      double dz        = sp_positions[i].Z() - closest.Z();
      double dist2     = use_y ? (dx * dx + dy * dy + dz * dz) : (dx * dx + dz * dz);
      if (dist2 > (22 * 22)) { // more than 22 cm from the hodoscope hit
        return -5;             // Track does not point to the hodoscope hit
      }
    }
  }

  double chi2 = -kBig;
  if (!fMinimizer)
    fMinimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
  fMinimizer->Clear(); // reset variables/state from the previous fit (reused across calls)
  fMinimizer->SetMaxFunctionCalls(1000);
  fMinimizer->SetTolerance(1e-6);
  ROOT::Math::Minimizer *minimizer = fMinimizer;
  ROOT::Math::Functor f(
      [=](const double *params) {
        double y0         = params[0];
        double z0         = params[1];
        double theta      = params[2]; // in radians
        double phi        = params[3]; // in radians
        double chi2_local = 0;
        double sx         = TMath::Sin(theta) * TMath::Cos(phi);
        double sy         = TMath::Sin(theta) * TMath::Sin(phi);
        double sz         = TMath::Cos(theta);
        for (int i = 0; i < nPoints; i++) {
          // closest approach of the track to the point
          double t =
              ((sp_positions[i].X() - x_ref) * sx + (sp_positions[i].Y() - y0) * sy + (sp_positions[i].Z() - z0) * sz);
          double x_closest = x_ref + t * sx;
          double y_closest = y0 + t * sy;
          double z_closest = z0 + t * sz;
          double dx        = sp_positions[i].X() - x_closest;
          double dy        = sp_positions[i].Y() - y_closest;
          double dz        = sp_positions[i].Z() - z_closest;
          double dist2     = dx * dx + dz * dz;
          if (i >= 2) {
            // Hodoscope hit: no chi-square penalty if within half the paddle width
            double paddle_width = 22; // in cm, TODO: get actual paddle width from database
            if (dist2 < (paddle_width / 2.0) * (paddle_width / 2.0)) {
              dist2 = 0;
            }
          }
          // use_y = false drops the (dy)^2 term so only the x-z residual is fit
          chi2_local += (dist2 + (use_y ? dy * dy : 0.0)) / (sp_resolutions[i] * sp_resolutions[i]);
        }
        return chi2_local;
      },
      4);
  minimizer->SetFunction(f);

  minimizer->SetVariable(0, "y0", anchor[0], 0.1);
  minimizer->SetVariable(1, "z0", anchor[1], 0.1);
  minimizer->SetLimitedVariable(2, "theta", dir[0], 0.1, fThetaMin,
                                fThetaMax); // Limit theta to avoid unphysical solutions (tracks going backwards)
  minimizer->SetLimitedVariable(3, "phi", dir[1], 0.1, fPhiMin,
                                fPhiMax); // Limit phi as the GEMs are only on the left side of the target
  minimizer->Minimize();
  if (minimizer->Status() != 0) {
    return -6; // Fit did not converge (minimizer is reused, do not delete)
  }
  const double *best_params = minimizer->X();
  anchor[0]                 = best_params[0]; // y0
  anchor[1]                 = best_params[1]; // z0
  dir[0]                    = best_params[2]; // theta in radians
  dir[1]                    = best_params[3]; // phi in radians
  chi2                      = minimizer->MinValue();
  // minimizer is reused across calls; do not delete here

  return chi2;
}
//_____________________________________________________________________________
Double_t THcLADKine::FitProj(const TVector3 &vtx, const TVector3 &seed, const TVector3 &tdir,
                             const std::vector<GEM1DMeas> &gem, const std::vector<TVector3> &hodo_pts,
                             double hodo_sigma, bool zero_hodo_xz) {
  // One-parameter projective fit. The track passes through the fixed vertex; its
  // direction is d(a) = (seed + a*tdir).Unit(), i.e. the seed tilted by a single
  // transverse angle 'a' along tdir. Each GEM cluster contributes a plane-crossing
  // residual along its true lab axis; each hodo point contributes the transverse
  // (tdir) distance between the track and the point. A few Gauss-Newton steps
  // handle non-small deflections. Returns chi-square (>= 0), or -1 if there is no
  // residual to fit.
  int nres = (int)gem.size() + (int)hodo_pts.size();
  if (nres < 1 || seed.Mag() < 1e-9 || tdir.Mag() < 1e-9)
    return -1;

  const double eps = 1e-4; // finite-difference step in the transverse angle (rad-ish)

  auto residuals = [&](double a, std::vector<double> &r, std::vector<double> &sig) {
    TVector3 d = (seed + a * tdir).Unit();
    r.clear();
    sig.clear();
    for (const auto &m : gem) {
      double ndotd = m.normal.Dot(d);
      if (fabs(ndotd) < 1e-6) { // track nearly parallel to the plane
        r.push_back(1e3);
        sig.push_back(fSigma_GEM);
        continue;
      }
      double t   = m.normal.Dot(m.origin - vtx) / ndotd;
      TVector3 P = vtx + t * d; // track/plane crossing
      r.push_back(m.axisHat.Dot(P - m.origin) - m.meas);
      sig.push_back(fSigma_GEM);
    }
    for (const auto &H : hodo_pts) {
      double t    = (H - vtx).Dot(d);
      TVector3 Pc = vtx + t * d; // closest point on the track to H
      double res  = tdir.Dot(Pc - H);
      if (zero_hodo_xz && fabs(res) < fHodoXZZeroHW)
        res = 0.0;
      r.push_back(res);
      sig.push_back(hodo_sigma);
    }
  };

  double a = 0.0;
  std::vector<double> r, sig, rp, sp;
  for (int iter = 0; iter < 4; iter++) {
    residuals(a, r, sig);
    residuals(a + eps, rp, sp);
    double num = 0.0, den = 0.0;
    for (size_t i = 0; i < r.size(); i++) {
      double g = (rp[i] - r[i]) / eps; // d r_i / d a
      double w = 1.0 / (sig[i] * sig[i]);
      num += g * r[i] * w;
      den += g * g * w;
    }
    if (den < 1e-12)
      break; // no sensitivity to a (e.g. a single zeroed hodo residual)
    a -= num / den;
  }

  residuals(a, r, sig);
  double chi2 = 0.0;
  for (size_t i = 0; i < r.size(); i++)
    chi2 += (r[i] * r[i]) / (sig[i] * sig[i]);
  return chi2;
}

//_____________________________________________________________________________
void THcLADKine::Do1DClusterTracking() {
  // Projective 1D-cluster tracking. For each hodoscope good hit we require the
  // target vertex and at least one hodo plane hit, then fit target + hodo + GEM
  // clusters in two INDEPENDENT projections (x-z from the V strips, y from the U
  // strips) for three GEM sets: front only (GEM0), back only (GEM1), and both.
  // The nine chi-squares (xz/y/combined x GEM0/GEM1/both) are stored on the good
  // hit. Hypotheses with no usable cluster keep the -1 sentinel, so the many
  // hodo hits that have no real track simply carry -1.
  if (fGEM == nullptr || fHodoscope == nullptr)
    return;

  Int_t nLayers = fGEM->GetNLayers();
  if (nLayers < 2)
    return;
  Int_t layFront = nLayers - 2; // nearer the target (GEM0, z ~ 79 cm)
  Int_t layBack  = nLayers - 1; // farther from the target (GEM1, z ~ 96 cm)

  // Vertex-constrained pass: nothing to do without a vertex.
  if (!(fVertexModule && fVertexModule->HasVertex()))
    return;
  TVector3 vtx = fVertexModule->GetVertex();
  if (fNfixed_z > 0 && fFixed_z) { // snap z to the nearest fixed foil, like the main pass
    double best = 1e30, bz = vtx.Z();
    for (Int_t j = 0; j < fNfixed_z; j++) {
      double dd = fabs(vtx.Z() - fFixed_z[j]);
      if (dd < best) {
        best = dd;
        bz   = fFixed_z[j];
      }
    }
    vtx.SetZ(bz);
  }

  // Collect per-(layer,axis) clusters, capped at fMax1DClusterPerSlot (highest ADC).
  const Int_t maxPerSlot = (fMax1DClusterPerSlot > 0) ? fMax1DClusterPerSlot : std::numeric_limits<Int_t>::max();
  auto collect = [&](Int_t layer, Int_t axis) {
    std::vector<GEM1DMeas> v;
    for (const auto &m : fGEM->Get1DMeas(layer))
      if (m.axis == axis)
        v.push_back(m);
    std::sort(v.begin(), v.end(), [](const GEM1DMeas &a, const GEM1DMeas &b) { return a.adc > b.adc; });
    if ((Int_t)v.size() > maxPerSlot)
      v.resize(maxPerSlot);
    return v;
  };
  std::vector<GEM1DMeas> frontV = collect(layFront, LADGEM::kVaxis); // x-z strips, front
  std::vector<GEM1DMeas> backV  = collect(layBack, LADGEM::kVaxis);  // x-z strips, back
  std::vector<GEM1DMeas> frontU = collect(layFront, LADGEM::kUaxis); // y strips, front
  std::vector<GEM1DMeas> backU  = collect(layBack, LADGEM::kUaxis);  // y strips, back

  const Int_t maxFits = (fMax1DCandidates > 0) ? fMax1DCandidates : std::numeric_limits<Int_t>::max();
  Int_t nFits         = 0;

  // Run one projection (x-z or y) for one GEM set, trying every cluster (or
  // cluster pair for "both") and keeping the smallest chi-square. Returns -1 if
  // no fit could be formed.
  // bestPosA/B (when non-null) receive the lab-frame measurement coordinate of
  // the winning cluster on layer A (front) / layer B (back): the x-z projection
  // reports the lab x, the y projection reports the lab y. Set to -1000 when no
  // fit is formed. The measured lab point of a cluster is origin + meas*axisHat.
  // bestAdcA/B (when non-null) likewise receive the winning cluster's ADC sum.
  auto fitProjection = [&](const std::vector<GEM1DMeas> &clA, const std::vector<GEM1DMeas> &clB, bool both, bool isXZ,
                           const std::vector<TVector3> &hpts, double *bestPosA = nullptr,
                           double *bestPosB = nullptr, double *bestAdcA = nullptr,
                           double *bestAdcB = nullptr) -> double {
    double sigHodo = isXZ ? fSigma_Hodo_x : fSigma_Hodo_y;
    bool zeroXZ    = isXZ; // only the coarse x-z (paddle) hodo residual is optionally zeroed
    double bestchi = -1.0;
    auto coordOf   = [&](const GEM1DMeas &m) {
      TVector3 P = m.origin + m.meas * m.axisHat; // lab point of the cluster along its axis
      return isXZ ? P.X() : P.Y();
    };
    auto tryOne = [&](const std::vector<GEM1DMeas> &gm) {
      if (nFits >= maxFits || gm.empty())
        return;
      nFits++;
      TVector3 meanO(0, 0, 0);
      for (const auto &m : gm)
        meanO += m.origin;
      meanO *= 1.0 / gm.size();
      TVector3 seed = meanO - vtx;
      if (seed.Mag() < 1e-9)
        return;
      seed = seed.Unit();
      TVector3 yhat(0, 1, 0);
      TVector3 ey = yhat - (yhat.Dot(seed)) * seed; // lab-y made transverse to the seed
      if (ey.Mag() < 1e-6)
        return; // seed ~ parallel to lab y (never happens for LAD)
      ey            = ey.Unit();
      TVector3 exz  = seed.Cross(ey).Unit(); // transverse direction in the x-z plane
      TVector3 tdir = isXZ ? exz : ey;
      double chi    = FitProj(vtx, seed, tdir, gm, hpts, sigHodo, zeroXZ);
      if (chi >= 0 && (bestchi < 0 || chi < bestchi)) {
        bestchi = chi;
        // gm[0] is the layer-A (front) cluster; for "both", gm[1] is layer-B (back).
        if (bestPosA)
          *bestPosA = coordOf(gm[0]);
        if (bestPosB)
          *bestPosB = (gm.size() > 1) ? coordOf(gm[1]) : -1000.0;
        if (bestAdcA)
          *bestAdcA = gm[0].adc;
        if (bestAdcB)
          *bestAdcB = (gm.size() > 1) ? gm[1].adc : -1000.0;
      }
    };
    if (bestPosA)
      *bestPosA = -1000.0;
    if (bestPosB)
      *bestPosB = -1000.0;
    if (bestAdcA)
      *bestAdcA = -1000.0;
    if (bestAdcB)
      *bestAdcB = -1000.0;
    if (!both) {
      for (const auto &c : clA)
        tryOne({c});
    } else {
      for (const auto &a : clA)
        for (const auto &b : clB)
          tryOne({a, b});
    }
    return bestchi;
  };

  TClonesArray *ladHits = fHodoscope->GetLADGoodHits();
  Int_t nHodo           = (ladHits ? ladHits->GetLast() + 1 : 0);

  std::vector<GEM1DMeas> empty;
  for (int ih = 0; ih < nHodo; ih++) {
    THcGoodLADHit *gh = static_cast<THcGoodLADHit *>(ladHits->At(ih));
    if (gh == nullptr)
      continue;

    // Hodo point(s): at least one is mandatory; use both if the good hit spans
    // both hodoscope planes.
    std::vector<TVector3> hpts;
    if (gh->GetPlaneHit0() >= 0 && gh->GetPlaneHit0() < 999)
      hpts.push_back(fHodoscope->GetHitPositionLab(gh->GetPlaneHit0(), gh->GetPaddleHit0(), gh->GetHitYPosHit0()));
    if (gh->GetPlaneHit1() >= 0 && gh->GetPlaneHit1() < 999)
      hpts.push_back(fHodoscope->GetHitPositionLab(gh->GetPlaneHit1(), gh->GetPaddleHit1(), gh->GetHitYPosHit1()));
    if (hpts.empty())
      continue; // hodo hit mandatory

    // Winning-cluster lab positions (px/py) and ADC sums (ax/ay) from the
    // single-GEM projections: x from the x-z (V) fit, y from the y (U) fit, per
    // layer. GEM0 = front, GEM1 = back.
    double px0 = -1000., px1 = -1000., py0 = -1000., py1 = -1000.;
    double ax0 = -1000., ax1 = -1000., ay0 = -1000., ay1 = -1000.;
    // x-z projection (V strips)
    double xz0 = frontV.empty() ? -1.0 : fitProjection(frontV, empty, false, true, hpts, &px0, nullptr, &ax0);
    double xz1 = backV.empty() ? -1.0 : fitProjection(backV, empty, false, true, hpts, &px1, nullptr, &ax1);
    double xzB = (frontV.empty() || backV.empty()) ? -1.0 : fitProjection(frontV, backV, true, true, hpts);
    // y projection (U strips)
    double y0 = frontU.empty() ? -1.0 : fitProjection(frontU, empty, false, false, hpts, &py0, nullptr, &ay0);
    double y1 = backU.empty() ? -1.0 : fitProjection(backU, empty, false, false, hpts, &py1, nullptr, &ay1);
    double yB = (frontU.empty() || backU.empty()) ? -1.0 : fitProjection(frontU, backU, true, false, hpts);

    gh->SetTrkChiSqr_1D_xz_GEM0(xz0);
    gh->SetTrkChiSqr_1D_xz_GEM1(xz1);
    gh->SetTrkChiSqr_1D_xz_GEMboth(xzB);
    gh->SetTrkChiSqr_1D_y_GEM0(y0);
    gh->SetTrkChiSqr_1D_y_GEM1(y1);
    gh->SetTrkChiSqr_1D_y_GEMboth(yB);
    // Combined = like-to-like sum, only when both projections exist for that GEM set.
    gh->SetTrkChiSqr_1D_GEM0((xz0 >= 0 && y0 >= 0) ? xz0 + y0 : -1.0);
    gh->SetTrkChiSqr_1D_GEM1((xz1 >= 0 && y1 >= 0) ? xz1 + y1 : -1.0);
    gh->SetTrkChiSqr_1D_GEMboth((xzB >= 0 && yB >= 0) ? xzB + yB : -1.0);
    // Degrees of freedom of each projective fit = (#GEM clusters + #hodo points)
    // - 1 free param. Single GEM -> nH; both GEMs -> nH + 1. The combined (xz+y)
    // dof is the sum of the two projections'. -1 mirrors a missing chi-square.
    const int nH    = (int)hpts.size();
    const int dof_s = nH;     // one GEM cluster + nH hodo - 1
    const int dof_b = nH + 1; // two GEM clusters + nH hodo - 1
    gh->SetTrkNdof_1D_xz_GEM0(xz0 >= 0 ? dof_s : -1);
    gh->SetTrkNdof_1D_xz_GEM1(xz1 >= 0 ? dof_s : -1);
    gh->SetTrkNdof_1D_xz_GEMboth(xzB >= 0 ? dof_b : -1);
    gh->SetTrkNdof_1D_y_GEM0(y0 >= 0 ? dof_s : -1);
    gh->SetTrkNdof_1D_y_GEM1(y1 >= 0 ? dof_s : -1);
    gh->SetTrkNdof_1D_y_GEMboth(yB >= 0 ? dof_b : -1);
    gh->SetTrkNdof_1D_GEM0((xz0 >= 0 && y0 >= 0) ? 2 * dof_s : -1);
    gh->SetTrkNdof_1D_GEM1((xz1 >= 0 && y1 >= 0) ? 2 * dof_s : -1);
    gh->SetTrkNdof_1D_GEMboth((xzB >= 0 && yB >= 0) ? 2 * dof_b : -1);
    // Winning-cluster lab (x,y) per GEM layer (x from the x-z fit, y from the y
    // fit). Parallels the 2D trk.x1/y1 (GEM0) and trk.x2/y2 (GEM1). -1000 = none.
    gh->SetTrk1DX0(px0);
    gh->SetTrk1DY0(py0);
    gh->SetTrk1DX1(px1);
    gh->SetTrk1DY1(py1);
    // ADC sums of those same winning clusters (adcx from the V/x-z cluster, adcy
    // from the U/y cluster), per GEM layer. -1000 = no cluster.
    gh->SetTrk1DAdcX0(ax0);
    gh->SetTrk1DAdcY0(ay0);
    gh->SetTrk1DAdcX1(ax1);
    gh->SetTrk1DAdcY1(ay1);
  }
}
