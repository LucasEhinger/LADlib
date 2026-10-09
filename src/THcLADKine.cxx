
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
#include <cmath>

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
  fVertexZMax     = 20.0; // cm
  fRejectBadZ     = 1;
  fZTof           = 0.0;
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

  fDoNoVertexTracking = 1; // run the no-vertex tracking pass by default
  fDoXZTracking       = 0; // x-z (no-y) tracking off by default
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
  fZTof           = 0.0;
  vertex.SetXYZ(0, 0, 0);
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
  fVertexZMax          = 20.0; // cm
  fRejectBadZ          = 1;
  fTrk_dtCut           = 10.0;
  fDoNoVertexTracking  = 1; // default: no-vertex tracking on
  fDoXZTracking        = 0; // default: x-z (no-y) tracking off

  cout << "Reading LAD Kinematics parameters from database..." << endl;

  DBRequest list[] = {{"d0_cut_wVertex", &fD0Cut_wVertex, kDouble, 0, 1},
                      {"d0_cut_noVertex", &fD0Cut_noVertex, kDouble, 0, 1},
                      {"nfixed_z", &fNfixed_z, kInt, 0, 1},
                      {"global_time_offset", &fglobal_time_offset, kDouble, 0, 1},
                      {"vertex_zmax", &fVertexZMax, kDouble, 0, 1},
                      {"vertex_bad_z_reject", &fRejectBadZ, kInt, 0, 1},
                      {"trk_dt_cut", &fTrk_dtCut, kDouble, 0, 1},
                      {"do_noVertex_tracking", &fDoNoVertexTracking, kInt, 0, 1},
                      {"do_xz_tracking", &fDoXZTracking, kInt, 0, 1},
                      {"_rf_period", &rf_period, kDouble, 0, 1},
                      {"proton_front_edep_cut", &fFrontPlaneEdepCut, kDouble, 0, 1},
                      {"proton_back_dt_min", &fBackPlaneDtMin, kDouble, 0, 1},
                      {"proton_back_dt_max", &fBackPlaneDtMax, kDouble, 0, 1},
                      {"proton_back_y_diag", &fBackPlaneYDiag, kDouble, 0, 1},
                      {"proton_back_x_diag", &fBackPlaneXDiag, kDouble, 0, 1},
                      {0}};
  gHcParms->LoadParmValues((DBRequest *)&list, prefix);

  // Per-spectrometer LAD ToF offset. The LAD hit time is referenced to the vertex
  // time of this spectrometer, and the SHMS and HMS vertex times differ, so each
  // needs its own offset to put the photon peak at tof = L/c. Optional: when the
  // key for this spectrometer is absent, the shared global_time_offset is used.
  Double_t offset_shms = kBig, offset_hms = kBig;
  DBRequest list_offsets[] = {{"global_time_offset_shms", &offset_shms, kDouble, 0, 1},
                              {"global_time_offset_hms", &offset_hms, kDouble, 0, 1},
                              {0}};
  gHcParms->LoadParmValues((DBRequest *)&list_offsets, prefix);
  const char spec_prefix = tolower(fSpecName[0]);
  const Double_t spec_offset = (spec_prefix == 'p') ? offset_shms : (spec_prefix == 'h') ? offset_hms : kBig;
  if (spec_offset < 0.1 * kBig)
    fglobal_time_offset = spec_offset;
  cout << "  " << GetName() << ": LAD ToF offset = " << fglobal_time_offset << " ns"
       << (spec_offset < 0.1 * kBig ? "" : " (shared global_time_offset)") << endl;

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
    const int cat = hit_index_for_best_chisq / 3, sub = hit_index_for_best_chisq % 3;
    track->SetChisq(bestchisq[cat][sub]);
    track->SetAngles(best_gemdir[cat][sub][0], best_gemdir[cat][sub][1]);
    track->SetProjVertex(vertex.X(), vertex.Y(), best_gemdir[cat][sub][2]);
    Double_t d0 = fabs(vertex.Z() - best_gemdir[cat][sub][2]);
    track->SetD0(d0);

    // Track cuts first, so that a track they reject cannot claim a hodo hit below.
    bool trackOK = track->GetGoodD0() ? (d0 > 0 && d0 < fD0Cut_wVertex) : (d0 > 0 && d0 < fD0Cut_noVertex);
    if (track->GetChisq() < 0 || track->GetdT() > fTrk_dtCut)
      trackOK = false;

    track->SetHasHodoHit(cat);
    track->SetBestHodoHit(nullptr);
    if (cat >= 1) {
      THcGoodLADHit *besthit = static_cast<THcGoodLADHit *>(LADHits_unfiltered->At(bestHodoHitIndex[cat][sub]));
      goodhit_n++;
      if (besthit == nullptr) {
        trackOK = false;
      } else if (goodhit_n > MAXGOODHITS) {
        cout << "THcLADKine::Process: Exceeded maximum number of good hits (" << goodhit_n << "/" << MAXGOODHITS
             << "), skipping hit." << endl;
        trackOK = false;
      } else {
        const int nhodo = (cat == 2 && sub == 0) ? 2 : 1;
        // dof = 3*(#GEM + #hodo) - 3 params (full 3D residual, vertex-constrained)
        const int ndof = 3 * (2 + nhodo) - 3;
        track->SetHasHodoHit(nhodo);
        track->SetBestHodoHit(besthit);
        // Several tracks can pick the same hit as their best; the hit keeps the good one with the lowest
        // chi2/ndof (it used to keep the last one in the loop, including tracks rejected by the cuts above).
        if (trackOK && BetterTrack(track->GetChisq(), ndof, besthit->GetTrkChiSqr(), besthit->GetTrkNdof())) {
          besthit->SetTrackID(track->GetTrackID());
          besthit->SetTrkChiSqr(track->GetChisq());
          besthit->SetTrkNdof(ndof);
        }
      }
    }
    isGoodTrack[iTrack] = trackOK;
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
    bool good     = kTRUE;
    Double_t d0nv = track->GetD0_noTrackVertex();
    if (!(d0nv > 0 && d0nv < fD0Cut_noVertex))
      good = kFALSE;
    if (bestchisq[cat][sub] < 0)
      good = kFALSE;
    if (track->GetdT() > fTrk_dtCut)
      good = kFALSE;
    track->SetIsGoodTrack_noTrackVertex(good);

    // dof = 3*(#GEM + #hodo) - 4 params (full 3D residual, free line, no vertex)
    const int ndof_nv = 3 * (2 + nplanes) - 4;
    // per-hit association from the no-vertex fit: the hit keeps the good track with the lowest chi2/ndof
    if (besthit && good &&
        BetterTrack(track->GetChisq_noTrackVertex(), ndof_nv, besthit->GetTrkChiSqr_noTrackVertex(),
                    besthit->GetTrkNdof_noTrackVertex())) {
      besthit->SetTrackID_noTrackVertex(track->GetTrackID());
      besthit->SetTrkChiSqr_noTrackVertex(track->GetChisq_noTrackVertex());
      besthit->SetTrkNdof_noTrackVertex(ndof_nv);
    }
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
    const int cat = hit_index_for_best_chisq / 3, sub = hit_index_for_best_chisq % 3;
    track->SetChisq_xz(bestchisq[cat][sub]);
    track->SetAngles_xz(best_gemdir[cat][sub][0], best_gemdir[cat][sub][1]);
    track->SetProjVertex_xz(vertex.X(), vertex.Y(), best_gemdir[cat][sub][2]);
    Double_t d0_xz = fabs(vertex.Z() - best_gemdir[cat][sub][2]);
    track->SetD0_xz(d0_xz);

    // Track cuts first, so that a track they reject cannot claim a hodo hit below.
    bool trackOK = track->GetGoodD0() ? (d0_xz > 0 && d0_xz < fD0Cut_wVertex) : (d0_xz > 0 && d0_xz < fD0Cut_noVertex);
    if (track->GetChisq_xz() < 0 || track->GetdT() > fTrk_dtCut)
      trackOK = false;

    track->SetHasHodoHit_xz(cat);
    track->SetBestHodoHit_xz(nullptr);
    if (cat >= 1) {
      THcGoodLADHit *besthit = static_cast<THcGoodLADHit *>(LADHits_unfiltered->At(bestHodoHitIndex[cat][sub]));
      if (besthit == nullptr) {
        trackOK = false;
      } else {
        const int nhodo = (cat == 2 && sub == 0) ? 2 : 1;
        // dof = 2*(#GEM + #hodo) - 3 params (x-z residual only, vertex-constrained)
        const int ndof = 2 * (2 + nhodo) - 3;
        track->SetHasHodoHit_xz(nhodo);
        track->SetBestHodoHit_xz(besthit);
        // the hit keeps the good track with the lowest chi2/ndof (see the vertex-constrained pass)
        if (trackOK && BetterTrack(track->GetChisq_xz(), ndof, besthit->GetTrkChiSqr_xz(), besthit->GetTrkNdof_xz())) {
          besthit->SetTrackID_xz(track->GetTrackID());
          besthit->SetTrkChiSqr_xz(track->GetChisq_xz());
          besthit->SetTrkNdof_xz(ndof);
        }
      }
    }
    isGoodTrack_xz[iTrack] = trackOK;
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
    bool good     = kTRUE;
    Double_t d0nv = track->GetD0_noTrackVertex_xz();
    if (!(d0nv > 0 && d0nv < fD0Cut_noVertex))
      good = kFALSE;
    if (bestchisq[cat][sub] < 0)
      good = kFALSE;
    if (track->GetdT() > fTrk_dtCut)
      good = kFALSE;
    track->SetIsGoodTrack_noTrackVertex_xz(good);

    // dof = 2*(#GEM + #hodo) - 4 params (x-z residual only, free line, no vertex)
    const int ndof_nv = 2 * (2 + nplanes) - 4;
    // per-hit association from the no-vertex x-z fit: the hit keeps the good track with the lowest chi2/ndof
    if (besthit && good &&
        BetterTrack(track->GetChisq_noTrackVertex_xz(), ndof_nv, besthit->GetTrkChiSqr_noTrackVertex_xz(),
                    besthit->GetTrkNdof_noTrackVertex_xz())) {
      besthit->SetTrackID_noTrackVertex_xz(track->GetTrackID());
      besthit->SetTrkChiSqr_noTrackVertex_xz(track->GetChisq_noTrackVertex_xz());
      besthit->SetTrkNdof_noTrackVertex_xz(ndof_nv);
    }
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

    for (int s = 0; s < 2; s++) {
      const Int_t plane = s ? goodhit->GetPlaneHit1() : goodhit->GetPlaneHit0();
      if (plane < 0 || plane >= 999)
        continue;
      const Double_t t     = s ? goodhit->GetHitTimeHit1() : goodhit->GetHitTimeHit0();
      Double_t tof         = CalculateToF(t);
      Double_t tof_rfcorr  = CalculateTOFRFcorr(t);
      goodhit->SetHitTOF(s, tof);
      goodhit->SetHitTOFRFcorr(s, tof_rfcorr);
      // Path from the vertex to the hit (paddle centre, measured y, plane position): no GEM track needed.
      const Int_t paddle = s ? goodhit->GetPaddleHit1() : goodhit->GetPaddleHit0();
      const Double_t y   = s ? goodhit->GetHitYPosHit1() : goodhit->GetHitYPosHit0();
      const Double_t path = (fHodoscope->GetHitPositionLab(plane, paddle, y) - vertex).Mag(); // cm
      goodhit->SetHitPath(s, path);
      if (path > 0 && fabs(tof) < 1e9)
        goodhit->SetHitTOFperM(s, tof / (0.01 * path));
      if (path > 0 && fabs(tof_rfcorr) < 1e9) {
        goodhit->SetHitTOFperMRFcorr(s, tof_rfcorr / (0.01 * path));
        if (tof_rfcorr != 0)
          goodhit->SetHitBeta(s, path / (tof_rfcorr * 29.9792458));
      }
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
  // Electron interaction time at the vertex (fTVertex) and the RF-bucket time at
  // the target centre (fTVertex_RFcorr), from the golden track.

  // Vertex z used by the ToF corrections: this spectrometer's reaction point, or
  // the target centre (z = 0) when there is none. A reaction point outside
  // |z| < vertex_zmax is a junk reconstruction (some reach |z| of metres, i.e.
  // ToF shifts of many ns, z/c = 3.3 ns per metre): with vertex_bad_z_reject = 1
  // the ToF of that event is not computed (fZTof = kBig, ToF = kBig); with 0 the
  // target centre is used. Rejecting gave the narrower photon peak (SHMS: ~9% of
  // LAD hits have such a vertex, and their timing is worse).
  fZTof = 0.0;
  if (fVertexModule && fVertexModule->HasVertex()) {
    const Double_t z = fVertexModule->GetVertex().Z();
    if (fVertexZMax <= 0 || std::fabs(z) < fVertexZMax)
      fZTof = z;
    else if (fRejectBadZ)
      fZTof = kBig;
  }
  const Double_t zcorr = (fZTof < 0.1 * kBig) ? fZTof : 0.0; // z used in the time arithmetic

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

  // The central path length is measured from the target centre. A vertex at z on
  // the beam line is z*cos(theta_e) closer to the spectrometer (theta_e = electron
  // lab angle), so the electron path is shorter by that much. Without this term
  // fTVertex is early by z*cos(theta_e)/c (0.32 ns for z = 10 cm at 13 deg).
  const TVector3 &pv = fTrack->GetPvect();
  const Double_t costh_e =
      (pv.Mag() > 0. && pv.Mag() < 1e10) ? pv.CosTheta() : std::cos(fSpectro ? fSpectro->GetThetaGeo() : 0.);
  DeltaPathLength -= zcorr * costh_e;

  Double_t beta_calc      = elec_P / sqrt(elec_P * elec_P + elecMass * elecMass);
  Double_t PathLengthCorr = DeltaPathLength / (lightSpeed * beta_calc);

  fTVertex = fPtime - PathLengthCorr;

  // Time the electron's bunch passed the target centre: the beam reaches z a time
  // z/c after z = 0. This is what the RF correction snaps to the nearest bucket,
  // so fTVertex_RFcorr is the bucket time at the target centre (CalculateTOFRFcorr
  // adds the z/c back).
  const Double_t t_bunch0 = fTVertex - zcorr / this->lightSpeed;

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
    double tmp      = std::remainder(t_bunch0 - fRFTime + rf_offset[fRFTimeIndex], rf_period);
    fTVertex_RFcorr = t_bunch0 - tmp;
  } else {
    fTVertex_RFcorr = t_bunch0; // no RF offsets available
  }

  // cout << "Offset" << rf_offset[fRFTimeIndex] << " Offset Index " << fRFTimeIndex << " RF Time " << fRFTime
  //      << " TVertex before RF corr " << fTVertex << " TVertex after RF corr " << fTVertex_RFcorr << endl;
  return;
}
//_____________________________________________________________________________
Double_t THcLADKine::CalculateToF(Double_t t_raw) {
  if (fZTof >= 0.1 * kBig)
    return kBig; // vertex outside |z| < vertex_zmax (vertex_bad_z_reject = 1)

  Double_t tof = t_raw - fTVertex + fglobal_time_offset;

  return tof;
}
Double_t THcLADKine::CalculateTOFRFcorr(Double_t t_raw) {
  if (fZTof >= 0.1 * kBig)
    return kBig; // vertex outside |z| < vertex_zmax (vertex_bad_z_reject = 1)

  Double_t tof = t_raw - fTVertex_RFcorr + fglobal_time_offset;

  // fTVertex_RFcorr is the bucket time at the target centre; the interaction at z
  // happened z/c later (fZTof: within |z| < vertex_zmax, or 0 without a vertex).
  tof -= fZTof / lightSpeed;

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
        {"t_vertex", "Electron interaction time at the vertex (path length from the vertex z)", "fTVertex"},
        {"t_vertex_RFcorr", "RF-bucket time at the target centre (RF corrected vertex time)", "fTVertex_RFcorr"},
        {"z_tof", "Vertex z used for the ToF corrections (cm); 0 if no vertex, kBig if outside vertex_zmax and rejected", "fZTof"},
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
