#ifndef ROOT_THcLADKine
#define ROOT_THcLADKine
//////////////////////////////////////////////////////////////////////////

#include "THaSpectrometer.h"
#include "THaTrack.h"
#include "THcGoodLADHit.h"
#include "THcLADGEM.h"
#include "THcLADHodoscope.h"
#include "THcPrimaryKine.h"
#include "THcReactionPoint.h"
#include "THcTrigDet.h"

class THcReactionPoint;
class THcLADHodoscope;
class THcLADGEM;

namespace ROOT {
namespace Math {
class Minimizer;
}
} // namespace ROOT

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
  THcLADGEM *fGEM;
  THcLADHodoscope *fHodoscope;
  TString fVertexModuleName;
  THcReactionPoint *fVertexModule;
  THaTrack *fTrack; // Golden Track
  THcTrigDet *fTrigDet;
  ROOT::Math::Minimizer *fMinimizer; //! reusable Minuit2 minimizer for track fitting (transient)

  Int_t MAXGOODHITS = 500;
  Int_t goodhit_n;
  TClonesArray *fGoodLADHits;
  TClonesArray *fGEMTracks;
  Double_t fD0Cut_wVertex;
  Double_t fD0Cut_noVertex;
  Double_t fMax_dTrk_horiz_match;
  Double_t fMax_dTrk_vert_match;
  Double_t fTrk_dtCut;
  Int_t fNfixed_z;
  Double_t *fFixed_z;
  Double_t fglobal_time_offset;
  Double_t fTVertex;
  Double_t fRFTime;
  Double_t fTVertex_RFcorr;
  Int_t n_rf_offsets;
  Double_t *rf_offset;
  Double_t rf_period;
  TVector3 vertex;
  Double_t lightSpeed = 29.9792458; // Speed of light in cm/ns

  Double_t fZCellMin; // default -15.0 cm
  Double_t fZCellMax; // default +15.0 cm
  Double_t fThetaMin; // default 60.0 deg
  Double_t fThetaMax; // default 170.0 deg
  Double_t fPhiMin;   // default -50.0 deg
  Double_t fPhiMax;   // default +50.0 deg

  Double_t fchisq_cut[2]; // chisq difference between 1 and 2 hodo hit track fits, used to determine if we accept tracks
                          // with only 1 hodo hit (if chisq_2hit - chisq_1hit > fchisq_cut[0]), or if we have no hodo
                          // hits (if chisq_1hit - chisq_0hit < fchisq_cut[1])

  Double_t fFrontPlaneEdepCut; // MeV threshold for front-plane proton ID (default 100)
  Double_t fBackPlaneDtMin;    // dt (ns) of vertical left edge of back-plane proton cut (default 2.8)
  Double_t fBackPlaneDtMax;    // dt (ns) of vertical right edge; hits beyond here are not proton (default 10)
  Double_t fBackPlaneYDiag;    // edep (MeV) at the top of the diagonal segment (default 35)
  Double_t fBackPlaneXDiag;    // dt (ns) where the diagonal meets y=0 (default 5)

  Double_t fSigma_GEM; // GEM resolution in cm, used for track fitting, should be set based on detector performance
  Double_t
      fSigma_Hodo; // Hodoscope resolution in cm, used for track fitting, should be set based on detector performance
  // Hodoscope resolutions used ONLY by the 1D-cluster projective tracking pass
  // (DB keys lsigma_hodo_x / lsigma_hodo_y / lhodo_xz_zero_hw). The x-z (paddle)
  // direction is coarse (~ paddle width / sqrt(12)); the along-paddle y from
  // timing is finer. fHodoXZZeroHW zeroes the x-z hodo residual within this
  // half-width (cm); 0 disables it (default).
  Double_t fSigma_Hodo_x;  // cm, x-z / paddle-direction resolution (default 6.0)
  Double_t fSigma_Hodo_y;  // cm, along-paddle (y) resolution (default 3.0)
  Double_t fHodoXZZeroHW;  // cm, half-width for zeroing the x-z hodo residual (default 0.0 = off)

  // Tracking-mode toggles (DB keys ldo_noVertex_tracking / ldo_xz_tracking /
  // ldo_1Dcluster_tracking). Stored as Int_t (0/1) because the DB loader writes kInt.
  Int_t fDoNoVertexTracking;   // 1 = run the no-vertex tracking pass (default 1)
  Int_t fDoXZTracking;         // 1 = run the x-z (no-y) tracking passes (default 0)
  Int_t fDo1DClusterTracking;  // 1 = run the 1D-cluster (unpaired U/V) tracking pass (default 0)

  // 1D-cluster tracking combinatorics limits (DB keys lmax_1Dcluster_per_slot /
  // lmax_1Dcluster_candidates). PerSlot caps clusters kept per (layer,axis);
  // Candidates caps the total projective fits attempted per event (safety).
  Int_t fMax1DClusterPerSlot;  // max clusters kept per (layer,axis) slot, highest ADC first (default 3)
  Int_t fMax1DCandidates;      // max projective fits attempted per event (default 2000)

  virtual Int_t DefineVariables(EMode mode = kDefine);
  void CalculateTVertex();
  Double_t CalculateToF(Double_t t_raw);
  Double_t CalculateTOFRFcorr(Double_t t_raw);
  // use_y = true  -> full 3D perpendicular residual (x, y, z) as before.
  // use_y = false -> the (dy)^2 term is dropped, so only the x and z residuals
  //                  enter the chi-square ("x-z tracking").
  Double_t FitTrack(TVector3 vertex, std::vector<TVector3> sp_positions, std::vector<double> sp_resolutions,
                    double dir[3], bool use_y = true);
  // No-vertex variant: fit a free 3D line (4 DOF) through the GEM + hodoscope
  // space points, WITHOUT constraining the track to originate at the target
  // vertex. The line is anchored at a fixed x-plane x = sp_positions[0].X()
  // (= GEM1 x), which is single-valued over the LAD acceptance. On input
  // dir[0]=theta seed, dir[1]=phi seed and anchor[0]=y0 seed, anchor[1]=z0 seed
  // (the line's y,z where x = x_ref). On output the fitted values are written
  // back. Returns the chi-square, or a negative code on failure. use_y = false
  // drops the (dy)^2 term so only x and z residuals enter the chi-square.
  Double_t FitTrack_noTrackVertex(std::vector<TVector3> sp_positions, std::vector<double> sp_resolutions, double dir[3],
                                  double anchor[2], bool use_y = true);

  // One-parameter projective fit used by the 1D-cluster tracking pass.
  // The track passes through the fixed vertex; its direction is the seed tilted
  // by a single angle along the transverse direction tdir (the orthogonal
  // transverse angle stays at the seed, giving two independent projective fits
  // -- an x-z fit along the in-plane V axis and a y fit along ~lab-y). Each GEM
  // measurement contributes a plane-crossing residual along its true lab axis;
  // each hodoscope point contributes the transverse (tdir) distance, zeroed
  // within fHodoXZZeroHW for the x-z projection. Returns the chi-square (>= 0),
  // or -1 if under-determined. A few Gauss-Newton iterations handle non-small
  // deflections.
  Double_t FitProj(const TVector3 &vtx, const TVector3 &seed, const TVector3 &tdir,
                   const std::vector<GEM1DMeas> &gem, const std::vector<TVector3> &hodo_pts, double hodo_sigma,
                   bool zero_hodo_xz);
  // Driver for the 1D-cluster projective tracking pass: for each hodoscope good
  // hit, fit target+hodo+{front, back, both} GEM clusters in the x-z and y
  // projections and store the chi-squares on the THcGoodLADHit.
  void Do1DClusterTracking();
  void MakeProtonCut(TClonesArray *hits);

  ClassDef(THcLADKine, 0)
};

#endif /* ROOT_THcLADKine */