#include "THcLADHodoscope.h"
#include "THaCutList.h"
#include "THaEvData.h"
#include "THaGlobals.h"
#include "THaSpectrometer.h"
#include "THaTrack.h"
#include "THcDetectorMap.h"
#include "THcGlobals.h"
#include "THcHitList.h"
#include "THcLADSpectrometer.h"
#include "THcParmList.h"
#include "THcSignalHit.h"
#include "VarDef.h"
#include "VarType.h"

//_________________________________________________________________
THcLADHodoscope::THcLADHodoscope(const char *name, const char *description, THaApparatus *apparatus)
    : THaNonTrackingDetector(name, description, apparatus) {
  // Constructor
  fNPlanes = 1;
}

//_________________________________________________________________
THcLADHodoscope::~THcLADHodoscope() {
  // Destructor

  for (int ip = 0; ip < fNPlanes; ip++)
    delete fPlanes[ip];
  delete[] fPlanes;

  delete[] fNPaddle;
  fNPaddle = NULL;
  delete[] fTdcOffset;
  fTdcOffset = NULL;
  delete[] fAdcTdcOffset;
  fAdcTdcOffset = NULL;
  delete[] fHodoSlop;
  fHodoSlop = NULL;

  delete[] fHodoBtmAdcTimeWindowMin;
  fHodoBtmAdcTimeWindowMin = NULL;
  delete[] fHodoBtmAdcTimeWindowMax;
  fHodoBtmAdcTimeWindowMax = NULL;
  delete[] fHodoTopAdcTimeWindowMin;
  fHodoTopAdcTimeWindowMin = NULL;
  delete[] fHodoTopAdcTimeWindowMax;
  fHodoTopAdcTimeWindowMax = NULL;

  delete[] fPlaneNames;
  fPlaneNames = NULL;

  delete[] fHodoVelLight;
  fHodoVelLight = NULL;

  // Time walk
  delete[] fHodoVelFit;
  fHodoVelFit = NULL;
  delete[] fHodoCableFit;
  fHodoCableFit = NULL;
  delete[] fHodo_LCoeff;
  fHodo_LCoeff = NULL;
  delete[] fHodoTop_c1;
  fHodoTop_c1 = NULL;
  delete[] fHodoBtm_c1;
  fHodoBtm_c1 = NULL;
  delete[] fHodoTop_c2;
  fHodoTop_c2 = NULL;
  delete[] fHodoBtm_c2;
  fHodoBtm_c2 = NULL;
}

//_________________________________________________________________
void THcLADHodoscope::Clear(Option_t *opt) {

  /*! \brief Clears variables
   *
   *  Called by  THcHodoscope::Decode
   *
   */
  // fTimeHist_StartTime_Sigma=  kBig;
  // fTimeHist_StartTime_Peak=  kBig;
  // fTimeHist_StartTime_NumPeaks=  0;
  // fTimeHist_StartTime_Hits=  kBig;
  // fTimeHist_FpTime_Sigma=  kBig;
  // fTimeHist_FpTime_Peak=  kBig;
  // fTimeHist_FpTime_NumPeaks=  0;
  // fTimeHist_FpTime_Hits=  kBig;

  // fBeta = 0.0;
  // fBetaNoTrk = 0.0;
  // fBetaNoTrkChiSq = 0.0;
  fStartTime = -1000.;
  // fADCStartTime  = -1000.;
  // fOffsetTime  = kBig;
  fFPTimeAll = -1000.;
  // fGoodStartTime = kFALSE;
  // fGoodScinHits = 0;

  if (*opt != 'I') {
    for (Int_t ip = 0; ip < fNPlanes; ip++) {
      fPlanes[ip]->Clear();
      fFPTime[ip] = 0.;
      // fPlaneCenter[ip]=0.;
      // fPlaneSpacing[ip]=0.;
      // for (UInt_t iPaddle = 0; iPaddle < fNPaddle[ip]; ++iPaddle) {
      //   fScinHitPaddle[ip][iPaddle] = 0;
      // }
    }
  }
  fdEdX.clear();
  fNScinHit.clear();
  // fNClust.clear();
  // fClustSize.clear();
  // fClustPos.clear();
  // fNCluster.clear();
  // fClusterSize.clear();
  // fClusterXPos.clear();
  // fClusterYPos.clear();
  // fThreeScin.clear();
  // fGoodScinHitsX.clear();
  fGoodFlags.clear();

  // Hodo good hit variables
  goodhit_plane.clear();
  goodhit_paddle.clear();
  goodhit_track_id.clear();
  goodhit_beta.clear();
  goodhit_delta_pos_trans.clear();
  goodhit_delta_pos_long.clear();
  goodhit_hit_time.clear();
  goodhit_matching_hit_index.clear();
  goodhit_hit_theta.clear();
  goodhit_hit_phi.clear();
  goodhit_hit_edep.clear();
}
//_________________________________________________________________

void THcLADHodoscope::Setup(const char *name, const char *description) {
  /**
     Create the scintillator plane objects for the hodoscope.

     Uses the Xhodo_num_planes and Xhodo_plane_names to get the number of
     planes and their names.

  */
  if (IsZombie())
    return;

  // fDebug = 1;  // Keep this at one while we're working on the code

  char prefix[2];

  prefix[0] = tolower(GetApparatus()->GetName()[0]);
  // Fix to prevent different param files for SHMS & HMS LAD hodoscopes
  prefix[0] = 'l';
  prefix[1] = '\0';

  TString temp(prefix[0]);

  TString histname = temp + "_timehist";
  hTime            = new TH1F(histname, "", 400, 0, 200);

  string planenamelist;
  DBRequest listextra[] = {{"ladhodo_num_planes", &fNPlanes, kInt},
                           {"ladhodo_plane_names", &planenamelist, kString},
                           {"ladhodo_tdcrefcut", &fTDC_RefTimeCut, kInt, 0, 1},
                           {"ladhodo_adcrefcut", &fADC_RefTimeCut, kInt, 0, 1},
                           {0}};
  // fNPlanes = 4; 		// Default if not defined
  fTDC_RefTimeCut = 0; // Minimum allowed reference times
  fADC_RefTimeCut = 0;
  gHcParms->LoadParmValues((DBRequest *)&listextra, prefix);

  cout << "Plane Name List : " << planenamelist << endl;

  vector<string> plane_names = Podd::vsplit(planenamelist);
  // Plane names
  if (plane_names.size() != (UInt_t)fNPlanes) {
    cout << "ERROR: Number of planes " << fNPlanes << " doesn't agree with number of plane names " << plane_names.size()
         << endl;
    // Should quit.  Is there an official way to quit?
  }

  fPlaneNames = new char *[fNPlanes];
  for (Int_t i = 0; i < fNPlanes; i++) {
    fPlaneNames[i] = new char[plane_names[i].length() + 1];
    strcpy(fPlaneNames[i], plane_names[i].c_str());
  }

  // Probably shouldn't assume that description is defined
  char *desc = new char[strlen(description) + 100];
  fPlanes    = new THcLADHodoPlane *[fNPlanes];
  for (Int_t i = 0; i < fNPlanes; i++) {
    strcpy(desc, description);
    strcat(desc, " Plane ");
    strcat(desc, fPlaneNames[i]);
    fPlanes[i] = new THcLADHodoPlane(fPlaneNames[i], desc, i + 1, this); // Number planes starting from zero!!
    cout << "Created Scintillator Plane " << fPlaneNames[i] << ", " << desc << endl;
  }
  // Save the nominal particle mass
  // This throws an error if it's not a part of the LADspectrometer
  // THcLADSpectrometer *app = dynamic_cast<THcLADSpectrometer *>(GetApparatus());
  // fPartMass               = app->GetParticleMass();
  // fBetaNominal            = app->GetBetaAtPcentral();

  delete[] desc;
}

//_________________________________________________________________
THaAnalysisObject::EStatus THcLADHodoscope::Init(const TDatime &date) {
  Setup(GetName(), GetTitle());

  char EngineDID[] = "xLADSCIN"; // LADSCIN avoids confusion with spectrometer scintillators
  EngineDID[0]     = toupper(GetApparatus()->GetName()[0]);
  if (gHcDetectorMap->FillMap(fDetMap, EngineDID) < 0) {
    static const char *const here = "Init()";
    Error(Here(here), "Error filling detectormap for %s.", EngineDID);
    return kInitError;
  }

  // Should probably put this in ReadDatabase as we will know the
  // maximum number of hits after setting up the detector map
  // But it needs to happen before the sub detectors are initialized
  // so that they can get the pointer to the hitlist.
  cout << " Hodo tdc ref time cut = " << fTDC_RefTimeCut << " " << fADC_RefTimeCut << endl;

  InitHitList(fDetMap, "THcRawHodoHit", fDetMap->GetTotNumChan() + 1, fTDC_RefTimeCut, fADC_RefTimeCut);

  EStatus status;
  if ((status = THaNonTrackingDetector::Init(date)))
    return fStatus = status;

  for (Int_t ip = 0; ip < fNPlanes; ip++) {
    if ((status = fPlanes[ip]->Init(date))) {
      return fStatus = status;
    }
  }

  fNScinHits     = new Int_t[fNPlanes];
  fGoodPlaneTime = new Bool_t[fNPlanes];
  fNPlaneTime    = new Int_t[fNPlanes];
  fSumPlaneTime  = new Double_t[fNPlanes];

  //  Double_t  fHitCnt4 = 0., fHitCnt3 = 0.;

  // Int_t m = 0;
  // fScinHit = new Double_t*[fNPlanes];
  // for ( m = 0; m < fNPlanes; m++ ){
  //   fScinHit[m] = new Double_t[fNPaddle[0]];
  // }

  // for (int ip=0; ip<fNPlanes; ++ip) {
  //   fScinHitPaddle.emplace_back(fNPaddle[ip], 0);
  // }

  fPresentP        = 0;
  THaVar *vpresent = gHaVars->Find(Form("%s.present", GetApparatus()->GetName()));
  if (vpresent) {
    fPresentP = (Bool_t *)vpresent->GetValuePointer();
  }

  return kOK;
}

//_________________________________________________________________
Int_t THcLADHodoscope::End(THaRunBase *run) {
  // Do we really need this function?

  return 0;
}

//_________________________________________________________________
Int_t THcLADHodoscope::DefineVariables(EMode mode) {

  // Initialize variables to read out
  RVarDef vars[] = {{"goodhit_plane", "Good hit plane", "goodhit_plane"},
                    {"goodhit_paddle", "Good hit paddle", "goodhit_paddle"},
                    {"goodhit_track_id", "Good hit track ID", "goodhit_track_id"},
                    {"goodhit_beta", "Good hit beta", "goodhit_beta"},
                    {"goodhit_delta_pos_trans", "Good hit delta position transverse", "goodhit_delta_pos_trans"},
                    {"goodhit_delta_pos_long", "Good hit delta position longitudinal", "goodhit_delta_pos_long"},
                    {"goodhit_hit_time", "Good hit time", "goodhit_hit_time"},
                    {"goodhit_matching_hit_index", "Good hit matching hit index", "goodhit_matching_hit_index"},
                    {"goodhit_hit_theta", "Good hit theta", "goodhit_hit_theta"},
                    {"goodhit_hit_phi", "Good hit phi", "goodhit_hit_phi"},
                    {"goodhit_hit_edep", "Good hit energy deposition", "goodhit_hit_edep"},
                    {0}};

  return DefineVarsFromList(vars, mode);
}

//_________________________________________________________________
Int_t THcLADHodoscope::ReadDatabase(const TDatime &date) {

  cout << "THcLADHodoscope::ReadDatabase()" << endl;
  char prefix[2];
  prefix[0] = tolower(GetApparatus()->GetName()[0]); // "lad"
  // Fix to prevent different param files for SHMS & HMS LAD hodoscopes
  prefix[0] = 'l';
  prefix[1] = '\0';

  // since we define each hodoscope as a separate detector
  // we will need to use the detector name to load parameters
  // for each detector -- to be updated
  fNPaddle = new Int_t[fNPlanes];
  fFPTime  = new Double_t[fNPlanes];
  for (int ip = 0; ip < fNPlanes; ip++) {
    string parname    = "ladhodo_" + string(fPlanes[ip]->GetName()) + "_nr";
    DBRequest list2[] = {{parname.c_str(), &fNPaddle[ip], kInt}, {0}};

    gHcParms->LoadParmValues((DBRequest *)&list2, prefix);
  }

  // for all planes
  fTdcOffset    = new Int_t[fNPlanes];
  fAdcTdcOffset = new Double_t[fNPlanes];
  fHodoSlop     = new Double_t[fNPlanes];
  for (int ip = 0; ip < fNPlanes; ip++) {
    fTdcOffset[ip]    = 0;
    fAdcTdcOffset[ip] = 0;
    fHodoSlop[ip]     = 0;
  }

  // for all elements
  fMaxHodoScin             = fNPaddle[0] * fNPlanes;
  fHodoTopAdcTimeWindowMin = new Double_t[fMaxHodoScin];
  fHodoTopAdcTimeWindowMax = new Double_t[fMaxHodoScin];
  fHodoBtmAdcTimeWindowMin = new Double_t[fMaxHodoScin];
  fHodoBtmAdcTimeWindowMax = new Double_t[fMaxHodoScin];

  fHodoVelLight = new Double_t[fMaxHodoScin];

  // Time walk
  fHodoVelFit   = new Double_t[fMaxHodoScin];
  fHodoCableFit = new Double_t[fMaxHodoScin];
  fHodo_LCoeff  = new Double_t[fMaxHodoScin];
  fHodoTop_c1   = new Double_t[fMaxHodoScin];
  fHodoBtm_c1   = new Double_t[fMaxHodoScin];
  fHodoTop_c2   = new Double_t[fMaxHodoScin];
  fHodoBtm_c2   = new Double_t[fMaxHodoScin];

  for (int ii = 0; ii < fMaxHodoScin; ii++) {
    fHodoTopAdcTimeWindowMin[ii] = -1000.;
    fHodoTopAdcTimeWindowMax[ii] = 1000.;
    fHodoBtmAdcTimeWindowMin[ii] = -1000.;
    fHodoBtmAdcTimeWindowMax[ii] = 1000.;
  }

  DBRequest list3[] = {{"ladcosmicflag", &fCosmicFlag, kInt, 0, 1},
                       {"ladNumPlanesBetaCalc", &fNumPlanesBetaCalc, kInt, 0, 1},
                       {"ladhodo_tdc_min", &fScinTdcMin, kDouble},
                       {"ladhodo_tdc_max", &fScinTdcMax, kDouble},
                       {"ladtof_tolerance", &fTofTolerance, kDouble, 0, 1},
                       {"ladhodo_tdc_offset", fTdcOffset, kInt, (UInt_t)fNPlanes, 1},
                       {"ladhodo_adc_tdc_offset", fAdcTdcOffset, kDouble, (UInt_t)fNPlanes, 1},
                       {"ladhodo_tdc_to_time", &fScinTdcToTime, kDouble},
                       {"ladhodo_TopAdcTimeWindowMin", fHodoTopAdcTimeWindowMin, kDouble, (UInt_t)fMaxHodoScin, 1},
                       {"ladhodo_TopAdcTimeWindowMax", fHodoTopAdcTimeWindowMax, kDouble, (UInt_t)fMaxHodoScin, 1},
                       {"ladhodo_BtmAdcTimeWindowMin", fHodoBtmAdcTimeWindowMin, kDouble, (UInt_t)fMaxHodoScin, 1},
                       {"ladhodo_BtmAdcTimeWindowMax", fHodoBtmAdcTimeWindowMax, kDouble, (UInt_t)fMaxHodoScin, 1},
                       {0}};

  fCosmicFlag        = 0;
  fNumPlanesBetaCalc = 2;
  fTofTolerance      = 3.0;
  fScinTdcMin        = 0;
  fScinTdcMax        = 0;
  fScinTdcToTime     = 0;

  gHcParms->LoadParmValues((DBRequest *)&list3, prefix);

  DBRequest list5[] = {{"is_mc", &fIsMC, kInt, 0, 1}, {0}};
  fIsMC             = 0;
  gHcParms->LoadParmValues((DBRequest *)&list5, "");

  DBRequest list[] = {{"ladhodo_vel_light", &fHodoVelLight[0], kDouble, (UInt_t)fMaxHodoScin, 1}, {0}};
  gHcParms->LoadParmValues((DBRequest *)&list, prefix);

  DBRequest list4[] = {{"ladhodo_velFit", &fHodoVelFit[0], kDouble, (UInt_t)fMaxHodoScin, 1},
                       {"ladhodo_cableFit", &fHodoCableFit[0], kDouble, (UInt_t)fMaxHodoScin, 1},
                       {"ladhodo_LCoeff", &fHodo_LCoeff[0], kDouble, (UInt_t)fMaxHodoScin, 1},
                       {"lad_c1_Top", &fHodoTop_c1[0], kDouble, (UInt_t)fMaxHodoScin, 1},
                       {"lad_c1_Btm", &fHodoBtm_c1[0], kDouble, (UInt_t)fMaxHodoScin, 1},
                       {"lad_c2_Top", &fHodoTop_c2[0], kDouble, (UInt_t)fMaxHodoScin, 1},
                       {"lad_c2_Btm", &fHodoBtm_c2[0], kDouble, (UInt_t)fMaxHodoScin, 1},
                       {0}};

  // fTdc_Thrs = 1.0;
  // Set Default Values if NOT defined in param file
  for (int i = 0; i < fMaxHodoScin; i++) {

    // Turn OFF Time-Walk Correction if param file NOT found
    fHodoTop_c1[i] = 0.0;
    fHodoTop_c2[i] = 0.0;
    fHodoBtm_c1[i] = 0.0;
    fHodoBtm_c2[i] = 0.0;
  }
  for (int i = 0; i < fMaxHodoScin; i++) {
    // Set scin Velocity/Cable to default
    fHodoCableFit[i] = 0.0;
    fHodoVelFit[i]   = 15.0;
    // set time coeff between paddles to default
    fHodo_LCoeff[i] = 0.0;
  }

  gHcParms->LoadParmValues((DBRequest *)&list4, prefix);

  return kOK;
}

//_________________________________________________________________
Int_t THcLADHodoscope::Decode(const THaEvData &evdata) {

  // Decode raw data and pass it to hitlist
  // Read raw data -- THcHitList::DecodeToHitList
  // Processing hitlist for each plane -- THcLADHodoPlane::ProcessHits

  Bool_t present = kTRUE;
  if (fPresentP) {
    present = *fPresentP;
  }

  if (fIsMC) {
    fNSA          = 1;
    fNSB          = 0;
    fNPED         = 1;
    fHaveFADCInfo = true;
  }

  fNHits = DecodeToHitList(evdata, !present);

  // To analyze pedestal events -- Must define "Pedestal_event" cut in the cuts .def file
  // do we want to do this? or calculate pedestal for each event (using the first # of samples, e.g)
  // keeping it for now
  if (gHaCuts->Result("Pedestal_event")) {
    Int_t nexthit = 0;
    for (Int_t ip = 0; ip < fNPlanes; ip++) {
      nexthit = fPlanes[ip]->AccumulatePedestals(fRawHitList, nexthit);
    }
    fAnalyzePedestals = 1; // Analyze pedestals first normal events
    return (0);
  }
  if (fAnalyzePedestals) {
    for (Int_t ip = 0; ip < fNPlanes; ip++) {
      fPlanes[ip]->CalculatePedestals();
    }
    fAnalyzePedestals = 0; // Don't analyze pedestals next event
  }

  Int_t nexthit = 0;
  for (Int_t iplane = 0; iplane < fNPlanes; iplane++) {

    nexthit = fPlanes[iplane]->ProcessHits(fRawHitList, nexthit);
  }

  return fNHits;
}

//_________________________________________________________________
Int_t THcLADHodoscope::CoarseProcess(TClonesArray &tracks) {

  // FIXME. (Potentially) temp fix to prevent error when tracks is nullptr

  fSpectro = dynamic_cast<THaApparatus *>(GetApparatus());
  fGEM     = dynamic_cast<THcLADGEM *>(fSpectro->GetDetector("gem"));
  if (!fGEM) {
    return 0;
  }
  TClonesArray *fGEMTracks = fGEM->GetTracks();

  if (!fGEMTracks) {
    return 0;
  }

  // TClonesArray* gemTracks = spectrometer->GetGEM()->GetTracks();

  Int_t ntracks = fGEMTracks->GetLast() + 1; // Number of reconstructed tracks

  if (ntracks > 0) {
    vector<Double_t> nPmtHit(ntracks);
    vector<Double_t> timeAtFP(ntracks);
    fdEdX.reserve(ntracks);
    fGoodFlags.reserve(ntracks);
    // Loop over all tracks and get corrected time, tof, beta...
    // TODO: initialize fNumPlanesBetaCalc
    for (Int_t itrack = 0; itrack < ntracks; itrack++) { // Line 133
      nPmtHit[itrack]  = 0;
      timeAtFP[itrack] = 0;

      THcLADGEMTrack *theTrack = dynamic_cast<THcLADGEMTrack *>(fGEMTracks->At(itrack));
      if (!theTrack)
        return -1;

      // Calculate Theta, Phi, X and Y (intercepts with z=0 plane) from THcLADGEMTrack
      // FIXME. GEMs currently use m, whereas LAD uses cm
      TVector3 sp1(theTrack->GetX1() * 100, theTrack->GetY1() * 100, theTrack->GetZ1() * 100);
      TVector3 sp2(theTrack->GetX2() * 100, theTrack->GetY2() * 100, theTrack->GetZ2() * 100);

      Double_t track_dz = sp2.Z() - sp1.Z();
      if (track_dz == 0) {
        return 0;
      }
      Double_t track_dx = sp2.X() - sp1.X();
      Double_t track_dy = sp2.Y() - sp1.Y();

      Double_t track_theta = TMath::ATan2(TMath::Sqrt(track_dx * track_dx + track_dy * track_dy), track_dz);
      Double_t track_phi   = TMath::ATan2(track_dy, track_dx);

      Double_t track_x0 = sp1.X() - sp1.Z() * track_dx / track_dz;
      Double_t track_y0 = sp1.Y() - sp1.Z() * track_dy / track_dz;

      for (Int_t ip = 0; ip < fNumPlanesBetaCalc; ip++) {
        fGoodPlaneTime[ip] = kFALSE;
        fNScinHits[ip]     = 0;
        fNPlaneTime[ip]    = 0;
        fSumPlaneTime[ip]  = 0.;
      }

      std::vector<Double_t> dedx_temp;
      std::vector<std::vector<GoodFlags>> goodflagstmp1;
      goodflagstmp1.reserve(fNumPlanesBetaCalc);
#if __cplusplus >= 201103L
      fdEdX.push_back(std::move(dedx_temp)); // Create array of dedx per hit
      fGoodFlags.push_back(std::move(goodflagstmp1));
#else
      fdEdX.push_back(dedx_temp); // Create array of dedx per hit
      fGoodFlags.push_back(goodflagstmp1);
#endif
      Int_t nFPTime      = 0;
      Double_t betaChiSq = -3;
      Double_t beta      = 0;
      //      timeAtFP[itrack] = 0.;
      Double_t sumFPTime = 0.; // Line 138
      fNScinHit.push_back(0);

      hTime->Reset();
      fTOFCalc.clear();  // SAW - Can we
      fTOFPInfo.clear(); // SAW - combine these two?
      Int_t ihhit = 0;   // Hit # overall

      // Loop over all scintillator planes
      for (Int_t ip = 0; ip < fNPlanes; ip++) {

        std::vector<GoodFlags> goodflagstmp2;
        goodflagstmp2.reserve(fNScinHits[ip]);
#if __cplusplus >= 201103L
        fGoodFlags[itrack].push_back(std::move(goodflagstmp2));
#else
        fGoodFlags[itrack].push_back(goodflagstmp2);
#endif
        fNScinHits[ip]         = fPlanes[ip]->GetNScinHits();
        TClonesArray *hodoHits = fPlanes[ip]->GetHits();

        Double_t zPos       = fPlanes[ip]->GetZpos();
        Double_t dzPos      = fPlanes[ip]->GetDzpos();
        Double_t planeTheta = fPlanes[ip]->GetTheta();

        // Loop over hits with in a single plane
        for (Int_t iphit = 0; iphit < fNScinHits[ip]; iphit++) {
          // iphit is hit # within a plane
          THcLADHodoHit *hit = (THcLADHodoHit *)hodoHits->At(iphit);

          fTOFPInfo.emplace_back();

          fTOFPInfo[ihhit].scin_pos_time = 0.0;
          fTOFPInfo[ihhit].scin_neg_time = 0.0;
          fTOFPInfo[ihhit].hit           = hit;
          fTOFPInfo[ihhit].planeIndex    = ip;
          fTOFPInfo[ihhit].hitNumInPlane = iphit;
          fTOFPInfo[ihhit].onTrack       = kFALSE;

          Int_t paddle       = hit->GetPaddleNumber() - 1;
          Double_t zposition = zPos; // TODO: fix this (lad won't have this offset)

          Double_t track_TrnsCoord, track_LongCoord;
          track_TrnsCoord = track_x0 * TMath::Cos(planeTheta - TMath::Pi() / 2);
          track_LongCoord = track_y0;

          Double_t scinTrnsCoord, scinLongCoord;
          // x & y cooridnates are with respect to the central angle pointing to the scintillator plane
          scinTrnsCoord = track_TrnsCoord + TMath::Tan(track_theta - planeTheta) * (zposition); // Line 183

          scinLongCoord =
              track_LongCoord + TMath::Tan(track_phi) / TMath::Cos(track_theta - planeTheta) * (zposition); // Line 184
          // TODO: need to modify getX and getX to get the correct x and y coordinates

          fTOFPInfo[ihhit].scinTrnsCoord = scinTrnsCoord;
          fTOFPInfo[ihhit].scinLongCoord = scinLongCoord;

          Double_t scinCenter = fPlanes[ip]->GetPosCenter(paddle) + fPlanes[ip]->GetPosOffset();

          // Index to access the 2d arrays of paddle/scintillator properties
          Int_t fPIndex = GetScinIndex(ip, paddle);
          Double_t betatrack =
              0; // Todo. FixMe. Set to zero for now.
                 //  theTrack->GetP() / TMath::Sqrt(theTrack->GetP() * theTrack->GetP() + fPartMass * fPartMass);

          if (TMath::Abs(scinCenter - scinTrnsCoord) <
              (fPlanes[ip]->GetSize() * 0.5 + fPlanes[ip]->GetHodoSlop())) { // Line 293

            goodhit_plane.push_back(ip);
            goodhit_paddle.push_back(paddle);
            goodhit_track_id.push_back(itrack);
            goodhit_delta_pos_trans.push_back(scinCenter - scinTrnsCoord);
            goodhit_delta_pos_long.push_back(scinLongCoord - hit->GetCalcPosition());
            goodhit_hit_time.push_back(hit->GetScinCorrectedTime());
            goodhit_hit_theta.push_back(track_theta);
            goodhit_hit_phi.push_back(track_phi);
            goodhit_hit_edep.push_back(TMath::Sqrt(TMath::Max(0., hit->GetTopADC() * hit->GetBtmADC())));

            // Calculate beta
            TVector3 hit_vertex(0, 0, 0);
            double xhit = track_x0 + zposition * TMath::Cos(planeTheta - TMath::Pi() / 2) -
                          scinTrnsCoord * TMath::Sin(planeTheta - TMath::Pi() / 2);
            double yhit = scinLongCoord + track_y0;
            double zhit = -zposition * TMath::Sin(planeTheta - TMath::Pi() / 2) -
                          scinTrnsCoord * TMath::Cos(planeTheta - TMath::Pi() / 2);
            TVector3 hit_momentum(xhit, yhit, zhit);
            TVector3 hit_end(scinTrnsCoord, scinLongCoord, zPos);
            double path_length = (hit_vertex - hit_end).Mag();
            double time        = hit->GetScinCorrectedTime(); // todo. Subtract reference time
            double beta        = path_length / (time * 29.9792458);
            goodhit_beta.push_back(beta);

            // LHE. Everything from here can probably be deleted. It's just left over from the spectrometer hodoscope
            // code (keep z-correction?).

            fTOFPInfo[ihhit].onTrack = kTRUE;
            Double_t zcor =
                zposition / (29.979 * betatrack) * TMath::Sqrt(1. + track_theta * track_theta + track_phi * track_phi);
            fTOFPInfo[ihhit].zcor = zcor;
            if (fCosmicFlag) {
              Double_t zcor =
                  -zposition / (29.979 * 1.0) * TMath::Sqrt(1. + track_theta * track_theta + track_phi * track_phi);
              fTOFPInfo[ihhit].zcor = zcor;
            }
            Double_t tdc_pos = hit->GetTopTDC();
            Double_t tdc_neg = hit->GetBtmTDC();
            //
            if ((tdc_pos >= fScinTdcMin && tdc_pos <= fScinTdcMax) &&
                (tdc_neg >= fScinTdcMin && tdc_neg <= fScinTdcMax)) {
              fTOFPInfo[ihhit].scin_pos_time = hit->GetTopCorrectedTime();
              Double_t timep                 = hit->GetTopCorrectedTime() - zcor;
              fTOFPInfo[ihhit].time_pos      = timep;
              hTime->Fill(timep);
              fTOFPInfo[ihhit].scin_neg_time = hit->GetBtmCorrectedTime();
              Double_t timen                 = hit->GetBtmCorrectedTime() - zcor;
              fTOFPInfo[ihhit].time_neg      = timen;
              hTime->Fill(timen);
            } else {
              // Pretty sure this can be simplified to not recalculate time walk
              //  if (fTrackBetaIncludeSinglePmtHits == 1) {
              //    if (tdc_pos >= fScinTdcMin && tdc_pos <= fScinTdcMax) {
              //      Double_t adc_pos       = hit->GetPosADC();
              //      Double_t adcamp_pos    = hit->GetPosADCpeak();
              //      Double_t pathp         = fPlanes[ip]->GetPosLeft() - scinLongCoord;
              //      fTOFPInfo[ihhit].pathp = pathp;
              //      Double_t timep         = tdc_pos * fScinTdcToTime;

              //     // Double_t tw_corr_pos = fHodoPos_c1[fPIndex]/pow(adcamp_pos/fTdc_Thrs,fHodoPos_c2[fPIndex]) -
              //     // fHodoPos_c1[fPIndex]/pow(200./fTdc_Thrs, fHodoPos_c2[fPIndex]);
              //     Double_t tw_corr_pos = 0.;
              //     pathp                = scinLongCoord;
              //     if (adcamp_pos > 0)
              //       tw_corr_pos = 1. / pow(adcamp_pos / fTdc_Thrs, fHodoPos_c2[fPIndex]) -
              //                     1. / pow(200. / fTdc_Thrs, fHodoPos_c2[fPIndex]);
              //     timep += -tw_corr_pos + fHodo_LCoeff[fPIndex] + pathp / fHodoVelFit[fPIndex];

              //     fTOFPInfo[ihhit].scin_pos_time = timep;
              //     timep -= zcor;
              //     fTOFPInfo[ihhit].time_pos = timep;

              //     hTime->Fill(timep);
              //   }
              //   if (tdc_neg >= fScinTdcMin && tdc_neg <= fScinTdcMax) {
              //     Double_t adc_neg       = hit->GetNegADC();
              //     Double_t adcamp_neg    = hit->GetNegADCpeak();
              //     Double_t pathn         = scinLongCoord - fPlanes[ip]->GetPosRight();
              //     fTOFPInfo[ihhit].pathn = pathn;
              //     Double_t timen         = tdc_neg * fScinTdcToTime;
              //     if (fTofUsingInvAdc) {
              //       timen -= fHodoNegInvAdcOffset[fPIndex] + pathn / fHodoNegInvAdcLinear[fPIndex] +
              //                fHodoNegInvAdcAdc[fPIndex] / TMath::Sqrt(TMath::Max(20.0 * .020, adc_neg));
              //     } else {
              //       pathn                = scinLongCoord;
              //       Double_t tw_corr_neg = 0;
              //       if (adcamp_neg > 0)
              //         tw_corr_neg = 1. / pow(adcamp_neg / fTdc_Thrs, fHodoNeg_c2[fPIndex]) -
              //                       1. / pow(200. / fTdc_Thrs, fHodoNeg_c2[fPIndex]);
              //       timen += -tw_corr_neg - 2 * fHodoCableFit[fPIndex] + fHodo_LCoeff[fPIndex] -
              //                pathn / fHodoVelFit[fPIndex];
              //     }
              //     fTOFPInfo[ihhit].scin_neg_time = timen;
              //     timen -= zcor;
              //     fTOFPInfo[ihhit].time_neg = timen;
              //     hTime->Fill(timen);
              //   }
              // } // new fTrackBetaIncludeSinglePmtHits
            } // matches else
          } // condition for cenetr on a paddle
          ihhit++;
        } // First loop over hits in a plane <---------

        //-----------------------------------------------------------------------------------------------
        //------------- First large loop over scintillator hits ends here --------------------
        //-----------------------------------------------------------------------------------------------
      }
      Int_t nhits = ihhit;

      Double_t TimePeak = DetermineTimePeak(2);
      if (TimePeak > 0) {

        for (Int_t ih = 0; ih < nhits; ih++) { // loop over all scintillator hits
          if ((fTOFPInfo[ih].time_pos > (TimePeak - fTofTolerance)) &&
              (fTOFPInfo[ih].time_pos < (TimePeak + fTofTolerance))) {
            fTOFPInfo[ih].keep_pos = kTRUE;
          }
          if ((fTOFPInfo[ih].time_neg > (TimePeak - fTofTolerance)) &&
              (fTOFPInfo[ih].time_neg < (TimePeak + fTofTolerance))) {
            fTOFPInfo[ih].keep_neg = kTRUE;
          }
        }
      }

      //---------------------------------------------------------------------------------------------
      // ---------------------- Second loop over scint hits in a plane -----------------------------
      //---------------------------------------------------------------------------------------------

      fdEdX[itrack].reserve(nhits);
      fTOFCalc.reserve(nhits);
      for (Int_t ih = 0; ih < nhits; ih++) {
        THcLADHodoHit *hit = fTOFPInfo[ih].hit;
        Int_t iphit        = fTOFPInfo[ih].hitNumInPlane;
        Int_t ip           = fTOFPInfo[ih].planeIndex;
        //         fDumpOut << " looping over hits = " << ih << " plane = " << ip+1 << endl;
        // Flags are used by THcHodoEff
        fGoodFlags[itrack][ip].reserve(nhits);
        fGoodFlags[itrack][ip].push_back(GoodFlags());
        assert(iphit >= 0 && (size_t)iphit < fGoodFlags[itrack][ip].size());
        fGoodFlags[itrack][ip][iphit].onTrack      = kFALSE;
        fGoodFlags[itrack][ip][iphit].goodScinTime = kFALSE;
        fGoodFlags[itrack][ip][iphit].goodTdcNeg   = kFALSE;
        fGoodFlags[itrack][ip][iphit].goodTdcPos   = kFALSE;

        fTOFCalc.emplace_back();
        // Do we set back to false for each track, or just once per event?
        assert(ih >= 0 && (size_t)ih < fTOFCalc.size());
        fTOFCalc[ih].good_scin_time = kFALSE;
        // These need a track index too to calculate efficiencies
        fTOFCalc[ih].good_tdc_pos = kFALSE;
        fTOFCalc[ih].good_tdc_neg = kFALSE;
        fTOFCalc[ih].pindex       = ip;

        Int_t paddle              = hit->GetPaddleNumber() - 1;
        fTOFCalc[ih].hit_paddle   = paddle;
        fTOFCalc[ih].good_raw_pad = paddle;

        //	Double_t scinCenter = fPlanes[ip]->GetPosCenter(paddle) + fPlanes[ip]->GetPosOffset();
        //	Double_t scinTrnsCoord = fTOFPInfo[ih].scinTrnsCoord;
        //	Double_t scinLongCoord = fTOFPInfo[ih].scinLongCoord;

        Int_t fPIndex = GetScinIndex(ip, paddle);

        if (fTOFPInfo[ih].onTrack) {
          fGoodFlags[itrack][ip][iphit].onTrack = kTRUE;
          if (fTOFPInfo[ih].keep_pos) { // 301
            fTOFCalc[ih].good_tdc_pos                = kTRUE;
            fGoodFlags[itrack][ip][iphit].goodTdcPos = kTRUE;
          }
          if (fTOFPInfo[ih].keep_neg) { //
            fTOFCalc[ih].good_tdc_neg                = kTRUE;
            fGoodFlags[itrack][ip][iphit].goodTdcNeg = kTRUE;
          }
          // ** Calculate ave time for scin and error.
          if (fTOFCalc[ih].good_tdc_pos) {
            if (fTOFCalc[ih].good_tdc_neg) {
              fTOFCalc[ih].scin_time    = (fTOFPInfo[ih].scin_pos_time + fTOFPInfo[ih].scin_neg_time) / 2.;
              fTOFCalc[ih].scin_time_fp = (fTOFPInfo[ih].time_pos + fTOFPInfo[ih].time_neg) / 2.;

              fTOFCalc[ih].scin_sigma = 1;
              //  TMath::Sqrt(fHodoSigmaPos[fPIndex] * fHodoSigmaPos[fPIndex] +
              //             fHodoSigmaNeg[fPIndex] * fHodoSigmaNeg[fPIndex]) /
              // 2.; TODO: fix with actual sigma

              fTOFCalc[ih].good_scin_time                = kTRUE;
              fGoodFlags[itrack][ip][iphit].goodScinTime = kTRUE;
            } else {
              fTOFCalc[ih].scin_time    = fTOFPInfo[ih].scin_pos_time;
              fTOFCalc[ih].scin_time_fp = fTOFPInfo[ih].time_pos;

              fTOFCalc[ih].scin_sigma = 1; // fHodoSigmaPos[fPIndex]; TODO: fix with actual sigma

              fTOFCalc[ih].good_scin_time                = kTRUE;
              fGoodFlags[itrack][ip][iphit].goodScinTime = kTRUE;
            }
          } else {
            if (fTOFCalc[ih].good_tdc_neg) {
              fTOFCalc[ih].scin_time    = fTOFPInfo[ih].scin_neg_time;
              fTOFCalc[ih].scin_time_fp = fTOFPInfo[ih].time_neg;

              fTOFCalc[ih].scin_sigma = 1; // fHodoSigmaNeg[fPIndex]; TODO: fix with actual sigma

              fTOFCalc[ih].good_scin_time                = kTRUE;
              fGoodFlags[itrack][ip][iphit].goodScinTime = kTRUE;
            }
          } // In h_tof.f this includes the following if condition for time at focal plane
            // // because it is written in FORTRAN code

          // c     Get time at focal plane
          if (fTOFCalc[ih].good_scin_time) {

            // scin_time_fp doesn't need to be an array
            // Is this any different than the average of time_pos and time_neg?
            //	    Double_t scin_time_fp = ( fTOFPInfo[ih].time_pos +
            //				      fTOFPInfo[ih].time_neg ) / 2.;
            Double_t scin_time_fp = fTOFCalc[ih].scin_time_fp;

            sumFPTime = sumFPTime + scin_time_fp;
            nFPTime++;

            fSumPlaneTime[ip] = fSumPlaneTime[ip] + scin_time_fp;
            fNPlaneTime[ip]++;
            fNScinHit[itrack]++;

            if ((fTOFCalc[ih].good_tdc_pos) && (fTOFCalc[ih].good_tdc_neg)) {
              nPmtHit[itrack] = nPmtHit[itrack] + 2;
            } else {
              nPmtHit[itrack] = nPmtHit[itrack] + 1;
            }

            fdEdX[itrack].push_back(0.0);
            assert(fNScinHit[itrack] > 0 && (size_t)fNScinHit[itrack] < fdEdX[itrack].size() + 1);

            // --------------------------------------------------------------------------------------------
            if (fTOFCalc[ih].good_tdc_pos) {
              if (fTOFCalc[ih].good_tdc_neg) {
                fdEdX[itrack][fNScinHit[itrack] - 1] = TMath::Sqrt(TMath::Max(0., hit->GetTopADC() * hit->GetBtmADC()));
              } else {
                fdEdX[itrack][fNScinHit[itrack] - 1] = TMath::Max(0., hit->GetTopADC());
              }
            } else {
              if (fTOFCalc[ih].good_tdc_neg) {
                fdEdX[itrack][fNScinHit[itrack] - 1] = TMath::Max(0., hit->GetBtmADC());
              } else {
                fdEdX[itrack][fNScinHit[itrack] - 1] = 0.0;
              }
            }
            // --------------------------------------------------------------------------------------------

          } // time at focal plane condition
        } // on track condition

        // ** See if there are any good time measurements in the plane.
        if (fTOFCalc[ih].good_scin_time) {
          fGoodPlaneTime[ip] = kTRUE;
          fTOFCalc[ih].dedx  = fdEdX[itrack][fNScinHit[itrack] - 1];
        } else {
          fTOFCalc[ih].dedx = 0.0;
        }

      } // Second loop over hits of a scintillator plane ends here

      // TODO FIXME. Uncommenting these lines.
      //  theTrack->SetGoodPlane3(fGoodPlaneTime[2] ? 1 : 0);
      //  if (fNumPlanesBetaCalc == 4)
      //    theTrack->SetGoodPlane4(fGoodPlaneTime[3] ? 1 : 0);
      //  End Uncomment

      //
      //------------------------------------------------------------------------------
      //------------------------------------------------------------------------------
      //------------------------------------------------------------------------------
      //------------------------------------------------------------------------------
      //------------------------------------------------------------------------------
      //------------------------------------------------------------------------------
      //------------------------------------------------------------------------------
      //------------------------------------------------------------------------------

      // * * Fit beta if there are enough time measurements (one upper, one lower)
      // From h_tof_fit
      if (((fGoodPlaneTime[0]) || (fGoodPlaneTime[1])) && ((fGoodPlaneTime[2]) || (fGoodPlaneTime[3]))) {

        Double_t sumW  = 0.;
        Double_t sumT  = 0.;
        Double_t sumZ  = 0.;
        Double_t sumZZ = 0.;
        Double_t sumTZ = 0.;

        for (Int_t ih = 0; ih < nhits; ih++) {
          Int_t ip = fTOFPInfo[ih].planeIndex;

          if (fTOFCalc[ih].good_scin_time) {

            Double_t scinWeight = 1 / (fTOFCalc[ih].scin_sigma * fTOFCalc[ih].scin_sigma);
            Double_t zPosition  = (fPlanes[ip]->GetZpos() + (fTOFCalc[ih].hit_paddle % 2) * fPlanes[ip]->GetDzpos());

            sumW += scinWeight;
            sumT += scinWeight * fTOFCalc[ih].scin_time;
            sumZ += scinWeight * zPosition;
            sumZZ += scinWeight * (zPosition * zPosition);
            sumTZ += scinWeight * zPosition * fTOFCalc[ih].scin_time;

          } // condition of good scin time
        } // loop over hits

        Double_t tmp      = sumW * sumZZ - sumZ * sumZ;
        Double_t t0       = (sumT * sumZZ - sumZ * sumTZ) / tmp;
        Double_t tmpDenom = sumW * sumTZ - sumZ * sumT;

        if (TMath::Abs(tmpDenom) > (1 / 10000000000.0)) {

          beta      = tmp / tmpDenom;
          betaChiSq = 0.;

          for (Int_t ih = 0; ih < nhits; ih++) {
            Int_t ip = fTOFPInfo[ih].planeIndex;

            if (fTOFCalc[ih].good_scin_time) {

              Double_t zPosition = (fPlanes[ip]->GetZpos() + (fTOFCalc[ih].hit_paddle % 2) * fPlanes[ip]->GetDzpos());
              Double_t timeDif   = (fTOFCalc[ih].scin_time - t0);
              betaChiSq += ((zPosition / beta - timeDif) * (zPosition / beta - timeDif)) /
                           (fTOFCalc[ih].scin_sigma * fTOFCalc[ih].scin_sigma);

            } // condition for good scin time
          } // loop over hits

          Double_t pathNorm = TMath::Sqrt(1. + track_theta * track_theta + track_phi * track_phi);
          // Take angle into account
          beta = beta / pathNorm;
          beta = beta / 29.979; // velocity / c

        } // condition for fTmpDenom
        else {
          beta      = 0.;
          betaChiSq = -2.;
        } // else condition for fTmpDenom
      } else {
        beta      = 0.;
        betaChiSq = -1;
      }

      if (nFPTime != 0) {
        timeAtFP[itrack] = (sumFPTime / nFPTime);
      }
      //
      // ---------------------------------------------------------------------------

      Double_t FPTimeSum   = 0.0;
      Int_t nFPTimeSum     = 0;
      Int_t nGoodPlanesHit = 0;
      for (Int_t ip = 0; ip < fNumPlanesBetaCalc; ip++) {
        if (fNPlaneTime[ip] != 0) {
          nGoodPlanesHit++;
          fFPTime[ip] = (fSumPlaneTime[ip] / fNPlaneTime[ip]);
          FPTimeSum += fSumPlaneTime[ip];
          nFPTimeSum += fNPlaneTime[ip];
        } else {
          fFPTime[ip] = 1000. * (ip + 1);
        }
      }
      Double_t fptime = -2000;
      fptime          = fStartTime;
      if (nGoodPlanesHit >= 3)
        fptime = FPTimeSum / nFPTimeSum;
      fFPTimeAll    = fptime;
      Double_t dedx = 0.0;
      for (auto &ih : fTOFCalc) {
        if (ih.good_scin_time) {
          dedx = ih.dedx;
          break;
        }
      }

      // TODO: Fixme. Uncommented the following lines. Not currently saving results
      //  theTrack->SetDedx(dedx);
      //  theTrack->SetFPTime(fptime);
      //  theTrack->SetBeta(beta);
      //  theTrack->SetBetaChi2(betaChiSq);
      //  theTrack->SetNPMT(nPmtHit[itrack]);

    } // Main loop over tracks ends here.

  } // If condition for at least one track

  // OriginalTrackEffTest();
  // TrackEffTest();
  // //
  // CalcCluster();

  return 0;
}

//_________________________________________________________________
Int_t THcLADHodoscope::FineProcess(TClonesArray &tracks) { return 0; }

ClassImp(THcLADHodoscope)
