//*-- Author :    Ole Hansen (ole@jlab.org)    9-Dec-2011
//*-- Edited by Lucas Ehinger (ehingerl@jlab.org) for LAD experiment in Hall C. Aug-2024

/////////////////////////////////////////////////////////////////////
//
//   LADSimDecoder
//
//   Decoder for LAD simulation data
//
//   Interprets event buffer from input as LADSimEvent objects
//   (containing digitized simulation data) and unpacks them into
//   crateslot arrays for low-level decoding by detectors.
//
/////////////////////////////////////////////////////////////////////

#include "LADSimDecoder.h"
#include "LADSimDataDecoder.h"
#include "THaBenchmark.h"
#include "THaCrateMap.h"
#include "THaSlotData.h"
#include "VarDef.h"

#include "TDatabasePDG.h"
#include "TError.h"
#include "THaCrateMap.h"
#include "THaDetMap.h"
#include "THaDetector.h"
#include "THaVarList.h"
#include "TMath.h"
#include "TRandom.h"
#include "TSystem.h"
#include "Textvars.h"

#include "TList.h"
#include "TObject.h"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>

using namespace std;
using namespace Podd;

class THaAnalysisObject;

ClassImp(LADSimDecoder) // Implements LADSimDecoder

    //-----------------------------------------------------------------------------
    LADSimDecoder::LADSimDecoder() // : fCheckedForEnabledDetectors(false), fTreeIsSet(false)
{
  // Constructor
  DefineVariables();
  fDetectors.clear();

  gSystem->Load("libEG.so"); // for TDatabasePDG
  // Get MPD encoder for GEMs
  // FIXME: a bit of a kludge...
  // we shouldn't have to do that to initialize all encoders... shall we?
  fDecoderMPD = dynamic_cast<LADSimSADCEncoder *>(LADSimDataDecoder::GetEncoderByName("mpd"));

  fIsInit = false;
}

//-----------------------------------------------------------------------------
LADSimDecoder::~LADSimDecoder() {
  // h1_sizeHCal->Write();
  // h1_sizeGEMs->Write();
  // DefineVariables( THaAnalysisObject::kDelete );
  // h1_sizeHCal->Delete();
  // h1_sizeGEMs->Delete();
}

Int_t LADSimDecoder::Init() {

  Int_t status = THaEvData::Init();

  fDetectors.clear();

  // status += DefineVariables();
  SetDetectors();

  fIsInit = true;

  return status;
}

//-----------------------------------------------------------------------------
Int_t LADSimDecoder::DefineVariables(THaAnalysisObject::EMode mode) {
  // Define global variables for the MC quantities. Extends the base
  // class method.

  const char *const here = "LADSimDecoder::DefineVariables";

  if (mode == THaAnalysisObject::kDefine && fIsSetup)
    return THaAnalysisObject::kOK;

  SimDecoder::DefineVariables(mode);

  cout << "Read LADSimDecoder variables " << endl;

  RVarDef vars[] = {// simc variables
                    // {"simc_sigma", "MC cross section from SIMC gen.", "fSigma_simc"},
                    // {"simc_Weight", "MC cross section weight from SIMC gen.", "fWeight_simc"},
                    // {"simc_Q2", "MC Q2 from SIMC gen.", "fQ2_simc"},
                    // {"simc_xbj", "MC xbj from SIMC gen.", "fXbj_simc"},
                    // {"simc_nu", "MC nu from SIMC gen.", "fNu_simc"},
                    // {"simc_W", "MC W from SIMC gen.", "fW_simc"},
                    // {"simc_epsilon", "MC epsilon from SIMC gen.", "fEpsilon_simc"},
                    // {"simc_Ebeam", "MC Ebeam from SIMC gen.", "fEbeam_simc"},
                    // {"simc_p_e", "MC e momentum from SIMC gen.", "fEp_simc"},
                    // {"simc_theta_e", "MC e polar angle from SIMC gen.", "fEtheta_simc"},
                    // {"simc_phi_e", "MC e azimuthal angle from SIMC gen.", "fEphi_simc"},
                    // {"simc_px_e", "MC e mom. x componant from SIMC gen.", "fEPx_simc"},
                    // {"simc_py_e", "MC e mom. y componant from SIMC gen.", "fEPy_simc"},
                    // {"simc_pz_e", "MC e mom. z componant from SIMC gen.", "fEPz_simc"},
                    // {"simc_fnucl", "MC final-state nucleon type from SIMC gen.", "fFnucl_simc"},
                    // {"simc_p_n", "MC nucleon mom. from SIMC gen.", "fNp_simc"},
                    // {"simc_theta_n", "MC nucleon polar angle from SIMC gen.", "fNtheta_simc"},
                    // {"simc_phi_n", "MC nucleon azimuthal angle from SIMC gen.", "fNphi_simc"},
                    // {"simc_px_n", "MC nucleon mom. x componant from SIMC gen.", "fNPx_simc"},
                    // {"simc_py_n", "MC nucleon mom. y componant from SIMC gen.", "fNPy_simc"},
                    // {"simc_pz_n", "MC nucleon mom. z componant from SIMC gen.", "fNPz_simc"},
                    // {"simc_vx", "MC vertex x co-ordinate from SIMC gen.", "fVx_simc"},
                    // {"simc_vy", "MC vertex y co-ordinate from SIMC gen.", "fVy_simc"},
                    // {"simc_vz", "MC vertex z co-ordinate from SIMC gen.", "fVz_simc"},
                    // {"simc_veE", "MC scattered e- energy at vertex from SIMC gen.", "fVeE_simc"},
                    // {"simc_vetheta", "MC scattered e- theta at vertex from SIMC gen.", "fVetheta_simc"},
                    // ** ^^ **
                    {"mc_sigma", "MC cross section", "fSigma"},
                    {"mc_omega", "MC phase spece generation", "fOmega"},
                    {"mc_epx", "MC electron momentum x", "fEPx"},
                    {"mc_epy", "MC electron momentum y", "fEPy"},
                    {"mc_epz", "MC electron momentum z", "fEPz"},
                    {"mc_npx", "MC nucleon momentum x", "fNPx"},
                    {"mc_npy", "MC nucleon momentum y", "fNPy"},
                    {"mc_npz", "MC nucleon momentum z", "fNPz"},
                    {"mc_vx", "MC vertex x", "fVx"},
                    {"mc_vy", "MC vertex y", "fVy"},
                    {"mc_vz", "MC vertex z", "fVz"},
                    {"mc_ep", "MC Initial momentum of the final state electron in GeV", "fEp"},
                    {"mc_np", "MC Initial momentum of the final state nucleon in GeV", "fNp"},
                    {"mc_nucl", "MC Initial (struck) nucleon type: 1 = proton, 0 = neutron", "fNucl"},
                    {"mc_fnucl", "MC Final-state (detected) nucleon type: 1 = proton, 0 = neutron", "fFnucl"},
                    {"nbbtracks", "number of BB MC tracks", "fNBBtracks"},
                    {"bbtrack_nhits", "BB MC track hit mult", "fBBtrack_Nhits"},
                    {"bbtrack_tid", "BB MC track TID", "fBBtrack_TID"},
                    {"bbtrack_pid", "BB MC track PID", "fBBtrack_PID"},
                    {"bbtrack_mid", "BB MC track MID", "fBBtrack_MID"},
                    {"bbtrack_p", "BB MC track momentum", "fBBtrack_P"},
                    {"bbtrack_x", "BB MC track transport X position", "fBBtrack_X"},
                    {"bbtrack_y", "BB MC track transport Y position", "fBBtrack_Y"},
                    {"bbtrack_dx", "BB MC track transport dX slope", "fBBtrack_dX"},
                    {"bbtrack_dy", "BB MC track transport dY slope", "fBBtrack_dY"},
                    {"nbbgemhits", "number of BBGEM MC hits", "fNBBGEMhits"},
                    {"bbgemhit_plane", "BBGEM MC hit plane", "fBBGEMhit_plane"},
                    {"bbgemhit_tid", "BBGEM MC hit TID", "fBBGEMhit_TID"},
                    {"bbgemhit_pid", "BBGEM MC hit PID", "fBBGEMhit_PID"},
                    {"bbgemhit_mid", "BBGEM MC hit MID", "fBBGEMhit_MID"},
                    {"bbgemhit_edep", "BBGEM MC hit edep", "fBBGEMhit_edep"},
                    {"bbgemhit_x", "BBGEM MC hit transport X", "fBBGEMhit_x"},
                    {"bbgemhit_y", "BBGEM MC hit transport Y", "fBBGEMhit_y"},
                    {"bbps_esum", "BBPS total energy sum", "fBBPS_esum"},
                    {"bbsh_esum", "BBSH total energy sum", "fBBSH_esum"},
                    {"bbgemhit_ptridx", "Primary track index for BBGEM SD", "fBBGEMhit_ptridx"},
                    {"bbgemhit_sdtridx", "SD track index for BBGEM SD", "fBBGEMhit_sdtridx"},
                    {"bbgemtrack_ptridx", "Primary track index for BBGEM Track SD", "fBBGEMtrack_ptridx"},
                    {"bbgemtrack_sdtridx", "SD track index for BBGEM Track SD", "fBBGEMtrack_sdtridx"},
                    {"bbhodohit_ptridx", "Primary track index for BBHodo SD", "fBBHODOhit_ptridx"},
                    {"bbhodohit_sdtridx", "SD track index for BBHodo SD", "fBBHODOhit_sdtridx"},
                    {"bbpshit_ptridx", "Primary track index for BBPSTF1 SD", "fBBPSTF1hit_ptridx"},
                    {"bbpshit_sdtridx", "SD track index for BBPSTF1 SD", "fBBPSTF1hit_sdtridx"},
                    {"bbshhit_ptridx", "Primary track index for BBSHTF1 SD", "fBBSHTF1hit_ptridx"},
                    {"bbshhit_sdtridx", "SD track index for BBSHTF1 SD", "fBBSHTF1hit_sdtridx"},
                    {"hcalhit_ptridx", "Primary track index for HCalScint SD", "fHCALhit_ptridx"},
                    {"hcalhit_sdtridx", "SD track index for HCalScint SD", "fHCALhit_sdtridx"},
                    {"ptrack_ntracks", "Primary track ntracks", "fPTrack_ntacks"},
                    {"ptrack_tid", "Primary track TID", "fPTrack_TID"},
                    {"ptrack_pid", "Primary track PID", "fPTrack_PID"},
                    {"ptrack_posx", "Primary track posx", "fPTrack_posx"},
                    {"ptrack_posy", "Primary track posy", "fPTrack_posy"},
                    {"ptrack_posz", "Primary track posz", "fPTrack_posz"},
                    {"ptrack_momx", "Primary track momx", "fPTrack_momx"},
                    {"ptrack_momy", "Primary track momy", "fPTrack_momy"},
                    {"ptrack_momz", "Primary track momz", "fPTrack_momz"},
                    {"ptrack_polx", "Primary track polx", "fPTrack_polx"},
                    {"ptrack_poly", "Primary track poly", "fPTrack_poly"},
                    {"ptrack_polz", "Primary track polz", "fPTrack_polz"},
                    {"ptrack_etot", "Primary track Etot", "fPTrack_Etot"},
                    {"ptrack_t", "Primary track T", "fPTrack_T"},
                    {"sdtrack_ntracks", "SD track ntracks", "fSDTrack_ntacks"},
                    {"sdtrack_tid", "SD track TID", "fSDTrack_TID"},
                    {"sdtrack_mid", "SD track MID", "fSDTrack_MID"},
                    {"sdtrack_pid", "SD track PID", "fSDTrack_PID"},
                    {"sdtrack_posx", "SD track posx", "fSDTrack_posx"},
                    {"sdtrack_posy", "SD track posy", "fSDTrack_posy"},
                    {"sdtrack_posz", "SD track posz", "fSDTrack_posz"},
                    {"sdtrack_momx", "SD track momx", "fSDTrack_momx"},
                    {"sdtrack_momy", "SD track momy", "fSDTrack_momy"},
                    {"sdtrack_momz", "SD track momz", "fSDTrack_momz"},
                    {"sdtrack_polx", "SD track polx", "fSDTrack_polx"},
                    {"sdtrack_poly", "SD track poly", "fSDTrack_poly"},
                    {"sdtrack_polz", "SD track polz", "fSDTrack_polz"},
                    {"sdtrack_etot", "SD track Etot", "fSDTrack_Etot"},
                    {"sdtrack_t", "SD track T", "fSDTrack_T"},
                    {"sdtrack_vx", "SD track vx", "fSDTrack_vx"},
                    {"sdtrack_vy", "SD track vy", "fSDTrack_vy"},
                    {"sdtrack_vz", "SD track vz", "fSDTrack_vz"},
                    {"sdtrack_vnx", "SD track vnx", "fSDTrack_vnx"},
                    {"sdtrack_vny", "SD track vny", "fSDTrack_vny"},
                    {"sdtrack_vnz", "SD track vnz", "fSDTrack_vnz"},
                    {"sdtrack_vEkin", "SD track vEkin", "fSDTrack_vEkin"},
                    {0}};

  return THaAnalysisObject::DefineVarsFromList(vars, THaAnalysisObject::kRVarDef, mode, "", this, Podd::MC_PREFIX,
                                               here);
}

//-----------------------------------------------------------------------------
void LADSimDecoder::Clear(Option_t *opt) {
  // Clear track and plane data

  SimDecoder::Clear(opt); // clears fMCCherHits, fMCCherClus

  // fPMTMap.clear();
}

//-----------------------------------------------------------------------------
#if ANALYZER_VERSION_CODE >= ANALYZER_VERSION(1, 6, 0)
int LADSimDecoder::LoadEvent(const UInt_t *evbuffer)
#else
int LADSimDecoder::LoadEvent(const Int_t *evbuffer)
#endif
{
  // Wrapper around DoLoadEvent so we can conveniently stop the benchmark
  // counter in case of errors

  if (!fIsInit)
    Init();

  int ret = -1;
  if (sizeof(evbuffer) != 0) {
    ret = DoLoadEvent(evbuffer);
  }

  if (fDoBench)
    fBench->Stop("physics_decode");

  return ret;
}

//-----------------------------------------------------------------------------
#if ANALYZER_VERSION_CODE >= ANALYZER_VERSION(1, 6, 0)
Int_t LADSimDecoder::DoLoadEvent(const UInt_t *evbuffer)
#else
Int_t LADSimDecoder::DoLoadEvent(const Int_t *evbuffer)
#endif
{
  // Fill crateslot structures with Monte Carlo event data in 'evbuffer'
  static const char *const here = "LADSimDecoder::LoadEvent";

#if ANALYZER_VERSION_CODE < ANALYZER_VERSION(1, 6, 0)
  Bool_t fNeedInit = fgNeedInit;
#endif
  assert(fMap || fNeedInit);

  // Local copy of evbuffer pointer, used in GetMCHitInfo
  buffer = evbuffer;

  if (fDebug > 2)
    std::cout << "Processing " << here << std::endl;

  const LADSimEvent *simEvent = reinterpret_cast<const LADSimEvent *>(buffer);
  // add a check here!!!

  // simc variables

  // fSigma_simc   = simEvent->->simc_sigma;
  // fWeight_simc  = simEvent->Tlad->simc_Weight;
  // fQ2_simc      = simEvent->Tlad->simc_Q2;
  // fXbj_simc     = simEvent->Tlad->simc_xbj;
  // fNu_simc      = simEvent->Tlad->simc_nu;
  // fW_simc       = simEvent->Tlad->simc_W;
  // fEpsilon_simc = simEvent->Tlad->simc_epsilon;
  // fEbeam_simc   = simEvent->Tlad->simc_Ebeam;
  // fEp_simc      = simEvent->Tlad->simc_p_e;
  // fEtheta_simc  = simEvent->Tlad->simc_theta_e;
  // fEphi_simc    = simEvent->Tlad->simc_phi_e;
  // fEPx_simc     = simEvent->Tlad->simc_px_e;
  // fEPy_simc     = simEvent->Tlad->simc_py_e;
  // fEPz_simc     = simEvent->Tlad->simc_pz_e;
  // fFnucl_simc   = simEvent->Tlad->simc_fnucl;
  // fNp_simc      = simEvent->Tlad->simc_p_n;
  // fNtheta_simc  = simEvent->Tlad->simc_theta_n;
  // fNphi_simc    = simEvent->Tlad->simc_phi_n;
  // fNPx_simc     = simEvent->Tlad->simc_px_n;
  // fNPy_simc     = simEvent->Tlad->simc_py_n;
  // fNPz_simc     = simEvent->Tlad->simc_pz_n;
  // fVx_simc      = simEvent->Tlad->simc_vx;
  // fVy_simc      = simEvent->Tlad->simc_vy;
  // fVz_simc      = simEvent->Tlad->simc_vz;
  // fVeE_simc     = simEvent->Tlad->simc_veE;
  // fVetheta_simc = simEvent->Tlad->simc_vetheta;
  // g4LAD variables

  // TODO: fixme. Actually load the events we need. Currently commenting out everything except for hodo
  //   fSigma          = simEvent->Tlad->ev_sigma;
  //   fOmega          = simEvent->Tlad->ev_solang;
  //   fEPx            = simEvent->Tlad->ev_epx;
  //   fEPy            = simEvent->Tlad->ev_epy;
  //   fEPz            = simEvent->Tlad->ev_epz;
  //   fNPx            = simEvent->Tlad->ev_npx;
  //   fNPy            = simEvent->Tlad->ev_npy;
  //   fNPz            = simEvent->Tlad->ev_npz;
  //   fVx             = simEvent->Tlad->ev_vx;
  //   fVy             = simEvent->Tlad->ev_vy;
  //   fVz             = simEvent->Tlad->ev_vz;
  //   fEp             = simEvent->Tlad->ev_ep;
  //   fNp             = simEvent->Tlad->ev_np;
  //   fNucl           = simEvent->Tlad->ev_nucl;
  //   fFnucl          = simEvent->Tlad->ev_fnucl;
  //   fNBBtracks      = simEvent->Tlad->Earm_BBGEM_Track_ntracks;
  //   fBBtrack_Nhits  = *(simEvent->Tlad->Earm_BBGEM_Track_NumHits);
  //   fBBtrack_TID    = *(simEvent->Tlad->Earm_BBGEM_Track_TID);
  //   fBBtrack_PID    = *(simEvent->Tlad->Earm_BBGEM_Track_PID);
  //   fBBtrack_MID    = *(simEvent->Tlad->Earm_BBGEM_Track_MID);
  //   fBBtrack_P      = *(simEvent->Tlad->Earm_BBGEM_Track_P);
  //   fBBtrack_X      = *(simEvent->Tlad->Earm_BBGEM_Track_X);
  //   fBBtrack_Y      = *(simEvent->Tlad->Earm_BBGEM_Track_Y);
  //   fBBtrack_dX     = *(simEvent->Tlad->Earm_BBGEM_Track_Xp);
  //   fBBtrack_dY     = *(simEvent->Tlad->Earm_BBGEM_Track_Yp);
  //   fNBBGEMhits     = simEvent->Tlad->Earm_BBGEM_hit_nhits;
  //   fBBGEMhit_plane = *(simEvent->Tlad->Earm_BBGEM_hit_plane);
  //   fBBGEMhit_TID   = *(simEvent->Tlad->Earm_BBGEM_hit_trid);
  //   fBBGEMhit_PID   = *(simEvent->Tlad->Earm_BBGEM_hit_pid);
  //   fBBGEMhit_MID   = *(simEvent->Tlad->Earm_BBGEM_hit_mid);
  //   fBBGEMhit_edep  = *(simEvent->Tlad->Earm_BBGEM_hit_edep);
  //   fBBGEMhit_x     = *(simEvent->Tlad->Earm_BBGEM_hit_tx);
  //   fBBGEMhit_y     = *(simEvent->Tlad->Earm_BBGEM_hit_ty);
  //   // fBBPS_esum          = simEvent->Tlad->Earm_BBPSTF1_det_esum;
  //   // fBBSH_esum          = simEvent->Tlad->Earm_BBSHTF1_det_esum;
  //   fBBGEMhit_ptridx    = *(simEvent->Tlad->Earm_BBGEM_hit_ptridx);
  //   fBBGEMhit_sdtridx   = *(simEvent->Tlad->Earm_BBGEM_hit_sdtridx);
  //   fBBGEMtrack_ptridx  = *(simEvent->Tlad->Earm_BBGEM_Track_ptridx);
  //   fBBGEMtrack_sdtridx = *(simEvent->Tlad->Earm_BBGEM_Track_sdtridx);
  //   fBBHODOhit_ptridx   = *(simEvent->Tlad->Earm_BBHodoScint_hit_ptridx);
  //   fBBHODOhit_sdtridx  = *(simEvent->Tlad->Earm_BBHodoScint_hit_sdtridx);

  Int_t ret = HED_OK;
  if (first_decode || fNeedInit) {
    // fMap->print();
    if ((ret = init_cmap()) != HED_OK)
      return ret;
#if ANALYZER_VERSION_CODE < ANALYZER_VERSION(1, 6, 0)
    if ((ret = init_slotdata(fMap)) != HED_OK)
#else
    if ((ret = init_slotdata()) != HED_OK)
#endif
      return ret;
    first_decode = false;
  }

  if (fDoBench)
    fBench->Begin("clearEvent");
  Clear();
  for (unsigned short i : fSlotClear)
    crateslot[i]->clearEvent();
  if (fDoBench)
    fBench->Stop("clearEvent");

  // FIXME: needed?
  evscaler     = 0;
  event_length = 0;

  // event_type = 1;//not smart to set event_type to 1 automatically...
  event_type = 0; // event_type set to 0 by default
  // only set it to 1 if there is some signal in at least one detector...
  event_num = simEvent->EvtID; //++;
  // int recent_event = event_num; // no longer used

  // Event weight
  // fWeight = simEvent->Tlad->ev_sigma * simEvent->Tlad->ev_solang; //TODO: Actually use the real event weight
  fWeight = 1.0;

  //
  if (fDoBench)
    fBench->Begin("physics_decode");

  // Bool_t newclus;
  // Int_t crate, slot, chan,lchan;

  std::vector<std::map<Decoder::THaSlotData *, std::vector<UInt_t>>> detmaps;
  detmaps.resize(fDetectors.size());

  for (size_t d = 0; d < fDetectors.size(); d++) {
    if (fDebug > 2)
      cout << fDetectors[d] << endl;
    LoadDetector(detmaps[d], fDetectors[d], simEvent);
  }

  // Now call LoadSlot for the different detectors
  for (size_t d = 0; d < fDetectors.size(); d++) {
    // int size_det = 0;
    if (fDebug > 2)
      cout << " " << fDetectors[d] << endl;
    for (std::map<Decoder::THaSlotData *, std::vector<UInt_t>>::iterator it = detmaps[d].begin();
         it != detmaps[d].end(); ++it) {
      if (it->first->GetModule() == 0) {
        if (fDebug > 2) {
          std::cout << "No data available for detector " << fDetectors[d] << std::endl;
        }
      } else {
        event_type = 1;
        // if there is data in at least one detector, event_type set to 1
        if (fDebug > 2) {
          std::cout << "load crate/slot: " << it->first->getCrate() << "/" << it->first->getSlot() << " it->second = {";
          for (size_t k = 0; k < it->second.size(); k++)
            std::cout << it->second[k] << " ; ";
          std::cout << " } " << std::endl;
        }
        it->first->GetModule()->LoadSlot(it->first, it->second.data(), 0, it->second.size());
      }
      // cout << fDetectors[d].c_str() << " " << it->first->getCrate() << " " << it->first->getSlot() << " " <<
      // it->second.size() << endl; size_det+=it->second.size();
    }
    // if(strcmp(fDetectors[d].c_str(), "bb.gem")==0)h1_sizeGEMs->Fill(size_det);
  }

  return HED_OK;
}

// Utilities
/*
Int_t LADSimDecoder::RetrieveDetMapParam(const char* detname,
                                          int& chanperslot, int& slotpercrate,
                                          int& firstcrate, int& firstslot)
{
  // chanperslot = ((TDetInfo &)fManager->GetDetInfo("hcal")).ChanPerSlot();
  // slotpercrate = ((TDetInfo &)fManager->GetDetInfo("hcal")).SlotPerCrate();
  // firstslot = ((TDetInfo &)fManager->GetDetInfo("hcal")).FirstSlot();
  // firstcrate = ((TDetInfo &)fManager->GetDetInfo("hcal")).FirstCrate();
  TDetInfo detinfo = fManager->GetDetInfo(detname);
  chanperslot = detinfo.ChanPerSlot();
  slotpercrate = detinfo.SlotPerCrate();
  firstslot = detinfo.FirstSlot();
  firstcrate = detinfo.FirstCrate();
}
*/

Int_t LADSimDecoder::LoadDetector(std::map<Decoder::THaSlotData *, std::vector<UInt_t>> &map,
                                  const std::string &detname, const LADSimEvent *simev) {
  if (fDebug > 1)
    std::cout << "LADSimDecoder::LoadDectector(" << detname << ")" << std::endl;
  // int detid = detinfo.DetUniqueId();
  Int_t crate, slot;
  // unsigned int nwords = 0;
  unsigned short chan = 0; //, data_type = 0, chan_mult = 0;
  int lchan;
  int mod, apvnum;
  // SimEncoder::mpd_data tmp_mpd;
  // UInt_t* mpd_hdr = new UInt_t[2];
  std::vector<UInt_t> strips;
  std::vector<UInt_t> samps;
  std::vector<UInt_t> times;

  bool loadevt = false;
  // int cur_apv = -1;

  Decoder::THaSlotData *sldat = 0;
  // This should be *general* and work for *every* subsystem
  //  Loop over all raw data in this event
  // UInt_t j = 0;
  // FIXME: we don't want that, I just set it up this way for the sake of going forward
  // Simple fix (might not be ideal): do "if(detname=="xyz")"
  // cout << detname.c_str() << endl;
  int row, col;

  if (strcmp(detname.c_str(), "L.hod") == 0) {
    // cout << " ouh " << detname.c_str() << " " << simev->Tlad->Earm_BBHodoScint_hit_nhits << " " <<
    // simev->Tlad->Earm_BBHodo_dighit_nchan << endl;
    //  cout << simev->Tlad->Earm_BBHodo_dighit_chan->size() << " "
    //  	 << simev->Tlad->Earm_BBHodo_dighit_adc->size() << " "
    //  	 << simev->Tlad->Earm_BBHodo_dighit_tdc_l->size() << " "
    //  	 << simev->Tlad->Earm_BBHodo_dighit_tdc_t->size() << endl;
    /*
    ChanToROC(detname, 180, crate, slot, chan);
    cout << crate << " " << slot << " " << chan << endl;
    if( crate >= 0 || slot >=  0 ) {
      sldat = crateslot[idx(crate,slot)].get();
    }
    std::vector<UInt_t> *myev = &(map[sldat]);
    myev->push_back(LADSimDataDecoder::EncodeHeader(1, chan, 2));
    myev->push_back(0);
    */
    int ntdc = 0;
    assert(simev->Tlad->b_LAD_Hodo_dighit_nchan);
    for (int j = 0; j < simev->Tlad->LAD_Hodo_dighit_nchan; j++) {
      ntdc  = 0;
      lchan = simev->Tlad->LAD_Hodo_dighit_chan->at(j); 
      ChanToROC(detname, lchan, crate, slot, chan);

      if (crate >= 0 || slot >= 0) {
        sldat = crateslot[idx(crate, slot)].get();
      }
      if (simev->Tlad->LAD_Hodo_dighit_tdc_l->at(j) > -1000000)
        ntdc++;
      if (simev->Tlad->LAD_Hodo_dighit_tdc_t->at(j) > -1000000)
        ntdc++;

      if (ntdc) {
        std::vector<UInt_t> *myev = &(map[sldat]);

        // Only save leading tdc here.
        // myev->push_back(LADSimDataDecoder::EncodeHeader(1, chan, ntdc));

        // if (simev->Tlad->LAD_Hodo_dighit_tdc_l->at(j) > -1000000)
        //   myev->push_back(simev->Tlad->LAD_Hodo_dighit_tdc_l->at(j));
        // if (simev->Tlad->LAD_Hodo_dighit_tdc_t->at(j) > -1000000) {
        //   uint tdc = simev->Tlad->LAD_Hodo_dighit_tdc_t->at(j) | (1 << 31);
        //   //   // cout << tdc << endl;
        //   myev->push_back(tdc);
        // }

        myev->push_back(LADSimDataDecoder::EncodeHeader(1, chan, 1));
        myev->push_back(simev->Tlad->LAD_Hodo_dighit_tdc_l->at(j));
        // if (simev->Tlad->LAD_Hodo_dighit_tdc_l->at(j) > -1000000)
        //   myev->push_back(simev->Tlad->LAD_Hodo_dighit_tdc_l->at(j));
        // else
        //   myev->push_back(0);

        ChanToROC(detname, lchan, crate, slot, chan); //+91 ??? that might be the trick
        if (crate >= 0 || slot >= 0) {
          sldat = crateslot[idx(crate, slot + 2)].get();
        }
        myev = &(map[sldat]);

        myev->push_back(LADSimDataDecoder::EncodeHeader(5, chan, 3));
        // if (simev->Tlad->Earm_BBHodo_dighit_adc->at(j) > -1000000) {
        myev->push_back(simev->Tlad->LAD_Hodo_dighit_adc->at(j));
        myev->push_back(simev->Tlad->LAD_Hodo_dighit_adc_amp->at(j));
        myev->push_back(simev->Tlad->LAD_Hodo_dighit_adc_time->at(j));
        // }

        if (fDebug > 2) {
          std::cout << " j = " << j << " my ev = {";
          for (size_t k = 0; k < myev->size(); k++)
            std::cout << myev->at(k) << " ; ";
          std::cout << " } " << std::endl;
        }
      }
    }
  }

  if (strcmp(detname.c_str(), "bb.gem") == 0) {
    // cout << fPx << " " << fPy << " " << fPz << "   " << fVz << endl;
    samps.clear();
    strips.clear();
    // cout << " ouh " << detname.c_str() << " " << simev->Tlad->Earm_BBGEM_dighit_nstrips << endl;
    assert(simev->Tlad->b_Earm_BBGEM_dighit_nstrips);
    for (int j = 0; j < simev->Tlad->Earm_BBGEM_dighit_nstrips; j++) {
      loadevt = false;
      mod     = simev->Tlad->Earm_BBGEM_dighit_module->at(j);
      lchan   = simev->Tlad->Earm_BBGEM_dighit_strip->at(j);
      apvnum  = APVnum(detname, mod, lchan, crate, slot, chan);

      if (simev->Tlad->Earm_BBGEM_dighit_samp->at(j) >= 0) {
        strips.push_back(chan);
        samps.push_back(simev->Tlad->Earm_BBGEM_dighit_adc->at(j));
      }

      if (fDebug > 3)
        cout << " mod " << mod << " lchan " << lchan << " crate " << crate << " slot " << slot << " apvnum " << apvnum
             << " chan " << chan << " samp " << simev->Tlad->Earm_BBGEM_dighit_samp->at(j) << " adc "
             << simev->Tlad->Earm_BBGEM_dighit_adc->at(j) << endl;
      // if(mod>=26 && simev->Tlad->Earm_BBGEM_dighit_samp->at(j)==5)cout << mod << " " << lchan << " " << apvnum <<
      // endl;

      if (j == simev->Tlad->Earm_BBGEM_dighit_nstrips - 1) {
        loadevt = true;
      } else if (mod != simev->Tlad->Earm_BBGEM_dighit_module->at(j + 1) ||
                 // fabs(lchan-simev->Tlad->Earm_BBGEM_dighit_strip->at(j+1))>=128
                 floor(simev->Tlad->Earm_BBGEM_dighit_strip->at(j + 1) / 128) != floor(lchan / 128)) {
        loadevt = true;
      }

      if (loadevt) {
        if (crate >= 0 || slot >= 0) {
          sldat = crateslot[idx(crate, slot)].get();
        }
        std::vector<UInt_t> *myev = &(map[sldat]);

        if (!samps.empty()) {
          // myev->push_back(LADSimDataDecoder::EncodeHeader(5, apvnum, samps.size()));
          // I think I'm onto something here, but I also need to transmit strip num
          myev->push_back(LADSimDataDecoder::EncodeHeader(9, apvnum, samps.size()));
          for (int k = 0; k < (int)samps.size(); k++) {
            // cout << " " << samps[k];
            myev->push_back(strips[k] * 8192 + samps[k]); // strips[k]<< 13 | samps[k]);
          }
          // for(int l = 0; l<myev->size();l++)cout << myev->at(l) << " ";
          // cout << endl;
        }
        // cout << endl;

        samps.clear();
        strips.clear();
      }
    }
  }

  return HED_OK;
}

/*
void LADSimDecoder::SetDetMapParam(const std::string detname, int cps, int spc, int fs, int fc)
{
  fChansPerSlotDetMap[detname] = cps;
  fSlotsPerCrateDetMap[detname] = spc;
  fFirstSlotDetMap[detname] = fs;
  fFirstCrateDetMap[detname] = fc;
}
*/

void LADSimDecoder::CheckForEnabledDetectors() {
  // fDetectors = fManager->GetAllDetInfo();
  if (fDebug > 0) {
    for (size_t i = 0; i < fDetectors.size(); i++) {
      std::cout << "Found detector: " << fDetectors[i].c_str() << endl;
      //<< ", ID: " << fDetectors[i].DetUniqueId() << std::endl;
    }
  }
  fCheckedForEnabledDetectors = true;
}

/*
void LADSimDecoder::SetTree(TTree *t)
{
  if(t==0)return;
  fTree = new digsim_tree(t);
  if(fTree==0)return;
  fTreeIsSet = true;
}
*/

void LADSimDecoder::SetDetectors() {
  // std::cout << "[LADSimDecoder::SetDetectors()]: rundate = ";

  // TDatime rundate = gHaRun->GetDate(); //will this work? answer appears to be NO

  // If the following works, we should be gold:
  // TDatime rundate;
  // rundate.Set(GetRunTime()); // GetRunTime() gives the run time as a UNIX time
  TDatime rundate(124, 1, 1, 0, 0, 0); // FIXME: this is a hack, but it should work for now

  rundate.Print();

  TIter aiter(gHaApps);
  THaApparatus *app = 0;
  while ((app = (THaApparatus *)aiter())) {
    TList *listdet = app->GetDetectors();
    TIter diter(listdet);
    TObject *det = 0;
    while ((det = (TObject *)diter())) {
      cout << "Setting det " << app->GetName() << "." << det->GetName() << " into LADSimDecoder" << endl;
      
      // AddDetector(Form("%s.%s",app->GetName(), det->GetName()),
      // 	    (app->GetDetector(det->GetName()))->GetInitDate());

      // AddDetector(Form("%s.%s", app->GetName(), det->GetName()), rundate);
      string tmp = "L.hod"; //TODO: FIXME. This is hard coded
      AddDetector(tmp, rundate); // FIXME: this is a hack, but it should work for now
    }
  }
}

Int_t LADSimDecoder::AddDetector(string detname, TDatime date) {
  fDetectors.push_back(detname);
  return ReadDetectorDB(detname, date);
}

Int_t LADSimDecoder::ReadDetectorDB(std::string detname, TDatime date) {
  // EPAF: in here the det name is the "full" det name i.e. including the spectro name
  // std::string path = std::string(std::getenv("LAD")) + "/DB/";//TODO: FIXME
  // if (std::getenv("DB_DIR")) {
  //   path = std::string(std::getenv("DB_DIR")) + "/";
  // }
  // const string &fileName = path + "db_" + detname + ".dat";

  // const string fileName = "DBASE/LAD/general.param"; // TODO: FIXME. This is hard coded
  const string fileName = "MAPS/LAD/DETEC/HODO/lhodo_mc.map"; // TODO: FIXME. This is hard coded

  const string prefix = detname + ".";
  // First, open the common db file and parse info there, later, the
  // digitization specific db can be used to override any values
  // FILE* file  = Podd::OpenDBFile(fileName.c_str(), date);

  std::cout << "Calling ReadDetectorDB for detector " << detname << ", Date = ";
  date.Print();

  // FILE *file = Podd::OpenDBFile( detname.c_str(), date, "LADSimDecoder::ReadDetectorDB()",
  // 				 "r", 2 );

  FILE *file = Podd::OpenDBFile(fileName.c_str(), date);
  // FILE *file = Podd::OpenDBFile(detname.c_str(), date);

  if (!file)
    return THaAnalysisObject::kFileError;

  std::vector<int> detmap, chanmap, detmap_adc;
  uint nchan, nlogchan = 0, chanmapstart = 0;

  // int cps, spc, fs, fc;

  bool isgem  = (detname.find("gem") != std::string::npos);
  int apv_num = -1, mpd = -1, mod = 0, axis = -1;
  // int pos = -1;

  DBRequest request[] = {{"nchan", &nchan, kInt, 0, false},               //
                         {"nlog_chan", &nlogchan, kInt, 0, true},         // <- optional
                         {"detmap", &detmap, kIntV, 0, false},            //
                         {"chanmap", &chanmap, kIntV, 0, true},           // <- optional
                         {"chanmap_start", &chanmapstart, kInt, 0, true}, // <- optional
                         {"detmap_adc", &detmap_adc, kIntV, 0, true},     // <- optional
                         /*
                         {"first_crate", &fc, kInt, 0, true},// <- optional
                         {"first_slot", &fs, kInt, 0, true},//  <- optional
                         {"chan_per_slot", &cps, kInt, 0, true},//  <- optional
                         {"slot_per_crate", &spc, kInt, 0, true},//  <- optional
                         */
                         {0}};
  Int_t err;
  int nparam_mod = 5;
  if (isgem) { // gem detectors
    nparam_mod = 9;
  }
  int crate, slot, ch_lo, ch_hi, ch_ref, ch_count = 0, ch_map = 0;

  if (isgem) { // it's easier if gems are their own thing
    /*
    std::string chambers;
    DBRequest req_chambers[] = {
      {"chambers", &chambers, kString, 0, false}, //
      { 0 }
    };
    err = THaAnalysisObject::LoadDB(file, date, req_chambers, prefix.c_str());

    //cout << " prefix " << prefix.c_str() << " err " << err << " chambers " << chambers.c_str() << " size ? " <<
    chambers.size() << endl;

    std::vector<std::string> chambers_names;
    if(err==0)chambers_names = vsplit(chambers);

    if(!chambers_names.empty()){
      for (std::vector<std::string>::iterator it = chambers_names.begin() ; it != chambers_names.end(); ++it){
    */
    std::string modules;
    // std::string pref_cham = prefix+(*it)+".";
    // std::string pref_cham = prefix;//+".";
    // cout << "prefix chamber "  << pref_cham.c_str() << endl;
    DBRequest req_modules[] = {{"modules", &modules, kString, 0, false}, //
                               {0}};
    err                     = THaAnalysisObject::LoadDB(file, date, req_modules, prefix.c_str());

    if (err)
      return THaAnalysisObject::kInitError;

    // cout << " prefix " << pref_cham.c_str() << " err " << err << " modules " << modules.c_str() << " size ? " <<
    // modules.size() << endl;

    std::vector<std::string> modules_names;
    if (err == 0)
      modules_names = vsplit(modules);
    if (!modules_names.empty()) {
      for (std::vector<std::string>::iterator jt = modules_names.begin(); jt != modules_names.end(); ++jt) {
        std::string pref_mod = prefix + (*jt) + ".";

        DBRequest request_gem[] = {{"chanmap", &chanmap, kIntV, 0, false}, {0}};
        err += THaAnalysisObject::LoadDB(file, date, request_gem, pref_mod.c_str());

        fInvGEMDetMap[detname].resize(
            fInvGEMDetMap[detname].size() +
            2); // increments the size of this container by two. But it never gets initialized prior to now! We have to
                // trust that the size is zero to start with. Is this safe? Probably not...

        // This resizes the fInvGEMDetMap[detname][2*module+axis] to the total size of the decode map
        for (int m = 0; m < 2; m++)
          (fInvGEMDetMap[detname])[mod * 2 + m].resize(chanmap.size() / nparam_mod);

        // std::cout << "(detname, mod, nparam_mod)=(" << detname << ", " << mod << ", " << nparam_mod
        // 	  << ")" << std::endl;

        //	int nparam_mod = 9;
        int ax_prev = 0;
        int n_ax = 0, n_ax_x = 0, n_ax_y = 0;
        for (size_t k = 0; k < chanmap.size(); k += nparam_mod) {
          // for(int m = 0; m<nparam_mod; m++)std::cout << chanmap[k+m] << " ";
          // std::cout << std::endl;
          crate   = chanmap[k];
          slot    = chanmap[k + 1];
          mpd     = chanmap[k + 2];
          apv_num = mpd << 4 | chanmap[k + 4]; //
          // pos = chanmap[k+6];
          axis = chanmap[k + 8];
          if (axis == 0)
            n_ax_x++;
          if (axis == 1)
            n_ax_y++;
          if (ax_prev != axis) {
            n_ax    = 0;
            ax_prev = axis;
          }
          ch_lo = 128 * n_ax;
          ch_hi = 128 * (n_ax + 1) - 1;
          // mod*2+axis???
          // std::cout << mod << " " << mod*2+axis << " " << fInvGEMDetMap[detname].size() << " " << mpd << " " <<
          // chanmap[k+4] << " " << apv_num << " " << n_ax << endl;
          (fInvGEMDetMap[detname])[mod * 2 + axis][n_ax] = gemstripinfo(crate, slot, apv_num);
          n_ax++;
        }
        (fInvGEMDetMap[detname])[mod * 2 + 0].resize(n_ax_x);
        (fInvGEMDetMap[detname])[mod * 2 + 1].resize(n_ax_y);
        /*
        std::string planeconfig;
        //cout << "prefix module "  << pref_mod.c_str() << endl;
        DBRequest req_planeconfig[] = {
          {"planeconfig", &planeconfig, kString, 0, false}, //
          { 0 }
        };
        err+= THaAnalysisObject::LoadDB(file, date, req_planeconfig, pref_mod.c_str());
        //cout << " prefix " << pref_mod.c_str() << " err " << err << " planeconfig " << planeconfig.c_str() << " size ?
        " << planeconfig.size() << endl; std::vector<std::string> plane_readouts; if(err==0)plane_readouts =
        vsplit(planeconfig); if(!plane_readouts.empty()){ for (std::vector<std::string>::iterator kt =
        plane_readouts.begin() ; kt != plane_readouts.end(); ++kt){
            //mod++;
            std::string pref_ro = pref_mod+(*kt)+".";
            ch_count = 0;
            fInvGEMDetMap[detname].resize(fInvGEMDetMap[detname].size()+1);

            if(fDebug>=2)cout << fInvGEMDetMap[detname].size() << " module number " << mod << endl;

            err+= THaAnalysisObject::LoadDB(file, date, request, pref_ro.c_str());
            //if(nlogchan==0)
            nlogchan = nchan;
            if(fDebug>=2)cout << " prefix " << pref_ro.c_str() << " err " << err << endl;

            if(err==0)fInvGEMDetMap[detname][mod].resize(nchan);

            //int nparam_mod = 5;
            for(size_t k = 0; k < detmap.size(); k+=nparam_mod) {
              crate  = detmap[k];
              slot   = detmap[k+1];
              ch_lo  = detmap[k+2];
              ch_hi  = detmap[k+3];

              for(int i = ch_lo; i<=ch_hi; i++, ch_count++){
                if(i%128==0){
                  apv_num++;
                  if(fDebug>=3)cout << crate << " " << slot << " " << i << " " << apv_num << endl;
                }
                if(ch_count>nlogchan){
                  std::cout << " <0> number of channels defined in detmap ( >= " << ch_count << ") exceeds logical
        number of channels = " << nlogchan << std::endl; return THaAnalysisObject::kInitError;
                }
                (fInvGEMDetMap[detname])[mod][ch_count]=gemstripinfo(crate, slot, i, apv_num);
              }

            }
          }//end loop on kt
        }//end if !plane_readouts
        */
        mod++;
      } // end loop on jt
    } // end if !modules_names
    //  }//end loop on it
    //}//end if !chambers_names

  } else { // not gem
    err = THaAnalysisObject::LoadDB(file, date, request, prefix.c_str());
    //}
    // Could close the common file already
    fclose(file);
    if (nlogchan == 0)
      nlogchan = nchan;

    if (err)
      return THaAnalysisObject::kInitError;

    // fNChanDet[detname] = nchan;
    // fChanMapStartDet[detname] = chanmapstart;
    (fInvDetMap[detname]).resize(nlogchan + 1); // for ref
    // if(detmap[4]==-1)nparam_mod = 5;
    for (size_t k = 0; k < detmap.size(); k += nparam_mod) {
      crate  = detmap[k];
      slot   = detmap[k + 1];
      ch_lo  = detmap[k + 2];
      ch_hi  = detmap[k + 3];
      ch_ref = detmap[k + 4];
      if (ch_ref == -1) {
        (fInvDetMap[detname])[nlogchan] = detchaninfo(crate, slot, ch_lo);
        continue;
      }
      /*
        if(detname.find("hodo")!=std::string::npos)
        cout << " crate " << crate << " slot " << slot
        << " ch_lo " << ch_lo << " ch_hi " << ch_hi << endl;
      */
      if (chanmap.empty()) {
        for (int i = ch_lo; i <= ch_hi; i++, ch_count++) {
          /*
          if(isgem && i%128==0){
            apv_num++;
            cout << crate << " " << slot << " " << i << " " << apv_num << endl;
          }
          */
          if (ch_count > (int)nlogchan) {
            std::cout << " <1> number of channels defined in detmap ( >= " << ch_count
                      << ") exceeds logical number of channels = " << nlogchan << std::endl;
            return THaAnalysisObject::kInitError;
          }

          (fInvDetMap[detname])[ch_count] = detchaninfo(crate, slot, i);
          // cout << "ch_count " << ch_count << " crate " << crate << " slot " << slot << " i " << i << "
          // &(fInvDetMap[detname]).at(ch_count) " << &(fInvDetMap[detname]).at(ch_count) << endl ;
          /*
            if(detname.find("hodo")!=std::string::npos){
            cout << " crate " << crate << " slot " << slot
            << " i " << i << " ch_count " << ch_count << endl;
            cout << &(fInvDetMap.at(detname)).at(ch_count) << endl;
            }
          */
        }
      } else {
        int chan_offset = 1;
        if (detname.find("hod") != std::string::npos)
          chan_offset = 0;
        for (int i = ch_lo; i <= ch_hi; i++, ch_map++) {
          if (ch_count > (int)nlogchan) {
            std::cout << " <2> number of channels defined in detmap ( >= " << ch_count
                      << ") exceeds logical number of channels = " << nlogchan << std::endl;
            return THaAnalysisObject::kInitError;
          }
          if (fDebug >= 2)
            std::cout << " i = " << i << ", crate = " << crate << ", slot = " << slot << ", ch_count = " << ch_count
                      << " chan = " << chanmap[ch_map] - chan_offset << " (+" << nchan << ") " << std::endl;
          if (chanmap[ch_map] >= 0) {
            if (ch_count < (int)nchan) {
              (fInvDetMap[detname])[chanmap[ch_map] - chan_offset] = detchaninfo(crate, slot, i);
              if (fDebug >= 3)
                std::cout << chanmap[ch_map] - chan_offset << " "
                          << &(fInvDetMap.at(detname)).at(chanmap[ch_map] - chan_offset) << std::endl;
            } else {
              (fInvDetMap[detname])[chanmap[ch_map] + nchan - chan_offset] = detchaninfo(crate, slot, i);
              if (fDebug >= 3)
                std::cout << &(fInvDetMap.at(detname)).at(chanmap[ch_map] + nchan - chan_offset) << std::endl;
            }
            ch_count++;
          }
        }
      }
    }
  } // end else (if isgem)
  /*
  fChansPerSlotDetMap[detname] = cps;
  fSlotsPerCrateDetMap[detname] = spc;
  fFirstSlotDetMap[detname] = fs;
  fFirstCrateDetMap[detname] = fc;
  */

  return (THaAnalysisObject::kOK);
}

//-----------------------------------------------------------------------------
// static inline
void LADSimDecoder::ChanToROC(const std::string &detname, Int_t h_chan, Int_t &crate, Int_t &slot,
                              UShort_t &chan) const {
  // Convert location parameters (row, col, chan) of the given Channel
  // to hardware channel (crate,slot,chan)
  // The (crate,slot,chan) assignment must match the detmap definition in
  // the database!  See TreeSearch/dbconvert.cxx
  // In the case of GRINCH/RICH:
  // crate = GTP; slot = VETROC; chan = PMT. (NINOs are "transparent", in a similar way to the MPDs)

  /*
  int CPS = fChansPerSlotDetMap.at(detname);
  int SPC = fSlotsPerCrateDetMap.at(detname);
  int FS = fFirstSlotDetMap.at(detname);
  int FC = fFirstCrateDetMap.at(detname);

  //div_t d = div( h_chan, fManager->GetChanPerSlot() );
  div_t d = div( h_chan, CPS );
  slot = d.quot;
  chan = d.rem;

  d = div( slot, SPC );
  crate = d.quot+FC;
  slot  = d.rem+FS;
  */
  if ((size_t)h_chan >= fInvDetMap.at(detname).size())
    std::cout << " " << detname << " " << h_chan << " " << &fInvDetMap.at(detname) << std::endl;
  assert((size_t)h_chan < fInvDetMap.at(detname).size());

  if (fDebug > 3) {

    std::cout << &(fInvDetMap.at(detname)).at(h_chan) << std::endl;
  }
  crate = ((fInvDetMap.at(detname)).at(h_chan)).crate;
  slot  = ((fInvDetMap.at(detname)).at(h_chan)).slot;
  chan  = ((fInvDetMap.at(detname)).at(h_chan)).chan;

}

int LADSimDecoder::APVnum(const std::string &detname, Int_t mod, Int_t h_chan, Int_t &crate, Int_t &slot,
                          UShort_t &chan) const {
  chan  = h_chan % 128;
  int n = (h_chan - chan) / 128;

  // std::cout << "(detname, mod, h_chan, chan, n )= (" << detname << ", " << mod << ", "
  // 	    << h_chan << ", " << chan << ", " << n << ")" << std::endl;

  assert((size_t)mod < fInvGEMDetMap.at(detname).size());
  assert((size_t)n < (fInvGEMDetMap.at(detname)[mod]).size());

  // if( mod>fInvGEMDetMap.at(detname).size() ){
  //   std::err << "ERROR: map size =  " << " for detector " << detname
  // 	     << " is smaller than module size =  " << mod << std:: endl;
  //   //return -1;
  // }else if( n>(fInvGEMDetMap.at(detname)[mod]).size() ){
  //   return -1;
  // }
  // if((fInvGEMDetMap.at(detname))[mod][n].chan_lo<=h_chan &&
  // hchan <= (fInvGEMDetMap.at(detname))[mod][n].chan_hi){
  crate = ((fInvGEMDetMap.at(detname))[mod][n]).crate;
  slot  = ((fInvGEMDetMap.at(detname))[mod][n]).slot;
  return ((fInvGEMDetMap.at(detname))[mod][n]).apvnum;
  //}else{
  // return -1;
  //}
}

/*
//-----------------------------------------------------------------------------
static inline
Int_t MakeROCKey( Int_t crate, Int_t slot, Int_t chan )
{
  return chan;// +
  //fManager->GetChanPerSlot()*( slot + fManager->GetSlotPerCrate()*crate );
}

//-----------------------------------------------------------------------------
Int_t LADSimDecoder::ChanFromROC( Int_t crate, Int_t slot, Int_t chan ) const
{
  // Return index of digitized strip correspomding to hardware channel
  // (crate,slot,chan)

  if( fPMTMap.empty() )
    return -1;

  PMTMap_t::const_iterator found = fPMTMap.find( MakeROCKey(crate,slot,chan) );
  if( found == fPMTMap.end() )
    return -1;

  return found->second;
}
*/
