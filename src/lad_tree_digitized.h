//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Jun  7 11:45:28 2021 by ROOT version 6.14/04
// from TTree T/Geant4 SBS Simulation
// found on file: /volatile/halla/sbs/efuchey/digitized_files/simdigtest_2.root
//////////////////////////////////////////////////////////

#ifndef lad_tree_digitized_h
#define lad_tree_digitized_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"

class lad_tree_digitized {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Double_t        ev_count;
   Double_t        ev_rate;
   Double_t        ev_solang;
   Double_t        ev_sigma;
   Double_t        ev_W2;
   Double_t        ev_xbj;
   Double_t        ev_Q2;
   Double_t        ev_th;
   Double_t        ev_ph;
   Double_t        ev_Aperp;
   Double_t        ev_Apar;
   Double_t        ev_Pt;
   Double_t        ev_Pl;
   Double_t        ev_vx;
   Double_t        ev_vy;
   Double_t        ev_vz;
   Double_t        ev_ep;
   Double_t        ev_np;
   Double_t        ev_epx;
   Double_t        ev_epy;
   Double_t        ev_epz;
   Double_t        ev_npx;
   Double_t        ev_npy;
   Double_t        ev_npz;
   Double_t        ev_nth;
   Double_t        ev_nph;
   Double_t        ev_pmperp;
   Double_t        ev_pmpar;
   Double_t        ev_pmparsm;
   Double_t        ev_z;
   Double_t        ev_phperp;
   Double_t        ev_phih;
   Double_t        ev_phiS;
   Double_t        ev_thetaS;
   Double_t        ev_MX2;
   Double_t        ev_Sx;
   Double_t        ev_Sy;
   Double_t        ev_Sz;
   Double_t        ev_s;
   Double_t        ev_t;
   Double_t        ev_u;
   Double_t        ev_costhetaCM;
   Double_t        ev_Egamma;
   Int_t           ev_nucl;
   Int_t           ev_fnucl;
   Int_t           ev_hadr;
   Int_t           ev_earmaccept;
   Int_t           ev_harmaccept;
   Int_t           Earm_BBGEM_hit_nhits;
   std::vector<int>     *Earm_BBGEM_hit_plane;
   std::vector<int>     *Earm_BBGEM_hit_strip;
   std::vector<double>  *Earm_BBGEM_hit_x;
   std::vector<double>  *Earm_BBGEM_hit_y;
   std::vector<double>  *Earm_BBGEM_hit_z;
   std::vector<double>  *Earm_BBGEM_hit_polx;
   std::vector<double>  *Earm_BBGEM_hit_poly;
   std::vector<double>  *Earm_BBGEM_hit_polz;
   std::vector<double>  *Earm_BBGEM_hit_t;
   std::vector<double>  *Earm_BBGEM_hit_trms;
   std::vector<double>  *Earm_BBGEM_hit_tmin;
   std::vector<double>  *Earm_BBGEM_hit_tmax;
   std::vector<double>  *Earm_BBGEM_hit_tx;
   std::vector<double>  *Earm_BBGEM_hit_ty;
   std::vector<double>  *Earm_BBGEM_hit_xin;
   std::vector<double>  *Earm_BBGEM_hit_yin;
   std::vector<double>  *Earm_BBGEM_hit_zin;
   std::vector<double>  *Earm_BBGEM_hit_xout;
   std::vector<double>  *Earm_BBGEM_hit_yout;
   std::vector<double>  *Earm_BBGEM_hit_zout;
   std::vector<double>  *Earm_BBGEM_hit_txp;
   std::vector<double>  *Earm_BBGEM_hit_typ;
   std::vector<double>  *Earm_BBGEM_hit_xg;
   std::vector<double>  *Earm_BBGEM_hit_yg;
   std::vector<double>  *Earm_BBGEM_hit_zg;
   std::vector<int>     *Earm_BBGEM_hit_trid;
   std::vector<int>     *Earm_BBGEM_hit_mid;
   std::vector<int>     *Earm_BBGEM_hit_pid;
   std::vector<double>  *Earm_BBGEM_hit_vx;
   std::vector<double>  *Earm_BBGEM_hit_vy;
   std::vector<double>  *Earm_BBGEM_hit_vz;
   std::vector<double>  *Earm_BBGEM_hit_p;
   std::vector<double>  *Earm_BBGEM_hit_edep;
   std::vector<double>  *Earm_BBGEM_hit_beta;
   std::vector<int>     *Earm_BBGEM_hit_otridx;
   std::vector<int>     *Earm_BBGEM_hit_ptridx;
   std::vector<int>     *Earm_BBGEM_hit_sdtridx;
   Int_t           Earm_BBGEM_Track_ntracks;
   std::vector<int>     *Earm_BBGEM_Track_TID;
   std::vector<int>     *Earm_BBGEM_Track_PID;
   std::vector<int>     *Earm_BBGEM_Track_MID;
   std::vector<int>     *Earm_BBGEM_Track_NumHits;
   std::vector<int>     *Earm_BBGEM_Track_NumPlanes;
   std::vector<int>     *Earm_BBGEM_Track_NDF;
   std::vector<double>  *Earm_BBGEM_Track_Chi2fit;
   std::vector<double>  *Earm_BBGEM_Track_Chi2true;
   std::vector<double>  *Earm_BBGEM_Track_X;
   std::vector<double>  *Earm_BBGEM_Track_Y;
   std::vector<double>  *Earm_BBGEM_Track_Xp;
   std::vector<double>  *Earm_BBGEM_Track_Yp;
   std::vector<double>  *Earm_BBGEM_Track_T;
   std::vector<double>  *Earm_BBGEM_Track_P;
   std::vector<double>  *Earm_BBGEM_Track_Sx;
   std::vector<double>  *Earm_BBGEM_Track_Sy;
   std::vector<double>  *Earm_BBGEM_Track_Sz;
   std::vector<double>  *Earm_BBGEM_Track_Xfit;
   std::vector<double>  *Earm_BBGEM_Track_Yfit;
   std::vector<double>  *Earm_BBGEM_Track_Xpfit;
   std::vector<double>  *Earm_BBGEM_Track_Ypfit;
   std::vector<int>     *Earm_BBGEM_Track_otridx;
   std::vector<int>     *Earm_BBGEM_Track_ptridx;
   std::vector<int>     *Earm_BBGEM_Track_sdtridx;
   Double_t        Earm_BBHodoScint_det_esum;
   Int_t           Earm_BBHodoScint_hit_nhits;
   std::vector<int>     *Earm_BBHodoScint_hit_row;
   std::vector<int>     *Earm_BBHodoScint_hit_col;
   std::vector<int>     *Earm_BBHodoScint_hit_cell;
   std::vector<int>     *Earm_BBHodoScint_hit_plane;
   std::vector<int>     *Earm_BBHodoScint_hit_wire;
   std::vector<double>  *Earm_BBHodoScint_hit_xcell;
   std::vector<double>  *Earm_BBHodoScint_hit_ycell;
   std::vector<double>  *Earm_BBHodoScint_hit_zcell;
   std::vector<double>  *Earm_BBHodoScint_hit_xcellg;
   std::vector<double>  *Earm_BBHodoScint_hit_ycellg;
   std::vector<double>  *Earm_BBHodoScint_hit_zcellg;
   std::vector<double>  *Earm_BBHodoScint_hit_xhit;
   std::vector<double>  *Earm_BBHodoScint_hit_yhit;
   std::vector<double>  *Earm_BBHodoScint_hit_zhit;
   std::vector<double>  *Earm_BBHodoScint_hit_xhitg;
   std::vector<double>  *Earm_BBHodoScint_hit_yhitg;
   std::vector<double>  *Earm_BBHodoScint_hit_zhitg;
   std::vector<double>  *Earm_BBHodoScint_hit_sumedep;
   std::vector<double>  *Earm_BBHodoScint_hit_tavg;
   std::vector<double>  *Earm_BBHodoScint_hit_trms;
   std::vector<double>  *Earm_BBHodoScint_hit_tmin;
   std::vector<double>  *Earm_BBHodoScint_hit_tmax;
   std::vector<int>     *Earm_BBHodoScint_hit_otridx;
   std::vector<int>     *Earm_BBHodoScint_hit_ptridx;
   std::vector<int>     *Earm_BBHodoScint_hit_sdtridx;
   Int_t           OTrack_ntracks;
   std::vector<int>     *OTrack_TID;
   std::vector<int>     *OTrack_MID;
   std::vector<int>     *OTrack_PID;
   std::vector<double>  *OTrack_posx;
   std::vector<double>  *OTrack_posy;
   std::vector<double>  *OTrack_posz;
   std::vector<double>  *OTrack_momx;
   std::vector<double>  *OTrack_momy;
   std::vector<double>  *OTrack_momz;
   std::vector<double>  *OTrack_polx;
   std::vector<double>  *OTrack_poly;
   std::vector<double>  *OTrack_polz;
   std::vector<double>  *OTrack_Etot;
   std::vector<double>  *OTrack_T;
   Int_t           PTrack_ntracks;
   std::vector<int>     *PTrack_TID;
   std::vector<int>     *PTrack_PID;
   std::vector<double>  *PTrack_posx;
   std::vector<double>  *PTrack_posy;
   std::vector<double>  *PTrack_posz;
   std::vector<double>  *PTrack_momx;
   std::vector<double>  *PTrack_momy;
   std::vector<double>  *PTrack_momz;
   std::vector<double>  *PTrack_polx;
   std::vector<double>  *PTrack_poly;
   std::vector<double>  *PTrack_polz;
   std::vector<double>  *PTrack_Etot;
   std::vector<double>  *PTrack_T;
   Int_t           SDTrack_ntracks;
   std::vector<int>     *SDTrack_TID;
   std::vector<int>     *SDTrack_MID;
   std::vector<int>     *SDTrack_PID;
   std::vector<double>  *SDTrack_posx;
   std::vector<double>  *SDTrack_posy;
   std::vector<double>  *SDTrack_posz;
   std::vector<double>  *SDTrack_momx;
   std::vector<double>  *SDTrack_momy;
   std::vector<double>  *SDTrack_momz;
   std::vector<double>  *SDTrack_polx;
   std::vector<double>  *SDTrack_poly;
   std::vector<double>  *SDTrack_polz;
   std::vector<double>  *SDTrack_Etot;
   std::vector<double>  *SDTrack_T;
   std::vector<double>  *SDTrack_vx;
   std::vector<double>  *SDTrack_vy;
   std::vector<double>  *SDTrack_vz;
   std::vector<double>  *SDTrack_vnx;
   std::vector<double>  *SDTrack_vny;
   std::vector<double>  *SDTrack_vnz;
   std::vector<double>  *SDTrack_vEkin;
   Int_t           LAD_Hodo_dighit_nchan;
   std::vector<int>     *LAD_Hodo_dighit_chan;
   std::vector<int>     *LAD_Hodo_dighit_adc;
   std::vector<int>     *LAD_Hodo_dighit_tdc_l;
   std::vector<int>     *LAD_Hodo_dighit_tdc_t;
   Int_t           Earm_BBGEM_dighit_nstrips;
   std::vector<int>     *Earm_BBGEM_dighit_module;
   std::vector<int>     *Earm_BBGEM_dighit_strip;
   std::vector<int>     *Earm_BBGEM_dighit_adc;
   std::vector<int>     *Earm_BBGEM_dighit_samp;

   // List of branches
   TBranch        *b_simc_sigma;   //!
   TBranch        *b_simc_Weight;   //!
   TBranch        *b_simc_Q2;   //!
   TBranch        *b_simc_xbj;   //!
   TBranch        *b_simc_nu;   //!
   TBranch        *b_simc_W;   //!
   TBranch        *b_simc_epsilon;   //!
   TBranch        *b_simc_Ebeam;   //!
   TBranch        *b_simc_p_e;   //!
   TBranch        *b_simc_theta_e;   //!
   TBranch        *b_simc_phi_e;   //!
   TBranch        *b_simc_px_e;   //!
   TBranch        *b_simc_py_e;   //!
   TBranch        *b_simc_pz_e;   //!
   TBranch        *b_simc_fnucl;   //!
   TBranch        *b_simc_p_n;   //!
   TBranch        *b_simc_theta_n;   //!
   TBranch        *b_simc_phi_n;   //!
   TBranch        *b_simc_px_n;   //!
   TBranch        *b_simc_py_n;   //!
   TBranch        *b_simc_pz_n;   //!
   TBranch        *b_simc_vx;   //!
   TBranch        *b_simc_vy;   //!
   TBranch        *b_simc_vz;   //!
   TBranch        *b_simc_veE;   //!
   TBranch        *b_simc_vetheta;   //!
   TBranch        *b_ev;   //!
   TBranch        *b_Earm_BBGEM_hit_nhits;   //!
   TBranch        *b_Earm_BBGEM_hit_plane;   //!
   TBranch        *b_Earm_BBGEM_hit_strip;   //!
   TBranch        *b_Earm_BBGEM_hit_x;   //!
   TBranch        *b_Earm_BBGEM_hit_y;   //!
   TBranch        *b_Earm_BBGEM_hit_z;   //!
   TBranch        *b_Earm_BBGEM_hit_polx;   //!
   TBranch        *b_Earm_BBGEM_hit_poly;   //!
   TBranch        *b_Earm_BBGEM_hit_polz;   //!
   TBranch        *b_Earm_BBGEM_hit_t;   //!
   TBranch        *b_Earm_BBGEM_hit_trms;   //!
   TBranch        *b_Earm_BBGEM_hit_tmin;   //!
   TBranch        *b_Earm_BBGEM_hit_tmax;   //!
   TBranch        *b_Earm_BBGEM_hit_tx;   //!
   TBranch        *b_Earm_BBGEM_hit_ty;   //!
   TBranch        *b_Earm_BBGEM_hit_xin;   //!
   TBranch        *b_Earm_BBGEM_hit_yin;   //!
   TBranch        *b_Earm_BBGEM_hit_zin;   //!
   TBranch        *b_Earm_BBGEM_hit_xout;   //!
   TBranch        *b_Earm_BBGEM_hit_yout;   //!
   TBranch        *b_Earm_BBGEM_hit_zout;   //!
   TBranch        *b_Earm_BBGEM_hit_txp;   //!
   TBranch        *b_Earm_BBGEM_hit_typ;   //!
   TBranch        *b_Earm_BBGEM_hit_xg;   //!
   TBranch        *b_Earm_BBGEM_hit_yg;   //!
   TBranch        *b_Earm_BBGEM_hit_zg;   //!
   TBranch        *b_Earm_BBGEM_hit_trid;   //!
   TBranch        *b_Earm_BBGEM_hit_mid;   //!
   TBranch        *b_Earm_BBGEM_hit_pid;   //!
   TBranch        *b_Earm_BBGEM_hit_vx;   //!
   TBranch        *b_Earm_BBGEM_hit_vy;   //!
   TBranch        *b_Earm_BBGEM_hit_vz;   //!
   TBranch        *b_Earm_BBGEM_hit_p;   //!
   TBranch        *b_Earm_BBGEM_hit_edep;   //!
   TBranch        *b_Earm_BBGEM_hit_beta;   //!
   TBranch        *b_Earm_BBGEM_hit_otridx;   //!
   TBranch        *b_Earm_BBGEM_hit_ptridx;   //!
   TBranch        *b_Earm_BBGEM_hit_sdtridx;   //!
   TBranch        *b_Earm_BBGEM_Track_ntracks;   //!
   TBranch        *b_Earm_BBGEM_Track_TID;   //!
   TBranch        *b_Earm_BBGEM_Track_PID;   //!
   TBranch        *b_Earm_BBGEM_Track_MID;   //!
   TBranch        *b_Earm_BBGEM_Track_NumHits;   //!
   TBranch        *b_Earm_BBGEM_Track_NumPlanes;   //!
   TBranch        *b_Earm_BBGEM_Track_NDF;   //!
   TBranch        *b_Earm_BBGEM_Track_Chi2fit;   //!
   TBranch        *b_Earm_BBGEM_Track_Chi2true;   //!
   TBranch        *b_Earm_BBGEM_Track_X;   //!
   TBranch        *b_Earm_BBGEM_Track_Y;   //!
   TBranch        *b_Earm_BBGEM_Track_Xp;   //!
   TBranch        *b_Earm_BBGEM_Track_Yp;   //!
   TBranch        *b_Earm_BBGEM_Track_T;   //!
   TBranch        *b_Earm_BBGEM_Track_P;   //!
   TBranch        *b_Earm_BBGEM_Track_Sx;   //!
   TBranch        *b_Earm_BBGEM_Track_Sy;   //!
   TBranch        *b_Earm_BBGEM_Track_Sz;   //!
   TBranch        *b_Earm_BBGEM_Track_Xfit;   //!
   TBranch        *b_Earm_BBGEM_Track_Yfit;   //!
   TBranch        *b_Earm_BBGEM_Track_Xpfit;   //!
   TBranch        *b_Earm_BBGEM_Track_Ypfit;   //!
   TBranch        *b_Earm_BBGEM_Track_otridx;   //!
   TBranch        *b_Earm_BBGEM_Track_ptridx;   //!
   TBranch        *b_Earm_BBGEM_Track_sdtridx;   //!
   TBranch        *b_Earm_BBHodoScint_det_esum;   //!
   TBranch        *b_Earm_BBHodoScint_hit_nhits;   //!
   TBranch        *b_Earm_BBHodoScint_hit_row;   //!
   TBranch        *b_Earm_BBHodoScint_hit_col;   //!
   TBranch        *b_Earm_BBHodoScint_hit_cell;   //!
   TBranch        *b_Earm_BBHodoScint_hit_plane;   //!
   TBranch        *b_Earm_BBHodoScint_hit_wire;   //!
   TBranch        *b_Earm_BBHodoScint_hit_xcell;   //!
   TBranch        *b_Earm_BBHodoScint_hit_ycell;   //!
   TBranch        *b_Earm_BBHodoScint_hit_zcell;   //!
   TBranch        *b_Earm_BBHodoScint_hit_xcellg;   //!
   TBranch        *b_Earm_BBHodoScint_hit_ycellg;   //!
   TBranch        *b_Earm_BBHodoScint_hit_zcellg;   //!
   TBranch        *b_Earm_BBHodoScint_hit_xhit;   //!
   TBranch        *b_Earm_BBHodoScint_hit_yhit;   //!
   TBranch        *b_Earm_BBHodoScint_hit_zhit;   //!
   TBranch        *b_Earm_BBHodoScint_hit_xhitg;   //!
   TBranch        *b_Earm_BBHodoScint_hit_yhitg;   //!
   TBranch        *b_Earm_BBHodoScint_hit_zhitg;   //!
   TBranch        *b_Earm_BBHodoScint_hit_sumedep;   //!
   TBranch        *b_Earm_BBHodoScint_hit_tavg;   //!
   TBranch        *b_Earm_BBHodoScint_hit_trms;   //!
   TBranch        *b_Earm_BBHodoScint_hit_tmin;   //!
   TBranch        *b_Earm_BBHodoScint_hit_tmax;   //!
   TBranch        *b_Earm_BBHodoScint_hit_otridx;   //!
   TBranch        *b_Earm_BBHodoScint_hit_ptridx;   //!
   TBranch        *b_Earm_BBHodoScint_hit_sdtridx;   //!
   TBranch        *b_OTrack_ntracks;   //!
   TBranch        *b_OTrack_TID;   //!
   TBranch        *b_OTrack_MID;   //!
   TBranch        *b_OTrack_PID;   //!
   TBranch        *b_OTrack_posx;   //!
   TBranch        *b_OTrack_posy;   //!
   TBranch        *b_OTrack_posz;   //!
   TBranch        *b_OTrack_momx;   //!
   TBranch        *b_OTrack_momy;   //!
   TBranch        *b_OTrack_momz;   //!
   TBranch        *b_OTrack_polx;   //!
   TBranch        *b_OTrack_poly;   //!
   TBranch        *b_OTrack_polz;   //!
   TBranch        *b_OTrack_Etot;   //!
   TBranch        *b_OTrack_T;   //!
   TBranch        *b_PTrack_ntracks;   //!
   TBranch        *b_PTrack_TID;   //!
   TBranch        *b_PTrack_PID;   //!
   TBranch        *b_PTrack_posx;   //!
   TBranch        *b_PTrack_posy;   //!
   TBranch        *b_PTrack_posz;   //!
   TBranch        *b_PTrack_momx;   //!
   TBranch        *b_PTrack_momy;   //!
   TBranch        *b_PTrack_momz;   //!
   TBranch        *b_PTrack_polx;   //!
   TBranch        *b_PTrack_poly;   //!
   TBranch        *b_PTrack_polz;   //!
   TBranch        *b_PTrack_Etot;   //!
   TBranch        *b_PTrack_T;   //!
   TBranch        *b_SDTrack_ntracks;   //!
   TBranch        *b_SDTrack_TID;   //!
   TBranch        *b_SDTrack_MID;   //!
   TBranch        *b_SDTrack_PID;   //!
   TBranch        *b_SDTrack_posx;   //!
   TBranch        *b_SDTrack_posy;   //!
   TBranch        *b_SDTrack_posz;   //!
   TBranch        *b_SDTrack_momx;   //!
   TBranch        *b_SDTrack_momy;   //!
   TBranch        *b_SDTrack_momz;   //!
   TBranch        *b_SDTrack_polx;   //!
   TBranch        *b_SDTrack_poly;   //!
   TBranch        *b_SDTrack_polz;   //!
   TBranch        *b_SDTrack_Etot;   //!
   TBranch        *b_SDTrack_T;   //!
   TBranch        *b_SDTrack_vx;   //!
   TBranch        *b_SDTrack_vy;   //!
   TBranch        *b_SDTrack_vz;   //!
   TBranch        *b_SDTrack_vnx;   //!
   TBranch        *b_SDTrack_vny;   //!
   TBranch        *b_SDTrack_vnz;   //!
   TBranch        *b_SDTrack_vEkin;   //!

   TBranch        *b_LAD_Hodo_dighit_nchan;   //!
   TBranch        *b_LAD_Hodo_dighit_chan;   //!
   TBranch        *b_LAD_Hodo_dighit_adc;   //!
   TBranch        *b_LAD_Hodo_dighit_tdc_l;   //!
   TBranch        *b_LAD_Hodo_dighit_tdc_t;   //!
   TBranch        *b_Earm_BBGEM_dighit_nstrips;   //!
   TBranch        *b_Earm_BBGEM_dighit_module;   //!
   TBranch        *b_Earm_BBGEM_dighit_strip;   //!
   TBranch        *b_Earm_BBGEM_dighit_adc;   //!
   TBranch        *b_Earm_BBGEM_dighit_samp;   //!

   lad_tree_digitized(TTree *tree=0);
   virtual ~lad_tree_digitized();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};
#endif

#ifdef lad_tree_digitized_cxx
lad_tree_digitized::lad_tree_digitized(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/volatile/halla/sbs/efuchey/digitized_files/simdigtest_2.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/volatile/halla/sbs/efuchey/digitized_files/simdigtest_2.root");
      }
      f->GetObject("T",tree);

   }
   Init(tree);
}

lad_tree_digitized::~lad_tree_digitized()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t lad_tree_digitized::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t lad_tree_digitized::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void lad_tree_digitized::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   Earm_BBGEM_hit_plane = 0;
   Earm_BBGEM_hit_strip = 0;
   Earm_BBGEM_hit_x = 0;
   Earm_BBGEM_hit_y = 0;
   Earm_BBGEM_hit_z = 0;
   Earm_BBGEM_hit_polx = 0;
   Earm_BBGEM_hit_poly = 0;
   Earm_BBGEM_hit_polz = 0;
   Earm_BBGEM_hit_t = 0;
   Earm_BBGEM_hit_trms = 0;
   Earm_BBGEM_hit_tmin = 0;
   Earm_BBGEM_hit_tmax = 0;
   Earm_BBGEM_hit_tx = 0;
   Earm_BBGEM_hit_ty = 0;
   Earm_BBGEM_hit_xin = 0;
   Earm_BBGEM_hit_yin = 0;
   Earm_BBGEM_hit_zin = 0;
   Earm_BBGEM_hit_xout = 0;
   Earm_BBGEM_hit_yout = 0;
   Earm_BBGEM_hit_zout = 0;
   Earm_BBGEM_hit_txp = 0;
   Earm_BBGEM_hit_typ = 0;
   Earm_BBGEM_hit_xg = 0;
   Earm_BBGEM_hit_yg = 0;
   Earm_BBGEM_hit_zg = 0;
   Earm_BBGEM_hit_trid = 0;
   Earm_BBGEM_hit_mid = 0;
   Earm_BBGEM_hit_pid = 0;
   Earm_BBGEM_hit_vx = 0;
   Earm_BBGEM_hit_vy = 0;
   Earm_BBGEM_hit_vz = 0;
   Earm_BBGEM_hit_p = 0;
   Earm_BBGEM_hit_edep = 0;
   Earm_BBGEM_hit_beta = 0;
   Earm_BBGEM_hit_otridx = 0;
   Earm_BBGEM_hit_ptridx = 0;
   Earm_BBGEM_hit_sdtridx = 0;
   Earm_BBGEM_Track_TID = 0;
   Earm_BBGEM_Track_PID = 0;
   Earm_BBGEM_Track_MID = 0;
   Earm_BBGEM_Track_NumHits = 0;
   Earm_BBGEM_Track_NumPlanes = 0;
   Earm_BBGEM_Track_NDF = 0;
   Earm_BBGEM_Track_Chi2fit = 0;
   Earm_BBGEM_Track_Chi2true = 0;
   Earm_BBGEM_Track_X = 0;
   Earm_BBGEM_Track_Y = 0;
   Earm_BBGEM_Track_Xp = 0;
   Earm_BBGEM_Track_Yp = 0;
   Earm_BBGEM_Track_T = 0;
   Earm_BBGEM_Track_P = 0;
   Earm_BBGEM_Track_Sx = 0;
   Earm_BBGEM_Track_Sy = 0;
   Earm_BBGEM_Track_Sz = 0;
   Earm_BBGEM_Track_Xfit = 0;
   Earm_BBGEM_Track_Yfit = 0;
   Earm_BBGEM_Track_Xpfit = 0;
   Earm_BBGEM_Track_Ypfit = 0;
   Earm_BBGEM_Track_otridx = 0;
   Earm_BBGEM_Track_ptridx = 0;
   Earm_BBGEM_Track_sdtridx = 0;
   Earm_BBHodoScint_hit_row = 0;
   Earm_BBHodoScint_hit_col = 0;
   Earm_BBHodoScint_hit_cell = 0;
   Earm_BBHodoScint_hit_plane = 0;
   Earm_BBHodoScint_hit_wire = 0;
   Earm_BBHodoScint_hit_xcell = 0;
   Earm_BBHodoScint_hit_ycell = 0;
   Earm_BBHodoScint_hit_zcell = 0;
   Earm_BBHodoScint_hit_xcellg = 0;
   Earm_BBHodoScint_hit_ycellg = 0;
   Earm_BBHodoScint_hit_zcellg = 0;
   Earm_BBHodoScint_hit_xhit = 0;
   Earm_BBHodoScint_hit_yhit = 0;
   Earm_BBHodoScint_hit_zhit = 0;
   Earm_BBHodoScint_hit_xhitg = 0;
   Earm_BBHodoScint_hit_yhitg = 0;
   Earm_BBHodoScint_hit_zhitg = 0;
   Earm_BBHodoScint_hit_sumedep = 0;
   Earm_BBHodoScint_hit_tavg = 0;
   Earm_BBHodoScint_hit_trms = 0;
   Earm_BBHodoScint_hit_tmin = 0;
   Earm_BBHodoScint_hit_tmax = 0;
   Earm_BBHodoScint_hit_otridx = 0;
   Earm_BBHodoScint_hit_ptridx = 0;
   Earm_BBHodoScint_hit_sdtridx = 0;
   OTrack_TID = 0;
   OTrack_MID = 0;
   OTrack_PID = 0;
   OTrack_posx = 0;
   OTrack_posy = 0;
   OTrack_posz = 0;
   OTrack_momx = 0;
   OTrack_momy = 0;
   OTrack_momz = 0;
   OTrack_polx = 0;
   OTrack_poly = 0;
   OTrack_polz = 0;
   OTrack_Etot = 0;
   OTrack_T = 0;
   PTrack_TID = 0;
   PTrack_PID = 0;
   PTrack_posx = 0;
   PTrack_posy = 0;
   PTrack_posz = 0;
   PTrack_momx = 0;
   PTrack_momy = 0;
   PTrack_momz = 0;
   PTrack_polx = 0;
   PTrack_poly = 0;
   PTrack_polz = 0;
   PTrack_Etot = 0;
   PTrack_T = 0;
   SDTrack_TID = 0;
   SDTrack_MID = 0;
   SDTrack_PID = 0;
   SDTrack_posx = 0;
   SDTrack_posy = 0;
   SDTrack_posz = 0;
   SDTrack_momx = 0;
   SDTrack_momy = 0;
   SDTrack_momz = 0;
   SDTrack_polx = 0;
   SDTrack_poly = 0;
   SDTrack_polz = 0;
   SDTrack_Etot = 0;
   SDTrack_T = 0;
   SDTrack_vx = 0;
   SDTrack_vy = 0;
   SDTrack_vz = 0;
   SDTrack_vnx = 0;
   SDTrack_vny = 0;
   SDTrack_vnz = 0;
   SDTrack_vEkin = 0;
   LAD_Hodo_dighit_chan = 0;
   LAD_Hodo_dighit_adc = 0;
   LAD_Hodo_dighit_tdc_l = 0;
   LAD_Hodo_dighit_tdc_t = 0;
   Earm_BBGEM_dighit_module = 0;
   Earm_BBGEM_dighit_strip = 0;
   Earm_BBGEM_dighit_adc = 0;
   Earm_BBGEM_dighit_samp = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("ev", &ev_count, &b_ev);
   fChain->SetBranchAddress("Earm.BBGEM.hit.nhits", &Earm_BBGEM_hit_nhits, &b_Earm_BBGEM_hit_nhits);
   fChain->SetBranchAddress("Earm.BBGEM.hit.plane", &Earm_BBGEM_hit_plane, &b_Earm_BBGEM_hit_plane);
   fChain->SetBranchAddress("Earm.BBGEM.hit.strip", &Earm_BBGEM_hit_strip, &b_Earm_BBGEM_hit_strip);
   fChain->SetBranchAddress("Earm.BBGEM.hit.x", &Earm_BBGEM_hit_x, &b_Earm_BBGEM_hit_x);
   fChain->SetBranchAddress("Earm.BBGEM.hit.y", &Earm_BBGEM_hit_y, &b_Earm_BBGEM_hit_y);
   fChain->SetBranchAddress("Earm.BBGEM.hit.z", &Earm_BBGEM_hit_z, &b_Earm_BBGEM_hit_z);
   fChain->SetBranchAddress("Earm.BBGEM.hit.polx", &Earm_BBGEM_hit_polx, &b_Earm_BBGEM_hit_polx);
   fChain->SetBranchAddress("Earm.BBGEM.hit.poly", &Earm_BBGEM_hit_poly, &b_Earm_BBGEM_hit_poly);
   fChain->SetBranchAddress("Earm.BBGEM.hit.polz", &Earm_BBGEM_hit_polz, &b_Earm_BBGEM_hit_polz);
   fChain->SetBranchAddress("Earm.BBGEM.hit.t", &Earm_BBGEM_hit_t, &b_Earm_BBGEM_hit_t);
   fChain->SetBranchAddress("Earm.BBGEM.hit.trms", &Earm_BBGEM_hit_trms, &b_Earm_BBGEM_hit_trms);
   fChain->SetBranchAddress("Earm.BBGEM.hit.tmin", &Earm_BBGEM_hit_tmin, &b_Earm_BBGEM_hit_tmin);
   fChain->SetBranchAddress("Earm.BBGEM.hit.tmax", &Earm_BBGEM_hit_tmax, &b_Earm_BBGEM_hit_tmax);
   fChain->SetBranchAddress("Earm.BBGEM.hit.tx", &Earm_BBGEM_hit_tx, &b_Earm_BBGEM_hit_tx);
   fChain->SetBranchAddress("Earm.BBGEM.hit.ty", &Earm_BBGEM_hit_ty, &b_Earm_BBGEM_hit_ty);
   fChain->SetBranchAddress("Earm.BBGEM.hit.xin", &Earm_BBGEM_hit_xin, &b_Earm_BBGEM_hit_xin);
   fChain->SetBranchAddress("Earm.BBGEM.hit.yin", &Earm_BBGEM_hit_yin, &b_Earm_BBGEM_hit_yin);
   fChain->SetBranchAddress("Earm.BBGEM.hit.zin", &Earm_BBGEM_hit_zin, &b_Earm_BBGEM_hit_zin);
   fChain->SetBranchAddress("Earm.BBGEM.hit.xout", &Earm_BBGEM_hit_xout, &b_Earm_BBGEM_hit_xout);
   fChain->SetBranchAddress("Earm.BBGEM.hit.yout", &Earm_BBGEM_hit_yout, &b_Earm_BBGEM_hit_yout);
   fChain->SetBranchAddress("Earm.BBGEM.hit.zout", &Earm_BBGEM_hit_zout, &b_Earm_BBGEM_hit_zout);
   fChain->SetBranchAddress("Earm.BBGEM.hit.txp", &Earm_BBGEM_hit_txp, &b_Earm_BBGEM_hit_txp);
   fChain->SetBranchAddress("Earm.BBGEM.hit.typ", &Earm_BBGEM_hit_typ, &b_Earm_BBGEM_hit_typ);
   fChain->SetBranchAddress("Earm.BBGEM.hit.xg", &Earm_BBGEM_hit_xg, &b_Earm_BBGEM_hit_xg);
   fChain->SetBranchAddress("Earm.BBGEM.hit.yg", &Earm_BBGEM_hit_yg, &b_Earm_BBGEM_hit_yg);
   fChain->SetBranchAddress("Earm.BBGEM.hit.zg", &Earm_BBGEM_hit_zg, &b_Earm_BBGEM_hit_zg);
   fChain->SetBranchAddress("Earm.BBGEM.hit.trid", &Earm_BBGEM_hit_trid, &b_Earm_BBGEM_hit_trid);
   fChain->SetBranchAddress("Earm.BBGEM.hit.mid", &Earm_BBGEM_hit_mid, &b_Earm_BBGEM_hit_mid);
   fChain->SetBranchAddress("Earm.BBGEM.hit.pid", &Earm_BBGEM_hit_pid, &b_Earm_BBGEM_hit_pid);
   fChain->SetBranchAddress("Earm.BBGEM.hit.vx", &Earm_BBGEM_hit_vx, &b_Earm_BBGEM_hit_vx);
   fChain->SetBranchAddress("Earm.BBGEM.hit.vy", &Earm_BBGEM_hit_vy, &b_Earm_BBGEM_hit_vy);
   fChain->SetBranchAddress("Earm.BBGEM.hit.vz", &Earm_BBGEM_hit_vz, &b_Earm_BBGEM_hit_vz);
   fChain->SetBranchAddress("Earm.BBGEM.hit.p", &Earm_BBGEM_hit_p, &b_Earm_BBGEM_hit_p);
   fChain->SetBranchAddress("Earm.BBGEM.hit.edep", &Earm_BBGEM_hit_edep, &b_Earm_BBGEM_hit_edep);
   fChain->SetBranchAddress("Earm.BBGEM.hit.beta", &Earm_BBGEM_hit_beta, &b_Earm_BBGEM_hit_beta);
   fChain->SetBranchAddress("Earm.BBGEM.hit.otridx", &Earm_BBGEM_hit_otridx, &b_Earm_BBGEM_hit_otridx);
   fChain->SetBranchAddress("Earm.BBGEM.hit.ptridx", &Earm_BBGEM_hit_ptridx, &b_Earm_BBGEM_hit_ptridx);
   fChain->SetBranchAddress("Earm.BBGEM.hit.sdtridx", &Earm_BBGEM_hit_sdtridx, &b_Earm_BBGEM_hit_sdtridx);
   fChain->SetBranchAddress("Earm.BBGEM.Track.ntracks", &Earm_BBGEM_Track_ntracks, &b_Earm_BBGEM_Track_ntracks);
   fChain->SetBranchAddress("Earm.BBGEM.Track.TID", &Earm_BBGEM_Track_TID, &b_Earm_BBGEM_Track_TID);
   fChain->SetBranchAddress("Earm.BBGEM.Track.PID", &Earm_BBGEM_Track_PID, &b_Earm_BBGEM_Track_PID);
   fChain->SetBranchAddress("Earm.BBGEM.Track.MID", &Earm_BBGEM_Track_MID, &b_Earm_BBGEM_Track_MID);
   fChain->SetBranchAddress("Earm.BBGEM.Track.NumHits", &Earm_BBGEM_Track_NumHits, &b_Earm_BBGEM_Track_NumHits);
   fChain->SetBranchAddress("Earm.BBGEM.Track.NumPlanes", &Earm_BBGEM_Track_NumPlanes, &b_Earm_BBGEM_Track_NumPlanes);
   fChain->SetBranchAddress("Earm.BBGEM.Track.NDF", &Earm_BBGEM_Track_NDF, &b_Earm_BBGEM_Track_NDF);
   fChain->SetBranchAddress("Earm.BBGEM.Track.Chi2fit", &Earm_BBGEM_Track_Chi2fit, &b_Earm_BBGEM_Track_Chi2fit);
   fChain->SetBranchAddress("Earm.BBGEM.Track.Chi2true", &Earm_BBGEM_Track_Chi2true, &b_Earm_BBGEM_Track_Chi2true);
   fChain->SetBranchAddress("Earm.BBGEM.Track.X", &Earm_BBGEM_Track_X, &b_Earm_BBGEM_Track_X);
   fChain->SetBranchAddress("Earm.BBGEM.Track.Y", &Earm_BBGEM_Track_Y, &b_Earm_BBGEM_Track_Y);
   fChain->SetBranchAddress("Earm.BBGEM.Track.Xp", &Earm_BBGEM_Track_Xp, &b_Earm_BBGEM_Track_Xp);
   fChain->SetBranchAddress("Earm.BBGEM.Track.Yp", &Earm_BBGEM_Track_Yp, &b_Earm_BBGEM_Track_Yp);
   fChain->SetBranchAddress("Earm.BBGEM.Track.T", &Earm_BBGEM_Track_T, &b_Earm_BBGEM_Track_T);
   fChain->SetBranchAddress("Earm.BBGEM.Track.P", &Earm_BBGEM_Track_P, &b_Earm_BBGEM_Track_P);
   fChain->SetBranchAddress("Earm.BBGEM.Track.Sx", &Earm_BBGEM_Track_Sx, &b_Earm_BBGEM_Track_Sx);
   fChain->SetBranchAddress("Earm.BBGEM.Track.Sy", &Earm_BBGEM_Track_Sy, &b_Earm_BBGEM_Track_Sy);
   fChain->SetBranchAddress("Earm.BBGEM.Track.Sz", &Earm_BBGEM_Track_Sz, &b_Earm_BBGEM_Track_Sz);
   fChain->SetBranchAddress("Earm.BBGEM.Track.Xfit", &Earm_BBGEM_Track_Xfit, &b_Earm_BBGEM_Track_Xfit);
   fChain->SetBranchAddress("Earm.BBGEM.Track.Yfit", &Earm_BBGEM_Track_Yfit, &b_Earm_BBGEM_Track_Yfit);
   fChain->SetBranchAddress("Earm.BBGEM.Track.Xpfit", &Earm_BBGEM_Track_Xpfit, &b_Earm_BBGEM_Track_Xpfit);
   fChain->SetBranchAddress("Earm.BBGEM.Track.Ypfit", &Earm_BBGEM_Track_Ypfit, &b_Earm_BBGEM_Track_Ypfit);
   fChain->SetBranchAddress("Earm.BBGEM.Track.otridx", &Earm_BBGEM_Track_otridx, &b_Earm_BBGEM_Track_otridx);
   fChain->SetBranchAddress("Earm.BBGEM.Track.ptridx", &Earm_BBGEM_Track_ptridx, &b_Earm_BBGEM_Track_ptridx);
   fChain->SetBranchAddress("Earm.BBGEM.Track.sdtridx", &Earm_BBGEM_Track_sdtridx, &b_Earm_BBGEM_Track_sdtridx);
   fChain->SetBranchAddress("LAD.HodoScint.det.esum", &Earm_BBHodoScint_det_esum, &b_Earm_BBHodoScint_det_esum);
   fChain->SetBranchAddress("LAD.HodoScint.hit.nhits", &Earm_BBHodoScint_hit_nhits, &b_Earm_BBHodoScint_hit_nhits);
   fChain->SetBranchAddress("LAD.HodoScint.hit.row", &Earm_BBHodoScint_hit_row, &b_Earm_BBHodoScint_hit_row);
   fChain->SetBranchAddress("LAD.HodoScint.hit.col", &Earm_BBHodoScint_hit_col, &b_Earm_BBHodoScint_hit_col);
   fChain->SetBranchAddress("LAD.HodoScint.hit.cell", &Earm_BBHodoScint_hit_cell, &b_Earm_BBHodoScint_hit_cell);
   fChain->SetBranchAddress("LAD.HodoScint.hit.plane", &Earm_BBHodoScint_hit_plane, &b_Earm_BBHodoScint_hit_plane);
   fChain->SetBranchAddress("LAD.HodoScint.hit.wire", &Earm_BBHodoScint_hit_wire, &b_Earm_BBHodoScint_hit_wire);
   fChain->SetBranchAddress("LAD.HodoScint.hit.xcell", &Earm_BBHodoScint_hit_xcell, &b_Earm_BBHodoScint_hit_xcell);
   fChain->SetBranchAddress("LAD.HodoScint.hit.ycell", &Earm_BBHodoScint_hit_ycell, &b_Earm_BBHodoScint_hit_ycell);
   fChain->SetBranchAddress("LAD.HodoScint.hit.zcell", &Earm_BBHodoScint_hit_zcell, &b_Earm_BBHodoScint_hit_zcell);
   fChain->SetBranchAddress("LAD.HodoScint.hit.xcellg", &Earm_BBHodoScint_hit_xcellg, &b_Earm_BBHodoScint_hit_xcellg);
   fChain->SetBranchAddress("LAD.HodoScint.hit.ycellg", &Earm_BBHodoScint_hit_ycellg, &b_Earm_BBHodoScint_hit_ycellg);
   fChain->SetBranchAddress("LAD.HodoScint.hit.zcellg", &Earm_BBHodoScint_hit_zcellg, &b_Earm_BBHodoScint_hit_zcellg);
   fChain->SetBranchAddress("LAD.HodoScint.hit.xhit", &Earm_BBHodoScint_hit_xhit, &b_Earm_BBHodoScint_hit_xhit);
   fChain->SetBranchAddress("LAD.HodoScint.hit.yhit", &Earm_BBHodoScint_hit_yhit, &b_Earm_BBHodoScint_hit_yhit);
   fChain->SetBranchAddress("LAD.HodoScint.hit.zhit", &Earm_BBHodoScint_hit_zhit, &b_Earm_BBHodoScint_hit_zhit);
   fChain->SetBranchAddress("LAD.HodoScint.hit.xhitg", &Earm_BBHodoScint_hit_xhitg, &b_Earm_BBHodoScint_hit_xhitg);
   fChain->SetBranchAddress("LAD.HodoScint.hit.yhitg", &Earm_BBHodoScint_hit_yhitg, &b_Earm_BBHodoScint_hit_yhitg);
   fChain->SetBranchAddress("LAD.HodoScint.hit.zhitg", &Earm_BBHodoScint_hit_zhitg, &b_Earm_BBHodoScint_hit_zhitg);
   fChain->SetBranchAddress("LAD.HodoScint.hit.sumedep", &Earm_BBHodoScint_hit_sumedep, &b_Earm_BBHodoScint_hit_sumedep);
   fChain->SetBranchAddress("LAD.HodoScint.hit.tavg", &Earm_BBHodoScint_hit_tavg, &b_Earm_BBHodoScint_hit_tavg);
   fChain->SetBranchAddress("LAD.HodoScint.hit.trms", &Earm_BBHodoScint_hit_trms, &b_Earm_BBHodoScint_hit_trms);
   fChain->SetBranchAddress("LAD.HodoScint.hit.tmin", &Earm_BBHodoScint_hit_tmin, &b_Earm_BBHodoScint_hit_tmin);
   fChain->SetBranchAddress("LAD.HodoScint.hit.tmax", &Earm_BBHodoScint_hit_tmax, &b_Earm_BBHodoScint_hit_tmax);
   fChain->SetBranchAddress("LAD.HodoScint.hit.otridx", &Earm_BBHodoScint_hit_otridx, &b_Earm_BBHodoScint_hit_otridx);
   fChain->SetBranchAddress("LAD.HodoScint.hit.ptridx", &Earm_BBHodoScint_hit_ptridx, &b_Earm_BBHodoScint_hit_ptridx);
   fChain->SetBranchAddress("LAD.HodoScint.hit.sdtridx", &Earm_BBHodoScint_hit_sdtridx, &b_Earm_BBHodoScint_hit_sdtridx);
   
   fChain->SetBranchAddress("OTrack.ntracks", &OTrack_ntracks, &b_OTrack_ntracks);
   fChain->SetBranchAddress("OTrack.TID", &OTrack_TID, &b_OTrack_TID);
   fChain->SetBranchAddress("OTrack.MID", &OTrack_MID, &b_OTrack_MID);
   fChain->SetBranchAddress("OTrack.PID", &OTrack_PID, &b_OTrack_PID);
   fChain->SetBranchAddress("OTrack.posx", &OTrack_posx, &b_OTrack_posx);
   fChain->SetBranchAddress("OTrack.posy", &OTrack_posy, &b_OTrack_posy);
   fChain->SetBranchAddress("OTrack.posz", &OTrack_posz, &b_OTrack_posz);
   fChain->SetBranchAddress("OTrack.momx", &OTrack_momx, &b_OTrack_momx);
   fChain->SetBranchAddress("OTrack.momy", &OTrack_momy, &b_OTrack_momy);
   fChain->SetBranchAddress("OTrack.momz", &OTrack_momz, &b_OTrack_momz);
   fChain->SetBranchAddress("OTrack.polx", &OTrack_polx, &b_OTrack_polx);
   fChain->SetBranchAddress("OTrack.poly", &OTrack_poly, &b_OTrack_poly);
   fChain->SetBranchAddress("OTrack.polz", &OTrack_polz, &b_OTrack_polz);
   fChain->SetBranchAddress("OTrack.Etot", &OTrack_Etot, &b_OTrack_Etot);
   fChain->SetBranchAddress("OTrack.T", &OTrack_T, &b_OTrack_T);
   fChain->SetBranchAddress("PTrack.ntracks", &PTrack_ntracks, &b_PTrack_ntracks);
   fChain->SetBranchAddress("PTrack.TID", &PTrack_TID, &b_PTrack_TID);
   fChain->SetBranchAddress("PTrack.PID", &PTrack_PID, &b_PTrack_PID);
   fChain->SetBranchAddress("PTrack.posx", &PTrack_posx, &b_PTrack_posx);
   fChain->SetBranchAddress("PTrack.posy", &PTrack_posy, &b_PTrack_posy);
   fChain->SetBranchAddress("PTrack.posz", &PTrack_posz, &b_PTrack_posz);
   fChain->SetBranchAddress("PTrack.momx", &PTrack_momx, &b_PTrack_momx);
   fChain->SetBranchAddress("PTrack.momy", &PTrack_momy, &b_PTrack_momy);
   fChain->SetBranchAddress("PTrack.momz", &PTrack_momz, &b_PTrack_momz);
   fChain->SetBranchAddress("PTrack.polx", &PTrack_polx, &b_PTrack_polx);
   fChain->SetBranchAddress("PTrack.poly", &PTrack_poly, &b_PTrack_poly);
   fChain->SetBranchAddress("PTrack.polz", &PTrack_polz, &b_PTrack_polz);
   fChain->SetBranchAddress("PTrack.Etot", &PTrack_Etot, &b_PTrack_Etot);
   fChain->SetBranchAddress("PTrack.T", &PTrack_T, &b_PTrack_T);
   fChain->SetBranchAddress("SDTrack.ntracks", &SDTrack_ntracks, &b_SDTrack_ntracks);
   fChain->SetBranchAddress("SDTrack.TID", &SDTrack_TID, &b_SDTrack_TID);
   fChain->SetBranchAddress("SDTrack.MID", &SDTrack_MID, &b_SDTrack_MID);
   fChain->SetBranchAddress("SDTrack.PID", &SDTrack_PID, &b_SDTrack_PID);
   fChain->SetBranchAddress("SDTrack.posx", &SDTrack_posx, &b_SDTrack_posx);
   fChain->SetBranchAddress("SDTrack.posy", &SDTrack_posy, &b_SDTrack_posy);
   fChain->SetBranchAddress("SDTrack.posz", &SDTrack_posz, &b_SDTrack_posz);
   fChain->SetBranchAddress("SDTrack.momx", &SDTrack_momx, &b_SDTrack_momx);
   fChain->SetBranchAddress("SDTrack.momy", &SDTrack_momy, &b_SDTrack_momy);
   fChain->SetBranchAddress("SDTrack.momz", &SDTrack_momz, &b_SDTrack_momz);
   fChain->SetBranchAddress("SDTrack.polx", &SDTrack_polx, &b_SDTrack_polx);
   fChain->SetBranchAddress("SDTrack.poly", &SDTrack_poly, &b_SDTrack_poly);
   fChain->SetBranchAddress("SDTrack.polz", &SDTrack_polz, &b_SDTrack_polz);
   fChain->SetBranchAddress("SDTrack.Etot", &SDTrack_Etot, &b_SDTrack_Etot);
   fChain->SetBranchAddress("SDTrack.T", &SDTrack_T, &b_SDTrack_T);
   fChain->SetBranchAddress("SDTrack.vx", &SDTrack_vx, &b_SDTrack_vx);
   fChain->SetBranchAddress("SDTrack.vy", &SDTrack_vy, &b_SDTrack_vy);
   fChain->SetBranchAddress("SDTrack.vz", &SDTrack_vz, &b_SDTrack_vz);
   fChain->SetBranchAddress("SDTrack.vnx", &SDTrack_vnx, &b_SDTrack_vnx);
   fChain->SetBranchAddress("SDTrack.vny", &SDTrack_vny, &b_SDTrack_vny);
   fChain->SetBranchAddress("SDTrack.vnz", &SDTrack_vnz, &b_SDTrack_vnz);
   fChain->SetBranchAddress("SDTrack.vEkin", &SDTrack_vEkin, &b_SDTrack_vEkin);

   fChain->SetBranchAddress("LAD.Hodo.dighit.nchan", &LAD_Hodo_dighit_nchan, &b_LAD_Hodo_dighit_nchan);
   fChain->SetBranchAddress("LAD.Hodo.dighit.chan", &LAD_Hodo_dighit_chan, &b_LAD_Hodo_dighit_chan);
   fChain->SetBranchAddress("LAD.Hodo.dighit.adc", &LAD_Hodo_dighit_adc, &b_LAD_Hodo_dighit_adc);
   fChain->SetBranchAddress("LAD.Hodo.dighit.tdc_l", &LAD_Hodo_dighit_tdc_l, &b_LAD_Hodo_dighit_tdc_l);
   fChain->SetBranchAddress("LAD.Hodo.dighit.tdc_t", &LAD_Hodo_dighit_tdc_t, &b_LAD_Hodo_dighit_tdc_t);
   fChain->SetBranchAddress("Earm.BBGEM.dighit.nstrips", &Earm_BBGEM_dighit_nstrips, &b_Earm_BBGEM_dighit_nstrips);
   fChain->SetBranchAddress("Earm.BBGEM.dighit.module", &Earm_BBGEM_dighit_module, &b_Earm_BBGEM_dighit_module);
   fChain->SetBranchAddress("Earm.BBGEM.dighit.strip", &Earm_BBGEM_dighit_strip, &b_Earm_BBGEM_dighit_strip);
   fChain->SetBranchAddress("Earm.BBGEM.dighit.adc", &Earm_BBGEM_dighit_adc, &b_Earm_BBGEM_dighit_adc);
   fChain->SetBranchAddress("Earm.BBGEM.dighit.samp", &Earm_BBGEM_dighit_samp, &b_Earm_BBGEM_dighit_samp);
   Notify();
}

Bool_t lad_tree_digitized::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void lad_tree_digitized::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t lad_tree_digitized::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef lad_tree_digitized_cxx
