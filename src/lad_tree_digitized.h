#ifndef lad_tree_digitized_h
#define lad_tree_digitized_h

#include <TChain.h>
#include <TFile.h>
#include <TROOT.h>

// Header file for the classes stored in the TTree if any.
#include "vector"

class lad_tree_digitized {
public:
  TTree *fChain;  //! pointer to the analyzed TTree or TChain
  Int_t fCurrent; //! current Tree number in a TChain

  // Fixed size dimensions of array or collections stored in the TTree if any.

  // Declaration of leaf types
  Double_t ev_count;
  // LAD GEM True
//   Int_t Earm_BBGEM_hit_nhits;
//   std::vector<int> *Earm_BBGEM_hit_plane;
//   std::vector<int> *Earm_BBGEM_hit_strip;
//   std::vector<double> *Earm_BBGEM_hit_x;
//   std::vector<double> *Earm_BBGEM_hit_y;
//   std::vector<double> *Earm_BBGEM_hit_z;
//   std::vector<double> *Earm_BBGEM_hit_polx;
//   std::vector<double> *Earm_BBGEM_hit_poly;
//   std::vector<double> *Earm_BBGEM_hit_polz;
//   std::vector<double> *Earm_BBGEM_hit_t;
//   std::vector<double> *Earm_BBGEM_hit_trms;
//   std::vector<double> *Earm_BBGEM_hit_tmin;
//   std::vector<double> *Earm_BBGEM_hit_tmax;
//   std::vector<double> *Earm_BBGEM_hit_tx;
//   std::vector<double> *Earm_BBGEM_hit_ty;
//   std::vector<double> *Earm_BBGEM_hit_xin;
//   std::vector<double> *Earm_BBGEM_hit_yin;
//   std::vector<double> *Earm_BBGEM_hit_zin;
//   std::vector<double> *Earm_BBGEM_hit_xout;
//   std::vector<double> *Earm_BBGEM_hit_yout;
//   std::vector<double> *Earm_BBGEM_hit_zout;
//   std::vector<double> *Earm_BBGEM_hit_txp;
//   std::vector<double> *Earm_BBGEM_hit_typ;
//   std::vector<double> *Earm_BBGEM_hit_xg;
//   std::vector<double> *Earm_BBGEM_hit_yg;
//   std::vector<double> *Earm_BBGEM_hit_zg;
//   std::vector<int> *Earm_BBGEM_hit_trid;
//   std::vector<int> *Earm_BBGEM_hit_mid;
//   std::vector<int> *Earm_BBGEM_hit_pid;
//   std::vector<double> *Earm_BBGEM_hit_vx;
//   std::vector<double> *Earm_BBGEM_hit_vy;
//   std::vector<double> *Earm_BBGEM_hit_vz;
//   std::vector<double> *Earm_BBGEM_hit_p;
//   std::vector<double> *Earm_BBGEM_hit_edep;
//   std::vector<double> *Earm_BBGEM_hit_beta;
//   std::vector<int> *Earm_BBGEM_hit_otridx;
//   std::vector<int> *Earm_BBGEM_hit_ptridx;
//   std::vector<int> *Earm_BBGEM_hit_sdtridx;
//   Int_t Earm_BBGEM_Track_ntracks;
//   std::vector<int> *Earm_BBGEM_Track_TID;
//   std::vector<int> *Earm_BBGEM_Track_PID;
//   std::vector<int> *Earm_BBGEM_Track_MID;
//   std::vector<int> *Earm_BBGEM_Track_NumHits;
//   std::vector<int> *Earm_BBGEM_Track_NumPlanes;
//   std::vector<int> *Earm_BBGEM_Track_NDF;
//   std::vector<double> *Earm_BBGEM_Track_Chi2fit;
//   std::vector<double> *Earm_BBGEM_Track_Chi2true;
//   std::vector<double> *Earm_BBGEM_Track_X;
//   std::vector<double> *Earm_BBGEM_Track_Y;
//   std::vector<double> *Earm_BBGEM_Track_Xp;
//   std::vector<double> *Earm_BBGEM_Track_Yp;
//   std::vector<double> *Earm_BBGEM_Track_T;
//   std::vector<double> *Earm_BBGEM_Track_P;
//   std::vector<double> *Earm_BBGEM_Track_Sx;
//   std::vector<double> *Earm_BBGEM_Track_Sy;
//   std::vector<double> *Earm_BBGEM_Track_Sz;
//   std::vector<double> *Earm_BBGEM_Track_Xfit;
//   std::vector<double> *Earm_BBGEM_Track_Yfit;
//   std::vector<double> *Earm_BBGEM_Track_Xpfit;
//   std::vector<double> *Earm_BBGEM_Track_Ypfit;
//   std::vector<int> *Earm_BBGEM_Track_otridx;
//   std::vector<int> *Earm_BBGEM_Track_ptridx;
//   std::vector<int> *Earm_BBGEM_Track_sdtridx;

   // LAD Hodo True
//   Double_t Earm_BBHodoScint_det_esum;
//   Int_t Earm_BBHodoScint_hit_nhits;
//   std::vector<int> *Earm_BBHodoScint_hit_row;
//   std::vector<int> *Earm_BBHodoScint_hit_col;
//   std::vector<int> *Earm_BBHodoScint_hit_cell;
//   std::vector<int> *Earm_BBHodoScint_hit_plane;
//   std::vector<int> *Earm_BBHodoScint_hit_wire;
//   std::vector<double> *Earm_BBHodoScint_hit_xcell;
//   std::vector<double> *Earm_BBHodoScint_hit_ycell;
//   std::vector<double> *Earm_BBHodoScint_hit_zcell;
//   std::vector<double> *Earm_BBHodoScint_hit_xcellg;
//   std::vector<double> *Earm_BBHodoScint_hit_ycellg;
//   std::vector<double> *Earm_BBHodoScint_hit_zcellg;
//   std::vector<double> *Earm_BBHodoScint_hit_xhit;
//   std::vector<double> *Earm_BBHodoScint_hit_yhit;
//   std::vector<double> *Earm_BBHodoScint_hit_zhit;
//   std::vector<double> *Earm_BBHodoScint_hit_xhitg;
//   std::vector<double> *Earm_BBHodoScint_hit_yhitg;
//   std::vector<double> *Earm_BBHodoScint_hit_zhitg;
//   std::vector<double> *Earm_BBHodoScint_hit_sumedep;
//   std::vector<double> *Earm_BBHodoScint_hit_tavg;
//   std::vector<double> *Earm_BBHodoScint_hit_trms;
//   std::vector<double> *Earm_BBHodoScint_hit_tmin;
//   std::vector<double> *Earm_BBHodoScint_hit_tmax;
//   std::vector<int> *Earm_BBHodoScint_hit_otridx;
//   std::vector<int> *Earm_BBHodoScint_hit_ptridx;
//   std::vector<int> *Earm_BBHodoScint_hit_sdtridx;

  // LAD Hodo Digitized
  Int_t LAD_Hodo_dighit_nchan;
  std::vector<int> *LAD_Hodo_dighit_chan;
  std::vector<int> *LAD_Hodo_dighit_adc;
  std::vector<int> *LAD_Hodo_dighit_tdc_l;
  std::vector<int> *LAD_Hodo_dighit_tdc_t;

  // LAD GEM Digitized
  Int_t Earm_BBGEM_dighit_nstrips;
  std::vector<int> *Earm_BBGEM_dighit_module;
  std::vector<int> *Earm_BBGEM_dighit_strip;
  std::vector<int> *Earm_BBGEM_dighit_adc;
  std::vector<int> *Earm_BBGEM_dighit_samp;

  // List of branches
  TBranch *b_ev; //!

  // LAD GEM True
  // TBranch        *b_Earm_BBGEM_hit_nhits;   //!
  // TBranch        *b_Earm_BBGEM_hit_plane;   //!
  // TBranch        *b_Earm_BBGEM_hit_strip;   //!
  // TBranch        *b_Earm_BBGEM_hit_x;   //!
  // TBranch        *b_Earm_BBGEM_hit_y;   //!
  // TBranch        *b_Earm_BBGEM_hit_z;   //!
  // TBranch        *b_Earm_BBGEM_hit_polx;   //!
  // TBranch        *b_Earm_BBGEM_hit_poly;   //!
  // TBranch        *b_Earm_BBGEM_hit_polz;   //!
  // TBranch        *b_Earm_BBGEM_hit_t;   //!
  // TBranch        *b_Earm_BBGEM_hit_trms;   //!
  // TBranch        *b_Earm_BBGEM_hit_tmin;   //!
  // TBranch        *b_Earm_BBGEM_hit_tmax;   //!
  // TBranch        *b_Earm_BBGEM_hit_tx;   //!
  // TBranch        *b_Earm_BBGEM_hit_ty;   //!
  // TBranch        *b_Earm_BBGEM_hit_xin;   //!
  // TBranch        *b_Earm_BBGEM_hit_yin;   //!
  // TBranch        *b_Earm_BBGEM_hit_zin;   //!
  // TBranch        *b_Earm_BBGEM_hit_xout;   //!
  // TBranch        *b_Earm_BBGEM_hit_yout;   //!
  // TBranch        *b_Earm_BBGEM_hit_zout;   //!
  // TBranch        *b_Earm_BBGEM_hit_txp;   //!
  // TBranch        *b_Earm_BBGEM_hit_typ;   //!
  // TBranch        *b_Earm_BBGEM_hit_xg;   //!
  // TBranch        *b_Earm_BBGEM_hit_yg;   //!
  // TBranch        *b_Earm_BBGEM_hit_zg;   //!
  // TBranch        *b_Earm_BBGEM_hit_trid;   //!
  // TBranch        *b_Earm_BBGEM_hit_mid;   //!
  // TBranch        *b_Earm_BBGEM_hit_pid;   //!
  // TBranch        *b_Earm_BBGEM_hit_vx;   //!
  // TBranch        *b_Earm_BBGEM_hit_vy;   //!
  // TBranch        *b_Earm_BBGEM_hit_vz;   //!
  // TBranch        *b_Earm_BBGEM_hit_p;   //!
  // TBranch        *b_Earm_BBGEM_hit_edep;   //!
  // TBranch        *b_Earm_BBGEM_hit_beta;   //!
  // TBranch        *b_Earm_BBGEM_hit_otridx;   //!
  // TBranch        *b_Earm_BBGEM_hit_ptridx;   //!
  // TBranch        *b_Earm_BBGEM_hit_sdtridx;   //!
  // TBranch        *b_Earm_BBGEM_Track_ntracks;   //!
  // TBranch        *b_Earm_BBGEM_Track_TID;   //!
  // TBranch        *b_Earm_BBGEM_Track_PID;   //!
  // TBranch        *b_Earm_BBGEM_Track_MID;   //!
  // TBranch        *b_Earm_BBGEM_Track_NumHits;   //!
  // TBranch        *b_Earm_BBGEM_Track_NumPlanes;   //!
  // TBranch        *b_Earm_BBGEM_Track_NDF;   //!
  // TBranch        *b_Earm_BBGEM_Track_Chi2fit;   //!
  // TBranch        *b_Earm_BBGEM_Track_Chi2true;   //!
  // TBranch        *b_Earm_BBGEM_Track_X;   //!
  // TBranch        *b_Earm_BBGEM_Track_Y;   //!
  // TBranch        *b_Earm_BBGEM_Track_Xp;   //!
  // TBranch        *b_Earm_BBGEM_Track_Yp;   //!
  // TBranch        *b_Earm_BBGEM_Track_T;   //!
  // TBranch        *b_Earm_BBGEM_Track_P;   //!
  // TBranch        *b_Earm_BBGEM_Track_Sx;   //!
  // TBranch        *b_Earm_BBGEM_Track_Sy;   //!
  // TBranch        *b_Earm_BBGEM_Track_Sz;   //!
  // TBranch        *b_Earm_BBGEM_Track_Xfit;   //!
  // TBranch        *b_Earm_BBGEM_Track_Yfit;   //!
  // TBranch        *b_Earm_BBGEM_Track_Xpfit;   //!
  // TBranch        *b_Earm_BBGEM_Track_Ypfit;   //!
  // TBranch        *b_Earm_BBGEM_Track_otridx;   //!
  // TBranch        *b_Earm_BBGEM_Track_ptridx;   //!
  // TBranch        *b_Earm_BBGEM_Track_sdtridx;   //!

  // LAD Hodo True
  //  TBranch        *b_Earm_BBHodoScint_det_esum;   //!
  //  TBranch        *b_Earm_BBHodoScint_hit_nhits;   //!
  //  TBranch        *b_Earm_BBHodoScint_hit_row;   //!
  //  TBranch        *b_Earm_BBHodoScint_hit_col;   //!
  //  TBranch        *b_Earm_BBHodoScint_hit_cell;   //!
  //  TBranch        *b_Earm_BBHodoScint_hit_plane;   //!
  //  TBranch        *b_Earm_BBHodoScint_hit_wire;   //!
  //  TBranch        *b_Earm_BBHodoScint_hit_xcell;   //!
  //  TBranch        *b_Earm_BBHodoScint_hit_ycell;   //!
  //  TBranch        *b_Earm_BBHodoScint_hit_zcell;   //!
  //  TBranch        *b_Earm_BBHodoScint_hit_xcellg;   //!
  //  TBranch        *b_Earm_BBHodoScint_hit_ycellg;   //!
  //  TBranch        *b_Earm_BBHodoScint_hit_zcellg;   //!
  //  TBranch        *b_Earm_BBHodoScint_hit_xhit;   //!
  //  TBranch        *b_Earm_BBHodoScint_hit_yhit;   //!
  //  TBranch        *b_Earm_BBHodoScint_hit_zhit;   //!
  //  TBranch        *b_Earm_BBHodoScint_hit_xhitg;   //!
  //  TBranch        *b_Earm_BBHodoScint_hit_yhitg;   //!
  //  TBranch        *b_Earm_BBHodoScint_hit_zhitg;   //!
  //  TBranch        *b_Earm_BBHodoScint_hit_sumedep;   //!
  //  TBranch        *b_Earm_BBHodoScint_hit_tavg;   //!
  //  TBranch        *b_Earm_BBHodoScint_hit_trms;   //!
  //  TBranch        *b_Earm_BBHodoScint_hit_tmin;   //!
  //  TBranch        *b_Earm_BBHodoScint_hit_tmax;   //!
  //  TBranch        *b_Earm_BBHodoScint_hit_otridx;   //!
  //  TBranch        *b_Earm_BBHodoScint_hit_ptridx;   //!
  //  TBranch        *b_Earm_BBHodoScint_hit_sdtridx;   //!

  TBranch *b_LAD_Hodo_dighit_nchan; //!
  TBranch *b_LAD_Hodo_dighit_chan;  //!
  TBranch *b_LAD_Hodo_dighit_adc;   //!
  TBranch *b_LAD_Hodo_dighit_tdc_l; //!
  TBranch *b_LAD_Hodo_dighit_tdc_t; //!

  TBranch *b_Earm_BBGEM_dighit_nstrips; //!
  TBranch *b_Earm_BBGEM_dighit_module;  //!
  TBranch *b_Earm_BBGEM_dighit_strip;   //!
  TBranch *b_Earm_BBGEM_dighit_adc;     //!
  TBranch *b_Earm_BBGEM_dighit_samp;    //!

  lad_tree_digitized(TTree *tree = 0);
  virtual ~lad_tree_digitized();
  virtual Int_t Cut(Long64_t entry);
  virtual Int_t GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void Init(TTree *tree);
  virtual void Loop();
  virtual Bool_t Notify();
  virtual void Show(Long64_t entry = -1);
};
#endif

#ifdef lad_tree_digitized_cxx
lad_tree_digitized::lad_tree_digitized(TTree *tree) : fChain(0) {
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  if (tree == 0) {
    TFile *f =
        (TFile *)gROOT->GetListOfFiles()->FindObject("/volatile/halla/sbs/efuchey/digitized_files/simdigtest_2.root");
    if (!f || !f->IsOpen()) {
      f = new TFile("/volatile/halla/sbs/efuchey/digitized_files/simdigtest_2.root");
    }
    f->GetObject("T", tree);
  }
  Init(tree);
}

lad_tree_digitized::~lad_tree_digitized() {
  if (!fChain)
    return;
  delete fChain->GetCurrentFile();
}

Int_t lad_tree_digitized::GetEntry(Long64_t entry) {
  // Read contents of entry.
  if (!fChain)
    return 0;
  return fChain->GetEntry(entry);
}
Long64_t lad_tree_digitized::LoadTree(Long64_t entry) {
  // Set the environment to read one entry
  if (!fChain)
    return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0)
    return centry;
  if (fChain->GetTreeNumber() != fCurrent) {
    fCurrent = fChain->GetTreeNumber();
    Notify();
  }
  return centry;
}

void lad_tree_digitized::Init(TTree *tree) {
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  // Set object pointer
  // LAD GEM True
  // Earm_BBGEM_hit_plane = 0;
  // Earm_BBGEM_hit_strip = 0;
  // Earm_BBGEM_hit_x = 0;
  // Earm_BBGEM_hit_y = 0;
  // Earm_BBGEM_hit_z = 0;
  // Earm_BBGEM_hit_polx = 0;
  // Earm_BBGEM_hit_poly = 0;
  // Earm_BBGEM_hit_polz = 0;
  // Earm_BBGEM_hit_t = 0;
  // Earm_BBGEM_hit_trms = 0;
  // Earm_BBGEM_hit_tmin = 0;
  // Earm_BBGEM_hit_tmax = 0;
  // Earm_BBGEM_hit_tx = 0;
  // Earm_BBGEM_hit_ty = 0;
  // Earm_BBGEM_hit_xin = 0;
  // Earm_BBGEM_hit_yin = 0;
  // Earm_BBGEM_hit_zin = 0;
  // Earm_BBGEM_hit_xout = 0;
  // Earm_BBGEM_hit_yout = 0;
  // Earm_BBGEM_hit_zout = 0;
  // Earm_BBGEM_hit_txp = 0;
  // Earm_BBGEM_hit_typ = 0;
  // Earm_BBGEM_hit_xg = 0;
  // Earm_BBGEM_hit_yg = 0;
  // Earm_BBGEM_hit_zg = 0;
  // Earm_BBGEM_hit_trid = 0;
  // Earm_BBGEM_hit_mid = 0;
  // Earm_BBGEM_hit_pid = 0;
  // Earm_BBGEM_hit_vx = 0;
  // Earm_BBGEM_hit_vy = 0;
  // Earm_BBGEM_hit_vz = 0;
  // Earm_BBGEM_hit_p = 0;
  // Earm_BBGEM_hit_edep = 0;
  // Earm_BBGEM_hit_beta = 0;
  // Earm_BBGEM_hit_otridx = 0;
  // Earm_BBGEM_hit_ptridx = 0;
  // Earm_BBGEM_hit_sdtridx = 0;
  // Earm_BBGEM_Track_TID = 0;
  // Earm_BBGEM_Track_PID = 0;
  // Earm_BBGEM_Track_MID = 0;
  // Earm_BBGEM_Track_NumHits = 0;
  // Earm_BBGEM_Track_NumPlanes = 0;
  // Earm_BBGEM_Track_NDF = 0;
  // Earm_BBGEM_Track_Chi2fit = 0;
  // Earm_BBGEM_Track_Chi2true = 0;
  // Earm_BBGEM_Track_X = 0;
  // Earm_BBGEM_Track_Y = 0;
  // Earm_BBGEM_Track_Xp = 0;
  // Earm_BBGEM_Track_Yp = 0;
  // Earm_BBGEM_Track_T = 0;
  // Earm_BBGEM_Track_P = 0;
  // Earm_BBGEM_Track_Sx = 0;
  // Earm_BBGEM_Track_Sy = 0;
  // Earm_BBGEM_Track_Sz = 0;
  // Earm_BBGEM_Track_Xfit = 0;
  // Earm_BBGEM_Track_Yfit = 0;
  // Earm_BBGEM_Track_Xpfit = 0;
  // Earm_BBGEM_Track_Ypfit = 0;
  // Earm_BBGEM_Track_otridx = 0;
  // Earm_BBGEM_Track_ptridx = 0;
  // Earm_BBGEM_Track_sdtridx = 0;

  // LAD Hodo True
  // Earm_BBHodoScint_hit_row = 0;
  // Earm_BBHodoScint_hit_col = 0;
  // Earm_BBHodoScint_hit_cell = 0;
  // Earm_BBHodoScint_hit_plane = 0;
  // Earm_BBHodoScint_hit_wire = 0;
  // Earm_BBHodoScint_hit_xcell = 0;
  // Earm_BBHodoScint_hit_ycell = 0;
  // Earm_BBHodoScint_hit_zcell = 0;
  // Earm_BBHodoScint_hit_xcellg = 0;
  // Earm_BBHodoScint_hit_ycellg = 0;
  // Earm_BBHodoScint_hit_zcellg = 0;
  // Earm_BBHodoScint_hit_xhit = 0;
  // Earm_BBHodoScint_hit_yhit = 0;
  // Earm_BBHodoScint_hit_zhit = 0;
  // Earm_BBHodoScint_hit_xhitg = 0;
  // Earm_BBHodoScint_hit_yhitg = 0;
  // Earm_BBHodoScint_hit_zhitg = 0;
  // Earm_BBHodoScint_hit_sumedep = 0;
  // Earm_BBHodoScint_hit_tavg = 0;
  // Earm_BBHodoScint_hit_trms = 0;
  // Earm_BBHodoScint_hit_tmin = 0;
  // Earm_BBHodoScint_hit_tmax = 0;
  // Earm_BBHodoScint_hit_otridx = 0;
  // Earm_BBHodoScint_hit_ptridx = 0;
  // Earm_BBHodoScint_hit_sdtridx = 0;

  // LAD Hodo Digitized
  LAD_Hodo_dighit_chan  = 0;
  LAD_Hodo_dighit_adc   = 0;
  LAD_Hodo_dighit_tdc_l = 0;
  LAD_Hodo_dighit_tdc_t = 0;

  // LAD GEM Digitized
  // Earm_BBGEM_dighit_module = 0;
  // Earm_BBGEM_dighit_strip = 0;
  // Earm_BBGEM_dighit_adc = 0;
  // Earm_BBGEM_dighit_samp = 0;

  // Set branch addresses and branch pointers
  if (!tree)
    return;
  fChain   = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);

  fChain->SetBranchAddress("ev", &ev_count, &b_ev);
  // LAD GEM True
  // fChain->SetBranchAddress("Earm.BBGEM.hit.nhits", &Earm_BBGEM_hit_nhits, &b_Earm_BBGEM_hit_nhits);
  // fChain->SetBranchAddress("Earm.BBGEM.hit.plane", &Earm_BBGEM_hit_plane, &b_Earm_BBGEM_hit_plane);
  // fChain->SetBranchAddress("Earm.BBGEM.hit.strip", &Earm_BBGEM_hit_strip, &b_Earm_BBGEM_hit_strip);
  // fChain->SetBranchAddress("Earm.BBGEM.hit.x", &Earm_BBGEM_hit_x, &b_Earm_BBGEM_hit_x);
  // fChain->SetBranchAddress("Earm.BBGEM.hit.y", &Earm_BBGEM_hit_y, &b_Earm_BBGEM_hit_y);
  // fChain->SetBranchAddress("Earm.BBGEM.hit.z", &Earm_BBGEM_hit_z, &b_Earm_BBGEM_hit_z);
  // fChain->SetBranchAddress("Earm.BBGEM.hit.polx", &Earm_BBGEM_hit_polx, &b_Earm_BBGEM_hit_polx);
  // fChain->SetBranchAddress("Earm.BBGEM.hit.poly", &Earm_BBGEM_hit_poly, &b_Earm_BBGEM_hit_poly);
  // fChain->SetBranchAddress("Earm.BBGEM.hit.polz", &Earm_BBGEM_hit_polz, &b_Earm_BBGEM_hit_polz);
  // fChain->SetBranchAddress("Earm.BBGEM.hit.t", &Earm_BBGEM_hit_t, &b_Earm_BBGEM_hit_t);
  // fChain->SetBranchAddress("Earm.BBGEM.hit.trms", &Earm_BBGEM_hit_trms, &b_Earm_BBGEM_hit_trms);
  // fChain->SetBranchAddress("Earm.BBGEM.hit.tmin", &Earm_BBGEM_hit_tmin, &b_Earm_BBGEM_hit_tmin);
  // fChain->SetBranchAddress("Earm.BBGEM.hit.tmax", &Earm_BBGEM_hit_tmax, &b_Earm_BBGEM_hit_tmax);
  // fChain->SetBranchAddress("Earm.BBGEM.hit.tx", &Earm_BBGEM_hit_tx, &b_Earm_BBGEM_hit_tx);
  // fChain->SetBranchAddress("Earm.BBGEM.hit.ty", &Earm_BBGEM_hit_ty, &b_Earm_BBGEM_hit_ty);
  // fChain->SetBranchAddress("Earm.BBGEM.hit.xin", &Earm_BBGEM_hit_xin, &b_Earm_BBGEM_hit_xin);
  // fChain->SetBranchAddress("Earm.BBGEM.hit.yin", &Earm_BBGEM_hit_yin, &b_Earm_BBGEM_hit_yin);
  // fChain->SetBranchAddress("Earm.BBGEM.hit.zin", &Earm_BBGEM_hit_zin, &b_Earm_BBGEM_hit_zin);
  // fChain->SetBranchAddress("Earm.BBGEM.hit.xout", &Earm_BBGEM_hit_xout, &b_Earm_BBGEM_hit_xout);
  // fChain->SetBranchAddress("Earm.BBGEM.hit.yout", &Earm_BBGEM_hit_yout, &b_Earm_BBGEM_hit_yout);
  // fChain->SetBranchAddress("Earm.BBGEM.hit.zout", &Earm_BBGEM_hit_zout, &b_Earm_BBGEM_hit_zout);
  // fChain->SetBranchAddress("Earm.BBGEM.hit.txp", &Earm_BBGEM_hit_txp, &b_Earm_BBGEM_hit_txp);
  // fChain->SetBranchAddress("Earm.BBGEM.hit.typ", &Earm_BBGEM_hit_typ, &b_Earm_BBGEM_hit_typ);
  // fChain->SetBranchAddress("Earm.BBGEM.hit.xg", &Earm_BBGEM_hit_xg, &b_Earm_BBGEM_hit_xg);
  // fChain->SetBranchAddress("Earm.BBGEM.hit.yg", &Earm_BBGEM_hit_yg, &b_Earm_BBGEM_hit_yg);
  // fChain->SetBranchAddress("Earm.BBGEM.hit.zg", &Earm_BBGEM_hit_zg, &b_Earm_BBGEM_hit_zg);
  // fChain->SetBranchAddress("Earm.BBGEM.hit.trid", &Earm_BBGEM_hit_trid, &b_Earm_BBGEM_hit_trid);
  // fChain->SetBranchAddress("Earm.BBGEM.hit.mid", &Earm_BBGEM_hit_mid, &b_Earm_BBGEM_hit_mid);
  // fChain->SetBranchAddress("Earm.BBGEM.hit.pid", &Earm_BBGEM_hit_pid, &b_Earm_BBGEM_hit_pid);
  // fChain->SetBranchAddress("Earm.BBGEM.hit.vx", &Earm_BBGEM_hit_vx, &b_Earm_BBGEM_hit_vx);
  // fChain->SetBranchAddress("Earm.BBGEM.hit.vy", &Earm_BBGEM_hit_vy, &b_Earm_BBGEM_hit_vy);
  // fChain->SetBranchAddress("Earm.BBGEM.hit.vz", &Earm_BBGEM_hit_vz, &b_Earm_BBGEM_hit_vz);
  // fChain->SetBranchAddress("Earm.BBGEM.hit.p", &Earm_BBGEM_hit_p, &b_Earm_BBGEM_hit_p);
  // fChain->SetBranchAddress("Earm.BBGEM.hit.edep", &Earm_BBGEM_hit_edep, &b_Earm_BBGEM_hit_edep);
  // fChain->SetBranchAddress("Earm.BBGEM.hit.beta", &Earm_BBGEM_hit_beta, &b_Earm_BBGEM_hit_beta);
  // fChain->SetBranchAddress("Earm.BBGEM.hit.otridx", &Earm_BBGEM_hit_otridx, &b_Earm_BBGEM_hit_otridx);
  // fChain->SetBranchAddress("Earm.BBGEM.hit.ptridx", &Earm_BBGEM_hit_ptridx, &b_Earm_BBGEM_hit_ptridx);
  // fChain->SetBranchAddress("Earm.BBGEM.hit.sdtridx", &Earm_BBGEM_hit_sdtridx, &b_Earm_BBGEM_hit_sdtridx);
  // fChain->SetBranchAddress("Earm.BBGEM.Track.ntracks", &Earm_BBGEM_Track_ntracks, &b_Earm_BBGEM_Track_ntracks);
  // fChain->SetBranchAddress("Earm.BBGEM.Track.TID", &Earm_BBGEM_Track_TID, &b_Earm_BBGEM_Track_TID);
  // fChain->SetBranchAddress("Earm.BBGEM.Track.PID", &Earm_BBGEM_Track_PID, &b_Earm_BBGEM_Track_PID);
  // fChain->SetBranchAddress("Earm.BBGEM.Track.MID", &Earm_BBGEM_Track_MID, &b_Earm_BBGEM_Track_MID);
  // fChain->SetBranchAddress("Earm.BBGEM.Track.NumHits", &Earm_BBGEM_Track_NumHits, &b_Earm_BBGEM_Track_NumHits);
  // fChain->SetBranchAddress("Earm.BBGEM.Track.NumPlanes", &Earm_BBGEM_Track_NumPlanes, &b_Earm_BBGEM_Track_NumPlanes);
  // fChain->SetBranchAddress("Earm.BBGEM.Track.NDF", &Earm_BBGEM_Track_NDF, &b_Earm_BBGEM_Track_NDF);
  // fChain->SetBranchAddress("Earm.BBGEM.Track.Chi2fit", &Earm_BBGEM_Track_Chi2fit, &b_Earm_BBGEM_Track_Chi2fit);
  // fChain->SetBranchAddress("Earm.BBGEM.Track.Chi2true", &Earm_BBGEM_Track_Chi2true, &b_Earm_BBGEM_Track_Chi2true);
  // fChain->SetBranchAddress("Earm.BBGEM.Track.X", &Earm_BBGEM_Track_X, &b_Earm_BBGEM_Track_X);
  // fChain->SetBranchAddress("Earm.BBGEM.Track.Y", &Earm_BBGEM_Track_Y, &b_Earm_BBGEM_Track_Y);
  // fChain->SetBranchAddress("Earm.BBGEM.Track.Xp", &Earm_BBGEM_Track_Xp, &b_Earm_BBGEM_Track_Xp);
  // fChain->SetBranchAddress("Earm.BBGEM.Track.Yp", &Earm_BBGEM_Track_Yp, &b_Earm_BBGEM_Track_Yp);
  // fChain->SetBranchAddress("Earm.BBGEM.Track.T", &Earm_BBGEM_Track_T, &b_Earm_BBGEM_Track_T);
  // fChain->SetBranchAddress("Earm.BBGEM.Track.P", &Earm_BBGEM_Track_P, &b_Earm_BBGEM_Track_P);
  // fChain->SetBranchAddress("Earm.BBGEM.Track.Sx", &Earm_BBGEM_Track_Sx, &b_Earm_BBGEM_Track_Sx);
  // fChain->SetBranchAddress("Earm.BBGEM.Track.Sy", &Earm_BBGEM_Track_Sy, &b_Earm_BBGEM_Track_Sy);
  // fChain->SetBranchAddress("Earm.BBGEM.Track.Sz", &Earm_BBGEM_Track_Sz, &b_Earm_BBGEM_Track_Sz);
  // fChain->SetBranchAddress("Earm.BBGEM.Track.Xfit", &Earm_BBGEM_Track_Xfit, &b_Earm_BBGEM_Track_Xfit);
  // fChain->SetBranchAddress("Earm.BBGEM.Track.Yfit", &Earm_BBGEM_Track_Yfit, &b_Earm_BBGEM_Track_Yfit);
  // fChain->SetBranchAddress("Earm.BBGEM.Track.Xpfit", &Earm_BBGEM_Track_Xpfit, &b_Earm_BBGEM_Track_Xpfit);
  // fChain->SetBranchAddress("Earm.BBGEM.Track.Ypfit", &Earm_BBGEM_Track_Ypfit, &b_Earm_BBGEM_Track_Ypfit);
  // fChain->SetBranchAddress("Earm.BBGEM.Track.otridx", &Earm_BBGEM_Track_otridx, &b_Earm_BBGEM_Track_otridx);
  // fChain->SetBranchAddress("Earm.BBGEM.Track.ptridx", &Earm_BBGEM_Track_ptridx, &b_Earm_BBGEM_Track_ptridx);
  // fChain->SetBranchAddress("Earm.BBGEM.Track.sdtridx", &Earm_BBGEM_Track_sdtridx, &b_Earm_BBGEM_Track_sdtridx);
  // LAD Hodo True
  // fChain->SetBranchAddress("LAD.HodoScint.det.esum", &Earm_BBHodoScint_det_esum, &b_Earm_BBHodoScint_det_esum);
  // fChain->SetBranchAddress("LAD.HodoScint.hit.nhits", &Earm_BBHodoScint_hit_nhits, &b_Earm_BBHodoScint_hit_nhits);
  // fChain->SetBranchAddress("LAD.HodoScint.hit.row", &Earm_BBHodoScint_hit_row, &b_Earm_BBHodoScint_hit_row);
  // fChain->SetBranchAddress("LAD.HodoScint.hit.col", &Earm_BBHodoScint_hit_col, &b_Earm_BBHodoScint_hit_col);
  // fChain->SetBranchAddress("LAD.HodoScint.hit.cell", &Earm_BBHodoScint_hit_cell, &b_Earm_BBHodoScint_hit_cell);
  // fChain->SetBranchAddress("LAD.HodoScint.hit.plane", &Earm_BBHodoScint_hit_plane, &b_Earm_BBHodoScint_hit_plane);
  // fChain->SetBranchAddress("LAD.HodoScint.hit.wire", &Earm_BBHodoScint_hit_wire, &b_Earm_BBHodoScint_hit_wire);
  // fChain->SetBranchAddress("LAD.HodoScint.hit.xcell", &Earm_BBHodoScint_hit_xcell, &b_Earm_BBHodoScint_hit_xcell);
  // fChain->SetBranchAddress("LAD.HodoScint.hit.ycell", &Earm_BBHodoScint_hit_ycell, &b_Earm_BBHodoScint_hit_ycell);
  // fChain->SetBranchAddress("LAD.HodoScint.hit.zcell", &Earm_BBHodoScint_hit_zcell, &b_Earm_BBHodoScint_hit_zcell);
  // fChain->SetBranchAddress("LAD.HodoScint.hit.xcellg", &Earm_BBHodoScint_hit_xcellg, &b_Earm_BBHodoScint_hit_xcellg);
  // fChain->SetBranchAddress("LAD.HodoScint.hit.ycellg", &Earm_BBHodoScint_hit_ycellg, &b_Earm_BBHodoScint_hit_ycellg);
  // fChain->SetBranchAddress("LAD.HodoScint.hit.zcellg", &Earm_BBHodoScint_hit_zcellg, &b_Earm_BBHodoScint_hit_zcellg);
  // fChain->SetBranchAddress("LAD.HodoScint.hit.xhit", &Earm_BBHodoScint_hit_xhit, &b_Earm_BBHodoScint_hit_xhit);
  // fChain->SetBranchAddress("LAD.HodoScint.hit.yhit", &Earm_BBHodoScint_hit_yhit, &b_Earm_BBHodoScint_hit_yhit);
  // fChain->SetBranchAddress("LAD.HodoScint.hit.zhit", &Earm_BBHodoScint_hit_zhit, &b_Earm_BBHodoScint_hit_zhit);
  // fChain->SetBranchAddress("LAD.HodoScint.hit.xhitg", &Earm_BBHodoScint_hit_xhitg, &b_Earm_BBHodoScint_hit_xhitg);
  // fChain->SetBranchAddress("LAD.HodoScint.hit.yhitg", &Earm_BBHodoScint_hit_yhitg, &b_Earm_BBHodoScint_hit_yhitg);
  // fChain->SetBranchAddress("LAD.HodoScint.hit.zhitg", &Earm_BBHodoScint_hit_zhitg, &b_Earm_BBHodoScint_hit_zhitg);
  // fChain->SetBranchAddress("LAD.HodoScint.hit.sumedep", &Earm_BBHodoScint_hit_sumedep,
  // &b_Earm_BBHodoScint_hit_sumedep); fChain->SetBranchAddress("LAD.HodoScint.hit.tavg", &Earm_BBHodoScint_hit_tavg,
  // &b_Earm_BBHodoScint_hit_tavg); fChain->SetBranchAddress("LAD.HodoScint.hit.trms", &Earm_BBHodoScint_hit_trms,
  // &b_Earm_BBHodoScint_hit_trms); fChain->SetBranchAddress("LAD.HodoScint.hit.tmin", &Earm_BBHodoScint_hit_tmin,
  // &b_Earm_BBHodoScint_hit_tmin); fChain->SetBranchAddress("LAD.HodoScint.hit.tmax", &Earm_BBHodoScint_hit_tmax,
  // &b_Earm_BBHodoScint_hit_tmax); fChain->SetBranchAddress("LAD.HodoScint.hit.otridx", &Earm_BBHodoScint_hit_otridx,
  // &b_Earm_BBHodoScint_hit_otridx); fChain->SetBranchAddress("LAD.HodoScint.hit.ptridx", &Earm_BBHodoScint_hit_ptridx,
  // &b_Earm_BBHodoScint_hit_ptridx); fChain->SetBranchAddress("LAD.HodoScint.hit.sdtridx",
  // &Earm_BBHodoScint_hit_sdtridx, &b_Earm_BBHodoScint_hit_sdtridx);

  // LAD Hodo Digitized
  fChain->SetBranchAddress("LAD.Hodo.dighit.nchan", &LAD_Hodo_dighit_nchan, &b_LAD_Hodo_dighit_nchan);
  fChain->SetBranchAddress("LAD.Hodo.dighit.chan", &LAD_Hodo_dighit_chan, &b_LAD_Hodo_dighit_chan);
  fChain->SetBranchAddress("LAD.Hodo.dighit.adc", &LAD_Hodo_dighit_adc, &b_LAD_Hodo_dighit_adc);
  fChain->SetBranchAddress("LAD.Hodo.dighit.tdc_l", &LAD_Hodo_dighit_tdc_l, &b_LAD_Hodo_dighit_tdc_l);
  fChain->SetBranchAddress("LAD.Hodo.dighit.tdc_t", &LAD_Hodo_dighit_tdc_t, &b_LAD_Hodo_dighit_tdc_t);
  // LAD GEM Digitized
  //  fChain->SetBranchAddress("Earm.BBGEM.dighit.nstrips", &Earm_BBGEM_dighit_nstrips, &b_Earm_BBGEM_dighit_nstrips);
  //  fChain->SetBranchAddress("Earm.BBGEM.dighit.module", &Earm_BBGEM_dighit_module, &b_Earm_BBGEM_dighit_module);
  //  fChain->SetBranchAddress("Earm.BBGEM.dighit.strip", &Earm_BBGEM_dighit_strip, &b_Earm_BBGEM_dighit_strip);
  //  fChain->SetBranchAddress("Earm.BBGEM.dighit.adc", &Earm_BBGEM_dighit_adc, &b_Earm_BBGEM_dighit_adc);
  //  fChain->SetBranchAddress("Earm.BBGEM.dighit.samp", &Earm_BBGEM_dighit_samp, &b_Earm_BBGEM_dighit_samp);
  Notify();
}

Bool_t lad_tree_digitized::Notify() {
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}

void lad_tree_digitized::Show(Long64_t entry) {
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain)
    return;
  fChain->Show(entry);
}
Int_t lad_tree_digitized::Cut(Long64_t entry) {
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}
#endif // #ifdef lad_tree_digitized_cxx
