#ifndef THcLADGEMModule_h
#define THcLADGEMModule_h

#include "THaSubDetector.h"
#include "TClonesArray.h"
#include "THcLADGEMCluster.h"

#include <vector>
#include <set>
#include <map>
#include <array>
#include <deque>

namespace LADGEM {
  enum GEMaxis_t { kUaxis=0, kVaxis }; // U=X, V=Y
  enum APVmap_t { kINFN=0, kUVA_XY, kUVA_UV, kMC };
}

struct mpdmap_t {
  UInt_t crate;
  UInt_t slot;
  UInt_t mpd_id;
  UInt_t gem_id;
  UInt_t adc_id;
  UInt_t i2c;
  UInt_t pos;
  UInt_t invert;
  UInt_t axis; //needed to add axis to the decode map
  UInt_t index;
};

class THcLADGEMModule : public THaSubDetector {
 public:

  THcLADGEMModule( const char* name, const char* description="", Int_t imod=0,
		   THaDetectorBase* parent = nullptr );
  virtual ~THcLADGEMModule();
  
  virtual EStatus Init( const TDatime& run_time );
  virtual void    Clear( Option_t* opt="" );
  virtual Int_t   Decode( const THaEvData& );
  virtual Int_t   CoarseProcess( TClonesArray& tracks );
  virtual Int_t   FineProcess( TClonesArray& tracks );

  std::vector<THcLADGEMCluster> GetClusters(int axis){
    if( axis == LADGEM::kUaxis )
      return fClustersU;
    else
      return fClustersV;
  }

  Int_t GetNClusters(int axis) {
    if( axis == LADGEM::kUaxis )
      return fClustersU.size();
    else
      return fClustersV.size();
  }

  Int_t GetModuleID() { return fModuleID; }
  Int_t GetLayerNum() { return fLayer; }
  Int_t GetStripNumber( UInt_t rawstrip, UInt_t pos, UInt_t invert );

  std::vector<mpdmap_t> fMPDmap;
  std::vector<Int_t> fChanMapData;

  //Pedestal means and RMS values for all channels:
  std::vector<Double_t> fPedestalU, fPedRMSU; 
  std::vector<Double_t> fPedestalV, fPedRMSV;

  //Optional database parameters for monitoring common-mode fluctuations in pedestal runs and/or full readout events:
  std::vector<Double_t> fCommonModeMeanU; 
  std::vector<Double_t> fCommonModeMeanV;
  std::vector<Double_t> fCommonModeRMSU;
  std::vector<Double_t> fCommonModeRMSV;

  Int_t GetNStripsHitU() { return fNstrips_hitU; }
  Int_t GetNStripsHitV() { return fNstrips_hitV; }
  Int_t GetNStripsHit()  { return fNstrips_hit; }
  Int_t GetN2DHits()     { return fN2Dhits; }
  
 protected:
  
  //Constant, module-specific parameters:
  Int_t    fModuleID;
  UShort_t fLayer;
  UInt_t   fNstripsU; // Total number of strips in this module along the generic "U" axis
  UInt_t   fNstripsV; // Total number of strips in this module along the generic "V" axis

  bool fIsDecoded;

  std::string fChanMapFileName;

  std::vector<std::deque<Double_t> > fCommonModeResultContainer_by_APV;
  std::vector<Double_t> fCommonModeRollingAverage_by_APV;
  std::vector<Double_t> fCommonModeRollingRMS_by_APV;
  std::vector<UInt_t> fNeventsRollingAverage_by_APV;

  std::vector<std::deque<Double_t> > fCMbiasResultContainer_by_APV;
  std::vector<Double_t> fCommonModeOnlineBiasRollingAverage_by_APV;
  std::vector<Double_t> fCommonModeOnlineBiasRollingRMS_by_APV;
  std::vector<UInt_t> fNeventsOnlineBias_by_APV;

  LADGEM::APVmap_t fAPVmapping; //choose APV channel --> strip mapping; there are only three possible values supported for now (see SBSGEM::APVmap_t)

  std::array<std::vector<UInt_t>, 4> APVMAP;

  std::vector<double> fT0_by_APV; //T0 for each APV card (really for each MPD)
  double fTref_coarse;      //T0 of the "reference" channel (time stamp coarse - T0 of the "reference" APV)
  std::vector<double> fTcoarse_by_APV; //Time stamp coarse - T0 of this channel - Tref for each channel
  std::vector<UInt_t> fTfine_by_APV; //"Fine time stamp" by APV:
  std::vector<UInt_t> fEventCount_by_APV; //MPD event counter by APV:
  std::vector<double> fTimeStamp_ns_by_APV; //Coarse time stamp - T0 + fine time stamp % 6

  // Flags
  Bool_t   fPedestalMode;
  Bool_t   fSubtractPedBeforeCommonMode;
  Double_t fZeroSuppressRMS;
  Bool_t   fZeroSuppress;
  Bool_t   fNegSignalStudy;

  Bool_t fOnlineZeroSuppression; //this MIGHT be redundant with fZeroSuppress (or not)
  Int_t fCODA_BUILD_ALL_SAMPLES;
  Int_t fCODA_CM_ENABLED;

  Int_t fCommonModeFlag;    //default = 0 = sorting method, 1 = Danning method, 2 = histogramming method, 3 = "online" Danning-method
  Int_t fCommonModeOnlFlag; //default = 3 = Danning method during GMn, 4 = Danning method during GEn
  Int_t fPedSubFlag; //default = 0 (pedestal subtraction NOT done for full readout events). 
  Double_t fCommonModeDanningMethod_NsigmaCut;
  Int_t fSuppressFirstLast;  // Suppress strips peaking in first or last time sample:
  Int_t fUseStripTimingCuts; // Apply strip timing cuts:

  Double_t fStripTau; //time constant for strip timing fit
  Double_t fDeconv_weights[3]; //

  //Make these fixed-size arrays defined per strip axis direction, this will require some painful code changes but will improve efficiency and S/N
  // FIXME: these variables are hard-coded, and mostly used in one specific function. No point to set the default values at constructor.
  // Should set the values in DB file and load them -- SP
  Double_t fStripMaxTcut_central[2], fStripMaxTcut_width[2], fStripMaxTcut_sigma[2]; // Strip timing cuts for local maximum used to seed cluster
  Double_t fStripMaxTcut_central_deconv[2], fStripMaxTcut_width_deconv[2], fStripMaxTcut_sigma_deconv[2]; //Strip timing cuts based on deconvoluted strip time
  Double_t fStripMaxTcut_central_fit[2], fStripMaxTcut_width_fit[2], fStripMaxTcut_sigma_fit[2]; //Strip timing cuts based on strip "fit" time. 
  Double_t fStripAddTcut_width; //Time cut for adding strips to a cluster
  Double_t fStripAddCorrCoeffCut; //cut on correlation coefficient for adding neighboring strips to a cluster.

  bool fUseTSchi2cut;
  std::vector<double> fGoodStrip_TSfrac_mean;  //should have same dimension as number of APV25 time samples:
  std::vector<double> fGoodStrip_TSfrac_sigma; //
  Double_t fStripTSchi2Cut;

  Bool_t fMeasureCommonMode; //Default = false; use full-readout events to measure the common-mode mean and rms in real time
  UInt_t fNeventsCommonModeLookBack; // number of events to use for rolling average; default = 100

  Bool_t fCorrectCommonMode; //Default = false; use rolling common-mode mean and calculated common-mode values to detect a condition where negative pulses bias the common-mode down and attempt to correct it.
  UInt_t fCorrectCommonModeMinStrips;
  Double_t fCorrectCommonMode_Nsigma;

  // Strip cuts:
  //Number of strips on low and high side to reject for common-mode calculation:
  Int_t fCommonModeNstripRejectHigh; //default = 28;
  Int_t fCommonModeNstripRejectLow; //default = 28;
  Int_t fCommonModeNumIterations; //number of iterations for Danning Method: default = 3
  Int_t fCommonModeMinStripsInRange; //Minimum strips in range for Danning Method: default = 10;
  Double_t fCommonModeBinWidth_Nsigma; //Bin width for "histogramming-method" common-mode calculation, in units of the RMS in a particular APV, default = 2
  Double_t fCommonModeScanRange_Nsigma; //Scan range for common-mode histogramming method calculation, in units of RMS, +/- Nsigma about the mean, default = 4
  Double_t fCommonModeStepSize_Nsigma; //Step size in units of RMS, default = 1/5:

  //same as in MPDModule: unavoidable to duplicate this definition unless someone smarter than I cam
  //can come up with a way to do so:
  UInt_t fChan_CM_flags;  
  UInt_t fChan_TimeStamp_low;
  UInt_t fChan_TimeStamp_high;
  UInt_t fChan_MPD_EventCount;
  UInt_t fChan_MPD_Debug;

  Double_t fTrigTime; //trigger time; to be decoded once by parent class
  Double_t fMaxTrigTimeCorrection; //Maximum (absolute) correction to be applied to strip times based on trigger time (Default = 25 ns)
  Double_t fTrigTimeSlope; //Slope of GEM time versus trig time correlation (default = 1)

  UChar_t fN_APV25_CHAN;     // default: 128
  UChar_t fN_MPD_TIME_SAMP;  // default: 6 number of MPD time samples
  UShort_t fMPDMAP_ROW_SIZE; // default: 9

  Double_t fSamplePeriod; //for timing calculations: default = 24 ns for Hall A .

 //variables defining rectangular track search region constraint (NOTE: these will change event-to-event, they are NOT constant!)
  
  std::vector<Double_t> fxcmin, fxcmax;
  std::vector<Double_t> fycmin, fycmax;

  //Arrays to temporarily hold raw data from ONE APV card:
  std::vector<UInt_t> fStripAPV;
  std::vector<UInt_t> fRawStripAPV;
  std::vector<Int_t>  fRawADC_APV;
  std::vector<Double_t> fPedSubADC_APV;
  std::vector<Double_t> fCommonModeSubtractedADC_APV;

  // Online calculated common-mode array
  std::vector<Double_t> fCM_online; //size equal to fN_MPD_TIME_SAMP

  //BASIC DECODED STRIP HIT INFO:
  //By the time the information is populated here, the ADC values are already assumed to be pedestal/common-mode subtracted and/or zero-suppressed as appropriate:
  Int_t fNstrips_hit; //total Number of strips fired (after common-mode subtraction and zero suppression)
  Int_t fNstrips_hit_pos; //total Number of strips fired after positive zero suppression
  Int_t fNstrips_hit_neg; //total Number of strips fired after negative zero suppression
  Int_t fNdecoded_ADCsamples; //= fNstrips_hit * fN_MPD_TIME_SAMP
  UInt_t fNstrips_hitU; //total number of U strips fired
  UInt_t fNstrips_hitV; //total number of V strips fired
  UInt_t fNstrips_hitU_neg; //total number of U strips fired negative
  UInt_t fNstrips_hitV_neg; //total number of V strips fired negative

  // Number of strips passing basic zero suppression thresholds:
  UInt_t fNstrips_keep;
  UInt_t fNstrips_keepU;
  UInt_t fNstrips_keepV;
  //Number of strips passing "local max" thresholds:
  UInt_t fNstrips_keep_lmax;
  UInt_t fNstrips_keep_lmaxU;
  UInt_t fNstrips_keep_lmaxV;


  std::vector<UInt_t> fStrip;  //Strip index of hit (these could be "U" or "V" generalized X and Y), assumed to run from 0..N-1
 std::vector<LADGEM::GEMaxis_t>  fAxis;  //We just made our enumerated type that has two possible values, makes the code more readable (maybe)
  std::vector<std::vector<Double_t> > fADCsamples; //2D array of ADC samples by hit: Outer index runs over hits; inner index runs over ADC samples
  std::vector<std::vector<Int_t> > fRawADCsamples; //2D array of raw (non-baseline-subtracted) ADC values.
  std::vector<std::vector<Double_t> > fADCsamples_deconv; //"Deconvoluted" ADC samples

  std::vector<Double_t> fADCsums;
  std::vector<Double_t> fADCsumsDeconv; //deconvoluted strip ADC sums
  std::vector<Double_t> fStripADCavg;
  std::vector<UInt_t> fStripIsU; // is this a U strip? 0/1
  std::vector<UInt_t> fStripIsV; // is this a V strip? 0/1
  std::vector<UInt_t> fStripOnTrack; //Is this strip on any track?
  std::vector<UInt_t> fStripIsNeg; //Is this strip negative?
  std::vector<UInt_t> fStripIsNegU; //Is this strip negative?
  std::vector<UInt_t> fStripIsNegV; //Is this strip negative?
  std::vector<UInt_t> fStripIsNegOnTrack; //Is this strip negative and on a track?
  std::vector<UInt_t> fStripIsNegOnTrackU; //Is this strip negative and on a track?
  std::vector<UInt_t> fStripIsNegOnTrackV; //Is this strip negative and on a track?
  std::vector<UInt_t> fStripRaw; //Raw strip numbers on track?
  std::vector<UInt_t> fStripEvent; //strip raw info
  std::vector<UInt_t> fStripCrate; //strip raw info
  std::vector<UInt_t> fStripMPD; //strip raw info
  std::vector<UInt_t> fStripADC_ID; //strip raw info
  std::vector<Int_t> fStripTrackIndex; // If this strip is included in a cluster that ends up on a good track, we want to record the index in the track array of the track that contains this strip.
  std::vector<bool> fKeepStrip; //keep this strip?
  //std::vector<bool> fNegStrip; //Does this strip pass negative zero suppression? - neg pulse study
  //std::vector<Int_t> fStripKeep; //Strip passes timing cuts (and part of a cluster)?
  std::vector<UInt_t> fMaxSamp; //APV25 time sample with maximum ADC;
  std::vector<UInt_t> fMaxSampDeconv; //Sample with largest deconvoluted ADC value
  std::vector<UInt_t> fMaxSampDeconvCombo; //last sample of two-sample combination with largest deconvoluted ADC
  std::vector<Double_t> fADCmax; //largest ADC sample on the strip:
  std::vector<Double_t> fADCmaxDeconv; //Largest deconvoluted ADC sample
  std::vector<Double_t> fADCmaxDeconvCombo; //Largest combination of two deconvoluted samples
  std::vector<Double_t> fTmeanDeconv; //mean deconvoluted strip time
  std::vector<Double_t> fTmean; //ADC-weighted mean strip time:
  std::vector<Double_t> fTsigma; //ADC-weighted RMS deviation from the mean
  std::vector<Double_t> fStripTfit; //Dumb strip fit:
  std::vector<Double_t> fStripTdiff;  //strip time diff wrt max strip in cluster
  std::vector<Double_t> fStripTSchi2; //strip time-sample chi2 wrt "good" pulse shape
  std::vector<Double_t> fStripTSprob; //p-value associated with strip TS chi2
  std::vector<Double_t> fStripCorrCoeff; //strip correlation coeff wrt max strip in cluster 
  std::vector<Double_t> fTcorr; //Strip time with all applicable corrections; e.g., trigger time, etc.
  std::vector<UInt_t> fStrip_ENABLE_CM; //Flag to indicate whether CM was done online or offline for this strip
  std::vector<UInt_t> fStrip_CM_GOOD; //Flag to indicate whether online CM succeeded
  std::vector<UInt_t> fStrip_BUILD_ALL_SAMPLES; //Flag to indicate whether online zero suppression was enabled 

  // FIXME:
  // Cluster related variables and containers are initially here
  // We want to store them separately using GEMHit Object -- SP
  Int_t fClusteringFlag; //Which quantities to use for clustering?
  Int_t fDeconvolutionFlag; //reject strips failing basic deconvolution criteria?

  Double_t fThresholdSample; //Threshold on the (gain-matched and pedestal-subtracted) max. sample on a strip to keep that strip for clustering 
  Double_t fThresholdStripSum; //Threshold on the sum of (pedestal-subtracted) ADC samples on a strip
  Double_t fThresholdClusterSum; //Threshold on the sum of (pedestal-subtracted) ADC samples 

  Double_t fThresholdSampleDeconv; //Threshold on the max. deconbo
  Double_t fThresholdDeconvADCMaxCombo; //Threshold on the maximal sum of two consecutive deconvoluted ADC samples;
  Double_t fThresholdClusterSumDeconv;

  Double_t fADCratioSigma; // sigma of ADCV/ADCU-1
  Double_t fADCasymCut,fADCasymSigma;       // Filtering criterion for ADC X/Y (or U/V) asymmetry
  Double_t fTimeCutUVdiff,fTimeCutUVsigma;    // Filtering criterion for ADC X/Y (or U/V) time difference (this is distinct from any timing cuts relative to
  //trigger or reference timing at the individual strip level
  Double_t fTimeCutUVdiffDeconv,fTimeCutUVsigmaDeconv;
  Double_t fTimeCutUVdiffFit,fTimeCutUVsigmaFit;

  Double_t fHitTimeMean[2], fHitTimeSigma[2];
  Double_t fHitTimeMeanDeconv[2], fHitTimeSigmaDeconv[2];
  Double_t fHitTimeMeanFit[2], fHitTimeSigmaFit[2];
  Double_t fSigmaHitTimeAverageCorrected;

  //Parameters controlling cluster splitting and insignificant peak elimination based on "peak prominence" calculation
  Double_t fThresh_2ndMax_nsigma;   //Number of sigmas above noise level for minimum peak prominence in splitting overlapping clusters
  Double_t fThresh_2ndMax_fraction; //Peak prominence threshold as a fraction of peak height for splitting overlapping clusters

  UShort_t fMaxNeighborsU_totalcharge; //Only strips within +/- fMaxNeighborsU of the peak can be added to a cluster for total charge calculation
  UShort_t fMaxNeighborsV_totalcharge; //Only strips within +/- fMaxNeighborsV of the peak can be added to a cluster for total charge calculation

  //Only strips within these limits around the peak can be used for hit position reconstruction
  UShort_t fMaxNeighborsU_hitpos; 
  UShort_t fMaxNeighborsV_hitpos; 

  std::vector<Double_t> fADCsamples1D; //1D array to hold ADC samples; should end up with dimension fNstrips_hit*fN_MPD_TIME_SAMP
  std::vector<Int_t> fRawADCsamples1D;
  std::vector<Double_t> fADCsamplesDeconv1D; //1D array of deconvoluted ADC samples

  Double_t fRMS_ConversionFactor; // = sqrt(fN_MPD_TIME_SAMP);

  UInt_t fNAPVs_U; //Number of APV cards per module along "U" strip direction; this is typically 8, 10, or 12, but could be larger for U/V GEMs
  UInt_t fNAPVs_V; //Number of APV cards per module along "V" strip direction; 
  std::vector<Double_t> fUgain; // Internal "gain match" coefficients for U strips by APV card;
  std::vector<Double_t> fVgain; // Internal "gain match" coefficients for V strips by APV card;
  Double_t fModuleGain; // Module gain relative to some "Target" ADC value:

  std::vector<Double_t> fCMbiasU;
  std::vector<Double_t> fCMbiasV;
  
  double fCommonModeRange_nsigma; //default = 5


  // Clustering parameters

  Double_t fSigma_hitshape; // controll hit shape for cluster-splitting algorithm
  
  //GEOMETRICAL PARAMETERS:
  Double_t fUStripPitch;    //strip pitch along U, will virtually always be 0.4 mm
  Double_t fVStripPitch;    //strip pitch along V, will virtually always be 0.4 mm
  Double_t fUStripOffset;   //position of first U strip along the direction it measures:
  Double_t fVStripOffset;   //position of first V sttrip alogn the direction it measures:
  Double_t fUAngle;         //Angle between U strips and "X" axis of TRANSPORT coordinates;
  Double_t fVAngle;         //Angle between V strips and "X" axis of TRANSPORT coordinates;
  Double_t fPxU;            //U Strip X projection = cos( UAngle );
  Double_t fPyU;            //U Strip Y projection = sin( UAngle );
  Double_t fPxV;            //V Strip X projection = cos( VAngle );
  Double_t fPyV;            //V Strip Y projection = sin( VAngle );
  Double_t fCenter[3];      //Position center of the module in the local detector coord  

  Bool_t fIsMC;

  bool fMakeEfficiencyPlots;
  bool fEfficiencyInitialized;
  bool fMakeCommonModePlots; //diagnostic plots for offline common-mode stuff: default = false;
  bool fCommonModePlotsInitialized;
  bool fCommonModePlots_DBoverride;


  virtual Int_t ReadDatabase(const TDatime& date );
  virtual Int_t DefineVariables( EMode mode = kDefine );

  // Helper functions
  Double_t StripTSchi2(int hitindex);
  TVector2 UVtoXY( TVector2 UV );
  TVector2 XYtoUV( TVector2 XY );
  Double_t CorrCoeff(int nsamples, const std::vector<double> &Usamples, const std::vector<double> &Vsamples, int firstsample=0 );

  double   GetCommonMode( UInt_t isamp, Int_t flag, const mpdmap_t &apvinfo, UInt_t nhits=128 ); //default to "sorting" method:
  double   GetCommonModeCorrection( int isamp, const mpdmap_t &apvinfo, UInt_t &ngood, const UInt_t &nhits=128, bool fullreadout=false, Int_t flag=0 );
  double   FitStripTime( int striphitindex, double RMS=20.0 ); // "dumb" fit method 
  double   CalcFitTime( const std::vector<Double_t> &samples, double RMS=20.0 );
  void     UpdateRollingAverage( int iapv, double val, std::vector<std::deque<Double_t> > &RC, std::vector<Double_t> &AVG, std::vector<Double_t> &RMS, std::vector<UInt_t> &Nevt ); 
  void     CalcDeconvolutedSamples( const std::vector<Double_t> &ADCs, std::vector<Double_t> &DeconvADCs );
  void     SetTriggerTime( Double_t ttrig );

  void     InitAPVMAP();
  Int_t    GetChannelMap(const char* prefix, const TDatime& date);

  void FindClusters1D(LADGEM::GEMaxis_t axis);
  void Find2DHits();

  std::vector<THcLADGEMCluster> fClustersU;
  std::vector<THcLADGEMCluster> fClustersV;

  Int_t fNClus;
  Int_t fMAX2DHITS;
  Int_t fN2Dhits;

  THaDetectorBase* fParent;

  ClassDef(THcLADGEMModule,0)
    
};

#endif
