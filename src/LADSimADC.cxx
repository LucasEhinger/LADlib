//////////////////////////////////////////////////////////////////
//
//   LADSimADC

#include "LADSimADC.h"
#include "THaEvData.h"
#include "THaSlotData.h"
#include "TMath.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring> // for memset
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>
#include <stdint.h>
#include <string>
#include <unistd.h>

using namespace std;

// #define DEBUG
// #define WITH_DEBUG

namespace Decoder {

Module::TypeIter_t LADSimADC::fgThisType =
    // DoRegister( ModuleType( "Decoder::LADSimADC" , 50250 ));
    DoRegister(ModuleType("Decoder::LADSimADC", -250));
Module::TypeIter_t LADSimADC::fgType1 = DoRegister(ModuleType("Decoder::LADSimADC", -1881));
Module::TypeIter_t LADSimADC::fgType2 = DoRegister(ModuleType("Decoder::LADSimADC", -792));
Module::TypeIter_t LADSimADC::fgType3 = DoRegister(ModuleType("Decoder::LADSimADC", -3561));
// Int_t modid[4] = {-250, -1881, -792, -3561};

LADSimADC::LADSimADC() { sadc_data.resize(NADCCHAN); }

LADSimADC::LADSimADC(Int_t crate, Int_t slot)
    : PipeliningModule(crate, slot) // EPAF: I think we need that don't we ?
{
  sadc_data.resize(NADCCHAN);
  IsInit = kFALSE;
  Init();
}

LADSimADC::~LADSimADC() {
#if defined DEBUG && defined WITH_DEBUG
  // delete fDebugFile; fDebugFile = 0;
#endif
}

/*
Bool_t LADSimADC::HasCapability(Decoder::EModuleType type) {
  return ( type == kSampleADC || type == kPulseIntegral || type == kPulseTime
      || type == kPulsePeak || type == kPulsePedestal || type == kCoarseTime || type == kFineTime);
} */

// Clear all data vectors
void LADSimADC::ClearDataVectors() {
  // Clear all data objects
  assert(sadc_data.size() == NADCCHAN); // Initialization error in constructor
  for (uint32_t i = 0; i < NADCCHAN; i++) {
    sadc_data[i].integral = 0;
    sadc_data[i].peak_amp = 0;
    sadc_data[i].samples.clear();
  }
}

void LADSimADC::Clear(const Option_t *opt) {
  // Clear event-by-event data
  ClearDataVectors();
}

void LADSimADC::Init() {
  Clear();
  IsInit = kTRUE;
  fName  = "LADSimADC (Simple JLab Flash ADC Simulated Module)";
}

void LADSimADC::CheckDecoderStatus() const { std::cout << "LADSimADC has been called" << std::endl; }

Bool_t LADSimADC::HasCapability(Decoder::EModuleType type) {
  return (type == kSampleADC || type == kPulseIntegral || type == kPulseTime || type == kPulsePeak ||
          type == kPulsePedestal);
}

UInt_t LADSimADC::GetData(Decoder::EModuleType emode, UInt_t chan, UInt_t ievent) const {
  switch (emode) {
  case kSampleADC:
    return 0;
    return GetPulseSamplesData(chan, ievent);
  case kPulseIntegral:
    return GetPulseIntegralData(chan, ievent);
  case kPulseTime:
    return GetPulseTimeData(chan, ievent);
  case kPulsePeak:
    return GetPulsePeakData(chan, ievent);
  case kPulsePedestal:
    return GetPulsePedestalData(chan, ievent);
  }
  return 0;
}

UInt_t LADSimADC::GetPulsePeakData(UInt_t chan, UInt_t ievent) const {
  if (ievent >= sadc_data[chan].samples.size()) {
    std::cerr << "Error: LADSimADC::GetPulsePeakData: ievent " << ievent << " out of range for channel " << chan
              << std::endl;
    return 0;
  }
  return sadc_data[chan].peak_amp;
}

UInt_t LADSimADC::GetPulsePedestalData(UInt_t chan, UInt_t ievent) const {
  if (ievent >= sadc_data[chan].samples.size()) {
    std::cerr << "Error: LADSimADC::GetPulsePedestalData: ievent " << ievent << " out of range for channel " << chan
              << std::endl;
    return 0;
  }
  return 0; // need to return pedestal
}

UInt_t LADSimADC::GetPulseSamplesData(UInt_t chan, UInt_t ievent) const {
  if (ievent >= sadc_data[chan].samples.size()) {
    std::cerr << "Error: LADSimADC::GetPulseSamplesData: ievent " << ievent << " out of range for channel " << chan
              << std::endl;
    return 0;
  }
  return sadc_data[chan].samples[ievent];
}

UInt_t LADSimADC::GetPulseIntegralData(UInt_t chan, UInt_t ievent) const {
  if (ievent >= sadc_data[chan].samples.size()) {
    std::cerr << "Error: LADSimADC::GetPulseIntegralData: ievent " << ievent << " out of range for channel " << chan
              << std::endl;
    return 0;
  }
  return sadc_data[chan].integral;
}

UInt_t LADSimADC::GetPulseTimeData(UInt_t chan, UInt_t ievent) const {
  if (ievent >= sadc_data[chan].samples.size()) {
    std::cerr << "Error: LADSimADC::GetPulseTimeData: ievent " << ievent << " out of range for channel " << chan
              << std::endl;
    return 0;
  }
  return 0; // need to return time
}

UInt_t LADSimADC::GetNumEvents(Decoder::EModuleType emode, UInt_t chan) const {
  switch (emode) {
  case kSampleADC:
    return sadc_data[chan].samples.size();
  case kPulseIntegral:
    return sadc_data[chan].samples.size();
  case kPulseTime:
    return sadc_data[chan].samples.size();
  case kPulsePeak:
    return sadc_data[chan].samples.size();
  case kPulsePedestal:
    return sadc_data[chan].samples.size();
  }
  return 0;
}

UInt_t LADSimADC::LoadSlot(THaSlotData *sldat, const UInt_t *evbuffer, const UInt_t *pstop) {
  Clear();
  unsigned int nwords = 0;
  unsigned short chan = 0, type;
  UInt_t raw_buff, strip;
  bool printed  = false;
  bool is_first = true;
  // std::cout << "LADSimADC load crate/slot: " << sldat->getCrate() << "/" << sldat->getSlot() << std::endl;
  while (evbuffer < pstop) {
    // First, decode the header
    LADSimDataDecoder::DecodeHeader(*evbuffer++, type, chan, nwords);
    // std::cout << type << " " << sldat->getCrate() << " " << sldat->getSlot() << " " << chan << " " << nwords << endl;
    LADSimDataDecoder *enc = LADSimDataDecoder::GetEncoder(type);
    if (!enc) {
      std::cerr << "Could not find ADC decoder of type: " << type << ", is_first: " << is_first << std::endl;
    } else {
      if (!enc->IsADC() && !enc->IsMPD()) {
        std::cerr << "Encoder " << enc->GetName() << " of type " << type << " is not an ADC nor an MPD!" << std::endl;
      } else if (nwords > 0) {
        // if(enc->IsFADC() || enc->IsMPD()) { // FADC with samples
        if (enc->IsSADC()) { // FADC with samples
          SimEncoder::sadc_data tmp_sadc_data;
          // enc->DecodeFADC(tmp_sadc_data,evbuffer,nwords);
          enc->DecodeSADC(tmp_sadc_data, evbuffer, nwords);
          // std::cout << tmp_sadc_data.samples.size() << std::endl;
          // std::cout << chan << " " << sadc_data[chan].samples.size() << std::endl;
                      //LHE: temp fix to add pulse amp and integral. Should later fix to run as FADC
          sadc_data[chan].integral = tmp_sadc_data.integral;
          sadc_data[chan].peak_amp = tmp_sadc_data.peak_amp;
          for (size_t i = 0; i < tmp_sadc_data.samples.size(); i++) {
            raw_buff = tmp_sadc_data.samples[i];
            // std::cout << i << " " << tmp_sadc_data.samples[i] << endl;
            sadc_data[chan].samples.push_back(tmp_sadc_data.samples[i]);
            // std::cout << i << " " << sadc_data[chan].samples.size() << " " << raw_buff << endl;


            sldat->loadData("adc", chan, raw_buff, raw_buff);
          }
        } else if (enc->IsMPD()) {
          SimEncoder::mpd_data tmp_mpd_data;
          enc->DecodeMPD(tmp_mpd_data, evbuffer, nwords);
          for (size_t i = 0; i < tmp_mpd_data.samples.size(); i++) {
            raw_buff = tmp_mpd_data.samples[i];
            strip    = tmp_mpd_data.strips[i];
            // std::cout << i << " " << chan << " " << tmp_mpd_data.strips[i] << " " << tmp_mpd_data.samples[i] << endl;
            sldat->loadData("adc", chan, raw_buff, strip);
          }
        } else if (enc->IsADC()) { // Integral of ADC
          SimEncoder::adc_data tmp_adc_data;
          enc->DecodeADC(tmp_adc_data, evbuffer, nwords);
          raw_buff = tmp_adc_data.integral;
          // std::cout << " raw_buff " << raw_buff << endl;
          sadc_data[chan].integral = raw_buff; // LHE: Added 9/17/24 to add pulse amp compatability
          // sadc_data[chan].integral = tmp_adc_data.integral;
          // sadc_data[chan].peak_amp = tmp_adc_data.peak_amp;
          sldat->loadData("adc", chan, raw_buff, raw_buff);
        }
      }
    }
    evbuffer += nwords; // Skip ahead the number of words processed
    is_first = false;
  }
  if (printed)
    std::cerr << std::endl;
  return 0;
}

UInt_t LADSimADC::LoadSlot(THaSlotData *sldat, const UInt_t *evbuffer, UInt_t pos, UInt_t len) {
  return LoadSlot(sldat, evbuffer + pos, evbuffer + pos + len);
  // return LADSimADC::LoadSlot(sldat,evbuffer,len);
}

} // namespace Decoder

ClassImp(Decoder::LADSimADC)
