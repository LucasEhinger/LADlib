//////////////////////////////////////////////////////////////////
//
//   LADSimTDC

#include "LADSimTDC.h"
#include "LADSimDataDecoder.h"
#include "THaEvData.h"
#include "THaSlotData.h"
#include "TMath.h"

#include <unistd.h>
#include <stdint.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>  // for memset
#include <iostream>
#include <string>
#include <sstream>
#include <cmath>
#include <iomanip>
#include <numeric>
#include <cassert>

using namespace std;

//#define DEBUG
//#define WITH_DEBUG

namespace Decoder {

  Module::TypeIter_t LADSimTDC::fgThisType =
    //DoRegister( ModuleType( "Decoder::LADSimTDC" , 53204 ));
    DoRegister( ModuleType( "Decoder::LADSimTDC" , -3204 ));
  
  Module::TypeIter_t LADSimTDC::fgType0 =
    DoRegister( ModuleType( "Decoder::LADSimTDC" , -6401 ));
  Module::TypeIter_t LADSimTDC::fgType1 =
    DoRegister( ModuleType( "Decoder::LADSimTDC" , -1877 ));
  Module::TypeIter_t LADSimTDC::fgType2 =
    DoRegister( ModuleType( "Decoder::LADSimTDC" , -1190 ));
  Module::TypeIter_t LADSimTDC::fgType3 =
    DoRegister( ModuleType( "Decoder::LADSimTDC" , -3201 ));

  LADSimTDC::LADSimTDC()
  {
    tdc_data.resize(NTDCCHAN);
  }

  LADSimTDC::LADSimTDC(Int_t crate, Int_t slot)
    : PipeliningModule(crate, slot)//EPAF: I think we need that don't we ?
  {
    tdc_data.resize(NTDCCHAN);
    IsInit = kFALSE;
    Init();
  }

  LADSimTDC::~LADSimTDC() {
#if defined DEBUG && defined WITH_DEBUG
    // delete fDebugFile; fDebugFile = 0;
#endif
  }

  /*
  Bool_t LADSimTDC::HasCapability(Decoder::EModuleType type) {
    return ( type == kSampleADC || type == kPulseIntegral || type == kPulseTime
        || type == kPulsePeak || type == kPulsePedestal || type == kCoarseTime || type == kFineTime);
  } */

  // Clear all data vectors
  void LADSimTDC::ClearDataVectors() {
    // Clear all data objects
    assert(tdc_data.size() == NTDCCHAN);  // Initialization error in constructor
    for (uint32_t i = 0; i < NTDCCHAN; i++) {
      tdc_data[i].lead_time.clear();
      tdc_data[i].trail_time.clear();
    }
  }

  void LADSimTDC::Clear( const Option_t* opt) {
    // Clear event-by-event data
    ClearDataVectors();
  }

  void LADSimTDC::Init() {
    Clear();
    IsInit = kTRUE;
    fName = "LADSimTDC Generic TDC Module";
  }

  void LADSimTDC::CheckDecoderStatus() const {
    std::cout << "LADSimTDC has been called" << std::endl;
  }

  UInt_t LADSimTDC::LoadSlot(THaSlotData *sldat, const UInt_t *evbuffer,
      const UInt_t *pstop) {
    Clear();
    unsigned int nwords = 0;
    unsigned short chan = 0, type = 0;
    UInt_t raw_buff;
    SimEncoder::tdc_data tmp_tdc_data;
    //std::cout << "LADSimTDC load crate/slot: " << sldat->getCrate() << "/" << sldat->getSlot() << std::endl;
    while(evbuffer < pstop) {
      // First, decode the header
      chan = type = nwords = 0;
      LADSimDataDecoder::DecodeHeader(*evbuffer++,type,chan,nwords);
      //std::cout << " nwords " << nwords << " chan " << chan << " type " << type << std::endl; 
      LADSimDataDecoder *enc = LADSimDataDecoder::GetEncoder(type);
      if(!enc) {
        std::cerr << "Could not find TDC decoder of type: " << type
          << std::endl;
      } else {
        if(!enc->IsTDC()) {
          std::cerr << "Encoder " << enc->GetName() << " of type " << type
            << " is not a TDC!" << std::endl;
        } else if ( nwords > 0 ) {
	  enc->DecodeTDC(tmp_tdc_data,evbuffer,nwords);
          //std::cerr << "Got TDC encoder for type: " << type
	  //<< ", name: " << enc->GetName() << std::endl;
	  //cout << "time array size ? " << tmp_tdc_data.time.size() << endl;
          for(size_t i = 0; i < tmp_tdc_data.time.size(); i++ ) {
            raw_buff = tmp_tdc_data.getTime(i);
	    //cout << " chan " << chan << " raw_buff " << raw_buff << "; edge ? " <<  tmp_tdc_data.getEdge(i) << endl;
            if(tmp_tdc_data.getEdge(i)) { // Trail
              tdc_data[chan].trail_time.push_back(raw_buff);
            } else { // Lead
	      tdc_data[chan].lead_time.push_back(raw_buff);
            }
            // TODO: Figure out what to do with the edge information
            // I'd imagine we need to distinguish it somehow!
	    sldat->loadData("tdc",chan,raw_buff,tmp_tdc_data.getEdge(i));
          }
          tmp_tdc_data.time.clear(); // Clear it to prepare for next read
        }
      }
      evbuffer += nwords; // Skip ahead in the buffer
    }
   return 0;
  }

  UInt_t LADSimTDC::LoadSlot( THaSlotData *sldat, const UInt_t *evbuffer,
                              UInt_t pos, UInt_t len) {
    return LoadSlot(sldat,evbuffer+pos,evbuffer+pos+len);
    //return LADSimTDC::LoadSlot(sldat,evbuffer,len);
  }
  
}

ClassImp(Decoder::LADSimTDC)
