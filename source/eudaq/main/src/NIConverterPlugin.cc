#include "eudaq/DataConverterPlugin.hh"
#include "eudaq/Exception.hh"
#include "eudaq/RawDataEvent.hh"

#if USE_LCIO
#  include "IMPL/LCEventImpl.h"
#  include "IMPL/TrackerDataImpl.h"
#  include "IMPL/LCCollectionVec.h"
#  include "UTIL/CellIDEncoder.h"
#  include "lcio.h"
#endif

#include <iostream>
#include <ostream>
#include <string>
#include <vector>
#include <memory>
#include <iomanip>
#include <algorithm>

#define GET(d, i) getlittleendian<unsigned>(&(d)[(i)*4])

namespace eudaq {
  static const int dbg = 0; // 0=off, 1=structure, 2=structure+data
  static const int PIVOTPIXELOFFSET = 64;
  
  class NIConverterPlugin : public DataConverterPlugin {
    typedef std::vector<unsigned char> datavect;
    typedef std::vector<unsigned char>::const_iterator datait;
  public:
    virtual void Initialize(const Event & bore) {
      m_boards = from_string(bore.GetTag("BOARDS"), 0);
      m_ids.clear();
      for (unsigned i = 0; i < m_boards; ++i) {
        unsigned id = from_string(bore.GetTag("ID" + to_string(i)), i);
        m_ids.push_back(id);
      }
    }
    
    void DecodeFrame(lcio::TrackerDataImpl* zsFrame, size_t len, datait it, int frame, unsigned pivotpixel) const {
      std::vector<unsigned short> vec;
      for (size_t i = 0; i < len; ++i) {
        unsigned v = GET(it, i);
        vec.push_back(v & 0xffff);
        vec.push_back(v >> 16);
      }
      
      unsigned npixels = 0;
      for (size_t i = 0; i < vec.size(); ++i) {
        //  std::cout << "  " << i << " : " << hexdec(vec[i]) << std::endl;
        if (i == vec.size() - 1) break;
        unsigned numstates = vec[i] & 0xf;
        unsigned row = vec[i]>>4 & 0x7ff;
        if (numstates+1 > vec.size()-i) {
            // Ignoring bad line
            //std::cout << "Ignoring bad line " << row << " (too many states)" << std::endl;
            break;
        }
        if (dbg>1) std::cout << "Hit line " << (vec[i] & 0x8000 ? "* " : ". ") << row
                             << ", states " << numstates << ":";
        bool pivot = (row >= (pivotpixel / 16));
        for (unsigned s = 0; s < numstates; ++s) {
          unsigned v = vec.at(++i);
          unsigned column = v>>2 & 0x7ff;
          unsigned num = v & 3;
          if (dbg>1) std::cout << (s ? "," : " ") << column;
          if (dbg>1) if ((v&3) > 0) std::cout << "-" << (column + num);
          for (unsigned j = 0; j < num+1; ++j) {
            zsFrame->chargeValues().push_back( column+j );
            zsFrame->chargeValues().push_back( row );
            zsFrame->chargeValues().push_back( 1 );       
          }
          npixels += num + 1;
        }
        if (dbg>1) std::cout << std::endl;
      }
      if (dbg) std::cout << "Total pixels " << frame << " = " << npixels << std::endl;
    }
    
    virtual unsigned GetTriggerID(Event const & ev) const {
      const RawDataEvent & rawev = dynamic_cast<const RawDataEvent &>(ev);
      if (rawev.NumBlocks() < 1 || rawev.GetBlock(0).size() < 8) return (unsigned)-1;
      return GET(rawev.GetBlock(0), 1) >> 16;
    }
       
#if USE_LCIO 
    virtual void GetLCIORunHeader(lcio::LCRunHeader & /*header*/, eudaq::Event const & /*bore*/) const;
    virtual bool GetLCIOSubEvent(lcio::LCEvent & result, const Event & source) const;
#endif
  
  protected:
    
  private:
    NIConverterPlugin() : DataConverterPlugin("NI"), m_boards(0) {}
    unsigned m_boards;
    std::vector<int> m_ids;
    static NIConverterPlugin const m_instance;
  };
  
  NIConverterPlugin const NIConverterPlugin::m_instance;
  
#if USE_LCIO 
  void NIConverterPlugin::GetLCIORunHeader(lcio::LCRunHeader & header, eudaq::Event const & /*bore*/) const { }

  bool NIConverterPlugin::GetLCIOSubEvent(lcio::LCEvent & result, const Event & source) const {
    if (source.IsBORE()) {
      // shouldn't happen 
      return true;
    } else if (source.IsEORE()) {
      // nothing to do
      return true;
    }
    
    // prepare the collections for the zs ones
    LCCollectionVec * zs2DataCollection;
    bool zs2DataCollectionExists = false;
    
    try {
      zs2DataCollection = static_cast< LCCollectionVec* > ( result.getCollection( "zsdata_m26" ) );
      zs2DataCollectionExists = true;
    } catch ( lcio::DataNotAvailableException& e ) {
      zs2DataCollection = new LCCollectionVec( lcio::LCIO::TRACKERDATA );
    }
    
    // set the proper cell encoder
    CellIDEncoder<TrackerDataImpl> zs2DataEncoder( "sensorID:6,sparsePixelType:5", zs2DataCollection  );
     
    // to understand if we have problem with de-syncronisation, let
    // me prepare a vector of size_t to contain the
    // pivot pixel position
    std::vector<size_t > pivotPixelPosVec;
    
    
    // If we get here it must be a data event
    const RawDataEvent & rawev = dynamic_cast<const RawDataEvent &>(source);
    if (rawev.NumBlocks() != 2 || rawev.GetBlock(0).size() < 20 ||
          rawev.GetBlock(1).size() < 20) {
        std::cout << "Ignoring bad event " + to_string(source.GetEventNumber()) << std::endl;
        return false;
    }
    const datavect & data0 = rawev.GetBlock(0);
    const datavect & data1 = rawev.GetBlock(1);
    unsigned header0 = GET(data0, 0);
    unsigned header1 = GET(data1, 0);
    unsigned tluid = GetTriggerID(source);
    if (dbg) std::cout << "TLU id = " << hexdec(tluid, 4) << std::endl;
    unsigned pivot = GET(data0, 1) & 0xffff;
    if (dbg) std::cout << "Pivot = " << hexdec(pivot, 4) << std::endl;
    datait it0 = data0.begin() + 8;
    datait it1 = data1.begin() + 8;
    unsigned board = 0;
    while (it0 < data0.end() && it1 < data1.end()) {
      unsigned id = board;
      if (id < m_ids.size()) id = m_ids[id];
      if (dbg) std::cout << "Sensor " << board << ", id = " << id << std::endl;
      if (dbg) std::cout << "Mimosa_header0 = " << hexdec(header0) << std::endl;
      if (dbg) std::cout << "Mimosa_header1 = " << hexdec(header1) << std::endl;
      if (it0 + 2 >= data0.end()) {
        std::cout << "Trailing rubbish in first frame" << std::endl;
        break;
      }
      if (it1 + 2 >= data1.end()) {
        std::cout << "Trailing rubbish in second frame" << std::endl;
        break;
      }
      if (dbg) std::cout << "Mimosa_framecount0 = " << hexdec(GET(it0, 0)) << std::endl;
      if (dbg) std::cout << "Mimosa_framecount1 = " << hexdec(GET(it1, 0)) << std::endl;
      unsigned len0 = GET(it0, 1);
      if (dbg) std::cout << "Mimosa_wordcount0 = " << hexdec(len0 & 0xffff, 4)
                         << ", " << hexdec(len0 >> 16, 4) << std::endl;
      if ((len0 & 0xffff) != (len0 >> 16)) {
        std::cout << "Mismatched lengths decoding first frame (" +
                   to_string(len0 & 0xffff) + ", " + to_string(len0 >> 16) + ")" << std::endl;
        len0 = std::max(len0 & 0xffff, len0 >> 16);
      }
      len0 &= 0xffff;
      unsigned len1 = GET(it1, 1);
      if (dbg) std::cout << "Mimosa_wordcount1 = " << hexdec(len1 & 0xffff, 4)
                         << ", " << hexdec(len1 >> 16, 4) << std::endl;
      if ((len1 & 0xffff) != (len1 >> 16)) {
        std::cout << "Mismatched lengths decoding second frame (" +
                   to_string(len1 & 0xffff) + ", " + to_string(len1 >> 16) + ")" << std::endl;
        len1 = std::max(len1 & 0xffff, len1 >> 16);
      }
      len1 &= 0xffff;
      if (it0 + len0*4 + 12 > data0.end()) {
        std::cout << "Bad length in first frame" << std::endl;
        break;
      }
      if (it1 + len1*4 + 12 > data1.end()) {
        std::cout << "Bad length in second frame" << std::endl;
        break;
      }

      unsigned pivotpixel = (9216 + pivot + PIVOTPIXELOFFSET) % 9216;
      pivotPixelPosVec.push_back( pivotpixel );
       
      // Prepare a new lcio::TrackerData for the ZS data
      lcio::TrackerDataImpl* zsFrame =  new lcio::TrackerDataImpl;
      zs2DataEncoder["sensorID"] = id;
      zs2DataEncoder["sparsePixelType"] = 0;
      zs2DataEncoder.setCellID( zsFrame );
      
      // Fill ZS data into new lcio::TrackerData object
      DecodeFrame(zsFrame, len0, it0+8, 0, pivotpixel);
      DecodeFrame(zsFrame, len1, it1+8, 1, pivotpixel);
      
      // Now add the TrackerData to the collection
      zs2DataCollection->push_back( zsFrame );
      
      if (dbg) std::cout << "Mimosa_trailer0 = " << hexdec(GET(it0, len0+2)) << std::endl;
      if (dbg) std::cout << "Mimosa_trailer1 = " << hexdec(GET(it1, len1+2)) << std::endl;
      it0 += len0*4 + 16;
      it1 += len1*4 + 16;
      if (it0 <= data0.end()) header0 = GET(it0, -1);
      if (it1 <= data1.end()) header1 = GET(it1, -1);
      ++board;
    }

    
    // check if all the boards where running in synchronous mode or
    // not. Remember that the last pivot pixel entry is the one of the
    // master board.
    bool outOfSyncFlag = false;
    std::vector<size_t >::iterator masterBoardPivotAddress = pivotPixelPosVec.end() - 1;
    std::vector<size_t >::iterator slaveBoardPivotAddress  = pivotPixelPosVec.begin();
    while ( slaveBoardPivotAddress < masterBoardPivotAddress ) {
      if ( *slaveBoardPivotAddress - *masterBoardPivotAddress >= 2 ) {
        outOfSyncFlag = true;

        // we don't need to continue looping over all boards if one of
        // them is already out of sync
        break;
      }
      ++slaveBoardPivotAddress;
    }
    if ( outOfSyncFlag ) {

      if ( result.getEventNumber()  < 20 ) {
        // in this case we have the responsibility to tell the user that
        // the event was out of sync
        std::cout << "Event number " << result.getEventNumber() << " seems to be out of sync" << std::endl;
        std::vector<size_t >::iterator masterBoardPivotAddress = pivotPixelPosVec.end() - 1;
        std::vector<size_t >::iterator slaveBoardPivotAddress  = pivotPixelPosVec.begin();
        while ( slaveBoardPivotAddress < masterBoardPivotAddress ) {
          // print out all the slave boards first
          std::cout << " --> Board (S) " <<  std::setw(3) << setiosflags( std::ios::right )
                    << slaveBoardPivotAddress - pivotPixelPosVec.begin() << resetiosflags( std::ios::right )
                    << " = " << std::setw(15) << setiosflags( std::ios::right )
                    << (*slaveBoardPivotAddress) << resetiosflags( std::ios::right )
                    << " (" << std::setw(15) << setiosflags( std::ios::right )
                    << (signed) (*masterBoardPivotAddress) - (signed) (*slaveBoardPivotAddress) << resetiosflags( std::ios::right)
                    << ")" << std::endl;
          ++slaveBoardPivotAddress;
        }
        // print out also the master. It is impossible that the master
        // is out of sync with respect to itself, but for completeness...
        std::cout  << " --> Board (M) "  <<  std::setw(3) << setiosflags( std::ios::right )
                   << slaveBoardPivotAddress - pivotPixelPosVec.begin() << resetiosflags( std::ios::right )
                   << " = " << std::setw(15) << setiosflags( std::ios::right )
                   << (*slaveBoardPivotAddress) << resetiosflags( std::ios::right )
                   << " (" << std::setw(15)  << setiosflags( std::ios::right )
                   << (signed) (*masterBoardPivotAddress) - (signed) (*slaveBoardPivotAddress) << resetiosflags( std::ios::right)
                   << ")" << std::endl;

      } else if ( result.getEventNumber()  == 20 ) {
        // if the number of consecutive warnings is equal to the maximum
        // allowed, don't bother the user anymore with this message,
        // because it's very likely the run was taken unsynchronized on
        // purpose
        std::cout << "The maximum number of unsychronized events has been reached." << std::endl
                  << "Assuming the run was taken in asynchronous mode" << std::endl;
      }
    }
     
    if ( !zs2DataCollectionExists && ( zs2DataCollection->size() != 0 )) {
      result.addCollection( zs2DataCollection, "zsdata_m26" );
    }
     
    return true;
  }
  
#endif


} //namespace eudaq
