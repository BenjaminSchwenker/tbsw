#include "eudaq/DataConverterPlugin.hh"
#include "eudaq/StandardEvent.hh"
#include "eudaq/Utils.hh"
#include <iostream>
#include <fstream>
#include <string>
#include <stdint.h>

using namespace std;

struct mPixel  {
    uint8_t x; 
    uint8_t y; 
    uint8_t val;
    bool pivot; 
};


struct mHeader  {
    unsigned char hw1; 
    unsigned char hw2; 
    unsigned char hw3;
    unsigned char hw4; 
    unsigned char hw5; 
};

// All LCIO-specific parts are put in conditional compilation blocks
// so that the other parts may still be used if LCIO is not available.

#if USE_LCIO
#  include "IMPL/LCEventImpl.h"
#  include "IMPL/TrackerRawDataImpl.h"
#  include "IMPL/TrackerDataImpl.h"
#  include "IMPL/LCCollectionVec.h"
#  include "IMPL/LCGenericObjectImpl.h"
#  include "UTIL/CellIDEncoder.h"
#  include "lcio.h"
#endif


namespace eudaq {

  static const char* EVENT_TYPE = "GENERIC";

  class GenericConverterPlugin : public DataConverterPlugin {

    public:

      virtual bool GetStandardSubEvent(StandardEvent & sev, const Event & ev) const {	    
        if (ev.IsBORE()) {
            return true;
        } else if (ev.IsEORE()) {
            return true;
        }
        /*
        const RawDataEvent & ev_raw = dynamic_cast<const RawDataEvent &>(ev);
        if (ev_raw.NumBlocks() != 1)
            return false;
        
        const std::vector<unsigned char> & data = ev_raw.GetBlock(0);
    
        pb::StandardEvent event;
        event.ParseFromArray( data.data(), data.size());
        
        //cout << event.DebugString() << endl;

        for (int i = 0; i < event.plane_size(); i++) {
            const pb::StandardPlane& pb_plane = event.plane(i);

            StandardPlane plane(pb_plane.id(), pb_plane.type(), pb_plane.sensor());
            plane.SetTLUEvent(pb_plane.tluevent());
            
            int flags = 0;
            const pb::Flags& pb_flags = pb_plane.flags();

            if(pb_flags.zs())
                flags |= StandardPlane::FLAG_ZS;
            if(pb_flags.needcds())
                flags |= StandardPlane::FLAG_NEEDCDS;
            if(pb_flags.negative())
                flags |= StandardPlane::FLAG_NEGATIVE;
            if(pb_flags.accumulate())
                flags |= StandardPlane::FLAG_ACCUMULATE;
            if(pb_flags.withpivot())
                flags |= StandardPlane::FLAG_WITHPIVOT;
            if(pb_flags.withsubmat())
                flags |= StandardPlane::FLAG_WITHSUBMAT;
            if(pb_flags.diffcoords())
                flags |= StandardPlane::FLAG_DIFFCOORDS;

            //TODO: Add somethign smarted then this for raw data (need protocol change?)

            if (pb_flags.zs())
                plane.SetSizeZS(pb_plane.xsize(), pb_plane.ysize(), 0, pb_plane.frame_size(), flags);
            else
                plane.SetSizeRaw(pb_plane.xsize(), pb_plane.ysize(), pb_plane.frame_size(), flags);
            
            for (int iframe = 0; iframe < pb_plane.frame_size(); iframe++) {
                const pb::Frame& pb_frame = pb_plane.frame(iframe);
                
                for (int i = 0; i < pb_frame.pixel_size(); i++) {
                    const pb::Pixel& pb_pix = pb_frame.pixel(i);
                    plane.PushPixel(pb_pix.x(), pb_pix.y(), pb_pix.val(), pb_pix.pivot(), iframe);
                }
            }
            
            sev.AddPlane(plane);        
        }
        */
        return true;


      }

#if USE_LCIO

    virtual void GetLCIORunHeader(lcio::LCRunHeader & /*header*/, eudaq::Event const & /*bore*/, eudaq::Configuration const & /*conf*/) const;
    virtual bool GetLCIOSubEvent(lcio::LCEvent & result, const Event & source) const;

    
#endif



    private:

      GenericConverterPlugin()
        : DataConverterPlugin(EVENT_TYPE)
      {}

      static GenericConverterPlugin m_instance;
  };

  GenericConverterPlugin GenericConverterPlugin::m_instance;

#if USE_LCIO

void GenericConverterPlugin::GetLCIORunHeader(lcio::LCRunHeader & header, eudaq::Event const & /*bore*/, eudaq::Configuration const & /*conf*/) const { }


bool GenericConverterPlugin::GetLCIOSubEvent(lcio::LCEvent & result, const Event & source) const {
  
  if (source.IsBORE()) {
    // shouldn't happen 
    return true;
  } else if (source.IsEORE()) {
    // nothing to do
    return true;
  }

  // If we get here it must be a data event
  
  //-----------------------------------------------
  // Decode event data to a DEPFETEvent format 
  
  
   
  
  const RawDataEvent & source_raw = dynamic_cast<const RawDataEvent &>(source);
  if (source_raw.NumBlocks() != 1)
    return false;
        
  const std::vector<unsigned char> & data = source_raw.GetBlock(0);
   
 
  std::cout << "new generic data block has size " << data.size() << std::endl; 
  //std::cout << "uchar has size " << sizeof(unsigned char) << std::endl; 
  //std::cout << "pixel has size " << sizeof(struct mPixel) << std::endl; 
  //std::cout << "max pixels " << sizeof(unsigned char)*data.size()/sizeof(struct mPixel) << std::endl; 
  /*
  long int rc = 0; 
  long int npixels=0; 

  // Skipping the header 

  struct mHeader GROUP_Header = *(struct mHeader *) &data[rc];
  rc += sizeof(struct mHeader); 

  // Loop over pixel data 

  npixels = sizeof(unsigned char)*data.size()/sizeof(struct mPixel);
  
  for (int i = 0; i < npixels; i++) {
    
    struct mPixel gdata = *(struct mPixel *) &data[rc];
    rc += sizeof(struct mPixel);     

    std::cout << "  generic col "  << gdata.x << std::endl; 
    std::cout << "  generic row "  << gdata.y << std::endl;
    std::cout << "  generic val "  << gdata.val << std::endl; 
    std::cout << "  generic pivot "  << gdata.pivot << std::endl;  

    
    
  }
  */      
  
   
  
  
  //-----------------------------------------------
  // Decode event data to a LCIO format 
    
  LCCollectionVec * ZSDataCollection = new LCCollectionVec(LCIO::TRACKERDATA);
  CellIDEncoder<TrackerDataImpl> ZSDataEncoder( "sensorID:6,sparsePixelType:5" , ZSDataCollection );   
  
  /*       
  for (size_t iframe=0;iframe<m_event.size();iframe++) {
    
    // Get read fully decoded data  
    const depfet::DEPFETADCValues& data = m_event[iframe];
     
    // Prepare a TrackerData to store zs pixels
    TrackerDataImpl* zspixels = new TrackerDataImpl;
      
    // Set description for zspixels 
    ZSDataEncoder["sensorID"] = data.getModuleNr();
    ZSDataEncoder["sparsePixelType"] = 0;
    ZSDataEncoder.setCellID( zspixels );
       
    int nPixel = (int) data.at(0).size(); 
           
    for ( int iPixel = 0; iPixel < nPixel; iPixel++ ) {
      int val = data.at(2).at(iPixel);
      int col = data.at(1).at(iPixel);   
      int row = data.at(0).at(iPixel);  
      
      cout << " col  " << col << endl;    
      cout << " row  " << row << endl;   
      cout << " val  " << val << endl;   
      cout << " cm  " << cm << endl;   
       
      // Store pixel data int EUTelescope format 
      zspixels->chargeValues().push_back( col );
      zspixels->chargeValues().push_back( row );
      zspixels->chargeValues().push_back( val );   
    }
      
    // Add event to LCIO collection 
    ZSDataCollection->push_back (zspixels);
  }
  */
  
  result.addCollection( ZSDataCollection, "zsdata_generic" );    
    
  return true;
}
  
#endif



} // namespace eudaq
