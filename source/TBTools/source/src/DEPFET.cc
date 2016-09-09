#include "DEPFET.h"

using namespace depfet;


const char *   DEPFET::MATRIXDEFAULTENCODING    = "sensorID:17,uMax:12,vMax:12"; 
const char *   DEPFET::ZSDATADEFAULTENCODING    = "sensorID:17,sparsePixelType:5";
const char *   DEPFET::ZSCLUSTERDEFAULTENCODING = "sensorID:17,clusterID:5,sparsePixelType:5,quality:2";


namespace depfet {

  /*
  
  PixelQuality operator&(PixelQuality a, PixelQuality b) {
    return PixelQuality( static_cast<int>(a) & static_cast<int>(b) );
  }

  PixelQuality operator|(PixelQuality a, PixelQuality b) {
    return PixelQuality( static_cast<int>(a) | static_cast<int>(b) );
  }
  
  std::ostream& operator<<(std::ostream& os, const PixelQuality& quality) {
    if ( quality == kGoodPixel ) {
      os << "kGoodPixel (" << static_cast<int > (quality) << ")";
      return os;
    }
    bool moreThanOne = false;
    if ( quality & kDeadPixel ) {
      if ( moreThanOne ) os << ", ";
      os << "kDeadPixel";
      moreThanOne = true;
    }
    if ( quality & kNoisyPixel ) {
      if ( moreThanOne ) os << ", ";
      os << "kNoisyPixel";
      moreThanOne = true;
    }
    if ( quality & kHotPixel ) {
      if ( moreThanOne ) os << ", ";
      os << "kHotPixel";
      moreThanOne = true;
    }
    if ( quality & kBrigthPixel ) {
      if ( moreThanOne ) os << ", ";
      os << "kBrigthPixel";
      moreThanOne = true;
    }
    if ( quality & kBadChannel ) {
      if ( moreThanOne ) os << ", ";
      os << "kBadChannel";
      moreThanOne = true;
    }
    if ( quality & kBadPixel ) {
      if ( moreThanOne ) os << ", ";
      os << "kBadPixel";
      moreThanOne = true;
    }
    os << " (" << static_cast<int> (quality) << ")";
    return os;
  }

  */

} // end namespace depfet
