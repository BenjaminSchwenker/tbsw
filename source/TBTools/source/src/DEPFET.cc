#include "DEPFET.h"

using namespace depfet;

const int      DEPFET::GOODPIXEL           =  0;
const int      DEPFET::BADPIXEL            =  1;
const int      DEPFET::HITPIXEL            = -1;
const int      DEPFET::MISSINGPIXEL        =  2;
const int      DEPFET::FIRINGPIXEL         =  3;
const char *   DEPFET::MATRIXDEFAULTENCODING    = "sensorID:16,xMin:12,xMax:12,yMin:12,yMax:12";
const char *   DEPFET::ZSDATADEFAULTENCODING    = "sensorID:16,sparsePixelType:5";
const char *   DEPFET::CLUSTERDEFAULTENCODING   = "sensorID:16,clusterID:8,xSeed:12,ySeed:12,xCluSize:5,yCluSize:5,quality:5";
const char *   DEPFET::PULSEDEFAULTENCODING     = "sensorID:16,clusterID:8,xSeed:12,ySeed:12,xCluSize:5,yCluSize:5,type:5";
const char *   DEPFET::ZSCLUSTERDEFAULTENCODING = "sensorID:16,clusterID:8,sparsePixelType:5,quality:5";


namespace depfet {

  ClusterQuality operator&(ClusterQuality a, ClusterQuality b) {
    return ClusterQuality( static_cast<int>(a) & static_cast<int>(b) );
  }

  ClusterQuality operator|(ClusterQuality a, ClusterQuality b) {
    return ClusterQuality( static_cast<int>(a) | static_cast<int>(b) );
  }

  std::ostream& operator<<(std::ostream& os, const ClusterQuality& quality) {
    if ( quality == kGoodCluster ) {
      os << "kGoodCluster (" << static_cast<int > (quality) << ")";
      return os;
    }
    bool moreThanOne = false;
    if ( quality & kIncompleteCluster ) {
      if ( moreThanOne ) os << ", ";
      os << "kIncompleteCluster";
      moreThanOne = true;
    }
    if ( quality & kBorderCluster ) {
      if ( moreThanOne ) os << ", ";
      os << "kBorderCluster";
      moreThanOne = true;
    }
    if ( quality & kMergedCluster ) {
      if ( moreThanOne ) os << ", ";
      os << "kMergedCluster";
      moreThanOne = true;
    }
    os << " (" << static_cast<int> (quality) << ")";
    return os;
  }


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

} // end namespace depfet
