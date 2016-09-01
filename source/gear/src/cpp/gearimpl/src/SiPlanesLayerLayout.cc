#include "gearimpl/SiPlanesLayerLayoutImpl.h"
#include <math.h>

namespace gear{

  void SiPlanesLayerLayoutImpl::addLayer(
			      // layer
			      int layerID,	 
			      double layerPositionX, double layerPositionY, double layerPositionZ,
			      double layerSizeX, double layerSizeY, double layerThickness,
			      double layerRadLength,
			      // sensitive
			      int sensitiveID, int sensitivePixType, 
			      double sensitivePositionX, double sensitivePositionY, double sensitivePositionZ,
			      double sensitiveSizeX, double sensitiveSizeY, double sensitiveThickness,
			      int sensitiveNpixelX, int sensitiveNpixelY,
			      double sensitivePitchX,double sensitivePitchY,
			      double sensitiveResolutionX, double sensitiveResolutionY,
                              double sensitiveEulerAlpha,
                              double sensitiveEulerBeta,
                              double sensitiveEulerGamma,
			      double sensitiveRotation1,
			      double sensitiveRotation2,
			      double sensitiveRotation3,
			      double sensitiveRotation4,
			      double sensitiveRadLength )
  {

    Layer lL ;
    SensLayer sL ;
    lL.ID         = layerID ;
    lL.PositionX  = layerPositionX ;
    lL.PositionY  = layerPositionY ;
    lL.PositionZ  = layerPositionZ ;
    lL.SizeX      = layerSizeX ;
    lL.SizeY      = layerSizeY ;
    lL.Thickness  = layerThickness ;
    lL.RadLength  = layerRadLength ;

    sL.ID         = sensitiveID ;
    sL.PixType    = sensitivePixType ;
    sL.PositionX  = sensitivePositionX ;
    sL.PositionY  = sensitivePositionY ;
    sL.PositionZ  = sensitivePositionZ ;
    sL.SizeX      = sensitiveSizeX ;
    sL.SizeY      = sensitiveSizeY ;
    sL.Thickness  = sensitiveThickness ;
    sL.NpixelX    = sensitiveNpixelX;
    sL.NpixelY    = sensitiveNpixelY;
    sL.PitchX     = sensitivePitchX;
    sL.PitchY     = sensitivePitchY;
    sL.ResolutionX = sensitiveResolutionX;
    sL.ResolutionY = sensitiveResolutionY;
    sL.EulerAlpha = sensitiveEulerAlpha;
    sL.EulerBeta  = sensitiveEulerBeta;
    sL.EulerGamma = sensitiveEulerGamma;
    sL.Rotation1  = sensitiveRotation1;
    sL.Rotation2  = sensitiveRotation2;
    sL.Rotation3  = sensitiveRotation3;
    sL.Rotation4  = sensitiveRotation4;
    sL.RadLength  = sensitiveRadLength ;

    _lVec.push_back( lL ) ;
    _sVec.push_back( sL ) ;

  }

  void SiPlanesLayerLayoutImpl::addDUT(
			      // layer
			      int dutID,	 
			      double dutPositionX, double dutPositionY, double dutPositionZ,
			      double dutSizeX, double dutSizeY, double dutThickness,
			      double dutRadLength,
			      // sensitive
			      int dutsensitiveID, int dutsensitivePixType, 
			      double dutsensitivePositionX, double dutsensitivePositionY, double dutsensitivePositionZ,
			      double dutsensitiveSizeX, double dutsensitiveSizeY, double dutsensitiveThickness,
			      int dutsensitiveNpixelX, int dutsensitiveNpixelY,
			      double dutsensitivePitchX,double dutsensitivePitchY,
			      double dutsensitiveResolutionX, double dutsensitiveResolutionY,
                              double dutsensitiveEulerAlpha,
                              double dutsensitiveEulerBeta,
                              double dutsensitiveEulerGamma, 
			      double dutsensitiveRotation1,
			      double dutsensitiveRotation2,
			      double dutsensitiveRotation3,
			      double dutsensitiveRotation4,
			      double dutsensitiveRadLength)
  {

    Layer lDUT ;
    SensLayer sDUT ;
    lDUT.ID         = dutID ;
    lDUT.PositionX  = dutPositionX ;
    lDUT.PositionY  = dutPositionY ;
    lDUT.PositionZ  = dutPositionZ ;
    lDUT.SizeX      = dutSizeX ;
    lDUT.SizeY      = dutSizeY ;
    lDUT.Thickness  = dutThickness ;
    lDUT.RadLength  = dutRadLength ;

    sDUT.ID         = dutsensitiveID ;
    sDUT.PixType    = dutsensitivePixType ;
    sDUT.PositionX  = dutsensitivePositionX ;
    sDUT.PositionY  = dutsensitivePositionY ;
    sDUT.PositionZ  = dutsensitivePositionZ ;
    sDUT.SizeX      = dutsensitiveSizeX ;
    sDUT.SizeY      = dutsensitiveSizeY ;
    sDUT.Thickness  = dutsensitiveThickness ;
    sDUT.NpixelX    = dutsensitiveNpixelX;
    sDUT.NpixelY    = dutsensitiveNpixelY;
    sDUT.PitchX     = dutsensitivePitchX;
    sDUT.PitchY     = dutsensitivePitchY;
    sDUT.ResolutionX = dutsensitiveResolutionX;
    sDUT.ResolutionY = dutsensitiveResolutionY;
    sDUT.EulerAlpha = dutsensitiveEulerAlpha;
    sDUT.EulerBeta  = dutsensitiveEulerBeta;
    sDUT.EulerGamma = dutsensitiveEulerGamma;
    sDUT.Rotation1  = dutsensitiveRotation1;
    sDUT.Rotation2  = dutsensitiveRotation2;
    sDUT.Rotation3  = dutsensitiveRotation3;
    sDUT.Rotation4  = dutsensitiveRotation4;
    sDUT.RadLength  = dutsensitiveRadLength ;

    _lDut = lDUT ;
    _sDut = sDUT  ;

  }

} //namespace
