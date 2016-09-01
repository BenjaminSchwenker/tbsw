// -*- C++ -*-
#ifndef GEAR_SiPlanesParametersImpl_H
#define GEAR_SiPlanesParametersImpl_H 1

#include "gear/SiPlanesParameters.h"
#include "gearimpl/GearParametersImpl.h"
#include "gearimpl/SiPlanesLayerLayoutImpl.h"


namespace gear {

class SiPlanesLayerLayout;

  /** Abstract description of layers in pixel beam telescope with or without DUT.
   * 
   *  @author T Klimkovich, DESY
   *  @version $Id: 
   */
class SiPlanesParametersImpl : public GearParametersImpl, public SiPlanesParameters {

public: 
  /** C'tor  
   *  @param siplanesID             ID of SiPlanes detector setup   
   *  @param siplanesType           the type of the SiPlanes detector: TelescopeWithDUT, TelescopeWithoutDUT
   *  @param siplanesNumber         the number of Si planes
   */
  SiPlanesParametersImpl(int siplanesID, int siplanesType, int siplanesNumber) ;

  // Destructor.
  virtual ~SiPlanesParametersImpl() { /* nop */; }
  
  // static constants enums( for the type of layout )
  
  static const int TelescopeWithDUT    = 1 ;
  static const int TelescopeWithoutDUT = 2 ;

  /** Adding a Layer to the SiPlanes detector
   *
   * @param layerID            ID of nonsensitive volume of telescope plane
   * @param layerPositionX     x position of nonsensitive volume of telescope plane (mm)
   * @param layerPositionY     y position of nonsensitive volume of telescope plane (mm)
   * @param layerPositionZ     z position of nonsensitive volume of telescope plane (mm)
   * @param layerSizeX         size in x direction of nonsensitive volume of telescope plane (mm)
   * @param layerSizeY         size in y direction of nonsensitive volume of telescope plane (mm)
   * @param layerThickness     the thickness of nonsensitive volume of telescope plane (mm)
   * @param layerRadLenght     the radiation lenght of nonsensitive volume of telescope plane (mm)
   * @param sensitiveID        ID of sensitive volume of telescope plane
   * @param sensitivePixType   PixType of sensitive volume of telescope plane
   * @param sensitivePositionX x position of sensitive volume of telescope plane (mm)
   * @param sensitivePositionY y position of sensitive volume of telescope plane (mm)
   * @param sensitivePositionZ z position of sensitive volume of telescope plane (mm)
   * @param sensitiveSizeX     size in x direction of sensitive volume of telescope plane (mm)
   * @param sensitiveSizeY     size in y direction of sensitive volume of telescope plane (mm)
   * @param sensitiveThickness the thickness of sensitive volume of telescope plane (mm)
   * @param sensitiveNpixelX   number of pixels in x direction of the sensitive area of telescope plane
   * @param sensitiveNpixelY   number of pixels in y direction of the sensitive area of telescope plane
   * @param sensitivePitchX    x size of pitch of sensitive area of telescope plane (mm)
   * @param sensitivePitchY    y size of pitch of sensitive area of telescope plane (mm)
   * @param sensitiveResolutionX intrinsic resolution in X of sensitive area of telescope plane (mm) 
   * @param sensitiveResolutionY intrinsic resolution in Y of sensitive area of telescope plane (mm)
   * @param sensitiveEulerAlpha  Euler alpha angle of sensitive area of telescope plane (deg)
   * @param sensitiveEulerBeta   Euler beta angle of sensitive area of telescope plane (deg)
   * @param sensitiveEulerGamma  Euler gamma angle of sensitive area of telescope plane (deg)
   * @param sensitiveRotation1 = cos(theta): element (11) of the rotation matrix of sensitive area of telescope plane
   * @param sensitiveRotation2 = -sin(theta): element (12) of the rotation matrix of sensitive area of telescope plane
   * @param sensitiveRotation3 = sin(theta): element (21) of the rotation matrix of sensitive area of telescope plane
   * @param sensitiveRotation4 = cos(theta): element (22) of the rotation matrix of sensitive area of telescope plane
   * @param sensitiveRadLenght the radiation lenght of sensitive area of telescope plane (mm)
   */
  virtual void addLayer(int layerID,
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
			double sensitiveRadLength)
  {
    _layer.addLayer( layerID,
                     layerPositionX, layerPositionY, layerPositionZ,
                     layerSizeX, layerSizeY, layerThickness,
		     layerRadLength,
		     sensitiveID, sensitivePixType,
                     sensitivePositionX, sensitivePositionY, sensitivePositionZ,
		     sensitiveSizeX, sensitiveSizeY, sensitiveThickness,
		     sensitiveNpixelX, sensitiveNpixelY,
		     sensitivePitchX, sensitivePitchY,
		     sensitiveResolutionX, sensitiveResolutionY,
                     sensitiveEulerAlpha,
                     sensitiveEulerBeta,
                     sensitiveEulerGamma,
		     sensitiveRotation1,
		     sensitiveRotation2,
		     sensitiveRotation3,
		     sensitiveRotation4,
		     sensitiveRadLength ) ;
    return ;
  }
  
  /** Adding DUT to the SiPlanes detector
   *
   * @param dutID                  ID of nonsensitive volume of the DUT
   * @param dutPositionX           x position of nonsensitive volume of the DUT plane (mm)
   * @param dutPositionY           y position of nonsensitive volume of the DUT plane (mm)
   * @param dutPositionZ           z position of nonsensitive volume of the DUT plane (mm)    
   * @param dutSizeX               size in x direction of nonsensitive volume of the DUT (mm)
   * @param dutSizeY               size in y direction of nonsensitive volume of the DUT (mm)
   * @param dutThickness           the thickness of nonsensitive volume of the DUT layer (mm)    
   * @param dutRadLenght           the radiation lenght of nonsensitive volume of the DUT layer (mm)           
   * @param dutsensitiveID         ID of nonsensitive volume of the DUT
   * @param dutsensitivePixType    PixType of nonsensitive volume of the DUT
   * @param dutsensitivePositionX  x position of sensitive volume of the DUT plane (mm)
   * @param dutsensitivePositionY  y position of sensitive volume of the DUT plane (mm)
   * @param dutsensitivePositionZ  z position of sensitive volume of the DUT plane (mm)    
   * @param dutsensitiveSizeX      size in x direction of sensitive volume of the DUT (mm)
   * @param dutsensitiveSizeY      size in y direction of sensitive volume of the DUT (mm)
   * @param dutsensitiveThickness  the thickness of sensitive volume of the DUT layer (mm)
   * @param dutsensitiveNpixelX    number of pixels in x direction of sensitive area of the DUT
   * @param dutsensitiveNpixelY    number of pixels in y direction of sensitive area of the DUT
   * @param dutsensitivePitchX     x size of pitch of sensitive area of the DUT (mm)
   * @param dutsensitivePitchY     y size of pitch of sensitive area of the DUT(mm)
   * @param dutsensitiveResolutionX intrinsic resolution in X of sensitive area of the DUT(mm), not relevant, just for completeness
   * @param dutsensitiveResolutionY intrinsic resolution in Y of sensitive area of the DUT(mm), not relevant, just for completeness
   * @param dutsensitiveEulerAlpha  Euler alpha angle of sensitive area of the DUT (deg)
   * @param dutsensitiveEulerBeta   Euler beta angle of sensitive area of the DUT (deg)
   * @param dutsensitiveEulerGamma  Euler gamma angle of sensitive area of the DUT (deg)
   * @param dutsensitiveRotation1  = cos(theta): element (11) of the rotation matrix of sensitive area of the DUT
   * @param dutsensitiveRotation2  = -sin(theta): element (12) of the rotation matrix of sensitive area of the DUT
   * @param dutsensitiveRotation3  = sin(theta): element (21) of the rotation matrix of sensitive area of the DUT
   * @param dutsensitiveRotation4  = cos(theta): element (22) of the rotation matrix of sensitive area of the DUT
   * @param dutsensitiveRadLenght  radiation lenght of sensitive area of the DUT layer (mm)                         
   */


  virtual void addDUT(int dutID,
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
    _layer.addDUT( dutID,
		   dutPositionX, dutPositionY, dutPositionZ,
		   dutSizeX, dutSizeY, dutThickness,
		   dutRadLength,
		   dutsensitiveID, dutsensitivePixType,
		   dutsensitivePositionX, dutsensitivePositionY, dutsensitivePositionZ,
		   dutsensitiveSizeX, dutsensitiveSizeY, dutsensitiveThickness,
		   dutsensitiveNpixelX, dutsensitiveNpixelY,
		   dutsensitivePitchX, dutsensitivePitchY,
		   dutsensitiveResolutionX, dutsensitiveResolutionY, 
                   dutsensitiveEulerAlpha, 
                   dutsensitiveEulerBeta, 
                   dutsensitiveEulerGamma,
		   dutsensitiveRotation1,
		   dutsensitiveRotation2,
		   dutsensitiveRotation3,
		   dutsensitiveRotation4,
		   dutsensitiveRadLength ) ;
    return ;
  }
  
      
  /** Returns the layer layout of SiPlanes detector 
   */
  virtual const SiPlanesLayerLayout & getSiPlanesLayerLayout() const { return _layer ; }
  

  /** Returns the ID of SiPlanes detector setup
   */

  virtual int getSiPlanesID() const { return _siplanesID ; }

  /** Returns the type of SiPlanes detector: SiPlanesParametersImpl.TelescopeWithDUT 
   *  or SiPlanesParametersImpl.TelescopeWithoutDUT
   */
  virtual int getSiPlanesType() const { return _siplanesType ; }

  /** Returns the number of Si planes
   */
  virtual int getSiPlanesNumber() const { return _siplanesNumber ; }
  
protected:
  
  SiPlanesLayerLayoutImpl _layer ;
  
  int _siplanesID;

  int _siplanesType ;
  
  int _siplanesNumber ;

private:

}; // class

} // namespace gear

#endif /* ifndef GEAR_SIPLANESPARAMETERS_H */
