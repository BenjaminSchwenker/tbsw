// -*- C++ -*-
#ifndef GEAR_SiPlanesParametersImpl_H
#define GEAR_SiPlanesParametersImpl_H 1

#include "gear/SiPlanesParameters.h"
#include "gearimpl/GearParametersImpl.h"
#include "gearimpl/SiPlanesLayerLayoutImpl.h"


namespace gear {

class SiPlanesLayerLayout;

  /** Abstract description of layers in pixel beam telescope.
   * 
   *  @author T Klimkovich, DESY
   *  @version $Id: 
   */
class SiPlanesParametersImpl : public GearParametersImpl, public SiPlanesParameters {

public: 
  /** C'tor  
   *  @param siplanesID             ID of SiPlanes detector setup   
   *  @param siplanesNumber         the number of Si planes
   */
  SiPlanesParametersImpl(int siplanesID, int siplanesNumber) ;

  // Destructor.
  virtual ~SiPlanesParametersImpl() { /* nop */; }
  
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
   
  /** Returns the layer layout of SiPlanes detector 
   */
  virtual const SiPlanesLayerLayout & getSiPlanesLayerLayout() const { return _layer ; }
  
  /** Returns the ID of SiPlanes detector setup
   */
  virtual int getSiPlanesID() const { return _siplanesID ; }
  
  /** Returns the number of Si planes
   */
  virtual int getSiPlanesNumber() const { return _siplanesNumber ; }
  
protected:
  
  SiPlanesLayerLayoutImpl _layer ;
  
  int _siplanesID;
  
  int _siplanesNumber ;

private:

}; // class

} // namespace gear

#endif /* ifndef GEAR_SIPLANESPARAMETERS_H */
