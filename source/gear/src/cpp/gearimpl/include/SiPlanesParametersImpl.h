// -*- C++ -*-
#ifndef GEAR_SiPlanesParametersImpl_H
#define GEAR_SiPlanesParametersImpl_H 1

#include "gear/SiPlanesParameters.h"
#include "gearimpl/GearParametersImpl.h"
#include "gearimpl/SiPlanesLayerLayoutImpl.h"

#include <vector>
#include <tuple>


namespace gear {

class SiPlanesLayerLayout;

  /** Abstract description of layers in pixel beam telescope.
   * 
   *  @author T Klimkovich, DESY
   *  @author B. Schwenker, Uni Göttingen
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
   * @param sensitiveID          ID of sensitive volume of telescope plane
   * @param sensitivePositionX   x position of sensitive volume of telescope plane (mm)
   * @param sensitivePositionY   y position of sensitive volume of telescope plane (mm)
   * @param sensitivePositionZ   z position of sensitive volume of telescope plane (mm)
   * @param sensitiveThickness   the thickness of sensitive volume of telescope plane (mm)
   * @param sensitiveRadLenght   the radiation lenght of sensitive volume of telescope plane (mm)
   * @param sensitiveAtomicNum   the atomic number Z of sensitive volume of telescope plane 
   * @param sensitiveAtomicMass  the atomic mass A of sensitive volume of telescope plane  
   * @param sensitiveEulerAlpha  Euler alpha angle of sensitive area of telescope plane (deg)
   * @param sensitiveEulerBeta   Euler beta angle of sensitive area of telescope plane (deg)
   * @param sensitiveEulerGamma  Euler gamma angle of sensitive area of telescope plane (deg)
   * @param sensitiveRotation1   = cos(theta): element (11) of the rotation matrix of sensitive area of telescope plane
   * @param sensitiveRotation2   = -sin(theta): element (12) of the rotation matrix of sensitive area of telescope plane
   * @param sensitiveRotation3   = sin(theta): element (21) of the rotation matrix of sensitive area of telescope plane
   * @param sensitiveRotation4   = cos(theta): element (22) of the rotation matrix of sensitive area of telescope plane
   * @param sensitiveUCells      numbering and pitches for pixel cells along sensitive u axis   
   * @param sensitiveVCells      numbering and pitches for pixel cells along sensitive v axis   
   * @param layerSizeX          size in x direction of nonsensitive volume of telescope plane (mm)
   * @param layerSizeY          size in y direction of nonsensitive volume of telescope plane (mm)
   * @param layerThickness      the thickness of nonsensitive volume of telescope plane (mm)
   * @param layerRadLenght      the radiation lenght of nonsensitive volume of telescope plane (mm)
   * @param layerAtomicNum      the atomic number Z of ladder volume  of telescope plane 
   * @param layerAtomicMass     the atomic mass A of ladder volume of telescope plane  
   */
  virtual void addLayer( 
            int sensitiveID, 
			double sensitivePositionX, double sensitivePositionY, double sensitivePositionZ,
			double sensitiveThickness, double sensitiveRadLength,
            double sensitiveAtomicNum, double sensitiveAtomicMass,
            double sensitiveEulerAlpha,
            double sensitiveEulerBeta,
            double sensitiveEulerGamma,
			double sensitiveRotation1,
			double sensitiveRotation2,
			double sensitiveRotation3,
			double sensitiveRotation4,
            std::vector< std::tuple<int,int,double> > sensitiveUCells,
            std::vector< std::tuple<int,int,double> > sensitiveVCells, 
			double layerSizeX, double layerSizeY, 
            double layerThickness, double layerRadLength,
			double layerAtomicNum, double layerAtomicMass
			)
  {
    _layer.addLayer(
            sensitiveID, 
			sensitivePositionX, sensitivePositionY, sensitivePositionZ,
			sensitiveThickness, sensitiveRadLength,
            sensitiveAtomicNum, sensitiveAtomicMass,
            sensitiveEulerAlpha,
            sensitiveEulerBeta,
            sensitiveEulerGamma,
			sensitiveRotation1,
			sensitiveRotation2,
			sensitiveRotation3,
			sensitiveRotation4,
            sensitiveUCells,
            sensitiveVCells, 
			layerSizeX, layerSizeY, 
            layerThickness, layerRadLength,
			layerAtomicNum, layerAtomicMass ) ;
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
