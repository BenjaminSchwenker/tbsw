// -*- C++ -*-
#ifndef GEAR_SiPlanesLayerLayoutImpl_H
#define GEAR_SiPlanesLayerLayoutImpl_H 1

#include "gear/SiPlanesLayerLayout.h"


namespace gear {

/** Abstract description of layers in pixel beam telescope. 
 * 
 * @author T. Klimkovich, DESY
 * @author B. Schwenker, Uni Göttingen
 */

class SiPlanesLayerLayoutImpl : public SiPlanesLayerLayout {

public: 
  
  /** Helper class for layer properties */
  struct SensLayer {
    int ID ;
    double PositionX ;
    double PositionY ;
    double PositionZ ;
    double Thickness ;
    double RadLength ;
    double AtomicNumber ;
    double AtomicMass ;     
    double EulerAlpha ; 
    double EulerBeta ; 
    double EulerGamma ; 
    double Rotation1 ;
    double Rotation2 ;
    double Rotation3 ;
    double Rotation4 ;
    std::vector< std::tuple<int,int,double> > uCells ; 
    std::vector< std::tuple<int,int,double> > vCells ; 
  } ;
  
  struct Layer {
    double SizeU ;
    double SizeV ;
    double Thickness ;
    double RadLength ;
    double AtomicNumber ;
    double AtomicMass ;
  } ;
  
  typedef std::vector<SensLayer> SensLayerVec ;
  typedef std::vector<Layer> LayerVec ;
  
  // Destructor.
  virtual ~SiPlanesLayerLayoutImpl() { /* nop */; }
  
  virtual int getNLayers() const { return _sensVec.size() ; }
  
  virtual int getSensitiveID(int layerIndex) const { return _sensVec.at( layerIndex ).ID  ; }
  virtual double getSensitivePositionX(int layerIndex) const { return _sensVec.at( layerIndex ).PositionX  ; }
  virtual double getSensitivePositionY(int layerIndex) const { return _sensVec.at( layerIndex ).PositionY  ; }
  virtual double getSensitivePositionZ(int layerIndex) const { return _sensVec.at( layerIndex ).PositionZ  ; }
  virtual double getSensitiveThickness(int layerIndex) const { return _sensVec.at( layerIndex ).Thickness  ; }
  virtual double getSensitiveRadLength(int layerIndex) const { return _sensVec.at( layerIndex ).RadLength  ; }
  virtual double getSensitiveAtomicNumber(int layerIndex) const { return _sensVec.at( layerIndex ).AtomicNumber  ; } 
  virtual double getSensitiveAtomicMass(int layerIndex) const { return _sensVec.at( layerIndex ).AtomicMass  ; }   
  virtual double getSensitiveRotationAlpha(int layerIndex) const { return _sensVec.at( layerIndex ).EulerAlpha  ; }
  virtual double getSensitiveRotationBeta(int layerIndex) const { return _sensVec.at( layerIndex ).EulerBeta  ; } 
  virtual double getSensitiveRotationGamma(int layerIndex) const { return _sensVec.at( layerIndex ).EulerGamma  ; }  
  virtual double getSensitiveRotation1(int layerIndex) const { return _sensVec.at( layerIndex ).Rotation1  ; }
  virtual double getSensitiveRotation2(int layerIndex) const { return _sensVec.at( layerIndex ).Rotation2  ; }
  virtual double getSensitiveRotation3(int layerIndex) const { return _sensVec.at( layerIndex ).Rotation3  ; }
  virtual double getSensitiveRotation4(int layerIndex) const { return _sensVec.at( layerIndex ).Rotation4  ; }
  virtual std::vector< std::tuple<int,int,double> > getSensitiveUCells(int layerIndex) const { return _sensVec.at( layerIndex ).uCells  ; }
  virtual std::vector< std::tuple<int,int,double> > getSensitiveVCells(int layerIndex) const { return _sensVec.at( layerIndex ).vCells  ; }
  
  virtual double getLayerSizeU(int layerIndex) const { return _lVec.at( layerIndex ).SizeU  ; }
  virtual double getLayerSizeV(int layerIndex) const { return _lVec.at( layerIndex ).SizeV  ; }
  virtual double getLayerThickness(int layerIndex) const { return _lVec.at( layerIndex ).Thickness  ; }
  virtual double getLayerRadLength(int layerIndex) const { return _lVec.at( layerIndex ).RadLength  ; }
  virtual double getLayerAtomicNumber(int layerIndex) const { return _lVec.at( layerIndex ).AtomicNumber  ; } 
  virtual double getLayerAtomicMass(int layerIndex) const { return _lVec.at( layerIndex ).AtomicMass  ; }   

  /** Add a new layer at the given position
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
			double layerSizeU, double layerSizeV, 
            double layerThickness, double layerRadLength,
			double layerAtomicNum, double layerAtomicMass
       );
 
protected:
  
  // Support layer 
  LayerVec _lVec ;
  // Sensitive layer
  SensLayerVec _sensVec ;
 
}; // class
} // namespace gear
#endif /* ifndef GEAR_SIPLANESLAYERLAYOUT_H */
