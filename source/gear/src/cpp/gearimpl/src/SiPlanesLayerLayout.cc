#include "gearimpl/SiPlanesLayerLayoutImpl.h"
#include <math.h>

namespace gear{

  void SiPlanesLayerLayoutImpl::addLayer(
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
			double layerAtomicNum, double layerAtomicMass)
  {
    SensLayer sL ;
    sL.ID         = sensitiveID ;
    sL.PositionX  = sensitivePositionX ;
    sL.PositionY  = sensitivePositionY ;
    sL.PositionZ  = sensitivePositionZ ;
    sL.Thickness  = sensitiveThickness ;
    sL.RadLength  = sensitiveRadLength ;
    sL.AtomicNumber = sensitiveAtomicNum; 
    sL.AtomicMass = sensitiveAtomicMass;
    sL.EulerAlpha = sensitiveEulerAlpha;
    sL.EulerBeta  = sensitiveEulerBeta;
    sL.EulerGamma = sensitiveEulerGamma;
    sL.Rotation1  = sensitiveRotation1;
    sL.Rotation2  = sensitiveRotation2;
    sL.Rotation3  = sensitiveRotation3;
    sL.Rotation4  = sensitiveRotation4;
    sL.uCells = sensitiveUCells;
    sL.vCells = sensitiveVCells;

    Layer lL ;
    lL.SizeU      = layerSizeU ;
    lL.SizeV      = layerSizeV ;
    lL.Thickness  = layerThickness ;
    lL.RadLength  = layerRadLength ;
    lL.AtomicNumber = layerAtomicNum; 
    lL.AtomicMass = layerAtomicMass;
    
    _lVec.push_back( lL ) ;
    _sensVec.push_back( sL ) ;

  }

  
} //namespace
