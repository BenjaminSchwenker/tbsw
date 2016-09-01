// -*- C++ -*-
#ifndef GEAR_VXDLayerLayoutImpl_H
#define GEAR_VXDLayerLayoutImpl_H 1

#include "gear/VXDLayerLayout.h"
#include <vector>
#include <math.h>

namespace gear {

/** Abstract description of layers in a Vertex detector. <br>
 *  This assumes a symmetric layout of ladders, arranged in layers equidistant 
 *  to IP.  <br>
 *  The sensitive area is assumed to be inside the ladders but can be positioned independently.
 * 
 * @author R. Lippe, DESY
 * @version $Id: 
 */

class VXDLayerLayoutImpl : public VXDLayerLayout {


  friend class VXDParametersImpl ;

public: 
  
  /** Helper class for layer properties */
  struct Layer {
    int    NLadders ;
    double internalPhi0 ;
    double Distance ;
    double Offset ;
    double Thickness ;
    double Length ;
    double Width ;
    double RadLength ;
  } ;

  typedef std::vector<Layer> LayerVec ;

  // Destructor.
  virtual ~VXDLayerLayoutImpl() { /* nop */; }
  
  /** The total number of layers.
   */
  virtual int getNLayers() const { return _lVec.size() ; }
  
  /** The number of ladders in the layer layerIndex - layer indexing starts at 0
   *  for the layer closest to IP.
   */
  virtual int getNLadders(int layerIndex) const { return _lVec.at( layerIndex ).NLadders  ; }

  /** Azimuthal angle of the (outward pointing) normal of the first ladder.
   *  Phi0==0 corresponds to the first ladder's normal coinciding (if offset==0) with the x-axis.
   *  The layerIndex starts at 0 for the layer closest to IP.
   */
  virtual double getPhi0(int layerIndex) const { return M_PI/2. - _lVec.at( layerIndex ).internalPhi0  ; }
  
  /** The radiation length in the support structure ladders of layer layerIndex - layer indexing starts at 0
   *  for the layer closest to IP.
   */
  virtual double getLadderRadLength(int layerIndex) const { return _lVec.at( layerIndex ).RadLength  ; }
  
  /** The distance of ladders in layer layerIndex from the IP - layer indexing starts at 0
   *  for the layer closest to IP.
   */
  virtual double getLadderDistance(int layerIndex) const { return _lVec.at( layerIndex ).Distance  ; }

  /** The thickness in mm of the ladders in layerIndex - layer indexing starting at 0
   *  for the layer closest to IP.
   */
  virtual double getLadderThickness(int layerIndex) const { return _lVec.at( layerIndex ).Thickness  ; }

  /** The offset of the ladder in mm defines the shift of the ladder in the direction of increasing phi
   *  perpendicular to the ladder's normal. For example if the first ladder is at phi0==0 then the offset 
   *  defines an upward shift of the ladder (parallel to the y-axis).  
   *  Layer indexing starts at 0 for the layer closest to IP.
   *  @see getPhi0
   *  @see getSensitiveOffset
   */
  virtual double getLadderOffset(int layerIndex) const { return  _lVec.at( layerIndex ).Offset  ; }

  /** The width of the ladder in layer in mm for ladders in layer layerIndex -
   *  layer indexing starting at 0 from the layer closest to IP.
   */
  virtual double getLadderWidth(int layerIndex) const { return _lVec.at( layerIndex ).Width  ; }

  /** The length of the ladder in z direction in mm for ladders in layer layerIndex -
   *  layer indexing starting at 0 from the layer closest to IP.
   */
  virtual double getLadderLength(int layerIndex) const { return _lVec.at( layerIndex ).Length  ; }
  
  /** The radiation length in sensitive volumes in layer layerIndex - layer indexing starts at 0
   *  for the layer closest to IP.
   */
  virtual double getSensitiveRadLength(int layerIndex) const { return _sVec.at( layerIndex ).RadLength  ; }

  /** The distance of sensitive area in ladders in layer layerIndex from the IP.
   */
  virtual double getSensitiveDistance(int layerIndex) const { return _sVec.at( layerIndex ).Distance  ; }

  /** The thickness in mm of the sensitive area in ladders in layer layerIndex.
   */
  virtual double getSensitiveThickness(int layerIndex) const { return _sVec.at( layerIndex ).Thickness  ; }

  /** Same as getLadderOffset() except for the sensitive part of the ladder.
   * @see getLadderOffset
   */
  virtual double getSensitiveOffset(int layerIndex) const { return _sVec.at( layerIndex ).Offset  ; }

  /** The width of the sensitive area in ladders in layer layerIndex in mm.
   */
  virtual double getSensitiveWidth(int layerIndex) const { return _sVec.at( layerIndex ).Width  ; }

  /** The length of the sensitive area in ladders in z direction in mm for ladders in layer layerIndex.
   */
  virtual double getSensitiveLength(int layerIndex) const { return _sVec.at( layerIndex ).Length  ; }

  /** Add a new layer at the given positon
   */
  virtual void addLayer(int nLadders, double phi0, 
			// ladder
			double ladderDistance, double ladderOffset, double ladderThickness,
			double ladderLength, double ladderWidth, double ladderRadLength,
			// sensitive
			double sensitiveDistance, double sensitiveOffset, double sensitiveThickness,
			double sensitiveLength, double sensitiveWidth, double sensitiveRadLength);
 
  /** returns maximum radius for a given layer
   */
  virtual double getMaxRadius( int layerIndex , bool sensitive=false ) const ; 
  
  /** returns starting phi for first ladder in layer layerIndex (on side facing IP)
   */
  virtual double getStartInnerPhi( int layerIndex , bool sensitive=false ) const ;

  /** returns ending phi for first ladder in layer layerIndex (on side facing IP)
   */
  virtual double getEndInnerPhi( int layerIndex , bool sensitive=false ) const ;

  /** returns starting phi for first ladder in layer layerIndex (on side away from IP)
   */
  virtual double getStartOuterPhi( int layerIndex , bool sensitive=false ) const ;

  /** returns ending phi for first ladder in layer layerIndex (on side away from IP)
   */
  virtual double getEndOuterPhi( int layerIndex , bool sensitive=false ) const ;
  
  /** returns thickness as viewed under the angle phi only for the first ladder
   *  in layer layerIndex.
   *  returns -1 if phi is not covered by ladder
   */
  virtual double getThicknessForPhi( int layerIndex, double phi , bool sensitive=false ) const ;
  
  
protected:

  double getInternalPhi0(int layerIndex) const { return _lVec.at( layerIndex ).internalPhi0  ; }


  // Ladder
  LayerVec _lVec ;
  // Sensitive
  LayerVec _sVec ;
    
private:


}; // class
} // namespace gear
#endif /* ifndef GEAR_VXDLAYERLAYOUT_H */
