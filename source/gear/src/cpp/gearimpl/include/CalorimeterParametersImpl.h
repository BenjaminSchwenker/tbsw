#ifndef GEAR_CalorimeterParametersImpl_H
#define GEAR_CalorimeterParametersImpl_H 1

#include <vector>

#include "gear/CalorimeterParameters.h"
#include "gearimpl/GearParametersImpl.h"
#include "gearimpl/LayerLayoutImpl.h"

namespace gear {

class LayerLayoutImpl ;

/** Implementation of 
 * 
 * @author F. Gaede, DESY
 * @version $Id: CalorimeterParametersImpl.h,v 1.1.1.1 2005-09-13 14:41:03 gaede Exp $
 */

  class CalorimeterParametersImpl :  public GearParametersImpl,  public CalorimeterParameters {
    
  public: 

      /** C'tor for a cylindrical (octagonal) BARREL calorimeter. 
       */
      CalorimeterParametersImpl(double rMin, double zMax, int symOrder=8, double phi0=0.0)  ;


      /** C'tor for a cylindrical (octagonal) ENDCAP calorimeter. 
       */
      CalorimeterParametersImpl(double rMin, double rMax, double zMin, int symOrder=2, double phi0=0.0) ; 


    /// Destructor.
    virtual ~CalorimeterParametersImpl() { /* nop */; }
    
    
    /** The layer layout in the calorimeter*/
    virtual const LayerLayout & getLayerLayout() const  { return _layout ; }
    
    /** The type of the calorimeter layout: CalorimeterParametersImpl.BARREL or
     *  CalorimeterParametersImpl.ENDCAP.
     */
    virtual int getLayoutType() const { return _type ; }
    
    /** Extent of the calorimeter in the r-z-plane  [ rmin, rmax, zmin, zmax ] in mm.
     */
    virtual const std::vector<double>& getExtent() const ;
    

    /** The order of the rotational symmetry <br>
     *   8 for an octagonal barrel calorimeter<br>
     *   2 for an endcap calorimeter<br>
     *   1 for a standalone prototype<br>
     *   0 for an idealized cylindrical calorimeter.
     */
    virtual int getSymmetryOrder() const { return _sym ; }
    
    /** The azimuthal angle of the first module if symmetryOrder > 0.
     */
    virtual double getPhi0() const { return _phi ; }
    

    virtual LayerLayoutImpl& layerLayout() { return _layout ; }

  protected:
    
    int _type ;
    int _sym ;
    double _phi ;
    LayerLayoutImpl _layout ;
    mutable std::vector<double> _extent ;
    
  }; // class
} // namespace gear
#endif /* ifndef GEAR_CalorimeterParametersImpl_H */
