#ifndef GEAR_TPCParametersImpl_H
#define GEAR_TPCParametersImpl_H 1

#include "gear/TPCParameters.h"

#include "gearimpl/GearParametersImpl.h"

namespace gear {

class PadRowLayout2D;

/**The geometry properties of
 * a TPC like detector needed for reconstruction code. <br>
 * This assumes a row based layout of the pads, either in a cartesian or in a 
 * polar angle coordinate system.  <br>
 * Most of this class is based on ideas by Dean Karlen.
 * 
 * @author F. Gaede, DESY
 * @version $Id: TPCParametersImpl.h,v 1.1.1.1 2005-09-13 14:41:03 gaede Exp $
 */

  class TPCParametersImpl :  public GearParametersImpl,  public TPCParameters {
    
  public: 
    
    TPCParametersImpl() ;

//     /** Copy constructor for base parameters */
//     TPCParametersImpl( GearParameters& gearParameters ) ;

    /// Destructor.
    virtual ~TPCParametersImpl() ;
    
    /** The pad layout in the endcaps */
    virtual const PadRowLayout2D & getPadLayout() const;

    /** The maximum drift length in the TPC in mm.
     */
    virtual double getMaxDriftLength() const;

    /** The electron drift velocity in the TPC in mm/s.
     */
    virtual double getDriftVelocity() const;

    /** The readout frequency in Hz.
     */
    virtual double getReadoutFrequency() const;



  virtual void setPadLayout( PadRowLayout2D * padLayout ) ;

  virtual void  setMaxDriftLength( double maxDriftLength ) ;

  virtual void  setDriftVelocity( double driftVelocity ) ;

  virtual void  setReadoutFrequency( double readoutFrequency ) ;
  



protected:

  PadRowLayout2D * _padLayout ;

  double _maxDriftLength ;
  double _driftVelocity ;
  double _readoutFrequency ;
 

}; // class

} // namespace gear

#endif /* ifndef GEAR_TPCParametersImpl_H */
