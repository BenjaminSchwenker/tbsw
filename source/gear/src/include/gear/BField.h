
#ifndef GEAR_BFIELD_H
#define GEAR_BFIELD_H 1

#include "GEAR.h"
#include "gear/GearParameters.h"

namespace gear {

  /** Global B field map.
   *
   *  @author F.Gaede, DESY
   *  @version $Id: 
   */
	
class BField : virtual public GearParameters {

public: 
    /// Destructor.
    virtual ~BField() { /* nop */; }

    /** Returns the B field vector in Tesla at  given point.
     */	
    virtual Vector3D at(const Vector3D &point) const = 0;
}; // class
} // namespace gear
#endif /* ifndef GEAR_BFIELD_H */
