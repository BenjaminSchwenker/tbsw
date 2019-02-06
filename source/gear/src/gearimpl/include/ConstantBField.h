#ifndef GEAR_ConstantBField_h
#define GEAR_ConstantBField_h 

#include "gear/BField.h"
#include "gearimpl/GearParametersImpl.h"

namespace gear {

  /** Global B field map implementation with constant field.
   *
   *  @author F.Gaede, DESY
   *  @version $Id: 
   */
	
class ConstantBField : public GearParametersImpl , public BField {
protected: 	
    Vector3D _b ;
public:		

    ConstantBField(Vector3D b) : _b(b) {} 

	/// Destructor.
    virtual ~ConstantBField() { /* nop */; }

      //extends GearParameters {
    /** Returns the B field vector in Tesla at  given point.
     */	
    virtual Vector3D at(Vector3D point) const { return _b ; }
}; // class
} // namespace gear
#endif /* ifndef GEAR_BFIELD_H */
