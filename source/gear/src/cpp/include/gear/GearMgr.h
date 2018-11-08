
#ifndef GEAR_GEARMGR_H
#define GEAR_GEARMGR_H 1


#include <string>
#include <vector>

#include "GEAR.h"
#include "gear/GearParameters.h"

namespace gear {

class BField;
class SiPlanesParameters;


/** Abstract interface for a manager class that returns the Gear classes for the 
 *  relevant subdetectors.
 * 
 *
 * @author F. Gaede, DESY
 * @version $Id: GearMgr.aid,v 1.9 2008-10-22 15:10:46 engels Exp $
 */
class GearMgr {

public: 
    /// Destructor.
    virtual ~GearMgr() { /* nop */; }

   /** The unique detector name - typically the model name used in the simulation program
    */
    virtual const std::string & getDetectorName() const = 0;

    /** Get named parameters for key. This can be used to describe a subdetector that is not 
     *  yet forseen in the Gear API.
     * 
     *  @throws UnknownParameterException
     */
    virtual const GearParameters & getGearParameters(const std::string & key) const throw (UnknownParameterException, std::exception )  = 0;

    /** Get the B field map.
     *
     *  @throws UnknownParameterException
     */
    virtual const BField & getBField() const throw (UnknownParameterException, std::exception )  = 0;

    /** Get the SiPlanes parameters.
     *
     *  @throws UnknownParameterException
     */
    virtual const SiPlanesParameters & getSiPlanesParameters() const throw (UnknownParameterException, std::exception )  = 0;


    /** Keys of all GearParameters. 
     */ 
    virtual const std::vector<std::string>  & getGearParameterKeys() const = 0;

    /** Set detector name.
     */
    virtual void setDetectorName(const std::string & name) = 0;

    /** Set named parameters for key. This can be used to describe a subdetector that is not 
     *  yet forseen in the Gear API.
     */
    virtual void setGearParameters(const std::string & key, GearParameters * gearParameters) = 0;

    /** Set the BField.
     */
    virtual void setBField(BField * bField) = 0;

    /** Set the SiPlanesParameters.
     */
    virtual void setSiPlanesParameters(SiPlanesParameters * siplanesParameters) = 0;
    
}; // class
} // namespace gear
#endif /* ifndef GEAR_GEARMGR_H */
