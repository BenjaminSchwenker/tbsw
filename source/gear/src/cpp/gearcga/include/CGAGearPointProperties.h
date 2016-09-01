#ifndef GEAR_CGAGEARPOINTPROPERTIES_H
#define GEAR_CGAGEARPOINTPROPERTIES_H 1

#include <string>
#include <vector>

#include "gear/GEAR.h"
#include "gear/GearPointProperties.h"

namespace gear {

/** CGA implementation of the abstract interface that returns the (material) 
 *  properties of a given point in world coordinates.  
 * @author G. Musat, Ecole Polytechnique
 * @version $Id: CGAGearPointProperties.h,v 1.2 2007-03-07 15:37:45 gaede Exp $
 */
class CGAGearPointProperties: public GearPointProperties {

public: 
    CGAGearPointProperties(std::string steer, std::string model,
                        std::string setup, std::string host, std::string user,
                        std::string password);

    /// Destructor.
    virtual ~CGAGearPointProperties() { /* nop */; }

    //vec getPosition returns the current position (if ID is known)
    /** The cellID of the the sensitive detector at pos.
     */
    virtual long64 getCellID(const Vector3D & pos) const throw (NotImplementedException, std::exception );

    // 	bool getMaterial  (true if defined material, false if problem)
    /** Name of material at pos.  
     */
    virtual const std::string & getMaterialName(const Vector3D & pos) const throw (NotImplementedException, std::exception );

    /** Density in kg/m^3 at pos.  
     */
    virtual double getDensity(const Vector3D & pos) const throw (NotImplementedException, std::exception );

//     /** Name of material at pos.  
//      */
//     public double getState( const Vector3D& pos) const throws NotImplementedException ;
    /** Temperature in K of material at pos.  
     */
    virtual double getTemperature(const Vector3D & pos) const throw (NotImplementedException, std::exception );

    /** Pressure in P at pos.  
     */
    virtual double getPressure(const Vector3D & pos) const throw (NotImplementedException, std::exception );

    /** Radiation length of material in mm at pos.  
     */
    virtual double getRadlen(const Vector3D & pos) const throw (NotImplementedException, std::exception );

    /** Interaction length of material in mm at pos.  
     */
    virtual double getIntlen(const Vector3D & pos) const throw (NotImplementedException, std::exception );

    /** Position in local coordinate 
     */
    virtual Vector3D getLocalPosition(const Vector3D & pos) const throw (NotImplementedException, std::exception );

    /** The magnetic field vector at pos in [Tesla].  
     */
    virtual Vector3D getB(const Vector3D & pos) const throw (NotImplementedException, std::exception );

    /** The electric field vector at pos in [V/m].  
     */
    virtual Vector3D getE(const Vector3D & pos) const throw (NotImplementedException, std::exception );

    /** Names of (geant4) logical volumes in heirarchy starting at given pos ending with the world volume. 
     */
    virtual std::vector<std::string>  getListOfLogicalVolumes(const Vector3D & pos) const throw (NotImplementedException, std::exception );

    /** Names of (geant4) physical volumes in heirarchy starting at given pos ending with the world volume. 
     */
    virtual std::vector<std::string>  getListOfPhysicalVolumes(const Vector3D & pos) const throw (NotImplementedException, std::exception );

    /** Names of (geant4) region that contains the given pos.
     */
    virtual std::string getRegion(const Vector3D & pos) const throw (NotImplementedException, std::exception );

    /** True if region that contains pos is defined as a tracker.
     */
    virtual bool isTracker(const Vector3D & pos) const throw (NotImplementedException, std::exception );

    /** True if region that contains pos is defined as a calorimeter.
     */
    virtual bool isCalorimeter(const Vector3D & pos) const throw (NotImplementedException, std::exception );
}; // class
} // namespace gear
#endif /* ifndef GEAR_CGAGEARPOINTPROPERTIES_H */
