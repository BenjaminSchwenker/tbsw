#include <string>
#include <vector>
#include <iostream>

#include "gear/GEAR.h"
#include "gearcga/CGAGearPointProperties.h"
#include "gearcga/CGAGeometryInitializer.h"
#include "CGADefs.h"

namespace gear {

    CGAGearPointProperties::CGAGearPointProperties(std::string steer,
                        std::string model, std::string setup, std::string host,
                        std::string user, std::string password) {
                                                                                
        CGAGeometryInitializer::GetCGAGeometryInitializer(steer, model,
        setup, host, user, password);
                                                                                
    }

    //vec getPosition returns the current position (if ID is known)
    /** The cellID of the the sensitive detector at pos.
     */
    long64 CGAGearPointProperties::getCellID(const Vector3D & pos) const throw (NotImplementedException, std::exception ) {
	
	double x, y, z;
	int flag;
	x = pos[0];
	y = pos[1];
	z = pos[2];
	cell_ids cellIDs = CGAGetCellId(x, y, z, flag, 0.0, 0.0, 1.0);

	char text[2048];
	sprintf(text, "CGAGearPointProperties::getCellID: Point (%f, %f, %f) ",
			 x , y , z);
	if(cellIDs.id0 < 0) {
		throw gear::Exception((const std::string)(text) + 
			" is outside of world volume." );
	}
	else if(flag == -1) {
		throw gear::Exception((const std::string)(text) + 
				" is outside of sensitive volume");
	}
	else if(flag == 0)
		 std::cout << text << 
			 " is in guard-ring. Returning CellID of nearest cell."
			 << std::endl;
	
	int SHIFT_A_64 = 0;
	int SHIFT_B_64 = 32;

	unsigned long long MASK_A_64 = 0x00000000FFFFFFFFLL;
	unsigned long long MASK_B_64 = 0xFFFFFFFF00000000LL;

	long64 result = (unsigned long long) (
                ((((unsigned long long)(cellIDs.id0))<<SHIFT_A_64)&MASK_A_64)|
                ((((unsigned long long)(cellIDs.id1))<<SHIFT_B_64)&MASK_B_64)
					     );

	return result;
    }

    // 	bool getMaterial  (true if defined material, false if problem)
    /** Name of material at pos.  
     */
    const std::string & CGAGearPointProperties::getMaterialName(const Vector3D & pos) const throw (NotImplementedException, std::exception ) {

	static std::string matName("unknown");

	double position[3];
	position[0] = pos[0];
	position[1] = pos[1];
	position[2] = pos[2];

	matName = CGAGetMaterialName(position);

	return matName;

    }
    /** Density in kg/m^3 at pos.  
     */
    double CGAGearPointProperties::getDensity(const Vector3D & pos) const throw (NotImplementedException, std::exception ) {

	double position[3];
        position[0] = pos[0];
        position[1] = pos[1];
        position[2] = pos[2];

	return CGAGetDensity(position);
    }

//     /** Name of material at pos.  
//      */
//     public double getState( const Vector3D& pos) const throws NotImplementedException ;
    /** Temperature in K of material at pos.  
     */
    double CGAGearPointProperties::getTemperature(const Vector3D & pos) const throw (NotImplementedException, std::exception ) {

	double position[3];
        position[0] = pos[0];
        position[1] = pos[1];
        position[2] = pos[2];

	return CGAGetTemperature(position);
    }

    /** Pressure in P at pos.  
     */
    double CGAGearPointProperties::getPressure(const Vector3D & pos) const throw (NotImplementedException, std::exception ) {

	double position[3];
        position[0] = pos[0];
        position[1] = pos[1];
        position[2] = pos[2];

	return CGAGetPressure(position);

    }

    /** Radiation length of material in mm at pos.  
     */
    double CGAGearPointProperties::getRadlen(const Vector3D & pos) const throw (NotImplementedException, std::exception ) {

	double position[3];
        position[0] = pos[0];
        position[1] = pos[1];
        position[2] = pos[2];

	return CGAGetRadLen(position);
    }

    /** Interaction length of material in mm at pos.  
     */
    double CGAGearPointProperties::getIntlen(const Vector3D & pos) const throw (NotImplementedException, std::exception ) {

	double position[3];
        position[0] = pos[0];
        position[1] = pos[1];
        position[2] = pos[2];

	return CGAGetIntLen(position);
    }

    /** Position in local coordinate 
     */
    Vector3D CGAGearPointProperties::getLocalPosition(const Vector3D & pos) const throw (NotImplementedException, std::exception ) {

	double position[3];
        position[0] = pos[0];
        position[1] = pos[1];
        position[2] = pos[2];

	std::vector<double> result = CGAGetLocalPosition(position);

	return Vector3D(result[0], result[1], result[2]);
    }

    /** The magnetic field vector at pos in [Tesla].  
     */
    Vector3D CGAGearPointProperties::getB(const Vector3D & pos) const throw (NotImplementedException, std::exception ) {

	double position[3];
        position[0] = pos[0];
        position[1] = pos[1];
        position[2] = pos[2];

	std::vector<double> result = CGAGetB(position);

	return Vector3D(result[0], result[1], result[2]);

    }

    /** The electric field vector at pos in [V/m].  
     */
    Vector3D CGAGearPointProperties::getE(const Vector3D & pos) const throw (NotImplementedException, std::exception ) {

	double position[3];
        position[0] = pos[0];
        position[1] = pos[1];
        position[2] = pos[2];

	std::vector<double> result = CGAGetE(position);

	return Vector3D(result[0], result[1], result[2]);
    }

    /** Names of (geant4) logical volumes in heirarchy starting at given pos ending with the world volume. 
     */
    std::vector<std::string>  CGAGearPointProperties::getListOfLogicalVolumes(const Vector3D & pos) const throw (NotImplementedException, std::exception ) {

	double position[3];
        position[0] = pos[0];
        position[1] = pos[1];
        position[2] = pos[2];

	return CGAGetListOfLogicalVolumes(position);
    }
    /** Names of (geant4) physical volumes in heirarchy starting at given pos ending with the world volume. 
     */
    std::vector<std::string>  CGAGearPointProperties::getListOfPhysicalVolumes(const Vector3D & pos) const throw (NotImplementedException, std::exception ) {

	double position[3];
        position[0] = pos[0];
        position[1] = pos[1];
        position[2] = pos[2];

	return CGAGetListOfPhysicalVolumes(position);
    }

    /** Names of (geant4) region that contains the given pos.
     */
    std::string CGAGearPointProperties::getRegion(const Vector3D & pos) const throw (NotImplementedException, std::exception ) {

	double position[3];
        position[0] = pos[0];
        position[1] = pos[1];
        position[2] = pos[2];

	return CGAGetRegionName(position);
    }

    /** True if region that contains pos is defined as a tracker.
     */
    bool CGAGearPointProperties::isTracker(const Vector3D & pos) const throw (NotImplementedException, std::exception ) {

	double position[3];
        position[0] = pos[0];
        position[1] = pos[1];
        position[2] = pos[2];

	return CGAIsTracker(position);
    }

    /** True if region that contains pos is defined as a calorimeter.
     */
    bool CGAGearPointProperties::isCalorimeter(const Vector3D & pos) const throw (NotImplementedException, std::exception ) {

	double position[3];
        position[0] = pos[0];
        position[1] = pos[1];
        position[2] = pos[2];

	return CGAIsCalorimeter(position);
    }

} // namespace gear
