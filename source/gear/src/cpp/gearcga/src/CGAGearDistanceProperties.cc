#include <string>
#include <vector>

#include "gear/GEAR.h"
#include "CGADefs.h"
#include "gearcga/CGAGearDistanceProperties.h"
#include "gearcga/CGAGeometryInitializer.h"

namespace gear {

    CGAGearDistanceProperties::CGAGearDistanceProperties(std::string steer, 
			std::string model, std::string setup, std::string host,
			std::string user, std::string password) {
	
	CGAGeometryInitializer::GetCGAGeometryInitializer(steer, model,
        setup, host, user, password);

    }
                                                                                
    void CGAGearDistanceProperties::beamOn(const Vector3D & p0, const Vector3D & p1) const{

	double initial[3], final[3], direction[3];
	// p0 and p1 are in mm; CGABeamOn requires cm
	for(unsigned int i=0; i<3; i++) {
		initial[i] = (p0[i])*0.1;
		final[i]   = (p1[i])*0.1;
		direction[i] = final[i] - initial[i];
	}

	char particle[]="geantino";
	CGABeamOn(initial, final, direction, particle, 20, 1);
    }

    /** List of matrial names along the distance between [p0,p1] .
     */
    const std::vector<std::string> & CGAGearDistanceProperties::getMaterialNames(const Vector3D & p0, const Vector3D & p1) const throw (NotImplementedException, std::exception ) {

	std::vector<std::string> volNames;
        static std::vector<std::string> matNames;
	std::vector<double> distance;
        std::vector<double> x;
	std::vector<double> y;
	std::vector<double> z;
        std::vector<double> nbX0;
	std::vector<double> nInterLen;
	
	beamOn(p0, p1);

	CGAGetStepsForJava(volNames,matNames,distance,x,y,z,nbX0,nInterLen);
                                                                                
	return matNames;
    }

    /** List of matrial thicknesses in mm along the distance between [p0,p1] - runs parallel to the array
     *  returned by  getMaterialNames().
     */
    const std::vector<double>  & CGAGearDistanceProperties::getMaterialThicknesses(const Vector3D & p0, const Vector3D & p1) const throw (NotImplementedException, std::exception ) {

	std::vector<std::string> volNames;
        std::vector<std::string> matNames; 
	static std::vector<double> distance;
        std::vector<double> x;
	std::vector<double> y;
	std::vector<double> z;
        std::vector<double> nbX0;
	std::vector<double> nInterLen;
	
	beamOn(p0, p1);

	CGAGetStepsForJava(volNames,matNames,distance,x,y,z,nbX0,nInterLen);
                                                                                
	return distance;
    }

    /** The number of radiation lengths along the distance between [p0,p1] .
     */
    double CGAGearDistanceProperties::getNRadlen(const Vector3D & p0, const Vector3D & p1) const throw (NotImplementedException, std::exception ) {

	std::vector<std::string> volNames;
        std::vector<std::string> matNames; 
	std::vector<double> distance;
        std::vector<double> x;
	std::vector<double> y;
	std::vector<double> z;
        std::vector<double> nbX0;
	std::vector<double> nInterLen;
	
	beamOn(p0, p1);

	CGAGetStepsForJava(volNames,matNames,distance,x,y,z,nbX0,nInterLen);

	double nRadLen=0;
	for(unsigned int i=0; i < nbX0.size(); i++)
		nRadLen += nbX0[i];
                                                                             
	return nRadLen;
    }


    /** The number of interaction lengths along the distance between [p0,p1] .
     */
    double CGAGearDistanceProperties::getNIntlen(const Vector3D & p0, const Vector3D & p1) const throw (NotImplementedException, std::exception ) {

	std::vector<std::string> volNames;
        std::vector<std::string> matNames; 
	std::vector<double> distance;
        std::vector<double> x;
	std::vector<double> y;
	std::vector<double> z;
        std::vector<double> nbX0;
	std::vector<double> nbIL;
	
	beamOn(p0, p1);

	CGAGetStepsForJava(volNames,matNames,distance,x,y,z,nbX0,nbIL);
                                                                                
	double nInterLen = 0;
	for(unsigned int i=0; i<nbIL.size(); i++)
		nInterLen += nbIL[i];

	return nInterLen;
    }


    /** The integrated magnetic field along  the distance between [p0,p1] in Tesla*mm.  
     */
    double CGAGearDistanceProperties::getBdL(const Vector3D & p0, const Vector3D & p1) const throw (NotImplementedException, std::exception ) {

	double start[3], end[3];
	for(int i=0; i<3; i++) {
		start[i] = p0[i];
		end[i]   = p1[i];
	}

	return CGAGetBdl(start, end);
    }

    /** The integrated electric field along  the distance between [p0,p1] in  mVolt.  
     */
    double CGAGearDistanceProperties::getEdL(const Vector3D & p0, const Vector3D & p1) const throw (NotImplementedException, std::exception ) {

	double start[3], end[3];
	for(int i=0; i<3; i++) {
		start[i] = p0[i];
		end[i]   = p1[i];
	}

	return CGAGetEdl(start, end);
    }
} // namespace gear
