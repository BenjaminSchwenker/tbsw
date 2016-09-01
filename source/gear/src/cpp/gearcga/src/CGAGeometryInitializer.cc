#include <string>
#include <vector>

#include "gear/GEAR.h"
#include "gearcga/CGAGeometryInitializer.h"
#include "CGADefs.h"

namespace gear {

/** CGA Geometry Initializer class
 * @author G. Musat, Ecole Polytechnique
 * @version $Id: CGAGeometryInitializer.cc,v 1.1 2006-05-29 13:03:18 musat Exp $
 */
    CGAGeometryInitializer * CGAGeometryInitializer::theInitializer = NULL;

    CGAGeometryInitializer * CGAGeometryInitializer::GetCGAGeometryInitializer(
	std::string steeringFile, std::string model,
        std::string setup, std::string host, std::string user,
        std::string password) {

	if(theInitializer == NULL)
		theInitializer = new CGAGeometryInitializer(steeringFile,
			model, setup, host, user, password);
	return theInitializer;
    }

    CGAGeometryInitializer::CGAGeometryInitializer(std::string steeringFile, 
	std::string model, std::string setup, std::string host, 
	std::string user, std::string password) {

	CGAInit(steeringFile.c_str(), model.c_str(), setup.c_str(), 
		host.c_str(), user.c_str(), password.c_str());

    }

} // namespace gear
