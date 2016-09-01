
#include "gearimpl/Util.h"
#include "gearxml/GearXML.h"
#include "gear/GearMgr.h"
#include "gear/GEAR.h"
#include "gear/SiPlanesParameters.h"

#ifdef CGA
#include "gearcga/CGAGearDistanceProperties.h"
#include "gearcga/CGAGearPointProperties.h"
#endif

#include <iostream>
#include <assert.h>

#include <exception>
#include <cstdlib>

#include <sstream>
#include <fstream>

#include <math.h>

#ifdef GEAR_USE_AIDA  
#include "AIDA/AIDA.h"
using namespace AIDA ;
#endif



using namespace gear ;

void testSiPlanes( const SiPlanesParameters& siplanesParams ) ;

void gear_unexpected(){

  try {

    throw ;

  } catch( std::exception& e) {

    std::cout << " A runtime error has occured : "
              << e.what()
              << std::endl
              << " the program will have to be terminated - sorry." << std::endl ;
    exit(1) ;
  }
}



/** Testprogram for gear classes.
 * 
 */

int main(int argc, char**argv){
  
  
  std::set_unexpected( gear_unexpected ) ;
  std::set_terminate( gear_unexpected ) ;
  
  if( argc < 2 ) {
    std::cout << " testgear:  Testprogram for SiPlanes detector in gear. " << std::endl 
	      << " usage:     testgear GearOutput.xml  " << std::endl ;
    exit(1) ;
   }
  
  std::string fileName( argv[1] ) ;
  
  GearXML gearXML( fileName ) ;
  
  GearMgr* gearMgr = gearXML.createGearMgr() ;
  
  
  std::cout << " testSiPlanesgear - instantiated  GearMgr from file " << fileName 
	    << std::endl ;

  const SiPlanesParameters& vp = 
    dynamic_cast<const SiPlanesParameters&>( gearMgr->getSiPlanesParameters() ) ;
 

  std::cout << " SiPlanes parameters : " << std::endl 
	    << vp 
	    << std::endl ;

}
