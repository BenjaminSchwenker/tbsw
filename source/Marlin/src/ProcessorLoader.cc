#ifndef MARLIN_NO_DLL
#include "marlin/ProcessorLoader.h"

#include <dlfcn.h>
#include <iostream>
#include <cstdlib>

namespace marlin{


ProcessorLoader::ProcessorLoader(
        lcio::StringVec::const_iterator first, 
        lcio::StringVec::const_iterator last ) {


    _loadError=false;
    lcio::StringVec::const_iterator current = first ;

    while( current != last ){

        //std::string libName( "lib" + *current + ".so" ) ;
        std::string libName( *current ) ;
        std::cout << "<!-- Loading shared library : " << libName << " -->" << std::endl ;
        
        //void* libPointer  = dlopen( libName.c_str() , RTLD_NOW) ;
        //void* libPointer  = dlopen( libName.c_str() , RTLD_LAZY ) ;
        //void* libPointer  = dlopen( libName.c_str() , RTLD_NOW | RTLD_GLOBAL) ;
        void* libPointer  = dlopen( libName.c_str() , RTLD_LAZY | RTLD_GLOBAL) ;

        if( libPointer == 0 ){
            std::cout << std::endl << "<!-- ERROR loading shared library : " << libName << std::endl
                      << "    ->    "   << dlerror() << " -->" << std::endl ;
            _loadError=true;
        }
        else{
            _libs.push_back( libPointer ) ;
        }
        ++current ;
    }
}


ProcessorLoader::~ProcessorLoader() {


  char * s = getenv("MARLIN_DEBUG" ) ;
  
  // do not unload processors if $MARLIN_DEBUG is set to 1
  // useful for debugging with valgrind (https://jira.slac.stanford.edu/browse/MAR-45)
  if( std::string("1").compare( (s?s:"0") ) == 0 ) {
    std::cout << std::endl << "<!-- MARLIN_DEBUG=1 set in your environment - skip unloading processors --> " << std::endl ;
  } else {
    for( LibVec::iterator it = _libs.begin() ; it != _libs.end() ; ++it ) {
        dlclose( *it ) ;
    }
  }
}


} // namespace marlin
#endif //MARLIN_NO_DLL
