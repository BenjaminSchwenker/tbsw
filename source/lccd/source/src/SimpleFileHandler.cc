#include "lccd/SimpleFileHandler.hh"

#include "lcio.h"

#include "IO/LCReader.h"

#include <iostream>

using namespace lcio ;

namespace lccd {


  SimpleFileHandler::SimpleFileHandler(const std::string& fileName, 
				       const std::string& name,
				       const std::string& inputCollection ) : 
    ConditionsHandlerBase(name) ,
    _fileName( fileName ) ,
    _inputCollection( inputCollection ),
    _firstCall( true )  {
    
    // use the same name for the input collection
    // and the data set name 
    if( _inputCollection.size() == 0 ) 
      _inputCollection = name ;

    _col = 0 ;

  }
  
  SimpleFileHandler::~SimpleFileHandler() {
//     std::cout << " in SimpleFileHandler::~SimpleFileHandler() " << std::endl ;
    delete _col ;
  }
  
  
  void SimpleFileHandler::update( LCCDTimeStamp timestamp ) {
    
//     if( timestamp >= _validTill ) { // true for first call only ! 
    if( _firstCall ) {

      _firstCall = false ;
      

//       std::cout << "SimpleFileHandler::update: reading constants from file for time stamp : " 
// 		<< timestamp << std::endl ;
      
      
      LCReader* lcReader = LCFactory::getInstance()->createLCReader() ;
      
      lcReader->open( _fileName ) ;
      
      LCEvent* evt = lcReader->readNextEvent() ;
      
      if( evt == 0 ) 
	throw Exception("SimpleFileHandler::update: No event found in file: "+_fileName) ;
      
      if( _col != 0 ) 
	delete _col ;

      _col = evt->takeCollection( _inputCollection ) ;
      
      // set infinite validity range
      _validSince =  LCCDMinusInf ;
      _validTill = LCCDPlusInf ;
      
      notifyListeners() ;
      
      lcReader->close() ;
      delete lcReader ;
      
    }
  }
  
}
