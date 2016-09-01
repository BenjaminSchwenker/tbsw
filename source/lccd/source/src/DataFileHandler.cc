#include "lccd/DataFileHandler.hh"

#include "lcio.h"

#include "IO/LCReader.h"

#include <iostream>

using namespace lcio ;

namespace lccd {


  DataFileHandler::DataFileHandler(const std::string& name ) : 
    ConditionsHandlerBase(name) {
    
    _col = 0 ;
  }
  
  DataFileHandler::~DataFileHandler() {

    delete _col ;
  }
  

  
  void DataFileHandler::updateEvent( LCEvent* evt ) {
    
    
    try{

       LCCollection* col = evt->takeCollection( name() ) ;

       if( _col != 0 )  delete _col ;
       
       _col = col ;

       _validSince =  evt->getTimeStamp()  ;
       _validTill = LCCDPlusInf ;

       notifyListeners() ;
    }
    catch( DataNotAvailableException ) { // no new data -> simply add current collection

      evt->addCollection( currentCollection() , name()  ) ;
      evt->takeCollection( name() ) ;

    }

  }
  
}
