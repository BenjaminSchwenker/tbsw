#include "lccd/LCConditionsMgr.hh"

#include "lcio.h"
#include <iostream>

namespace lccd{

  LCConditionsMgr* LCConditionsMgr::_me  = 0 ;
  

  LCConditionsMgr* LCConditionsMgr::instance(){
    if( _me == 0) 
      _me = new LCConditionsMgr ;
    return _me ;
  }
  
  LCConditionsMgr::LCConditionsMgr(){
  }
  

  void LCConditionsMgr::registerChangeListener( IConditionsChangeListener* cl, const std::string&  name) {

    HandlerMap::iterator it =  _map.find( name ) ;

    if( it != _map.end() ) {
      it->second->registerChangeListener( cl ) ;
    }      
    else
      throw lcio::Exception(" LCConditionsMgr::registerChangeListener: no handler "
			    " for name: " + name ) ; 
  }
  
  void LCConditionsMgr::removeHandler(const std::string&  name ) {

    HandlerMap::iterator it =  _map.find( name ) ;

    if( it != _map.end() ) {
      delete it->second ;
      _map.erase( it )  ;
    }
  }
  
  void LCConditionsMgr::updateEvent( lcio::LCEvent* evt ) {
    
    for( HandlerMap::iterator it =  _map.begin() ; it != _map.end() ; it++ ){
      
      try {
	it->second->updateEvent( evt ) ;
      }
      catch(lcio::Exception e ){
	std::cout << " LCConditionsMgr::updateEvent: an exception occured: " << e.what() << std::endl ;
      }
    }
  }
  
  void LCConditionsMgr::update( LCCDTimeStamp timestamp ) {

    for( HandlerMap::iterator it =  _map.begin() ; it != _map.end() ; it++ ){
      
      try {
      it->second->update( timestamp ) ;
      }
      catch(lcio::Exception e ){
	std::cout << " LCConditionsMgr::update: an exception occured: " << e.what() << std::endl ;
      }
    }

  }

  IConditionsHandler* LCConditionsMgr::getHandler( const std::string&  name ) {

    HandlerMap::iterator it =  _map.find( name ) ;

    if( it != _map.end() ) 
      return it->second ;
    else 
      return 0 ;  
  }
  
  
  void LCConditionsMgr::registerHandler( const std::string&  name, IConditionsHandler* handler) {
    
    HandlerMap::iterator it =  _map.find( name ) ;

    if( it != _map.end() ) {

      throw lcio::Exception(" LCConditionsMgr::registerHandler: a handler allready exists"
			    " for name " + name ) ; 
    }

    _map[ name ] = handler  ;
  }
  

  void LCConditionsMgr::clear() {

    for( HandlerMap::iterator it =  _map.begin() ; it != _map.end() ; it++ ){

      delete it->second ;
    }

    _map.clear() ;
  }


}

