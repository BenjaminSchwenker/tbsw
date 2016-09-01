#include "gearimpl/GearParametersImpl.h"

#include <algorithm>

namespace gear{
  
  GearParametersImpl::~GearParametersImpl() { /* nop */; }
  
  int GearParametersImpl::getIntVal(const std::string & key) const throw (UnknownParameterException, std::exception ) { 
    
    IntMap::const_iterator it = _intMap.find( key ) ;
    if( it == _intMap.end() )
      throw UnknownParameterException( key ) ;
    return it->second ;
  }
  
  double GearParametersImpl::getDoubleVal(const std::string & key) const throw (UnknownParameterException, std::exception ) {
    
    DoubleMap::const_iterator it = _doubleMap.find( key ) ;
    if( it == _doubleMap.end() )
      throw UnknownParameterException( key ) ;
    return it->second ;
  }
  
  const std::string & GearParametersImpl::getStringVal(const std::string & key) const 
    throw (UnknownParameterException, std::exception ) { 

    StringMap::const_iterator it = _stringMap.find( key ) ;
    if( it == _stringMap.end() )
      throw UnknownParameterException( key ) ;
    return it->second ;
  }

  const std::vector<int>  & GearParametersImpl::getIntVals(const std::string & key) const 
    throw (UnknownParameterException, std::exception ) { 

    IntVecMap::const_iterator it = _intVecMap.find( key ) ;
    if( it == _intVecMap.end() )
      throw UnknownParameterException( key ) ;
    return it->second ;
  }

  const std::vector<double>  & GearParametersImpl::getDoubleVals(const std::string & key) const 
    throw (UnknownParameterException, std::exception ) { 

    DoubleVecMap::const_iterator it = _doubleVecMap.find( key ) ;
    if( it == _doubleVecMap.end() )
      throw UnknownParameterException( key ) ;
    return it->second ;
  }

  const std::vector<std::string>  & GearParametersImpl::getStringVals(const std::string & key) const 
    throw (UnknownParameterException, std::exception ) { 

    StringVecMap::const_iterator it = _stringVecMap.find( key ) ;
    if( it == _stringVecMap.end() )
      throw UnknownParameterException( key ) ;
    return it->second ;
  }

  void GearParametersImpl::setIntVal(const std::string & key , int val ) {
    _intMap[ key ] = val ;
  }
    
  void GearParametersImpl::setDoubleVal(const std::string & key, double val ) {
    _doubleMap[ key ] = val ;
  }

    
  void GearParametersImpl::setStringVal(const std::string & key , const std::string & val) {
    _stringMap[ key ] = val ;
  }

    
  void GearParametersImpl::setIntVals(const std::string & key, const std::vector<int>& vals) {
    _intVecMap[ key ] = vals ;
  }


  void GearParametersImpl::setDoubleVals(const std::string & key, const std::vector<double>& vals) {
    _doubleVecMap[ key ] = vals ;
  }

    
  void GearParametersImpl::setStringVals(const std::string & key, const std::vector<std::string>& vals) {
    _stringVecMap[ key ] = vals ;
  }

  const std::vector<std::string>& GearParametersImpl::getIntKeys() const {

    _intKeys.clear() ;
    _intKeys.reserve( _intMap.size() ) ;

    for( IntMap::const_iterator it = _intMap.begin() ; it != _intMap.end() ; ++it ){
      _intKeys.push_back( it->first ) ;
    }
    return _intKeys ;
  }
  
  const std::vector<std::string>& GearParametersImpl::getDoubleKeys() const {   

    _doubleKeys.clear() ;
    _doubleKeys.reserve( _doubleMap.size() ) ;
    
    for( DoubleMap::const_iterator it = _doubleMap.begin() ; it != _doubleMap.end() ; ++it ){
      _doubleKeys.push_back( it->first ) ;
    }
    return _doubleKeys ;
  } 
  


  const std::vector<std::string>& GearParametersImpl::getStringKeys() const { 

    _stringKeys.clear() ;
    _stringKeys.reserve( _stringMap.size() ) ;

    for( StringMap::const_iterator it = _stringMap.begin() ; it != _stringMap.end() ; ++it ){
      _stringKeys.push_back( it->first ) ;
    }
    return _stringKeys ;
  } 
  
  const std::vector<std::string>& GearParametersImpl::getIntVecKeys() const { 

    _intVecKeys.clear() ;
    _intVecKeys.reserve( _intVecMap.size() ) ;

    for( IntVecMap::const_iterator it = _intVecMap.begin() ; it != _intVecMap.end() ; ++it ){
      _intVecKeys.push_back( it->first ) ;
    }
    return _intVecKeys ;
  } 
  
  const std::vector<std::string>& GearParametersImpl::getDoubleVecKeys() const { 

    _doubleVecKeys.clear() ;
    _doubleVecKeys.reserve( _doubleVecMap.size() ) ;
    
    for( DoubleVecMap::const_iterator it = _doubleVecMap.begin() ; it != _doubleVecMap.end() ; ++it ){
      _doubleVecKeys.push_back( it->first ) ;
    }
    return _doubleVecKeys ;
  }


  const std::vector<std::string>& GearParametersImpl::getStringVecKeys() const { 

    _stringVecKeys.clear() ;
    _stringVecKeys.reserve( _stringVecMap.size() ) ;
    
    for( StringVecMap::const_iterator it = _stringVecMap.begin() ; it != _stringVecMap.end() ; ++it ){
      _stringVecKeys.push_back( it->first ) ;
    }
    return _stringVecKeys ;
  } 
  
}
