
#include "gearimpl/GearMgrImpl.h"
#include "gearimpl/GearParametersImpl.h"
#include "gear/SiPlanesParameters.h"
#include "gear/GearPointProperties.h"
#include "gear/GearDistanceProperties.h"
#include "gear/BField.h"

namespace gear{


  GearMgrImpl::GearMgrImpl() :
    _siplanesParameters(0) ,
    _pointProperties(0) ,
    _distanceProperties(0) ,
    _bField(0),
    _detectorName(""){
  }
  
  GearMgrImpl::~GearMgrImpl() {
    
    // clean up all parameters
    if( _siplanesParameters ) delete _siplanesParameters ;
    if( _pointProperties ) delete _pointProperties ;
    if( _distanceProperties ) delete _distanceProperties ;
    if( _bField  ) delete  _bField ;
    
    
    ParameterMap::iterator it_end = _map.end() ;

    for( ParameterMap::iterator it = _map.begin() ; it != it_end ; ++ it ) {
      delete it->second ;
    }
    
  }
  

  const std::string& GearMgrImpl::getDetectorName() const    

    throw (UnknownParameterException, std::exception ) { 

    if( _detectorName.size() == 0 )
      throw UnknownParameterException( "No DetectorName set ") ;


    return _detectorName ; 
  }


  const GearParameters & GearMgrImpl::getGearParameters(const std::string & key) const 
    
    throw (UnknownParameterException, std::exception ) {
    
    ParameterMap::const_iterator it = _map.find( key ) ;
    if( it == _map.end() || it->second == 0 )
      throw UnknownParameterException( "No parameters set for : " + key ) ;
    return * it->second ;
    
  }   
  
  const BField & GearMgrImpl::getBField() const
    throw (UnknownParameterException, std::exception ) {
    
    if( _bField == 0 )
      throw UnknownParameterException( "No BField set ") ;

    return  *_bField ;

  }


  const SiPlanesParameters & GearMgrImpl::getSiPlanesParameters() const
    throw (UnknownParameterException, std::exception ) {

    if( _siplanesParameters == 0 )
      throw UnknownParameterException( "No SiPlanesParameters set ") ;

    return *_siplanesParameters ;

  }

  
  const GearPointProperties & GearMgrImpl::getPointProperties() const 
    throw (NotImplementedException, std::exception ) {

    if( _pointProperties == 0 )
      throw UnknownParameterException( "No PointProperties set or implemented ") ;
    
    return  *_pointProperties ;
  }
  
  
  
  const GearDistanceProperties & GearMgrImpl::getDistanceProperties() const 
    throw (NotImplementedException, std::exception ) {

    if( _distanceProperties == 0 )
      throw UnknownParameterException( "No DistanceProperties set or implemented ") ;
    
    return  *_distanceProperties ;

  }

  void GearMgrImpl::setGearParameters( const std::string & key, GearParameters* parameters ) {

    if( parameters == 0 )   // don't allow  null pointers 
      return  ;    
    
    ParameterMap::iterator it = _map.find( key ) ;

    if( it != _map.end() ) {
      
      delete it->second ;
      it->second = parameters ;
      
    } else {
      
      _map[ key ] = parameters ; 
      
    }
    
    
  }
  
  

  void GearMgrImpl::setBField( BField* b){
    
    _bField = b ;
  }


  void GearMgrImpl::setSiPlanesParameters( SiPlanesParameters* siplanesParameters ) {

    _siplanesParameters = siplanesParameters ;
  }

 

  void GearMgrImpl::setPointProperties( GearPointProperties* pointProperties) {
    
    _pointProperties = pointProperties ;
  }

  void GearMgrImpl::setDistanceProperties( GearDistanceProperties* distanceProperties) {
 
   _distanceProperties = distanceProperties ;
  }
  
  const std::vector<std::string>& GearMgrImpl::getGearParameterKeys() const {

    _keys.clear() ;
    _keys.reserve( _map.size() ) ;

    for( ParameterMap::const_iterator it = _map.begin() ; it != _map.end() ; ++it ){
      _keys.push_back( it->first ) ;
    }
    return _keys ;
  }

}
