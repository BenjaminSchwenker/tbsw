#include "gearimpl/CalorimeterParametersImpl.h"

namespace gear{


    CalorimeterParametersImpl::CalorimeterParametersImpl(double rMin, double zMax, int symOrder, double phi0) :
	_type( CalorimeterParameters::BARREL) ,
	_sym( symOrder ) ,
	_phi( phi0 ) 
	
    {
	_extent.resize(4) ;
	_extent[0] = rMin ;
	_extent[1] = rMin ; // dynamically computed from layers
	_extent[2] = 0.0 ;
	_extent[3] = zMax ;
    }



    CalorimeterParametersImpl::CalorimeterParametersImpl(double rMin, double rMax, double zMin, int symOrder, double phi0) :
	_type( CalorimeterParameters::ENDCAP) ,
	_sym( symOrder ) ,
	_phi( phi0 ) 
	
    {
	_extent.resize(4) ;
	_extent[0] = rMin ;
	_extent[1] = rMax ; 
	_extent[2] = zMin ;
	_extent[3] = zMin ;// dynamically computed from layers
    }

  const std::vector<double>& CalorimeterParametersImpl::getExtent() const { 

    int nLayer = _layout.getNLayers() ;


    if(  nLayer > 0 ) {

      if( _type == CalorimeterParameters::BARREL ) { 

	// the first layer starts at rMin and the last layer ends at its distance from the origin plus its thickness 
	
	_extent[1] = _layout.getDistance( nLayer-1 )  +  _layout.getThickness( nLayer-1 )  ;
	
      } else {
	
	_extent[3] = _layout.getDistance( nLayer-1 )  +  _layout.getThickness( nLayer-1 )  ;
      }
    }

    return _extent ; 
  }
  
  
}
