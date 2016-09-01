#include "gearimpl/LayerLayoutImpl.h"


namespace gear{



  void LayerLayoutImpl::addLayer(double thickness, double cellSize0, double cellSize1, double absorberThickness ) {

    Layer l ;

    if( _vec.size() == 0 ) {

      l.Distance = 0.0 ;
      
    } else {

      l.Distance = _vec.back().Distance + _vec.back().Thickness ;
    }

    l.Thickness = thickness ;
    l.AbsorberThickness = absorberThickness ;
    l.CellSize0 = cellSize0 ;
    l.CellSize1 = cellSize1 ;
    _vec.push_back( l ) ;
    
  } 

  
  void LayerLayoutImpl::positionLayer(double distance, double thickness, double cellSize0, double cellSize1, double absorberThickness ) {
    
    Layer l ;

    l.Distance = distance ;
    l.Thickness = thickness ;
    l.AbsorberThickness = absorberThickness ;
    l.CellSize0 = cellSize0 ;
    l.CellSize1 = cellSize1 ;
    _vec.push_back( l ) ;
       
  }

}
