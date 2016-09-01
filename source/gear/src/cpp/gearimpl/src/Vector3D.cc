#include "gearimpl/Vector3D.h"
//#include <iostream>

namespace gear{


std::ostream & operator << (std::ostream & os, const Vector3D &v){
  
  os << "  ( " << v[0] 
     << ", " << v[1]
     << ", " << v[2]
//   os << "  ( " << v.x() 
//      << ", " << v.y()
//      << ", " << v.z()
     << " ) -  [ phi: " << v.phi()
     << " , rho: " << v.rho() << " ] "  
     << "  [ theta: " << v.theta()
     << " , r: " << v.r() << " ] "  
     << std::endl ;

  return os ;
}

}
