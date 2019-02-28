#ifndef GEAR_Util_H
#define GEAR_Util_H 1


#include "gear/GearParameters.h"
#include "gear/BField.h"
#include "gear/SiPlanesParameters.h"
#include "gear/GearMgr.h"
#include "gear/GEAR.h"

#include <iostream>

namespace gear {
    

  std::ostream& operator<< (  std::ostream& s,  const GearMgr& m ) ;


  std::ostream& operator<< (  std::ostream& s,  const GearParameters& p ) ;
    

  std::ostream& operator<< (  std::ostream& s,  const BField& b ) ;


  std::ostream& operator<< (  std::ostream& s,  const SiPlanesParameters& p ) ;


} // namespace gear

#endif /* ifndef GEAR_Util_H */
