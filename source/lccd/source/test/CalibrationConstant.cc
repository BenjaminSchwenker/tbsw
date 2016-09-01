#include "CalibrationConstant.hh"


void CalibrationConstant::print(  std::ostream& os ){

  os << " cellId: ["  << std::hex  << getCellID() << "]" <<  std::dec 
     << ", offset: " <<  getOffset() 
     << ", gain: "   <<  getGain() 
     << std::endl ;
  
}
