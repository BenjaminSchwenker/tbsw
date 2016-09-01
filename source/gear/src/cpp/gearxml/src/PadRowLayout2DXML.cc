#include "gearxml/PadRowLayout2DXML.h"

#include "gearxml/FixedPadSizeDiskLayoutXML.h"
#include "gearxml/FixedPadAngleDiskLayoutXML.h"
#include "gearxml/RectangularPadRowLayoutXML.h"

#include "gearimpl/FixedPadSizeDiskLayout.h"
#include "gearimpl/FixedPadAngleDiskLayout.h"
#include "gearimpl/RectangularPadRowLayout.h"

#include <typeinfo>
#include <iostream>

namespace gear {
  
  PadRowLayout2DXML* PadRowLayout2DXML::getHandler( const std::string& typeName ) {

    return PadRowLayout2DXML::map()[ typeName ] ;
  }

  PadRowLayout2DXML* PadRowLayout2DXML::getHandler(const PadRowLayout2D* l ) {
    
    std::string typeName("UNKNOWN") ;
    
    if( dynamic_cast<const FixedPadSizeDiskLayout*>(l) != 0 ) typeName = "FixedPadSizeDiskLayout" ;
    if( dynamic_cast<const FixedPadAngleDiskLayout*>(l) != 0 ) typeName = "FixedPadAngleDiskLayout" ;
    if( dynamic_cast<const RectangularPadRowLayout*>(l) != 0 ) typeName = "RectangularPadRowLayout" ;
    
    return getHandler( typeName ) ;
  }
  

  PadRowLayout2DXML::Map& PadRowLayout2DXML::map() {
    
    static  PadRowLayout2DXML::Map theMap ;
    static FixedPadSizeDiskLayoutXML  aFixedPadSizeDiskLayoutXML ;
    static FixedPadAngleDiskLayoutXML  aFixedPadAngleDiskLayoutXML ;
    static RectangularPadRowLayoutXML aRectangularPadRowLayoutXML ;

    if( theMap.size() == 0 ) {
      
      theMap["FixedPadSizeDiskLayout"] = & aFixedPadSizeDiskLayoutXML ;
      theMap["FixedPadAngleDiskLayout"] = & aFixedPadAngleDiskLayoutXML ;
      theMap["RectangularPadRowLayout"] = & aRectangularPadRowLayoutXML ;

    }
    return theMap ;
  }

  
//   std::string PadRowLayout2DXML::typeName( const PadRowLayout2D* l ) {
    
//     if( dynamic_cast<const FixedPadSizeDiskLayout*>(l) != 0 ) return  "FixedPadSizeDiskLayout" ;
//     else  
//       return "UNKNOWN" ;
//   }
}

