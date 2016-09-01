#ifndef GEAR_FixedPadAngleDiskLayoutXML_H
#define GEAR_FixedPadAngleDiskLayoutXML_H 1


#include "gearxml/PadRowLayout2DXML.h"

#include <string>
#include <map>

namespace gear {
  
   
  /** Abstract XML handler for FixedPadAngleDiskLayoutXML.
   * 
   * @author F. Gaede, DESY
   * @version $Id: FixedPadAngleDiskLayoutXML.h,v 1.1 2008-10-30 14:27:06 gaede Exp $
   */
  class FixedPadAngleDiskLayoutXML : public PadRowLayout2DXML {
    
  public: 
    
    /** Creates an XML node for the given  PadRowLayout2D.
     */
    virtual TiXmlElement toXML( const PadRowLayout2D* layout ) const ;
    
    
    /** Creates the appropriate PadRowLayout2D subclass from the given
     *  XML element (node) 
     */
    virtual PadRowLayout2D* fromXML( const TiXmlElement* xmlElement) const ;

  }; // class
} // namespace gear

#endif /* ifndef GEAR_FixedPadAngleDiskLayoutXML_H */
