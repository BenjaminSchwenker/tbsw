#ifndef GEAR_RectangularPadRowLayoutXML_H
#define GEAR_RectangularPadRowLayoutXML_H 1


#include "gearxml/PadRowLayout2DXML.h"

#include <string>
#include <map>

namespace gear {
  
   
  /** Abstract XML handler for RectangularPadRowLayoutXML.
   * 
   * @author F. Gaede, DESY
   * @version $Id: RectangularPadRowLayoutXML.h,v 1.1 2006-05-12 12:53:55 gaede Exp $
   */
  class RectangularPadRowLayoutXML : public PadRowLayout2DXML {
    
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

#endif /* ifndef GEAR_RectangularPadRowLayoutXML_H */
