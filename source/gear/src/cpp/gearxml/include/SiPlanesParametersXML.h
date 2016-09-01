#ifndef GEAR_SiPlanesParametersXML_H
#define GEAR_SiPlanesParametersXML_H 1

#include "gearxml/XMLHandler.h"

#include <string>

namespace gear {

  class TiXmlNode ;
  
  /** XML handler for SiPlanesParameters.
   *
   * @author T. Klimkovich, DESY
   * @version $Id $
   */
  class SiPlanesParametersXML : public XMLHandler {

  public:
    
    /** Creates an XML node for the given parameters
     */
    virtual TiXmlElement toXML( const GearParameters & parameters ) const ;

    /** Creates the appropriate GearParameters subclass from the given
     *  XML element (node)
     */
    virtual GearParameters* fromXML( const TiXmlElement* xmlElement, GearMgr* gearMgr=0) const ;


  }; //class

} // namespace

#endif
