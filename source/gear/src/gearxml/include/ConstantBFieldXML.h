#ifndef GEAR_ConstantBFieldXML_H
#define GEAR_ConstantBFieldXML_H 1


#include "gearxml/XMLHandler.h"
#include "gearimpl/ConstantBField.h"

#include <string>


namespace gear {
  
  
  class TiXmlNode ;
  
  
  /** XML handler for ConstantBField
   * 
   * @author F. Gaede, DESY
   * @version $Id: ConstantBFieldXML.h,v 1.1 2007-04-28 21:45:45 gaede Exp $
   */
  class ConstantBFieldXML : public XMLHandler {
    
  public: 
    
    /** Creates an XML node for the given parameters 
     */
    virtual TiXmlElement toXML( const GearParameters & parameters ) const ;
    
    
    /** Creates the appropriate GearParameters subclass from the given
     *  XML element (node) 
     */
    virtual GearParameters* fromXML( const TiXmlElement* xmlElement, GearMgr* gearMgr=0) const ;
    

  protected:
    
    //  std::string getAttribute(const  TiXmlNode* node , const std::string& name ) const ;    

    }; // class
  
} // namespace gear

#endif /* ifndef GEAR_ConstantBFieldXML_H */
