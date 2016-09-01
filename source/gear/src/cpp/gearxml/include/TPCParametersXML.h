#ifndef GEAR_TPCParametersXML_H
#define GEAR_TPCParametersXML_H 1


#include "gearxml/XMLHandler.h"

#include <string>


namespace gear {
  
  
  class TiXmlNode ;
  
  
  /** XML handler for TPCParameters.
   * 
   * @author F. Gaede, DESY
   * @version $Id: TPCParametersXML.h,v 1.2 2005-11-25 16:08:15 gaede Exp $
   */
  class TPCParametersXML : public XMLHandler {
    
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

#endif /* ifndef GEAR_TPCParametersXML_H */
