#ifndef GEAR_GearParametersXML_H
#define GEAR_GearParametersXML_H 1


#include "gearxml/XMLHandler.h"

#include "gearimpl/GearParametersImpl.h"

#include <string>


namespace gear {
  
  
  class TiXmlNode ;
  
  
  /** XML handler for GearParameters.
   * 
   * @author F. Gaede, DESY
   * @version $Id: GearParametersXML.h,v 1.2 2005-11-25 16:08:15 gaede Exp $
   */
  class GearParametersXML : public XMLHandler {
    
  public: 
    
    /** Creates an XML node for the given parameters 
     */
    virtual TiXmlElement toXML( const GearParameters& parameters ) const ;
    
    
    /** Creates the appropriate GearParameters subclass from the given
     *  XML element (node) 
     */
    virtual GearParameters* fromXML( const TiXmlElement* xmlElement, GearMgr* gearMgr=0) const ;
    

    /** Static helper function that can be used by other subclass handlers to read parameters 
     */
    static void setParametersFromXML( const TiXmlElement* xmlElement, GearParametersImpl* gearParams );


    /** Static helper function that can be used by other subclass handlers to create XML for parameters 
     */
    static void getXMLForParameters( TiXmlElement* xmlElement, const GearParameters* gearParams );


  protected:
    
    //  std::string getAttribute(const  TiXmlNode* node , const std::string& name ) const ;    

    
    
  }; // class
  
} // namespace gear

#endif /* ifndef GEAR_GearParametersXML_H */
