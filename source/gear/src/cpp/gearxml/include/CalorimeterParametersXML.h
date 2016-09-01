#ifndef GEAR_CalorimeterParametersXML_H
#define GEAR_CalorimeterParametersXML_H 1


#include "gearxml/XMLHandler.h"

#include <string>


namespace gear {
  
  
  class TiXmlNode ;
  
  
  /** XML handler for CalorimeterParameters.
   * 
   * @author F. Gaede, DESY
   * @version $Id: CalorimeterParametersXML.h,v 1.4 2006-04-13 10:16:20 lippe Exp $
   */
  class CalorimeterParametersXML : public XMLHandler {
    
  public: 
    
    /** Creates an XML node for the given parameters 
     *  @author R.Lippe, DESY
     */
    virtual TiXmlElement toXML( const GearParameters & parameters ) const ;
    
    
    /** Creates the appropriate GearParameters subclass from the given
     *  XML element (node) 
     */
    virtual GearParameters* fromXML( const TiXmlElement* xmlElement, GearMgr* gearMgr=0) const ;
    

  protected:
    
    //  std::string getAttribute(const  TiXmlNode* node , const std::string& name ) const ;    

  private:

    /** Private function that compares two double values and returns true if they are almost equal
     *  precision to be equal is relativly given in EPSILON
     */
    bool isEqual( const double valueOne, const double valueTwo ) const ;

    }; // class
  
} // namespace gear

#endif /* ifndef GEAR_CalorimeterParametersXML_H */
