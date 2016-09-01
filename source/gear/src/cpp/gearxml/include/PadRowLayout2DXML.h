#ifndef GEAR_PadRowLayout2DXML_H
#define GEAR_PadRowLayout2DXML_H 1


//#include "gearxml/XMLHandler.h"
//#include "gear/PadRowLayout2D.h"

#include <string>
#include <map>

namespace gear {
  
  class PadRowLayout2D ;
  class TiXmlElement ;
  class GearMgr ;

  class PadRowLayout2DXML ;

   
  /** Abstract XML handler for PadRowLayout2DXML.
   * 
   * @author F. Gaede, DESY
   * @version $Id: PadRowLayout2DXML.h,v 1.1 2006-05-12 12:53:55 gaede Exp $
   */
  class PadRowLayout2DXML {
    
    typedef std::map< std::string, PadRowLayout2DXML* > Map ;


  public: 
    
    /** Creates an XML node for the given  PadRowLayout2D.
     */
    virtual TiXmlElement toXML( const PadRowLayout2D* layout ) const=0 ;
    
    
    /** Creates the appropriate PadRowLayout2D subclass from the given
     *  XML element (node) 
     */
    virtual PadRowLayout2D* fromXML( const TiXmlElement* xmlElement) const=0 ;


    /** Get handler for given type */
    static PadRowLayout2DXML* getHandler(const std::string& typeName ) ;

    /** Get handler for given type of PadRowLayout2D */
    static PadRowLayout2DXML* getHandler(const PadRowLayout2D* l ) ;


    /** The map with available xml handler objects (PadRowLayout2DXML) 
     */
    static Map& map() ;

//     /** Helper function for getting the proper key in the handler map */
//     static std::string typeName( const PadRowLayout2D* l ) ;

  }; // class
} // namespace gear

#endif /* ifndef GEAR_PadRowLayout2DXML_H */
