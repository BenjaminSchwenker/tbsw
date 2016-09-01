#ifndef GEAR_XMLHandlerMgr_H
#define GEAR_XMLHandlerMgr_H 1


#include <map>

#include "gearxml/XMLHandler.h"


namespace gear {
  
  typedef std::map< std::string, XMLHandler* >   XMLHandlerMap ;
  
  
  /** Manager for  XML handlers of GEAR classes.
   * 
   * @author F. Gaede, DESY
   * @version $Id: XMLHandlerMgr.h,v 1.2 2007-11-06 09:06:00 gaede Exp $
   */
  
  class XMLHandlerMgr {
    
    
  public: 
    
    static XMLHandlerMgr* instance() ;
    
    
    /** XMLHandler for given type. 
     */
    const XMLHandler* getHandler( const std::string& type )  ;
    
    
    /** Register XMLHandler for type. 
     */
    const void setHandler( const std::string& type , XMLHandler* handler )  ;
   
    ~XMLHandlerMgr() ;

  protected:

    XMLHandlerMap _map ;

  private:
    XMLHandlerMgr() ;
    //    static XMLHandlerMgr* _me ;

  }; // class
  
} // namespace gear

#endif /* ifndef GEAR_XMLHandlerMgr_H */
