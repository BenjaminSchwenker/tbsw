#ifndef ConditionsMap_h
#define ConditionsMap_h

#include "lcio.h"
#include "EVENT/LCCollection.h"
#include "Exceptions.h"

#include "lccd/IConditionsChangeListener.hh"

#include <map>
#include <typeinfo>
#include <sstream>

#include <iostream>

namespace lccd {
  
  /** Template class for maps of conditions data. The class LCCONDOBJECT
   *  has to have a constructor that takes an LCObject* as argument.
   *  Typically this is the case for all LCGenericObject subclasses.
   * 
   *  @author F.Gaede, DESY
   *  @version $Id: ConditionsMap.hh,v 1.9 2008-11-04 10:32:20 meyern Exp $
   */
  template <class KEY, class LCCONDOBJECT>
  class ConditionsMap : public lccd::IConditionsChangeListener {

  public:
    
    typedef typename std::map< KEY, LCCONDOBJECT >::iterator MapIter ;

    // PMF: pointer to member function that returns the key
    typedef KEY (LCCONDOBJECT::*PMF)() ;  

    // PMF: pointer to member function that returns the key if it is
    // declared const
    typedef KEY (LCCONDOBJECT::*PMFC)() const ;  

    /** Constructor: provide the pointer to the member fuction that
        returns the key */
    ConditionsMap( PMF pmf ) : _pmf(pmf) { _isPMFconst = false; }

    /** Constructor: provide the pointer to the member fuction that
        returns the key if it is declared const */
    ConditionsMap( PMFC pmfc ) : _pmfc(pmfc) { _isPMFconst = true; }


    /** The d'tor.*/
    virtual ~ConditionsMap() {} 


    /** The map.*/
    const std::map< KEY, LCCONDOBJECT> & map() { return _map ; }
    
    
    /** Returns a reference to the conditions object for key 
     */
    LCCONDOBJECT& find( KEY key) {

      MapIter it = _map.find( key ) ;

      if( it == _map.end() ){
	std::stringstream err ;
	err << "ConditionsMap::find: no entry for key: " << key ;
	throw lcio::Exception( err.str() ) ;
      }
      else
	return it->second ;
    }


    /** Repopulate the map with new conditions data.
     */ 
    void conditionsChanged( lcio::LCCollection* col ) {
      
      _map.clear() ;
      
      // add all elements in the collection to the map
      for( int i=0; i< col->getNumberOfElements() ; i++ ){
	
	// this requires that the c'tor exists in LCCONDOBJECT !
	LCCONDOBJECT obj( col->getElementAt( i ) ) ;
	
	
        if(_isPMFconst == false) {

          _map.insert( std::make_pair( (obj.*_pmf)() , obj)  ) ;
          
        } else if(_isPMFconst == true) {

          _map.insert( std::make_pair( (obj.*_pmfc)() , obj)  ) ;

        }

	// call member function for key  ^^^^^^^  ,e.g.
	// _map.insert( std::make_pair( obj.getCellID() , obj)  ) ;
      }

//       print( std::cout ) ;
    }

    /** Print the conditions map to the specified output stream.<br>
     *  Debug method - prints the key and the id of the stored objects.
     */
    
    void print(  std::ostream& os ) {
      
      os << "ConditionsMap"  
//    	 << std::typeid(keyType).name() << "," 
//   	 << typeid( typename LCCONDOBJECT ).name() 
	 << std::endl ;
      
      //   for( std::map< KEY, LCCONDOBJECT >::iterator it = _map.begin() ;
      for( MapIter it = _map.begin() ;
	   it != _map.end() ;
	   it++ ){
	
	os << "  key: " <<  it->first  << " - [" << it->second.id() << "]" 
	   << " - " << typeid( it->second ).name() 
	   << std::endl ;
      }
      
    }
    

  protected:
    
    /** No default c'tor */
    ConditionsMap() { }
 
    std::map< KEY, LCCONDOBJECT > _map ;
    PMF _pmf ;
    PMFC _pmfc ;

    bool _isPMFconst;
};

}  //end namespace

#endif // ConditionsMap_h
