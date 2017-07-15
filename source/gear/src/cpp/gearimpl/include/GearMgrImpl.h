#ifndef GEAR_GearMgrImpl_H
#define GEAR_GearMgrImpl_H 1

#include <string>

#include "gear/GearMgr.h"
#include "gear/GearParameters.h"


namespace gear {
    
    
  /** Manager class that returns the Gear classes for the  relevant subdetectors.
   * 
   *  Based on ideas discussed at the 2004 Argonne Simulation Workshop as summarized by T.Behnke.
   *
   * @author F. Gaede, DESY
   * @version $Id: GearMgrImpl.h,v 1.11 2008-10-22 15:10:46 engels Exp $
   */
  class GearMgrImpl : public GearMgr {
	
    typedef std::map< std::string ,  GearParameters* > ParameterMap ;

  public: 

    // C'tor 
    GearMgrImpl() ;
    
    /// Destructor.
    virtual ~GearMgrImpl() ;
	
   /** The unique detector name - typically the model name used in the simulation program
    */
    virtual const std::string& getDetectorName() const  throw (UnknownParameterException, std::exception ) ;


    /** Get named parameters for key. This can be used to describe a subdetector that is not 
     *  yet forseen in the Gear API.
     * 
     *  @throws UnknownParameterException
     */
    virtual const GearParameters & getGearParameters(const std::string & key) const 
      throw (UnknownParameterException, std::exception )  ;

    /** Get the BField.
     */
    virtual const BField & getBField() const 
      throw (UnknownParameterException, std::exception ) ;

    
    /** Get the SiPlanes parameters.
     *
     *  @throws UnknownParameterException
     */
    virtual const SiPlanesParameters & getSiPlanesParameters() const
      throw (UnknownParameterException, std::exception )  ;

    
   /** Get the point properties object.
     * 
     *  @throws NotImplementedException
     */
    virtual const GearPointProperties & getPointProperties() const 
      throw (NotImplementedException, std::exception ) ;
    
    
    /** Get the distance properties object.
     * 
     *  @throws NotImplementedException
     */
    virtual const GearDistanceProperties & getDistanceProperties() const 
      throw (NotImplementedException, std::exception ) ;


    /** Keys of all GearParameters. 
     */ 
    virtual const std::vector<std::string>  & getGearParameterKeys() const ;


    virtual void setDetectorName(const std::string& name) { _detectorName = name ; }

    /** Set the GearParameters for the given key - overwrites any 
     *  existing entries.
     */
    virtual void setGearParameters( const std::string & key, GearParameters* parameters ) ;
  
     /** Set the SiPlanesParameters.
     */
    virtual void setSiPlanesParameters( SiPlanesParameters * siplanesParameters ) ;


    /** Set the point properties object 
     */
    virtual void  setPointProperties( GearPointProperties* pointProperties) ; 

    /** Set the distance properties object 
     */
    virtual void  setDistanceProperties( GearDistanceProperties* distanceProperties) ; 
    
    /** Set the b field object
     */
    virtual void setBField( BField* bField ) ;
    
    
    
    
  protected:
    
    ParameterMap _map ;
    SiPlanesParameters* _siplanesParameters ;
    GearPointProperties*  _pointProperties ;
    GearDistanceProperties*  _distanceProperties ;
    BField* _bField ;
    std::string _detectorName ;

    mutable StringVec _keys ;

  }; // class
} // namespace gear
#endif /* ifndef GEAR_GEARMGR_H */
