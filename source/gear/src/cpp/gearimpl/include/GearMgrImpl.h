#ifndef GEAR_GearMgrImpl_H
#define GEAR_GearMgrImpl_H 1

#include <string>

#include "gear/GearMgr.h"
#include "gear/GearParameters.h"
#include "gear/TPCParameters.h"


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

    /** Get the TPCParameters.
     */
    virtual const TPCParameters & getTPCParameters() const
      throw (UnknownParameterException, std::exception ) ;
    
     /** Get the Ecal barrel parameters.
     *
     *  @throws UnknownParameterException
     */
    virtual const CalorimeterParameters & getEcalBarrelParameters() const 
	throw (UnknownParameterException, std::exception ) ;

    /** Get the Ecal endcap parameters.
     *
     *  @throws UnknownParameterException
     */
    virtual const CalorimeterParameters & getEcalEndcapParameters() const 
	throw (UnknownParameterException, std::exception )  ;

    /** Get the Ecal plug parameters.
     *
     *  @throws UnknownParameterException
     */
    virtual const CalorimeterParameters & getEcalPlugParameters() const 
	throw (UnknownParameterException, std::exception )  ;

      /** Get the Yoke barrel parameters.
     *
     *  @throws UnknownParameterException
     */
    virtual const CalorimeterParameters & getYokeBarrelParameters() const 
	throw (UnknownParameterException, std::exception ) ;

    /** Get the Yoke endcap parameters.
     *
     *  @throws UnknownParameterException
     */
    virtual const CalorimeterParameters & getYokeEndcapParameters() const 
	throw (UnknownParameterException, std::exception )  ;

    /** Get the Yoke plug parameters.
     *
     *  @throws UnknownParameterException
     */
    virtual const CalorimeterParameters & getYokePlugParameters() const 
	throw (UnknownParameterException, std::exception )  ;

   /** Get the Hcal barrel parameters.
     *
     *  @throws UnknownParameterException
     */
    virtual const CalorimeterParameters & getHcalBarrelParameters() const 
	throw (UnknownParameterException, std::exception )  ;


    /** Get the Hcal endcap parameters.
     *
     *  @throws UnknownParameterException
     */
    virtual const CalorimeterParameters & getHcalEndcapParameters() const 
	throw (UnknownParameterException, std::exception )  ;


    /** Get the Hcal ring parameters.
     *
     *  @throws UnknownParameterException
     */
    virtual const CalorimeterParameters & getHcalRingParameters() const 
	throw (UnknownParameterException, std::exception )  ;


    /** Get the Lcal parameters.
     *
     *  @throws UnknownParameterException
     */
    virtual const CalorimeterParameters & getLcalParameters() const 
      throw (UnknownParameterException, std::exception ) ;
    
    /** Get the LHcal parameters.
     *
     *  @throws UnknownParameterException
     */
    virtual const CalorimeterParameters & getLHcalParameters() const 
      throw (UnknownParameterException, std::exception ) ;
 
    /** Get the BeamCal parameters.
     *
     *  @throws UnknownParameterException
     */
    virtual const CalorimeterParameters & getBeamCalParameters() const 
      throw (UnknownParameterException, std::exception ) ;
 
    /** Get the VXD parameters.
     *
     *  @throws UnknownParameterException
     */
    virtual const VXDParameters & getVXDParameters() const
      throw (UnknownParameterException, std::exception )  ;

    /** Get the SiPlanes parameters.
     *
     *  @throws UnknownParameterException
     */
    virtual const SiPlanesParameters & getSiPlanesParameters() const
      throw (UnknownParameterException, std::exception )  ;

    /** Get the TBSi parameters. 
      * 
      *  @throws UnknownParameterException 
      */ 
     virtual const TBSiParameters & getTBSiParameters() const 
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
    
    /** Set the TPCParameters.
     */
    virtual void setTPCParameters( TPCParameters* tpcParameters ) ;


    /** Set the EcalBarrelParameters.
     */
    virtual void setEcalBarrelParameters( CalorimeterParameters* ecalBarrelParameters ) ;

    /** Set the EcalEndcapParameters.
     */
    virtual void setEcalEndcapParameters( CalorimeterParameters* ecalEndcapParameters ) ;
    
    /** Set the EcalPlugParameters.
     */
    virtual void setEcalPlugParameters( CalorimeterParameters* ecalPlugParameters ) ;
    
    /** Set the YokeBarrelParameters.
     */
    virtual void setYokeBarrelParameters( CalorimeterParameters* yokeBarrelParameters ) ;

    /** Set the YokeEndcapParameters.
     */
    virtual void setYokeEndcapParameters( CalorimeterParameters* yokeEndcapParameters ) ;
    
    /** Set the YokePlugParameters.
     */
    virtual void setYokePlugParameters( CalorimeterParameters* yokePlugParameters ) ;
 
	
	/** Set the HcalBarrelParameters.
     */
    virtual void setHcalBarrelParameters( CalorimeterParameters* hcalBarrelParameters ) ;

    /** Set the HcalEndcapParameters.
     */
    virtual void setHcalEndcapParameters( CalorimeterParameters* hcalEndcapParameters ) ;

    /** Set the HcalRingParameters.
     */
    virtual void setHcalRingParameters( CalorimeterParameters* hcalRingParameters ) ;

    /** Set the LcalParameters.
     */
    virtual void setLcalParameters(CalorimeterParameters * lcalParameters) ;

    /** Set the LHcalParameters.
     */
    virtual void setLHcalParameters(CalorimeterParameters * lhcalParameters) ;

    /** Set the BeamCalParameters.
     */
    virtual void setBeamCalParameters(CalorimeterParameters * beamcalParameters) ;

     /** Set the VXDParameters.
     */
    virtual void setVXDParameters( VXDParameters * vxdParameters ) ;


     /** Set the SiPlanesParameters.
     */
    virtual void setSiPlanesParameters( SiPlanesParameters * siplanesParameters ) ;

    /** Set TBSiParameters. 
      */
    virtual void setTBSiParameters( TBSiParameters * TBSiParameters ) ;

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
    TPCParameters* _tpcParameters ;
    CalorimeterParameters* _ecalBarrelParameters ;
    CalorimeterParameters* _ecalEndcapParameters ;
    CalorimeterParameters* _ecalPlugParameters ;
    CalorimeterParameters* _yokeBarrelParameters ;
    CalorimeterParameters* _yokeEndcapParameters ;
    CalorimeterParameters* _yokePlugParameters ;
    CalorimeterParameters* _hcalBarrelParameters ;
    CalorimeterParameters* _hcalEndcapParameters ;
    CalorimeterParameters* _hcalRingParameters ;
    CalorimeterParameters* _lcalParameters ;
    CalorimeterParameters* _lhcalParameters ;
    CalorimeterParameters* _beamcalParameters ;
    VXDParameters* _vxdParameters ;
    SiPlanesParameters* _siplanesParameters ;
    TBSiParameters* _tbsiParameters ;
    GearPointProperties*  _pointProperties ;
    GearDistanceProperties*  _distanceProperties ;
    BField* _bField ;
    std::string _detectorName ;

    mutable StringVec _keys ;

  }; // class
} // namespace gear
#endif /* ifndef GEAR_GEARMGR_H */
