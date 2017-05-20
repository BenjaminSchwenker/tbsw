#ifndef PixelCluster_H
#define PixelCluster_H


// Include LCIO header files 
#include "lcio.h"
#include <IMPL/TrackerDataImpl.h>

#include <string>
#include <ostream>
#include <sstream>

namespace depfet {
  
  //! RawDigits
  /*! This class represents raw measurements from the pixel sensor
   * 
   *  @Author Benjamin Schwenker, Universitaet Goettingen
   *  <mailto:benjamin.schwenker@phys.uni-goettingen.de>
   */ 
  
  class RawDigit
  {
    public: 
     
     /** Default constructor */
     RawDigit(unsigned short iU=0, unsigned short iV=0, unsigned short charge=1) : m_cellIDU(iU), m_cellIDV(iV), m_charge(charge)
     {}
     
     /** Sorting RawDigits */
     bool operator < (const RawDigit& other) const
     {
       if ( m_cellIDV < other.m_cellIDV ) 
         return true;  
       else if (  other.m_cellIDV < m_cellIDV)
         return false;
       else 
	     return m_cellIDU < other.m_cellIDU;  
     }
     
     unsigned short m_cellIDU;
     unsigned short m_cellIDV;
     unsigned short m_charge;
  };


  //! PixelCluster
  /*! This class represents a sorted pixel cluster 
   * 
   *  @Author Benjamin Schwenker, Universitaet Goettingen
   *  <mailto:benjamin.schwenker@phys.uni-goettingen.de>
   */ 
  
  class PixelCluster  {
  
  public:

    /** Default constructor */
    PixelCluster():
      m_sensorID(0), m_clsCharge(0), m_seedCharge(0),
      m_clsSize(0), m_uSize(0), m_vSize(0), m_uStart(0), m_vStart(0)
    {}
    
    /** Constructor */
    PixelCluster(lcio::TrackerData * Digits, unsigned short sensorID = 0) ;
    
    //! Destructor
    virtual ~PixelCluster() { /* NOOP */ ; }
    
    /** Convert to string */
    operator std::string() const;
    
    /** Get the sensor ID.
     * @return ID of the sensor.
     */
    unsigned short getSensorID() const { return m_sensorID; }
    
    /** Get collected charge.
     * @return charge collected in the cluster.
     */
    unsigned short getCharge() const { return m_clsCharge; }
    
    /** Get seed charge.
     * @return seed charge of the cluster.
     */
    unsigned short getSeedCharge() const { return m_seedCharge; }
    
    /** Get cluster size.
     * @return number of pixels contributing to the cluster.
     */
    unsigned short getSize() const { return m_clsSize; }
    
    /** Get cluster size in u direction.
     * @return number of uCells contributing to the cluster.
     */
    unsigned short getUSize() const { return m_uSize; }
    
    /** Get cluster size in v direction.
     * @return number of vCells contributing to the cluster.
     */
    unsigned short getVSize() const { return m_vSize; }
    
    /** Get cluster start cell in u direction.
     * @return first uCell contributing to the cluster.
     */
    unsigned short getUStart() const { return m_uStart; }
    
    /** Get cluster start cell in v direction.
     * @return first vCell contributing to the cluster.
     */
    unsigned short getVStart() const { return m_vStart; }
    
    /** Get cluster label string  
     * @return cluster label used for cluster calibration 
     */
    std::string getLabel(int scale = 1) const;

    /** Get cluster label string 
     * @return cluster label used for cluster calibration 
     */
    std::string getLabel(const std::vector<int>& jumps) const; 
     
    /** Get cluster type string. 
     * @return type string
     */
    std::string getType() const;
    
    /** Get sorted vector of raw digits 
     */
    const std::vector<RawDigit>& getRawDigits() const { return m_sortedDigits; }
    
  protected:
     
    // SensorID 
    unsigned short m_sensorID;  
    // Cluster charge charge in ADC units 
    unsigned short m_clsCharge;  
    // Seed charge in ADC units 
    unsigned short m_seedCharge; 
    // Cluster size in digits 
    unsigned short m_clsSize;   
    // Cluster size in ucells 
    unsigned short m_uSize;  
    // Cluster size in vcells     
    unsigned short m_vSize; 
    // Start ucell of the cluster      
    unsigned short m_uStart;   
    // Start vcell of the cluster   
    unsigned short m_vStart;     
    // Sorted vector of raw digits 
    std::vector<RawDigit> m_sortedDigits;
  };
}
#endif
