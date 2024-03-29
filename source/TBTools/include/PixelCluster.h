#ifndef PixelCluster_H
#define PixelCluster_H

#include "Det.h"

// Include LCIO header files 
#include "lcio.h"
#include <IMPL/TrackerDataImpl.h>

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
     RawDigit(unsigned short iU=0, unsigned short iV=0, unsigned short charge=0, unsigned short time=0) : m_cellIDU(iU), m_cellIDV(iV), m_charge(charge), m_time(time)
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
     unsigned short m_time;
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
      m_clsSize(0), m_uSize(0), m_vSize(0), m_uStart(0), m_vStart(0), m_uSeed(0), m_vSeed(0), m_minTime(0), m_maxTime(0)
    {}
    
    /** Constructor */
    PixelCluster(lcio::TrackerData * Digits, unsigned short sensorID = 0) ;
    
    //! Destructor
    ~PixelCluster() { /* NOOP */ ; }
    
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
    
    /** Get seed pixel cell in u direction.
     * @return uCell of seed pixel.
     */
    short getUSeed() const { return m_uSeed; }
    
    /** Get seed pixel cell in v direction.
     * @return vCell of seed pixel.
     */
    short getVSeed() const { return m_vSeed; }

    /** Get min time.
     * @return minimum time over digit times.
     */
    short getMinTime() const { return m_minTime; }

    /** Get max time.
     * @return maximum time over digit times.
     */
    short getMaxTime() const { return m_maxTime; }
    
    /** Get sorted vector of raw digits 
     */
    const std::vector<RawDigit>& getRawDigits() const { return m_sortedDigits; }
    
    /** Compute center of gravity hit position and covariance matrix
    */
    void getCenterOfGravity(const Det& Sensor, double& u, double& v, double& cov_u, double& cov_v, double& cov_uv) const;

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
    // Seed pixel uCell address   
    short m_uSeed;
    // Seed pixel vCell address  
    short m_vSeed;
    // Smallest time of hit   
    unsigned short m_minTime;
    // Largest time of hit   
    unsigned short m_maxTime;

    // Sorted vector of raw digits 
    std::vector<RawDigit> m_sortedDigits;
  };
}
#endif
