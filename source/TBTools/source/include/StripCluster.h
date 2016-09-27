#ifndef StripCluster_H
#define StripCluster_H


// Include LCIO header files 
#include "lcio.h"
#include <IMPL/TrackerDataImpl.h>


namespace depfet {
  
  //! StripCluster
  /*! This class stores all information about reconstructed strip clusters
   * 
   *  @Author Benjamin Schwenker, Universitaet Goettingen
   *  <mailto:benjamin.schwenker@phys.uni-goettingen.de>
   */ 
  
  class StripCluster  {
  
  public:

    /** Default constructor */
    StripCluster():
      m_sensorID(0), m_uClsCharge(0), m_uSeedCharge(0), m_vClsCharge(0), m_vSeedCharge(0),
      m_uSize(0), m_vSize(0), m_uStart(0), m_vStart(0), m_clsType(0) 
    {}
    
    /**  Constructor
     * @param sensorID The DAQ ID 
     * @param clsCharge The cluster charge
     * @param seedCharge The charge of the cluster seed
     * @param clsSize Size of the cluster in pixels
     * @param uSize Number of uCells contributing to the cluster
     * @param vSize Number of vCells contributing to the cluster
     * @param uStart First uCell contributing to the cluster
     * @param vStart First vCell contributing to the cluster
     */
    StripCluster(unsigned short sensorID, float uSeedCharge, float uClsCharge,  float vSeedCharge, float vClsCharge,  
               unsigned short uSize, unsigned short vSize, unsigned short uStart, unsigned short vStart, unsigned short clsType):
      m_sensorID(sensorID), m_uClsCharge(uClsCharge), m_uSeedCharge(uSeedCharge), m_vClsCharge(uClsCharge), m_vSeedCharge(uSeedCharge),  
      m_uSize(uSize), m_vSize(vSize), m_uStart(uStart), m_vStart(vStart), m_clsType(clsType)
    {}
    
    /** Constructor */
    StripCluster(lcio::LCObjectVec& DigitVec, unsigned short sensorID = 0, unsigned short clsType=0) ;
    
    //! Destructor
    virtual ~StripCluster() { /* NOOP */ ; }
    
    /** Get the sensor ID.
     * @return ID of the sensor.
     */
    unsigned short getSensorID() const { return m_sensorID; }
    
    /** Get collected charge.
     * @return charge collected in the cluster.
     */
    float getUCharge() const { return m_uClsCharge; }

    /** Get collected charge.
     * @return charge collected in the cluster.
     */
    float getVCharge() const { return m_vClsCharge; }

    /** Get seed charge.
     * @return seed charge of the cluster.
     */
    float getUSeedCharge() const { return m_uSeedCharge; }

    /** Get seed charge.
     * @return seed charge of the cluster.
     */
    float getVSeedCharge() const { return m_vSeedCharge; }
    
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

    /** Get cluster type
     * @return the cluster type.
     */
    unsigned short getClusterType() const { return m_clsType; }
        
  protected:
      
    unsigned short m_sensorID; // SensorID
    float m_uClsCharge;        // Deposited charge 
    float m_uSeedCharge;        // Cluster seed charge 
    float m_vClsCharge;        // Deposited charge 
    float m_vSeedCharge;        // Cluster seed charge  
    unsigned short m_uSize;    // Cluster size in ucells
    unsigned short m_vSize;    // Cluster size in vcells  
    unsigned short m_uStart;   // Start ucell of the cluster 
    unsigned short m_vStart;   // Start vcell of the cluster 
    unsigned short m_clsType;  // Cluster type   
    
    
  };
 
}
#endif
