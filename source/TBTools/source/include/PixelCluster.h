#ifndef PixelCluster_H
#define PixelCluster_H


// Include LCIO header files 
#include "lcio.h"
#include <IMPL/TrackerDataImpl.h>


namespace depfet {
  
  //! PixelCluster
  /*! This class stores all information about reconstructed pixel clusters
   * 
   *  @Author Benjamin Schwenker, Universitaet Goettingen
   *  <mailto:benjamin.schwenker@phys.uni-goettingen.de>
   */ 
  
  class PixelCluster  {
  
  public:

    /** Default constructor */
    PixelCluster():
      m_sensorID(0), m_clsCharge(0), m_seedCharge(0),
      m_clsSize(0), m_uSize(0), m_vSize(0), m_uStart(0), m_vStart(0), m_clsType(0) 
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
    PixelCluster(unsigned short sensorID, float seedCharge, float clsCharge, unsigned short clsSize, 
               unsigned short uSize, unsigned short vSize, unsigned short uStart, unsigned short vStart, unsigned short clsType):
      m_sensorID(sensorID), m_clsCharge(clsCharge), m_seedCharge(seedCharge),  
      m_clsSize(clsSize), m_uSize(uSize), m_vSize(vSize), m_uStart(uStart), m_vStart(vStart), m_clsType(clsType)
    {}
    
    /** Constructor */
    PixelCluster(lcio::TrackerData * Digits, unsigned short sensorID = 0, unsigned short clsType=0) ;
    
    //! Destructor
    virtual ~PixelCluster() { /* NOOP */ ; }

    /** Get the sensor ID.
     * @return ID of the sensor.
     */
    unsigned short getSensorID() const { return m_sensorID; }

    /** Get collected charge.
     * @return charge collected in the cluster.
     */
    float getCharge() const { return m_clsCharge; }

    /** Get seed charge.
     * @return seed charge of the cluster.
     */
    float getSeedCharge() const { return m_seedCharge; }
    
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

    /** Get cluster type
     * @return the cluster type.
     */
    unsigned short getClusterType() const { return m_clsType; }
        
  protected:
      
    unsigned short m_sensorID; // SensorID
    float m_clsCharge;         // Deposited charge 
    float m_seedCharge;        // Cluster seed charge 
    unsigned short m_clsSize;  // Cluster size in pixels 
    unsigned short m_uSize;    // Cluster size in ucells
    unsigned short m_vSize;    // Cluster size in vcells  
    unsigned short m_uStart;   // Start ucell of the cluster 
    unsigned short m_vStart;   // Start vcell of the cluster 
    unsigned short m_clsType;  // Cluster type   
    
    
  };
 
}
#endif
