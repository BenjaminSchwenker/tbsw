#ifndef PolyPixelCluster_H
#define PolyPixelCluster_H

#include "Det.h"

// Include LCIO header files 
#include "lcio.h"
#include <IMPL/TrackerDataImpl.h>

#include <string>
#include <ostream>
#include <sstream>

namespace depfet {
  
  //! PolyRawDigit
  /*! This class represents raw measurements from the pixel sensor
   *  
   *  A PolyRawDigit stores the pixel charge as well as the position of 
   *  of the pixel center in sensor uv coordinates and the pixelType. 
   *  
   *  @Author Benjamin Schwenker, Universitaet Goettingen
   *  <mailto:benjamin.schwenker@phys.uni-goettingen.de>
   */ 
  
  class PolyRawDigit
  {
    public: 
     
     /** Default constructor */
     PolyRawDigit(float u=0, float v=0, unsigned short charge=1, int type=0) : m_cellPosU(u), m_cellPosV(v), m_charge(charge), m_pixelType(type)
     {}
     
     /** Sorting PolyRawDigits */
     bool operator < (const PolyRawDigit& other) const
     {
       if ( m_cellPosV < other.m_cellPosV ) 
         return true;  
       else if (  other.m_cellPosV < m_cellPosV)
         return false;
       else 
	     return m_cellPosU < other.m_cellPosU;  
     }
     
     float m_cellPosU;
     float m_cellPosV;
     unsigned short m_charge;
     int m_pixelType;
  };


  //! PolyPixelCluster
  /*! This class represents a sorted poly pixel cluster 
   * 
   *  @Author Benjamin Schwenker, Universitaet Goettingen
   *  <mailto:benjamin.schwenker@phys.uni-goettingen.de>
   */ 
  
  class PolyPixelCluster  {
  
  public:

    /** Default constructor */
    PolyPixelCluster():
      m_sensorID(0), m_clsCharge(0), m_seedCharge(0),
      m_clsSize(0), m_uSize(0), m_vSize(0), m_uStart(0), m_vStart(0)
    {}
    
    /** Constructor */
    PolyPixelCluster(lcio::TrackerData * Digits, const Det& Sensor) ;
    
    //! Destructor
    virtual ~PolyPixelCluster() { /* NOOP */ ; }
    
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
    
    /** Get cluster shape string 
     * @return shape string containing all features of cluster used for position reconstruction 
     */
    std::string getShape(int pixeltype = 1, int vCellPeriod = 1, int uCellPeriod = 1, int etaBin=0) const; 
     
    /** Get cluster type string. 
     * @return type string, same as shape but without eta information
     */
    std::string getType(int pixeltype = 1, int vCellPeriod = 1, int uCellPeriod = 1) const;

    /** Get cluster shape string
     * @return shape string containing all features of cluster used for position reconstruction
     */
    std::string getEtaBinString(int etaBin=0) const;

    /** Get sorted vector of raw digits 
     */
    const std::vector<PolyRawDigit>& getRawDigits() const { return m_sortedDigits; }

    /** Get cluster eta
    */
    double computeEta(double thetaU, double thetaV) const; 
    
    /** Get cluster eta bin  
    */
    int computeEtaBin(double eta, const std::vector<double>& etaBinEdges) const; 

    /** Get index of head pixel in cluster
    */
    int getHeadPixelIndex(double thetaU, double thetaV) const; 
    
    /** Get index of tail pixel in cluster
    */
    int getTailPixelIndex(double thetaU, double thetaV) const;
    
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
    
    
    // Sorted vector of raw digits 
    std::vector<PolyRawDigit> m_sortedDigits;
  };
}
#endif
