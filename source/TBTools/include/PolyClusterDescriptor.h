#ifndef PolyClusterDescriptor_H
#define PolyClusterDescriptor_H

#include "Det.h"
#include "PixelCluster.h"


#include <string>


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
     PolyRawDigit(float u=0, float v=0, float charge=1, int type=0) : m_cellPosU(u), m_cellPosV(v), m_charge(charge), m_pixelType(type)
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
     float m_charge;
     int m_pixelType;
  };


  //! PolyClusterDescriptor
  /*! This class represents a sorted poly pixel cluster 
   * 
   *  @Author Benjamin Schwenker, Universitaet Goettingen
   *  <mailto:benjamin.schwenker@phys.uni-goettingen.de>
   */ 
  
  class PolyClusterDescriptor  {
  
  public:

    /** Default constructor */
    PolyClusterDescriptor():
      m_originU(0), m_originV(0), m_uStart(0), m_vStart(0), m_lowerRight(0), m_upperLeft(0), m_upperRight(0),  m_lowerLeft(0)
    {}
    
    /** Constructor */
    PolyClusterDescriptor(const PixelCluster& Cluster, const Det& Sensor) ;

    /** Get origin coordinate
     */
    float getOriginU() const { return  m_originU; }

    /** Get origin coordinate
     */
    float getOriginV() const { return  m_originV; }
    
    /** Get cluster shape string 
     * @return shape string containing all features of cluster used for position reconstruction 
     */
    std::string getShape(int vCellPeriod = 1, int uCellPeriod = 1, int etaBin=0) const; 
     
    /** Get cluster type string. 
     * @return type string, same as shape but without eta information
     */
    std::string getType(int vCellPeriod = 1, int uCellPeriod = 1) const;

    /** Get cluster shape string
     * @return shape string containing all features of cluster used for position reconstruction
     */
    std::string getEtaBinString(int etaBin=0) const;
    
    /** Get cluster eta
    */
    double computeEta(double thetaU, double thetaV, int etaIndex=0) const; 
    
    /** Get cluster eta bin  
    */
    static int computeEtaBin(double eta, const std::vector<double>& etaBinEdges); 

  private:
    
    /** Get indices of head pixels in cluster
    */
    std::vector<int> getHeadPixels(double thetaU, double thetaV) const; 
    
    /** Get indices of tail pixels in cluster
    */
    std::vector<int> getTailPixels(double thetaU, double thetaV) const;
    
    // Origin of cluster coordinates
    float m_originU;   
    // Origin of cluster coordinates
    float m_originV;  
    // Start ucell of the cluster      
    unsigned short m_uStart;   
    // Start vcell of the cluster   
    unsigned short m_vStart;   
    // Index of digit at lower right corner of cluster
    int m_lowerRight;
    // Index of digit at upper left corner of cluster
    int m_upperLeft; 
    // Index of digit at upper right corner of cluster
    int m_upperRight;
    // Index of digit at lower left corner of cluster
    int m_lowerLeft;
      
    // Sorted vector of raw digits 
    std::vector<PolyRawDigit> m_sortedDigits;
  };
}
#endif
