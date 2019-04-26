#ifndef PolyClusterDescriptor_H
#define PolyClusterDescriptor_H

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
      m_originU(0), m_originV(0)
    {}
    
    /** Constructor */
    PolyClusterDescriptor(lcio::TrackerData * Digits, const Det& Sensor) ;

    /** Get origin coordinate
     */
    float getOriginU() const { return  m_originU; }

    /** Get origin coordinate
     */
    float getOriginV() const { return  m_originV; }
    
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
    
  protected:
    
    // Origin of cluster coordinates
    float m_originU;   
    // Origin of cluster coordinates
    float m_originV;   
    
    // Sorted vector of raw digits 
    std::vector<PolyRawDigit> m_sortedDigits;
  };
}
#endif
