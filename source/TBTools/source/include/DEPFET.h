#ifndef DEPFET_NAMESPACE_H
#define DEPFET_NAMESPACE_H

//! The depfet namespace.
/*! This namespace is used in order not to pollute neither the lcio,
 *  nor the Marlin, not the standard namespaces. It contains all
 *  classes defined by the DEPFET collaboration in order to
 *  develop both the analysis/reconstruction software.
 *
 *  @author Benjamin Schwenker, Göttingen University <mailto:benjamin.schwenker@phys.uni-goettingen.de>
 */

namespace depfet {}
#endif

#ifndef DEPFET_H
#define DEPFET_H

// system includes <>
#include <iostream>
#include <string>
#include <vector>

namespace depfet
{

  //! Global constants used in the Eutelescope package
  /*!
   * This class has only static data members used only to define global
   * constant to be used within the DEPFET package. Please add here
   * whatever constant you want to use.  A typical useful of this class
   * is to define name of collection to be retrieved/saved from/to
   * files.
   *
   */
  
  class DEPFET
  {
  
  public:
    //! Default destructor.
    /*! This is the default destructor, but it is actually a NO-OP
     *  since there is nothing to be destroyed.
     */
    virtual ~ DEPFET ()  {;  }
    
  public:
    
    
    // Encoding strings
    
    //! Default encoding for full data 
    /*! 
     *  "sensorID:17,uMax:12,vMax:12"
     */
    static const char * MATRIXDEFAULTENCODING;
         
    //! Default encoding for zero suppress data
    /*!
     *  "sensorID:5,sparsePixelType:5"
     */
    static const char * ZSDATADEFAULTENCODING;
    
    //! Cluster default encoding
    /*! 
     *  "sensorID:5,clusterID:8,sparsePixelType:5,quality:5"
     */
    static const char * ZSCLUSTERDEFAULTENCODING;
    
  }; // End class
  
  

  //! Cluster quality enum
  /*! This enum can be attached to a LCIO class describing a cluster
   *  or it can be inserted into the CellID describing the cluste
   */
   
  /*
  enum ClusterQuality {
    kGoodCluster       = 0,
    kIncompleteCluster = 1L << 0,
    kBorderCluster     = 1L << 1,
    kMergedCluster     = 1L << 2
  };
  */
  
  //! Cluster quality bit-wise AND operator
  /*! This is a convenience operator used to identify the reason of a
   *  non good quality clusters. Bad quality clusters may be so for more
   *  than one reason simultaneously. This operator is used in the
   *  identification of such reasons.
   *
   *  @param a A cluster quality value
   *  @param b Another cluster quality value
   *  @return the bit wise and among @a a and @a b
   */
  //ClusterQuality operator&(ClusterQuality a, ClusterQuality b);
  
  //! Cluster quality bit-wise OR operator
  /*! This is a crucial operator for ClusterQuality since, during the
   *  cluster search processor, a cluster maybe flagged with one or
   *  more than one "bad" qualities. For this reason, using this
   *  operator can allow to flag the same cluster with more than one
   *  bad qualities.
   *
   *  @param a A cluster quality value
   *  @param b Another cluster quality value
   *  @return the bit wise or among @a a and @a b
   */
  //ClusterQuality operator|(ClusterQuality a, ClusterQuality b);
  
  //! Cluster quality operator <<
  /*! This operator can be used to stream out the value of a cluster
   *  quality enum. Both the numerical and the textual values are shown.
   *
   *  @param os The input output stream
   *  @param quality The variable to the be stream out
   *  @return The input output stream
   */
  //std::ostream& operator<<(std::ostream& os, const ClusterQuality & quality);  
  
  //! Pixel quality enum
  /*! This enum can be attached to a LCIO class describing a pixel
   *  or it can be inserted into the CellID describing the pixel
   *  collection. It is a five bit flag, that can be used to
   *  discriminate among different pixel qualities. This is because
   *  not all pixels can be considered to be at the same quality levels. 
   *
   *  Here a description of all allowed value of cluster quality and
   *  their meaning:
   *
   *  \li <b>kGoodPixel</b>: this flag is used to identify pixels
   *  having no problem at all.
   *
   *  \li <b>kDeadPixel</b>: this flag is used to identify
   *  pixels with very low noise.
   *
   *  \li <b>kNoisyPixel</b>: this flag is used to identify 
   *  pixel with very high noise.
   *
   *  \li <b>kHotPixel</b>: this flag is used to label pixels 
   *  having a too high mean firing frequency 
   * 
   *  \li <b>kBrigthPixel</b>: this flag is used to label pixels 
   *  having strange pedestal.  
   *  
   *  \li <b>kBadChannel</b>: this flag is used to label pixels 
   *  connected to bad DCDB channel  
   * 
   *  \li <b>kBadPixel</b>: this flag is used to label pixels 
   *  considered bad for any other reason 
   * 
   *  There are still two "not assigned" bits that can be used in the
   *  future to mark other different kind of bad quality pixels.
   *
   */

  enum PixelQuality {
    kGoodPixel         = 0,
    kDeadPixel         = 1L << 0,
    kNoisyPixel        = 1L << 1,
    kHotPixel          = 1L << 2,
    kBrigthPixel       = 1L << 3,
    kBadChannel        = 1L << 4,
    kBadPixel          = 1L << 5
    
  };
  
  //! Pixel quality bit-wise AND operator
  /*! This is a convenience operator used to identify the reason of a
   *  non good quality pixel. Bad quality pixels may be so for more
   *  than one reason simultaneously. This operator is used in the
   *  identification of such reasons.
   *
   *  @param a A pixel quality value
   *  @param b Another pixel quality value
   *  @return the bit wise and among @a a and @a b
   */
  PixelQuality operator&(PixelQuality a, PixelQuality b);

  //! Pixel quality bit-wise OR operator
  /*! This is a crucial operator for PixelQuality since, during the
   *  quality tests, a pixel maybe flagged with one or
   *  more than one "bad" qualities. For this reason, using this
   *  operator can allow to flag the same pixel with more than one
   *  bad qualities.
   *
   *  @param a A pixel quality value
   *  @param b Another pixel quality value
   *  @return the bit wise or among @a a and @a b
   */
  PixelQuality operator|(PixelQuality a, PixelQuality b);

  //! Pixel quality operator <<
  /*! This operator can be used to stream out the value of a pixel
   *  quality enum. Both the numerical and the textual values are shown.
   *
   *  @param os The input output stream
   *  @param quality The variable to the be stream out
   *  @return The input output stream
   */
  std::ostream& operator<<(std::ostream& os, const PixelQuality & quality);

  //! Sparse pixel type enum
  /*! This enumerator is used to define the sparsified pixel type.
   *
   */
  enum SparsePixelType {
    kBaseSparsePixel   = 0,
    kSimpleSparsePixel = 1,
    // add here your implementation
    kUnknownPixelType       = 31
  };

}

#endif
