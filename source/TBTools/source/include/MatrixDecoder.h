#ifndef MATRIXDECODER_H
#define MATRIXDECODER_H 1

// Include LCIO
#include <UTIL/CellIDDecoder.h>

// Include Marlin
#include <marlin/Global.h>
#include <streamlog/streamlog.h>

// system includes <>
#include <string>
#include <iostream>
#include <iomanip>
#include <cstdlib> 

namespace depfet {
  

  //! Pixel matrix decoder.
  /*! This helper class is used to find the x and y coordinates from a
   *  pixel index and vice versa. 
   *
   *  In the Tracker(Raw)Data object all the pixels are lined up into
   *  a 1D vector and the 2 coordinates can be obtained using the this
   *  decoder class, once the number of pixels along x and y are
   *  properly set. 
   *  
   *
   */ 
  class MatrixDecoder {
    
  public:
    //! Constructor with two integer numbers
    /*! This is the basic constructor where the number of pixels along
     *  the directions is given directly as parameter of the constructor
     *
     *  @param xNoOfPixel The number of pixels along X
     *  @param yNoOfPixel The number of pixels along Y
     *
     */
    MatrixDecoder(int xNoOfPixel, int yNoOfPixel);

    //! Constructor with four integer numbers
    /*! This is the complete version of the previous constructor
     *  representing the most general case. In the case the full
     *  matrix is divided into submatrix and each of them is
     *  considered a separated item being stored into a different
     *  Trakere(Raw)Data, so the first pixel will have coordinates
     *  different from (0,0). Those two numbers should be then passed
     *  to the MatrixDecoder in order to properly work. 
     *
     *  Remember that according to the general numbering convention,
     *  pixels are always counted starting from zero and that the
     *  pixel in the frame of reference origin has coordinates (0, 0).
     * 
     *  @param xNoOfPixel The number of pixels along X
     *  @param yNoOfPixel The number of pixels along Y
     *  @param xMin       X coordinate of the first pixel
     *  @param yMin       Y coordinate of the first pixel
     */
    MatrixDecoder(int xNoOfPixel, int yNoOfPixel, int xMin, int yMin);

    //! Constructor with a cell id decoder
    /*! The MatrixDecoder is used to decode data included into a
     *  Tracker(Raw)Data, so it is very likely that a CellIDDecoder
     *  for the corresponding collection already exists and can be
     *  used to get these information.
     *
     *  The current implementation is using a template for both the
     *  CellIDDecoder and the rawdata matrix. In case the decoder does
     *  not contain the required fields, the lcio::exception is caught
     *  and the program execution is stopped.
     *
     *  @param decoder A reference to the cell id decoder 
     *  @param rawData A pointer to the tracker (raw) data needed
     *  for the decoding
     *
     *
     */
    template <class T>
    MatrixDecoder(UTIL::CellIDDecoder<T >& decoder, T * rawData) {
  
      _xNoOfPixel = decoder(rawData)["xMax"] - decoder(rawData)["xMin"] + 1;
      _yNoOfPixel = decoder(rawData)["yMax"] - decoder(rawData)["yMin"] + 1;
      _xMin       = decoder(rawData)["xMin"];
      _yMin       = decoder(rawData)["yMin"];

      if ( _xNoOfPixel <= 0 ) {
         streamlog_out (ERROR3) << "MatrixDecoder: xNoOfPixel has to be positive" << std::endl;
         exit(-1); 
      }
      if ( _yNoOfPixel <= 0 ) {
        streamlog_out (ERROR3) << "MatrixDecoder: yNoOfPixel has to be positive" << std::endl;
        exit(-1); 
      }
    }

    //! Finds the index position having the two coordinates
    /*! Since the relation between the position index and the x, y
     *  coordinate is 1 to 1, also the inverse function can be defined
     *
     *  @param x The x coordinate of the pixel you want to know the index
     *  @param y The y coordinate of the pixel you want to know the index
     *  @return the position index corresponding to x and y
     */
    int getIndexFromXY(int x, int y) const;

    //! Finds the Y coordinate having the index
    /*! Pixels in the Tracker(Raw)Data vector are pushed back
     *  row-wise. This means that to obtain the Y coordinate is just
     *  enough to make the integer division between @c index and the
     *  number of pixels along x (@c _xNoOfPixel). To this number the
     *  eventual offset along the y axis (@c _yMin) is added.
     * 
     *  @param index The position for which the Y coordinate is
     *  requested
     *  @return The Y coordinate corresponding to the @c index
     */
    int getYFromIndex(int index) const;

    //! Finds the X coordinate having the index
    /*! Pixels in the Tracker(Raw)Data vector are pushed back
     *  row-wise. To obtain the X coordinate, first we need to count
     *  how many full row index corresponds to. This is obtained
     *  using the % binary operator between index and _xNoOfPixel. The
     *  eventual x offset is also added.
     */
    int getXFromIndex(int index) const;
     
    //! Finds both coordinates having the index
    /*! Utility to get both coordinates having the index. 
     *  
     *  @param index The pixel index
     *  @param x A reference to the corresponding X coordinate
     *  @param y A reference to the corresponding Y coordinate
     */
    void getXYFromIndex(int index, int& x, int& y) const ;

    //! Returns the minimum value of X
    /*! 
     *  @return The minimum value of X
     */ 
    inline int getMinX() const { return _xMin ; }

    //! Returns the minimum value of Y
    /*! 
     *  @return The minimum value of Y
     */ 
    inline int getMinY() const { return _yMin ; }

   //! Returns the maximum value of X
    /*! 
     *  @return The minimum value of X
     */ 
    inline int getMaxX() const { return _xMin + _xNoOfPixel - 1 ; }

    //! Returns the minimum value of Y
    /*! 
     *  @return The minimum value of Y
     */ 
    inline int getMaxY() const { return _yMin + _yNoOfPixel - 1; }

    

  private:
    
    //! The number of pixels along x
    int _xNoOfPixel;

    //! The number of pixels along y
    int _yNoOfPixel;

    //! The origin of the x axis
    /*! Especially in the case of matrix divided into sub-matrix, the
     *  minimum value of the x axis can be different from 0
     */ 
    int _xMin;

    //! The origin of the y axis
    /*! Especially in the case of matrix divided into sub-matrix, the
     *  minimum value of the y axis can be different from 0
     */ 
    int _yMin;


  };
  
}
#endif
