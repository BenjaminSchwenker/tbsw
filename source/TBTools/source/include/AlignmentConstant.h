#ifndef ALIGNMENTCONSTAN_H
#define ALIGNMENTCONSTAN_H

#define ALIGN_CONST_MAX_SIZE 20

// lcio includes <.h>
#include <lcio.h>
#include <IMPL/LCGenericObjectImpl.h>

// system includes <>
#include <string>

namespace depfet {

  //! Alignment constant for the DEPFET TB package
  /*! The aim of this class is to store into a LCGenericObject the
   *  alignment contants. 
   */ 
  class AlignmentConstant : public IMPL::LCGenericObjectImpl {

  public:
    //! Default constructor
    /*! This is the default constructor of the alignment constants
     *  
     */
    AlignmentConstant(); 

    //! Constructor with all the needed parameters
    /*! This constructor is the recommended one because it set all the
     *  needed parameters in one go. 
     *
     *  @param sensorID The sensor number as from the GEAR file
     *  @param xOff Offset of the sensor along the x direction
     *  @param yOff Offset of the sensor along the y direction
     *  @param zOff Offset of the sensor along the z direction
     *  @param alpha Angle of the sensor along x
     *  @param beta Angle of the sensor along y
     *  @param gamma Angle of the sensor along z
     *  @param alphaErr Error on the angle of the sensor along x
     *  @param betaErr Error on the angle of the sensor along y
     *  @param gammaErr Error on the angle of the sensor along z
     *  @param xOffErr Error on the offset of the sensor along the x direction
     *  @param yOffErr Error on the offset of the sensor along the y direction
     *  @param zOffErr Error on the offset of the sensor along the z direction
     *
     */ 
    AlignmentConstant( int sensorID, 
			    double xOff,   double yOff,   double zOff,
			    double alpha, double beta, double gamma,
			    double xOffErr,   double yOffErr,   double zOffErr,
			    double alphaErr, double betaErr, double gammaErr );
      
    //! Default destructor
    virtual ~AlignmentConstant() { /* NO-OP */ ; }

    //! Set the sensor id
    void setSensorID( int id ) ;

    //! Set the x offset in mm 
    void setXOffset( double off ) ;

    //! Set the y offset in mm 
    void setYOffset( double off ) ;

    //! Set the z offset in mm 
    void setZOffset( double off ) ;

    //! Set the angle around x
    void setAlpha( double theta );

    //! Set the angle around y
    void setBeta( double theta );

    //! Set the angle around z
    void setGamma( double theta );

    //! Set the error of the offset along x
    void setXOffsetError( double err ) ;

    //! Set the error of the offset along y
    void setYOffsetError( double err ) ;

    //! Set the error of the offset along z
    void setZOffsetError( double err ) ;

    //! Set the error of the angle around x
    void setAlphaError( double err ) ;

    //! Set the error of the angle around y
    void setBetaError( double err ) ;

    //! Set the error of the angle around z
    void setGammaError( double err ) ;

    //! Get the sensor ID
    int getSensorID() const ;
    
    //! Get the offset along x
    double getXOffset() const;

    //! Get the offset along y
    double getYOffset() const;

    //! Get the offset along z
    double getZOffset() const;

    //! Get the theta along x
    double getAlpha() const;

    //! Get the theta along y
    double getBeta() const;

    //! Get the theta along z
    double getGamma() const;

    //! Get the error of the offset along x
    double getXOffsetError() const;

    //! Get the error of the offset along y
    double getYOffsetError() const;

    //! Get the error of the offset along z
    double getZOffsetError() const;

    //! Get the error of the theta along x
    double getAlphaError() const;

    //! Get the error of the theta along y
    double getBetaError() const;

    //! Get the error of the theta along z
    double getGammaError() const;

    //! Print the output
    /*! This method is used to print out the constant 
     *
     *  @param os The input output stream
     */
    virtual void print(std::ostream& os) const ;

    //! Overload of operator<<
    /*! This friend function is the overload of the operator << for
     *  the AlignmentConstant
     *
     *  @param os The input output stream as modified by the print
     *  method
     *  @param c The alignment constant
     *  @return The output stream
     *
     */ 
    friend std::ostream& operator<< (std::ostream& os, const AlignmentConstant & c) { c.print(os); return os; }

  };

  

}

#endif
