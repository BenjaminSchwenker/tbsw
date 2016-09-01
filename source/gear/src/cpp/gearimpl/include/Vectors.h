#ifndef GEAR_Vectors_H
#define GEAR_Vectors_H 1

#include <stdexcept> 
#include "Vector3D.h"

namespace gear {


  /** Trivial N-dimensional vector. 
   *  Access to coordinates is provided through the [] operator.
   *  Could be replaced or make use of a  a more elaborate implementation, 
   *  e.g. from CLHEP or root::MathCore in the future.
   * 
   * @author F. Gaede, DESY
   * @version $Id: Vectors.h,v 1.5 2007-03-16 15:43:17 gaede Exp $
   */
  
  template <int N, typename Float_T=double>
  class VectorND_T {
    
  public: 

    VectorND_T() {} ;

    /** The i-th coordinate */
    inline Float_T& operator[](unsigned i) {  
      
      if( i > N-1 ) throw std::out_of_range( "VectorND_T::operator[]" ) ;
      
      return _c[i]  ; 
    } 
    inline Float_T operator[](unsigned i) const {  
      
      if( i > N-1 ) throw std::out_of_range( "VectorND_T::operator[]" ) ;
      
      return _c[i]  ; 
    } 
    
  protected:
    
    Float_T _c[N] ; 
    
  }; // class
  

  /** Specialization for 2D */
  template <typename Float_T=double>
  struct Vector2D_T : public VectorND_T<2,Float_T> {
    typedef VectorND_T<2,Float_T> Base ;
    Vector2D_T() { 
      Base::_c[0] = 0. ; 
      Base::_c[1] = 0. ; 
    }
    Vector2D_T( Float_T c0, Float_T c1 ) {
      Base::_c[0] = c0 ; 
      Base::_c[1] = c1 ; 
    } 
  } ;
  
  /** Specialization for 3D */
  template <typename Float_T=double>
  struct Vector3D_T : public VectorND_T<3,Float_T> {
    typedef VectorND_T<3,Float_T> Base ;
    Vector3D_T() { 
      Base::_c[0] = 0. ; 
      Base::_c[1] = 0. ; 
      Base::_c[2] = 0. ;
    }
    Vector3D_T( Float_T c0, Float_T c1, Float_T c2 ) { 
      Base::_c[0] = c0 ;
      Base::_c[1] = c1 ; 
      Base::_c[2] = c2 ; 
    } 
  } ;
  


  typedef Vector2D_T<double> Vector2D ;
  
//   typedef Vector3D_T<double> Vector3D ;
  
  

//   /** Defined for backward compatibility */
//   struct Point2D : public Vector2D{
    
//     Point2D() : Vector2D( 0., 0. ) {}
//     Point2D(double d0, double d1 ) : Vector2D( d0, d1 ) {}

//   } ;

//   /** Defined for backward compatibility */
//   struct Point3D : public Vector3D{
    
//     Point3D() : Vector3D( 0., 0., 0.) {}
//     Point3D(double d0, double d1, double d2 ) : Vector3D( d0, d1,d2 ) {}

//   } ;


} // namespace gear

#endif /* ifndef GEAR_Vectors_H */

