#ifndef Vector3D_h
#define Vector3D_h 1

#include <math.h>
#include <iostream> 


namespace gear {

/** Simple three dimensional vector providing the components for cartesian, 
 *  cylindrical and spherical coordinate  systems - internal reperesentation is
 *  cartesian.
 *  @author F. Gaede, DESY
 *  @version $Id: Vector3D.h,v 1.1 2007-03-07 15:37:45 gaede Exp $
 */

class Vector3D{
  
public:
  
  /** Default c'tor - zero vector */
  Vector3D() : _x(0.0),_y(0.0),_z(0.0) {}
  
  
  /** Templated c'tor - allows to have overloaded c'tors for different coordinates */
  template <class T>
  Vector3D( double x,double y, double z , T(&)() ) ;
  

  /** Default corrdinate system for initialization is cartesian */
  Vector3D( double x,double y, double z ) :
    _x(x),
    _y(y),
    _z(z) {
  }
  
  
  template <class T>
  /**Copy c'tor for three vectors from other packages - requires T::x(),T::y(), T::z().
   */
  Vector3D( const T& t) :
    _x( t.x() ) , 
    _y( t.y() ) , 
    _z( t.z() ){
  }
  
  /** Cartesian x coordinate */
  inline double x() const { return  _x ; }

  /** Cartesian y coordinate */
  inline double y() const { return  _y ; }

  /** Cartesian cartesian  z coordinate */
  inline double z() { return  _z ; }

  /** Assign to  cartesian x coordinate */
  inline double x() { return  _x ; }

  /**  Assign to  cartesian y coordinate */
  inline double y() { return  _y ; }

  /**  Assign to cartesian z coordinate */
  inline double z() const { return  _z ; }
  

  /** Accessing x,y,z with bracket operator */
  inline double operator[](int i) const {
    switch(i) {
    case 0: return _x ; break ;
    case 1: return _y ; break ;
    case 2: return _z ; break ;
    }
    return 0.0 ;
  }
  /** Accessing x,y,z with bracket operator for assignment */
  inline double& operator[](int i) {
    switch(i) {
    case 0: return _x ; break ;
    case 1: return _y ; break ;
    case 2: return _z ; break ;
    }
    static double dummy(0.0) ;
    return dummy ;
  }
  
  /** Azimuthal angle - cylindrical and spherical */
  inline double phi() const {
    
    return _x == 0.0 && _y == 0.0 ? 0.0 : atan2(_y,_x);
  }
  
  /** Transversal component - cylindrical 'r' */
  inline double rho() const {
    
    return trans() ;
  }
  
  /** Transversal component */
  inline double trans() const {
    
    return sqrt( _x*_x + _y*_y ) ; 
  }

  /** Transversal component squared */
  inline double trans2() const {
    
    return  _x*_x + _y*_y  ;
  }
  
  /** Spherical r/magnitude */
  inline double r() const {
    
    return sqrt( _x*_x + _y*_y + _z*_z ) ; 
  }
  
  
  /** Spherical r/magnitude, squared */
  inline double r2() const {
    
    return  _x*_x + _y*_y + _z*_z  ; 
  }
  
  /** Polar angle - spherical */
  inline double theta() const {
    
    return _x == 0.0 && _y == 0.0 && _z == 0.0 ? 0.0 : atan2( rho(),_z) ;
  }
  
  /** Scalar product */
  inline double dot( const Vector3D& v) { 
    return _x * v.x() + _y * v.y() + _z * v.z() ;
  }
  
  /** Vector product */
  inline Vector3D cross( const Vector3D& v) { 

    return Vector3D( _y * v.z() - _z * v.y() ,
		     _z * v.x() - _x * v.z() ,
		     _x * v.y() - _y * v.x() )  ;
  }

  /** Parallel unit vector */
  inline Vector3D unit() { 

    double n = r() ;
    return Vector3D( _x / n , _y / n , _z / n ) ;
  }


  /** Implicit templated conversion to anything that has a c'tor T(x,y,z) 
   *  and accessor functions x(),y(),z(). For safety the result is checked which 
   *  causes a small performance penalty.
   *  @see to()
   *  
   */
  template <class T>
  inline operator T() { 

    T t( _x, _y , _z ) ;

    assert( t.x()== _x && t.y()== _y && t.z()== _z ) ;

    return t ;

//     return T( _x, _y, _z ) ; 
  } 


  /** Explicit, unchecked conversion to anything that has a c'tor T(x,y,z). 
   *  Example: 
   *  CLHEP::Vector3D clhv = v.to< CLHEP::Vector3D>() ;
   *  @see operator T()
   */
  template <class T>
  inline T to() { return T( _x, _y, _z ) ; } 


protected:
  
  double _x,_y,_z ;
  

  // helper classes and function to allow 
  // different c'tors selected at compile time
public:
  
  struct Cartesian   { } ;
  struct Cylindrical { } ;
  struct Spherical   { } ;
  
  static Cartesian   cartesian()  { return Cartesian() ; }
  static Cylindrical cylindrical(){ return Cylindrical() ;}
  static Spherical   spherical()  { return Spherical() ; } 

} ;

/** Addition of two vectors */
inline Vector3D operator+(  const Vector3D& a, const Vector3D& b ) { 
  
  return Vector3D( a.x() + b.x()  , a.y() + b.y(), a.z() + b.z()  ) ;
}
/** Subtraction of two vectors */
inline Vector3D operator-(  const Vector3D& a, const Vector3D& b ) { 
  
  return Vector3D( a.x() - b.x()  , a.y() - b.y(), a.z() - b.z()  ) ;
}

/** Multiplication with scalar */
inline Vector3D operator*( double s , const Vector3D& v ) { 
  
  return Vector3D( s * v.x()  , s * v.y()  ,  s * v.z() ) ;
}


// template specializations for constructors of  different coordinate systems

/** Cartesian c'tor  - example: <br> 
 *  Vector3D  v( x, y, c , Vector3D::cartesian ) ;
 */
template <>
inline Vector3D::Vector3D( double x,double y, double z, Vector3D::Cartesian (&)() ) : 
  _x(x),
  _y(y),
  _z(z) {
}

/** Cylindrical c'tor  - example: <br> 
 *  Vector3D  v( rho, phi, z , Vector3D::cylindrical ) ;
 */
template <>
inline Vector3D::Vector3D( double rho,double phi, double z, Vector3D::Cylindrical (&)() ) : _z(z)  {
  
  _x = rho * cos( phi ) ;
  _y = rho * sin( phi ) ;
}


/** Spherical c'tor  - example: <br> 
 *  Vector3D  v( r, phi, theta , Vector3D::spherical ) ;
 */
template <>
inline Vector3D::Vector3D( double r,double phi, double theta, Vector3D::Spherical (&)() ) {
  double rst =  r * sin( theta ) ;
  _x = rst * cos( phi ) ;
  _y = rst * sin( phi ) ;
  _z = r * cos( theta ) ;
}



/** Output operator */
  std::ostream & operator << (std::ostream & os, const Vector3D &v) ;



// --- optionally one could have additional vectors with internal representations
// --- in other coordinates, e.g. cylindrical :

// /** Native cylindrical Vector 3D - internal data is rho, phi, z 
//  */
// class Vector3DCylindrical {
  
  
// public:
  
//   /** Conversion c'tor */
//   Vector3DCylindrical( const Vector3D& v ) : 
//     _rho( v.rho() ), 
//     _phi( v.phi() ), 
//     _z(   v.z()   ) { 
//   }

//   template <class T>
//   Vector3DCylindrical( double rho, double phi, double z , T(&)() ) ;
  
//   Vector3DCylindrical( double rho, double phi, double z ) : 
//     _rho(rho),_phi(phi),_z(z){ 
//   }
  
    
//   inline double x() const { return  _rho * cos( _phi ) ; }
  
//   inline double y() const { return  _rho * sin( _phi )  ; }
  
//   inline double z() const { return  _z ; }
  
//   inline double phi() const { return _phi ; }
  
//   inline double rho() const { return _rho ; }
    
//   inline double r() const {
    
//     return sqrt( _rho*_rho + _z*_z ) ; 
//   }
  
//   inline double theta() const {
    
//     return _rho == 0.0 && _z == 0.0 ? 0.0 : std::atan2( _rho ,_z);
//   }
  

//   /** Conversion operator */
//   operator Vector3D() {  return Vector3D( x() , y() , z() ) ; }

// protected:
  
//   double _rho,_phi,_z ;
  
// } ;

// // template specializations for constructors of  different coordinate systems
// template <>
// Vector3DCylindrical::Vector3DCylindrical( double x,double y, double z, 
// 					  Vector3D::Cartesian (&)() ) : _z(z) { 
  
//   _rho =  sqrt( x*x + y*y ) ;

//   _phi = ( x == 0.0 && y == 0.0 ? 0.0 : std::atan2(y,x) ) ;

// }

// template <>
// Vector3DCylindrical::Vector3DCylindrical( double rho,double phi, double z,  
// 					  Vector3D::Cylindrical (&)() ) :
//   _rho(rho),_phi(phi),_z(z){ 
// }


// template <>
// Vector3DCylindrical::Vector3DCylindrical( double r,double phi, double theta, 
// 					  Vector3D::Spherical (&)() ) : _phi( phi) {
  
//   _rho =  r * sin( theta ) ;
//   _z =    r * cos( theta ) ;
// }


} // namespace
#endif
