#include "gearimpl/FixedPadAngleDiskLayout.h"

#include <cmath>
#include <math.h>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace gear {
  
  FixedPadAngleDiskLayout::~FixedPadAngleDiskLayout()
  {
      cleanup();
  }
      
  void FixedPadAngleDiskLayout::cleanup()  
  { 
    for( unsigned i=0; i<_padIndices.size(); ++i ){
      delete _padIndices[i] ;
    }
  }

    FixedPadAngleDiskLayout::FixedPadAngleDiskLayout( double   rMin, double   rMax, int nRow,
						      double phiMin, double phiMax, int nPadsInRow )
	: _rMin( rMin ),
	  _rMax( rMax ),
	  _nRow( nRow ) ,
	  _phiMin( phiMin ),
	  _phiMax( phiMax ),
	  _nPadsInRow( nPadsInRow ) 
  {
    
    _rowHeight =  ( rMax - rMin ) / _nRow ;
    _padAngle  =  (phiMax - phiMin) / nPadsInRow;

//    // minimal sanity check:
//    if ( _padHeight <= 0.0  || _padWidth <= 0.0 ) {
//
//      std::stringstream s ;
//      s << " FixedPadAngleDiskLayout: missing or inconsistent parameters"
//	<< " - padHeight: " << _padHeight 
//	<< " - padWidth : " << _padWidth ;
//
//      throw ParseException( s.str() ) ;
//    }

    _extent.resize(4) ;
    _extent[0] = rMin ;
    _extent[1] = rMax ;
    _extent[2] = phiMin ;
    _extent[3] = phiMax ;

//    // compute the row properties
//
//    int nr = (int) std::floor ( ( rMax - rMin ) /  _padHeight ) ;
//
//    _nRow = ( nRow > 0 ) ? nRow : nr ;
//
//    if( nRow > nr ) {
//
//      std::cout << " WARNING: FixedPadAngleDiskLayout() :  cannot place  " <<  nRow 
//		<< " rows of height " <<  _padHeight
//		<< " mm between rMin = " << rMin << " mm and rMax = " << rMax 
//		<< " mm -  will use nRow = " << nr << " !! " 
//		<< std::endl ;
//	
//      _nRow = nr ;
//
//    } 

    _padIndices.resize( nRow ) ;
    // initialise all pointers to 0
    // this thing is only filled when accessed by getPadsInRow() the first time
    // to allocate the memory only when needed
    for (std::vector< std::vector<int>* >::iterator rowIter = _padIndices.begin() ; 
	 rowIter < _padIndices.end(); rowIter++)
    {
	*rowIter=0;
    }


    _nPad  = nRow * nPadsInRow ;

//    for( int i = 0 ; i < _nRow ; i++ ) {
//
//      Row row ;
//      
//      row.RCenter  =  rMin + ( i * _rowHeight )  + .5 * _rowHeight ;
//
//      double u = 2. * row.RCenter * M_PI  ;
//
//      row.NPad = (int) std::floor( u / (_padWidth + padGap )  ) ; 
//
//      row.PhiPad = 2. * M_PI  / row.NPad ;
//
////       std::cout << "  row, row.NPad , row.PhiPad : " << i << " , " << row.NPad << " , " <<  row.PhiPad 
//// 		<< std::endl; 
//
//      _nPad += row.NPad ;
//
//      _rows.push_back( row ) ;
//    }
      }

  FixedPadAngleDiskLayout::FixedPadAngleDiskLayout( const FixedPadAngleDiskLayout & right) 
  {
      copy_and_assign(right);
  }

  FixedPadAngleDiskLayout & FixedPadAngleDiskLayout::operator = ( const FixedPadAngleDiskLayout & right) 
  {
      cleanup();
      copy_and_assign(right);
      return *this;
  }

  void FixedPadAngleDiskLayout::copy_and_assign(const  FixedPadAngleDiskLayout & right)
  {
      _rMin = right._rMin;
      _rMax = right._rMax;
      _nRow = right._nRow;
      _rowHeight = right._rowHeight;
      _phiMin = right._phiMin;
      _phiMax = right._phiMax;

      _nPad = right._nPad;
      _nPadsInRow = right._nPadsInRow;

      _extent = right._extent;

      // just set the _padIndices to an vector of empty entries.
      // In every instance, also the copied ones, it is filled only at first
      // usage to save memory
      // just like in the constructor
      _padIndices.resize( _nRow ) ;
      for (std::vector< std::vector<int>* >::iterator rowIter = _padIndices.begin() ; 
	   rowIter < _padIndices.end(); rowIter++)
      {
	  *rowIter=0;
      }
      
     
  }
  
  Vector2D FixedPadAngleDiskLayout::getPadCenter(int padIndex)  const {
    
    int rowNum =  getRowNumber( padIndex ) ;
    int padNum =  getPadNumber( padIndex ) ;

    double r = (static_cast<double>(rowNum)+0.5) * _rowHeight + _rMin;

    double phi =  ( static_cast<double>(padNum) + 0.5 ) * _padAngle + _phiMin ; 
    
    return Vector2D( r , phi ) ;
  }

  const std::vector<int>& FixedPadAngleDiskLayout::getPadsInRow(int rowNumber) const {

    static std::vector<int> empty ;
    try 
    {
      if( _padIndices.at(rowNumber )  == 0 ) 
      {
	  _padIndices[ rowNumber ]  = new std::vector<int>( _nPadsInRow )  ;
	
	for(int i = 0 ; i < _nPadsInRow ; i++)
	{
	  _padIndices[rowNumber]->operator[](i)  = getPadIndex( rowNumber , i  ) ;
	}
	
      }
    } catch( std::out_of_range& r) {
      
      //       std::cout << " FixedPadAngleDiskLayout::getPadsInRow : no row " << rowNumber << std::endl ;
      return empty ;

    }

    return *_padIndices[ rowNumber ] ;

  }

  int FixedPadAngleDiskLayout::getPadIndex(int rowNum, int padNum) const {

      if( rowNum  > _nRow - 1 )
      {
	  throw std::out_of_range(" FixedPadAngleDiskLayout::getPadIndex row number too large !");
      }
    
    if( padNum > _nPadsInRow - 1 ) {
      
      std::stringstream sstr ;
      
      sstr << "FixedPadAngleDiskLayout::getPadIndex: pad number too large: "
	   << padNum << ", only " <<   _nPadsInRow << " pads in row " << rowNum ;

      throw std::out_of_range( sstr.str() );
    }
    
    return  (rowNum << 16 ) | ( 0x0000ffff & padNum ) ;


  }

  int FixedPadAngleDiskLayout::getNearestPad(double r, double phi) const 
  {

    if( r < 0.0 ) 
	throw  std::out_of_range("FixedPadAngleDiskLayout::getNearestPad: radius must not be negative in polar coordinates!" );

    // wrap phi to +- 2*pi, which is the max range for phiMin and phiMax
    while( phi < -2*M_PI ) {  phi += 2. * M_PI  ; }
    while( phi > 2*M_PI ) {  phi -= 2. * M_PI  ; }

    if  ( (phi <= _phiMax) && (phi >= _phiMin) )
    {
	// This is easy, the pad number in the row is the one corresponding to the angle
	int padNumber = static_cast<int>( (phi - _phiMin) / _padAngle );
	
	int rowNum = r < _rMin ? 0 : (int) std::floor( ( r - _rMin ) / _rowHeight  ) ;
    
	if( rowNum >= _nRow  )       
	    rowNum = _nRow -1 ;

	return getPadIndex( rowNum , padNumber ) ;
    }

    // the difficult part has not been imlemented yet
    throw gear::NotImplementedException("FixedPadAngleDiskLayout::getNearestPad: coordinates out of plane extend are not implemented yet");

    // calculate pca to straight line with phi = _phiMin

    // determine r_intersect for x = r_intersect * cos(_phiMin) = x(r, phi) - m * sin(_phiMin)
    //                           y = r_intersect * sin(_phiMin) = y(r, phi) + m * cos(_phiMin)
//    double dist_right_edge = -1.;
//    if ( r_intersect >= rMin &&  r_intersect <= rMin)
//	// possible candidate for closest distance
//	dist_right_edge = m;
//    else
//
//    double dist_left_edge = -1.;
//    if ( r_intersect >= rMax &&  r_intersect <= rMax)
//	// possible candidate for closest distance
//	dist_left_edge = m;

    return 0;
//    return getPadIndex( rowNum , padNum ) ;
				  
  }

  int FixedPadAngleDiskLayout::getRightNeighbour(int padIndex) const {
    
    int pn = getPadNumber( padIndex ) + 1  ;
    int rn = getRowNumber( padIndex)  ;
    
    // anything larger than _nPadsInRow would have caused an out_of_range exception in
    // getPadNumber
    if(pn == _nPadsInRow)
    {
    // take into account wrap arround in case of full circle (incl. rounding errors)
	if (fabs(_phiMax - _phiMin - 2*M_PI) < 1e-6 )
	    pn = 0 ;
	else
	    throw Exception("RectangularPadRowLayout::getRightNeighbour: no right neighbour pad !");
    }

    return getPadIndex(  rn , pn ) ;
  }

  int FixedPadAngleDiskLayout::getLeftNeighbour(int padIndex) const {

    int pn = getPadNumber( padIndex ) - 1  ;
    int rn = getRowNumber( padIndex)  ;

    // anything larger than _nPadsInRow would have caused an out_of_range exception in
    // getPadNumber
    if (pn == _nPadsInRow)
    {
    // take into account wrap arround in case of full circle (incl. rounding errors)
	if (fabs(_phiMax - _phiMin - 2*M_PI) < 1e-6 )
	    pn = 0 ;
	else
	    throw Exception("RectangularPadRowLayout::getRightNeighbour: no right neighbour pad !");
    }

    return getPadIndex(  rn , pn ) ;
  }
  

  bool FixedPadAngleDiskLayout::isInsidePad(double r, double phi, int padIndex) const {
      
//     double rCenter =  0. ;
//     int rowNumber = getRowNumber( padIndex ) ;
//     try {  rCenter = _rows.at( rowNumber ).RCenter ; 
//     } 
//     catch(std::out_of_range){
//       return false ;
//     }
      Vector2D padCenter = getPadCenter( padIndex);

      if( (r  < padCenter[0] - _rowHeight/2.) || (r > padCenter[0] + _rowHeight/2.) ||
	  (phi < padCenter[1]- _padAngle/2.)  || (phi < padCenter[1]+_padAngle/2. ) )
	  return false ;
      else
	  return true;
  }

  bool FixedPadAngleDiskLayout::isInsidePad(double r, double phi) const {

      if( (r < _rMin) || (r > _rMax) || (phi < _phiMin) || (phi > _phiMax) )
	  return false ;
      else
	  return true ;
  }
  
} // namespace

