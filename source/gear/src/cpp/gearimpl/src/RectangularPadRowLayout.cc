#include "gearimpl/RectangularPadRowLayout.h"

#include <cmath>
#include <math.h>
#include <iostream>
#include <stdexcept>
#include <sstream>
#include <assert.h>
#include <cstdlib>

namespace gear {
  
    
  
  RectangularPadRowLayout::RectangularPadRowLayout( double xMin, double xMax , double yMin ) :

    _nRow(0) ,
    _nPad(0),
    _repeatRows(0)  {
    
    _extent.resize(4) ;
    _extent[0] = xMin ;
    _extent[1] = xMax ;
    _extent[2] = yMin ;
    _extent[3] = yMin ;
    
  }
  
  
  void RectangularPadRowLayout::addRow( int nRow, int nPad , double padWidth , double padHeight , 
				   double rowHeight , double leftOffset, 
				   double rightOffset ) {
    

    // a few sanity checks:
    assert( nRow > 0 ) ;
    assert( nPad > 0 ) ;
    assert( ! ( padWidth < 0. )  ) ;
    assert( ! ( padHeight < 0. )  ) ;
    assert( ! ( rowHeight < 0. )  ) ;
    assert( ! ( leftOffset < 0. )  ) ;
    assert( ! ( rightOffset < 0. )  ) ;


    double xMin = _extent[0] + leftOffset ;
    double xMax = _extent[1] - rightOffset ;

    double rowWidth = xMax - xMin ;

    if( nPad * padWidth   >  rowWidth ) {

      std::stringstream sstr ;
      sstr << "RectangularPadRowLayout::addRow: can't fit " 
	   << nPad << " pads in row with width " <<  rowWidth 
	   <<  " ( offsets left/right: " << leftOffset << " / "   << rightOffset << " ) " ; 

      throw Exception( sstr.str() ) ;
     }

    Row row ; 

    _nPad += abs( nRow ) * abs( nPad )  ;
    
    row.WidthPerPad = rowWidth / nPad ;

    if( rowHeight == 0.0 ) 
      rowHeight = padHeight ;
    
    if( rowHeight <  padHeight ){ 

      std::stringstream sstr ;
      sstr << "RectangularPadRowLayout::addRow: row height ( " 
	   << rowHeight << " ) can't be smaller than pad height ( " << padHeight << " ) !" ;

      throw Exception( sstr.str() ) ;
    }


    row.NPad = nPad  ;
    row.PadWidth = padWidth ;
    row.PadHeight = padHeight ;
    
    row.Height =  rowHeight ;

    row.LeftOffset =  leftOffset ;
    row.RightOffset = rightOffset  ;
    

    for( int i = 0 ; i < nRow ; i++ ) {

      row.Center = _extent[3] + rowHeight / 2. ;

      _extent[3]+= rowHeight ;

      _rows.push_back( row ) ;

      _padIndices.push_back( 0 ) ;
    }

    _nRows.push_back( nRow ) ; // keep track of nRows as they are added

  }

  void RectangularPadRowLayout::repeatRows(unsigned count) {

    if( _repeatRows > 0  ) // ignore all but first call
      return ;

    _repeatRows = count ;
    
    int nRow = _rows.size() ;
    
    for(unsigned i=0 ; i< count-1 ; ++i ){
      for(int j=0 ; j< nRow; ++j ){
	
	const Row& r = _rows[j] ;
	
	addRow( 1 , r.NPad , r.PadWidth ,r.PadHeight , 
		r.Height , r.LeftOffset, r.RightOffset ) ;
	
      }
    }
  }
  
  int RectangularPadRowLayout::getNRows() const {
    return _rows.size() ;
  }

  double RectangularPadRowLayout::getRowHeight(int rowNumber) const { 
    return _rows.at(rowNumber).Height ; 
  }
  
  double RectangularPadRowLayout::getPadHeight(int padIndex) const {
    
    int rowNum =  getRowNumber( padIndex ) ;

    return _rows[rowNum].PadHeight ; 
  }
  
  
  double RectangularPadRowLayout::getPadWidth(int padIndex) const {

    int rowNum =  getRowNumber( padIndex ) ;

    return _rows[rowNum].PadWidth ; 
  } 


  Vector2D RectangularPadRowLayout::getPadCenter(int padIndex)  const {
    
    int rowNum =  getRowNumber( padIndex ) ;
    int padNum =  getPadNumber( padIndex ) ;
    
    const Row& row =  _rows[ rowNum ] ;

    double x = _extent[0] + row.LeftOffset + padNum * row.WidthPerPad + row.WidthPerPad / 2. ;  

    double y = row.Center ; 

    return Vector2D( x , y  ) ;
  }
  

  RectangularPadRowLayout::~RectangularPadRowLayout() { 
    for( unsigned i=0; i<_padIndices.size(); ++i ){
      delete _padIndices[i] ;
    }
  }

  const std::vector<int>& RectangularPadRowLayout::getPadsInRow(int rowNumber) const {
    
    static std::vector<int> empty ;
    
    try {
      
      if( _padIndices.at(rowNumber )  == 0 ) {
	
	int nPad = _rows.at(rowNumber).NPad ;
	
	_padIndices[ rowNumber ]  = new std::vector<int>( nPad )  ;
	
	for(int i = 0 ; i < nPad ; i++){
	  
	  _padIndices[rowNumber]->operator[](i)  = getPadIndex( rowNumber , i  ) ;
	}
	
      }
    } catch( std::out_of_range& r) {
      
      //       std::cout << " RectangularPadRowLayout::getPadsInRow : no row " << rowNumber << std::endl ;
      return empty ;
      
    }
    
    return *_padIndices[ rowNumber ] ;
  }
  

  int RectangularPadRowLayout::getRowNumber(int padIndex)  const { 

    int rn = ( 0xffff0000 & padIndex ) >> 16 ; 

    if( rn < 0 ||  rn > int(_rows.size() - 1) ){
      
      std::stringstream sstr ;
      
      sstr << "RectangularPadRowLayout::getRowNumber: illegal rownumber: " 
	   << rn << " for padIndex  0x" << std::hex << padIndex << " nRows: " << _rows.size() << std::dec ;
      
      throw Exception( sstr.str() ) ;
    }

    return rn ;
  } 


  int RectangularPadRowLayout::getPadIndex(int rowNum, int padNum) const { 


    if( (unsigned) rowNum > _rows.size() - 1 ) {
      
      throw std::out_of_range(" RectangularPadRowLayout::getPadIndex row number too large !");
    }
    
    if( padNum > _rows[rowNum].NPad - 1 ) {
      
      std::stringstream sstr ;
      
      sstr << "RectangularPadRowLayout::getPadIndex: pad number too large: "
	   << padNum << " only " <<   _rows[rowNum].NPad << " pads in row " << rowNum ;

      throw std::out_of_range( sstr.str() );
    }
    
    return  (rowNum << 16 ) | ( 0x0000ffff & padNum ) ;

  }

int RectangularPadRowLayout::getNearestPadOld(double x, double y) const {
    
    //  move from outside of pad plane to edge of pad plane
    if(  x < _extent[0] ) x =  _extent[0] ;
    if(  x > _extent[1] ) x =  _extent[1] ;
    if(  y < _extent[2] ) y =  _extent[2] ;
    if(  y > _extent[3] ) y =  _extent[3] ;
    
    int rowNum = 0  ; 
    
    double dMin = 2. * ( _extent[3] -  _extent[2] ) ;
    
    // find closest row - NB row heights might vary and rows can have offsets !
    // thus we need to loop over all rows

    for( unsigned i=0 ; i < _rows.size() ; ++i ) {
      
      const Row& row =  _rows[i] ;

      double xmin =  _extent[0] +  row.LeftOffset ;
      double xmax =  _extent[1] -  row.RightOffset ;
      double ymin =  row.Center -  row.PadHeight / 2. ;
      double ymax =  row.Center +  row.PadHeight / 2. ;

      Vector2D  p( x, y ) ;

      double d ;

      if( ( d = distanceToBox( p , xmin, ymin, xmax, ymax ) ) < dMin ){

	dMin = d ; 
	rowNum = i ;


	if( dMin < 0.0 ) break ; // if distance is negative we are already inside a pad
      }
    }
    
    int padNum = 0 ;
    
    // now that we have the row we need the closest pad
    const Row& row =  _rows[rowNum] ;
    
    if( x < (  _extent[0] + row.LeftOffset  ) ){
      
      padNum = 0 ;  
        
    } else if( x >= (  _extent[1] - row.RightOffset ) ) {
      
      padNum = row.NPad - 1 ; 

    } else {

      padNum =  (int) ( ( x - (  _extent[0] + row.LeftOffset ) )  /  row.WidthPerPad ) ;  

    }

    return getPadIndex( rowNum , padNum ) ;
  }


int RectangularPadRowLayout::getNearestPad(double x, double y) const {

    //  move from outside of pad plane to edge of pad plane
    if(  x < _extent[0] ) x =  _extent[0] ;
    if(  x > _extent[1] ) x =  _extent[1] ;
    if(  y < _extent[2] ) y =  _extent[2] ;
    if(  y > _extent[3] ) y =  _extent[3] ;

	// make a best guess at the starting row by assuming all rows are the same height
    int row_number = (int)( round( (_rows.size() - 1) * (y - _extent[2]) / (_extent[3] - _extent[2]) ) );

    // do a binary search for the nearest row
	while( true ) {

		const Row &row = _rows[row_number];
		const double y_displacement = y - row.Center;
		
		//std::cout << "Checking row " << row_number << ", centre: " << row.Center << ", displacement: " << y_displacement << "mm" << std::endl;

		if( fabs(y_displacement) <= (row.PadHeight/2) ) {

			// this is the correct row
			break;
		}
		else if( fabs(y_displacement) <= (row.Height/2) ) {
			
			// check the next row of pads, it could be closer to them then this row
			const int next_row_number = (y_displacement < 0)? row_number - 1: row_number + 1;
			const Row &next_row = _rows[next_row_number];

			const double y_distance_from_pads = fabs(y_displacement) - (row.PadHeight/2);
			const double y_distance_from_next_pads = fabs(y - next_row.Center) - (next_row.PadHeight/2);

			row_number = (y_distance_from_next_pads < y_distance_from_pads)? next_row_number: row_number;
			break;
		}
		else {
			// it wasn't in this row... find the next best guess
			int remaining_rows;
			double remaining_distance;
			double offset;

			if( y_displacement < 0 ) {
				
				remaining_rows = row_number;
				remaining_distance = row.Center - _extent[2] - (row.Height/2.);
				offset = y_displacement + (row.Height/2);
			}
			else {
				
				remaining_rows = _rows.size() - row_number - 1;
				remaining_distance = _extent[3] - row.Center - (row.Height/2.);
				offset = y_displacement - (row.Height/2);
			}

			const double param = remaining_rows * (offset / remaining_distance);
			const int delta_row = (int)((param < 0)? floor(param): ceil(param));

			row_number += delta_row;
		}
		
		/*std::cout << "Hit Enter to continue... \n";
		std::cin.get();*/
	}

	// we have the row number... find the pad number
    int pad_number = 0 ;
    
    // now that we have the row we need the closest pad
    const Row& row =  _rows[row_number] ;
    
    if( x <= (  _extent[0] + row.LeftOffset  ) ) {
      
      pad_number = 0 ;  
        
    } else if( x >= (  _extent[1] - row.RightOffset ) ) {
      
      pad_number = row.NPad - 1;

    } else {

      pad_number =  (int) ( ( x - (  _extent[0] + row.LeftOffset ) )  /  row.WidthPerPad ) ;  

    }

    return getPadIndex( row_number , pad_number ) ;
  }

  
  int RectangularPadRowLayout::getRightNeighbour(int padIndex) const {
    
    int pn = getPadNumber( padIndex ) + 1  ;
    int rn = getRowNumber( padIndex)  ;
    
    int nPad = _rows.at(rn).NPad ;
    
    if( pn > nPad-1 ){ 
      throw Exception("RectangularPadRowLayout::getRightNeighbour: no right neighbour pad !");
//       throw std::out_of_range("RectangularPadRowLayout::getRightNeighbour: no right neighbour pad !");
    }
    
    return getPadIndex(  rn , pn ) ;
  }

  int RectangularPadRowLayout::getLeftNeighbour(int padIndex) const {

    int pn = getPadNumber( padIndex ) - 1  ;
    int rn = getRowNumber( padIndex)  ;
    
    if( pn < 0  ){ 
//       std::cerr << " hgoing top throw exception: getLeftNeighbour: no left neighbour pad !" << std::endl ;
      throw Exception("RectangularPadRowLayout::getLeftNeighbour: no left neighbour pad !");
//       throw std::out_of_range("RectangularPadRowLayout::getLeftNeighbour: no left neighbour pad !");
    }
    
    return getPadIndex(  rn , pn ) ;
  }


  bool RectangularPadRowLayout::isInsidePad(double x, double y, int padIndex) const {
      
    int pn = getPadNumber( padIndex ) ;
    int rn = getRowNumber( padIndex) ;

    const Row& row =  _rows[rn] ;

    // ---                                                                       ----  insensitive gap --------
    double pXMin = _extent[0] + row.LeftOffset +  pn    * row.WidthPerPad + (  row.WidthPerPad - row.PadWidth )/ 2.;
    double pXMax = _extent[0] + row.LeftOffset + (pn+1) * row.WidthPerPad - (  row.WidthPerPad - row.PadWidth )/ 2.;
    double pYMin = row.Center -  row.PadHeight / 2. ;
    double pYMax = row.Center +  row.PadHeight / 2. ;
    
    return ( pXMin <= x && x <= pXMax  &&
	     pYMin <= y && y <= pYMax ) ;
  }

  bool RectangularPadRowLayout::isInsidePad(double x, double y) const {

    //  outside of pad plane
    if( x < _extent[0] || x > _extent[1] || 
	y < _extent[2] || y > _extent[3] )
      
      return false ;
    
    return  isInsidePad( x , y , getNearestPad( x, y ) )  ;
  }
  
  double RectangularPadRowLayout::distanceToBox( const Vector2D& p , 
						 double xmin, double ymin, 
						 double xmax, double ymax ) const {
    
    // depending on the points position wrt. the box
    // the distance is computed in a different way:
    //          xmin        xmax
    //       x   |     x     |      x
    //     ------|===========|--------- ymax
    //       x   ||         ||      x
    //     ------|===========|--------- ymin
    //       x   |     x     |      x


    double x = p[0] ;
    double y = p[1] ;
    
    double d ;
    
    if( x < xmin ){
      
      if( y < ymin ) {
	
	d = sqrt( (xmin-x) * (xmin-x) + (ymin-y) * (ymin-y) ) ;
	
      } else if ( y > ymax ) {
	
	d =  sqrt( (xmin-x) * (xmin-x) + (ymax-y) * (ymax-y) ) ;
	
      } else { // ( ymin <= y && y <= ymax )
	
	d = xmin - x ;
      }
      
    } else if(  x > xmax ) {
      
      if( y < ymin ) {
	
	d = sqrt( (xmax-x) * (xmax-x) + (ymin-y) * (ymin-y) ) ;
	
      } else if ( y > ymax ) {
	
	d = sqrt( (xmax-x) * (xmax-x) + (ymax-y) * (ymax-y) ) ;
	
      } else { // ( ymin <= y && y <= ymax )
	
	d = x - xmax ;
      }
      
    } else { // ( xmin <= x && x <= xmax )
      
      if( y < ymin ) {
	
	d = ymin - y ; 
	
      } else if ( y > ymax ) {
	
	d = y - ymax ;
	
      } else { // ( ymin <= y && y <= ymax )
	
	d = - 1.  ;  // inside the box !
      }
    }
    return d ;
  }
  

} // namespace

