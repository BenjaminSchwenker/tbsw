#ifndef RectangularPadRowLayout_h
#define RectangularPadRowLayout_h 1

#include <vector>
#include "gear/PadRowLayout2D.h"


namespace gear {



/** Implementation of PadRowLayout2D for a rectangular row based layout where
 *  all pads in a given row are equal and have rectangular shape. 
 *  However rows can vary in height and layout.
 *  The pads with the given height and width are distributed in the rows symmetrically
 *  where an optional uninstrumented left or right offset is taken into account (this allows, e.g. 
 *  staggering of rows). Rows are specified/filled from the bottom up. 
 *  Uninstrumented areas on the plane will occur if a row height is given that exceeds the pad height 
 *  and/or if the pads don't fill the complete row (minus the left/right offset).<br>
 *  This code is in parts based on ideas and an extension of original implementations 
 *  by M. Weber RWTH Aachen and J. McGeachie UVic, Ca..
 * 
 * @author F. Gaede, DESY
 * @version $Id: RectangularPadRowLayout.h,v 1.4 2008-06-04 16:04:27 gaede Exp $
 */
  class RectangularPadRowLayout : public PadRowLayout2D {
    
  public:

    /** Internal helper class for RectangularPadRowLayout */
    struct Row{
      int NPad ;
      double PadWidth ;
      double PadHeight ;
      double Height ;
      double LeftOffset ;
      double RightOffset ;
      double Center ;
      double WidthPerPad ; // pad + insensitive ( gap ) 
    } ;
    
    
  protected:
    int _nRow ;
    int _nPad ;
    int _repeatRows ;
    std::vector<Row> _rows ;
    std::vector<double> _extent ;
    mutable std::vector< std::vector<int>* > _padIndices ;
    std::vector<unsigned> _nRows ;  // helper vector to keep track of equal rows

  public: 

    /** Construct the empty RectangularPadRowLayout with the width and x position
     *  specified through  xMin and xMax and an optional yMin - yMax will be computed from 
     *  the rows that are added with addRow().
     */
    RectangularPadRowLayout( double xMin, double xMax , double yMin=0.0) ;
    
    /// Destructor.
    virtual ~RectangularPadRowLayout() ; 
    

    /** Add nRow rows with the given parameters.
     */
    virtual void addRow(  int nRow, int nPad , double padWidth , double padHeight , 
		    double rowHeight =0.0, double leftOffset =0.0 , 
		    double rightOffset =0.0 );
    

    /** Repeat the current rows 'count' times - this allows to easily repeat a pattern of
     *  several rows, e.g. 2 rows offset to each other by half a padWidth (staggering).
     *  Can only be called once per layout, i.e. no multiple rpeat patterns are allowed.
     */
    virtual void repeatRows(unsigned count) ;

    /** Returns the number number for which a given row pattern has been repaeted 
     */
    int getRepeatRowCount() const { return _repeatRows ; } 


    /** The type of the row layout: PadRowLayout2D::CARTESIAN.
     */
    virtual int getPadLayoutType() const { return PadRowLayout2D::CARTESIAN ; } 
    
    /** The shape of the pads: PadRowLayout2D::RECTANGLE (i.e. keystone).
     */ 
    virtual int getPadShape() const { return PadRowLayout2D::RECTANGLE ; }
    
    /** The total number of pads in the TPC.
     */
    virtual int getNPads() const { return _nPad ; }
    
    /** The number of rows.
     */
    virtual int getNRows() const ;
    
    /** The row height in mm.
     */
    virtual double getRowHeight(int rowNumber) const ;

    /** The height of the pad in mm.
     */
    virtual double getPadHeight(int padIndex) const ;

    /** The width of the pad in radians NB! 
     */
    virtual double getPadWidth(int padIndex) const ;

    /** The center of the pad in 2d coordinates, (x,y) or (r,phi).
     */
    virtual Vector2D getPadCenter(int padIndex)  const;


    /** Indices of all pads in row rowNumber (row indices start from 0 at the
     * bottom (CARTESIAN) or at the center (POLAR)). 
     */
    virtual const std::vector<int>& getPadsInRow(int rowNumber) const ;

    /** Extent of the sensitive plane - [xmin,xmax,ymin,ymax] CARTESIAN or 
     *	[rmin,rmax,phimin,phimax] POLAR.
     */
    virtual const std::vector<double>& getPlaneExtent()  const { return _extent ; }

    /** The number of the row that contains the pad at padIndex - numbering starts at y==0 (bottom).
     */
    virtual int getRowNumber(int padIndex)  const ;

    /** The pad number (column) within the row  - numbering starts at x==0 (left).
     */
    virtual int getPadNumber(int padIndex) const { return ( 0x0000ffff & padIndex ) ; }

    /** Create a padIndex for the given row and pad ( column ) number 
     */
    virtual int getPadIndex(int rowNum, int padNum)  const;

    /** The index of the pad nearest to the given point in 2d coordinates
     * (x,y,) or (r,phi).
     */
    virtual int getNearestPad(double x, double y)  const;

	int getNearestPadOld(double x, double y) const;

    /** The index of the right neighbour pad.
     */
    virtual int getRightNeighbour(int padIndex)  const;

    /** The index of the left neighbour pad.
     */
    virtual int getLeftNeighbour(int padIndex)  const;

    /** True if coordinate (x,y) is within the given pad.
     */
    virtual bool isInsidePad(double x, double y, int padIndex)  const;

    /** True if coordinate (x,y) is within any pad.
     */
    virtual bool isInsidePad(double x, double y)  const;


    /** Helper method to identify equal row numbers in this layout (as they have been added).*/
    const std::vector<unsigned>& equalRowNumbers() const { return _nRows ; }

    /** Helper method with all row data.*/
    const std::vector<Row>& rows() const { return _rows ; }


  protected:
    
    /** Helper function for finding the nearest pad - returns the distance for points outside 
     *  the box and -1. for points inside.*/

    double distanceToBox( const Vector2D& p , double xMin, double yMin, double xMax, double yMax ) const ;

  }; // class

} // namespace gear
#endif // ifndef RectangularPadRowLayout_h
