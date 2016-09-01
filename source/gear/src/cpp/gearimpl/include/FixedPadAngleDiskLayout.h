#ifndef FixedPadAngleDiskLayout_h
#define FixedPadAngleDiskLayout_h 1

#include <vector>
#include "gear/PadRowLayout2D.h"


namespace gear {



/** Implementation of PadRowLayout2D for a disk with fixed angled keystone pads.
 *  The number of pads is distrubuted on a wedge shaped segment starting 
 *  with the right edge of the first pad in throw at phiMin and the left edge
 *  of the last pad in the row.
 *  All rows have the same height, starting at rmin for the lower edge of the
 *  first row to rmax for the upper edge of the last row.
 *
 *  Row height and pad height are identical, and pad gap is 0
 *  (the only thing that matters is the effective pitch)
 * 
 * @author M. Killenberg, University of Bonn
 * @version $$
 */
  class FixedPadAngleDiskLayout : public PadRowLayout2D {
    
  public:

//    /** Internal helper class for FixedPadAngleDiskLayout */
//    struct Row{
//      int NPad ;
//      double RCenter ;
//      double PhiPad ;
//    } ;
    

  protected:
    double _rMin ;
    double _rMax ;
    int _nRow ;
    double _rowHeight ;
    double _phiMin;
    double _phiMax;
//    double _padWidth ;
//    double _padHeight ;
//    double _padGap ;
    double _padAngle;
    int _nPadsInRow;
    int _nPad ;
//    std::vector<Row> _rows ;
    std::vector<double> _extent ;
    mutable std::vector< std::vector<int>* > _padIndices ;

  /** function to copy all internal variables, incl. the objects
   *  pointed to and owned by the ModularTPCParameters.
   *  Used by constructor and assigment operator to avoid code duplication
   */
  void copy_and_assign(const  FixedPadAngleDiskLayout & );

  /** function to delete all the objects
   *  pointed to and owned by the GearMgr.
   *  Used by desctructor and assigment operator to avoid code duplication
   */
  void cleanup();
    

  public: 

    /** Construct the FixedPadAngleDiskLayout from the given parameters rMin, rMax,
     *  padHeight and PadWidth where the maximum number of pads are fitted to the plane.
     *  @param nRow :  limits the number of rows - if given 
     *  @param padGap : additional gap between pads in r-phi in mm 
     */
    
  FixedPadAngleDiskLayout( double   rMin, double   rMax, int nRow,
			   double phiMin, double phiMax, int nPadsInRow ) ;

    /** The copy constructor.
     *	Needed because _padIndices allocates memory dynamically
     */
    FixedPadAngleDiskLayout( const FixedPadAngleDiskLayout &);
    
    /// The assignment operator
    FixedPadAngleDiskLayout & operator = ( const FixedPadAngleDiskLayout &);

    /// Destructor.
    virtual ~FixedPadAngleDiskLayout() ; 
    

    /** The gap width in mm that was given in the C'tor. */
    virtual double getPadGap() const { return 0 ; } 

    /** The fixed width of the pads in mm. 
     */
    virtual double getFixedPadAngle() const { return _padAngle ; }
    

    /** The type of the row layout: PadRowLayout2D::POLAR.
     */
    virtual int getPadLayoutType() const { return PadRowLayout2D::POLAR ; } 
    
    /** The shape of the pads: PadRowLayout2D::RECTANGLE (i.e. keystone).
     */ 
    virtual int getPadShape() const { return PadRowLayout2D::RECTANGLE ; }
    
    /** The total number of pads in the TPC.
     */
    virtual int getNPads() const { return _nPad ; }
    
    /** The number of rows.
     */
    virtual int getNRows() const { return _nRow ; }
    
    /** The row height in mm.
     */
    virtual double getRowHeight(int rowNumber) const { return _rowHeight ; }

    /** The height of the pad in mm.
     */
    virtual double getPadHeight(int padIndex) const { return _rowHeight ; }

    /** The width of the pad in radians NB! 
     */
    virtual double getPadWidth(int padIndex) const { return _padAngle;}

    /** The center of the pad in 2d coordinates, (x,y) or (r,phi).
     */
    virtual Vector2D getPadCenter(int padIndex)  const;

    virtual int getNPadsInRow() const {return _nPadsInRow; }

    /** Indices of all pads in row rowNumber (row indices start from 0 at the
     * bottom (CARTESIAN) or at the center (POLAR)). 
     */
    virtual const std::vector<int>& getPadsInRow(int rowNumber) const ;

    /** Extent of the sensitive plane - [xmin,xmax,ymin,ymax] CARTESIAN or 
     *	[rmin,rmax,phimin,phimax] POLAR.
     */
    virtual const std::vector<double>& getPlaneExtent()  const { return _extent ; }

    /** The number of the row that contains the pad at padIndex - numbering starts at r/y==0.
     */
    virtual int getRowNumber(int padIndex)  const { return ( 0xffff0000 & padIndex ) >> 16 ; } 

    /** The pad number (column) within the row  - numbering starts at phi/x =.
     */
    virtual int getPadNumber(int padIndex) const { return ( 0x0000ffff & padIndex ) ; }

    /** Create a padIndex for the given row and pad ( column ) number 
     */
    virtual int getPadIndex(int rowNum, int padNum)  const;

    /** The index of the pad nearest to the given point in 2d coordinates
     * (x,y,) or (r,phi).
     */
    virtual int getNearestPad(double c0, double c1)  const;

    /** The index of the right neighbour pad.
     */
    virtual int getRightNeighbour(int padIndex)  const;

    /** The index of the left neighbour pad.
     */
    virtual int getLeftNeighbour(int padIndex)  const;

    /** True if coordinate (c0,c1) is within the given pad.
     */
    virtual bool isInsidePad(double c0, double c1, int padIndex)  const;

    /** True if coordinate (c0,c1) is within any pad.
     */
    virtual bool isInsidePad(double c0, double c1)  const;

  }; // class

} // namespace gear
#endif // ifndef FixedPadAngleDiskLayout_h
