#include "gearxml/RectangularPadRowLayoutXML.h"

#include "gearxml/XMLHandler.h"

#include "gearxml/tinyxml.h"
#include "gearimpl/RectangularPadRowLayout.h"
#include "gear/GearMgr.h"

namespace gear {
  

  TiXmlElement RectangularPadRowLayoutXML::toXML( const PadRowLayout2D* layout ) const {

    // Check if type is RectangularPadRowLayout
    
    const RectangularPadRowLayout* rectangularPadLayout=dynamic_cast<const RectangularPadRowLayout*>(layout);
    
    if (rectangularPadLayout==NULL) {
      
      throw ParseException("RectangularPadRowLayoutXML::toXML wrong type !");
    }
    
    // append data to PadRowLayout2D
    TiXmlElement padRowLayout2DXML("PadRowLayout2D");    
    
    padRowLayout2DXML.SetAttribute("type","RectangularPadRowLayout"); 
    
    padRowLayout2DXML.SetDoubleAttribute( "xMin", rectangularPadLayout->getPlaneExtent()[0] );
    padRowLayout2DXML.SetDoubleAttribute( "xMax", rectangularPadLayout->getPlaneExtent()[1] );
    padRowLayout2DXML.SetDoubleAttribute( "yMin", rectangularPadLayout->getPlaneExtent()[2] );

    int repeatRowCount =  rectangularPadLayout->getRepeatRowCount() ;

    if( repeatRowCount > 1 )
      padRowLayout2DXML.SetDoubleAttribute( "RepeatRows", repeatRowCount );


    const std::vector<unsigned>& rowNumbers = rectangularPadLayout->equalRowNumbers() ;
    int theRow = 0 ;

    int rowNum  = rowNumbers.size() ; 
    
    if( repeatRowCount > 1  ) 
      rowNum /=  repeatRowCount ;

    
    
    for(int i=0 ; i< rowNum ; ++i){
      
      TiXmlElement rowElement("row");    
      
      const RectangularPadRowLayout::Row& row = rectangularPadLayout->rows()[ theRow ]   ;
      
      int nRow =  rowNumbers[i]  ;
      
      if( nRow > 1 ) 
	rowElement.SetAttribute("repeat", nRow) ;
      
      rowElement.SetAttribute("nPad", row.NPad ) ;

      rowElement.SetDoubleAttribute("padWidth"    , row.PadWidth ) ; 
      rowElement.SetDoubleAttribute("padHeight"   , row.PadHeight ) ; 
      rowElement.SetDoubleAttribute("rowHeight"   , row.Height ) ; 
      rowElement.SetDoubleAttribute("leftOffset"  , row.LeftOffset ) ; 
      rowElement.SetDoubleAttribute("rightOffset" , row.RightOffset ) ; 
      
      padRowLayout2DXML.InsertEndChild( rowElement ) ;

      theRow += nRow ;
    }
    
    return padRowLayout2DXML ;
  }
  
    
  PadRowLayout2D* RectangularPadRowLayoutXML::fromXML( const TiXmlElement* xmlElement ) const {

 
    double xMin =  atof(  getXMLAttribute( xmlElement , "xMin" ) .c_str() ) ;
    double xMax =  atof(  getXMLAttribute( xmlElement , "xMax" ) .c_str() ) ;
    double yMin =  atof(  getXMLAttribute( xmlElement , "yMin" ) .c_str() ) ;
    
    int repeatRowCount = atoi(  getOptionalXMLAttribute( xmlElement , "repeatRows" , "1" ) .c_str() ) ;
    
    RectangularPadRowLayout* padRowLayout = new RectangularPadRowLayout( xMin, xMax, yMin ) ;

    // ----- rows ------------
    int rowId = 0 ;

    const TiXmlNode* row = 0 ;
    while( ( row = xmlElement->IterateChildren( "row", row ) )  != 0  ){

      int    nRow             = atoi(  getOptionalXMLAttribute( row , "repeat" , "1" ) .c_str() ) ;

      int    nPad             = atoi(  getXMLAttribute( row , "nPad"   ) .c_str() ) ;

      double padWidth    = atof(  getXMLAttribute( row , "padWidth"  ) .c_str() ) ;
      double padHeight   = atof(  getXMLAttribute( row , "padHeight" ) .c_str() ) ;

      double rowHeight   = atof(  getOptionalXMLAttribute( row , "rowHeight"  , "0.0" ).c_str() ) ;
      double leftOffset  = atof(  getOptionalXMLAttribute( row , "leftOffset" , "0.0" ).c_str() ) ;
      double rightOffset = atof(  getOptionalXMLAttribute( row , "rightOffset", "0.0" ).c_str() ) ;
      
      padRowLayout->addRow( nRow, nPad, padWidth, padHeight, rowHeight, leftOffset, rightOffset ) ;
	
      rowId+= nRow ;
    }

    padRowLayout->repeatRows( repeatRowCount ) ;
    
    return padRowLayout ;
  }

}

