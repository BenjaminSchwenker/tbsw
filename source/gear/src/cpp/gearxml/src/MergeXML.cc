#include "gear/GEAR.h"
#include "gearxml/MergeXML.h"

#include <sstream>

namespace gear{


  bool MergeXML::setFile1( const std::string& fileName ) {

    _file1 = new TiXmlDocument( fileName ) ;
    
    return _file1->LoadFile( fileName ) ;
  }


  bool MergeXML::setFile2( const std::string& fileName ) {
    
    _file2 = new TiXmlDocument( fileName ) ;
    
    return _file2->LoadFile( fileName ) ;
  }


  void MergeXML::correctDominantFile() {
    
    if ( _domFile == 1 ) {
      _fileSwap = new TiXmlDocument ;
      *_fileSwap = *_file1 ;
      *_file1 = *_file2 ;
      *_file2 = *_fileSwap ;
    }
  }


  bool MergeXML::mergeFiles (const std::string& fileName) {

    // debug information
    //std::cout << "mergeFiles called for "<< fileName << std::endl ;

    // variables
    const std::string gear = "gear" ;
    const std::string detectors = "detectors" ;

    // correct the files for dominancy
    correctDominantFile() ;

    // copy all dominant information in final file
    _fileMerged = new TiXmlDocument( fileName ) ;
    *_fileMerged = *_file2 ;
    
    // check if files are geartype
    if ( not ( _file1->FirstChild( gear ) && _fileMerged->FirstChild( gear )
	       && (_file1->FirstChild( gear ))->FirstChild( detectors ) 
	       && (_fileMerged->FirstChild( gear ))->FirstChild( detectors )
	       )) {
      
      // give out information and quit
      std::cout << "Files are not of type gear.\nAborted." <<std::endl ;
      return false;
    }

    // go to specified nodes
    TiXmlNode* nodeGear1 = _file1->FirstChild( gear ) ;
    TiXmlNode* nodeDetectors1 = nodeGear1->FirstChild( detectors ) ;

    TiXmlNode* nodeGearMerged = _fileMerged->FirstChild( gear ) ;
    TiXmlNode* nodeDetectorsMerged = nodeGearMerged->FirstChild( detectors ) ;
    

    // loop through all nodes
    mergeNode( nodeDetectors1, nodeDetectorsMerged ) ;
    
    // save file
    return _fileMerged->SaveFile( fileName ) ;
  }


  // compare nodes
  void MergeXML::mergeNode( TiXmlNode* node, TiXmlNode* domNode ) {

    // during the loop count occurences of same Name
    int occurence = 0 ;
    std::string lastName ;


    // debug information
//     if ( node->ToComment() ) std::cout << "node is Comment." << std::endl;
//     if ( node->ToText() )  std::cout << "node is Text." << std::endl;
//     if ( node->ToElement() )  std::cout << "node is Element." << std::endl;
//     if ( node->ToUnknown() )  std::cout << "node is Unknown." << std::endl;
    
    // skip if node is comment or pure text
    if (!( !( node->ToText() ) || !( node->ToComment() ) || !( node->ToUnknown() ))  ) {
      // debug information
      //std::cout << "skipping node "<< getNodeName(node)  <<" b/c it is text or comment" << std::endl ;
      return ;
    }

    // loop through all child nodes
    TiXmlNode* loopNode = 0 ;
    while (( loopNode = ( node->IterateChildren( loopNode ) )))
      {
	// debug information
	//std::cout << " loopNode " << loopNode -> Value() << std::endl ;

	// check how often nodes appear
	int countRec = countNodes( loopNode, node ) ;
	int countDom = countNodes( loopNode, domNode ) ;

	// debug information
	//std::cout << " counter in rec: " << countRec << "   counter in dom: " << countDom << std::endl ;

	// if multiple occurences and number does not match throw
	if ( ((countRec > 1) || (countDom >1)) &&  ( countDom != countRec ) ) {

	    // throw 
	  std::stringstream s ;
	  s << "ERROR. Counter for " << loopNode->Value() << " returns " << countRec
	    << " while counter for domNode " << domNode->Value() << " returns " 
	    << countDom << "." << std::endl ;

	  throw Exception( s.str() ) ;

	    return ;
	}
	

	// count occurences of same name
	if ( lastName == getNodeName( loopNode ) ) {
	     occurence ++ ;
	}
	else {
	  occurence = 1 ;
	}
	lastName = getNodeName( loopNode ) ;

	// get the appropriate child note from dominant node
	TiXmlNode* domChildNode = getChildNode( loopNode, domNode, occurence ) ;

	if (domChildNode == NULL ) {

	  // add to domNode
	  domNode->InsertEndChild( *loopNode ) ;

	  // debug information
	  //std::cout << " Added node " << lastName << std::endl ;
	}
	else {
	  // check if nodes are elements
	  if (loopNode->ToElement() && domChildNode->ToElement() ) {

	    // merge the attributes of these elements
	    mergeAttribute( (loopNode->ToElement()), (domChildNode->ToElement()) ) ;
	  }
 
	  // start mergeNode for subordniate nodes
	  if ( hasAttribute(loopNode) && !((loopNode->FirstChild())->ToText()) ) {
	    mergeNode( loopNode, domChildNode ) ;
	  }
	    //else {
	    // std::cout << "skipping to check childnodes because hasAttribtue: "
	    //	      << hasAttribute(loopNode) << "    ToText: " 
	    //      << !((loopNode->FirstChild())->ToText()) 
	    //      << "." << " true = " << true << std::endl ;
	    //  }

	}
      } // close Loop
  }


  // returns name for Node
  std::string MergeXML::getNodeName( TiXmlNode* node ) {

    // if name has to be taken than get name - otherwise take value of node  
    std::string strName ;
     
    if ( hasAttribute(node) ) { 
      strName = (node->ToElement())->Attribute( _name ) ;
    }
    else {
      strName = node->Value() ;
    }

    return strName ;
  }


  // check if node has attribute name
  bool MergeXML::hasAttribute( TiXmlNode* node ) {
    
    return ( (node->ToElement()) && ((node->ToElement())->Attribute(_name)) ) ;
  }


  // counts how often nodes appear
  int MergeXML::countNodes( TiXmlNode* childNode, TiXmlNode* parentNode ) {

    // counter for appearances
    int counter = 0 ;

    // get Name from Node either as value or as attribute name
    std::string strName = getNodeName( childNode ) ;

    // loop through all nodes and count
    TiXmlNode* loopNode = parentNode->FirstChild() ;
    while ( loopNode ) {
      
      // debug information
      // int i = 0 ; i++ ;
      // std::cout << "Iterating " << i << std::endl ;

      std::string loopName = getNodeName( loopNode ) ;

      // if we have to look for name
      if ( ( loopName == strName ) && ( hasAttribute(childNode) == hasAttribute(loopNode) ) ) {
	// count up
	counter ++ ;
      }
      
      loopNode = loopNode->NextSibling() ;
	
    } // close loop	
    return counter ;
  }
 

  // get ChildElement
  TiXmlNode* MergeXML::getChildNode( TiXmlNode* node, TiXmlNode* domNode, int getOccurence ) {

    // debug information
    //std::cout << "   getChildElement called for " << node->Value() << std::endl ;

    // get attribute name
    std::string searchName = getNodeName( node ) ;
      
    // debug information
    //std::cout << "     searching " << searchName << std::endl ;

    // get first child element
    // debug information
    //std::cout << "     Starting to loopElements " << node->Value() << std::endl ;

    TiXmlNode* loopNode = domNode->FirstChild() ;

    // the occurences of the same element are numbered
    int occurence = 0 ;

    while (loopNode )
      {

	std::string loopName = getNodeName( loopNode ) ;
	// debug information
	//std::cout << "       search " << searchName << " in " << loopName << std::endl ;
	
	// if found return
	// if names and attributes match
	if( (loopName == searchName) && ( hasAttribute(node) == hasAttribute(loopNode) ) ) {
	  occurence ++ ;
	  // debug information
	  //std::cout << "      found "<< occurence << " of " << getOccurence << std::endl ;
	  if ( occurence == getOccurence ) {
	    // debug information
	    //std::cout << "       returning in occurence " << occurence << std::endl ;
	    return loopNode ;
	  }
	}
	
	// continue loop
	loopNode = loopNode->NextSibling() ;

      } // close loop
    // debug information
    //std::cout << "     returning NULL" << std::endl ;
    return NULL ;
  }
    

  // compare Attributes
  void MergeXML::mergeAttribute( TiXmlElement* element, TiXmlElement* domElement ) {

    // loop through all attributs of given node
    TiXmlAttribute* loopAttribute = element->FirstAttribute() ;
    while ( loopAttribute ) 
      {
	// which Attribut
	std::string attName = loopAttribute->Name() ;
	
	// debug information
	//std::cout << "     loopAttribut " << attName ;

	// get same attribute in dominant element
	const char* domAttribute = domElement->Attribute( attName ) ;

	if( domAttribute == NULL ) {

	  // add to domElement 
	  domElement->SetAttribute( attName, loopAttribute->Value() ) ;
	  
	  // debug information
	  //std::cout << " added." << std::endl ;
	}
	else {
	  
	  // debug information
	  //std::cout << " exists." << std::endl ;
	}

	// go to next Attribute
	loopAttribute = loopAttribute->Next() ;
      }
  }

}  //namespace gear







