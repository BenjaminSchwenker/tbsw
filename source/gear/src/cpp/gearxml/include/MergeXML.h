#ifndef GEAR_MergeXML_H
#define GEAR_MergeXML_H

#include <string>
#include "gearxml/tinyxml.h"

#define NAME_STR "name"

namespace gear {

  class TiXmlDocument ;

  /**class to merge two xml files.
   *
   *@author R. Lippe, DESY
   *@version $Id:
   */
  class MergeXML {

  public:

    /**Sets first filename
     *returns true if file can be read.
     */
    bool setFile1( const std::string& fileName ) ;

    /**Sets second filename
     *returns true if file can be read.
     */
    bool setFile2( const std::string& fileName ) ;

    /**determines what file is dominant.
     *if nothing is set, file2 will overrule file1.
     *argument is only 1 or 2
     */
    void makeDominant(const int domFileNr )
      {
	if (domFileNr==1 || domFileNr==2) { 
	  _domFile = domFileNr ;
	}
      }

    /**merges files
     */
    bool mergeFiles (const std::string& fileName ) ;

    /**constructor
     */
    MergeXML() : 
      _domFile(2),
      _name(NAME_STR) {
    }
    

  private:

    //checks whether file can be read in
    bool checkXML( const std::string& fileName ) ;

    //get files correct: dominant will be _file2 and recessiv will be _file1
    void correctDominantFile() ;

    //compare nodes
    void mergeNode( TiXmlNode* node, TiXmlNode* domNode ) ;

    // merge childsNodes where domNode overrules Node
    void mergeElement( TiXmlNode* node, TiXmlNode* domNode ) ;

    // merge childElements where domNode overrules Node
    void mergeAttribute( TiXmlElement* element, TiXmlElement* domElement ) ;

    // count the number of nodes that look like childNode in parentNode
    int countNodes( TiXmlNode* childNode, TiXmlNode* domNode ) ;

    // get chilElement with same Attribute name
    TiXmlNode* getChildNode( TiXmlNode* node, TiXmlNode* domNode, int getOccurence = 1 ) ;

    // get the Name of the Node - either as value of the attribute name or as node->Value() 
    std::string getNodeName( TiXmlNode* node ) ;

    // check if the Node has the required attribute
    bool hasAttribute( TiXmlNode* node ) ;

    //file1
    TiXmlDocument* _file1 ;

    //file2
    TiXmlDocument* _file2 ;

    //fileSwap - if dominant changes
    TiXmlDocument* _fileSwap ;

    //fileMerged - final document
    TiXmlDocument* _fileMerged ;

    //dominant file - default file2
    int _domFile ;

    // unique attribute as identifier for nodes
    const char* _name ;
    
  }; //class

} //namespace
    
#endif
