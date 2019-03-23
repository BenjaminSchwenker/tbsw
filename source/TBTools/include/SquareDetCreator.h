#ifndef SQUAREDETCREATOR_H
#define SQUAREDETCREATOR_H

#include "DetCreatorBase.h"

namespace depfet {
     
  /** Class SquareDetCreator
   *  
   *  The creator for are pixel detector with square checkerboard pixels. 
   *  
   *  @Author B. Schwenker, University of Göttingen (March 2019)
   *  <mailto:benjamin.schwenker@phys.uni-goettingen.de>
   */
  class SquareDetCreator : public DetCreatorBase {
   public:
    
    /** Default Constructor */
    SquareDetCreator() {}
    /** Default Destructor */
    ~SquareDetCreator() {}
    
    /**
     * Function to actually create the detector
     */
    void create(const TiXmlNode* content, std::vector<Det*>& dets) override;
  };
  
} //end of namespace depfet

#endif /* SQUAREDETCREATOR_H */
