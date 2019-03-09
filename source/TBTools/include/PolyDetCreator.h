#ifndef POLYDETCREATOR_H
#define POLYDETCREATOR_H

#include "DetCreatorBase.h"

namespace depfet {
     
  /** Class PolyDetCreator
   *  
   *  The creator for are pixel detector with polygonal pixels. 
   *  
   *  @Author H. C. Beck, University of Göttingen (March 2019)
   *  <mailto:helge-christoph.beck@phys.uni-goettingen.de>
   */
  class PolyDetCreator : public DetCreatorBase {
   public:
    
    /** Default Constructor */
    PolyDetCreator() {}
    /** Default Destructor */
    ~PolyDetCreator() {}
    
    /**
     * Function to actually create the detector
     */
    void create(const TiXmlNode* content, std::vector<Det*>& dets) override;
  };
  
} //end of namespace depfet

#endif /* POLYDETCREATOR_H */
