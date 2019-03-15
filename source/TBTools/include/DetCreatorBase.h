#ifndef DETCREATORBASE_H
#define DETCREATORBASE_H

#include <vector>

class TiXmlNode;

namespace depfet {
  
  class Det;
     
  /** Class DetCreatorBase
   *  
   *  Pure virtual base class for all detector creators
   *  
   *  @Author B. Schwenker, University of Göttingen (March 2019)
   *  <mailto:benjamin.schwenker@phys.uni-goettingen.de>
   */
  class DetCreatorBase {
   public:
      
    /** Default Constructor */
    DetCreatorBase() {}
    /** Default Destructor */
    virtual ~DetCreatorBase() {}
       
    /**
     * Function to actually create the detector, has to be overridden by derived classes
     */
    virtual void create(const TiXmlNode* content, std::vector<Det*>& dets) = 0;
    
  };
  
} //end of namespace depfet

#endif /* DETCREATORBASE_H */
