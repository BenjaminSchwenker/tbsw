#ifndef DETCREATORMANAGER_H
#define DETCREATORMANAGER_H

#include <string>
#include <map>

namespace depfet {
  
  class DetCreatorBase;

  /** Class DetCreatorManager
   *  
   *  Class to manage all detector creators and provide factory access 
   *
   *  @Author B. Schwenker, University of Göttingen (March 2019)
   *  <mailto:benjamin.schwenker@phys.uni-goettingen.de>
   */
  class DetCreatorManager {
   public:
    /** Typedef for a factory function */
    typedef DetCreatorBase* DetCreatorFactory();
    
    /**
     * Register a new creator by providing a name and a pointer to a factory
     * for this kind of creator.
     */
    static void registerCreatorFactory(const std::string& name, DetCreatorFactory* factory);

    /**
     * Return a new instance of a creator with the given name. 
     *
     * Returns 0 if no creator with the given name is registered. Ownership
     * of the creator is transferred to the caller who is responsible of
     * freeing the creator.
     */
    static DetCreatorBase* getCreator(const std::string& name);
    
   protected:
    /** singleton, hide constructor */
    DetCreatorManager() {}
    /** singleton, hide copy constructor */
    DetCreatorManager(const DetCreatorManager&) = delete;
    /** singleton, hide assignment operator */
    void operator=(const DetCreatorManager&) = delete;
    /** getter for the singleton instance */
    static DetCreatorManager& getInstance();
    /** Static map to hold all registered factories */
    std::map<std::string, DetCreatorFactory*> m_creatorFactories;
  };  
} //end of namespace depfet

#endif /* DETCREATORMANAGER_H_ */
