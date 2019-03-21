#ifndef DETCREATORFACTORY_H
#define DETCREATORFACTORY_H

#include "DetCreatorManager.h"

namespace depfet {
  
  class DetCreatorBase;

  /**
   * Very simple class to provide an easy way to register creators with the
   * CreatorManager. When defining a new creator, add a
   *
   * CreatorFactory<Classname> Classname_factory("CreatorName");
   *
   * or similiar to the source file to automatically provide the needed
   * factory function and register the creator with the CreatorManager
   * 
   * @Author B. Schwenker, University of Göttingen (March 2019)
   * <mailto:benjamin.schwenker@phys.uni-goettingen.de>
   */
  template<class T> struct DetCreatorFactory {
    /** Constructor to register the DetCreator with the DetCreatorManager */
    explicit DetCreatorFactory(const std::string& name)
    {
      DetCreatorManager::registerCreatorFactory(name, factory);
    }
    /**
     * Static factory function to return a new instance of the given Creator
     * class. Ownership of the object goes to the caller who is responsible
     * of freeing the creator once it is done
     */
    static DetCreatorBase* factory()
    {
      return new T();
    }
  };
  
} //end of namespace depfet

#endif /* DETCREATORFACTORY_H */
