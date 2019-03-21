#include "DetCreatorManager.h"
#include <memory>
#include <streamlog/streamlog.h>


using namespace std;

namespace depfet {

  DetCreatorManager& DetCreatorManager::getInstance() 
  {
    static unique_ptr<DetCreatorManager> instance(new DetCreatorManager());
    return *instance;
  }

  void DetCreatorManager::registerCreatorFactory(const std::string& name, DetCreatorFactory* factory)
  {
    getInstance().m_creatorFactories[name] = factory;
  }
  
  DetCreatorBase* DetCreatorManager::getCreator(const string& name)
  {
    DetCreatorManager& instance = getInstance();
    map<string, DetCreatorFactory*>::const_iterator it = instance.m_creatorFactories.find(name);
    if (it == instance.m_creatorFactories.end()) {
      streamlog_out (ERROR3) << "Could not find a detector creator named" << name << endl;
      return 0;
    }
    return it->second();
  }
}
