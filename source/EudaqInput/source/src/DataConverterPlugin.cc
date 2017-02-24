#include "DataConverterPlugin.hh"
#include "PluginManager.hh"

#include <iostream>

namespace eudaqinput {

  unsigned DataConverterPlugin::GetTriggerID(eudaqinput::Event const &) const {
    return (unsigned)-1;
  }

  DataConverterPlugin::DataConverterPlugin(std::string subtype)
    : m_eventtype(make_pair(Event::str2id("_RAW"), subtype))
  {
    //std::cout << "DEBUG: Registering DataConverterPlugin: " << Event::id2str(m_eventtype.first) << ":" << m_eventtype.second << std::endl;
    PluginManager::GetInstance().RegisterPlugin(this);
  }

  DataConverterPlugin::DataConverterPlugin(unsigned type, std::string subtype)
    : m_eventtype(make_pair(type, subtype))
  {
    //std::cout << "DEBUG: Registering DataConverterPlugin: " << Event::id2str(m_eventtype.first) << ":" << m_eventtype.second << std::endl;
    PluginManager::GetInstance().RegisterPlugin(this);
  }

}//namespace eudaqinput
