#include "PluginManager.hh"
#include "Exception.hh"



#include "lcio.h"
#include "IMPL/LCEventImpl.h"
#include "IMPL/LCRunHeaderImpl.h"


//#include <iostream>
#include <string>
using namespace std;

namespace eudaqinput {

  PluginManager & PluginManager::GetInstance() {
    // the only one static instance of the plugin manager is in the getInstance function
    // like this it is ensured that the instance is created before it is used
    static PluginManager manager;
    return manager;
  }
  
  void PluginManager::RegisterPlugin(DataConverterPlugin * plugin) {
    m_pluginmap[plugin->GetEventType()] = plugin;
  }

  DataConverterPlugin & PluginManager::GetPlugin(const Event & event) {
    return GetPlugin(std::make_pair(event.get_id(), event.GetSubType()));
  }

  DataConverterPlugin & PluginManager::GetPlugin(PluginManager::t_eventid eventtype) {
    std::map<t_eventid, DataConverterPlugin *>::iterator pluginiter
      = m_pluginmap.find(eventtype);

    if (pluginiter == m_pluginmap.end()) {
      EUDAQINPUT_THROWX(FileReadException, "PluginManager::GetPlugin(): Unkown event type "+Event::id2str(eventtype.first)+":"+eventtype.second);
    }

    return *pluginiter->second;
  }

  void PluginManager::Initialize(const DetectorEvent & dev) {
    for (size_t i = 0; i < dev.NumEvents(); ++i) {
      const eudaqinput::Event & subev = *dev.GetEvent(i);
      GetInstance().GetPlugin(subev).Initialize(subev);
    }
  }

  unsigned PluginManager::GetTriggerID(const Event & ev) {
    return GetInstance().GetPlugin(ev).GetTriggerID(ev);
  }


  lcio::LCRunHeader * PluginManager::GetLCRunHeader(const DetectorEvent & bore) {
    IMPL::LCRunHeaderImpl * lcHeader = new IMPL::LCRunHeaderImpl;
    lcHeader->setRunNumber(bore.GetRunNumber());
    lcHeader->setDetectorName("EUTelescope");
    
    for (size_t i = 0; i < bore.NumEvents(); ++i) {
      const eudaqinput::Event & subev = *bore.GetEvent(i);
      GetInstance().GetPlugin(subev).GetLCIORunHeader(*lcHeader, subev);
    }
    return lcHeader;
  }


  lcio::LCEvent * PluginManager::ConvertToLCIO(const DetectorEvent & dev) {
    lcio::LCEventImpl * event = new lcio::LCEventImpl;
    event->setEventNumber(dev.GetEventNumber());
    event->setRunNumber(dev.GetRunNumber());
    event->setTimeStamp(dev.GetTimestamp());

    for (size_t i = 0; i < dev.NumEvents(); ++i) {
      ConvertLCIOSubEvent(*event, *dev.GetEvent(i));
    }

    return event;
  }


  
  
  void PluginManager::ConvertLCIOSubEvent(lcio::LCEvent & dest, const Event & source) {
    
      GetInstance().GetPlugin(source).GetLCIOSubEvent(dest, source);
    
  }

}//namespace eudaqinput
