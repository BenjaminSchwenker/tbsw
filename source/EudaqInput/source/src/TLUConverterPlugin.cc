#include "DataConverterPlugin.hh"
#include "TLUEvent.hh"
#include <string>


#include <IMPL/LCEventImpl.h>
#include <lcio.h>


namespace eudaqinput{

  class TLUConverterPlugin : public DataConverterPlugin {
  public:
    TLUConverterPlugin() : DataConverterPlugin(Event::str2id("_TLU")) {}

    
    virtual unsigned GetTriggerID(const eudaqinput::Event & ev) const {
      return ev.GetEventNumber();
    }



    virtual bool GetLCIOSubEvent(lcio::LCEvent & result, const eudaqinput::Event & source) const {
      dynamic_cast<lcio::LCEventImpl &>(result).setTimeStamp(source.GetTimestamp());
      dynamic_cast<lcio::LCEventImpl &>(result).setEventNumber(source.GetEventNumber());
      dynamic_cast<lcio::LCEventImpl &>(result).setRunNumber(source.GetRunNumber());
      return true;
    }


  private:
    static TLUConverterPlugin const m_instance;
  };

  TLUConverterPlugin const TLUConverterPlugin::m_instance;

}//namespace eudaqinput
