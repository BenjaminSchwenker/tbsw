#include "FileReader.hh"
#include "FileNamer.hh"
#include "Event.hh"
#include <list>

namespace eudaqinput {
  
  
  
  namespace {

    static bool ReadEvent(FileDeserializer & des, int ver, eudaqinput::Event * & ev, size_t skip = 0) {
      if (!des.HasData()) {
        return false;
      }
      if (ver < 2) {
        for (size_t i = 0; i <= skip; ++i) {
          if (!des.HasData()) break;
          ev = EventFactory::Create(des);
        }
      } else {
        BufferSerializer buf;
        for (size_t i = 0; i <= skip; ++i) {
          if (!des.HasData()) break;
          des.read(buf);
        }
        ev = eudaqinput::EventFactory::Create(buf);
      }
      return true;
    }
  }
  
  FileReader::FileReader(const std::string & file, const std::string & filepattern, bool synctriggerid)
    : m_filename(FileNamer(filepattern).Set('X', ".raw").SetReplace('R', file)),
      m_des(m_filename),
      m_ev(EventFactory::Create(m_des)),
      m_ver(1)
  {
  }
  
  FileReader::~FileReader()
  {  
  }
  
  bool FileReader::NextEvent(size_t skip) {
    eudaqinput::Event * ev = 0;
    bool result = ReadEvent(m_des, m_ver, ev, skip);
    if (ev) m_ev = ev;
    return result;
  }

  unsigned FileReader::RunNumber() const {
    return m_ev->GetRunNumber();
  }

  const Event & FileReader::GetEvent() const {
    return *m_ev;
  }

  const DetectorEvent & FileReader::GetDetectorEvent() const {
    return dynamic_cast<const DetectorEvent &>(*m_ev);
  }

}
