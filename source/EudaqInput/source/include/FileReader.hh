#ifndef EUDAQINPUT_INCLUDED_FileReader
#define EUDAQINPUT_INCLUDED_FileReader

#include "FileSerializer.hh"
#include "DetectorEvent.hh"
#include "counted_ptr.hh"
#include <string>

namespace eudaqinput {

  class FileReader {
  public:
    FileReader(const std::string & filename, const std::string & filepattern = "", bool synctriggerid = false);
    ~FileReader();
    bool NextEvent(size_t skip = 0);
    std::string Filename() const { return m_filename; }
    unsigned RunNumber() const;
    const eudaqinput::Event & GetEvent() const;
    const DetectorEvent & Event() const { return GetDetectorEvent(); } // for backward compatibility
    const DetectorEvent & GetDetectorEvent() const;
    void Interrupt() { m_des.Interrupt(); }
    
  private:
    std::string m_filename;
    FileDeserializer m_des;
    counted_ptr<eudaqinput::Event> m_ev;
    unsigned m_ver;
  };

}

#endif // EUDAQINPUT_INCLUDED_FileReader
