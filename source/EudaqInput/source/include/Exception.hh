#ifndef EUDAQINPUT_INCLUDED_Exception
#define EUDAQINPUT_INCLUDED_Exception

#include "Utils.hh"
#include <exception>
#include <string>

#ifndef EUDAQINPUT_FUNC
# define EUDAQINPUT_FUNC ""
#endif

#define EUDAQINPUT_THROWX(exc, msg) throw ::eudaqinput::InitException(exc(msg), __FILE__, __LINE__, EUDAQINPUT_FUNC)


#define EUDAQINPUT_EXCEPTIONX(name, base) \
  class name : public base {         \
  public:                            \
  name(const std::string & msg)      \
    : base(msg) {}                   \
  }

#define EUDAQINPUT_EXCEPTION(name) EUDAQINPUT_EXCEPTIONX(name, ::eudaqinput::Exception)

namespace eudaqinput {

  class Exception : public std::exception {
  public:
    Exception(const std::string & msg);
    const char * what() const throw() {
      if (m_text.length() == 0) make_text();
      return m_text.c_str();
    }
    // This shouldn't really be const, but it must be callable on temporary objects...
    const Exception & SetLocation(const std::string & file = "",
                                  unsigned line = 0,
                                  const std::string & func = "") const;
    virtual ~Exception() throw() {
    }
  protected:
    std::string m_msg;
  private:
    void make_text() const;
    mutable std::string m_text;
    mutable std::string m_file, m_func;
    mutable unsigned m_line;
  };

 

  namespace {
    void do_log(const Exception &) {
    }
    
  }

  template <typename T>
  const T & InitException(const T & e, const std::string & file, int line = 0, const std::string func = "") {
    e.SetLocation(file, line, func);
    do_log(e); 
    return e;
  }

  // Some useful predefined exceptions
  EUDAQINPUT_EXCEPTION(FileNotFoundException);
  EUDAQINPUT_EXCEPTION(FileExistsException);
  EUDAQINPUT_EXCEPTION(FileNotWritableException);
  EUDAQINPUT_EXCEPTION(FileReadException);
  EUDAQINPUT_EXCEPTION(FileWriteException);
  EUDAQINPUT_EXCEPTION(FileFormatException);
  EUDAQINPUT_EXCEPTION(CommunicationException);
  EUDAQINPUT_EXCEPTIONX(BusError, CommunicationException);

}

#endif // EUDAQINPUT_INCLUDED_Exception
