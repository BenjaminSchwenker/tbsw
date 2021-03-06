#ifndef EUDAQINPUT_INCLUDED_BufferSerializer
#define EUDAQINPUT_INCLUDED_BufferSerializer

#include <deque>
#include <algorithm>
#include <iostream>
#include "Serializer.hh"
#include "Exception.hh"

namespace eudaqinput {

  class BufferSerializer : public Serializer, public Deserializer, public Serializable {
  public:
    BufferSerializer() : m_offset(0) {}
    template<typename InIt>
    BufferSerializer(InIt first, InIt last) : m_data(first, last), m_offset(0) {
    }
    BufferSerializer(Deserializer &);
    void clear() { m_data.clear(); m_offset = 0; }
    const unsigned char & operator [] (size_t i) const { return m_data[i]; }
    size_t size() const { return m_data.size(); }
    virtual bool HasData() { return m_data.size() != 0; }
    virtual void Serialize(Serializer &) const;
  private:
    virtual void Serialize(const unsigned char * data, size_t len);
    virtual void Deserialize(unsigned char * data, size_t len);
    std::vector<unsigned char> m_data;
    size_t m_offset;
  };

}

#endif // EUDAQINPUT_INCLUDED_BufferSerializer
