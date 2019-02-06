#ifndef EUDAQINPUT_INCLUDED_Serializable
#define EUDAQINPUT_INCLUDED_Serializable


namespace eudaqinput {

  class Serializer;

  class Serializable {
  public:
    virtual void Serialize(Serializer &) const = 0;
    virtual ~Serializable() {}

  };

}

#endif // EUDAQINPUT_INCLUDED_Serializable
