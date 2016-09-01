#ifndef DEPFET_S3ACONVERTER_H
#define DEPFET_S3ACONVERTER_H

#include "RawData.h"
#include "ADCValues.h"

namespace depfet {

  struct S3AConverter {
    void operator()(const RawData& rawData, ADCValues& adcValues);
  };

}
#endif
