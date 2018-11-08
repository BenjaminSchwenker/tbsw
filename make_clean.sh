#!/bin/bash



############################################################
# Clean gear
rm -rf source/gear/GEARConfig.cmake
rm -rf source/gear/build/
rm -rf source/gear/include/
rm -rf source/gear/lib/
rm -rf source/gear/src/cpp/include/gearimpl
rm -rf source/gear/src/cpp/include/gearxml

############################################################
# Clean lcio
rm -rf source/lcio/LCIOConfig.cmake
rm -rf source/lcio/bin/
rm -rf source/lcio/build/
rm -rf source/lcio/include/
rm -rf source/lcio/lib/

############################################################
# Clean Marlin  
rm -rf source/Marlin/MarlinConfig.cmake
rm -rf source/Marlin/bin/
rm -rf source/Marlin/build/
rm -rf source/Marlin/lib/
rm -rf source/Marlin/streamlogConfig.cmake

############################################################
# Clean EudaqInput
rm -rf source/EudaqInput/EudaqInputConfig.cmake
rm -rf source/EudaqInput/build/

############################################################
# Clean TBTools  
rm -rf source/TBTools/TBToolsConfig.cmake
rm -rf source/TBTools/build/
rm -rf source/TBTools/source/include/AlignEventDict_rdict.pcm
rm -rf source/TBTools/source/src/AlignEventDict.C

############################################################
# Clean TBReco  
rm -rf source/TBReco/TBRecoConfig.cmake
rm -rf source/TBReco/build/

############################################################
# Remove init scripts
rm -rf init_tbsw.sh
rm -rf workspace/init_tbsw.sh
