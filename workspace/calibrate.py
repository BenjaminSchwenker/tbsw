#!/usr/bin/python3

from tbsw import *
import getopt, sys

if __name__ == '__main__':
  
  steerfiles = ''
  ifile = ''
  caltag = 'default'
  path = [ 
           'hotpixelkiller.xml' , 
           #'cluster-calibration-mc.xml',
           'correlator.xml' ,
           'kalmanalign-iteration-1.xml',
           'kalmanalign-iteration-2.xml',
           'kalmanalign-iteration-2.xml',
           'kalmanalign-iteration-2.xml',
           'telescope-dqm.xml',
           #'cluster-calibration-tb.xml',
         ]
  
  try:
    opts, args = getopt.getopt(sys.argv[1:],"hs:i:c:",["steerfiles=", "ifile=", "caltag="])
  except getopt.GetoptError: 
    print ('For help, please execute: reco.py -h ')   
    sys.exit()   
  
  for opt, arg in opts:
    if opt == '-h':
      print ('reco.py -s <steerfiles>  -i <ifile> -c <caltag>')
      sys.exit()
    elif opt in ("-s", "--steerfiles"):
      steerfiles = arg
    elif opt in ("-i", "--ifile"):
      ifile = arg
    elif opt in ("-c", "--caltag"):
      caltag = arg
   
  name = os.path.splitext(os.path.basename(ifile))[0] + '-' + caltag + '-cal' 
  CalObj = Calibration(steerfiles=steerfiles, name=name )
  CalObj.calibrate(path=path,ifile=ifile,caltag=caltag)  
