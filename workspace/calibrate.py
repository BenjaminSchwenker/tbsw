#!/usr/bin/python3

import tbsw
import getopt 
import sys

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
  
  if steerfiles == '':
    print ('missing option: -s path/to/steerfiles')
    sys.exit()  
   
  if ifile == '':
    print ('missing option: -i path/to/inputfilename')
    sys.exit()  
  
  tbsw.calibrate(steerfiles=steerfiles, path=path, caltag=caltag, ifile=ifile) 

