#!/usr/bin/python3

import tbsw
import getopt 
import sys


if __name__ == '__main__':
  
  steerfiles = ''
  ifile = ''
  caltag = 'default'
  path = []

  try:
    opts, args = getopt.getopt(sys.argv[1:],"hs:x:i:c:",["steerfiles=", "xmlfile=", "ifile=", "caltag="])
  except getopt.GetoptError: 
    print ('For help, please execute: reco.py -h ')   
    sys.exit()   

  for opt, arg in opts:
    if opt == '-h':
      print ('reco.py -s <steerfiles>  -x <xmlfile> -i <ifile> -c <caltag>')
      sys.exit()
    elif opt in ("-s", "--steerfiles"):
      steerfiles = arg
    elif opt in ("-x", "--xmlfile"):
      path = [arg]
    elif opt in ("-i", "--ifile"):
      ifile = arg
    elif opt in ("-c", "--caltag"):
      caltag = arg
  
  if steerfiles == '':
    print ('missing option: -s path/to/steerfiles')
    sys.exit()  
  
  if path == []:
    print ('missing option: -x steerfile.xml')
    sys.exit()   

  if ifile == '':
    print ('missing option: -i path/to/inputfilename')
    sys.exit()  
   
  tbsw.reconstruct(steerfiles=steerfiles, path=path, caltag=caltag, ifile=ifile)
