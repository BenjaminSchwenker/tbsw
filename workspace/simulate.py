#!/usr/bin/python3

import tbsw
import getopt 
import sys

if __name__ == '__main__':
  
  steerfiles = ''
  ofile = 'mcdata.slcio'
  caltag = 'default'
  path = []
  
  try:
    opts, args = getopt.getopt(sys.argv[1:],"hs:x:o:c:",["steerfiles=", "xmlfile=", "ofile=", "caltag="])
  except getopt.GetoptError:    
    print ('For help, please execute: simulate.py -h ')
    sys.exit()
  
  for opt, arg in opts:
    if opt == '-h':
      print ('simulate.py -s <steerfiles>  -x <xmlfile> -o <ofile> -c <caltag>')
      sys.exit()
    elif opt in ("-s", "--steerfiles"):
      steerfiles = arg
    elif opt in ("-x", "--xmlfile"):
      path = [arg]
    elif opt in ("-o", "--ofile"):
      ofile = arg
    elif opt in ("-c", "--caltag"):
      caltag = arg
  
  if steerfiles == '':
    print ('missing option: -s path/to/steerfiles')
    sys.exit(2)  
  
  if path == []:
    print ('missing option: -x steerfile.xml')
    sys.exit(2)  

  tbsw.simulate(steerfiles=steerfiles, path=path, caltag=caltag, ofile=ofile)

