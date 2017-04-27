#!/usr/bin/python3

from tbsw import *

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
  
  name = os.path.splitext(os.path.basename(ifile))[0] + '-' + caltag + '-reco' 
  RecObj = Reconstruction(steerfiles=steerfiles, name=name )
  RecObj.reconstruct(path=path,ifile=ifile,caltag=caltag) 
