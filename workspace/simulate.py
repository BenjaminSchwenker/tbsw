#!/usr/bin/python3

from tbsw import *

if __name__ == '__main__':
  
  steerfiles = ''
  ofile = 'mcdata.slcio'
  caltag = 'default'
  path = []
  
  try:
    opts, args = getopt.getopt(sys.argv[1:],"hs:x:o:c:",["steerfiles=", "xmlfile=", "ofile=", "caltag="])
  except getopt.GetoptError:    
    print ('For help, please execute: scan.py -h ')
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
     
    name = os.path.splitext(os.path.basename(ofile))[0] + '-' + caltag + '-sim' 
    SimObj = Simulation(steerfiles=steerfiles, name=name )
    SimObj.simulate(path=path, ofile=ofile, caltag=caltag)  
    
  

