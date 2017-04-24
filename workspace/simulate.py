#!/usr/bin/python3

import tbsw

if __name__ == '__main__':
  
  steerfiles = ''
  path = ['reco.xml']
  ofile = 'mcdata.slcio'
  caltag = 'default'
  
  try:
    opts, args = getopt.getopt(sys.argv[1:],"hs:x:o:c:",["steerfiles=", "xmlfile=", "ofile=", "caltag="])
  except getopt.GetoptError:    
    sys.exit(2)
  
  for opt, arg in opts:
    if opt == '-h':
      print ('simulate.py -s <steerfiles>  -x <xmlfile> -o <outputfile> -c <caltag>')
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
  
  if xmlfile == '':
    print ('missing option: -x steerfile.xml')
    sys.exit(2)  

  tbsw.simulate(steerfiles, path, caltag, ofile)

