from tbsw import *
import xml.etree.ElementTree
import getopt, sys


steerfiles = 'steering-files/x0-sim/'
caltag = 'test'  
filename = 'mc.slcio'


try:
  opts, args = getopt.getopt(sys.argv[1:],"hs:f:c:",["steerfiles=", "filename=", "caltag=",])
except getopt.GetoptError:    
  print ('For help, please execute: scan.py -h ')
  sys.exit()
  
for opt, arg in opts:
  if opt in ("-h"):
    print ('scan.py -s <steerfiles>  -f <filename> -c <caltag>')
    sys.exit()
  elif opt in ("-s", "--steerfiles"):
    steerfiles = arg
  elif opt in ("-f", "--filename"):
    filename = arg
  elif opt in ("-c", "--caltag"):
    caltag = arg


calpath = [ 
           'hotpixelkiller.xml' , 
           'cluster-calibration-mc.xml',
           'correlator.xml' ,
           'kalmanalign-iteration-1.xml',
           'kalmanalign-iteration-2.xml',
           'kalmanalign-iteration-2.xml',
           'kalmanalign-iteration-2.xml',
           'telescope-dqm.xml',
           #'cluster-calibration-tb.xml',
         ]

name = os.path.splitext(os.path.basename(filename))[0] + '-' + caltag  

softscales = [ 6, ]


for scale in softscales:
  
  SimObj = Simulation(steerfiles=steerfiles, name=name + '-sim' )
  SimObj.simulate(path=['simulation.xml'], ofile=filename, caltag=None)  
   
  CalObj = Calibration(steerfiles=steerfiles, name=name + '-cal')
  
  # set the desired scale
  xmlfile = CalObj.get_name('cluster-calibration-mc.xml')
  tree = xml.etree.ElementTree.parse(xmlfile)
  root = tree.getroot() 
  
  for proc in root.findall('processor'):
    if proc.get('name') == 'M26ClusterDBCreator':
      for param in proc.findall('parameter'):
        if param.get('name') == 'SoftScale':
          param.set('value', str(scale)) 
  
  tree.write(xmlfile)    
  
  CalObj.calibrate(path=calpath,ifile=filename,caltag=caltag)  
   
  #RecObj = Reconstruction(steerfiles=steerfiles, name=name + '-reco' )
  #RecObj.reconstruct(path=['reco.xml'],ifile=filename,caltag=caltag) 
  
