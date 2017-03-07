#!/usr/bin/python3

import sys, getopt 
import os
import shutil
import subprocess
import glob
from multiprocessing import Pool
from itertools import repeat


def reco_run(params):
  
  # unpack arguments
  (rawfile,xmlfile,caltag) = params

  # remember current working dir 
  fullpath = os.getcwd() 
   
  store_lcio = False
  
  # get runtag
  runtag = os.path.splitext(os.path.basename(rawfile))[0]
                    
  print ('[Print] Starting to process file ' + runtag + ' ...')

  # create tmp dir
  if not os.path.isdir(fullpath+'/tmp-runs'):
    os.mkdir(fullpath+'/tmp-runs')
  tmpdir = os.path.join(fullpath+'/tmp-runs',runtag)  
  
  if os.path.isdir(tmpdir):
    shutil.rmtree(tmpdir)
       
  os.mkdir(tmpdir)
  
  # copy lcio file that is to be processed
  shutil.copy(rawfile, tmpdir+'/inputfilename')
  # copy steering file 
  shutil.copy(xmlfile, tmpdir+'/reco.xml')
  
  # check that calibration files exist
  caldir = caldir = fullpath+'/cal-files/'+caltag   
  if not os.path.isdir(caldir):
    print ('[Print] warning: no calibration data found')
    return  
  else: 
    # copy calibration files 
    shutil.copytree(caldir,tmpdir+'/cal-files')
  
  # run reco in tmp dir 
  os.chdir(tmpdir)
                       
  #subprocess.call('/$MARLIN/bin/Marlin reco.xml > log.txt 2>&1', shell=True)
  subprocess.call('/$MARLIN/bin/Marlin reco.xml', shell=True)
  
  # clean up inputfile
  os.remove('inputfilename')

  # clean up tmp files 
  for tmpfile in glob.glob('tmp*'):
    os.remove(tmpfile)

  # store dqm files 
  if not os.path.isdir(fullpath+'/root-files'):
    os.mkdir(fullpath+'/root-files')
  
  for dqmfile in glob.glob('*.root'): 
    name = os.path.splitext(os.path.basename(dqmfile))[0]
    shutil.move(dqmfile, fullpath+'/root-files/'+name+'-'+runtag+'-'+cal-tag+'.root')  
  
  os.chdir(fullpath)
  #shutil.rmtree(tmpdir)

  print ('[Print] Processing file ' + runtag + ' done!')
	                        

if __name__ == '__main__':
  
  inputfiles = []
  xmlfile = ''
  caltag='dummy'
    
  try:
    opts, args = getopt.getopt(sys.argv[1:],"hi:x:c:",["ifile=","xmlfile=","caltag="])
  except getopt.GetoptError:
    print ('reco.py -i <inputfiles>  -x <xmlfile> -c <caltag>')
    print ('-c is optional and defaults to: ' + caltag )
    sys.exit(2)
  for opt, arg in opts:
    if opt == '-h':
      print ('reco.py -i <inputfiles> -x <xmlfile> -c <caltag>')
      print ('-c is optional and defaults to: ' + caltag )
      sys.exit()
    elif opt in ("-i", "--ifile"):
      inputfiles = glob.glob(arg)
    elif opt in ("-x", "--xmlfile"):
      xmlfile = arg
    elif opt in ("-c", "--caltag"):
      caltag = arg
  
  if inputfiles == []:
    print ('missing option: -i path/to/inputfiles')
    sys.exit(2)  
  
  if xmlfile == '':
    print ('missing option: -x path/to/steeringfilename.xml')
    sys.exit(2)  
  
  # pack params for multiprocessing
  params = [(rawfile,xmlfile,caltag) for rawfile in inputfiles]
  
  p = Pool()
  p.map(reco_run, params )	


 
