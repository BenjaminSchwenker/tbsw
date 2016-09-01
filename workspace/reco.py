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
  (rawfile,xmlfile) = params

  # remember current working dir 
  fullpath = os.getcwd() 

  store_lcio = False
  
  # get runtag
  runtag = os.path.splitext(os.path.basename(rawfile))[0]
                    
  print ('[Print] Starting to process file ' + runtag + ' ...')

  # create tmp dir
  tmpdir = os.path.join(fullpath+'/tmp-runs',runtag)  
  
  if os.path.isdir(tmpdir):
    shutil.rmtree(tmpdir)
       
  os.mkdir(tmpdir)
  
  # and populat it with needed files 
  shutil.copy(xmlfile, tmpdir+'/reco.xml')
  shutil.copy(rawfile, tmpdir+'/tmp-rawdata.slcio')
  shutil.copytree('cal-files',tmpdir+'/cal-files')
  shutil.copytree('gear-files',tmpdir+'/gear-files')
  
  # run reco in tmp dir 
  os.chdir(tmpdir)
                       
  subprocess.call('/$MARLIN/bin/Marlin reco.xml > /dev/null 2>&1', shell=True)
  

  # clean up tmp files 
  for tmpfile in glob.glob('tmp*'):
    os.remove(tmpfile)

  # store dqm files 
  for dqmfile in glob.glob('*.root'): 
    name = os.path.splitext(os.path.basename(dqmfile))[0]
    shutil.move(dqmfile, fullpath+'/dqm-files/'+name+'-'+runtag+'.root')  
  
  # store processed lcio files
  if store_lcio == True: 
    for slciofile in glob.glob('*.slcio'):    
      name = os.path.splitext(os.path.basename(slciofile))[0]
      shutil.move(slciofile, fullpath+'/track-files/'+name+'-'+runtag+'.slcio')  
  
  os.chdir(fullpath)
  shutil.rmtree(tmpdir)

  print ('[Print] Processing file ' + runtag + ' done!')
	                        

if __name__ == '__main__':
  
  inputfiles = []
  xmlfile = ''
    
  try:
    opts, args = getopt.getopt(sys.argv[1:],"hi:x:",["ifile=","xmlfile="])
  except getopt.GetoptError:
    print ('reco.py -i <inputfiles>  -x <xmlfile>')
    sys.exit(2)
  for opt, arg in opts:
    if opt == '-h':
      print ('reco.py -i <inputfiles> -x <xmlfile>')
      sys.exit()
    elif opt in ("-i", "--ifile"):
      inputfiles = glob.glob(arg)
    elif opt in ("-x", "--xmlfile"):
      xmlfile = arg
  
  if inputfiles == []:
    print ('missing option: -i inputfile')
    sys.exit(2)  
  
  if xmlfile == '':
    print ('missing option: -x xmlfile')
    sys.exit(2)  
  
  # pack params for multiprocessing
  params = [(rawfile,xmlfile) for rawfile in inputfiles]
  
  p = Pool()
  p.map(reco_run, params )	


 
