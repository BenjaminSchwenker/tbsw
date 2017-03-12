#!/usr/bin/python3

import os
import shutil
import subprocess
import sys, getopt 
import glob

if __name__ == '__main__':
  
  gearfile = ''
  xmlfile = ''
  ofile = 'mcdata.slcio'
  dbfile = ''
  
  try:
    opts, args = getopt.getopt(sys.argv[1:],"hg:x:o:d:",["gearfile=", "xmlfile=", "ofile=", "dbfile="])
  except getopt.GetoptError:
    
    sys.exit(2)
  for opt, arg in opts:
    if opt == '-h':
      print ('simulate.py -g <gearfile>  -x <xmlfile> -o <outputfile> -d <alignDB>')
      sys.exit()
    elif opt in ("-g", "--gearfile"):
      gearfile = arg
    elif opt in ("-x", "--xmlfile"):
      xmlfile = arg
    elif opt in ("-o", "--ofile"):
      ofile = arg
    elif opt in ("-d", "--dbfile"):
      dbfile = arg
  
  if gearfile == '':
    print ('missing option: -g path/to/gearfile.xml')
    sys.exit(2)  
  
  if xmlfile == '':
    print ('missing option: -x path/to/steerfile.xml')
    sys.exit(2)  
    
  # remember current working dir 
  fullpath = os.getcwd() 
   
  # get runtag
  runtag = os.path.splitext(os.path.basename(ofile))[0]
   
  # create tmp dir if it does not exist
  tmpdir = os.path.join(fullpath,'tmp-runs')
  if not os.path.isdir(tmpdir):
    os.mkdir(tmpdir)
          
  tmpdir = os.path.join(tmpdir,runtag)  
  
  if os.path.isdir(tmpdir):
    shutil.rmtree(tmpdir)
       
  os.mkdir(tmpdir)
  
  # populate tmpdir with all needed files 
  shutil.copy(gearfile, tmpdir+'/gear.xml')
  shutil.copy(xmlfile, tmpdir+'/sim.xml')
  if os.path.isfile(dbfile): 
    shutil.copy(dbfile, tmpdir+'/eudet-alignmentDB.slcio')  
  
  # run simulation in tmp dir  
  os.chdir(tmpdir)
  
  #subprocess.call('/$MARLIN/bin/Marlin sim.xml > log-sim.txt 2>&1', shell=True)
  subprocess.call('/$MARLIN/bin/Marlin sim.xml', shell=True)
  print ('[Print] Simulation done ...')
  
  def move_over(src, dest):
    if os.path.exists(dest):
      os.remove(dest)
    if os.path.isfile(src):
      shutil.move(src, dest)
  
  src = 'outputfile.slcio'
  dest = os.path.join(fullpath, ofile)	
  move_over(src, dest)      	
  	                       
  os.chdir(fullpath)
                  
  
  
     	           
                


 	


