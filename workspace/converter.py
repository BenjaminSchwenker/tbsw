#!/usr/bin/python3

import os
import shutil
import subprocess
import glob
import sys, getopt 
from multiprocessing import Pool


def convert_run(params):
  
  # unpack arguments
  (rawfile,outputdir) = params
  
  # remember current working dir 
  fullpath = os.getcwd() 
  
  # make sure we have gets an absolute path
  if not os.path.isabs(outputdir): 
    outputdir = os.path.join(fullpath, outputdir)
     
  # make sure we have gets an absolute path
  if not os.path.isabs(rawfile): 
    rawfile = os.path.join(fullpath, rawfile)
    
  # can only run it here (FIXME)  
  os.chdir(os.environ['EUDAQ']+'/bin')
     
  # run converter
  subprocess.call('./Converter.exe ' + rawfile, shell=True)
  
  # get runtag
  runtag = os.path.splitext(os.path.basename(rawfile))[0]
    
  # move lcio file  
  outfile = os.path.join(outputdir,runtag+'.slcio')
  shutil.move(runtag+'.slcio', outfile)  
  
  # go back to initial dir  
  os.chdir(fullpath)


if __name__ == '__main__':
  
  inputfiles=[]
  outputdir=''
  
  try:
    opts, args = getopt.getopt(sys.argv[1:],"hi:o:",["ifile=","odir="])
  except getopt.GetoptError:
    print ('converter.py -i <inputfile>  -o <outputdir>')
    sys.exit(2)
  for opt, arg in opts:
    if opt == '-h':
      print ('converter.py -i <inputfile>  -o <outputdir>')
      sys.exit()
    elif opt in ("-i", "--ifile"):
      inputfiles = glob.glob(arg)
    elif opt in ("-o", "--odir"):
      outputdir = arg
    
  if inputfiles == []:
    print ('missing option: inputfile')
    sys.exit(2)  
  
  if outputdir == '':
    print ('missing option: outputdir')
    sys.exit(2)  

  # make sure outputdir exists
  if not os.path.isdir(outputdir):
    os.mkdir(outputdir)
  
  # pack params for multiprocessing
  params = [(rawfile,outputdir) for rawfile in inputfiles]
  
  p = Pool()
  p.map(convert_run, params )	
   
  
