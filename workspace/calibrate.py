#!/usr/bin/python3

import os
import shutil
import subprocess
import sys, getopt 
import glob

if __name__ == '__main__':
  
  rawfile = ''
  xmlpath = ''
  caltag='dummy'
  
  try:
    opts, args = getopt.getopt(sys.argv[1:],"hi:x:c:",["ifile=","xmlfile=","caltag="])
  except getopt.GetoptError:
    print ('calibrate.py -i <inputfile>  -x <xmlfile> -c <caltag>')
    sys.exit(2)
  for opt, arg in opts:
    if opt == '-h':
      print ('calibrate.py -i <inputfile> -x <xmlfile> -c <caltag>')
      sys.exit()
    elif opt in ("-i", "--ifile"):
      rawfile = arg
    elif opt in ("-x", "--xmlfile"):
      xmlpath = arg
    elif opt in ("-c", "--caltag"):
      caltag = arg
  
  if rawfile == '':
    print ('missing option: -i inputfile')
    sys.exit(2)  
  
  if xmlpath == '':
    print ('missing option: -x path/to/xmlfiles')
    sys.exit(2)  
    
  # remember current working dir 
  fullpath = os.getcwd() 
  
  # get runtag
  runtag = os.path.splitext(os.path.basename(rawfile))[0]
  
  # create tmp dir
  if not os.path.isdir(fullpath+'/tmp-runs'):
    os.mkdir(fullpath+'/tmp-runs')
  
  tmpdir = os.path.join(fullpath+'/tmp-runs',runtag)
  tmpdir = tmpdir + '-' + caltag  
  
  if os.path.isdir(tmpdir):
    shutil.rmtree(tmpdir)
       
  os.mkdir(tmpdir)
  
  # copy steer files 
  for name in os.listdir(xmlpath):
    srcname = os.path.join(xmlpath, name)
    shutil.copy(srcname, tmpdir)
  
  # copy rawdata 
  shutil.copy(rawfile, tmpdir+'/tmp-rawdata.slcio')
  
  # run calibration in tmp dir  
  os.chdir(tmpdir)
  
  subprocess.call('/$MARLIN/bin/Marlin hotpixelkiller.xml > log-hotpixel.txt 2>&1', shell=True)
  print ('[Print] HotPixelKiller done ...')
  		  
  subprocess.call('/$MARLIN/bin/Marlin correlator.xml > log-correlator.txt 2>&1', shell=True)    
  print ('[Print] Correlator done ...')           
                  		              
  subprocess.call('/$MARLIN/bin/Marlin kalmanalign-it1.xml > log-align-it1.txt 2>&1', shell=True)     
  print ('[Print] Alignment first iteration done ...')               
               	               
  subprocess.call('/$MARLIN/bin/Marlin kalmanalign-it2.xml > log-align-it2.txt 2>&1', shell=True)   
  print ('[Print] Alignment second iteration done ...')          
		               
  subprocess.call('/$MARLIN/bin/Marlin kalmanalign-it2.xml > log-align-it3.txt 2>&1', shell=True)  
  print ('[Print] Alignment third iteration done ...')     
		               
  subprocess.call('/$MARLIN/bin/Marlin kalmanalign-it2.xml > log-align-it4.txt 2>&1', shell=True)   
  print ('[Print] Alignment fourth iteration done ...')   
	               
  subprocess.call('/$MARLIN/bin/Marlin telescope-dqm.xml > log-dqm.txt 2>&1', shell=True)     
  print ('[Print] TelescopeDQM done.')   	

  # clean up tmp files 
  for tmpfile in glob.glob('tmp*'):
    os.remove(tmpfile)
     
  # create new calibration tag    
  if not os.path.isdir(fullpath+'/cal-files'):
    os.mkdir(fullpath+'/cal-files')   
  caldir = fullpath+'/cal-files/'+caltag  
  
  if os.path.isdir(caldir):
    shutil.rmtree(caldir)
       
  os.mkdir(caldir)		

  # save all files with calibration data to common folder   
  for dbfile in glob.glob('*.slcio'): 
    shutil.copy(dbfile, os.path.join(caldir,dbfile))  
  # treat the gear file as a calibration file as well
  shutil.copy('gear.xml', os.path.join(caldir,'gear.xml'))  
  	                       
  os.chdir(fullpath)
                  
  print ('[Print] All calibration files written to ' + caldir)   
               
		           
                


 	


