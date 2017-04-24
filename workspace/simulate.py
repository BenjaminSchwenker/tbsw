#!/usr/bin/python3

import utility

if __name__ == '__main__':
  
  gearfile = ''
  xmlfile = ''
  ofile = 'mcdata.slcio'
  caltag='default'
  
  try:
    opts, args = getopt.getopt(sys.argv[1:],"hg:x:o:c:",["gearfile=", "xmlfile=", "ofile=", "caltag="])
  except getopt.GetoptError:
    
    sys.exit(2)
  for opt, arg in opts:
    if opt == '-h':
      print ('simulate.py -g <gearfile>  -x <xmlfile> -o <outputfile> -c <caltag>')
      sys.exit()
    elif opt in ("-g", "--gearfile"):
      gearfile = arg
    elif opt in ("-x", "--xmlfile"):
      xmlfile = arg
    elif opt in ("-o", "--ofile"):
      ofile = arg
    elif opt in ("-c", "--caltag"):
      caltag = arg
  
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

  # create tmp dir
  if not os.path.isdir(fullpath+'/tmp-runs'):
    os.mkdir(fullpath+'/tmp-runs')
  tmpdir = os.path.join(fullpath+'/tmp-runs',runtag)  
  tmpdir = tmpdir + '-' + caltag + '-sim' 
    
  if os.path.isdir(tmpdir):
    shutil.rmtree(tmpdir)
       
  os.mkdir(tmpdir)
  
  # populate tmpdir with all needed files 
  shutil.copy(gearfile, tmpdir+'/gear.xml')
  shutil.copy(xmlfile, tmpdir+'/sim.xml')
  
  # check that calibration files exist
  caldir = caldir = fullpath+'/cal-files/'+caltag   
  if not os.path.isdir(caldir):
    print ('[Print] warning: caltag not found')
    return  
  else: 
    # copy calibration files 
    shutil.copytree(caldir,tmpdir+'/cal-files')  
  
  # run simulation in tmp dir  
  os.chdir(tmpdir)
  
  subprocess.call('/$MARLIN/bin/Marlin sim.xml > log-sim.txt 2>&1', shell=True)
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
                  
  
  
     	           
                


 	


