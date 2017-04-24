#!/usr/bin/python3

"""
Utility code to use TBSW in pyhton scripts 

:author: benjamin.schwenker@phys.uni-goettinge.de  
"""

import os
import shutil
import subprocess
import sys, getopt 
import glob


def simulate(steerfiles, marlinpath, caltag, ofile):
  """
  Creates a lcio file containing N events. The events contain collections for 
  simulated tracks measured signals. 
  
  :@steerfile:  name of folder containing all steering files  
  :@marlinpath: list containing Marlin xml files that will be executed 
  :@caltag:     name of calibration tag (optional)
  :@ofile:      name of output lcio file
  
  :author: benjamin.schwenker@phys.uni-goettinge.de  
  """

  return None


def reconstruct(steerfiles, marlinpath, caltag, ifile):
  """
  Reconstructs an input file with raw data using a caltag for calibration. 
  
  :@steerfile:  name of folder containing all steering files 
  :@marlinpath: list containing Marlin xml files that will be executed 
  :@caltag:     name of calibration tag (optional)
  :@ifile:      name of input file with raw data
  
  :author: benjamin.schwenker@phys.uni-goettinge.de  
  """   

  return None


def calibrate(steerfiles, marlinpath, caltag, ifile):
  """
  Calibrate raw data from a test beam experiment. 
  
  :@steerfile:  name of folder containing all steering files 
  :@marlinpath: list containing Marlin xml files that will be executed 
  :@caltag:     name of calibration tag (optional)
  :@ifile:      name of input file with raw data
  
  :author: benjamin.schwenker@phys.uni-goettinge.de  
  """     

  return None

