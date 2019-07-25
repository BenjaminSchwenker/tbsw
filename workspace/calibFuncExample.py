"""
  Script generating an example GainCalibrationDB file 
  with the structure as it is used for the 
  PixelChargeCalibrator processor.

  Calibration function and its parameters are stored in 
  folders in the TFile. The folder names are build by:
  d + sensorID (i.e. d1, d21, ...)
  (d for detector)

  Sensors which do not have calibrations provided 
  do not get calibrated.

  The function parameters are  binned in pixel of the sensor (column, row).

  author: Helge Christoph Beck
  email: helge-christoph.beck@phys.uni-goettingen.de
"""

from ROOT import TFile, TF1, TH2F

if __name__ == '__main__':
  
  # parameters to provide
  cols = 80
  rows = 336

  sensorIDList = [21, 22]
  
  baseCalibFuncName = "calibFunc"
  baseCalibParaName = "para"

  gainCalibrationFileName = "calibrationDBFile.root"

  # generating file structure
  foutfile = TFile(gainCalibrationFileName, "RECREATE")
  
  for sensorID in sensorIDList:
    foutfile.cd()
    detDir = foutfile.mkdir("d" + str(sensorID)) # creating folder for sensor
    detDir.cd()
    funcName = baseCalibFuncName # function name is the same for all sensors in the file, differentiated by the folders
    fucalib = TF1(funcName, "pol2", -1.0, 100.0) # creating calibration function, use whatever you need as function and range
    # setting base parameters could be needed for non standard functions. Consult the TF1 reference.
    fucalib.Draw()
    fucalib.Write()

    nparFunc = fucalib.GetNpar()

    for par in range(nparFunc):
      histoParaName = baseCalibParaName + "_" + str(par) # same name for all sensors but in different folder
      hpara = TH2F(histoParaName, "", cols, 0, cols, rows, 0, rows) # creating histogram for parameters

      weight = par + 1 # get your parameters from a dedicated calibration probably
      for c in range(cols):
        for r in range(rows):
          hpara.Fill(c, r, weight)
      hpara.Write()

  foutfile.Close()
