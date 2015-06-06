"""
Jython utilities for Teapot Dome data set.
Author: Xinming Wu, Colorado School of Mines
Version: 2015.06.01
"""
from common import *

#############################################################################
# Internal constants

_datdir = "../../../data/seis/swt/csm/"
_welldir = "../../../data/seis/swt/csm/welllogs/"
#############################################################################
# Setup

def setupForSubset(name):
  """
  Setup for a specified directory includes:
    seismic directory
    samplings s1,s2,s3
  Example: setupForSubset("s1")
  """
  global seismicDir
  global welllogDir
  global s1,s2,s3
  global n1,n2,n3
  if name=="seismict":
    print "setupForSubset: seismict"
    seismicDir = _datdir+"seismict/"
    n1,n2,n3 = 1501,357,161
    d1,d2,d3 = 1.0,1.0,1.0 
    d1,d2,d3 = 0.002,0.025,0.025 # (s,km,km)
    f1,f2,f3 = 0.000,0.000,0.000
    s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n3,d3,f3)

  elif name=="subt":
    print "setupForSubset: subt"
    seismicDir = _datdir+"seismict/subt/"
    n1,n2,n3 = 1025,240,80
    d1,d2,d3 = 1.0,1.0,1.0 
    d1,d2,d3 = 0.002,0.025,0.025 # (s,km,km)
    f1,f2,f3 = 0.000,d2*29,d3*46
    s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n3,d3,f3)

  elif name=="seismicz":
    print "setupForSubset: seismicz"
    seismicDir = _datdir+"seismicz/"
    n1,n2,n3 = 2762,357,161
    d1,d2,d3 = 1.0,1.0,1.0 
    #d1,d2,d3 = 0.002,0.025,0.025 # (km,km,km)
    f1,f2,f3 = 0.000,0.000,0.000
    s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n3,d3,f3)

  elif name=="subw":
    print "setupForSubset: subset of well logs"
    welllogDir = _welldir+"subw/"
    #n3 number of wells
    #n2 number of log curves
    #n1 depth samples of the longest well
    n1,n2,n3 = 10558,2,12
    d1,d2,d3 = 0.000152,1.0,1.0 
    f1,f2,f3 = 0.246278,1.0,1.0
    s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n3,d3,f3)
  else:
    print "unrecognized subset:",name
    System.exit

def getSamplings():
  return s1,s2,s3

def getSeismicDir():
  return seismicDir

#############################################################################
# read/write images

def readImage(basename):
  """ 
  Reads an image from a file with specified basename
  """
  fileName = seismicDir+basename+".dat"
  image = zerofloat(n1,n2,n3)
  ais = ArrayInputStream(fileName)
  ais.readFloats(image)
  ais.close()
  return image


def writeImage(basename,image):
  """ 
  Writes an image to a file with specified basename
  """
  fileName = seismicDir+basename+".dat"
  aos = ArrayOutputStream(fileName)
  aos.writeFloats(image)
  aos.close()
  return image

#############################################################################
# read/write logs
def getLogs():
  wlName = _welldir+"tpwa.dat"
  wldata = WellLog.Data.readBinary(wlName)
  logs = wldata.getAll()
  welllogs=[]
  for log in logs:
    if(log!=None and log.getCurve("den")!=None and log.getCurve("vel")!=None):
      model = SynSeis.getModel(log)
      if(model!=None and model.nz>2000 and model.x1<0.8):
        welllogs.append(log)
  return welllogs

def writeLogs(basename,logs):
  wldata = WellLog.Data()
  for log in logs:
    wldata.add(log)
  writeLogData(basename,wldata)

def writeLogData(basename,wldata):
  fileName = welllogDir+basename+".dat"
  wldata.writeBinary(fileName)

def readLogData(basename):
  fileName = welllogDir+basename+".dat"
  return WellLog.Data.readBinary(fileName)

def writeLogDataToArray(wldata):
  dz = 0.000152
  zf,zl =[],[]
  curves=["den","vel"]
  for curve in curves:
    for log in wldata.getLogsWith(curve):
      f,z,y,x = log.getSamples(curve)
      nz = len(z)
      zf.append(z[0])
      zl.append(z[nz-1])
  nz = round((max(zl)-min(zf))/dz)
  fz = round(min(zf)*1000000.0)/1000000.0
  nullValue = -999.25
  nl = 12 # number of well logs
  nc = 2  # number of log curves
  wa = fillfloat(nullValue,nz,nc,nl)
  for ic, curve in enumerate(curves):
    logs = wldata.getLogsWith(curve)
    for il, log in enumerate(logs):
      f,z,y,x = log.getSamples(curve)
      n = len(f)
      j = round((z[0]-fz)/dz)
      copy(n,0,f,j,wa[il][ic]) 
  fileName = welllogDir+"wx"+".dat"
  aos = ArrayOutputStream(fileName)
  aos.writeFloats(wa)
  aos.close()
  print nz
  print dz
  print fz

def readLogArray(basename):
  """ 
  Reads an array of well logs from a file with specified basename
  """
  fileName = welllogDir+basename+".dat"
  image = zerofloat(n1,n2,n3)
  ais = ArrayInputStream(fileName)
  ais.readFloats(image)
  ais.close()
  return image


from org.python.util import PythonObjectInputStream
def readObject(name):
  fis = FileInputStream(seismicDir+name+".dat")
  ois = PythonObjectInputStream(fis)
  obj = ois.readObject()
  ois.close()
  return obj
def writeObject(name,obj):
  fos = FileOutputStream(seismicDir+name+".dat")
  oos = ObjectOutputStream(fos)
  oos.writeObject(obj)
  oos.close()
