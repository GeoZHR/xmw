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
  global sz,sl,sc
  global nz,nl,nc
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
    f1,f2,f3 = 0.0,0.0,0.0
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
    #n3 number of log curves
    #n2 number of wells
    #n1 depth samples of the longest well
    #s1,s2,s3 = sz,sl,sc
  else:
    print "unrecognized subset:",name
    System.exit

def getSamplings():
  return s1,s2,s3

def getWellSamplings():
  return sz,sl,sc

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

def writeLogsToArray(basename,logs):
  dz = 0
  nl = len(logs) # number of well logs
  fzs = zerofloat(nl)
  lzs = zerofloat(nl)
  for il, log in enumerate(logs):
    z  = log.z
    nz = len(z)
    fzs[il] = z[0]
    lzs[il] = z[nz-1]
    dz = z[1]-z[0]
    print z[0]
  fz = min(fzs)
  nz = round((max(lzs)-fz)/dz)+1
  nullValue = -999.25
  nc = 2  # number of log curves
  wa = fillfloat(nullValue,nz,nl,nc)
  for il, log in enumerate(logs):
    d = log.d
    v = log.v
    n = len(d)
    j = round((log.z[0]-fz)/dz)
    copy(n,0,d,j,wa[0][il]) 
    copy(n,0,v,j,wa[1][il]) 
  fileName = welllogDir+basename+".dat"
  aos = ArrayOutputStream(fileName)
  aos.writeFloats(wa)
  aos.close()
  dz = dz*0.0003048 # ft to km
  fz = fz*0.0003048 # ft to km
  return sz,sl,sc,wa
def logsToArray(basename,logs):
  fk = 0.0003048
  dz = 0.5*fk
  fzs,lzs =[],[]
  curves=["den","vel"]
  for curve in curves:
    for log in wldata.getLogsWith(curve):
      f,z,y,x = log.getSamples(curve)
      nz = len(z)
      fzs.append(z[0])
      lzs.append(z[nz-1])
  fz = min(fzs)
  nz = round((max(zl)-fz)/dz)+1
  nullValue = -999.25
  nl = 12 # number of well logs
  nc = 2  # number of log curves
  wa = fillfloat(nullValue,nz,nl,nc)
  for ic, curve in enumerate(curves):
    logs = wldata.getLogsWith(curve)
    for il, log in enumerate(logs):
      f,z,y,x = log.getSamples(curve)
      n = len(f)
      j = round((z[0]-fz)/dz)
      copy(n,0,f,j,wa[ic][il]) 

  return wa

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
