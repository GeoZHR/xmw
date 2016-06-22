"""
Jython utilities for Teapot Dome data set.
Author: Xinming Wu, Colorado School of Mines
Version: 2015.06.01
"""
from common import *

#############################################################################
# Internal constants

_datdir = "../../../data/seis/cfd/"
_welldir = "../../../data/seis/cfd/wells/"
#############################################################################
# Setup


def setupForSubset(name):
  """
  Setup for a specified directory includes:
    seismic directory
    samplings s1,s2,s3
  Example: setupForSubset("s1")
  """
  # Location of w281
  _x2,_x3=242557.93,393670.74 #(ft)
  fk = 0.0003048 # 1 ft = 0.0003048 km
  _x2 = _x2*fk #1 ft = 0.0003048 km
  _x3 = _x3*fk #1 ft = 0.0003048 km
  _i2,_i3 = 165, 81 # inline and crossline position
  global _d2,_d3 # shifts to the seismic
  global seismicDir
  global welllogDir
  global s1,s2,s3
  global n1,n2,n3
  global sz,sl,sc
  global nz,nl,nc
  if name=="cfd2007":
    print "setupForSubset: cranfield 2007"
    seismicDir = _datdir+"cfd2007/"
    n1,n2,n3 = 950,243,222
    d1,d2,d3 = 0.002, 0.025146,  0.025146 # (s,km,km)
    f1,f2,f3 = 0.600,70.873811,115.947758 # (s,km,km)
    #d1,d2,d3 = 1.0,1.0,1.0 
    #f1,f2,f3 = 0.0,0.0,0.0 
    _d2 = f2+_i2*d2-_x2
    _d3 = f3+_i3*d3-_x3
    s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n3,d3,f3)
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

def readImageL(basename):
  """ 
  Reads an image from a file with specified basename
  """
  fileName = seismicDir+basename+".dat"
  image = zerofloat(n1,n2,n3)
  ais = ArrayInputStream(fileName,ByteOrder.LITTLE_ENDIAN)
  ais.readFloats(image)
  ais.close()
  return image


def readImage1(n1,basename):
  """ 
  Reads an image from a file with specified basename
  """
  fileName = seismicDir+basename+".dat"
  image = zerofloat(n1)
  ais = ArrayInputStream(fileName)
  ais.readFloats(image)
  ais.close()
  return image

def readImage2(n1,n2,basename):
  """ 
  Reads an image from a file with specified basename
  """
  fileName = seismicDir+basename+".dat"
  image = zerofloat(n1,n2)
  ais = ArrayInputStream(fileName,ByteOrder.LITTLE_ENDIAN)
  ais.readFloats(image)
  ais.close()
  return image

def readLog(n1,n2,fileName):
  """ 
  Reads an image from a file with specified basename
  """
  image = zerofloat(n1,n2)
  ais = ArrayInputStream(fileName)
  ais.readFloats(image)
  ais.close()
  return image



def writeImage(basename,image):
  """ 
  Writes an image to a file with specified basename
  """
  fileName = seismicDir+basename+".dat"
  aos = ArrayOutputStream(fileName,ByteOrder.LITTLE_ENDIAN)
  aos.writeFloats(image)
  aos.close()
  return image

#############################################################################
# read/write logs
def getLog242():
  n1,n2=8842,3
  #w271.dat  n1,n2=9123,3
  #w281.dat  n1,n2=8729,3
  fk = 0.0003048 # 1 ft = 0.0003048 km
  x2,x3=236264.72,395566.00 #(ft)
  x2 = x2*fk+_d2
  x3 = x3*fk+_d3
  wlName = _welldir+"w242.dat"
  ws = readLog(n1,n2,wlName)
  mul(ws[1],fk,ws[1])
  return x2,x3,ws
def getLog271():
  n1,n2=9123,3
  #w281.dat  n1,n2=8729,3
  fk = 0.0003048 # 1 ft = 0.0003048 km
  x2,x3=242565.70,393699.70 #(ft)
  x2 = x2*fk+_d2
  x3 = x3*fk+_d3
  wlName = _welldir+"w271.dat"
  ws = readLog(n1,n2,wlName)
  mul(ws[1],fk,ws[1])
  return x2,x3,ws
def getLog281():
  n1,n2=8729,3
  fk = 0.0003048 # 1 ft = 0.0003048 km
  x2,x3=242557.93,393670.74 #(ft)
  x2 = x2*fk+_d2
  x3 = x3*fk+_d3
  wlName = _welldir+"w281.dat"
  ws = readLog(n1,n2,wlName)
  mul(ws[1],fk,ws[1])
  return x2,x3,ws

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
  #return welllogs
  logs = []
  for il in range(len(welllogs)):
    if abs(il-3)>0:
      logs.append(welllogs[il])
  return logs

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

#############################################################################
# read/write fault skins

def skinName(basename,index):
  return basename+("%03i"%(index))
def skinIndex(basename,fileName):
  assert fileName.startswith(basename)
  i = len(basename)
  return int(fileName[i:i+3])

def listAllSkinFiles(basename):
  """ Lists all skins with specified basename, sorted by index. """
  fileNames = []
  for fileName in File(seismicDir).list():
    if fileName.startswith(basename):
      fileNames.append(fileName)
  fileNames.sort()
  return fileNames

def removeAllSkinFiles(basename):
  """ Removes all skins with specified basename. """
  fileNames = listAllSkinFiles(basename)
  for fileName in fileNames:
    File(seismicDir+fileName).delete()

def readSkin(basename,index):
  """ Reads one skin with specified basename and index. """
  return FaultSkin.readFromFile(seismicDir+skinName(basename,index)+".dat")

def readSkins(basename):
  """ Reads all skins with specified basename. """
  fileNames = []
  for fileName in File(seismicDir).list():
    if fileName.startswith(basename):
      fileNames.append(fileName)
  fileNames.sort()
  skins = []
  for iskin,fileName in enumerate(fileNames):
    index = skinIndex(basename,fileName)
    skin = readSkin(basename,index)
    skins.append(skin)
  return skins

def writeSkin(basename,index,skin):
  """ Writes one skin with specified basename and index. """
  FaultSkin.writeToFile(seismicDir+skinName(basename,index)+".dat",skin)

def writeSkins(basename,skins):
  """ Writes all skins with specified basename. """
  for index,skin in enumerate(skins):
    writeSkin(basename,index,skin)

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
