"""
Jython utilities for fake image processing.
Author: Dave Hale, Colorado School of Mines
Version: 2014.06.17
"""
from common import *

#############################################################################
# Internal constants

_datdir = "../../../data/seis/avo/"

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
  global s1,s2,s3
  global n1,n2,n3
  if name=="fake":
    """ fake image """
    print "setupForSubset: fake"
    seismicDir = _datdir+"fake/20171215-NiuLP-NewTest/"
    n1,n2,n3 = 121,152,153
    n1,n2,n3 = 273,1001,153
    s1,s2,s3 = Sampling(n1),Sampling(n2),Sampling(n3)
  elif name=="tp":
    """ subset of Teapot dome """
    print "setupForSubset: tp"
    seismicDir = _datdir+"tp/"
    n1,n2,n3 = 251,357,143
    d1,d2,d3 = 1.0,1.0,1.0 
    f1,f2,f3 = 0.0,0.0,0.0 # = 0.000,0.000,0.000
    s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n3,d3,f3)
  elif name=="pnz":
    """ subset of Teapot dome """
    print "setupForSubset: pnz"
    seismicDir = _datdir+"pnz/"
    n1,n2,n3 = 300,450,450
    d1,d2,d3 = 1.0,1.0,1.0 
    f1,f2,f3 = 0.0,0.0,0.0 # = 0.000,0.000,0.000
    s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n3,d3,f3)
  elif name=="f3d":
    """ subset of F3 """
    print "setupForSubset: f3d"
    seismicDir = _datdir+"f3d/"
    #n1,n2,n3 = 2121,945,645 #rxf
    n1,n2,n3 = 1940,945,645
    d1,d2,d3 = 1.0,1.0,1.0 
    f1,f2,f3 = 0.0,0.0,0.0 # = 0.000,0.000,0.000
    s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n3,d3,f3)
  elif name=="f3dSub":
    """ subset of F3Sub """
    print "setupForSubset: f3dSub"
    seismicDir = _datdir+"f3d/sub/"
    #n1,n2,n3 = 2897,400,500
    n1,n2,n3 = 1547,350,500
    d1,d2,d3 = 1.0,1.0,1.0 
    f1,f2,f3 = 0.0,0.0,0.0 # = 0.000,0.000,0.000
    s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n3,d3,f3)
  elif name=="ch":
    """ subset with channels """
    print "setupForSubset: ch"
    seismicDir = _datdir+"ch/"
    #n1,n2,n3 = 2897,400,500
    n1,n2,n3 = 240,880,500
    d1,d2,d3 = 1.0,1.0,1.0 
    f1,f2,f3 = 0.0,0.0,0.0 # = 0.000,0.000,0.000
    s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n3,d3,f3)
  elif name=="f3dFaultSub":
    """ subset of F3Sub """
    print "setupForSubset: f3dSub"
    seismicDir = _datdir+"f3d/sub/"
    #n1,n2,n3 = 2897,400,500
    n1,n2,n3 = 194,350,500
    d1,d2,d3 = 1.0,1.0,1.0 
    f1,f2,f3 = 0.0,0.0,0.0 # = 0.000,0.000,0.000
    s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n3,d3,f3)
  elif name=="tp1":
    print "setupForSubset: subt"
    seismicDir = _datdir+"tp1/"
    n1,n2,n3 = 1025,240,80
    d1,d2,d3 = 1.0,1.0,1.0 
    f1,f2,f3 = 0.0,0.0,0.0
    #d1,d2,d3 = 0.002,0.025,0.025 # (s,km,km)
    #f1,f2,f3 = 0.000,d2*29,d3*46
    s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n3,d3,f3)

def getSamplings():
  return s1,s2,s3

def getSeismicDir():
  return seismicDir

#############################################################################
# read/write images
def uncName(basename,index):
  return basename+("%04i"%(index))

def readTxtLogs(name):
  fileName = seismicDir+name+".txt"
  tr = TxtReader()
  rvps = tr.readLogs(fileName)
  return rvps


def readImage(name):
  """ 
  Reads an image from a file with specified name.
  name: base name of image file; e.g., "tpsz"
  """
  fileName = seismicDir+name+".dat"
  image = zerofloat(n1,n2,n3)
  ais = ArrayInputStream(fileName)
  ais.readFloats(image)
  ais.close()
  return image

def readImage3D(m1,m2,m3,name):
  """ 
  Reads an image from a file with specified name.
  name: base name of image file; e.g., "tpsz"
  """
  fileName = seismicDir+name+".dat"
  image = zerofloat(m1,m2,m3)
  ais = ArrayInputStream(fileName,ByteOrder.LITTLE_ENDIAN)
  ais.readFloats(image)
  ais.close()
  return image

def readImage2D(n1,n2,name):
  """ 
  Reads an image from a file with specified name.
  name: base name of image file; e.g., "tpsz"
  """
  fileName = seismicDir+name+".dat"
  image = zerofloat(n1,n2)
  ais = ArrayInputStream(fileName)
  ais.readFloats(image)
  ais.close()
  return image

def readImage2DL(n1,n2,name):
  """ 
  Reads an image from a file with specified name.
  name: base name of image file; e.g., "tpsz"
  """
  fileName = seismicDir+name+".dat"
  image = zerofloat(n1,n2)
  ais = ArrayInputStream(fileName,ByteOrder.LITTLE_ENDIAN)
  ais.readFloats(image)
  ais.close()
  return image

def writeImage(name,image):
  """ 
  Writes an image to a file with specified name.
  name: base name of image file; e.g., "tpgp"
  image: the image
  """
  fileName = seismicDir+name+".dat"
  aos = ArrayOutputStream(fileName,ByteOrder.LITTLE_ENDIAN)
  aos.writeFloats(image)
  aos.close()
  return image

def writeImageL(name,image):
  """ 
  Writes an image to a file with specified name.
  name: base name of image file; e.g., "tpgp"
  image: the image
  """
  fileName = seismicDir+name+".dat"
  aos = ArrayOutputStream(fileName,ByteOrder.LITTLE_ENDIAN)
  aos.writeFloats(image)
  aos.close()
  return image
def readUnc(basename,index):
  """ Reads one unconformity with specified basename and index. """
  fileName = seismicDir+uncName(basename,index)+".dat"
  unc = zerofloat(n2,n3)
  ais = ArrayInputStream(fileName)
  ais.readFloats(unc)
  ais.close()
  return unc

def readUncs(basename):
  """ Reads all unconformities with specified basename. """
  fileNames = []
  for fileName in File(seismicDir).list():
    if fileName.startswith(basename):
      fileNames.append(fileName)
  fileNames.sort()
  uncs = []
  for iskin,fileName in enumerate(fileNames):
    index = skinIndex(basename,fileName)
    unc = readUnc(basename,index)
    uncs.append(unc)
  return uncs

def writeUnc(basename,index,unc):
  """ Writes one unconformity with specified basename and index. """
  fileName = seismicDir+uncName(basename,index)+".dat"
  aos = ArrayOutputStream(fileName)
  aos.writeFloats(unc)
  aos.close()
  return unc

def writeUncs(basename,uncs):
  """ Writes all unconformities with specified basename. """
  for index,unc in enumerate(uncs):
    writeUnc(basename,index,unc)


#############################################################################
# read/write fault skins

def skinName(basename,index):
  return basename+("%04i"%(index))
def skinIndex(basename,fileName):
  assert fileName.startswith(basename)
  i = len(basename)
  return int(fileName[i:i+4])

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
