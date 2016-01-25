"""
Jython utilities for fake image processing.
Author: Dave Hale, Colorado School of Mines
Version: 2014.06.17
"""
from common import *

#############################################################################
# Internal constants

_datdir = "../../../data/seis/slt/"

#############################################################################
# Setup
def setupForSubset(name):
  """
  Setup for a specified directory includes:
    seismic directory
    samplings s1,s2,s3
  Example: setupForSubset("s1")
  """
  global pngDir
  global seismicDir
  global s1,s2,s3
  global n1,n2,n3
  global sz,sl,sc
  global nz,nl,nc
  global pngDir
  if name=="3d":
    print "setupForSubset: 3d"
    seismicDir = _datdir+"3d/"
    pngDir = "../../../png/slt/3d/"
    n1,n2,n3 = 242,611,591
    d1,d2,d3 = 1.0,1.0,1.0 
    #d1,d2,d3 = 0.002,0.025,0.025 # (s,km,km)
    f1,f2,f3 = 0.000,0.000,0.000
    #f1,f2,f3 = 220,340,0
    s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n3,d3,f3)
  elif name=="3dSub":
    print "setupForSubset: 3dSub"
    seismicDir = _datdir+"3d/sub/"
    pngDir = "../../../png/slt/3d/sub/"
    n1,n2,n3 = 242,211,241
    d1,d2,d3 = 1.0,1.0,1.0 
    #d1,d2,d3 = 0.002,0.025,0.025 # (s,km,km)
    f1,f2,f3 = 0.000,0.000,0.000
    #f1,f2,f3 = 220,740,350
    s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n3,d3,f3)
  elif name=="2dSub1":
    print "setupForSubset: 2d"
    seismicDir = _datdir+"2d/sub1/"
    pngDir = "../../../png/slt/2d/sub1/"
    n1,n2 = 162,461  #copy(n1-80,n2-150,80,150)
    d1,d2 = 0.004,0.025 # (s,km,km)
    f1,f2 = 0.004+80*d1,150*d2
    #d1,d2 = 1.00,1.0
    #f1,f2 = 0.00,0.0
    s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n2,d2,f2)
  elif name=="2dSub2":
    print "setupForSubset: 2d"
    seismicDir = _datdir+"2d/sub2/"
    pngDir = "../../../png/slt/2d/sub2/"
    #n1,n2 = 296,1941
    n1,n2 = 162,461  #copy(n1-80,n2-150,80,150)
    d1,d2 = 0.004,0.025 # (s,km,km)
    f1,f2 = 0.004+80*d1,0*d2
    #f1,f2 = 220,340
    s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n2,d2,f2)
  else:
    print "unrecognized subset:",name
    System.exit

def getSamplings():
  return s1,s2,s3

def getSeismicDir():
  return seismicDir

def getPngDir():
  return pngDir

#############################################################################
# read/write images

def readImage2d(name):
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

def writeImage(name,image):
  """ 
  Writes an image to a file with specified name.
  name: base name of image file; e.g., "tpgp"
  image: the image
  """
  fileName = seismicDir+name+".dat"
  aos = ArrayOutputStream(fileName)
  aos.writeFloats(image)
  aos.close()
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
