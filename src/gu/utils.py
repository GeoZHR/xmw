"""
Jython utilities for fake image processing.
Author: Dave Hale, Colorado School of Mines
Version: 2014.06.17
"""
from common import *

#############################################################################
# Internal constants

_datdir = "../../../data/seis/gu/"

#############################################################################
# Setup

seismicDir = _datdir
horizonDir = _datdir+"shb5_hor/"
n1,n2,n3 = 1400,950,677
d1,d2,d3 = 0.001,1,1
f1,f2,f3 = 4,388,1508
s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n3,d3,f3)

def getSamplings():
  return s1,s2,s3

def getSeismicDir():
  return seismicDir

def getHorizonDir():
  return horizonDir

#############################################################################
# read/write images

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
# read horizon
def readHorizon(name):
  """ 
  Reads a horizon from a file with specified name.
  """
  fileName = horizonDir+name+".dat"
  image = zerofloat(n2,n3)
  ais = ArrayInputStream(fileName)#,ByteOrder.LITTLE_ENDIAN)
  ais.readFloats(image)
  ais.close()
  return image
def writeHorizon(name):
  """ 
  Reads a horizon from a file with specified name.
  """
  fileName = horizonDir+name+".dat"
  aos = ArrayOutputStream(fileName)#,ByteOrder.LITTLE_ENDIAN)
  aos.writeFloats(image)
  aos.close()
  return image

def writeHorizon(name,image):
  """ 
  Reads a horizon from a file with specified name.
  """
  fileName = horizonDir+name+".dat"
  aos = ArrayOutputStream(fileName)#,ByteOrder.LITTLE_ENDIAN)
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
