"""
Jython utilities for edge detection.
Author: Xinming Wu, Colorado School of Mines
Version: 2016.01.28
"""
from common import *

#############################################################################
# Internal constants

_datdir = "../../../data/seis/spv/"
_pngdir = "../../../png/spv/"

#############################################################################
# Setup

def setupForSubset(name):
  """
  Setup for a specified directory includes:
    seismic directory
    samplings s1,s2
  Example: setupForSubset("pnz")
  """
  global pngDir
  global seismicDir
  global s1,s2,s3
  global n1,n2,n3
  if name=="fault2d":
    """ 2d fault """
    print "setupForSubset: 2d fault"
    pngDir = _pngdir+"2d/f3d/"
    seismicDir = _datdir+"fault/2d/"
    n1,n2 = 651,601 #fxnwc
    n1,n2 = 300,1200 #gx238 crf dataset
    n1,n2 = 400,801 #cylde200
    n1,n2 = 222,440 #f3d75s
    n1,n2 = 300,550 #gx396
    n1,n2 = 400,420 #ep75
    n1,n2 = 380,591 #ep56
    d1,d2 = 1,1 # (s,km/s)
    f1,f2 = 0,0
    s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n2,d2,f2)
  elif name=="f3d2d":
    """ 2d fault """
    print "setupForSubset: 2d fault"
    pngDir = _pngdir+"fault/2d/f3d/"
    seismicDir = _datdir+"fault/2d/f3d/"
    n1,n2 = 222,440 #f3d75s
    d1,d2 = 1,1 # (s,km/s)
    f1,f2 = 0,0
    s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n2,d2,f2)
  elif name=="crf2d":
    """ costa rica field """
    print "setupForSubset: crf"
    pngDir = _pngdir+"2d/crf/"
    seismicDir = _datdir+"fault/2d/crf/"
    n1,n2 = 210,825
    n1,n2 = 180,500
    d1,d2 = 1,1 # (s,km/s)
    f1,f2 = 0,0
    s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n2,d2,f2)
  elif name=="campos2d":
    """ costa rica field """
    print "setupForSubset: crf"
    pngDir = _pngdir+"fault/2d/campos/"
    seismicDir = _datdir+"fault/2d/campos/"
    n1,n2 = 300,550
    d1,d2 = 1,1 # (s,km/s)
    f1,f2 = 0,0
    s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n2,d2,f2)
  elif name=="tj2d":
    """ costa rica field """
    print "setupForSubset: crf"
    pngDir = _pngdir+"fault/2d/campos/"
    seismicDir = "../../../data/seis/mhe/tj/"
    n1,n2 = 187,801
    d1,d2 = 1,1 # (s,km/s)
    f1,f2 = 0,0
    s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n2,d2,f2)
  elif name=="edge":
    """ edge """
    print "setupForSubset: edge"
    pngDir = _pngdir+"fault/2d/"
    seismicDir = _datdir+"edge/"
    n1,n2 = 400,801 #cylde200
    n1,n2 = 300,1200 #gx238 crf dataset
    n1,n2 = 651,601 #fxnwc
    n1,n2 = 222,440 #f3d75s
    n1,n2 = 380,591 #ep56
    n1,n2 = 1052,698 #teste
    n1,n2 = 321,481 #test
    #n1,n2 = 321,481 #135069.dat
    d1,d2 = 1,1 # (s,km/s)
    f1,f2 = 0,0
    s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n2,d2,f2)
  elif name=="bsds500":
    """ bsds500 """
    print "setupForSubset: bsds500"
    pngDir = _pngdir+"BSDS500/data/images/test/"
    seismicDir = _datdir+"BSDS500/data/images/"
    n1,n2 = 321,481 #test
    d1,d2 = 1,1 # (s,km/s)
    f1,f2 = 0,0
    s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n2,d2,f2)
  elif name=="campos":
    """ campos """
    print "setupForSubset: campos"
    pngDir = _pngdir+"fault/3d/campos/"
    seismicDir = _datdir+"fault/3d/campos/"
    n1,n2,n3 = 300,600,550 #ep56
    d1,d2,d3 = 1,1,1 # (s,km/s)
    f1,f2,f3 = 0,0,0
    s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n3,d3,f3)
  elif name=="crf":
    """ costa rica field """
    print "setupForSubset: crf"
    seismicDir = _datdir+"fault/3d/crf/"
    pngDir = _pngdir+"3d/crf/"
    n1,n2,n3 = 210,920,825 #ep56
    d1,d2,d3 = 1,1,1 # (s,km/s)
    f1,f2,f3 = 0,0,0
    s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n3,d3,f3)
  elif name=="clyde":
    """ costa rica field """
    print "setupForSubset: clyde"
    seismicDir = _datdir+"fault/3d/clyde/"
    pngDir = _pngdir+"3d/clyde/"
    n1,n2,n3 = 400,801,300 #ep56
    d1,d2,d3 = 1,1,1 # (s,km/s)
    f1,f2,f3 = 0,0,0
    s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n3,d3,f3)
  elif name=="f3d":
    print "setupForSubset: f3d"
    pngDir = _pngdir+"3d/f3d/"
    seismicDir = _datdir+"fault/3d/f3d/"
    n1,n2,n3 = 65,380,591
    n1,n2,n3 = 100,400,420
    d1,d2,d3 = 1.0,1.0,1.0 
    #j1,j2,j3 = 344,0,0
    #d1,d2,d3 = 0.004,0.025,0.024999 # (s,km,km)
    #f1,f2,f3 = 0.472,0.0,0.0 # = 0.000,0.000,0.000
    f1,f2,f3 = 0.0,0.0,0.0 # = 0.000,0.000,0.000
    s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n3,d3,f3)
  elif name=="xue":
    """ xue """
    print "setupForSubset: xue"
    pngDir = _pngdir+"xue/"
    seismicDir = _datdir+"xue/"
    n1,n2 = 1251,7
    d1,d2 = 1.0,1.0 
    f1,f2 = 0.0,0.0
    s1,s2 = Sampling(n1,d1,f1),Sampling(n2,d2,f2)
    n3,d3,f3 = 1,1,1
    s3 = Sampling(n3,d3,f3)
  elif name=="scan":
    """ xue """
    print "setupForSubset: scan"
    pngDir = _pngdir+"scan/"
    seismicDir = _datdir+"scan/"
    n1,n2 = 1000,101
    d1,d2 = 1.0,1.0 
    f1,f2 = 0.0,0.0
    s1,s2 = Sampling(n1,d1,f1),Sampling(n2,d2,f2)
    n3,d3,f3 = 1,1,1
    s3 = Sampling(n3,d3,f3)
  elif name=="fd2":
    """ fd2 """
    print "setupForSubset: fd2"
    pngDir = _pngdir+"fd/"
    seismicDir = _datdir+"fd/"
    n1,n2 = 442,951
    d1,d2 = 1.0,1.0 
    f1,f2 = 0.0,0.0
    s1,s2 = Sampling(n1,d1,f1),Sampling(n2,d2,f2)
    n3,d3,f3 = 1,1,1
    s3 = Sampling(n3,d3,f3)
  elif name=="dgb":
    """ dgb """
    print "setupForSubset: dgb"
    pngDir = _pngdir+"dgb/"
    seismicDir = _datdir+"dgb/"
    n1,n2 = 301,300
    d1,d2 = 1.0,1.0 
    f1,f2 = 0.0,0.0
    s1,s2 = Sampling(n1,d1,f1),Sampling(n2,d2,f2)
    n3,d3,f3 = 1,1,1
    s3 = Sampling(n3,d3,f3)
  elif name=="channel":
    """ nwc channel (slice 120) """
    print "setupForSubset: channel"
    pngDir = _pngdir+"channel/"
    seismicDir = _datdir+"channel/"
    n1,n2 = 601,301
    d1,d2 = 1.0,1.0 
    f1,f2 = 0.0,0.0 # = 0.000,0.000,0.000
    s1,s2 = Sampling(n1,d1,f1),Sampling(n2,d2,f2)
    n3,d3,f3 = 1,1,1
    s3 = Sampling(n3,d3,f3)
  elif name=="surface":
    print "setupForSubset: surface"
    seismicDir = _datdir+"surface/"
    pngDir = _pngdir+"surface/"
    n1,n2,n3 = 242,611,591
    n1,n2,n3 = 55,254,137
    d1,d2,d3 = 1.0,1.0,1.0 
    f1,f2,f3 = 0.000,0.000,0.000
    s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n3,d3,f3)
  elif name=="semblance3d":
    print "setupForSubset: 3D semblance"
    seismicDir = _datdir+"semblance/3d/"
    pngDir = _pngdir+"semblance/3d/"
    n1,n2,n3 = 1000,201,250
    d1,d2,d3 = 1.0,1.0,1.0 
    f1,f2,f3 = 0.000,0.000,0.000
    s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n3,d3,f3)
  elif name=="env":
    print "setupForSubset: 3D env"
    seismicDir = _datdir+"env/"
    pngDir = _pngdir+"env/"
    n1,n2,n3 = 1001,1028,91
    n1,n2,n3 = 91,1028,1001
    d1,d2,d3 = 1.0,1.0,1.0 
    f1,f2,f3 = 0.000,0.000,0.000
    s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n3,d3,f3)

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
def readImageChannels(basename):
  """ 
  Reads three channels of a color image
  """
  fileName = seismicDir+basename+".jpg"
  il = ImageLoader()
  image = il.readThreeChannels(fileName)
  return image
def readColorImage(basename):
  """ 
  Reads three channels of a color image
  """
  fileName = seismicDir+basename+".jpg"
  il = ImageLoader()
  image = il.readColorImage(fileName)
  return image

def readImage(basename):
  """ 
  Reads an image from a file with specified basename
  """
  fileName = seismicDir+basename+".dat"
  image = zerofloat(n1,n2)
  ais = ArrayInputStream(fileName)
  ais.readFloats(image)
  ais.close()
  return image
def readImageX(m1,m2,basename):
  """ 
  Reads an image from a file with specified basename
  """
  fileName = seismicDir+basename+".dat"
  image = zerofloat(m1,m2)
  ais = ArrayInputStream(fileName)
  ais.readFloats(image)
  ais.close()
  return image

def readImage3D(basename):
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

def readImage2L(n1,n2,basename):
  """ 
  Reads an image from a file with specified basename
  """
  fileName = seismicDir+basename+".dat"
  image = zerofloat(n1,n2)
  ais = ArrayInputStream(fileName,ByteOrder.LITTLE_ENDIAN)
  ais.readFloats(image)
  ais.close()
  return image


def readImage1D(basename):
  """ 
  Reads an image from a file with specified basename
  """
  fileName = seismicDir+basename+".dat"
  image = zerofloat(n1)
  ais = ArrayInputStream(fileName)
  ais.readFloats(image)
  ais.close()
  return image

def readImage1L(basename):
  """ 
  Reads an image from a file with specified basename
  """
  fileName = seismicDir+basename+".dat"
  image = zerofloat(n1)
  ais = ArrayInputStream(fileName,ByteOrder.LITTLE_ENDIAN)
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

def writeImageL(basename,image):
  """ 
  Writes an image to a file with specified basename
  """
  fileName = seismicDir+basename+".dat"
  aos = ArrayOutputStream(fileName,ByteOrder.LITTLE_ENDIAN)
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
