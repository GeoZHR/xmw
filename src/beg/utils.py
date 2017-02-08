"""
Jython utilities for fake image processing.
Author: Dave Hale, Colorado School of Mines
Version: 2014.06.17
"""
from common import *

#############################################################################
# Internal constants

_datdir = "../../../data/seis/beg/"

#############################################################################
# Setup
def setupForSubset(name):
  """
  Setup for a specified directory includes:
    seismic directory
    samplings s1,s2,s3
  Example: setupForSubset("hongliu")
  """
  global seismicDir
  global welllogDir
  global s1,s2,s3
  global n1,n2,n3
  global sz,sl,sc
  global nz,nl,nc
  if name=="hongliu":
    print "setupForSubset: hongliu"
    seismicDir = _datdir+"hongliu/"
    n1,n2,n3 = 751,201,201
    d1,d2,d3 = 1.0,1.0,1.0 
    #d1,d2,d3 = 0.002,0.025,0.025 # (s,km,km)
    f1,f2,f3 = 0.000,400,1500
    s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n3,d3,f3)
  elif name=="australia":
    print "setupForDataset: australia"
    seismicDir = _datdir+"xavier/australia/"
    n1,n2,n3 = 1640,1065,1905
    d1,d2,d3 = 1.0,1.0,1.0 
    f1,f2,f3 = 0.0,0.0,0.0
    #d1,d2,d3 = 0.004,1,2 # (s,tracet,trace)
    #f1,f2,f3 = 128,0,0 (sample,trace,trace)
    s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n3,d3,f3)
  elif name=="bahamas":
    print "setupForDataset: bahamas"
    seismicDir = _datdir+"xavier/bahamas/"
    #n1,n2,n3 = 200,2359,4320 #gs
    #n1,n2,n3 = 120,2359,4320 #gxc
    n1,n2,n3 = 381,2359,4320 #gx
    d1,d2,d3 = 1.0,1.0,1.0 
    f1,f2,f3 = 0.0,0.0,0.0
    #f1,f2,f3 = 100,2000,6400 (sample,trace,trace)
    s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n3,d3,f3)
  elif name=="wasson":
    print "setupForDataset: bahamas"
    seismicDir = _datdir+"xavier/wasson/"
    n1,n2,n3 = 381,2359,4320
    d1,d2,d3 = 1.0,1.0,1.0 
    f1,f2,f3 = 0.0,0.0,0.0
    s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n3,d3,f3)
  elif name=="mason":
    print "setupForDataset: mason"
    seismicDir = _datdir+"mason/"
    n1,n2,n3 = 55,250,200
    d1,d2,d3 = 1.0,1.0,1.0 
    f1,f2,f3 = 0.0,0.0,0.0
    s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n3,d3,f3)
  elif name=="nathan":
    print "setupForDataset: nathan"
    seismicDir = _datdir+"nathan/"
    n1,n2,n3 = 601,4974,660
    d1,d2,d3 = 1.0,1.0,1.0 
    f1,f2,f3 = 0.0,0.0,0.0
    #d1,d2,d3 = 0.002,0.025,0.025 # (s,km,km)
    #f1,f2,f3 = 0.000,5400,10744
    s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n3,d3,f3)
  elif name=="nathanSub1":
    print "setupForSubset: nathanSub1"
    seismicDir = _datdir+"nathan/sub1/"
    n1,n2,n3 = 400,800,550
    d1,d2,d3 = 1.0,1.0,1.0 
    f1,f2,f3 = 0.0,0.0,0.0
    #d1,d2,d3 = 0.002,0.025,0.025 # (s,km,km)
    #f1,f2,f3 = 0.000,5400,10744
    s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n3,d3,f3)
  elif name=="nathanSub2":
    print "setupForSubset: nathanSub2"
    seismicDir = _datdir+"nathan/sub2/"
    n1,n2,n3 = 400,800,550
    d1,d2,d3 = 1.0,1.0,1.0 
    f1,f2,f3 = 0.0,0.0,0.0
    #d1,d2,d3 = 0.002,0.025,0.025 # (s,km,km)
    #f1,f2,f3 = 0.000,5400,10744
    s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n3,d3,f3)
  elif name=="nathanSub3":
    print "setupForSubset: nathanSub3"
    seismicDir = _datdir+"nathan/sub3/"
    n1,n2,n3 = 601,3675,550  #fx = copy(n1,3675,550,0,1100,60,fx)
    d1,d2,d3 = 1.0,1.0,1.0 
    f1,f2,f3 = 0.0,0.0,0.0
    #d1,d2,d3 = 0.002,0.025,0.025 # (s,km,km)
    #f1,f2,f3 = 0.000,5400,10744
    s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n3,d3,f3)
  elif name=="nathanSub5":
    print "setupForSubset: nathanSub5"
    seismicDir = _datdir+"nathan/sub5/"
    n1,n2,n3 = 250,600,450  #fx = copy(n1,3675,550,0,1100,60,fx)
    d1,d2,d3 = 1.0,1.0,1.0 
    f1,f2,f3 = 0.0,0.0,0.0
    #d1,d2,d3 = 0.002,0.025,0.025 # (s,km,km)
    #f1,f2,f3 = 0.000,5400,10744
    s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n3,d3,f3)
  elif name=="jake":
    print "setupForSubset: jake"
    seismicDir = _datdir+"jake/"
    n1,n2,n3 = 211,1656,1001
    d1,d2,d3 = 1.0,1.0,1.0 
    f1,f2,f3 = 0.0,0.0,0.0
    #d1,d2,d3 = 0.002,0.025,0.025 # (s,km,km)
    #f1,f2,f3 = 0.000,5400,10744
    s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n3,d3,f3)
  elif name=="jakeSub2":
    print "setupForSubset: jakeSub2"
    seismicDir = _datdir+"jake/sub2/"
    n1,n2,n3 = 426,1856,830
    d1,d2,d3 = 1.0,1.0,1.0 
    f1,f2,f3 = 0.0,0.0,0.0
    #d1,d2,d3 = 0.002,0.025,0.025 # (s,km,km)
    #f1,f2,f3 = 0.000,5400,10744
    s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n3,d3,f3)
  elif name=="bp":
    print "setupForSubset: bp"
    seismicDir = _datdir+"bp/"
    n1,n2,n3 = 501,2110,2818
    d1,d2,d3 = 1.0,1.0,1.0 
    #d1,d2,d3 = 0.002,0.025,0.025 # (s,km,km)
    f1,f2,f3 = 0.000,5400,10744
    s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n3,d3,f3)
  elif name=="bpSub":
    print "setupForSubset: bp subset"
    seismicDir = _datdir+"bp/sub/"
    n1,n2,n3 = 191,1005,1251
    d1,d2,d3 = 1.0,1.0,1.0 
    #d1,d2,d3 = 0.002,0.025,0.025 # (s,km,km)
    #f1,f2,f3 = 0.000,5400,11400
    f1,f2,f3 = 0.000,0.000,0.000
    s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n3,d3,f3)
  elif name=="bpSub1":
    print "setupForSubset: bp subet"
    seismicDir = _datdir+"bp/sub1/"
    n1,n2,n3 = 191,1005,1251
    d1,d2,d3 = 1.0,1.0,1.0 
    #d1,d2,d3 = 0.002,0.025,0.025 # (s,km,km)
    #f1,f2,f3 = 0.000,5400,11400
    f1,f2,f3 = 0.000,0.000,0.000
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
def readHorizon(name):
  """ 
  Reads an image from a file with specified name.
  name: base name of image file; e.g., "tpsz"
  """
  fileName = seismicDir+name+".dat"
  image = zerofloat(n2,n3)
  ais = ArrayInputStream(fileName)
  ais.readFloats(image)
  ais.close()
  return image
def readHorizons(ns,name):
  """ 
  Reads an image from a file with specified basename
  """
  fileName = seismicDir+name+".dat"
  image = zerofloat(n2,n3,ns)
  ais = ArrayInputStream(fileName)
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

def readImageL(name):
  """ 
  Reads an image from a file with specified name.
  name: base name of image file; e.g., "tpsz"
  """
  fileName = seismicDir+name+".dat"
  image = zerofloat(n1,n2,n3)
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

def writeAscii(name,image):
  """ 
  Writes an image to a file with specified name.
  name: base name of image file; e.g., "tpgp"
  image: the image
  """
  fileName = seismicDir+name+".xyz"
  hp = Helpler()
  hp.writeAsciHorizon(name,image)
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

#############################################################################
# read/write fault skins

def skinName(basename,index):
  return basename+("%05i"%(index))
def skinIndex(basename,fileName):
  assert fileName.startswith(basename)
  i = len(basename)
  return int(fileName[i:i+5])

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
