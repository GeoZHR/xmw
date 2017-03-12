"""
Jython utilities for fake image processing.
Author: Dave Hale, Colorado School of Mines
Version: 2014.06.17
"""
from common import *

#############################################################################
# Internal constants

_datdir = "../../../data/seis/slt/"
_datdir = "../../../data/seis/fls/"

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
  global s1,s2
  global n1,n2
  global sz,sl,sc
  global nz,nl,nc
  global pngDir
  if name=="seam2dSub1":
    print "setupForSubset: seam 3d sub1"
    seismicDir = _datdir+"sub1/"
    pngDir = "../../../png/slt/3d/"
    n1,n2= 580,1169
    d1,d2 = 1.0,1.0
    #d1,d2,d3 = 0.002,0.025,0.025 # (s,km,km)
    f1,f2 = 0.000,0.000
    #f1,f2,f3 = 220,340,0
    s1,s2 = Sampling(n1,d1,f1),Sampling(n2,d2,f2)
  elif name=="seam2d":
    print "setupForSubset: seam 2d"
    seismicDir = _datdir+"seam/2d/"
    pngDir = "../../../png/fls/seam/2d/"
    n1,n2 = 400,400 #gs = copy(400,400,50,150,gx)
    d1,d2 = 1.0,1.0
    #d1,d2,d3 = 0.002,0.025,0.025 # (s,km,km)
    f1,f2 = 0.000,0.000
    #f1,f2,f3 = 220,740,350
    s1,s2 = Sampling(n1,d1,f1),Sampling(n2,d2,f2)
  elif name=="seam2dSub2":
    print "setupForSubset: seam 2d sub2"
    seismicDir = _datdir+"sub2/"
    pngDir = "../../../png/fls/seam/3d/"
    #n1,n2,n3 = 751,1169,1002
    n1,n2 = 580,1002
    d1,d2 = 1.0,1.0
    #d1,d2,d3 = 0.002,0.025,0.025 # (s,km,km)
    f1,f2 = 0.000,0.000
    #f1,f2,f3 = 220,740,350
    s1,s2 = Sampling(n1,d1,f1),Sampling(n2,d2,f2)
  elif name=="pik":
    print "setupForSubset: pik"
    seismicDir = _datdir+"pik/"
    pngDir = "../../../png/fls/pik/"
    #n1,n2 = 1800,5516 #full xline slice 1205
    #n1,n2 = 1000,1600 #gs f1,f2=100,1600
    n1,n2 = 1000,101 #full inline slice 2538
    n1,n2 = 240,480 #full inline slice 2538
    d1,d2 = 1.0,1.0
    f1,f2 = 0.000,0.000
    s1,s2 = Sampling(n1,d1,f1),Sampling(n2,d2,f2)
  elif name=="bag2dTwo":
    print "setupForSubset: bag 2d"
    seismicDir = _datdir+"bag/2d/sub2/"
    pngDir = "../../../png/fls/bag/2d/sub2/"
    #n1,n2 = 1800,5516 #full xline slice 1205
    #n1,n2 = 1000,1600 #gs f1,f2=100,1600
    n1,n2 = 1800,2111 #full inline slice 2538
    n1,n2 = 850,1000 #fs=copy(850,1000,100,0,fx)full inline slice 2538
    d1,d2 = 1.0,1.0
    f1,f2 = 0.000,0.000
    #f1,f2 = 100,1600
    s1,s2 = Sampling(n1,d1,f1),Sampling(n2,d2,f2)
  elif name=="bag2dOne":
    print "setupForSubset: bag 2d"
    seismicDir = _datdir+"bag/2d/sub1/"
    pngDir = "../../../png/fls/bag/2d/sub1/"
    n1,n2 = 280,260 #fs=copy(850,1000,100,0,fx)full inline slice 2538
    d1,d2 = 1.0,1.0
    f1,f2 = 0.000,0.000
    #f1,f2 = 100,1600
    s1,s2 = Sampling(n1,d1,f1),Sampling(n2,d2,f2)
  elif name=="bag2d":
    print "setupForSubset: bag 2d"
    seismicDir = _datdir+"bag/2d/"
    pngDir = "../../../png/fls/bag/2d/"
    #n1,n2 = 1800,5516 #full xline slice 1205
    #n1,n2 = 1000,1600 #gs f1,f2=100,1600
    n1,n2 = 950,1129 #gs4370 f1.f2=150,757
    n1,n2 = 1800,2111 #full inline slice 4370
    d1,d2 = 1.0,1.0
    f1,f2 = 0.000,0.000
    #f1,f2 = 100,1600
    s1,s2 = Sampling(n1,d1,f1),Sampling(n2,d2,f2)
  elif name=="seam3dSub":
    print "setupForSubset: seam 3d"
    seismicDir = _datdir+"3dSub/"
    pngDir = "../../../png/fls/seam/3d/"
    #n1,n2,n3 = 751,1169,1002
    n1,n2,n3 = 580,1000,466
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
  return s1,s2

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
  image = zerofloat(n1,n2)
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
  aos = ArrayOutputStream(fileName)
  aos.writeFloats(image)
  aos.close()
  return image

#############################################################################
# read/write fault skins


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
