"""
Jython utilities for channel enhancing.
Author: Xinming Wu, Colorado School of Mines
Version: 2016.01.28
"""
from common import *

#############################################################################
# Internal constants

_datdir = "../../../data/seis/scl/"
_pngdir = "../../../png/scl/"

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
  global s1,s2
  global n1,n2
  if name=="pnz":
    """ pnz horizon """
    print "setupForSubset: pnz"
    pngDir = _pngdir+"pnz/"
    seismicDir = _datdir+"pnz/"
    n1,n2 = 879,752
    d1,d2 = 1.0,1.0
    #d1,d2,d3 = 0.002,0.008,0.008 # (s,km,km)
    f1,f2 = 0.0,0.0
    s1,s2 = Sampling(n1,d1,f1),Sampling(n2,d2,f2)
  elif name=="parihaka":
    """ parihaka horizon """
    print "setupForSubset: parihaka"
    pngDir = _pngdir+"parihaka/"
    seismicDir = _datdir+"parihaka/"
    n1,n2 = 946,830
    d1,d2 = 1.0,1.0 
    #d1,d2,d3 = 0.002,0.008,0.008 # (s,km,km)
    f1,f2 = 0.0,0.0 # = 0.000,0.000,0.000
    s1,s2 = Sampling(n1,d1,f1),Sampling(n2,d2,f2)
  elif name=="vessel":
    """ nwc horizon """
    print "setupForSubset: vessel"
    pngDir = _pngdir+"vessel/"
    seismicDir = _datdir+"vessel/"
    n1,n2 = 512,512
    d1,d2 = 1.0,1.0 
    #d1,d2,d3 = 0.002,0.008,0.008 # (s,km,km)
    f1,f2 = 0.0,0.0 # = 0.000,0.000,0.000
    s1,s2 = Sampling(n1,d1,f1),Sampling(n2,d2,f2)

  elif name=="nwc":
    """ nwc horizon """
    print "setupForSubset: nwc"
    pngDir = _pngdir+"nwc/"
    seismicDir = _datdir+"nwc/"
    n1,n2 = 601,401
    d1,d2 = 1.0,1.0 
    #d1,d2,d3 = 0.002,0.008,0.008 # (s,km,km)
    f1,f2 = 0.0,0.0 # = 0.000,0.000,0.000
    s1,s2 = Sampling(n1,d1,f1),Sampling(n2,d2,f2)

  elif name=="fake":
    """ fake """
    print "setupForSubset: fake"
    pngDir = _pngdir+"fake/"
    seismicDir = _datdir+"fake/"
    n1,n2 = 200,200
    d1,d2 = 1.0,1.0 
    #d1,d2,d3 = 0.002,0.008,0.008 # (s,km,km)
    f1,f2 = 0.0,0.0 # = 0.000,0.000,0.000
    s1,s2 = Sampling(n1,d1,f1),Sampling(n2,d2,f2)
  elif name=="tccs":
    """ tccs """
    print "setupForSubset: tccs"
    pngDir = _pngdir+"tccs/"
    seismicDir = _datdir+"tccs/"
    n1,n2 = 344,826
    d1,d2 = 1.0,1.0 
    #d1,d2,d3 = 0.002,0.008,0.008 # (s,km,km)
    f1,f2 = 0.0,0.0 # = 0.000,0.000,0.000
    s1,s2 = Sampling(n1,d1,f1),Sampling(n2,d2,f2)


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

def writeImage(basename,image):
  """ 
  Writes an image to a file with specified basename
  """
  fileName = seismicDir+basename+".dat"
  aos = ArrayOutputStream(fileName)
  aos.writeFloats(image)
  aos.close()
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
