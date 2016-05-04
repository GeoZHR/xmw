"""
Jython utilities for fake image processing.
Author: Xinming Wu, Colorado School of Mines
Version: 2016.04.29
"""
from common import *

#############################################################################
# Internal constants

_datdir = "../../../data/seis/"

#############################################################################
# Setup
def setupForSubset(name):
  """
  Setup for a specified directory includes:
    seismic directory
    samplings s1,s2
  Example: setupForSubset("oregan")
  """
  global seismicDir
  global s1,s2
  global n1,n2
  if name=="oregan":
    """ oregan image """
    print "setupForSubset: oregan"
    seismicDir = _datdir+"oregan/"
    n1,n2 = 2000,60201 # f1=1200,f2=0
    s1,s2 = Sampling(n1),Sampling(n2)
  elif name=="fake":
    """ oregan image """
    print "setupForSubset: oreganSub"
    seismicDir = _datdir+"oregan/fake/"
    n1,n2 = 121,152
    s1,s2 = Sampling(n1),Sampling(n2)
  elif name=="oreganSub":
    """ oregan image """
    print "setupForSubset: oreganSub"
    seismicDir = _datdir+"oregan/sub/"
    n1,n2 = 2351,15051
    s1,s2 = Sampling(n1),Sampling(n2)
  elif name=="oreganSub1":
    """ oregan image """
    print "setupForSubset: oreganSub"
    seismicDir = _datdir+"oregan/sub1/"
    n1,n2 = 1200,3000
    s1,s2 = Sampling(n1),Sampling(n2)

def getSamplings():
  return s1,s2

def getSeismicDir():
  return seismicDir

#############################################################################
# read/write images

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
