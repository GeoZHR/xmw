import sys

from java.awt import *
from java.io import *
from java.lang import *
from javax.swing import *

from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.io import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.util import *
from edu.mines.jtk.util.ArrayMath import *

from hvc import *

seismicDir = "../../../data/seis/hvc/2d/tpd/"
ffile = "tp73"
ffile = "fx126s"
ffile = "gxsb1"
#s1 = Sampling(251,0.004,0.500)
#s2 = Sampling(357,0.025,0.000)
s1 = Sampling(251,1.0,0.0)
s2 = Sampling(357,1.0,0.0)
s1 = Sampling(202,1.0,0.0)
s2 = Sampling(400,1.0,0.0)
s1 = Sampling(210,1.0,0.0)
s2 = Sampling(2000,1.0,0.0)

n1,n2 = s1.count,s2.count
d1,d2 = s1.delta,s2.delta

def main(args):
  #flattenWithSlopes()
  flattenWithSlopesAndCorrelations()
def flattenWithSlopes():
  f = readImage(ffile)
  f = gain(f)
  sig1,sig2 = 8.0,2.0 # good for Teapot Dome image tp73
  pmax = 10.0
  lsf = LocalSlopeFinder(sig1,sig2,pmax)
  p2 = zerofloat(n1,n2)
  wp = zerofloat(n1,n2)
  lsf.findSlopes(f,p2,wp)
  wp = pow(wp,14)
  fl = Flattener2Dw()
  fl.setWeight1(0.00)
  fl.setIterations(0.01,200)
  fl.setSmoothings(8.0,2.0)
  #ks,dt=None,None
  fm = fl.getMappingsFromSlopes(s1,s2,p2,wp)
  g = fm.flatten(f)
  plot(s1,s2,f,title="Input",png="f")
  plot(s1,s2,g,title="Flattened",png="g")
  plot(s1,s2,f,u=fm.u1,title="Input",png="f")
  plot(s1,s2,fm.u1,cmin=10,cmax=250,cmap=ColorMap.JET,title="RGT",png="g")
def flattenWithSlopesAndCorrelations():
  f = readImage(ffile)
  f = gain(f)
  sig1,sig2 = 8.0,2.0 # good for Teapot Dome image tp73
  pmax = 10.0
  lsf = LocalSlopeFinder(sig1,sig2,pmax)
  p2 = zerofloat(n1,n2)
  wp = zerofloat(n1,n2)
  lsf.findSlopes(f,p2,wp)
  wp = pow(wp,14)
  ps = []
  for ik in range(0,n2,5):
    ps.append(ik)
  gcf = GlobalCorrelationFinder(-35,35)
  gcf.setStrainMax(0.3)
  dks = gcf.findCorrelations(150,ps,f)
  ns = len(dks)
  ks = zeroint(3,ns)
  dt = zerofloat(n1,ns)
  gcf.unpack(dks,ks,dt)
  fl = Flattener2Dw()
  fl.setWeight1(0.0)
  fl.setIterations(0.01,500)
  fl.setSmoothings(8.0,2.0)
  #ks,dt=None,None
  fm = fl.getMappingsFromSlopesAndCorrelations(s1,s2,0.001,p2,wp,ks,dt)
  g = fm.flatten(f)
  plot(s1,s2,f,title="Input",png="f")
  plot(s1,s2,g,title="Flattened",png="g")
  plot(s1,s2,f,u=fm.u1,title="Input",png="f")
  plot(s1,s2,fm.u1,cmin=10,cmax=250,cmap=ColorMap.JET,title="RGT",png="g")


def flip2(f):
  n2 = len(f)
  for i2 in range(n2/2):
    fi2 = f[i2]
    f[i2] = f[n2-1-i2]
    f[n2-1-i2] = fi2

def gain(x):
  n2 = len(x)
  n1 = len(x[0])
  g = mul(x,x) 
  ref = RecursiveExponentialFilter(50.0)
  ref.apply1(g,g)
  y = zerofloat(n1,n2)
  div(x,sqrt(g),y)
  return y

#############################################################################
# graphics

gray = ColorMap.GRAY
jet = ColorMap.JET
pngDir = None
#pngDir = "./png/"
def plot(s1,s2,x,u=None,c=None,cmap=ColorMap.GRAY,clab=None,cmin=-2,cmax=2,
         title=None,png=None):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  if title:
    sp.setTitle(title)
  sp.addColorBar(clab)
  sp.setSize(600,900)
  sp.setSize(800,400)
  sp.plotPanel.setColorBarWidthMinimum(80)
  pv = sp.addPixels(s1,s2,x)
  pv.setColorModel(cmap)
  pv.setInterpolation(PixelsView.Interpolation.LINEAR)
  if cmin<cmax:
    pv.setClips(cmin,cmax)
  if u:
    cv = sp.addContours(s1,s2,u)
    cv.setContours(40)
    cv.setLineColor(Color.YELLOW)
  if c:
    k1s,k2s = c
    pv = PointsView(k1s,k2s)
    pv.setOrientation(PointsView.Orientation.X1DOWN_X2RIGHT)
    pv.setLineColor(Color.WHITE)
    pv.setLineWidth(2)
    pv.setMarkStyle(PointsView.Mark.HOLLOW_CIRCLE)
    pv.setMarkColor(Color.WHITE)
    pv.setMarkSize(10)
    sp.add(pv)
  if pngDir and png:
    sp.paintToPng(300,3.333,pngDir+png+".png")

#############################################################################
# utilities

def readImage(name):
  fileName = seismicDir+name+".dat"
  n1,n2 = s1.count,s2.count
  image = zerofloat(n1,n2)
  ais = ArrayInputStream(fileName)
  ais.readFloats(image)
  ais.close()
  return image

def writeImage(name,image):
  fileName = seismicDir+name+".dat"
  aos = ArrayOutputStream(fileName)
  aos.writeFloats(image)
  aos.close()
  return image

#############################################################################
# Run the function main on the Swing thread
import sys
class _RunMain(Runnable):
  def __init__(self,main):
    self.main = main
  def run(self):
    self.main(sys.argv)
def run(main):
  SwingUtilities.invokeLater(_RunMain(main)) 
run(main)
