import sys

from java.awt import *
from java.io import *
from java.nio import *
from java.lang import *
from javax.swing import *

from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.io import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.util import *
from edu.mines.jtk.util.ArrayMath import *

from hv import *
from util import *

seismicDir = "../../../data/seis/tpd/"
ffile = "tp73"
s1 = Sampling(251,1.0,0.0)
s2 = Sampling(357,1.0,0.0)
sp1 = Sampling(251,0.004,0.500)
sp2 = Sampling(357,0.025,0.000)

n1,n2 = s1.count,s2.count
d1,d2 = s1.delta,s2.delta
f1,f2 = sp1.first,sp2.first

def main(args):
  #flatten()
  goRgt()

def goRgt():
  #f = FakeData.seismic2d2011A(n1,n2,30)
  f = readImage(ffile)
  f = gain(f)
  plot(sp1,sp2,f,clab="Amplitude",vlabel="Time (s)",hlabel="Inline (km)",png="f")
  #sigma = 1.0 # good for fake data
  sigma1,sigma2=4.0,2.0 # good for Teapot Dome image tp73
  lof = LocalOrientFilter(sigma1,sigma2)
  u1 = zerofloat(n1,n2)
  u2 = zerofloat(n1,n2)
  wp = zerofloat(n1,n2)
  lof.applyForNormalLinear(f,u1,u2,wp)
  wp = pow(wp,8)
  rgt = RelativeGeologicTime2()
  rgt.setWeight1(0.0)
  rgt.setIterations(0.01,1000)
  rgt.setSmoothings(4.0,8.0)
  u = rgt.applyForRgt(u1,u2,wp)
  mul(u,sp1.delta,u)
  add(u,sp1.first,u)
  plot(sp1,sp2,u,cmap=jet,clab="RGT",vlabel="Time (s)",hlabel="Inline (km)",png="rgt")
  plot(sp1,sp2,f,u=u,clab="RGT",vlabel="Time (s)",hlabel="Inline (km)",png="hs")

def flatten():
  #f = FakeData.seismic2d2011A(n1,n2,30)
  f = readImage(ffile)
  f = gain(f)
  plot(sp1,sp2,f,clab="Amplitude",vlabel="Time (s)",hlabel="Inline (km)",png="f")
  #sigma = 1.0 # good for fake data
  sigma1,sigma2=4.0,1.0 # good for Teapot Dome image tp73
  pmax = 10.0
  lsf = LocalSlopeFinder(sigma1,sigma2,pmax)
  p2 = zerofloat(n1,n2)
  wp = zerofloat(n1,n2)
  lsf.findSlopes(f,p2,wp)
  p2 = mul(d1/d2,p2)
  wp = pow(wp,8)
  fl = Flattener2()
  fl.setWeight1(0.02)
  fl.setIterations(0.01,1000)
  fl.setSmoothings(4.0,8.0)
  fm = fl.getMappingsFromSlopes(s1,s2,p2,wp)
  u = fm.u1
  mul(u,sp1.delta,u)
  add(u,sp1.first,u)
  plot(sp1,sp2,u,cmap=jet,clab="RGT",vlabel="Time (s)",hlabel="Inline (km)",png="rgt")
  plot(sp1,sp2,f,u=u,clab="RGT",vlabel="Time (s)",hlabel="Inline (km)",png="hs")

def flip2(f):
  n2 = len(f)
  for i2 in range(n2/2):
    fi2 = f[i2]
    f[i2] = f[n2-1-i2]
    f[n2-1-i2] = fi2

def faults(f):
  fs = FaultSemblance()
  p2 = fs.slopes(f)
  snd = fs.semblanceNumDen(p2,f)
  fs = FaultScanner2(20,snd,FaultScanner2.Smoother.SHEAR)
  (fl,ft) = fs.scan(-20,20)
  #(fl,ft) = fs.thin((fl,ft))
  return fl
def gain(x):
  g = mul(x,x) 
  ref = RecursiveExponentialFilter(20.0)
  ref.apply1(g,g)
  y = zerofloat(n1,n2)
  div(x,sqrt(g),y)
  return y

#############################################################################
# graphics

gray = ColorMap.GRAY
jet = ColorMap.JET
pngDir = "../../../png/hv/2d/"
pngDir = None
def plot(s1,s2,x,u=None,c=None,cmap=ColorMap.GRAY,clab=None,vlabel=None,hlabel=None,
  cmin=0,cmax=0,title=None,png=None):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  #if title:
  #  sp.setTitle(title)
  #cb=sp.addColorBar(clab)
  pv = sp.addPixels(s1,s2,x)
  sp.setHLimits(s2.first,s2.last)
  sp.setVLimits(s1.first,s1.last)
  pv.setColorModel(cmap)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  if cmin<cmax:
    pv.setClips(cmin,cmax)
  if u:
    cv = sp.addContours(s1,s2,u)
    cv.setContours(50)
    cv.setColorModel(jet)
    cv.setLineWidth(2.0)
    #cv.setLineColor(Color.YELLOW)
  if c:
    k1s,k2s = c
    d1,d2= s1.delta,s2.delta
    f1,f2= s1.first,s2.first
    nc = len(k1s)
    tp1 = PointsView.Mark.HOLLOW_CIRCLE
    tp2 = PointsView.Mark.HOLLOW_SQUARE
    tp3 = PointsView.Mark.PLUS
    tps = [tp1,tp2,tp3]
    for ic in range(nc):
      np = len(k1s[ic])
      k1p = zerofloat(np)
      k2p = zerofloat(np)
      for ip in range(np):
        k1p[ip] = k1s[ic][ip]*d1
        k1p[ip] += f1
        k2p[ip] = k2s[ic][ip]*d2
        k2p[ip] += f2
      pv = PointsView(k1p,k2p)
      pv.setOrientation(PointsView.Orientation.X1DOWN_X2RIGHT)
      pv.setLineStyle(PointsView.Line.NONE)
      pv.setMarkStyle(tps[ic])
      pv.setMarkColor(Color.YELLOW)
      pv.setMarkSize(15)
      pv.setLineWidth(6)
      sp.add(pv)
  sp.setVLabel(vlabel)
  sp.setHLabel(hlabel)
  sp.setSize(500,1000)
  sp.setFontSizeForPrint(8.0,200)
  sp.plotPanel.setColorBarWidthMinimum(45)
  if pngDir and png:
    sp.paintToPng(720,2.2222,pngDir+png+".png")

#############################################################################
# utilities

def readImage(name):
  fileName = seismicDir+name+".dat"
  n1,n2 = s1.count,s2.count
  image = zerofloat(n1,n2)
  ais = ArrayInputStream(fileName,ByteOrder.BIG_ENDIAN)
  ais.readFloats(image)
  ais.close()
  return image

def writeImage(name,image):
  fileName = seismicDir+name+".dat"
  aos = ArrayOutputStream(fileName,ByteOrder.BIG_ENDIAN)
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
