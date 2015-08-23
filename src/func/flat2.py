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

from func import *

seismicDir = "../../../data/seis/func/2d/"
ffile = "f3d148"
s1 = Sampling(120,1.0,0.0)
s2 = Sampling(260,1.0,0.0)
n1,n2 = s1.count,s2.count
d1,d2 = s1.delta,s2.delta

def main(args):
  flatten()

def flatten():
  f = readImage(ffile)
  sigma1 = 2.0 
  sigma2 = 2.0 
  pmax = 5.0
  lof = LocalOrientFilter(sigma1,sigma1)
  u1 = zerofloat(n1,n2)
  u2 = zerofloat(n1,n2)
  wp = zerofloat(n1,n2)
  lof.applyForNormalLinear(f,u1,u2,wp)
  wp = pow(wp,8)
  fl = Rgt2()
  fl.setWeight1(0.0)
  fl.setIterations(0.01,1000)
  fl.setSmoothings(6.0,6.0)
  rgf = RecursiveGaussianFilter(1.0)
  sa = fillfloat(1.0,n1,n2)
  ut = fillfloat(0.0,n1,n2)
  for k in range(50):
    fl.setInitial(ut)
    v1 = mul(u1,sa)
    v2 = mul(u2,sa)
    fl.updateRgtFromNormals(v1,v2,wp,ut)
    sa = fl.applyForG1(ut)
    #rgf.apply1X(ut,sa)
  g = fl.flatten(ut,f)
  plot(s1,s2,f,title="Input",png="f")
  plot(s1,s2,g,title="Flattened",png="g")
  plot(s1,s2,ut,cmap=jet,title="RGT",cmin=5,cmax=n1-5,png="rgt")
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

#############################################################################
# graphics

gray = ColorMap.GRAY
jet = ColorMap.JET
pngDir = None
#pngDir = "./png/"
def plot(s1,s2,x,u=None,cmap=ColorMap.GRAY,clab=None,cmin=0,cmax=0,
         title=None,png=None):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  if title:
    sp.setTitle(title)
  sp.addColorBar(clab)
  sp.setSize(900,500)
  sp.plotPanel.setColorBarWidthMinimum(80)
  pv = sp.addPixels(s1,s2,x)
  pv.setColorModel(cmap)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  if cmin<cmax:
    pv.setClips(cmin,cmax)
  if u:
    cv = sp.addContours(s1,s2,u)
    cv.setContours(50)
    cv.setLineColor(Color.YELLOW)
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
