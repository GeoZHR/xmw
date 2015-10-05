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

from he import *
from util import *

seismicDir = "../../../data/seis/javad/"
pngDir = "../../../png/javad/"
ffile = "event"
s1 = Sampling(1100,1.0,0.0)
s2 = Sampling(201,1.0,0.0)
n1,n2 = s1.count,s2.count
d1,d2 = s1.delta,s2.delta
f1,f2 = s1.first,s2.first

def main(args):
  goEvent()
def goEvent():
  gx = readImage(ffile)
  pmax = 10.0
  sigma1,sigma2=4.0,2.0  
  lsf = LocalSlopeFinder(sigma1,sigma2,pmax)
  p2 = zerofloat(n1,n2)
  wp = zerofloat(n1,n2)
  k1,k2=[784,699,821],[65,119,187]
  k1,k2=[784],[65]
  lsf.findSlopes(gx,p2,wp)
  he2 = HorizonExtractor2()
  he2.setWeight(0.0)
  he2.setSmoothing(6.0)
  he2.setExternalIterations(10)
  cs = [k1],[k2]
  cv = he2.curveInitialization(n2,n1-1,k1,k2)
  plot(s1,s2,gx,c=cs,ps=cv,png="initial")
  he2.curveUpdateFromSlopes(wp,p2,k1,k2,cv)
  plot(s1,s2,gx,c=cs,ps=cv,png="updat10")
  #plot(s1,s2,wp)

def flip2(f):
  n2 = len(f)
  for i2 in range(n2/2):
    fi2 = f[i2]
    f[i2] = f[n2-1-i2]
    f[n2-1-i2] = fi2

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

def plot(s1,s2,x,u=None,c=None,ps=None, 
        cmap=ColorMap.GRAY,clab=None,vlabel=None,hlabel=None,
        cmin=0,cmax=0,title=None,png=None):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  #if title:
  #  sp.setTitle(title)
  cb=sp.addColorBar(clab)
  pv = sp.addPixels(s1,s2,x)
  pv.setColorModel(cmap)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  if cmin<cmax:
    pv.setClips(cmin,cmax)
  if u:
    n1t = n1-6
    n2t = n2-6
    fxt = 3*dx+fx
    ftt = 3*dt+ft
    s1t = Sampling(n1t,dt,ftt)
    s2t = Sampling(n2t,dx,fxt)
    ut = copy(n1t,n2t,2,2,u)
    cv = sp.addContours(s1t,s2t,ut)
    cv.setClips(0.52,1.42)
    cv.setColorModel(jet)
    cv.setContours(30)
    cv.setLineWidth(3.0)
    #cv.setLineColor(Color.YELLOW)
  if ps:
    px = zerofloat(n2)
    for i2 in range(n2):
      px[i2] = i2
    pvc=sp.addPoints(ps,px)
    pvc.setLineWidth(3.0)
    pvc.setLineColor(Color.RED)
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
      pv.setMarkStyle(tps[2])
      pv.setMarkColor(Color.GREEN)
      pv.setMarkSize(12)
      pv.setLineWidth(5)
      sp.add(pv)
  sp.setVLabel(vlabel)
  sp.setHLabel(hlabel)
  sp.setSize(600+150,650)
  #sp.setFontSizeForSlide(0.5,0.9,16.0/9.0)
  sp.plotPanel.setColorBarWidthMinimum(150)
  if pngDir and png:
    sp.paintToPng(720,2.2222,pngDir+png+".png")

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
