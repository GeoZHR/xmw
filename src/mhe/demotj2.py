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

from mhe import *

seismicDir = "../../../data/seis/mhe/tj/"
ffile = "tp73"
ffile = "lulia561"
ffile = "image"
s1 = Sampling(251,1.0,0.0)
s2 = Sampling(357,1.0,0.0)
s1 = Sampling(876,1.0,0.0)
s2 = Sampling(1131,1.0,0.0)
s1 = Sampling(187,1.0,0.0)
s2 = Sampling(801,1.0,0.0)

n1,n2 = s1.count,s2.count
d1,d2 = s1.delta,s2.delta

def main(args):
  goF3dPick()
  #goTpdPick()
  #goLuliaPick()

def goF3dPick():
  f = readImageL(ffile)
  #fs = copy(242,n2,220,0,f)
  #writeImage(ffile,fs)
  f = gain(f)
  fs = zerofloat(n1,n2)
  rgf = RecursiveExponentialFilter(1)
  rgf.apply1(f,fs)
  p2 = zerofloat(n1,n2)
  wp = zerofloat(n1,n2)
  sigma1,sigma2=8.0,2.0
  pmax = 5
  lsf = LocalSlopeFinder(sigma1,sigma2,pmax)
  lsf.findSlopes(f,p2,wp) # estimate slopes and linearity
  wp = pow(wp,8)
  he2 = MultigridHorizon2()
  he2.setSmoothing(12)
  he2.setWeight(0.001)
  he2.setExternalIterations(20)
  # pick a horizon with 2 control points
  k11,k12=[24,44],[580,730]  
  cv1 = he2.curveInitialization(n2,n1-1,k11,k12)
  cv2 = copy(cv1)
  cv1 = he2.curveUpdateFromSlopes(10,35,60,f,wp,p2,k11,k12,cv1);
  cv2 = he2.curveUpdateFromSlopes(wp,p2,k11,k12,cv2);
  plot(s1,s2,f,h=[cv1],k2=k12,k1=k11,w1=500,w2=900)
  #plot(s1,s2,f,h=cv1,k2=k12,k1=k11,w1=500,w2=900)
  #plot(s1,s2,f,h=[cv2],k2=k12,k1=k11,w1=500,w2=900)

def goTpdPick():
  f = readImage(ffile)
  p2 = zerofloat(n1,n2)
  wp = zerofloat(n1,n2)
  sigma1,sigma2=8.0,2.0
  pmax = 5
  lsf = LocalSlopeFinder(sigma1,sigma2,pmax)
  lsf.findSlopes(f,p2,wp) # estimate slopes and linearity
  wp = pow(wp,8)
  he2 = MultigridHorizon2()
  he2.setSmoothing(12)
  he2.setWeight(0.005)
  he2.setExternalIterations(20)
  # pick a horizon with 2 control points
  c2 = 70
  gcf = GlobalCorrelationFinder(-35,35)
  gcf.setStrainMax(0.25)
  pc = gcf.getTraceIndexes(5,60,n2,n2,[c2],0.001)
  #us = gcf.findSlopes(pc,f)
  cv1s = []
  cv2s = []
  k1,k2=[],[]
  for i1 in range(10,240,15):
    k1.append(i1)
    k2.append(c2)
  k1[3] = k1[3]-5
  k1[4] = k1[4]-2
  k1[5] = k1[5]-2
  k1[6] = k1[6]-2
  k1[8] = k1[8]-5
  k1[9] = k1[9]-3
  k1[10] = k1[10]-5
  k1[11] = k1[11]-3
  k1[12] = k1[12]+2
  for ip in range(len(k1)):
    k1i = [k1[ip]]
    k2i = [k2[ip]]
    cv1 = he2.curveInitialization(n2,n1-1,k1i,k2i)
    cv2 = copy(cv1)
    cv1 = he2.curveUpdateFromSlopes(f,wp,p2,pc,k1i,k2i,cv1);
    cv2 = he2.curveUpdateFromSlopes(wp,p2,k1i,k2i,cv2);
    cv1s.append(cv1)
    cv2s.append(cv2)
  plot(s1,s2,f,h=cv1s,k2=k2,k1=k1)
  plot(s1,s2,f,h=cv2s,k2=k2,k1=k1)
  '''
  # pick another horizon with 5 control points
  k21,k22=[522,355,478,378,175],[800,4200,7250,8050,9800]
  cv2 = he2.curveInitialization(n2,n1-1,k21,k22)
  he2.curveUpdateFromSlopes(wp,p2,k21,k22,cv2);
  plot(s1,s2,f,h=[cv2],k2=k22,k1=k21,w1=800,w2=1400,png="horizonTwo")
  '''
def goLuliaPick():
  f = readImage(ffile)
  fs = zerofloat(n1,n2)
  rgf = RecursiveExponentialFilter(1)
  rgf.apply1(f,fs)
  p2 = zerofloat(n1,n2)
  wp = zerofloat(n1,n2)
  sigma1,sigma2=12.0,2.0
  pmax = 5
  lsf = LocalSlopeFinder(sigma1,sigma2,pmax)
  lsf.findSlopes(f,p2,wp) # estimate slopes and linearity
  wp = pow(wp,6)
  he2 = MultigridHorizon2()
  he2.setSmoothing(12)
  he2.setWeight(0.000)
  he2.setExternalIterations(20)
  # pick a horizon with 2 control points
  k11,k12=[425,280,189],[65,260,900]  
  k11,k12=[320,172],[60,785]  
  gcf = GlobalCorrelationFinder(-65,65)
  gcf.setStrainMax(0.25)
  pc = gcf.getTraceIndexes(10,60,n2,n2,k12,0.01)
  #us = gcf.findSlopes(pc,fs)
  cv1 = he2.curveInitialization(n2,n1-1,k11,k12)
  cv2 = copy(cv1)
  cv1 = he2.curveUpdateFromSlopes(f,wp,p2,pc,k11,k12,cv1);
  cv2 = he2.curveUpdateFromSlopes(wp,p2,k11,k12,cv2);
  plot(s1,s2,f,h=[cv1],k2=k12,k1=k11,w1=900,w2=900)
  plot(s1,s2,f,h=[cv2],k2=k12,k1=k11,w1=900,w2=900)

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
pngDir = None
pngDir = "./"
def plot(s1,s2,x,h=None,k2=None,k1=None,w1=1000,w2=500,cmap=ColorMap.GRAY,
    clab=None,vlabel=None,hlabel=None,
  cmin=0,cmax=0,title=None,png=None):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  #if title:
  #  sp.setTitle(title)
  #cb=sp.addColorBar(clab)
  pv = sp.addPixels(s1,s2,x)
  sp.setHLimits(0,n2-1)
  sp.setVLimits(0,n1-1)
  sp.setHLabel("Inline (sample)")
  sp.setVLabel("Time (sample)")
  pv.setColorModel(cmap)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  if cmin<cmax:
    pv.setClips(cmin,cmax)
  if h:
    for hi in h:
      x = rampfloat(0,1,n2)
      pv = sp.addPoints(hi,x)
      pv.setLineColor(Color.YELLOW)
      pv.setLineWidth(2)
      sp.add(pv)
  if k2 and k1:
    np = len(k2)
    for ip in range(np):
      pv = PointsView([k1[ip]],[k2[ip]])
      pv.setOrientation(PointsView.Orientation.X1DOWN_X2RIGHT)
      pv.setLineStyle(PointsView.Line.NONE)
      pv.setMarkStyle(PointsView.Mark.HOLLOW_CIRCLE)
      pv.setMarkColor(Color.RED)
      pv.setMarkSize(5)
      pv.setLineWidth(3)
      sp.add(pv)
  sp.setSize(w2,w1)
  sp.setFontSize(16)
  sp.plotPanel.setColorBarWidthMinimum(45)
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

def readImageL(name):
  fileName = seismicDir+name+".dat"
  n1,n2 = s1.count,s2.count
  image = zerofloat(n1,n2)
  ais = ArrayInputStream(fileName,ByteOrder.LITTLE_ENDIAN)
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
