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
_dataDir = "../../../data/seis/mhe/2d/"
_pngDir = "../../../png/mhe/2d/"


def main(args):
  #goCrfPick()
  #goF3dPick()
  #goCurtPick()
  #goTpdPick()
  goSlopeDisp()

def goSlopeDisp():
  ffile = "gx596Sub"
  global s1,s2
  global n1,n2
  global seismicDir,pngDir
  pngDir = _pngDir+"crf/"
  seismicDir = _dataDir+"crf/"
  s1 = Sampling(130,1.0,0.0)
  s2 = Sampling(400,1.0,0.0)
  n1,n2 = s1.count,s2.count
  f = readImage(ffile)
  f = gain(f)
  v1 = zerofloat(n1,n2)
  v2 = zerofloat(n1,n2)
  lof = LocalOrientFilter(4,2)
  lof.apply(f,None,None,None,v1,v2,None,None,None)
  plot(s1,s2,f,v1=v1,v2=v2,w1=400,w2=650,png="crfSlope")
def goCrfPick():
  ffile = "gx596Sub"
  global s1,s2
  global n1,n2
  global seismicDir,pngDir
  pngDir = _pngDir+"crf/"
  seismicDir = _dataDir+"crf/"
  s1 = Sampling(130,1.0,0.0)
  s2 = Sampling(400,1.0,0.0)
  n1,n2 = s1.count,s2.count
  f = readImage(ffile)
  #f = copy(200,400,170,900,f)
  #writeImage("gx596Sub",f)
  f = gain(f)
  fs = zerofloat(n1,n2)
  p2 = zerofloat(n1,n2)
  wp = zerofloat(n1,n2)
  sigma1,sigma2=4.0,2.0
  pmax = 5
  lsf = LocalSlopeFinder(sigma1,sigma2,pmax)
  lsf.findSlopes(f,p2,wp) # estimate slopes and linearity
  wp = pow(wp,8)
  he2 = MultigridHorizon2()
  he2.setSmoothing(12)
  he2.setWeight(0.001)
  he2.setExternalIterations(20)
  # pick a horizon with 2 control points
  k1,k2=[53],[160]  
  k1,k2=[62],[180]  
  cv1 = he2.curveInitialization(n2,n1-1,k1,k2)
  cv2 = copy(cv1)
  cv1 = he2.curveUpdateFromSlopes(10,35,60,f,wp,p2,k1,k2,cv1);
  cv2 = he2.curveUpdateFromSlopes(wp,p2,k1,k2,cv2);
  #plot(s1,s2,f,w1=400,w2=650,png="f3dSub1")
  cv1[5]= cv1[19]
  colors = [Color.BLUE,Color.CYAN,Color.GREEN,
              Color.ORANGE,Color.MAGENTA,Color.RED]
  for ic in range(6):
    fs = he2.getCorrelationTraces(40,f,cv1[ic])
    sv = Sampling(81,1,-40)
    plot1s(sv,fs,hcolor=colors[ic],png="iter"+str(ic))
  plot(s1,s2,f,w1=400,w2=650,png="crfSub1")
  plot(s1,s2,f,hs=cv1,k2=k2,k1=k1,w1=400,w2=650,png="crfSub1m")
  plot(s1,s2,f,hs=cv2,k2=k2,k1=k1,w1=400,w2=650,png="crfSub1s")

def goF3dPick():
  ffile = "f3d178Sub2"
  global s1,s2
  global n1,n2
  global seismicDir,pngDir
  pngDir = _pngDir+"f3d/"
  seismicDir = _dataDir+"f3d/"
  s1 = Sampling(242,1.0,0.0)
  s2 = Sampling(951,1.0,0.0)
  n1,n2 = s1.count,s2.count
  f = readImage(ffile)
  f = gain(f)
  fs = zerofloat(n1,n2)
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
  k11,k12=[189],[806]  
  k21,k22=[170],[880]  
  k31,k32=[128,96],[442,805]  
  k41,k42=[42 ],[330]  
  k51,k52=[65],[756]  
  k51,k52=[59],[753]  
  k1s = [k11,k21,k31,k41,k51]
  k2s = [k12,k22,k32,k42,k52]
  cv1s,cv2s,cv3s = [],[],[]
  p1s,p2s = [],[]
  np = len(k1s)
  for ip in range(np):
    k1 = k1s[ip]
    k2 = k2s[ip]
    cv1i = he2.curveInitialization(n2,n1-1,k1,k2)
    cv2i = copy(cv1i)
    cv1i = he2.curveUpdateFromSlopes(18,35,60,f,wp,p2,k1,k2,cv1i);
    cv2i = he2.curveUpdateFromSlopes(wp,p2,k1,k2,cv2i);
    cv3i = he2.dipPredict1(k1[0],k2[0],p2)
    cv1s.append(cv1i)
    cv2s.append(cv2i)
    cv3s.append(cv3i)
    for ik in range(len(k1)):
      p1s.append(k1[ik])
      p2s.append(k2[ik])
  plot(s1,s2,f,w1=400,w2=650,png="f3d")
  plot(s1,s2,f,h=cv1s,k2=p2s,k1=p1s,w1=400,w2=650,png="f3dm")
  plot(s1,s2,f,h=cv2s,k2=p2s,k1=p1s,w1=400,w2=650,png="f3ds")
  p1s[3] = -100
  v1 = zerofloat(n1,n2)
  v2 = zerofloat(n1,n2)
  lof = LocalOrientFilter(8,2)
  lof.apply(f,None,None,None,v1,v2,None,None,None)
  plot(s1,s2,f,h=cv3s,k2=p2s,k1=p1s,w1=400,w2=650,png="f3dp")
  plot(s1,s2,f,v1=v1,v2=v2,w1=400,w2=650,png="f3dSlope")

def goCurtPick():
  ffile = "gxsb1"
  global s1,s2
  global n1,n2
  global seismicDir,pngDir
  pngDir = _pngDir+"curt/"
  seismicDir = _dataDir+"curt/"
  s1 = Sampling(210,1.0,0.0)
  s2 = Sampling(2900,1.0,0.0)
  n1,n2 = s1.count,s2.count
  f = readImage(ffile)
  f = gain(f)
  fs = zerofloat(n1,n2)
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
  k11,k12=[39],[500]  
  k21,k22=[53],[500]  
  k31,k32=[76],[500]  
  k41,k42=[105],[500]  
  k51,k52=[130],[488]  
  k61,k62=[148],[460]  
  k71,k72=[188],[500]  
  k1s = [k11,k21,k31,k41,k51,k61,k71]
  k2s = [k12,k22,k32,k42,k52,k62,k72]
  cv1s,cv2s,cv3s = [],[],[]
  p1s,p2s = [],[]
  np = len(k1s)
  for ip in range(np):
    k1 = k1s[ip]
    k2 = k2s[ip]
    cv1i = he2.curveInitialization(n2,n1-1,k1,k2)
    cv2i = copy(cv1i)
    cv1i = he2.curveUpdateFromSlopes(10,40,200,f,wp,p2,k1,k2,cv1i);
    cv2i = he2.curveUpdateFromSlopes(wp,p2,k1,k2,cv2i);
    cv3i = he2.dipPredict1(k1[0],k2[0],p2)
    cv1s.append(cv1i)
    cv2s.append(cv2i)
    cv3s.append(cv3i)
    for ik in range(len(k1)):
      p1s.append(k1[ik])
      p2s.append(k2[ik])
  plot(s1,s2,f,w1=400,w2=1050,png="curt")
  plot(s1,s2,f,h=cv1s,k2=p2s,k1=p1s,w1=400,w2=1050,png="crutm")
  plot(s1,s2,f,h=cv2s,k2=p2s,k1=p1s,w1=400,w2=1050,png="curts")
  #p1s[3] = -100
  plot(s1,s2,f,h=cv3s,k2=p2s,k1=p1s,w1=400,w2=1050,png="curtp")


def goTpdPick():
  ffile = "tp73"
  global s1,s2
  global n1,n2
  global seismicDir,pngDir
  pngDir = _pngDir+"tpd/"
  seismicDir = _dataDir+"tpd/"
  s1 = Sampling(240,1.0,0.0)
  s2 = Sampling(357,1.0,0.0)
  n1,n2 = s1.count,s2.count
  f = readImage(ffile)
  p2 = zerofloat(n1,n2)
  wp = zerofloat(n1,n2)
  sigma1,sigma2=8.0,2.0
  pmax = 5
  lsf = LocalSlopeFinder(sigma1,sigma2,pmax)
  lsf.findSlopes(f,p2,wp) # estimate slopes and linearity
  wp = pow(wp,8)
  he2 = MultigridHorizon2()
  he2.setSmoothing(6)
  he2.setWeight(0.001)
  he2.setExternalIterations(20)
  # pick a horizon with 2 control points
  #us = gcf.findSlopes(pc,f)
  cv1s = []
  cv2s = []
  cv3s = []
  k1,k2=[],[]
  for i1 in range(10,230,15):
    k1.append(i1)
    k2.append(70)
  k1[3] = k1[3]-5
  k1[4] = k1[4]-2
  k1[5] = k1[5]-2
  k1[6] = k1[6]-2
  k1[8] = k1[8]-2
  k1[9] = k1[9]-3
  k1[10] = k1[10]-5
  k1[11] = k1[11]-4
  k1[12] = k1[12]+2
  for ip in range(len(k1)):
    k1i = [k1[ip]]
    k2i = [k2[ip]]
    cv1 = he2.curveInitialization(n2,n1-1,k1i,k2i)
    cv2 = copy(cv1)
    cv1 = he2.curveUpdateFromSlopes(10,30,60,f,wp,p2,k1i,k2i,cv1);
    cv2 = he2.curveUpdateFromSlopes(wp,p2,k1i,k2i,cv2);
    cv3 = he2.dipPredict1(k1[ip],k2[ip],p2)
    cv1s.append(cv1)
    cv2s.append(cv2)
    cv3s.append(cv3)
  plot(s1,s2,f,w1=500,w2=320,png="seis")
  plot(s1,s2,f,h=cv1s,k2=k2,k1=k1,w1=500,w2=320,png="tpdm")
  plot(s1,s2,f,h=cv2s,k2=k2,k1=k1,w1=500,w2=320,png="tpds")
  plot(s1,s2,f,h=cv3s,k2=k2,k1=k1,w1=500,w2=320,png="tpdp")

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

def plot1s(s1,ys,vmin=None,vmax=None,color=Color.BLACK,hcolor=Color.CYAN,
  hlabel="Inline (trace number)",vlabel="Time (sample)",png=None):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sf = 1.0
  yf = sf
  sp.setVLimits(-40,40)
  if vmin and vmax:
    sp.setVLimits(vmin,vmax)
  sp.setHLimits(0,n2-1)
  sp.setHInterval(20)
  sp.setVInterval(20)
  y = zerofloat(n2)
  x = rampfloat(0,1,n2)
  pv2 = sp.addPoints(y,x)
  pv2.setLineColor(hcolor)
  pv2.setLineWidth(3.0)
  for il,y in enumerate(ys):
    ya = sum(y)/len(y)
    y = sub(y,ya)
    y = mul(y,3)
    y = add(y,yf)
    pv1 = sp.addPoints(s1,y)
    pv1.setLineColor(color)
    if(il==9):
      pv1.setLineColor(Color.RED)
    pv1.setLineWidth(2.0)
    yf = yf+20
  sp.setSize(650,246)
  sp.setHLabel(hlabel)
  sp.setVLabel(vlabel)
  sp.setFontSize(14) #for print
  #sp.setFontSize(30) #for slides
  if png and pngDir:
    sp.paintToPng(300,7.0,pngDir+png+".png")

def plot(s1,s2,x,h=None,hs=None,k2=None,k1=None,v1=None,v2=None,
         w1=1000,w2=500,cmap=ColorMap.GRAY,cmin=0,cmax=0,title=None,png=None):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  #if title:
  #  sp.setTitle(title)
  #cb=sp.addColorBar(clab)
  pv = sp.addPixels(s1,s2,x)
  sp.setHLimits(0,n2-1)
  sp.setVLimits(0,n1-1)
  sp.setHLabel("Inline (trace number)")
  sp.setVLabel("Time (sample)")
  pv.setColorModel(cmap)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  pv.setClips(-2,2)
  if (v1 and v2):
    x1 = zerofloat(2)
    x2 = zerofloat(2)
    dx1,dx2 = 5, 12
    scale = 4
    d1,d2=1,1
    f1,f2=0,0
    for i2 in range(0,n2,dx2):
      for i1 in range(0,n1,dx1):
        x2[0] = (i2-v2[i2][i1]*scale)*d2+f2
        x2[1] = (i2+v2[i2][i1]*scale)*d2+f2
        x1[0] = (i1-v1[i2][i1]*scale)*d1+f1
        x1[1] = (i1+v1[i2][i1]*scale)*d1+f1
        pvu = sp.addPoints(x1,x2)
        pvu.setLineWidth(4)
        pvu.setLineColor(Color.CYAN)

  if h:
    for hi in h:
      x = rampfloat(0,1,n2)
      pv = sp.addPoints(hi,x)
      pv.setLineColor(Color.YELLOW)
      pv.setLineWidth(3.0)
      sp.add(pv)

  if hs:
    colors = [Color.BLUE,Color.CYAN,Color.GREEN,
              Color.ORANGE,Color.MAGENTA,Color.RED]
    for k in range(6):
      x = rampfloat(0,1,n2)
      pv = sp.addPoints(hs[k],x)
      pv.setLineColor(colors[k])
      pv.setLineWidth(4.0)
      if(k<5):
        pv.setLineStyle(PointsView.Line.SOLID)
      sp.add(pv)
      k = k+1
  if k2 and k1:
    np = len(k2)
    for ip in range(np):
      pv = PointsView([k1[ip]],[k2[ip]])
      pv.setOrientation(PointsView.Orientation.X1DOWN_X2RIGHT)
      pv.setLineStyle(PointsView.Line.NONE)
      pv.setMarkStyle(PointsView.Mark.HOLLOW_CIRCLE)
      pv.setMarkColor(Color.RED)
      pv.setMarkSize(12)
      pv.setLineWidth(6)
      pv.setMarkSize(8)
      pv.setLineWidth(4)
      sp.add(pv)
  sp.setSize(w2,w1)
  sp.setFontSize(14)
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
