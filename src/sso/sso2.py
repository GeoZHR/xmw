import sys

from java.awt import *
from java.io import *
from java.nio import *
from java.lang import *
from javax.swing import *
from java.util import *


from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.io import *
from edu.mines.jtk.interp import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.util import *
from edu.mines.jtk.util.ArrayMath import *

from he import *
from sso import *

pngDir = "../../../png/sso/2d/poseidon/"
pngDir = None

seismicDir = "../../../data/seis/sso/2d/mexico/"
seismicDir = "../../../data/seis/sso/2d/poseidon/"
fxfile = "fx"
fxfile = "gx"
p2lfile = "p2l"
ellfile = "ell"
p2sfile = "p2s"
elsfile = "els"
f1,f2 = 0,0
d1,d2 = 1,1
n1,n2 = 301,920 # mexico fx
n1,n2 = 120,923 # poseidon fx
n1,n2 = 100,923 # poseidon gx
n1,n2 = 501,501 # fake
d1,d2 = 1.0,1.0
f1,f2 = 0.0,0.0
s1 = Sampling(n1,d1,f1)
s2 = Sampling(n2,d2,f2)
def main(args):
  #goLof()
  #goLoe()
  goFakeLoe()
  #goHorizonLof()
  #goHorizonLoe()

def goResample():
  fx = readImage(fxfile)
  fs = zerofloat(n1,600)
  c2 = Sampling(600,1.5,0)
  si = SincInterpolator()
  for i2 in range(600):
    x2 = i2*1.5
    for i1 in range(n1):
      fs[i2][i1]=si.interpolate(s1,s2,fx,i1,x2)
  writeImage("fs",fs)
def goLof():
  fx = readImage(fxfile)
  fx = gain(fx)
  addNoise(0.2,fx)
  sig1,sig2=8,2
  el = zerofloat(n1,n2)
  u1 = zerofloat(n1,n2)
  u2 = zerofloat(n1,n2)
  lof = LocalOrientFilter(sig1,sig2)
  lof.applyForNormalLinear(fx,u1,u2,el)
  ets = lof.applyForTensors(fx)
  ets.setEigenvalues(0.0001,1.0)
  gx = zerofloat(n1,n2)
  lsf = LocalSmoothingFilter()
  lsf.apply(ets,30,fx,gx)

  p2 = mul(-1,div(u2,u1))
  writeImage(p2lfile,p2)
  writeImage(ellfile,el)
  plot(s1,s2,sub(fx,gx),cmin=-0.5,cmax=0.5,cint=0.2,clab="Amplitude",png="dif1")
  plot(s1,s2,fx,cmin=-1,cmax=1,cint=1.0,clab="Amplitude",png="seis")
  plot(s1,s2,gx,cmin=-1,cmax=1,cint=1.0,clab="Amplitude",png="smooth1")
  fs = copy(54,n2,33,0,fx)
  ps = copy(54,n2,33,0,p2)
  c1 = Sampling(54,1,33)
  plot(c1,s2,fs,ps,cmin=-0.6,cmax=0.6,cmap=jetFill(0.5),cint=0.2,
       clab="Crosline slope (samples/trace)",png="slope1")


def goFakeLoe():
  fx,px=FakeData.seismicAndSlopes2d2014A(0.0)
  sig1,sig2=23,1
  el = zerofloat(n1,n2)
  u1 = zerofloat(n1,n2)
  u2 = zerofloat(n1,n2)
  ps = zerofloat(n1,n2)
  lof = LocalOrientFilter(sig1,sig2)
  lof.applyForNormalLinear(fx,u1,u2,el)
  pl = mul(-1,div(u2,u1))
  ets = lof.applyForTensors(fx)
  loe = LocalOrientEstimator(ets,10)
  #loe.setGradientSmoothing(3)
  loe.setEigenvalues(0.0001,0.5)
  loe.applyForSlopeLinear(5,fx,ps,el)
  '''
  ets = loe.applyForTensors(fx)
  loe = LocalOrientEstimator(ets,10)
  loe.setEigenvalues(0.001,0.5)
  loe.applyForSlopeLinear(5,fx,ps,el)
  ets = loe.applyForTensors(fx)
  loe = LocalOrientEstimator(ets,10)
  loe.setEigenvalues(0.001,0.5)
  loe.applyForSlopeLinear(5,fx,ps,el)
  '''
  '''
  plot(s1,s2,fx,px,cmin=-4.0,cmax=4.0,cmap=jetFill(0.5),cint=1.0,
       clab="True slopes (samples/trace)",png="slope2")
  '''
  plot(s1,s2,fx,abs(sub(pl,px)),cmap=jetFill(1.0),cmin=0.0,cmax=0.5,cint=0.1,
       clab="Slopes (samples/trace)",png="slope2")
  plot(s1,s2,fx,abs(sub(ps,px)),cmap=jetFill(1.0),cmin=0.0,cmax=0.5,cint=0.1,
       clab="Improved slopes (samples/trace)",png="slope2")

def goLoe():
  fx = readImage(fxfile)
  fx = gain(fx)
  addNoise(0.2,fx)
  sig1,sig2=8,2
  el = zerofloat(n1,n2)
  u1 = zerofloat(n1,n2)
  u2 = zerofloat(n1,n2)
  #lof = LocalOrientFilterP(sig1,sig2)
  lof = LocalOrientFilter(sig1,sig2)
  ets = lof.applyForTensors(fx)
  loe = LocalOrientEstimator(ets,10)
  loe.setEigenvalues(1,0.05)
  loe.applyForNormalLinear(fx,u1,u2,el)
  ets = loe.applyForTensors(fx)
  ets.setEigenvalues(0.0001,1.0)
  gx = zerofloat(n1,n2)
  lsf = LocalSmoothingFilter()
  lsf.apply(ets,30,fx,gx)

  p2 = mul(-1,div(u2,u1))
  writeImage(p2sfile,p2)
  writeImage(elsfile,el)
  plot(s1,s2,sub(fx,gx),cmin=-0.5,cmax=0.5,cint=0.2,clab="Amplitude",png="dif2")
  plot(s1,s2,gx,cmin=-1,cmax=1,cint=1.0,clab="Amplitude",png="smooth2")
  fs = copy(54,n2,33,0,fx)
  ps = copy(54,n2,33,0,p2)
  c1 = Sampling(54,1,33)
  plot(c1,s2,fs,ps,cmin=-0.6,cmax=0.6,cmap=jetFill(0.5),cint=0.2,
       clab="Crossline slope (samples/trace)",png="slope2")

def goHorizonLof():
  fx = readImage(fxfile)
  fx = gain(fx)
  addNoise(0.2,fx)
  p2 = readImage(p2lfile)
  el = readImage(ellfile)
  k2 = 200
  k1 = rampfloat(28,1,30)
  k2 = 430
  k1 = rampfloat(40,1,32)
  k1 = rampfloat(46,1,20)
  k2 = rampfloat(410,5,20)

  hv = HorizonVolume2()
  hc = hv.applyForHorizonVolume(k1,k2,el,p2)
  fs = copy(54,n2,33,0,fx)
  c1 = Sampling(54,1,33)
  plot(c1,s2,fs,cmin=-1.5,cmax=1.5,clab="Amplitude",png="seisSub")
  plot(c1,s2,fs,hv=hc,cmin=-1.5,cmax=1.5,clab="Amplitude",png="horizon1")


def goHorizonLoe():
  fx = readImage(fxfile)
  fx = gain(fx)
  addNoise(0.2,fx)
  p2 = readImage(p2sfile)
  el = readImage(elsfile)
  k2 = 200
  k1 = rampfloat(28,1,30)
  k2 = 430
  k1 = rampfloat(40,1,32)
  k1 = rampfloat(46,1,20)
  k2 = rampfloat(410,5,20)
  hv = HorizonVolume2()
  hc = hv.applyForHorizonVolume(k1,k2,el,p2)
  fs = copy(54,n2,33,0,fx)
  c1 = Sampling(54,1,33)
  plot(c1,s2,fs,hv=hc,cmin=-1.5,cmax=1.5,clab="Amplitude",png="horizon2")
 

def normalize(ss):
  sub(ss,min(ss),ss)
  div(ss,max(ss),ss)
  
def gain(x):
  g = mul(x,x) 
  ref = RecursiveExponentialFilter(100.0)
  ref.apply1(g,g)
  y = zerofloat(n1,n2)
  div(x,sqrt(g),y)
  return y

def addNoise(nrms, fx):
  r = Random(1);
  gx = mul(2.0,sub(randfloat(r,n1,n2),0.5));
  rgf = RecursiveGaussianFilter(2.0);
  rgf.apply2X(gx,gx)
  gx = mul(gx,nrms*rms(fx)/rms(gx));
  return add(fx,gx,fx);
def rms(fx):
  fs = sum(mul(fx,fx))
  return sqrt(fs/n1/n2)


#############################################################################
# graphics

gray = ColorMap.GRAY
jet = ColorMap.JET
def jetFill(alpha):
  return ColorMap.setAlpha(ColorMap.JET,alpha)
def jetFillExceptMin(alpha):
  a = fillfloat(alpha,256)
  a[0] = 0.0
  return ColorMap.setAlpha(ColorMap.JET,a)
def jetRamp(alpha):
  return ColorMap.setAlpha(ColorMap.JET,rampfloat(0.0,alpha/256,256))
def grayRamp(alpha):
  return ColorMap.setAlpha(ColorMap.GRAY,rampfloat(0.0,alpha/256,256))

def plot(s1,s2,f,g=None,hv=None,k1=None,k2=None,cmap=None,cmin=None,cmax=None,
        cint=None,clab=None,png=None): 
  f1 = s1.getFirst()
  f2 = s2.getFirst()
  d1 = s1.getDelta()
  d2 = s2.getDelta()
  n1 = s1.getCount()
  print cmin
  print cmax
  orientation = PlotPanel.Orientation.X1DOWN_X2RIGHT;
  panel = PlotPanel(1,1,orientation)#,PlotPanel.AxesPlacement.NONE)
  #panel.setVInterval(0.1)
  #panel.setHInterval(1.0)
  panel.setHLabel("Crossline (traces)")
  panel.setVLabel("Samples")
  pxv = panel.addPixels(0,0,s1,s2,f);
  pxv.setColorModel(ColorMap.GRAY)
  #pxv.setInterpolation(PixelsView.Interpolation.NEAREST)
  if g:
    pxv.setClips(-1,1)
  else:
    if cmin and cmax:
      pxv.setClips(cmin,cmax)
  if g:
    pv = panel.addPixels(s1,s2,g)
    pv.setInterpolation(PixelsView.Interpolation.NEAREST)
    pv.setColorModel(cmap)
    pv.setClips(cmin,cmax)
    if cmin and cmax:
      pv.setClips(cmin,cmax)
  if hv:
    nh = len(hv)
    mp = ColorMap(0,nh,ColorMap.JET)
    px = rampfloat(f2,d2,n2)
    for ih in range(0,nh,1):
      pvh=panel.addPoints(hv[ih],px)
      pvh.setLineWidth(1.0)
      pvh.setLineColor(mp.getColor(ih))
  cb = panel.addColorBar();
  if cint:
    cb.setInterval(cint)
  if clab:
    cb.setLabel(clab)
  panel.setColorBarWidthMinimum(50)
  moc = panel.getMosaic();
  frame = PlotFrame(panel);
  frame.setDefaultCloseOperation(PlotFrame.EXIT_ON_CLOSE);
  #frame.setTitle("normal vectors")
  frame.setVisible(True);
  #frame.setSize(1020,700) #for f3d
  frame.setSize(920,290) #for poseidon
  frame.setFontSize(13)
  if pngDir and png:
    frame.paintToPng(720,3.333,pngDir+png+".png")

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
