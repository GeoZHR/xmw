import sys

from java.awt import *
from java.io import *
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

from util import *
from sso import *
from ad import *

pngDir = None
pngDir = "../../../png/ad/fed2/"

seismicDir = "../../../data/seis/ad/fed/"
fxfile = "tp73"
fxfile = "clyde200"
fxfile = "f3d75s"
f1,f2 = 0,0
d1,d2 = 1,1
n1,n2 = 251,357
n1,n2 = 400,801
n1,n2 = 200,200
n1,n2 = 421,421
n1,n2 = 222,440 # f3d75s
d1,d2 = 0.004,0.025
f1,f2 = 0.004+d1*240,0.000
s1 = Sampling(n1,d1,f1)
s2 = Sampling(n2,d2,f2)

def main(args):
  #goGaussian()
  #goDiffusivity()
  goLinearDiffusion()
  #goNonlinearDiffusion()
  #goHilbert()
  #goGaussianD()
  #goFaultSmooth()
  #goSlope()
  #goSlopeX()
def goSlope():
  gx,px=FakeData.seismicAndSlopes2d2014A(0.0)
  gx = copy(n1,n2,0,0,gx)
  px = copy(n1,n2,0,0,px)
  sig1,sig2=10,1
  lof = LocalOrientFilter(sig1,sig2)
  u1 = zerofloat(n1,n2)
  u2 = zerofloat(n1,n2)
  el = zerofloat(n1,n2)
  lof.applyForNormalLinear(gx,u1,u2,el)
  ps = mul(-1,div(u2,u1))
  pmin = min(px)
  pmax = max(px)
  clip(pmin,pmax,ps)
  plot(gx)
  plot(gx,px,cmin=-1,cmax=1,cmap=jetFill(1.0),cint=1.0)
  plot(gx,ps,cmin=-1,cmax=1,cmap=jetFill(1.0),cint=1.0)
  plot(gx,abs(sub(ps,px)),cmin=0.001,cmax=0.5,cmap=jetFill(1.0),cint=0.2)

def goSlopeX():
  gx,px=FakeData.seismicAndSlopes2d2014A(0.0)
  gx = copy(n1,n2,0,0,gx)
  px = copy(n1,n2,0,0,px)

  sig1,sig2=10,1
  e1 = zerofloat(n1,n2)
  e2 = zerofloat(n1,n2)
  u1 = zerofloat(n1,n2)
  u2 = zerofloat(n1,n2)
  lof = LocalOrientFilter(sig1,sig2)
  lof.applyForNormalLinear(gx,u1,u2,e1)
  ets = lof.applyForTensors(gx)
  sof = StratigraphicOrientFilter(sig1,sig2)
  sof.applyForNormalLinear(ets,gx,u1,u2,e2)
  ps = mul(-1,div(u2,u1))
  pmin = min(px)
  pmax = max(px)
  clip(pmin,pmax,ps)
  plot(gx)
  plot(gx,px,cmin=-1,cmax=1,cmap=jetFill(1.0),cint=1.0)
  plot(gx,ps,cmin=-1,cmax=1,cmap=jetFill(1.0),cint=1.0)
  plot(gx,abs(sub(ps,px)),cmin=0.001,cmax=0.5,cmap=jetFill(1.0),cint=0.2)


def goFaultSmooth():
  fx = readImage(fxfile)
  fx = gain(fx)
  sig1,sig2=3,1
  lof = LocalOrientFilterP(sig1,sig2)
  u1 = zerofloat(n1,n2)
  u2 = zerofloat(n1,n2)
  el = zerofloat(n1,n2)
  lof.applyForNormalLinear(fx,u1,u2,el)
  el = sub(1,el)
  los = LocalOrientFilterS(sig1,sig2)
  ets = los.applyForTensors(pow(el,0.5))
  ets.setEigenvalues(0.0001,1.0)
  sig = 32
  cycle,limit=3,0.5
  fed = FastExplicitDiffusion()
  fed.setCycles(cycle,limit)
  es = fed.apply(sig,ets,el)
  plot(fx)
  es = sub(es,min(es))
  es = div(es,max(es))
  plot(fx,el,cmin=0.1,cmax=0.2,cmap=jetRamp(1.0),cint=0.2,png="wsn")
  plot(fx,es,cmin=0.3,cmax=0.6,cmap=jetRamp(1.0),cint=0.2,png="wsn")


def goGaussianD():
  fx = readImage(fxfile)
  fx = gain(fx)
  fx = zerofloat(n1,n2)
  for i2 in range(50,150,1):
    fx[i2][i2] = 20

  '''
  for i2 in range(200,300,1):
    fx[i2][i2-150] = 10
  '''
  fx[n2/2][n1/2] = 20
  g1 = zerofloat(n1,n2)
  g2 = zerofloat(n1,n2)
  rgf = RecursiveGaussianFilterP(1.0)
  rgf.apply10(fx,g1)
  rgf.apply01(fx,g2)
  plot(g2,cmin=-1,cmax=1,cint=1.0)
  plot(g1,cmin=-1,cmax=1,cint=1.0)
  plot(fx,cmin=-1,cmax=1,cint=1.0)


def goGaussian():
  gx = zerofloat(n1,n2)
  fx = readImage(fxfile)
  fx = gain(fx)
  rgf = RecursiveGaussianFilter(16)
  rgf.apply00(fx,gx)
  plot(gx,cmin=-1,cmax=1,cint=1.0,png="gxg")

def goLinearDiffusion():
  fx = readImage(fxfile)
  fx = gain(fx)
  sig1,sig2=4,2
  lof = LocalOrientFilterP(sig1,sig2)
  ets = lof.applyForTensors(fx)
  ets.setEigenvalues(0.001,1.0)
  sig = 8
  cycle,limit=3,0.5
  fed = FastExplicitDiffusion()
  fed.setCycles(cycle,limit)
  gx = fed.apply(sig,ets,fx)
  plot(sub(fx,gx),cmin=-0.5,cmax=0.5,cint=0.2,label="Amplitude",png="fgl")
  plot(gx,cmin=-1,cmax=1,cint=1.0,label="Amplitude",png="gxl")
  plot(fx,cmin=-1,cmax=1,cint=1.0,label="Amplitude",png="fx")
def goDiffusivity():
  fx = readImage(fxfile)
  fx = gain(fx)
  sig1,sig2=4,2
  lof = LocalOrientFilterP(sig1,sig2)
  ets = lof.applyForTensors(fx)
  ets.setEigenvalues(0.001,1.0)
  sig = 8
  lbd = 0.12
  cycle,limit=3,0.5
  fed = FastExplicitDiffusion()
  fed.setCycles(cycle,limit)
  fw = zerofloat(n1,n2)
  fws = zerofloat(n1,n2)
  fed.applyForWeightsP0(lbd,ets,fx,fw)
  fed.applyForWeightsP1(lbd,ets,fx,fws)
  fws = pow(fws,4)
  plot(fx,sub(1,fw),cmin=0.1,cmax=1,cmap=jetFillExceptMin(1.0),cint=0.2,
       label="Diffusivity",png="fw")
  plot(fx,sub(1,fws),cmin=0.1,cmax=1,cmap=jetFillExceptMin(1.0),cint=0.2,
       label="Diffusivity",png="fws")

def goNonlinearDiffusion():
  fx = readImage(fxfile)
  fx = gain(fx)
  sig1,sig2=4,2
  lof = LocalOrientFilterP(sig1,sig2)
  ets = lof.applyForTensors(fx)
  sig = 8
  lbd = 0.12
  cycle,limit=3,0.5
  fed = FastExplicitDiffusion()
  fed.setCycles(cycle,limit)
  gx = fed.apply(sig,lbd,ets,fx)
  wp = zerofloat(n1,n2)
  ws = zerofloat(n1,n2)
  fed.applyForWeightsP(lbd,ets,gx,wp)
  fed.applyForWeightsP(lbd,ets,fx,ws)
  plot(sub(fx,gx),cmin=-0.5,cmax=0.5,cint=0.2,label="Amplitude",png="fgn")
  plot(fx,sub(1,ws),cmin=0.1,cmax=1,cmap=jetFillExceptMin(1.0),cint=0.2,
       label="Diffusivity",png="wsn")
  plot(gx,sub(1,wp),cmin=0.1,cmax=1,cmap=jetFillExceptMin(1.0),cint=0.2,
       label="Diffusivity",png="wpn")
  plot(gx,cmin=-1,cmax=1,cint=1.0,label="Amplitude",png="gxn")
  plot(fx,cmin=-1,cmax=1,cint=1.0,label="Amplitude",png="fx")
  plot(fx,t=ets,cmin=-1,cmax=1,cint=1.0,label="Amplitude",png="fxet")

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
  rgf.apply10(gx,gx)
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
  return ColorMap.setAlpha(ColorMap.HUE_BLUE_TO_RED,rampfloat(0.0,alpha/256,256))
def grayRamp(alpha):
  return ColorMap.setAlpha(ColorMap.GRAY,rampfloat(0.0,alpha/256,256))

def plot(f,g=None,t=None,cmap=None,cmin=None,cmax=None,cint=None,
        label=None,neareast=False,png=None): 
  orientation = PlotPanel.Orientation.X1DOWN_X2RIGHT;
  panel = PlotPanel(1,1,orientation)#,PlotPanel.AxesPlacement.NONE)
  panel.setVInterval(0.1)
  panel.setHInterval(1.0)
  panel.setHLabel("Inline (km)")
  panel.setVLabel("Time (s)")
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
    if cmin and cmax:
      pv.setClips(cmin,cmax)
  cb = panel.addColorBar();
  if cint:
    cb.setInterval(cint)
  if label:
    cb.setLabel(label)
  if t:
    tv = TensorsView(s1,s2,t)
    tv.setOrientation(TensorsView.Orientation.X1DOWN_X2RIGHT)
    tv.setLineColor(Color.YELLOW)
    tv.setLineWidth(5)
    tv.setEllipsesDisplayed(12)
    tv.setScale(0.9)
    tile = panel.getTile(0,0)
    tile.addTiledView(tv)
  panel.setColorBarWidthMinimum(100)
  moc = panel.getMosaic();
  frame = PlotFrame(panel);
  frame.setDefaultCloseOperation(PlotFrame.EXIT_ON_CLOSE);
  #frame.setTitle("normal vectors")
  frame.setVisible(True);
  frame.setSize(1050,700) #for f3d
  frame.setFontSize(30)
  if pngDir and png:
    frame.paintToPng(300,3.333,pngDir+png+".png")

def plot2(f,s1,s2,g=None,cmin=None,cmax=None,label=None,png=None,et=None):
  n2 = len(f)
  n1 = len(f[0])
  panel = panel2Teapot()
  panel = PlotPanel(1,1,PlotPanel.Orientation.X1DOWN_X2RIGHT)
  panel.setHInterval(1.0)
  panel.setVInterval(0.1)
  panel.setHLabel("Inline (km)")
  panel.setVLabel("Time (s)")
  if label:
    panel.addColorBar(label)
  else:
    panel.addColorBar()
  panel.setColorBarWidthMinimum(180)
  pv = panel.addPixels(s1,s2,f)
  pv.setInterpolation(PixelsView.Interpolation.LINEAR)
  pv.setColorModel(ColorMap.GRAY)
  #pv.setClips(min(f),max(f))
  if g:
    alpha = 0.8
  else:
    g = zerofloat(s1.count,s2.count)
    alpha = 0.0
  pv = panel.addPixels(s1,s2,g)
  #pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  #pv.setColorModel(jetFillExceptMin(1.0))
  pv.setColorModel(jetRamp(1.0))
  if cmin and cmax:
    pv.setClips(cmin,cmax)
  frame2Teapot(panel,png)
def panel2Teapot():
  panel = PlotPanel(1,1,
    PlotPanel.Orientation.X1DOWN_X2RIGHT,
    PlotPanel.AxesPlacement.NONE)
  return panel
def frame2Teapot(panel,png=None):
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  #frame.setFontSizeForPrint(8,240)
  #frame.setSize(1240,774)
  frame.setFontSizeForSlide(1.0,0.8)
  frame.setSize(880,700)
  frame.setVisible(True)
  if png and pngDir:
    frame.paintToPng(400,3.2,pngDir+"/"+png+".png")
  return frame
def makePointSets(cmap,f,x1,x2):
  sets = {}
  for i in range(len(f)):
    if f[i] in sets:
      points = sets[f[i]]
      points[0].append(f[i])
      points[1].append(x1[i])
      points[2].append(x2[i])
    else:
      points = [[f[i]],[x1[i]],[x2[i]]] # lists of f, x1, x2
      sets[f[i]] = points
  ns = len(sets)
  fs = zerofloat(1,ns)
  x1s = zerofloat(1,ns)
  x2s = zerofloat(1,ns)
  il = 0
  for points in sets:
    fl = sets[points][0]
    x1l = sets[points][1]
    x2l = sets[points][2]
    nl = len(fl)
    fs[il] = zerofloat(nl)
    x1s[il] = zerofloat(nl)
    x2s[il] = zerofloat(nl)
    copy(fl,fs[il])
    copy(x1l,x1s[il])
    copy(x2l,x2s[il])
    il += 1
  return fs,x1s,x2s

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
