"""
Demo of salt likelihoods and salt boundaries
Author: Xinming Wu, Colorado School of Mines
Version: 2015.12.19
"""
from utils import *

setupForSubset("2dSub1")
s1,s2,s3 = getSamplings()
n1,n2,n3 = s1.count,s2.count,s3.count

# Names and descriptions of image files used below.
gxfile = "gx" # input image 
elfile = "el" # linearity image
slfile = "sl" # salt likelihoods
st1file = "st1" # thinned salt likelihoods
st2file = "st2" # thinned salt likelihoods
sffile = "sf"   # salt indicator

pngDir = None
pngDir = getPngDir() #"../../../png/slt/2d/sub1/"

# Processing begins here. When experimenting with one part of this demo, we
# can comment out earlier parts that have already written results to files.
def main(args):
  #goLinear()       # compute lilearity
  #goSaltLike()     # compute salt likelihoods
  #goSaltSurfer()   # compute salt boundaries
  goSaltSurferC()   # compute salt boundaries with constraints
  #goTest()

def goTest():
  u1 = zerofloat(n1,n2)
  g1 = zerofloat(n1,n2)
  g2 = zerofloat(n1,n2)
  u1[n2/2][n1/2] = 100
  rgf = RecursiveGaussianFilterP(1.0)
  rgf.apply00(u1,u1)
  m = 4
  sigma = sqrt(m*(m+1)/3)
  alpha = m*(m+1)/6
  print sigma
  print alpha
  rgf = RecursiveGaussianFilterP(sigma)
  rgf.apply00(u1,g1)
  lsf = LocalSmoothingFilter()
  lsf.apply(alpha,u1,g2)
  plot2(g1,s1,s2,cmin=min(g1),cmax=max(g1))
  plot2(g2,s1,s2,cmin=min(g2),cmax=max(g2))

def goLinear():
  gx = readImage2d(gxfile)
  u1 = zerofloat(n1,n2)
  u2 = zerofloat(n1,n2)
  el1 = zerofloat(n1,n2)
  el2 = zerofloat(n1,n2)
  sig1,sig2 = 60,60
  lof = LocalOrientFilterP(sig1,sig2)
  lof.applyForNormalLinear(gx,u1,u2,el1)
  ets = lof.applyForTensors(gx)
  ets.setEigenvalues(0.01,1.0)
  ss = SaltScanner()
  el2 = ss.applyForLinear(4,ets,gx)
  plot2(gx,s1,s2)
  plot2(gx,s1,s2,el1,cmin=0.1,cmax=0.9,cmap=jetRamp(1.0),
        label="Linearity",png="el1")
  plot2(gx,s1,s2,el2,cmin=0.1,cmax=0.9,cmap=jetRamp(1.0),
        label="Linearity",png="el2")

def goSaltLike():
  gx = readImage2d(gxfile)
  u1 = zerofloat(n1,n2)
  u2 = zerofloat(n1,n2)
  el = zerofloat(n1,n2)
  lof = LocalOrientFilterP(4,2)
  ets = lof.applyForTensors(gx)
  ets.setEigenvalues(0.01,1.0)
  ss = SaltScanner()
  el = ss.applyForLinear(80,ets,gx)
  lof = LocalOrientFilterP(4,4)
  lof.applyForNormal(gx,u1,u2)
  sl = ss.saltLikelihood(9,el,u1,u2)
  ssf = SaltSurfer()
  st1,st2 = ssf.thin(0.3,sl,u1,u2)
  writeImage(elfile,el)
  writeImage(slfile,sl)
  writeImage(st1file,st1)
  writeImage(st2file,st2)
  plot2(gx,s1,s2,label="Amplitude",png="gx")
  plot2(gx,s1,s2,el,cmin=0.1,cmax=0.9,cmap=jetRamp(1.0),
        label="Linearity",png="el")
  plot2(gx,s1,s2,sl,cmin=0.2,cmax=0.9,cmap=jetRamp(1.0),
        label="Salt likelihood",png="sl")
  plot2(gx,s1,s2,st2,cmin=0.2,cmax=0.9,cmap=jetRamp(1.0),
        label="Salt likelihood", png="st")
def goSaltSurfer():
  gx = readImage2d(gxfile)
  el = readImage2d(elfile)
  sl = readImage2d(slfile)
  st1 = readImage2d(st1file)
  st2 = readImage2d(st2file)
  u1 = zerofloat(n1,n2)
  u2 = zerofloat(n1,n2)
  g1 = zerofloat(n1,n2)
  g2 = zerofloat(n1,n2)
  lof = LocalOrientFilterP(8,2)
  lof.applyForNormal(el,u1,u2)
  ss = SaltSurfer()
  mul(u1,sl,g1)
  mul(u2,sl,g2)
  sps = ScreenPoissonSurfer()
  sps.setSmoothings(20,20,20)
  mk = pow(st2,7)
  sf = sps.saltIndicator(mk,g1,g2)
  writeImage(sffile,sf)
  print min(sf)
  print max(sf)
  plot2(gx,s1,s2,st2,cmin=0.1,cmax=1,cmap=jetRamp(1.0),
        label="Salt likelihood")
  plot2(gx,s1,s2,st2,u=sf,cmin=0.1,cmax=1,cmap=jetRamp(1.0),
        label="Salt likelihood",png="sb")
  plot2(gx,s1,s2,sf,cmin=-15,cmax=15,cmap=bwrRamp(1.0),
        label="Salt indicator",png="sf")
  plot2(gx,s1,s2,sf,u=sf,cmin=-15,cmax=15,cmap=bwrRamp(1.0),
        label="Salt indicator",png="sfb")

def goSaltSurferC():
  gx = readImage2d(gxfile)
  el = readImage2d(elfile)
  sl = readImage2d(slfile)
  st1 = readImage2d(st1file)
  st2 = readImage2d(st2file)
  u1 = zerofloat(n1,n2)
  u2 = zerofloat(n1,n2)
  g1 = zerofloat(n1,n2)
  g2 = zerofloat(n1,n2)
  lof = LocalOrientFilterP(8,2)
  lof.applyForNormal(el,u1,u2)
  ss = SaltSurfer()
  mul(u1,sl,g1)
  mul(u2,sl,g2)
  spc = ScreenPoissonSurferC()
  spc.setSmoothings(20,20,20)
  k1 = [116, 90, 98, 73, 80, 97, 85]#106, 90, 97, 78, 72, 80, 86, 97, 74]
  k2 = [139,162,173,190,213,235,269]#118,162,173,186,198,213,221,235,276]
  mk = pow(st1,12)
  sf = spc.saltIndicator(k1,k2,mk,g1,g2)
  writeImage(sffile,sf)
  print min(sf)
  print max(sf)
  plot2(gx,s1,s2,st2,cmin=0.1,cmax=1,cmap=jetRamp(1.0),
        label="Salt likelihood")
  plot2(gx,s1,s2,sf,cmin=-15,cmax=15,
        cmap=bwrRamp(1.0),label="Salt indicator",png="sfc")
  plot2(gx,s1,s2,st2,u=sf,k1=k1,k2=k2,cmin=0.1,cmax=1,cmap=jetRamp(1.0),
        label="Salt likelihood",png="stbc")
  plot2(gx,s1,s2,sf,u=sf,k1=k1,k2=k2,cmin=-15,cmax=15,
        cmap=bwrRamp(1.0),label="Salt indicator",png="sfbc")


def gain(x):
  g = mul(x,x) 
  ref = RecursiveExponentialFilter(10.0)
  ref.apply1(g,g)
  y = zerofloat(n1,n2)
  div(x,sqrt(g),y)
  return y

#############################################################################
# graphics

gray = ColorMap.GRAY
jet = ColorMap.JET
def jetFillExceptMin(alpha):
  a = fillfloat(alpha,256)
  a[0] = 0.0
  return ColorMap.setAlpha(ColorMap.JET,a)

def bwrFillExceptMin(alpha):
  a = fillfloat(alpha,256)
  a[0] = 0.0
  return ColorMap.setAlpha(ColorMap.BLUE_WHITE_RED,a)

def jetRamp(alpha):
  return ColorMap.setAlpha(ColorMap.JET,rampfloat(0.0,alpha/256,256))

def bwrRamp(alpha):
  return ColorMap.setAlpha(ColorMap.BLUE_WHITE_RED,rampfloat(0.0,alpha/256,256))

def grayRamp(alpha):
  return ColorMap.setAlpha(ColorMap.GRAY,rampfloat(0.0,alpha/256,256))


def plot2(f,s1,s2,g=None,u=None,k1=None,k2=None,
          cmin=None,cmax=None,cmap=None,label=None,png=None):
  n2 = len(f)
  n1 = len(f[0])
  f1,f2 = s1.getFirst(),s2.getFirst()
  d1,d2 = s1.getDelta(),s2.getDelta()
  panel = panel2Teapot()
  #panel.setHInterval(2.0)
  #panel.setVInterval(0.2)
  panel.setHLabel("Inline (km)")
  panel.setVLabel("Time (s)")
  #panel.setHInterval(100.0)
  #panel.setVInterval(100.0)
  #panel.setHLabel("Pixel")
  #panel.setVLabel("Pixel")
  if label:
    panel.addColorBar(label)
  else:
    panel.addColorBar()
  panel.setColorBarWidthMinimum(130)
  pv = panel.addPixels(s1,s2,f)
  pv.setInterpolation(PixelsView.Interpolation.LINEAR)
  pv.setColorModel(ColorMap.GRAY)
  if g:
    pv = panel.addPixels(s1,s2,g)
    pv.setInterpolation(PixelsView.Interpolation.NEAREST)
    pv.setColorModel(cmap)
  if cmin and cmax:
    pv.setClips(cmin,cmax)
  if u:
    p1,p2=2,4
    fc1 = f1+p1*d1
    fc2 = f2+p2*d2
    mc1 = n1-2*p1
    mc2 = n2-2*p2
    uc = copy(mc1,mc2,p1,p2,u)
    c1 = Sampling(mc1,d1,fc1)
    c2 = Sampling(mc2,d2,fc2)
    cv = panel.addContours(c1,c2,uc)
    cv.setContours([0])
    cv.setLineWidth(6.0)
    cv.setLineColor(Color.MAGENTA)
  if k1 and k2:
    for ip in range(len(k1)):
      k1[ip] = f1+k1[ip]*d1
      k2[ip] = f2+k2[ip]*d2
    pv = panel.addPoints(k1,k2)
    pv.setLineStyle(PointsView.Line.NONE)
    pv.setMarkStyle(PointsView.Mark.PLUS)
    pv.setMarkColor(Color.BLACK)
    pv.setMarkSize(14)
    pv.setLineWidth(6)
  frame2Teapot(panel,png)
def panel2Teapot():
  panel = PlotPanel(1,1,
    PlotPanel.Orientation.X1DOWN_X2RIGHT)#,PlotPanel.AxesPlacement.NONE)
  return panel
def frame2Teapot(panel,png=None):
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  #frame.setFontSizeForPrint(8,240)
  #frame.setSize(1240,774)
  #frame.setFontSizeForSlide(1.0,0.9)
  frame.setFontSize(36)
  frame.setSize(1080,650)
  frame.setVisible(True)
  if png and pngDir:
    frame.paintToPng(400,3.2,pngDir+"/"+png+".png")
  return frame

#############################################################################
# utilities

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
