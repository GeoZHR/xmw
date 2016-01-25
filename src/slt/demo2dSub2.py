from utils import *

setupForSubset("2dSub2")
s1,s2,s3 = getSamplings()
n1,n2,n3 = s1.count,s2.count,s3.count

pngDir = None
pngDir = getPngDir() #"../../../png/slt/2d/sub1/"
gxfile = "gx" # input image (maybe after bilateral filtering)
elfile = "el"
slfile = "sl"
stfile = "st"
sffile = "sf"


def main(args):
  goSaltLike()
  goSaltSurfer()
  #goSaltSurferC()

def goSaltLike():
  gx = readImage2d(gxfile)
  u1 = zerofloat(n1,n2)
  u2 = zerofloat(n1,n2)
  el = zerofloat(n1,n2)
  lof = LocalOrientFilterP(4,1)
  ets = lof.applyForTensors(gx)
  ets.setEigenvalues(0.01,1.0)
  ss = SaltScanner()
  el = ss.applyForLinear(200,ets,gx)
  lof = LocalOrientFilterP(4,4)
  lof.applyForNormal(gx,u1,u2)
  sl = ss.saltLikelihood(15,el,u1,u2)
  ssf = SaltSurfer()
  st1,st2 = ssf.thin(0.3,sl,u1,u2)
  writeImage(elfile,el)
  writeImage(slfile,sl)
  writeImage(stfile,st1)
  plot2(gx,s1,s2,label="Amplitude",png="gx")
  plot2(gx,s1,s2,el,cmin=0.1,cmax=1,cmap=jetRamp(1.0),
        label="Linearity",png="el")
  plot2(gx,s1,s2,sl,cmin=0.2,cmax=1,cmap=jetRamp(1.0),
        label="Salt likelihood",png="sl")
  plot2(gx,s1,s2,st2,cmin=0.1,cmax=1,cmap=jetRamp(1.0),
        label="Thinned salt likelihood", png="st")

def goSaltSurfer():
  gx = readImage2d(gxfile)
  el = readImage2d(elfile)
  sl = readImage2d(slfile)
  st = readImage2d(stfile)
  u1 = zerofloat(n1,n2)
  u2 = zerofloat(n1,n2)
  g1 = zerofloat(n1,n2)
  g2 = zerofloat(n1,n2)
  lof = LocalOrientFilterP(8,2)
  lof.applyForNormal(gx,u1,u2)
  ss = SaltSurfer()
  mul(u1,sl,g1)
  mul(u2,sl,g2)
  mk = pow(st,7) 
  sps = ScreenPoissonSurfer()
  sps.setSmoothings(20,20,20)
  sf = sps.saltIndicator(mk,g1,g2)
  writeImage(sffile,sf)
  print min(sf)
  print max(sf)
  plot2(gx,s1,s2,st,cmin=0.1,cmax=1,cmap=jetRamp(1.0),
        label="Thinned salt likelihood")
  plot2(gx,s1,s2,sf,u=sf,cmin=-18,cmax=18,cmap=bwrRamp(1.0),
        label="Salt indicator",png="sf")

def goSaltSurferC():
  gx = readImage2d(gxfile)
  el = readImage2d(elfile)
  sl = readImage2d(slfile)
  st = readImage2d(stfile)
  u1 = zerofloat(n1,n2)
  u2 = zerofloat(n1,n2)
  g1 = zerofloat(n1,n2)
  g2 = zerofloat(n1,n2)
  lof = LocalOrientFilterP(2,1)
  lof.applyForNormal(el,u1,u2)
  ss = SaltSurfer()
  mul(u1,sl,g1)
  mul(u2,sl,g2)
  spc = ScreenPoissonSurferC()
  spc.setSmoothings(20,20,20)
  k1 = [116,106, 90, 97, 78, 72, 80, 86, 97, 74]
  k2 = [139,118,162,173,186,198,213,221,235,276]
  sf = spc.saltIndicator(k1,k2,st,g1,g2)
  writeImage(sffile,sf)
  print min(sf)
  print max(sf)
  plot2(gx,s1,s2,st,cmin=0.1,cmax=1,cmap=jetRamp(1.0),
        label="Thinned salt likelihood")
  plot2(gx,s1,s2,sf,u=sf,k1=k1,k2=k2,cmin=-15,cmax=15,
        cmap=bwrRamp(1.0),label="Salt indicator",png="sfc")


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
  panel.setHInterval(2.0)
  panel.setVInterval(0.2)
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
    cv.setLineColor(Color.CYAN)
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
